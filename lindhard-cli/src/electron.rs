//! Electron runs: an input with an `[electron]` table
//! ([`lindhard::input::electron`]). Builds the elastic and inelastic
//! cross-section tables of every material the target uses, runs
//! [`lindhard::electron::transport`] with the full electron tally
//! ([`lindhard::tally::FullElectronTally`]) and writes
//! `electron_summary.json` plus CSV files for the spectra, the deposition
//! grids and the tables. Layout: `docs/cli.md`, section "Electron runs".
//!
//! Nothing here is physics: every model is the library's, chosen and
//! parameterised by the input, and every choice and data provenance is
//! written to the summary so a result can be reproduced from its own header.

use std::fmt::Write as _;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use lindhard::electron::data::{CrossSectionTable, ShellChannelTable};
use lindhard::electron::inelastic::table::mean_loss_ev;
use lindhard::electron::transport::{LayerTables, RunMetadata, Transport};
use lindhard::input::electron::{
    assemble, DataFile, ElectronInput, ResolvedElectron, ResolvedElectronMaterial,
};
use lindhard::input::ModelInfo;
use lindhard::material::MaterialSpec;
use lindhard::tally::{
    fit_psf, ChargeCells, ElectronReport, FullElectronTally, Histogram, PsfFitOptions, PsfModel,
    PsfReport,
};
use serde::Serialize;

use crate::output::{density, software, Format, Software, NM};
use crate::table_cache::{CacheFile, TableCache, TableKind, TableSource};

/// Name of the electron summary format.
pub const FORMAT_NAME: &str = "lindhard-electron-summary";
/// Bumped only on a breaking change (a key removed or changed in meaning).
pub const FORMAT_VERSION: u32 = 1;

pub const SUMMARY_FILE: &str = "electron_summary.json";
pub const SPECTRA_FILE: &str = "electron_escape_spectra.csv";
pub const CARTESIAN_FILE: &str = "electron_deposition_cartesian.csv";
pub const CYLINDRICAL_FILE: &str = "electron_deposition_cylindrical.csv";
pub const CHARGE_CARTESIAN_FILE: &str = "electron_charge_cartesian.csv";
pub const CHARGE_CYLINDRICAL_FILE: &str = "electron_charge_cylindrical.csv";
pub const PSF_PROFILE_FILE: &str = "electron_psf_profile.csv";
pub const PSF_PARAMETERS_FILE: &str = "electron_psf_parameters.csv";
pub const DEPTH_DOSE_FILE: &str = "electron_depth_dose.csv";
pub const TRANSMISSION_FILE: &str = "electron_interface_transmission.csv";
pub const EXPOSURE_FILE: &str = "electron_exposure_histogram.csv";
pub const TABLES_FILE: &str = "electron_tables.csv";
pub use crate::output::DAMAGE_EXPORT_FILE;

/// Histories per work chunk of the parallel driver. Fixed, so the summation
/// order, and therefore every output bit, does not depend on the thread
/// count; recorded in `physics.transport.chunk_size`. Electron histories
/// (with their secondaries) are long, so chunks are smaller than the ion
/// engine's.
pub const CHUNK_SIZE: u64 = 16;

/// The two tables of one material.
#[derive(Debug)]
pub struct MaterialTables {
    /// Elastic table.
    pub elastic: CrossSectionTable,
    /// Inelastic table.
    pub inelastic: CrossSectionTable,
    /// Where the elastic table came from.
    pub elastic_origin: TableOrigin,
    /// Where the inelastic table came from.
    pub inelastic_origin: TableOrigin,
}

/// Where a table came from: built by the run or read from the table cache,
/// and the cache file when a cache is in use.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct TableOrigin {
    /// Built or read.
    pub source: TableSource,
    /// The cache file (`None` without `--table-cache`).
    pub cache: Option<CacheFile>,
}

impl TableOrigin {
    fn built() -> Self {
        Self {
            source: TableSource::Built,
            cache: None,
        }
    }
}

/// One table of one material: from `cache` when given, else built.
fn table(
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
    kind: TableKind,
    cache: Option<&TableCache>,
) -> Result<(CrossSectionTable, TableOrigin)> {
    let build = || -> Result<CrossSectionTable> {
        Ok(match kind {
            TableKind::Elastic => assemble::elastic_table(r, &m.material)?,
            TableKind::Inelastic => assemble::inelastic_table(r, m)?,
        })
    };
    let field = format!("electron.materials.{}", m.name);
    let what = format!("{field}: the {} table", kind.label());
    match cache {
        None => Ok((
            build().with_context(|| format!("{what}: building"))?,
            TableOrigin::built(),
        )),
        Some(c) => {
            let t = c
                .get_or_build(r, m, kind, || build().context("building"))
                .with_context(|| what.clone())?;
            Ok((
                t.table,
                TableOrigin {
                    source: t.source,
                    cache: t.cache,
                },
            ))
        }
    }
}

/// Thread count and timings (the only thread-dependent part of the output).
#[derive(Debug, Clone, Copy, Serialize)]
pub struct ElectronRunInfo {
    pub threads: usize,
    pub table_build_s: f64,
    pub transport_s: f64,
    pub histories_per_s: f64,
}

/// Everything a finished electron run produced.
#[derive(Debug)]
pub struct ElectronSimulation {
    /// Tables per material, in [`ResolvedElectron::materials`] order.
    pub tables: Vec<MaterialTables>,
    /// The run configuration as the engine recorded it.
    pub metadata: RunMetadata,
    /// The tally's report.
    pub report: ElectronReport,
    /// The PSF profile and its fits (with `tally.psf`).
    pub psf: Option<PsfOutcome>,
    /// The per-voxel deposition (with `tally.voxel`).
    pub voxel: Option<lindhard::tally::VoxelReport>,
    /// Threads and timings.
    pub info: ElectronRunInfo,
}

/// The radial profile of `tally.psf` with the fits that succeeded, and the
/// error of each fit that did not.
#[derive(Debug)]
pub struct PsfOutcome {
    /// The profile and the successful fits, in input order.
    pub report: PsfReport,
    /// The models whose fit failed, with the error.
    pub errors: Vec<(PsfModel, String)>,
}

/// Fits each requested model to the profile. A failing fit is recorded and
/// does not stop the others or the run.
fn psf_outcome(r: &ResolvedElectron, report: &ElectronReport) -> Option<PsfOutcome> {
    let profile = report.deposition.psf.clone()?;
    let options = PsfFitOptions {
        normalization: r.psf_normalization,
        ..Default::default()
    };
    let mut fits = Vec::new();
    let mut errors = Vec::new();
    for &m in &r.psf_fits {
        match fit_psf(&profile, m, &options) {
            Ok(f) => fits.push(f),
            Err(e) => errors.push((m, e.to_string())),
        }
    }
    Some(PsfOutcome {
        report: PsfReport { profile, fits },
        errors,
    })
}

/// The pure parts of the assembly live in the library
/// ([`lindhard::input::electron::assemble`], #750) so that the library's own
/// benches, examples and tests build the tables `lindhard run` builds. These
/// two are part of the table cache key and are used by `lindhard-py` and the
/// CLI examples, so they keep their `lindhard_cli::electron` path.
pub use lindhard::input::electron::assemble::{inelastic_axis, rate_refinement};

/// The cross-section tables of every material of `r` (or read from `cache`)
/// and the transport engine built on them, with the insulator and inner-shell
/// channels the input asks for, plus the seconds the tables took. The tables
/// are built on `pool`.
pub fn build_transport(
    r: &ResolvedElectron,
    pool: &rayon::ThreadPool,
    cache: Option<&TableCache>,
) -> Result<(Transport, Vec<MaterialTables>, f64)> {
    let t0 = Instant::now();
    // The table cache and the per-material `anyhow` context are the CLI's; the
    // builds themselves are the library's.
    let built: Vec<(LayerTables, TableOrigin, TableOrigin)> = pool.install(|| {
        r.materials
            .iter()
            .map(|m| {
                let (elastic, elastic_origin) = table(r, m, TableKind::Elastic, cache)?;
                let (inelastic, inelastic_origin) = table(r, m, TableKind::Inelastic, cache)?;
                Ok((
                    LayerTables { elastic, inelastic },
                    elastic_origin,
                    inelastic_origin,
                ))
            })
            .collect::<Result<_>>()
    })?;
    let shells: Vec<Vec<ShellChannelTable>> = pool.install(|| {
        r.materials
            .iter()
            .map(|m| {
                assemble::shell_tables(r, m).with_context(|| {
                    format!(
                        "electron.materials.{}.inner_shells: building the channels",
                        m.name
                    )
                })
            })
            .collect::<Result<_>>()
    })?;
    let table_build_s = t0.elapsed().as_secs_f64();

    let mut layer_tables = Vec::with_capacity(built.len());
    let mut origins = Vec::with_capacity(built.len());
    for (t, elastic_origin, inelastic_origin) in built {
        layer_tables.push(t);
        origins.push((elastic_origin, inelastic_origin));
    }
    let assembled = assemble::MaterialTables::new(r, layer_tables, shells)?;
    let transport = assemble::transport(r, &assembled)?;
    let tables = assembled
        .tables
        .into_iter()
        .zip(origins)
        .map(|(t, (elastic_origin, inelastic_origin))| MaterialTables {
            elastic: t.elastic,
            inelastic: t.inelastic,
            elastic_origin,
            inelastic_origin,
        })
        .collect();
    Ok((transport, tables, table_build_s))
}

/// Build the tables (or read them from `cache`), run every history on
/// `threads` workers (`None`: all cores) and report. Neither the thread count
/// nor the cache changes the results: a cached table is the built one, bit
/// for bit (the cache form round-trips every `f64` exactly).
pub fn simulate_electron(
    r: &ResolvedElectron,
    threads: Option<usize>,
    cache: Option<&TableCache>,
) -> Result<ElectronSimulation> {
    if threads == Some(0) {
        bail!("threads must be at least 1");
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.unwrap_or(0))
        .build()
        .context("building the thread pool")?;

    let (transport, tables, table_build_s) = build_transport(r, &pool, cache)?;
    let mut proto = FullElectronTally::new(&transport, r.tally.clone())
        .context("setting up the electron tally")?;
    if let Some(g) = r.voxel_grid {
        proto = proto.with_voxel_grid(g);
    }

    let t1 = Instant::now();
    let run = pool
        .install(|| {
            transport.run_beam(
                r.input.run.seed,
                r.input.run.histories,
                CHUNK_SIZE,
                &r.primary,
                &r.beam,
                || proto.clone(),
            )
        })
        .context("electron transport failed")?;
    let transport_s = t1.elapsed().as_secs_f64();
    let report = run.tally.report();
    let psf = psf_outcome(r, &report);
    let voxel = run.tally.voxel_report();
    Ok(ElectronSimulation {
        tables,
        metadata: run.metadata,
        report,
        psf,
        voxel,
        info: ElectronRunInfo {
            threads: pool.current_num_threads(),
            table_build_s,
            transport_s,
            histories_per_s: r.input.run.histories as f64 / transport_s.max(f64::MIN_POSITIVE),
        },
    })
}

// ---------------------------------------------------------------------------
// electron_summary.json

#[derive(Serialize)]
struct LayerOut {
    index: usize,
    source: String,
    front_nm: f64,
    back_nm: Option<f64>,
    atom_density_per_cm3: f64,
    material: MaterialSpec,
}

#[derive(Serialize)]
struct ElfOut<'a> {
    #[serde(flatten)]
    file: &'a DataFile,
    material: &'a str,
    provenance: &'a str,
    energy_min_ev: f64,
    energy_max_ev: f64,
    points: usize,
}

#[derive(Serialize)]
struct TableOut<'a> {
    model: &'a str,
    material: &'a str,
    provenance: &'a str,
    format_version: u32,
    energy_min_ev: f64,
    energy_max_ev: f64,
    energies: usize,
    probabilities: usize,
    #[serde(flatten)]
    origin: &'a TableOrigin,
}

fn table_out<'a>(t: &'a CrossSectionTable, origin: &'a TableOrigin) -> TableOut<'a> {
    let e = t.energy_ev();
    TableOut {
        model: t.model(),
        material: t.material(),
        provenance: t.provenance(),
        format_version: t.format_version(),
        energy_min_ev: e[0],
        energy_max_ev: e[e.len() - 1],
        energies: e.len(),
        probabilities: t.probability().len(),
        origin,
    }
}

/// `physics.materials[].inner_shells`: the subshells of every element that
/// are channels, and the binding-energy table they came from.
#[derive(Serialize)]
struct InnerShellsOut<'a> {
    subshells: Vec<&'static str>,
    /// `"embedded EADL2017"` or `"file"`.
    binding_source: &'static str,
    /// The `binding` file (path, SHA-256); absent for the embedded table.
    #[serde(skip_serializing_if = "Option::is_none")]
    binding_file: Option<&'a DataFile>,
    binding_provenance: &'a str,
}

#[derive(Serialize)]
struct MaterialOut<'a> {
    name: &'a str,
    optical_elf: ElfOut<'a>,
    band: Option<&'a lindhard::electron::boundary::BandStructure>,
    phonon: Option<&'a lindhard::electron::phonon::FrohlichPhonon>,
    polaron: Option<&'a lindhard::electron::phonon::PolaronTrapping>,
    #[serde(skip_serializing_if = "Option::is_none")]
    inner_shells: Option<InnerShellsOut<'a>>,
    elastic_table: TableOut<'a>,
    inelastic_table: TableOut<'a>,
}

#[derive(Serialize)]
struct Physics<'a> {
    models: Vec<ModelInfo>,
    transport: &'a RunMetadata,
    target: Vec<LayerOut>,
    materials: Vec<MaterialOut<'a>>,
}

#[derive(Serialize)]
struct Files {
    escape_spectra: &'static str,
    tables: &'static str,
    deposition_cartesian: Option<&'static str>,
    deposition_cylindrical: Option<&'static str>,
    charge_cartesian: Option<&'static str>,
    charge_cylindrical: Option<&'static str>,
    psf_profile: Option<&'static str>,
    psf_parameters: Option<&'static str>,
    depth_dose: Option<&'static str>,
    interface_transmission: Option<&'static str>,
    exposure_histogram: Option<&'static str>,
    damage_export: Option<&'static str>,
}

#[derive(Serialize)]
struct Summary<'a> {
    format: Format,
    software: Software,
    input: &'a ElectronInput,
    physics: Physics<'a>,
    results: serde_json::Value,
    files: Files,
    run: ElectronRunInfo,
}

/// The report as summary data: the histograms and grid arrays, which go to
/// CSV files, are replaced by their totals.
fn results_json(
    report: &ElectronReport,
    psf: Option<&PsfOutcome>,
) -> serde_json::Result<serde_json::Value> {
    let mut v = serde_json::to_value(report)?;
    // The per-bin profile goes to the CSV; keep the totals and the fits. A
    // null `deposition.psf` stays, so the existing summary key is not removed.
    if let Some(d) = v["deposition"].as_object_mut() {
        if d.get("psf").is_some_and(|p| !p.is_null()) {
            d.remove("psf");
        }
    }
    if let (Some(o), Some(p)) = (v.as_object_mut(), psf) {
        o.insert("psf".into(), psf_json(p)?);
    }
    // The per-electron histograms go to the CSV; keep the moments and the
    // compound-Poisson summaries.
    if let Some(vols) = v["deposition"]["exposure"]["volumes"].as_array_mut() {
        for vol in vols {
            if let Some(h) = vol["histogram"].as_object_mut() {
                h.remove("counts");
            }
        }
    }
    for face in ["front", "back"] {
        let f = &mut v[face];
        if let Some(o) = f.as_object_mut() {
            o.remove("energy_histogram");
        }
        for class in ["slow", "fast"] {
            if let Some(o) = f[class].as_object_mut() {
                o.remove("polar_histogram");
            }
        }
    }
    // The per-bin depth profile and the spectra go to CSV files.
    if let Some(d) = v["dose"].as_object_mut() {
        if d.get("depth").is_some_and(|x| !x.is_null()) {
            d.remove("depth");
        }
        if let Some(t) = d.get_mut("transmission").and_then(|t| t.as_array_mut()) {
            for e in t {
                if let Some(o) = e.as_object_mut() {
                    o.remove("spectrum");
                }
            }
        }
    }
    for grid in ["cartesian", "cylindrical"] {
        let g = &mut v["deposition"][grid];
        if let Some(o) = g.as_object_mut() {
            let inside: f64 = o
                .remove("energy_ev")
                .and_then(|a| {
                    a.as_array()
                        .map(|a| a.iter().filter_map(|x| x.as_f64()).sum())
                })
                .unwrap_or(0.0);
            o.insert("inside_ev".into(), serde_json::json!(inside));
            // The per-cell channel maps go to the CSV too.
            o.remove("ionisation_ev");
            o.remove("lattice_ev");
        }
    }
    // The per-cell charge goes to the CSV files; keep the totals inside the
    // grid and outside it.
    for grid in ["cartesian", "cylindrical"] {
        let Some(cells) = v["charge"][grid]["cells"].as_object_mut() else {
            continue;
        };
        let (mut stopped, mut created) = (0u64, 0u64);
        for c in cells
            .remove("counts")
            .and_then(|a| a.as_array().cloned())
            .unwrap_or_default()
        {
            for k in ["stopped_primary", "stopped_secondary"] {
                stopped += c[k].as_u64().unwrap_or(0);
            }
            for k in ["created_by_primary", "created_by_secondary"] {
                created += c[k].as_u64().unwrap_or(0);
            }
        }
        cells.remove("net_sum_sq");
        cells.insert("inside_stopped".into(), serde_json::json!(stopped));
        cells.insert("inside_created".into(), serde_json::json!(created));
    }
    Ok(v)
}

/// `results.damage`: per layer, the ionisation and lattice energy of the
/// deposition (eV, summed over histories and per primary), the electron-hole
/// pairs per incident electron where the layer's material has an `E_pair`
/// (otherwise the key is absent), and the knock-on threshold report at the
/// beam energy.
///
/// The ionisation channel includes the residual energy of electrons that fell
/// below the cutoff (a stated convention, see `lindhard::tally::electron`).
/// Elastic nuclear recoil energy is not tallied by the electron engine.
/// `T_max` uses the beam energy as given in the input (for a step-barrier run
/// that is the vacuum energy). An element without `E_d` is reported as
/// `"not set"`.
fn damage_json(r: &ResolvedElectron, report: &ElectronReport) -> serde_json::Value {
    use serde_json::json;
    let d = &report.deposition;
    let h = report.histories.max(1) as f64;
    let beam = r.input.electron.beam.energy_ev;
    let mut layers = Vec::new();
    let mut pairs_total = Some(0.0);
    for (i, &mi) in r.layer_material.iter().enumerate() {
        let m = &r.materials[mi];
        let ion = d.per_layer_ionisation_ev.get(i).copied().unwrap_or(0.0);
        let lat = d.per_layer_lattice_ev.get(i).copied().unwrap_or(0.0);
        let mut o = json!({
            "layer": i,
            "material": m.name,
            "ionisation_ev": ion,
            "lattice_ev": lat,
            "ionisation_ev_per_primary": ion / h,
            "lattice_ev_per_primary": lat / h,
        });
        match &m.e_pair {
            Some(p) => {
                o["e_pair_ev"] = json!(p.ev);
                o["e_pair_source"] = json!(p.source);
                o["pairs_per_primary"] = json!(ion / h / p.ev);
                if let Some(t) = pairs_total.as_mut() {
                    *t += ion / h / p.ev;
                }
            }
            None => {
                // A layer with ionisation energy and no E_pair has no
                // defined pair count, so there is no total either.
                if ion > 0.0 {
                    pairs_total = None;
                }
            }
        }
        let k = lindhard::electron::knock_on::layer_report(i, r.stack.layers()[i].material(), beam);
        o["knock_on"] = json!(k
            .elements
            .iter()
            .map(|e| json!({
                "element": e.symbol,
                "z": e.z,
                "t_max_ev": e.t_max_ev,
                "e_d_ev": e.e_d_ev.map_or(json!("not set"), |v| json!(v)),
                "can_displace": e.can_displace.map_or(json!("not set"), |v| json!(v)),
                "threshold_beam_ev": e.threshold_beam_ev.map_or(json!("not set"), |v| json!(v)),
            }))
            .collect::<Vec<_>>());
        layers.push(o);
    }
    let mut out = json!({
        "beam_energy_ev": beam,
        "residual_channel": "ionisation",
        "layers": layers,
    });
    if let Some(t) = pairs_total {
        if r.layer_material
            .iter()
            .any(|&mi| r.materials[mi].e_pair.is_some())
        {
            out["pairs_per_primary"] = json!(t);
        }
    }
    out
}

/// `results.psf`: the profile totals and, per fit, its parameters and quality.
/// Per-bin arrays and residuals are in the CSV files.
fn psf_json(p: &PsfOutcome) -> serde_json::Result<serde_json::Value> {
    let pr = &p.report.profile;
    let mut fits = Vec::new();
    for f in &p.report.fits {
        fits.push(serde_json::json!({
            "model": f.model,
            "source": f.source,
            "normalization": f.normalization,
            "parameter_names": f.parameter_names,
            "values": f.values,
            "std_errors": f.std_errors,
            "reduced_chi2": f.reduced_chi2,
            "dof": f.dof,
            "converged": f.converged,
        }));
    }
    let mut errors = Vec::new();
    for (m, e) in &p.errors {
        errors.push(serde_json::json!({ "model": m, "error": e }));
    }
    Ok(serde_json::json!({
        "histories": pr.histories,
        "depth_lo_m": pr.depth_lo_m,
        "depth_hi_m": pr.depth_hi_m,
        "bins": pr.len(),
        "total_ev": pr.total_ev,
        "total_std_err_ev": pr.total_std_err_ev,
        "beyond_ev": pr.beyond_ev,
        "beyond_std_err_ev": pr.beyond_std_err_ev,
        "fits": fits,
        "fit_errors": errors,
    }))
}

/// `electron_psf_profile.csv` (with `tally.psf`).
pub fn psf_profile_csv(sim: &ElectronSimulation) -> Option<String> {
    sim.psf.as_ref().map(|p| p.report.profile_csv())
}

/// `electron_exposure_histogram.csv` (with `tally.exposure`).
pub fn exposure_histogram_csv(sim: &ElectronSimulation) -> Option<String> {
    sim.report
        .deposition
        .exposure
        .as_ref()
        .map(|e| e.histogram_csv())
}

/// `electron_psf_parameters.csv` (with `tally.psf`).
pub fn psf_parameters_csv(sim: &ElectronSimulation) -> Option<String> {
    sim.psf.as_ref().map(|p| p.report.parameters_csv())
}

/// The `electron_summary.json` text (pretty-printed, trailing newline).
pub fn summary_json(r: &ResolvedElectron, sim: &ElectronSimulation) -> Result<String> {
    let target = r
        .stack
        .layers()
        .iter()
        .zip(&r.layers)
        .enumerate()
        .map(|(index, (g, l))| LayerOut {
            index,
            source: l.source.clone(),
            front_nm: g.front_m() / NM,
            back_nm: g.back_m().is_finite().then(|| g.back_m() / NM),
            atom_density_per_cm3: g.material().atom_number_density() * 1e-6,
            material: MaterialSpec::from(g.material().clone()),
        })
        .collect();
    let materials = r
        .materials
        .iter()
        .zip(&sim.tables)
        .map(|(m, t)| {
            let (lo, hi) = m.optical_elf.energy_range_ev();
            MaterialOut {
                name: &m.name,
                optical_elf: ElfOut {
                    file: &m.optical_elf_file,
                    material: m.optical_elf.material(),
                    provenance: m.optical_elf.provenance(),
                    energy_min_ev: lo,
                    energy_max_ev: hi,
                    points: m.optical_elf.energy_ev().len(),
                },
                band: m.band.as_ref(),
                phonon: m.channels.phonon.as_ref(),
                polaron: m.channels.polaron.as_ref(),
                inner_shells: m.inner_shells.as_ref().map(|i| InnerShellsOut {
                    subshells: i.subshells.iter().map(|s| s.label()).collect(),
                    binding_source: if i.binding_file.is_some() {
                        "file"
                    } else {
                        "embedded EADL2017"
                    },
                    binding_file: i.binding_file.as_ref(),
                    binding_provenance: i.binding.provenance(),
                }),
                elastic_table: table_out(&t.elastic, &t.elastic_origin),
                inelastic_table: table_out(&t.inelastic, &t.inelastic_origin),
            }
        })
        .collect();
    let files = Files {
        escape_spectra: SPECTRA_FILE,
        tables: TABLES_FILE,
        deposition_cartesian: r.tally.cartesian.map(|_| CARTESIAN_FILE),
        deposition_cylindrical: r.tally.cylindrical.map(|_| CYLINDRICAL_FILE),
        charge_cartesian: r.tally.cartesian.map(|_| CHARGE_CARTESIAN_FILE),
        charge_cylindrical: r.tally.cylindrical.map(|_| CHARGE_CYLINDRICAL_FILE),
        psf_profile: sim.psf.as_ref().map(|_| PSF_PROFILE_FILE),
        psf_parameters: sim.psf.as_ref().map(|_| PSF_PARAMETERS_FILE),
        depth_dose: r.tally.dose.and_then(|d| d.depth).map(|_| DEPTH_DOSE_FILE),
        // A single-layer target has no interface to cross.
        interface_transmission: transmission_csv(&sim.report).map(|_| TRANSMISSION_FILE),
        exposure_histogram: r.tally.exposure.as_ref().map(|_| EXPOSURE_FILE),
        damage_export: sim.voxel.as_ref().map(|_| DAMAGE_EXPORT_FILE),
    };
    let summary = Summary {
        format: Format {
            name: FORMAT_NAME,
            version: FORMAT_VERSION,
        },
        software: software(),
        input: &r.input,
        physics: Physics {
            models: r.models(),
            transport: &sim.metadata,
            target,
            materials,
        },
        results: {
            let mut v = results_json(&sim.report, sim.psf.as_ref())?;
            v["damage"] = damage_json(r, &sim.report);
            v
        },
        files,
        run: sim.info,
    };
    let mut text = serde_json::to_string_pretty(&summary)?;
    text.push('\n');
    Ok(text)
}

// ---------------------------------------------------------------------------
// CSV files

/// `damage_export.csv`: the electron voxel tally in the versioned exchange
/// format (engine `electron`), or `None` without `[electron.tally.voxel]`.
pub fn damage_export_csv(r: &ResolvedElectron, sim: &ElectronSimulation) -> Option<String> {
    let report = sim.voxel.as_ref()?;
    Some(crate::output::export_csv(
        report,
        "electron",
        r.input.run.seed,
        &r.input,
    ))
}

/// One row per bin (edges scaled by `scale`), then underflow and overflow
/// rows with empty densities. `prefix` starts every row.
fn hist_rows(out: &mut String, prefix: &str, h: &Histogram, histories: u64, scale: f64) {
    let b = &h.binning;
    let w = b.width() * scale;
    for (i, &c) in h.counts.iter().enumerate() {
        writeln!(
            out,
            "{prefix}{:?},{:?},{c},{:?}",
            b.edge(i) * scale,
            b.edge(i + 1) * scale,
            density(c, histories, w)
        )
        .unwrap();
    }
    writeln!(out, "{prefix}-inf,{:?},{},", b.lo * scale, h.underflow).unwrap();
    writeln!(out, "{prefix}{:?},inf,{},", b.hi * scale, h.overflow).unwrap();
}

/// `electron_escape_spectra.csv`: per face, the energy spectrum of all
/// escaping electrons (eV) and the polar-angle spectrum of each class
/// (degrees from the outward normal). Densities are per primary per eV or per
/// degree.
pub fn spectra_csv(report: &ElectronReport) -> String {
    let n = report.histories;
    let deg = 180.0 / std::f64::consts::PI;
    let mut out = String::from("face,spectrum,class,lo,hi,count,per_primary_per_unit\n");
    for (face, f) in [("front", &report.front), ("back", &report.back)] {
        hist_rows(
            &mut out,
            &format!("{face},energy_ev,all,"),
            &f.energy_histogram,
            n,
            1.0,
        );
        for (class, c) in [("slow", &f.slow), ("fast", &f.fast)] {
            hist_rows(
                &mut out,
                &format!("{face},polar_deg,{class},"),
                &c.polar_histogram,
                n,
                deg,
            );
        }
    }
    out
}

/// `electron_deposition_cartesian.csv`: deposited energy per voxel (eV, and
/// eV per primary per nm³), in [`lindhard::tally::CartesianGrid::index`]
/// order.
pub fn cartesian_csv(report: &ElectronReport) -> Option<String> {
    let d = report.deposition.cartesian.as_ref()?;
    let g = &d.grid;
    let v = g.voxel_volume_m3() / (NM * NM * NM);
    let n = report.histories.max(1) as f64;
    let mut out = String::from(
        "ix,iy,iz,x_lo_nm,x_hi_nm,y_lo_nm,y_hi_nm,z_lo_nm,z_hi_nm,energy_ev,ev_per_primary_per_nm3,ionisation_ev,lattice_ev\n",
    );
    for ix in 0..g.x.bins {
        for iy in 0..g.y.bins {
            for iz in 0..g.z.bins {
                let e = d.at(ix, iy, iz);
                let k = g.index(ix, iy, iz);
                writeln!(
                    out,
                    "{ix},{iy},{iz},{:?},{:?},{:?},{:?},{:?},{:?},{e:?},{:?},{:?},{:?}",
                    g.x.edge(ix) / NM,
                    g.x.edge(ix + 1) / NM,
                    g.y.edge(iy) / NM,
                    g.y.edge(iy + 1) / NM,
                    g.z.edge(iz) / NM,
                    g.z.edge(iz + 1) / NM,
                    e / (n * v),
                    d.ionisation_ev.get(k).copied().unwrap_or(0.0),
                    d.lattice_ev.get(k).copied().unwrap_or(0.0)
                )
                .unwrap();
            }
        }
    }
    Some(out)
}

/// `electron_deposition_cylindrical.csv`: deposited energy per `(r, depth)`
/// cell (eV, and eV per primary per nm³).
pub fn cylindrical_csv(report: &ElectronReport) -> Option<String> {
    let d = report.deposition.cylindrical.as_ref()?;
    let g = &d.grid;
    let n = report.histories.max(1) as f64;
    let mut out = String::from(
        "ir,ix,r_lo_nm,r_hi_nm,depth_lo_nm,depth_hi_nm,energy_ev,ev_per_primary_per_nm3,ionisation_ev,lattice_ev\n",
    );
    for ir in 0..g.r.bins {
        let v = g.cell_volume_m3(ir) / (NM * NM * NM);
        for ix in 0..g.depth.bins {
            let e = d.at(ir, ix);
            let k = g.index(ir, ix);
            writeln!(
                out,
                "{ir},{ix},{:?},{:?},{:?},{:?},{e:?},{:?},{:?},{:?}",
                g.r.edge(ir) / NM,
                g.r.edge(ir + 1) / NM,
                g.depth.edge(ix) / NM,
                g.depth.edge(ix + 1) / NM,
                e / (n * v),
                d.ionisation_ev.get(k).copied().unwrap_or(0.0),
                d.lattice_ev.get(k).copied().unwrap_or(0.0)
            )
            .unwrap();
        }
    }
    Some(out)
}

/// `electron_depth_dose.csv` (with `tally.dose.depth_bins_per_layer`): per depth
/// bin, the energy per primary (eV, and eV per primary per nm) and the dose per
/// unit incident fluence (Gy m², empty in a semi-infinite layer), with
/// standard errors.
pub fn depth_dose_csv(report: &ElectronReport) -> Option<String> {
    let d = report.dose.as_ref()?.depth.as_ref()?;
    let mut out = String::from(
        "bin,layer,x_lo_nm,x_hi_nm,energy_ev,energy_ev_std_err,ev_per_nm,ev_per_nm_std_err,\
         dose_gy_m2,dose_gy_m2_std_err\n",
    );
    for i in 0..d.energy_ev.len() {
        let (e, p) = (&d.energy_ev[i], &d.ev_per_nm[i]);
        let (g, ge) = match &d.dose_gy_m2[i] {
            Some(g) => (format!("{:?}", g.mean), format!("{:?}", g.std_err)),
            None => (String::new(), String::new()),
        };
        writeln!(
            out,
            "{i},{},{:?},{:?},{:?},{:?},{:?},{:?},{g},{ge}",
            d.layer_of_bin[i],
            d.edges_m[i] / NM,
            d.edges_m[i + 1] / NM,
            e.mean,
            e.std_err,
            p.mean,
            p.std_err
        )
        .unwrap();
    }
    Some(out)
}

/// `electron_interface_transmission.csv` (with `tally.dose.transmission`): one
/// row per interface, direction, particle class and energy bin of the
/// spectrum, with the crossing totals repeated on each row.
pub fn transmission_csv(report: &ElectronReport) -> Option<String> {
    let dose = report.dose.as_ref()?;
    if dose.transmission.is_empty() {
        return None;
    }
    let mut out = String::from(
        "interface,from_layer,to_layer,class,count_per_primary,count_std_err,\
         energy_ev_per_primary,energy_std_err,energy_fraction,bin,e_lo_ev,e_hi_ev,crossings\n",
    );
    for t in &dose.transmission {
        let b = &t.spectrum.binning;
        for (i, c) in t.spectrum.counts.iter().enumerate() {
            writeln!(
                out,
                "{},{},{},{},{:?},{:?},{:?},{:?},{:?},{i},{:?},{:?},{c}",
                t.interface,
                t.from_layer,
                t.to_layer,
                if t.secondary { "secondary" } else { "primary" },
                t.count.mean,
                t.count.std_err,
                t.energy_ev.mean,
                t.energy_ev.std_err,
                t.energy_fraction,
                b.edge(i),
                b.edge(i + 1)
            )
            .unwrap();
        }
    }
    Some(out)
}

/// The charge columns shared by both charge CSV files: the four counts, the
/// net electron number (stopped minus created), and the net per incident
/// electron with its standard error.
const CHARGE_COLUMNS: &str = "stopped_primary,stopped_secondary,created_by_primary,created_by_secondary,net,net_per_primary,net_per_primary_se";

fn charge_row(out: &mut String, cells: &ChargeCells, i: usize, histories: u64) {
    let c = &cells.counts[i];
    let (mean, se) = cells.net_per_primary(i, histories);
    write!(
        out,
        "{},{},{},{},{},{mean:?},{se:?}",
        c.stopped_primary,
        c.stopped_secondary,
        c.created_by_primary,
        c.created_by_secondary,
        c.net()
    )
    .unwrap();
}

/// `electron_charge_cartesian.csv`: electrons stopped and holes created per
/// voxel, by generation (see [`lindhard::tally::ChargeReport`] for the
/// convention), in [`lindhard::tally::CartesianGrid::index`] order.
pub fn charge_cartesian_csv(report: &ElectronReport) -> Option<String> {
    let d = report.charge.cartesian.as_ref()?;
    let g = &d.grid;
    let mut out =
        format!("ix,iy,iz,x_lo_nm,x_hi_nm,y_lo_nm,y_hi_nm,z_lo_nm,z_hi_nm,{CHARGE_COLUMNS}\n");
    for ix in 0..g.x.bins {
        for iy in 0..g.y.bins {
            for iz in 0..g.z.bins {
                write!(
                    out,
                    "{ix},{iy},{iz},{:?},{:?},{:?},{:?},{:?},{:?},",
                    g.x.edge(ix) / NM,
                    g.x.edge(ix + 1) / NM,
                    g.y.edge(iy) / NM,
                    g.y.edge(iy + 1) / NM,
                    g.z.edge(iz) / NM,
                    g.z.edge(iz + 1) / NM,
                )
                .unwrap();
                charge_row(&mut out, &d.cells, g.index(ix, iy, iz), report.histories);
                out.push('\n');
            }
        }
    }
    Some(out)
}

/// `electron_charge_cylindrical.csv`: electrons stopped and holes created per
/// `(r, depth)` cell. Summed over `ir`, the rows give the depth profile.
pub fn charge_cylindrical_csv(report: &ElectronReport) -> Option<String> {
    let d = report.charge.cylindrical.as_ref()?;
    let g = &d.grid;
    let mut out = format!("ir,ix,r_lo_nm,r_hi_nm,depth_lo_nm,depth_hi_nm,{CHARGE_COLUMNS}\n");
    for ir in 0..g.r.bins {
        for ix in 0..g.depth.bins {
            write!(
                out,
                "{ir},{ix},{:?},{:?},{:?},{:?},",
                g.r.edge(ir) / NM,
                g.r.edge(ir + 1) / NM,
                g.depth.edge(ix) / NM,
                g.depth.edge(ix + 1) / NM,
            )
            .unwrap();
            charge_row(&mut out, &d.cells, g.index(ir, ix), report.histories);
            out.push('\n');
        }
    }
    Some(out)
}

/// A CSV text field per RFC 4180: quoted, with embedded quotes doubled, when
/// it contains a comma, double quote, CR or LF; otherwise unchanged.
fn csv_field(s: &str) -> std::borrow::Cow<'_, str> {
    if s.contains([',', '"', '\r', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\"")).into()
    } else {
        s.into()
    }
}

/// `electron_tables.csv`: the inverse mean free paths of the tables the run
/// used, per material and grid energy, and the inelastic mean loss and
/// stopping power implied by the stored loss distribution.
pub fn tables_csv(r: &ResolvedElectron, sim: &ElectronSimulation) -> String {
    let mut out = String::from(
        "material,energy_ev,elastic_inverse_mfp_per_nm,inelastic_inverse_mfp_per_nm,\
         inelastic_mean_loss_ev,inelastic_stopping_ev_per_nm\n",
    );
    for (m, t) in r.materials.iter().zip(&sim.tables) {
        // One line per inelastic row. The elastic table has the rows of
        // electron.tables, which the inelastic table contains; the rows the
        // rate refinement adds (#339) have no elastic value.
        let elastic = t.elastic.energy_ev();
        for (i, &e) in t.inelastic.energy_ev().iter().enumerate() {
            let el = elastic
                .binary_search_by(|x| x.total_cmp(&e))
                .ok()
                .map(|j| t.elastic.inverse_mfp_per_m()[j] * NM);
            let inv = t.inelastic.inverse_mfp_per_m()[i];
            let mean = t
                .inelastic
                .quantiles(i)
                .map_or(0.0, |q| mean_loss_ev(t.inelastic.probability(), q));
            let el = el.map_or(String::new(), |x| format!("{x:?}"));
            writeln!(
                out,
                "{},{e:?},{el},{:?},{mean:?},{:?}",
                csv_field(&m.name),
                inv * NM,
                inv * mean * NM
            )
            .unwrap();
        }
    }
    out
}

#[cfg(test)]
mod csv_tests {
    use super::csv_field;

    #[test]
    fn simple_names_are_unchanged_and_special_ones_quoted() {
        assert_eq!(csv_field("Si"), "Si");
        assert_eq!(csv_field("SiO2 film"), "SiO2 film");
        assert_eq!(csv_field("a,b"), "\"a,b\"");
        assert_eq!(csv_field("a\"b"), "\"a\"\"b\"");
        assert_eq!(csv_field("a\r\nb\nc"), "\"a\r\nb\nc\"");
    }
}
