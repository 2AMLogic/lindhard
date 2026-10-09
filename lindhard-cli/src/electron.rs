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
use lindhard::electron::data::CrossSectionTable;
use lindhard::electron::elastic::table::{
    combine, default_probability_grid, AtomicElastic, PotentialSource, SalvatDhfsTable,
    ThomasFermiYukawa, DEFAULT_REFINE_TOLERANCE,
};
use lindhard::electron::elastic::{SalvatDhfs, SolverOptions};
use lindhard::electron::inelastic::table::{
    build_inelastic_table_for_model, mean_loss_ev, InelasticTableOptions,
};
use lindhard::electron::inelastic::PennInelastic;
use lindhard::electron::transport::{LayerTables, RunMetadata, Transport};
use lindhard::input::electron::{
    DataFile, ElectronInput, PotentialChoice, ResolvedElectron, ResolvedElectronMaterial,
};
use lindhard::input::ModelInfo;
use lindhard::material::{Material, MaterialSpec};
use lindhard::tally::{
    fit_psf, ElectronReport, FullElectronTally, Histogram, PsfFitOptions, PsfModel, PsfReport,
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
pub const PSF_PROFILE_FILE: &str = "electron_psf_profile.csv";
pub const PSF_PARAMETERS_FILE: &str = "electron_psf_parameters.csv";
pub const TABLES_FILE: &str = "electron_tables.csv";

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
    let build = || match kind {
        TableKind::Elastic => elastic_table(r, &m.material),
        TableKind::Inelastic => inelastic_table(r, m),
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

/// The elastic table of one material.
///
/// Each element's potential doubles as the electron density the optional
/// exchange and correlation-polarization corrections need
/// ([`AtomicElastic::compute_corrected`]): the Yukawa stand-in's own Poisson
/// density, or the DHFS Poisson density of Salvat et al. (1987) Eq. (12)
/// (`lindhard::electron::elastic::corrections`). Without corrections this is
/// the plain partial-wave table of `build_elastic_table`.
fn elastic_table(r: &ResolvedElectron, m: &Material) -> Result<CrossSectionTable> {
    let grid = &r.table_energy_ev;
    let solver = SolverOptions::default();
    let desc = match r.elastic.potential {
        PotentialChoice::SalvatDhfs => SalvatDhfsTable.description(),
        PotentialChoice::ThomasFermiYukawa => ThomasFermiYukawa.description(),
    };
    let desc = format!("{desc}{}", r.elastic.corrections_description());
    let mut atoms = Vec::new();
    for c in m.components() {
        if c.atom_fraction() <= 0.0 {
            continue;
        }
        let z = c.z();
        let corrections = r.elastic.corrections(z);
        atoms.push(match r.elastic.potential {
            PotentialChoice::SalvatDhfs => {
                let p = SalvatDhfs::for_element(u32::from(z))?;
                AtomicElastic::compute_corrected(z, &p, &p, &desc, grid, &corrections, solver)?
            }
            PotentialChoice::ThomasFermiYukawa => {
                let y = ThomasFermiYukawa::yukawa(z)?;
                AtomicElastic::compute_corrected(z, &y, &y, &desc, grid, &corrections, solver)?
            }
        });
    }
    Ok(combine(
        m,
        &atoms,
        &default_probability_grid(),
        Some(DEFAULT_REFINE_TOLERANCE),
    )?)
}

/// The inelastic table of one material.
fn inelastic_table(
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
) -> Result<CrossSectionTable> {
    let model = PennInelastic::try_new(r.inelastic, m.optical_elf.clone())?
        .with_fermi_energy_ev(r.inelastic_fermi_ev)?;
    Ok(build_inelastic_table_for_model(
        &model,
        &m.material,
        &InelasticTableOptions::new(r.table_energy_ev.clone()),
    )?)
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

    let t0 = Instant::now();
    let tables: Vec<MaterialTables> = pool.install(|| {
        r.materials
            .iter()
            .map(|m| {
                let (elastic, elastic_origin) = table(r, m, TableKind::Elastic, cache)?;
                let (inelastic, inelastic_origin) = table(r, m, TableKind::Inelastic, cache)?;
                Ok(MaterialTables {
                    elastic,
                    inelastic,
                    elastic_origin,
                    inelastic_origin,
                })
            })
            .collect::<Result<_>>()
    })?;
    let table_build_s = t0.elapsed().as_secs_f64();

    let layer_tables: Vec<LayerTables> = r
        .layer_material
        .iter()
        .map(|&i| LayerTables {
            elastic: tables[i].elastic.clone(),
            inelastic: tables[i].inelastic.clone(),
        })
        .collect();
    let bands: Option<Vec<_>> = r
        .layer_material
        .iter()
        .map(|&i| r.materials[i].band.clone())
        .collect();
    let mut transport = match bands {
        Some(b) => Transport::with_band_structures(r.stack.clone(), layer_tables, b, r.config),
        None => Transport::new(r.stack.clone(), layer_tables, r.config),
    }
    .context("setting up the electron transport")?;
    for (layer, &i) in r.layer_material.iter().enumerate() {
        let c = &r.materials[i].channels;
        if c.any() {
            transport = transport.with_insulator_channels(layer, c.clone())?;
        }
    }
    let proto =
        FullElectronTally::new(&transport, r.tally).context("setting up the electron tally")?;

    let t1 = Instant::now();
    let run = pool
        .install(|| {
            transport.run(
                r.input.run.seed,
                r.input.run.histories,
                CHUNK_SIZE,
                &r.primary,
                || proto.clone(),
            )
        })
        .context("electron transport failed")?;
    let transport_s = t1.elapsed().as_secs_f64();
    let report = run.tally.report();
    let psf = psf_outcome(r, &report);
    Ok(ElectronSimulation {
        tables,
        metadata: run.metadata,
        report,
        psf,
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

#[derive(Serialize)]
struct MaterialOut<'a> {
    name: &'a str,
    optical_elf: ElfOut<'a>,
    band: Option<&'a lindhard::electron::boundary::BandStructure>,
    phonon: Option<&'a lindhard::electron::phonon::FrohlichPhonon>,
    polaron: Option<&'a lindhard::electron::phonon::PolaronTrapping>,
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
    psf_profile: Option<&'static str>,
    psf_parameters: Option<&'static str>,
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
        }
    }
    Ok(v)
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
        psf_profile: sim.psf.as_ref().map(|_| PSF_PROFILE_FILE),
        psf_parameters: sim.psf.as_ref().map(|_| PSF_PARAMETERS_FILE),
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
        results: results_json(&sim.report, sim.psf.as_ref())?,
        files,
        run: sim.info,
    };
    let mut text = serde_json::to_string_pretty(&summary)?;
    text.push('\n');
    Ok(text)
}

// ---------------------------------------------------------------------------
// CSV files

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
        "ix,iy,iz,x_lo_nm,x_hi_nm,y_lo_nm,y_hi_nm,z_lo_nm,z_hi_nm,energy_ev,ev_per_primary_per_nm3\n",
    );
    for ix in 0..g.x.bins {
        for iy in 0..g.y.bins {
            for iz in 0..g.z.bins {
                let e = d.at(ix, iy, iz);
                writeln!(
                    out,
                    "{ix},{iy},{iz},{:?},{:?},{:?},{:?},{:?},{:?},{e:?},{:?}",
                    g.x.edge(ix) / NM,
                    g.x.edge(ix + 1) / NM,
                    g.y.edge(iy) / NM,
                    g.y.edge(iy + 1) / NM,
                    g.z.edge(iz) / NM,
                    g.z.edge(iz + 1) / NM,
                    e / (n * v)
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
        "ir,ix,r_lo_nm,r_hi_nm,depth_lo_nm,depth_hi_nm,energy_ev,ev_per_primary_per_nm3\n",
    );
    for ir in 0..g.r.bins {
        let v = g.cell_volume_m3(ir) / (NM * NM * NM);
        for ix in 0..g.depth.bins {
            let e = d.at(ir, ix);
            writeln!(
                out,
                "{ir},{ix},{:?},{:?},{:?},{:?},{e:?},{:?}",
                g.r.edge(ir) / NM,
                g.r.edge(ir + 1) / NM,
                g.depth.edge(ix) / NM,
                g.depth.edge(ix + 1) / NM,
                e / (n * v)
            )
            .unwrap();
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
        // Both tables share the grid of electron.tables.
        for (i, &e) in t.inelastic.energy_ev().iter().enumerate() {
            let el = t
                .elastic
                .energy_ev()
                .get(i)
                .filter(|&&x| x == e)
                .map(|_| t.elastic.inverse_mfp_per_m()[i] * NM);
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

#[cfg(test)]
mod elastic_tests {
    use super::*;
    use std::path::Path;

    /// The CLI example on a two-point grid with the given potential and
    /// exchange setting.
    fn resolved(potential: PotentialChoice, exchange: bool) -> ResolvedElectron {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples/electron");
        let text = std::fs::read_to_string(dir.join("e_10keV_si.toml")).unwrap();
        let mut r = ElectronInput::from_toml_str(&text)
            .unwrap()
            .resolve_in(&dir)
            .unwrap();
        r.elastic.potential = potential;
        r.elastic.exchange = exchange;
        r.table_energy_ev = vec![100.0, 1000.0];
        r
    }

    fn table(potential: PotentialChoice, exchange: bool) -> CrossSectionTable {
        let r = resolved(potential, exchange);
        elastic_table(&r, &r.materials[0].material).unwrap()
    }

    /// The corrections reach the DHFS table as they reach the stand-in's
    /// (#169: they were dropped for `salvat-dhfs`).
    #[test]
    fn exchange_applies_to_both_potentials() {
        for p in [
            PotentialChoice::SalvatDhfs,
            PotentialChoice::ThomasFermiYukawa,
        ] {
            let (off, on) = (table(p, false), table(p, true));
            assert!(!off.model().contains("exchange"), "{p:?}: {}", off.model());
            assert!(on.model().contains("exchange"), "{p:?}: {}", on.model());
            assert_ne!(
                off.inverse_mfp_per_m(),
                on.inverse_mfp_per_m(),
                "{p:?}: exchange left the table unchanged"
            );
        }
        let dhfs = table(PotentialChoice::SalvatDhfs, false);
        assert!(dhfs.model().contains("Salvat"), "{}", dhfs.model());
    }
}
