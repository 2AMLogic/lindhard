//! `summary.json` and the CSV profiles.
//!
//! Layout and compatibility rules are documented in `docs/cli.md`. In short:
//! everything except the top-level `run` object is a pure function of the
//! input (and the binary), so it is byte-identical at any thread count; new
//! results are added as new keys under `results` and new entries under
//! `files`, and existing keys keep their meaning. `format.version` changes
//! only when an existing key changes or disappears.

use std::fmt::Write as _;

use lindhard::input::TuningReport;
use lindhard::input::{ModelInfo, Resolved};
use lindhard::ion::bca::{CrystalMetadata, EnergyBudget, MeanFreePath};
use lindhard::ion::scattering::ScatteringTable;
use lindhard::material::MaterialSpec;
use lindhard::tally::{
    CascadeDefects, DualPearsonFit, Histogram, IonReport, MomentSummary, NrtDamage, PearsonIv,
};
use serde::Serialize;

use crate::dynamic::DynamicSimulation;
use crate::tally::CliTally;

/// Name and version of the summary format.
pub const FORMAT_NAME: &str = "lindhard-summary";
/// Name of the `dynamic_summary.json` format (same version counter rules).
pub const DYNAMIC_FORMAT_NAME: &str = "lindhard-dynamic-summary";
/// Bumped only on a breaking change (a key removed or changed in meaning).
pub const FORMAT_VERSION: u32 = 1;

pub const SUMMARY_FILE: &str = "summary.json";
pub const DEPTH_FILE: &str = "depth_profile.csv";
pub const IONS_FILE: &str = "ions.csv";
pub const LATERAL_FILE: &str = "lateral_profile.csv";
pub const DAMAGE_FILE: &str = "damage_profile.csv";
pub const ESCAPES_FILE: &str = "escape_spectra.csv";

pub const DYNAMIC_SUMMARY_FILE: &str = "dynamic_summary.json";
pub const DYNAMIC_STEPS_FILE: &str = "dynamic_steps.csv";
pub const DYNAMIC_COMPOSITION_FILE: &str = "dynamic_composition.csv";

pub const NM: f64 = 1e-9;

#[derive(Debug, Serialize)]
pub struct Format {
    pub name: &'static str,
    pub version: u32,
}

#[derive(Debug, Serialize)]
pub struct Software {
    pub name: &'static str,
    pub version: &'static str,
    pub git_describe: &'static str,
}

#[derive(Serialize)]
struct Engine {
    primary_cutoff_ev: f64,
    recoil_cutoff_ev: f64,
    follow_recoils: bool,
    primary_surface_binding_ev: f64,
    free_path: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    min_cm_angle_rad: Option<f64>,
    weak_collisions: u8,
    electronic_loss: &'static str,
    seed: u64,
    chunk_size: u64,
}

#[derive(Serialize)]
struct Table {
    eps_min: f64,
    eps_max: f64,
    beta_min: f64,
    beta_max: f64,
    per_decade: usize,
    max_abs_error_rad: f64,
    max_rel_error: f64,
}

#[derive(Serialize)]
struct LayerOut {
    index: usize,
    source: String,
    front_nm: f64,
    /// `null` for a semi-infinite substrate.
    back_nm: Option<f64>,
    atom_density_per_cm3: f64,
    material: MaterialSpec,
}

/// A user stopping table, identified by file, for provenance. User tables are
/// user data: the record is the path as given, where it was read from, the
/// SHA-256 of the file's bytes and the table's own provenance string.
#[derive(Serialize)]
struct StoppingTableOut {
    path: String,
    resolved_path: String,
    sha256: String,
    provenance: String,
    ion_z: u8,
    ion_mass_amu: f64,
    target_z: u8,
    energy_min_ev: f64,
    energy_max_ev: f64,
}

#[derive(Serialize)]
struct Physics {
    models: Vec<ModelInfo>,
    /// Present only when the input has a `[stopping]` table.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    stopping_tables: Vec<StoppingTableOut>,
    engine: Engine,
    scattering_table: Table,
    target: Vec<LayerOut>,
    /// Present only when `[physics] tuning` names a set: the phenomenological
    /// calibration applied, with original and effective energies.
    #[serde(skip_serializing_if = "Option::is_none")]
    tuning: Option<TuningReport>,
    /// Present only for a run with `[[crystal]]`: what the engine reports for
    /// each crystal (lattice, orientation, search and thermal parameters,
    /// `electronic_constants_unverified` where it applies). Added without a
    /// `format.version` bump (docs/cli.md, "Compatibility and extension").
    #[serde(skip_serializing_if = "Vec::is_empty")]
    crystal: Vec<CrystalMetadata>,
    /// Present only for a run with `[beam.divergence]`: the resolved
    /// distribution, width (radians and degrees), incidence policy and
    /// random-stream segment. Added without a `format.version` bump
    /// (docs/cli.md, "Compatibility and extension").
    #[serde(skip_serializing_if = "Option::is_none")]
    beam_divergence: Option<lindhard::ion::bca::DivergenceMetadata>,
}

#[derive(Serialize)]
struct Primaries {
    stopped: u64,
    backscattered: u64,
    transmitted: u64,
    /// Mean depth at rest of stopped primaries, nm (`null` if none stopped).
    /// The same quantity as `results.range.depth.mean_nm` (the projected
    /// range), kept under its original key.
    stopped_depth_mean_nm: Option<f64>,
    /// Standard deviation of that depth, nm (`results.range.depth.std_dev_nm`).
    stopped_depth_std_nm: Option<f64>,
}

#[derive(Serialize)]
struct Recoils {
    displaced: u64,
    sputtered: u64,
    transmitted: u64,
}

#[derive(Serialize)]
struct Yields {
    backscattered_per_ion: f64,
    transmitted_per_ion: f64,
    sputtered_per_ion: f64,
}

#[derive(Serialize)]
struct Budget {
    incident: f64,
    electronic_nonlocal: f64,
    electronic_local: f64,
    lattice: f64,
    surface_barrier: f64,
    backscattered: f64,
    sputtered: f64,
    transmitted: f64,
    rest: f64,
    max_relative_residual: f64,
}

#[derive(Serialize)]
struct Results {
    histories: u64,
    primaries: Primaries,
    recoils: Recoils,
    yields: Yields,
    /// Mean per primary ion, eV.
    energy_budget_ev_per_ion: Budget,
    range: Range,
    damage: Damage,
    sputtering: Sputtering,
    escapes: Escapes,
}

/// Moments of a length distribution, in nm.
#[derive(Serialize)]
struct MomentsNm {
    n: u64,
    mean_nm: f64,
    std_dev_nm: f64,
    skewness: f64,
    kurtosis: f64,
    mean_std_err_nm: f64,
    std_dev_std_err_nm: f64,
    skewness_std_err: f64,
    kurtosis_std_err: f64,
}

fn moments_nm(m: &MomentSummary) -> MomentsNm {
    MomentsNm {
        n: m.n,
        mean_nm: m.mean / NM,
        std_dev_nm: m.std_dev / NM,
        skewness: m.skewness,
        kurtosis: m.kurtosis,
        mean_std_err_nm: m.mean_std_err / NM,
        std_dev_std_err_nm: m.std_dev_std_err / NM,
        skewness_std_err: m.skewness_std_err,
        kurtosis_std_err: m.kurtosis_std_err,
    }
}

/// A Pearson IV density: its four moments and Heinrich's parameters.
#[derive(Serialize)]
struct PearsonOut {
    mean_nm: f64,
    std_dev_nm: f64,
    skewness: f64,
    kurtosis: f64,
    m: f64,
    nu: f64,
    a_nm: f64,
    lambda_nm: f64,
}

fn pearson_out(p: &PearsonIv) -> PearsonOut {
    PearsonOut {
        mean_nm: p.mean / NM,
        std_dev_nm: p.std_dev / NM,
        skewness: p.skewness,
        kurtosis: p.kurtosis,
        m: p.m,
        nu: p.nu,
        a_nm: p.a / NM,
        lambda_nm: p.lambda / NM,
    }
}

#[derive(Serialize)]
struct DualPearsonOut {
    head_fraction: f64,
    head: PearsonOut,
    tail: PearsonOut,
    chi_square: f64,
    single_chi_square: Option<f64>,
    bins: usize,
    evaluations: u64,
}

fn dual_pearson_out(f: &DualPearsonFit) -> DualPearsonOut {
    DualPearsonOut {
        head_fraction: f.profile.head_fraction,
        head: pearson_out(&f.profile.head),
        tail: pearson_out(&f.profile.tail),
        chi_square: f.chi_square,
        single_chi_square: f.single_chi_square,
        bins: f.bins,
        evaluations: f.evaluations,
    }
}

/// Depth statistics of the beam particles at rest in one layer.
#[derive(Serialize)]
struct LayerRangeOut {
    layer: usize,
    stopped: u64,
    depth: Option<MomentsNm>,
}

/// `results.range`: where the beam particles came to rest. Lengths in nm.
#[derive(Serialize)]
struct Range {
    stopped: u64,
    /// `mean_nm` is the projected range `Rp`, `std_dev_nm` the straggle
    /// `dRp`; `null` with fewer than two stopped primaries.
    depth: Option<MomentsNm>,
    /// Pearson IV with the depth moments, when they are in the type IV
    /// region; otherwise `null` and the reason is in `pearson_iv_error`.
    pearson_iv: Option<PearsonOut>,
    pearson_iv_error: Option<String>,
    /// Present only with `tally.dual_pearson = true`.
    #[serde(skip_serializing_if = "Option::is_none")]
    dual_pearson: Option<DualPearsonOut>,
    #[serde(skip_serializing_if = "Option::is_none")]
    dual_pearson_error: Option<String>,
    lateral_y: Option<MomentsNm>,
    lateral_z: Option<MomentsNm>,
    radial: Option<MomentsNm>,
    layers: Vec<LayerRangeOut>,
}

/// Per-ion damage numbers.
#[derive(Serialize)]
struct DamagePerIon {
    pka: f64,
    pka_energy_ev: f64,
    damage_energy_ev: f64,
    nrt_displacements: f64,
    kinchin_pease_displacements: f64,
    displacements: f64,
    replacements: f64,
    vacancies: f64,
    interstitials: f64,
}

#[derive(Serialize)]
struct LayerDamageOut {
    layer: usize,
    nrt: NrtDamage,
    cascade: CascadeDefects,
}

/// `results.damage`: displacement-model estimates (`nrt`) and defects counted
/// in the simulated cascades (`cascade`), kept apart because they are not the
/// same quantity (see `lindhard::ion::damage`). Totals over all histories;
/// energies in eV.
#[derive(Serialize)]
struct Damage {
    nrt: NrtDamage,
    cascade: CascadeDefects,
    per_ion: DamagePerIon,
    layers: Vec<LayerDamageOut>,
}

/// Atoms leaving through one face, of one element.
#[derive(Serialize)]
struct SputteredElement {
    z: u8,
    symbol: &'static str,
    count: u64,
    per_ion: f64,
    mean_energy_ev: f64,
}

/// `results.sputtering`: target atoms leaving through the front face, per
/// incident ion, by element.
#[derive(Serialize)]
struct Sputtering {
    yield_per_ion: f64,
    by_element: Vec<SputteredElement>,
}

#[derive(Serialize)]
struct FaceOut {
    count: u64,
    per_ion: f64,
    mean_energy_ev: f64,
}

#[derive(Serialize)]
struct SpeciesOut {
    species: usize,
    z: u8,
    symbol: &'static str,
    beam: bool,
    front: FaceOut,
    back: FaceOut,
}

/// `results.escapes`: every species leaving each face (the spectra are in
/// `escape_spectra.csv`).
#[derive(Serialize)]
struct Escapes {
    backscatter_coefficient: f64,
    transmission_coefficient: f64,
    energy_reflection_coefficient: f64,
    species: Vec<SpeciesOut>,
}

pub fn symbol(z: u8) -> &'static str {
    lindhard::elements::element(z).map_or("?", |e| e.symbol)
}

fn range_out(r: &IonReport) -> Range {
    let g = &r.range;
    Range {
        stopped: g.stopped,
        depth: g.depth.as_ref().map(moments_nm),
        pearson_iv: g.pearson_iv.as_ref().map(pearson_out),
        pearson_iv_error: g.pearson_iv_error.as_ref().map(ToString::to_string),
        dual_pearson: g.dual_pearson.as_ref().map(dual_pearson_out),
        dual_pearson_error: g.dual_pearson_error.as_ref().map(ToString::to_string),
        lateral_y: g.lateral_y.as_ref().map(moments_nm),
        lateral_z: g.lateral_z.as_ref().map(moments_nm),
        radial: g.radial.as_ref().map(moments_nm),
        layers: g
            .layers
            .iter()
            .map(|l| LayerRangeOut {
                layer: l.layer,
                stopped: l.stopped,
                depth: l.depth.as_ref().map(moments_nm),
            })
            .collect(),
    }
}

fn damage_out(r: &IonReport) -> Damage {
    let d = &r.damage;
    let n = (r.histories as f64).max(1.0);
    Damage {
        nrt: d.nrt,
        cascade: d.cascade,
        per_ion: DamagePerIon {
            pka: d.nrt.pka_count as f64 / n,
            pka_energy_ev: d.nrt.pka_energy_ev / n,
            damage_energy_ev: d.nrt.damage_energy_ev / n,
            nrt_displacements: d.nrt.nrt_displacements / n,
            kinchin_pease_displacements: d.nrt.kinchin_pease_displacements / n,
            displacements: d.cascade.displacements as f64 / n,
            replacements: d.cascade.replacements as f64 / n,
            vacancies: d.cascade.vacancies as f64 / n,
            interstitials: d.cascade.interstitials as f64 / n,
        },
        layers: d
            .layers
            .iter()
            .map(|l| LayerDamageOut {
                layer: l.layer,
                nrt: l.nrt,
                cascade: l.cascade,
            })
            .collect(),
    }
}

fn sputtering_out(r: &IonReport) -> Sputtering {
    Sputtering {
        yield_per_ion: r.escapes.sputter_yield,
        by_element: r
            .escapes
            .species
            .iter()
            .filter(|s| !s.beam)
            .map(|s| SputteredElement {
                z: s.z,
                symbol: symbol(s.z),
                count: s.front.count,
                per_ion: s.front.per_ion,
                mean_energy_ev: s.front.mean_energy_ev,
            })
            .collect(),
    }
}

fn escapes_out(r: &IonReport) -> Escapes {
    let e = &r.escapes;
    let face = |f: &lindhard::tally::FaceEscape| FaceOut {
        count: f.count,
        per_ion: f.per_ion,
        mean_energy_ev: f.mean_energy_ev,
    };
    Escapes {
        backscatter_coefficient: e.backscatter_coefficient,
        transmission_coefficient: e.transmission_coefficient,
        energy_reflection_coefficient: e.energy_reflection_coefficient,
        species: e
            .species
            .iter()
            .map(|s| SpeciesOut {
                species: s.species,
                z: s.z,
                symbol: symbol(s.z),
                beam: s.beam,
                front: face(&s.front),
                back: face(&s.back),
            })
            .collect(),
    }
}

#[derive(Serialize)]
struct Files {
    depth_profile: &'static str,
    ions: Option<&'static str>,
    lateral_profile: &'static str,
    damage_profile: &'static str,
    escape_spectra: &'static str,
}

/// The only nondeterministic part of the summary.
#[derive(Debug, Serialize, Clone, Copy)]
pub struct RunInfo {
    pub threads: usize,
    pub table_build_s: f64,
    pub transport_s: f64,
    pub ions_per_s: f64,
}

#[derive(Serialize)]
struct Summary<'a> {
    format: Format,
    software: Software,
    input: &'a lindhard::input::Input,
    physics: Physics,
    results: Results,
    files: Files,
    run: RunInfo,
}

pub fn software() -> Software {
    Software {
        name: "lindhard",
        version: lindhard::VERSION,
        git_describe: env!("LINDHARD_GIT_DESCRIBE"),
    }
}

fn budget_per_ion(b: &EnergyBudget, n: f64, max_rel: f64) -> Budget {
    Budget {
        incident: b.incident / n,
        electronic_nonlocal: b.electronic_nonlocal / n,
        electronic_local: b.electronic_local / n,
        lattice: b.lattice / n,
        surface_barrier: b.surface_barrier / n,
        backscattered: b.backscattered / n,
        sputtered: b.sputtered / n,
        transmitted: b.transmitted / n,
        rest: b.rest / n,
        max_relative_residual: max_rel,
    }
}

/// The `summary.json` text (pretty-printed, trailing newline).
pub fn summary_json(
    r: &Resolved,
    table: &ScatteringTable,
    t: &CliTally,
    report: &IonReport,
    crystals: &[CrystalMetadata],
    run: RunInfo,
) -> serde_json::Result<String> {
    let c = &r.config;
    let spec = table.spec();
    let s = &t.summary;
    let n = s.histories as f64;
    let stopped = s.primaries_stopped > 0;
    // The projected range and straggle come from the range tally's
    // accumulator, so these keys equal `results.range.depth` exactly. With a
    // single stopped primary that tally has no moments; use the summary's.
    let (depth_mean, depth_std) = match &report.range.depth {
        Some(d) => (d.mean, d.std_dev),
        None => (s.mean_depth(), s.depth_std()),
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
            stopping_tables: r
                .stopping_tables
                .iter()
                .map(|l| {
                    let (lo, hi) = l.table.energy_range_ev();
                    StoppingTableOut {
                        path: l.path.clone(),
                        resolved_path: l.resolved_path.display().to_string(),
                        sha256: l.sha256.clone(),
                        provenance: l.table.provenance().to_string(),
                        ion_z: l.table.ion_z(),
                        ion_mass_amu: l.table.ion_mass_amu(),
                        target_z: l.table.target_z(),
                        energy_min_ev: lo,
                        energy_max_ev: hi,
                    }
                })
                .collect(),
            engine: Engine {
                primary_cutoff_ev: c.primary_cutoff_ev,
                recoil_cutoff_ev: c.recoil_cutoff_ev,
                follow_recoils: c.follow_recoils,
                primary_surface_binding_ev: c.primary_surface_binding_ev,
                free_path: match c.mean_free_path {
                    MeanFreePath::Constant => "constant",
                    MeanFreePath::EnergyDependent { .. } => "energy-dependent",
                },
                min_cm_angle_rad: match c.mean_free_path {
                    MeanFreePath::Constant => None,
                    MeanFreePath::EnergyDependent { min_cm_angle_rad } => Some(min_cm_angle_rad),
                },
                weak_collisions: c.weak_collisions,
                electronic_loss: match c.electronic {
                    lindhard::ion::bca::ElectronicLoss::NonLocal => "nonlocal",
                    lindhard::ion::bca::ElectronicLoss::EquipartitionLsOr => "equipartition-ls-or",
                },
                seed: c.seed,
                chunk_size: c.chunk_size,
            },
            scattering_table: Table {
                eps_min: spec.eps_min,
                eps_max: spec.eps_max,
                beta_min: spec.beta_min,
                beta_max: spec.beta_max,
                per_decade: spec.per_decade,
                max_abs_error_rad: table.max_abs_error(),
                max_rel_error: table.max_rel_error(),
            },
            target: r
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
                .collect(),
            tuning: r.tuning.clone(),
            crystal: crystals.to_vec(),
            beam_divergence: lindhard::ion::bca::DivergenceMetadata::new(&r.divergence),
        },
        results: Results {
            histories: s.histories,
            primaries: Primaries {
                stopped: s.primaries_stopped,
                backscattered: s.backscattered,
                transmitted: s.transmitted,
                stopped_depth_mean_nm: stopped.then(|| depth_mean / NM),
                stopped_depth_std_nm: stopped.then(|| depth_std / NM),
            },
            recoils: Recoils {
                displaced: s.recoils,
                sputtered: s.sputtered,
                transmitted: s.recoils_transmitted,
            },
            yields: Yields {
                backscattered_per_ion: s.backscattered as f64 / n,
                transmitted_per_ion: s.transmitted as f64 / n,
                sputtered_per_ion: s.sputtered as f64 / n,
            },
            energy_budget_ev_per_ion: budget_per_ion(&s.budget, n, s.max_relative_residual),
            range: range_out(report),
            damage: damage_out(report),
            sputtering: sputtering_out(report),
            escapes: escapes_out(report),
        },
        files: Files {
            depth_profile: DEPTH_FILE,
            ions: t.per_ion.then_some(IONS_FILE),
            lateral_profile: LATERAL_FILE,
            damage_profile: DAMAGE_FILE,
            escape_spectra: ESCAPES_FILE,
        },
        run,
    };
    let mut text = serde_json::to_string_pretty(&summary)?;
    text.push('\n');
    Ok(text)
}

/// `depth_profile.csv`: stopped primaries per depth bin. The last bin also
/// holds everything deeper; its upper edge is written as `inf` and its
/// density is left empty.
pub fn depth_csv(t: &CliTally) -> String {
    let s = &t.summary;
    let w = s.bin_width_m / NM;
    let n = s.histories as f64;
    let last = s.depth_hist.len() - 1;
    let mut out = String::from("depth_lo_nm,depth_hi_nm,stopped_primaries,fraction_per_nm\n");
    for (i, &c) in s.depth_hist.iter().enumerate() {
        let lo = i as f64 * w;
        if i == last {
            writeln!(out, "{lo:?},inf,{c},").unwrap();
        } else {
            let hi = (i + 1) as f64 * w;
            let f = c as f64 / (n * w);
            writeln!(out, "{lo:?},{hi:?},{c},{f:?}").unwrap();
        }
    }
    out
}

/// `ions.csv`: final state of every primary, by history index.
pub fn ions_csv(t: &CliTally) -> String {
    let mut out = String::from("index,fate,x_nm,y_nm,z_nm,energy_ev,dir_x,dir_y,dir_z,layer\n");
    for f in &t.finals {
        writeln!(
            out,
            "{},{},{:?},{:?},{:?},{:?},{:?},{:?},{:?},{}",
            f.index,
            f.fate.as_str(),
            f.pos[0] / NM,
            f.pos[1] / NM,
            f.pos[2] / NM,
            f.energy_ev,
            f.dir[0],
            f.dir[1],
            f.dir[2],
            f.layer
        )
        .unwrap();
    }
    out
}

/// Per-ion, per-nm density of a histogram bin.
pub fn density(count: u64, histories: u64, width_nm: f64) -> f64 {
    count as f64 / (histories.max(1) as f64 * width_nm)
}

/// `lateral_profile.csv`: stopped primaries by lateral position `y`, `z` and
/// radial distance (long format). Rows of each quantity are in bin order,
/// followed by its `underflow` and `overflow` rows (edges `-inf`/`inf`).
pub fn lateral_csv(r: &IonReport) -> String {
    let g = &r.range;
    let mut out = String::from("quantity,lo_nm,hi_nm,count,per_ion_per_nm\n");
    for (name, h) in [
        ("y", &g.lateral_y_histogram),
        ("z", &g.lateral_z_histogram),
        ("radial", &g.radial_histogram),
    ] {
        hist_rows(&mut out, &format!("{name},"), h, r.histories, NM, 1.0);
    }
    out
}

/// Append one row per bin of `h` (edges divided by `unit`) then underflow and
/// overflow rows. `prefix` starts every row.
fn hist_rows(out: &mut String, prefix: &str, h: &Histogram, histories: u64, unit: f64, scale: f64) {
    let b = &h.binning;
    let w = b.width() / unit * scale;
    for (i, &c) in h.counts.iter().enumerate() {
        let (lo, hi) = (b.edge(i) / unit * scale, b.edge(i + 1) / unit * scale);
        writeln!(
            out,
            "{prefix}{lo:?},{hi:?},{c},{:?}",
            density(c, histories, w)
        )
        .unwrap();
    }
    writeln!(
        out,
        "{prefix}-inf,{:?},{},",
        b.lo / unit * scale,
        h.underflow
    )
    .unwrap();
    writeln!(out, "{prefix}{:?},inf,{},", b.hi / unit * scale, h.overflow).unwrap();
}

/// `damage_profile.csv`: simulated cascade defects by depth. `vacancies` is
/// displacements minus replacements. Densities are per incident ion per nm.
/// The last row is everything deeper than the grid.
pub fn damage_csv(r: &IonReport) -> String {
    let d = &r.damage;
    let n = r.histories;
    let b = &d.vacancy_histogram.binning;
    let w = b.width() / NM;
    let mut out = String::from(
        "depth_lo_nm,depth_hi_nm,vacancies,interstitials,replacements,\
         vacancies_per_ion_per_nm,interstitials_per_ion_per_nm,replacements_per_ion_per_nm\n",
    );
    for i in 0..b.bins {
        let (v, it, rp) = (
            d.vacancy_histogram.counts[i],
            d.interstitial_histogram.counts[i],
            d.replacement_histogram.counts[i],
        );
        writeln!(
            out,
            "{:?},{:?},{v},{it},{rp},{:?},{:?},{:?}",
            b.edge(i) / NM,
            b.edge(i + 1) / NM,
            density(v, n, w),
            density(it, n, w),
            density(rp, n, w)
        )
        .unwrap();
    }
    writeln!(
        out,
        "{:?},inf,{},{},{},,,",
        b.hi / NM,
        d.vacancy_histogram.overflow,
        d.interstitial_histogram.overflow,
        d.replacement_histogram.overflow
    )
    .unwrap();
    out
}

/// `escape_spectra.csv`: energy and polar-angle spectra of every species
/// leaving each face (long format, in species, face, spectrum order). Energy
/// bins are in eV and angle bins in degrees from the outward surface normal;
/// `per_ion_per_unit` is per incident ion per eV or per degree. Each spectrum
/// is followed by its `underflow` and `overflow` rows.
pub fn escapes_csv(r: &IonReport) -> String {
    let mut out =
        String::from("species_z,symbol,beam,face,spectrum,lo,hi,count,per_ion_per_unit\n");
    for s in &r.escapes.species {
        for (face, f) in [("front", &s.front), ("back", &s.back)] {
            let prefix = |spec: &str| format!("{},{},{},{face},{spec},", s.z, symbol(s.z), s.beam);
            hist_rows(
                &mut out,
                &prefix("energy_ev"),
                &f.energy_histogram,
                r.histories,
                1.0,
                1.0,
            );
            hist_rows(
                &mut out,
                &prefix("polar_deg"),
                &f.polar_histogram,
                r.histories,
                1.0,
                180.0 / std::f64::consts::PI,
            );
        }
    }
    out
}

/// `dynamic_steps.csv`: one row per accepted fluence step (step 0 is the
/// initial target), with cumulative yields. Layout in `docs/cli.md`.
pub fn dynamic_steps_csv(d: &DynamicSimulation) -> String {
    let mut out = String::from(
        "step,first_index,ions,ions_done,fluence_cm2,attempts,max_change,clamped,\
         removed_slabs,n_slabs,surface_nm,thickness_nm,cum_backscattered,cum_transmitted,\
         cum_stopped_in_target,cum_stopped_in_substrate,cum_sputtered,sputter_yield,\
         cum_recoils_transmitted",
    );
    for &z in &d.species {
        write!(out, ",cum_sputtered_{}", symbol(z)).unwrap();
    }
    out.push('\n');
    for s in &d.steps {
        let y = &s.cumulative;
        let thickness: f64 = s.slabs.iter().map(|b| b.back_m - b.front_m).sum();
        let sputtered = y.sputtered_total();
        let per_ion = if y.histories > 0 {
            format!("{:?}", sputtered as f64 / y.histories as f64)
        } else {
            String::new()
        };
        write!(
            out,
            "{},{},{},{},{:?},{},{:?},{},{},{},{:?},{:?},{},{},{},{},{},{},{}",
            s.step,
            s.first_index,
            s.ions,
            s.ions_done,
            s.fluence_m2 * 1e-4,
            s.attempts,
            s.max_change,
            s.clamped,
            s.removed_slabs,
            s.slabs.len(),
            s.recession_m / NM,
            thickness / NM,
            y.backscattered,
            y.transmitted,
            y.primaries_in_slabs,
            y.primaries_in_substrate,
            sputtered,
            per_ion,
            y.recoils_transmitted,
        )
        .unwrap();
        for &z in &d.species {
            write!(out, ",{}", y.sputtered.get(&z).copied().unwrap_or(0)).unwrap();
        }
        out.push('\n');
    }
    out
}

/// `dynamic_composition.csv`: the slab profile after every step (step 0 is the
/// initial target), one row per step and slab.
pub fn dynamic_composition_csv(d: &DynamicSimulation) -> String {
    let mut out = String::from("step,slab,front_nm,back_nm,thickness_nm");
    for &z in &d.species {
        write!(out, ",atoms_per_cm2_{}", symbol(z)).unwrap();
    }
    for &z in &d.species {
        write!(out, ",fraction_{}", symbol(z)).unwrap();
    }
    out.push('\n');
    for s in &d.steps {
        for (i, b) in s.slabs.iter().enumerate() {
            write!(
                out,
                "{},{},{:?},{:?},{:?}",
                s.step,
                i,
                b.front_m / NM,
                b.back_m / NM,
                (b.back_m - b.front_m) / NM
            )
            .unwrap();
            let get = |z: u8| {
                b.inventory
                    .iter()
                    .find(|&&(zz, _)| zz == z)
                    .map_or(0.0, |&(_, a)| a)
            };
            let total: f64 = b.inventory.iter().map(|&(_, a)| a).sum();
            for &z in &d.species {
                write!(out, ",{:?}", get(z) * 1e-4).unwrap();
            }
            for &z in &d.species {
                write!(out, ",{:?}", get(z) / total).unwrap();
            }
            out.push('\n');
        }
    }
    out
}

#[derive(Serialize)]
struct DynamicFiles {
    steps: &'static str,
    composition: &'static str,
}

#[derive(Serialize)]
struct DynamicTotals {
    ions: u64,
    steps: u64,
    rejected_attempts: u64,
    fluence_cm2: f64,
    n_slabs_initial: usize,
    n_slabs_final: usize,
    thickness_initial_nm: f64,
    thickness_final_nm: f64,
    backscattered_per_ion: f64,
    transmitted_per_ion: f64,
    sputtered_per_ion: f64,
    /// Present only with erosion on, so the summary of a run without it is
    /// byte-identical to what it was before erosion existed.
    #[serde(skip_serializing_if = "Option::is_none")]
    recession_nm: Option<f64>,
}

#[derive(Serialize)]
struct DynamicSummary<'a> {
    format: Format,
    software: Software,
    input: &'a lindhard::input::Input,
    models: Vec<ModelInfo>,
    species: Vec<&'static str>,
    totals: DynamicTotals,
    files: DynamicFiles,
    run: RunInfo,
}

/// `dynamic_summary.json`: the echoed input, the models, run totals and the
/// file list. Everything except the `run` object is a pure function of the
/// input, so it is byte-identical at any thread count.
pub fn dynamic_summary_json(r: &Resolved, d: &DynamicSimulation) -> serde_json::Result<String> {
    let first = d.steps.first().expect("step 0");
    let last = d.steps.last().expect("step 0");
    let thickness = |s: &crate::dynamic::StepRow| -> f64 {
        s.slabs.iter().map(|b| b.back_m - b.front_m).sum::<f64>() / NM
    };
    let n = last.cumulative.histories.max(1) as f64;
    let summary = DynamicSummary {
        format: Format {
            name: DYNAMIC_FORMAT_NAME,
            version: FORMAT_VERSION,
        },
        software: software(),
        input: &r.input,
        models: r.models(),
        species: d.species.iter().map(|&z| symbol(z)).collect(),
        totals: DynamicTotals {
            ions: last.ions_done,
            steps: last.step,
            rejected_attempts: d.rejected_attempts,
            fluence_cm2: last.fluence_m2 * 1e-4,
            n_slabs_initial: first.slabs.len(),
            n_slabs_final: last.slabs.len(),
            thickness_initial_nm: thickness(first),
            thickness_final_nm: thickness(last),
            backscattered_per_ion: last.cumulative.backscattered as f64 / n,
            transmitted_per_ion: last.cumulative.transmitted as f64 / n,
            sputtered_per_ion: last.cumulative.sputtered_total() as f64 / n,
            recession_nm: r
                .input
                .dynamic
                .as_ref()
                .is_some_and(|d| d.erosion)
                .then(|| last.recession_m / NM),
        },
        files: DynamicFiles {
            steps: DYNAMIC_STEPS_FILE,
            composition: DYNAMIC_COMPOSITION_FILE,
        },
        run: d.info,
    };
    let mut s = serde_json::to_string_pretty(&summary)?;
    s.push('\n');
    Ok(s)
}
