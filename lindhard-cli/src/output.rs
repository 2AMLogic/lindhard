//! `summary.json` and the CSV profiles.
//!
//! Layout and compatibility rules are documented in `docs/cli.md`. In short:
//! everything except the top-level `run` object is a pure function of the
//! input (and the binary), so it is byte-identical at any thread count; new
//! results are added as new keys under `results` and new entries under
//! `files`, and existing keys keep their meaning. `format.version` changes
//! only when an existing key changes or disappears.

use std::fmt::Write as _;

use lindhard::input::{ModelInfo, Resolved};
use lindhard::ion::bca::{EnergyBudget, MeanFreePath};
use lindhard::ion::scattering::ScatteringTable;
use lindhard::material::MaterialSpec;
use serde::Serialize;

use crate::tally::CliTally;

/// Name and version of the summary format.
pub const FORMAT_NAME: &str = "lindhard-summary";
/// Bumped only on a breaking change (a key removed or changed in meaning).
pub const FORMAT_VERSION: u32 = 1;

pub const SUMMARY_FILE: &str = "summary.json";
pub const DEPTH_FILE: &str = "depth_profile.csv";
pub const IONS_FILE: &str = "ions.csv";

const NM: f64 = 1e-9;

#[derive(Serialize)]
pub struct Format {
    pub name: &'static str,
    pub version: u32,
}

#[derive(Serialize)]
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

#[derive(Serialize)]
struct Physics {
    models: Vec<ModelInfo>,
    engine: Engine,
    scattering_table: Table,
    target: Vec<LayerOut>,
}

#[derive(Serialize)]
struct Primaries {
    stopped: u64,
    backscattered: u64,
    transmitted: u64,
    /// Mean depth at rest of stopped primaries, nm (`null` if none stopped).
    stopped_depth_mean_nm: Option<f64>,
    /// Standard deviation of that depth, nm.
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
}

#[derive(Serialize)]
struct Files {
    depth_profile: &'static str,
    ions: Option<&'static str>,
}

/// The only nondeterministic part of the summary.
#[derive(Serialize)]
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
    run: RunInfo,
) -> serde_json::Result<String> {
    let c = &r.config;
    let spec = table.spec();
    let s = &t.summary;
    let n = s.histories as f64;
    let stopped = s.primaries_stopped > 0;
    let summary = Summary {
        format: Format {
            name: FORMAT_NAME,
            version: FORMAT_VERSION,
        },
        software: software(),
        input: &r.input,
        physics: Physics {
            models: r.models(),
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
        },
        results: Results {
            histories: s.histories,
            primaries: Primaries {
                stopped: s.primaries_stopped,
                backscattered: s.backscattered,
                transmitted: s.transmitted,
                stopped_depth_mean_nm: stopped.then(|| s.mean_depth() / NM),
                stopped_depth_std_nm: stopped.then(|| s.depth_std() / NM),
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
        },
        files: Files {
            depth_profile: DEPTH_FILE,
            ions: t.per_ion.then_some(IONS_FILE),
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
