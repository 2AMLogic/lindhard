//! Shared fixture of the electron benchmarks (`benches/electron.rs`), the
//! thread-scaling example (`examples/electron_scaling.rs`) and its
//! determinism test (`tests/electron_bench_determinism.rs`).
//!
//! The problems are the matched problems of #150
//! (`validation/oracles/electron_problems.json`): electrons at normal incidence
//! into bulk Si and bulk Cu at 1, 5 and 20 keV, with the physics lindhard runs
//! there (`"lindhard"` block of that file): Mott elastic scattering on the
//! Thomas-Fermi Yukawa stand-in with the Furness-McCarthy exchange correction,
//! single-pole Penn inelastic, Kieft-Bosch secondaries, the step barrier with
//! quantum transmission and refraction, cutoff 0.01 eV above the vacuum level,
//! tables from 10 eV at 20 points per decade, the escape spectra and a
//! cylindrical deposition grid. Nothing here is physics: the input goes
//! through [`lindhard::input::electron`] and the tables are built the way
//! `lindhard-cli` builds them (`lindhard-cli/src/electron.rs`).
//!
//! **Material data.** No optical or band data of real materials is in this
//! tree (`docs/data-provenance.md`), so by default the fixture writes a
//! SYNTHETIC input to a temporary directory: the synthetic Drude plasmon ELF
//! of the other electron benches (`E_p = 20 eV`, width 5 eV; not optical data
//! of any material) and the synthetic band parameters of the electron test
//! fixtures (Si: insulator `W_v = 10`, `E_g = 2`, `chi = 3` eV, as
//! `examples/electron/e_10keV_si.toml`; Cu: metal `E_F = 5`, `Phi = 4` eV, as
//! `tests/electron_secondaries.rs`). Composition and density are the element
//! table's. The timings then show the cost of the engine on these problems'
//! energies, geometry and settings, not on the materials' real cross sections.
//!
//! With `LINDHARD_BENCH_ELECTRON_INPUTS=<dir>`, the fixture instead reads
//! `<dir>/<problem>/input.toml`, the matched input that
//! `validation/oracles/bench_electron.py --write-inputs` writes from the
//! user's cstool clone (outside this tree, into the gitignored
//! `validation/oracle-runs/`), so the bench runs on the same material data as
//! the Nebula comparison.

#![allow(dead_code)]

use std::error::Error;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use lindhard::electron::data::CrossSectionTable;
use lindhard::electron::elastic::table::{
    build_elastic_table, combine, default_probability_grid, AtomicElastic, ElasticTableOptions,
    PotentialSource, SalvatDhfsTable, ThomasFermiYukawa, DEFAULT_REFINE_TOLERANCE,
};
use lindhard::electron::elastic::SolverOptions;
use lindhard::electron::inelastic::table::{
    build_inelastic_table_for_model, EnergyAxis, InelasticTableOptions,
};
use lindhard::electron::inelastic::{DrudeLorentz, DrudeLorentzOscillator, PennInelastic};
use lindhard::electron::transport::{LayerTables, Transport};
use lindhard::input::electron::{
    ElectronInput, PotentialChoice, ResolvedElectron, ResolvedElectronMaterial,
};
use lindhard::material::Material;
use lindhard::tally::{ElectronReport, FullElectronTally};
use sha2::{Digest, Sha256};

pub(crate) type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

/// Histories per work chunk, as `lindhard-cli` (`CHUNK_SIZE` in
/// `lindhard-cli/src/electron.rs`). Fixed, so results do not depend on the
/// thread count.
pub(crate) const CHUNK_SIZE: u64 = 16;

/// Seed of every benchmark run.
pub(crate) const SEED: u64 = 1;

/// Environment variable naming a directory of matched inputs.
pub(crate) const INPUTS_ENV: &str = "LINDHARD_BENCH_ELECTRON_INPUTS";

/// One benchmark problem: the #150 problem id, the target element and the
/// beam energy, and the histories one Criterion iteration runs (chosen so an
/// iteration takes of the order of a second on one thread).
#[derive(Debug, Clone, Copy)]
pub(crate) struct ElectronProblem {
    pub id: &'static str,
    pub symbol: &'static str,
    pub energy_ev: f64,
    pub bench_histories: u64,
}

pub(crate) const PROBLEMS: [ElectronProblem; 6] = [
    ElectronProblem {
        id: "e_1keV_si",
        symbol: "Si",
        energy_ev: 1.0e3,
        bench_histories: 1000,
    },
    ElectronProblem {
        id: "e_5keV_si",
        symbol: "Si",
        energy_ev: 5.0e3,
        bench_histories: 200,
    },
    ElectronProblem {
        id: "e_20keV_si",
        symbol: "Si",
        energy_ev: 2.0e4,
        bench_histories: 50,
    },
    ElectronProblem {
        id: "e_1keV_cu",
        symbol: "Cu",
        energy_ev: 1.0e3,
        bench_histories: 1000,
    },
    ElectronProblem {
        id: "e_5keV_cu",
        symbol: "Cu",
        energy_ev: 5.0e3,
        bench_histories: 200,
    },
    ElectronProblem {
        id: "e_20keV_cu",
        symbol: "Cu",
        energy_ev: 2.0e4,
        bench_histories: 50,
    },
];

pub(crate) fn problem(id: &str) -> Option<ElectronProblem> {
    PROBLEMS.iter().copied().find(|p| p.id == id)
}

/// Settings of the synthetic input; [`Settings::matched`] is the #150 one.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Settings {
    pub points_per_decade: f64,
    /// Outer edge (nm) and bins of the radial deposition grid.
    pub rmax_nm: f64,
    pub rbins: usize,
}

impl Settings {
    /// `validation/oracles/electron_problems.json`: 20 points per decade, a
    /// radial grid out to `E0 / (1 eV)` nm in 10^4 bins.
    pub(crate) fn matched(p: &ElectronProblem) -> Self {
        Self {
            points_per_decade: 20.0,
            rmax_nm: p.energy_ev,
            rbins: 10_000,
        }
    }
}

/// The synthetic Drude plasmon ELF as a data file (`DrudeLorentz`, the
/// parameters of `benches/penn.rs` and the CLI example; not optical data of
/// any material), sampled at 80 log-spaced energies from 0.5 eV to 30 keV.
fn synthetic_elf_toml() -> Result<String> {
    let elf = DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])?.to_optical_elf(
        "synthetic Drude plasmon (electron benchmarks)",
        0.5,
        3.0e4,
        80,
    )?;
    let list = |v: &[f64]| {
        let items: Vec<String> = v.iter().map(|x| format!("{x:?}")).collect();
        format!("[{}]", items.join(", "))
    };
    Ok(format!(
        "material = \"synthetic Drude plasmon (electron benchmarks)\"\n\
         provenance = \"synthetic Drude-Lorentz ELF, not physical data: (A = 400 eV^2, E = 20 eV, \
         gamma = 5 eV); 80 log-spaced samples from 0.5 to 30000 eV (lindhard/benches/electron_common)\"\n\
         energy_ev = {}\nelf = {}\n",
        list(elf.energy_ev()),
        list(elf.elf_values())
    ))
}

/// Synthetic band parameters (see the module docs); not data of the element.
fn synthetic_band(symbol: &str) -> &'static str {
    match symbol {
        "Cu" => {
            "{ kind = \"metal\", fermi_ev = 5.0, work_function_ev = 4.0, provenance = \"SYNTHETIC \
             round numbers of the electron test fixtures (lindhard/benches/electron_common), not band \
             parameters of Cu\" }"
        }
        _ => {
            "{ kind = \"insulator\", valence_band_width_ev = 10.0, band_gap_ev = 2.0, affinity_ev = 3.0, \
             provenance = \"SYNTHETIC round numbers of the CLI example (lindhard/benches/electron_common), \
             not band parameters of Si\" }"
        }
    }
}

/// The synthetic input of a problem (see the module docs).
pub(crate) fn synthetic_input_toml(p: &ElectronProblem, s: &Settings, histories: u64) -> String {
    let mut t = String::new();
    let _ = write!(
        t,
        "# Written by lindhard/benches/electron_common: SYNTHETIC electron data.\n\
         [electron.beam]\nenergy_ev = {e:?}\n\
         [electron.transport]\ncutoff_ev = 0.01\ncutoff_reference = \"vacuum-level\"\n\
         secondaries = \"kieft-bosch\"\nboundary = \"step-barrier\"\n\
         [electron.elastic]\nmodel = \"mott\"\npotential = \"thomas-fermi-yukawa\"\nexchange = true\n\
         [electron.inelastic]\nmodel = \"penn-single-pole\"\n\
         [electron.tables]\nmin_energy_ev = 10.0\npoints_per_decade = {ppd:?}\n\
         [electron.materials.{sym}]\noptical_elf = \"synthetic_elf.toml\"\nband = {band}\n\
         [electron.tally]\nse_bse_split_ev = 50.0\n\
         [electron.tally.cylindrical]\nr = {{ lo_nm = 0.0, hi_nm = {rmax:?}, bins = {rbins} }}\n\
         depth = {{ lo_nm = 0.0, hi_nm = 100000000.0, bins = 1 }}\n\
         [target]\nsubstrate = \"{sym}\"\n\
         [run]\nhistories = {histories}\nseed = {SEED}\n",
        e = p.energy_ev,
        ppd = s.points_per_decade,
        sym = p.symbol,
        band = synthetic_band(p.symbol),
        rmax = s.rmax_nm,
        rbins = s.rbins,
    );
    t
}

/// Write the synthetic input of a problem into `dir` and return its path.
pub(crate) fn write_synthetic(
    dir: &Path,
    p: &ElectronProblem,
    s: &Settings,
    histories: u64,
) -> Result<PathBuf> {
    std::fs::create_dir_all(dir)?;
    std::fs::write(dir.join("synthetic_elf.toml"), synthetic_elf_toml()?)?;
    let path = dir.join(format!("{}.toml", p.id));
    std::fs::write(&path, synthetic_input_toml(p, s, histories))?;
    Ok(path)
}

/// Read and validate an electron input file.
pub(crate) fn resolve(path: &Path) -> Result<ResolvedElectron> {
    let text = std::fs::read_to_string(path)?;
    let input = ElectronInput::from_toml_str(&text)?;
    let base = path.parent().unwrap_or(Path::new("."));
    Ok(input.resolve_in(base)?)
}

/// The input of a problem: the matched one from [`INPUTS_ENV`] if set, else
/// the synthetic one written under the system temporary directory. Returns
/// the resolved input and whether it is the matched one.
pub(crate) fn problem_input(
    p: &ElectronProblem,
    histories: u64,
) -> Result<(ResolvedElectron, bool)> {
    if let Some(dir) = std::env::var_os(INPUTS_ENV) {
        let path = Path::new(&dir).join(p.id).join("input.toml");
        return Ok((resolve(&path)?, true));
    }
    let dir = std::env::temp_dir().join(format!("lindhard-bench-electron-{}", std::process::id()));
    let path = write_synthetic(&dir, p, &Settings::matched(p), histories)?;
    Ok((resolve(&path)?, false))
}

/// The elastic table of one material, as `lindhard-cli` builds it.
pub(crate) fn elastic_table(r: &ResolvedElectron, m: &Material) -> Result<CrossSectionTable> {
    let grid = &r.table_energy_ev;
    match r.elastic.potential {
        PotentialChoice::SalvatDhfs => {
            let opts = ElasticTableOptions {
                energy_ev: grid.clone(),
                ..ElasticTableOptions::default()
            };
            Ok(build_elastic_table(m, &SalvatDhfsTable, &opts)?)
        }
        PotentialChoice::ThomasFermiYukawa => {
            let desc = format!(
                "{}{}",
                ThomasFermiYukawa.description(),
                r.elastic.corrections_description()
            );
            let mut atoms = Vec::new();
            for c in m.components() {
                if c.atom_fraction() <= 0.0 {
                    continue;
                }
                let z = c.z();
                let y = ThomasFermiYukawa::yukawa(z)?;
                atoms.push(AtomicElastic::compute_corrected(
                    z,
                    &y,
                    &y,
                    &desc,
                    grid,
                    &r.elastic.corrections(z),
                    SolverOptions::default(),
                )?);
            }
            Ok(combine(
                m,
                &atoms,
                &default_probability_grid(),
                Some(DEFAULT_REFINE_TOLERANCE),
            )?)
        }
    }
}

/// The inelastic table of one material, as `lindhard-cli` builds it: on the
/// band-bottom axis with the band's minimum excitation energy as the model's
/// Fermi energy for a material with a band, on the model's own axis with
/// `[electron.inelastic] fermi_energy_ev` otherwise (#241).
pub(crate) fn inelastic_table(
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
) -> Result<CrossSectionTable> {
    let (fermi_ev, axis) = match &m.band {
        Some(b) => (b.min_excitation_ev(), EnergyAxis::BandBottom),
        None => (r.inelastic_fermi_ev, EnergyAxis::ModelFermiLevel),
    };
    let model = PennInelastic::try_new(r.inelastic, m.optical_elf.clone())?
        .with_fermi_energy_ev(fermi_ev)?;
    Ok(build_inelastic_table_for_model(
        &model,
        &m.material,
        &InelasticTableOptions::new(r.table_energy_ev.clone()).with_axis(axis),
    )?)
}

/// Both tables of every material, in `r.materials` order.
pub(crate) fn tables(r: &ResolvedElectron) -> Result<Vec<LayerTables>> {
    r.materials
        .iter()
        .map(|m| {
            Ok(LayerTables {
                elastic: elastic_table(r, &m.material)?,
                inelastic: inelastic_table(r, m)?,
            })
        })
        .collect()
}

/// The transport and a fresh full tally, set up as `lindhard-cli` does.
pub(crate) fn transport(
    r: &ResolvedElectron,
    tables: &[LayerTables],
) -> Result<(Transport, FullElectronTally)> {
    let layer_tables: Vec<LayerTables> = r
        .layer_material
        .iter()
        .map(|&i| tables[i].clone())
        .collect();
    let bands: Option<Vec<_>> = r
        .layer_material
        .iter()
        .map(|&i| r.materials[i].band.clone())
        .collect();
    let mut t = match bands {
        Some(b) => Transport::with_band_structures(r.stack.clone(), layer_tables, b, r.config),
        None => Transport::new(r.stack.clone(), layer_tables, r.config),
    }?;
    for (layer, &i) in r.layer_material.iter().enumerate() {
        let c = &r.materials[i].channels;
        if c.any() {
            t = t.with_insulator_channels(layer, c.clone())?;
        }
    }
    let proto = FullElectronTally::new(&t, r.tally)?;
    Ok((t, proto))
}

/// Run `histories` primaries on the current rayon pool and report.
pub(crate) fn run(
    r: &ResolvedElectron,
    t: &Transport,
    proto: &FullElectronTally,
    histories: u64,
) -> Result<ElectronReport> {
    let run = t.run(SEED, histories, CHUNK_SIZE, &r.primary, || proto.clone())?;
    Ok(run.tally.report())
}

/// SHA-256 of the report's JSON text: equal digests mean bit-identical
/// reports (`serde_json` writes every `f64` so that it reads back exactly).
pub(crate) fn report_digest(report: &ElectronReport) -> Result<String> {
    let text = serde_json::to_string(report)?;
    Ok(Sha256::digest(text.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect())
}
