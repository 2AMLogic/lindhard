//! Helpers shared by the CLI integration tests. Cargo does not build this
//! directory as its own test target; each test binary pulls it in with
//! `mod common;`.
//!
//! The `env!` lookups below are expanded in whichever test binary includes
//! this module. `CARGO_TARGET_TMPDIR` is one workspace-wide directory, so
//! [`scratch`] also joins `CARGO_CRATE_NAME` and every test binary's
//! scratch directories stay separate under parallel runners;
//! `CARGO_BIN_EXE_lindhard` is the `lindhard` binary Cargo built for it.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::Command;

use lindhard::input::Input;
use lindhard::ion::bca::{Bca, CrystalTarget};
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;
use lindhard::tally::IonTally;
use lindhard_cli::output;
use lindhard_cli::output::RunInfo;
use lindhard_cli::sim::Simulation;
use lindhard_cli::tally::CliTally;
use serde_json::Value;

/// A fresh scratch directory `name` under this test binary's corner of the
/// Cargo-managed temporary directory, emptied first so a rerun does not see
/// the last run's files.
pub(crate) fn scratch(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("lindhard-cli-tests")
        .join(env!("CARGO_CRATE_NAME"))
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// `text` parsed and resolved through the shared input path.
pub(crate) fn resolve(text: &str) -> lindhard::input::Resolved {
    Input::from_toml_str(text).unwrap().resolve().unwrap()
}

/// The silicon `(001)` cut with the `[010]` reference direction, at the
/// 7 deg tilt and 22 deg azimuth of the CLI test inputs.
pub(crate) fn si_001_target() -> CrystalTarget {
    let lattice = Lattice::silicon();
    let orientation = Orientation::new(
        &lattice,
        [0, 0, 1],
        [0, 1, 0],
        7f64.to_radians(),
        22f64.to_radians(),
        0.0,
    )
    .unwrap();
    CrystalTarget::new(lattice, orientation)
}

/// The engine set up by hand, mirroring `simulate`: the ZBL scattering
/// table of `r`, its stopping model, a `Bca` customized by `setup`
/// (divergence, crystal partners; `|bca| bca` for the plain amorphous
/// engine), the ion and CLI tallies, and a rayon pool of `threads`
/// workers. Timings are zeroed, so only `threads` distinguishes the
/// result from a CLI run.
pub(crate) fn direct_simulation(
    r: &lindhard::input::Resolved,
    threads: usize,
    setup: impl FnOnce(Bca<'_>) -> Bca<'_>,
) -> Simulation {
    let pot = Potential::new(r.screening, f64::from(r.beam.ion.z()), 14.0)
        .with_length(r.screening_length);
    let table = ScatteringTable::build(&pot, &r.table_spec);
    let stopping = r.stopping_model();
    let bca = setup(Bca::new(r.beam, r.geometry(), r.config, &*stopping, &table).unwrap());
    let crystals = bca.crystal_metadata();
    let proto = IonTally::new(
        r.geometry(),
        &bca.species_z(),
        r.ion_tally_config().unwrap(),
    )
    .unwrap();
    let t = &r.input.tally;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let tally = pool
        .install(|| {
            bca.run(|| {
                CliTally::new(
                    t.depth_bin_nm * 1e-9,
                    t.depth_bins,
                    t.per_ion,
                    proto.clone(),
                )
            })
        })
        .unwrap();
    let report = tally.ion.report(t.dual_pearson);
    Simulation {
        tally,
        table,
        report,
        crystals,
        info: RunInfo {
            threads,
            table_build_s: 0.0,
            transport_s: 0.0,
            ions_per_s: 0.0,
        },
    }
}

/// `lindhard check` on `text`, written to `scratch(name)/in.toml`:
/// (success, stdout followed by stderr).
pub(crate) fn check(name: &str, text: &str) -> (bool, String) {
    let p = scratch(name).join("in.toml");
    std::fs::write(&p, text).unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_lindhard"))
        .args(["check", p.to_str().unwrap()])
        .output()
        .unwrap();
    (
        o.status.success(),
        format!(
            "{}{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        ),
    )
}

/// The summary JSON of `s` with the `run` block dropped: that block holds
/// timings, which no two runs share.
pub(crate) fn summary(r: &lindhard::input::Resolved, s: &Simulation) -> Value {
    let text = output::summary_json(r, &s.table, &s.tally, &s.report, &s.crystals, s.info).unwrap();
    let mut v: Value = serde_json::from_str(&text).unwrap();
    v.as_object_mut().unwrap().remove("run");
    v
}
