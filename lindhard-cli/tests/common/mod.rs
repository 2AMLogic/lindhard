//! Helpers shared by the CLI integration tests. Cargo does not build this
//! directory as its own test target; each test binary pulls it in with
//! `mod common;`.
//!
//! The `env!` lookups below are expanded in whichever test binary includes
//! this module, so each one gets its own `CARGO_TARGET_TMPDIR` and the
//! `lindhard` binary Cargo built for it.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::Command;

use lindhard_cli::output;
use lindhard_cli::sim::Simulation;
use serde_json::Value;

/// A fresh scratch directory `name` under the test target's temporary
/// directory, emptied first so a rerun does not see the last run's files.
pub(crate) fn scratch(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("lindhard-cli-tests")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
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
