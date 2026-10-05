//! Level-1 validation harness (`docs/validation.md`): the analytic and
//! internal checks of the scattering, stopping, BCA and tally work, gathered
//! into one target that prints one summary table.
//!
//! Runs under `cargo test` (it is a `harness = false` test target, so it
//! prints its table and exits non-zero if any check fails its tolerance):
//!
//! ```text
//! cargo test -p lindhard --test validation            # CI ("quick" statistics)
//! LINDHARD_VALIDATION=full cargo test --release -p lindhard --test validation
//! validation/run.sh                                   # full, and rewrite docs/validation.md
//! ```
//!
//! Environment:
//!
//! * `LINDHARD_VALIDATION=full`: more histories for the Monte Carlo checks
//!   (tolerances that contain a statistical term tighten with them).
//! * `LINDHARD_VALIDATION_OUT=<path>`: also write the table as Markdown.
//!
//! A first positional argument filters checks by substring of their id
//! (like libtest; every group still runs, only matching rows are reported and
//! asserted, and a filter that matches no id fails); `--list` lists the
//! groups. Rows are `pass`, `FAIL` or `info`
//! (reported, not asserted, with the reason in the note). Every tolerance is
//! the one achieved and justified in the PR that introduced the model, not a
//! band widened to make a check pass.

mod damage;
mod engine;
mod lss;
mod report;
mod scattering;
mod stopping;

use std::process::ExitCode;
use std::time::Instant;

use lindhard::input::TABLE_SPEC;
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::ScatteringTable;

use report::{Check, Status};

type Group = (&'static str, fn(&ScatteringTable, bool) -> Vec<Check>);

const GROUPS: [Group; 5] = [
    ("scatter", |t, _| scattering::checks(t)),
    ("stopping", |_, _| stopping::checks()),
    ("range", lss::checks),
    ("engine", engine::checks),
    ("damage", damage::checks),
];

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--list") {
        for (g, _) in GROUPS {
            println!("{g}: test");
        }
        return ExitCode::SUCCESS;
    }
    // libtest-style flags (--nocapture, --test-threads=..., --ignored) are
    // accepted and ignored; the first non-flag argument is a filter.
    let filter = args.iter().find(|a| !a.starts_with('-')).cloned();

    let full = std::env::var("LINDHARD_VALIDATION").is_ok_and(|v| v == "full");
    let t0 = Instant::now();
    // The engine's angle table, exactly as the CLI builds it (ZBL universal;
    // the angle is a function of the screening function only).
    let table = ScatteringTable::build(
        &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
        &TABLE_SPEC,
    );
    let mut rows = Vec::new();
    for (g, f) in GROUPS {
        let t = Instant::now();
        let mut r = f(&table, !full);
        if let Some(p) = filter.as_deref() {
            r.retain(|c| c.id.contains(p));
        }
        eprintln!(
            "validation: {g}: {} checks in {:.1} s",
            r.len(),
            t.elapsed().as_secs_f64()
        );
        rows.extend(r);
    }
    if rows.is_empty() {
        eprintln!(
            "validation: no check id matches the filter {:?}",
            filter.as_deref().unwrap_or("")
        );
        return ExitCode::FAILURE;
    }

    println!(
        "\nLevel-1 validation ({} statistics), {:.1} s\n",
        if full { "full" } else { "quick" },
        t0.elapsed().as_secs_f64()
    );
    print!("{}", report::text_table(&rows));
    let failed: Vec<_> = rows.iter().filter(|r| r.status == Status::Fail).collect();
    let passed = rows.iter().filter(|r| r.status == Status::Pass).count();
    let info = rows.len() - passed - failed.len();
    println!(
        "\n{passed} passed, {} failed, {info} reported (not asserted)",
        failed.len()
    );

    if let Ok(path) = std::env::var("LINDHARD_VALIDATION_OUT") {
        let md = format!(
            "Level 1, `{}` statistics, `lindhard` {}.\n\n{}",
            if full { "full" } else { "quick" },
            lindhard::VERSION,
            report::markdown_table(&rows)
        );
        if let Err(e) = std::fs::write(&path, md) {
            eprintln!("validation: cannot write {path}: {e}");
            return ExitCode::FAILURE;
        }
    }

    if failed.is_empty() {
        ExitCode::SUCCESS
    } else {
        for r in failed {
            eprintln!(
                "FAILED {}: {} = {} (tolerance {})",
                r.id, r.what, r.value, r.tolerance
            );
        }
        ExitCode::FAILURE
    }
}
