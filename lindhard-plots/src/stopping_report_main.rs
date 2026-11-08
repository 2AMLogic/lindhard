//! `lindhard-stopping-report`: write the stopping-residual page of the book, or
//! check it.
//!
//! ```text
//! cargo run -p lindhard-plots --bin lindhard-stopping-report             # write book/src/models/stopping-residuals.md
//! cargo run -p lindhard-plots --bin lindhard-stopping-report -- --check  # exit 1 if the committed page is stale or missing
//! ```
//!
//! `--input` (a directory of datasets) and `--output` override the default
//! paths; the defaults are resolved against the workspace root.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use lindhard_plots::{compare, stopping_report, CheckOutcome};

const DEFAULT_INPUT: &str = "validation/data/stopping";
const DEFAULT_OUTPUT: &str = "book/src/models/stopping-residuals.md";
const USAGE: &str =
    "usage: lindhard-stopping-report [--check] [--input <dataset dir>] [--output <page.md>]";

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("the crate lives in a workspace subdirectory")
        .to_path_buf()
}

fn run() -> Result<bool> {
    let root = workspace_root();
    let (mut check, mut input, mut output) =
        (false, root.join(DEFAULT_INPUT), root.join(DEFAULT_OUTPUT));
    let mut label = DEFAULT_INPUT.to_string();
    let mut it = std::env::args_os().skip(1);
    while let Some(arg) = it.next() {
        match arg.to_str() {
            Some("--check") => check = true,
            Some("--input") => {
                input = it.next().context(USAGE)?.into();
                label = input.display().to_string();
            }
            Some("--output") => output = it.next().context(USAGE)?.into(),
            Some("-h" | "--help") => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            _ => bail!("unexpected argument {arg:?}\n{USAGE}"),
        }
    }
    let page = stopping_report::render(&input, &label)?;

    if !check {
        std::fs::write(&output, &page).with_context(|| format!("writing {}", output.display()))?;
        println!("wrote {} ({} bytes)", output.display(), page.len());
        return Ok(true);
    }
    let committed = match std::fs::read(&output) {
        Ok(bytes) => Some(bytes),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => None,
        Err(e) => return Err(e).with_context(|| format!("reading {}", output.display())),
    };
    let regenerate = "regenerate it with `cargo run -p lindhard-plots --bin \
                      lindhard-stopping-report` and commit it";
    match compare(committed.as_deref(), page.as_bytes()) {
        CheckOutcome::Matches => {
            println!("{} is up to date", output.display());
            Ok(true)
        }
        CheckOutcome::Missing => {
            eprintln!("{} is missing; {regenerate}", output.display());
            Ok(false)
        }
        CheckOutcome::Differs {
            first_difference,
            committed_len,
            rendered_len,
        } => {
            eprintln!(
                "{} differs from a fresh render of {} (first difference at byte \
                 {first_difference}; committed {committed_len} bytes, rendered {rendered_len} \
                 bytes); {regenerate}",
                output.display(),
                input.display()
            );
            Ok(false)
        }
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::from(2)
        }
    }
}
