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
//!
//! The flags, the write-or-compare step and the exit codes are the shared
//! driver in `lindhard_plots`; only the paths, the usage and regenerate lines
//! and the render call are here. The page records where its datasets came
//! from, so the render call takes the input label the driver parsed.

use std::process::ExitCode;

use anyhow::Result;
use lindhard_plots::{exit_code, parse_args, stopping_report, write_or_check};

const DEFAULT_INPUT: &str = "validation/data/stopping";
const DEFAULT_OUTPUT: &str = "book/src/models/stopping-residuals.md";
const USAGE: &str =
    "usage: lindhard-stopping-report [--check] [--input <dataset dir>] [--output <page.md>]";
const REGENERATE: &str = "cargo run -p lindhard-plots --bin lindhard-stopping-report";

fn run() -> Result<bool> {
    let args = parse_args(USAGE, DEFAULT_INPUT, DEFAULT_OUTPUT)?;
    let page = stopping_report::render(&args.input, &args.input_label)?;
    write_or_check(&args, page.as_bytes(), REGENERATE)
}

fn main() -> ExitCode {
    exit_code(run())
}
