//! `lindhard-plots`: render the committed validation figures, or check them.
//!
//! ```text
//! cargo run -p lindhard-plots                # write docs/img/backscatter-eta.svg
//! cargo run -p lindhard-plots -- --check     # exit 1 if the committed SVG is stale or missing
//! ```
//!
//! `--input` and `--output` override the default paths (both relative to the
//! current directory when given; the defaults are resolved against the
//! workspace root, so the command works from anywhere in the checkout).
//!
//! The flags, the write-or-compare step and the exit codes are the shared
//! driver in `lindhard_plots`; only the paths, the usage and regenerate lines
//! and the render call are here.

use std::process::ExitCode;

use anyhow::{Context, Result};
use lindhard_plots::{backscatter, exit_code, parse_args, write_or_check, Args};

const DEFAULT_INPUT: &str = "validation/experiments/backscatter_results.json";
const DEFAULT_OUTPUT: &str = "docs/img/backscatter-eta.svg";
const USAGE: &str =
    "usage: lindhard-plots [--check] [--input <results.json>] [--output <figure.svg>]";
const REGENERATE: &str = "cargo run -p lindhard-plots";

fn render(args: &Args) -> Result<String> {
    let json = std::fs::read_to_string(&args.input)
        .with_context(|| format!("reading {}", args.input.display()))?;
    let results =
        backscatter::parse(&json).with_context(|| format!("in {}", args.input.display()))?;
    Ok(backscatter::render_svg(&results))
}

fn run() -> Result<bool> {
    let args = parse_args(USAGE, DEFAULT_INPUT, DEFAULT_OUTPUT)?;
    let svg = render(&args)?;
    write_or_check(&args, svg.as_bytes(), REGENERATE)
}

fn main() -> ExitCode {
    exit_code(run())
}
