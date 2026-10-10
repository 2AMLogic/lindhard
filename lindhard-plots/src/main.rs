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

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use lindhard_plots::{backscatter, compare, CheckOutcome};

const DEFAULT_INPUT: &str = "validation/experiments/backscatter_results.json";
const DEFAULT_OUTPUT: &str = "docs/img/backscatter-eta.svg";
const USAGE: &str =
    "usage: lindhard-plots [--check] [--input <results.json>] [--output <figure.svg>]";

#[derive(Debug)]
struct Args {
    check: bool,
    input: PathBuf,
    output: PathBuf,
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("the crate lives in a workspace subdirectory")
        .to_path_buf()
}

fn parse_args() -> Result<Args> {
    let root = workspace_root();
    let mut args = Args {
        check: false,
        input: root.join(DEFAULT_INPUT),
        output: root.join(DEFAULT_OUTPUT),
    };
    let mut it = std::env::args_os().skip(1);
    while let Some(arg) = it.next() {
        match arg.to_str() {
            Some("--check") => args.check = true,
            Some("--input") => args.input = it.next().context(USAGE)?.into(),
            Some("--output") => args.output = it.next().context(USAGE)?.into(),
            Some("-h" | "--help") => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            _ => bail!("unexpected argument {arg:?}\n{USAGE}"),
        }
    }
    Ok(args)
}

fn run() -> Result<bool> {
    let args = parse_args()?;
    let json = std::fs::read_to_string(&args.input)
        .with_context(|| format!("reading {}", args.input.display()))?;
    let results =
        backscatter::parse(&json).with_context(|| format!("in {}", args.input.display()))?;
    let svg = backscatter::render_svg(&results);

    if !args.check {
        std::fs::write(&args.output, &svg)
            .with_context(|| format!("writing {}", args.output.display()))?;
        println!("wrote {} ({} bytes)", args.output.display(), svg.len());
        return Ok(true);
    }

    let committed = match std::fs::read(&args.output) {
        Ok(bytes) => Some(bytes),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => None,
        Err(e) => return Err(e).with_context(|| format!("reading {}", args.output.display())),
    };
    match compare(committed.as_deref(), svg.as_bytes()) {
        CheckOutcome::Matches => {
            println!("{} is up to date", args.output.display());
            Ok(true)
        }
        CheckOutcome::Missing => {
            eprintln!(
                "{} is missing; regenerate it with `cargo run -p lindhard-plots` and commit it",
                args.output.display()
            );
            Ok(false)
        }
        CheckOutcome::Differs {
            first_difference,
            committed_len,
            rendered_len,
        } => {
            eprintln!(
                "{} differs from a fresh render of {} (first difference at byte {first_difference}; \
                 committed {committed_len} bytes, rendered {rendered_len} bytes); regenerate it with \
                 `cargo run -p lindhard-plots` and commit it",
                args.output.display(),
                args.input.display()
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
