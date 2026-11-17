//! Validation figures for `docs/`, and the stopping-residual page of the book,
//! rendered from committed results.
//!
//! Every generated file is a pure function of committed inputs: no simulation
//! runs, nothing is read from the environment, and the output is text that is
//! byte-identical from run to run, so CI can re-render it and compare with the
//! committed copy (the `--check` mode of the two binaries), the way
//! `validation/update_docs.py --check` does for the tables.
//!
//! The two binaries differ only in their default paths, their usage and
//! regenerate lines, and the render call; everything else — the flags
//! ([`parse_args`]), writing or comparing the result ([`write_or_check`]) and
//! the exit code ([`exit_code`]) — is the shared driver below.
//!
//! Plotting is done with [rizzma](https://docs.rs/rizzma), with its default
//! features off; only its SVG backend is used. Text is drawn as glyph outlines
//! from the font rizzma embeds, so the output does not depend on fonts
//! installed on the machine.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};

pub mod backscatter;
pub mod stopping_report;

/// Outcome of comparing a freshly rendered figure with the file on disk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckOutcome {
    /// The file exists and matches the rendering byte for byte.
    Matches,
    /// The file does not exist.
    Missing,
    /// The file exists but differs; the first differing byte offset, and the
    /// two lengths.
    Differs {
        /// Byte offset of the first difference (the shorter length when one
        /// is a prefix of the other).
        first_difference: usize,
        /// Length of the committed file in bytes.
        committed_len: usize,
        /// Length of the fresh rendering in bytes.
        rendered_len: usize,
    },
}

/// Compare a rendering with the committed bytes (`None` when the file is
/// missing).
#[must_use]
pub fn compare(committed: Option<&[u8]>, rendered: &[u8]) -> CheckOutcome {
    let Some(committed) = committed else {
        return CheckOutcome::Missing;
    };
    if committed == rendered {
        return CheckOutcome::Matches;
    }
    let first_difference = committed
        .iter()
        .zip(rendered)
        .position(|(a, b)| a != b)
        .unwrap_or_else(|| committed.len().min(rendered.len()));
    CheckOutcome::Differs {
        first_difference,
        committed_len: committed.len(),
        rendered_len: rendered.len(),
    }
}

/// The workspace root: the parent of this crate's directory, so that the
/// default input and output paths work from anywhere in the checkout.
fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("the crate lives in a workspace subdirectory")
        .to_path_buf()
}

/// The command line of a generator binary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Args {
    /// Compare the committed file with a fresh rendering instead of writing
    /// it (`--check`).
    pub check: bool,
    /// What to render from: a results file or a dataset directory.
    pub input: PathBuf,
    /// The generated file to write, or to compare against.
    pub output: PathBuf,
    /// How the input is named in the generated file: the default relative
    /// path, or exactly what `--input` was given. A generator that records
    /// its input in its output uses this rather than [`Args::input`], whose
    /// default is an absolute path specific to the machine.
    pub input_label: String,
}

/// Parse `--check`, `--input`, `--output` and `-h`/`--help` from the process
/// arguments. The defaults are resolved against the workspace root (the
/// parent of this crate's directory); a path given on the command line is
/// used as it stands, so it is relative to the current directory.
///
/// `-h`/`--help` prints `usage` and exits 0.
///
/// # Errors
///
/// An unexpected argument, or `--input`/`--output` with no value after it;
/// either error quotes `usage`.
pub fn parse_args(usage: &'static str, default_input: &str, default_output: &str) -> Result<Args> {
    let root = workspace_root();
    let mut args = Args {
        check: false,
        input: root.join(default_input),
        output: root.join(default_output),
        input_label: default_input.to_string(),
    };
    let mut it = std::env::args_os().skip(1);
    while let Some(arg) = it.next() {
        match arg.to_str() {
            Some("--check") => args.check = true,
            Some("--input") => {
                args.input = it.next().context(usage)?.into();
                args.input_label = args.input.display().to_string();
            }
            Some("--output") => args.output = it.next().context(usage)?.into(),
            Some("-h" | "--help") => {
                println!("{usage}");
                std::process::exit(0);
            }
            _ => bail!("unexpected argument {arg:?}\n{usage}"),
        }
    }
    Ok(args)
}

/// Write `rendered` to `args.output`, or, under `--check`, compare it with the
/// committed file and report what differs.
///
/// `regenerate_command` is the command that rewrites the file, named in both
/// failure messages. Returns whether the file is now up to date: `true` when
/// it was written or already matched, `false` when it is missing or stale
/// (which [`exit_code`] turns into exit 1).
///
/// # Errors
///
/// The output file cannot be written, or cannot be read for a reason other
/// than not existing.
pub fn write_or_check(args: &Args, rendered: &[u8], regenerate_command: &str) -> Result<bool> {
    let output = &args.output;
    if !args.check {
        std::fs::write(output, rendered)
            .with_context(|| format!("writing {}", output.display()))?;
        println!("wrote {} ({} bytes)", output.display(), rendered.len());
        return Ok(true);
    }

    let committed = match std::fs::read(output) {
        Ok(bytes) => Some(bytes),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => None,
        Err(e) => return Err(e).with_context(|| format!("reading {}", output.display())),
    };
    match compare(committed.as_deref(), rendered) {
        CheckOutcome::Matches => {
            println!("{} is up to date", output.display());
            Ok(true)
        }
        CheckOutcome::Missing => {
            eprintln!(
                "{} is missing; regenerate it with `{regenerate_command}` and commit it",
                output.display()
            );
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
                 bytes); regenerate it with `{regenerate_command}` and commit it",
                output.display(),
                args.input.display()
            );
            Ok(false)
        }
    }
}

/// The exit code of a generator binary: 0 written or up to date, 1 missing or
/// stale, 2 error (printed to stderr).
#[must_use]
pub fn exit_code(outcome: Result<bool>) -> ExitCode {
    match outcome {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::from(2)
        }
    }
}
