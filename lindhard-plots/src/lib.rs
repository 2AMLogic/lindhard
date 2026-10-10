//! Validation figures for `docs/`, rendered from committed results.
//!
//! Every figure is a pure function of a committed results file: no simulation
//! runs, nothing is read from the environment, and the output is SVG text that
//! is byte-identical from run to run, so CI can re-render it and compare with
//! the committed copy (the `--check` mode of the `lindhard-plots` binary), the
//! way `validation/update_docs.py --check` does for the tables.
//!
//! Plotting is done with [rizzma](https://docs.rs/rizzma), with its default
//! features off; only its SVG backend is used. Text is drawn as glyph outlines
//! from the font rizzma embeds, so the output does not depend on fonts
//! installed on the machine.

pub mod backscatter;

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
