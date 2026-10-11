//! Data the electron engine consumes, as validated types.
//!
//! Three kinds of data live here, each with a mandatory provenance string:
//!
//! - [`OpticalElf`]: a tabulated optical energy-loss function
//!   `ELF(ħω) = Im[-1/ε(ω)]` of one material at zero momentum transfer, the
//!   input of dielectric-function inelastic models (Ritchie, Phys. Rev. 106,
//!   874 (1957); Penn, Phys. Rev. B 35, 482 (1987)).
//! - [`SubshellBindingTable`]: atomic subshell binding energies and ground-state
//!   occupancies, for inner-shell ionization. It can be read from the
//!   ENDF-6 File 28 (atomic relaxation) format in which EADL is distributed
//!   ([`SubshellBindingTable::from_endf6_mf28_str`]).
//! - [`CrossSectionTable`]: a precomputed cross-section cache on an
//!   incident-energy grid: the inverse mean free path and the inverse CDF of
//!   one sampled variable (elastic polar angle or inelastic energy loss).
//!   [`ShellChannelTable`] wraps an energy-loss table of one inner-shell
//!   ionisation channel with the shell's identity, binding energy and the
//!   provenance of that energy.
//!
//! # Validation
//!
//! Every type has private fields and a checked constructor. The serde
//! implementations go through the same checks (`#[serde(try_from = ...)]`), so
//! a value read from disk by any serde format is valid by construction; a
//! derived deserializer cannot bypass the invariants.
//!
//! A provenance (and every identity string: material, model) must be
//! non-blank after trimming and must not contain control characters (it is a
//! one-line citation that is copied into run metadata). Missing or blank
//! provenance is [`ElectronDataError::MissingProvenance`], as for user stopping
//! tables (`ion::stopping::table`): data without an origin is not admitted
//! (`docs/data-provenance.md`).
//!
//! # Units
//!
//! SI internally (`crate::units`), with the unit spelled in each name at the
//! API boundary: energies in eV (`_ev`), inverse mean free paths in m⁻¹
//! (`_per_m`), angles in radians.
//!
//! The EADL2017 binding-energy table is committed with this module
//! ([`SubshellBindingTable::eadl2017`]). No optical data is; see
//! `docs/data-provenance.md` for the sources and for which redistribution
//! terms are still open.

mod binding;
mod endf6;
mod optical;
mod table;

pub use binding::{
    AtomBindings, ShellBinding, Subshell, SubshellBindingTable, OCCUPANCY_SUM_TOLERANCE,
};
pub use optical::OpticalElf;
pub use table::{
    CrossSectionTable, CrossSectionTableParts, SamplingAxis, ShellChannelTable,
    SHELL_CHANNEL_FORMAT_VERSION,
};

use std::path::Path;

/// The cross-section cache format this build reads and writes.
///
/// Bump it whenever the meaning or layout of a [`CrossSectionTable`] file
/// changes; a reader rejects every other version with
/// [`ElectronDataError::UnsupportedVersion`] instead of guessing.
pub const CACHE_FORMAT_VERSION: u32 = 1;

/// Errors from building, reading or querying electron data.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ElectronDataError {
    /// No provenance, or one that is blank after trimming.
    #[error("data has no provenance; a citation or description of its origin is required")]
    MissingProvenance,
    /// A required identity string (material, model) is missing or blank.
    #[error("data has no {0}; it is required")]
    MissingIdentity(&'static str),
    /// A text field is present but unusable.
    #[error("malformed {field}: {reason}")]
    MalformedText {
        /// Field name.
        field: &'static str,
        /// What is wrong with it.
        reason: String,
    },
    /// A numeric field or the shape of the data violates an invariant.
    #[error("invalid {what}: {reason}")]
    Invalid {
        /// Which quantity.
        what: &'static str,
        /// What is wrong with it.
        reason: String,
    },
    /// A query outside the tabulated range. Nothing is extrapolated.
    #[error("{quantity} {value} outside table range [{min}, {max}]")]
    OutOfRange {
        /// Which quantity was queried.
        quantity: &'static str,
        /// Queried value.
        value: f64,
        /// Table minimum.
        min: f64,
        /// Table maximum.
        max: f64,
    },
    /// A grid index past the end of the table.
    #[error("grid index {index} out of range for a table of {len} points")]
    IndexOutOfRange {
        /// Requested index.
        index: usize,
        /// Number of grid points.
        len: usize,
    },
    /// A cache file declares a format version this build does not read.
    #[error(
        "unsupported cross-section cache format version {}; this build reads version {supported}",
        found.map_or_else(|| "(none)".to_string(), |v| v.to_string())
    )]
    UnsupportedVersion {
        /// Version found in the file (`None` if absent).
        found: Option<u32>,
        /// Version this build supports ([`CACHE_FORMAT_VERSION`]).
        supported: u32,
    },
    /// Sampling was requested at an energy whose inverse mean free path is
    /// zero: there is no interaction there and no distribution to sample.
    #[error(
        "no interaction at {energy_ev} eV: the inverse mean free path is zero, nothing to sample"
    )]
    ZeroRate {
        /// Grid energy, eV.
        energy_ev: f64,
    },
    /// A file could not be parsed.
    #[error("parse error: {0}")]
    Parse(String),
    /// A file could not be read or written.
    #[error("I/O error: {0}")]
    Io(String),
}

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Check a required one-line text field: trimmed, non-blank, no control
/// characters. A missing or blank `provenance` is
/// [`ElectronDataError::MissingProvenance`]; any other missing field is
/// [`ElectronDataError::MissingIdentity`].
fn checked_text(field: &'static str, value: Option<&str>) -> Result<String> {
    let text = value.unwrap_or_default().trim();
    if text.is_empty() {
        return Err(if field == "provenance" {
            ElectronDataError::MissingProvenance
        } else {
            ElectronDataError::MissingIdentity(field)
        });
    }
    if let Some(c) = text.chars().find(|c| c.is_control()) {
        return Err(ElectronDataError::MalformedText {
            field,
            reason: format!("contains the control character {c:?}; use one line of text"),
        });
    }
    Ok(text.to_string())
}

fn invalid(what: &'static str, reason: impl Into<String>) -> ElectronDataError {
    ElectronDataError::Invalid {
        what,
        reason: reason.into(),
    }
}

/// An energy grid: at least two points, finite, positive, strictly increasing.
fn check_energy_grid(what: &'static str, grid: &[f64]) -> Result<()> {
    if grid.len() < 2 {
        return Err(invalid(what, "needs at least 2 points"));
    }
    if let Some(x) = grid.iter().find(|x| !(x.is_finite() && **x > 0.0)) {
        return Err(invalid(
            what,
            format!("values must be finite and positive, got {x}"),
        ));
    }
    if let Some(w) = grid.windows(2).find(|w| w[1] <= w[0]) {
        return Err(invalid(
            what,
            format!(
                "must be strictly increasing (no duplicates), got {} then {}",
                w[0], w[1]
            ),
        ));
    }
    Ok(())
}

/// Finite and non-negative values (zero allowed).
fn check_nonnegative(what: &'static str, values: &[f64]) -> Result<()> {
    if let Some(x) = values.iter().find(|x| !(x.is_finite() && **x >= 0.0)) {
        return Err(invalid(
            what,
            format!("values must be finite and non-negative, got {x}"),
        ));
    }
    Ok(())
}

/// Piecewise-linear interpolation on a validated, strictly increasing grid,
/// written as a convex combination so that knots are reproduced exactly and
/// non-negative ordinates give non-negative results. `x` must be in range.
fn lerp_on_grid(grid: &[f64], values: &[f64], x: f64) -> f64 {
    let n = grid.len();
    let i = grid.partition_point(|&g| g <= x).clamp(1, n - 1);
    let (x0, x1) = (grid[i - 1], grid[i]);
    let t = (x - x0) / (x1 - x0);
    values[i - 1] * (1.0 - t) + values[i] * t
}

fn check_in_range(quantity: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if !value.is_finite() || value < min || value > max {
        return Err(ElectronDataError::OutOfRange {
            quantity,
            value,
            min,
            max,
        });
    }
    Ok(())
}

fn parse_error(e: impl std::fmt::Display) -> ElectronDataError {
    ElectronDataError::Parse(e.to_string())
}

fn read_text(path: &Path) -> Result<String> {
    std::fs::read_to_string(path)
        .map_err(|e| ElectronDataError::Io(format!("{}: {e}", path.display())))
}
