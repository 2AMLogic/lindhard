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

use serde::de::IgnoredAny;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
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

/// The committed EADL2017 binding-energy table (see
/// [`SubshellBindingTable::eadl2017`]); its header carries the credit.
const EADL2017_TOML: &str = include_str!("eadl2017_binding.toml");

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

// ---------------------------------------------------------------------------
// Optical energy-loss function
// ---------------------------------------------------------------------------

/// A tabulated optical energy-loss function `ELF(ħω) = Im[-1/ε(ω)]` of one
/// material, at zero momentum transfer.
///
/// The ELF is the quantity a dielectric-function inelastic model extends into
/// momentum transfer (Ritchie, Phys. Rev. 106, 874 (1957); Penn, Phys. Rev. B
/// 35, 482 (1987)). It is dimensionless.
///
/// Invariants (checked by [`OpticalElf::new`] and by every loader):
/// non-blank `material` and `provenance`; at least two samples; photon
/// energies `ħω` finite, positive and strictly increasing; ELF values finite
/// and non-negative (zero is allowed, e.g. inside a band gap), one per energy.
///
/// **Interpolation** is piecewise linear in `(ħω, ELF)`. It reproduces the
/// tabulated values at the knots, never goes negative, and handles zero-valued
/// ELF (where log-log interpolation would fail). Queries outside
/// `[ħω_min, ħω_max]` are an [`ElectronDataError::OutOfRange`] error: nothing
/// is extrapolated here. Extrapolation to high energy (for example by a sum
/// rule or atomic photoabsorption data) belongs to a later physical model.
///
/// TOML form (all four keys required, no others accepted):
///
/// ```toml
/// material = "Si"
/// provenance = "Author, Journal vol, page (year), Table N"
/// energy_ev = [1.0, 10.0, 16.7, 100.0]
/// elf = [0.0, 0.5, 4.0, 0.1]
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawOpticalElf", into = "RawOpticalElf")]
pub struct OpticalElf {
    material: String,
    provenance: String,
    energy_ev: Vec<f64>,
    elf: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawOpticalElf {
    material: Option<String>,
    provenance: Option<String>,
    energy_ev: Vec<f64>,
    elf: Vec<f64>,
}

impl TryFrom<RawOpticalElf> for OpticalElf {
    type Error = ElectronDataError;
    fn try_from(raw: RawOpticalElf) -> Result<Self> {
        let provenance = checked_text("provenance", raw.provenance.as_deref())?;
        let material = checked_text("material", raw.material.as_deref())?;
        Self::new(material, provenance, raw.energy_ev, raw.elf)
    }
}

impl From<OpticalElf> for RawOpticalElf {
    fn from(e: OpticalElf) -> Self {
        Self {
            material: Some(e.material),
            provenance: Some(e.provenance),
            energy_ev: e.energy_ev,
            elf: e.elf,
        }
    }
}

impl OpticalElf {
    /// Build a validated ELF table. See the type docs for the invariants.
    pub fn new(
        material: impl Into<String>,
        provenance: impl Into<String>,
        energy_ev: Vec<f64>,
        elf: Vec<f64>,
    ) -> Result<Self> {
        let provenance = checked_text("provenance", Some(&provenance.into()))?;
        let material = checked_text("material", Some(&material.into()))?;
        check_energy_grid("photon energy grid", &energy_ev)?;
        if elf.len() != energy_ev.len() {
            return Err(invalid(
                "ELF values",
                format!(
                    "need one value per energy: {} energies, {} values",
                    energy_ev.len(),
                    elf.len()
                ),
            ));
        }
        check_nonnegative("ELF values", &elf)?;
        Ok(Self {
            material,
            provenance,
            energy_ev,
            elf,
        })
    }

    /// Parse and validate the TOML form (see the type docs).
    pub fn from_toml_str(text: &str) -> Result<Self> {
        let raw: RawOpticalElf = toml::from_str(text).map_err(parse_error)?;
        Self::try_from(raw)
    }

    /// Load the TOML form from a file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_toml_str(&read_text(path.as_ref())?)
    }

    /// Material identity, as declared by the data.
    pub fn material(&self) -> &str {
        &self.material
    }

    /// Provenance string (to be copied into run metadata).
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Tabulated photon energies `ħω`, eV.
    pub fn energy_ev(&self) -> &[f64] {
        &self.energy_ev
    }

    /// Tabulated ELF values (dimensionless), one per energy.
    pub fn elf_values(&self) -> &[f64] {
        &self.elf
    }

    /// The tabulated range `(min, max)` of `ħω`, eV.
    pub fn energy_range_ev(&self) -> (f64, f64) {
        (self.energy_ev[0], self.energy_ev[self.energy_ev.len() - 1])
    }

    /// `ELF(ħω)` by piecewise-linear interpolation; an error outside the
    /// tabulated range.
    pub fn elf(&self, energy_ev: f64) -> Result<f64> {
        let (lo, hi) = self.energy_range_ev();
        check_in_range("photon energy (eV)", energy_ev, lo, hi)?;
        Ok(lerp_on_grid(&self.energy_ev, &self.elf, energy_ev))
    }
}

// ---------------------------------------------------------------------------
// Subshell binding energies
// ---------------------------------------------------------------------------

/// Subshell labels and orbitals by ENDF subshell designator 1..=39, from the
/// ENDF-6 Formats Manual (Trkov, Herman, Brown (eds.), CSEWG Document
/// ENDF-102, report BNL-203218-2018-INRE (2018)), section 28.2 and Appendix
/// B.1: designator `n` is the subshell whose photo- or electro-atomic
/// cross section has `MT = 533 + n` (MT 534 = K ... MT 572 = Q3). The second
/// field is `2j`, so the subshell holds at most `2j + 1` electrons.
const SUBSHELLS: [(&str, &str, u8); 39] = [
    ("K", "1s1/2", 1),
    ("L1", "2s1/2", 1),
    ("L2", "2p1/2", 1),
    ("L3", "2p3/2", 3),
    ("M1", "3s1/2", 1),
    ("M2", "3p1/2", 1),
    ("M3", "3p3/2", 3),
    ("M4", "3d3/2", 3),
    ("M5", "3d5/2", 5),
    ("N1", "4s1/2", 1),
    ("N2", "4p1/2", 1),
    ("N3", "4p3/2", 3),
    ("N4", "4d3/2", 3),
    ("N5", "4d5/2", 5),
    ("N6", "4f5/2", 5),
    ("N7", "4f7/2", 7),
    ("O1", "5s1/2", 1),
    ("O2", "5p1/2", 1),
    ("O3", "5p3/2", 3),
    ("O4", "5d3/2", 3),
    ("O5", "5d5/2", 5),
    ("O6", "5f5/2", 5),
    ("O7", "5f7/2", 7),
    ("O8", "5g7/2", 7),
    ("O9", "5g9/2", 9),
    ("P1", "6s1/2", 1),
    ("P2", "6p1/2", 1),
    ("P3", "6p3/2", 3),
    ("P4", "6d3/2", 3),
    ("P5", "6d5/2", 5),
    ("P6", "6f5/2", 5),
    ("P7", "6f7/2", 7),
    ("P8", "6g7/2", 7),
    ("P9", "6g9/2", 9),
    ("P10", "6h9/2", 9),
    ("P11", "6h11/2", 11),
    ("Q1", "7s1/2", 1),
    ("Q2", "7p1/2", 1),
    ("Q3", "7p3/2", 3),
];

/// One atomic subshell, identified by its ENDF subshell designator
/// (1 = K (1s1/2), 2 = L1 (2s1/2), 3 = L2 (2p1/2), 4 = L3 (2p3/2), ...,
/// 39 = Q3 (7p3/2); see [`Subshell::label`] and [`Subshell::orbital`]).
///
/// Serialized as the bare designator number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "u8", into = "u8")]
pub struct Subshell(u8);

impl Subshell {
    /// The K shell (designator 1).
    pub const K: Subshell = Subshell(1);

    /// Largest designator defined by the ENDF-6 manual (Q3).
    pub const MAX_DESIGNATOR: u8 = 39;

    /// A subshell from its ENDF designator, 1..=39.
    pub fn new(designator: u8) -> Result<Self> {
        if (1..=Self::MAX_DESIGNATOR).contains(&designator) {
            Ok(Self(designator))
        } else {
            Err(invalid(
                "subshell designator",
                format!("{designator} is not in 1..={}", Self::MAX_DESIGNATOR),
            ))
        }
    }

    /// Look a subshell up by its label (`"K"`, `"L1"`, ..., `"Q3"`).
    pub fn from_label(label: &str) -> Option<Self> {
        SUBSHELLS
            .iter()
            .position(|s| s.0 == label)
            .map(|i| Self(i as u8 + 1))
    }

    /// The ENDF designator.
    pub fn designator(self) -> u8 {
        self.0
    }

    /// X-ray notation label, e.g. `"L3"`.
    pub fn label(self) -> &'static str {
        SUBSHELLS[usize::from(self.0) - 1].0
    }

    /// Orbital, e.g. `"2p3/2"`.
    pub fn orbital(self) -> &'static str {
        SUBSHELLS[usize::from(self.0) - 1].1
    }

    /// Maximum number of electrons, `2j + 1`.
    pub fn capacity(self) -> u8 {
        SUBSHELLS[usize::from(self.0) - 1].2 + 1
    }
}

impl TryFrom<u8> for Subshell {
    type Error = ElectronDataError;
    fn try_from(d: u8) -> Result<Self> {
        Self::new(d)
    }
}

impl From<Subshell> for u8 {
    fn from(s: Subshell) -> u8 {
        s.0
    }
}

/// The binding energy and ground-state occupancy of one occupied subshell.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawShellBinding", into = "RawShellBinding")]
pub struct ShellBinding {
    subshell: Subshell,
    binding_energy_ev: f64,
    occupancy: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawShellBinding {
    subshell: Subshell,
    binding_energy_ev: f64,
    occupancy: f64,
}

impl TryFrom<RawShellBinding> for ShellBinding {
    type Error = ElectronDataError;
    fn try_from(r: RawShellBinding) -> Result<Self> {
        Self::new(r.subshell, r.binding_energy_ev, r.occupancy)
    }
}

impl From<ShellBinding> for RawShellBinding {
    fn from(s: ShellBinding) -> Self {
        Self {
            subshell: s.subshell,
            binding_energy_ev: s.binding_energy_ev,
            occupancy: s.occupancy,
        }
    }
}

impl ShellBinding {
    /// An occupied subshell: binding energy finite and positive (eV), and
    /// occupancy finite, positive and at most the subshell capacity `2j + 1`.
    /// Occupancies may be fractional: EADL spreads the electrons of an open
    /// shell over its `j`-split subshells in proportion to their capacity.
    pub fn new(subshell: Subshell, binding_energy_ev: f64, occupancy: f64) -> Result<Self> {
        if !(binding_energy_ev.is_finite() && binding_energy_ev > 0.0) {
            return Err(invalid(
                "binding energy",
                format!(
                    "{} binding energy must be finite and positive, got {binding_energy_ev} eV",
                    subshell.label()
                ),
            ));
        }
        let cap = f64::from(subshell.capacity());
        if !(occupancy.is_finite() && occupancy > 0.0 && occupancy <= cap) {
            return Err(invalid(
                "occupancy",
                format!(
                    "{} occupancy must be in (0, {cap}], got {occupancy}",
                    subshell.label()
                ),
            ));
        }
        Ok(Self {
            subshell,
            binding_energy_ev,
            occupancy,
        })
    }

    /// The subshell.
    pub fn subshell(&self) -> Subshell {
        self.subshell
    }

    /// Binding energy, eV.
    pub fn binding_energy_ev(&self) -> f64 {
        self.binding_energy_ev
    }

    /// Number of electrons in the subshell in the neutral ground state.
    pub fn occupancy(&self) -> f64 {
        self.occupancy
    }
}

/// Absolute tolerance on `Σ occupancy = Z`. Source tables quote fractional
/// occupancies to two or three digits (e.g. 0.67 + 1.33).
pub const OCCUPANCY_SUM_TOLERANCE: f64 = 1.0e-3;

/// The occupied subshells of one neutral atom.
///
/// **Absent shells are not listed.** A subshell that is unoccupied in the
/// neutral ground state has no entry; [`AtomBindings::binding_energy_ev`]
/// returns `None` for it. There is no zero or sentinel binding energy.
///
/// Invariants: `1 <= z <= 92` (the crate's element table); at least one
/// subshell; subshells strictly increasing by designator (so no duplicates);
/// occupancies summing to `z` within [`OCCUPANCY_SUM_TOLERANCE`] (a neutral
/// atom).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawAtomBindings", into = "RawAtomBindings")]
pub struct AtomBindings {
    z: u8,
    shells: Vec<ShellBinding>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawAtomBindings {
    z: u8,
    shells: Vec<ShellBinding>,
}

impl TryFrom<RawAtomBindings> for AtomBindings {
    type Error = ElectronDataError;
    fn try_from(r: RawAtomBindings) -> Result<Self> {
        Self::new(r.z, r.shells)
    }
}

impl From<AtomBindings> for RawAtomBindings {
    fn from(a: AtomBindings) -> Self {
        Self {
            z: a.z,
            shells: a.shells,
        }
    }
}

impl AtomBindings {
    /// Build a validated atom. Shells must be given in increasing designator
    /// order (the order of the ENDF format).
    pub fn new(z: u8, shells: Vec<ShellBinding>) -> Result<Self> {
        if crate::elements::element(z).is_none() {
            return Err(invalid(
                "atomic number",
                format!("Z = {z} is not in 1..={}", crate::elements::NUM_ELEMENTS),
            ));
        }
        if shells.is_empty() {
            return Err(invalid("subshells", format!("Z = {z} has no subshells")));
        }
        for w in shells.windows(2) {
            let (a, b) = (w[0].subshell, w[1].subshell);
            if a == b {
                return Err(invalid(
                    "subshells",
                    format!("Z = {z} lists subshell {} twice", a.label()),
                ));
            }
            if b < a {
                return Err(invalid(
                    "subshells",
                    format!(
                        "Z = {z}: subshells must be in increasing designator order, got {} after {}",
                        b.label(),
                        a.label()
                    ),
                ));
            }
        }
        let electrons: f64 = shells.iter().map(|s| s.occupancy).sum();
        if (electrons - f64::from(z)).abs() > OCCUPANCY_SUM_TOLERANCE {
            return Err(invalid(
                "occupancy",
                format!("Z = {z}: occupancies sum to {electrons}, not {z} (neutral atom)"),
            ));
        }
        Ok(Self { z, shells })
    }

    /// Atomic number.
    pub fn z(&self) -> u8 {
        self.z
    }

    /// Occupied subshells, in increasing designator order.
    pub fn shells(&self) -> &[ShellBinding] {
        &self.shells
    }

    /// The entry for a subshell, or `None` if it is unoccupied.
    pub fn shell(&self, subshell: Subshell) -> Option<&ShellBinding> {
        self.shells
            .binary_search_by_key(&subshell, |s| s.subshell)
            .ok()
            .map(|i| &self.shells[i])
    }

    /// Binding energy of a subshell, eV, or `None` if it is unoccupied.
    pub fn binding_energy_ev(&self, subshell: Subshell) -> Option<f64> {
        self.shell(subshell).map(|s| s.binding_energy_ev)
    }
}

/// Subshell binding energies for a set of elements, with provenance.
///
/// Invariants: non-blank provenance; at least one atom; atoms strictly
/// increasing in `Z` (no duplicates). Coverage of a required range of `Z` is
/// checked separately by [`SubshellBindingTable::require_coverage`], because
/// a table for a few elements is legitimate.
///
/// No table is committed to this repository: the redistribution terms of the
/// candidate source (EADL, as distributed in EPICS2017) are an open question
/// (`docs/data-provenance.md`). A user supplies a local copy and reads it with
/// [`SubshellBindingTable::from_endf6_mf28_file`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawSubshellBindingTable", into = "RawSubshellBindingTable")]
pub struct SubshellBindingTable {
    provenance: String,
    atoms: Vec<AtomBindings>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSubshellBindingTable {
    provenance: Option<String>,
    atoms: Vec<AtomBindings>,
}

impl TryFrom<RawSubshellBindingTable> for SubshellBindingTable {
    type Error = ElectronDataError;
    fn try_from(r: RawSubshellBindingTable) -> Result<Self> {
        let provenance = checked_text("provenance", r.provenance.as_deref())?;
        Self::new(provenance, r.atoms)
    }
}

impl From<SubshellBindingTable> for RawSubshellBindingTable {
    fn from(t: SubshellBindingTable) -> Self {
        Self {
            provenance: Some(t.provenance),
            atoms: t.atoms,
        }
    }
}

impl SubshellBindingTable {
    /// Build a validated table; atoms in increasing `Z`.
    pub fn new(provenance: impl Into<String>, atoms: Vec<AtomBindings>) -> Result<Self> {
        let provenance = checked_text("provenance", Some(&provenance.into()))?;
        if atoms.is_empty() {
            return Err(invalid("binding table", "has no atoms"));
        }
        if let Some(w) = atoms.windows(2).find(|w| w[1].z <= w[0].z) {
            return Err(invalid(
                "binding table",
                format!(
                    "atoms must be in strictly increasing Z (no duplicates), got Z = {} after Z = {}",
                    w[1].z, w[0].z
                ),
            ));
        }
        Ok(Self { provenance, atoms })
    }

    /// Provenance string (to be copied into run metadata).
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// All atoms, in increasing `Z`.
    pub fn atoms(&self) -> &[AtomBindings] {
        &self.atoms
    }

    /// The atom with this `Z`, if tabulated.
    pub fn atom(&self, z: u8) -> Option<&AtomBindings> {
        self.atoms
            .binary_search_by_key(&z, |a| a.z)
            .ok()
            .map(|i| &self.atoms[i])
    }

    /// The `Z` in `range` that the table lacks.
    pub fn missing(&self, range: std::ops::RangeInclusive<u8>) -> Vec<u8> {
        range.filter(|&z| self.atom(z).is_none()).collect()
    }

    /// An error naming the missing elements unless every `Z` in `range` is
    /// tabulated.
    pub fn require_coverage(&self, range: std::ops::RangeInclusive<u8>) -> Result<()> {
        let (lo, hi) = (*range.start(), *range.end());
        let missing = self.missing(range);
        if missing.is_empty() {
            Ok(())
        } else {
            Err(invalid(
                "binding table coverage",
                format!("Z = {lo}..={hi} required; missing {missing:?}"),
            ))
        }
    }

    /// Parse and validate the TOML form:
    ///
    /// ```toml
    /// provenance = "..."
    /// [[atoms]]
    /// z = 1
    /// shells = [{ subshell = 1, binding_energy_ev = 13.6, occupancy = 1.0 }]
    /// ```
    pub fn from_toml_str(text: &str) -> Result<Self> {
        let raw: RawSubshellBindingTable = toml::from_str(text).map_err(parse_error)?;
        Self::try_from(raw)
    }

    /// Serialize to the TOML form.
    pub fn to_toml_string(&self) -> Result<String> {
        toml::to_string(self).map_err(|e| ElectronDataError::Parse(e.to_string()))
    }

    /// The committed EADL2017 table, Z = 1..92 (EADL as distributed in
    /// EPICS2017; D. E. Cullen, IAEA-NDS-224 Rev. 1 (April 2018), issued by the
    /// IAEA Nuclear Data Section and the NNDC). The data and the credit
    /// notice are in `eadl2017_binding.toml`, embedded in the library; the
    /// provenance row is in `docs/data-provenance.md`.
    ///
    /// # Panics
    ///
    /// Never in a correct build: the embedded file is validated by the test
    /// suite.
    pub fn eadl2017() -> Self {
        Self::from_toml_str(EADL2017_TOML).expect("the embedded EADL2017 table is valid")
    }

    /// Read subshell binding energies and occupancies from ENDF-6 File 28
    /// (MF = 28, MT = 533, atomic relaxation data), the format in which EADL
    /// is distributed (ENDF-6 Formats Manual, CSEWG Document ENDF-102,
    /// BNL-203218-2018-INRE (2018), section 28.2).
    ///
    /// For each material the section is a HEAD record `[ZA, AWR, 0, 0, NSS,
    /// 0]` followed by NSS LIST records `[SUBI, 0, 0, 0, NW, NTR / EBI, ELN,
    /// 0, 0, 0, 0, (SUBJ, SUBK, ETR, FTR, 0, 0) x NTR]` with `NW = 6 (1 +
    /// NTR)`. Only `SUBI` (designator), `EBI` (binding energy, eV) and `ELN`
    /// (electrons in the neutral atom) are kept; the transition data are
    /// skipped. Lines of other files and sections are ignored. Materials with
    /// `Z` beyond the crate's element table (92) are skipped; every other
    /// record must be well formed, or the whole read fails.
    ///
    /// `provenance` is required (non-blank): name the file, its version and
    /// where it was retrieved, so the run metadata records it.
    pub fn from_endf6_mf28_str(text: &str, provenance: &str) -> Result<Self> {
        let provenance = checked_text("provenance", Some(provenance))?;
        let atoms = endf6::read_mf28(text)?;
        Self::new(provenance, atoms)
    }

    /// [`SubshellBindingTable::from_endf6_mf28_str`] on a file.
    pub fn from_endf6_mf28_file(path: impl AsRef<Path>, provenance: &str) -> Result<Self> {
        Self::from_endf6_mf28_str(&read_text(path.as_ref())?, provenance)
    }
}

/// A minimal reader for the ENDF-6 records used by File 28.
///
/// Written from the ENDF-6 Formats Manual (ENDF-102, 2018): an 80-column
/// card has six 11-column data fields (columns 1-66), then MAT (67-70), MF
/// (71-72) and MT (73-75). Floating-point fields may omit the `E` of the
/// exponent (`1.234567+2`); the EPICS files also use Fortran `D` exponents.
/// A blank field is zero.
mod endf6 {
    use super::{AtomBindings, ElectronDataError, Result, ShellBinding, Subshell};

    struct Card<'a> {
        line_no: usize,
        fields: &'a str,
        mat: i64,
    }

    fn parse_err(line_no: usize, reason: impl std::fmt::Display) -> ElectronDataError {
        ElectronDataError::Parse(format!("ENDF line {line_no}: {reason}"))
    }

    /// An ENDF floating-point field: blank is zero, `D` exponents are
    /// accepted, and the `E` of the exponent may be omitted.
    pub(super) fn parse_float(raw: &str) -> Option<f64> {
        let raw = raw.trim();
        if raw.is_empty() {
            return Some(0.0);
        }
        let s = raw.replace(['D', 'd'], "E");
        if let Ok(v) = s.parse::<f64>() {
            return Some(v);
        }
        // E-less exponent: the last sign that is not the leading one.
        let p = s[1..].rfind(['+', '-'])? + 1;
        let (m, e) = s.split_at(p);
        format!("{m}E{e}").parse::<f64>().ok()
    }

    fn float(card: &Card, k: usize) -> Result<f64> {
        let raw = &card.fields[11 * k..11 * (k + 1)];
        parse_float(raw).ok_or_else(|| {
            parse_err(
                card.line_no,
                format!("field {} is not a number: {:?}", k + 1, raw.trim()),
            )
        })
    }

    fn int(card: &Card, k: usize) -> Result<i64> {
        let raw = card.fields[11 * k..11 * (k + 1)].trim();
        if raw.is_empty() {
            return Ok(0);
        }
        raw.parse::<i64>().map_err(|_| {
            parse_err(
                card.line_no,
                format!("field {} is not an integer: {raw:?}", k + 1),
            )
        })
    }

    /// The MF = 28, MT = 533 cards of a file, in order.
    fn mf28_cards(text: &str) -> Result<Vec<Card<'_>>> {
        let mut cards = Vec::new();
        for (i, line) in text.lines().enumerate() {
            let line_no = i + 1;
            let line = line.trim_end_matches('\r');
            if line.len() < 75 {
                continue; // too short to carry MF/MT: not a File 28 card
            }
            if !line.is_ascii() {
                return Err(parse_err(line_no, "non-ASCII characters"));
            }
            let (mf, mt) = (line[70..72].trim(), line[72..75].trim());
            if mf != "28" || mt != "533" {
                continue;
            }
            let mat = line[66..70]
                .trim()
                .parse::<i64>()
                .map_err(|_| parse_err(line_no, "MAT field is not an integer"))?;
            cards.push(Card {
                line_no,
                fields: &line[..66],
                mat,
            });
        }
        Ok(cards)
    }

    pub(super) fn read_mf28(text: &str) -> Result<Vec<AtomBindings>> {
        let cards = mf28_cards(text)?;
        if cards.is_empty() {
            return Err(ElectronDataError::Parse(
                "no ENDF-6 MF=28 MT=533 records found".into(),
            ));
        }
        let mut atoms = Vec::new();
        let mut start = 0;
        while start < cards.len() {
            let mat = cards[start].mat;
            let end = start
                + cards[start..]
                    .iter()
                    .position(|c| c.mat != mat)
                    .unwrap_or(cards.len() - start);
            if let Some(atom) = read_section(&cards[start..end])? {
                atoms.push(atom);
            }
            start = end;
        }
        Ok(atoms)
    }

    /// One material's section: HEAD then NSS subshell LIST records.
    fn read_section(cards: &[Card]) -> Result<Option<AtomBindings>> {
        let head = &cards[0];
        let za = float(head, 0)?;
        if !(za.is_finite() && za >= 1000.0 && za.fract() == 0.0 && za % 1000.0 == 0.0) {
            return Err(parse_err(
                head.line_no,
                format!("ZA = {za} is not an elemental ZA (1000 Z)"),
            ));
        }
        let z = za / 1000.0;
        let nss = int(head, 4)?;
        if nss < 1 {
            return Err(parse_err(head.line_no, format!("NSS = {nss}, need >= 1")));
        }
        let mut next = 1;
        let mut shells = Vec::new();
        for _ in 0..nss {
            let Some(h) = cards.get(next) else {
                return Err(parse_err(
                    cards[cards.len() - 1].line_no,
                    "section ends before all NSS subshells were read",
                ));
            };
            let subi = float(h, 0)?;
            let nw = int(h, 4)?;
            let ntr = int(h, 5)?;
            if ntr < 0 || nw != 6 * (1 + ntr) {
                return Err(parse_err(
                    h.line_no,
                    format!("NW = {nw} with NTR = {ntr}; expected NW = 6 (1 + NTR)"),
                ));
            }
            if !(subi.fract() == 0.0 && (1.0..=f64::from(Subshell::MAX_DESIGNATOR)).contains(&subi))
            {
                return Err(parse_err(
                    h.line_no,
                    format!("SUBI = {subi} is not a subshell designator"),
                ));
            }
            let data_lines = usize::try_from(nw / 6).unwrap_or(0);
            let Some(first) = cards.get(next + 1) else {
                return Err(parse_err(h.line_no, "LIST record has no data line"));
            };
            if cards.len() < next + 1 + data_lines {
                return Err(parse_err(h.line_no, "LIST record is truncated"));
            }
            let (ebi, eln) = (float(first, 0)?, float(first, 1)?);
            // Validate the skipped transition lines as numbers too, so a
            // corrupt record is not silently accepted.
            for c in &cards[next + 1..next + 1 + data_lines] {
                for k in 0..6 {
                    float(c, k)?;
                }
            }
            next += 1 + data_lines;
            shells.push((subi as u8, ebi, eln, h.line_no));
        }
        if next != cards.len() {
            return Err(parse_err(
                cards[next].line_no,
                "unexpected records after the last subshell",
            ));
        }
        if z > crate::elements::NUM_ELEMENTS as f64 {
            return Ok(None);
        }
        let shells = shells
            .into_iter()
            .map(|(d, ebi, eln, line_no)| {
                ShellBinding::new(Subshell::new(d)?, ebi, eln).map_err(|e| parse_err(line_no, e))
            })
            .collect::<Result<Vec<_>>>()?;
        AtomBindings::new(z as u8, shells)
            .map(Some)
            .map_err(|e| parse_err(head.line_no, e))
    }
}

// ---------------------------------------------------------------------------
// Cross-section cache
// ---------------------------------------------------------------------------

/// What a [`CrossSectionTable`] samples. Elastic and inelastic tables are
/// distinct types of data and are never interchangeable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SamplingAxis {
    /// Elastic scattering: the polar deflection angle `θ` in radians,
    /// `0 <= θ <= π`.
    ElasticPolarAngle,
    /// Inelastic scattering: the energy loss `W` in eV, `0 <= W <= E` for an
    /// incident energy `E`.
    InelasticEnergyLoss,
}

impl SamplingAxis {
    /// Unit of the sampled variable, as written in the cache file.
    pub fn unit(self) -> &'static str {
        match self {
            SamplingAxis::ElasticPolarAngle => "rad",
            SamplingAxis::InelasticEnergyLoss => "eV",
        }
    }
}

/// The unvalidated contents of a [`CrossSectionTable`], for
/// [`CrossSectionTable::new`]. See [`CrossSectionTable`] for the meaning and
/// invariants of each field.
#[derive(Debug, Clone, PartialEq)]
pub struct CrossSectionTableParts {
    /// Identity of the model that produced the table (name and version).
    pub model: String,
    /// Material identity.
    pub material: String,
    /// Provenance: the inputs of the model and their origin.
    pub provenance: String,
    /// The sampled variable.
    pub axis: SamplingAxis,
    /// Incident electron energies, eV.
    pub energy_ev: Vec<f64>,
    /// Inverse mean free path at each energy, m⁻¹.
    pub inverse_mfp_per_m: Vec<f64>,
    /// Cumulative-probability grid of the inverse CDFs.
    pub probability: Vec<f64>,
    /// Inverse-CDF ordinates, one row per energy (empty for zero-rate rows).
    pub quantiles: Vec<Vec<f64>>,
}

/// A precomputed cross-section table for one material and one process,
/// cacheable on disk (format version [`CACHE_FORMAT_VERSION`]).
///
/// # Contents and dimensions
///
/// With `N` incident energies and `M` probability points:
///
/// - `energy_ev[N]`: incident electron energies, eV; finite, positive,
///   strictly increasing, `N >= 2`.
/// - `inverse_mfp_per_m[N]`: inverse mean free path `λ⁻¹(E) = n σ(E)` for the
///   process, in m⁻¹; finite and non-negative. Zero means no interaction at
///   that energy (for example below a threshold).
/// - `probability[M]`: the cumulative-probability axis shared by every row;
///   finite, strictly increasing, starting at exactly 0 and ending at exactly
///   1, `M >= 2`.
/// - `quantiles[N][..]`: **inverse-CDF ordinates** (the cache stores the
///   inverse CDF, not the forward CDF). `quantiles[i][j]` is the value `x` of
///   the sampled variable at which the cumulative distribution at energy
///   `energy_ev[i]` reaches `probability[j]`. A row is non-decreasing and lies
///   in the axis's domain (`[0, π]` rad for the elastic polar angle, `[0,
///   energy_ev[i]]` eV for the inelastic energy loss). A row whose inverse
///   mean free path is zero must be **empty**: there is no distribution to
///   store, and [`CrossSectionTable::inverse_cdf`] returns
///   [`ElectronDataError::ZeroRate`] for it. Every other row has exactly `M`
///   entries.
///
/// The file also records the model identity, the material, the provenance,
/// the sampling axis and its unit (`ordinate_unit`, which must match the
/// axis), and the format version.
///
/// # Sampling
///
/// [`CrossSectionTable::inverse_cdf`] evaluates a stored row at a
/// probability `u` by linear interpolation between neighbouring probability
/// points (exact at the stored points). Interpolation between incident
/// energies is left to the transport code, which owns that choice.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawCacheIn", into = "RawCacheOut")]
pub struct CrossSectionTable {
    parts: CrossSectionTableParts,
}

/// Deserialization form. Every field is optional so that the format version
/// is checked first, whatever else the file holds; unknown keys are collected
/// and rejected after the version check.
#[derive(Debug, Deserialize)]
struct RawCacheIn {
    format_version: Option<u32>,
    model: Option<String>,
    material: Option<String>,
    provenance: Option<String>,
    axis: Option<SamplingAxis>,
    ordinate_unit: Option<String>,
    energy_ev: Option<Vec<f64>>,
    inverse_mfp_per_m: Option<Vec<f64>>,
    probability: Option<Vec<f64>>,
    quantiles: Option<Vec<Vec<f64>>>,
    #[serde(flatten)]
    unknown: BTreeMap<String, IgnoredAny>,
}

#[derive(Debug, Clone, Serialize)]
struct RawCacheOut {
    format_version: u32,
    model: String,
    material: String,
    provenance: String,
    axis: SamplingAxis,
    ordinate_unit: String,
    energy_ev: Vec<f64>,
    inverse_mfp_per_m: Vec<f64>,
    probability: Vec<f64>,
    quantiles: Vec<Vec<f64>>,
}

fn required<T>(field: &'static str, v: Option<T>) -> Result<T> {
    v.ok_or_else(|| ElectronDataError::Parse(format!("missing field `{field}`")))
}

impl TryFrom<RawCacheIn> for CrossSectionTable {
    type Error = ElectronDataError;
    fn try_from(r: RawCacheIn) -> Result<Self> {
        if r.format_version != Some(CACHE_FORMAT_VERSION) {
            return Err(ElectronDataError::UnsupportedVersion {
                found: r.format_version,
                supported: CACHE_FORMAT_VERSION,
            });
        }
        if let Some(k) = r.unknown.keys().next() {
            return Err(ElectronDataError::Parse(format!(
                "unknown field `{k}` in a version {CACHE_FORMAT_VERSION} cache"
            )));
        }
        let provenance = checked_text("provenance", r.provenance.as_deref())?;
        let model = checked_text("model", r.model.as_deref())?;
        let material = checked_text("material", r.material.as_deref())?;
        let axis = required("axis", r.axis)?;
        let unit = checked_text("ordinate_unit", r.ordinate_unit.as_deref())?;
        if unit != axis.unit() {
            return Err(ElectronDataError::MalformedText {
                field: "ordinate_unit",
                reason: format!(
                    "{unit:?} does not match axis {axis:?} (unit {:?})",
                    axis.unit()
                ),
            });
        }
        Self::new(CrossSectionTableParts {
            model,
            material,
            provenance,
            axis,
            energy_ev: required("energy_ev", r.energy_ev)?,
            inverse_mfp_per_m: required("inverse_mfp_per_m", r.inverse_mfp_per_m)?,
            probability: required("probability", r.probability)?,
            quantiles: required("quantiles", r.quantiles)?,
        })
    }
}

impl From<CrossSectionTable> for RawCacheOut {
    fn from(t: CrossSectionTable) -> Self {
        let p = t.parts;
        Self {
            format_version: CACHE_FORMAT_VERSION,
            ordinate_unit: p.axis.unit().to_string(),
            model: p.model,
            material: p.material,
            provenance: p.provenance,
            axis: p.axis,
            energy_ev: p.energy_ev,
            inverse_mfp_per_m: p.inverse_mfp_per_m,
            probability: p.probability,
            quantiles: p.quantiles,
        }
    }
}

impl CrossSectionTable {
    /// Build a validated table. See the type docs for the invariants.
    pub fn new(parts: CrossSectionTableParts) -> Result<Self> {
        let provenance = checked_text("provenance", Some(&parts.provenance))?;
        let model = checked_text("model", Some(&parts.model))?;
        let material = checked_text("material", Some(&parts.material))?;
        let n = parts.energy_ev.len();
        check_energy_grid("incident energy grid", &parts.energy_ev)?;
        if parts.inverse_mfp_per_m.len() != n {
            return Err(invalid(
                "inverse mean free path",
                format!(
                    "need one value per energy: {n} energies, {} values",
                    parts.inverse_mfp_per_m.len()
                ),
            ));
        }
        check_nonnegative("inverse mean free path", &parts.inverse_mfp_per_m)?;
        let p = &parts.probability;
        if p.len() < 2 {
            return Err(invalid("probability axis", "needs at least 2 points"));
        }
        if p.iter().any(|x| !x.is_finite()) {
            return Err(invalid("probability axis", "values must be finite"));
        }
        if p[0] != 0.0 || p[p.len() - 1] != 1.0 {
            return Err(invalid(
                "probability axis",
                format!(
                    "must start at 0 and end at 1, got {} .. {}",
                    p[0],
                    p[p.len() - 1]
                ),
            ));
        }
        if let Some(w) = p.windows(2).find(|w| w[1] <= w[0]) {
            return Err(invalid(
                "probability axis",
                format!("must be strictly increasing, got {} then {}", w[0], w[1]),
            ));
        }
        if parts.quantiles.len() != n {
            return Err(invalid(
                "inverse-CDF rows",
                format!(
                    "need one row per energy: {n} energies, {} rows",
                    parts.quantiles.len()
                ),
            ));
        }
        for (i, row) in parts.quantiles.iter().enumerate() {
            let e = parts.energy_ev[i];
            if parts.inverse_mfp_per_m[i] == 0.0 {
                if !row.is_empty() {
                    return Err(invalid(
                        "inverse-CDF rows",
                        format!(
                            "row {i} ({e} eV) has zero inverse mean free path and must be empty"
                        ),
                    ));
                }
                continue;
            }
            if row.len() != p.len() {
                return Err(invalid(
                    "inverse-CDF rows",
                    format!(
                        "row {i} ({e} eV) has {} entries; the probability axis has {}",
                        row.len(),
                        p.len()
                    ),
                ));
            }
            let max = match parts.axis {
                SamplingAxis::ElasticPolarAngle => std::f64::consts::PI,
                SamplingAxis::InelasticEnergyLoss => e,
            };
            if let Some(x) = row
                .iter()
                .find(|x| !(x.is_finite() && **x >= 0.0 && **x <= max))
            {
                return Err(invalid(
                    "inverse-CDF rows",
                    format!(
                        "row {i} ({e} eV): {x} is outside [0, {max}] {}",
                        parts.axis.unit()
                    ),
                ));
            }
            if let Some(w) = row.windows(2).find(|w| w[1] < w[0]) {
                return Err(invalid(
                    "inverse-CDF rows",
                    format!(
                        "row {i} ({e} eV) must be non-decreasing, got {} then {}",
                        w[0], w[1]
                    ),
                ));
            }
        }
        Ok(Self {
            parts: CrossSectionTableParts {
                model,
                material,
                provenance,
                ..parts
            },
        })
    }

    /// Parse and validate the TOML cache form.
    pub fn from_toml_str(text: &str) -> Result<Self> {
        // Parse the raw form, then convert, so the caller gets the typed
        // error (an unsupported version, a missing provenance) rather than
        // the message toml wraps a conversion error in.
        let raw: RawCacheIn = toml::from_str(text).map_err(parse_error)?;
        Self::try_from(raw)
    }

    /// Serialize to the TOML cache form.
    pub fn to_toml_string(&self) -> Result<String> {
        toml::to_string(self).map_err(|e| ElectronDataError::Parse(e.to_string()))
    }

    /// Load a TOML cache file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_toml_str(&read_text(path.as_ref())?)
    }

    /// Write a TOML cache file.
    pub fn write_toml_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        std::fs::write(path, self.to_toml_string()?)
            .map_err(|e| ElectronDataError::Io(format!("{}: {e}", path.display())))
    }

    /// Format version of this table ([`CACHE_FORMAT_VERSION`]).
    pub fn format_version(&self) -> u32 {
        CACHE_FORMAT_VERSION
    }

    /// Model identity.
    pub fn model(&self) -> &str {
        &self.parts.model
    }

    /// Material identity.
    pub fn material(&self) -> &str {
        &self.parts.material
    }

    /// Provenance string.
    pub fn provenance(&self) -> &str {
        &self.parts.provenance
    }

    /// The sampled variable.
    pub fn axis(&self) -> SamplingAxis {
        self.parts.axis
    }

    /// Incident energies, eV.
    pub fn energy_ev(&self) -> &[f64] {
        &self.parts.energy_ev
    }

    /// Inverse mean free paths, m⁻¹, one per energy.
    pub fn inverse_mfp_per_m(&self) -> &[f64] {
        &self.parts.inverse_mfp_per_m
    }

    /// The cumulative-probability axis.
    pub fn probability(&self) -> &[f64] {
        &self.parts.probability
    }

    /// The inverse-CDF ordinates at grid energy `i`, or `None` for a
    /// zero-rate row (or an index past the end).
    pub fn quantiles(&self, i: usize) -> Option<&[f64]> {
        self.parts
            .quantiles
            .get(i)
            .filter(|r| !r.is_empty())
            .map(Vec::as_slice)
    }

    /// The inverse CDF of grid row `i` at cumulative probability `u` in
    /// `[0, 1]`, in the axis's unit. A zero-rate row is
    /// [`ElectronDataError::ZeroRate`]: the caller must not sample an
    /// interaction that cannot happen.
    pub fn inverse_cdf(&self, i: usize, u: f64) -> Result<f64> {
        let n = self.parts.energy_ev.len();
        if i >= n {
            return Err(ElectronDataError::IndexOutOfRange { index: i, len: n });
        }
        check_in_range("cumulative probability", u, 0.0, 1.0)?;
        let Some(row) = self.quantiles(i) else {
            return Err(ElectronDataError::ZeroRate {
                energy_ev: self.parts.energy_ev[i],
            });
        };
        Ok(lerp_on_grid(&self.parts.probability, row, u))
    }

    /// The validated contents.
    pub fn parts(&self) -> &CrossSectionTableParts {
        &self.parts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subshell_table_matches_endf_mt_numbering() {
        // Designator n is MT 533 + n: K = 534, L1 = 535, ..., Q3 = 572.
        assert_eq!(Subshell::K.label(), "K");
        assert_eq!(Subshell::new(4).unwrap().label(), "L3");
        assert_eq!(Subshell::new(4).unwrap().orbital(), "2p3/2");
        assert_eq!(Subshell::new(39).unwrap().label(), "Q3");
        assert_eq!(Subshell::new(16).unwrap().capacity(), 8); // 4f7/2
        assert_eq!(Subshell::from_label("P11").unwrap().designator(), 36);
        assert!(Subshell::new(0).is_err());
        assert!(Subshell::new(40).is_err());
        // Every capacity is 2j + 1 with j read from the orbital string.
        for d in 1..=Subshell::MAX_DESIGNATOR {
            let s = Subshell::new(d).unwrap();
            let j2: u8 = s.orbital().split('/').next().unwrap()[2..].parse().unwrap();
            assert_eq!(s.capacity(), j2 + 1, "{}", s.label());
        }
    }

    #[test]
    fn endf_floats_with_and_without_exponent_letter() {
        let cases = [
            (" 1.500000+3", Some(1500.0)),
            ("1.0000D+11", Some(1.0e11)),
            ("13.6000000", Some(13.6)),
            (" .999241400", Some(0.9992414)),
            ("-2.5-1", Some(-0.25)),
            ("", Some(0.0)),
            ("7", Some(7.0)),
            ("1.2.3", None),
        ];
        for (raw, want) in cases {
            assert_eq!(endf6::parse_float(raw), want, "{raw:?}");
        }
    }
}
