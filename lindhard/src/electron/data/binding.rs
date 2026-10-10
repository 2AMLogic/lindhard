use super::endf6;
use super::{checked_text, invalid, parse_error, read_text, ElectronDataError, Result};
use serde::{Deserialize, Serialize};
use std::path::Path;

/// The committed EADL2017 binding-energy table (see
/// [`SubshellBindingTable::eadl2017`]); its header carries the credit.
const EADL2017_TOML: &str = include_str!("../eadl2017_binding.toml");

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
}
