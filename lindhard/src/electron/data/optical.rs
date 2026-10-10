use super::{
    check_energy_grid, check_in_range, check_nonnegative, checked_text, invalid, lerp_on_grid,
    parse_error, read_text, ElectronDataError, Result,
};
use serde::{Deserialize, Serialize};
use std::path::Path;

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
