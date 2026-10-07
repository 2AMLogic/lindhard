//! Electronic stopping of ions in matter, from published formulas only.
//!
//! # Conventions
//!
//! * The stopping quantity returned by every model is the **stopping cross
//!   section per target atom**, `S_e = (1/N) (-dE/dx)`, in **J m² per atom**
//!   (SI). [`to_ev_1e15_cm2`] converts to the customary eV·10⁻¹⁵ cm²/atom.
//! * Energies are the ion's kinetic energy in the laboratory frame, in eV
//!   (`_ev` suffix at the API boundary, as in [`crate::units`]).
//! * Reduced quantities (`ε`, `ρ`) are the Lindhard-Scharff ones, defined in
//!   [`lindhard_scharff`].
//! * Every model is a pure function of `(ion, target Z, energy)`: no RNG, no
//!   state, so results are bitwise reproducible across threads.
//!
//! Validity ranges per model are in `docs/stopping-models.md` and are exposed
//! at run time by [`ElectronicStopping::validity`].
//!
//! No tabulated SRIM or ICRU data enters this module. Terms whose published
//! coefficients are only available from such tables are omitted, and the
//! omission is listed in `docs/stopping-models.md`.

pub mod bethe;
pub mod bragg;
pub mod dataset;
pub mod lindhard_scharff;
pub mod mix;
pub mod oen_robinson;
pub mod straggling;
pub mod table;

use crate::constants::{ATOMIC_MASS_UNIT, FINE_STRUCTURE, SPEED_OF_LIGHT};
use crate::elements::element;
use crate::units::J_PER_EV;

/// Errors from stopping models, tables and Bragg sums.
#[derive(Debug, Clone, thiserror::Error, PartialEq)]
pub enum StoppingError {
    /// Energy was not finite and positive.
    #[error("energy must be finite and positive, got {0} eV")]
    InvalidEnergy(f64),
    /// Atomic number outside the element table (1..=92).
    #[error("unknown element Z = {0}")]
    UnknownElement(u8),
    /// The model's formula is not meaningful at this energy (for example the
    /// Bethe logarithm bracket is not positive far below its validity range).
    #[error("{model} is not applicable at {energy_ev} eV for this ion")]
    NotApplicable {
        /// Model name.
        model: &'static str,
        /// Energy, eV.
        energy_ev: f64,
    },
    /// A model parameter (for example a mean excitation energy) was not finite,
    /// or not positive where it must be.
    #[error("invalid model parameter {name} = {value}")]
    InvalidParameter {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },
    /// A user table has no (or an empty) provenance field.
    #[error(
        "stopping table has no provenance; a citation or description of its origin is required"
    )]
    MissingProvenance,
    /// A user table is malformed.
    #[error("invalid stopping table: {0}")]
    InvalidTable(String),
    /// A line of a stopping dataset ([`dataset`]) is malformed or fails
    /// validation.
    #[error("stopping dataset, line {line}: {reason}")]
    InvalidDataset {
        /// 1-based line number in the source text.
        line: usize,
        /// What is wrong with the line.
        reason: String,
    },
    /// A user table was queried outside its energy range.
    #[error("energy {energy_ev} eV outside table range [{min_ev}, {max_ev}] eV")]
    OutOfTableRange {
        /// Queried energy, eV.
        energy_ev: f64,
        /// Table minimum, eV.
        min_ev: f64,
        /// Table maximum, eV.
        max_ev: f64,
    },
    /// A user table was queried for a different ion or target than it holds.
    #[error("table is for Z1={table_ion}, Z2={table_target}; asked for Z1={ion}, Z2={target}")]
    TableMismatch {
        /// Table ion Z.
        table_ion: u8,
        /// Table target Z.
        table_target: u8,
        /// Requested ion Z.
        ion: u8,
        /// Requested target Z.
        target: u8,
    },
    /// Projectile mass was not finite and positive.
    #[error("projectile mass must be finite and positive, got {0} u")]
    InvalidMass(f64),
    /// A user table was queried for a projectile mass different from the one
    /// it was declared for (a different isotope).
    #[error("table is for an ion of mass {table_mass_amu} u; asked for {mass_amu} u")]
    TableMassMismatch {
        /// Mass the table was declared for, u.
        table_mass_amu: f64,
        /// Requested projectile mass, u.
        mass_amu: f64,
    },
    /// A material error while summing a compound.
    #[error(transparent)]
    Material(#[from] crate::material::MaterialError),
}

/// A projectile: atomic number and mass.
#[derive(Debug, Clone, Copy, PartialEq)]
///
/// Fields are private so that every `Ion` satisfies its invariants: a known
/// element and a finite, positive mass.
pub struct Ion {
    z: u8,
    mass_amu: f64,
}

impl Ion {
    /// An ion of element `z` with the standard atomic weight as its mass.
    pub fn new(z: u8) -> Result<Self, StoppingError> {
        let e = element(z).ok_or(StoppingError::UnknownElement(z))?;
        Ok(Self {
            z,
            mass_amu: e.atomic_weight,
        })
    }

    /// An ion with an explicit (isotopic) mass in amu.
    ///
    /// The mass must be finite and positive.
    pub fn with_mass(z: u8, mass_amu: f64) -> Result<Self, StoppingError> {
        element(z).ok_or(StoppingError::UnknownElement(z))?;
        check_mass(mass_amu)?;
        Ok(Self { z, mass_amu })
    }

    /// Atomic number `Z1`.
    pub fn z(&self) -> u8 {
        self.z
    }

    /// Mass in unified atomic mass units.
    pub fn mass_amu(&self) -> f64 {
        self.mass_amu
    }

    /// A proton (CODATA proton mass).
    pub fn proton() -> Self {
        Self {
            z: 1,
            mass_amu: crate::constants::PROTON_MASS / ATOMIC_MASS_UNIT,
        }
    }

    /// Mass in kg.
    pub fn mass_kg(&self) -> f64 {
        self.mass_amu * ATOMIC_MASS_UNIT
    }
}

/// Advisory energy range over which a model is meant to be used, eV.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValidityRange {
    /// Lower bound, eV (0 if none).
    pub min_energy_ev: f64,
    /// Upper bound, eV (`INFINITY` if none).
    pub max_energy_ev: f64,
}

impl ValidityRange {
    /// True if `energy_ev` is inside the range.
    pub fn contains(&self, energy_ev: f64) -> bool {
        energy_ev >= self.min_energy_ev && energy_ev <= self.max_energy_ev
    }
}

/// An electronic stopping model.
pub trait ElectronicStopping {
    /// Short stable model name (used in output metadata).
    fn name(&self) -> &'static str;

    /// Stopping cross section per atom of element `target_z`, J m².
    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError>;

    /// Advisory validity range for this ion.
    fn validity(&self, ion: &Ion) -> ValidityRange;
}

/// Convert J m² per atom to eV·10⁻¹⁵ cm² per atom (1e-15 cm² = 1e-19 m²).
pub fn to_ev_1e15_cm2(s_j_m2: f64) -> f64 {
    s_j_m2 / J_PER_EV / 1.0e-19
}

/// Convert eV·10⁻¹⁵ cm² per atom to J m² per atom.
pub fn from_ev_1e15_cm2(s: f64) -> f64 {
    s * 1.0e-19 * J_PER_EV
}

/// Bohr velocity `v0 = α c`, m/s.
pub fn bohr_velocity() -> f64 {
    FINE_STRUCTURE * SPEED_OF_LIGHT
}

/// Ion speed from kinetic energy, non-relativistic `sqrt(2E/M)`, m/s.
pub fn speed_nonrel(ion: &Ion, energy_ev: f64) -> f64 {
    (2.0 * energy_ev * J_PER_EV / ion.mass_kg()).sqrt()
}

pub(crate) fn check_mass(mass_amu: f64) -> Result<(), StoppingError> {
    if mass_amu.is_finite() && mass_amu > 0.0 {
        Ok(())
    } else {
        Err(StoppingError::InvalidMass(mass_amu))
    }
}

pub(crate) fn check_energy(energy_ev: f64) -> Result<(), StoppingError> {
    if energy_ev.is_finite() && energy_ev > 0.0 {
        Ok(())
    } else {
        Err(StoppingError::InvalidEnergy(energy_ev))
    }
}

pub(crate) fn target(z: u8) -> Result<&'static crate::elements::Element, StoppingError> {
    element(z).ok_or(StoppingError::UnknownElement(z))
}

/// Energy per amu (eV/u) at which the ion speed equals `n_v0 * v0 * Z1^(2/3)`:
/// the Lindhard-Scharff velocity scale that separates the low-energy and
/// Bethe regimes.
pub(crate) fn velocity_scale_energy_per_amu_ev(z1: u8, n_v0: f64) -> f64 {
    let v = n_v0 * bohr_velocity() * f64::from(z1).powf(2.0 / 3.0);
    0.5 * ATOMIC_MASS_UNIT * v * v / J_PER_EV
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn with_mass_rejects_invalid_mass() {
        for m in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                matches!(Ion::with_mass(1, m), Err(StoppingError::InvalidMass(_))),
                "mass {m}"
            );
        }
        assert!(Ion::with_mass(1, 2.014).is_ok());
    }

    #[test]
    fn unit_round_trip() {
        let s = 12.5;
        assert!((to_ev_1e15_cm2(from_ev_1e15_cm2(s)) - s).abs() < 1e-12);
    }

    #[test]
    fn models_are_bitwise_deterministic_across_threads() {
        let ion = Ion::new(15).unwrap();
        let m = bethe::BetheBloch::new();
        let r0 = m.stopping(&ion, 14, 1.2e8).unwrap();
        let hs: Vec<_> = (0..4)
            .map(|_| {
                std::thread::spawn(move || {
                    let ion = Ion::new(15).unwrap();
                    bethe::BetheBloch::new().stopping(&ion, 14, 1.2e8).unwrap()
                })
            })
            .collect();
        for h in hs {
            assert_eq!(h.join().unwrap().to_bits(), r0.to_bits());
        }
        assert_eq!(m.stopping(&ion, 14, 1.2e8).unwrap().to_bits(), r0.to_bits());
    }
}
