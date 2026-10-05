//! Electronic energy-loss straggling.
//!
//! Bohr (1948), Mat. Fys. Medd. Dan. Vid. Selsk. 18 (8): for a thick-target
//! free-electron estimate the variance of the energy loss per unit areal
//! density is
//!
//! `Ω_B² = 4π z² Z2 (e²/(4π ε0))²` (J² m² per target atom),
//!
//! independent of energy in the non-relativistic regime. A relativistic factor
//! `(1 − β²/2) / (1 − β²)` (Bethe and Livingston 1937; Fano 1963) is included
//! in [`StragglingModel::BohrRelativistic`].
//!
//! # Not implemented
//!
//! The Chu and Yang-O'Connor-Wang corrections, which reduce Bohr's value at
//! intermediate energies, are published as Hartree-Fock-Slater tables and as
//! fits to data (declined, see `docs/stopping-models.md`), and the Lindhard-
//! Scharff low-velocity correction was not verified against the original. They
//! are omitted. Bohr is therefore an upper-end estimate below about
//! 1 MeV/u; see `docs/stopping-models.md`.

use super::bethe::EffectiveCharge;
use super::{check_energy, target, Ion, StoppingError};
use crate::constants::COULOMB_E2;
use crate::material::Material;
use std::f64::consts::PI;

/// Straggling model choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StragglingModel {
    /// Bohr, non-relativistic.
    #[default]
    Bohr,
    /// Bohr with the relativistic factor.
    BohrRelativistic,
}

/// Energy-loss variance per target atom of element `target_z`, J² m².
pub fn variance_per_atom(
    model: StragglingModel,
    charge: EffectiveCharge,
    ion: &Ion,
    target_z: u8,
    energy_ev: f64,
) -> Result<f64, StoppingError> {
    check_energy(energy_ev)?;
    target(target_z)?;
    let mc2 = ion.mass_kg() * crate::constants::SPEED_OF_LIGHT.powi(2);
    let gamma = 1.0 + energy_ev * crate::units::J_PER_EV / mc2;
    let beta2 = 1.0 - 1.0 / (gamma * gamma);
    let z = charge.charge(ion.z(), beta2.sqrt());
    let base = 4.0 * PI * z * z * f64::from(target_z) * COULOMB_E2 * COULOMB_E2;
    Ok(match model {
        StragglingModel::Bohr => base,
        StragglingModel::BohrRelativistic => base * (1.0 - 0.5 * beta2) / (1.0 - beta2),
    })
}

/// Bragg-additive variance per average atom of a material, J² m²
/// (`Σ_j x_j Ω_j²`).
pub fn variance_per_atom_material(
    model: StragglingModel,
    charge: EffectiveCharge,
    ion: &Ion,
    material: &Material,
    energy_ev: f64,
) -> Result<f64, StoppingError> {
    let mut s = 0.0;
    for c in material.components() {
        s += c.atom_fraction() * variance_per_atom(model, charge, ion, c.z(), energy_ev)?;
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bohr_scales_as_z1_squared_z2() {
        let f = |z1, z2| {
            variance_per_atom(
                StragglingModel::Bohr,
                EffectiveCharge::Bare,
                &Ion::new(z1).unwrap(),
                z2,
                1e5,
            )
            .unwrap()
        };
        assert!((f(2, 14) / f(1, 14) - 4.0).abs() < 1e-12);
        assert!((f(1, 28) / f(1, 14) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn bohr_magnitude_for_proton_in_si() {
        // 4 pi * 14 * (2.307e-28 J m)^2 = 9.36e-54 J^2 m^2
        let v = variance_per_atom(
            StragglingModel::Bohr,
            EffectiveCharge::Bare,
            &Ion::proton(),
            14,
            1e6,
        )
        .unwrap();
        assert!((v / 9.36e-54 - 1.0).abs() < 5e-3, "{v}");
    }

    #[test]
    fn relativistic_factor_exceeds_one() {
        let g = |m| variance_per_atom(m, EffectiveCharge::Bare, &Ion::proton(), 14, 5e8).unwrap();
        assert!(g(StragglingModel::BohrRelativistic) > g(StragglingModel::Bohr));
    }
}
