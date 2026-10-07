//! Nonrelativistic Bethe stopping, the high-energy reference for the
//! dielectric models.
//!
//! The Bethe equation as given by the Particle Data Group, S. Navas et al.,
//! Phys. Rev. D 110, 030001 (2024), review 34 "Passage of particles through
//! matter" (D. E. Groom and S. R. Klein, rev. August 2023), eq. (34.5):
//!
//! ```text
//! -dE/dx = K z² (Z/A) (1/β²) [ (1/2) ln(2 m_e c² β² γ² W_max / I²) - β² - δ(βγ)/2 ].
//! ```
//!
//! Here it is used **nonrelativistically and for a unit charge** (`z = 1`,
//! `γ -> 1`, the `-β²` and density-effect terms dropped, `m_e c² β² = m_e v²
//! = 2E` for an electron of kinetic energy `E`), and as a linear stopping
//! power for an electron density `n` (the `K (Z/A) ρ` of the mass stopping
//! power is `4π n r_e² m_e c²`):
//!
//! ```text
//! S(E) = (4π n (e²/4πε₀)² / (m_e v²)) · (1/2) ln(2 m_e v² W_max / I²)
//!      = (2π n (e²/4πε₀)² / E) · (1/2) ln(4 E W_max / I²).
//! ```
//!
//! `W_max` is the largest energy loss the model allows. The PDG gives eq.
//! (34.5) for heavy particles, where `W_max` is the kinematic maximum
//! transfer to a free electron; for an incident electron the recommended
//! stopping formula differs (identical-particle exchange). The single-pole
//! Penn model has no exchange and allows losses up to the full kinetic energy
//! (`W_max = E`), so its high-energy limit is the expression above with
//! `W_max = E`, `S = (2π n (e²/4πε₀)² / E) ln(2E/I)`. That is the comparison
//! this module exists for; it is not a recommended electron stopping power.

use crate::constants::{COULOMB_E2, ELEMENTARY_CHARGE};

/// Nonrelativistic Bethe stopping power of a unit-charge particle of
/// electron mass (see the module docs), eV/m.
///
/// `electron_density_per_m3` is the density of the electrons the mean
/// excitation energy `mean_excitation_ev` refers to, `energy_ev` the kinetic
/// energy and `max_loss_ev` the largest energy loss `W_max`. Returns `NaN`
/// unless all four are finite and positive.
pub fn stopping_power_ev_per_m(
    electron_density_per_m3: f64,
    mean_excitation_ev: f64,
    energy_ev: f64,
    max_loss_ev: f64,
) -> f64 {
    let ok = [
        electron_density_per_m3,
        mean_excitation_ev,
        energy_ev,
        max_loss_ev,
    ]
    .iter()
    .all(|v| v.is_finite() && *v > 0.0);
    if !ok {
        return f64::NAN;
    }
    let e_j = energy_ev * ELEMENTARY_CHARGE;
    let log =
        0.5 * (4.0 * energy_ev * max_loss_ev / (mean_excitation_ev * mean_excitation_ev)).ln();
    // J/m -> eV/m
    2.0 * std::f64::consts::PI * electron_density_per_m3 * COULOMB_E2 * COULOMB_E2 / e_j * log
        / ELEMENTARY_CHARGE
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_atomic_unit_form() {
        // In Hartree atomic units S = (Ω²/(2E)) ln(2E/I) with Ω² = 4π n,
        // for W_max = E. Check the SI evaluation against it.
        let a0 = crate::constants::BOHR_RADIUS;
        let eh = crate::constants::HARTREE_ENERGY / ELEMENTARY_CHARGE;
        let (n_au, e_au, i_au) = (0.05, 300.0, 0.9);
        let want_au = 4.0 * std::f64::consts::PI * n_au / (2.0 * e_au) * (2.0 * e_au / i_au).ln();
        let got = stopping_power_ev_per_m(n_au / a0.powi(3), i_au * eh, e_au * eh, e_au * eh);
        let got_au = got * a0 / eh;
        // e²/(4πε₀) = 1 hartree bohr holds to the CODATA rounding, ~1e-10.
        assert!(
            (got_au / want_au - 1.0).abs() < 1e-9,
            "{got_au} vs {want_au}"
        );
        assert!(stopping_power_ev_per_m(1e29, 0.0, 1e3, 1e3).is_nan());
    }
}
