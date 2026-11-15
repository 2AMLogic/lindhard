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
//! # Use in transport
//!
//! The amorphous BCA engine samples electronic straggling when it is enabled
//! with `[physics] straggling = "bohr"` (or `"bohr-relativistic"`; input key
//! `StragglingChoice`, engine setting `BcaConfig::straggling`). Every free
//! flight's nonlocal loss is then drawn from a Gaussian centred on the mean
//! loss `N S_e s` it would take anyway, with standard deviation
//! `sqrt(n Omega^2 s)`: the Bragg-added variance of this module
//! (`variance_per_atom_material`) scaled by the flight's areal density
//! `n s`. The draw comes from the history's own transport random stream,
//! so results stay deterministic per `(seed, index)` at any thread count.
//! The effective charge is `EffectiveCharge::Bare`, the form in which the
//! Bohr expression is written above; the engine's mean loss applies no
//! charge correction either.
//!
//! Because the Gaussian is unbounded below, a draw can exceed the ion's
//! energy or come out negative (an energy gain). The sampled loss is
//! truncated by clamping to `[0, E]`, the same ceiling the mean loss is
//! clamped to. This biases the mean loss upward when the standard deviation
//! is comparable to the mean loss per flight — thin layers, low energy —
//! which is also the regime where Bohr's variance is an overestimate
//! ("Not implemented" below), so the two errors reinforce rather than
//! cancel. The bias is not cosmetic at short flights: with the constant
//! free path (flights of `N^(-1/3)`, 0.21 nm in carbon) the per-flight
//! standard deviation at MeV energies is several times the per-flight mean
//! loss, and the sampled transmitted variance of a 1 MeV proton through
//! 100 nm of carbon then measures about 0.42 of the Bohr value while the
//! mean loss inflates about 2.5x — 5.8 keV to 14.9 keV of electronic loss
//! per history, energy conservation delivering the rest into the clamped
//! draws (measured 2026-11-14 while building this feature). The input
//! resolution warns when `[physics] straggling` is combined with the
//! constant free path, and when the beam energy sits below Bohr's
//! `v >> v0 Z1^(2/3)` floor. With flights tens of
//! nanometres long (`free_path = "energy-dependent"` at a small
//! `min_cm_angle_deg`) the truncation is inactive and the same foil
//! reproduces the Bohr variance to about two percent;
//! `lindhard/tests/straggling.rs` pins that configuration.
//!
//! Turning the option on shifts every subsequent transport draw of a
//! history, so a run with straggling is a different trajectory sequence
//! than a run without it, not a per-history perturbation of one; comparisons
//! between the two must be statistical over many histories. A default
//! (`"none"`) run makes no draw at all, so its results are bit-identical to
//! a build without the option. Before 2026-11-14 this module was
//! library-only: nothing sampled it (checked with
//! `grep -rnE 'variance_per_atom|StragglingModel|straggling::' --include=*.rs .`
//! piped through `grep -v 'stopping/straggling.rs'`, which printed nothing).
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
    /// Bohr with the relativistic factor `(1 - β²/2) / (1 - β²)`. The same
    /// correction appears in J. F. Ziegler, J. P. Biersack, M. D. Ziegler,
    /// *SRIM: The Stopping and Range of Ions in Matter* (SRIM Co., 2015),
    /// Eq. 7-47, p. 7-25, written there with `β = c²/v²`.
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
