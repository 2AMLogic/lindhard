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
//! The ion BCA engine samples electronic straggling when it is enabled with
//! `[physics] straggling = "bohr"` (or `"bohr-relativistic"`; input key
//! `StragglingChoice`, engine setting `BcaConfig::straggling`). Every free
//! flight still takes its **mean** nonlocal loss `N S_e s` unchanged, so
//! layer attribution and deposition profiles keep their deterministic
//! backbone; the flight also adds to that particle's pending straggling
//! segment both its mean loss and its Bohr variance `n Omega^2 s`, with the
//! variance Bragg-added over the flight's own material
//! (`variance_per_atom_material`). Because Bohr's variance is linear in path
//! length, the variance accumulated over any stretch of flights is exactly
//! the Bohr variance of that stretch. When the accumulated mean loss has
//! grown to `k = 4` times the accumulated standard deviation, the segment
//! closes and **one** deviation is drawn from a Gaussian of the segment's
//! accumulated standard deviation, from the history's own transport stream
//! (so results stay deterministic per `(seed, index)` at any thread count).
//! At `k = 4` a standard normal draw falls below zero with probability
//! `Phi(-4) ~ 3.2e-5` (about one segment in 31,000), so the safety clamp of
//! the draw to `[-mu_acc, E]` — which keeps the segment's sampled loss
//! non-negative and within the ion's energy — is essentially inactive, and
//! the mean loss stays unbiased to that truncation probability.
//!
//! The segment, not the flight, is the sampling unit because Bohr's variance
//! is a many-collisions result: a Gaussian describes the loss over a path
//! long enough that it sits many standard deviations above zero, which a
//! single free flight is not (the constant free path flies `lambda =
//! N^(-1/3)`, 0.21 nm in carbon, where a flight's mean loss is a small
//! fraction of its Bohr standard deviation; a crystal region flies the
//! lattice, `search_length_nm` defaulting to the lattice constant, 0.543 nm
//! in Si, and has no long-flight option). Sampling per flight and truncating
//! each draw to `[0, E]` left-censored the Gaussian and inflated the mean
//! electronic loss — measured, before 2026-11-15, at about 2.5x on a 1 MeV
//! proton through 100 nm of carbon with the constant free path (variance
//! depressed to about 0.42 of the Bohr value) and about 3.7x through crystal
//! Si on lattice flights; the per-segment cadence removes that bias, which
//! is why those numbers are stated as history, not as current behaviour.
//!
//! A segment also closes — with the same draw and clamp — when the
//! particle's transport ends: it stops below its cutoff, escapes through any
//! face, or reflects at a surface. A partial segment closed at the end of a
//! transport can sit closer to zero than `k` standard deviations, so there
//! the clamp is not necessarily inactive: for a history made of many
//! segments it affects only the last, partial one, while for a history that
//! *is* a single partial segment — an ion stopped in a layer thinner than
//! the closing distance — it affects the history's only draw and biases its
//! mean loss upward. Layer boundaries do **not** close a segment:
//! each flight already contributes its own material's variance, so a
//! segment spanning a boundary carries Bohr's total variance across it. The
//! closing deviation is tallied on the flight that closes the segment — the
//! fluctuation is localized at the segment's end rather than spread along
//! it — a second-order effect on the deposition profile next to the mean
//! bias the cadence removes. The deviation draw changes the draw cadence of
//! a straggling run, so its numbers move relative to per-flight sampling;
//! the option is opt-in, so no released behaviour is affected.
//!
//! The effective charge is `EffectiveCharge::Bare`, the form in which the
//! Bohr expression is written above; the engine's mean loss applies no
//! charge correction either. The variance stays pinned to that charge
//! regardless of the caller's stopping model: a `bethe-bloch` mean loss
//! configured with a non-bare effective charge
//! (`BetheBloch::effective_charge`) samples against a variance
//! computed for a bare projectile, so the pairing is not self-consistent —
//! the mean loss carries the charge correction, the variance does not.
//! That is a known limitation, recorded here rather than fixed: threading
//! the mean loss's charge into the variance would change sampled physics.
//! Because the bare-charge variance does not depend on energy, the BCA
//! engine hoists it per layer beside its mean-loss coefficient and reads
//! the precomputed value while accumulating a segment's variance in the
//! amorphous free-flight path; the relativistic variant keeps its beta²
//! dependence and is evaluated per flight.
//!
//! Turning the option on shifts every subsequent transport draw of a
//! history, so a run with straggling is a different trajectory sequence
//! than a run without it, not a per-history perturbation of one;
//! comparisons between the two must be statistical over many histories. A
//! default (`"none"`) run makes no draw at all, so its results are
//! bit-identical to a build without the option. Before 2026-11-14 this
//! module was library-only: nothing sampled it (checked with
//! `grep -rnE 'variance_per_atom|StragglingModel|straggling::' --include=*.rs .`
//! piped through `grep -v 'stopping/straggling.rs'`, which printed nothing).
//!
//! `lindhard/tests/straggling.rs` pins the sampling: the transmitted-energy
//! variance of a 1 MeV proton through 100 nm of carbon recovers the Bohr
//! value with the constant free path (as with long flights), the mean
//! nonlocal loss agrees with a `straggling = "none"` run in amorphous C, Si
//! and a crystal Si foil, and energy conservation and thread-count
//! determinism hold with the new cadence.
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
