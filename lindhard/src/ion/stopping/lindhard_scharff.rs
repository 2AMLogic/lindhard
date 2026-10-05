//! Lindhard-Scharff electronic stopping.
//!
//! J. Lindhard and M. Scharff, Phys. Rev. 124, 128 (1961): at ion velocities
//! below `v0 Z1^(2/3)` the electronic stopping is proportional to velocity.
//!
//! # Reduced units (Lindhard, Scharff, Schiott, Mat. Fys. Medd. 33 (14), 1963)
//!
//! * screening length `a = 0.8853 a0 (Z1^(2/3) + Z2^(2/3))^(-1/2)`
//!   (Thomas-Fermi, universal-length form; `a0` is the Bohr radius),
//! * reduced energy `ε = E a M2 / (Z1 Z2 e² (M1 + M2))`, with `e²` meaning
//!   `e²/(4π ε0)`,
//! * reduced path length `ρ = N x 4π a² M1 M2 / (M1 + M2)²`,
//! * reduced stopping `(dε/dρ)_e = k_L ε^(1/2)`, with
//!   `k_L = ξ_e 0.0793 Z1^(1/2) Z2^(1/2) (A1+A2)^(3/2)
//!   / ((Z1^(2/3)+Z2^(2/3))^(3/4) A1^(3/2) A2^(1/2))`.
//!
//! `ξ_e = Z1^(1/6)` is the Lindhard-Scharff (1961) velocity-proportionality
//! factor (`ξ_e ≈ Z1^(1/6)`, order 1 to 2). **Deviation from the issue text:**
//! the formula quoted there omits `ξ_e`; without it the Z1 dependence does not
//! agree with the dimensional LS form below (off by `Z1^(1/6)`), so it is
//! included here. The 0.0793 is the rounded coefficient; a test checks the
//! whole expression against the dimensional form
//! `S_e = Z1^(1/6) 8π e² a0 Z1 Z2 / (Z1^(2/3)+Z2^(2/3))^(3/2) · v/v0` to 1 %.
//!
//! A per-element multiplicative correction `f(Z2)` may be attached
//! ([`LindhardScharff::with_correction`]); it defaults to 1.

use super::{check_energy, target, ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::constants::{BOHR_RADIUS, COULOMB_E2};
use crate::units::J_PER_EV;
use std::f64::consts::PI;

/// Lindhard-Scharff screening length `a`, m.
pub fn screening_length(z1: u8, z2: u8) -> f64 {
    0.8853 * BOHR_RADIUS / (f64::from(z1).powf(2.0 / 3.0) + f64::from(z2).powf(2.0 / 3.0)).sqrt()
}

/// Reduced energy `ε` for ion energy `energy_ev` on a target of mass `m2_amu`.
pub fn reduced_energy(ion: &Ion, z2: u8, m2_amu: f64, energy_ev: f64) -> f64 {
    let a = screening_length(ion.z, z2);
    energy_ev * J_PER_EV * a * m2_amu
        / (f64::from(ion.z) * f64::from(z2) * COULOMB_E2 * (ion.mass_amu + m2_amu))
}

/// Factor converting reduced stopping `dε/dρ` to a cross section in J m².
pub fn reduced_to_si_factor(ion: &Ion, z2: u8, m2_amu: f64) -> f64 {
    let a = screening_length(ion.z, z2);
    4.0 * PI * a * f64::from(ion.z) * f64::from(z2) * COULOMB_E2 * ion.mass_amu
        / (ion.mass_amu + m2_amu)
}

/// The Lindhard-Scharff coefficient `k_L` including `ξ_e = Z1^(1/6)` (see
/// module docs).
pub fn k_l(ion: &Ion, z2: u8, m2_amu: f64) -> f64 {
    let (z1, z2f) = (f64::from(ion.z), f64::from(z2));
    let (a1, a2) = (ion.mass_amu, m2_amu);
    0.0793 * z1.powf(1.0 / 6.0) * z1.sqrt() * z2f.sqrt() * (a1 + a2).powf(1.5)
        / ((z1.powf(2.0 / 3.0) + z2f.powf(2.0 / 3.0)).powf(0.75) * a1.powf(1.5) * a2.sqrt())
}

/// Lindhard-Scharff model with optional per-element correction factors.
#[derive(Debug, Clone, Default)]
pub struct LindhardScharff {
    corrections: Vec<(u8, f64)>,
}

impl LindhardScharff {
    /// Plain LS (all corrections 1).
    pub fn new() -> Self {
        Self::default()
    }

    /// Multiply the stopping in target element `z2` by `factor`.
    pub fn with_correction(mut self, z2: u8, factor: f64) -> Self {
        self.corrections.retain(|(z, _)| *z != z2);
        self.corrections.push((z2, factor));
        self
    }

    fn correction(&self, z2: u8) -> f64 {
        self.corrections
            .iter()
            .find(|(z, _)| *z == z2)
            .map_or(1.0, |(_, f)| *f)
    }

    /// Reduced stopping `dε/dρ` at reduced energy `eps` (without correction).
    pub fn reduced_stopping(&self, ion: &Ion, z2: u8, m2_amu: f64, eps: f64) -> f64 {
        k_l(ion, z2, m2_amu) * eps.sqrt()
    }
}

impl ElectronicStopping for LindhardScharff {
    fn name(&self) -> &'static str {
        "lindhard-scharff"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        check_energy(energy_ev)?;
        let m2 = target(target_z)?.atomic_weight;
        let eps = reduced_energy(ion, target_z, m2, energy_ev);
        Ok(self.correction(target_z)
            * self.reduced_stopping(ion, target_z, m2, eps)
            * reduced_to_si_factor(ion, target_z, m2))
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        // v < v0 Z1^(2/3): E/A below about 25 keV Z1^(4/3).
        ValidityRange {
            min_energy_ev: 0.0,
            max_energy_ev: ion.mass_amu * super::velocity_scale_energy_per_amu_ev(ion.z, 1.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elements::element;

    #[test]
    fn reduced_stopping_is_k_sqrt_eps() {
        let ion = Ion::new(15).unwrap();
        let m = LindhardScharff::new();
        let m2 = element(14).unwrap().atomic_weight;
        let f = reduced_to_si_factor(&ion, 14, m2);
        let mut ratio0 = None;
        for e in [1.0e2, 1.0e3, 1.0e4, 1.0e5] {
            let s = m.stopping(&ion, 14, e).unwrap();
            let eps = reduced_energy(&ion, 14, m2, e);
            let r = s / f / eps.sqrt();
            let r0 = *ratio0.get_or_insert(r);
            assert!((r / r0 - 1.0).abs() < 1e-12);
            assert!((r / k_l(&ion, 14, m2) - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn matches_dimensional_form_within_one_percent() {
        for (z1, z2) in [(1u8, 14u8), (15, 14), (11, 29), (33, 14), (5, 79)] {
            let ion = Ion::new(z1).unwrap();
            let e = 2.0e3;
            let s = LindhardScharff::new().stopping(&ion, z2, e).unwrap();
            let (z1f, z2f) = (f64::from(z1), f64::from(z2));
            let v = super::super::speed_nonrel(&ion, e) / super::super::bohr_velocity();
            let dim = z1f.powf(1.0 / 6.0) * 8.0 * PI * COULOMB_E2 * BOHR_RADIUS * z1f * z2f
                / (z1f.powf(2.0 / 3.0) + z2f.powf(2.0 / 3.0)).powf(1.5)
                * v;
            assert!((s / dim - 1.0).abs() < 0.01, "{z1} {z2}: {}", s / dim);
        }
    }

    #[test]
    fn correction_factor_scales() {
        let ion = Ion::new(5).unwrap();
        let base = LindhardScharff::new().stopping(&ion, 14, 1e4).unwrap();
        let c = LindhardScharff::new()
            .with_correction(14, 1.25)
            .stopping(&ion, 14, 1e4)
            .unwrap();
        assert!((c / base - 1.25).abs() < 1e-14);
    }

    #[test]
    fn rejects_bad_energy() {
        let ion = Ion::new(5).unwrap();
        assert!(LindhardScharff::new().stopping(&ion, 14, 0.0).is_err());
        assert!(LindhardScharff::new().stopping(&ion, 0, 1.0).is_err());
    }
}
