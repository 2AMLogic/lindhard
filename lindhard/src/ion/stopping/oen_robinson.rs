//! Oen-Robinson local (impact-parameter dependent) electronic energy loss.
//!
//! O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132, 647 (1976).
//!
//! The nonlocal part is Lindhard-Scharff with the same `k_L`. The local model
//! assigns a collision with closest approach `r_min` the loss
//!
//! `ΔE_e(E, r_min) = S_LS(E) · c² / (2π a²) · exp(-c r_min / a)`, `c = 0.3`,
//!
//! with `a` the Firsov screening length `0.8853 a0 (Z1^(1/2) + Z2^(1/2))^(-2/3)`.
//! The prefactor normalises `∫ ΔE_e 2π p dp = S_LS` when `r_min` is identified
//! with the impact parameter `p`; that identification is the approximation made
//! in the original paper, so the impact-averaged stopping equals the LS value.
//!
//! The constants `c = 0.3` and the choice of length were entered from the
//! literature as remembered by the contributor and have **not** been verified
//! against the 1976 paper; see `docs/data-provenance.md` and
//! [`OR_CONSTANTS_UNVERIFIED`].

use super::lindhard_scharff::LindhardScharff;
use super::{ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::constants::BOHR_RADIUS;
use std::f64::consts::PI;

/// Whether the constants of the local model ([`OR_EXPONENT_COEFFICIENT`] and
/// the length of [`local_length`]) are still **not verified** against Oen and
/// Robinson, Nucl. Instrum. Methods 132, 647 (1976). `true` while the
/// Oen-Robinson row of `docs/data-provenance.md` says "not verified".
///
/// This is the one place to change when the constants are checked against the
/// paper: run metadata that reports the flag
/// ([`crate::ion::bca::CrystalMetadata::electronic_constants_unverified`])
/// follows it.
pub const OR_CONSTANTS_UNVERIFIED: bool = true;

/// Exponent coefficient `c` in `exp(-c r_min / a)`.
pub const OR_EXPONENT_COEFFICIENT: f64 = 0.3;

/// Firsov-type screening length used by the local model, m.
pub fn local_length(z1: u8, z2: u8) -> f64 {
    0.8853 * BOHR_RADIUS / (f64::from(z1).sqrt() + f64::from(z2).sqrt()).powf(2.0 / 3.0)
}

/// Oen-Robinson model. `stopping` is the impact-averaged value; use
/// [`OenRobinson::local_loss`] for the per-collision loss.
#[derive(Debug, Clone, Default)]
pub struct OenRobinson {
    ls: LindhardScharff,
}

impl OenRobinson {
    /// New model with plain LS nonlocal coefficient.
    pub fn new() -> Self {
        Self::default()
    }

    /// Local electronic energy loss, J, for a collision with closest approach
    /// `r_min_m` (metres), which must be finite and non-negative.
    pub fn local_loss(
        &self,
        ion: &Ion,
        target_z: u8,
        energy_ev: f64,
        r_min_m: f64,
    ) -> Result<f64, StoppingError> {
        let s = self.ls.stopping(ion, target_z, energy_ev)?;
        if !(r_min_m.is_finite() && r_min_m >= 0.0) {
            return Err(StoppingError::InvalidParameter {
                name: "r_min_m",
                value: r_min_m,
            });
        }
        let a = local_length(ion.z(), target_z);
        let c = OR_EXPONENT_COEFFICIENT;
        Ok(s * c * c / (2.0 * PI * a * a) * (-c * r_min_m / a).exp())
    }
}

impl ElectronicStopping for OenRobinson {
    fn name(&self) -> &'static str {
        "oen-robinson"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        self.ls.stopping(ion, target_z, energy_ev)
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        self.ls.validity(ion)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn impact_average_recovers_stopping() {
        // Integrate local_loss * 2 pi p dp over p numerically (r_min = p).
        let ion = Ion::new(15).unwrap();
        let m = OenRobinson::new();
        let s = m.stopping(&ion, 14, 1.0e4).unwrap();
        let a = local_length(15, 14);
        let (n, pmax) = (200_000, 60.0 * a);
        let dp = pmax / f64::from(n);
        let mut sum = 0.0;
        for i in 0..n {
            let p = (f64::from(i) + 0.5) * dp;
            sum += m.local_loss(&ion, 14, 1.0e4, p).unwrap() * 2.0 * PI * p * dp;
        }
        assert!((sum / s - 1.0).abs() < 1e-4, "{}", sum / s);
    }

    #[test]
    fn loss_decreases_with_r_min() {
        let ion = Ion::new(15).unwrap();
        let m = OenRobinson::new();
        let a = m.local_loss(&ion, 14, 1e4, 1e-11).unwrap();
        let b = m.local_loss(&ion, 14, 1e4, 1e-10).unwrap();
        assert!(a > b && b > 0.0);
    }

    #[test]
    fn rejects_bad_closest_approach() {
        let ion = Ion::proton();
        let m = OenRobinson::new();
        for bad in [-1e-10, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                m.local_loss(&ion, 14, 1e4, bad),
                Err(StoppingError::InvalidParameter {
                    name: "r_min_m",
                    ..
                })
            ));
        }
    }

    #[test]
    fn zero_closest_approach_is_the_maximum_loss() {
        let ion = Ion::proton();
        let m = OenRobinson::new();
        let s = m.stopping(&ion, 14, 1e4).unwrap();
        let a = local_length(1, 14);
        let c = OR_EXPONENT_COEFFICIENT;
        let l0 = m.local_loss(&ion, 14, 1e4, 0.0).unwrap();
        assert!((l0 / (s * c * c / (2.0 * PI * a * a)) - 1.0).abs() < 1e-14);
    }
}
