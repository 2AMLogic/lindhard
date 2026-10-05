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
//! against the 1976 paper; see `docs/data-provenance.md`.

use super::lindhard_scharff::LindhardScharff;
use super::{ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::constants::BOHR_RADIUS;
use std::f64::consts::PI;

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
    /// `r_min_m` (metres).
    pub fn local_loss(
        &self,
        ion: &Ion,
        target_z: u8,
        energy_ev: f64,
        r_min_m: f64,
    ) -> Result<f64, StoppingError> {
        let s = self.ls.stopping(ion, target_z, energy_ev)?;
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
}
