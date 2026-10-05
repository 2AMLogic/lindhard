//! Interpolation joining the low-energy and Bethe regimes.
//!
//! Biersack and Haggmark, Nucl. Instrum. Methods 174, 257 (1980), and Ziegler,
//! Biersack and Littmark, *The Stopping and Range of Ions in Solids* (1985),
//! join a low-velocity stopping `S_low` and a high-velocity stopping `S_high`
//! by the harmonic (reciprocal-sum) form
//!
//! `1/S = 1/S_low + 1/S_high`.
//!
//! That is the form used here. Where those works fit `S_low` and `S_high` with
//! tabulated per-element coefficients (Tier C data, not used), this module
//! uses formulas only:
//!
//! * `S_low` = Lindhard-Scharff ([`super::lindhard_scharff`]);
//! * `S_high` = the Bethe-Bloch prefactor times the softened bracket
//!   `ln(1 + e^B)`, with `B` the Bethe-Bloch bracket
//!   ([`super::bethe::BetheBloch::bracket`]). `ln(1 + e^B) → B` for `B ≫ 1`, so
//!   it equals Bethe-Bloch at high energy, and stays positive where the Bethe
//!   bracket changes sign. This regularisation is our own device in the spirit
//!   of the `ln(1 + B/E + CE)` form of the literature; it carries no fitted
//!   coefficient.
//!
//! The default high-energy part uses the bare charge `z = Z1` and omits the
//! Bloch term. Both choices keep `S_high` bounded below as `E → 0` (it tends to
//! a constant), which makes the joined result tend to LS there. The Barkas
//! effective charge (it scales `S_high` down like `E` at low energy) and the
//! Bloch term (it makes `S_high` fall like `√E`, comparable to the LS branch
//! and so spoiling the LS limit) are available through [`BiersackVarelas::new`].
//! For heavy ions the default is an overestimate of the true stopping between
//! about 0.1 and 1 MeV/u, where the projectile is not yet stripped.
//!
//! The result is smooth (no crossover switch), tends to `S_low` at low energy
//! and to Bethe-Bloch at high energy, and has a single maximum (tested for the
//! cases in the unit tests).

use super::bethe::BetheBloch;
use super::lindhard_scharff::LindhardScharff;
use super::{check_energy, target, ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::constants::{COULOMB_E2, ELECTRON_REST_ENERGY};
use std::f64::consts::PI;

/// Harmonic joining of Lindhard-Scharff and softened Bethe-Bloch.
#[derive(Debug, Clone)]
pub struct BiersackVarelas {
    low: LindhardScharff,
    high: BetheBloch,
}

impl Default for BiersackVarelas {
    fn default() -> Self {
        Self {
            low: LindhardScharff::new(),
            high: BetheBloch::new().without_bloch(),
        }
    }
}

impl BiersackVarelas {
    /// Join the given low-energy model and Bethe-Bloch configuration.
    pub fn new(low: LindhardScharff, high: BetheBloch) -> Self {
        Self { low, high }
    }

    /// The softened high-energy stopping, J m².
    pub fn high_energy_part(
        &self,
        ion: &Ion,
        target_z: u8,
        energy_ev: f64,
    ) -> Result<f64, StoppingError> {
        check_energy(energy_ev)?;
        target(target_z)?;
        let b = self.high.bracket(ion, target_z, energy_ev)?;
        // softplus, written to avoid overflow for large b.
        let soft = if b > 30.0 { b } else { b.exp().ln_1p() };
        let mc2_ion = ion.mass_kg() * crate::constants::SPEED_OF_LIGHT.powi(2);
        let gamma = 1.0 + energy_ev * crate::units::J_PER_EV / mc2_ion;
        let beta2 = 1.0 - 1.0 / (gamma * gamma);
        let z = self.high.effective_charge.charge(ion.z, beta2.sqrt());
        Ok(
            4.0 * PI * COULOMB_E2 * COULOMB_E2 * z * z * f64::from(target_z)
                / (ELECTRON_REST_ENERGY * beta2)
                * soft,
        )
    }
}

impl ElectronicStopping for BiersackVarelas {
    fn name(&self) -> &'static str {
        "biersack-varelas"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        let lo = self.low.stopping(ion, target_z, energy_ev)?;
        let hi = self.high_energy_part(ion, target_z, energy_ev)?;
        Ok(lo * hi / (lo + hi))
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        ValidityRange {
            min_energy_ev: 0.0,
            max_energy_ev: ion.mass_amu * 1.0e9,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> Vec<f64> {
        // 10 eV .. 1 GeV, 40 points per decade.
        (0..=320)
            .map(|i| 10f64.powf(1.0 + f64::from(i) / 40.0))
            .collect()
    }

    fn check_unimodal_and_continuous(ion: Ion, z2: u8) {
        let m = BiersackVarelas::default();
        let s: Vec<f64> = grid()
            .iter()
            .map(|&e| m.stopping(&ion, z2, e).unwrap())
            .collect();
        let mut peak = 0;
        for i in 1..s.len() {
            if s[i] > s[peak] {
                peak = i;
            }
        }
        assert!(peak > 0 && peak < s.len() - 1, "peak at edge");
        for i in 1..=peak {
            assert!(s[i] > s[i - 1], "not increasing at {i}");
        }
        for i in peak + 1..s.len() {
            assert!(s[i] < s[i - 1], "not decreasing at {i}");
        }
        // continuity: 40 points/decade => neighbouring ratio bounded
        for i in 1..s.len() {
            assert!((s[i] / s[i - 1] - 1.0).abs() < 0.2, "jump at {i}");
        }
    }

    #[test]
    fn proton_in_si_smooth_and_unimodal() {
        check_unimodal_and_continuous(Ion::proton(), 14);
    }

    #[test]
    fn phosphorus_in_si_smooth_and_unimodal() {
        check_unimodal_and_continuous(Ion::new(15).unwrap(), 14);
    }

    #[test]
    fn limits() {
        let m = BiersackVarelas::default();
        let ion = Ion::proton();
        let e = 1.0e1;
        let ls = LindhardScharff::new().stopping(&ion, 14, e).unwrap();
        let r = m.stopping(&ion, 14, e).unwrap() / ls;
        assert!((r - 1.0).abs() < 0.05, "{r}");
        let e = 1.0e9;
        let bb = BetheBloch::new()
            .without_bloch()
            .stopping(&ion, 14, e)
            .unwrap();
        assert!((m.stopping(&ion, 14, e).unwrap() / bb - 1.0).abs() < 0.05);
    }
}
