//! Bethe-Bloch electronic stopping (relativistic) with Bloch and optional
//! user-supplied shell/density corrections, and a named effective-charge model.
//!
//! # Formula
//!
//! Per target atom (Bethe 1930/1932; Fano, Ann. Rev. Nucl. Sci. 13, 1 (1963);
//! PDG review "Passage of particles through matter", the same equation):
//!
//! `S = 4π (e²)² z² Z2 / (m c² β²) ·
//!      [ ½ ln(2 m c² β² γ² W_max / I²) − β² − δ/2 − C/Z2 + L1 + L2 ]`
//!
//! with `e²` meaning `e²/(4π ε0)`, `z` the projectile effective charge,
//! `W_max = 2 m c² β²γ² / (1 + 2γ m/M + (m/M)²)`, `I` the mean excitation
//! energy, `δ` the density-effect correction and `C` the shell correction.
//!
//! * `L2` (Bloch 1933, Ann. Phys. 408, 285): the exact series
//!   `L2 = −y² Σ_{n≥1} 1 / (n (n² + y²))`, `y = z α / β`.
//! * `I` defaults to the Bloch rule `I = 10 eV · Z2` (Bloch 1933). This is a
//!   rough estimate, good to about 20 % in `I` (about 2 % in `S`); pass a
//!   measured `I` (cited) with [`BetheBloch::with_mean_excitation_ev`].
//! * `δ` and `C` default to 0 (see below) and can be set by the caller as
//!   constants if they come from a cited source.
//!
//! # Omitted terms
//!
//! The Barkas term `L1` (Ashley-Ritchie-Brandt 1972) needs a tabulated
//! function, and the usual shell-correction and density-effect parameter sets
//! are published as tables (ICRU 49/37 and Sternheimer-Berger-Seltzer), which
//! this project may not ingest. They are not implemented. See
//! `docs/stopping-models.md`; expect errors of a few percent in the range
//! 1 to 10 MeV/u until they are added.
//!
//! # Effective charge
//!
//! [`EffectiveCharge::BarkasEmpirical`]: `z_eff = Z1 [1 − exp(−125 β / Z1^(2/3))]`
//! (W. H. Barkas, *Nuclear Research Emulsions* I, Academic Press (1963);
//! quoted in Northcliffe, Ann. Rev. Nucl. Sci. 13, 67 (1963)).

use super::{check_energy, target, ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::constants::{COULOMB_E2, ELECTRON_REST_ENERGY, FINE_STRUCTURE};
use crate::units::J_PER_EV;
use std::f64::consts::PI;

/// Projectile effective-charge model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EffectiveCharge {
    /// `z = Z1` (fully stripped); right for protons and alphas at high energy.
    #[default]
    Bare,
    /// Barkas empirical form, see the module docs.
    BarkasEmpirical,
}

impl EffectiveCharge {
    /// Effective charge number at speed `beta`.
    pub fn charge(self, z1: u8, beta: f64) -> f64 {
        let z = f64::from(z1);
        match self {
            Self::Bare => z,
            Self::BarkasEmpirical => z * (1.0 - (-125.0 * beta / z.powf(2.0 / 3.0)).exp()),
        }
    }
}

/// Bloch correction `L2(y)`, `y = z α / β` (see module docs).
pub fn bloch_term(y: f64) -> f64 {
    let y2 = y * y;
    const N: u32 = 20_000;
    let mut sum = 0.0;
    // Sum small terms first for accuracy.
    for n in (1..=N).rev() {
        let n = f64::from(n);
        sum += 1.0 / (n * (n * n + y2));
    }
    // Tail: sum_{n>N} 1/(n^3) ~ 1/(2 N^2) (y << N).
    let nn = f64::from(N);
    sum += 1.0 / (2.0 * nn * nn);
    -y2 * sum
}

/// Bethe-Bloch model.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BetheBloch {
    /// Effective-charge model.
    pub effective_charge: EffectiveCharge,
    /// Mean excitation energy override, eV.
    pub mean_excitation_ev: Option<f64>,
    /// Shell correction `C / Z2` (dimensionless, added as `−C/Z2`); default 0.
    pub shell_over_z: f64,
    /// Density-effect `δ` (dimensionless); default 0.
    pub density_delta: f64,
    /// Include the Bloch term `L2` (default: true via [`BetheBloch::new`]).
    pub no_bloch: bool,
}

impl BetheBloch {
    /// Bare-charge Bethe-Bloch with Bloch term, Bloch-rule `I`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the effective-charge model.
    pub fn with_effective_charge(mut self, m: EffectiveCharge) -> Self {
        self.effective_charge = m;
        self
    }

    /// Set the mean excitation energy, eV (cite its source at the call site).
    pub fn with_mean_excitation_ev(mut self, i_ev: f64) -> Self {
        self.mean_excitation_ev = Some(i_ev);
        self
    }

    /// Drop the Bloch term (pure Bethe).
    pub fn without_bloch(mut self) -> Self {
        self.no_bloch = true;
        self
    }

    fn kinematics(ion: &Ion, energy_ev: f64) -> (f64, f64, f64) {
        let mc2 = ion.mass_kg() * crate::constants::SPEED_OF_LIGHT.powi(2);
        let gamma = 1.0 + energy_ev * J_PER_EV / mc2;
        let beta2 = 1.0 - 1.0 / (gamma * gamma);
        (gamma, beta2, mc2)
    }

    /// Bracket `[...]` of the formula (dimensionless).
    pub fn bracket(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        check_energy(energy_ev)?;
        target(target_z)?;
        let (gamma, beta2, mc2_ion) = Self::kinematics(ion, energy_ev);
        let beta = beta2.sqrt();
        let i_j = self
            .mean_excitation_ev
            .unwrap_or(10.0 * f64::from(target_z))
            * J_PER_EV;
        let m = ELECTRON_REST_ENERGY;
        let ratio = m / mc2_ion;
        let wmax = 2.0 * m * beta2 * gamma * gamma / (1.0 + 2.0 * gamma * ratio + ratio * ratio);
        let z = self.effective_charge.charge(ion.z, beta);
        let mut b = 0.5 * (2.0 * m * beta2 * gamma * gamma * wmax / (i_j * i_j)).ln()
            - beta2
            - 0.5 * self.density_delta
            - self.shell_over_z;
        if !self.no_bloch {
            b += bloch_term(z * FINE_STRUCTURE / beta);
        }
        Ok(b)
    }
}

impl ElectronicStopping for BetheBloch {
    fn name(&self) -> &'static str {
        "bethe-bloch"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        let b = self.bracket(ion, target_z, energy_ev)?;
        if b <= 0.0 {
            return Err(StoppingError::NotApplicable {
                model: "bethe-bloch",
                energy_ev,
            });
        }
        let (_, beta2, _) = Self::kinematics(ion, energy_ev);
        let z = self.effective_charge.charge(ion.z, beta2.sqrt());
        Ok(
            4.0 * PI * COULOMB_E2 * COULOMB_E2 * z * z * f64::from(target_z)
                / (ELECTRON_REST_ENERGY * beta2)
                * b,
        )
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        // v >= 3 v0 Z1^(2/3) on the low side; no density effect on the high
        // side, so stop at 1 GeV/u.
        ValidityRange {
            min_energy_ev: ion.mass_amu * super::velocity_scale_energy_per_amu_ev(ion.z, 3.0),
            max_energy_ev: ion.mass_amu * 1.0e9,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{AVOGADRO, BOHR_RADIUS, PROTON_MASS, SPEED_OF_LIGHT};

    /// Hand calculation for protons in Si (Z2 = 14, I = 10 eV * 14 = 140 eV by
    /// the Bloch rule), no Bloch/shell/density terms:
    ///
    ///   gamma = 1 + T / (m_p c^2); beta^2 = 1 - 1/gamma^2
    ///   S = 4 pi (e^2)^2 Z2 / (m c^2 beta^2) [ ln(2 m c^2 beta^2 gamma^2 / I) - beta^2 ]
    ///
    /// with W_max replaced by its heavy-particle limit. The expected values in
    /// eV 1e-15 cm^2 were computed independently from this arithmetic with the
    /// CODATA constants (T = 1, 10, 100 MeV).
    #[test]
    fn proton_in_si_hand_calculation() {
        let m = BetheBloch::new()
            .with_mean_excitation_ev(140.0)
            .without_bloch();
        let ion = Ion::proton();
        let expected = [
            (1.0e6, 9.201_572_04),
            (1.0e7, 1.712_108_53),
            (1.0e8, 0.281_049_35),
        ];
        for (e, want) in expected {
            let got = super::super::to_ev_1e15_cm2(m.stopping(&ion, 14, e).unwrap());
            // W_max exact vs heavy limit differs by < 1e-4 in S.
            assert!((got / want - 1.0).abs() < 2.0e-4, "{e}: {got} vs {want}");
        }
    }

    #[test]
    fn independent_formula_check() {
        // Same hand formula, written in MeV cm^2 / g via K = 4 pi N_A r_e^2 m c^2.
        let re = FINE_STRUCTURE * FINE_STRUCTURE * BOHR_RADIUS;
        let mc2_mev = ELECTRON_REST_ENERGY / J_PER_EV / 1e6;
        let k = 4.0 * PI * AVOGADRO * (re * 100.0).powi(2) * mc2_mev; // MeV cm^2/mol
        let t = 10.0;
        let g = 1.0 + t / (PROTON_MASS * SPEED_OF_LIGHT.powi(2) / J_PER_EV / 1e6);
        let b2 = 1.0 - 1.0 / (g * g);
        let l = (2.0 * mc2_mev * 1e6 * b2 * g * g / 140.0).ln() - b2;
        let per_g = k * 14.0 / 28.0855 / b2 * l; // MeV cm^2 / g
        let m = BetheBloch::new()
            .with_mean_excitation_ev(140.0)
            .without_bloch();
        let s = super::super::to_ev_1e15_cm2(m.stopping(&Ion::proton(), 14, t * 1e6).unwrap());
        let per_g_model = s * 1e-15 * 1e-6 * AVOGADRO / 28.0855;
        assert!((per_g_model / per_g - 1.0).abs() < 2e-4);
    }

    #[test]
    fn bloch_term_small_y_matches_zeta3() {
        let y = 0.01;
        assert!((bloch_term(y) / (-1.202_056_903_16 * y * y) - 1.0).abs() < 1e-3);
    }

    #[test]
    fn bloch_reduces_stopping() {
        let ion = Ion::new(6).unwrap();
        let a = BetheBloch::new().stopping(&ion, 14, 1.2e8).unwrap();
        let b = BetheBloch::new()
            .without_bloch()
            .stopping(&ion, 14, 1.2e8)
            .unwrap();
        assert!(a < b);
    }

    #[test]
    fn effective_charge_limits() {
        assert!((EffectiveCharge::BarkasEmpirical.charge(15, 0.9) / 15.0 - 1.0).abs() < 1e-6);
        assert!(EffectiveCharge::BarkasEmpirical.charge(15, 0.01) < 5.0);
    }

    #[test]
    fn not_applicable_at_low_energy() {
        let r = BetheBloch::new().stopping(&Ion::proton(), 14, 10.0);
        assert!(matches!(r, Err(StoppingError::NotApplicable { .. })));
    }
}
