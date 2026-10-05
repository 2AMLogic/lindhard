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
//! * `δ` and `C` default to 0 and can be set by the caller as constants if
//!   they come from a cited source. `δ` can instead be computed with the
//!   opt-in single-oscillator model, see [`density_effect_single_oscillator`]
//!   and [`BetheBloch::with_density_effect`].
//!
//! # Omitted terms
//!
//! The Barkas term `L1` (Ashley-Ritchie-Brandt 1972) needs a tabulated
//! function, and the usual shell-correction parameter sets and the
//! multi-oscillator density-effect parameter sets are published as tables
//! (ICRU 49/37 and Sternheimer-Berger-Seltzer), which this project may not
//! ingest. They are not implemented. See `docs/stopping-models.md` for the
//! literature findings; expect errors of a few percent in the range
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

/// Plasma energy `ħ ω_p = ħ sqrt(n e² / (ε0 m))` in eV for an electron
/// density `n` in m^-3 (the free-electron plasma frequency).
pub fn plasma_energy_ev(electron_density_m3: f64) -> f64 {
    let wp2 = 4.0 * PI * electron_density_m3 * COULOMB_E2 / crate::constants::ELECTRON_MASS;
    crate::constants::HBAR * wp2.sqrt() / J_PER_EV
}

/// Density-effect correction `δ` of a single-oscillator dielectric model.
///
/// Takes the squared product `(βγ)²`, the target plasma energy `ħω_p` and the
/// mean excitation energy `I` (both eV). Returns 0 below the threshold
/// `(βγ)² = (I/ħω_p)²` and otherwise
///
/// `δ = ln((βγ)² (ħω_p/I)²) − 1 + (I/ħω_p)² / (βγ)²`.
///
/// # Source and status
///
/// This is the Fermi (Phys. Rev. 57, 485 (1940)) / Sternheimer (Phys. Rev. 88,
/// 851 (1952)) dielectric formulation, which in the form given by Fano (Ann.
/// Rev. Nucl. Sci. 13, 1 (1963)) reads
/// `δ = Σ f_i ln(1 + L²/ω_i²) − L² (1−β²)/(β² ω_p²)`, with `L` fixed by
/// `Σ f_i ω_p² / (ω_i² + L²) = 1/β² − 1`. We specialise it to ONE oscillator
/// (`f = 1`, `ω_0 = I/ħ`), where the constraint solves in closed form
/// (`L² = ω_p² (βγ)² − ω_0²`), giving the expression above. The reduction was
/// done by us and is not quoted from a paper; equation numbers of the sources
/// are not cited because they were not checked. It is verified in the tests by
/// its limits: `δ → 0` continuously (value and slope) at threshold, and
/// `δ → 2 ln(βγ ħω_p/I) − 1` at large `βγ`, the standard high-energy form.
///
/// # Validity
///
/// Exact in the high-energy limit when `I` is the true mean excitation
/// energy. Near and above threshold a single oscillator is a crude stand-in
/// for the real oscillator spectrum, and it contains no conductor
/// (free-electron) terms, so the transition region is only approximate. It
/// is zero below `βγ = I/ħω_p`, which is about 5 to 6 for Si (a proton of
/// about 5 GeV). Real materials switch the effect on much earlier (around
/// `βγ ≈ 1.5` for Si in the published parameterisations), because their
/// oscillator spectrum is spread; this model therefore UNDER-estimates `δ`
/// through the transition and, in particular, is identically 0 over the whole
/// range `BetheBloch::validity` advertises (up to 1 GeV/u). Its value is the
/// correct asymptote for ultra-relativistic ions only.
pub fn density_effect_single_oscillator(
    beta_gamma_sq: f64,
    plasma_energy_ev: f64,
    mean_excitation_ev: f64,
) -> f64 {
    let a = (mean_excitation_ev / plasma_energy_ev).powi(2);
    if beta_gamma_sq <= a {
        return 0.0;
    }
    (beta_gamma_sq / a).ln() - 1.0 + a / beta_gamma_sq
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
    /// Target electron density (m^-3) for the opt-in single-oscillator
    /// density effect; `None` (default) leaves it off. Added to
    /// `density_delta`.
    pub density_effect_electron_density_m3: Option<f64>,
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

    /// Enable the single-oscillator density effect
    /// ([`density_effect_single_oscillator`]) for a target of electron density
    /// `n_e` (m^-3). The density effect needs a meaningful `I`: set a cited
    /// one with [`BetheBloch::with_mean_excitation_ev`]. Opt-in; relevant only
    /// at relativistic energies.
    pub fn with_density_effect(mut self, electron_density_m3: f64) -> Self {
        self.density_effect_electron_density_m3 = Some(electron_density_m3);
        self
    }

    /// Drop the Bloch term (pure Bethe).
    pub fn without_bloch(mut self) -> Self {
        self.no_bloch = true;
        self
    }

    /// `(gamma, beta^2, M c^2)` for kinetic energy `energy_ev`.
    ///
    /// With `tau = T / (M c^2)`, `beta^2 = tau (tau + 2) / (1 + tau)^2`, which
    /// stays accurate when `tau` is below the rounding error of `1 + tau`
    /// (the form `1 - 1/gamma^2` cancels to exactly 0 there).
    pub(super) fn kinematics(ion: &Ion, energy_ev: f64) -> (f64, f64, f64) {
        let mc2 = ion.mass_kg() * crate::constants::SPEED_OF_LIGHT.powi(2);
        let tau = energy_ev * J_PER_EV / mc2;
        let gamma = 1.0 + tau;
        let beta2 = tau * (tau + 2.0) / (gamma * gamma);
        (gamma, beta2, mc2)
    }

    /// Reject non-physical parameters. Checked on every evaluation so that
    /// direct assignment to the public fields is covered too.
    fn validate(&self) -> Result<(), StoppingError> {
        if let Some(i) = self.mean_excitation_ev {
            if !(i.is_finite() && i > 0.0) {
                return Err(StoppingError::InvalidParameter {
                    name: "mean_excitation_ev",
                    value: i,
                });
            }
        }
        if let Some(n) = self.density_effect_electron_density_m3 {
            if !(n.is_finite() && n > 0.0) {
                return Err(StoppingError::InvalidParameter {
                    name: "density_effect_electron_density_m3",
                    value: n,
                });
            }
        }
        for (name, value) in [
            ("shell_over_z", self.shell_over_z),
            ("density_delta", self.density_delta),
        ] {
            if !value.is_finite() {
                return Err(StoppingError::InvalidParameter { name, value });
            }
        }
        Ok(())
    }

    /// Bracket `[...]` of the formula (dimensionless).
    pub fn bracket(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        check_energy(energy_ev)?;
        target(target_z)?;
        self.validate()?;
        let (gamma, beta2, mc2_ion) = Self::kinematics(ion, energy_ev);
        let beta = beta2.sqrt();
        let i_j = self
            .mean_excitation_ev
            .unwrap_or(10.0 * f64::from(target_z))
            * J_PER_EV;
        let m = ELECTRON_REST_ENERGY;
        let ratio = m / mc2_ion;
        let wmax = 2.0 * m * beta2 * gamma * gamma / (1.0 + 2.0 * gamma * ratio + ratio * ratio);
        let z = self.effective_charge.charge(ion.z(), beta);
        let delta = self.density_delta
            + self.density_effect_electron_density_m3.map_or(0.0, |n| {
                density_effect_single_oscillator(
                    beta2 * gamma * gamma,
                    plasma_energy_ev(n),
                    i_j / J_PER_EV,
                )
            });
        let mut b = 0.5 * (2.0 * m * beta2 * gamma * gamma * wmax / (i_j * i_j)).ln()
            - beta2
            - 0.5 * delta
            - self.shell_over_z;
        if !self.no_bloch {
            b += bloch_term(z * FINE_STRUCTURE / beta);
        }
        if b.is_nan() {
            return Err(StoppingError::NotApplicable {
                model: "bethe-bloch",
                energy_ev,
            });
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
        let (_, beta2, _) = Self::kinematics(ion, energy_ev);
        let z = self.effective_charge.charge(ion.z(), beta2.sqrt());
        let s = 4.0 * PI * COULOMB_E2 * COULOMB_E2 * z * z * f64::from(target_z)
            / (ELECTRON_REST_ENERGY * beta2)
            * b;
        if !(b > 0.0 && s.is_finite()) {
            return Err(StoppingError::NotApplicable {
                model: "bethe-bloch",
                energy_ev,
            });
        }
        Ok(s)
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        // v >= 3 v0 Z1^(2/3) on the low side; no density effect on the high
        // side, so stop at 1 GeV/u.
        ValidityRange {
            min_energy_ev: ion.mass_amu() * super::velocity_scale_energy_per_amu_ev(ion.z(), 3.0),
            max_energy_ev: ion.mass_amu() * 1.0e9,
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
    fn plasma_energy_of_aluminium() {
        // n_e = 13 rho N_A / A, rho = 2.699 g/cm^3, A = 26.98 (elementary
        // facts); the free-electron plasma energy is about 32.9 eV.
        let n = 13.0 * 2.699e3 * AVOGADRO / 26.98e-3;
        assert!((plasma_energy_ev(n) / 32.86 - 1.0).abs() < 3e-3);
    }

    #[test]
    fn density_effect_limits() {
        let (ehp, i): (f64, f64) = (30.0, 150.0);
        let a = (i / ehp).powi(2);
        // zero at and below threshold, continuous in value and slope
        assert_eq!(density_effect_single_oscillator(0.5 * a, ehp, i), 0.0);
        let eps = 1e-4 * a;
        let above = density_effect_single_oscillator(a + eps, ehp, i);
        assert!((0.0..1e-6).contains(&above));
        // high-energy form 2 ln(bg ehp / I) - 1
        for bg2 in [1e6, 1e8] {
            let want = (bg2 * (ehp / i).powi(2)).ln() - 1.0;
            let got = density_effect_single_oscillator(bg2, ehp, i);
            assert!((got - want).abs() < 2.0 * a / bg2 + 1e-12, "{got} {want}");
        }
        // monotonic
        let mut prev = 0.0;
        for k in 0..60 {
            let d = density_effect_single_oscillator(a * 10f64.powf(f64::from(k) * 0.1), ehp, i);
            assert!(d >= prev);
            prev = d;
        }
    }

    #[test]
    fn density_effect_is_opt_in_and_reduces_stopping() {
        let n = 14.0 * 2.33e3 * AVOGADRO / 28.0855e-3; // Si
        let base = BetheBloch::new().with_mean_excitation_ev(173.0);
        let on = base.with_density_effect(n);
        let p = Ion::proton();
        // 10 MeV and 1 GeV: below this model's threshold (beta gamma about
        // 5.6 for Si), identical
        for e in [1.0e7, 1.0e9] {
            assert_eq!(
                base.stopping(&p, 14, e).unwrap(),
                on.stopping(&p, 14, e).unwrap()
            );
        }
        // 100 GeV: reduced
        let (a, b) = (
            base.stopping(&p, 14, 1.0e11).unwrap(),
            on.stopping(&p, 14, 1.0e11).unwrap(),
        );
        assert!(b < a);
        let mut bad = base;
        bad.density_effect_electron_density_m3 = Some(-1.0);
        assert!(matches!(
            bad.stopping(&p, 14, 1.0e9),
            Err(StoppingError::InvalidParameter { .. })
        ));
    }

    #[test]
    fn effective_charge_limits() {
        assert!((EffectiveCharge::BarkasEmpirical.charge(15, 0.9) / 15.0 - 1.0).abs() < 1e-6);
        assert!(EffectiveCharge::BarkasEmpirical.charge(15, 0.01) < 5.0);
    }

    #[test]
    fn tiny_energy_is_a_defined_error_not_nan() {
        let u = Ion::with_mass(92, 238.028_91).unwrap();
        for e in [1.0e-6, 1.0e-30, 1.0e-300] {
            for m in [BetheBloch::new(), BetheBloch::new().without_bloch()] {
                let r = m.stopping(&u, 14, e);
                assert!(
                    matches!(r, Err(StoppingError::NotApplicable { .. })),
                    "{e}: {r:?}"
                );
                assert!(!m.bracket(&u, 14, e).is_ok_and(f64::is_nan));
            }
        }
    }

    #[test]
    fn beta2_is_stable_at_tiny_tau() {
        let u = Ion::with_mass(92, 238.028_91).unwrap();
        let (_, beta2, _) = BetheBloch::kinematics(&u, 1.0e-6);
        assert!(beta2 > 0.0 && beta2.is_finite());
        // non-relativistic limit beta^2 = 2 T / (M c^2)
        let mc2_ev = u.mass_kg() * crate::constants::SPEED_OF_LIGHT.powi(2) / J_PER_EV;
        assert!((beta2 / (2.0e-6 / mc2_ev) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let ion = Ion::proton();
        for i in [0.0, -140.0, f64::NAN, f64::INFINITY] {
            let mut m = BetheBloch::new().without_bloch();
            m.mean_excitation_ev = Some(i);
            assert!(matches!(
                m.stopping(&ion, 14, 1.0e7),
                Err(StoppingError::InvalidParameter {
                    name: "mean_excitation_ev",
                    ..
                })
            ));
            assert!(matches!(
                BetheBloch::new()
                    .with_mean_excitation_ev(i)
                    .bracket(&ion, 14, 1.0e7),
                Err(StoppingError::InvalidParameter { .. })
            ));
        }
        for bad in [f64::NAN, f64::INFINITY] {
            let mut m = BetheBloch::new();
            m.shell_over_z = bad;
            assert!(matches!(
                m.stopping(&ion, 14, 1.0e7),
                Err(StoppingError::InvalidParameter {
                    name: "shell_over_z",
                    ..
                })
            ));
            let mut m = BetheBloch::new();
            m.density_delta = bad;
            assert!(matches!(
                m.stopping(&ion, 14, 1.0e7),
                Err(StoppingError::InvalidParameter {
                    name: "density_delta",
                    ..
                })
            ));
        }
    }

    #[test]
    fn not_applicable_at_low_energy() {
        let r = BetheBloch::new().stopping(&Ion::proton(), 14, 10.0);
        assert!(matches!(r, Err(StoppingError::NotApplicable { .. })));
    }
}
