//! An analytic Drude-Lorentz energy-loss function with closed-form sum rules.
//!
//! This is a **synthetic** model ELF for tests, examples and benchmarks: its
//! parameters are chosen by the caller, not fitted to any material, and its
//! sum rules are known exactly, so a numerical sum-rule evaluation can be
//! checked against them.
//!
//! # Model
//!
//! Each oscillator `i` has a strength `A_i` (eV²), a resonance energy `E_i`
//! (eV) and a width `γ_i` (eV). The inverse dielectric function is
//!
//! ```text
//! 1/ε(W) = 1 + Σ_i A_i / (W² - E_i² + i γ_i W),
//! ```
//!
//! so the energy-loss function is
//!
//! ```text
//! ELF(W) = Im[-1/ε(W)] = Σ_i A_i γ_i W / ((W² - E_i²)² + γ_i² W²),   W > 0.
//! ```
//!
//! One term with `A = E_p²` and `E_i = E_p` is the Mermin-type optical ELF of
//! Shinotsuka et al., Surf. Interface Anal. 49, 238 (2017),
//! doi:10.1002/sia.6123, eq. (16) (`a_i = 1`); there it describes a free-electron
//! plasmon of energy `E_p` and width `γ`.
//!
//! # Closed-form sum rules
//!
//! The two sum rules of the same paper, eqs. (18) and (19), are, for the full
//! range `0 < W < ∞` and in energy units,
//!
//! ```text
//! ∫_0^∞ W ELF(W) dW         = (π/2) Σ_i A_i           (f-sum)
//! (2/π) ∫_0^∞ ELF(W) / W dW = Σ_i A_i / E_i²          (P_eff, without the n(0)⁻² term)
//! ```
//!
//! Both follow from the two elementary integrals (by residues: the integrands
//! are even in `W`, and `(W² - E²)² + γ² W² = |W² - E² + iγW|²` vanishes at
//! `W = ±sqrt(E² - γ²/4) ± iγ/2`)
//!
//! ```text
//! ∫_0^∞ γ W² / ((W² - E²)² + γ² W²) dW = π/2,
//! ∫_0^∞ γ    / ((W² - E²)² + γ² W²) dW = π / (2 E²).
//! ```
//!
//! The derivation is ours; no external source for these two integrals was
//! consulted. The unit tests check both against an independent numerical
//! quadrature. With `A = E_p²` and `E_i = E_p` (a free-electron plasmon), the
//! f-sum gives the plasma energy `E_p` back and `P_eff = 1`, the values the
//! paper states for a complete ELF (`P_eff -> 1` as `ΔE_max -> ∞`, text below
//! eq. (19)).

use crate::electron::data::{ElectronDataError, OpticalElf};

/// One Drude-Lorentz oscillator (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DrudeLorentzOscillator {
    /// Strength `A`, eV². `A = E²` is a free-electron plasmon of energy `E`.
    pub strength_ev2: f64,
    /// Resonance energy `E`, eV.
    pub energy_ev: f64,
    /// Width `γ`, eV.
    pub width_ev: f64,
}

impl DrudeLorentzOscillator {
    /// A free-electron plasmon of energy `plasmon_ev` and width `width_ev`:
    /// strength `plasmon_ev²`, resonance at `plasmon_ev`.
    pub fn plasmon(plasmon_ev: f64, width_ev: f64) -> Self {
        Self {
            strength_ev2: plasmon_ev * plasmon_ev,
            energy_ev: plasmon_ev,
            width_ev,
        }
    }

    fn elf(&self, w: f64) -> f64 {
        let d = w * w - self.energy_ev * self.energy_ev;
        self.strength_ev2 * self.width_ev * w / (d * d + self.width_ev * self.width_ev * w * w)
    }
}

/// A sum of Drude-Lorentz oscillators: an analytic, synthetic ELF with
/// closed-form sum rules (see the module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct DrudeLorentz {
    oscillators: Vec<DrudeLorentzOscillator>,
}

impl DrudeLorentz {
    /// At least one oscillator; every strength, energy and width finite and
    /// positive.
    pub fn new(oscillators: Vec<DrudeLorentzOscillator>) -> Result<Self, ElectronDataError> {
        if oscillators.is_empty() {
            return Err(ElectronDataError::Invalid {
                what: "Drude-Lorentz model",
                reason: "needs at least one oscillator".into(),
            });
        }
        for o in &oscillators {
            for (name, v) in [
                ("strength", o.strength_ev2),
                ("energy", o.energy_ev),
                ("width", o.width_ev),
            ] {
                if !(v.is_finite() && v > 0.0) {
                    return Err(ElectronDataError::Invalid {
                        what: "Drude-Lorentz oscillator",
                        reason: format!("{name} must be finite and positive, got {v}"),
                    });
                }
            }
        }
        Ok(Self { oscillators })
    }

    /// The oscillators.
    pub fn oscillators(&self) -> &[DrudeLorentzOscillator] {
        &self.oscillators
    }

    /// `ELF(W)` at energy loss `W` (eV); zero for `W <= 0`.
    pub fn elf(&self, energy_ev: f64) -> f64 {
        if energy_ev <= 0.0 {
            return 0.0;
        }
        self.oscillators.iter().map(|o| o.elf(energy_ev)).sum()
    }

    /// `∫_0^∞ W ELF(W) dW = (π/2) Σ A_i`, eV² (closed form).
    pub fn f_sum_ev2(&self) -> f64 {
        std::f64::consts::FRAC_PI_2 * self.oscillators.iter().map(|o| o.strength_ev2).sum::<f64>()
    }

    /// The plasma energy `ħΩ_p = sqrt(Σ A_i)`, eV, defined by the f-sum
    /// `∫ W ELF dW = (π/2) (ħΩ_p)²` (closed form).
    pub fn plasma_energy_ev(&self) -> f64 {
        self.oscillators
            .iter()
            .map(|o| o.strength_ev2)
            .sum::<f64>()
            .sqrt()
    }

    /// `P_eff = (2/π) ∫_0^∞ ELF(W)/W dW = Σ A_i / E_i²` (closed form).
    pub fn p_eff(&self) -> f64 {
        self.oscillators
            .iter()
            .map(|o| o.strength_ev2 / (o.energy_ev * o.energy_ev))
            .sum()
    }

    /// Sample the model on `n` log-spaced energies from `lo_ev` to `hi_ev`
    /// into an [`OpticalElf`] (interpolated linearly between the samples).
    /// The provenance says that the table is synthetic and lists the
    /// parameters.
    pub fn to_optical_elf(
        &self,
        material: &str,
        lo_ev: f64,
        hi_ev: f64,
        n: usize,
    ) -> Result<OpticalElf, ElectronDataError> {
        if !(lo_ev.is_finite() && hi_ev.is_finite() && lo_ev > 0.0 && hi_ev > lo_ev && n >= 2) {
            return Err(ElectronDataError::Invalid {
                what: "Drude-Lorentz sampling grid",
                reason: format!("need 0 < lo < hi and n >= 2, got {lo_ev}, {hi_ev}, {n}"),
            });
        }
        let ratio = hi_ev / lo_ev;
        let energy: Vec<f64> = (0..n)
            .map(|i| {
                if i == n - 1 {
                    hi_ev
                } else {
                    lo_ev * ratio.powf(i as f64 / (n - 1) as f64)
                }
            })
            .collect();
        let elf = energy.iter().map(|&w| self.elf(w)).collect();
        let params: Vec<String> = self
            .oscillators
            .iter()
            .map(|o| {
                format!(
                    "(A = {} eV^2, E = {} eV, gamma = {} eV)",
                    o.strength_ev2, o.energy_ev, o.width_ev
                )
            })
            .collect();
        let provenance = format!(
            "synthetic Drude-Lorentz ELF, not physical data: {}; {n} log-spaced samples from {lo_ev} to {hi_ev} eV",
            params.join(", ")
        );
        OpticalElf::new(material, provenance, energy, elf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electron::inelastic::quadrature::{integrate_segments, GaussLegendre};

    /// `∫_0^∞ f(W) dW` by the substitution `W = e^u`, `u` in [-30, 30] (the
    /// integrands below fall off as `W` and `W^-2` or faster at the ends).
    fn numeric(f: impl Fn(f64) -> f64, peak: f64) -> f64 {
        let gl = GaussLegendre::new(10);
        let breaks: Vec<f64> = [-30.0, peak.ln() - 1.0, peak.ln(), peak.ln() + 1.0, 30.0].to_vec();
        integrate_segments(&gl, &mut |u: f64| [f(u.exp()) * u.exp()], &breaks, 1e-12)[0]
    }

    #[test]
    fn closed_form_integrals_match_quadrature() {
        for (e, g) in [(20.0, 5.0), (16.7, 0.5), (5.0, 12.0), (100.0, 300.0)] {
            let den = |w: f64| (w * w - e * e).powi(2) + g * g * w * w;
            let first = numeric(|w| g * w * w / den(w), e);
            let second = numeric(|w| g / den(w), e);
            assert!(
                (first / std::f64::consts::FRAC_PI_2 - 1.0).abs() < 1e-10,
                "{e} {g}: {first}"
            );
            assert!(
                (second / (std::f64::consts::PI / (2.0 * e * e)) - 1.0).abs() < 1e-10,
                "{e} {g}: {second}"
            );
        }
    }

    #[test]
    fn sum_rules_of_a_two_oscillator_model() {
        let m = DrudeLorentz::new(vec![
            DrudeLorentzOscillator::plasmon(16.0, 4.0),
            DrudeLorentzOscillator {
                strength_ev2: 300.0,
                energy_ev: 110.0,
                width_ev: 60.0,
            },
        ])
        .unwrap();
        let f = numeric(|w| w * m.elf(w), 16.0);
        let p = numeric(|w| m.elf(w) / w, 16.0) * 2.0 / std::f64::consts::PI;
        assert!((f / m.f_sum_ev2() - 1.0).abs() < 1e-10);
        assert!((p / m.p_eff() - 1.0).abs() < 1e-10);
        assert!((m.plasma_energy_ev() - (256.0_f64 + 300.0_f64).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(DrudeLorentz::new(vec![]).is_err());
        assert!(DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(10.0, 0.0)]).is_err());
        assert!(DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(f64::NAN, 1.0)]).is_err());
        let m = DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(10.0, 1.0)]).unwrap();
        assert!(m.to_optical_elf("x", 1.0, 1.0, 10).is_err());
        assert!(m.to_optical_elf("x", 1.0, 10.0, 1).is_err());
        let t = m.to_optical_elf("x", 0.1, 1e3, 50).unwrap();
        assert!(t.provenance().contains("synthetic"));
        assert_eq!(t.energy_range_ev(), (0.1, 1e3));
    }
}
