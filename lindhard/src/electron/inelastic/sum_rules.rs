//! Sum rules of an optical energy-loss function.
//!
//! Two checks of the internal consistency of an optical ELF, from Shinotsuka,
//! Da, Tanuma, Yoshikawa, Powell, Penn, Surf. Interface Anal. 49, 238 (2017),
//! doi:10.1002/sia.6123 (open copy: PMC5524379), in Hartree atomic units
//! there and in eV here:
//!
//! - **f-sum** (eq. (18)): the effective number of electrons per target unit
//!   (atom or molecule) taking part in losses up to `W`,
//!   `N_eff(W) = (2 / (π (ħΩ_a)²)) ∫_0^W W' ELF(W') dW'`, with
//!   `(ħΩ_a)² = ħ² n_a e² / (ε₀ m_e)` for the target-unit density `n_a`.
//!   Without a target density, the report gives the effective electron
//!   density `n_eff(W) = (ε₀ m_e / ħ²) (2/π) ∫_0^W W' ELF dW'` instead
//!   (`N_eff = n_eff / n_a`).
//! - **Perfect-screening (Kramers-Kronig) sum** (eq. (19)):
//!   `P_eff = (2/π) ∫ ELF(W)/W dW`. Eq. (19) adds `n(0)⁻²`, the inverse square
//!   of the static refractive index, which is zero for a metal and is not
//!   known from the ELF alone; the report gives the integral only, and the
//!   caller adds `n(0)⁻²` for an insulator. For a complete ELF `P_eff -> 1`
//!   (text below eq. (19)).
//!
//! The report also gives the mean excitation energy `I` from the same ELF,
//! `ln I = ∫ W ELF ln W dW / ∫ W ELF dW` (the definition the Bethe limit uses;
//! see [`super::penn`]).
//!
//! # Integration
//!
//! The integrals run over the tabulated range only (nothing is extrapolated;
//! `OpticalElf` interpolates linearly and the ELF is taken as zero outside the
//! table). They are exact for the piecewise-linear interpolant: the f-sum by
//! Simpson's rule (exact for the quadratic `W ELF(W)` on each segment), `P_eff`
//! in closed form (`∫ (α + βW)/W dW = α ln(W₁/W₀) + β (W₁ - W₀)`), and the
//! `ln W`-weighted integral by an 8-point Gauss-Legendre rule per segment
//! (the integrand is smooth on each segment; relative error below 1e-12).

use std::fmt;

use super::quadrature::GaussLegendre;
use crate::constants::{ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY};
use crate::electron::data::{ElectronDataError, OpticalElf};

/// Electron density (m⁻³) whose free-electron plasma energy squared is
/// `x_ev2` (eV²): `n = (ħΩ)² ε₀ m_e / (ħ² e²)` with `(ħΩ)²` in J².
fn density_from_plasma_energy_sq(x_ev2: f64) -> f64 {
    let x_j2 = x_ev2 * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE;
    x_j2 * VACUUM_PERMITTIVITY * ELECTRON_MASS
        / (HBAR * HBAR * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE)
}

/// The sum rules of one [`OpticalElf`]. See the module docs for the
/// definitions. Built by [`SumRuleReport::new`]; printed by its `Display`
/// implementation.
#[derive(Debug, Clone, PartialEq)]
pub struct SumRuleReport {
    /// Material identity of the ELF.
    pub material: String,
    /// Provenance of the ELF.
    pub provenance: String,
    /// The integration range, the tabulated range of the ELF, eV.
    pub energy_range_ev: (f64, f64),
    /// `∫ W ELF(W) dW` over the table, eV².
    pub f_sum_ev2: f64,
    /// The plasma energy of all the electrons the ELF accounts for,
    /// `ħΩ_p = sqrt((2/π) ∫ W ELF dW)`, eV.
    pub plasma_energy_ev: f64,
    /// The effective electron density `n_eff` of the whole table, m⁻³.
    pub electron_density_per_m3: f64,
    /// `P_eff = (2/π) ∫ ELF(W)/W dW` over the table (no `n(0)⁻²` term).
    pub p_eff: f64,
    /// Mean excitation energy `I`, eV: `ln I = ∫ W ELF ln W dW / ∫ W ELF dW`.
    /// `NaN` if the f-sum is zero.
    pub mean_excitation_energy_ev: f64,
    /// The target-unit (atom or molecule) density used for `N_eff`, m⁻³, if
    /// one was given.
    pub target_density_per_m3: Option<f64>,
    /// The tabulated energies `W_k`, eV.
    pub energy_ev: Vec<f64>,
    /// `n_eff(W_k)`, the cumulative effective electron density up to each
    /// tabulated energy, m⁻³.
    pub cumulative_electron_density_per_m3: Vec<f64>,
    /// The tabulated ELF values, for partial segments.
    elf_values: Vec<f64>,
}

impl SumRuleReport {
    /// The sum rules of `elf`, without a target density (`N_eff` is then not
    /// available, `n_eff` is).
    pub fn new(elf: &OpticalElf) -> Self {
        let w = elf.energy_ev();
        let e = elf.elf_values();
        let gl = GaussLegendre::new(8);
        let mut f_cum = Vec::with_capacity(w.len());
        f_cum.push(0.0);
        let (mut f_sum, mut p_sum, mut ln_sum) = (0.0, 0.0, 0.0);
        for k in 0..w.len() - 1 {
            let (w0, w1, e0, e1) = (w[k], w[k + 1], e[k], e[k + 1]);
            let h = w1 - w0;
            let wm = 0.5 * (w0 + w1);
            let em = 0.5 * (e0 + e1);
            f_sum += h / 6.0 * (w0 * e0 + 4.0 * wm * em + w1 * e1);
            f_cum.push(f_sum);
            let beta = (e1 - e0) / h;
            p_sum += (e0 - beta * w0) * (w1 / w0).ln() + (e1 - e0);
            ln_sum += gl.integrate(&mut |x: f64| [x * (e0 + beta * (x - w0)) * x.ln()], w0, w1)[0];
        }
        let two_over_pi = 2.0 / std::f64::consts::PI;
        let mean_excitation_energy_ev = if f_sum > 0.0 {
            (ln_sum / f_sum).exp()
        } else {
            f64::NAN
        };
        Self {
            material: elf.material().to_string(),
            provenance: elf.provenance().to_string(),
            energy_range_ev: elf.energy_range_ev(),
            f_sum_ev2: f_sum,
            plasma_energy_ev: (two_over_pi * f_sum).sqrt(),
            electron_density_per_m3: density_from_plasma_energy_sq(two_over_pi * f_sum),
            p_eff: two_over_pi * p_sum,
            mean_excitation_energy_ev,
            target_density_per_m3: None,
            energy_ev: w.to_vec(),
            cumulative_electron_density_per_m3: f_cum
                .iter()
                .map(|f| density_from_plasma_energy_sq(two_over_pi * f))
                .collect(),
            elf_values: e.to_vec(),
        }
    }

    /// The sum rules of `elf` with the density of target units (atoms or
    /// molecules, m⁻³), so that `N_eff` per unit is available. The density
    /// must be finite and positive.
    pub fn with_target_density(
        elf: &OpticalElf,
        target_density_per_m3: f64,
    ) -> Result<Self, ElectronDataError> {
        if !(target_density_per_m3.is_finite() && target_density_per_m3 > 0.0) {
            return Err(ElectronDataError::Invalid {
                what: "target density",
                reason: format!("must be finite and positive, got {target_density_per_m3} m^-3"),
            });
        }
        let mut r = Self::new(elf);
        r.target_density_per_m3 = Some(target_density_per_m3);
        Ok(r)
    }

    /// `n_eff(W)`, m⁻³: the effective electron density for losses up to `W`
    /// (eV), exact for the interpolant. Zero below the table and the whole
    /// table's value above it.
    pub fn electron_density_up_to(&self, w_ev: f64) -> f64 {
        let w = &self.energy_ev;
        let c = &self.cumulative_electron_density_per_m3;
        if w_ev.is_nan() || w_ev <= w[0] {
            return 0.0;
        }
        if w_ev >= w[w.len() - 1] {
            return c[c.len() - 1];
        }
        let i = w.partition_point(|&x| x <= w_ev).clamp(1, w.len() - 1);
        let (w0, w1) = (w[i - 1], w[i]);
        let (e0, e1) = (self.elf_values[i - 1], self.elf_values[i]);
        // ∫_{w0}^{x} W ELF(W) dW with the segment's linear ELF, by Simpson's
        // rule (exact for the quadratic integrand).
        let x = w_ev;
        let ex = e0 + (e1 - e0) * (x - w0) / (w1 - w0);
        let (xm, em) = (0.5 * (w0 + x), 0.5 * (e0 + ex));
        let part = (x - w0) / 6.0 * (w0 * e0 + 4.0 * xm * em + x * ex);
        c[i - 1] + density_from_plasma_energy_sq(2.0 / std::f64::consts::PI * part)
    }

    /// `N_eff(W)`: effective electrons per target unit for losses up to `W`
    /// (eV), if a target density was given.
    pub fn effective_electrons_up_to(&self, w_ev: f64) -> Option<f64> {
        self.target_density_per_m3
            .map(|n| self.electron_density_up_to(w_ev) / n)
    }

    /// `N_eff` of the whole table, if a target density was given.
    pub fn effective_electrons(&self) -> Option<f64> {
        self.target_density_per_m3
            .map(|n| self.electron_density_per_m3 / n)
    }
}

impl fmt::Display for SumRuleReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (lo, hi) = self.energy_range_ev;
        writeln!(f, "Sum rules of the optical ELF of {}", self.material)?;
        writeln!(f, "  provenance:               {}", self.provenance)?;
        writeln!(f, "  integration range:        {lo} to {hi} eV")?;
        writeln!(f, "  f-sum, int W ELF dW:      {:.6e} eV^2", self.f_sum_ev2)?;
        writeln!(
            f,
            "  plasma energy (f-sum):    {:.6} eV",
            self.plasma_energy_ev
        )?;
        writeln!(
            f,
            "  effective electron density n_eff: {:.6e} m^-3",
            self.electron_density_per_m3
        )?;
        if let (Some(n), Some(neff)) = (self.target_density_per_m3, self.effective_electrons()) {
            writeln!(f, "  target-unit density:      {n:.6e} m^-3")?;
            writeln!(f, "  N_eff per target unit:    {neff:.6}")?;
        }
        writeln!(
            f,
            "  P_eff = (2/pi) int ELF/W dW: {:.6}  (add n(0)^-2 for an insulator)",
            self.p_eff
        )?;
        writeln!(
            f,
            "  mean excitation energy I: {:.6} eV",
            self.mean_excitation_energy_ev
        )
    }
}
