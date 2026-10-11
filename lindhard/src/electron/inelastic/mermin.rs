//! The Mermin dielectric function of a free-electron gas, and the Mermin-ELF
//! inelastic model built from a fit of oscillators to an optical ELF.
//!
//! # Sources and what was and was not read
//!
//! Mermin's relaxation-time dielectric function is N. D. Mermin, "Lindhard
//! dielectric function in the relaxation-time approximation", Phys. Rev. B 1,
//! 2362 (1970), doi:10.1103/PhysRevB.1.2362. That paper is closed access and
//! **was not opened**, so its own equation numbers are not known here. The
//! formula is taken from the account in P. de Vera, S. Taioli, P. E.
//! Trevisanutto, M. Dapor, I. Abril, S. Simonucci and R. Garcia-Molina,
//! Int. J. Mol. Sci. 23, 6121 (2022), doi:10.3390/ijms23116121 (open access,
//! PMC9181504, "dV2022", read 2026-10-07), section 2.1.1, whose eq. (4)
//! cites Mermin (1970) as its reference 33:
//!
//! ```text
//!                     (1 + iγ/ω) [ε_L(q, ω + iγ) - 1]
//! ε_M(q, ω) = 1 + ---------------------------------------------   (dV2022 eq. (4))
//!                  1 + (iγ/ω) [ε_L(q, ω + iγ) - 1] / [ε_L(q, 0) - 1]
//! ```
//!
//! (energies `ħω`, `ħγ`; `ε_L` the Lindhard function,
//! [`super::lindhard_gas`]). The Mermin-ELF (MELF) of the MELF-GOS method is
//! (dV2022 eq. (3), without the threshold step `Θ(E - E_th)`)
//!
//! ```text
//! Im[-1/ε(q, ω)] = Σ_j (A_j/(ħω_j)²) Im[-1/ε_M(q, ω; ω_j, γ_j)],
//! ```
//!
//! where `ε_M(q, ω; ω_j, γ_j)` uses the Lindhard function of a gas of plasma
//! frequency `ω_j`, and `A_j`, `ω_j`, `γ_j` are fitted to the optical ELF by
//! dV2022 eq. (5) (see [`super::mermin_fit`]); at `q = 0` the Mermin-type
//! ELF is the Drude-Lorentz one ([`super::drude`]).
//!
//! Checks made here instead of a comparison with Mermin's own paper (all in
//! the unit tests):
//!
//! * `γ -> 0` gives the Lindhard function to 1e-8 (the analytic
//!   continuation of the Lindhard function to complex frequency is checked
//!   against the real-frequency form of [`super::lindhard_gas`]);
//! * `q -> 0` gives the Drude function `1 - ω_p²/(ω² + iγω)`;
//! * `ω -> 0` gives the static Lindhard function `ε_L(q, 0)` (the
//!   local-equilibrium property that motivates the formula);
//! * the f-sum rule `∫ ω Im[-1/ε_M] dω = (π/2) ω_p²` at finite `q`. We did
//!   not find it stated in the sources read, so it is a numerical check of the
//!   transcription, not an input.
//!
//! The independent algebraic form of Latyshev and Yushkanov, arXiv:1212.6260
//! (eq. (5.5″), read as described in [`super::lindhard_gas`]), is the same
//! function with `x + iy` for `ω + iγ`.
//!
//! # Complex-frequency Lindhard function
//!
//! `ε_L(q, ω) = 1 + (k_TF²/q²) [1/2 + (A(z - u) + A(z + u))/(8z)]` with
//! `A(s) = (1 - s²) Log((s + 1)/(s - 1))` (principal logarithm of the ratio,
//! analytic off `[-1, 1]` and vanishing at infinity) and `u = (ω + iγ)/(q
//! k_F)`. For `γ -> 0+` this gives the real-frequency Lindhard function with
//! the electron-hole continuum as its imaginary part. When `|z ± u| > 4` the
//! bracket is evaluated from the same series as the real-frequency code
//! (every term of which is `O(1/|s|^m)`), avoiding cancellation.
//!
//! # The model
//!
//! [`MerminPenn`] uses the MELF loss function in the DIIMFP of the same
//! nonrelativistic kinematics as the Penn models ([`super::penn`],
//! S2017 eqs. (2)-(3)):
//!
//! ```text
//! p(T, ω) = (1/(π T')) ∫_{q-}^{q+} (dq/q) Im[-1/ε(q, ω)],
//! ```
//!
//! and integrates it over `ω` (in `ln ω`, from `1e-8 T` to `T`) for the
//! inverse IMFP and the stopping power. This is our own evaluation of a
//! published construction; the omission of the threshold step is a choice
//! documented in `docs/data-provenance.md` (a fitted Drude tail extends below
//! any real band gap).
//!
//! # Comparison with the single-pole model (synthetic Drude fixture)
//!
//! For the synthetic Drude plasmon of the tests (`E_p = 20 eV`, `γ = 5 eV`,
//! 240 tabulated energies; the fit returns it to 1e-13), ratios Mermin /
//! single-pole of the stopping power `S` and of the inverse IMFP `1/λ`
//! (default tolerances, release build; outputs of this code, not reference
//! values):
//!
//! ```text
//!    E / eV     S ratio   1/λ ratio
//!      30.0      1.9222     2.0818
//!     100.0      0.9832     1.0235
//!    1000.0      1.0002     1.0173
//!   10000.0      0.9991     1.0121
//!   50000.0      0.9988     1.0100
//! ```
//!
//! The stopping powers agree to 2e-3 from 1 keV (both obey the f-sum);
//! below about 100 eV the models differ strongly, as the full Penn model does
//! from the single pole (see [`super::full_penn`]). Nothing was tuned.

use super::drude::DrudeLorentzOscillator;
use super::lindhard_gas::LindhardGas;
use super::mermin_fit::{fit_mermin_oscillators, MerminFit, MerminFitOptions};
use super::penn::{hartree_ev, InelasticPoint};
use super::quadrature::{integrate_segments, GaussLegendre};
use crate::constants::BOHR_RADIUS;
use crate::electron::data::{ElectronDataError, OpticalElf};
use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Sub};

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Default relative tolerance of the integrals of the Mermin model.
pub const DEFAULT_MERMIN_TOLERANCE: f64 = 1.0e-5;

/// The lower limit of the loss integrals of the inverse IMFP and the
/// stopping power, as a fraction of the kinetic energy `T`: they integrate
/// `ω` over `[LOWEST_LOSS_FRACTION · T, T]` (module docs, "The model"). A
/// numerical lower limit of the integration in `ln ω`, not a physical
/// threshold (the fitted Drude-Lorentz ELF has none, module docs). The
/// inelastic table opens the sampled loss window of a Mermin
/// row at the same limit ([`super::table`], "Rows below the ELF table"), so
/// the stored rate and the sampled losses cover one window (#342).
pub(crate) const LOWEST_LOSS_FRACTION: f64 = 1.0e-8;

const GL_ORDER: usize = 5;

// ---- a minimal complex number (no dependency) ----

#[derive(Debug, Clone, Copy, PartialEq)]
struct C {
    re: f64,
    im: f64,
}

impl C {
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }
    fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    fn ln(self) -> Self {
        Self::new(self.abs().ln(), self.im.atan2(self.re))
    }
    fn scale(self, k: f64) -> Self {
        Self::new(self.re * k, self.im * k)
    }
}

impl Add for C {
    type Output = C;
    fn add(self, o: C) -> C {
        C::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for C {
    type Output = C;
    fn sub(self, o: C) -> C {
        C::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for C {
    type Output = C;
    fn mul(self, o: C) -> C {
        C::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}
impl Div for C {
    type Output = C;
    fn div(self, o: C) -> C {
        let d = o.norm_sq();
        C::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

/// `A(s) = (1 - s²) Log((s + 1)/(s - 1))`, zero at `s = ±1`.
fn a_complex(s: C) -> C {
    let one = C::real(1.0);
    let d = one - s * s;
    if d.re == 0.0 && d.im == 0.0 {
        return C::real(0.0);
    }
    d * ((s + one) / (s - one)).ln()
}

/// The bracket `1/2 + (A(z - u) + A(z + u))/(8z)` for complex `u`.
fn bracket(z: f64, u: C) -> C {
    let hi = u + C::real(z);
    let lo = u - C::real(z);
    if hi.abs() > 4.0 && lo.abs() > 4.0 {
        return bracket_series(hi, lo);
    }
    let half = C::real(0.5);
    half + (a_complex(C::real(z) - u) + a_complex(C::real(z) + u)).scale(1.0 / (8.0 * z))
}

/// `-(1/4) Σ_k c_k S_m/(u² - z²)^m`, `m = 2k + 1`, `c_k = 4/((2k+1)(2k+3))`,
/// `S_m = Σ_j hi^(m-1-j) lo^j` (`hi = u + z`, `lo = u - z`), from the
/// expansion of `A` for `|s| > 1` (see `lindhard_gas::shape`).
fn bracket_series(hi: C, lo: C) -> C {
    let d = hi * lo;
    let mut sum = C::real(0.0);
    let mut s_m = C::real(1.0); // S_1
    let mut lo_pow = lo; // lo^m
    let mut d_pow = d; // d^m
    for m in 1..=400usize {
        if m % 2 == 1 {
            let k = (m / 2) as f64;
            let ck = 4.0 / ((2.0 * k + 1.0) * (2.0 * k + 3.0));
            let term = (s_m / d_pow).scale(ck);
            sum = sum + term;
            if term.abs() <= 1e-17 * sum.abs() {
                break;
            }
        }
        // S_{m+1} = hi S_m + lo^m
        s_m = hi * s_m + lo_pow;
        lo_pow = lo_pow * lo;
        d_pow = d_pow * d;
    }
    sum.scale(-0.25)
}

/// The Lindhard dielectric function at complex frequency, atomic units.
fn lindhard_complex(gas: &LindhardGas, q: f64, w: C) -> C {
    let kf = gas.fermi_wavenumber();
    let z = q / (2.0 * kf);
    let u = w.scale(1.0 / (q * kf));
    C::real(1.0) + bracket(z, u).scale(gas.thomas_fermi_wavenumber_sq() / (q * q))
}

/// The Mermin dielectric function of a free-electron gas with a relaxation
/// energy `ħγ` (see the module docs). Atomic units (Hartree, bohr).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MerminGas {
    gas: LindhardGas,
    gamma: f64,
}

impl MerminGas {
    /// The gas `gas` with relaxation energy `gamma >= 0` (Hartree).
    /// `gamma = 0` is the Lindhard function.
    pub fn new(gas: LindhardGas, gamma: f64) -> Result<Self> {
        if !(gamma.is_finite() && gamma >= 0.0) {
            return Err(ElectronDataError::Invalid {
                what: "Mermin relaxation energy",
                reason: format!("must be finite and non-negative, got {gamma} Hartree"),
            });
        }
        Ok(Self { gas, gamma })
    }

    /// The underlying Lindhard gas.
    pub fn gas(&self) -> &LindhardGas {
        &self.gas
    }

    /// The relaxation energy `ħγ`, Hartree.
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// `ε_M(q, ω)` as `(Re, Im)` for `q > 0`, `ω > 0` (dV2022 eq. (4)).
    pub fn epsilon(&self, q: f64, w: f64) -> (f64, f64) {
        if self.gamma == 0.0 {
            return self.gas.epsilon(q, w);
        }
        let l_w = lindhard_complex(&self.gas, q, C::new(w, self.gamma)) - C::real(1.0);
        let z = q / (2.0 * self.gas.fermi_wavenumber());
        let l_0 = bracket(z, C::real(0.0)).re * self.gas.thomas_fermi_wavenumber_sq() / (q * q);
        let ig = C::new(0.0, self.gamma / w);
        let num = (C::real(1.0) + ig) * l_w;
        let den = C::real(1.0) + ig * l_w.scale(1.0 / l_0);
        let eps = C::real(1.0) + num / den;
        (eps.re, eps.im)
    }

    /// The loss function `Im[-1/ε_M(q, ω)] = Im ε/|ε|²`.
    pub fn loss(&self, q: f64, w: f64) -> f64 {
        let (re, im) = self.epsilon(q, w);
        im / (re * re + im * im)
    }
}

#[derive(Clone)]
struct Term {
    /// `A_j/(ħω_j)²`, dimensionless.
    weight: f64,
    gas: MerminGas,
}

/// The Mermin-ELF inelastic model: a fit of Mermin oscillators to one
/// optical ELF, extended to finite momentum by the Mermin function, and the
/// DIIMFP, IMFP and stopping power built from it. See the module docs.
#[derive(Clone)]
pub struct MerminPenn {
    elf: OpticalElf,
    fit: MerminFit,
    terms: Vec<Term>,
    fermi: f64,
    rel_tol: f64,
    gl: GaussLegendre,
}

impl std::fmt::Debug for MerminPenn {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MerminPenn")
            .field("material", &self.elf.material())
            .field("oscillators", &self.fit.oscillators)
            .field("fermi_hartree", &self.fermi)
            .field("rel_tol", &self.rel_tol)
            .finish()
    }
}

fn check_energy(what: &'static str, v: f64) -> Result<()> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(ElectronDataError::Invalid {
            what,
            reason: format!("must be finite and positive, got {v} eV"),
        })
    }
}

impl MerminPenn {
    /// The model of `fit`, a fit to `elf` (Fermi energy zero, default
    /// tolerance). Oscillators of zero amplitude are carried in the fit
    /// report but do not contribute.
    pub fn new(elf: OpticalElf, fit: MerminFit) -> Result<Self> {
        let h = hartree_ev();
        let mut terms = Vec::new();
        for o in &fit.oscillators {
            check_oscillator(o)?;
            if o.strength_ev2 == 0.0 {
                continue;
            }
            let wp = o.energy_ev / h;
            terms.push(Term {
                weight: o.strength_ev2 / (o.energy_ev * o.energy_ev),
                gas: MerminGas::new(LindhardGas::from_plasma_frequency(wp), o.width_ev / h)?,
            });
        }
        if terms.is_empty() {
            return Err(ElectronDataError::Invalid {
                what: "Mermin fit",
                reason: "every oscillator has zero amplitude".into(),
            });
        }
        Ok(Self {
            elf,
            fit,
            terms,
            fermi: 0.0,
            rel_tol: DEFAULT_MERMIN_TOLERANCE,
            gl: GaussLegendre::new(GL_ORDER),
        })
    }

    /// Fit `elf` with `options` ([`fit_mermin_oscillators`]) and build the
    /// model from the result.
    pub fn fit(elf: OpticalElf, options: &MerminFitOptions) -> Result<Self> {
        let fit = fit_mermin_oscillators(&elf, options)?;
        Self::new(elf, fit)
    }

    /// Set the Fermi energy (eV, finite and non-negative), as for the Penn
    /// models: energies are kinetic energies above the Fermi level and the
    /// kinematics use `T' = T + E_F`.
    pub fn with_fermi_energy_ev(mut self, fermi_ev: f64) -> Result<Self> {
        if !(fermi_ev.is_finite() && fermi_ev >= 0.0) {
            return Err(ElectronDataError::Invalid {
                what: "Fermi energy",
                reason: format!("must be finite and non-negative, got {fermi_ev} eV"),
            });
        }
        self.fermi = fermi_ev / hartree_ev();
        Ok(self)
    }

    /// Set the relative integration tolerance, in `[1e-9, 1e-2]`.
    pub fn with_relative_tolerance(mut self, rel_tol: f64) -> Result<Self> {
        if !(1e-9..=1e-2).contains(&rel_tol) {
            return Err(ElectronDataError::Invalid {
                what: "integration tolerance",
                reason: format!("must be in [1e-9, 1e-2], got {rel_tol}"),
            });
        }
        self.rel_tol = rel_tol;
        Ok(self)
    }

    /// The optical ELF the oscillators were fitted to.
    pub fn optical_elf(&self) -> &OpticalElf {
        &self.elf
    }

    /// The fit: parameters, residuals, f-sum and `P_eff`.
    pub fn fit_report(&self) -> &MerminFit {
        &self.fit
    }

    /// The Fermi energy, eV.
    pub fn fermi_energy_ev(&self) -> f64 {
        self.fermi * hartree_ev()
    }

    /// The relative integration tolerance.
    pub fn relative_tolerance(&self) -> f64 {
        self.rel_tol
    }

    fn loss_au(&self, q: f64, w: f64) -> f64 {
        self.terms.iter().map(|t| t.weight * t.gas.loss(q, w)).sum()
    }

    /// `Im[-1/ε(q, ω)]` of the Mermin-ELF model at momentum transfer `q`
    /// (m⁻¹) and energy loss `ω` (eV). At `q = 0` it is the fitted
    /// Drude-Lorentz ELF.
    pub fn loss_function(&self, q_per_m: f64, loss_ev: f64) -> f64 {
        if !(q_per_m.is_finite() && q_per_m >= 0.0 && loss_ev.is_finite() && loss_ev > 0.0) {
            return 0.0;
        }
        if q_per_m == 0.0 {
            return self.fit.model_elf(loss_ev);
        }
        self.loss_au(q_per_m * BOHR_RADIUS, loss_ev / hartree_ev())
    }

    /// The DIIMFP `p(T, ω)` (S2017 eq. (2), nonrelativistic) in m⁻¹ eV⁻¹ at
    /// kinetic energy `energy_ev` above the Fermi level and loss `loss_ev`.
    /// Zero for `ω <= 0` or `ω > T`.
    pub fn diimfp_per_m_ev(&self, energy_ev: f64, loss_ev: f64) -> Result<f64> {
        check_energy("electron energy", energy_ev)?;
        if !loss_ev.is_finite() {
            return Err(ElectronDataError::Invalid {
                what: "energy loss",
                reason: format!("must be finite, got {loss_ev} eV"),
            });
        }
        let h = hartree_ev();
        Ok(self.diimfp_au(energy_ev / h, loss_ev / h) / (BOHR_RADIUS * h))
    }

    fn diimfp_au(&self, t: f64, w: f64) -> f64 {
        if !(w > 0.0 && w <= t) {
            return 0.0;
        }
        let tp = t + self.fermi;
        let k = (2.0 * tp).sqrt();
        let s = (2.0 * (tp - w)).max(0.0).sqrt();
        let (q_minus, q_plus) = (2.0 * w / (k + s), k + s);
        if q_plus.is_nan() || q_plus <= q_minus {
            return 0.0;
        }
        let (lo, hi) = (q_minus.ln(), q_plus.ln());
        let mut breaks: Vec<f64> = (0..=8).map(|i| lo + (hi - lo) * i as f64 / 8.0).collect();
        // the edges of each gas's electron-hole continuum, q = ω/k_F ± ...
        for t in &self.terms {
            let kf = t.gas.gas().fermi_wavenumber();
            for q in [2.0 * kf, w / kf] {
                if q > q_minus && q < q_plus {
                    breaks.push(q.ln());
                }
            }
        }
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        let integral = integrate_segments(
            &self.gl,
            &mut |u: f64| [self.loss_au(u.exp(), w)],
            &breaks,
            self.rel_tol,
        )[0];
        integral / (PI * tp)
    }

    /// `(λ⁻¹, S)` at kinetic energy `t` (Hartree), atomic units.
    fn imfp_and_stopping_au(&self, t: f64) -> (f64, f64) {
        if t.is_nan() || t <= 0.0 {
            return (0.0, 0.0);
        }
        let (lo, hi) = ((LOWEST_LOSS_FRACTION * t).ln(), t.ln());
        let n = 24;
        let mut breaks: Vec<f64> = (0..=n)
            .map(|i| lo + (hi - lo) * i as f64 / n as f64)
            .collect();
        for term in &self.terms {
            let e = term.gas.gas().plasma_frequency();
            if e > lo.exp() && e < t {
                breaks.push(e.ln());
            }
        }
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        let r = integrate_segments(
            &self.gl,
            &mut |u: f64| {
                let w = u.exp();
                let p = self.diimfp_au(t, w);
                // d ω = ω d ln ω
                [p * w, p * w * w]
            },
            &breaks,
            self.rel_tol,
        );
        (r[0], r[1])
    }

    /// The inverse IMFP (m⁻¹) and the stopping power (eV/m) at kinetic
    /// energy `energy_ev` above the Fermi level.
    pub fn imfp_and_stopping(&self, energy_ev: f64) -> Result<InelasticPoint> {
        check_energy("electron energy", energy_ev)?;
        let h = hartree_ev();
        let (inv, s) = self.imfp_and_stopping_au(energy_ev / h);
        Ok(InelasticPoint {
            energy_ev,
            inverse_imfp_per_m: inv / BOHR_RADIUS,
            stopping_ev_per_m: s * h / BOHR_RADIUS,
        })
    }

    /// The inelastic mean free path `λ(E)`, m.
    pub fn imfp_m(&self, energy_ev: f64) -> Result<f64> {
        Ok(self.imfp_and_stopping(energy_ev)?.imfp_m())
    }

    /// The stopping power `S(E)`, eV/m.
    pub fn stopping_power_ev_per_m(&self, energy_ev: f64) -> Result<f64> {
        Ok(self.imfp_and_stopping(energy_ev)?.stopping_ev_per_m)
    }

    /// IMFP and stopping power on a grid of energies, in order. Sequential
    /// and deterministic.
    pub fn tabulate(&self, energies_ev: &[f64]) -> Result<Vec<InelasticPoint>> {
        energies_ev
            .iter()
            .map(|&e| self.imfp_and_stopping(e))
            .collect()
    }
}

/// Amplitude finite and non-negative, energy and width finite and positive.
pub(crate) fn check_oscillator(o: &DrudeLorentzOscillator) -> Result<()> {
    let bad = |name: &str, v: f64| ElectronDataError::Invalid {
        what: "Mermin oscillator",
        reason: format!("{name} is invalid: {v}"),
    };
    if !(o.strength_ev2.is_finite() && o.strength_ev2 >= 0.0) {
        return Err(bad("amplitude (must be finite, >= 0)", o.strength_ev2));
    }
    if !(o.energy_ev.is_finite() && o.energy_ev > 0.0) {
        return Err(bad("energy (must be finite, > 0)", o.energy_ev));
    }
    if !(o.width_ev.is_finite() && o.width_ev > 0.0) {
        return Err(bad("width (must be finite, > 0)", o.width_ev));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gas(kf: f64) -> LindhardGas {
        LindhardGas::from_fermi_wavenumber(kf)
    }

    #[test]
    fn complex_bracket_series_matches_direct_where_both_hold() {
        for (z, u) in [
            (0.3, C::new(5.0, 0.4)),
            (1.5, C::new(7.0, 2.0)),
            (0.1, C::new(4.5, 0.1)),
            (2.0, C::new(0.5, 5.0)),
            (3.0, C::new(0.0, 5.0)),
        ] {
            let (hi, lo) = (u + C::real(z), u - C::real(z));
            assert!(hi.abs() > 4.0 && lo.abs() > 4.0, "{z} {u:?}");
            let direct = C::real(0.5)
                + (a_complex(C::real(z) - u) + a_complex(C::real(z) + u)).scale(1.0 / (8.0 * z));
            let series = bracket_series(hi, lo);
            let diff = (direct - series).abs();
            assert!(diff <= 1e-9 * series.abs().max(1e-3), "{z} {u:?} {diff:e}");
        }
    }

    #[test]
    fn complex_lindhard_reduces_to_the_real_frequency_form() {
        // tiny positive imaginary frequency -> the real-frequency Lindhard
        // function, including the electron-hole imaginary part
        for kf in [0.4, 0.9] {
            let g = gas(kf);
            for (q, w) in [(0.3, 0.2), (0.8, 0.5), (1.5, 1.4), (0.2, 1.0), (2.5, 0.3)] {
                let c = lindhard_complex(&g, q, C::new(w, 1e-13));
                let (re, im) = g.epsilon(q, w);
                assert!(
                    (c.re - re).abs() <= 1e-9 * re.abs().max(1.0),
                    "{kf} {q} {w}: {c:?} vs {re} {im}"
                );
                assert!(
                    (c.im - im).abs() <= 1e-9 * im.abs().max(1.0),
                    "{kf} {q} {w}: {c:?} vs {re} {im}"
                );
            }
        }
    }

    #[test]
    fn mermin_reduces_to_lindhard_as_gamma_goes_to_zero_to_1e8() {
        for kf in [0.4, 0.9] {
            let g = gas(kf);
            for (q, w) in [(0.3, 0.2), (0.8, 0.5), (1.5, 1.4), (0.2, 1.0), (2.5, 0.3)] {
                let (lre, lim) = g.epsilon(q, w);
                // exactly zero: the identity
                let m0 = MerminGas::new(g, 0.0).unwrap().epsilon(q, w);
                assert_eq!(m0, (lre, lim));
                for gamma in [1e-11, 1e-13] {
                    let (re, im) = MerminGas::new(g, gamma).unwrap().epsilon(q, w);
                    assert!(
                        (re - lre).abs() <= 1e-8 * lre.abs().max(1.0),
                        "{kf} {q} {w}"
                    );
                    assert!(
                        (im - lim).abs() <= 1e-8 * lim.abs().max(1.0),
                        "{kf} {q} {w}"
                    );
                }
            }
        }
    }

    #[test]
    fn mermin_is_drude_at_small_q() {
        // ε_M(q -> 0, ω) = 1 - ω_p²/(ω² + iγω)
        let g = LindhardGas::from_plasma_frequency(0.5);
        let gamma = 0.12;
        let m = MerminGas::new(g, gamma).unwrap();
        let wp2 = 0.25;
        for w in [0.1, 0.5, 1.3] {
            let (re, im) = m.epsilon(1e-3, w);
            let d = C::new(w * w, gamma * w);
            let drude = C::real(1.0) - C::real(wp2) / d;
            assert!((re - drude.re).abs() < 1e-3, "{w} {re} {}", drude.re);
            assert!((im - drude.im).abs() < 1e-3, "{w} {im} {}", drude.im);
        }
    }

    #[test]
    fn mermin_static_limit_is_the_static_lindhard_function() {
        let g = gas(0.7);
        let m = MerminGas::new(g, 0.05).unwrap();
        for q in [0.3, 0.9, 2.0] {
            let (re, im) = m.epsilon(q, 1e-9);
            let (lre, _) = g.epsilon(q, 0.0);
            assert!((re - lre).abs() <= 1e-6 * lre, "{q} {re} {lre}");
            assert!(im.abs() <= 1e-6 * lre, "{q} {im}");
        }
    }

    #[test]
    fn mermin_obeys_the_f_sum_rule_at_finite_q() {
        // ∫ ω Im[-1/ε_M] dω = (π/2) ω_p²; the tail above the cutoff is
        // ω_p² γ/ω_c (ε - 1 ≈ -ω_p²/(ω(ω + iγ)) there)
        let g = LindhardGas::from_plasma_frequency(0.6);
        let wp2 = g.plasma_frequency().powi(2);
        let gl = GaussLegendre::new(8);
        for (gamma, q) in [(0.08, 0.3), (0.15, 0.8), (0.05, 1.2)] {
            let m = MerminGas::new(g, gamma).unwrap();
            let w_c = 4000.0;
            let breaks: Vec<f64> = (0..=60)
                .map(|i| (1e-6f64).ln() + ((w_c / 1e-6f64).ln()) * i as f64 / 60.0)
                .collect();
            let r = integrate_segments(
                &gl,
                &mut |u: f64| {
                    let w = u.exp();
                    [w * w * m.loss(q, w)]
                },
                &breaks,
                1e-9,
            )[0];
            let total = r + wp2 * gamma / w_c;
            let want = 0.5 * PI * wp2;
            assert!(
                ((total - want) / want).abs() < 2e-4,
                "gamma {gamma} q {q}: {total} vs {want}"
            );
        }
    }

    #[test]
    fn invalid_relaxation_energy_is_rejected() {
        let g = gas(0.5);
        assert!(MerminGas::new(g, -1.0).is_err());
        assert!(MerminGas::new(g, f64::NAN).is_err());
    }
}
