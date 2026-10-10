//! Fitting Mermin oscillators to an optical energy-loss function (the
//! MELF-GOS fit).
//!
//! # What is fitted
//!
//! In the MELF-GOS method the optical ELF is represented by a weighted sum of
//! Mermin-type energy-loss functions, whose parameters `A_j`, `ħω_j`, `ħγ_j`
//! (intensity, position, width) are "determined by fitting the outer-shell
//! ELF contribution to the available experimental optical spectrum" through
//! the `k = 0` form, in which the Mermin ELF is the Drude-Lorentz one:
//!
//! ```text
//! ELF(E) = Σ_j A_j E ħγ_j / ([(ħω_j)² - E²]² + [E ħγ_j]²)
//! ```
//!
//! (P. de Vera et al., Int. J. Mol. Sci. 23, 6121 (2022),
//! doi:10.3390/ijms23116121, PMC9181504, section 2.1.1, eqs. (3) and (5),
//! read 2026-10-07; see [`super::mermin`]). The same paper says the
//! consistency of the fit "is checked by fulfilling the Kramers-Kronig and
//! f-sum rules"; this module reports the f-sum and `P_eff` of the fit next to
//! those of the data ([`MerminFit`]). It does not say how the least-squares
//! problem is solved, and **the algorithm below is our own**, built from
//! textbook pieces that were not opened for this implementation:
//!
//! * the amplitudes enter linearly, so for fixed `(ω_j, γ_j)` they are the
//!   solution of a **non-negative least-squares** problem (Lawson-Hanson
//!   active-set algorithm, C. L. Lawson and R. J. Hanson, *Solving Least
//!   Squares Problems*, 1974, with each subproblem solved by Householder QR);
//!   amplitudes are therefore non-negative by construction;
//! * the energies and widths are found by **Levenberg-Marquardt** iterations on
//!   `(ln ω_j, ln γ_j)` of the amplitude-eliminated residual (variable
//!   projection, G. H. Golub and V. Pereyra, 1973), with a central-difference
//!   Jacobian.
//!
//! Every step is a fixed sequence of floating-point operations on the input
//! table: the fit is deterministic and does not depend on the thread count.
//!
//! # Limits
//!
//! The fit is a local optimisation. The amplitudes of the start are ignored
//! (they are solved for); the energies and widths of the start decide which
//! local minimum is reached. The default start places the oscillators at the
//! tallest local maxima of the tabulated ELF (completed by log-spaced
//! energies), which is a heuristic; for real data the caller should inspect
//! [`MerminFit::residuals`] and, if needed, supply a start. The threshold
//! step `Θ(E - E_th)` of dV2022 eq. (3) is not fitted.
//!
//! The whole tabulated ELF is fitted, inner-shell edges included, where
//! dV2022 fits Mermin oscillators to the outer-shell ELF only (their eqs. (1)
//! and (3); inner shells by GOS, eq. (2)) and checks the fit against the
//! Kramers-Kronig and f-sum rules (after eq. (5)). The default of three
//! oscillators is coarse for real data: for the committed Cu ELF it carries
//! 57 % of the table's f-sum and 67 % of its `P_eff`, and the Mermin IMFP is
//! 18 to 34 % longer than the single-pole one on the table (#300;
//! `docs/validation.md`, "Cu: the Mermin IMFP and the default oscillator
//! fit"; the default is #306). The widths are bounded below by the local knot
//! spacing (see [`fit_mermin_oscillators`]); without that bound, fits of the
//! Cu ELF with 10 or 16 oscillators placed an oscillator far narrower than
//! the knot spacing between two knots, which the residuals do not see and the
//! sum rules do (#307). The floor does not stop the same fits from placing
//! a very wide oscillator whose f-sum lies above the table (#311). Read
//! [`MerminFit::f_sum_ev2`] against [`MerminFit::data_f_sum_ev2`] in any
//! case.
//!
//! # Closed-form sum rules of the fit
//!
//! For the Drude-Lorentz sum the f-sum is `(π/2) Σ A_j` and
//! `P_eff = Σ A_j/(ħω_j)²` ([`super::drude`]; the derivation there is ours, by
//! residues, with a numerical cross-check), against the numerical sums of the
//! data table over its tabulated range ([`super::sum_rules`]).

// Index loops are the clearest form of the small dense linear algebra below.
#![allow(clippy::needless_range_loop)]

use super::drude::DrudeLorentzOscillator;
use super::mermin::check_oscillator;
use super::sum_rules::SumRuleReport;
use crate::electron::data::{ElectronDataError, OpticalElf};
use std::f64::consts::FRAC_PI_2;
use std::fmt;

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Largest number of oscillators accepted.
pub const MAX_OSCILLATORS: usize = 16;

/// How the residuals of the tabulated points are weighted in the fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FitWeighting {
    /// Every tabulated point counts equally (absolute residuals).
    Uniform,
    /// Residuals are divided by `ELF_i + floor · max ELF`, so that the
    /// tails of the ELF (decades below the peak) are not ignored. `floor`
    /// must be in `(0, 1]`.
    Relative {
        /// Fraction of the maximum ELF added to the divisor.
        floor: f64,
    },
}

/// Settings of [`fit_mermin_oscillators`].
#[derive(Debug, Clone, PartialEq)]
pub struct MerminFitOptions {
    /// Number of oscillators, `1..=MAX_OSCILLATORS`. Ignored if `start` is
    /// given (its length is used).
    pub n_oscillators: usize,
    /// Starting energies and widths (amplitudes are ignored), or `None` for
    /// the default start described in the module docs.
    pub start: Option<Vec<DrudeLorentzOscillator>>,
    /// The weighting of the residuals.
    pub weighting: FitWeighting,
    /// Maximum number of Levenberg-Marquardt iterations.
    pub max_iterations: usize,
}

impl Default for MerminFitOptions {
    fn default() -> Self {
        Self {
            n_oscillators: 3,
            start: None,
            weighting: FitWeighting::Relative { floor: 1e-2 },
            max_iterations: 500,
        }
    }
}

/// The result of a fit: parameters, residuals and the sum rules. See the
/// module docs.
#[derive(Debug, Clone, PartialEq)]
pub struct MerminFit {
    /// The fitted oscillators (`strength_ev2` is the amplitude `A_j`, eV²,
    /// `>= 0`; `energy_ev` is `ħω_j`; `width_ev` is `ħγ_j`), sorted by energy.
    pub oscillators: Vec<DrudeLorentzOscillator>,
    /// The tabulated energies, eV.
    pub energy_ev: Vec<f64>,
    /// The tabulated ELF values.
    pub data: Vec<f64>,
    /// The fitted ELF at the tabulated energies.
    pub model: Vec<f64>,
    /// `model - data` at the tabulated energies.
    pub residuals: Vec<f64>,
    /// The weighting that was minimised.
    pub weighting: FitWeighting,
    /// Root-mean-square of the weighted residuals (the minimised quantity).
    pub weighted_rms: f64,
    /// Root-mean-square of `model - data`.
    pub rms_residual: f64,
    /// Largest `|model - data|`.
    pub max_abs_residual: f64,
    /// Levenberg-Marquardt iterations used.
    pub iterations: usize,
    /// Whether a convergence criterion was met before `max_iterations`.
    pub converged: bool,
    /// `(π/2) Σ A_j`, eV²: the f-sum `∫ E ELF dE` of the fit over `(0, ∞)`.
    pub f_sum_ev2: f64,
    /// `Σ A_j/(ħω_j)²`: `P_eff` of the fit over `(0, ∞)` (no `n(0)⁻²` term).
    pub p_eff: f64,
    /// `∫ E ELF dE` of the data over its tabulated range, eV².
    pub data_f_sum_ev2: f64,
    /// `P_eff` of the data over its tabulated range.
    pub data_p_eff: f64,
}

fn drude_basis(w: f64, e: f64, g: f64) -> f64 {
    let d = w * w - e * e;
    g * w / (d * d + g * g * w * w)
}

impl MerminFit {
    /// The fitted ELF at `energy_ev` (the Drude-Lorentz sum of the
    /// oscillators; zero for `energy_ev <= 0`).
    pub fn model_elf(&self, energy_ev: f64) -> f64 {
        if energy_ev <= 0.0 {
            return 0.0;
        }
        self.oscillators
            .iter()
            .map(|o| o.strength_ev2 * drude_basis(energy_ev, o.energy_ev, o.width_ev))
            .sum()
    }

    /// The plasma energy of the fit, `sqrt(Σ A_j)`, eV.
    pub fn plasma_energy_ev(&self) -> f64 {
        self.oscillators
            .iter()
            .map(|o| o.strength_ev2)
            .sum::<f64>()
            .sqrt()
    }

    /// One line to record in metadata: the number of oscillators and their
    /// parameters.
    pub fn summary(&self) -> String {
        let osc: Vec<String> = self
            .oscillators
            .iter()
            .map(|o| {
                format!(
                    "(A={:.6e} eV^2, E={:.6e} eV, gamma={:.6e} eV)",
                    o.strength_ev2, o.energy_ev, o.width_ev
                )
            })
            .collect();
        format!("{} oscillators {}", self.oscillators.len(), osc.join(" "))
    }
}

impl fmt::Display for MerminFit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Mermin-ELF fit: {} oscillators, {} iterations, converged: {}",
            self.oscillators.len(),
            self.iterations,
            self.converged
        )?;
        for (j, o) in self.oscillators.iter().enumerate() {
            writeln!(
                f,
                "  {j}: A = {:.6e} eV^2, E = {:.6e} eV, gamma = {:.6e} eV",
                o.strength_ev2, o.energy_ev, o.width_ev
            )?;
        }
        writeln!(
            f,
            "  residual rms {:.3e} (weighted {:.3e}), max |r| {:.3e}",
            self.rms_residual, self.weighted_rms, self.max_abs_residual
        )?;
        writeln!(
            f,
            "  f-sum: fit {:.6e} eV^2, data {:.6e} eV^2 (ratio {:.6})",
            self.f_sum_ev2,
            self.data_f_sum_ev2,
            self.f_sum_ev2 / self.data_f_sum_ev2
        )?;
        write!(
            f,
            "  P_eff: fit {:.6}, data {:.6}",
            self.p_eff, self.data_p_eff
        )
    }
}

// ---- small dense linear algebra ----

/// Least squares `min |A x - y|` for the columns `cols` (each of length
/// `m`), by Householder QR. Assumes full column rank; returns `None` when a
/// column is numerically dependent.
fn lstsq(cols: &[&[f64]], y: &[f64]) -> Option<Vec<f64>> {
    let n = cols.len();
    let m = y.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if n > m {
        return None;
    }
    let mut a: Vec<Vec<f64>> = cols.iter().map(|c| c.to_vec()).collect();
    let mut b = y.to_vec();
    let mut diag = vec![0.0; n];
    for k in 0..n {
        let norm = a[k][k..].iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-13 {
            return None;
        }
        let alpha = if a[k][k] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = a[k][k..].to_vec();
        v[0] -= alpha;
        let vnorm2: f64 = v.iter().map(|x| x * x).sum();
        diag[k] = alpha;
        if vnorm2 > 0.0 {
            for j in (k + 1)..n {
                let dot: f64 = v.iter().zip(&a[j][k..]).map(|(p, q)| p * q).sum();
                let s = 2.0 * dot / vnorm2;
                for (i, vi) in v.iter().enumerate() {
                    a[j][k + i] -= s * vi;
                }
            }
            let dot: f64 = v.iter().zip(&b[k..]).map(|(p, q)| p * q).sum();
            let s = 2.0 * dot / vnorm2;
            for (i, vi) in v.iter().enumerate() {
                b[k + i] -= s * vi;
            }
        }
    }
    // back substitution with R (upper triangle: R[k][j] = a[j][k], diag[k])
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let mut s = b[k];
        for j in (k + 1)..n {
            s -= a[j][k] * x[j];
        }
        x[k] = s / diag[k];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Non-negative least squares, `min |A x - y|` with `x >= 0`: the
/// Lawson-Hanson active-set algorithm. The columns should be scaled to unit
/// norm.
fn nnls(cols: &[Vec<f64>], y: &[f64]) -> Vec<f64> {
    let n = cols.len();
    let mut x = vec![0.0; n];
    let mut passive = vec![false; n];
    let ynorm = y.iter().map(|v| v * v).sum::<f64>().sqrt();
    let tol = 1e-13 * ynorm.max(f64::MIN_POSITIVE);
    let gradient = |x: &[f64]| -> Vec<f64> {
        let mut r = y.to_vec();
        for (c, xc) in cols.iter().zip(x) {
            if *xc != 0.0 {
                for (ri, ci) in r.iter_mut().zip(c) {
                    *ri -= xc * ci;
                }
            }
        }
        cols.iter()
            .map(|c| c.iter().zip(&r).map(|(p, q)| p * q).sum::<f64>())
            .collect()
    };
    for _ in 0..(3 * n + 3) {
        let w = gradient(&x);
        let Some((t, wt)) = (0..n)
            .filter(|&j| !passive[j])
            .map(|j| (j, w[j]))
            .max_by(|a, b| a.1.total_cmp(&b.1))
        else {
            break;
        };
        if wt <= tol {
            break;
        }
        passive[t] = true;
        for _ in 0..(3 * n + 3) {
            let idx: Vec<usize> = (0..n).filter(|&j| passive[j]).collect();
            let sub: Vec<&[f64]> = idx.iter().map(|&j| cols[j].as_slice()).collect();
            let Some(s) = lstsq(&sub, y) else {
                // dependent column: drop the one just added and stop
                passive[t] = false;
                return x;
            };
            if s.iter().all(|&v| v > 0.0) {
                x = vec![0.0; n];
                for (&j, &v) in idx.iter().zip(&s) {
                    x[j] = v;
                }
                break;
            }
            let mut alpha = f64::INFINITY;
            for (&j, &sj) in idx.iter().zip(&s) {
                if sj <= 0.0 {
                    let a = x[j] / (x[j] - sj);
                    if a < alpha {
                        alpha = a;
                    }
                }
            }
            for (&j, &sj) in idx.iter().zip(&s) {
                x[j] += alpha * (sj - x[j]);
            }
            for &j in &idx {
                if x[j] <= 1e-14 {
                    x[j] = 0.0;
                    passive[j] = false;
                }
            }
        }
    }
    x
}

/// Solve the small symmetric positive system `M x = r` by Gaussian
/// elimination with partial pivoting.
fn solve(mut m: Vec<Vec<f64>>, mut r: Vec<f64>) -> Option<Vec<f64>> {
    let n = r.len();
    for k in 0..n {
        let p = (k..n).max_by(|&a, &b| m[a][k].abs().total_cmp(&m[b][k].abs()))?;
        if m[p][k].abs() < 1e-300 {
            return None;
        }
        m.swap(k, p);
        r.swap(k, p);
        for i in (k + 1)..n {
            let f = m[i][k] / m[k][k];
            for j in k..n {
                m[i][j] -= f * m[k][j];
            }
            r[i] -= f * r[k];
        }
    }
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let mut s = r[k];
        for j in (k + 1)..n {
            s -= m[k][j] * x[j];
        }
        x[k] = s / m[k][k];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

// ---- the fit ----

/// The local knot spacing of the sorted grid `w` at the energy `e`, eV: the
/// larger of the two gaps next to the knot nearest to `e` (the single gap of
/// an end knot; ties go to the lower knot). The lower bound of the fitted
/// widths (see [`fit_mermin_oscillators`]); a pure function of the knots.
fn local_knot_spacing(w: &[f64], e: f64) -> f64 {
    let n = w.len();
    if n < 2 {
        return 0.0;
    }
    let k = match w.binary_search_by(|x| x.total_cmp(&e)) {
        Ok(k) => k,
        Err(0) => 0,
        Err(i) if i >= n => n - 1,
        Err(i) => {
            if e - w[i - 1] <= w[i] - e {
                i - 1
            } else {
                i
            }
        }
    };
    let left = if k > 0 { w[k] - w[k - 1] } else { 0.0 };
    let right = if k + 1 < n { w[k + 1] - w[k] } else { 0.0 };
    left.max(right)
}

struct Problem<'a> {
    w: &'a [f64],
    y: Vec<f64>,
    row_weight: Vec<f64>,
}

impl Problem<'_> {
    /// The energy and width (eV) of the parameters `(ln E, ln γ)`, with the
    /// width raised to the local knot spacing at `E` if it is below it (the
    /// width floor of [`fit_mermin_oscillators`]). Above the floor the width
    /// is `exp(ln γ)` unchanged, bit for bit.
    fn energy_width(&self, ln_e: f64, ln_g: f64) -> (f64, f64) {
        let e = ln_e.exp();
        (e, ln_g.exp().max(local_knot_spacing(self.w, e)))
    }

    /// Project the stored width parameters `ln γ_j` of `theta` onto the width
    /// floor: raise each one that is below `ln` of the local knot spacing at
    /// its `E_j` up to it. Returns whether any parameter changed (if none is
    /// below the floor, `theta` is untouched, bit for bit). Without this, a
    /// `ln γ_j` below the floor would not enter the residual (see
    /// [`Self::energy_width`]), its Jacobian column would be zero and the
    /// fit could never move it again.
    fn raise_widths_to_floor(&self, theta: &mut [f64]) -> bool {
        let mut changed = false;
        for j in 0..theta.len() / 2 {
            let ln_floor = local_knot_spacing(self.w, theta[2 * j].exp()).ln();
            if theta[2 * j + 1] < ln_floor {
                theta[2 * j + 1] = ln_floor;
                changed = true;
            }
        }
        changed
    }

    /// The amplitudes `A` (eV², `>= 0`) and the weighted residual vector for
    /// the nonlinear parameters `theta = (ln E_j, ln γ_j)`.
    fn project(&self, theta: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = theta.len() / 2;
        let yw: Vec<f64> = self
            .y
            .iter()
            .zip(&self.row_weight)
            .map(|(y, s)| y * s)
            .collect();
        let mut cols = Vec::with_capacity(n);
        let mut norms = Vec::with_capacity(n);
        for j in 0..n {
            let (e, g) = self.energy_width(theta[2 * j], theta[2 * j + 1]);
            let mut c: Vec<f64> = self
                .w
                .iter()
                .zip(&self.row_weight)
                .map(|(&w, &s)| s * drude_basis(w, e, g))
                .collect();
            let norm = c.iter().map(|v| v * v).sum::<f64>().sqrt();
            let norm = if norm > 0.0 && norm.is_finite() {
                norm
            } else {
                1.0
            };
            for v in &mut c {
                *v /= norm;
            }
            cols.push(c);
            norms.push(norm);
        }
        let coef = nnls(&cols, &yw);
        let mut r: Vec<f64> = yw.iter().map(|v| -v).collect();
        for (c, k) in cols.iter().zip(&coef) {
            for (ri, ci) in r.iter_mut().zip(c) {
                *ri += k * ci;
            }
        }
        let amps = coef.iter().zip(&norms).map(|(k, n)| k / n).collect();
        (amps, r)
    }
}

fn default_start(elf: &OpticalElf, n: usize) -> Vec<DrudeLorentzOscillator> {
    let w = elf.energy_ev();
    let e = elf.elf_values();
    let mut peaks: Vec<(f64, f64)> = Vec::new();
    for i in 0..w.len() {
        let left = i == 0 || e[i] > e[i - 1];
        let right = i + 1 == w.len() || e[i] >= e[i + 1];
        if left && right && e[i] > 0.0 {
            peaks.push((e[i], w[i]));
        }
    }
    peaks.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.total_cmp(&b.1)));
    let mut energies: Vec<f64> = peaks.iter().take(n).map(|p| p.1).collect();
    // complete with log-spaced energies
    let (lo, hi) = (w[0].ln(), w[w.len() - 1].ln());
    let mut k = 1usize;
    while energies.len() < n {
        let cand = (lo + (hi - lo) * k as f64 / (n as f64 + 1.0)).exp();
        if !energies.iter().any(|&x| (x / cand).ln().abs() < 0.05) {
            energies.push(cand);
        }
        k += 1;
        if k > 10 * n + 10 {
            energies.push(cand * (1.0 + 0.1 * k as f64));
        }
    }
    energies.sort_by(f64::total_cmp);
    energies
        .into_iter()
        .map(|e| DrudeLorentzOscillator {
            strength_ev2: 0.0,
            energy_ev: e,
            width_ev: 0.25 * e,
        })
        .collect()
}

/// Fit a sum of Mermin oscillators to `elf` (see the module docs).
///
/// # Width floor
///
/// Each width `ħγ_j` is bounded below by the local knot spacing of the table
/// at `ħω_j`: the larger of the two gaps next to the tabulated energy nearest
/// to `ħω_j` (the single gap of an end knot). The stored parameter `ln γ_j`
/// is projected onto the bound (raised to `ln` of the floor) at the start and
/// after every accepted step, so a floored width keeps a nonzero Jacobian
/// column and can move off the floor when the data want it wider. The floor
/// is also applied inside every residual evaluation (trial and Jacobian
/// points) and to the reported oscillators, so all three agree. A fit in
/// which no width ever falls below the floor is unchanged bit for bit; of
/// the committed ELFs at the default options only C is affected (its start
/// width 0.05 eV at the 0.2 eV peak is raised to the 0.1 eV knot spacing,
/// and the fit ends within about 1e-7 of the unbounded one, in 146
/// iterations instead of 140). The reason: the residuals are evaluated at
/// the knots only, and a
/// Lorentzian narrower than the gap it sits in can fall between two knots,
/// where no residual constrains it while its amplitude, and with it the
/// f-sum `(π/2) A_j`, is free (#307: the 10- and 16-oscillator fits of the
/// committed Cu ELF reached 6771 and 16467 times the table's f-sum this way).
/// Such a width is not determined by the tabulated values. **This bound is
/// our own choice**: de Vera et al. (2022), section 2.1.1, give no fitting
/// algorithm, only the Kramers-Kronig and f-sum consistency check, which
/// [`MerminFit`] reports. The floor is a pure function of the knots and adds
/// no dependence on the thread count.
///
/// The floor is a step function of `ħω_j`: it jumps where the nearest knot
/// changes, halfway between two knots. For a width at the floor the residual
/// is therefore discontinuous in `ln ħω_j` there, and a central difference
/// straddling the jump gives a large Jacobian entry. The fit accepts only
/// strict decreases of the cost, so it does not move uphill across a jump,
/// but it can stop at one: if no damping gives a decrease, it reports
/// `converged = true` at a point that is a minimum only on one side of the
/// jump. Likewise, at a width on the floor where the data want a narrower
/// line, `converged = true` means a minimum subject to the bound.
///
/// The floor does not bound the f-sum against very wide oscillators: with
/// 10 or 16 oscillators the committed Cu ELF still gets one of width
/// 4e5 to 8e5 eV whose weight lies above the table, where the default
/// relative weighting barely sees it (f-sum 188 and 106 times the table's,
/// unconverged; #311).
pub fn fit_mermin_oscillators(elf: &OpticalElf, options: &MerminFitOptions) -> Result<MerminFit> {
    let w = elf.energy_ev();
    let y = elf.elf_values();
    let start = match &options.start {
        Some(s) => s.clone(),
        None => {
            if !(1..=MAX_OSCILLATORS).contains(&options.n_oscillators) {
                return Err(ElectronDataError::Invalid {
                    what: "Mermin fit",
                    reason: format!(
                        "number of oscillators must be in 1..={MAX_OSCILLATORS}, got {}",
                        options.n_oscillators
                    ),
                });
            }
            default_start(elf, options.n_oscillators)
        }
    };
    let n = start.len();
    if !(1..=MAX_OSCILLATORS).contains(&n) {
        return Err(ElectronDataError::Invalid {
            what: "Mermin fit",
            reason: format!("number of oscillators must be in 1..={MAX_OSCILLATORS}, got {n}"),
        });
    }
    for o in &start {
        let mut o = *o;
        o.strength_ev2 = 0.0; // amplitudes of the start are not used
        check_oscillator(&o)?;
    }
    if w.len() < 2 * n + 1 {
        return Err(ElectronDataError::Invalid {
            what: "Mermin fit",
            reason: format!(
                "{} tabulated points cannot determine {n} oscillators (need at least {})",
                w.len(),
                2 * n + 1
            ),
        });
    }
    let ymax = y.iter().cloned().fold(0.0_f64, f64::max);
    if ymax.is_nan() || ymax <= 0.0 {
        return Err(ElectronDataError::Invalid {
            what: "Mermin fit",
            reason: "the ELF is zero everywhere".into(),
        });
    }
    let row_weight: Vec<f64> = match options.weighting {
        FitWeighting::Uniform => vec![1.0; w.len()],
        FitWeighting::Relative { floor } => {
            if !(floor > 0.0 && floor <= 1.0) {
                return Err(ElectronDataError::Invalid {
                    what: "Mermin fit weighting",
                    reason: format!("relative floor must be in (0, 1], got {floor}"),
                });
            }
            y.iter().map(|v| 1.0 / (v + floor * ymax)).collect()
        }
    };
    let problem = Problem {
        w,
        y: y.to_vec(),
        row_weight,
    };

    let mut theta: Vec<f64> = start
        .iter()
        .flat_map(|o| [o.energy_ev.ln(), o.width_ev.ln()])
        .collect();
    // the width floor is enforced on the stored parameters (a projection),
    // at the start and after every accepted step
    problem.raise_widths_to_floor(&mut theta);
    let np = theta.len();
    let (mut amps, mut r) = problem.project(&theta);
    let mut cost = 0.5 * r.iter().map(|v| v * v).sum::<f64>();
    let scale = 0.5
        * problem
            .y
            .iter()
            .zip(&problem.row_weight)
            .map(|(y, s)| (y * s).powi(2))
            .sum::<f64>();
    let mut lambda = 1e-3;
    let mut iterations = 0;
    let mut converged = false;
    let h = 1e-6;

    while iterations < options.max_iterations {
        if cost <= 1e-28 * scale {
            converged = true;
            break;
        }
        iterations += 1;
        // central-difference Jacobian of the projected residual
        let mut jac: Vec<Vec<f64>> = Vec::with_capacity(np); // column-major
        for p in 0..np {
            let mut tp = theta.clone();
            tp[p] += h;
            let (_, rp) = problem.project(&tp);
            tp[p] -= 2.0 * h;
            let (_, rm) = problem.project(&tp);
            jac.push(
                rp.iter()
                    .zip(&rm)
                    .map(|(a, b)| (a - b) / (2.0 * h))
                    .collect(),
            );
        }
        let jtj: Vec<Vec<f64>> = (0..np)
            .map(|a| {
                (0..np)
                    .map(|b| jac[a].iter().zip(&jac[b]).map(|(p, q)| p * q).sum())
                    .collect()
            })
            .collect();
        let g: Vec<f64> = (0..np)
            .map(|a| jac[a].iter().zip(&r).map(|(p, q)| p * q).sum())
            .collect();
        let gmax = g.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        if gmax <= 1e-14 * scale.sqrt().max(f64::MIN_POSITIVE) {
            converged = true;
            break;
        }
        // inner loop: increase damping until the step reduces the cost
        let mut accepted = false;
        for _ in 0..40 {
            let mut m = jtj.clone();
            for (a, row) in m.iter_mut().enumerate() {
                row[a] += lambda * jtj[a][a].max(1e-12 * gmax.max(1e-300));
            }
            let rhs: Vec<f64> = g.iter().map(|v| -v).collect();
            if let Some(mut delta) = solve(m, rhs) {
                for d in &mut delta {
                    *d = d.clamp(-1.0, 1.0);
                }
                let trial: Vec<f64> = theta.iter().zip(&delta).map(|(t, d)| t + d).collect();
                let (a2, r2) = problem.project(&trial);
                let c2 = 0.5 * r2.iter().map(|v| v * v).sum::<f64>();
                if c2.is_finite() && c2 < cost {
                    let step = delta.iter().fold(0.0_f64, |m, d| m.max(d.abs()));
                    let drop = (cost - c2) / cost;
                    theta = trial;
                    amps = a2;
                    r = r2;
                    cost = c2;
                    if problem.raise_widths_to_floor(&mut theta) {
                        // the residual at the projected point (equal to the
                        // trial's up to rounding: the trial already used the
                        // floor for these widths)
                        (amps, r) = problem.project(&theta);
                        cost = 0.5 * r.iter().map(|v| v * v).sum::<f64>();
                    }
                    lambda = (lambda * 0.1).max(1e-15);
                    accepted = true;
                    if step < 1e-13 || drop < 1e-15 {
                        converged = true;
                    }
                    break;
                }
            }
            lambda *= 10.0;
        }
        if !accepted {
            // no downhill step at any damping: at a (local) minimum
            converged = true;
            break;
        }
        if converged {
            break;
        }
    }

    let mut oscillators: Vec<DrudeLorentzOscillator> = (0..n)
        .map(|j| {
            let (energy_ev, width_ev) = problem.energy_width(theta[2 * j], theta[2 * j + 1]);
            DrudeLorentzOscillator {
                strength_ev2: amps[j],
                energy_ev,
                width_ev,
            }
        })
        .collect();
    oscillators.sort_by(|a, b| {
        a.energy_ev
            .total_cmp(&b.energy_ev)
            .then(a.width_ev.total_cmp(&b.width_ev))
    });

    let mut fit = MerminFit {
        oscillators,
        energy_ev: w.to_vec(),
        data: y.to_vec(),
        model: Vec::new(),
        residuals: Vec::new(),
        weighting: options.weighting,
        weighted_rms: (r.iter().map(|v| v * v).sum::<f64>() / r.len() as f64).sqrt(),
        rms_residual: 0.0,
        max_abs_residual: 0.0,
        iterations,
        converged,
        f_sum_ev2: 0.0,
        p_eff: 0.0,
        data_f_sum_ev2: 0.0,
        data_p_eff: 0.0,
    };
    fit.model = w.iter().map(|&e| fit.model_elf(e)).collect();
    fit.residuals = fit.model.iter().zip(y).map(|(m, d)| m - d).collect();
    fit.rms_residual = (fit.residuals.iter().map(|v| v * v).sum::<f64>() / w.len() as f64).sqrt();
    fit.max_abs_residual = fit.residuals.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    fit.f_sum_ev2 = FRAC_PI_2 * fit.oscillators.iter().map(|o| o.strength_ev2).sum::<f64>();
    fit.p_eff = fit
        .oscillators
        .iter()
        .map(|o| o.strength_ev2 / (o.energy_ev * o.energy_ev))
        .sum();
    let rules = SumRuleReport::new(elf);
    fit.data_f_sum_ev2 = rules.f_sum_ev2;
    fit.data_p_eff = rules.p_eff;
    Ok(fit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electron::inelastic::drude::DrudeLorentz;

    fn osc(a: f64, e: f64, g: f64) -> DrudeLorentzOscillator {
        DrudeLorentzOscillator {
            strength_ev2: a,
            energy_ev: e,
            width_ev: g,
        }
    }

    #[test]
    fn lstsq_solves_an_overdetermined_system() {
        let c0 = [1.0, 1.0, 1.0, 1.0];
        let c1 = [0.0, 1.0, 2.0, 3.0];
        let y = [1.0, 3.0, 5.0, 7.0]; // 1 + 2x exactly
        let x = lstsq(&[&c0, &c1], &y).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-12 && (x[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn nnls_clips_negative_coefficients() {
        let c0 = vec![1.0, 0.0, 0.0];
        let c1 = vec![0.0, 1.0, 0.0];
        // y = 2 c0 - 3 c1: unconstrained solution has a negative coefficient
        let x = nnls(&[c0, c1], &[2.0, -3.0, 0.0]);
        assert!((x[0] - 2.0).abs() < 1e-12);
        assert_eq!(x[1], 0.0);
    }

    #[test]
    fn nnls_matches_least_squares_when_the_solution_is_positive() {
        let c0 = vec![1.0, 1.0, 0.5, 0.0];
        let c1 = vec![0.0, 1.0, 1.0, 1.0];
        let y: Vec<f64> = c0.iter().zip(&c1).map(|(a, b)| 1.5 * a + 0.7 * b).collect();
        let x = nnls(&[c0, c1], &y);
        assert!((x[0] - 1.5).abs() < 1e-12 && (x[1] - 0.7).abs() < 1e-12);
    }

    fn synthetic() -> (Vec<DrudeLorentzOscillator>, OpticalElf) {
        let truth = vec![
            osc(120.0, 9.0, 3.0),
            osc(400.0, 22.0, 7.0),
            osc(900.0, 55.0, 25.0),
        ];
        let elf = DrudeLorentz::new(truth.clone())
            .unwrap()
            .to_optical_elf("synthetic", 0.5, 2000.0, 400)
            .unwrap();
        (truth, elf)
    }

    #[test]
    fn recovers_three_known_oscillators_from_a_perturbed_start() {
        let (truth, elf) = synthetic();
        let start = vec![
            osc(1.0, 9.0 * 1.15, 3.0 * 0.8),
            osc(1.0, 22.0 * 0.88, 7.0 * 1.2),
            osc(1.0, 55.0 * 1.1, 25.0 * 0.85),
        ];
        let fit = fit_mermin_oscillators(
            &elf,
            &MerminFitOptions {
                start: Some(start),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(fit.converged, "{fit}");
        assert_eq!(fit.oscillators.len(), 3);
        for (got, want) in fit.oscillators.iter().zip(&truth) {
            for (name, g, w) in [
                ("A", got.strength_ev2, want.strength_ev2),
                ("E", got.energy_ev, want.energy_ev),
                ("gamma", got.width_ev, want.width_ev),
            ] {
                assert!(((g - w) / w).abs() < 1e-3, "{name}: {g} vs {w}\n{fit}");
            }
            assert!(got.strength_ev2 >= 0.0);
        }
        // the sum rules of the fit are the closed forms of the truth
        let want = DrudeLorentz::new(truth).unwrap();
        assert!(((fit.f_sum_ev2 - want.f_sum_ev2()) / want.f_sum_ev2()).abs() < 1e-3);
        assert!(((fit.p_eff - want.p_eff()) / want.p_eff()).abs() < 1e-3);
        assert!(fit.max_abs_residual < 1e-6, "{fit}");
        assert_eq!(fit.residuals.len(), elf.energy_ev().len());
    }

    #[test]
    fn default_start_also_recovers_the_synthetic_elf() {
        let (truth, elf) = synthetic();
        let fit = fit_mermin_oscillators(&elf, &MerminFitOptions::default()).unwrap();
        assert!(fit.rms_residual < 1e-4, "{fit}");
        let want = DrudeLorentz::new(truth).unwrap();
        assert!(((fit.f_sum_ev2 - want.f_sum_ev2()) / want.f_sum_ev2()).abs() < 1e-2);
    }

    #[test]
    fn amplitudes_stay_non_negative_for_a_surplus_of_oscillators() {
        let (_, elf) = synthetic();
        let fit = fit_mermin_oscillators(
            &elf,
            &MerminFitOptions {
                n_oscillators: 5,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(fit.oscillators.iter().all(|o| o.strength_ev2 >= 0.0));
        assert!(fit.rms_residual < 1e-3, "{fit}");
    }

    #[test]
    fn fit_is_bitwise_deterministic_across_threads() {
        let (_, elf) = synthetic();
        let opts = MerminFitOptions::default();
        let base = fit_mermin_oscillators(&elf, &opts).unwrap();
        let handles: Vec<_> = (0..4)
            .map(|_| {
                let (elf, opts) = (elf.clone(), opts.clone());
                std::thread::spawn(move || fit_mermin_oscillators(&elf, &opts).unwrap())
            })
            .collect();
        for h in handles {
            let other = h.join().unwrap();
            for (a, b) in base.oscillators.iter().zip(&other.oscillators) {
                assert_eq!(a.strength_ev2.to_bits(), b.strength_ev2.to_bits());
                assert_eq!(a.energy_ev.to_bits(), b.energy_ev.to_bits());
                assert_eq!(a.width_ev.to_bits(), b.width_ev.to_bits());
            }
            assert_eq!(base.iterations, other.iterations);
        }
    }

    #[test]
    fn invalid_options_are_rejected() {
        let (_, elf) = synthetic();
        for n in [0, MAX_OSCILLATORS + 1] {
            assert!(fit_mermin_oscillators(
                &elf,
                &MerminFitOptions {
                    n_oscillators: n,
                    ..Default::default()
                }
            )
            .is_err());
        }
        assert!(fit_mermin_oscillators(
            &elf,
            &MerminFitOptions {
                weighting: FitWeighting::Relative { floor: 0.0 },
                ..Default::default()
            }
        )
        .is_err());
        assert!(fit_mermin_oscillators(
            &elf,
            &MerminFitOptions {
                start: Some(vec![osc(1.0, -1.0, 1.0)]),
                ..Default::default()
            }
        )
        .is_err());
    }

    #[test]
    fn local_knot_spacing_is_the_larger_gap_next_to_the_nearest_knot() {
        // non-uniform grid: gaps 1, 2, 4, 8
        let w = [1.0, 2.0, 4.0, 8.0, 16.0];
        // interior knot and points nearest to it: max(2, 4)
        assert_eq!(local_knot_spacing(&w, 4.0), 4.0);
        assert_eq!(local_knot_spacing(&w, 3.5), 4.0);
        assert_eq!(local_knot_spacing(&w, 5.9), 4.0);
        // nearer to 8: max(4, 8)
        assert_eq!(local_knot_spacing(&w, 6.1), 8.0);
        // a tie goes to the lower knot (2: max(1, 2))
        assert_eq!(local_knot_spacing(&w, 3.0), 2.0);
        // end knots use their single gap, also outside the table
        assert_eq!(local_knot_spacing(&w, 1.0), 1.0);
        assert_eq!(local_knot_spacing(&w, 0.1), 1.0);
        assert_eq!(local_knot_spacing(&w, 16.0), 8.0);
        assert_eq!(local_knot_spacing(&w, 1e6), 8.0);
    }

    #[test]
    fn widths_are_not_fitted_below_the_knot_spacing() {
        // a start far narrower than the grid: the reported widths are at
        // least the local knot spacing
        let (_, elf) = synthetic();
        let w = elf.energy_ev();
        let start = vec![osc(1.0, 9.0, 1e-4), osc(1.0, 55.0, 1e-4)];
        let fit = fit_mermin_oscillators(
            &elf,
            &MerminFitOptions {
                start: Some(start),
                max_iterations: 5,
                ..Default::default()
            },
        )
        .unwrap();
        for o in &fit.oscillators {
            assert!(o.width_ev >= local_knot_spacing(w, o.energy_ev), "{fit}");
        }
    }

    #[test]
    fn data_sum_rules_are_reported() {
        let (truth, elf) = synthetic();
        let fit = fit_mermin_oscillators(&elf, &MerminFitOptions::default()).unwrap();
        let want = DrudeLorentz::new(truth).unwrap();
        // the table ends at 2 keV, so the data's sums are a little short
        assert!(fit.data_f_sum_ev2 > 0.9 * want.f_sum_ev2());
        assert!(fit.data_p_eff > 0.9 * want.p_eff());
        assert!(fit.summary().starts_with("3 oscillators"));
    }
}
