//! Per-pair multiplicative correction of an electronic-stopping model, fitted
//! to measured points, with a bootstrap uncertainty and held-out residuals.
//!
//! A [`Correction`] wraps a base [`ElectronicStopping`] model and multiplies
//! its stopping cross section by a smooth factor `g(E)` for **one** ion–target
//! pair. The factor is
//!
//! ```text
//! g(E) = exp( sum_{k=0}^{d} c_k T_k(x) ),   x = (2 ln E - ln E_lo - ln E_hi) / (ln E_hi - ln E_lo)
//! ```
//!
//! where `T_k` are Chebyshev polynomials of the first kind (J. C. Mason and
//! D. C. Handscomb, *Chebyshev Polynomials*, CRC Press (2003), ch. 2) and
//! `[E_lo, E_hi]` is the energy span of the training points. The exponential
//! keeps `g` positive. The coefficients minimise
//! `sum_i w_i (ln(S_i / B_i) - sum_k c_k T_k(x_i))^2`, where `S_i` is the
//! measured stopping, `B_i` the base model at the same point and
//! `w_i = 1 / (sigma_i / S_i)^2` when **every** training point states an
//! uncertainty `sigma_i > 0`, otherwise `w_i = 1`. The normal equations are
//! solved by Gaussian elimination with partial pivoting (G. H. Golub and
//! C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins (2013),
//! sec. 3.4 and 4.2). The form of the correction (a low-order polynomial of
//! `ln E` for the log of the ratio) is our own choice; no published fit is
//! reproduced.
//!
//! # Where the correction applies
//!
//! Only inside `[E_lo, E_hi]` of the training points, and only for the ion
//! (Z and mass) and target element it was fitted to. Anywhere else
//! [`Correction::stopping`] returns the base model unchanged, because a
//! polynomial extrapolated beyond its data has no physical basis.
//!
//! # Uncertainty
//!
//! The nonparametric bootstrap of B. Efron, Ann. Statist. 7, 1 (1979): the
//! training points are resampled with replacement `replicates` times, the
//! correction is refitted to each resample, and the pointwise band at energy
//! `E` is the central interval of `g_b(E)` over the replicates (linear
//! interpolation between order statistics). A resample whose fit is singular
//! (for example fewer distinct energies than coefficients) is discarded and
//! counted in [`CorrectionFit::failed_replicates`]. The band covers the
//! scatter of the points about the fitted form; it does not include any
//! systematic error of the measurement or any bias of the polynomial form.
//! If the measured points carry no stated uncertainty (sigma 0, as in tables
//! read off a smoothed curve), the band says nothing about measurement error.
//!
//! # Reproducibility
//!
//! The held-out split uses the stream `(seed, 0)` of [`crate::rng::stream`]
//! and replicate `b` uses `(seed, b + 1)`. Replicates run on the current rayon
//! pool and are collected in index order; every sum is sequential, so the
//! result is bit-identical for a given seed whatever the thread count.
//!
//! # Held-out points
//!
//! A deterministic, seeded subset of the interior points (never the lowest and
//! highest energy, so they stay inside the fitted span) is left out of the fit.
//! [`CorrectionFit::held_out_report`] gives the residuals of the base model and
//! of the corrected model on exactly those points, side by side, whether or not
//! the correction helps.

use super::dataset::{StoppingDataset, Target};
use super::{to_ev_1e15_cm2, ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::rng::stream;
use rand_core::Rng;
use rayon::prelude::*;

/// Name returned by [`Correction::name`].
pub const CORRECTED_MODEL_NAME: &str = "corrected";

/// Settings of [`fit_correction`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FitSettings {
    /// Degree `d` of the Chebyshev polynomial; `d + 1` coefficients.
    pub degree: usize,
    /// Number of interior points left out of the fit and reported as held out.
    pub held_out: usize,
    /// Number of bootstrap replicates (0 gives no band).
    pub replicates: usize,
    /// Seed of the split and of the bootstrap.
    pub seed: u64,
}

/// Statistics of `model - measured` over a set of points.
///
/// Signs follow [`super::dataset::PointResidual`]: positive when the model is
/// above the measurement. With no points the statistics are NaN.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResidualStats {
    /// Number of points.
    pub n: usize,
    /// Mean of `model / measured - 1`.
    pub mean_relative: f64,
    /// Root mean square of `model / measured - 1`.
    pub rms_relative: f64,
    /// Mean of `model - measured`, eV·10⁻¹⁵ cm².
    pub mean_absolute_ev_1e15_cm2: f64,
    /// Root mean square of `model - measured`, eV·10⁻¹⁵ cm².
    pub rms_absolute_ev_1e15_cm2: f64,
}

impl ResidualStats {
    /// Statistics of `(model, measured)` pairs, summed sequentially in order.
    pub fn of(pairs: &[(f64, f64)]) -> Self {
        let n = pairs.len();
        if n == 0 {
            return Self {
                n,
                mean_relative: f64::NAN,
                rms_relative: f64::NAN,
                mean_absolute_ev_1e15_cm2: f64::NAN,
                rms_absolute_ev_1e15_cm2: f64::NAN,
            };
        }
        let (mut sr, mut sr2, mut sa, mut sa2) = (0.0, 0.0, 0.0, 0.0);
        for &(m, s) in pairs {
            let r = m / s - 1.0;
            let a = m - s;
            sr += r;
            sr2 += r * r;
            sa += a;
            sa2 += a * a;
        }
        let nf = n as f64;
        Self {
            n,
            mean_relative: sr / nf,
            rms_relative: (sr2 / nf).sqrt(),
            mean_absolute_ev_1e15_cm2: sa / nf,
            rms_absolute_ev_1e15_cm2: (sa2 / nf).sqrt(),
        }
    }
}

/// One held-out point with the base and corrected model values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeldOutPoint {
    /// Index of the point in the dataset.
    pub index: usize,
    /// Projectile energy, eV.
    pub energy_ev: f64,
    /// Measured stopping, eV·10⁻¹⁵ cm².
    pub measured_ev_1e15_cm2: f64,
    /// Base model, eV·10⁻¹⁵ cm².
    pub base_ev_1e15_cm2: f64,
    /// Corrected model, eV·10⁻¹⁵ cm².
    pub corrected_ev_1e15_cm2: f64,
}

/// Residuals of the base and the corrected model on the held-out points.
#[derive(Debug, Clone, PartialEq)]
pub struct HeldOutReport {
    /// The held-out points, in dataset order.
    pub points: Vec<HeldOutPoint>,
    /// Residual statistics of the uncorrected base model on those points.
    pub base: ResidualStats,
    /// Residual statistics of the corrected model on the same points.
    pub corrected: ResidualStats,
}

impl HeldOutReport {
    /// True if the correction lowered the held-out relative rms. A correction
    /// that does not is still reported; it is never hidden.
    pub fn reduces_rms_relative(&self) -> bool {
        self.corrected.rms_relative < self.base.rms_relative
    }
}

/// A fitted correction: coefficients, bootstrap replicates and the split.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrectionFit {
    ion: Ion,
    target_z: u8,
    e_lo_ev: f64,
    e_hi_ev: f64,
    coefficients: Vec<f64>,
    replicates: Vec<Vec<f64>>,
    failed_replicates: usize,
    weighted: bool,
    settings: FitSettings,
    train: Vec<usize>,
    held_out: Vec<usize>,
    held_out_report: HeldOutReport,
    train_base: ResidualStats,
    train_corrected: ResidualStats,
}

impl CorrectionFit {
    /// The projectile the correction was fitted for.
    pub fn ion(&self) -> &Ion {
        &self.ion
    }

    /// Atomic number of the target element.
    pub fn target_z(&self) -> u8 {
        self.target_z
    }

    /// Energy span `(lowest, highest)` of the training points, eV; the
    /// correction applies inside it only (both ends included).
    pub fn span_ev(&self) -> (f64, f64) {
        (self.e_lo_ev, self.e_hi_ev)
    }

    /// Chebyshev coefficients `c_0..=c_d` of `ln g`.
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// The settings the fit was made with (including the seed).
    pub fn settings(&self) -> FitSettings {
        self.settings
    }

    /// Dataset indices of the training points, ascending.
    pub fn train_indices(&self) -> &[usize] {
        &self.train
    }

    /// Dataset indices of the held-out points, ascending.
    pub fn held_out_indices(&self) -> &[usize] {
        &self.held_out
    }

    /// True if the fit used weights `1 / (sigma / S)^2`, which requires every
    /// training point to state an uncertainty; otherwise it is unweighted.
    pub fn weighted(&self) -> bool {
        self.weighted
    }

    /// Bootstrap resamples discarded because their fit was singular.
    pub fn failed_replicates(&self) -> usize {
        self.failed_replicates
    }

    /// Number of bootstrap replicates kept.
    pub fn n_replicates(&self) -> usize {
        self.replicates.len()
    }

    /// Residuals of the base and corrected model on the held-out points.
    pub fn held_out_report(&self) -> &HeldOutReport {
        &self.held_out_report
    }

    /// Residual statistics of the base model on the training points.
    pub fn train_base(&self) -> ResidualStats {
        self.train_base
    }

    /// Residual statistics of the corrected model on the training points.
    pub fn train_corrected(&self) -> ResidualStats {
        self.train_corrected
    }

    /// True if `energy_ev` is inside the span the correction applies to.
    pub fn in_span(&self, energy_ev: f64) -> bool {
        energy_ev >= self.e_lo_ev && energy_ev <= self.e_hi_ev
    }

    /// The correction factor `g(E)` of the fit, or `None` outside the span.
    pub fn factor(&self, energy_ev: f64) -> Option<f64> {
        self.in_span(energy_ev)
            .then(|| eval_factor(&self.coefficients, self.e_lo_ev, self.e_hi_ev, energy_ev))
    }

    /// Pointwise bootstrap band of `g(E)`: the central interval holding
    /// `level` (for example 0.68) of the replicates, as `(low, high)`.
    ///
    /// `None` outside the span, with no replicates, or for `level` outside
    /// `(0, 1)`.
    pub fn band(&self, energy_ev: f64, level: f64) -> Option<(f64, f64)> {
        if !self.in_span(energy_ev) || self.replicates.is_empty() || !(level > 0.0 && level < 1.0) {
            return None;
        }
        let mut g: Vec<f64> = self
            .replicates
            .iter()
            .map(|c| eval_factor(c, self.e_lo_ev, self.e_hi_ev, energy_ev))
            .collect();
        g.sort_by(f64::total_cmp);
        let tail = 0.5 * (1.0 - level);
        Some((quantile(&g, tail), quantile(&g, 1.0 - tail)))
    }
}

/// An [`ElectronicStopping`] model multiplied by a fitted [`CorrectionFit`].
///
/// For the fitted ion and target element and an energy inside the fitted span
/// the base value is multiplied by `g(E)`; in every other case the base model
/// is returned unchanged.
#[derive(Debug, Clone)]
pub struct Correction<M> {
    base: M,
    fit: CorrectionFit,
}

impl<M: ElectronicStopping> Correction<M> {
    /// Wrap `base` with `fit`.
    pub fn new(base: M, fit: CorrectionFit) -> Self {
        Self { base, fit }
    }

    /// The wrapped base model.
    pub fn base(&self) -> &M {
        &self.base
    }

    /// The fit.
    pub fn fit(&self) -> &CorrectionFit {
        &self.fit
    }
}

impl<M: ElectronicStopping> ElectronicStopping for Correction<M> {
    fn name(&self) -> &'static str {
        CORRECTED_MODEL_NAME
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        let s = self.base.stopping(ion, target_z, energy_ev)?;
        let same_pair = ion.z() == self.fit.ion.z()
            && ion.mass_amu().to_bits() == self.fit.ion.mass_amu().to_bits()
            && target_z == self.fit.target_z;
        match self.fit.factor(energy_ev) {
            Some(g) if same_pair => Ok(s * g),
            _ => Ok(s),
        }
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        self.base.validity(ion)
    }
}

fn invalid(msg: impl Into<String>) -> StoppingError {
    StoppingError::InvalidFit(msg.into())
}

/// Chebyshev abscissa of `energy_ev` for the span `[lo, hi]`.
fn abscissa(lo: f64, hi: f64, energy_ev: f64) -> f64 {
    let (a, b) = (lo.ln(), hi.ln());
    ((2.0 * energy_ev.ln() - a - b) / (b - a)).clamp(-1.0, 1.0)
}

fn eval_factor(coefficients: &[f64], lo: f64, hi: f64, energy_ev: f64) -> f64 {
    let x = abscissa(lo, hi, energy_ev);
    let (mut t_prev, mut t) = (1.0, x);
    let mut sum = coefficients[0];
    if coefficients.len() > 1 {
        sum += coefficients[1] * x;
    }
    for &c in coefficients.iter().skip(2) {
        let t_next = 2.0 * x * t - t_prev;
        sum += c * t_next;
        t_prev = t;
        t = t_next;
    }
    sum.exp()
}

/// Quantile `q` in `[0, 1]` of a sorted slice, linear between order statistics.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    let pos = q * (sorted.len() - 1) as f64;
    let i = pos.floor() as usize;
    let j = (i + 1).min(sorted.len() - 1);
    let f = pos - i as f64;
    sorted[i] + f * (sorted[j] - sorted[i])
}

/// One training point in the fit: abscissa data and weight.
#[derive(Clone, Copy)]
struct Sample {
    energy_ev: f64,
    log_ratio: f64,
    weight: f64,
}

/// Weighted least squares of `log_ratio` on `T_0..=T_d(x)`; `None` if singular.
fn solve(samples: &[Sample], degree: usize, lo: f64, hi: f64) -> Option<Vec<f64>> {
    let m = degree + 1;
    let mut a = vec![vec![0.0; m + 1]; m];
    let mut basis = vec![0.0; m];
    for s in samples {
        let x = abscissa(lo, hi, s.energy_ev);
        basis[0] = 1.0;
        if m > 1 {
            basis[1] = x;
        }
        for k in 2..m {
            basis[k] = 2.0 * x * basis[k - 1] - basis[k - 2];
        }
        for r in 0..m {
            for c in 0..m {
                a[r][c] += s.weight * basis[r] * basis[c];
            }
            a[r][m] += s.weight * basis[r] * s.log_ratio;
        }
    }
    let scale = (0..m).map(|i| a[i][i]).fold(0.0, f64::max);
    for col in 0..m {
        let mut piv = col;
        for r in col + 1..m {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs().is_nan() || a[piv][col].abs() <= 1e-12 * scale {
            return None;
        }
        a.swap(col, piv);
        let pivot_row = a[col].clone();
        for row in a.iter_mut().skip(col + 1) {
            let f = row[col] / pivot_row[col];
            for (x, p) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *x -= f * p;
            }
        }
    }
    let mut c = vec![0.0; m];
    for r in (0..m).rev() {
        let mut v = a[r][m];
        for k in r + 1..m {
            v -= a[r][k] * c[k];
        }
        c[r] = v / a[r][r];
    }
    c.iter().all(|v| v.is_finite()).then_some(c)
}

/// Uniform index in `0..n` from one 64-bit draw (multiply-shift; the bias is
/// below `n / 2^64`).
fn draw_index(rng: &mut impl Rng, n: usize) -> usize {
    ((u128::from(rng.next_u64()) * n as u128) >> 64) as usize
}

/// Fit a [`Correction`] of `base` to a single-system element-target dataset.
///
/// The dataset must hold points of one ion (same Z and mass) on one element
/// target, at no fewer than `degree + 2 + held_out` points with at least two
/// distinct energies, and `base` must return a positive value at every point.
/// Anything else is an [`StoppingError::InvalidFit`] (or the base model's own
/// error), so no point is dropped silently.
pub fn fit_correction(
    base: &dyn ElectronicStopping,
    dataset: &StoppingDataset,
    settings: FitSettings,
) -> Result<CorrectionFit, StoppingError> {
    let pts = dataset.points();
    let n = pts.len();
    let first = pts
        .first()
        .ok_or_else(|| invalid("dataset has no points"))?;
    let Target::Element(target_z) = *first.target() else {
        return Err(invalid("compound targets are not supported"));
    };
    let ion = *first.ion();
    for p in pts {
        if p.ion().z() != ion.z()
            || p.ion().mass_amu().to_bits() != ion.mass_amu().to_bits()
            || p.target().label() != first.target().label()
        {
            return Err(invalid(format!(
                "line {}: the dataset must hold a single ion-target system",
                p.line()
            )));
        }
    }
    let mut base_values = Vec::with_capacity(n);
    for p in pts {
        let b = to_ev_1e15_cm2(base.stopping(&ion, target_z, p.energy_ev())?);
        if !(b.is_finite() && b > 0.0) {
            return Err(invalid(format!(
                "line {}: base model value {b} is not positive",
                p.line()
            )));
        }
        base_values.push(b);
    }

    // Indices by energy; the lowest and highest stay in training.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        pts[i]
            .energy_ev()
            .total_cmp(&pts[j].energy_ev())
            .then(i.cmp(&j))
    });
    if n < settings.degree + 2 + settings.held_out {
        return Err(invalid(format!(
            "{n} points are too few for degree {} with {} held out",
            settings.degree, settings.held_out
        )));
    }
    let interior = &order[1..n - 1];
    if settings.held_out > interior.len() {
        return Err(invalid("more held-out points than interior points"));
    }
    let mut shuffled = interior.to_vec();
    let mut rng = stream(settings.seed, 0);
    for i in (1..shuffled.len()).rev() {
        shuffled.swap(i, draw_index(&mut rng, i + 1));
    }
    let mut held_out: Vec<usize> = shuffled[..settings.held_out].to_vec();
    held_out.sort_unstable();
    let train: Vec<usize> = (0..n)
        .filter(|i| held_out.binary_search(i).is_err())
        .collect();

    let e_lo = pts[order[0]].energy_ev();
    let e_hi = pts[order[n - 1]].energy_ev();
    if e_hi <= e_lo {
        return Err(invalid("the points span a single energy"));
    }
    let weighted = train.iter().all(|&i| pts[i].sigma_ev_1e15_cm2() > 0.0);
    let samples: Vec<Sample> = train
        .iter()
        .map(|&i| {
            let meas = pts[i].s_e_ev_1e15_cm2();
            let sigma_log = pts[i].sigma_ev_1e15_cm2() / meas;
            Sample {
                energy_ev: pts[i].energy_ev(),
                log_ratio: (meas / base_values[i]).ln(),
                weight: if weighted {
                    1.0 / (sigma_log * sigma_log)
                } else {
                    1.0
                },
            }
        })
        .collect();
    let coefficients = solve(&samples, settings.degree, e_lo, e_hi)
        .ok_or_else(|| invalid("the least-squares system is singular"))?;

    let outcomes: Vec<Option<Vec<f64>>> = (0..settings.replicates)
        .into_par_iter()
        .map(|b| {
            let mut rng = stream(settings.seed, b as u64 + 1);
            let resample: Vec<Sample> = (0..samples.len())
                .map(|_| samples[draw_index(&mut rng, samples.len())])
                .collect();
            solve(&resample, settings.degree, e_lo, e_hi)
        })
        .collect();
    let failed_replicates = outcomes.iter().filter(|o| o.is_none()).count();
    let replicates: Vec<Vec<f64>> = outcomes.into_iter().flatten().collect();

    let corrected_at =
        |i: usize| base_values[i] * eval_factor(&coefficients, e_lo, e_hi, pts[i].energy_ev());
    let pairs = |idx: &[usize], corrected: bool| -> Vec<(f64, f64)> {
        idx.iter()
            .map(|&i| {
                let m = if corrected {
                    corrected_at(i)
                } else {
                    base_values[i]
                };
                (m, pts[i].s_e_ev_1e15_cm2())
            })
            .collect()
    };
    let held_out_report = HeldOutReport {
        points: held_out
            .iter()
            .map(|&i| HeldOutPoint {
                index: i,
                energy_ev: pts[i].energy_ev(),
                measured_ev_1e15_cm2: pts[i].s_e_ev_1e15_cm2(),
                base_ev_1e15_cm2: base_values[i],
                corrected_ev_1e15_cm2: corrected_at(i),
            })
            .collect(),
        base: ResidualStats::of(&pairs(&held_out, false)),
        corrected: ResidualStats::of(&pairs(&held_out, true)),
    };

    let train_base = ResidualStats::of(&pairs(&train, false));
    let train_corrected = ResidualStats::of(&pairs(&train, true));
    Ok(CorrectionFit {
        ion,
        target_z,
        e_lo_ev: e_lo,
        e_hi_ev: e_hi,
        coefficients,
        replicates,
        failed_replicates,
        weighted,
        settings,
        train_base,
        train_corrected,
        train,
        held_out,
        held_out_report,
    })
}
