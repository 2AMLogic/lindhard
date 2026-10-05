//! Pearson type IV and dual-Pearson representations of range profiles.
//!
//! # Pearson IV
//!
//! The Pearson system (K. Pearson, "Contributions to the mathematical theory
//! of evolution. II. Skew variation in homogeneous material", Phil. Trans.
//! R. Soc. Lond. A 186 (1895) 343) is the family of densities with
//! `d ln f / dx = -(x - c) / (b0 + b1 x + b2 x^2)`; its type IV member is
//! fixed by the first four moments. W. K. Hofker, "Implantation of boron in
//! silicon", Philips Res. Repts. Suppl. No. 8 (1975), introduced it for ion
//! implantation profiles, characterised by `Rp`, `ΔRp`, `γ` and `β`.
//!
//! We use the explicit form and moment inversion of J. Heinrich, "A guide to
//! the Pearson type IV distribution", CDF/MEMO/STATISTICS/PUBLIC/6820 (2004):
//!
//! ```text
//! f(x) = k [1 + ((x - λ)/a)^2]^(-m) exp(-ν atan((x - λ)/a)),
//! r = 6 (β - γ^2 - 1) / (2β - 3γ^2 - 6),   m = 1 + r/2,
//! D = 16 (r - 1) - γ^2 (r - 2)^2,
//! ν = -r (r - 2) γ / D^(1/2),   a = σ D^(1/2) / 4,   λ = μ - (r - 2) γ σ / 4,
//! k = |Γ(m + iν/2) / Γ(m)|^2 / (a B(m - 1/2, 1/2)).
//! ```
//!
//! **Validity region.** The inversion needs `2β - 3γ^2 - 6 > 0` (then
//! `r > 3`, i.e. `m > 5/2`, so the fourth moment exists) and `D > 0` (real
//! `ν`). Together these are Pearson's type IV criterion `0 < κ < 1`; the
//! boundary `D = 0` (`κ = 1`, the type IV / VI border) is
//! `β_IV(γ) = (48 + 39γ^2 + 6 (γ^2 + 4)^(3/2)) / (32 - γ^2)` for `γ^2 < 32`,
//! which follows from solving `D = 0` for `β` (a test checks the
//! equivalence). The Gaussian point `γ = 0, β = 3` lies on the boundary, so a
//! nearly Gaussian sample can fall outside by noise; that is reported as
//! [`PearsonError::OutsideTypeIv`], never as NaN.
//!
//! The normalisation uses the real part of the complex log-gamma function by
//! the Stirling series (M. Abramowitz and I. A. Stegun, *Handbook of
//! Mathematical Functions* (1964), 6.1.40) after upward recurrence.
//!
//! # Dual Pearson
//!
//! A. F. Tasch, H. Shin, C. Park, J. Alvis and S. Novak, "An improved
//! approach to accurately model shallow B and BF2 implants in silicon",
//! J. Electrochem. Soc. 136 (1989) 810, represent a profile as
//! `R f_1(x) + (1 - R) f_2(x)`: two Pearson IV densities with their own four
//! moments and a fraction `R`. [`DualPearson::fit`] fits the nine parameters
//! to a depth histogram; the method is described there.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use super::hist::Histogram;

/// Errors from building a Pearson representation.
#[derive(Debug, Clone, PartialEq, Error, Serialize, Deserialize)]
pub enum PearsonError {
    /// A moment is not finite, or the standard deviation is not positive.
    #[error("invalid moments: mean = {mean}, std_dev = {std_dev}, skewness = {skewness}, kurtosis = {kurtosis}")]
    InvalidMoments {
        /// Mean.
        mean: f64,
        /// Standard deviation.
        std_dev: f64,
        /// Skewness.
        skewness: f64,
        /// Kurtosis.
        kurtosis: f64,
    },
    /// `(γ, β)` is outside the Pearson type IV region.
    #[error("moments outside the Pearson IV region: skewness = {skewness}, kurtosis = {kurtosis} (needs kurtosis > {min_kurtosis:?})")]
    OutsideTypeIv {
        /// Skewness `γ`.
        skewness: f64,
        /// Kurtosis `β`.
        kurtosis: f64,
        /// The boundary `β_IV(γ)`; `None` when `γ^2 >= 32` (no type IV
        /// distribution has that skewness).
        min_kurtosis: Option<f64>,
    },
    /// The histogram has too few counts or bins to fit.
    #[error("histogram unsuitable for a fit: {0}")]
    BadHistogram(String),
}

/// The type IV / VI boundary `β_IV(γ)`; `None` for `γ^2 >= 32`.
pub fn type_iv_min_kurtosis(skewness: f64) -> Option<f64> {
    let g2 = skewness * skewness;
    (g2 < 32.0).then(|| (48.0 + 39.0 * g2 + 6.0 * (g2 + 4.0).powf(1.5)) / (32.0 - g2))
}

/// A Pearson type IV density, stored both as the four moments it matches and
/// as Heinrich's parameters.
///
/// Units: `mean`, `std_dev`, `a` and `lambda` are in the unit of the samples
/// (metres for depth profiles); the rest are dimensionless. `log_norm` is
/// `ln k` with `k` in inverse sample units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PearsonIv {
    /// Mean (`Rp` for a depth profile).
    pub mean: f64,
    /// Standard deviation (`ΔRp`).
    pub std_dev: f64,
    /// Skewness `γ`.
    pub skewness: f64,
    /// Kurtosis `β` (not excess).
    pub kurtosis: f64,
    /// Shape parameter `m > 5/2`.
    pub m: f64,
    /// Asymmetry parameter `ν`.
    pub nu: f64,
    /// Scale `a > 0`.
    pub a: f64,
    /// Location `λ`.
    pub lambda: f64,
    /// `ln k`, the log of the normalisation constant.
    pub log_norm: f64,
}

/// Real part of `ln Γ(x + iy)` for `x > 0` (Stirling series after recurrence;
/// Abramowitz and Stegun 6.1.40).
fn re_ln_gamma(x: f64, y: f64) -> f64 {
    let shift = if x < 12.0 {
        (12.0 - x).ceil() as usize
    } else {
        0
    };
    let mut correction = 0.0;
    for k in 0..shift {
        let xr = x + k as f64;
        correction += 0.5 * (xr * xr + y * y).ln();
    }
    let (wr, wi) = (x + shift as f64, y);
    let modulus2 = wr * wr + wi * wi;
    let ln_mod = 0.5 * modulus2.ln();
    let arg = wi.atan2(wr);
    // Re[(w - 1/2) ln w - w + ln(2 pi)/2]
    let mut s = (wr - 0.5) * ln_mod - wi * arg - wr + 0.5 * (2.0 * std::f64::consts::PI).ln();
    // Re of 1/(12 w) - 1/(360 w^3) + 1/(1260 w^5) - 1/(1680 w^7).
    let (ir, ii) = (wr / modulus2, -wi / modulus2); // 1/w
    let (i2r, i2i) = (ir * ir - ii * ii, 2.0 * ir * ii); // 1/w^2
    let mut pr = ir;
    let mut pi = ii;
    for c in [1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0, -1.0 / 1680.0] {
        s += c * pr;
        let t = pr * i2r - pi * i2i;
        pi = pr * i2i + pi * i2r;
        pr = t;
    }
    s - correction
}

impl PearsonIv {
    /// The Pearson IV density with the given mean, standard deviation,
    /// skewness and kurtosis (Heinrich 2004). Errors if a moment is not
    /// finite, `std_dev <= 0`, or `(γ, β)` is outside the type IV region.
    pub fn from_moments(
        mean: f64,
        std_dev: f64,
        skewness: f64,
        kurtosis: f64,
    ) -> Result<Self, PearsonError> {
        if !(mean.is_finite()
            && std_dev.is_finite()
            && std_dev > 0.0
            && skewness.is_finite()
            && kurtosis.is_finite())
        {
            return Err(PearsonError::InvalidMoments {
                mean,
                std_dev,
                skewness,
                kurtosis,
            });
        }
        let (g, b) = (skewness, kurtosis);
        let outside = || PearsonError::OutsideTypeIv {
            skewness: g,
            kurtosis: b,
            min_kurtosis: type_iv_min_kurtosis(g),
        };
        let den = 2.0 * b - 3.0 * g * g - 6.0;
        if den.is_nan() || den <= 0.0 {
            return Err(outside());
        }
        let r = 6.0 * (b - g * g - 1.0) / den;
        let d = 16.0 * (r - 1.0) - g * g * (r - 2.0) * (r - 2.0);
        if !(d > 0.0 && r > 3.0) {
            return Err(outside());
        }
        let m = 1.0 + 0.5 * r;
        let sd = d.sqrt();
        let nu = -r * (r - 2.0) * g / sd;
        let a = std_dev * sd / 4.0;
        let lambda = mean - (r - 2.0) * g * std_dev / 4.0;
        let ln_beta = re_ln_gamma(m - 0.5, 0.0) + re_ln_gamma(0.5, 0.0) - re_ln_gamma(m, 0.0);
        let log_norm = 2.0 * (re_ln_gamma(m, 0.5 * nu) - re_ln_gamma(m, 0.0)) - a.ln() - ln_beta;
        Ok(Self {
            mean,
            std_dev,
            skewness,
            kurtosis,
            m,
            nu,
            a,
            lambda,
            log_norm,
        })
    }

    /// Probability density at `x`, in inverse sample units.
    pub fn pdf(&self, x: f64) -> f64 {
        let t = (x - self.lambda) / self.a;
        (self.log_norm - self.m * t.mul_add(t, 1.0).ln() - self.nu * t.atan()).exp()
    }
}

/// Probability mass of each bin of `hist`'s binning under `pdf`, by
/// three-point Gauss-Legendre quadrature per bin.
fn bin_masses(hist: &Histogram, pdf: impl Fn(f64) -> f64) -> Vec<f64> {
    let b = &hist.binning;
    let h = b.width();
    let off = 0.5 * h * (0.6f64).sqrt();
    (0..b.bins)
        .map(|i| {
            let c = b.center(i);
            0.5 * h * ((5.0 / 9.0) * (pdf(c - off) + pdf(c + off)) + (8.0 / 9.0) * pdf(c))
        })
        .collect()
}

/// Neyman chi-square `Σ (c_i - N P_i)^2 / max(c_i, 1)` of a model with bin
/// masses `p` against the counts, `N` the histogram total.
fn chi_square(hist: &Histogram, p: &[f64]) -> f64 {
    let n = hist.total() as f64;
    let s: f64 = hist
        .counts
        .iter()
        .zip(p)
        .map(|(&c, &pi)| {
            let d = c as f64 - n * pi;
            d * d / (c as f64).max(1.0)
        })
        .sum();
    if s.is_finite() {
        s
    } else {
        f64::INFINITY
    }
}

/// A dual-Pearson profile `R f_head + (1 - R) f_tail` (Tasch et al. 1989).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DualPearson {
    /// Fraction `R` in the head component, in `[0, 1]`.
    pub head_fraction: f64,
    /// Head component.
    pub head: PearsonIv,
    /// Tail component.
    pub tail: PearsonIv,
}

/// The result of [`DualPearson::fit`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DualPearsonFit {
    /// The fitted profile.
    pub profile: DualPearson,
    /// Neyman chi-square of the fit over the in-range bins.
    pub chi_square: f64,
    /// The same statistic for the single Pearson IV from the histogram's
    /// moments, if those are in the type IV region.
    pub single_chi_square: Option<f64>,
    /// Number of in-range bins (data points of the fit).
    pub bins: usize,
    /// Objective evaluations used.
    pub evaluations: u64,
}

/// Map an unconstrained 4-vector to a component in the type IV region:
/// `(mean, ln σ, atanh(γ/4), ln(β - β_IV(γ)))`.
fn component(p: &[f64]) -> Option<PearsonIv> {
    let g = 4.0 * p[2].tanh();
    let b = type_iv_min_kurtosis(g)? + p[3].exp().max(1e-6);
    PearsonIv::from_moments(p[0], p[1].exp(), g, b).ok()
}

fn unconstrained(mean: f64, sd: f64, g: f64, b: f64) -> [f64; 4] {
    let g = g.clamp(-3.6, 3.6);
    let excess = (b - type_iv_min_kurtosis(g).expect("|g| < 4")).max(0.2);
    [mean, sd.ln(), (g / 4.0).atanh(), excess.ln()]
}

fn logistic(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Minimise `f` from `x0` with the Nelder-Mead simplex method (J. A. Nelder
/// and R. Mead, Comput. J. 7 (1965) 308) with the standard coefficients
/// (reflection 1, expansion 2, contraction 1/2, shrink 1/2) as stated by
/// J. C. Lagarias, J. A. Reeds, M. H. Wright and P. E. Wright, SIAM J. Optim.
/// 9 (1998) 112. Deterministic: fixed initial simplex, ties broken by vertex
/// index, fixed evaluation cap.
fn nelder_mead(
    f: &dyn Fn(&[f64]) -> f64,
    x0: &[f64],
    steps: &[f64],
    max_evals: u64,
    evals: &mut u64,
) -> (Vec<f64>, f64) {
    let n = x0.len();
    let used = std::cell::Cell::new(0u64);
    let eval = |x: &[f64]| {
        used.set(used.get() + 1);
        let v = f(x);
        if v.is_nan() {
            f64::INFINITY
        } else {
            v
        }
    };
    let mut pts: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    let v0 = eval(x0);
    pts.push((x0.to_vec(), v0));
    for i in 0..n {
        let mut x = x0.to_vec();
        x[i] += steps[i];
        let v = eval(&x);
        pts.push((x, v));
    }
    while used.get() < max_evals {
        // Stable sort keeps vertex order on ties.
        pts.sort_by(|a, b| a.1.total_cmp(&b.1));
        let (best, worst) = (pts[0].1, pts[n].1);
        if (worst - best).abs() <= 1e-12 * (best.abs() + 1e-12) {
            break;
        }
        let mut c = vec![0.0; n];
        for (x, _) in &pts[..n] {
            for (ci, xi) in c.iter_mut().zip(x) {
                *ci += xi / n as f64;
            }
        }
        let along = |t: f64| -> Vec<f64> {
            c.iter()
                .zip(&pts[n].0)
                .map(|(ci, wi)| ci + t * (ci - wi))
                .collect()
        };
        let xr = along(1.0);
        let fr = eval(&xr);
        if fr < pts[0].1 {
            let xe = along(2.0);
            let fe = eval(&xe);
            pts[n] = if fe < fr { (xe, fe) } else { (xr, fr) };
        } else if fr < pts[n - 1].1 {
            pts[n] = (xr, fr);
        } else {
            let (xc, fc) = if fr < pts[n].1 {
                let x = along(0.5);
                let v = eval(&x);
                (x, v)
            } else {
                let x = along(-0.5);
                let v = eval(&x);
                (x, v)
            };
            if fc < fr.min(pts[n].1) {
                pts[n] = (xc, fc);
            } else {
                let x_best = pts[0].0.clone();
                for (x, v) in pts.iter_mut().skip(1) {
                    for (xi, bi) in x.iter_mut().zip(&x_best) {
                        *xi = bi + 0.5 * (*xi - bi);
                    }
                    *v = eval(x);
                }
            }
        }
    }
    pts.sort_by(|a, b| a.1.total_cmp(&b.1));
    *evals += used.get();
    pts.swap_remove(0)
}

impl DualPearson {
    /// Probability density at `x`.
    pub fn pdf(&self, x: f64) -> f64 {
        self.head_fraction * self.head.pdf(x) + (1.0 - self.head_fraction) * self.tail.pdf(x)
    }

    /// Mean, standard deviation, skewness and kurtosis of the mixture, from
    /// the components' moments.
    pub fn moments(&self) -> (f64, f64, f64, f64) {
        let w = [self.head_fraction, 1.0 - self.head_fraction];
        let c = [&self.head, &self.tail];
        let mean = w[0] * c[0].mean + w[1] * c[1].mean;
        // Central moments of the mixture about `mean`: for each component,
        // E[(x - M)^p] = Σ_j C(p, j) μ_j d^(p-j) with d = mean_k - M.
        const BINOM: [[f64; 5]; 5] = [
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 3.0, 1.0, 0.0],
            [1.0, 4.0, 6.0, 4.0, 1.0],
        ];
        let mut mu = [0.0f64; 5];
        for (wk, ck) in w.iter().zip(c) {
            let s = ck.std_dev;
            let own = [
                1.0,
                0.0,
                s * s,
                ck.skewness * s.powi(3),
                ck.kurtosis * s.powi(4),
            ];
            let d = ck.mean - mean;
            for (p, slot) in mu.iter_mut().enumerate() {
                let acc: f64 = (0..=p)
                    .map(|j| BINOM[p][j] * own[j] * d.powi((p - j) as i32))
                    .sum();
                *slot += wk * acc;
            }
        }
        let sd = mu[2].sqrt();
        (mean, sd, mu[3] / (sd * sd * sd), mu[4] / (mu[2] * mu[2]))
    }

    /// Fit a dual-Pearson profile to a histogram (for example the depth
    /// histogram of stopped ions).
    ///
    /// **Method.** The bins are rescaled to `z = (x - μ) / σ` with the mean
    /// and standard deviation of the in-range bin centres. The nine
    /// parameters are mapped to unconstrained variables (the logit of `R`;
    /// for each component the mean, `ln σ`, `atanh(γ/4)` and
    /// `ln(β - β_IV(γ))`), so every trial point is a valid pair of Pearson IV
    /// densities with `|γ| < 4`. The objective is the Neyman chi-square
    /// `Σ (c_i - N P_i)^2 / max(c_i, 1)` over the in-range bins, with `P_i`
    /// the model mass of bin `i` by three-point Gauss-Legendre quadrature and
    /// `N` the histogram total (so mass outside the range counts against the
    /// model). It is minimised by Nelder-Mead from two fixed starting points
    /// (a head at the overall moments with `R = 0.9`, and a tail one standard
    /// deviation deeper or shallower), each restarted once from its best
    /// point, with a cap of 6000 evaluations per run. The lowest objective
    /// wins; if the single Pearson IV from the histogram's moments is better,
    /// it is returned as a degenerate dual (`R = 1`, both components equal).
    ///
    /// The result is a deterministic function of the histogram (no random
    /// numbers, fixed iteration order). It is a least-squares representation,
    /// not a unique decomposition: different parameter sets can describe the
    /// same profile almost equally well, and the split into "head" and
    /// "tail" is only physically meaningful when the histogram resolves it.
    pub fn fit(hist: &Histogram) -> Result<DualPearsonFit, PearsonError> {
        let in_range: u64 = hist.counts.iter().sum();
        let nonzero = hist.counts.iter().filter(|&&c| c > 0).count();
        if in_range < 100 || nonzero < 10 {
            return Err(PearsonError::BadHistogram(format!(
                "{in_range} in-range counts in {nonzero} non-empty bins (need >= 100 and >= 10)"
            )));
        }
        // Moments of the in-range bin centres.
        let b = &hist.binning;
        let mut mom = super::moments::Moments::new();
        for (i, &c) in hist.counts.iter().enumerate() {
            let mut one = super::moments::Moments::new();
            if c > 0 {
                // c identical samples at the bin centre.
                one.n = c;
                one.mean = b.center(i);
                mom.merge(&one);
            }
        }
        let s = mom
            .summary()
            .ok_or_else(|| PearsonError::BadHistogram("zero variance".into()))?;
        let (mu, sigma) = (s.mean, s.std_dev);

        // The objective in z units.
        let objective = |p: &[f64]| -> f64 {
            let (Some(h), Some(t)) = (component(&p[1..5]), component(&p[5..9])) else {
                return f64::INFINITY;
            };
            let r = logistic(p[0]);
            let masses = bin_masses(hist, |x| {
                let z = (x - mu) / sigma;
                (r * h.pdf(z) + (1.0 - r) * t.pdf(z)) / sigma
            });
            chi_square(hist, &masses)
        };

        let single = PearsonIv::from_moments(0.0, 1.0, s.skewness, s.kurtosis).ok();
        let single_chi = single
            .map(|p| chi_square(hist, &bin_masses(hist, |x| p.pdf((x - mu) / sigma) / sigma)));

        let head0 = unconstrained(0.0, 1.0, s.skewness, s.kurtosis);
        let steps = [0.5, 0.25, 0.2, 0.2, 0.5, 0.25, 0.2, 0.2, 0.5];
        let mut evals = 0;
        let mut best: Option<(Vec<f64>, f64)> = None;
        for tail_shift in [1.0, -1.0] {
            let tail0 = unconstrained(tail_shift, 1.5, 0.0, 4.0);
            let mut x0 = vec![(0.9f64 / 0.1).ln()];
            x0.extend_from_slice(&head0);
            x0.extend_from_slice(&tail0);
            let (x1, _) = nelder_mead(&objective, &x0, &steps, 6000, &mut evals);
            let (x2, v2) = nelder_mead(&objective, &x1, &steps, 6000, &mut evals);
            if best.as_ref().is_none_or(|(_, v)| v2 < *v) {
                best = Some((x2, v2));
            }
        }
        let (x, chi) = best.expect("two starts");
        let to_x = |c: PearsonIv| {
            PearsonIv::from_moments(
                mu + sigma * c.mean,
                sigma * c.std_dev,
                c.skewness,
                c.kurtosis,
            )
        };
        let fit = match (single, single_chi) {
            (Some(p), Some(sc)) if sc <= chi => {
                let p = to_x(p)?;
                DualPearsonFit {
                    profile: DualPearson {
                        head_fraction: 1.0,
                        head: p,
                        tail: p,
                    },
                    chi_square: sc,
                    single_chi_square: Some(sc),
                    bins: b.bins,
                    evaluations: evals,
                }
            }
            _ => {
                let h = component(&x[1..5]).expect("finite objective");
                let t = component(&x[5..9]).expect("finite objective");
                DualPearsonFit {
                    profile: DualPearson {
                        head_fraction: logistic(x[0]),
                        head: to_x(h)?,
                        tail: to_x(t)?,
                    },
                    chi_square: chi,
                    single_chi_square: single_chi,
                    bins: b.bins,
                    evaluations: evals,
                }
            }
        };
        Ok(fit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_gamma_matches_known_values() {
        let half = 0.5 * std::f64::consts::PI.ln();
        assert!((re_ln_gamma(0.5, 0.0) - half).abs() < 1e-12);
        assert!((re_ln_gamma(10.0, 0.0) - 362_880f64.ln()).abs() < 1e-12);
        assert!((re_ln_gamma(1.0, 0.0)).abs() < 1e-12);
        // |Γ(1 + iy)|^2 = π y / sinh(π y).
        for y in [0.3, 2.0, 15.0] {
            let pi_y = std::f64::consts::PI * y;
            let want = 0.5 * (pi_y / pi_y.sinh()).ln();
            assert!((re_ln_gamma(1.0, y) - want).abs() < 1e-12, "y = {y}");
        }
    }

    #[test]
    fn boundary_formula_is_d_equals_zero() {
        // gamma = 0: the boundary is the Gaussian point beta = 3.
        assert_eq!(type_iv_min_kurtosis(0.0), Some(3.0));
        assert!(PearsonIv::from_moments(0.0, 1.0, 0.0, 3.0).is_err());
        assert!(PearsonIv::from_moments(0.0, 1.0, 0.0, 3.001).is_ok());
        for g in [0.3, -1.0, 2.0, -4.5] {
            let b = type_iv_min_kurtosis(g).unwrap();
            let r = 6.0 * (b - g * g - 1.0) / (2.0 * b - 3.0 * g * g - 6.0);
            let d = 16.0 * (r - 1.0) - g * g * (r - 2.0) * (r - 2.0);
            assert!(d.abs() < 1e-9 * r * r, "g = {g}: D = {d}");
            assert!(PearsonIv::from_moments(0.0, 1.0, g, b * 0.999).is_err());
            assert!(PearsonIv::from_moments(0.0, 1.0, g, b * 1.001 + 1e-6).is_ok());
        }
        assert_eq!(type_iv_min_kurtosis(6.0), None);
    }

    #[test]
    fn symmetric_case_is_student_like() {
        // gamma = 0: nu = 0, lambda = mean.
        let p = PearsonIv::from_moments(2.0, 1.0, 0.0, 4.0).unwrap();
        assert_eq!(p.nu, 0.0);
        assert_eq!(p.lambda, 2.0);
        assert!((p.pdf(2.5) - p.pdf(1.5)).abs() < 1e-15);
    }
}
