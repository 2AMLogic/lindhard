//! One-pass, mergeable central moments with standard errors.
//!
//! # Accumulator
//!
//! [`Moments`] keeps the count `n`, the mean, and the centred power sums
//! `M_p = Σ (x_i - mean)^p` for `p = 2..=8`. Two accumulators `A`, `B` combine
//! with the arbitrary-order pairwise update of P. Pébay, "Formulas for robust,
//! one-pass parallel computation of covariances and arbitrary-order
//! statistical moments", Sandia report SAND2008-6212 (2008), which generalises
//! the second-order update of T. F. Chan, G. H. Golub and R. J. LeVeque,
//! "Updating formulae and a pairwise algorithm for computing sample
//! variances", Stanford report STAN-CS-79-773 (1979). With `δ = mean_B -
//! mean_A` and `n = n_A + n_B`,
//!
//! ```text
//! M_p = M_p^A + M_p^B
//!     + Σ_{k=1}^{p-2} C(p,k) δ^k [ (-n_B/n)^k M_{p-k}^A + (n_A/n)^k M_{p-k}^B ]
//!     + (n_A n_B δ / n)^p [ n_B^{1-p} - (-n_A)^{1-p} ].
//! ```
//!
//! Adding one sample is a merge with `n_B = 1`, `M^B = 0`. The update is exact
//! in real arithmetic; in floating point the result depends on the order of
//! merges, so it is **deterministic for a fixed order** (bit-identical) and
//! agrees across different splittings to rounding. [`crate::rng::run_particles`]
//! merges in chunk order, which fixes the order independently of the thread
//! count.
//!
//! # Conventions
//!
//! Central moments are the method-of-moments estimates `μ_p = M_p / n` (no
//! `n - 1` correction), as used for fitting Pearson distributions:
//! standard deviation `σ = μ_2^(1/2)`, skewness `γ = μ_3 / μ_2^(3/2)`,
//! kurtosis `β = μ_4 / μ_2^2` (**not** excess kurtosis: a Gaussian has
//! `β = 3`). For a depth distribution these are `Rp`, `ΔRp`, `γ`, `β`.
//!
//! # Standard errors
//!
//! Large-sample standard errors from the asymptotic covariance of sample
//! central moments, H. Cramér, *Mathematical Methods of Statistics*
//! (Princeton, 1946), §27.7:
//!
//! ```text
//! n Cov(m_j, m_k) = μ_{j+k} - μ_j μ_k - j μ_{j-1} μ_{k+1} - k μ_{j+1} μ_{k-1}
//!                 + j k μ_{j-1} μ_{k-1} μ_2,
//! ```
//!
//! propagated to `σ`, `γ` and `β` by the delta method (first-order Taylor
//! expansion), with the population moments replaced by their estimates up to
//! `μ_8`. The standard error of the mean is `σ / n^(1/2)`. No normality is
//! assumed; for a Gaussian these reduce to the familiar `(6/n)^(1/2)` for `γ`
//! and `(24/n)^(1/2)` for `β` (checked in a test). They are asymptotic, so
//! treat them as approximate for small `n` or heavy tails.

use serde::{Deserialize, Serialize};

/// Highest power sum kept.
const MAX_ORDER: usize = 8;

/// Mergeable accumulator of the count, mean and centred power sums up to
/// order 8. See the module docs for the update formula.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Moments {
    /// Number of samples.
    pub n: u64,
    /// Running mean.
    pub mean: f64,
    /// Centred power sums `M_2 .. M_8` (index `p - 2`).
    pub power_sums: [f64; MAX_ORDER - 1],
}

fn binom(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

impl Moments {
    /// Empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// `M_p` for `p` in `0..=8` (`M_0 = n`, `M_1 = 0`).
    fn m(&self, p: usize) -> f64 {
        match p {
            0 => self.n as f64,
            1 => 0.0,
            _ => self.power_sums[p - 2],
        }
    }

    /// Add one sample.
    pub fn push(&mut self, x: f64) {
        self.merge_parts(1, x, &[0.0; MAX_ORDER - 1]);
    }

    /// Fold `other` into `self` (Pébay 2008).
    pub fn merge(&mut self, other: &Moments) {
        self.merge_parts(other.n, other.mean, &other.power_sums);
    }

    fn merge_parts(&mut self, n_b: u64, mean_b: f64, sums_b: &[f64; MAX_ORDER - 1]) {
        if n_b == 0 {
            return;
        }
        if self.n == 0 {
            self.n = n_b;
            self.mean = mean_b;
            self.power_sums = *sums_b;
            return;
        }
        let (na, nb) = (self.n as f64, n_b as f64);
        let n = na + nb;
        let delta = mean_b - self.mean;
        let m_b = |p: usize| match p {
            0 => nb,
            1 => 0.0,
            _ => sums_b[p - 2],
        };
        let mut out = [0.0; MAX_ORDER - 1];
        for (idx, slot) in out.iter_mut().enumerate() {
            let p = idx + 2;
            let mut s = self.m(p) + m_b(p);
            for k in 1..=p - 2 {
                s += binom(p, k)
                    * delta.powi(k as i32)
                    * ((-nb / n).powi(k as i32) * self.m(p - k)
                        + (na / n).powi(k as i32) * m_b(p - k));
            }
            let pi = p as i32;
            s += (na * nb * delta / n).powi(pi) * (nb.powi(1 - pi) - (-na).powi(1 - pi));
            *slot = s;
        }
        self.power_sums = out;
        self.mean += delta * nb / n;
        self.n += n_b;
    }

    /// Central moment estimate `μ_p = M_p / n` (`μ_0 = 1`, `μ_1 = 0`).
    pub fn central(&self, p: usize) -> f64 {
        assert!(p <= MAX_ORDER, "central moments are kept up to order 8");
        match p {
            0 => 1.0,
            1 => 0.0,
            _ => self.power_sums[p - 2] / self.n as f64,
        }
    }

    /// Mean, standard deviation, skewness, kurtosis and their standard
    /// errors. `None` if fewer than two samples or zero variance.
    pub fn summary(&self) -> Option<MomentSummary> {
        if self.n < 2 {
            return None;
        }
        let mu = |p: usize| self.central(p);
        let (m2, m3, m4) = (mu(2), mu(3), mu(4));
        if m2.is_nan() || m2 <= 0.0 {
            return None;
        }
        let n = self.n as f64;
        // n Cov(m_j, m_k), Cramér (1946) §27.7.
        let cov = |j: usize, k: usize| {
            let (jf, kf) = (j as f64, k as f64);
            mu(j + k) - mu(j) * mu(k) - jf * mu(j - 1) * mu(k + 1) - kf * mu(j + 1) * mu(k - 1)
                + jf * kf * mu(j - 1) * mu(k - 1) * m2
        };
        let quad = |g: [(usize, f64); 2]| {
            let mut v = 0.0;
            for &(j, gj) in &g {
                for &(k, gk) in &g {
                    v += gj * gk * cov(j, k);
                }
            }
            (v.max(0.0) / n).sqrt()
        };
        let sd = m2.sqrt();
        let skew = m3 / (m2 * sd);
        let kurt = m4 / (m2 * m2);
        Some(MomentSummary {
            n: self.n,
            mean: self.mean,
            std_dev: sd,
            skewness: skew,
            kurtosis: kurt,
            mean_std_err: (m2 / n).sqrt(),
            std_dev_std_err: (cov(2, 2).max(0.0) / n).sqrt() / (2.0 * sd),
            skewness_std_err: quad([(2, -1.5 * m3 / (m2 * m2 * sd)), (3, 1.0 / (m2 * sd))]),
            kurtosis_std_err: quad([(2, -2.0 * m4 / (m2 * m2 * m2)), (4, 1.0 / (m2 * m2))]),
        })
    }
}

/// The first four moments of a sample, with standard errors.
///
/// Units: `mean`, `std_dev` and their errors are in the unit of the samples
/// (metres for depths and lateral positions); `skewness` and `kurtosis` are
/// dimensionless. For a depth distribution `mean` is the projected range
/// `Rp`, `std_dev` the straggle `ΔRp`, `skewness` `γ` and `kurtosis` `β`
/// (Gaussian `β = 3`). See the module docs for conventions.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MomentSummary {
    /// Number of samples.
    pub n: u64,
    /// Mean.
    pub mean: f64,
    /// Standard deviation `μ_2^(1/2)`.
    pub std_dev: f64,
    /// Skewness `γ = μ_3 / μ_2^(3/2)`.
    pub skewness: f64,
    /// Kurtosis `β = μ_4 / μ_2^2` (not excess).
    pub kurtosis: f64,
    /// Standard error of `mean`.
    pub mean_std_err: f64,
    /// Standard error of `std_dev`.
    pub std_dev_std_err: f64,
    /// Standard error of `skewness`.
    pub skewness_std_err: f64,
    /// Standard error of `kurtosis`.
    pub kurtosis_std_err: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_pass(xs: &[f64]) -> (f64, Vec<f64>) {
        let n = xs.len() as f64;
        let mean = xs.iter().sum::<f64>() / n;
        let sums = (2..=8)
            .map(|p| xs.iter().map(|x| (x - mean).powi(p)).sum::<f64>())
            .collect();
        (mean, sums)
    }

    #[test]
    fn push_matches_two_pass() {
        // A skewed deterministic sample.
        let xs: Vec<f64> = (0..500)
            .map(|i| ((i as f64) * 0.37).sin().exp() * 3.0 + 1.0)
            .collect();
        let mut m = Moments::new();
        xs.iter().for_each(|&x| m.push(x));
        let (mean, sums) = two_pass(&xs);
        assert!((m.mean - mean).abs() < 1e-13 * mean.abs());
        for (a, b) in m.power_sums.iter().zip(&sums) {
            assert!((a - b).abs() <= 1e-10 * b.abs().max(1e-12), "{a} vs {b}");
        }
    }

    #[test]
    fn merge_of_empty_is_identity() {
        let mut a = Moments::new();
        a.push(1.0);
        a.push(4.0);
        let before = a.clone();
        a.merge(&Moments::new());
        assert_eq!(a, before);
        let mut e = Moments::new();
        e.merge(&before);
        assert_eq!(e, before);
    }

    #[test]
    fn normal_theory_standard_errors() {
        // Exact Gaussian moments with sigma = 1: mu4 = 3, mu6 = 15, mu8 = 105.
        let m = Moments {
            n: 100,
            mean: 0.0,
            power_sums: [100.0, 0.0, 300.0, 0.0, 1500.0, 0.0, 10500.0],
        };
        let s = m.summary().unwrap();
        assert!((s.skewness_std_err - (6.0f64 / 100.0).sqrt()).abs() < 1e-14);
        assert!((s.kurtosis_std_err - (24.0f64 / 100.0).sqrt()).abs() < 1e-14);
        assert!((s.std_dev_std_err - (0.5f64 / 100.0).sqrt()).abs() < 1e-14);
        assert!((s.mean_std_err - 0.1).abs() < 1e-15);
    }

    #[test]
    fn degenerate_samples_have_no_summary() {
        let mut m = Moments::new();
        assert!(m.summary().is_none());
        m.push(2.0);
        assert!(m.summary().is_none());
        m.push(2.0);
        assert!(m.summary().is_none());
    }
}
