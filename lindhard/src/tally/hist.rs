//! Fixed-binning histograms with integer counts.
//!
//! Counts are `u64`, so merging two histograms is an exact, associative and
//! commutative sum: the merged result does not depend on how the samples were
//! split between workers.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors from building a [`Binning`].
#[derive(Debug, Clone, PartialEq, Error)]
pub enum BinningError {
    /// The edges are not finite with `hi > lo`, or there are no bins.
    #[error(
        "invalid binning: lo = {lo}, hi = {hi}, bins = {bins} (need finite lo < hi and bins > 0)"
    )]
    Invalid {
        /// Lower edge.
        lo: f64,
        /// Upper edge.
        hi: f64,
        /// Number of bins.
        bins: usize,
    },
}

/// Uniform binning of `[lo, hi)` into `bins` equal bins.
///
/// The unit of `lo` and `hi` is the unit of the quantity binned (metres for
/// depths, eV for energies, radians for angles); each use documents it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Binning {
    /// Lower edge of the first bin.
    pub lo: f64,
    /// Upper edge of the last bin.
    pub hi: f64,
    /// Number of bins.
    pub bins: usize,
}

impl Binning {
    /// Validated binning of `[lo, hi)` into `bins` bins.
    pub fn new(lo: f64, hi: f64, bins: usize) -> Result<Self, BinningError> {
        if lo.is_finite() && hi.is_finite() && hi > lo && bins > 0 {
            Ok(Self { lo, hi, bins })
        } else {
            Err(BinningError::Invalid { lo, hi, bins })
        }
    }

    /// Width of one bin.
    pub fn width(&self) -> f64 {
        (self.hi - self.lo) / self.bins as f64
    }

    /// Lower edge of bin `i`.
    pub fn edge(&self, i: usize) -> f64 {
        self.lo + (self.hi - self.lo) * (i as f64 / self.bins as f64)
    }

    /// Centre of bin `i`.
    pub fn center(&self, i: usize) -> f64 {
        0.5 * (self.edge(i) + self.edge(i + 1))
    }

    /// Where `x` falls: `Ok(bin)`, or `Err(Below)` / `Err(Above)`. NaN counts
    /// as above.
    pub fn locate(&self, x: f64) -> Result<usize, OutOfRange> {
        if x < self.lo {
            return Err(OutOfRange::Below);
        }
        if x.is_nan() || x >= self.hi {
            return Err(OutOfRange::Above);
        }
        let i = ((x - self.lo) / (self.hi - self.lo) * self.bins as f64) as usize;
        Ok(i.min(self.bins - 1))
    }
}

/// Which side of a [`Binning`] a value fell on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutOfRange {
    /// `x < lo`.
    Below,
    /// `x >= hi`, or NaN.
    Above,
}

/// Counts over a [`Binning`], with underflow and overflow.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Histogram {
    /// The binning.
    pub binning: Binning,
    /// Count per bin.
    pub counts: Vec<u64>,
    /// Samples below `binning.lo`.
    pub underflow: u64,
    /// Samples at or above `binning.hi` (and NaN samples).
    pub overflow: u64,
}

impl Histogram {
    /// Empty histogram.
    pub fn new(binning: Binning) -> Self {
        Self {
            binning,
            counts: vec![0; binning.bins],
            underflow: 0,
            overflow: 0,
        }
    }

    /// Count one sample.
    pub fn fill(&mut self, x: f64) {
        match self.binning.locate(x) {
            Ok(i) => self.counts[i] += 1,
            Err(OutOfRange::Below) => self.underflow += 1,
            Err(OutOfRange::Above) => self.overflow += 1,
        }
    }

    /// All samples, in range or not.
    pub fn total(&self) -> u64 {
        self.counts.iter().sum::<u64>() + self.underflow + self.overflow
    }

    /// Add `other`, which must have the same binning.
    ///
    /// # Panics
    /// If the binnings differ (merging tallies built with different settings
    /// is a programming error).
    pub fn merge(&mut self, other: &Histogram) {
        assert_eq!(
            self.binning, other.binning,
            "merging histograms with different binnings"
        );
        for (a, b) in self.counts.iter_mut().zip(&other.counts) {
            *a += b;
        }
        self.underflow += other.underflow;
        self.overflow += other.overflow;
    }

    /// Probability density per bin: `count / (total * width)`, so that the
    /// in-range densities integrate to the in-range fraction. Unit: inverse
    /// of the binned quantity's unit. All zeros if the histogram is empty.
    pub fn density(&self) -> Vec<f64> {
        let n = self.total();
        if n == 0 {
            return vec![0.0; self.counts.len()];
        }
        let norm = 1.0 / (n as f64 * self.binning.width());
        self.counts.iter().map(|&c| c as f64 * norm).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binning_rejects_bad_input() {
        assert!(Binning::new(0.0, 0.0, 1).is_err());
        assert!(Binning::new(1.0, 0.0, 1).is_err());
        assert!(Binning::new(0.0, 1.0, 0).is_err());
        assert!(Binning::new(f64::NAN, 1.0, 3).is_err());
        assert!(Binning::new(0.0, f64::INFINITY, 3).is_err());
    }

    #[test]
    fn fill_and_edges() {
        let mut h = Histogram::new(Binning::new(0.0, 10.0, 5).unwrap());
        for x in [-1.0, 0.0, 1.99, 2.0, 9.999, 10.0, f64::NAN] {
            h.fill(x);
        }
        assert_eq!(h.counts, vec![2, 1, 0, 0, 1]);
        assert_eq!((h.underflow, h.overflow, h.total()), (1, 2, 7));
        assert_eq!(h.binning.center(0), 1.0);
        let mut g = h.clone();
        g.merge(&h);
        assert_eq!(g.counts, vec![4, 2, 0, 0, 2]);
        assert_eq!(g.total(), 14);
    }
}
