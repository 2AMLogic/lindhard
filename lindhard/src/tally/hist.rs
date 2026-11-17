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

/// Energy (or any other `f64` weight) over a [`Binning`], with underflow and
/// overflow: the weighted companion of [`Histogram`].
///
/// Sums are plain `f64` additions, merged bin by bin in the order
/// [`crate::rng::run_particles`] calls `merge` (chunk order), so for a fixed
/// chunk size the result is bit-identical at any thread count; unlike the
/// integer counts of [`Histogram`] it is not invariant to the chunk size.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WeightedHistogram {
    /// The binning.
    pub binning: Binning,
    /// Summed weight per bin.
    pub sums: Vec<f64>,
    /// Weight below `binning.lo`.
    pub underflow: f64,
    /// Weight at or above `binning.hi` (and at NaN positions).
    pub overflow: f64,
}

impl WeightedHistogram {
    /// Empty histogram.
    pub fn new(binning: Binning) -> Self {
        Self {
            binning,
            sums: vec![0.0; binning.bins],
            underflow: 0.0,
            overflow: 0.0,
        }
    }

    /// Add weight `w` at position `x`.
    pub fn fill(&mut self, x: f64, w: f64) {
        match self.binning.locate(x) {
            Ok(i) => self.sums[i] += w,
            Err(OutOfRange::Below) => self.underflow += w,
            Err(OutOfRange::Above) => self.overflow += w,
        }
    }

    /// Spread weight `w` uniformly along the segment from `x0` to `x1`: each
    /// bin (and the underflow and overflow) gets `w` times the fraction of
    /// the segment length inside it. A zero-length or non-finite segment is
    /// filled at `x1`.
    pub fn fill_segment(&mut self, x0: f64, x1: f64, w: f64) {
        let (a, b) = if x0 <= x1 { (x0, x1) } else { (x1, x0) };
        let len = b - a;
        if !(len > 0.0 && len.is_finite()) {
            self.fill(x1, w);
            return;
        }
        let (lo, hi) = (self.binning.lo, self.binning.hi);
        // Fraction of the segment at positions below `x`.
        let below = |x: f64| ((x.clamp(a, b) - a) / len).clamp(0.0, 1.0);
        self.underflow += w * below(lo);
        self.overflow += w * (1.0 - below(hi));
        if b < lo || a >= hi {
            return;
        }
        let first = self.binning.locate(a.max(lo)).unwrap_or(0);
        let last = self
            .binning
            .locate(b.min(hi))
            .unwrap_or(self.binning.bins - 1);
        for i in first..=last {
            let f = below(self.binning.edge(i + 1)) - below(self.binning.edge(i));
            self.sums[i] += w * f;
        }
    }

    /// Everything, in range or not.
    pub fn total(&self) -> f64 {
        self.sums.iter().sum::<f64>() + self.underflow + self.overflow
    }

    /// Add `other`, which must have the same binning.
    ///
    /// # Panics
    /// If the binnings differ.
    pub fn merge(&mut self, other: &WeightedHistogram) {
        assert_eq!(
            self.binning, other.binning,
            "merging histograms with different binnings"
        );
        for (a, b) in self.sums.iter_mut().zip(&other.sums) {
            *a += b;
        }
        self.underflow += other.underflow;
        self.overflow += other.overflow;
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

    #[test]
    fn weighted_segment_conserves_weight() {
        let mut h = WeightedHistogram::new(Binning::new(0.0, 10.0, 5).unwrap());
        h.fill(3.0, 1.0);
        h.fill(-1.0, 2.0);
        h.fill(f64::NAN, 4.0);
        // 12 long, weight 1 per unit length: 2 below lo, 2 per bin, none above.
        h.fill_segment(-2.0, 10.0, 12.0);
        assert!((h.underflow - 4.0).abs() < 1e-12);
        assert!((h.sums[0] - 2.0).abs() < 1e-12);
        assert!((h.sums[1] - 3.0).abs() < 1e-12);
        assert_eq!(h.overflow, 4.0);
        // Reversed, spanning both edges, crossing every bin.
        let mut g = WeightedHistogram::new(h.binning);
        g.fill_segment(15.0, -5.0, 20.0);
        assert!((g.underflow - 5.0).abs() < 1e-12);
        assert!((g.overflow - 5.0).abs() < 1e-12);
        assert!(g.sums.iter().all(|s| (s - 2.0).abs() < 1e-12));
        assert!((g.total() - 20.0).abs() < 1e-12);
        // Degenerate segment is a point fill.
        let mut p = WeightedHistogram::new(h.binning);
        p.fill_segment(5.0, 5.0, 1.5);
        assert_eq!(p.sums[2], 1.5);
        let mut m = h.clone();
        m.merge(&g);
        assert!((m.total() - h.total() - g.total()).abs() < 1e-12);
    }
}
