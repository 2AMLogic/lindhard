//! Radial point-spread function (PSF) of a pencil beam: extraction from the
//! electron deposition tally and double- or triple-Gaussian fits.
//!
//! # The radial profile
//!
//! A [`RadialAccumulator`] bins the energy deposited by the electron
//! transport in a depth slab `depth_lo <= x < depth_hi` by the distance
//! `r = (y² + z²)^(1/2)` from the beam axis (the axis of
//! [`CylindricalGrid`](super::CylindricalGrid)). The radial bins are
//! log-spaced ([`LogRadialBinning`]) so that one profile spans nanometres to
//! the backscatter range; bin 0 is the inner disc `[0, r_min)`. Deposits at
//! `r >= r_max` in the slab are kept in [`RadialProfile::beyond_ev`], so the
//! profile sums to the energy deposited in the slab. The accumulator is fed
//! by [`FullElectronTally`](super::FullElectronTally) when
//! [`ElectronTallyConfig::psf`](super::ElectronTallyConfig::psf) is set, from
//! the same deposits (and with the same bookkeeping convention) as its
//! cylindrical `r`-`x` grid; it is an `r`-`x` grid restricted to one depth
//! slab, with log radial bins.
//!
//! **Per-bin errors from per-history accumulation.** Within a history the
//! deposits of every bin are summed into a scratch array; at the end of the
//! history the bin's history total `x_h` is added to `S = Σ_h x_h` and its
//! square to `Q = Σ_h x_h²`. The standard error of the bin total `S` over `N`
//! histories is then `(N/(N-1) (Q - S²/N))^(1/2)`, `N` times the unbiased
//! sample variance of `x_h` taken to the square root (the usual batch
//! estimator with batches of one history). Energy deposited by the
//! secondaries of a history is part of that history. With `N < 2` the error
//! is reported as 0, and a fit then has no usable bins.
//!
//! # The fitted forms
//!
//! **Double Gaussian.** The normalised double-Gaussian proximity function of
//! T. H. P. Chang, "Proximity effect in electron-beam lithography",
//! J. Vac. Sci. Technol. 12, 1271 (1975), doi:10.1116/1.568515, as printed in
//! Q. Mao, J. Zhu, X. Cheng and Z. Wang, Discover Nano 20, 84 (2025),
//! doi:10.1186/s11671-025-04264-0 (open access, CC BY-NC-ND 4.0; PMC12089565), eq. (1),
//! section "Traditional Gaussian model":
//!
//! ```text
//! f(r) = 1/(π (1 + η)) [ (1/α²) exp(-r²/α²) + (η/β²) exp(-r²/β²) ]
//! ```
//!
//! with `α` the forward-scattering range, `β` the backscattering range and
//! `η` the ratio of backscattered to forward-scattered energy (Mao et al.'s
//! words). Chang's paper itself is closed access and was not seen (see
//! `docs/data-provenance.md`); Mao et al. attribute the form to it.
//!
//! **Triple Gaussian.** The sum-of-Gaussians PSF of T. Rosa Figueiro,
//! *Process modeling for proximity effect correction in electron beam
//! lithography*, PhD thesis, Université Grenoble Alpes (2015), HAL
//! tel-01206934, eq. (72), p. 111 (PDF page 112):
//!
//! ```text
//! PSF(x) = 1/(1 + Σ η_i) [ 1/(π α²) exp(-x²/α²) + Σ_i η_i/(π β_i²) exp(-x²/β_i²) ]
//! ```
//!
//! with two terms in the sum (the "3 Gaussian PSF model" of its Fig. 61).
//! Here `β_1 = β`, `η_1 = η`, `β_2 = γ`, `η_2 = ν`:
//!
//! ```text
//! f(r) = 1/(π (1 + η + ν)) [ (1/α²) e^(-r²/α²) + (η/β²) e^(-r²/β²) + (ν/γ²) e^(-r²/γ²) ]
//! ```
//!
//! With one term in the sum, eq. (72) is Mao et al.'s eq. (1). Both are
//! normalised, `∫₀^∞ f(r) 2πr dr = 1`; the fitted function is `E f(r)`, eV/m²,
//! with `E` (eV) a free parameter, the energy the PSF represents.
//!
//! # How the fit is done
//!
//! The data are the bin energies themselves (not their logarithm), and the
//! model of bin `[a, b)` is the exact annulus integral of `E f(r)`. Each
//! Gaussian term integrates in closed form (our own one-line integral,
//! checked numerically in the tests):
//!
//! ```text
//! ∫_a^b (1/(π w²)) exp(-r²/w²) 2πr dr = exp(-a²/w²) - exp(-b²/w²)
//! ```
//!
//! so there is no bin-centre approximation, however wide the log bins are.
//! The fit minimises the weighted sum of squares
//! `χ² = Σ_i ((E_i - μ_i)/σ_i)²` with `σ_i` the per-history standard error of
//! bin `i`; bins with `σ_i = 0` (no deposit, or fewer than two histories)
//! carry no information and are left out. Large dynamic ranges are handled by
//! the weights: a bin's error scales with its content.
//!
//! **Normalisation.** By default ([`PsfNormalization::SlabTotal`]) `E` is
//! fixed to the energy deposited in the slab ([`RadialProfile::total_ev`],
//! deposits beyond `r_max` included) and only the shape is fitted, so the
//! fitted function integrates to the slab energy whatever the quality of the
//! fit: the forms are normalised, and their scale is the deposited energy.
//! With [`PsfNormalization::Free`] `E` is fitted too. That recovers `E` when
//! the form describes the profile, but a Gaussian sum fitted to a profile it
//! does not describe follows the best-measured bins, and its integral can
//! then miss the slab energy by tens of percent. One case is a point beam,
//! which deposits a finite fraction of its energy on the axis itself, before
//! the first deflection; no Gaussian of finite width reproduces that.
//!
//! **The algorithm is our own**, built from textbook pieces not opened for
//! this implementation: a start from a grid search over the widths (for
//! fixed widths the model is linear in the term amplitudes, solved by
//! weighted linear least squares; starts with a non-positive amplitude are
//! skipped), then Levenberg-Marquardt iterations (K. Levenberg, Q. Appl.
//! Math. 2, 164 (1944); D. W. Marquardt, J. SIAM 11, 431 (1963)) on the
//! logarithms of the parameters (which keeps them positive) with the
//! analytic Jacobian. The covariance is `(Jᵀ W J)⁻¹` in the natural
//! parameters at the minimum, with `W = diag(1/σ_i²)`; it is **not** scaled
//! by the reduced `χ²`, so a reduced `χ²` well above 1 signals that the form
//! does not describe the profile and that the errors are too small. The
//! labels follow the start: the grid orders the widths `α < β (< γ)`.
//!
//! # Determinism
//!
//! The accumulator merges field by field in the order
//! [`crate::rng::run_particles`] merges chunk tallies, so the profile is
//! bit-identical at any thread count for a fixed chunk size, and the fit is a
//! fixed sequence of floating-point operations on the merged profile.
//!
//! # Export
//!
//! [`PsfReport`] holds a profile and its fits, derives `Serialize` and
//! `Deserialize` (callers serialise it, as for the other reports), and writes
//! the profile with the fitted model per bin ([`PsfReport::profile_csv`]) and
//! the parameters with their errors ([`PsfReport::parameters_csv`]) as CSV. Lengths are in metres, energies in
//! eV, areal densities in eV/m².

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Errors of the PSF extraction and fit.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum PsfError {
    /// The radial binning is not `0 < r_min < r_max` (finite) with `bins > 0`.
    #[error("invalid radial binning: r_min = {r_min} m, r_max = {r_max} m, bins = {bins} (need finite 0 < r_min < r_max and bins > 0)")]
    Binning {
        /// Inner edge of the first log bin.
        r_min: f64,
        /// Outer edge of the last bin.
        r_max: f64,
        /// Number of log bins.
        bins: usize,
    },
    /// The depth slab is not finite with `hi > lo`.
    #[error("invalid depth slab [{lo}, {hi}) m: need finite lo < hi")]
    Slab {
        /// Lower depth.
        lo: f64,
        /// Upper depth.
        hi: f64,
    },
    /// A profile's arrays are inconsistent.
    #[error("invalid radial profile: {0}")]
    Profile(String),
    /// Fewer usable bins (with a positive error) than parameters plus one.
    #[error("{used} usable bins for {parameters} parameters: need at least parameters + 1")]
    TooFewBins {
        /// Bins with a positive error.
        used: usize,
        /// Free parameters.
        parameters: usize,
    },
    /// No start with positive amplitudes was found.
    #[error("no start with positive term amplitudes was found")]
    NoStart,
    /// The curvature matrix at the minimum is singular.
    #[error("the curvature matrix is singular; the parameters are not determined by the profile")]
    Singular,
    /// PSF parameters are not positive and finite, or do not match the model.
    #[error("invalid PSF parameters: {0}")]
    Parameters(String),
}

/// Log-spaced radial bins: bin 0 is the disc `[0, r_min)`, bins `1..=bins`
/// split `[r_min, r_max)` into `bins` bins of equal width in `ln r`. Metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LogRadialBinning {
    /// Inner edge of the first log bin (outer edge of the central disc), m.
    pub r_min_m: f64,
    /// Outer edge of the last bin, m.
    pub r_max_m: f64,
    /// Number of log bins (the central disc is one more).
    pub bins: usize,
}

impl LogRadialBinning {
    /// Validated binning.
    pub fn new(r_min_m: f64, r_max_m: f64, bins: usize) -> Result<Self, PsfError> {
        let b = Self {
            r_min_m,
            r_max_m,
            bins,
        };
        b.validate()?;
        Ok(b)
    }

    fn validate(&self) -> Result<(), PsfError> {
        if self.r_min_m.is_finite()
            && self.r_max_m.is_finite()
            && self.r_min_m > 0.0
            && self.r_max_m > self.r_min_m
            && self.bins > 0
        {
            Ok(())
        } else {
            Err(PsfError::Binning {
                r_min: self.r_min_m,
                r_max: self.r_max_m,
                bins: self.bins,
            })
        }
    }

    /// Number of bins including the central disc, `bins + 1`.
    pub fn len(&self) -> usize {
        self.bins + 1
    }

    /// Always false: there is at least the central disc.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Inner edge of bin `i` (`i <= bins + 1`), m: `0`, then
    /// `r_min (r_max/r_min)^((i-1)/bins)`, and exactly `r_max` for the last.
    pub fn edge(&self, i: usize) -> f64 {
        if i == 0 {
            0.0
        } else if i > self.bins {
            self.r_max_m
        } else {
            let t = (i - 1) as f64 / self.bins as f64;
            self.r_min_m * (self.r_max_m / self.r_min_m).powf(t)
        }
    }

    /// All `bins + 2` edges, m.
    pub fn edges(&self) -> Vec<f64> {
        (0..=self.len()).map(|i| self.edge(i)).collect()
    }

    /// The bin of radius `r`, or `None` for `r >= r_max` (and NaN).
    pub fn locate(&self, r: f64) -> Option<usize> {
        if r < self.r_min_m {
            return Some(0);
        }
        if r.is_nan() || r >= self.r_max_m {
            return None;
        }
        let t = (r / self.r_min_m).ln() / (self.r_max_m / self.r_min_m).ln();
        let i = (t * self.bins as f64) as usize;
        Some(1 + i.min(self.bins - 1))
    }
}

/// Where the radial profile is taken: log radial bins and a depth slab
/// `depth_lo_m <= x < depth_hi_m`, m.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PsfConfig {
    /// Radial bins.
    pub radial: LogRadialBinning,
    /// Lower depth of the slab, m.
    pub depth_lo_m: f64,
    /// Upper depth of the slab (excluded), m.
    pub depth_hi_m: f64,
}

impl PsfConfig {
    /// Validated configuration.
    pub fn new(
        radial: LogRadialBinning,
        depth_lo_m: f64,
        depth_hi_m: f64,
    ) -> Result<Self, PsfError> {
        let c = Self {
            radial,
            depth_lo_m,
            depth_hi_m,
        };
        c.validate()?;
        Ok(c)
    }

    /// Checks the binning and the slab (a deserialised value skips `new`).
    pub fn validate(&self) -> Result<(), PsfError> {
        self.radial.validate()?;
        if self.depth_lo_m.is_finite()
            && self.depth_hi_m.is_finite()
            && self.depth_hi_m > self.depth_lo_m
        {
            Ok(())
        } else {
            Err(PsfError::Slab {
                lo: self.depth_lo_m,
                hi: self.depth_hi_m,
            })
        }
    }
}

/// Per-history accumulation of a radial profile (see
/// [the module docs](self#the-radial-profile)). Call
/// [`RadialAccumulator::deposit`] for every deposit of a history and
/// [`RadialAccumulator::end_history`] once at its end.
#[derive(Debug, Clone, PartialEq)]
pub struct RadialAccumulator {
    config: PsfConfig,
    histories: u64,
    /// `Σ_h x_h` per bin; the last entry is the slab beyond `r_max`, then the
    /// slab total.
    sum: Vec<f64>,
    /// `Σ_h x_h²`, same layout.
    sum_sq: Vec<f64>,
    /// The current history's totals, same layout.
    scratch: Vec<f64>,
    /// Entries of `scratch` touched in the current history (bin entries
    /// only; the total is always folded).
    touched: Vec<usize>,
}

impl RadialAccumulator {
    /// Empty accumulator.
    pub fn new(config: PsfConfig) -> Result<Self, PsfError> {
        config.validate()?;
        let n = config.radial.len() + 2;
        Ok(Self {
            config,
            histories: 0,
            sum: vec![0.0; n],
            sum_sq: vec![0.0; n],
            scratch: vec![0.0; n],
            touched: Vec::new(),
        })
    }

    /// The configuration.
    pub fn config(&self) -> &PsfConfig {
        &self.config
    }

    /// Histories ended so far.
    pub fn histories(&self) -> u64 {
        self.histories
    }

    /// A deposit of `e` eV at `pos` (m) in the current history. Ignored
    /// outside the depth slab or for `e <= 0`.
    pub fn deposit(&mut self, pos: [f64; 3], e: f64) {
        let in_slab = pos[0] >= self.config.depth_lo_m && pos[0] < self.config.depth_hi_m;
        if !(e.is_finite() && e > 0.0 && in_slab) {
            return;
        }
        let n = self.config.radial.len();
        let i = self.config.radial.locate(pos[1].hypot(pos[2])).unwrap_or(n);
        if self.scratch[i] == 0.0 {
            self.touched.push(i);
        }
        self.scratch[i] += e;
        self.scratch[n + 1] += e;
    }

    /// Ends the current history: folds its bin totals and their squares.
    pub fn end_history(&mut self) {
        self.histories += 1;
        // Fold in increasing bin order, independent of the deposit order.
        self.touched.sort_unstable();
        for &i in &self.touched {
            let x = self.scratch[i];
            self.sum[i] += x;
            self.sum_sq[i] += x * x;
            self.scratch[i] = 0.0;
        }
        self.touched.clear();
        let t = self.scratch.len() - 1;
        let x = self.scratch[t];
        if x != 0.0 {
            self.sum[t] += x;
            self.sum_sq[t] += x * x;
            self.scratch[t] = 0.0;
        }
    }

    /// Field-wise merge of the totals of another accumulator.
    ///
    /// # Panics
    /// If the configurations differ.
    pub fn merge(&mut self, o: &RadialAccumulator) {
        assert!(
            self.config == o.config,
            "merging radial accumulators with different configurations"
        );
        self.histories += o.histories;
        for (a, b) in self.sum.iter_mut().zip(&o.sum) {
            *a += b;
        }
        for (a, b) in self.sum_sq.iter_mut().zip(&o.sum_sq) {
            *a += b;
        }
    }

    /// The profile with per-bin standard errors.
    pub fn profile(&self) -> RadialProfile {
        let n = self.config.radial.len();
        let se = |i: usize| std_err_of_sum(self.sum[i], self.sum_sq[i], self.histories);
        RadialProfile {
            histories: self.histories,
            depth_lo_m: self.config.depth_lo_m,
            depth_hi_m: self.config.depth_hi_m,
            edges_m: self.config.radial.edges()[..=n].to_vec(),
            energy_ev: self.sum[..n].to_vec(),
            std_err_ev: (0..n).map(se).collect(),
            beyond_ev: self.sum[n],
            beyond_std_err_ev: se(n),
            total_ev: self.sum[n + 1],
            total_std_err_ev: se(n + 1),
        }
    }
}

/// Standard error of `S = Σ_h x_h` over `n` histories from `S` and
/// `Q = Σ_h x_h²`: `(n/(n-1) (Q - S²/n))^(1/2)`; 0 for `n < 2`.
fn std_err_of_sum(s: f64, q: f64, n: u64) -> f64 {
    if n < 2 {
        return 0.0;
    }
    let nf = n as f64;
    (nf / (nf - 1.0) * (q - s * s / nf)).max(0.0).sqrt()
}

/// A radial energy-deposition profile in a depth slab, summed over all
/// histories, with per-bin standard errors. Bin `i` is the annulus
/// `edges_m[i] <= r < edges_m[i + 1]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RadialProfile {
    /// Histories accumulated.
    pub histories: u64,
    /// Lower depth of the slab, m.
    pub depth_lo_m: f64,
    /// Upper depth of the slab (excluded), m.
    pub depth_hi_m: f64,
    /// Bin edges, m; one more than bins, starting at 0.
    pub edges_m: Vec<f64>,
    /// Energy per bin, eV.
    pub energy_ev: Vec<f64>,
    /// Standard error of each bin's energy, eV.
    pub std_err_ev: Vec<f64>,
    /// Energy deposited in the slab beyond the last edge, eV.
    pub beyond_ev: f64,
    /// Its standard error, eV.
    pub beyond_std_err_ev: f64,
    /// Energy deposited in the slab, eV: the bins plus `beyond_ev`.
    pub total_ev: f64,
    /// Its standard error, eV (from the per-history slab totals).
    pub total_std_err_ev: f64,
}

impl RadialProfile {
    /// Relative tolerance of the check that `total_ev` equals the bin
    /// energies plus `beyond_ev` (see [`RadialProfile::validate`]).
    ///
    /// [`RadialAccumulator`] sums the same deposits into the bins and into the
    /// slab total, but in a different order, so the two agree only to
    /// floating-point rounding: a naive sum of `n` non-negative terms is off
    /// by at most about `n ε` relative (`ε ≈ 1.1e-16`), well below `1e-6` for
    /// any feasible number of histories.
    pub const TOTAL_REL_TOLERANCE: f64 = 1e-6;

    /// Checks that the profile is consistent: edges start at 0 and increase,
    /// one energy and one error per bin, all finite, errors non-negative; the
    /// depth slab finite with `depth_lo_m < depth_hi_m`; `beyond_ev`,
    /// `total_ev` and their errors finite and non-negative; and `total_ev`
    /// equal to the bin energies plus `beyond_ev` within
    /// [`RadialProfile::TOTAL_REL_TOLERANCE`] relative.
    pub fn validate(&self) -> Result<(), PsfError> {
        let bad = |s: &str| Err(PsfError::Profile(s.to_string()));
        let n = self.energy_ev.len();
        if n == 0 || self.edges_m.len() != n + 1 || self.std_err_ev.len() != n {
            return bad("need n >= 1 bins with n + 1 edges, n energies and n errors");
        }
        if self.edges_m[0] != 0.0 || !self.edges_m.windows(2).all(|w| w[1] > w[0]) {
            return bad("edges must start at 0 and increase");
        }
        if !self.edges_m.iter().all(|e| e.is_finite())
            || !self.energy_ev.iter().all(|e| e.is_finite())
            || !self.std_err_ev.iter().all(|e| e.is_finite() && *e >= 0.0)
        {
            return bad("edges, energies and errors must be finite, errors non-negative");
        }
        if !(self.depth_lo_m.is_finite()
            && self.depth_hi_m.is_finite()
            && self.depth_lo_m < self.depth_hi_m)
        {
            return bad("depth bounds must be finite with depth_lo_m < depth_hi_m");
        }
        let scalars = [
            self.beyond_ev,
            self.beyond_std_err_ev,
            self.total_ev,
            self.total_std_err_ev,
        ];
        if !scalars.iter().all(|v| v.is_finite() && *v >= 0.0) {
            return bad(
                "beyond_ev, beyond_std_err_ev, total_ev and total_std_err_ev \
                 must be finite and non-negative",
            );
        }
        let parts = self.energy_ev.iter().sum::<f64>() + self.beyond_ev;
        let scale =
            self.energy_ev.iter().map(|e| e.abs()).sum::<f64>() + self.beyond_ev + self.total_ev;
        if (self.total_ev - parts).abs() > Self::TOTAL_REL_TOLERANCE * scale {
            return bad("total_ev must equal the bin energies plus beyond_ev");
        }
        Ok(())
    }

    /// Number of bins.
    pub fn len(&self) -> usize {
        self.energy_ev.len()
    }

    /// Whether there are no bins.
    pub fn is_empty(&self) -> bool {
        self.energy_ev.is_empty()
    }

    /// Area of the annulus of bin `i`, m².
    pub fn area_m2(&self, i: usize) -> f64 {
        let (a, b) = (self.edges_m[i], self.edges_m[i + 1]);
        PI * (b * b - a * a)
    }

    /// Representative radius of bin `i`, m: the geometric mean of its edges,
    /// and half the outer edge for the central disc.
    pub fn center_m(&self, i: usize) -> f64 {
        let (a, b) = (self.edges_m[i], self.edges_m[i + 1]);
        if a > 0.0 {
            (a * b).sqrt()
        } else {
            0.5 * b
        }
    }

    /// Mean areal energy density of each bin, eV/m².
    pub fn density_ev_per_m2(&self) -> Vec<f64> {
        (0..self.len())
            .map(|i| self.energy_ev[i] / self.area_m2(i))
            .collect()
    }

    /// Standard error of each bin's areal density, eV/m².
    pub fn density_std_err_ev_per_m2(&self) -> Vec<f64> {
        (0..self.len())
            .map(|i| self.std_err_ev[i] / self.area_m2(i))
            .collect()
    }
}

/// Which proximity function is fitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PsfModel {
    /// `(α, β, η)`: Mao et al. (2025) eq. (1), after Chang (1975).
    DoubleGaussian,
    /// `(α, β, η, γ, ν)`: Rosa Figueiro (2015) eq. (72) with two backscatter
    /// terms.
    TripleGaussian,
}

impl PsfModel {
    /// Number of Gaussian terms.
    pub fn terms(&self) -> usize {
        match self {
            PsfModel::DoubleGaussian => 2,
            PsfModel::TripleGaussian => 3,
        }
    }

    /// Number of free parameters, `E` included.
    pub fn parameters(&self) -> usize {
        2 * self.terms()
    }

    /// Parameter names in [`PsfFit::values`] order.
    pub fn parameter_names(&self) -> Vec<&'static str> {
        match self {
            PsfModel::DoubleGaussian => vec!["total_ev", "alpha_m", "beta_m", "eta"],
            PsfModel::TripleGaussian => {
                vec!["total_ev", "alpha_m", "beta_m", "eta", "gamma_m", "nu"]
            }
        }
    }

    /// Source of the normalised form, as recorded in a [`PsfFit`].
    pub fn source(&self) -> &'static str {
        match self {
            PsfModel::DoubleGaussian => {
                "f(r) = [exp(-r^2/alpha^2)/alpha^2 + eta exp(-r^2/beta^2)/beta^2] / (pi (1 + eta)): \
                 Mao, Zhu, Cheng and Wang, Discover Nano 20, 84 (2025), doi:10.1186/s11671-025-04264-0, \
                 eq. (1), after Chang, J. Vac. Sci. Technol. 12, 1271 (1975), doi:10.1116/1.568515"
            }
            PsfModel::TripleGaussian => {
                "f(r) = [exp(-r^2/alpha^2)/alpha^2 + eta exp(-r^2/beta^2)/beta^2 + nu exp(-r^2/gamma^2)/gamma^2] \
                 / (pi (1 + eta + nu)): Rosa Figueiro, PhD thesis, Univ. Grenoble Alpes (2015), \
                 HAL tel-01206934, eq. (72) with two backscatter terms (beta_1 = beta, eta_1 = eta, \
                 beta_2 = gamma, eta_2 = nu)"
            }
        }
    }
}

/// A double- or triple-Gaussian PSF scaled to the energy `total_ev`:
/// `E f(r)`, eV/m² (see [the module docs](self#the-fitted-forms)). The
/// triple form has `gamma_m` and `nu`; the double form has neither.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GaussianPsf {
    /// The energy `E` the PSF integrates to, eV.
    pub total_ev: f64,
    /// Forward-scattering width `α`, m.
    pub alpha_m: f64,
    /// Backscattering width `β`, m.
    pub beta_m: f64,
    /// Ratio `η` of the `β` term to the `α` term.
    pub eta: f64,
    /// Third width `γ`, m (triple Gaussian only).
    pub gamma_m: Option<f64>,
    /// Ratio `ν` of the `γ` term to the `α` term (triple Gaussian only).
    pub nu: Option<f64>,
}

impl GaussianPsf {
    /// Double Gaussian `(E, α, β, η)`.
    pub fn double(total_ev: f64, alpha_m: f64, beta_m: f64, eta: f64) -> Result<Self, PsfError> {
        let p = Self {
            total_ev,
            alpha_m,
            beta_m,
            eta,
            gamma_m: None,
            nu: None,
        };
        p.validate()?;
        Ok(p)
    }

    /// Triple Gaussian `(E, α, β, η, γ, ν)`.
    pub fn triple(
        total_ev: f64,
        alpha_m: f64,
        beta_m: f64,
        eta: f64,
        gamma_m: f64,
        nu: f64,
    ) -> Result<Self, PsfError> {
        let p = Self {
            total_ev,
            alpha_m,
            beta_m,
            eta,
            gamma_m: Some(gamma_m),
            nu: Some(nu),
        };
        p.validate()?;
        Ok(p)
    }

    /// The model: triple if `gamma_m` is set.
    pub fn model(&self) -> PsfModel {
        if self.gamma_m.is_some() {
            PsfModel::TripleGaussian
        } else {
            PsfModel::DoubleGaussian
        }
    }

    /// Checks that every parameter is positive and finite and that `gamma_m`
    /// and `nu` are both set or both unset.
    pub fn validate(&self) -> Result<(), PsfError> {
        if self.gamma_m.is_some() != self.nu.is_some() {
            return Err(PsfError::Parameters(
                "gamma_m and nu must both be set or both unset".into(),
            ));
        }
        if self.values().iter().all(|v| v.is_finite() && *v > 0.0) {
            Ok(())
        } else {
            Err(PsfError::Parameters(format!(
                "every parameter must be positive and finite, got {:?}",
                self.values()
            )))
        }
    }

    /// Parameters in [`PsfModel::parameter_names`] order.
    pub fn values(&self) -> Vec<f64> {
        let mut v = vec![self.total_ev, self.alpha_m, self.beta_m, self.eta];
        if let (Some(g), Some(n)) = (self.gamma_m, self.nu) {
            v.extend([g, n]);
        }
        v
    }

    fn from_values(model: PsfModel, v: &[f64]) -> Self {
        Self {
            total_ev: v[0],
            alpha_m: v[1],
            beta_m: v[2],
            eta: v[3],
            gamma_m: (model == PsfModel::TripleGaussian).then(|| v[4]),
            nu: (model == PsfModel::TripleGaussian).then(|| v[5]),
        }
    }

    /// `(widths, weights)` with the `α` term's weight 1.
    fn terms(&self) -> (Vec<f64>, Vec<f64>) {
        let mut w = vec![self.alpha_m, self.beta_m];
        let mut c = vec![1.0, self.eta];
        if let (Some(g), Some(n)) = (self.gamma_m, self.nu) {
            w.push(g);
            c.push(n);
        }
        (w, c)
    }

    /// The normalised PSF `f(r)`, 1/m².
    pub fn normalized_per_m2(&self, r_m: f64) -> f64 {
        let (w, c) = self.terms();
        let sum_c: f64 = c.iter().sum();
        let s: f64 = w
            .iter()
            .zip(&c)
            .map(|(w, c)| c / (w * w) * (-(r_m * r_m) / (w * w)).exp())
            .sum();
        s / (PI * sum_c)
    }

    /// The areal energy density `E f(r)`, eV/m².
    pub fn density_ev_per_m2(&self, r_m: f64) -> f64 {
        self.total_ev * self.normalized_per_m2(r_m)
    }

    /// Energy in the annulus `[a, b)`, eV, by the closed-form integral of
    /// [the module docs](self#how-the-fit-is-done).
    pub fn annulus_energy_ev(&self, a_m: f64, b_m: f64) -> f64 {
        let (w, c) = self.terms();
        let sum_c: f64 = c.iter().sum();
        let s: f64 = w
            .iter()
            .zip(&c)
            .map(|(w, c)| c * gauss_annulus(a_m, b_m, *w))
            .sum();
        self.total_ev * s / sum_c
    }
}

/// `exp(-a²/w²) - exp(-b²/w²)`, the fraction of a normalised 2D Gaussian
/// term of width `w` in the annulus `[a, b)`, without cancellation for thin
/// annuli.
fn gauss_annulus(a: f64, b: f64, w: f64) -> f64 {
    let w2 = w * w;
    -(-(a * a) / w2).exp() * (-(b * b - a * a) / w2).exp_m1()
}

/// `d/dw` of [`gauss_annulus`]: `(2/w³)(a² e^(-a²/w²) - b² e^(-b²/w²))`.
fn gauss_annulus_dw(a: f64, b: f64, w: f64) -> f64 {
    let w2 = w * w;
    2.0 / (w2 * w) * (a * a * (-(a * a) / w2).exp() - b * b * (-(b * b) / w2).exp())
}

/// How the energy scale `E` of the fitted PSF is set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum PsfNormalization {
    /// `E` is fixed to the energy deposited in the slab,
    /// [`RadialProfile::total_ev`] (beyond `r_max` included), and only the
    /// shape `(α, β, η[, γ, ν])` is fitted. The fitted function then
    /// integrates to the slab energy. This is the default: the forms are
    /// normalised PSFs, and their scale is the deposited energy.
    #[default]
    SlabTotal,
    /// `E` is a free parameter. Only a form that describes the profile
    /// (reduced `χ²` near 1) then integrates to the slab energy; on a
    /// profile it does not describe, the weighted fit follows the
    /// best-measured bins and `E` can be far off (by a third for the double
    /// Gaussian on the transport profile of `lindhard/tests/psf.rs`).
    Free,
}

/// Settings of [`fit_psf`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PsfFitOptions {
    /// How `E` is set.
    pub normalization: PsfNormalization,
    /// Start of the Levenberg-Marquardt iterations; `None` runs the grid
    /// search. Its model must match the fitted one; its `total_ev` is
    /// replaced by the slab energy under [`PsfNormalization::SlabTotal`].
    pub start: Option<GaussianPsf>,
    /// Points of the log-spaced width grid of the start search, from the
    /// first log edge to half the outer edge.
    pub grid_points: usize,
    /// Maximum Levenberg-Marquardt iterations.
    pub max_iterations: usize,
}

impl Default for PsfFitOptions {
    fn default() -> Self {
        Self {
            normalization: PsfNormalization::SlabTotal,
            start: None,
            grid_points: 40,
            max_iterations: 500,
        }
    }
}

/// One bin of a fit: data, model and normalised residual.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PsfResidual {
    /// Inner edge, m.
    pub r_lo_m: f64,
    /// Outer edge, m.
    pub r_hi_m: f64,
    /// Bin energy, eV.
    pub observed_ev: f64,
    /// Fitted model energy in the bin, eV.
    pub model_ev: f64,
    /// Standard error of the bin energy, eV.
    pub std_err_ev: f64,
    /// `(observed - model) / std_err`; `None` for a bin left out of the fit
    /// (zero error).
    pub normalized: Option<f64>,
}

/// The result of [`fit_psf`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PsfFit {
    /// The fitted form.
    pub model: PsfModel,
    /// Where the normalised form comes from.
    pub source: String,
    /// How the fit was done.
    pub method: String,
    /// How `E` was set.
    pub normalization: PsfNormalization,
    /// The fitted PSF.
    pub psf: GaussianPsf,
    /// Parameter names, in the order of `values`, `std_errors` and
    /// `covariance`.
    pub parameter_names: Vec<String>,
    /// Fitted values (`E` first, fixed or fitted).
    pub values: Vec<f64>,
    /// Standard errors, the square roots of the covariance diagonal.
    pub std_errors: Vec<f64>,
    /// Covariance matrix, row-major, not scaled by the reduced `χ²`: `(Jᵀ W
    /// J)⁻¹` over the fitted parameters. Under
    /// [`PsfNormalization::SlabTotal`] the `E` row and column hold only the
    /// variance of the slab energy ([`RadialProfile::total_std_err_ev`]
    /// squared); its covariances with the shape are not estimated and are 0.
    pub covariance: Vec<Vec<f64>>,
    /// `χ²` at the minimum.
    pub chi2: f64,
    /// Degrees of freedom: bins used minus fitted parameters.
    pub dof: usize,
    /// `χ² / dof`.
    pub reduced_chi2: f64,
    /// Bins used (positive error).
    pub bins_used: usize,
    /// Levenberg-Marquardt iterations.
    pub iterations: usize,
    /// Whether the iterations stopped on a converged step rather than the
    /// iteration cap.
    pub converged: bool,
    /// Per-bin data, model and residuals (every bin of the profile).
    pub residuals: Vec<PsfResidual>,
}

impl PsfFit {
    /// Standard error of the named parameter.
    pub fn std_error(&self, name: &str) -> Option<f64> {
        self.parameter_names
            .iter()
            .position(|n| n == name)
            .map(|i| self.std_errors[i])
    }
}

/// How the fit is done, as recorded in a [`PsfFit`].
pub const PSF_FIT_METHOD: &str = "weighted least squares on the bin energies (not their logarithm), \
     model = exact annulus integral of E f(r), weights 1/sigma_i^2 from per-history errors, bins with \
     sigma_i = 0 left out; start from a log grid over the widths with linear least squares for the \
     amplitudes; Levenberg-Marquardt on the logarithms of the parameters with the analytic Jacobian; \
     covariance (J^T W J)^-1 over the fitted parameters in their natural units, not scaled by the \
     reduced chi^2; E fixed to the slab energy (default) or fitted, see normalization";

/// The data a fit uses: bins with a positive error.
struct FitData {
    a: Vec<f64>,
    b: Vec<f64>,
    y: Vec<f64>,
    sigma: Vec<f64>,
}

impl FitData {
    fn chi2(&self, model: PsfModel, p: &[f64]) -> f64 {
        let psf = GaussianPsf::from_values(model, p);
        (0..self.y.len())
            .map(|i| {
                let r = (self.y[i] - psf.annulus_energy_ev(self.a[i], self.b[i])) / self.sigma[i];
                r * r
            })
            .sum()
    }

    /// Weighted residuals `(y - μ)/σ` and the Jacobian `(∂μ/∂p)/σ`, row-major
    /// (one row per bin), in the natural parameters.
    fn residuals_and_jacobian(&self, model: PsfModel, p: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let k = model.terms();
        let np = model.parameters();
        let e = p[0];
        let mut w = vec![p[1]];
        let mut c = vec![1.0];
        for t in 1..k {
            w.push(p[2 * t]);
            c.push(p[2 * t + 1]);
        }
        let sum_c: f64 = c.iter().sum();
        let n = self.y.len();
        let mut r = vec![0.0; n];
        let mut jac = vec![0.0; n * np];
        for i in 0..n {
            let g: Vec<f64> = w
                .iter()
                .map(|w| gauss_annulus(self.a[i], self.b[i], *w))
                .collect();
            let m: f64 = g.iter().zip(&c).map(|(g, c)| g * c).sum::<f64>() / sum_c;
            let s = self.sigma[i];
            r[i] = (self.y[i] - e * m) / s;
            let row = &mut jac[i * np..(i + 1) * np];
            row[0] = m / s;
            for t in 0..k {
                let dw = gauss_annulus_dw(self.a[i], self.b[i], w[t]);
                let iw = if t == 0 { 1 } else { 2 * t };
                row[iw] = e * c[t] * dw / sum_c / s;
                if t > 0 {
                    row[2 * t + 1] = e * (g[t] - m) / sum_c / s;
                }
            }
        }
        (r, jac)
    }
}

/// Fits the double- or triple-Gaussian PSF to a radial profile. See
/// [the module docs](self#how-the-fit-is-done).
pub fn fit_psf(
    profile: &RadialProfile,
    model: PsfModel,
    options: &PsfFitOptions,
) -> Result<PsfFit, PsfError> {
    profile.validate()?;
    let np = model.parameters();
    // Indices of the fitted parameters; `E` (index 0) only if free.
    let free: Vec<usize> = match options.normalization {
        PsfNormalization::Free => (0..np).collect(),
        PsfNormalization::SlabTotal => (1..np).collect(),
    };
    let nf = free.len();
    let used: Vec<usize> = (0..profile.len())
        .filter(|&i| profile.std_err_ev[i] > 0.0)
        .collect();
    if used.len() < nf + 1 {
        return Err(PsfError::TooFewBins {
            used: used.len(),
            parameters: nf,
        });
    }
    if options.normalization == PsfNormalization::SlabTotal && profile.total_ev <= 0.0 {
        return Err(PsfError::Profile(
            "no energy in the slab to normalise the PSF to".into(),
        ));
    }
    let data = FitData {
        a: used.iter().map(|&i| profile.edges_m[i]).collect(),
        b: used.iter().map(|&i| profile.edges_m[i + 1]).collect(),
        y: used.iter().map(|&i| profile.energy_ev[i]).collect(),
        sigma: used.iter().map(|&i| profile.std_err_ev[i]).collect(),
    };
    let mut start = match options.start {
        Some(s) => {
            s.validate()?;
            if s.model() != model {
                return Err(PsfError::Parameters(
                    "the start's model differs from the fitted model".into(),
                ));
            }
            s.values()
        }
        None => grid_start(profile, &data, model, options.grid_points.max(2))?,
    };
    if options.normalization == PsfNormalization::SlabTotal {
        start[0] = profile.total_ev;
    }

    // Levenberg-Marquardt on θ = ln p over the fitted parameters.
    let mut p = start;
    let mut chi2 = data.chi2(model, &p);
    let mut lambda = 1e-3;
    let mut iterations = 0;
    let mut converged = false;
    while iterations < options.max_iterations {
        iterations += 1;
        let (r, jac) = data.residuals_and_jacobian(model, &p);
        // Chain rule to θ: ∂μ/∂θ_j = p_j ∂μ/∂p_j.
        let mut a = vec![0.0; nf * nf];
        let mut g = vec![0.0; nf];
        for (i, ri) in r.iter().enumerate() {
            let row = &jac[i * np..(i + 1) * np];
            for (j, &fj) in free.iter().enumerate() {
                let jj = row[fj] * p[fj];
                g[j] += jj * ri;
                for (k, &fk) in free.iter().enumerate().take(j + 1) {
                    a[j * nf + k] += jj * row[fk] * p[fk];
                }
            }
        }
        for j in 0..nf {
            for k in 0..j {
                a[k * nf + j] = a[j * nf + k];
            }
        }
        let mut improved = false;
        while lambda <= 1e16 {
            let mut m = a.clone();
            for j in 0..nf {
                m[j * nf + j] += lambda * a[j * nf + j].max(f64::MIN_POSITIVE);
            }
            let Some(delta) = solve(m, g.clone(), nf) else {
                lambda *= 10.0;
                continue;
            };
            let mut pt = p.clone();
            for (d, &fj) in delta.iter().zip(&free) {
                pt[fj] = (p[fj].ln() + d).exp();
            }
            let c2 = data.chi2(model, &pt);
            if pt.iter().all(|v| v.is_finite() && *v > 0.0) && c2.is_finite() && c2 < chi2 {
                let step = delta.iter().fold(0.0f64, |m, d| m.max(d.abs()));
                let rel = (chi2 - c2) / c2.max(f64::MIN_POSITIVE);
                p = pt;
                chi2 = c2;
                lambda = (lambda * 0.1).max(1e-12);
                improved = true;
                if step < 1e-10 || rel < 1e-13 {
                    converged = true;
                }
                break;
            }
            lambda *= 10.0;
        }
        if !improved {
            // No step lowers χ²: a minimum to working precision.
            converged = true;
        }
        if converged {
            break;
        }
    }

    // Covariance over the fitted parameters, in their natural units.
    let (_, jac) = data.residuals_and_jacobian(model, &p);
    let n = data.y.len();
    let mut a = vec![0.0; nf * nf];
    for i in 0..n {
        let row = &jac[i * np..(i + 1) * np];
        for (j, &fj) in free.iter().enumerate() {
            for (k, &fk) in free.iter().enumerate() {
                a[j * nf + k] += row[fj] * row[fk];
            }
        }
    }
    let inv = invert(a, nf).ok_or(PsfError::Singular)?;
    let mut covariance = vec![vec![0.0; np]; np];
    for (j, &fj) in free.iter().enumerate() {
        for (k, &fk) in free.iter().enumerate() {
            covariance[fj][fk] = inv[j * nf + k];
        }
    }
    if options.normalization == PsfNormalization::SlabTotal {
        covariance[0][0] = profile.total_std_err_ev * profile.total_std_err_ev;
    }
    // Every variance finite; those of the fitted parameters also positive
    // (the fixed `E` of `SlabTotal` carries the slab total's variance, which
    // is 0 for fewer than two histories).
    if !(0..np).all(|j| covariance[j][j].is_finite())
        || !free.iter().all(|&j| covariance[j][j] > 0.0)
    {
        return Err(PsfError::Singular);
    }
    let std_errors = (0..np).map(|j| covariance[j][j].sqrt()).collect();
    let psf = GaussianPsf::from_values(model, &p);
    let residuals = (0..profile.len())
        .map(|i| {
            let (lo, hi) = (profile.edges_m[i], profile.edges_m[i + 1]);
            let model_ev = psf.annulus_energy_ev(lo, hi);
            let s = profile.std_err_ev[i];
            PsfResidual {
                r_lo_m: lo,
                r_hi_m: hi,
                observed_ev: profile.energy_ev[i],
                model_ev,
                std_err_ev: s,
                normalized: (s > 0.0).then(|| (profile.energy_ev[i] - model_ev) / s),
            }
        })
        .collect();
    let dof = n - nf;
    Ok(PsfFit {
        model,
        source: model.source().to_string(),
        method: PSF_FIT_METHOD.to_string(),
        normalization: options.normalization,
        psf,
        parameter_names: model
            .parameter_names()
            .into_iter()
            .map(String::from)
            .collect(),
        values: p,
        std_errors,
        covariance,
        chi2,
        dof,
        reduced_chi2: chi2 / dof as f64,
        bins_used: n,
        iterations,
        converged,
        residuals,
    })
}

/// The start: the best `χ²` over a log grid of ordered widths, with the
/// amplitudes of each grid point from weighted linear least squares.
fn grid_start(
    profile: &RadialProfile,
    data: &FitData,
    model: PsfModel,
    points: usize,
) -> Result<Vec<f64>, PsfError> {
    let k = model.terms();
    // First positive edge to half the outer edge.
    let lo = profile.edges_m[1];
    let hi = (0.5 * profile.edges_m[profile.len()]).max(lo * 1.0001);
    let widths: Vec<f64> = (0..points)
        .map(|i| lo * (hi / lo).powf(i as f64 / (points - 1) as f64))
        .collect();
    let n = data.y.len();
    // Weighted basis columns g_i(w)/σ_i and data y_i/σ_i.
    let cols: Vec<Vec<f64>> = widths
        .iter()
        .map(|&w| {
            (0..n)
                .map(|i| gauss_annulus(data.a[i], data.b[i], w) / data.sigma[i])
                .collect()
        })
        .collect();
    let yw: Vec<f64> = (0..n).map(|i| data.y[i] / data.sigma[i]).collect();
    let mut best: Option<(f64, Vec<usize>, Vec<f64>)> = None;
    let mut idx = vec![0usize; k];
    let mut consider = |idx: &[usize]| {
        let mut m = vec![0.0; k * k];
        let mut rhs = vec![0.0; k];
        for i in 0..n {
            for j in 0..k {
                let cj = cols[idx[j]][i];
                rhs[j] += cj * yw[i];
                for l in 0..k {
                    m[j * k + l] += cj * cols[idx[l]][i];
                }
            }
        }
        let Some(amp) = solve(m, rhs, k) else {
            return;
        };
        if !amp.iter().all(|a| a.is_finite() && *a > 0.0) {
            return;
        }
        let chi2: f64 = (0..n)
            .map(|i| {
                let mu: f64 = (0..k).map(|j| amp[j] * cols[idx[j]][i]).sum();
                (yw[i] - mu).powi(2)
            })
            .sum();
        if best.as_ref().is_none_or(|b| chi2 < b.0) {
            best = Some((chi2, idx.to_vec(), amp));
        }
    };
    // Every strictly increasing index tuple, in lexicographic order.
    fn next(idx: &mut [usize], points: usize) -> bool {
        let k = idx.len();
        let mut j = k;
        while j > 0 {
            j -= 1;
            if idx[j] < points - (k - j) {
                idx[j] += 1;
                for l in j + 1..k {
                    idx[l] = idx[l - 1] + 1;
                }
                return true;
            }
        }
        false
    }
    if points < k {
        return Err(PsfError::NoStart);
    }
    for (j, v) in idx.iter_mut().enumerate() {
        *v = j;
    }
    loop {
        consider(&idx);
        if !next(&mut idx, points) {
            break;
        }
    }
    let (_, idx, amp) = best.ok_or(PsfError::NoStart)?;
    // μ = Σ_t A_t g(w_t) with A_t = E c_t / Σc and c_0 = 1.
    let e: f64 = amp.iter().sum();
    let mut v = vec![e, widths[idx[0]]];
    for t in 1..k {
        v.push(widths[idx[t]]);
        v.push(amp[t] / amp[0]);
    }
    Ok(v)
}

/// Solves the `n x n` system `m x = r` (row-major) by Gaussian elimination
/// with partial pivoting; `None` if singular.
fn solve(mut m: Vec<f64>, mut r: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| m[i * n + c].abs().total_cmp(&m[j * n + c].abs()))?;
        if m[piv * n + c] == 0.0 || !m[piv * n + c].is_finite() {
            return None;
        }
        if piv != c {
            for k in 0..n {
                m.swap(c * n + k, piv * n + k);
            }
            r.swap(c, piv);
        }
        for i in c + 1..n {
            let f = m[i * n + c] / m[c * n + c];
            for k in c..n {
                m[i * n + k] -= f * m[c * n + k];
            }
            r[i] -= f * r[c];
        }
    }
    let mut x = vec![0.0; n];
    for c in (0..n).rev() {
        let s: f64 = (c + 1..n).map(|k| m[c * n + k] * x[k]).sum();
        x[c] = (r[c] - s) / m[c * n + c];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Inverse of an `n x n` matrix (row-major), column by column with
/// [`solve`]; `None` if singular.
fn invert(m: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; n * n];
    for c in 0..n {
        let mut e = vec![0.0; n];
        e[c] = 1.0;
        let x = solve(m.clone(), e, n)?;
        for r in 0..n {
            inv[r * n + c] = x[r];
        }
    }
    Some(inv)
}

/// A radial profile and its fits, for export.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PsfReport {
    /// The profile.
    pub profile: RadialProfile,
    /// Its fits.
    pub fits: Vec<PsfFit>,
}

impl PsfReport {
    /// The profile as CSV, one row per bin, with the model energy of each fit
    /// (`model_<model>_ev`) and the areal densities. Columns: `r_lo_m`,
    /// `r_hi_m`, `r_center_m`, `area_m2`, `energy_ev`, `std_err_ev`,
    /// `density_ev_per_m2`, `density_std_err_ev_per_m2`, then one per fit.
    pub fn profile_csv(&self) -> String {
        let p = &self.profile;
        let mut s = String::from(
            "r_lo_m,r_hi_m,r_center_m,area_m2,energy_ev,std_err_ev,density_ev_per_m2,density_std_err_ev_per_m2",
        );
        for f in &self.fits {
            s.push_str(&format!(",model_{}_ev", model_tag(f.model)));
        }
        s.push('\n');
        let d = p.density_ev_per_m2();
        let de = p.density_std_err_ev_per_m2();
        for i in 0..p.len() {
            s.push_str(&format!(
                "{},{},{},{},{},{},{},{}",
                p.edges_m[i],
                p.edges_m[i + 1],
                p.center_m(i),
                p.area_m2(i),
                p.energy_ev[i],
                p.std_err_ev[i],
                d[i],
                de[i]
            ));
            for f in &self.fits {
                s.push_str(&format!(",{}", f.residuals[i].model_ev));
            }
            s.push('\n');
        }
        s
    }

    /// The fitted parameters as CSV: `model,parameter,value,std_error`, then
    /// one `chi2`, `dof`, `reduced_chi2`, `converged` (1 or 0) and
    /// `iterations` row per fit (with an empty error). The last two are the
    /// fit's own diagnostics: a `converged` of 0 means the iteration budget
    /// ran out and the parameters are the last iterate.
    pub fn parameters_csv(&self) -> String {
        let mut s = String::from("model,parameter,value,std_error\n");
        for f in &self.fits {
            let tag = model_tag(f.model);
            for ((n, v), e) in f.parameter_names.iter().zip(&f.values).zip(&f.std_errors) {
                s.push_str(&format!("{tag},{n},{v},{e}\n"));
            }
            s.push_str(&format!("{tag},chi2,{},\n", f.chi2));
            s.push_str(&format!("{tag},dof,{},\n", f.dof));
            s.push_str(&format!("{tag},reduced_chi2,{},\n", f.reduced_chi2));
            s.push_str(&format!("{tag},converged,{},\n", u8::from(f.converged)));
            s.push_str(&format!("{tag},iterations,{},\n", f.iterations));
        }
        s
    }
}

fn model_tag(m: PsfModel) -> &'static str {
    match m {
        PsfModel::DoubleGaussian => "double_gaussian",
        PsfModel::TripleGaussian => "triple_gaussian",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_binning_edges_and_locate_agree() {
        let b = LogRadialBinning::new(1e-9, 1e-5, 40).unwrap();
        let e = b.edges();
        assert_eq!(e.len(), 42);
        assert_eq!(e[0], 0.0);
        assert_eq!(e[1], 1e-9);
        assert_eq!(e[41], 1e-5);
        assert!(e.windows(2).all(|w| w[1] > w[0]));
        assert_eq!(b.locate(0.0), Some(0));
        assert_eq!(b.locate(0.5e-9), Some(0));
        assert_eq!(b.locate(1e-5), None);
        assert_eq!(b.locate(f64::NAN), None);
        for i in 1..=40 {
            let mid = (e[i] * e[i + 1]).sqrt();
            assert_eq!(b.locate(mid), Some(i), "bin {i}");
        }
        assert!(LogRadialBinning::new(0.0, 1.0, 3).is_err());
        assert!(LogRadialBinning::new(1.0, 1.0, 3).is_err());
        assert!(LogRadialBinning::new(1.0, 2.0, 0).is_err());
        assert!(PsfConfig::new(b, 1.0, 1.0).is_err());
        assert!(PsfConfig::new(b, 0.0, f64::INFINITY).is_err());
    }

    #[test]
    fn annulus_integral_matches_quadrature() {
        let p = GaussianPsf::triple(7.0, 3e-9, 2e-7, 0.7, 4e-6, 0.2).unwrap();
        for (a, b) in [(0.0, 1e-9), (2e-9, 9e-9), (1e-7, 5e-7), (3e-6, 1e-5)] {
            // Simpson's rule on E f(r) 2πr dr.
            let n = 20_000;
            let h = (b - a) / n as f64;
            let f = |r: f64| p.density_ev_per_m2(r) * 2.0 * PI * r;
            let mut s = f(a) + f(b);
            for i in 1..n {
                s += f(a + i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
            }
            let q = s * h / 3.0;
            let c = p.annulus_energy_ev(a, b);
            assert!(
                (q - c).abs() <= 1e-9 * c.abs().max(1e-300),
                "[{a}, {b}): {q} vs {c}"
            );
        }
        // Normalised: the whole plane holds E.
        assert!((p.annulus_energy_ev(0.0, 1.0) - 7.0).abs() < 1e-12);
        let d = GaussianPsf::double(1.0, 1e-9, 1e-7, 0.5).unwrap();
        assert!((d.annulus_energy_ev(0.0, 1.0) - 1.0).abs() < 1e-15);
        // f(0) = (1/α² + η/β²) / (π (1 + η)).
        let f0 = (1.0 / 1e-18 + 0.5 / 1e-14) / (PI * 1.5);
        assert!((d.normalized_per_m2(0.0) - f0).abs() < 1e-12 * f0);
    }

    #[test]
    fn analytic_jacobian_matches_finite_differences() {
        for model in [PsfModel::DoubleGaussian, PsfModel::TripleGaussian] {
            let p: Vec<f64> = match model {
                PsfModel::DoubleGaussian => vec![100.0, 5e-9, 3e-7, 0.6],
                PsfModel::TripleGaussian => vec![100.0, 5e-9, 3e-7, 0.6, 4e-6, 0.3],
            };
            let b = LogRadialBinning::new(1e-10, 1e-4, 30).unwrap();
            let e = b.edges();
            let data = FitData {
                a: e[..e.len() - 1].to_vec(),
                b: e[1..].to_vec(),
                // y = 0, so the residuals are -μ without cancellation.
                y: vec![0.0; e.len() - 1],
                sigma: vec![1.0; e.len() - 1],
            };
            let np = model.parameters();
            let (r0, jac) = data.residuals_and_jacobian(model, &p);
            for j in 0..np {
                let col_max = (0..r0.len()).fold(0.0f64, |m, i| m.max(jac[i * np + j].abs()));
                let h = 1e-6 * p[j];
                let mut pp = p.clone();
                pp[j] += h;
                let mut pm = p.clone();
                pm[j] -= h;
                let (rp, _) = data.residuals_and_jacobian(model, &pp);
                let (rm, _) = data.residuals_and_jacobian(model, &pm);
                for i in 0..r0.len() {
                    let fd = -(rp[i] - rm[i]) / (2.0 * h);
                    let an = jac[i * np + j];
                    assert!(
                        (fd - an).abs() <= 1e-6 * an.abs() + 1e-8 * col_max,
                        "{model:?} p{j} bin {i}: {fd} vs {an}"
                    );
                }
            }
        }
    }

    #[test]
    fn accumulator_errors_follow_per_history_totals() {
        let b = LogRadialBinning::new(1.0, 100.0, 2).unwrap();
        let cfg = PsfConfig::new(b, 0.0, 1.0).unwrap();
        let mut acc = RadialAccumulator::new(cfg).unwrap();
        // Bin 0: r < 1. Bin 1: [1, 10). Bin 2: [10, 100).
        let hist: [&[(f64, f64)]; 3] = [
            &[(0.5, 2.0), (0.2, 1.0), (5.0, 3.0)],
            &[(5.0, 1.0)],
            &[(50.0, 4.0), (500.0, 1.0)],
        ];
        for h in hist {
            for &(r, e) in h {
                acc.deposit([0.5, r, 0.0], e);
            }
            // Outside the slab: ignored.
            acc.deposit([2.0, 0.0, 0.0], 100.0);
            acc.end_history();
        }
        let p = acc.profile();
        assert_eq!(p.histories, 3);
        assert_eq!(p.energy_ev, vec![3.0, 4.0, 4.0]);
        assert_eq!(p.beyond_ev, 1.0);
        assert_eq!(p.total_ev, 12.0);
        // Bin 0: per-history totals 3, 0, 0.
        let se = |xs: [f64; 3]| {
            let m = xs.iter().sum::<f64>() / 3.0;
            let v = xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / 2.0;
            (3.0 * v).sqrt()
        };
        assert!((p.std_err_ev[0] - se([3.0, 0.0, 0.0])).abs() < 1e-12);
        assert!((p.std_err_ev[1] - se([3.0, 1.0, 0.0])).abs() < 1e-12);
        assert!((p.total_std_err_ev - se([6.0, 1.0, 5.0])).abs() < 1e-12);
        p.validate().unwrap();
        // Merging two halves equals one accumulator.
        let mut a = RadialAccumulator::new(cfg).unwrap();
        let mut c = RadialAccumulator::new(cfg).unwrap();
        for (k, h) in hist.iter().enumerate() {
            let t = if k == 0 { &mut a } else { &mut c };
            for &(r, e) in *h {
                t.deposit([0.5, r, 0.0], e);
            }
            t.end_history();
        }
        a.merge(&c);
        assert_eq!(a.profile(), acc.profile());
    }

    #[test]
    fn solver_and_inverse() {
        let m = vec![4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let x = solve(m.clone(), vec![1.0, 2.0, 3.0], 3).unwrap();
        let r: Vec<f64> = (0..3)
            .map(|i| (0..3).map(|k| m[i * 3 + k] * x[k]).sum())
            .collect();
        for (a, b) in r.iter().zip([1.0, 2.0, 3.0]) {
            assert!((a - b).abs() < 1e-12);
        }
        let inv = invert(m.clone(), 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let s: f64 = (0..3).map(|k| m[i * 3 + k] * inv[k * 3 + j]).sum();
                assert!((s - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
        assert!(solve(vec![1.0, 2.0, 2.0, 4.0], vec![1.0, 1.0], 2).is_none());
    }
}
