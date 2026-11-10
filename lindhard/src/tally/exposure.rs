//! Exposure statistics: the distribution of the energy a small volume
//! receives from a finite number of incident electrons.
//!
//! A tally reports means over many histories. A small volume exposed to a few
//! electrons receives an energy that fluctuates, both from the number of
//! electrons and from the wide spread of the energy each one leaves there.
//!
//! # The per-electron distribution `X`
//!
//! [`ExposureAccumulator`] follows each primary history (the primary and all
//! its secondaries) and sums the energy it deposits inside every configured
//! [`ExposureVolume`]. That sum, `X` (eV), is binned into an
//! [`ExposureHistogram`]: an exact atom at `X = 0` (most histories deposit
//! nothing in a small volume), `bins` equal-width bins of width
//! `bin_width_ev` on `(0, bins * bin_width_ev]` (a bin is closed on its
//! upper edge: bin `j`, `j = 1..=bins`, holds `(j-1) w < X <= j w`), and the
//! running sums `Σ X^k` for `k = 1..=4`. A history with `X` beyond the last
//! edge is counted in `overflow_count` and clamped into the last bin; the
//! sums stay exact.
//!
//! # Beam profile and the lateral-uniformity assumption
//!
//! The beam is a [`BeamProfile`]: a point, a Gaussian of standard deviation
//! `sigma` per lateral axis, or a uniform disc. The landing offset `d` of each
//! electron is drawn from the counter-based stream
//! `rng::stream(seed, history index)` at the start of its history, and the
//! deposits are taken relative to the landing point (`pos - start - d`
//! compared with the volume centre, i.e. the volume is shifted by `-d`).
//! The stored `X` is therefore already averaged over the beam and the
//! electrons stay independent and identically distributed.
//!
//! **This is valid only while the target is laterally uniform** (a layered
//! stack, as everywhere in the electron engine): then the history of an
//! electron landing at `d` is the history of one landing at the axis, shifted
//! by `d`. Laterally varying targets and beams are out of scope here.
//!
//! # Compound Poisson statistics
//!
//! With `N ~ Poisson(λ)` electrons incident (`λ` the mean number), the energy
//! in the volume is `S = X_1 + … + X_N`, `X_i` independent copies of `X`.
//! Conditioning on `N` (law of total expectation and variance) gives
//!
//! ```text
//! E[S]   = E[N] E[X]                    = λ E[X]
//! Var[S] = E[N] Var[X] + Var[N] E[X]²   = λ E[X²]
//! ```
//!
//! because `E[N] = Var[N] = λ`. These are the textbook compound-Poisson
//! results (see W. Feller, *An Introduction to Probability Theory and Its
//! Applications*, vol. 1, 3rd ed., Wiley 1968, the chapter on compound
//! distributions; and N. R. Campbell, Proc. Camb. Phil. Soc. 15, 117 (1909),
//! for the mean and variance of such shot-noise sums). The derivation above
//! is self-contained; the chapter-level citation was not checked against the
//! page. The moments use the exact stored sums, with no binning error.
//!
//! The full distribution of `S` uses Poisson thinning: the histories that
//! deposit something number `K ~ Poisson(λ (1 - p0))`, `p0` the mass at
//! `X = 0`, and `S` is the sum of `K` draws from the conditional nonzero
//! distribution `f`. So
//!
//! ```text
//! P(S = s) = Σ_k Poisson(k; λ(1-p0)) · f^{*k}(s)
//! ```
//!
//! with `f^{*k}` the `k`-fold convolution of the binned `f`, built
//! iteratively, and the `k = 0` term the exact atom at `S = 0`. The sum is
//! truncated where the Poisson tail is below `tail_tolerance`; the dropped
//! mass is reported.
//!
//! **Binning bound.** Representing each bin by its lower edge, centre or
//! upper edge gives three distributions of `S`; the first is stochastically
//! smaller and the last larger than the true one (up to the clamped
//! overflow, reported). Every quantile and threshold probability is reported
//! as `low`, `central` (bin centres) and `high`: the true value lies
//! between `low` and `high` for quantiles and for `P(S > b)`, and the same
//! holds for `P(S < a)` with the roles of the two ends swapped (the field
//! names refer to the representation, see [`Bracket`]). The width of the
//! bracket on a quantile is at most `K_max * bin_width_ev`, `K_max` the
//! largest number of depositing electrons kept.
//!
//! Units: energies in eV, lengths in metres.

use crate::electron::transport::{disc_offset, finite_positive, gaussian_offset};
use crate::rng::{stream, uniform};
use serde::{Deserialize, Serialize};

/// Largest `bins * k_max * bins` multiply-add count the convolution is
/// allowed (about a second).
const MAX_CONVOLUTION_WORK: f64 = 2.0e9;

/// Errors of the exposure configuration.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ExposureError {
    /// A scalar of the configuration is out of range.
    #[error("exposure: {0}")]
    Invalid(String),
}

fn invalid<T>(msg: impl Into<String>) -> Result<T, ExposureError> {
    Err(ExposureError::Invalid(msg.into()))
}

/// Lateral profile of the beam.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum BeamProfile {
    /// All electrons land on the axis.
    Point,
    /// Landing offset Gaussian with standard deviation `sigma_m` per axis.
    Gaussian {
        /// Standard deviation per lateral axis, m.
        sigma_m: f64,
    },
    /// Landing offset uniform on a disc.
    UniformDisc {
        /// Disc radius, m.
        radius_m: f64,
    },
}

impl BeamProfile {
    /// The landing offset `(y, z)` of history `index`, drawn from
    /// `rng::stream(seed, index)`.
    pub fn offset(&self, seed: u64, index: u64) -> [f64; 2] {
        match *self {
            BeamProfile::Point => [0.0, 0.0],
            BeamProfile::Gaussian { sigma_m } => {
                let mut rng = stream(seed, index);
                let (u1, u2) = (uniform(&mut rng), uniform(&mut rng));
                gaussian_offset(sigma_m, u1, u2)
            }
            BeamProfile::UniformDisc { radius_m } => {
                let mut rng = stream(seed, index);
                let (u1, u2) = (uniform(&mut rng), uniform(&mut rng));
                disc_offset(radius_m, u1, u2)
            }
        }
    }

    fn validate(&self) -> Result<(), ExposureError> {
        match *self {
            BeamProfile::Point => Ok(()),
            BeamProfile::Gaussian { sigma_m } if finite_positive(sigma_m) => Ok(()),
            BeamProfile::UniformDisc { radius_m } if finite_positive(radius_m) => Ok(()),
            _ => invalid("beam sigma or radius must be positive and finite"),
        }
    }
}

/// Lateral shape of a target volume, centred on the volume's offset.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum VolumeShape {
    /// A cylinder along the depth axis.
    Cylinder {
        /// Radius, m.
        radius_m: f64,
    },
    /// A box with the given half widths along `y` and `z`.
    Box {
        /// Half width along `y`, m.
        half_y_m: f64,
        /// Half width along `z`, m.
        half_z_m: f64,
    },
}

/// A target volume: a lateral shape over a depth range `depth_lo_m <= x <
/// depth_hi_m`, centred at lateral offset `(offset_y_m, offset_z_m)` from the
/// beam axis. Deposits exist only inside the target, so a volume extending
/// beyond the film collects only the part inside it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExposureVolume {
    /// Label in the outputs.
    pub name: String,
    /// Lateral shape.
    pub shape: VolumeShape,
    /// Lower depth, m.
    pub depth_lo_m: f64,
    /// Upper depth (excluded), m.
    pub depth_hi_m: f64,
    /// Lateral centre `y`, m.
    pub offset_y_m: f64,
    /// Lateral centre `z`, m.
    pub offset_z_m: f64,
}

impl ExposureVolume {
    fn validate(&self) -> Result<(), ExposureError> {
        let ok = |v: f64| v.is_finite();
        let shape_ok = match self.shape {
            VolumeShape::Cylinder { radius_m } => ok(radius_m) && radius_m > 0.0,
            VolumeShape::Box { half_y_m, half_z_m } => {
                ok(half_y_m) && half_y_m > 0.0 && ok(half_z_m) && half_z_m > 0.0
            }
        };
        if !shape_ok {
            return invalid(format!("volume '{}': size must be positive", self.name));
        }
        if !(ok(self.depth_lo_m)
            && ok(self.depth_hi_m)
            && self.depth_hi_m > self.depth_lo_m
            && ok(self.offset_y_m)
            && ok(self.offset_z_m))
        {
            return invalid(format!(
                "volume '{}': need finite offsets and depth_hi > depth_lo",
                self.name
            ));
        }
        Ok(())
    }

    /// Whether `pos` (depth, then the lateral position relative to the
    /// landing point shifted into the volume's frame) lies inside.
    fn contains(&self, x: f64, dy: f64, dz: f64) -> bool {
        if !(x >= self.depth_lo_m && x < self.depth_hi_m) {
            return false;
        }
        let (y, z) = (dy - self.offset_y_m, dz - self.offset_z_m);
        match self.shape {
            VolumeShape::Cylinder { radius_m } => y.hypot(z) < radius_m,
            VolumeShape::Box { half_y_m, half_z_m } => y.abs() < half_y_m && z.abs() < half_z_m,
        }
    }
}

/// Settings of the exposure tally.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExposureConfig {
    /// The volumes (a neighbour is simply another volume at its offset).
    pub volumes: Vec<ExposureVolume>,
    /// Beam profile.
    pub beam: BeamProfile,
    /// Seed of the landing-offset stream (independent of the transport's).
    pub seed: u64,
    /// Bins of the histogram of `X`.
    pub bins: usize,
    /// Upper edge of the histogram, eV.
    pub x_max_ev: f64,
    /// Mean number of incident electrons `λ` of the exposure.
    pub lambda: f64,
    /// Probabilities `p` in `(0, 1)` whose quantiles of `S` are reported.
    pub quantiles: Vec<f64>,
    /// Thresholds `a`, eV, for `P(S < a)`.
    pub below_ev: Vec<f64>,
    /// Thresholds `b`, eV, for `P(S > b)`.
    pub above_ev: Vec<f64>,
    /// Poisson tail mass dropped from the convolution.
    pub tail_tolerance: f64,
}

impl ExposureConfig {
    /// Checks the configuration.
    pub fn validate(&self) -> Result<(), ExposureError> {
        if self.volumes.is_empty() {
            return invalid("at least one volume is required");
        }
        for v in &self.volumes {
            v.validate()?;
        }
        self.beam.validate()?;
        if self.bins == 0 || !(self.x_max_ev.is_finite() && self.x_max_ev > 0.0) {
            return invalid("bins must be at least 1 and x_max_ev positive");
        }
        if !(self.lambda.is_finite() && self.lambda > 0.0) {
            return invalid("lambda must be positive");
        }
        if self.quantiles.iter().any(|p| !(*p > 0.0 && *p < 1.0)) {
            return invalid("quantiles must lie in (0, 1)");
        }
        if self
            .below_ev
            .iter()
            .chain(&self.above_ev)
            .any(|t| !t.is_finite())
        {
            return invalid("thresholds must be finite");
        }
        if !(self.tail_tolerance > 0.0 && self.tail_tolerance < 1.0) {
            return invalid("tail_tolerance must lie in (0, 1)");
        }
        Ok(())
    }
}

/// The histogram of the per-electron deposit `X` in one volume, with exact
/// moment sums. Mergeable; see [the module docs](self).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExposureHistogram {
    /// Primary histories tallied.
    pub histories: u64,
    /// Histories with `X = 0` exactly.
    pub zero_count: u64,
    /// Bin width, eV.
    pub bin_width_ev: f64,
    /// Counts of `(j-1) w < X <= j w`, `j = 1..=bins` (index `j-1`).
    pub counts: Vec<u64>,
    /// Histories with `X` beyond the last edge (also counted in the last bin).
    pub overflow_count: u64,
    /// `Σ X`, eV.
    pub sum_x: f64,
    /// `Σ X²`, eV².
    pub sum_x2: f64,
    /// `Σ X³`, eV³.
    pub sum_x3: f64,
    /// `Σ X⁴`, eV⁴.
    pub sum_x4: f64,
}

impl ExposureHistogram {
    fn new(bins: usize, x_max_ev: f64) -> Self {
        Self {
            histories: 0,
            zero_count: 0,
            bin_width_ev: x_max_ev / bins as f64,
            counts: vec![0; bins],
            overflow_count: 0,
            sum_x: 0.0,
            sum_x2: 0.0,
            sum_x3: 0.0,
            sum_x4: 0.0,
        }
    }

    fn push(&mut self, x: f64) {
        self.histories += 1;
        if x <= 0.0 {
            self.zero_count += 1;
            return;
        }
        let bins = self.counts.len();
        let j = (x / self.bin_width_ev).ceil() as usize; // 1-based
        if j > bins {
            self.overflow_count += 1;
        }
        self.counts[j.clamp(1, bins) - 1] += 1;
        let x2 = x * x;
        self.sum_x += x;
        self.sum_x2 += x2;
        self.sum_x3 += x2 * x;
        self.sum_x4 += x2 * x2;
    }

    fn merge(&mut self, o: &Self) {
        self.histories += o.histories;
        self.zero_count += o.zero_count;
        self.overflow_count += o.overflow_count;
        for (a, b) in self.counts.iter_mut().zip(&o.counts) {
            *a += b;
        }
        self.sum_x += o.sum_x;
        self.sum_x2 += o.sum_x2;
        self.sum_x3 += o.sum_x3;
        self.sum_x4 += o.sum_x4;
    }

    /// `E[X]`, eV (0 with no histories).
    pub fn mean_ev(&self) -> f64 {
        self.sum_x / self.histories.max(1) as f64
    }

    /// `E[X²]`, eV².
    pub fn second_moment_ev2(&self) -> f64 {
        self.sum_x2 / self.histories.max(1) as f64
    }

    /// Standard error of `E[X]`: `sqrt((E[X²] - E[X]²)/(n-1))`; 0 for `n < 2`.
    pub fn mean_std_err_ev(&self) -> f64 {
        let n = self.histories as f64;
        if self.histories < 2 {
            return 0.0;
        }
        let var = (self.second_moment_ev2() - self.mean_ev().powi(2)).max(0.0);
        (var * n / (n - 1.0) / n).sqrt()
    }

    /// Standard error of `E[X²]` from `E[X⁴] - E[X²]²`.
    pub fn second_moment_std_err_ev2(&self) -> f64 {
        let n = self.histories as f64;
        if self.histories < 2 {
            return 0.0;
        }
        let m2 = self.second_moment_ev2();
        let var = (self.sum_x4 / n - m2 * m2).max(0.0);
        (var * n / (n - 1.0) / n).sqrt()
    }

    /// Fraction of histories with `X = 0`.
    pub fn zero_fraction(&self) -> f64 {
        self.zero_count as f64 / self.histories.max(1) as f64
    }

    /// The compound-Poisson summary of `S = X_1 + … + X_N`, `N ~ Poisson(λ)`,
    /// for any `λ`, from this histogram alone (no rerun). See
    /// [the module docs](self#compound-poisson-statistics).
    pub fn compound_poisson(
        &self,
        lambda: f64,
        request: &CompoundRequest,
    ) -> Result<CompoundPoisson, ExposureError> {
        if !(lambda.is_finite() && lambda > 0.0) {
            return invalid("lambda must be positive");
        }
        if self.histories == 0 {
            return invalid("no histories");
        }
        let mean = lambda * self.mean_ev();
        let variance = lambda * self.second_moment_ev2();
        let mut out = CompoundPoisson {
            lambda,
            mean_ev: mean,
            mean_std_err_ev: lambda * self.mean_std_err_ev(),
            variance_ev2: variance,
            variance_std_err_ev2: lambda * self.second_moment_std_err_ev2(),
            p_zero: 0.0,
            k_max: 0,
            truncated_mass: 0.0,
            bin_width_ev: self.bin_width_ev,
            overflow_fraction: self.overflow_count as f64 / self.histories as f64,
            distribution: None,
            note: None,
        };
        let p0 = self.zero_fraction();
        let nz_total: u64 = self.counts.iter().sum();
        // Atom at S = 0: no electron deposits anything.
        out.p_zero = (-lambda * (1.0 - p0)).exp();
        if nz_total == 0 {
            out.k_max = 0;
            out.distribution = Some(distribution_atom_at_zero(request));
            return Ok(out);
        }
        let mu = lambda * (1.0 - p0); // mean number of depositing electrons
        let (poisson, k_max, trunc) = poisson_pmf(mu, request.tail_tolerance);
        out.k_max = k_max;
        out.truncated_mass = trunc;
        let bins = self.counts.len();
        let work = (k_max as f64).powi(2) * (bins as f64).powi(2) / 2.0;
        if work > MAX_CONVOLUTION_WORK {
            out.note = Some(format!(
                "distribution not computed: lambda too large for the convolution \
                 (k_max {k_max}, {bins} bins); mean and variance are exact"
            ));
            return Ok(out);
        }
        let f: Vec<f64> = self
            .counts
            .iter()
            .map(|&c| c as f64 / nz_total as f64)
            .collect();
        // Mass on the integer grid m = Σ j_i, in half-units of the bin width:
        // lower edge 2(m-k), centre 2m-k, upper 2m. Index by 2m (so 0..=2 k_max bins).
        let size = 2 * k_max * bins + 1;
        let mut lo = vec![0.0; size];
        let mut mid = vec![0.0; size];
        let mut hi = vec![0.0; size];
        lo[0] += poisson[0];
        mid[0] += poisson[0];
        hi[0] += poisson[0];
        let mut g = vec![1.0]; // f^{*0}, support m = 0
        for k in 1..=k_max {
            // g <- g * f; support of g is m in [k-1, (k-1) bins] -> [k, k bins].
            let mut next = vec![0.0; g.len() + bins];
            for (m, &gm) in g.iter().enumerate() {
                if gm == 0.0 {
                    continue;
                }
                for (j, &fj) in f.iter().enumerate() {
                    next[m + j + 1] += gm * fj;
                }
            }
            g = next;
            // g index is m directly (shifted by the +1 per factor above).
            let pk = poisson[k];
            for (m, &gm) in g.iter().enumerate() {
                if gm == 0.0 {
                    continue;
                }
                let w = pk * gm;
                // Values in half-units: lower 2(m-k), centre 2m-k, upper 2m.
                lo[2 * (m - k)] += w;
                mid[2 * m - k] += w;
                hi[2 * m] += w;
            }
        }
        out.distribution = Some(summarise(self.bin_width_ev, &lo, &mid, &hi, request));
        Ok(out)
    }
}

/// Which compound-Poisson outputs to compute.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CompoundRequest {
    /// Probabilities whose quantiles of `S` are wanted.
    pub quantiles: Vec<f64>,
    /// Thresholds `a` for `P(S < a)`, eV.
    pub below_ev: Vec<f64>,
    /// Thresholds `b` for `P(S > b)`, eV.
    pub above_ev: Vec<f64>,
    /// Dropped Poisson tail mass.
    pub tail_tolerance: f64,
}

/// A value bracketed by the binning representations of `X` (see
/// [the module docs](self#compound-poisson-statistics)). `low` is the
/// result with every bin at its lower edge, `high` at its upper edge,
/// `central` at the bin centres. For a quantile and for `P(S > b)` the true
/// value lies in `[low, high]`; for `P(S < a)` it lies in `[high, low]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Bracket {
    /// Lower-edge representation.
    pub low: f64,
    /// Bin-centre representation.
    pub central: f64,
    /// Upper-edge representation.
    pub high: f64,
}

/// A reported quantile of `S`, eV.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuantileResult {
    /// The probability.
    pub p: f64,
    /// The quantile, eV: the least `s` with `P(S <= s) >= p`.
    pub energy_ev: Bracket,
}

/// A reported threshold probability.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThresholdResult {
    /// The threshold, eV.
    pub threshold_ev: f64,
    /// The probability.
    pub probability: Bracket,
}

/// The discrete distribution results of a compound-Poisson sum.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CompoundDistribution {
    /// Quantiles of `S`.
    pub quantiles: Vec<QuantileResult>,
    /// `P(S < a)`.
    pub below: Vec<ThresholdResult>,
    /// `P(S > b)`.
    pub above: Vec<ThresholdResult>,
}

/// Compound-Poisson summary of the energy `S` in a volume.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CompoundPoisson {
    /// Mean number of incident electrons.
    pub lambda: f64,
    /// `E[S] = λ E[X]`, eV.
    pub mean_ev: f64,
    /// Standard error of `mean_ev`, eV.
    pub mean_std_err_ev: f64,
    /// `Var[S] = λ E[X²]`, eV².
    pub variance_ev2: f64,
    /// Standard error of `variance_ev2`, eV².
    pub variance_std_err_ev2: f64,
    /// `P(S = 0)`, exact.
    pub p_zero: f64,
    /// Largest number of depositing electrons kept in the convolution.
    pub k_max: usize,
    /// Poisson mass beyond `k_max` that was dropped.
    pub truncated_mass: f64,
    /// Bin width of `X`, eV (the grid of the representation).
    pub bin_width_ev: f64,
    /// Fraction of histories whose `X` was clamped into the last bin.
    pub overflow_fraction: f64,
    /// Quantiles and threshold probabilities, if computed.
    pub distribution: Option<CompoundDistribution>,
    /// Why `distribution` is absent, if it is.
    pub note: Option<String>,
}

fn distribution_atom_at_zero(r: &CompoundRequest) -> CompoundDistribution {
    let b = |v: f64| Bracket {
        low: v,
        central: v,
        high: v,
    };
    CompoundDistribution {
        quantiles: r
            .quantiles
            .iter()
            .map(|&p| QuantileResult {
                p,
                energy_ev: b(0.0),
            })
            .collect(),
        below: r
            .below_ev
            .iter()
            .map(|&t| ThresholdResult {
                threshold_ev: t,
                probability: b(if t > 0.0 { 1.0 } else { 0.0 }),
            })
            .collect(),
        above: r
            .above_ev
            .iter()
            .map(|&t| ThresholdResult {
                threshold_ev: t,
                probability: b(if t < 0.0 { 1.0 } else { 0.0 }),
            })
            .collect(),
    }
}

/// Poisson pmf `P(k)` for `k = 0..=k_max`, where `k_max` is the smallest
/// count whose upper tail is below `tol`; returns also the dropped mass.
fn poisson_pmf(mu: f64, tol: f64) -> (Vec<f64>, usize, f64) {
    let mut p = vec![(-mu).exp()];
    let mut cum = p[0];
    let mut k = 0usize;
    // exp(-mu) underflows for very large mu; the work check rejects those.
    while 1.0 - cum > tol && k < 1_000_000 {
        k += 1;
        let next = p[k - 1] * mu / k as f64;
        p.push(next);
        cum += next;
        if next == 0.0 && k as f64 > mu {
            break;
        }
    }
    (p, k, (1.0 - cum).max(0.0))
}

/// Quantiles and thresholds from the three half-unit-indexed pmfs.
fn summarise(
    w: f64,
    lo: &[f64],
    mid: &[f64],
    hi: &[f64],
    r: &CompoundRequest,
) -> CompoundDistribution {
    // Entry i of each array is the mass at i * (w / 2).
    let half = w / 2.0;
    let cdf = |a: &[f64]| -> Vec<f64> {
        let mut c = 0.0;
        a.iter()
            .map(|&m| {
                c += m;
                c
            })
            .collect()
    };
    let (cl, cm, ch) = (cdf(lo), cdf(mid), cdf(hi));
    let quantile = |c: &[f64], p: f64| -> f64 {
        let i = c.iter().position(|&v| v >= p).unwrap_or(c.len() - 1);
        i as f64 * half
    };
    // P(S <= t) and P(S < t): mass strictly below t is the cdf at the last
    // grid point < t.
    let below = |c: &[f64], t: f64| -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        // Grid points i*half < t  <=>  i < t/half.
        let n = (t / half).ceil() as usize; // count of grid points strictly below t
        if n == 0 {
            0.0
        } else {
            c[(n - 1).min(c.len() - 1)]
        }
    };
    let le = |c: &[f64], t: f64| -> f64 {
        if t < 0.0 {
            return 0.0;
        }
        let i = (t / half).floor() as usize;
        c[i.min(c.len() - 1)]
    };
    CompoundDistribution {
        quantiles: r
            .quantiles
            .iter()
            .map(|&p| QuantileResult {
                p,
                energy_ev: Bracket {
                    low: quantile(&cl, p),
                    central: quantile(&cm, p),
                    high: quantile(&ch, p),
                },
            })
            .collect(),
        below: r
            .below_ev
            .iter()
            .map(|&t| ThresholdResult {
                threshold_ev: t,
                probability: Bracket {
                    low: below(&cl, t),
                    central: below(&cm, t),
                    high: below(&ch, t),
                },
            })
            .collect(),
        above: r
            .above_ev
            .iter()
            .map(|&t| ThresholdResult {
                threshold_ev: t,
                probability: Bracket {
                    low: 1.0 - le(&cl, t),
                    central: 1.0 - le(&cm, t),
                    high: 1.0 - le(&ch, t),
                },
            })
            .collect(),
    }
}

/// Per-history accumulation of the deposit in each configured volume.
/// Call [`begin_history`](Self::begin_history), then
/// [`deposit`](Self::deposit) for every deposit, then
/// [`end_history`](Self::end_history).
#[derive(Debug, Clone, PartialEq)]
pub struct ExposureAccumulator {
    config: ExposureConfig,
    hists: Vec<ExposureHistogram>,
    /// Energy of the current history per volume.
    scratch: Vec<f64>,
    /// Landing offset and start of the current history.
    shift: [f64; 2],
    origin: [f64; 2],
}

impl ExposureAccumulator {
    /// Empty accumulator.
    pub fn new(config: ExposureConfig) -> Result<Self, ExposureError> {
        config.validate()?;
        let n = config.volumes.len();
        let hists = (0..n)
            .map(|_| ExposureHistogram::new(config.bins, config.x_max_ev))
            .collect();
        Ok(Self {
            config,
            hists,
            scratch: vec![0.0; n],
            shift: [0.0; 2],
            origin: [0.0; 2],
        })
    }

    /// The configuration.
    pub fn config(&self) -> &ExposureConfig {
        &self.config
    }

    /// Starts history `index` of a primary starting at lateral position
    /// `start_yz` (m).
    pub fn begin_history(&mut self, index: u64, start_yz: [f64; 2]) {
        self.shift = self.config.beam.offset(self.config.seed, index);
        self.origin = start_yz;
        self.scratch.iter_mut().for_each(|s| *s = 0.0);
    }

    /// A deposit of `e` eV at `pos` (m) in the current history.
    pub fn deposit(&mut self, pos: [f64; 3], e: f64) {
        if !(e.is_finite() && e > 0.0) {
            return;
        }
        let dy = pos[1] - self.origin[0] + self.shift[0];
        let dz = pos[2] - self.origin[1] + self.shift[1];
        for (v, s) in self.config.volumes.iter().zip(&mut self.scratch) {
            if v.contains(pos[0], dy, dz) {
                *s += e;
            }
        }
    }

    /// Ends the current history.
    pub fn end_history(&mut self) {
        for (h, s) in self.hists.iter_mut().zip(&mut self.scratch) {
            h.push(*s);
            *s = 0.0;
        }
    }

    /// Field-wise merge of another accumulator's totals.
    ///
    /// # Panics
    /// If the configurations differ.
    pub fn merge(&mut self, o: &Self) {
        assert!(
            self.config == o.config,
            "merging exposure accumulators with different configurations"
        );
        for (a, b) in self.hists.iter_mut().zip(&o.hists) {
            a.merge(b);
        }
    }

    /// The per-volume results.
    pub fn report(&self) -> ExposureReport {
        let request = CompoundRequest {
            quantiles: self.config.quantiles.clone(),
            below_ev: self.config.below_ev.clone(),
            above_ev: self.config.above_ev.clone(),
            tail_tolerance: self.config.tail_tolerance,
        };
        let volumes = self
            .config
            .volumes
            .iter()
            .zip(&self.hists)
            .map(|(v, h)| {
                let (compound, error) = match h.compound_poisson(self.config.lambda, &request) {
                    Ok(c) => (Some(c), None),
                    Err(e) => (None, Some(e.to_string())),
                };
                ExposureVolumeReport {
                    volume: v.clone(),
                    mean_ev: h.mean_ev(),
                    mean_std_err_ev: h.mean_std_err_ev(),
                    second_moment_ev2: h.second_moment_ev2(),
                    second_moment_std_err_ev2: h.second_moment_std_err_ev2(),
                    zero_fraction: h.zero_fraction(),
                    histogram: h.clone(),
                    compound,
                    compound_error: error,
                }
            })
            .collect();
        ExposureReport {
            lambda: self.config.lambda,
            beam: self.config.beam,
            volumes,
        }
    }
}

/// Results for one volume.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExposureVolumeReport {
    /// The volume.
    pub volume: ExposureVolume,
    /// `E[X]`, eV.
    pub mean_ev: f64,
    /// Its standard error, eV.
    pub mean_std_err_ev: f64,
    /// `E[X²]`, eV².
    pub second_moment_ev2: f64,
    /// Its standard error, eV².
    pub second_moment_std_err_ev2: f64,
    /// `P(X = 0)`.
    pub zero_fraction: f64,
    /// The per-electron distribution.
    pub histogram: ExposureHistogram,
    /// The compound-Poisson summary at the configured `λ`.
    pub compound: Option<CompoundPoisson>,
    /// Why `compound` is absent, if it is.
    pub compound_error: Option<String>,
}

/// Exposure results of a run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExposureReport {
    /// Mean number of incident electrons.
    pub lambda: f64,
    /// Beam profile the landing offsets were drawn from.
    pub beam: BeamProfile,
    /// Per volume, in configuration order.
    pub volumes: Vec<ExposureVolumeReport>,
}

impl ExposureReport {
    /// CSV of the per-electron histograms: `volume,bin,lo_ev,hi_ev,count,
    /// fraction`, with bin 0 the exact `X = 0` atom (`lo = hi = 0`).
    pub fn histogram_csv(&self) -> String {
        let mut s = String::from("volume,bin,lo_ev,hi_ev,count,fraction\n");
        for v in &self.volumes {
            let h = &v.histogram;
            let n = h.histories.max(1) as f64;
            s += &format!(
                "{},0,0,0,{},{}\n",
                v.volume.name,
                h.zero_count,
                h.zero_count as f64 / n
            );
            for (i, &c) in h.counts.iter().enumerate() {
                s += &format!(
                    "{},{},{},{},{},{}\n",
                    v.volume.name,
                    i + 1,
                    i as f64 * h.bin_width_ev,
                    (i + 1) as f64 * h.bin_width_ev,
                    c,
                    c as f64 / n
                );
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn req() -> CompoundRequest {
        CompoundRequest {
            quantiles: vec![0.5],
            below_ev: vec![],
            above_ev: vec![],
            tail_tolerance: 1e-14,
        }
    }

    /// X = 1 eV with probability p (a single bin of width 1): then
    /// S ~ Poisson(λ p), and the three representations are all exact.
    #[test]
    fn thinned_poisson_matches_the_closed_form() {
        let mut h = ExposureHistogram::new(1, 1.0);
        for i in 0..100 {
            h.push(if i < 30 { 1.0 } else { 0.0 });
        }
        let lambda = 4.0;
        let r = CompoundRequest {
            quantiles: vec![0.5, 0.9],
            below_ev: vec![3.0],
            above_ev: vec![5.0],
            tail_tolerance: 1e-15,
        };
        let c = h.compound_poisson(lambda, &r).unwrap();
        let mu: f64 = lambda * 0.3;
        // Closed form.
        let pmf = |k: u32| (-mu).exp() * mu.powi(k as i32) / (1..=k).product::<u32>() as f64;
        let d = c.distribution.unwrap();
        let below3: f64 = (0..3).map(pmf).sum();
        let above5: f64 = 1.0 - (0..=5).map(pmf).sum::<f64>();
        // Upper edge representation is exact for X = 1 (bin (0,1]).
        assert!((d.below[0].probability.high - below3).abs() < 1e-12);
        assert!((d.above[0].probability.high - above5).abs() < 1e-12);
        assert!((c.p_zero - (-mu).exp()).abs() < 1e-15);
        assert!((c.mean_ev - lambda * 0.3).abs() < 1e-12);
        assert!((c.variance_ev2 - lambda * 0.3).abs() < 1e-12);
        // Median of Poisson(1.2) is 1.
        assert_eq!(d.quantiles[0].energy_ev.high, 1.0);
        assert_eq!(d.quantiles[0].energy_ev.low, 0.0);
    }

    #[test]
    fn empty_volume_is_an_atom_at_zero() {
        let mut h = ExposureHistogram::new(4, 4.0);
        for _ in 0..10 {
            h.push(0.0);
        }
        let mut r = req();
        r.below_ev = vec![1.0];
        r.above_ev = vec![0.0];
        let c = h.compound_poisson(7.0, &r).unwrap();
        assert_eq!(c.mean_ev, 0.0);
        let d = c.distribution.unwrap();
        assert_eq!(d.quantiles[0].energy_ev.central, 0.0);
        assert_eq!(d.below[0].probability.central, 1.0);
        assert_eq!(d.above[0].probability.central, 0.0);
    }

    #[test]
    fn brackets_order_and_moments_are_consistent() {
        let mut h = ExposureHistogram::new(10, 10.0);
        for i in 0..1000u32 {
            h.push(if i % 3 == 0 {
                0.0
            } else {
                (i % 17) as f64 * 0.5
            });
        }
        let mut r = req();
        r.quantiles = vec![0.1, 0.5, 0.9];
        r.above_ev = vec![8.0];
        let c = h.compound_poisson(3.0, &r).unwrap();
        let d = c.distribution.unwrap();
        for q in &d.quantiles {
            let b = q.energy_ev;
            assert!(b.low <= b.central && b.central <= b.high);
        }
        let b = d.above[0].probability;
        assert!(b.low <= b.central && b.central <= b.high);
        // Bin-centre mean of S is close to the exact mean to within K_max bins.
        assert!(c.truncated_mass <= 1e-14);
    }

    #[test]
    fn config_validation_rejects_bad_scalars() {
        let v = ExposureVolume {
            name: "a".into(),
            shape: VolumeShape::Cylinder { radius_m: 1e-8 },
            depth_lo_m: 0.0,
            depth_hi_m: 1e-8,
            offset_y_m: 0.0,
            offset_z_m: 0.0,
        };
        let ok = ExposureConfig {
            volumes: vec![v],
            beam: BeamProfile::Point,
            seed: 1,
            bins: 10,
            x_max_ev: 100.0,
            lambda: 1.0,
            quantiles: vec![0.5],
            below_ev: vec![],
            above_ev: vec![],
            tail_tolerance: 1e-12,
        };
        assert!(ok.validate().is_ok());
        let mut b = ok.clone();
        b.lambda = 0.0;
        assert!(b.validate().is_err());
        let mut b = ok.clone();
        b.volumes.clear();
        assert!(b.validate().is_err());
        let mut b = ok;
        b.quantiles = vec![1.0];
        assert!(b.validate().is_err());
    }

    #[test]
    fn beam_offsets_are_deterministic_and_on_the_disc() {
        let b = BeamProfile::UniformDisc { radius_m: 2.0 };
        assert_eq!(b.offset(5, 9), b.offset(5, 9));
        for i in 0..200 {
            let [y, z] = b.offset(5, i);
            assert!(y.hypot(z) <= 2.0);
        }
        assert_eq!(BeamProfile::Point.offset(1, 1), [0.0, 0.0]);
    }
}
