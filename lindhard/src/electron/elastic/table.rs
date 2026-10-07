//! Energy-grid elastic tables: the radial-Dirac solver of the parent module
//! run over a log energy grid for every element of a [`Material`], stored as
//! an [`electron::data::CrossSectionTable`](crate::electron::data::CrossSectionTable)
//! with [`SamplingAxis::ElasticPolarAngle`] (issue #90, step 2 of the Mott
//! elastic model).
//!
//! # Independent-atom approximation
//!
//! A compound or mixture is treated as a set of **independent atoms**: the
//! macroscopic differential cross section is the number-density-weighted sum
//! of the atomic ones, `n dσ/dΩ = Σ_z n_z dσ_z/dΩ`. Chemical bonding, the
//! solid-state potential and coherent (diffraction) effects between atoms are
//! ignored. This is an approximation; condensed-phase (muffin-tin) potentials
//! are follow-up work. Consequently
//!
//! * the inverse mean free path is `λ⁻¹(E) = Σ_z n_z σ_el,z(E)`, with `n_z`
//!   from [`Material::number_density_of`] and `σ_el,z` the solver's total
//!   elastic cross section ([`PartialWaves::sigma_el`]);
//! * the forward cumulative distribution of the polar angle is
//!   `F(θ) = Σ_z n_z S_z(θ) / Σ_z n_z S_z(π)`, with
//!   `S_z(θ) = 2π ∫_0^θ dσ_z/dΩ sin θ' dθ'`.
//!
//! The stored quantiles are the inverse of that `F`: they do not add
//! linearly between elements, the forward CDF does.
//!
//! # Angular distribution: exact Legendre form
//!
//! The Mott DCS `|f|² + |g|²` of a phase-shift set cut at `l_max` is a
//! polynomial of degree `2 l_max` in `x = cos θ` (the parent module documents
//! the amplitudes). Its Legendre coefficients
//! `c_L = (2L+1)/2 ∫ DCS(x) P_L(x) dx`, `L = 0..=2 l_max`, are therefore
//! obtained *exactly* (to rounding) with a `2 l_max + 2`-point Gauss-Legendre
//! rule, and the cumulative cross section follows in closed form from the
//! Legendre identity `(2L+1) P_L = P'_{L+1} - P'_{L-1}`. That identity is
//! DLMF 18.9.7 (`(n+λ) C_n^(λ) = λ (C_n^(λ+1) - C_{n-2}^(λ+1))`) at
//! `λ = 1/2`, with `C_n^(1/2) = P_n` (DLMF 18.7.9) and
//! `d/dx C_{n+1}^(1/2) = C_n^(3/2)` (DLMF 18.9.19); NIST Digital Library of
//! Mathematical Functions, <https://dlmf.nist.gov/18.9>, read 2026-10-07. The
//! unit tests also check it numerically against direct quadrature:
//!
//! ```text
//! S(θ) = 2π ∫_x^1 DCS dx' = 2π [ c_0 (1 - x) + Σ_{L>=1} c_L (P_{L-1}(x) - P_{L+1}(x)) / (2L+1) ]
//! ```
//!
//! with the complementary form (integral from `-1` to `x`) used for
//! cumulative probabilities above 1/2, so the backward tail keeps its
//! relative precision. `S(π) = 4π c_0`. No angular grid of the DCS is
//! interpolated anywhere: the only discretisation in a table is the
//! probability grid of the stored inverse CDF.
//!
//! # Inverse CDF
//!
//! At each energy the quantile `θ_j` with `F(θ_j) = u_j` is found by a
//! safeguarded Newton iteration (bisection fallback) on the exact `F`, after
//! bracketing on a fixed coarse angle grid. `u = 0` maps to `θ = 0` and
//! `u = 1` to `θ = π`. The probability grid is shared by all rows and must
//! resolve both tails: the Mott DCS is strongly forward peaked at high energy,
//! so most of the first transport cross section `σ_tr1 = σ_el ⟨1 - cos θ⟩`
//! sits in the large-angle tail at `1 - u` far below 1e-3. The starting grid
//! ([`default_probability_grid`]) is uniform in the logistic variable
//! `t = ln(u/(1-u))`, which is geometric in `u` near 0 and in `1 - u` near 1.
//!
//! A fixed grid of that kind is not enough. At intermediate energies the DCS
//! has deep diffraction minima, where `θ(u)` is steep, and linear
//! interpolation there biases `⟨1 - cos θ⟩`. Measured with the Thomas-Fermi
//! Yukawa stand-in on the default energy grid (Z = 1, 6, 14, 29, 79, 92), the
//! worst relative error of the recovered `σ_tr1` was 2e-2 for a logistic
//! step of 0.4, 3e-3 for a step of 0.1 and still 1.3e-3 for a step of 0.05
//! (1109 points). The grid is therefore **refined adaptively**: shared by all
//! rows of a table and split where any row needs it (see `refine` in the
//! source). With the defaults (step 0.4, tolerance
//! [`DEFAULT_REFINE_TOLERANCE`]) the same measurement gives a worst error of
//! 9.4e-5 or less, with 939 to 1604 probability points. The tests hold the
//! recovered `σ_tr1` to 1e-3 at every grid energy.
//!
//! # Potentials
//!
//! The atomic potential comes from a caller-supplied [`PotentialSource`]. The
//! Salvat et al. (1987) DHFS table is not in this tree (see the parent module
//! and `docs/data-provenance.md`), so [`SalvatDhfsTable`] always fails; the
//! tests and benchmarks use [`ThomasFermiYukawa`], a **stand-in** Yukawa
//! potential with the Thomas-Fermi length. The table's `model` and
//! `provenance` strings name the potential, so a stand-in table never reads
//! as DHFS.
//!
//! # Determinism
//!
//! Energies are solved in parallel with rayon and collected in grid order;
//! every per-energy computation is serial and in a fixed order (the solver's
//! own partial-wave parallelism is already order-preserving). Errors are
//! reported for the lowest failing energy, whatever the thread count. Tables
//! are bit-identical on any number of threads (tested).

use super::{
    gauss_legendre, ElasticError, ElasticSolver, PartialWaves, SalvatDhfs, ScreenedPotential,
    SolverOptions, Yukawa, BOHR2_TO_M2,
};
use crate::electron::data::{
    CrossSectionTable, CrossSectionTableParts, ElectronDataError, SamplingAxis,
};
use crate::elements::element;
use crate::material::{Material, MaterialError};
use rayon::prelude::*;
use std::f64::consts::PI;

/// Lowest energy of the default grid, eV.
pub const DEFAULT_MIN_ENERGY_EV: f64 = 10.0;
/// Highest energy of the default grid, eV.
pub const DEFAULT_MAX_ENERGY_EV: f64 = 5.0e4;
/// Minimum number of energy points per decade of the default grid.
pub const DEFAULT_POINTS_PER_DECADE: f64 = 20.0;
/// Default probability tails: the first and last interior probability points
/// are `DEFAULT_PROBABILITY_TAIL` and `1 - DEFAULT_PROBABILITY_TAIL`.
pub const DEFAULT_PROBABILITY_TAIL: f64 = 1.0e-12;
/// Default step of the probability grid in the logistic variable
/// `t = ln(u/(1-u))`. Chosen by measurement, see the module docs.
pub const DEFAULT_LOGISTIC_STEP: f64 = 0.4;

/// Default target of the probability-grid refinement, relative to
/// `⟨1 - cos θ⟩` (see [`ElasticTableOptions::refine_tolerance`]).
pub const DEFAULT_REFINE_TOLERANCE: f64 = 5.0e-4;

/// Errors from the table builder.
#[derive(Debug, thiserror::Error)]
pub enum ElasticTableError {
    /// An input violates an invariant.
    #[error("invalid {what}: {reason}")]
    Invalid {
        /// Which quantity.
        what: &'static str,
        /// What is wrong with it.
        reason: String,
    },
    /// The potential source has no potential for this element.
    #[error("no screened potential for Z={z}: {source}")]
    Potential {
        /// Atomic number.
        z: u8,
        /// The source's error.
        source: ElasticError,
    },
    /// The partial-wave solver failed.
    #[error("elastic solve failed for Z={z} at {energy_ev} eV: {source}")]
    Solve {
        /// Atomic number.
        z: u8,
        /// Kinetic energy, eV.
        energy_ev: f64,
        /// The solver's error.
        source: ElasticError,
    },
    /// A material lookup failed.
    #[error(transparent)]
    Material(#[from] MaterialError),
    /// The assembled table failed validation.
    #[error(transparent)]
    Data(#[from] ElectronDataError),
}

fn invalid(what: &'static str, reason: impl Into<String>) -> ElasticTableError {
    ElasticTableError::Invalid {
        what,
        reason: reason.into(),
    }
}

// ---------------------------------------------------------------------------
// Potential sources

/// Supplies the atomic screened potential of each element.
pub trait PotentialSource: Sync {
    /// The potential of element `z`.
    fn potential(&self, z: u8) -> Result<Box<dyn ScreenedPotential>, ElasticError>;

    /// One line naming the potential model and its origin. It is written into
    /// the table's `model` and `provenance` strings.
    fn description(&self) -> String;
}

/// **Stand-in** potential: a Yukawa (single exponential) screened Coulomb
/// potential `V = -(Z/r) exp(-r/a)` with the Thomas-Fermi length
/// `a = 0.8853 a0 Z^(-1/3)`
/// ([`thomas_fermi_constant`](crate::ion::potential::thomas_fermi_constant),
/// Firsov 1958 / Lindhard, Scharff & Schiott 1963, see
/// `docs/data-provenance.md`). It is not an atomic DHFS potential; it lets the
/// table machinery be built and tested while the Salvat et al. (1987)
/// coefficients are missing.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ThomasFermiYukawa;

impl PotentialSource for ThomasFermiYukawa {
    fn potential(&self, z: u8) -> Result<Box<dyn ScreenedPotential>, ElasticError> {
        let zf = f64::from(z);
        if z == 0 {
            return Err(ElasticError::Invalid {
                what: "atomic number",
                reason: "0".into(),
            });
        }
        let a = crate::ion::potential::thomas_fermi_constant() / zf.cbrt();
        Ok(Box::new(Yukawa::new(zf, a)?))
    }

    fn description(&self) -> String {
        "STAND-IN Yukawa potential -(Z/r) exp(-r/a) with Thomas-Fermi length \
         a = 0.8853 a0 Z^(-1/3); not a DHFS atomic potential"
            .into()
    }
}

/// The Salvat et al. (1987) DHFS potentials via [`SalvatDhfs::for_element`].
/// **Every call fails in this build**: the coefficient table is a documented
/// gap (`docs/data-provenance.md`).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SalvatDhfsTable;

impl PotentialSource for SalvatDhfsTable {
    fn potential(&self, z: u8) -> Result<Box<dyn ScreenedPotential>, ElasticError> {
        Ok(Box::new(SalvatDhfs::for_element(u32::from(z))?))
    }

    fn description(&self) -> String {
        "Salvat, Martinez, Mayol & Parellada, Phys. Rev. A 36, 467 (1987) DHFS analytic \
         screening"
            .into()
    }
}

// ---------------------------------------------------------------------------
// Grids

/// A log-spaced energy grid from `min_ev` to `max_ev` (both included exactly)
/// with at least `points_per_decade` points per decade.
pub fn log_energy_grid(
    min_ev: f64,
    max_ev: f64,
    points_per_decade: f64,
) -> Result<Vec<f64>, ElasticTableError> {
    if !(min_ev.is_finite() && max_ev.is_finite() && min_ev > 0.0 && max_ev > min_ev) {
        return Err(invalid(
            "energy range",
            format!("{min_ev} .. {max_ev} eV (need 0 < min < max, finite)"),
        ));
    }
    if !(points_per_decade.is_finite() && points_per_decade > 0.0) {
        return Err(invalid(
            "points per decade",
            format!("{points_per_decade} (need finite, > 0)"),
        ));
    }
    let decades = (max_ev / min_ev).log10();
    let intervals = ((decades * points_per_decade).ceil() as usize).max(1);
    let step = decades / intervals as f64;
    let mut grid: Vec<f64> = (0..=intervals)
        .map(|i| min_ev * 10f64.powf(step * i as f64))
        .collect();
    grid[0] = min_ev;
    grid[intervals] = max_ev;
    Ok(grid)
}

/// The default energy grid: [`DEFAULT_MIN_ENERGY_EV`] to
/// [`DEFAULT_MAX_ENERGY_EV`] with at least [`DEFAULT_POINTS_PER_DECADE`]
/// points per decade.
pub fn default_energy_grid() -> Vec<f64> {
    log_energy_grid(
        DEFAULT_MIN_ENERGY_EV,
        DEFAULT_MAX_ENERGY_EV,
        DEFAULT_POINTS_PER_DECADE,
    )
    .expect("valid default grid")
}

/// A probability grid uniform in the logistic variable `t = ln(u/(1-u))`
/// between `u = tail` and `u = 1 - tail` with step at most `step`, plus the
/// end points 0 and 1.
pub fn logistic_probability_grid(tail: f64, step: f64) -> Result<Vec<f64>, ElasticTableError> {
    if !(tail.is_finite() && tail > 1.0e-15 && tail < 0.5) {
        return Err(invalid(
            "probability tail",
            format!("{tail} (need 1e-15 < tail < 0.5)"),
        ));
    }
    if !(step.is_finite() && step >= 1.0e-3) {
        return Err(invalid(
            "logistic step",
            format!("{step} (need finite, >= 1e-3)"),
        ));
    }
    let t_max = ((1.0 - tail) / tail).ln();
    let n = ((2.0 * t_max / step).ceil() as usize).max(1);
    let h = 2.0 * t_max / n as f64;
    let mut grid = Vec::with_capacity(n + 3);
    grid.push(0.0);
    for i in 0..=n {
        let t = -t_max + h * i as f64;
        // Each branch keeps the small side (u or 1-u) at full relative precision.
        let u = if t <= 0.0 {
            let e = t.exp();
            e / (1.0 + e)
        } else {
            1.0 - 1.0 / (1.0 + t.exp())
        };
        if u > *grid.last().expect("non-empty") && u < 1.0 {
            grid.push(u);
        }
    }
    grid.push(1.0);
    Ok(grid)
}

/// The default probability grid: [`logistic_probability_grid`] with
/// [`DEFAULT_PROBABILITY_TAIL`] and [`DEFAULT_LOGISTIC_STEP`].
pub fn default_probability_grid() -> Vec<f64> {
    logistic_probability_grid(DEFAULT_PROBABILITY_TAIL, DEFAULT_LOGISTIC_STEP)
        .expect("valid default grid")
}

/// Grids and solver settings of a table build.
#[derive(Debug, Clone, PartialEq)]
pub struct ElasticTableOptions {
    /// Incident kinetic energies, eV (strictly increasing, >= 2 points).
    pub energy_ev: Vec<f64>,
    /// Cumulative-probability grid of the stored inverse CDFs: the starting
    /// grid when `refine_tolerance` is set, the final grid otherwise.
    pub probability: Vec<f64>,
    /// Target of the adaptive refinement of the probability grid, as a
    /// fraction of `⟨1 - cos θ⟩` at each energy (`None`: no refinement).
    /// The grid is shared by all rows of one table.
    pub refine_tolerance: Option<f64>,
    /// Radial-Dirac solver settings, used unchanged at every energy.
    pub solver: SolverOptions,
}

impl Default for ElasticTableOptions {
    fn default() -> Self {
        Self {
            energy_ev: default_energy_grid(),
            probability: default_probability_grid(),
            refine_tolerance: Some(DEFAULT_REFINE_TOLERANCE),
            solver: SolverOptions::default(),
        }
    }
}

// ---------------------------------------------------------------------------
// One element at one energy

/// Partial-wave results of one element at one energy, with the DCS in exact
/// Legendre form (see the module docs). Cross sections in bohr².
#[derive(Debug, Clone, PartialEq)]
pub struct AtomicElasticRow {
    energy_ev: f64,
    l_max: usize,
    sigma_el: f64,
    sigma_tr1: f64,
    legendre: Vec<f64>,
}

/// Values of the Legendre series at one angle: `(DCS, S_lower, S_upper)`
/// without the `2π` factors (`S_lower = ∫_x^1 DCS`, `S_upper = ∫_{-1}^x DCS`).
fn legendre_sums(c: &[f64], theta: f64) -> (f64, f64, f64) {
    let x = theta.cos();
    // 1 - x and 1 + x at full relative precision near either end
    let (sh, ch) = (0.5 * theta).sin_cos();
    let (omx, opx) = (2.0 * sh * sh, 2.0 * ch * ch);
    let mut dcs = c[0];
    let mut lower = c[0] * omx;
    let mut upper = c[0] * opx;
    let (mut p_prev, mut p) = (1.0_f64, x);
    for (l, &cl) in c.iter().enumerate().skip(1) {
        let lf = l as f64;
        let p_next = ((2.0 * lf + 1.0) * x * p - lf * p_prev) / (lf + 1.0);
        dcs += cl * p;
        let term = cl * (p_prev - p_next) / (2.0 * lf + 1.0);
        lower += term;
        upper -= term;
        p_prev = p;
        p = p_next;
    }
    (dcs, lower, upper)
}

impl AtomicElasticRow {
    /// Build from converged partial waves at `energy_ev`.
    pub fn from_partial_waves(energy_ev: f64, pw: &PartialWaves) -> Self {
        let l_max = pw.l_max();
        let deg = 2 * l_max;
        let (xs, ws) = gauss_legendre(deg + 2);
        let mut c = vec![0.0; deg + 1];
        for (&x, &w) in xs.iter().zip(&ws) {
            let (f, g) = pw.amplitudes(x);
            let d = w * (f[0] * f[0] + f[1] * f[1] + g[0] * g[0] + g[1] * g[1]);
            let (mut p_prev, mut p) = (0.0_f64, 1.0_f64);
            for (l, cl) in c.iter_mut().enumerate() {
                *cl += d * p;
                let lf = l as f64;
                let p_next = ((2.0 * lf + 1.0) * x * p - lf * p_prev) / (lf + 1.0);
                p_prev = p;
                p = p_next;
            }
        }
        for (l, cl) in c.iter_mut().enumerate() {
            *cl *= (2 * l + 1) as f64 / 2.0;
        }
        Self {
            energy_ev,
            l_max,
            sigma_el: pw.sigma_el(),
            sigma_tr1: pw.sigma_tr1(),
            legendre: c,
        }
    }

    /// Solve at one energy.
    pub fn solve(
        pot: &dyn ScreenedPotential,
        energy_ev: f64,
        opts: SolverOptions,
    ) -> Result<Self, ElasticError> {
        let pw = ElasticSolver::new(pot, energy_ev, opts)?.partial_waves()?;
        Ok(Self::from_partial_waves(energy_ev, &pw))
    }

    /// Kinetic energy, eV.
    pub fn energy_ev(&self) -> f64 {
        self.energy_ev
    }

    /// Highest partial wave of the solve.
    pub fn l_max(&self) -> usize {
        self.l_max
    }

    /// The solver's total elastic cross section ([`PartialWaves::sigma_el`]),
    /// bohr².
    pub fn sigma_el(&self) -> f64 {
        self.sigma_el
    }

    /// The solver's first transport cross section
    /// ([`PartialWaves::sigma_tr1`]), bohr².
    pub fn sigma_tr1(&self) -> f64 {
        self.sigma_tr1
    }

    /// Legendre coefficients `c_L` of the DCS, bohr²/sr.
    pub fn legendre_coefficients(&self) -> &[f64] {
        &self.legendre
    }

    /// DCS at polar angle `theta` (rad) from the Legendre form, bohr²/sr.
    pub fn dcs(&self, theta: f64) -> f64 {
        legendre_sums(&self.legendre, theta).0
    }

    /// Cumulative cross section `S(θ) = 2π ∫_0^θ DCS sin θ' dθ'`, bohr².
    pub fn cumulative(&self, theta: f64) -> f64 {
        2.0 * PI * legendre_sums(&self.legendre, theta).1
    }

    /// `S(π) = 4π c_0`, the total of the Legendre form, bohr². Equal to
    /// [`Self::sigma_el`] up to the solver's quadrature-vs-phase-shift
    /// rounding.
    pub fn cumulative_total(&self) -> f64 {
        4.0 * PI * self.legendre[0]
    }
}

/// Partial-wave results of one element over an energy grid.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomicElastic {
    z: u8,
    potential: String,
    solver: SolverOptions,
    rows: Vec<AtomicElasticRow>,
}

impl AtomicElastic {
    /// Solve element `z` with potential `pot` (described by `description`) at
    /// every energy of `energy_ev`. Energies run in parallel; the result and
    /// any error (the one at the lowest failing energy) do not depend on the
    /// thread count.
    pub fn compute(
        z: u8,
        pot: &dyn ScreenedPotential,
        description: &str,
        energy_ev: &[f64],
        solver: SolverOptions,
    ) -> Result<Self, ElasticTableError> {
        check_grid(energy_ev)?;
        let results: Vec<Result<AtomicElasticRow, ElasticError>> = energy_ev
            .par_iter()
            .map(|&e| AtomicElasticRow::solve(pot, e, solver))
            .collect();
        let mut rows = Vec::with_capacity(results.len());
        for (r, &e) in results.into_iter().zip(energy_ev) {
            rows.push(r.map_err(|source| ElasticTableError::Solve {
                z,
                energy_ev: e,
                source,
            })?);
        }
        Ok(Self {
            z,
            potential: description.to_string(),
            solver,
            rows,
        })
    }

    /// Atomic number.
    pub fn z(&self) -> u8 {
        self.z
    }

    /// Description of the potential used.
    pub fn potential(&self) -> &str {
        &self.potential
    }

    /// Solver settings used.
    pub fn solver(&self) -> SolverOptions {
        self.solver
    }

    /// One row per grid energy.
    pub fn rows(&self) -> &[AtomicElasticRow] {
        &self.rows
    }
}

fn check_grid(energy_ev: &[f64]) -> Result<(), ElasticTableError> {
    if energy_ev.len() < 2 {
        return Err(invalid("energy grid", "needs at least 2 points"));
    }
    if energy_ev.iter().any(|e| !(e.is_finite() && *e > 0.0))
        || energy_ev.windows(2).any(|w| w[1] <= w[0])
    {
        return Err(invalid(
            "energy grid",
            "values must be finite, positive and strictly increasing",
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Mixture inverse CDF

/// Angles of the coarse bracketing table: 0, then geometric from 1e-9 rad to π.
fn bracket_angles() -> Vec<f64> {
    const N: usize = 400;
    let (lo, hi) = (1.0e-9_f64, PI);
    let r = (hi / lo).ln() / (N - 1) as f64;
    let mut v = Vec::with_capacity(N + 1);
    v.push(0.0);
    for i in 0..N {
        v.push(lo * (r * i as f64).exp());
    }
    v[N] = PI;
    v
}

/// The weighted sum of atomic Legendre series at one energy.
struct MixtureRow<'a> {
    parts: Vec<(f64, &'a [f64])>,
    total: f64,
}

impl MixtureRow<'_> {
    /// `(dW/dθ, W_lower, W_upper)` with `W = Σ n_z S_z` (without the common
    /// `2π`, which cancels in the CDF).
    fn eval(&self, theta: f64) -> (f64, f64, f64) {
        let (mut d, mut lo, mut up) = (0.0, 0.0, 0.0);
        for (n, c) in &self.parts {
            let (a, b, cc) = legendre_sums(c, theta);
            d += n * a;
            lo += n * b;
            up += n * cc;
        }
        (d * theta.sin(), lo, up)
    }

    /// Exact `⟨1 - cos θ⟩ = 1 - Σ n c_1 / (3 Σ n c_0)` of the mixture
    /// (`∫ x DCS dx = 2 c_1 / 3`).
    fn mean_one_minus_cos(&self) -> f64 {
        let c1: f64 = self
            .parts
            .iter()
            .map(|(n, c)| n * c.get(1).copied().unwrap_or(0.0))
            .sum();
        1.0 - (2.0 / 3.0) * c1 / self.total
    }

    /// Residual increasing in θ whose root is the quantile of `u`, and its
    /// derivative.
    fn residual(&self, u: f64, theta: f64) -> (f64, f64) {
        let (d, lo, up) = self.eval(theta);
        if u <= 0.5 {
            (lo - u * self.total, d)
        } else {
            ((1.0 - u) * self.total - up, d)
        }
    }

    fn quantiles(&self, probability: &[f64], angles: &[f64]) -> Vec<f64> {
        let m = probability.len();
        let coarse: Vec<(f64, f64)> = angles
            .iter()
            .map(|&t| {
                let (_, lo, up) = self.eval(t);
                (lo, up)
            })
            .collect();
        let mut out = Vec::with_capacity(m);
        let mut prev = 0.0_f64;
        for (j, &u) in probability.iter().enumerate() {
            let theta = if j == 0 {
                0.0
            } else if j == m - 1 {
                PI
            } else {
                // first coarse node whose residual is positive
                let k = if u <= 0.5 {
                    let target = u * self.total;
                    coarse.partition_point(|c| c.0 <= target)
                } else {
                    let target = (1.0 - u) * self.total;
                    coarse.partition_point(|c| c.1 >= target)
                };
                let hi = angles[k.clamp(1, angles.len() - 1)];
                let lo = angles[k.saturating_sub(1)].max(prev).min(hi);
                self.solve(u, lo, hi).max(prev)
            };
            prev = theta.clamp(0.0, PI);
            out.push(prev);
        }
        out
    }

    /// Safeguarded Newton on `[lo, hi]`.
    fn solve(&self, u: f64, mut lo: f64, mut hi: f64) -> f64 {
        let (flo, _) = self.residual(u, lo);
        let (fhi, _) = self.residual(u, hi);
        if flo > 0.0 || fhi < 0.0 {
            // rounding made the coarse bracket fail; fall back to the widest
            if flo > 0.0 {
                lo = 0.0;
            }
            if fhi < 0.0 {
                hi = PI;
            }
        }
        let mut x = 0.5 * (lo + hi);
        for _ in 0..200 {
            let (f, df) = self.residual(u, x);
            if f == 0.0 {
                return x;
            }
            if f < 0.0 {
                lo = x;
            } else {
                hi = x;
            }
            if hi - lo <= 4.0 * f64::EPSILON * hi {
                break;
            }
            let newton = x - f / df;
            x = if df > 0.0 && newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
            if x == lo || x == hi {
                break;
            }
        }
        x
    }
}

/// `⟨1 - cos θ⟩` of a stored inverse-CDF row, integrated **exactly** over the
/// piecewise-linear `θ(u)` that
/// [`CrossSectionTable::inverse_cdf`](crate::electron::data::CrossSectionTable::inverse_cdf)
/// interpolates. `σ_tr1 = σ_el ⟨1 - cos θ⟩`.
pub fn mean_one_minus_cos(probability: &[f64], quantiles: &[f64]) -> f64 {
    let mut s = 0.0;
    for (p, q) in probability.windows(2).zip(quantiles.windows(2)) {
        let du = p[1] - p[0];
        let m = 0.5 * (q[0] + q[1]);
        let h = 0.5 * (q[1] - q[0]);
        // (1/Δθ) ∫ (1 - cos θ) dθ = 1 - cos(m) sinc(h)
        //                        = 2 sin²(m/2) + cos(m) (1 - sinc h)
        let one_minus_sinc = if h.abs() < 1.0e-3 {
            let h2 = h * h;
            h2 / 6.0 * (1.0 - h2 / 20.0)
        } else {
            1.0 - h.sin() / h
        };
        s += du * (2.0 * (0.5 * m).sin().powi(2) + m.cos() * one_minus_sinc);
    }
    s
}

/// A one-line identity of the material: its name, if set, and its atom
/// fractions and density.
fn material_identity(material: &Material) -> String {
    let comp: Vec<String> = material
        .components()
        .iter()
        .map(|c| {
            let sym = element(c.z()).map_or("?", |e| e.symbol);
            format!("{sym}:{}", c.atom_fraction())
        })
        .collect();
    let base = format!(
        "{} at {} kg/m^3 (atom fractions)",
        comp.join(" "),
        material.mass_density()
    );
    match material.name() {
        Some(n) => format!("{n}: {base}"),
        None => base,
    }
}

/// Combine atomic results into the elastic table of `material` by
/// independent-atom additivity (module docs). `atoms` must hold every
/// component of `material` with a non-zero atom fraction, all on the same
/// energy grid, potential and solver settings. `probability` is the
/// probability grid, used as given when `refine_tolerance` is `None` and as
/// the starting grid of the adaptive refinement otherwise (see
/// [`ElasticTableOptions::refine_tolerance`]).
pub fn combine(
    material: &Material,
    atoms: &[AtomicElastic],
    probability: &[f64],
    refine_tolerance: Option<f64>,
) -> Result<CrossSectionTable, ElasticTableError> {
    let used: Vec<_> = material
        .components()
        .iter()
        .filter(|c| c.atom_fraction() > 0.0)
        .collect();
    let mut parts = Vec::with_capacity(used.len());
    for c in &used {
        let atom = atoms
            .iter()
            .find(|a| a.z == c.z())
            .ok_or_else(|| invalid("atomic results", format!("missing Z={}", c.z())))?;
        parts.push((material.number_density_of(c.z())?, atom));
    }
    let first = parts
        .first()
        .ok_or_else(|| invalid("material", "no component with a non-zero fraction"))?
        .1;
    let energy: Vec<f64> = first.rows.iter().map(|r| r.energy_ev).collect();
    for (_, a) in &parts {
        let same_grid = a.rows.len() == energy.len()
            && a.rows.iter().zip(&energy).all(|(r, e)| r.energy_ev == *e);
        if !same_grid || a.potential != first.potential || a.solver != first.solver {
            return Err(invalid(
                "atomic results",
                "all elements must share the energy grid, potential and solver settings",
            ));
        }
    }
    check_probability(probability)?;
    if let Some(tol) = refine_tolerance {
        if !(tol.is_finite() && tol > 0.0) {
            return Err(invalid(
                "refine tolerance",
                format!("{tol} (need finite, > 0)"),
            ));
        }
    }
    let inverse_mfp_per_m: Vec<f64> = (0..energy.len())
        .map(|i| {
            parts
                .iter()
                .map(|(n, a)| n * a.rows[i].sigma_el * BOHR2_TO_M2)
                .sum()
        })
        .collect();
    let mixes: Vec<MixtureRow> = (0..energy.len())
        .map(|i| MixtureRow {
            parts: parts
                .iter()
                .map(|(n, a)| (*n, a.rows[i].legendre.as_slice()))
                .collect(),
            total: parts
                .iter()
                .map(|(n, a)| n * 2.0 * a.rows[i].legendre[0])
                .sum(),
        })
        .collect();
    let angles = bracket_angles();
    let quantiles: Vec<Vec<f64>> = mixes
        .par_iter()
        .map(|m| m.quantiles(probability, &angles))
        .collect();
    let (probability, quantiles) = match refine_tolerance {
        Some(tol) => refine(&mixes, probability.to_vec(), quantiles, tol),
        None => (probability.to_vec(), quantiles),
    };
    let s = &first.solver;
    let model = format!(
        "lindhard {} electron::elastic::table: Mott DCS from radial-Dirac partial waves, \
         independent-atom additivity; potential: {}",
        env!("CARGO_PKG_VERSION"),
        first.potential
    );
    let provenance = format!(
        "computed by lindhard {} from published formulas (see electron::elastic docs); \
         potential: {}; solver: matching_threshold={} phase_tolerance={} consecutive={} \
         step_scale={} max_l={} light_speed={}; {} energies {}..{} eV; {} probability points",
        env!("CARGO_PKG_VERSION"),
        first.potential,
        s.matching_threshold,
        s.phase_tolerance,
        s.consecutive,
        s.step_scale,
        s.max_l,
        s.light_speed,
        energy.len(),
        energy[0],
        energy[energy.len() - 1],
        probability.len()
    );
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model,
        material: material_identity(material),
        provenance,
        axis: SamplingAxis::ElasticPolarAngle,
        energy_ev: energy,
        inverse_mfp_per_m,
        probability,
        quantiles,
    })?)
}

fn check_probability(p: &[f64]) -> Result<(), ElasticTableError> {
    let ok = p.len() >= 2
        && p[0] == 0.0
        && p[p.len() - 1] == 1.0
        && p.iter().all(|x| x.is_finite())
        && p.windows(2).all(|w| w[1] > w[0]);
    if ok {
        Ok(())
    } else {
        Err(invalid(
            "probability grid",
            "need >= 2 finite, strictly increasing points from exactly 0 to exactly 1",
        ))
    }
}

/// `ln(u/(1-u))`, with the small side of `u` and `1-u` at full precision.
fn logit(u: f64) -> f64 {
    if u <= 0.5 {
        (u / (1.0 - u)).ln()
    } else {
        -((1.0 - u) / u).ln()
    }
}

/// Inverse of [`logit`].
fn logistic(t: f64) -> f64 {
    if t <= 0.0 {
        let e = t.exp();
        e / (1.0 + e)
    } else {
        1.0 - 1.0 / (1.0 + t.exp())
    }
}

/// Probability segments below this width are judged by their absolute
/// contribution rather than by the pointwise error (see [`refine`]).
const REFINE_SMALL_SEGMENT: f64 = 1.0e-4;
/// Upper bound on the refined probability grid.
pub const MAX_PROBABILITY_POINTS: usize = 20_000;

/// Adaptive refinement of the shared probability grid.
///
/// Each interior segment `[u_s, u_{s+1}]` (the two end segments, which touch
/// `u = 0` and `u = 1`, are left alone) is tested at its midpoint in the
/// logistic variable: the exact quantile `θ_mid` is solved for at every
/// energy and compared with the linear interpolation the table will use.
/// With `δ = |cos θ_lin - cos θ_mid|` (the error of `1 - cos θ`) and `μ1` the
/// exact `⟨1 - cos θ⟩` of the row, the segment is split when
/// `δ Δu > tol μ1 max(Δu, 1e-4)` at any energy. Large segments are thus held
/// to a pointwise error `tol μ1`; tiny tail segments to an absolute
/// contribution `1e-4 tol μ1` each. The split inserts the already solved
/// `θ_mid` and only the two halves are re-tested. Refinement stops when no
/// segment is split, or at [`MAX_PROBABILITY_POINTS`]. The midpoint estimate
/// is a heuristic: the achieved recovery of `σ_tr1` is what the tests check.
fn refine(
    mixes: &[MixtureRow],
    mut grid: Vec<f64>,
    mut q: Vec<Vec<f64>>,
    tol: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mu1: Vec<f64> = mixes.iter().map(MixtureRow::mean_one_minus_cos).collect();
    // segments to test, by left index
    let mut check: Vec<usize> = (1..grid.len().saturating_sub(2)).collect();
    while !check.is_empty() && grid.len() < MAX_PROBABILITY_POINTS {
        let mids: Vec<f64> = check
            .iter()
            .map(|&s| logistic(0.5 * (logit(grid[s]) + logit(grid[s + 1]))))
            .collect();
        let per_row: Vec<(Vec<f64>, Vec<bool>)> = mixes
            .par_iter()
            .zip(q.par_iter())
            .zip(mu1.par_iter())
            .map(|((m, row), &mu)| {
                let mut th = Vec::with_capacity(check.len());
                let mut flag = Vec::with_capacity(check.len());
                for (&s, &um) in check.iter().zip(&mids) {
                    let (ua, ub) = (grid[s], grid[s + 1]);
                    let (qa, qb) = (row[s], row[s + 1]);
                    let t = m.solve(um, qa, qb).clamp(qa, qb);
                    let lin = qa + (um - ua) / (ub - ua) * (qb - qa);
                    let du = ub - ua;
                    let delta = (lin.cos() - t.cos()).abs();
                    th.push(t);
                    flag.push(delta * du > tol * mu * du.max(REFINE_SMALL_SEGMENT));
                }
                (th, flag)
            })
            .collect();
        let split: Vec<bool> = (0..check.len())
            .map(|k| per_row.iter().any(|r| r.1[k]))
            .collect();
        if !split.iter().any(|b| *b) {
            break;
        }
        // which check entry (if any) splits old segment s
        let mut at: Vec<Option<usize>> = vec![None; grid.len()];
        for (k, &s) in check.iter().enumerate() {
            if split[k] && mids[k] > grid[s] && mids[k] < grid[s + 1] {
                at[s] = Some(k);
            }
        }
        let mut new_grid = Vec::with_capacity(grid.len() + check.len());
        let mut new_q: Vec<Vec<f64>> = q
            .iter()
            .map(|_| Vec::with_capacity(grid.len() + check.len()))
            .collect();
        let mut new_check = Vec::new();
        for s in 0..grid.len() {
            new_grid.push(grid[s]);
            for (nr, r) in new_q.iter_mut().zip(&q) {
                nr.push(r[s]);
            }
            if let Some(k) = at[s] {
                new_check.push(new_grid.len() - 1);
                new_grid.push(mids[k]);
                for (nr, pr) in new_q.iter_mut().zip(&per_row) {
                    nr.push(pr.0[k]);
                }
                new_check.push(new_grid.len() - 1);
            }
        }
        grid = new_grid;
        q = new_q;
        check = new_check;
    }
    (grid, q)
}

/// Build the elastic table of `material` with potentials from `source`.
/// Components with a zero atom fraction are skipped (their potential is not
/// requested).
pub fn build_elastic_table(
    material: &Material,
    source: &dyn PotentialSource,
    opts: &ElasticTableOptions,
) -> Result<CrossSectionTable, ElasticTableError> {
    check_grid(&opts.energy_ev)?;
    let description = source.description();
    let mut atoms = Vec::new();
    for c in material.components() {
        if c.atom_fraction() <= 0.0 {
            continue;
        }
        let z = c.z();
        let pot = source
            .potential(z)
            .map_err(|source| ElasticTableError::Potential { z, source })?;
        atoms.push(AtomicElastic::compute(
            z,
            pot.as_ref(),
            &description,
            &opts.energy_ev,
            opts.solver,
        )?);
    }
    combine(material, &atoms, &opts.probability, opts.refine_tolerance)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legendre_form_reproduces_the_partial_wave_dcs() {
        let pot = Yukawa::new(29.0, 0.8853 / 29f64.cbrt()).unwrap();
        for e in [50.0, 2000.0] {
            let pw = ElasticSolver::new(&pot, e, SolverOptions::default())
                .unwrap()
                .partial_waves()
                .unwrap();
            let row = AtomicElasticRow::from_partial_waves(e, &pw);
            for t in [0.0, 0.01, 0.3, 1.0, 2.0, 3.0, PI] {
                let (a, b) = (row.dcs(t), pw.dcs(t));
                assert!((a - b).abs() <= 1e-10 * pw.dcs(0.0), "{e} {t}: {a} {b}");
            }
            let s = row.cumulative_total();
            assert!(
                (s - pw.sigma_el()).abs() < 1e-10 * s,
                "{s} {}",
                pw.sigma_el()
            );
            assert!((row.cumulative(PI) - s).abs() < 1e-12 * s);
            assert_eq!(row.cumulative(0.0), 0.0);
        }
    }

    #[test]
    fn closed_form_cumulative_matches_direct_quadrature() {
        let pot = Yukawa::new(14.0, 0.8853 / 14f64.cbrt()).unwrap();
        let pw = ElasticSolver::new(&pot, 300.0, SolverOptions::default())
            .unwrap()
            .partial_waves()
            .unwrap();
        let row = AtomicElasticRow::from_partial_waves(300.0, &pw);
        for theta in [0.05, 0.4, 1.3, 2.9] {
            // composite Simpson of 2 pi DCS sin over [0, theta]
            let n = 4000;
            let h = theta / n as f64;
            let f = |t: f64| 2.0 * PI * pw.dcs(t) * t.sin();
            let mut s = f(0.0) + f(theta);
            for k in 1..n {
                s += if k % 2 == 1 { 4.0 } else { 2.0 } * f(k as f64 * h);
            }
            let quad = s * h / 3.0;
            let exact = row.cumulative(theta);
            assert!(
                (exact - quad).abs() < 1e-9 * row.cumulative_total(),
                "{theta}: {exact} {quad}"
            );
        }
    }

    #[test]
    fn one_minus_cos_is_exact_for_linear_rows() {
        // theta(u) = pi u: <1 - cos> = 1 - (sin pi)/pi = 1
        let p: Vec<f64> = (0..=4).map(|i| i as f64 / 4.0).collect();
        let q: Vec<f64> = p.iter().map(|u| PI * u).collect();
        assert!((mean_one_minus_cos(&p, &q) - 1.0).abs() < 1e-15);
        // constant theta = 0.5
        let q = vec![0.5; 5];
        assert!((mean_one_minus_cos(&p, &q) - (1.0 - 0.5f64.cos())).abs() < 1e-15);
    }

    #[test]
    fn probability_grid_shape() {
        let g = default_probability_grid();
        assert_eq!(g[0], 0.0);
        assert_eq!(*g.last().unwrap(), 1.0);
        assert!(g.windows(2).all(|w| w[1] > w[0]));
        assert!((g[1] - DEFAULT_PROBABILITY_TAIL).abs() < 1e-20);
        assert!((1.0 - g[g.len() - 2] - DEFAULT_PROBABILITY_TAIL).abs() < 1e-15);
        assert!(logistic_probability_grid(0.0, 0.1).is_err());
        assert!(logistic_probability_grid(1e-6, 0.0).is_err());
    }
}
