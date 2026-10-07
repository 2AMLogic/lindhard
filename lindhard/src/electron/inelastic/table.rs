//! Inelastic sampling tables from the single-pole Penn DIIMFP (issue #94):
//! an [`electron::data::CrossSectionTable`](crate::electron::data::CrossSectionTable)
//! with [`SamplingAxis::InelasticEnergyLoss`] (the inverse mean free path and
//! the inverse CDF of the energy loss `W` at every energy of a caller-chosen
//! grid), and a sampler of the momentum transfer `q` given `(E, W)`, from
//! which the deflection angle follows.
//!
//! Both are built from [`SinglePolePenn`] only; no data and no equation beyond
//! the model of [`super::penn`] enter here (the elementary kinematics of a
//! binary collision aside, see below).
//!
//! # Energy-loss table
//!
//! At incident energy `E` the loss distribution is the DIIMFP at fixed `E`,
//! `p(E, W)` ([`SinglePolePenn::diimfp_per_m_ev`]), normalised:
//! `F(W) = ∫_0^W p dW' / ∫_0^E p dW'`. The table stores
//!
//! * `λ⁻¹(E)` from [`SinglePolePenn::imfp_and_stopping`] (the closed-form
//!   reduction of #18, not the integral of `p` computed here; the two agree
//!   to the model's tolerance, 1e-5 in the #18 tests), and
//! * the quantiles `W(u_j)` with `F(W(u_j)) = u_j` on a probability grid
//!   shared by all rows.
//!
//! **Density model.** `W` spans decades, so the integrand is handled in
//! `s = ln W`, where `p dW = p W ds`. `g(s) = p W` is sampled on adaptively
//! chosen nodes between the lowest tabulated ELF energy (below which `p = 0`)
//! and `E`, and taken to be *linear between nodes*. Starting from
//! [`INITIAL_LOSS_PANELS`] equal panels, a panel is split at its midpoint
//! while the Simpson-minus-trapezoid estimate of its error in either
//! `∫ g ds` or `∫ W g ds` exceeds `tol · I · h/H` (`I` the whole integral, `h`
//! the panel width, `H` the whole width; `tol` =
//! [`InelasticTableOptions::density_tolerance`]). The quantiles are then the
//! *exact* inverse of this piecewise-quadratic CDF (a quadratic per panel,
//! solved without cancellation), so the density model is the only
//! approximation of the loss distribution and its error is controlled by that
//! tolerance.
//!
//! **Probability grid.** The stored inverse CDF is interpolated linearly in
//! `u` ([`CrossSectionTable::inverse_cdf`]). Its first moment, which is what
//! sets the stopping power `S = λ⁻¹ ⟨W⟩` with `⟨W⟩ = ∫_0^1 W(u) du`, is
//! dominated by the upper tail of `W`, so the starting grid
//! ([`logistic_probability_grid`]) is uniform in `t = ln(u/(1-u))`, and it is
//! **refined adaptively**, shared by all rows: a segment is split at its
//! logistic midpoint when `δ Δu > tol μ max(Δu, 1e-4)` at any energy, with `δ`
//! the difference between the exact quantile there and the linear
//! interpolation and `μ` the row's mean loss (the criterion of the elastic
//! tables, [`crate::electron::elastic::table`]). The test
//! `stopping_power_from_the_loss_cdf_matches_the_model` holds
//! [`stopping_power_ev_per_m`] to 1e-3 of `S(E)` of #18 at every grid
//! energy.
//!
//! Rows where `λ⁻¹ = 0` (no loss allowed) are stored empty. The table covers
//! the direct single-pole model only: no exchange, no inner shells and no
//! surface or band-gap structure (see [`super::penn`]). A model with the
//! exchange correction enabled ([`SinglePolePenn::with_exchange`]) is
//! rejected by [`build_inelastic_table`] and [`MomentumTransferSampler`].
//!
//! # Momentum-transfer sampler
//!
//! Given `(E, W)`, the DIIMFP integrand in `q` is, in `u = ln q`
//! (S2017 eq. (2), as in [`super::penn`]), `Im[-1/ε(q, W)]` on
//! `q- <= q <= q+`, `q± = k ± k'`, `k = sqrt(2T')`, `k' = sqrt(2(T' - W))`
//! (Hartree units, `T' = E + E_F`), cut to the `q` where the plasma frequency
//! of the pole lies in the tabulated ELF range. [`MomentumTransferSampler`] is
//! a **direct sampler built per `(E, W)`**: the integrand is sampled on
//! adaptive nodes (breaks at the kinematic limits and at every ELF knot, 8
//! panels per smooth segment to start, the same refinement rule as above with
//! the zeroth moment only and tolerance [`DEFAULT_Q_TOLERANCE`]), linear in
//! `ln q` between nodes, and inverted exactly. The density used is therefore
//! the model's to within that tolerance, and the histogram test compares it
//! with an independent quadrature of [`SinglePolePenn::loss_function`] at 1e6
//! samples. It is not a table: building one costs the integrand evaluations of
//! one DIIMFP slice, so a transport loop should keep the sampler for repeated
//! draws at one `(E, W)` and expect that cost for each new pair. A tabulated
//! `(E, W, u)` form is left for later work.
//!
//! The deflection angle follows from momentum conservation of the primary,
//! `q² = k² + k'² - 2 k k' cos θ`, i.e.
//! `cos θ = (k² + k'² - q²)/(2 k k')` ([`MomentumTransferSampler::cos_theta`]),
//! clamped to `[-1, 1]`; `q = q-` is `θ = 0` and `q = q+` is `θ = π`. This is
//! the elementary kinematics of the loss `W` (the final wavenumber `k'` fixed
//! by energy conservation) with the model's own `q±`.
//!
//! # Determinism
//!
//! Rows are built in parallel with rayon and collected in grid order; each
//! row, and the probability refinement (a parallel map over rows followed by
//! a serial split decision), is a fixed sequence of floating-point
//! operations, so tables are bit-identical on any thread count (tested).
//! Samplers consume uniforms `u` supplied by the caller (for example from
//! [`crate::rng::stream`]) and keep no state of their own.

use super::penn::{hartree_ev, SinglePolePenn};
use crate::constants::BOHR_RADIUS;
use crate::electron::data::{
    CrossSectionTable, CrossSectionTableParts, ElectronDataError, SamplingAxis,
};
use crate::electron::elastic::table::{
    logistic, logistic_probability_grid, logit, material_identity,
};
use crate::material::Material;
use rayon::prelude::*;

/// Number of equal panels (in `ln W`) the loss density starts from.
pub const INITIAL_LOSS_PANELS: usize = 64;
/// Default tolerance of the loss-density model (see the module docs).
pub const DEFAULT_DENSITY_TOLERANCE: f64 = 1.0e-5;
/// Default target of the probability-grid refinement, relative to the mean
/// loss of each row.
pub const DEFAULT_REFINE_TOLERANCE: f64 = 1.0e-4;
/// Default tolerance of the momentum-transfer density model.
pub const DEFAULT_Q_TOLERANCE: f64 = 1.0e-6;
/// Default probability tail of the starting grid (see
/// [`logistic_probability_grid`]).
pub const DEFAULT_PROBABILITY_TAIL: f64 = 1.0e-12;
/// Default logistic step of the starting probability grid.
pub const DEFAULT_LOGISTIC_STEP: f64 = 0.5;
/// Upper bound on the number of nodes of one density model and on the size of
/// the refined probability grid.
pub const MAX_NODES: usize = 20_000;

/// Errors from the table builder and the samplers.
#[derive(Debug, thiserror::Error)]
pub enum InelasticTableError {
    /// An input violates an invariant.
    #[error("invalid {what}: {reason}")]
    Invalid {
        /// Which quantity.
        what: &'static str,
        /// What is wrong with it.
        reason: String,
    },
    /// The model reports a loss rate but the density has no support.
    #[error(
        "no loss distribution at {energy_ev} eV (inverse mean free path \
         {inverse_mfp_per_m} m^-1, density integral {integral})"
    )]
    EmptyDistribution {
        /// Incident energy, eV.
        energy_ev: f64,
        /// The model's inverse mean free path, m⁻¹.
        inverse_mfp_per_m: f64,
        /// The integral of the sampled density.
        integral: f64,
    },
    /// The model or the assembled table failed.
    #[error(transparent)]
    Data(#[from] ElectronDataError),
}

fn invalid(what: &'static str, reason: impl Into<String>) -> InelasticTableError {
    InelasticTableError::Invalid {
        what,
        reason: reason.into(),
    }
}

/// The table and the q sampler cover the direct (no-exchange) model only:
/// the q sampler integrates the direct integrand and the table provenance
/// does not record an exchange setting. Reject a model with the exchange
/// correction enabled rather than build something inconsistent.
fn check_no_exchange(penn: &SinglePolePenn) -> Result<(), InelasticTableError> {
    match penn.exchange() {
        Some(x) => Err(invalid(
            "Penn model",
            format!(
                "exchange correction (below {} eV) is not supported by the inelastic \
                 table or the momentum-transfer sampler",
                x.applies_below_ev()
            ),
        )),
        None => Ok(()),
    }
}

// ---------------------------------------------------------------------------
// Piecewise-linear density in one variable

/// A density `g(x)` that is linear between nodes, with its cumulative
/// integral; `x` is `ln W` or `ln q`.
#[derive(Debug, Clone, PartialEq)]
struct LinearDensity {
    x: Vec<f64>,
    g: Vec<f64>,
    /// Cumulative integral at each node.
    cum: Vec<f64>,
}

impl LinearDensity {
    fn from_nodes(x: Vec<f64>, g: Vec<f64>) -> Self {
        let mut cum = Vec::with_capacity(x.len());
        let mut acc = 0.0;
        cum.push(0.0);
        for i in 1..x.len() {
            acc += 0.5 * (x[i] - x[i - 1]) * (g[i] + g[i - 1]);
            cum.push(acc);
        }
        Self { x, g, cum }
    }

    fn total(&self) -> f64 {
        *self.cum.last().expect("non-empty")
    }

    /// `∫ e^x g dx` of the piecewise-linear model (closed form per panel).
    fn first_moment(&self) -> f64 {
        (1..self.x.len())
            .map(|i| panel_exp_moment(self.x[i - 1], self.x[i], self.g[i - 1], self.g[i]))
            .sum()
    }

    /// The inverse CDF: the `x` with cumulative probability `u`.
    fn quantile(&self, u: f64) -> f64 {
        let n = self.x.len();
        if u <= 0.0 {
            return self.x[0];
        }
        if u >= 1.0 {
            return self.x[n - 1];
        }
        let target = u * self.total();
        // First node with cum >= target.
        let j = self.cum.partition_point(|&c| c < target).clamp(1, n - 1);
        let (x0, x1) = (self.x[j - 1], self.x[j]);
        let (g0, g1) = (self.g[j - 1], self.g[j]);
        let h = x1 - x0;
        let rem = (target - self.cum[j - 1]).max(0.0);
        // g0 t + (g1 - g0) t² / (2h) = rem, t in [0, h].
        let disc = (g0 * g0 + 2.0 * (g1 - g0) * rem / h).max(0.0);
        let den = g0 + disc.sqrt();
        let t = if den > 0.0 { 2.0 * rem / den } else { 0.0 };
        x0 + t.clamp(0.0, h)
    }
}

/// `∫_{x0}^{x1} e^x g(x) dx` for `g` linear from `g0` to `g1`.
fn panel_exp_moment(x0: f64, x1: f64, g0: f64, g1: f64) -> f64 {
    // ∫ e^x (g0 + m (x - x0)) dx = e^{x0} [ g0 (E - 1) + m (h E - (E - 1)) ]
    // with h = x1 - x0 and E = e^h.
    let h = x1 - x0;
    if h <= 0.0 {
        return 0.0;
    }
    let m = (g1 - g0) / h;
    let em1 = h.exp_m1();
    x0.exp() * (g0 * em1 + m * (h * (em1 + 1.0) - em1))
}

/// Build the density model of `f` by adaptive panel splitting (module docs):
/// `breaks` are fixed nodes, each segment is first cut into `sub` equal
/// panels.
fn adaptive_density(
    f: &(impl Fn(f64) -> f64 + ?Sized),
    breaks: &[f64],
    sub: usize,
    tol: f64,
    first_moment: bool,
) -> LinearDensity {
    let mut x: Vec<f64> = Vec::new();
    for w in breaks.windows(2) {
        for i in 0..sub {
            x.push(w[0] + (w[1] - w[0]) * i as f64 / sub as f64);
        }
    }
    x.push(*breaks.last().expect("two breaks"));
    let mut g: Vec<f64> = x.iter().map(|&xi| f(xi)).collect();
    // Cached midpoint values of the panels, once evaluated.
    let mut mids: Vec<Option<f64>> = vec![None; x.len() - 1];
    let span = x[x.len() - 1] - x[0];
    loop {
        for (i, m) in mids.iter_mut().enumerate() {
            if m.is_none() {
                *m = Some(f(0.5 * (x[i] + x[i + 1])));
            }
        }
        let d = LinearDensity::from_nodes(x.clone(), g.clone());
        let i0 = d.total();
        let i1 = if first_moment { d.first_moment() } else { 0.0 };
        let mut split = vec![false; mids.len()];
        let mut any = false;
        for i in 0..mids.len() {
            let h = x[i + 1] - x[i];
            let gm = mids[i].expect("filled");
            let allowed = h / span * tol;
            let e0 = (2.0 * h / 3.0 * (gm - 0.5 * (g[i] + g[i + 1]))).abs();
            let mut bad = e0 > allowed * i0;
            if first_moment {
                let w = |xx: f64, gg: f64| xx.exp() * gg;
                let xm = 0.5 * (x[i] + x[i + 1]);
                let lin = 0.5 * (w(x[i], g[i]) + w(x[i + 1], g[i + 1]));
                let e1 = (2.0 * h / 3.0 * (w(xm, gm) - lin)).abs();
                bad |= e1 > allowed * i1;
            }
            if bad && x.len() + 1 < MAX_NODES {
                split[i] = true;
                any = true;
            }
        }
        if !any {
            return d;
        }
        let mut nx = Vec::with_capacity(x.len() + 16);
        let mut ng = Vec::with_capacity(x.len() + 16);
        let mut nm = Vec::with_capacity(mids.len() + 16);
        for i in 0..mids.len() {
            nx.push(x[i]);
            ng.push(g[i]);
            if split[i] {
                nx.push(0.5 * (x[i] + x[i + 1]));
                ng.push(mids[i].expect("filled"));
                nm.push(None);
                nm.push(None);
            } else {
                nm.push(mids[i]);
            }
        }
        nx.push(x[x.len() - 1]);
        ng.push(g[g.len() - 1]);
        x = nx;
        g = ng;
        mids = nm;
    }
}

// ---------------------------------------------------------------------------
// Energy-loss table

/// Grids and tolerances of a table build.
#[derive(Debug, Clone, PartialEq)]
pub struct InelasticTableOptions {
    /// Incident kinetic energies above the Fermi level, eV (strictly
    /// increasing, >= 2 points).
    pub energy_ev: Vec<f64>,
    /// Starting cumulative-probability grid of the stored inverse CDFs (from
    /// exactly 0 to exactly 1); the final grid when `refine_tolerance` is
    /// `None`.
    pub probability: Vec<f64>,
    /// Tolerance of the loss-density model ([`DEFAULT_DENSITY_TOLERANCE`]).
    pub density_tolerance: f64,
    /// Target of the probability-grid refinement, relative to the mean loss
    /// of each row (`None`: no refinement).
    pub refine_tolerance: Option<f64>,
}

impl InelasticTableOptions {
    /// The defaults on the given energy grid.
    pub fn new(energy_ev: Vec<f64>) -> Self {
        Self {
            energy_ev,
            probability: default_probability_grid(),
            density_tolerance: DEFAULT_DENSITY_TOLERANCE,
            refine_tolerance: Some(DEFAULT_REFINE_TOLERANCE),
        }
    }
}

/// The default starting probability grid: [`logistic_probability_grid`] with
/// [`DEFAULT_PROBABILITY_TAIL`] and [`DEFAULT_LOGISTIC_STEP`].
pub fn default_probability_grid() -> Vec<f64> {
    logistic_probability_grid(DEFAULT_PROBABILITY_TAIL, DEFAULT_LOGISTIC_STEP)
        .expect("valid default grid")
}

/// The loss density of one energy.
fn loss_density(
    penn: &SinglePolePenn,
    energy_ev: f64,
    inverse_mfp: f64,
    tol: f64,
) -> Result<LinearDensity, InelasticTableError> {
    let w_lo = penn.optical_elf().energy_ev()[0];
    let f = |s: f64| {
        let w = s.exp().min(energy_ev);
        // diimfp_per_m_ev only fails for non-finite input.
        penn.diimfp_per_m_ev(energy_ev, w).unwrap_or(0.0) * w
    };
    let d = adaptive_density(
        &f,
        &[w_lo.ln(), energy_ev.ln()],
        INITIAL_LOSS_PANELS,
        tol,
        true,
    );
    if !(d.total() > 0.0 && d.total().is_finite()) {
        return Err(InelasticTableError::EmptyDistribution {
            energy_ev,
            inverse_mfp_per_m: inverse_mfp,
            integral: d.total(),
        });
    }
    Ok(d)
}

/// The loss `W` at cumulative probability `u` of a row, clamped to
/// `[0, energy_ev]`.
fn loss_quantile(d: &LinearDensity, energy_ev: f64, u: f64) -> f64 {
    d.quantile(u).exp().min(energy_ev)
}

/// `∫_0^1 W(u) du` by the trapezoid rule of the table's linear
/// interpolation.
pub fn mean_loss_ev(probability: &[f64], quantiles: &[f64]) -> f64 {
    probability
        .windows(2)
        .zip(quantiles.windows(2))
        .map(|(p, q)| 0.5 * (p[1] - p[0]) * (q[0] + q[1]))
        .sum()
}

/// The stopping power `λ⁻¹ ⟨W⟩` (eV/m) of row `i` of an energy-loss table,
/// with `⟨W⟩` from [`mean_loss_ev`]; zero for a zero-rate row. `None` if `i`
/// is out of range.
pub fn stopping_power_ev_per_m(table: &CrossSectionTable, i: usize) -> Option<f64> {
    let inv = *table.inverse_mfp_per_m().get(i)?;
    Some(match table.quantiles(i) {
        Some(q) => inv * mean_loss_ev(table.probability(), q),
        None => 0.0,
    })
}

/// Adaptive refinement of the shared probability grid (module docs).
fn refine_probability(
    rows: &[Option<LinearDensity>],
    energy_ev: &[f64],
    mut grid: Vec<f64>,
    mut quantiles: Vec<Vec<f64>>,
    tol: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let means: Vec<f64> = quantiles
        .iter()
        .map(|q| {
            if q.is_empty() {
                0.0
            } else {
                mean_loss_ev(&grid, q)
            }
        })
        .collect();
    // The two segments touching 0 and 1 are left alone.
    let mut todo: Vec<bool> = (0..grid.len() - 1)
        .map(|s| s > 0 && s + 2 < grid.len())
        .collect();
    loop {
        let segs: Vec<usize> = (0..todo.len()).filter(|&s| todo[s]).collect();
        if segs.is_empty() || grid.len() + segs.len() > MAX_NODES {
            return (grid, quantiles);
        }
        let umid: Vec<f64> = segs
            .iter()
            .map(|&s| logistic(0.5 * (logit(grid[s]) + logit(grid[s + 1]))))
            .collect();
        // Exact quantiles at the midpoints, one row per task.
        let exact: Vec<Vec<f64>> = rows
            .par_iter()
            .zip(energy_ev.par_iter())
            .map(|(r, &e)| match r {
                Some(d) => umid.iter().map(|&u| loss_quantile(d, e, u)).collect(),
                None => Vec::new(),
            })
            .collect();
        let mut split = vec![false; segs.len()];
        for (k, &s) in segs.iter().enumerate() {
            let du = grid[s + 1] - grid[s];
            if !(umid[k] > grid[s] && umid[k] < grid[s + 1]) {
                continue; // grid exhausted to rounding
            }
            // The table interpolates linearly in u, at the logistic midpoint
            // as at any other point.
            let t = (umid[k] - grid[s]) / du;
            for (i, q) in quantiles.iter().enumerate() {
                if q.is_empty() {
                    continue;
                }
                let lin = q[s] + t * (q[s + 1] - q[s]);
                if (lin - exact[i][k]).abs() * du > tol * means[i] * du.max(1.0e-4) {
                    split[k] = true;
                    break;
                }
            }
        }
        if !split.iter().any(|&b| b) {
            return (grid, quantiles);
        }
        let mut ngrid = Vec::with_capacity(grid.len() + segs.len());
        let mut nq: Vec<Vec<f64>> = quantiles
            .iter()
            .map(|q| {
                Vec::with_capacity(if q.is_empty() {
                    0
                } else {
                    q.len() + segs.len()
                })
            })
            .collect();
        let mut ntodo = Vec::with_capacity(todo.len() + segs.len());
        let mut k = 0;
        for s in 0..todo.len() {
            ngrid.push(grid[s]);
            for (i, q) in quantiles.iter().enumerate() {
                if !q.is_empty() {
                    nq[i].push(q[s]);
                }
            }
            let mut halves = false;
            if todo[s] {
                if split[k] {
                    ngrid.push(umid[k]);
                    for (i, q) in nq.iter_mut().enumerate() {
                        if !quantiles[i].is_empty() {
                            q.push(exact[i][k]);
                        }
                    }
                    halves = true;
                }
                k += 1;
            }
            ntodo.push(halves);
            if halves {
                ntodo.push(true);
            }
        }
        ngrid.push(grid[grid.len() - 1]);
        for (i, q) in quantiles.iter().enumerate() {
            if !q.is_empty() {
                nq[i].push(q[q.len() - 1]);
            }
        }
        grid = ngrid;
        quantiles = nq;
        todo = ntodo;
    }
}

fn check_options(o: &InelasticTableOptions) -> Result<(), InelasticTableError> {
    let p = &o.probability;
    let ok = p.len() >= 2
        && p[0] == 0.0
        && p[p.len() - 1] == 1.0
        && p.iter().all(|x| x.is_finite())
        && p.windows(2).all(|w| w[1] > w[0]);
    if !ok {
        return Err(invalid(
            "probability grid",
            "need >= 2 finite, strictly increasing points from exactly 0 to exactly 1",
        ));
    }
    if !(o.density_tolerance.is_finite() && (1.0e-9..=1.0e-2).contains(&o.density_tolerance)) {
        return Err(invalid(
            "density tolerance",
            format!("{} (need 1e-9 ..= 1e-2)", o.density_tolerance),
        ));
    }
    if let Some(t) = o.refine_tolerance {
        if !(t.is_finite() && t > 0.0) {
            return Err(invalid(
                "refine tolerance",
                format!("{t} (need finite, > 0)"),
            ));
        }
    }
    Ok(())
}

/// Build the inelastic energy-loss table of `penn` on the grid of `options`.
///
/// `material` only labels the table (its identity string); the density that
/// sets the absolute rate is inside the optical ELF of `penn`
/// ([`SinglePolePenn::imfp_and_stopping`] returns `λ⁻¹` in m⁻¹ with no
/// further input). Energies are kinetic energies above the Fermi level, as
/// in the model.
pub fn build_inelastic_table(
    penn: &SinglePolePenn,
    material: &Material,
    options: &InelasticTableOptions,
) -> Result<CrossSectionTable, InelasticTableError> {
    check_no_exchange(penn)?;
    check_options(options)?;
    let energy = &options.energy_ev;
    // Rows in parallel, collected in grid order; the lowest failing energy
    // is the reported error whatever the thread count.
    type Row = Result<(f64, Option<LinearDensity>), InelasticTableError>;
    let results: Vec<Row> = energy
        .par_iter()
        .map(|&e| {
            let pt = penn.imfp_and_stopping(e)?;
            if pt.inverse_imfp_per_m > 0.0 {
                let d = loss_density(penn, e, pt.inverse_imfp_per_m, options.density_tolerance)?;
                Ok((pt.inverse_imfp_per_m, Some(d)))
            } else {
                Ok((0.0, None))
            }
        })
        .collect();
    let mut inverse_mfp_per_m = Vec::with_capacity(energy.len());
    let mut rows: Vec<Option<LinearDensity>> = Vec::with_capacity(energy.len());
    for r in results {
        let (inv, d) = r?;
        inverse_mfp_per_m.push(inv);
        rows.push(d);
    }
    let probability = &options.probability;
    let quantiles: Vec<Vec<f64>> = rows
        .par_iter()
        .zip(energy.par_iter())
        .map(|(r, &e)| match r {
            Some(d) => probability
                .iter()
                .map(|&u| loss_quantile(d, e, u))
                .collect(),
            None => Vec::new(),
        })
        .collect();
    let (probability, quantiles) = match options.refine_tolerance {
        Some(tol) => refine_probability(&rows, energy, probability.clone(), quantiles, tol),
        None => (probability.clone(), quantiles),
    };
    let elf = penn.optical_elf();
    let model = format!(
        "lindhard {} electron::inelastic::table: single-pole Penn DIIMFP (nonrelativistic), \
         W inverse CDF from an adaptive piecewise-linear density in ln W",
        env!("CARGO_PKG_VERSION")
    );
    let provenance = format!(
        "computed by lindhard {} from published formulas (see electron::inelastic docs); \
         optical ELF: {} ({}); Fermi energy {} eV; integration tolerance {}; \
         density tolerance {}; {} energies {}..{} eV; {} probability points",
        env!("CARGO_PKG_VERSION"),
        elf.material(),
        elf.provenance(),
        penn.fermi_energy_ev(),
        penn.relative_tolerance(),
        options.density_tolerance,
        energy.len(),
        energy[0],
        energy[energy.len() - 1],
        probability.len()
    );
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model,
        material: material_identity(material),
        provenance,
        axis: SamplingAxis::InelasticEnergyLoss,
        energy_ev: energy.clone(),
        inverse_mfp_per_m,
        probability,
        quantiles,
    })?)
}

// ---------------------------------------------------------------------------
// Momentum-transfer sampler

/// Samples the momentum transfer `q` of an inelastic collision with incident
/// kinetic energy `E` and loss `W` (module docs). Built once per `(E, W)`.
#[derive(Debug, Clone, PartialEq)]
pub struct MomentumTransferSampler {
    density: LinearDensity,
    /// Hartree-unit wavenumbers `k` and `k'`, and `T' = E + E_F`.
    k: f64,
    kp: f64,
    tp: f64,
}

impl MomentumTransferSampler {
    /// Build the sampler with the default tolerance [`DEFAULT_Q_TOLERANCE`].
    pub fn new(
        penn: &SinglePolePenn,
        energy_ev: f64,
        loss_ev: f64,
    ) -> Result<Self, InelasticTableError> {
        Self::with_tolerance(penn, energy_ev, loss_ev, DEFAULT_Q_TOLERANCE)
    }

    /// Build the sampler with density-model tolerance `tol` in
    /// `[1e-9, 1e-2]`. Fails if `0 < W <= E` does not hold or the slice has
    /// no support (`W` below the ELF table, or no `q` allowed).
    pub fn with_tolerance(
        penn: &SinglePolePenn,
        energy_ev: f64,
        loss_ev: f64,
        tol: f64,
    ) -> Result<Self, InelasticTableError> {
        check_no_exchange(penn)?;
        if !(energy_ev.is_finite() && energy_ev > 0.0) {
            return Err(invalid("electron energy", format!("{energy_ev} eV")));
        }
        if !(loss_ev.is_finite() && loss_ev > 0.0 && loss_ev <= energy_ev) {
            return Err(invalid(
                "energy loss",
                format!("{loss_ev} eV (need 0 < W <= E = {energy_ev} eV)"),
            ));
        }
        if !(tol.is_finite() && (1.0e-9..=1.0e-2).contains(&tol)) {
            return Err(invalid(
                "q tolerance",
                format!("{tol} (need 1e-9 ..= 1e-2)"),
            ));
        }
        let h = hartree_ev();
        let (t, w) = (energy_ev / h, loss_ev / h);
        let breaks = penn.q_slice_breaks_au(t, w).ok_or_else(|| {
            invalid(
                "momentum-transfer slice",
                format!("no allowed q at E = {energy_ev} eV, W = {loss_ev} eV"),
            )
        })?;
        let f = |u: f64| penn.q_integrand_au(w, u);
        let density = adaptive_density(&f, &breaks, 8, tol, false);
        if !(density.total() > 0.0 && density.total().is_finite()) {
            return Err(invalid(
                "momentum-transfer slice",
                format!("zero loss function at E = {energy_ev} eV, W = {loss_ev} eV"),
            ));
        }
        let tp = t + penn.fermi_energy_ev() / h;
        Ok(Self {
            density,
            k: (2.0 * tp).sqrt(),
            kp: (2.0 * (tp - w)).max(0.0).sqrt(),
            tp,
        })
    }

    /// The momentum transfer (as a wavenumber, m⁻¹) at cumulative
    /// probability `u` in `[0, 1]` (clamped).
    pub fn sample(&self, u: f64) -> f64 {
        self.density.quantile(u.clamp(0.0, 1.0)).exp() / BOHR_RADIUS
    }

    /// The support `(q_min, q_max)` of the sampler, m⁻¹: the kinematic
    /// limits `q±` cut to the tabulated ELF range.
    pub fn q_range_per_m(&self) -> (f64, f64) {
        let n = self.density.x.len();
        (
            self.density.x[0].exp() / BOHR_RADIUS,
            self.density.x[n - 1].exp() / BOHR_RADIUS,
        )
    }

    /// The kinematic limits `q- = k - k'` and `q+ = k + k'`, m⁻¹.
    pub fn kinematic_limits_per_m(&self) -> (f64, f64) {
        (
            (self.k - self.kp) / BOHR_RADIUS,
            (self.k + self.kp) / BOHR_RADIUS,
        )
    }

    /// The cosine of the primary's deflection angle for momentum transfer
    /// `q` (m⁻¹), `(k² + k'² - q²)/(2 k k')`, clamped to `[-1, 1]`.
    pub fn cos_theta(&self, q_per_m: f64) -> f64 {
        let q = q_per_m * BOHR_RADIUS;
        let c = (self.k * self.k + self.kp * self.kp - q * q) / (2.0 * self.k * self.kp);
        c.clamp(-1.0, 1.0)
    }

    /// The DIIMFP `p(E, W)` (m⁻¹ eV⁻¹) of the sampler's own density model:
    /// the integral of the density over `ln q`, divided by `π T'`. It agrees
    /// with [`SinglePolePenn::diimfp_per_m_ev`] to the model tolerance.
    pub fn diimfp_per_m_ev(&self) -> f64 {
        let h = hartree_ev();
        self.density.total() / (std::f64::consts::PI * self.tp) / (BOHR_RADIUS * h)
    }
}
