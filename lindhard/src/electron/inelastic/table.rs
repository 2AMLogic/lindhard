//! Inelastic sampling tables from the single-pole Penn DIIMFP (issue #94):
//! an [`electron::data::CrossSectionTable`](crate::electron::data::CrossSectionTable)
//! with [`SamplingAxis::InelasticEnergyLoss`] (the inverse mean free path and
//! the inverse CDF of the energy loss `W` at every energy of a caller-chosen
//! grid), and a sampler of the momentum transfer `q` given `(E, W)`, from
//! which the deflection angle follows.
//!
//! Both are built from [`SinglePolePenn`] only; no data and no equation beyond
//! the model of [`super::penn`] enter here (the elementary kinematics of a
//! binary collision aside, see below). [`build_inelastic_table_for_model`]
//! runs the same energy-loss table procedure on any of the three models of
//! [`super::model::PennInelastic`] (single pole, full Penn, Mermin-ELF),
//! using that model's own DIIMFP and inverse IMFP; nothing else changes.
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
//! chosen nodes between the lowest loss of the model's rate (the lowest
//! tabulated ELF energy, below which `p = 0`, for the Penn models; see "Rows
//! below the ELF table" for the Mermin model) and `E`, and taken to be
//! *linear between nodes*. Starting from
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
//! **Full Penn rows.** For the full Penn model the density is built on the
//! same nodes from the DIIMFP of a [`DiimfpGrid`] tabulated once for the
//! table's energies ([`FullPenn::diimfp_grid`], #256), not from the direct
//! nested integrals: its interpolation error, held to the model tolerance
//! (1e-4, relative or relative to the row's mean density), is a second
//! approximation of the loss distribution, next to the density tolerance.
//! It is also smooth in `W`, where the direct DIIMFP carries quadrature
//! noise at the model tolerance, which is above the default density
//! tolerance and would drive the panel refinement to [`MAX_NODES`].
//!
//! **Unresolved grid cells.** Where the grid could not hold its
//! interpolation to the model tolerance even at its narrowest cell (1e-4 in
//! `ln ω`), it leaves the cell to the direct model
//! ([`DiimfpGrid::unresolved_cells`]). The row integrand then switches model
//! at each edge of such a cell, a step of the order of the model tolerance,
//! and inside it carries the direct model's quadrature noise. The split test
//! above does not depend on the panel width for a step (both sides scale
//! with `h`), so a step above about `3 · tol` of the row's mean density
//! would be bisected down to floating-point resolution and then padded with
//! duplicate nodes up to [`MAX_NODES`], each a direct evaluation. So a panel
//! within one grid-cell width (1e-4) of an unresolved cell is not split
//! once it is no wider than that width: that is the finest scale at which
//! the grid tested the DIIMFP, and below it the split test measures the
//! noise and the step, not the shape. What such a panel can still miss is
//! bounded by its width times the step, `1e-4 · (model tolerance) · g`,
//! per cell edge, far below the model tolerance of the row. Rows that do
//! not come near an unresolved cell, and every model without a grid, are
//! refined exactly as before. [`loss_density_node_count`] reports the node
//! count of a row.
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
//! Rows where `λ⁻¹ = 0` (no loss allowed) are stored empty, and so are rows
//! with `T` at or below the lowest loss of the model's rate (module docs,
//! "Energy axis"). The table covers
//! the direct single-pole model only: no exchange and no surface or band-gap
//! structure (see [`super::penn`]). Inner shells get one table each, beside
//! the valence table, from [`build_shell_channel_tables`]. A model with the
//! exchange correction enabled ([`SinglePolePenn::with_exchange`]) is
//! rejected by [`build_inelastic_table`] and [`MomentumTransferSampler`].
//!
//! # Energy axis
//!
//! The models take `T`, the kinetic energy above the model's Fermi level,
//! with kinematics on `T' = T + E_F` and losses up to `T` (S2017 eqs.
//! (2)-(3), [`super::penn`]). By default ([`EnergyAxis::ModelFermiLevel`])
//! a table is stored on that axis: the row at `E` is the model's at `T = E`.
//!
//! [`EnergyAxis::BandBottom`] stores it on the band-bottom axis instead: the
//! row at `E` is the model's at `T = E - E_F`, so its kinematics are on
//! `T' = E` and its losses stop at `E - E_F`; rows with `T <= 0` are stored
//! empty. With `E_F` the band's, this is the convention of the table
//! compiler of Nebula, `compile_full_imfp_icdf` in
//! `cstool/dielectric_function/compile.py` (Nebula-simulator/cstool commit
//! `0c739eb3fcc3fe5297e74c601ac4a9546db596cf`, BSD-3-Clause): it evaluates
//! every row at the electron's kinetic energy `K` (the momenta of the
//! kinematic limits are those of `K` and `K - ω`) and keeps only the losses
//! `ω < K - F` (its "Fermi correction"), `F` the energy its caller
//! `compile_full_penn` (`apps/cstool.py`, same commit) passes, the band's
//! `get_min_excitation()` (`cstool/input_data/band_structure.py`): the Fermi
//! energy of a metal, the conduction-band bottom `W_v + E_g` of an
//! insulator or semiconductor. The model's own `E_F` is the `F` of the axis,
//! so there is one Fermi energy, not two. The table builder knows no band;
//! the caller sets the model's Fermi energy (`lindhard run` does so per
//! material, `docs/cli.md`).
//!
//! **Rows below the ELF table.** The sampled loss window of a row is
//! `(W_min, T]`, `W_min` the lowest loss the model's rate `λ⁻¹(T)` counts,
//! so that the stored rate and the sampled losses cover the same window:
//!
//! * single-pole and full Penn: `W_min = elf_min`, the lowest tabulated ELF
//!   energy. Their DIIMFP is zero below it ([`super::penn`],
//!   [`super::full_penn`]), so their rate counts no loss there;
//! * Mermin-ELF: `W_min = 1e-8 T`, the lower limit of the `ω` integral of
//!   its rate and stopping power ([`super::mermin`], the same constant). Its
//!   fitted Drude-Lorentz ELF has no cutoff at `elf_min` (the omission of the
//!   threshold step is a documented choice of that model), so the DIIMFP
//!   below `elf_min` is the fit's extension below the tabulated range and is
//!   part of both the rate and the sampled loss distribution (#342). Before
//!   #342 the window opened at `elf_min` for every model, and a Mermin row
//!   drew the losses its rate counts below `elf_min` as losses above it: for
//!   the Si ELF tabulated from 0.5 eV (default fit, model Fermi energy
//!   13.46 eV) the window `(0.5 eV, T]` held 28 % of the rate at `T` = 0.6
//!   eV, 72 % at 1 eV, 92 % at 2 eV and 98.4 % at 5 eV.
//!
//! A row with `0 < T <= W_min` has an empty window, so no loss can be drawn
//! there, and it is stored empty (inverse mean free path 0, no quantiles),
//! as a row with `T <= 0` is, on either axis (#340): without this the `ln W`
//! window `[ln W_min, ln T]` is reversed (or of zero width at `T = W_min`)
//! and the density integrates to a negative value (or 0). For the Penn
//! models these are the rows with `0 < T <= elf_min`. A Mermin row has
//! `W_min < T` at every `T > 0`, so it is never empty for this reason: its
//! rows with `0 < T <= elf_min`, stored empty from #340 to #342, carry the
//! model's rate and losses in `(1e-8 T, T]`.
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
//! Rows (and the profiles of a full-Penn [`DiimfpGrid`]) are built in
//! parallel with rayon and collected in grid order; each
//! row, and the probability refinement (a parallel map over rows followed by
//! a serial split decision), is a fixed sequence of floating-point
//! operations, so tables are bit-identical on any thread count (tested).
//! Samplers consume uniforms `u` supplied by the caller (for example from
//! [`crate::rng::stream`]) and keep no state of their own.

use super::full_penn::DiimfpGrid;
#[cfg(doc)]
use super::full_penn::FullPenn;
use super::inner_shell::ShellResolvedChannels;
use super::mermin::LOWEST_LOSS_FRACTION;
use super::model::PennInelastic;
use super::penn::{hartree_ev, SinglePolePenn};
use crate::constants::BOHR_RADIUS;
use crate::electron::data::{
    CrossSectionTable, CrossSectionTableParts, ElectronDataError, SamplingAxis, ShellChannelTable,
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
/// panels. `rough` lists ranges of `x` (ascending) where `f` is known not to
/// be smooth below the width `rough_width` (the unresolved cells of a
/// [`DiimfpGrid`], module docs): a panel within `rough_width` of one of them
/// is not split once it is no wider than `rough_width`. With no ranges the
/// refinement is the plain one.
fn adaptive_density(
    f: &(impl Fn(f64) -> f64 + ?Sized),
    breaks: &[f64],
    sub: usize,
    tol: f64,
    first_moment: bool,
    rough: &[(f64, f64)],
    rough_width: f64,
) -> LinearDensity {
    let at_resolution = |a: f64, b: f64| {
        b - a <= rough_width
            && rough
                .iter()
                .any(|&(lo, hi)| a <= hi + rough_width && b >= lo - rough_width)
    };
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
            if bad && x.len() + 1 < MAX_NODES && !at_resolution(x[i], x[i + 1]) {
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

/// The energy axis a table is stored on (module docs, "Energy axis").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EnergyAxis {
    /// The model's own axis: a row at energy `E` is the model's at `T = E`,
    /// the kinetic energy above the model's Fermi level, with losses up to
    /// `E`.
    #[default]
    ModelFermiLevel,
    /// The band-bottom axis: a row at energy `E` is the model's at
    /// `T = E - E_F`, `E_F` the model's Fermi energy, so the kinematics use
    /// `T' = T + E_F = E` and the losses stop at `E - E_F`. Rows with
    /// `T <= 0` are stored empty. This is the convention of cstool's
    /// `compile_full_imfp_icdf` (module docs).
    BandBottom,
}

/// Grids and tolerances of a table build.
#[derive(Debug, Clone, PartialEq)]
pub struct InelasticTableOptions {
    /// Incident kinetic energies on the axis of [`Self::axis`], eV (strictly
    /// increasing, >= 2 points): above the model's Fermi level by default,
    /// above the band bottom with [`EnergyAxis::BandBottom`].
    pub energy_ev: Vec<f64>,
    /// The axis `energy_ev` is on ([`EnergyAxis::ModelFermiLevel`] by
    /// default).
    pub axis: EnergyAxis,
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
            axis: EnergyAxis::ModelFermiLevel,
            probability: default_probability_grid(),
            density_tolerance: DEFAULT_DENSITY_TOLERANCE,
            refine_tolerance: Some(DEFAULT_REFINE_TOLERANCE),
        }
    }

    /// The same options with the table stored on `axis`.
    pub fn with_axis(mut self, axis: EnergyAxis) -> Self {
        self.axis = axis;
        self
    }
}

/// The default starting probability grid: [`logistic_probability_grid`] with
/// [`DEFAULT_PROBABILITY_TAIL`] and [`DEFAULT_LOGISTIC_STEP`].
pub fn default_probability_grid() -> Vec<f64> {
    logistic_probability_grid(DEFAULT_PROBABILITY_TAIL, DEFAULT_LOGISTIC_STEP)
        .expect("valid default grid")
}

/// What the table builder needs from an ELF-extension model: the lowest
/// tabulated ELF energy, the lowest loss of its rate integral, the DIIMFP
/// and the inverse IMFP. Implemented for [`SinglePolePenn`] and for every
/// model of [`PennInelastic`]; the table is built the same way for each
/// (module docs).
trait LossModel: Sync {
    fn elf_min_ev(&self) -> f64;
    /// The lowest loss (eV) the model's inverse IMFP at `energy_ev` counts,
    /// which opens the sampled loss window of the row (module docs, "Rows
    /// below the ELF table"). By default the lowest tabulated ELF energy,
    /// below which the DIIMFP is zero.
    fn loss_min_ev(&self, _energy_ev: f64) -> f64 {
        self.elf_min_ev()
    }
    fn fermi_energy_ev(&self) -> f64;
    fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError>;
    fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError>;
    /// A tabulation of the DIIMFP for the rows at `energy_ev`, if the model
    /// has one ([`FullPenn::diimfp_grid`]); the rows then read the DIIMFP
    /// from it wherever it covers the loss.
    fn diimfp_grid(&self, _energy_ev: &[f64]) -> Result<Option<DiimfpGrid>, ElectronDataError> {
        Ok(None)
    }
}

impl LossModel for SinglePolePenn {
    fn elf_min_ev(&self) -> f64 {
        self.optical_elf().energy_ev()[0]
    }
    fn fermi_energy_ev(&self) -> f64 {
        SinglePolePenn::fermi_energy_ev(self)
    }
    fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
        self.diimfp_per_m_ev(energy_ev, loss_ev)
    }
    fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
        Ok(self.imfp_and_stopping(energy_ev)?.inverse_imfp_per_m)
    }
}

impl LossModel for PennInelastic {
    fn elf_min_ev(&self) -> f64 {
        self.optical_elf().energy_ev()[0]
    }
    fn loss_min_ev(&self, energy_ev: f64) -> f64 {
        match self {
            // the lower limit of the Mermin rate integral (#342); the fitted
            // ELF has no cutoff at the lowest tabulated energy
            PennInelastic::Mermin(_) => LOWEST_LOSS_FRACTION * energy_ev,
            // the single-pole and full Penn DIIMFPs are zero below it
            PennInelastic::SinglePole(_) | PennInelastic::Full(_) => self.elf_min_ev(),
        }
    }
    fn fermi_energy_ev(&self) -> f64 {
        PennInelastic::fermi_energy_ev(self)
    }
    fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
        self.diimfp_per_m_ev(energy_ev, loss_ev)
    }
    fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
        Ok(self.imfp_and_stopping(energy_ev)?.inverse_imfp_per_m)
    }
    fn diimfp_grid(&self, energy_ev: &[f64]) -> Result<Option<DiimfpGrid>, ElectronDataError> {
        match self {
            PennInelastic::Full(m) => Ok(Some(m.diimfp_grid(energy_ev)?)),
            _ => Ok(None),
        }
    }
}

/// The DIIMFP of `penn` at `(energy_ev, loss_ev)`, from `grid` where it
/// covers the point; zero where the model fails (only for non-finite input).
fn row_diimfp<M: LossModel>(
    penn: &M,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    loss_ev: f64,
) -> f64 {
    grid.and_then(|g| g.diimfp_per_m_ev(energy_ev, loss_ev))
        .unwrap_or_else(|| penn.diimfp(energy_ev, loss_ev).unwrap_or(0.0))
}

/// The density model of the row at `energy_ev`, on `ln W` from the model's
/// lowest loss ([`LossModel::loss_min_ev`]) to `energy_ev` (module docs).
fn row_density<M: LossModel>(
    penn: &M,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    tol: f64,
) -> LinearDensity {
    let w_lo = penn.loss_min_ev(energy_ev);
    let f = |s: f64| {
        let w = s.exp().min(energy_ev);
        row_diimfp(penn, grid, energy_ev, w) * w
    };
    // The unresolved cells of the grid, on the row's `ln W` axis (module
    // docs, "Unresolved grid cells"); none without a grid.
    let (rough, rough_width) = match grid {
        Some(g) => (g.unresolved_ln_loss_ev(), g.min_cell_width()),
        None => (Vec::new(), 0.0),
    };
    adaptive_density(
        &f,
        &[w_lo.ln(), energy_ev.ln()],
        INITIAL_LOSS_PANELS,
        tol,
        true,
        &rough,
        rough_width,
    )
}

/// The number of nodes of the loss-density model of the row at `energy_ev`
/// (eV, on the model's own axis) of a table of `model` built with density
/// tolerance `density_tolerance`, reading the DIIMFP from `grid` where it
/// covers the loss as a table build does (pass the grid of
/// [`FullPenn::diimfp_grid`] for the table's energies, or `None`). The
/// count is bounded by [`MAX_NODES`]; a row near that bound has spent its
/// evaluations on something the split test cannot resolve. Each node not
/// covered by the grid costs one direct DIIMFP evaluation.
pub fn loss_density_node_count(
    model: &PennInelastic,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    density_tolerance: f64,
) -> usize {
    row_density(model, grid, energy_ev, density_tolerance)
        .x
        .len()
}

/// The loss density of one energy.
fn loss_density<M: LossModel>(
    penn: &M,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    inverse_mfp: f64,
    tol: f64,
) -> Result<LinearDensity, InelasticTableError> {
    let d = row_density(penn, grid, energy_ev, tol);
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
/// further input). Energies are on the axis of `options.axis`: kinetic
/// energies above the model's Fermi level by default, above the band bottom
/// with [`EnergyAxis::BandBottom`] (module docs, "Energy axis").
pub fn build_inelastic_table(
    penn: &SinglePolePenn,
    material: &Material,
    options: &InelasticTableOptions,
) -> Result<CrossSectionTable, InelasticTableError> {
    check_no_exchange(penn)?;
    check_options(options)?;
    let (inverse_mfp_per_m, probability, quantiles) = loss_rows(penn, options)?;
    let energy = &options.energy_ev;
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
    ) + &axis_provenance(penn, options);
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

/// Build the inelastic energy-loss table of any of the three ELF-extension
/// models ([`PennInelastic`]: single-pole Penn, full Penn or Mermin-ELF) on
/// the grid of `options`, by the same procedure as [`build_inelastic_table`]
/// (module docs), with that model's DIIMFP and inverse IMFP in place of the
/// single pole's. For [`PennInelastic::SinglePole`] the rows are those of
/// [`build_inelastic_table`]; only the `model` and `provenance` strings
/// differ, and they carry [`PennInelastic::model_identity`], so a table
/// names the algorithm it was built with.
///
/// The full Penn and Mermin models evaluate their DIIMFP by numerical
/// integration, so a table of either costs far more than a single-pole table
/// on the same grid. For the full model the rows read the DIIMFP from a
/// [`DiimfpGrid`] built once for the grid energies
/// ([`FullPenn::diimfp_grid`]: interpolation error held to the model
/// tolerance), which is what makes such a table affordable (#256); the
/// inverse mean free paths are still the model's own.
pub fn build_inelastic_table_for_model(
    model: &PennInelastic,
    material: &Material,
    options: &InelasticTableOptions,
) -> Result<CrossSectionTable, InelasticTableError> {
    if let PennInelastic::SinglePole(p) = model {
        check_no_exchange(p)?;
    }
    check_options(options)?;
    let (inverse_mfp_per_m, probability, quantiles) = loss_rows(model, options)?;
    let energy = &options.energy_ev;
    let elf = model.optical_elf();
    let identity = model.model_identity();
    // the full model's rows read the DIIMFP from its loss grid (#256)
    let via = match model {
        PennInelastic::Full(_) => " (from the T-independent loss grid FullPenn::diimfp_grid)",
        _ => "",
    };
    let model_text = format!(
        "lindhard {} electron::inelastic::table: {identity} DIIMFP{via}, \
         W inverse CDF from an adaptive piecewise-linear density in ln W",
        env!("CARGO_PKG_VERSION")
    );
    let provenance = format!(
        "computed by lindhard {} from published formulas (see electron::inelastic docs); \
         model: {identity}; optical ELF: {} ({}); density tolerance {}; \
         {} energies {}..{} eV; {} probability points",
        env!("CARGO_PKG_VERSION"),
        elf.material(),
        elf.provenance(),
        options.density_tolerance,
        energy.len(),
        energy[0],
        energy[energy.len() - 1],
        probability.len()
    ) + &axis_provenance(model, options);
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model: model_text,
        material: material_identity(material),
        provenance,
        axis: SamplingAxis::InelasticEnergyLoss,
        energy_ev: energy.clone(),
        inverse_mfp_per_m,
        probability,
        quantiles,
    })?)
}

/// The model's kinetic energy `T` of each row of `options`: the axis energy
/// itself on [`EnergyAxis::ModelFermiLevel`], `E - E_F` (`E_F` the model's
/// Fermi energy) on [`EnergyAxis::BandBottom`].
fn model_energies_ev<M: LossModel>(penn: &M, options: &InelasticTableOptions) -> Vec<f64> {
    match options.axis {
        EnergyAxis::ModelFermiLevel => options.energy_ev.clone(),
        EnergyAxis::BandBottom => {
            let ef = penn.fermi_energy_ev();
            options.energy_ev.iter().map(|&e| e - ef).collect()
        }
    }
}

/// The provenance text of the axis: empty on the model's own axis (so tables
/// built before the axis option existed keep their bytes), the axis and the
/// Fermi energy otherwise.
fn axis_provenance<M: LossModel>(penn: &M, options: &InelasticTableOptions) -> String {
    match options.axis {
        EnergyAxis::ModelFermiLevel => String::new(),
        EnergyAxis::BandBottom => format!(
            "; energy axis: band bottom (row at E is the model's at T = E - {} eV)",
            penn.fermi_energy_ev()
        ),
    }
}

/// The tables of a [`ShellResolvedChannels`]: the valence channel's and one
/// per inner shell, in the order of
/// [`ShellResolvedChannels::inner_shells`]
/// ([`build_shell_channel_tables`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ShellResolvedTables {
    /// The valence channel, as [`build_inelastic_table`] of
    /// [`ShellResolvedChannels::valence`].
    pub valence: CrossSectionTable,
    /// The inner-shell channels.
    pub shells: Vec<ShellChannelTable>,
}

/// Build the energy-loss tables of every channel of `channels` on the grid of
/// `options`: the valence table, and for inner shell `i` the table of
/// [`ShellResolvedChannels::shell_model`] (that shell's own optical ELF with
/// the valence model's settings), wrapped with the shell's identity, binding
/// energy and `binding_provenance` (the origin of the binding energies, for
/// example [`crate::electron::data::SubshellBindingTable::provenance`]).
///
/// Each table is built by [`build_inelastic_table`], one after the other in
/// channel order, so the result is bit-identical at any thread count. Its
/// rate is the channel's own `λ⁻¹_j(E)` and its rows the inverse CDF of that
/// channel's DIIMFP: the **channel-first** construction, in which a
/// transport picks a channel with probability `λ⁻¹_j / Σ λ⁻¹` and then draws
/// the loss from that channel's table (not the conditional channel choice of
/// [`ShellResolvedChannels::sample_channel`], which belongs with a loss drawn
/// from the total DIIMFP; the two are not mixed).
///
/// **Supported model.** Single-pole Penn per channel only, without the
/// exchange correction (rejected, as by [`build_inelastic_table`]); so the
/// binding energy does not enter a shell's DIIMFP, which is zero below the
/// first energy of the shell's ELF (at or above `B` by construction of
/// [`ShellResolvedChannels::new`]). Each shell row is floored at `B`, which
/// only removes the rounding of `exp(ln W)` at the edge. No full Penn or
/// Mermin valence model is wired to shell channels, and the valence ELF must
/// not already contain the shells (no partition rule is applied, see
/// [`super::inner_shell`]).
///
/// **Energy axis.** As for [`build_inelastic_table`], every table (valence
/// and shells) is on the axis of `options.axis`: kinetic energies above the
/// model's Fermi energy (zero by default) with
/// [`EnergyAxis::ModelFermiLevel`], above the band bottom with
/// [`EnergyAxis::BandBottom`], each row then the model's at `T = E - E_F`
/// with `E_F` the Fermi energy of the model (the shell models share the
/// valence model's settings, so its Fermi energy) (module docs, "Energy
/// axis"). See `crate::electron::transport`, "Energy reference of the
/// inelastic table", for how the transport reads them.
pub fn build_shell_channel_tables(
    channels: &ShellResolvedChannels,
    material: &Material,
    options: &InelasticTableOptions,
    binding_provenance: &str,
) -> Result<ShellResolvedTables, InelasticTableError> {
    let valence = build_inelastic_table(channels.valence(), material, options)?;
    let mut shells = Vec::with_capacity(channels.inner_shells().len());
    for (i, s) in channels.inner_shells().iter().enumerate() {
        let b = s.binding_energy_ev;
        let base = build_inelastic_table(channels.shell_model(i), material, options)?;
        let mut parts = base.parts().clone();
        let label = s.subshell.label();
        parts.model = format!("{} [inner-shell channel Z = {} {label}]", parts.model, s.z);
        parts.provenance = format!(
            "{}; inner-shell channel Z = {} {label}, binding energy {b} eV ({}), \
             losses floored at the binding energy",
            parts.provenance,
            s.z,
            binding_provenance.trim()
        );
        for row in &mut parts.quantiles {
            for w in row.iter_mut() {
                *w = w.max(b);
            }
        }
        let table = CrossSectionTable::new(parts)?;
        shells.push(ShellChannelTable::new(
            s.z,
            s.subshell,
            b,
            binding_provenance,
            table,
        )?);
    }
    Ok(ShellResolvedTables { valence, shells })
}

/// Inverse mean free paths, the (refined) probability grid and the quantile
/// rows of a table of `penn` on the grid of `options` (already checked).
type LossRows = (Vec<f64>, Vec<f64>, Vec<Vec<f64>>);

fn loss_rows<M: LossModel>(
    penn: &M,
    options: &InelasticTableOptions,
) -> Result<LossRows, InelasticTableError> {
    // The model's energy `T` of each row (module docs, "Energy axis").
    let energy = &model_energies_ev(penn, options);
    // The DIIMFP grid covers only the rows that get losses: on the
    // band-bottom axis a row at or below the Fermi level has none (below),
    // and the grid would reject its energy.
    let grid = match options.axis {
        EnergyAxis::ModelFermiLevel => penn.diimfp_grid(energy)?,
        EnergyAxis::BandBottom => {
            let lossy: Vec<f64> = energy.iter().copied().filter(|&e| e > 0.0).collect();
            penn.diimfp_grid(&lossy)?
        }
    };
    // Rows in parallel, collected in grid order; the lowest failing energy
    // is the reported error whatever the thread count.
    type Row = Result<(f64, Option<LinearDensity>), InelasticTableError>;
    let results: Vec<Row> = energy
        .par_iter()
        .map(|&e| {
            if options.axis == EnergyAxis::BandBottom && e <= 0.0 {
                // A band-bottom energy at or below the Fermi level: no loss
                // is allowed. (On the model's own axis the model rejects
                // such an energy.)
                return Ok((0.0, None));
            }
            if e <= penn.loss_min_ev(e) {
                // The sampled loss window (loss_min, T] is empty (#340): no
                // loss can be drawn, so the row is empty, whatever rate the
                // model reports below its lowest loss (module docs, "Rows
                // below the ELF table").
                return Ok((0.0, None));
            }
            let inv = penn.inverse_imfp(e)?;
            if inv > 0.0 {
                let d = loss_density(penn, grid.as_ref(), e, inv, options.density_tolerance)?;
                Ok((inv, Some(d)))
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
    Ok((inverse_mfp_per_m, probability, quantiles))
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
        let density = adaptive_density(&f, &breaks, 8, tol, false, &[], 0.0);
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A row-like integrand in `x = ln W` (a smooth loss peak) that, inside
    /// the cell `[c, c + 1.5e-4]`, is replaced by a second model differing
    /// by 1e-3 relative and carrying noise at 1e-4: the shape of a full-Penn
    /// row at an unresolved [`DiimfpGrid`] cell, where the grid gives way to
    /// the direct model (module docs, "Unresolved grid cells"). Synthetic,
    /// not a material.
    fn stepped_row(c: f64) -> impl Fn(f64) -> f64 {
        move |x: f64| {
            let smooth = (-(x - 3.0).powi(2)).exp();
            if (c..=c + 1.5e-4).contains(&x) {
                smooth * (1.0 + 1e-3 + 1e-4 * (1e7 * x).sin())
            } else {
                smooth
            }
        }
    }

    #[test]
    fn a_step_at_an_unresolved_cell_does_not_exhaust_the_nodes() {
        // A starting node (x = 3, the 25th of 64 panels on [0, 8]) falls in
        // the cell, at the peak, so the step is sampled and far above the
        // split threshold (about 3 tol of the mean).
        let c = 3.0 - 0.5e-4;
        let f = stepped_row(c);
        let breaks = [0.0, 8.0];
        let cells = [(c, c + 1.5e-4)];
        let d = adaptive_density(&f, &breaks, INITIAL_LOSS_PANELS, 1e-5, true, &cells, 1e-4);
        let smooth = adaptive_density(
            &|x: f64| (-(x - 3.0).powi(2)).exp(),
            &breaks,
            INITIAL_LOSS_PANELS,
            1e-5,
            true,
            &[],
            0.0,
        );
        let (n, n0) = (d.x.len(), smooth.x.len());
        eprintln!("nodes: {n} with the stepped cell declared, {n0} for the smooth row");
        // The cell costs a few bisection levels around its two edges, not
        // the ~MAX_NODES of a bisection to floating-point resolution.
        assert!(n < n0 + 100, "{n} nodes vs {n0} (MAX_NODES = {MAX_NODES})");
        // The density still matches the smooth row's integral to far below
        // the model tolerance: the step adds ~1.5e-7 of it.
        let rel = (d.total() / smooth.total() - 1.0).abs();
        assert!(rel < 1e-5, "{rel:e}");
        // No panel was bisected towards floating-point resolution.
        let narrowest = d.x.windows(2).map(|w| w[1] - w[0]).fold(f64::MAX, f64::min);
        assert!(narrowest > 1e-5, "{narrowest:e}");
    }

    #[test]
    fn without_unresolved_cells_the_refinement_is_unchanged() {
        // A smooth row passes no rough ranges and is refined exactly as
        // before; declaring a cell far from the row's structure changes no
        // node either, since no panel near it ever gets that narrow.
        let f = |x: f64| (-(x - 3.0).powi(2)).exp();
        let breaks = [0.0, 8.0];
        let plain = adaptive_density(&f, &breaks, INITIAL_LOSS_PANELS, 1e-5, true, &[], 0.0);
        let far = adaptive_density(
            &f,
            &breaks,
            INITIAL_LOSS_PANELS,
            1e-5,
            true,
            &[(7.5, 7.5 + 1.5e-4)],
            1e-4,
        );
        assert_eq!(plain, far);
    }

    /// The committed Si ELF (Yang et al. 2019, `validation/data/optical`,
    /// tabulated from 0.5 eV), or `None` in a packaged crate without it.
    fn si_elf() -> Option<crate::electron::data::OpticalElf> {
        let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../validation/data/optical/si_elf_yang2019.toml");
        if !p.is_file() {
            eprintln!("{} not found; Si check skipped", p.display());
            return None;
        }
        Some(crate::electron::data::OpticalElf::from_toml_file(p).unwrap())
    }

    #[test]
    fn mermin_rows_below_the_elf_table_open_at_the_rate_lower_limit() {
        // Issues #340 and #342: the Mermin model of Si on the band-bottom
        // axis (E_F = W_v + E_g = 13.46 eV, the band of
        // validation/experiments/se_yield.py) with rows at T = E - E_F inside
        // (0, elf_min]. The Mermin rate there is positive (its fitted ELF has
        // no lower cutoff). Before #340 these rows were refused (a reversed
        // ln W window); from #340 to #342 they were stored empty; now the
        // window opens at the lower limit of the rate integral, 1e-8 T, and
        // they carry the model's rate.
        let Some(elf) = si_elf() else { return };
        let ef = 13.46;
        let model = PennInelastic::try_new(super::super::model::PennAlgorithm::Mermin, elf)
            .unwrap()
            .with_fermi_energy_ev(ef)
            .unwrap();
        let elf_min = model.elf_min_ev();
        assert_eq!(elf_min, 0.5);
        let t_low = [0.2085, 0.2530, 0.4216];
        let mut grid: Vec<f64> = t_low.iter().map(|t| ef + t).collect();
        // the T = elf_min row, exactly on the band-bottom axis, and one above
        grid.push(ef + elf_min);
        grid.push(ef + 1.0);
        let opts = InelasticTableOptions::new(grid.clone()).with_axis(EnergyAxis::BandBottom);
        let (inv, _, q) = loss_rows(&model, &opts).expect("rows inside (0, elf_min] build");
        for (k, &e) in grid.iter().enumerate() {
            let t = e - ef;
            let w_min = model.loss_min_ev(t);
            assert_eq!(w_min, LOWEST_LOSS_FRACTION * t);
            assert_eq!(inv[k], model.inverse_imfp(t).unwrap(), "T = {t}");
            assert!(inv[k] > 0.0, "T = {t}");
            let row = &q[k];
            // the window is (1e-8 T, T]: the lowest quantile is its lower
            // edge, far below elf_min, and the highest is T
            assert!((row[0] / w_min - 1.0).abs() < 1e-12, "T = {t}: {}", row[0]);
            assert!(row.iter().all(|&w| w > 0.0 && w <= t), "T = {t}");
            assert!((row[row.len() - 1] - t).abs() <= 1e-12 * t, "T = {t}");
            assert!(row.windows(2).all(|w| w[1] >= w[0]), "T = {t}");
        }
    }

    /// A loss model whose rate counts losses down to `loss_min < elf_min`,
    /// with a flat DIIMFP on `(loss_min, T]`: synthetic.
    struct WideWindow;

    impl LossModel for WideWindow {
        fn elf_min_ev(&self) -> f64 {
            0.5
        }
        fn loss_min_ev(&self, energy_ev: f64) -> f64 {
            0.01 * energy_ev
        }
        fn fermi_energy_ev(&self) -> f64 {
            0.0
        }
        fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(if loss_ev > 0.0 && loss_ev <= energy_ev {
                1.0
            } else {
                0.0
            })
        }
        fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(energy_ev * 0.99)
        }
    }

    #[test]
    fn a_row_opens_at_the_model_lowest_loss() {
        // T = 0.3 lies in (loss_min, elf_min]: the row is built, not empty,
        // and both rows open at loss_min = 0.01 T, not at elf_min.
        let opts = InelasticTableOptions::new(vec![0.3, 2.0]);
        let (inv, prob, q) = loss_rows(&WideWindow, &opts).unwrap();
        assert_eq!(inv, vec![0.3 * 0.99, 2.0 * 0.99]);
        for (row, t) in q.iter().zip([0.3, 2.0]) {
            assert!((row[0] / (0.01 * t) - 1.0).abs() < 1e-12, "{}", row[0]);
            assert!((row[row.len() - 1] - t).abs() <= 1e-12 * t);
            // a flat DIIMFP: W(u) = loss_min + u (T - loss_min)
            for (&u, &w) in prob.iter().zip(row) {
                let want = 0.01 * t + u * 0.99 * t;
                assert!((w - want).abs() <= 1e-4 * t, "u {u}: {w} vs {want}");
            }
        }
    }

    /// A loss model with a positive rate and no DIIMFP anywhere: synthetic.
    struct RateWithoutDensity;

    impl LossModel for RateWithoutDensity {
        fn elf_min_ev(&self) -> f64 {
            0.5
        }
        fn fermi_energy_ev(&self) -> f64 {
            0.0
        }
        fn diimfp(&self, _: f64, _: f64) -> Result<f64, ElectronDataError> {
            Ok(0.0)
        }
        fn inverse_imfp(&self, _: f64) -> Result<f64, ElectronDataError> {
            Ok(1.0e6)
        }
    }

    #[test]
    fn a_positive_rate_without_density_above_the_elf_table_is_still_refused() {
        // Below or at elf_min (the default lowest loss, as for the Penn
        // models) the row is empty; above it, a positive rate with a density
        // that integrates to zero is a defect, refused.
        let opts = InelasticTableOptions::new(vec![0.25, 0.5]);
        let (inv, _, q) = loss_rows(&RateWithoutDensity, &opts).unwrap();
        assert_eq!(inv, vec![0.0, 0.0]);
        assert!(q.iter().all(|r| r.is_empty()));
        let opts = InelasticTableOptions::new(vec![0.25, 2.0]);
        match loss_rows(&RateWithoutDensity, &opts) {
            Err(InelasticTableError::EmptyDistribution { energy_ev, .. }) => {
                assert_eq!(energy_ev, 2.0)
            }
            other => panic!("expected EmptyDistribution, got {other:?}"),
        }
    }
}
