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
//! tabulated ELF energy, below which `p = 0`, for the single-pole model; see
//! "Rows below the ELF table" for the Mermin and full Penn models) and `E`,
//! and taken to be
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
//! **Floating-point resolution** (#345). Every model is a numerical integral
//! at its own tolerance, so its DIIMFP carries small jumps where the
//! adaptive quadrature changes its subdivision between two adjacent `W`.
//! The Mermin DIIMFP ([`super::mermin`], integration tolerance 1e-5, the
//! same as the default density tolerance) has them: the Si row at
//! `T = 1000 eV` of the committed Yang et al. (2019) ELF (model Fermi energy
//! 13.46 eV) has 18, of 0.4 to 7 times 1e-5 of the local density, between
//! `W` = 20 and 976 eV. A jump is not a shape the density can resolve: as
//! for the model switch at an unresolved cell, the split test of a panel
//! across it does not fall with the width, so the panel is bisected until
//! its midpoint rounds onto one of its ends. Such a panel (one ulp of
//! `ln W` wide) is not split: its midpoint is not a new node, and splitting
//! it would only add zero-width panels, up to [`MAX_NODES`] of them. This is
//! this project's choice of floor, with no tuned width: it is the point
//! where bisection can no longer place a node. A jump bisected to it costs
//! about 50 nodes (the levels from the starting panel to one ulp), and
//! leaves nothing of the density unresolved: the panel holding it has no
//! interior. Before #345 the zero-width panels were added until the row
//! reached [`MAX_NODES`]; they carry no integral (`h = 0`), so a row whose
//! other panels had converged is stored the same, with fewer evaluations
//! (the Si row above: 20 099 nodes before, 5 792 now).
//!
//! **Node bound.** A row whose refinement still has a failing panel when it
//! reaches [`MAX_NODES`] nodes is refused,
//! [`InelasticTableError::DensityNodeLimit`], rather than stored with that
//! panel unresolved (the convention of the rate refinement's row bound,
//! below).
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
//! * single-pole Penn: `W_min = elf_min`, the lowest tabulated ELF energy.
//!   Its DIIMFP is zero below it ([`super::penn`]), so its rate counts no
//!   loss there;
//! * full Penn: `W_min = 1e-6 T` ([`super::full_penn`], its
//!   `LOWEST_LOSS_FRACTION`). The DIIMFP of the full model is *not* zero
//!   below `elf_min`: each Lindhard gas of the expansion has an
//!   electron-hole continuum down to `ω -> 0` at `q > 0` (S2017 eq. (4)),
//!   and the model's rate and stopping power integrate every loss from 0. A
//!   window opened at `elf_min` (until #361) left the losses below it out of
//!   the sampled distribution while the stored rate counted them: on the Si
//!   ELF tabulated from 0.5 eV the share of `λ⁻¹` below 0.5 eV is 4.7e-3 at
//!   `T` = 10 eV, and the table's `S = λ⁻¹ ⟨W⟩` was biased high by the
//!   same amount. The lower limit is a numerical floor, not a model cutoff
//!   (a literal zero cannot bound the `ln W` axis): `p ∝ ω` at small `ω`,
//!   so the rate below `W` falls as `W²` and the stopping power as `W³`,
//!   and the omitted share of both is below the model tolerance by orders
//!   of magnitude (`tests/penn_full.rs` bounds it against the direct
//!   DIIMFP). The row's rate is the model's own; the sampled distribution
//!   and `⟨W⟩` cover `(1e-6 T, T]`, which the model's `0 < ω <= T` equals up
//!   to that omitted share. The DIIMFP grid reaches down to the floor of the
//!   lowest row energy, so these losses are read from it, not computed
//!   directly;
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
//! single-pole model these are the rows with `0 < T <= elf_min`. A Mermin
//! or full Penn row has `W_min < T` at every `T > 0`, so it is never empty
//! for this reason: its rows with `0 < T <= elf_min`, stored empty from #340
//! to #342 (Mermin) or #361 (full Penn), carry the model's rate and losses
//! in `(W_min, T]`.
//!
//! # Rate refinement
//!
//! On the band-bottom axis the rate `λ⁻¹(E)` is zero at `E_F` and climbs
//! by two orders of magnitude within about 20 eV above it, and the transport
//! interpolates it linearly between rows (`crate::electron::transport`,
//! "Free flight"). A log grid fine enough there is far finer than needed
//! everywhere else (#291, #339). [`InelasticTableOptions::rate_refinement`]
//! ([`RateRefinement`], off by default) adds rows to this one table instead:
//! every cell of the grid that overlaps `E - E_F` in
//! [`RateRefinement::window_ev`] ([`RATE_WINDOW_EV`] by default) is bisected
//! at its midpoint in `E` while the linear interpolation of the exact row
//! rate between the cell's ends differs from the exact rate at any of the
//! cell's test points by more than [`RateRefinement::tolerance`]
//! ([`RATE_TOLERANCE`]) of it. The test points are the cell's quartiles and
//! midpoint (`rate_test_points`) and, where those pass, the estimated
//! position of the largest error: the vertex of the parabola through the
//! largest of the sampled errors and its two neighbours (`rate_error_peak`).
//! Children that no longer overlap the window are not tested, so the
//! bisection never approaches the zero of the rate at `E_F`, and cells
//! outside the window are never split. Only the scalar rate is evaluated
//! while refining (the row rule of the table: zero for an empty row); each
//! level's test points are evaluated in parallel and merged in sorted order,
//! so the grid is the same on any thread count. The rows are then built as
//! any other, from the merged grid, which contains the given grid. A
//! refinement that would add more than [`RateRefinement::max_rows`] rows is
//! an error naming the cell that still fails (a bound against a rate the
//! bisection cannot resolve, not a tuned value), and so is a cell still to
//! test that can no longer be bisected because its midpoint rounds to an
//! end. A step in the rate splits one cell per level, so it reaches
//! rounding (about 50 levels) before any useful row bound: that second
//! error, not the row bound, is what stops it.
//!
//! **Rate onset** (#357). The single-pole rate is
//! identically zero for `E - E_F <= elf_min` ("Rows below the ELF table",
//! `W_min = elf_min`) and positive above it: a step at the onset
//! `E_F + elf_min`, which the transport's linear interpolation from the
//! zero row below it cannot follow at any row density. Where the onset lies
//! at or above the window start (an ELF tabulated from 2 eV or more, at the
//! default window), a cell is tested only if it overlaps the window and its
//! lower end is at or above the onset: the cell straddling the onset and
//! the cells below it are left as given, the interpolation across the
//! straddling cell being the model's own step. Testing them would bisect
//! towards the step until its cell is exhausted to rounding, the error
//! above, and so refuse an ELF that builds without the option. Where the
//! onset lies below the window start nothing changes: the bisection's
//! children that leave the window are not tested, so it never reaches the
//! step, as it never reaches the zero at `E_F`. Every committed ELF starts
//! at or below 0.5 eV (`validation/data/optical/`), so for none of them
//! does the onset change the refined grid. The Mermin rate is positive for
//! every `T > 0` (it counts losses down to `1e-8 T`), and so is the full
//! Penn rate (down to `1e-6 T`, #361: it was grouped with the single-pole
//! onset before), so those models have no onset and keep the full window.
//!
//! The test is a criterion, not a bound. A test of the midpoint alone
//! (#351) left cells whose error peaked off the midpoint above the
//! tolerance, and so did one of the quartiles and midpoint (the Al single
//! pole from 5 eV at 10 per decade, 1.0008 % at 6.1 eV, between the first
//! quartile and the midpoint of its cell); the parabolic vertex locates such
//! a peak without a tuned number of samples. The acceptance check is the dense
//! evaluation of the `rate` mode of
//! `lindhard-cli/examples/inelastic_low_energy.rs` (every 0.1 eV over the
//! window, `docs/validation.md`). The table's other grid (the elastic table
//! of a run, or a table without the option) is unchanged: the transport
//! reads each table on its own [`CrossSectionTable::energy_ev`].
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

mod build;
mod density;
mod refine;
mod sampler;

pub use build::{
    build_inelastic_table, build_inelastic_table_for_model, build_shell_channel_tables,
    default_probability_grid, InelasticTableOptions, ShellResolvedTables,
};
pub use density::{loss_density_node_count, mean_loss_ev, stopping_power_ev_per_m};
pub use refine::{
    RateRefinement, DEFAULT_RATE_REFINEMENT_MAX_ROWS, RATE_TOLERANCE, RATE_WINDOW_EV,
};
pub use sampler::MomentumTransferSampler;

#[cfg(doc)]
use super::full_penn::{DiimfpGrid, FullPenn};
use super::penn::SinglePolePenn;
use crate::electron::data::ElectronDataError;
#[cfg(doc)]
use crate::electron::data::{CrossSectionTable, SamplingAxis};
#[cfg(doc)]
use crate::electron::elastic::table::logistic_probability_grid;

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
    /// The loss density of a row still failed its split test when it reached
    /// [`MAX_NODES`] nodes (module docs, "Floating-point resolution").
    #[error(
        "the loss density at {energy_ev} eV is not resolved to its tolerance \
         within {max_nodes} nodes"
    )]
    DensityNodeLimit {
        /// Incident energy, eV.
        energy_ev: f64,
        /// The node bound, [`MAX_NODES`].
        max_nodes: usize,
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
