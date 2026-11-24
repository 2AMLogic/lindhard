//! The electron tally: energy deposition, emission yields and spectra, wired
//! to the electron transport hooks.
//!
//! [`FullElectronTally`] implements [`ElectronTally`]. Build one prototype
//! with [`FullElectronTally::new`] and pass clones to
//! [`Transport::run`](crate::electron::transport::Transport::run). The
//! merged tally turns into a plain-data, serializable [`ElectronReport`] with
//! [`FullElectronTally::report`].
//!
//! Units: positions and lengths in metres, energies in eV, angles in radians.
//! Totals are summed over all histories. Divide by [`ElectronReport::histories`]
//! for per-primary values (the yields are already per primary).
//!
//! # What is tallied
//!
//! - **Energy deposition.** An inelastic collision deposits the energy it
//!   leaves in the solid at the collision point: its whole loss `W` without a
//!   secondary model, and with one only the part no secondary carries away
//!   (see [the energy balance](#the-energy-balance)). An electron whose
//!   energy falls below the stopping threshold of its layer deposits what it
//!   has left at the point where it stopped. Both go into a per-layer total
//!   and into optional 3D grids: a
//!   Cartesian `x`-`y`-`z` grid ([`CartesianGrid`]) and a cylindrical `r`-`x`
//!   grid about the beam axis ([`CylindricalGrid`]). Energy deposited outside
//!   a grid is kept in that grid's `outside_ev`, so every grid sums to the
//!   total deposited energy. Placing the loss at the collision point and the
//!   sub-threshold remainder at the stopping point is this crate's bookkeeping
//!   convention for a point-collision transport loop, not a physics model.
//!   The same deposits optionally feed a radial profile in one depth slab,
//!   on log radial bins and with per-history errors
//!   ([`ElectronTallyConfig::psf`], [`super::psf`]), which is the input of
//!   the point-spread-function fits.
//! - **Ionisation and lattice channels.** Every deposit is also booked to
//!   one of two channels, per layer and per cell of both grids, see
//!   [the channel split](#the-ionisation-and-lattice-channels).
//! - **Emission.** Every electron the `escaped` hook reports is counted at its
//!   face, with its energy and its polar angle from the outward surface
//!   normal. It is classed as *slow* (`E < split`, the secondary-electron
//!   class) or *fast* (`E >= split`, the backscattered class), see
//!   [the split convention](#the-sebse-split). The front-face counts per
//!   primary are the backscatter yield `η` (fast) and the secondary yield `δ`
//!   (slow). The energy spectrum and the polar-angle spectrum (one per class)
//!   are histograms per face.
//! - **Generation-volume moments.** Energy-weighted mean and standard
//!   deviation of the deposition position along `x`, `y`, `z` and of the
//!   radial distance from the beam axis (the line `y = z = 0` through the
//!   entry point), plus the energy-weighted RMS radius. See
//!   [`GenerationVolume`].
//! - **Stopping points.** Unweighted moments of the depth and radial distance
//!   of electrons that fell below the stopping threshold (an electron range
//!   distribution), and the same moments for the primaries alone
//!   ([`PrimaryStoppingPoints`]: the penetration depth of stopped primaries).
//! - **Energy balance.** The incident energy plus the Fermi-sea source equals
//!   the deposited, escaped, trapped and barrier terms, see
//!   [the energy balance](#the-energy-balance) and [`ElectronEnergyBudget`].
//! - **Table coverage.** A numerical diagnostic, not a physical result: per
//!   layer and channel, how many rate evaluations of the transport read the
//!   cross-section table inside its energy grid and how many used the
//!   constant continuation below or above it
//!   ([`super::table_coverage`], [`ElectronReport::table_coverage`]). These
//!   are rate evaluations, not collisions.
//! - **Event-cap counts.** Another numerical diagnostic: how many secondary
//!   tracks the collision cap cut off, and how many primary histories had
//!   any electron (the primary or a secondary) cut off by it
//!   ([`EventCapCounts`], [`ElectronReport::event_caps`]). The primary's own
//!   cap stays in [`FateCounts::event_capped`]; zero there does not mean the
//!   secondary cascades completed.
//!
//! # The SE/BSE split
//!
//! An electron leaving the target with energy below 50 eV is by convention a
//! secondary electron and one at higher energy a backscattered electron:
//! H. Chen, Y. Zou, S. Mao, M. S. S. Khan, K. Tőkési and Z. J. Ding,
//! "Influence of energy loss function to the Monte Carlo simulated electron
//! backscattering coefficient", Sci. Rep. 12, 18201 (2022),
//! doi:10.1038/s41598-022-20466-3 (open access), Introduction ("These
//! electrons, with their energies higher than 50 eV, are defined as
//! backscattered electrons; while the lower energy electrons are the excited
//! secondary electrons", citing H. Niedrig, Scanning 1, 17 (1978)) and section
//! "Monte Carlo simulation" ("an escaped electron is counted either as a true
//! secondary electron (< 50 eV) or a backscattered electron (> 50 eV)"). The
//! class depends only on the energy, not on whether the electron is a primary
//! or a secondary.
//!
//! That paper leaves `E = 50 eV` itself unassigned. This tally puts it in the
//! fast (backscattered) class, as Nebula's detector flags do: `DETECTOR_LT50`
//! ("if energy < 50 eV") and `DETECTOR_GE50` ("if energy >= 50 eV") in
//! `source/core/cpu_material_manager.h`, Nebula commit `a50a8e8`
//! (<https://github.com/Nebula-simulator/Nebula>, BSD-3-Clause). Only the
//! convention is taken from there; no code is ported.
//!
//! The split is configurable ([`ElectronTallyConfig::se_bse_split_ev`],
//! default [`SE_BSE_SPLIT_EV`]) and is recorded, with its source, in
//! [`ElectronTallyMetadata`].
//!
//! # Event kinds and how later channels slot in
//!
//! The tally reads only the hooks of [`ElectronTally`] and does not depend on
//! how many electrons a history contains:
//!
//! - Emission is counted per `escaped` call, so secondaries that escape are
//!   counted with the primaries, as the convention requires.
//! - `stopped` is classified by energy against the stopping threshold of the
//!   electron's layer ([`Transport::stopping_thresholds_ev`]: the cutoff, or
//!   `U + cutoff` under
//!   [`CutoffReference::VacuumLevel`](crate::electron::transport::CutoffReference::VacuumLevel)):
//!   below it the electron stopped and its energy is a deposit, at or above
//!   it the electron is trapped (the loop calls `stopped` for
//!   [`Fate::Trapped`] too, and a trapped electron always has at least the
//!   threshold energy, because the loop checks the threshold after every
//!   energy change, a primary starts above it and a secondary below it is
//!   not created).
//! - `absorbed` (the back face under
//!   [`EscapeRule::FrontOnly`](crate::electron::transport::EscapeRule::FrontOnly))
//!   counts as trapped energy.
//! - [`Fate::EventCap`] has no hook of its own, so the energy of the last
//!   state seen for the capped electron is counted as trapped: at
//!   `end_secondary` for a secondary, at `end_history` for the primary. The
//!   last state is kept per electron; the primary's is set aside when the
//!   first secondary begins, so a secondary cannot overwrite it. Each capped
//!   secondary is also counted once in [`EventCapCounts::secondary_tracks`],
//!   and a history with any capped electron once in
//!   [`EventCapCounts::affected_histories`]; these counts do not touch the
//!   energy accounting.
//!
//! # The energy balance
//!
//! Inside a layer, kinetic energies are measured from the band bottom, and
//! with [`BoundaryModel::StepBarrier`](crate::electron::transport::BoundaryModel::StepBarrier)
//! an electron gains or loses the step `ΔU = U' - U` in inner potential at
//! every face it gets through ([`crate::electron::boundary`]; `U = E_F + Φ`
//! in a metal, Verduin's thesis Eq. 3.136, cited there). The primary's
//! incident energy is its energy where it starts, its vacuum energy with the
//! step barrier. Summed over all histories:
//!
//! ```text
//! incident + fermi_sea = deposited + escaped + trapped + barrier
//! ```
//!
//! - `barrier` is the sum of `-ΔU` over every face transmission (the
//!   `barrier` hook): kinetic energy taken by the potential steps. It is
//!   `-U` for a primary that enters (its vacuum energy becomes `E + U`
//!   inside), `+U` for any electron that leaves through a surface, and
//!   `-(U' - U)` at an interface. A primary reflected at entry changes
//!   nothing and is counted as escaped with its incident energy.
//! - With a secondary model, an inelastic event takes `W` from the electron
//!   and liberates, if anything, an electron of kinetic energy
//!   `E_SE = E_F + W - B` (Verduin Eq. 3.86, see
//!   [`crate::electron::secondary`]), which the tally follows like any other
//!   electron until it stops, escapes or is trapped. Per event
//!   ([`SecondaryEvent`]) `W = secondary + binding + deposited`, with
//!   `binding = B - E_F`. Only `deposited`, plus `binding` where it is
//!   positive (an initial state below the band bottom), stays in the solid
//!   at the event and is deposited there. Where `binding` is negative (a
//!   conduction electron, as in a metal, which already had `E_F - B` above
//!   the band bottom), the liberated electron brings `-binding` of its own:
//!   that is the `fermi_sea` source term. So `W + fermi_sea = secondary +
//!   deposit` event by event, and the secondary's energy is counted once,
//!   where it ends.
//!
//! - With the opt-in insulator channels of a layer
//!   ([`crate::electron::phonon`]), an LO-phonon emission gives `ħω` to the
//!   lattice (deposited where it happens), an absorption takes `ħω` from it
//!   (the `phonon_absorbed` source, added to the left side of the balance),
//!   and an electron trapped as a polaron deposits its whole remaining energy
//!   where it is trapped (the `polaron` deposit). The balance is then
//!   `incident + fermi_sea + phonon_absorbed = deposited + escaped + trapped +
//!   barrier`, with `deposited` including `phonon_emitted` and `polaron`.
//!   Polaron traps are not counted in the stopping points, which are
//!   electrons that fell below the stopping threshold.
//!
//! Without a secondary model and without the step barrier, `fermi_sea` and
//! `barrier` are zero and the balance is `incident = deposited + escaped +
//! trapped`.
//!
//! ## Deposits and the energy imparted
//!
//! Deposits (`deposited`, the per-layer totals and every grid cell) are
//! measured from the band bottom. With a secondary model in a layer with a
//! Fermi energy, a liberated conduction electron's kinetic energy includes
//! the `-binding` it already had (the `fermi_sea` source), and wherever that
//! electron stops, its whole remaining energy above the band bottom is
//! deposited. So deposits include Fermi-sea energy the beam never supplied,
//! and `deposited` can exceed the energy the beam imparted to the target.
//! Moving the sources of the balance above to the right side gives the net
//! energy imparted:
//!
//! ```text
//! incident - escaped = deposited + trapped + barrier - fermi_sea - phonon_absorbed
//! ```
//!
//! Without a secondary model and without the insulator channels the two
//! source terms are zero.
//!
//! # The ionisation and lattice channels
//!
//! The kind of a deposit is known at the call site and carried with it, not
//! reconstructed from totals. Two channels, with `ionisation + lattice =
//! deposited` per layer, per grid cell and in the budget (to floating-point
//! rounding; each is its own accumulator):
//!
//! - **Ionisation** (electronic excitation, the source of electron-hole
//!   pairs): the part of an inelastic event left in the solid (`W` without a
//!   secondary model; `deposited + max(binding, 0)` with one, from the
//!   `secondary` hook), and the **residual** energy of electrons that fell
//!   below the stopping threshold ([`ElectronEnergyBudget::residual_ev`]).
//! - **Lattice**: LO-phonon emission (`phonon_emitted`) and the whole
//!   remaining energy of an electron trapped as a polaron (`polaron`).
//!
//! The residual is the cutoff-dependent term. It is placed in the ionisation
//! channel by convention: a sub-cutoff electron is below the range where this
//! engine follows collisions, and its energy goes into further electronic
//! excitation as it thermalises. This is a bookkeeping convention, not a
//! computed split; the lattice channel is the lower bound of heating, and the
//! ionisation channel includes the residual and, like all deposits, Fermi-sea
//! energy (see [the energy imparted](#deposits-and-the-energy-imparted)).
//! Elastic nuclear recoil energy is not tallied by the electron engine and is
//! in neither channel.
//!
//! # Determinism
//!
//! Counts are integers, histograms have integer bins, and the grids, totals
//! and weighted moments are `f64` accumulators merged field by field in the
//! order [`crate::rng::run_particles`] calls `merge` (chunk order). For a
//! fixed `chunk_size` the report is bit-identical at any thread count
//! (`lindhard/tests/electron_tally.rs`). Nothing is stored in a hash map.

use serde::{Deserialize, Serialize};

use super::dose::{DoseAcc, DoseConfig, DoseError};
use super::exposure::{ExposureAccumulator, ExposureConfig, ExposureError};
use super::hist::Binning;
use super::moments::Moments;
use super::psf::{PsfConfig, PsfError, RadialAccumulator};
use super::table_coverage::TableCoverageTally;
use super::voxel::{VoxelRecord, VoxelReport, VoxelTallyGrid};
use crate::electron::secondary::{SecondaryEvent, SecondaryModel};
use crate::electron::transport::{
    Boundary, ElectronState, ElectronTally, Face, Fate, GridCoverage, PhononEvent, TableChannel,
    Transport,
};
use crate::geometry::Geometry;

mod acc;
mod charge;
mod report;

pub use charge::{
    net_mean_se, CartesianCharge, ChargeCells, ChargeCounts, ChargeReport, CylindricalCharge,
    EscapedCharge,
};
pub use report::{
    CartesianDeposition, CylindricalDeposition, DepositionReport, ElectronEnergyBudget,
    ElectronReport, ElectronTallyMetadata, EmissionClass, EventCapCounts, FaceEmission, FateCounts,
    GenerationVolume, PrimaryStoppingPoints, StoppingPoints, Yields,
};

use acc::{per, Budget, Channel, DepositKind, FaceAcc, GridAcc, VoxelAcc, WeightedMoments};
use charge::ChargeAcc;

/// Default SE/BSE energy split, eV: an escaping electron below this energy is
/// in the secondary (slow) class, at or above it in the backscattered (fast)
/// class. Chen et al., Sci. Rep. 12, 18201 (2022),
/// doi:10.1038/s41598-022-20466-3, Introduction and section "Monte Carlo
/// simulation"; the assignment of the boundary itself follows Nebula's
/// `DETECTOR_GE50` flag. See [the module docs](self#the-sebse-split).
pub const SE_BSE_SPLIT_EV: f64 = 50.0;

/// The source of the default split, as recorded in [`ElectronTallyMetadata`].
pub const SE_BSE_SPLIT_SOURCE: &str = "Chen, Zou, Mao, Khan, Tokesi and Ding, Sci. Rep. 12, 18201 \
     (2022), doi:10.1038/s41598-022-20466-3 (Introduction; section Monte Carlo simulation), citing \
     Niedrig, Scanning 1, 17 (1978); E = split is counted as backscattered, as Nebula's \
     DETECTOR_GE50 flag (source/core/cpu_material_manager.h, commit a50a8e8) does";

/// The classification rule, as recorded in [`ElectronTallyMetadata`].
pub const SE_BSE_SPLIT_RULE: &str =
    "slow (secondary) if E < split; fast (backscattered) if E >= split; E as reported by the escaped hook";

/// Errors from building a [`FullElectronTally`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ElectronTallyError {
    /// The SE/BSE split is not finite and non-negative.
    #[error("invalid SE/BSE split {0} eV: must be finite and non-negative")]
    Split(f64),
    /// The radial binning of a cylindrical grid starts below `r = 0`.
    #[error("the cylindrical grid's radial bins must start at r >= 0, got lo = {0} m")]
    NegativeRadius(f64),
    /// The escape polar-angle binning is not within `[0, π]`.
    #[error("the escape polar-angle bins must lie in [0, pi], got [{lo}, {hi})")]
    PolarRange {
        /// Lower edge.
        lo: f64,
        /// Upper edge.
        hi: f64,
    },
    /// A grid has more cells than fit in memory addressing.
    #[error("the {0} grid has too many cells")]
    GridTooLarge(&'static str),
    /// The radial-profile (PSF) configuration is invalid.
    #[error("invalid PSF configuration: {0}")]
    Psf(#[from] PsfError),
    /// The dose / depth-dose / transmission configuration is invalid.
    #[error(transparent)]
    Dose(#[from] DoseError),
    /// The exposure settings are invalid.
    #[error(transparent)]
    Exposure(#[from] ExposureError),
}

/// A Cartesian deposition grid over `x` (depth), `y` and `z`, m.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CartesianGrid {
    /// Depth bins, m.
    pub x: Binning,
    /// Lateral `y` bins, m.
    pub y: Binning,
    /// Lateral `z` bins, m.
    pub z: Binning,
}

impl CartesianGrid {
    /// Number of voxels, or `None` on overflow.
    pub fn voxels(&self) -> Option<usize> {
        self.x
            .bins
            .checked_mul(self.y.bins)?
            .checked_mul(self.z.bins)
    }

    /// Flat index of voxel `(ix, iy, iz)`: `(ix * ny + iy) * nz + iz`.
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        (ix * self.y.bins + iy) * self.z.bins + iz
    }

    /// Volume of one voxel, m³.
    pub fn voxel_volume_m3(&self) -> f64 {
        self.x.width() * self.y.width() * self.z.width()
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        let ix = self.x.locate(pos[0]).ok()?;
        let iy = self.y.locate(pos[1]).ok()?;
        let iz = self.z.locate(pos[2]).ok()?;
        Some(self.index(ix, iy, iz))
    }
}

/// A cylindrical deposition grid about the beam axis (the line `y = z = 0`
/// through the entry point, along `x`): radial distance `r = (y² + z²)^(1/2)`
/// and depth `x`, m.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CylindricalGrid {
    /// Radial bins, m; `r.lo` must be at least 0.
    pub r: Binning,
    /// Depth bins, m.
    pub depth: Binning,
}

impl CylindricalGrid {
    /// Number of cells, or `None` on overflow.
    pub fn cells(&self) -> Option<usize> {
        self.r.bins.checked_mul(self.depth.bins)
    }

    /// Flat index of cell `(ir, ix)`: `ir * n_depth + ix`.
    pub fn index(&self, ir: usize, ix: usize) -> usize {
        ir * self.depth.bins + ix
    }

    /// Volume of a cell in radial bin `ir` (an annulus times a depth bin), m³.
    pub fn cell_volume_m3(&self, ir: usize) -> f64 {
        let (a, b) = (self.r.edge(ir), self.r.edge(ir + 1));
        std::f64::consts::PI * (b * b - a * a) * self.depth.width()
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        let ir = self.r.locate(pos[1].hypot(pos[2])).ok()?;
        let ix = self.depth.locate(pos[0]).ok()?;
        Some(self.index(ir, ix))
    }
}

/// Settings of a [`FullElectronTally`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronTallyConfig {
    /// Optional Cartesian deposition grid.
    pub cartesian: Option<CartesianGrid>,
    /// Optional cylindrical `r`-`x` deposition grid.
    pub cylindrical: Option<CylindricalGrid>,
    /// Energy of escaping electrons, eV.
    pub escape_energy: Binning,
    /// Polar angle of escaping electrons from the outward surface normal, rad
    /// (`0` along the normal, `π/2` grazing). Must lie within `[0, π]`.
    pub escape_polar: Binning,
    /// SE/BSE split, eV: escaping electrons below it are slow (secondary),
    /// at or above it fast (backscattered). See
    /// [the module docs](self#the-sebse-split).
    pub se_bse_split_ev: f64,
    /// Optional radial profile in a depth slab, log radial bins, with
    /// per-history errors: the input of the PSF fits of [`super::psf`].
    #[serde(default)]
    pub psf: Option<PsfConfig>,
    /// Optional dose per layer, depth-dose profile and interface
    /// transmission ([`super::dose`]).
    #[serde(default)]
    pub dose: Option<DoseConfig>,
    /// Optional exposure statistics: the per-electron deposit in target
    /// volumes and the compound-Poisson distribution of a finite exposure
    /// ([`super::exposure`]).
    #[serde(default)]
    pub exposure: Option<ExposureConfig>,
}

impl ElectronTallyConfig {
    /// Spectra with the given binnings, no deposition grids, no PSF profile,
    /// and the [`SE_BSE_SPLIT_EV`] split.
    pub fn new(escape_energy: Binning, escape_polar: Binning) -> Self {
        Self {
            cartesian: None,
            cylindrical: None,
            escape_energy,
            escape_polar,
            se_bse_split_ev: SE_BSE_SPLIT_EV,
            psf: None,
            dose: None,
            exposure: None,
        }
    }
}

/// Deposition, emission and energy-balance tally for the electron transport
/// loop. See the module docs.
#[derive(Debug, Clone, PartialEq)]
pub struct FullElectronTally {
    config: ElectronTallyConfig,
    cutoff_ev: f64,
    /// Stopping threshold per layer, eV, from the transport.
    thresholds: Vec<f64>,
    /// Whether the transport runs a secondary model (then the `secondary`
    /// hook, not `inelastic`, carries the deposit of an inelastic event).
    secondaries: bool,
    histories: u64,
    fates: FateCounts,
    event_caps: EventCapCounts,
    budget: Budget,
    layer_deposit_ev: Vec<f64>,
    layer_ionisation_ev: Vec<f64>,
    layer_lattice_ev: Vec<f64>,
    cartesian: Option<GridAcc>,
    cylindrical: Option<GridAcc>,
    /// Radial profile in a depth slab with per-history errors.
    psf: Option<RadialAccumulator>,
    /// Dose, depth-dose and transmission accumulators.
    dose: Option<DoseAcc>,
    /// Per-electron deposit in the exposure volumes.
    exposure: Option<ExposureAccumulator>,
    /// Energy-weighted deposition position: `x`, `y`, `z`, `r`.
    generation: [WeightedMoments; 4],
    stop_depth: Moments,
    stop_radial: Moments,
    /// The subset of `stop_depth` / `stop_radial` that are primaries.
    primary_stop_depth: Moments,
    primary_stop_radial: Moments,
    front: FaceAcc,
    back: FaceAcc,
    side: FaceAcc,
    /// Energy deposited per voxel of the optional voxel grid.
    voxel: Option<VoxelAcc>,
    /// Where electrons stopped and holes were created.
    charge: ChargeAcc,
    /// Rate evaluations against each table's energy grid.
    table_coverage: TableCoverageTally,
    /// Per-history scratch, not part of the result and cleared at every
    /// history boundary: the last state seen of the current electron, for
    /// [`Fate::EventCap`] (`None` once it ended through a hook).
    last: Option<ElectronState>,
    /// The primary's last state, set aside when the first secondary begins.
    primary_last: Option<ElectronState>,
    /// Whether the history has reached its secondaries.
    in_secondary: bool,
    /// Whether a secondary of the current history was cut off by the event
    /// cap (so the history is counted once in
    /// [`EventCapCounts::affected_histories`]).
    history_capped: bool,
}

impl FullElectronTally {
    /// Empty tally for runs of `transport` (its per-layer stopping thresholds
    /// classify `stopped` events, its secondary model decides which hook
    /// carries an inelastic deposit, its stack sets the number of layers).
    pub fn new<G: Geometry>(
        transport: &Transport<G>,
        config: ElectronTallyConfig,
    ) -> Result<Self, ElectronTallyError> {
        if !(config.se_bse_split_ev.is_finite() && config.se_bse_split_ev >= 0.0) {
            return Err(ElectronTallyError::Split(config.se_bse_split_ev));
        }
        let p = config.escape_polar;
        if p.lo < 0.0 || p.hi > std::f64::consts::PI {
            return Err(ElectronTallyError::PolarRange { lo: p.lo, hi: p.hi });
        }
        let cartesian = match config.cartesian {
            None => None,
            Some(g) => Some(GridAcc::new(
                g.voxels()
                    .ok_or(ElectronTallyError::GridTooLarge("Cartesian"))?,
            )),
        };
        let cylindrical = match config.cylindrical {
            None => None,
            Some(g) => {
                if g.r.lo < 0.0 {
                    return Err(ElectronTallyError::NegativeRadius(g.r.lo));
                }
                Some(GridAcc::new(
                    g.cells()
                        .ok_or(ElectronTallyError::GridTooLarge("cylindrical"))?,
                ))
            }
        };
        let psf = config.psf.map(RadialAccumulator::new).transpose()?;
        let dose = match config.dose {
            None => None,
            Some(d) => {
                let stack = transport.geometry().as_stack().ok_or_else(|| {
                    DoseError("the dose tally needs a layered stack target".into())
                })?;
                Some(DoseAcc::new(stack, d, config.escape_energy)?)
            }
        };
        let exposure = config
            .exposure
            .clone()
            .map(ExposureAccumulator::new)
            .transpose()?;
        let (front, back, side) = (
            FaceAcc::new(&config),
            FaceAcc::new(&config),
            FaceAcc::new(&config),
        );
        let charge = ChargeAcc::new(
            transport.n_materials(),
            config.cartesian,
            config.cylindrical,
        );
        Ok(Self {
            config,
            cutoff_ev: transport.config().cutoff_ev,
            thresholds: transport.stopping_thresholds_ev().to_vec(),
            secondaries: transport.config().secondaries != SecondaryModel::Off,
            histories: 0,
            fates: FateCounts::default(),
            event_caps: EventCapCounts::default(),
            budget: Budget::default(),
            layer_deposit_ev: vec![0.0; transport.n_materials()],
            layer_ionisation_ev: vec![0.0; transport.n_materials()],
            layer_lattice_ev: vec![0.0; transport.n_materials()],
            cartesian,
            cylindrical,
            psf,
            dose,
            exposure,
            generation: [WeightedMoments::default(); 4],
            stop_depth: Moments::new(),
            stop_radial: Moments::new(),
            primary_stop_depth: Moments::new(),
            primary_stop_radial: Moments::new(),
            front,
            back,
            side,
            voxel: None,
            charge,
            table_coverage: TableCoverageTally::new(transport),
            last: None,
            primary_last: None,
            in_secondary: false,
            history_capped: false,
        })
    }

    /// The configuration.
    pub fn config(&self) -> &ElectronTallyConfig {
        &self.config
    }

    /// Primary histories tallied so far.
    pub fn histories(&self) -> u64 {
        self.histories
    }

    /// Also tally the deposited energy per voxel of `grid`, for
    /// [`FullElectronTally::voxel_report`]. The grid is independent of the
    /// target: use [`VoxelTallyGrid::from_voxel_grid`] to tally on the voxels
    /// of a [`crate::geometry::VoxelGrid`] target.
    pub fn with_voxel_grid(mut self, grid: VoxelTallyGrid) -> Self {
        self.voxel = Some(VoxelAcc::new(grid));
        self
    }

    /// The per-voxel deposition as a [`VoxelReport`] (CSV and VTK export of
    /// [`crate::tally::voxel`]), or `None` without
    /// [`FullElectronTally::with_voxel_grid`].
    ///
    /// Each voxel's energy is split by how it was deposited: `electronic_ev`
    /// is the energy left at inelastic events (the
    /// [`ElectronEnergyBudget::inelastic_ev`] channel), `lattice_ev` the
    /// LO-phonon emission and polaron-trapping energy, and `rest_ev` the
    /// energy of electrons that fell below the stopping threshold
    /// ([`ElectronEnergyBudget::residual_ev`]); their sum over all voxels
    /// and `outside` is [`ElectronEnergyBudget::deposited_ev`].
    /// `implanted` and `displacements` are zero: the electron tally records
    /// neither, and the per-species defect fields (`species_z`, `defects`,
    /// `outside_defects`) are empty: the electron tally has no species axis.
    /// Deposits are point events at the collision (as for
    /// [`CartesianGrid`]), so a flight crossing voxel faces is not split.
    pub fn voxel_report(&self) -> Option<VoxelReport> {
        self.voxel.as_ref().map(|v| VoxelReport {
            grid: v.grid,
            histories: self.histories,
            voxels: v
                .cells
                .iter()
                .map(|(&flat, &cell)| VoxelRecord {
                    index: v.grid.unflatten(flat),
                    cell,
                    defects: Vec::new(),
                    mass_kg: None,
                })
                .collect(),
            outside: v.outside,
            species_z: Vec::new(),
            outside_defects: Vec::new(),
        })
    }

    /// Every energy deposit goes through here. `kind` picks the voxel field
    /// and, through [`DepositKind::channel`], the ionisation or lattice
    /// channel of the layer and grid maps.
    fn deposit(&mut self, at: &ElectronState, e: f64, kind: DepositKind) {
        if e <= 0.0 {
            return;
        }
        let pos = at.pos;
        if let Some(v) = &mut self.voxel {
            v.add(pos, e, kind);
        }
        self.layer_deposit_ev[at.layer] += e;
        let channel = kind.channel();
        match channel {
            Channel::Ionisation => self.layer_ionisation_ev[at.layer] += e,
            Channel::Lattice => self.layer_lattice_ev[at.layer] += e,
        }
        if let (Some(acc), Some(g)) = (&mut self.cartesian, &self.config.cartesian) {
            acc.add(g.locate(pos), e, channel);
        }
        if let (Some(acc), Some(g)) = (&mut self.cylindrical, &self.config.cylindrical) {
            acc.add(g.locate(pos), e, channel);
        }
        if let Some(p) = &mut self.psf {
            p.deposit(pos, e);
        }
        if let Some(d) = &mut self.dose {
            d.deposit(at.layer, pos[0], e);
        }
        if let Some(x) = &mut self.exposure {
            x.deposit(pos, e);
        }
        let r = pos[1].hypot(pos[2]);
        for (m, v) in self.generation.iter_mut().zip([pos[0], pos[1], pos[2], r]) {
            m.push(v, e);
        }
    }

    /// The plain-data result.
    pub fn report(&self) -> ElectronReport {
        let h = self.histories;
        let front = self.front.report(h);
        let back = self.back.report(h);
        let yields = Yields {
            backscatter_eta: front.fast.per_primary,
            secondary_delta: front.slow.per_primary,
            total_sigma: per((front.fast.count + front.slow.count) as f64, h),
            transmitted_fast: back.fast.per_primary,
            transmitted_slow: back.slow.per_primary,
        };
        let [gx, gy, gz, gr] = &self.generation;
        let generation_volume = (gx.weight > 0.0).then(|| GenerationVolume {
            energy_ev: gx.weight,
            mean_m: [gx.mean, gy.mean, gz.mean],
            std_dev_m: [
                gx.variance().sqrt(),
                gy.variance().sqrt(),
                gz.variance().sqrt(),
            ],
            mean_radius_m: gr.mean,
            radius_std_dev_m: gr.variance().sqrt(),
            rms_radius_m: (gy.variance() + gy.mean * gy.mean + gz.variance() + gz.mean * gz.mean)
                .sqrt(),
        });
        let deposition =
            DepositionReport {
                per_layer_ev: self.layer_deposit_ev.clone(),
                per_layer_ionisation_ev: self.layer_ionisation_ev.clone(),
                per_layer_lattice_ev: self.layer_lattice_ev.clone(),
                cartesian: self
                    .config
                    .cartesian
                    .zip(self.cartesian.as_ref())
                    .map(|(grid, a)| CartesianDeposition {
                        grid,
                        energy_ev: a.energy_ev.clone(),
                        outside_ev: a.outside_ev,
                        ionisation_ev: a.ionisation_ev.clone(),
                        lattice_ev: a.lattice_ev.clone(),
                        outside_ionisation_ev: a.outside_ionisation_ev,
                        outside_lattice_ev: a.outside_lattice_ev,
                    }),
                cylindrical: self.config.cylindrical.zip(self.cylindrical.as_ref()).map(
                    |(grid, a)| CylindricalDeposition {
                        grid,
                        energy_ev: a.energy_ev.clone(),
                        outside_ev: a.outside_ev,
                        ionisation_ev: a.ionisation_ev.clone(),
                        lattice_ev: a.lattice_ev.clone(),
                        outside_ionisation_ev: a.outside_ionisation_ev,
                        outside_lattice_ev: a.outside_lattice_ev,
                    },
                ),
                psf: self.psf.as_ref().map(RadialAccumulator::profile),
                exposure: self.exposure.as_ref().map(ExposureAccumulator::report),
            };
        let dose = self.dose.as_ref().map(|d| d.report(self.budget.incident));
        ElectronReport {
            histories: h,
            metadata: ElectronTallyMetadata {
                se_bse_split_ev: self.config.se_bse_split_ev,
                se_bse_split_rule: SE_BSE_SPLIT_RULE.to_string(),
                se_bse_split_source: SE_BSE_SPLIT_SOURCE.to_string(),
                cutoff_ev: self.cutoff_ev,
                stopping_threshold_ev: self.thresholds.clone(),
                layers: self.layer_deposit_ev.len(),
                config: self.config.clone(),
            },
            fates: self.fates,
            event_caps: self.event_caps,
            budget: self.budget.report(),
            yields,
            front,
            back,
            side: (self.side.slow.count + self.side.fast.count > 0).then(|| self.side.report(h)),
            deposition,
            dose,
            generation_volume,
            charge: self.charge.report(self.histories),
            stopping_points: StoppingPoints {
                stopped: self.stop_depth.n,
                depth: self.stop_depth.summary(),
                radial: self.stop_radial.summary(),
                primaries: PrimaryStoppingPoints {
                    stopped: self.primary_stop_depth.n,
                    depth: self.primary_stop_depth.summary(),
                    radial: self.primary_stop_radial.summary(),
                },
            },
            table_coverage: self.table_coverage.layers().to_vec(),
        }
    }
}

impl ElectronTally for FullElectronTally {
    fn begin_history(&mut self, index: u64, start: &ElectronState) {
        if let Some(x) = &mut self.exposure {
            x.begin_history(index, [start.pos[1], start.pos[2]]);
        }
        self.budget.incident += start.energy_ev;
        self.last = Some(*start);
        self.primary_last = None;
        self.in_secondary = false;
        self.history_capped = false;
    }

    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _length_m: f64) {
        self.last = Some(*end);
    }

    fn table_lookup(&mut self, at: &ElectronState, channel: TableChannel, coverage: GridCoverage) {
        self.table_coverage.record(at.layer, channel, coverage);
    }

    fn elastic(&mut self, after: &ElectronState, _theta: f64) {
        self.last = Some(*after);
    }

    fn inelastic(&mut self, after: &ElectronState, w_ev: f64) {
        self.last = Some(*after);
        // With a secondary model the `secondary` hook that follows carries
        // the part of `w_ev` left in the solid.
        if !self.secondaries {
            self.budget.inelastic += w_ev;
            self.deposit(after, w_ev, DepositKind::Electronic);
        }
    }

    fn secondary(
        &mut self,
        primary: &ElectronState,
        event: &SecondaryEvent,
        created: Option<&ElectronState>,
    ) {
        self.last = Some(*primary);
        if event.liberated {
            // A hole is left where the electron was lifted out. A liberated
            // electron that is not followed was below the stopping threshold
            // and is dropped here, so it stops where its hole is.
            self.charge.create(primary, self.in_secondary);
            if created.is_none() {
                self.charge.stop(primary, true);
            }
        }
        // `W = secondary + binding + deposited` (module docs): the created
        // secondary is counted where it ends, a positive binding stays in
        // the solid, a negative one is energy the liberated electron brings.
        let local = event.deposited_ev + event.binding_ev.max(0.0);
        self.budget.inelastic += local;
        self.budget.fermi_sea += (-event.binding_ev).max(0.0);
        self.deposit(primary, local, DepositKind::Electronic);
    }

    fn phonon(
        &mut self,
        after: &ElectronState,
        kind: PhononEvent,
        hbar_omega_ev: f64,
        _theta: f64,
    ) {
        self.last = Some(*after);
        match kind {
            PhononEvent::Emission => {
                self.budget.phonon_emitted += hbar_omega_ev;
                self.deposit(after, hbar_omega_ev, DepositKind::Lattice);
            }
            PhononEvent::Absorption => self.budget.phonon_absorbed += hbar_omega_ev,
        }
    }

    fn polaron_trapped(&mut self, at: &ElectronState) {
        self.last = None;
        self.charge.stop(at, self.in_secondary);
        self.budget.polaron += at.energy_ev;
        self.deposit(at, at.energy_ev, DepositKind::Lattice);
    }

    fn interface(&mut self, at: &ElectronState, from_layer: usize, to_layer: usize) {
        self.last = Some(*at);
        if let Some(d) = &mut self.dose {
            d.crossing(from_layer, to_layer, at.energy_ev, self.in_secondary);
        }
    }

    fn barrier(&mut self, at: &ElectronState, _boundary: Boundary, delta_u_ev: f64) {
        self.last = Some(*at);
        self.budget.barrier -= delta_u_ev;
    }

    fn reflected(&mut self, at: &ElectronState, _boundary: Boundary) {
        self.last = Some(*at);
    }

    fn begin_secondary(&mut self, start: &ElectronState, _generation: u32) {
        if !self.in_secondary {
            // The primary has ended; keep its last state for `end_history`.
            self.primary_last = self.last.take();
            self.in_secondary = true;
        }
        self.last = Some(*start);
    }

    fn end_secondary(&mut self, fate: Fate) {
        if fate == Fate::EventCap {
            if let Some(s) = self.last {
                self.budget.event_cap += s.energy_ev;
            }
            if let Some(s) = self.last {
                self.charge.stop(&s, true);
            }
            self.event_caps.secondary_tracks += 1;
            self.history_capped = true;
        }
        self.last = None;
    }

    fn stopped(&mut self, at: &ElectronState) {
        self.last = None;
        // Below the threshold or trapped without an interaction, the
        // electron ends inside the target either way.
        self.charge.stop(at, self.in_secondary);
        let e = at.energy_ev;
        if e < self.thresholds[at.layer] {
            self.budget.residual += e;
            self.deposit(at, e, DepositKind::Rest);
            let r = at.pos[1].hypot(at.pos[2]);
            self.stop_depth.push(at.pos[0]);
            self.stop_radial.push(r);
            if !self.in_secondary {
                // The primary is followed to its end before any secondary
                // begins, so this is the primary stopping.
                self.primary_stop_depth.push(at.pos[0]);
                self.primary_stop_radial.push(r);
            }
        } else {
            self.budget.no_interaction += e;
        }
    }

    fn escaped(&mut self, at: &ElectronState, face: Face) {
        self.last = None;
        let e = at.energy_ev;
        let polar = at.dir[0].abs().min(1.0).acos();
        self.charge.escaped(face, self.in_secondary);
        let acc = match face {
            Face::Front => {
                self.budget.escaped_front += e;
                &mut self.front
            }
            Face::Back => {
                self.budget.escaped_back += e;
                &mut self.back
            }
            Face::Side => {
                self.budget.escaped_side += e;
                &mut self.side
            }
        };
        acc.energy_hist.fill(e);
        let class = if e < self.config.se_bse_split_ev {
            &mut acc.slow
        } else {
            &mut acc.fast
        };
        class.count += 1;
        class.energy_ev += e;
        class.polar_hist.fill(polar);
    }

    fn absorbed(&mut self, at: &ElectronState) {
        self.last = None;
        self.charge.stop(at, self.in_secondary);
        self.budget.absorbed += at.energy_ev;
    }

    fn end_history(&mut self, _index: u64, fate: Fate) {
        self.histories += 1;
        if let Some(p) = &mut self.psf {
            p.end_history();
        }
        if let Some(d) = &mut self.dose {
            d.end_history();
        }
        if let Some(x) = &mut self.exposure {
            x.end_history();
        }
        match fate {
            Fate::Stopped => self.fates.stopped += 1,
            Fate::Escaped(Face::Front) => self.fates.escaped_front += 1,
            Fate::Escaped(Face::Back) => self.fates.escaped_back += 1,
            Fate::Escaped(Face::Side) => self.fates.escaped_side += 1,
            Fate::Absorbed => self.fates.absorbed += 1,
            Fate::Trapped => self.fates.trapped += 1,
            Fate::PolaronTrapped => self.fates.polaron_trapped += 1,
            Fate::EventCap => {
                self.fates.event_capped += 1;
                let primary = if self.in_secondary {
                    self.primary_last
                } else {
                    self.last
                };
                if let Some(s) = primary {
                    self.budget.event_cap += s.energy_ev;
                    self.charge.stop(&s, false);
                }
            }
        }
        self.charge.end_history();
        if self.history_capped || fate == Fate::EventCap {
            // Once per history, however many of its electrons were capped.
            self.event_caps.affected_histories += 1;
        }
        self.last = None;
        self.primary_last = None;
        self.in_secondary = false;
        self.history_capped = false;
    }

    /// Field-wise merge.
    ///
    /// # Panics
    /// If `other` was built with a different configuration, cutoff,
    /// thresholds, secondary model or stack.
    fn merge(&mut self, o: Self) {
        assert!(
            self.config == o.config
                && self.cutoff_ev == o.cutoff_ev
                && self.thresholds == o.thresholds
                && self.secondaries == o.secondaries
                && self.layer_deposit_ev.len() == o.layer_deposit_ev.len(),
            "merging electron tallies built for different runs"
        );
        self.histories += o.histories;
        self.fates.merge(&o.fates);
        self.event_caps.merge(&o.event_caps);
        self.budget.merge(&o.budget);
        for (a, b) in self.layer_deposit_ev.iter_mut().zip(&o.layer_deposit_ev) {
            *a += b;
        }
        for (a, b) in self
            .layer_ionisation_ev
            .iter_mut()
            .zip(&o.layer_ionisation_ev)
        {
            *a += b;
        }
        for (a, b) in self.layer_lattice_ev.iter_mut().zip(&o.layer_lattice_ev) {
            *a += b;
        }
        if let (Some(a), Some(b)) = (&mut self.cartesian, &o.cartesian) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.cylindrical, &o.cylindrical) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.psf, &o.psf) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.dose, &o.dose) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.exposure, &o.exposure) {
            a.merge(b);
        }
        for (a, b) in self.generation.iter_mut().zip(&o.generation) {
            a.merge(b);
        }
        self.stop_depth.merge(&o.stop_depth);
        self.stop_radial.merge(&o.stop_radial);
        self.primary_stop_depth.merge(&o.primary_stop_depth);
        self.primary_stop_radial.merge(&o.primary_stop_radial);
        self.front.merge(&o.front);
        self.back.merge(&o.back);
        self.side.merge(&o.side);
        match (&mut self.voxel, &o.voxel) {
            (Some(a), Some(b)) => a.merge(b),
            (None, None) => {}
            _ => panic!("merging electron tallies built for different runs"),
        }
        self.charge.merge(&o.charge);
        self.table_coverage.merge_from(&o.table_coverage);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weighted_moments_merge_matches_direct() {
        let xs = [(1.0, 2.0), (3.0, 1.0), (-2.0, 0.5), (7.0, 4.0), (0.5, 3.0)];
        let mut all = WeightedMoments::default();
        for &(x, w) in &xs {
            all.push(x, w);
        }
        let mut a = WeightedMoments::default();
        let mut b = WeightedMoments::default();
        for &(x, w) in &xs[..2] {
            a.push(x, w);
        }
        for &(x, w) in &xs[2..] {
            b.push(x, w);
        }
        a.merge(&b);
        let wsum: f64 = xs.iter().map(|p| p.1).sum();
        let mean = xs.iter().map(|p| p.0 * p.1).sum::<f64>() / wsum;
        let var = xs.iter().map(|p| p.1 * (p.0 - mean).powi(2)).sum::<f64>() / wsum;
        for m in [all, a] {
            assert!((m.weight - wsum).abs() < 1e-12);
            assert!((m.mean - mean).abs() < 1e-12);
            assert!((m.variance() - var).abs() < 1e-12);
        }
        // Zero weights are ignored.
        let before = all;
        all.push(100.0, 0.0);
        assert_eq!(all, before);
    }

    #[test]
    fn grid_indices() {
        let g = CartesianGrid {
            x: Binning::new(0.0, 1.0, 2).unwrap(),
            y: Binning::new(-1.0, 1.0, 3).unwrap(),
            z: Binning::new(-1.0, 1.0, 4).unwrap(),
        };
        assert_eq!(g.voxels(), Some(24));
        assert_eq!(g.locate([0.75, 0.9, -0.9]), Some(g.index(1, 2, 0)));
        assert_eq!(g.index(1, 2, 0), 20);
        assert_eq!(g.locate([1.5, 0.0, 0.0]), None);
        let c = CylindricalGrid {
            r: Binning::new(0.0, 2.0, 2).unwrap(),
            depth: Binning::new(0.0, 1.0, 5).unwrap(),
        };
        assert_eq!(c.locate([0.3, 0.6, 0.8]), Some(c.index(1, 1)));
        let total: f64 = (0..2).map(|i| c.cell_volume_m3(i)).sum::<f64>() * 5.0;
        assert!((total - std::f64::consts::PI * 4.0).abs() < 1e-12);
    }
}
