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
//!   distribution).
//! - **Energy balance.** The incident energy plus the Fermi-sea source equals
//!   the deposited, escaped, trapped and barrier terms, see
//!   [the energy balance](#the-energy-balance) and [`ElectronEnergyBudget`].
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
//!   first secondary begins, so a secondary cannot overwrite it.
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
//! Without a secondary model and without the step barrier, `fermi_sea` and
//! `barrier` are zero and the balance is `incident = deposited + escaped +
//! trapped`.
//!
//! # Determinism
//!
//! Counts are integers, histograms have integer bins, and the grids, totals
//! and weighted moments are `f64` accumulators merged field by field in the
//! order [`crate::rng::run_particles`] calls `merge` (chunk order). For a
//! fixed `chunk_size` the report is bit-identical at any thread count
//! (`lindhard/tests/electron_tally.rs`). Nothing is stored in a hash map.

use serde::{Deserialize, Serialize};

use super::hist::{Binning, Histogram};
use super::moments::{MomentSummary, Moments};
use crate::electron::secondary::{SecondaryEvent, SecondaryModel};
use crate::electron::transport::{Boundary, ElectronState, ElectronTally, Face, Fate, Transport};

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
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
}

impl ElectronTallyConfig {
    /// Spectra with the given binnings, no deposition grids, and the
    /// [`SE_BSE_SPLIT_EV`] split.
    pub fn new(escape_energy: Binning, escape_polar: Binning) -> Self {
        Self {
            cartesian: None,
            cylindrical: None,
            escape_energy,
            escape_polar,
            se_bse_split_ev: SE_BSE_SPLIT_EV,
        }
    }
}

/// Energy-weighted mean and centred second moment of one coordinate,
/// mergeable. The pairwise update is that of T. F. Chan, G. H. Golub and
/// R. J. LeVeque, Stanford report STAN-CS-79-773 (1979) (the update also used
/// by [`super::moments`]), with the counts replaced by weights, which is the
/// weighted form of D. H. D. West, "Updating mean and variance estimates: an
/// improved method", Commun. ACM 22, 532 (1979), doi:10.1145/359146.359153.
/// With `W = W_a + W_b` and `δ = mean_b - mean_a`:
/// `mean = mean_a + δ W_b / W`, `M2 = M2_a + M2_b + δ² W_a W_b / W`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct WeightedMoments {
    weight: f64,
    mean: f64,
    m2: f64,
}

impl WeightedMoments {
    fn push(&mut self, x: f64, w: f64) {
        self.merge_parts(w, x, 0.0);
    }

    fn merge(&mut self, o: &WeightedMoments) {
        self.merge_parts(o.weight, o.mean, o.m2);
    }

    fn merge_parts(&mut self, w_b: f64, mean_b: f64, m2_b: f64) {
        if w_b <= 0.0 {
            return;
        }
        if self.weight <= 0.0 {
            *self = Self {
                weight: w_b,
                mean: mean_b,
                m2: m2_b,
            };
            return;
        }
        let w = self.weight + w_b;
        let d = mean_b - self.mean;
        self.m2 += m2_b + d * d * self.weight * w_b / w;
        self.mean += d * w_b / w;
        self.weight = w;
    }

    fn variance(&self) -> f64 {
        if self.weight > 0.0 {
            (self.m2 / self.weight).max(0.0)
        } else {
            0.0
        }
    }
}

/// Energy sums, eV.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Budget {
    incident: f64,
    inelastic: f64,
    residual: f64,
    escaped_front: f64,
    escaped_back: f64,
    no_interaction: f64,
    absorbed: f64,
    event_cap: f64,
    fermi_sea: f64,
    barrier: f64,
}

impl Budget {
    fn merge(&mut self, o: &Budget) {
        self.incident += o.incident;
        self.inelastic += o.inelastic;
        self.residual += o.residual;
        self.escaped_front += o.escaped_front;
        self.escaped_back += o.escaped_back;
        self.no_interaction += o.no_interaction;
        self.absorbed += o.absorbed;
        self.event_cap += o.event_cap;
        self.fermi_sea += o.fermi_sea;
        self.barrier += o.barrier;
    }

    fn report(&self) -> ElectronEnergyBudget {
        let deposited = self.inelastic + self.residual;
        let escaped = self.escaped_front + self.escaped_back;
        let trapped = self.no_interaction + self.absorbed + self.event_cap;
        let source = self.incident + self.fermi_sea;
        let sink = deposited + escaped + trapped + self.barrier;
        let relative_imbalance = if source > 0.0 {
            (sink - source).abs() / source
        } else {
            0.0
        };
        ElectronEnergyBudget {
            incident_ev: self.incident,
            deposited_ev: deposited,
            escaped_ev: escaped,
            trapped_ev: trapped,
            inelastic_ev: self.inelastic,
            residual_ev: self.residual,
            escaped_front_ev: self.escaped_front,
            escaped_back_ev: self.escaped_back,
            no_interaction_ev: self.no_interaction,
            absorbed_ev: self.absorbed,
            event_cap_ev: self.event_cap,
            fermi_sea_ev: self.fermi_sea,
            barrier_ev: self.barrier,
            relative_imbalance,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct GridAcc {
    energy_ev: Vec<f64>,
    outside_ev: f64,
}

impl GridAcc {
    fn new(cells: usize) -> Self {
        Self {
            energy_ev: vec![0.0; cells],
            outside_ev: 0.0,
        }
    }

    fn add(&mut self, cell: Option<usize>, e: f64) {
        match cell {
            Some(i) => self.energy_ev[i] += e,
            None => self.outside_ev += e,
        }
    }

    fn merge(&mut self, o: &GridAcc) {
        for (a, b) in self.energy_ev.iter_mut().zip(&o.energy_ev) {
            *a += b;
        }
        self.outside_ev += o.outside_ev;
    }
}

#[derive(Debug, Clone, PartialEq)]
struct ClassAcc {
    count: u64,
    energy_ev: f64,
    polar_hist: Histogram,
}

impl ClassAcc {
    fn new(c: &ElectronTallyConfig) -> Self {
        Self {
            count: 0,
            energy_ev: 0.0,
            polar_hist: Histogram::new(c.escape_polar),
        }
    }

    fn merge(&mut self, o: &ClassAcc) {
        self.count += o.count;
        self.energy_ev += o.energy_ev;
        self.polar_hist.merge(&o.polar_hist);
    }

    fn report(&self, histories: u64) -> EmissionClass {
        EmissionClass {
            count: self.count,
            per_primary: per(self.count as f64, histories),
            energy_ev: self.energy_ev,
            mean_energy_ev: per(self.energy_ev, self.count),
            polar_histogram: self.polar_hist.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct FaceAcc {
    slow: ClassAcc,
    fast: ClassAcc,
    energy_hist: Histogram,
}

impl FaceAcc {
    fn new(c: &ElectronTallyConfig) -> Self {
        Self {
            slow: ClassAcc::new(c),
            fast: ClassAcc::new(c),
            energy_hist: Histogram::new(c.escape_energy),
        }
    }

    fn merge(&mut self, o: &FaceAcc) {
        self.slow.merge(&o.slow);
        self.fast.merge(&o.fast);
        self.energy_hist.merge(&o.energy_hist);
    }

    fn report(&self, histories: u64) -> FaceEmission {
        let count = self.slow.count + self.fast.count;
        FaceEmission {
            count,
            per_primary: per(count as f64, histories),
            energy_ev: self.slow.energy_ev + self.fast.energy_ev,
            energy_histogram: self.energy_hist.clone(),
            slow: self.slow.report(histories),
            fast: self.fast.report(histories),
        }
    }
}

fn per(x: f64, n: u64) -> f64 {
    if n == 0 {
        0.0
    } else {
        x / n as f64
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
    budget: Budget,
    layer_deposit_ev: Vec<f64>,
    cartesian: Option<GridAcc>,
    cylindrical: Option<GridAcc>,
    /// Energy-weighted deposition position: `x`, `y`, `z`, `r`.
    generation: [WeightedMoments; 4],
    stop_depth: Moments,
    stop_radial: Moments,
    front: FaceAcc,
    back: FaceAcc,
    /// Per-history scratch, not part of the result and cleared at every
    /// history boundary: the last state seen of the current electron, for
    /// [`Fate::EventCap`] (`None` once it ended through a hook).
    last: Option<ElectronState>,
    /// The primary's last state, set aside when the first secondary begins.
    primary_last: Option<ElectronState>,
    /// Whether the history has reached its secondaries.
    in_secondary: bool,
}

impl FullElectronTally {
    /// Empty tally for runs of `transport` (its per-layer stopping thresholds
    /// classify `stopped` events, its secondary model decides which hook
    /// carries an inelastic deposit, its stack sets the number of layers).
    pub fn new(
        transport: &Transport,
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
        Ok(Self {
            config,
            cutoff_ev: transport.config().cutoff_ev,
            thresholds: transport.stopping_thresholds_ev().to_vec(),
            secondaries: transport.config().secondaries != SecondaryModel::Off,
            histories: 0,
            fates: FateCounts::default(),
            budget: Budget::default(),
            layer_deposit_ev: vec![0.0; transport.stack().layers().len()],
            cartesian,
            cylindrical,
            generation: [WeightedMoments::default(); 4],
            stop_depth: Moments::new(),
            stop_radial: Moments::new(),
            front: FaceAcc::new(&config),
            back: FaceAcc::new(&config),
            last: None,
            primary_last: None,
            in_secondary: false,
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

    /// Every energy deposit goes through here.
    fn deposit(&mut self, at: &ElectronState, e: f64) {
        if e <= 0.0 {
            return;
        }
        let pos = at.pos;
        self.layer_deposit_ev[at.layer] += e;
        if let (Some(acc), Some(g)) = (&mut self.cartesian, &self.config.cartesian) {
            acc.add(g.locate(pos), e);
        }
        if let (Some(acc), Some(g)) = (&mut self.cylindrical, &self.config.cylindrical) {
            acc.add(g.locate(pos), e);
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
                cartesian: self
                    .config
                    .cartesian
                    .zip(self.cartesian.as_ref())
                    .map(|(grid, a)| CartesianDeposition {
                        grid,
                        energy_ev: a.energy_ev.clone(),
                        outside_ev: a.outside_ev,
                    }),
                cylindrical: self.config.cylindrical.zip(self.cylindrical.as_ref()).map(
                    |(grid, a)| CylindricalDeposition {
                        grid,
                        energy_ev: a.energy_ev.clone(),
                        outside_ev: a.outside_ev,
                    },
                ),
            };
        ElectronReport {
            histories: h,
            metadata: ElectronTallyMetadata {
                se_bse_split_ev: self.config.se_bse_split_ev,
                se_bse_split_rule: SE_BSE_SPLIT_RULE.to_string(),
                se_bse_split_source: SE_BSE_SPLIT_SOURCE.to_string(),
                cutoff_ev: self.cutoff_ev,
                stopping_threshold_ev: self.thresholds.clone(),
                layers: self.layer_deposit_ev.len(),
                config: self.config,
            },
            fates: self.fates,
            budget: self.budget.report(),
            yields,
            front,
            back,
            deposition,
            generation_volume,
            stopping_points: StoppingPoints {
                stopped: self.stop_depth.n,
                depth: self.stop_depth.summary(),
                radial: self.stop_radial.summary(),
            },
        }
    }
}

impl ElectronTally for FullElectronTally {
    fn begin_history(&mut self, _index: u64, start: &ElectronState) {
        self.budget.incident += start.energy_ev;
        self.last = Some(*start);
        self.primary_last = None;
        self.in_secondary = false;
    }

    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _length_m: f64) {
        self.last = Some(*end);
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
            self.deposit(after, w_ev);
        }
    }

    fn secondary(
        &mut self,
        primary: &ElectronState,
        event: &SecondaryEvent,
        _created: Option<&ElectronState>,
    ) {
        self.last = Some(*primary);
        // `W = secondary + binding + deposited` (module docs): the created
        // secondary is counted where it ends, a positive binding stays in
        // the solid, a negative one is energy the liberated electron brings.
        let local = event.deposited_ev + event.binding_ev.max(0.0);
        self.budget.inelastic += local;
        self.budget.fermi_sea += (-event.binding_ev).max(0.0);
        self.deposit(primary, local);
    }

    fn interface(&mut self, at: &ElectronState, _from_layer: usize, _to_layer: usize) {
        self.last = Some(*at);
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
        }
        self.last = None;
    }

    fn stopped(&mut self, at: &ElectronState) {
        self.last = None;
        let e = at.energy_ev;
        if e < self.thresholds[at.layer] {
            self.budget.residual += e;
            self.deposit(at, e);
            self.stop_depth.push(at.pos[0]);
            self.stop_radial.push(at.pos[1].hypot(at.pos[2]));
        } else {
            self.budget.no_interaction += e;
        }
    }

    fn escaped(&mut self, at: &ElectronState, face: Face) {
        self.last = None;
        let e = at.energy_ev;
        let polar = at.dir[0].abs().min(1.0).acos();
        let acc = match face {
            Face::Front => {
                self.budget.escaped_front += e;
                &mut self.front
            }
            Face::Back => {
                self.budget.escaped_back += e;
                &mut self.back
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
        self.budget.absorbed += at.energy_ev;
    }

    fn end_history(&mut self, _index: u64, fate: Fate) {
        self.histories += 1;
        match fate {
            Fate::Stopped => self.fates.stopped += 1,
            Fate::Escaped(Face::Front) => self.fates.escaped_front += 1,
            Fate::Escaped(Face::Back) => self.fates.escaped_back += 1,
            Fate::Absorbed => self.fates.absorbed += 1,
            Fate::Trapped => self.fates.trapped += 1,
            Fate::EventCap => {
                self.fates.event_capped += 1;
                let primary = if self.in_secondary {
                    self.primary_last
                } else {
                    self.last
                };
                if let Some(s) = primary {
                    self.budget.event_cap += s.energy_ev;
                }
            }
        }
        self.last = None;
        self.primary_last = None;
        self.in_secondary = false;
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
        self.budget.merge(&o.budget);
        for (a, b) in self.layer_deposit_ev.iter_mut().zip(&o.layer_deposit_ev) {
            *a += b;
        }
        if let (Some(a), Some(b)) = (&mut self.cartesian, &o.cartesian) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.cylindrical, &o.cylindrical) {
            a.merge(b);
        }
        for (a, b) in self.generation.iter_mut().zip(&o.generation) {
            a.merge(b);
        }
        self.stop_depth.merge(&o.stop_depth);
        self.stop_radial.merge(&o.stop_radial);
        self.front.merge(&o.front);
        self.back.merge(&o.back);
    }
}

/// Everything a [`FullElectronTally`] measured, as plain data.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronReport {
    /// Primary histories.
    pub histories: u64,
    /// The tally settings and conventions.
    pub metadata: ElectronTallyMetadata,
    /// How the primary histories ended.
    pub fates: FateCounts,
    /// Summed energy balance of all histories, eV.
    pub budget: ElectronEnergyBudget,
    /// Emission yields per primary.
    pub yields: Yields,
    /// Electrons leaving the front face (`x = 0`).
    pub front: FaceEmission,
    /// Electrons leaving the back face of a finite stack.
    pub back: FaceEmission,
    /// Where the energy was deposited.
    pub deposition: DepositionReport,
    /// Energy-weighted moments of the deposition position; `None` if nothing
    /// was deposited.
    pub generation_volume: Option<GenerationVolume>,
    /// Where electrons fell below the stopping threshold.
    pub stopping_points: StoppingPoints,
}

/// The settings and conventions a report was made with.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronTallyMetadata {
    /// SE/BSE split, eV.
    pub se_bse_split_ev: f64,
    /// How an escaping electron is classed against the split.
    pub se_bse_split_rule: String,
    /// Where the split convention comes from.
    pub se_bse_split_source: String,
    /// Transport cutoff, eV, as configured (measured from the band bottom or
    /// the vacuum level, see the run's
    /// [`RunMetadata`](crate::electron::transport::RunMetadata)).
    pub cutoff_ev: f64,
    /// Stopping threshold per layer, eV: a `stopped` electron below its
    /// layer's threshold is a deposit, at or above it a trapped electron.
    pub stopping_threshold_ev: Vec<f64>,
    /// Layers in the stack.
    pub layers: usize,
    /// The full tally configuration (grids and spectrum binnings).
    pub config: ElectronTallyConfig,
}

/// How the primary histories ended (one per history, from its
/// [`Fate`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct FateCounts {
    /// Fell below the stopping threshold in the target.
    pub stopped: u64,
    /// Left through the front face.
    pub escaped_front: u64,
    /// Left through the back face.
    pub escaped_back: u64,
    /// Absorbed at the back face.
    pub absorbed: u64,
    /// No interaction available and no face to reach.
    pub trapped: u64,
    /// Cut off by the collision cap.
    pub event_capped: u64,
}

impl FateCounts {
    fn merge(&mut self, o: &FateCounts) {
        self.stopped += o.stopped;
        self.escaped_front += o.escaped_front;
        self.escaped_back += o.escaped_back;
        self.absorbed += o.absorbed;
        self.trapped += o.trapped;
        self.event_capped += o.event_capped;
    }
}

/// The energy balance, eV, summed over all histories and every electron of
/// them (see [the module docs](self#the-energy-balance)):
/// `incident + fermi_sea = deposited + escaped + trapped + barrier`, with
///
/// - `deposited = inelastic + residual`: energy left in the solid at
///   inelastic events, plus the energy electrons had left when they fell
///   below the stopping threshold;
/// - `escaped = escaped_front + escaped_back`;
/// - `trapped = no_interaction + absorbed + event_cap`: energy still carried
///   by electrons that ended in the target without falling below the
///   threshold (no interaction available, absorbed at the back face, or cut
///   off by the collision cap);
/// - `fermi_sea`: the kinetic energy liberated secondaries already had (a
///   source; zero without a secondary model);
/// - `barrier`: the kinetic energy taken by the potential steps at the faces,
///   the sum of `-ΔU` (signed; zero with transparent faces).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElectronEnergyBudget {
    /// Energy of the primaries.
    pub incident_ev: f64,
    /// Deposited in the target.
    pub deposited_ev: f64,
    /// Carried out of the target.
    pub escaped_ev: f64,
    /// Carried by electrons that ended in the target at or above the
    /// stopping threshold.
    pub trapped_ev: f64,
    /// Energy left in the solid at inelastic events (part of
    /// `deposited_ev`): the loss `W` without a secondary model; with one,
    /// [`SecondaryEvent::deposited_ev`] plus any positive
    /// [`SecondaryEvent::binding_ev`].
    pub inelastic_ev: f64,
    /// Remaining energy of electrons that fell below the stopping threshold
    /// (part of `deposited_ev`).
    pub residual_ev: f64,
    /// Out of the front face (part of `escaped_ev`).
    pub escaped_front_ev: f64,
    /// Out of the back face (part of `escaped_ev`).
    pub escaped_back_ev: f64,
    /// Electrons with no interaction and no face to reach (part of
    /// `trapped_ev`).
    pub no_interaction_ev: f64,
    /// Electrons absorbed at the back face (part of `trapped_ev`).
    pub absorbed_ev: f64,
    /// Electrons cut off by the collision cap (part of `trapped_ev`).
    pub event_cap_ev: f64,
    /// Energy the liberated secondaries brought with them, the sum of
    /// `-binding_ev` over events where it is negative (a source).
    pub fermi_sea_ev: f64,
    /// Kinetic energy taken by the potential steps, the sum of `-ΔU` over
    /// face transmissions (signed).
    pub barrier_ev: f64,
    /// `|deposited + escaped + trapped + barrier - incident - fermi_sea| /
    /// (incident + fermi_sea)` (0 with no incident energy).
    pub relative_imbalance: f64,
}

/// Emission yields per primary.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Yields {
    /// Backscatter yield `η`: fast electrons out of the front face per
    /// primary.
    pub backscatter_eta: f64,
    /// Secondary yield `δ`: slow electrons out of the front face per primary.
    pub secondary_delta: f64,
    /// Total yield `σ = η + δ`.
    pub total_sigma: f64,
    /// Fast electrons out of the back face per primary.
    pub transmitted_fast: f64,
    /// Slow electrons out of the back face per primary.
    pub transmitted_slow: f64,
}

/// Electrons leaving through one face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FaceEmission {
    /// Electrons.
    pub count: u64,
    /// Electrons per primary.
    pub per_primary: f64,
    /// Sum of their energies, eV.
    pub energy_ev: f64,
    /// Energy spectrum, eV, both classes.
    pub energy_histogram: Histogram,
    /// Slow electrons (`E < split`, the secondary class).
    pub slow: EmissionClass,
    /// Fast electrons (`E >= split`, the backscattered class).
    pub fast: EmissionClass,
}

/// One energy class of the electrons leaving through a face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EmissionClass {
    /// Electrons.
    pub count: u64,
    /// Electrons per primary.
    pub per_primary: f64,
    /// Sum of their energies, eV.
    pub energy_ev: f64,
    /// Mean energy, eV (0 if none).
    pub mean_energy_ev: f64,
    /// Polar angle from the outward surface normal, rad.
    pub polar_histogram: Histogram,
}

/// Deposited energy by layer and on the configured grids, eV.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DepositionReport {
    /// By the layer the deposit was made in.
    pub per_layer_ev: Vec<f64>,
    /// On the Cartesian grid, if configured.
    pub cartesian: Option<CartesianDeposition>,
    /// On the cylindrical grid, if configured.
    pub cylindrical: Option<CylindricalDeposition>,
}

/// Deposited energy on a [`CartesianGrid`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CartesianDeposition {
    /// The grid.
    pub grid: CartesianGrid,
    /// Energy per voxel, eV, at [`CartesianGrid::index`].
    pub energy_ev: Vec<f64>,
    /// Energy deposited outside the grid, eV.
    pub outside_ev: f64,
}

impl CartesianDeposition {
    /// Energy in voxel `(ix, iy, iz)`, eV.
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.energy_ev[self.grid.index(ix, iy, iz)]
    }

    /// Energy per volume in each voxel, eV/m³.
    pub fn density_ev_per_m3(&self) -> Vec<f64> {
        let v = self.grid.voxel_volume_m3();
        self.energy_ev.iter().map(|e| e / v).collect()
    }
}

/// Deposited energy on a [`CylindricalGrid`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CylindricalDeposition {
    /// The grid.
    pub grid: CylindricalGrid,
    /// Energy per cell, eV, at [`CylindricalGrid::index`].
    pub energy_ev: Vec<f64>,
    /// Energy deposited outside the grid, eV.
    pub outside_ev: f64,
}

impl CylindricalDeposition {
    /// Energy in cell `(ir, ix)`, eV.
    pub fn at(&self, ir: usize, ix: usize) -> f64 {
        self.energy_ev[self.grid.index(ir, ix)]
    }

    /// Energy per volume in each cell, eV/m³.
    pub fn density_ev_per_m3(&self) -> Vec<f64> {
        self.energy_ev
            .iter()
            .enumerate()
            .map(|(i, e)| e / self.grid.cell_volume_m3(i / self.grid.depth.bins))
            .collect()
    }
}

/// Energy-weighted moments of the deposition position (weights: the energy
/// of each deposit). Lengths in metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GenerationVolume {
    /// Total weight: the deposited energy, eV.
    pub energy_ev: f64,
    /// Weighted mean of `[x, y, z]`.
    pub mean_m: [f64; 3],
    /// Weighted standard deviation of `[x, y, z]`.
    pub std_dev_m: [f64; 3],
    /// Weighted mean distance from the beam axis.
    pub mean_radius_m: f64,
    /// Weighted standard deviation of the distance from the beam axis.
    pub radius_std_dev_m: f64,
    /// Weighted RMS distance from the beam axis, `<y² + z²>^(1/2)`.
    pub rms_radius_m: f64,
}

/// Where electrons fell below the stopping threshold (unweighted, one sample
/// each).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoppingPoints {
    /// Electrons that fell below the stopping threshold.
    pub stopped: u64,
    /// Depth moments, m; `None` with fewer than two.
    pub depth: Option<MomentSummary>,
    /// Moments of the distance from the beam axis, m.
    pub radial: Option<MomentSummary>,
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
