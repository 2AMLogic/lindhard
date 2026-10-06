//! What the BCA engine reports, and a small built-in summary tally.
//!
//! The engine calls the [`BcaTally`] hooks as events happen. Every hook has a
//! no-op default, so a tally implements only what it needs; the depth,
//! lateral and damage tallies of a later milestone plug in here without any
//! change to the engine. Tallies are merged in chunk order by
//! [`crate::rng::run_particles`], so a tally whose `merge` is a plain sum is
//! bit-identical at any thread count.

use super::Particle;

/// Which face of the target a particle left through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Face {
    /// The front surface, `x = 0` (backscattering, sputtering).
    Front,
    /// The back face of a finite stack (transmission).
    Back,
}

/// Where an amount of nuclear energy given to the lattice came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatticeDeposit {
    /// A transfer `T <= E_d` (or `T <= E_b`): no atom is displaced and all of
    /// `T` stays at the collision site.
    Subthreshold,
    /// The lattice binding energy `E_b` subtracted from a displaced recoil.
    Binding,
    /// The transfer of a weak collision ([`super::BcaConfig::weak_collisions`]).
    /// Weak collisions make no recoils, whatever `T` is (Moller and Eckstein,
    /// IPP 9/64 (1988), p. 26), so all of `T` stays at the collision site.
    Weak,
}

/// Electronic energy-loss channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElectronicChannel {
    /// Continuous (nonlocal) loss along a free-flight segment.
    NonLocal,
    /// Local loss at a collision (Oen-Robinson).
    Local,
}

/// Energy bookkeeping for one primary history, eV.
///
/// Every eV of the incident energy ends in exactly one field, so
/// [`EnergyBudget::residual`] is zero up to floating-point rounding. A
/// surface barrier `E_s` paid by an escaping particle is counted in
/// `surface_barrier`, not in the escaped energy (which is the energy outside
/// the target).
#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct EnergyBudget {
    /// Incident kinetic energy of the primary.
    pub incident: f64,
    /// Nonlocal electronic loss, all particles.
    pub electronic_nonlocal: f64,
    /// Local electronic loss at collisions, all particles.
    pub electronic_local: f64,
    /// Nuclear energy left in the lattice: subthreshold transfers, weak-
    /// collision transfers, and `E_b` of displaced recoils.
    pub lattice: f64,
    /// Work done against the planar surface barrier by escaping particles.
    pub surface_barrier: f64,
    /// Energy carried out of the front face by the primary.
    pub backscattered: f64,
    /// Energy carried out of the front face by recoils (sputtered atoms).
    pub sputtered: f64,
    /// Energy carried out of the back face (primary and recoils).
    pub transmitted: f64,
    /// Kinetic energy of particles when they stopped (below their cutoff, or
    /// recoils not followed).
    pub rest: f64,
}

impl EnergyBudget {
    /// `incident` minus the sum of all other fields.
    pub fn residual(&self) -> f64 {
        self.incident
            - (self.electronic_nonlocal
                + self.electronic_local
                + self.lattice
                + self.surface_barrier
                + self.backscattered
                + self.sputtered
                + self.transmitted
                + self.rest)
    }

    /// Field-wise sum (for merging).
    pub fn add(&mut self, o: &EnergyBudget) {
        self.incident += o.incident;
        self.electronic_nonlocal += o.electronic_nonlocal;
        self.electronic_local += o.electronic_local;
        self.lattice += o.lattice;
        self.surface_barrier += o.surface_barrier;
        self.backscattered += o.backscattered;
        self.sputtered += o.sputtered;
        self.transmitted += o.transmitted;
        self.rest += o.rest;
    }
}

/// Event hooks called by the engine. All default to no-ops except `merge`.
///
/// Positions are `[x, y, z]` in metres, `x` the depth. The `Particle` passed
/// is the state **after** the event (for example after the energy loss).
pub trait BcaTally: Send {
    /// A primary history starts.
    fn begin_history(&mut self, _index: u64) {}

    /// Electronic loss of `energy_ev`. For [`ElectronicChannel::NonLocal`]
    /// it was spread along the straight segment from `from` to `p.pos`; for
    /// [`ElectronicChannel::Local`] it was deposited at `p.pos`.
    fn electronic(
        &mut self,
        _p: &Particle,
        _from: [f64; 3],
        _channel: ElectronicChannel,
        _energy_ev: f64,
    ) {
    }

    /// Nuclear energy `energy_ev` left in the lattice at `at` in layer `layer`.
    fn lattice(&mut self, _at: [f64; 3], _layer: usize, _kind: LatticeDeposit, _energy_ev: f64) {}

    /// A target atom was displaced (transfer above `E_d`); `recoil` is its
    /// state at creation, before it is followed (or stopped, if below the
    /// recoil cutoff or if cascades are off).
    fn recoil(&mut self, _recoil: &Particle) {}

    /// A particle came to rest at `p.pos` with kinetic energy `p.energy_ev`.
    fn stopped(&mut self, _p: &Particle) {}

    /// A particle left the target through `face`; `p` is its state outside
    /// (after refraction by the surface barrier).
    fn escaped(&mut self, _p: &Particle, _face: Face) {}

    /// A primary history ended; `budget` is its energy bookkeeping.
    fn end_history(&mut self, _index: u64, _budget: &EnergyBudget) {}

    /// Fold `other` (a later chunk) into `self`.
    fn merge(&mut self, other: Self)
    where
        Self: Sized;
}

/// A minimal summary: counts, energy sums and a depth histogram of where
/// primaries stop. Enough to test the engine; the full tallies come later.
#[derive(Debug, Clone, PartialEq)]
pub struct SummaryTally {
    /// Primary histories run.
    pub histories: u64,
    /// Sum of the per-history energy budgets.
    pub budget: EnergyBudget,
    /// Largest `|residual| / incident` over all histories.
    pub max_relative_residual: f64,
    /// Primaries that came to rest in the target.
    pub primaries_stopped: u64,
    /// Primaries that left through the front face.
    pub backscattered: u64,
    /// Primaries that left through the back face.
    pub transmitted: u64,
    /// Recoils that left through the front face.
    pub sputtered: u64,
    /// Recoils that left through the back face.
    pub recoils_transmitted: u64,
    /// Displaced atoms (recoils created).
    pub recoils: u64,
    /// Sum of the rest depth of stopped primaries, m.
    pub depth_sum: f64,
    /// Sum of squared rest depths of stopped primaries, m².
    pub depth_sq_sum: f64,
    /// Width of a depth bin, m.
    pub bin_width_m: f64,
    /// Stopped-primary counts per depth bin; the last bin collects overflow.
    pub depth_hist: Vec<u64>,
}

impl SummaryTally {
    /// Empty tally with `n_bins` depth bins of width `bin_width_m`.
    pub fn new(bin_width_m: f64, n_bins: usize) -> Self {
        Self {
            histories: 0,
            budget: EnergyBudget::default(),
            max_relative_residual: 0.0,
            primaries_stopped: 0,
            backscattered: 0,
            transmitted: 0,
            sputtered: 0,
            recoils_transmitted: 0,
            recoils: 0,
            depth_sum: 0.0,
            depth_sq_sum: 0.0,
            bin_width_m,
            depth_hist: vec![0; n_bins.max(1)],
        }
    }

    /// Mean rest depth of stopped primaries (the projected range), m.
    pub fn mean_depth(&self) -> f64 {
        self.depth_sum / self.primaries_stopped as f64
    }

    /// Standard deviation of the rest depth of stopped primaries, m.
    pub fn depth_std(&self) -> f64 {
        let n = self.primaries_stopped as f64;
        let m = self.depth_sum / n;
        (self.depth_sq_sum / n - m * m).max(0.0).sqrt()
    }
}

impl BcaTally for SummaryTally {
    fn recoil(&mut self, _recoil: &Particle) {
        self.recoils += 1;
    }

    fn stopped(&mut self, p: &Particle) {
        if p.is_primary() {
            self.primaries_stopped += 1;
            let x = p.pos[0];
            self.depth_sum += x;
            self.depth_sq_sum += x * x;
            let last = self.depth_hist.len() - 1;
            let bin = ((x / self.bin_width_m).max(0.0) as usize).min(last);
            self.depth_hist[bin] += 1;
        }
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        match (p.is_primary(), face) {
            (true, Face::Front) => self.backscattered += 1,
            (true, Face::Back) => self.transmitted += 1,
            (false, Face::Front) => self.sputtered += 1,
            (false, Face::Back) => self.recoils_transmitted += 1,
        }
    }

    fn end_history(&mut self, _index: u64, budget: &EnergyBudget) {
        self.histories += 1;
        self.budget.add(budget);
        let r = (budget.residual() / budget.incident).abs();
        self.max_relative_residual = self.max_relative_residual.max(r);
    }

    fn merge(&mut self, o: Self) {
        self.histories += o.histories;
        self.budget.add(&o.budget);
        self.max_relative_residual = self.max_relative_residual.max(o.max_relative_residual);
        self.primaries_stopped += o.primaries_stopped;
        self.backscattered += o.backscattered;
        self.transmitted += o.transmitted;
        self.sputtered += o.sputtered;
        self.recoils_transmitted += o.recoils_transmitted;
        self.recoils += o.recoils;
        self.depth_sum += o.depth_sum;
        self.depth_sq_sum += o.depth_sq_sum;
        for (a, b) in self.depth_hist.iter_mut().zip(o.depth_hist) {
            *a += b;
        }
    }
}
