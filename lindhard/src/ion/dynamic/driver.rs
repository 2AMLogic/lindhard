//! The fluence stepping loop.
//!
//! [`DynamicRun`] drives a [`CompositionGrid`] through a sequence of fluence
//! steps: run `n` primaries on the current target, turn the transport events
//! into inventory deltas ([`InventoryTally`]), apply them, and repeat until
//! all primaries of the beam are used.
//!
//! # Fluence
//!
//! The beam's `count` is the total number of primaries of the run and
//! [`DynamicConfig::fluence_m2`] the total fluence they represent, so each
//! primary stands for `fluence_m2 / count` atoms/m². A step of `n` primaries
//! scales its integer atom counts by that, giving the inventory deltas in
//! atoms/m². (The delta of a step is a mean over the `n` histories; the
//! target is held fixed during the step, so the step must be small enough
//! that the composition does not change much within it. That is what the
//! adaptive policy bounds.)
//!
//! # Determinism
//!
//! The primaries of a run are numbered `0..count` globally and primary `i`
//! always draws from [`stream`](crate::rng::stream)`(seed, i)`. A step that
//! covers global primaries `start..start + n` runs exactly those histories
//! ([`Bca::run_range`]), so a run is reproducible end to end, a step never
//! replays indices of an earlier one, and the thread count never matters:
//! per-chunk tallies are integer counts merged in chunk order and the deltas
//! are applied in `(slab, Z)` order.
//!
//! # Adaptive steps
//!
//! With [`StepPolicy::Adaptive`] a step is accepted only if its largest
//! relative composition change is at most `max_change`, where the change of a
//! slab is `|delta_Z| / (total atoms/m² of the slab)` maximised over the
//! elements `Z` (a slab emptied completely has change 1). A step that is too
//! large is **discarded** (the grid is not touched), shrunk to
//! `floor(n * max_change / change)` primaries (at least `min_ions_per_step`
//! and fewer than before) and retried from the **same** first index, so
//! rejected attempts consume no indices: the sequence of accepted
//! `(first_index, n)` pairs is decided from merged step results only, hence
//! deterministic. A step that fails [`CompositionGrid::apply`] (for example
//! it would remove more of an element than a slab holds) is treated the same
//! way, halving `n`. At `min_ions_per_step` a step is accepted whatever its
//! change; if it still removes more of an element than a slab holds (the last
//! atoms of an eroding slab), that removal is capped at what the slab holds
//! and counted in [`StepRecord::clamped`] (the step then creates the atoms the
//! cap refused to remove, so inventory is not conserved for that step by up to
//! one step's worth of removal). A [`StepPolicy::Fixed`] run does not cap: an
//! `apply` error is returned, with the grid unchanged, and the cure is a
//! smaller step or an adaptive policy. After an accepted step with
//! change below `max_change / 2` the next step size doubles, up to
//! `max_ions_per_step`. These rules are this crate's own design, not taken
//! from a published code.
//!
//! # Coordinates
//!
//! The front surface is `x = 0` of the grid at every step ([parent
//! module](super), "Coordinates"). Per-step outputs report the thickness of
//! the finite slabs and the interface depths measured from the surface of that
//! moment. With [`DynamicConfig::erosion`] off the surface never moves in the
//! sample frame either, and the depths are the depths of the original frame.
//! With erosion on, the surface recedes ("Erosion" below) and a depth `x`
//! at some step is `x + R` in the original frame, with `R` the cumulative
//! recession ([`StepRecord::recession_total_m`]) at that step.
//!
//! # Erosion
//!
//! With [`DynamicConfig::erosion`] on, sputtered atoms leave the target from
//! its front, not from the slab where they were displaced, as in the
//! dynamic-composition codes TRIDYN (Moller and Eckstein, Nucl. Instrum.
//! Methods B 2 (1984) 814, and the later TRIDYN papers: the target is updated
//! by removing the sputtered atoms from the surface layers and letting the
//! remaining material move up) and SDTrimSP (Mutzke et al., IPP report
//! 2019-02, the sections on the dynamic target and sputter erosion). Only the
//! published descriptions were used. The update, per accepted step:
//!
//! 1. The default tally subtracts every recoil at its origin slab, so an
//!    escaped recoil is already a loss there. Erosion adds those losses back
//!    (the tally counts them per origin slab, as integers, from
//!    `Particle::origin_layer`) and instead removes the same number of atoms
//!    of the same element, `Y_Z` per element, from the front of the target:
//!    slab 0 first, then deeper slabs when the element is used up (counting
//!    what the slab holds after the step's other deltas). Nothing is removed
//!    twice and no atom is created or lost: for each `Z` the net change of
//!    [`CompositionGrid::total_inventory`](super::CompositionGrid::total_inventory)
//!    is the implanted atoms minus the sputtered ones, as with erosion off.
//!    Atoms that came from the substrate were never in the inventory and
//!    remove nothing. Counts are scaled to atoms/m² like all other deltas.
//! 2. The step recession `dR` is the thickness the removed inventory occupied
//!    under the grid's [`Relaxation`](super::Relaxation):
//!    `sum_Z removed_Z v_Z` (ideal mixing) or `sum_Z removed_Z / n_mix` (fixed
//!    number density). It is tracked explicitly and summed over accepted steps
//!    into `R`; it cannot be recovered from the thickness, which also grows
//!    with implantation.
//! 3. Removal is part of the deltas the step is checked and applied with, so
//!    `max_change` and the adaptive policy see it. Removal never exceeds what
//!    the slabs hold, so it is never the cause of the `clamped` cap.
//! 4. The grid is rebuilt from `x = 0` as usual, which is the re-anchoring.
//!
//! Erosion is a pure function of the step's merged integer tally and draws no
//! random numbers, so it keeps the determinism of the run. With erosion off
//! none of this code runs and results are unchanged.
//!
//! # Retention limit
//!
//! [`DynamicRun::with_retention`] sets, per implanted species, a saturation
//! limit: the largest atomic fraction `f` of that species a slab may hold. The
//! dynamic-composition codes TRIDYN (Moller and Eckstein, Nucl. Instrum.
//! Methods B 2 (1984) 814) and SDTrimSP (Mutzke et al., IPP report 2019-02)
//! describe a saturation concentration of the implanted species above which
//! further atoms are not retained and are released again; only those published
//! descriptions were used. The update below is this crate's own statement of
//! that idea, and no numeric default is supplied: the limit is always given by
//! the caller.
//!
//! Per accepted step, after the transport deltas (and erosion) are formed and
//! before they are checked or applied, each slab is looked at with the
//! composition it would have after the step. For a limited species `Z` with
//! amount `a_Z` and `o` atoms/m² of all other elements, the slab may hold at
//! most `a_max = f o / (1 - f)` of `Z` (that is `a_Z / (a_Z + o) <= f`). Any
//! excess `a_Z - a_max` is **re-emitted**: it is added to the step's deltas as
//! a removal from that slab and leaves the target, so it is not inventory any
//! more and is reported in [`StepRecord::re_emitted_m2`]. Species are taken in
//! ascending `Z` and each sees the amounts left by the earlier ones, so with
//! several limited species in one slab the fractions are met in that order
//! (removing one raises the fraction of the others slightly); with one limited
//! species, the usual case, the limit is met exactly. A slab holding only the
//! limited species has `o = 0` and loses all of it. Because the removal is part
//! of the deltas, `max_change` and the adaptive policy see it. It is a pure
//! function of the merged step tally and the grid, draws no random numbers and
//! so keeps the thread-count invariance. The re-emitted atoms are not placed
//! anywhere: this is a bookkeeping limit, not a model of diffusion or of where
//! the released gas goes. With no limits set none of this code runs and
//! results are unchanged.
//!
//! The atom balance of an element `Z` over a step is
//! `Δ total_inventory(Z) = transport net - re_emitted`, where the transport net
//! counts retained primaries and moved, sputtered and transmitted recoils per
//! the conventions of the [parent module](super).

use std::collections::BTreeMap;

use crate::ion::bca::{Bca, BcaConfig, BcaError, Beam};
use crate::ion::scattering::ScatteringTable;
use crate::ion::stopping::ElectronicStopping;

use super::adapter::{InventoryTally, Yields};
use super::{CompositionGrid, DynamicError, InventoryDelta};

/// How the primaries are divided into fluence steps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StepPolicy {
    /// Every step has `ions_per_step` primaries (the last may have fewer).
    Fixed {
        /// Primaries per step, at least 1.
        ions_per_step: u64,
    },
    /// Steps start at `max_ions_per_step`, shrink when the composition change
    /// exceeds `max_change` and grow back (see the [module docs](super)).
    Adaptive {
        /// Largest and initial step, at least `min_ions_per_step`.
        max_ions_per_step: u64,
        /// Smallest step, at least 1.
        min_ions_per_step: u64,
        /// Largest accepted relative composition change per step, `> 0`.
        max_change: f64,
    },
}

/// Fluence-loop settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DynamicConfig {
    /// Total fluence of the whole run, ions/m², represented by the beam's
    /// `count` primaries.
    pub fluence_m2: f64,
    /// Step sizes.
    pub policy: StepPolicy,
    /// Sputter erosion and surface recession (see the [`DynamicRun`] docs,
    /// "Erosion"). Off reproduces the fixed-front behaviour exactly.
    pub erosion: bool,
}

/// Errors from the fluence loop.
#[derive(Debug, thiserror::Error)]
pub enum DynamicRunError {
    /// A setting is out of range.
    #[error("invalid dynamic run setting: {0}")]
    Invalid(String),
    /// The grid rejected an update (or could not be turned into a stack) and
    /// no smaller step is allowed.
    #[error("composition update failed: {0}")]
    Grid(#[from] DynamicError),
    /// The transport engine could not be set up or failed.
    #[error("transport failed: {0}")]
    Transport(#[from] BcaError),
}

/// The step size a run starts with, after validating `policy`. Both the 1-D
/// and the 2-D run call this in their constructors, so the two runs reject
/// the same policies and start at the same size.
pub(super) fn initial_n(policy: StepPolicy) -> Result<u64, DynamicRunError> {
    match policy {
        StepPolicy::Fixed { ions_per_step } => {
            if ions_per_step == 0 {
                return Err(DynamicRunError::Invalid(
                    "ions_per_step must be at least 1".to_string(),
                ));
            }
            Ok(ions_per_step)
        }
        StepPolicy::Adaptive {
            max_ions_per_step,
            min_ions_per_step,
            max_change,
        } => {
            if min_ions_per_step == 0 || max_ions_per_step < min_ions_per_step {
                return Err(DynamicRunError::Invalid(
                    "need 1 <= min_ions_per_step <= max_ions_per_step".to_string(),
                ));
            }
            if !(max_change.is_finite() && max_change > 0.0) {
                return Err(DynamicRunError::Invalid(
                    "max_change must be finite and positive".to_string(),
                ));
            }
            Ok(max_ions_per_step)
        }
    }
}

/// The step size to retry a step whose composition change exceeded
/// `max_change`: `Some(smaller)` when an adaptive run can still shrink `n`,
/// `None` otherwise (a [`StepPolicy::Fixed`] policy, a run already at
/// `min_ions_per_step`, or a change within budget). Shared by both runs.
pub(super) fn shrink(policy: StepPolicy, n: u64, change: f64) -> Option<u64> {
    let StepPolicy::Adaptive {
        min_ions_per_step,
        max_change,
        ..
    } = policy
    else {
        return None;
    };
    if n <= min_ions_per_step {
        return None;
    }
    // The reject condition, exactly as both runs wrote it: `change >
    // max_change`, which is false for a NaN change (then nothing shrinks).
    let over_budget = change > max_change;
    if !over_budget {
        return None;
    }
    let scaled = (n as f64 * max_change / change).floor() as u64;
    // The upper clamp is reachable: with n = 3 and a `change` one ulp above
    // `max_change`, double rounding puts the scaled value at `n` itself.
    Some(scaled.clamp(min_ions_per_step, n - 1))
}

/// The step size to retry a failed grid update at: half of `n`, floored at
/// `min_ions_per_step`; `None` when no smaller step exists (a
/// [`StepPolicy::Fixed`] policy, or an adaptive run already at the floor).
/// Shared by both runs.
pub(super) fn halve(policy: StepPolicy, n: u64) -> Option<u64> {
    let StepPolicy::Adaptive {
        min_ions_per_step, ..
    } = policy
    else {
        return None;
    };
    (n > min_ions_per_step).then(|| (n / 2).max(min_ions_per_step))
}

/// The next step size after an accepted step: `n` doubled when the change
/// stayed below half the budget, capped at `max_ions_per_step`; `None` for a
/// [`StepPolicy::Fixed`] policy, whose step size never changes (the caller
/// leaves its `current_n` alone). Shared by both runs.
pub(super) fn grow(policy: StepPolicy, n: u64, change: f64) -> Option<u64> {
    let StepPolicy::Adaptive {
        max_ions_per_step,
        max_change,
        ..
    } = policy
    else {
        return None;
    };
    Some(if change < 0.5 * max_change {
        n.saturating_mul(2).min(max_ions_per_step)
    } else {
        n
    })
}

/// What one accepted step did.
#[derive(Debug, Clone, PartialEq)]
pub struct StepRecord {
    /// Step number, from 1.
    pub step: u64,
    /// Global index of the first primary of the step.
    pub first_index: u64,
    /// Primaries in the step.
    pub ions: u64,
    /// Attempts made at this `first_index` (1 if never rejected).
    pub attempts: u32,
    /// Largest relative composition change of the step (see the
    /// [module docs](super)).
    pub max_change: f64,
    /// `(slab, element)` removals that were capped at what the slab held (see
    /// the [module docs](super)); 0 for a step that needed no capping.
    pub clamped: u32,
    /// Event counts of the step.
    pub yields: Yields,
    /// Slabs (indices before the step) that emptied and were removed.
    pub removed_slabs: Vec<usize>,
    /// Surface recession of this step, m (0 with erosion off).
    pub recession_m: f64,
    /// Cumulative recession up to and including this step, m (0 with erosion
    /// off): add it to a depth of this step to get the depth in the original
    /// frame.
    pub recession_total_m: f64,
    /// Atoms/m² re-emitted by the retention limit in this step, by atomic
    /// number (see the module docs, "Retention limit"); empty with no
    /// limit set or nothing in excess.
    pub re_emitted_m2: BTreeMap<u8, f64>,
}

/// A fluence-stepped run. See the [module docs](super).
pub struct DynamicRun<'a> {
    grid: CompositionGrid,
    beam: Beam,
    config: BcaConfig,
    stopping: &'a (dyn ElectronicStopping + Sync),
    table: &'a ScatteringTable,
    cfg: DynamicConfig,
    next_index: u64,
    steps: u64,
    current_n: u64,
    cumulative: Yields,
    recession_m: f64,
    retention: BTreeMap<u8, f64>,
    re_emitted_m2: BTreeMap<u8, f64>,
}

// Manual: the borrowed stopping model is a trait object without `Debug`, and
// the scattering table is large.
impl std::fmt::Debug for DynamicRun<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DynamicRun")
            .field("beam", &self.beam)
            .field("config", &self.config)
            .field("cfg", &self.cfg)
            .field("next_index", &self.next_index)
            .field("steps", &self.steps)
            .field("current_n", &self.current_n)
            .field("cumulative", &self.cumulative)
            .field("recession_m", &self.recession_m)
            .field("retention", &self.retention)
            .finish_non_exhaustive()
    }
}

impl<'a> DynamicRun<'a> {
    /// Set up a run of `beam.count` primaries on `grid`.
    pub fn new(
        grid: CompositionGrid,
        beam: Beam,
        config: BcaConfig,
        stopping: &'a (dyn ElectronicStopping + Sync),
        table: &'a ScatteringTable,
        cfg: DynamicConfig,
    ) -> Result<Self, DynamicRunError> {
        let bad = |m: &str| Err(DynamicRunError::Invalid(m.to_string()));
        if beam.count == 0 {
            return bad("the beam needs at least one primary");
        }
        if !(cfg.fluence_m2.is_finite() && cfg.fluence_m2 > 0.0) {
            return bad("fluence must be finite and positive");
        }
        let current_n = initial_n(cfg.policy)?;
        Ok(Self {
            grid,
            beam,
            config,
            stopping,
            table,
            cfg,
            next_index: 0,
            steps: 0,
            current_n,
            cumulative: Yields::default(),
            recession_m: 0.0,
            retention: BTreeMap::new(),
            re_emitted_m2: BTreeMap::new(),
        })
    }

    /// Set the retention limits (see the module docs, "Retention
    /// limit"): for each `(Z, f)` the largest atomic fraction `f`, with
    /// `0 < f < 1`, of element `Z` a slab may hold. Replaces any earlier
    /// limits; an empty map switches the model off (the default).
    pub fn with_retention(
        mut self,
        limits: impl IntoIterator<Item = (u8, f64)>,
    ) -> Result<Self, DynamicRunError> {
        let mut map = BTreeMap::new();
        for (z, f) in limits {
            if !(f.is_finite() && f > 0.0 && f < 1.0) {
                return Err(DynamicRunError::Invalid(format!(
                    "the retention limit of Z = {z} must be an atomic fraction in (0, 1)"
                )));
            }
            if map.insert(z, f).is_some() {
                return Err(DynamicRunError::Invalid(format!(
                    "the retention limit of Z = {z} is given twice"
                )));
            }
        }
        self.retention = map;
        Ok(self)
    }

    /// Atoms/m² re-emitted by the retention limit over the accepted steps, by
    /// atomic number.
    pub fn re_emitted_m2(&self) -> &BTreeMap<u8, f64> {
        &self.re_emitted_m2
    }

    /// The current target.
    pub fn grid(&self) -> &CompositionGrid {
        &self.grid
    }

    /// Primaries used so far (accepted steps only).
    pub fn ions_done(&self) -> u64 {
        self.next_index
    }

    /// Total primaries of the run.
    pub fn ions_total(&self) -> u64 {
        self.beam.count
    }

    /// Fluence delivered so far, ions/m².
    pub fn fluence_done_m2(&self) -> f64 {
        self.next_index as f64 * self.fluence_per_ion_m2()
    }

    /// Fluence one primary stands for, ions/m².
    pub fn fluence_per_ion_m2(&self) -> f64 {
        self.cfg.fluence_m2 / self.beam.count as f64
    }

    /// Accepted steps so far.
    pub fn steps_done(&self) -> u64 {
        self.steps
    }

    /// Event counts summed over the accepted steps.
    pub fn cumulative(&self) -> &Yields {
        &self.cumulative
    }

    /// Cumulative surface recession of the accepted steps, m (0 with erosion
    /// off).
    pub fn recession_m(&self) -> f64 {
        self.recession_m
    }

    /// Run the next step. `Ok(None)` when all primaries are used. On error the
    /// grid and the run state are unchanged.
    pub fn step(&mut self) -> Result<Option<StepRecord>, DynamicRunError> {
        let remaining = self.beam.count - self.next_index;
        if remaining == 0 {
            return Ok(None);
        }
        let first = self.next_index;
        let mut n = self.current_n.min(remaining);
        let mut attempts = 0u32;
        loop {
            attempts += 1;
            let n_slabs = self.grid.n_slabs();
            let stack = self.grid.to_stack()?;
            let beam = Beam {
                count: n,
                ..self.beam
            };
            let bca = Bca::new(beam, &stack, self.config, self.stopping, self.table)?;
            let tally = bca.run_range(first, n, || InventoryTally::new(n_slabs))?;

            let per_ion = self.fluence_per_ion_m2();
            let (deltas, step_recession, leftover) = if self.cfg.erosion {
                self.erosion_deltas(&tally, per_ion)?
            } else {
                let d = tally
                    .counts()
                    .map(|((slab, z), c)| InventoryDelta {
                        slab,
                        z,
                        delta_atoms_m2: c as f64 * per_ion,
                    })
                    .collect();
                (d, 0.0, 0)
            };
            let (deltas, re_emitted) = self.retention_deltas(deltas);
            let change = self.max_change(&deltas);

            if let Some(smaller) = shrink(self.cfg.policy, n, change) {
                n = smaller;
                continue;
            }

            let mut clamped = leftover;
            let outcome = match self.grid.apply(&deltas) {
                Ok(o) => o,
                Err(e) => {
                    if let Some(smaller) = halve(self.cfg.policy, n) {
                        n = smaller;
                        continue;
                    }
                    // No smaller step: a Fixed run returns the error; the
                    // smallest adaptive step caps the removal below.
                    if !matches!(self.cfg.policy, StepPolicy::Adaptive { .. }) {
                        return Err(e.into());
                    }
                    // Smallest step: a removal beyond what a slab holds is
                    // capped at what it holds.
                    if !matches!(e, DynamicError::NegativeInventory { .. }) {
                        return Err(e.into());
                    }
                    let capped: Vec<InventoryDelta> = deltas
                        .iter()
                        .map(|d| {
                            let have = self
                                .grid
                                .inventory(d.slab)
                                .and_then(|v| v.iter().find(|&&(z, _)| z == d.z).map(|&(_, a)| a))
                                .unwrap_or(0.0);
                            if d.delta_atoms_m2 < -have {
                                clamped += 1;
                                InventoryDelta {
                                    delta_atoms_m2: -have,
                                    ..*d
                                }
                            } else {
                                *d
                            }
                        })
                        .collect();
                    self.grid.apply(&capped)?
                }
            };
            self.steps += 1;
            self.next_index += n;
            self.cumulative.add(tally.yields());
            self.recession_m += step_recession;
            for (&z, &a) in &re_emitted {
                *self.re_emitted_m2.entry(z).or_insert(0.0) += a;
            }
            if let Some(next) = grow(self.cfg.policy, n, change) {
                self.current_n = next;
            }
            return Ok(Some(StepRecord {
                step: self.steps,
                first_index: first,
                ions: n,
                attempts,
                max_change: change,
                clamped,
                yields: tally.yields().clone(),
                removed_slabs: outcome.removed_slabs,
                recession_m: step_recession,
                recession_total_m: self.recession_m,
                re_emitted_m2: re_emitted,
            }));
        }
    }

    /// The deltas of a step with erosion on, its recession in m, and the
    /// number of elements whose sputtered count exceeded what the slabs hold
    /// (that excess is not removed; the step reports it in `clamped`). See the
    /// `DynamicRun` docs, "Erosion".
    fn erosion_deltas(
        &self,
        tally: &InventoryTally,
        per_ion: f64,
    ) -> Result<(Vec<InventoryDelta>, f64, u32), DynamicRunError> {
        // Integer net counts with the origin-slab loss of sputtered atoms
        // cancelled; sputtered atoms per element to take from the front.
        let mut net: BTreeMap<(usize, u8), i64> = tally.counts().collect();
        let mut sputtered: BTreeMap<u8, u64> = BTreeMap::new();
        for ((slab, z), n) in tally.sputtered_from() {
            *net.entry((slab, z)).or_insert(0) += n as i64;
            *sputtered.entry(z).or_insert(0) += n;
        }
        let mut delta: BTreeMap<(usize, u8), f64> = net
            .into_iter()
            .map(|(k, c)| (k, c as f64 * per_ion))
            .collect();
        let inventories: Vec<Vec<(u8, f64)>> = (0..self.grid.n_slabs())
            .map(|i| self.grid.inventory(i).unwrap_or_default())
            .collect();
        let mut removed: BTreeMap<u8, f64> = BTreeMap::new();
        let mut leftover = 0u32;
        for (&z, &y) in &sputtered {
            let mut rem = y as f64 * per_ion;
            for (slab, inv) in inventories.iter().enumerate() {
                if rem <= 0.0 {
                    break;
                }
                let have = inv
                    .iter()
                    .find(|&&(zz, _)| zz == z)
                    .map_or(0.0, |&(_, a)| a);
                let d = delta.entry((slab, z)).or_insert(0.0);
                let avail = (have + *d).max(0.0);
                let take = avail.min(rem);
                if take > 0.0 {
                    *d -= take;
                    rem -= take;
                    *removed.entry(z).or_insert(0.0) += take;
                }
            }
            if rem > 0.0 {
                leftover += 1;
            }
        }
        let recession = if removed.is_empty() {
            0.0
        } else {
            self.grid.relaxation().thickness_m(&removed)?
        };
        let deltas = delta
            .into_iter()
            .filter(|&(_, d)| d != 0.0)
            .map(|((slab, z), d)| InventoryDelta {
                slab,
                z,
                delta_atoms_m2: d,
            })
            .collect();
        Ok((deltas, recession, leftover))
    }

    /// Add the removals the retention limits require to the step's `deltas`;
    /// returns them with the atoms/m² re-emitted per element. See the module
    /// docs, "Retention limit". Returns the deltas unchanged with no limits.
    fn retention_deltas(
        &self,
        deltas: Vec<InventoryDelta>,
    ) -> (Vec<InventoryDelta>, BTreeMap<u8, f64>) {
        let mut re_emitted = BTreeMap::new();
        if self.retention.is_empty() {
            return (deltas, re_emitted);
        }
        let mut delta: BTreeMap<(usize, u8), f64> = BTreeMap::new();
        for d in &deltas {
            *delta.entry((d.slab, d.z)).or_insert(0.0) += d.delta_atoms_m2;
        }
        for slab in 0..self.grid.n_slabs() {
            let inv = self.grid.inventory(slab).unwrap_or_default();
            let have = |z: u8| {
                inv.iter()
                    .find(|&&(zz, _)| zz == z)
                    .map_or(0.0, |&(_, a)| a)
            };
            // Composition after the step, in Z order.
            let mut after: BTreeMap<u8, f64> = inv.iter().map(|&(z, a)| (z, a)).collect();
            for (&(s, z), &d) in delta.range((slab, 0)..=(slab, u8::MAX)) {
                debug_assert_eq!(s, slab);
                *after.entry(z).or_insert(0.0) += d;
            }
            for (&z, &f) in &self.retention {
                let a = after.get(&z).copied().unwrap_or(0.0);
                let others: f64 = after
                    .iter()
                    .filter(|(&zz, _)| zz != z)
                    .map(|(_, &x)| x)
                    .sum();
                let a_max = f * others / (1.0 - f);
                if a <= a_max {
                    continue;
                }
                let excess = a - a_max;
                let d = delta.entry((slab, z)).or_insert(0.0);
                // A slab with nothing else loses the species exactly.
                *d = if a_max == 0.0 { -have(z) } else { *d - excess };
                after.insert(z, a_max);
                *re_emitted.entry(z).or_insert(0.0) += excess;
            }
        }
        let deltas = delta
            .into_iter()
            .filter(|&(_, d)| d != 0.0)
            .map(|((slab, z), d)| InventoryDelta {
                slab,
                z,
                delta_atoms_m2: d,
            })
            .collect();
        (deltas, re_emitted)
    }

    /// Largest `|delta_Z| / (atoms/m^2 of the slab)` over the deltas.
    fn max_change(&self, deltas: &[InventoryDelta]) -> f64 {
        let totals: Vec<f64> = (0..self.grid.n_slabs())
            .map(|i| {
                self.grid
                    .inventory(i)
                    .map_or(0.0, |v| v.iter().map(|&(_, a)| a).sum())
            })
            .collect();
        deltas
            .iter()
            .map(|d| d.delta_atoms_m2.abs() / totals[d.slab])
            .fold(0.0, f64::max)
    }

    /// Run all remaining steps, calling `on_step` after each accepted step.
    pub fn run_all(
        &mut self,
        mut on_step: impl FnMut(&StepRecord, &CompositionGrid),
    ) -> Result<(), DynamicRunError> {
        while let Some(rec) = self.step()? {
            on_step(&rec, &self.grid);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn adaptive(max_ions_per_step: u64, min_ions_per_step: u64, max_change: f64) -> StepPolicy {
        StepPolicy::Adaptive {
            max_ions_per_step,
            min_ions_per_step,
            max_change,
        }
    }

    fn invalid(policy: StepPolicy) -> String {
        let DynamicRunError::Invalid(m) = initial_n(policy).unwrap_err() else {
            panic!("expected Invalid");
        };
        m
    }

    #[test]
    fn initial_n_validates_and_starts_the_run() {
        assert_eq!(
            initial_n(StepPolicy::Fixed { ions_per_step: 7 }).unwrap(),
            7
        );
        assert_eq!(initial_n(adaptive(100, 10, 0.1)).unwrap(), 100);
    }

    #[test]
    fn initial_n_rejects_a_zero_fixed_step() {
        assert_eq!(
            invalid(StepPolicy::Fixed { ions_per_step: 0 }),
            "ions_per_step must be at least 1"
        );
    }

    #[test]
    fn initial_n_rejects_an_empty_adaptive_range() {
        assert_eq!(
            invalid(adaptive(10, 0, 0.1)),
            "need 1 <= min_ions_per_step <= max_ions_per_step"
        );
        assert_eq!(
            invalid(adaptive(10, 11, 0.1)),
            "need 1 <= min_ions_per_step <= max_ions_per_step"
        );
    }

    #[test]
    fn initial_n_rejects_a_non_finite_or_non_positive_max_change() {
        assert_eq!(
            invalid(adaptive(100, 1, f64::NAN)),
            "max_change must be finite and positive"
        );
        assert_eq!(
            invalid(adaptive(100, 1, f64::INFINITY)),
            "max_change must be finite and positive"
        );
        assert_eq!(
            invalid(adaptive(100, 1, 0.0)),
            "max_change must be finite and positive"
        );
        assert_eq!(
            invalid(adaptive(100, 1, -0.5)),
            "max_change must be finite and positive"
        );
    }

    #[test]
    fn shrink_leaves_a_within_budget_change_alone() {
        let p = adaptive(1000, 10, 0.1);
        assert_eq!(shrink(p, 1000, 0.1), None);
        assert_eq!(shrink(p, 1000, 0.05), None);
        // At the floor there is nothing smaller to retry.
        assert_eq!(shrink(p, 10, 100.0), None);
        // A non-finite change is never above the budget (`NaN > x` is
        // false), so nothing shrinks.
        assert_eq!(shrink(p, 1000, f64::NAN), None);
    }

    #[test]
    fn shrink_scales_and_clamps_at_min_ions_per_step() {
        let p = adaptive(1000, 5, 0.1);
        // 1000 * 0.1 / 1000.0 floors to 0, clamped up to the floor of 5.
        assert_eq!(shrink(p, 1000, 1000.0), Some(5));
    }

    #[test]
    fn shrink_clamps_at_n_minus_one() {
        let p = adaptive(100, 1, 0.1);
        // One ulp above max_change: double rounding floors 3 * 0.1 / change
        // to 3 itself, so without the clamp the "smaller" step would not
        // shrink at all.
        let one_ulp_up = f64::from_bits(0.1_f64.to_bits() + 1);
        assert_eq!(shrink(p, 3, one_ulp_up), Some(2));
    }

    #[test]
    fn shrink_and_halve_do_nothing_for_a_fixed_policy() {
        let p = StepPolicy::Fixed { ions_per_step: 4 };
        assert_eq!(shrink(p, 4, 100.0), None);
        assert_eq!(halve(p, 4), None);
    }

    #[test]
    fn halve_floors_at_min_ions_per_step() {
        let p = adaptive(1000, 3, 0.1);
        assert_eq!(halve(p, 7), Some(3));
        assert_eq!(halve(p, 5), Some(3));
        // Halving above the floor stays above it.
        let q = adaptive(1000, 2, 0.1);
        assert_eq!(halve(q, 6), Some(3));
        // Already at the floor: the caller falls back to the cap path.
        assert_eq!(halve(p, 3), None);
    }

    #[test]
    fn grow_doubles_and_caps_at_max_ions_per_step() {
        let p = adaptive(100, 10, 0.1);
        assert_eq!(grow(p, 40, 0.01), Some(80));
        assert_eq!(grow(p, 60, 0.01), Some(100));
        assert_eq!(grow(p, 100, 0.01), Some(100));
    }

    #[test]
    fn grow_keeps_the_size_when_the_change_used_half_the_budget() {
        let p = adaptive(100, 10, 0.1);
        assert_eq!(grow(p, 40, 0.05), Some(40));
        assert_eq!(grow(p, 40, 0.5), Some(40));
    }

    #[test]
    fn grow_leaves_a_fixed_policy_alone() {
        let p = StepPolicy::Fixed { ions_per_step: 4 };
        assert_eq!(grow(p, 2, 0.0), None);
    }
}
