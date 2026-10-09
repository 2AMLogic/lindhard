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
        let current_n = match cfg.policy {
            StepPolicy::Fixed { ions_per_step } => {
                if ions_per_step == 0 {
                    return bad("ions_per_step must be at least 1");
                }
                ions_per_step
            }
            StepPolicy::Adaptive {
                max_ions_per_step,
                min_ions_per_step,
                max_change,
            } => {
                if min_ions_per_step == 0 || max_ions_per_step < min_ions_per_step {
                    return bad("need 1 <= min_ions_per_step <= max_ions_per_step");
                }
                if !(max_change.is_finite() && max_change > 0.0) {
                    return bad("max_change must be finite and positive");
                }
                max_ions_per_step
            }
        };
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
        })
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
            let change = self.max_change(&deltas);

            let shrink = match self.cfg.policy {
                StepPolicy::Adaptive {
                    min_ions_per_step,
                    max_change,
                    ..
                } if n > min_ions_per_step => {
                    if change > max_change {
                        let scaled = (n as f64 * max_change / change).floor() as u64;
                        Some(scaled.clamp(min_ions_per_step, n - 1))
                    } else {
                        None
                    }
                }
                _ => None,
            };
            if let Some(smaller) = shrink {
                n = smaller;
                continue;
            }

            let mut clamped = leftover;
            let outcome = match self.grid.apply(&deltas) {
                Ok(o) => o,
                Err(e) => {
                    let StepPolicy::Adaptive {
                        min_ions_per_step, ..
                    } = self.cfg.policy
                    else {
                        return Err(e.into());
                    };
                    if n > min_ions_per_step {
                        n = (n / 2).max(min_ions_per_step);
                        continue;
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
            if let StepPolicy::Adaptive {
                max_ions_per_step,
                max_change,
                ..
            } = self.cfg.policy
            {
                self.current_n = if change < 0.5 * max_change {
                    n.saturating_mul(2).min(max_ions_per_step)
                } else {
                    n
                };
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
