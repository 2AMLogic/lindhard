//! The tally-to-delta adapter: turns the transport events of one fluence step
//! into per-slab, per-element atom **counts**, which the driver scales to
//! areal inventory deltas.
//!
//! The conventions are those of the [parent module](super) ("Conservation and
//! update semantics"):
//!
//! * a recoil is subtracted from the slab it is created in, at the
//!   [`BcaTally::recoil`] event, and added to the slab where it comes to rest
//!   at the [`BcaTally::stopped`] event. A recoil that escapes therefore
//!   shows up as a loss with no further bookkeeping, and the origin slab of
//!   an escaped atom never has to be recovered from the (unlabelled) event
//!   sequence, which is what the aggregate exit totals could not give;
//! * a primary that comes to rest in a slab adds one atom of its species
//!   there; a primary that escapes adds nothing;
//! * the substrate is an immutable reservoir: an atom displaced from it
//!   removes nothing, and any atom that comes to rest in it adds nothing (the
//!   substrate is not part of the inventory).
//!
//! Layer indices of the stack built by
//! [`CompositionGrid::to_stack`](super::CompositionGrid::to_stack) are the slab
//! indices of the grid (the substrate, if any, is the one layer after the last
//! slab), so the adapter uses `Particle::layer` as the slab index directly.
//!
//! Counts are integers, so the merged result does not depend on the order in
//! which chunk tallies are folded, and the driver additionally merges them in
//! chunk order.

use std::collections::BTreeMap;

use crate::ion::bca::{BcaTally, Face, Particle};

/// Event counts of a block of primaries (one fluence step, or a sum of steps).
///
/// All counts are exact integers.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Yields {
    /// Primary histories run.
    pub histories: u64,
    /// Displaced target atoms (recoils created).
    pub displaced: u64,
    /// Primaries that came to rest in a finite slab.
    pub primaries_in_slabs: u64,
    /// Primaries that came to rest in the substrate.
    pub primaries_in_substrate: u64,
    /// Primaries that left through the front face.
    pub backscattered: u64,
    /// Primaries that left through the back face.
    pub transmitted: u64,
    /// Recoils that left through the front face (sputtered atoms), by atomic
    /// number.
    pub sputtered: BTreeMap<u8, u64>,
    /// Recoils that left through the back face.
    pub recoils_transmitted: u64,
}

impl Yields {
    /// Total sputtered atoms (all elements).
    pub fn sputtered_total(&self) -> u64 {
        self.sputtered.values().sum()
    }

    /// Add `other` field by field.
    pub fn add(&mut self, other: &Yields) {
        self.histories += other.histories;
        self.displaced += other.displaced;
        self.primaries_in_slabs += other.primaries_in_slabs;
        self.primaries_in_substrate += other.primaries_in_substrate;
        self.backscattered += other.backscattered;
        self.transmitted += other.transmitted;
        for (&z, &n) in &other.sputtered {
            *self.sputtered.entry(z).or_insert(0) += n;
        }
        self.recoils_transmitted += other.recoils_transmitted;
    }
}

/// A [`BcaTally`] that records the net atom count change per `(slab, Z)` and
/// the [`Yields`] of a step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryTally {
    n_slabs: usize,
    counts: BTreeMap<(usize, u8), i64>,
    yields: Yields,
}

impl InventoryTally {
    /// An empty tally for a grid of `n_slabs` finite slabs (layer indices at
    /// or above `n_slabs` are the substrate).
    pub fn new(n_slabs: usize) -> Self {
        Self {
            n_slabs,
            counts: BTreeMap::new(),
            yields: Yields::default(),
        }
    }

    /// Net atom count change per `(slab, Z)`, in `(slab, Z)` order. Zero
    /// entries are omitted.
    pub fn counts(&self) -> impl Iterator<Item = ((usize, u8), i64)> + '_ {
        self.counts
            .iter()
            .filter(|(_, &n)| n != 0)
            .map(|(&k, &n)| (k, n))
    }

    /// The event counts.
    pub fn yields(&self) -> &Yields {
        &self.yields
    }

    fn add(&mut self, layer: usize, z: u8, n: i64) {
        if layer < self.n_slabs {
            *self.counts.entry((layer, z)).or_insert(0) += n;
        }
    }
}

impl BcaTally for InventoryTally {
    fn begin_history(&mut self, _index: u64) {
        self.yields.histories += 1;
    }

    fn recoil(&mut self, recoil: &Particle) {
        self.yields.displaced += 1;
        self.add(recoil.layer, recoil.z, -1);
    }

    fn stopped(&mut self, p: &Particle) {
        self.add(p.layer, p.z, 1);
        if p.is_primary() {
            if p.layer < self.n_slabs {
                self.yields.primaries_in_slabs += 1;
            } else {
                self.yields.primaries_in_substrate += 1;
            }
        }
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        match (p.is_primary(), face) {
            (true, Face::Front) => self.yields.backscattered += 1,
            (true, Face::Back) => self.yields.transmitted += 1,
            (false, Face::Front) => *self.yields.sputtered.entry(p.z).or_insert(0) += 1,
            (false, Face::Back) => self.yields.recoils_transmitted += 1,
            // A stack has no lateral faces.
            (_, Face::Side) => {}
        }
    }

    fn merge(&mut self, other: Self) {
        for (k, n) in other.counts {
            *self.counts.entry(k).or_insert(0) += n;
        }
        self.yields.add(&other.yields);
    }
}
