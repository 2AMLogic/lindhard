//! How often the electron transport read each cross-section table inside its
//! energy grid, and how often it used the constant continuation beyond an end.
//!
//! The transport interpolates the rate (and the inverse CDF) of each table
//! linearly in energy between the rows of the table's own grid, and holds the
//! first or last row constant outside it (`crate::electron::transport`, module
//! docs, "Free flight"). That continuation is a numerical choice, not physics:
//! a validation run whose electrons spend many flights beyond a grid end
//! measures the continuation as much as the model. The input resolver can
//! warn when a configured grid does not cover the expected energies; this
//! tally counts what actually happened, primaries and secondaries together.
//!
//! # What is counted
//!
//! One count per call of
//! [`ElectronTally::table_lookup`],
//! per layer and per channel (elastic, inelastic), split by where the energy
//! lay against that table's grid `[E_min, E_max]`:
//!
//! - `below`: `E < E_min`, the first row was used;
//! - `within`: `E_min <= E <= E_max`, endpoints included (an energy equal to
//!   an endpoint reads that row exactly, with no continuation);
//! - `above`: `E > E_max`, the last row was used.
//!
//! The transport evaluates both rates of the electron's layer once per pass
//! of its loop, before drawing a free path. So the counts are **rate
//! evaluations**, not collisions: a flight cut short at a layer face, a pass
//! after a face reflection and a flight with zero total rate are each one
//! evaluation of both tables and no collision, and a channel with zero rate
//! is evaluated (and counted) like any other. They are not fractions of
//! deposited energy or of path length either. The elastic and inelastic
//! counts of a layer are always equal in total; they differ in their split
//! when the two grids differ.
//!
//! # Determinism
//!
//! All counts are integers, merged by addition in chunk order, so they are
//! identical at any thread count for a fixed seed and chunk size. Counting
//! draws no random number and changes no state of the transport.

use serde::{Deserialize, Serialize};

use crate::electron::data::CrossSectionTable;
use crate::electron::transport::{
    ElectronState, ElectronTally, GridCoverage, TableChannel, Transport,
};

/// Rate evaluations of one table, by where the energy lay against its grid,
/// with the grid's bounds.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TableCoverageCounts {
    /// First energy of the table's grid, eV.
    pub energy_min_ev: f64,
    /// Last energy of the table's grid, eV.
    pub energy_max_ev: f64,
    /// Evaluations at `E < energy_min_ev` (first row continued).
    pub below: u64,
    /// Evaluations at `energy_min_ev <= E <= energy_max_ev`.
    pub within: u64,
    /// Evaluations at `E > energy_max_ev` (last row continued).
    pub above: u64,
}

impl TableCoverageCounts {
    /// No evaluation yet, with the bounds of `table`'s grid.
    pub fn new(table: &CrossSectionTable) -> Self {
        let g = table.energy_ev();
        Self {
            energy_min_ev: g[0],
            energy_max_ev: g[g.len() - 1],
            below: 0,
            within: 0,
            above: 0,
        }
    }

    /// All evaluations.
    pub fn total(&self) -> u64 {
        self.below + self.within + self.above
    }

    fn add(&mut self, c: GridCoverage) {
        match c {
            GridCoverage::Below => self.below += 1,
            GridCoverage::Within => self.within += 1,
            GridCoverage::Above => self.above += 1,
        }
    }

    fn merge(&mut self, o: &Self) {
        self.below += o.below;
        self.within += o.within;
        self.above += o.above;
    }
}

/// Rate evaluations of the two tables of one layer.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LayerTableCoverage {
    /// Layer index (0 = front).
    pub layer: usize,
    /// The elastic table.
    pub elastic: TableCoverageCounts,
    /// The inelastic table.
    pub inelastic: TableCoverageCounts,
}

/// Counts the rate evaluations of every table of a run against its grid (see
/// [the module docs](self)). Usable on its own as an [`ElectronTally`];
/// [`crate::tally::FullElectronTally`] carries one and puts its
/// [`TableCoverageTally::layers`] in
/// [`ElectronReport::table_coverage`](crate::tally::ElectronReport::table_coverage).
#[derive(Debug, Clone, PartialEq)]
pub struct TableCoverageTally {
    layers: Vec<LayerTableCoverage>,
}

impl TableCoverageTally {
    /// Empty counts for the tables of `transport`, one entry per layer.
    pub fn new(transport: &Transport) -> Self {
        Self {
            layers: transport
                .layer_tables()
                .iter()
                .enumerate()
                .map(|(layer, t)| LayerTableCoverage {
                    layer,
                    elastic: TableCoverageCounts::new(&t.elastic),
                    inelastic: TableCoverageCounts::new(&t.inelastic),
                })
                .collect(),
        }
    }

    /// The counts, by layer.
    pub fn layers(&self) -> &[LayerTableCoverage] {
        &self.layers
    }

    /// Count one evaluation.
    pub fn record(&mut self, layer: usize, channel: TableChannel, coverage: GridCoverage) {
        let l = &mut self.layers[layer];
        match channel {
            TableChannel::Elastic => l.elastic.add(coverage),
            TableChannel::Inelastic => l.inelastic.add(coverage),
        }
    }

    /// Add `o`'s counts (a later chunk) to these.
    ///
    /// # Panics
    /// If `o` was built for tables with a different number of layers or
    /// different grid bounds.
    pub fn merge_from(&mut self, o: &Self) {
        assert!(
            self.layers.len() == o.layers.len()
                && self.layers.iter().zip(&o.layers).all(|(a, b)| {
                    a.elastic.energy_min_ev == b.elastic.energy_min_ev
                        && a.elastic.energy_max_ev == b.elastic.energy_max_ev
                        && a.inelastic.energy_min_ev == b.inelastic.energy_min_ev
                        && a.inelastic.energy_max_ev == b.inelastic.energy_max_ev
                }),
            "merging table-coverage tallies built for different tables"
        );
        for (a, b) in self.layers.iter_mut().zip(&o.layers) {
            a.elastic.merge(&b.elastic);
            a.inelastic.merge(&b.inelastic);
        }
    }
}

impl ElectronTally for TableCoverageTally {
    fn table_lookup(&mut self, at: &ElectronState, channel: TableChannel, coverage: GridCoverage) {
        self.record(at.layer, channel, coverage);
    }

    fn merge(&mut self, other: Self) {
        self.merge_from(&other);
    }
}
