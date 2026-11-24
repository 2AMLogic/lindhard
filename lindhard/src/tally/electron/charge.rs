//! Deposited-charge bookkeeping of the electron tally: where electrons
//! stopped, where holes were left, and the escapes, per layer and per cell of
//! the optional deposition grids, with the per-history data for the standard
//! error of the net charge. [`ChargeAcc`] owns the bookkeeping; the
//! [`FullElectronTally`](super::FullElectronTally) delegates to it.

use serde::{Deserialize, Serialize};

use super::{CartesianGrid, CylindricalGrid};
use crate::electron::transport::{ElectronState, Face};

/// Electron counts at one place, split by the generation of the electron:
/// charge in units of the electron's, per [`ChargeReport`]'s convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct ChargeCounts {
    /// Primaries that terminated here (every terminal fate but escape).
    pub stopped_primary: u64,
    /// Secondaries that terminated here, including a liberated electron
    /// dropped at its creation point because it was below the stopping
    /// threshold.
    pub stopped_secondary: u64,
    /// Holes (liberated electrons) created here by collisions of a primary.
    pub created_by_primary: u64,
    /// Holes created here by collisions of a secondary.
    pub created_by_secondary: u64,
}

impl ChargeCounts {
    /// Electrons that terminated here.
    pub fn stopped(&self) -> u64 {
        self.stopped_primary + self.stopped_secondary
    }

    /// Electrons liberated here (holes left behind).
    pub fn created(&self) -> u64 {
        self.created_by_primary + self.created_by_secondary
    }

    /// Net electron number gained here: stopped minus created.
    pub fn net(&self) -> i64 {
        self.stopped() as i64 - self.created() as i64
    }

    fn merge(&mut self, o: &ChargeCounts) {
        self.stopped_primary += o.stopped_primary;
        self.stopped_secondary += o.stopped_secondary;
        self.created_by_primary += o.created_by_primary;
        self.created_by_secondary += o.created_by_secondary;
    }
}

/// Charge counts per cell of one grid, with the per-history net charge
/// squared and summed for the standard error. Everything is integer, so the
/// merge is exact and independent of the order of the parts.
#[derive(Debug, Clone, PartialEq)]
struct ChargeGrid {
    counts: Vec<ChargeCounts>,
    /// Sum over histories of the squared net charge the history left in the
    /// cell.
    net_sum_sq: Vec<u64>,
    outside: ChargeCounts,
    outside_net_sum_sq: u64,
    /// Per-history scratch: `(cell, net change)`; `usize::MAX` is outside.
    touched: Vec<(usize, i64)>,
}

impl ChargeGrid {
    fn new(cells: usize) -> Self {
        Self {
            counts: vec![ChargeCounts::default(); cells],
            net_sum_sq: vec![0; cells],
            outside: ChargeCounts::default(),
            outside_net_sum_sq: 0,
            touched: Vec::new(),
        }
    }

    /// An electron terminated in `cell`.
    fn stop(&mut self, cell: Option<usize>, secondary: bool) {
        let c = self.slot(cell, 1);
        if secondary {
            c.stopped_secondary += 1;
        } else {
            c.stopped_primary += 1;
        }
    }

    /// A hole was created in `cell` by an electron of the given generation.
    fn create(&mut self, cell: Option<usize>, by_secondary: bool) {
        let c = self.slot(cell, -1);
        if by_secondary {
            c.created_by_secondary += 1;
        } else {
            c.created_by_primary += 1;
        }
    }

    fn slot(&mut self, cell: Option<usize>, delta: i64) -> &mut ChargeCounts {
        self.touched.push((cell.unwrap_or(usize::MAX), delta));
        match cell {
            Some(i) => &mut self.counts[i],
            None => &mut self.outside,
        }
    }

    /// Fold the history's per-cell net charge into the squared sums.
    fn end_history(&mut self) {
        self.touched.sort_unstable_by_key(|t| t.0);
        let mut i = 0;
        while i < self.touched.len() {
            let cell = self.touched[i].0;
            let mut net = 0i64;
            while i < self.touched.len() && self.touched[i].0 == cell {
                net += self.touched[i].1;
                i += 1;
            }
            let sq = (net * net) as u64;
            if cell == usize::MAX {
                self.outside_net_sum_sq += sq;
            } else {
                self.net_sum_sq[cell] += sq;
            }
        }
        self.touched.clear();
    }

    fn merge(&mut self, o: &ChargeGrid) {
        for (a, b) in self.counts.iter_mut().zip(&o.counts) {
            a.merge(b);
        }
        for (a, b) in self.net_sum_sq.iter_mut().zip(&o.net_sum_sq) {
            *a += b;
        }
        self.outside.merge(&o.outside);
        self.outside_net_sum_sq += o.outside_net_sum_sq;
    }

    fn report(&self) -> ChargeCells {
        ChargeCells {
            counts: self.counts.clone(),
            net_sum_sq: self.net_sum_sq.clone(),
            outside: self.outside,
            outside_net_sum_sq: self.outside_net_sum_sq,
        }
    }
}

/// Charge bookkeeping of a tally: per layer, per grid cell, and escapes.
/// Each optional grid is carried with its counts so the bookkeeping needs
/// nothing but the transport events.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct ChargeAcc {
    layers: ChargeGrid,
    cartesian: Option<(CartesianGrid, ChargeGrid)>,
    cylindrical: Option<(CylindricalGrid, ChargeGrid)>,
    /// Escapes through each face: `[primary, secondary]`.
    escaped_front: [u64; 2],
    escaped_back: [u64; 2],
    escaped_side: [u64; 2],
}

impl ChargeAcc {
    /// Charge grids for a run over `n_layers` layers with the optional
    /// deposition grids of the tally's configuration.
    pub(super) fn new(
        n_layers: usize,
        cartesian: Option<CartesianGrid>,
        cylindrical: Option<CylindricalGrid>,
    ) -> Self {
        Self {
            layers: ChargeGrid::new(n_layers),
            // The grids were size-validated before this (a `GridTooLarge`
            // configuration never reaches a tally), so the fallback count is
            // unreachable.
            cartesian: cartesian.map(|g| (g, ChargeGrid::new(g.voxels().unwrap_or(0)))),
            cylindrical: cylindrical.map(|g| (g, ChargeGrid::new(g.cells().unwrap_or(0)))),
            escaped_front: [0; 2],
            escaped_back: [0; 2],
            escaped_side: [0; 2],
        }
    }

    /// An electron of the given generation terminated at `at`.
    pub(super) fn stop(&mut self, at: &ElectronState, secondary: bool) {
        self.layers.stop(Some(at.layer), secondary);
        if let Some((g, a)) = &mut self.cartesian {
            a.stop(g.locate(at.pos), secondary);
        }
        if let Some((g, a)) = &mut self.cylindrical {
            a.stop(g.locate(at.pos), secondary);
        }
    }

    /// A hole was left at `at` by a collision of an electron of the given
    /// generation.
    pub(super) fn create(&mut self, at: &ElectronState, by_secondary: bool) {
        self.layers.create(Some(at.layer), by_secondary);
        if let Some((g, a)) = &mut self.cartesian {
            a.create(g.locate(at.pos), by_secondary);
        }
        if let Some((g, a)) = &mut self.cylindrical {
            a.create(g.locate(at.pos), by_secondary);
        }
    }

    /// An electron of the given generation left through `face`.
    pub(super) fn escaped(&mut self, face: Face, secondary: bool) {
        let e = match face {
            Face::Front => &mut self.escaped_front,
            Face::Back => &mut self.escaped_back,
            Face::Side => &mut self.escaped_side,
        };
        e[usize::from(secondary)] += 1;
    }

    /// Fold the history's per-cell net charge into the squared sums.
    pub(super) fn end_history(&mut self) {
        self.layers.end_history();
        if let Some((_, a)) = &mut self.cartesian {
            a.end_history();
        }
        if let Some((_, a)) = &mut self.cylindrical {
            a.end_history();
        }
    }

    pub(super) fn merge(&mut self, o: &ChargeAcc) {
        self.layers.merge(&o.layers);
        if let (Some(a), Some(b)) = (&mut self.cartesian, &o.cartesian) {
            a.1.merge(&b.1);
        }
        if let (Some(a), Some(b)) = (&mut self.cylindrical, &o.cylindrical) {
            a.1.merge(&b.1);
        }
        for g in 0..2 {
            self.escaped_front[g] += o.escaped_front[g];
            self.escaped_back[g] += o.escaped_back[g];
            self.escaped_side[g] += o.escaped_side[g];
        }
    }

    pub(super) fn report(&self, histories: u64) -> ChargeReport {
        let mut total = ChargeCounts::default();
        for l in &self.layers.counts {
            total.merge(l);
        }
        let esc = |e: [u64; 2]| EscapedCharge {
            primary: e[0],
            secondary: e[1],
        };
        let escaped = self
            .escaped_front
            .iter()
            .chain(&self.escaped_back)
            .chain(&self.escaped_side)
            .sum::<u64>();
        ChargeReport {
            incident: histories,
            total,
            per_layer: self.layers.report(),
            escaped_front: esc(self.escaped_front),
            escaped_back: esc(self.escaped_back),
            escaped_side: esc(self.escaped_side),
            imbalance: histories as i64 + total.created() as i64
                - total.stopped() as i64
                - escaped as i64,
            cartesian: self.cartesian.as_ref().map(|(grid, a)| CartesianCharge {
                grid: *grid,
                cells: a.report(),
            }),
            cylindrical: self
                .cylindrical
                .as_ref()
                .map(|(grid, a)| CylindricalCharge {
                    grid: *grid,
                    cells: a.report(),
                }),
        }
    }
}

/// Where the electrons of a run came to rest and where holes were left, as
/// integer counts, and the balance against the incident electrons.
///
/// Charge is counted in units of the electron's, per electron number:
///
/// - every electron that ends inside the target counts `+1` where it ends:
///   every terminal fate except escape (stopped below the layer's threshold,
///   trapped without an interaction, absorbed at the back face, polaron
///   trapped, cut off by the collision cap);
/// - every liberated electron (a secondary-model event with
///   `liberated = true`) counts `-1` at the event, the hole it leaves, whether
///   or not the electron is followed. A liberated electron below the stopping
///   threshold is not followed, it is dropped and counted `+1` where it was
///   created, so its `-1` and `+1` cancel;
/// - an electron leaving through a face carries its charge out and is counted
///   in `escaped_front` / `escaped_back`.
///
/// The balance is exact in integers: `incident + created = stopped +
/// escaped`, so `imbalance` is zero on every run and the net electron number
/// deposited in the target, `stopped - created`, equals `incident - escaped`.
/// Generation: `primary` is the incident electron, `secondary` every electron
/// created by the secondary model. A secondary below the cutoff is not
/// followed and is therefore counted where it is dropped (see
/// [`ElectronTallyMetadata::cutoff_ev`](super::ElectronTallyMetadata::cutoff_ev)).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ChargeReport {
    /// Incident electrons (primary histories).
    pub incident: u64,
    /// Totals over the target.
    pub total: ChargeCounts,
    /// Per layer.
    pub per_layer: ChargeCells,
    /// Electrons that left through the front face.
    pub escaped_front: EscapedCharge,
    /// Electrons that left through the back face.
    pub escaped_back: EscapedCharge,
    /// Electrons that left through the lateral boundary of a non-slab
    /// geometry.
    pub escaped_side: EscapedCharge,
    /// `incident + created - stopped - escaped`; zero if every electron was
    /// counted exactly once.
    pub imbalance: i64,
    /// The Cartesian grid, if the run has one.
    pub cartesian: Option<CartesianCharge>,
    /// The cylindrical grid, if the run has one.
    pub cylindrical: Option<CylindricalCharge>,
}

/// Electrons that left through one face, by generation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EscapedCharge {
    /// Primaries.
    pub primary: u64,
    /// Secondaries.
    pub secondary: u64,
}

/// Charge counts per cell (or per layer) with the data for the standard
/// error of the net charge per primary.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChargeCells {
    /// Counts per cell.
    pub counts: Vec<ChargeCounts>,
    /// Per cell, the sum over histories of the square of the net charge
    /// (stopped minus created) the history left there.
    pub net_sum_sq: Vec<u64>,
    /// Counts outside the grid.
    pub outside: ChargeCounts,
    /// Squared net charge sum outside the grid.
    pub outside_net_sum_sq: u64,
}

impl ChargeCells {
    /// Net charge per incident electron in cell `i` and its standard error
    /// (of the mean over `histories`).
    pub fn net_per_primary(&self, i: usize, histories: u64) -> (f64, f64) {
        net_mean_se(self.counts[i].net(), self.net_sum_sq[i], histories)
    }
}

/// Mean and standard error of a per-history integer quantity from its sum
/// and sum of squares.
pub fn net_mean_se(sum: i64, sum_sq: u64, n: u64) -> (f64, f64) {
    if n == 0 {
        return (0.0, 0.0);
    }
    let nf = n as f64;
    let mean = sum as f64 / nf;
    if n == 1 {
        return (mean, 0.0);
    }
    let var = ((sum_sq as f64 - nf * mean * mean) / (nf - 1.0)).max(0.0);
    (mean, (var / nf).sqrt())
}

/// Charge on the Cartesian grid, in [`CartesianGrid::index`] order.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CartesianCharge {
    /// The grid.
    pub grid: CartesianGrid,
    /// Per voxel.
    pub cells: ChargeCells,
}

/// Charge on the cylindrical `(r, depth)` grid, in
/// [`CylindricalGrid::index`] order.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CylindricalCharge {
    /// The grid.
    pub grid: CylindricalGrid,
    /// Per cell.
    pub cells: ChargeCells,
}
