//! Internal non-charge accumulators of the electron tally: energy-weighted
//! moments, the energy budget, the emission classes per face, the deposition
//! grid scratch and the sparse voxel deposition. Nothing here is public; the
//! report types these fill live in [`super::report`] and the charge
//! bookkeeping in [`super::charge`].

use std::collections::BTreeMap;

use super::report::{ElectronEnergyBudget, EmissionClass, FaceEmission};
use super::ElectronTallyConfig;
use crate::tally::hist::Histogram;
use crate::tally::voxel::{VoxelCell, VoxelTallyGrid};

/// Energy-weighted mean and centred second moment of one coordinate,
/// mergeable. The pairwise update is that of T. F. Chan, G. H. Golub and
/// R. J. LeVeque, Stanford report STAN-CS-79-773 (1979) (the update also used
/// by [`moments`](crate::tally::moments)), with the counts replaced by weights, which is the
/// weighted form of D. H. D. West, "Updating mean and variance estimates: an
/// improved method", Commun. ACM 22, 532 (1979), doi:10.1145/359146.359153.
/// With `W = W_a + W_b` and `δ = mean_b - mean_a`:
/// `mean = mean_a + δ W_b / W`, `M2 = M2_a + M2_b + δ² W_a W_b / W`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub(super) struct WeightedMoments {
    pub(super) weight: f64,
    pub(super) mean: f64,
    pub(super) m2: f64,
}

impl WeightedMoments {
    pub(super) fn push(&mut self, x: f64, w: f64) {
        self.merge_parts(w, x, 0.0);
    }

    pub(super) fn merge(&mut self, o: &WeightedMoments) {
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

    pub(super) fn variance(&self) -> f64 {
        if self.weight > 0.0 {
            (self.m2 / self.weight).max(0.0)
        } else {
            0.0
        }
    }
}

/// Energy sums, eV.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub(super) struct Budget {
    pub(super) incident: f64,
    pub(super) inelastic: f64,
    pub(super) residual: f64,
    pub(super) escaped_front: f64,
    pub(super) escaped_back: f64,
    pub(super) escaped_side: f64,
    pub(super) no_interaction: f64,
    pub(super) absorbed: f64,
    pub(super) event_cap: f64,
    pub(super) fermi_sea: f64,
    pub(super) barrier: f64,
    pub(super) phonon_emitted: f64,
    pub(super) phonon_absorbed: f64,
    pub(super) polaron: f64,
}

impl Budget {
    pub(super) fn merge(&mut self, o: &Budget) {
        self.incident += o.incident;
        self.inelastic += o.inelastic;
        self.residual += o.residual;
        self.escaped_front += o.escaped_front;
        self.escaped_back += o.escaped_back;
        self.escaped_side += o.escaped_side;
        self.no_interaction += o.no_interaction;
        self.absorbed += o.absorbed;
        self.event_cap += o.event_cap;
        self.fermi_sea += o.fermi_sea;
        self.barrier += o.barrier;
        self.phonon_emitted += o.phonon_emitted;
        self.phonon_absorbed += o.phonon_absorbed;
        self.polaron += o.polaron;
    }

    pub(super) fn report(&self) -> ElectronEnergyBudget {
        let deposited = self.inelastic + self.residual + self.phonon_emitted + self.polaron;
        let escaped = self.escaped_front + self.escaped_back + self.escaped_side;
        let trapped = self.no_interaction + self.absorbed + self.event_cap;
        let source = self.incident + self.fermi_sea + self.phonon_absorbed;
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
            escaped_side_ev: self.escaped_side,
            no_interaction_ev: self.no_interaction,
            absorbed_ev: self.absorbed,
            event_cap_ev: self.event_cap,
            fermi_sea_ev: self.fermi_sea,
            barrier_ev: self.barrier,
            phonon_emitted_ev: self.phonon_emitted,
            phonon_absorbed_ev: self.phonon_absorbed,
            polaron_ev: self.polaron,
            ionisation_ev: self.inelastic + self.residual,
            lattice_ev: self.phonon_emitted + self.polaron,
            relative_imbalance,
        }
    }
}

/// The channel a deposit is booked to (see
/// [the module docs](super#the-ionisation-and-lattice-channels)).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Channel {
    Ionisation,
    Lattice,
}

/// One accumulator per channel; the plain total is kept alongside so the
/// existing total maps are unchanged.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct GridAcc {
    pub(super) energy_ev: Vec<f64>,
    pub(super) ionisation_ev: Vec<f64>,
    pub(super) lattice_ev: Vec<f64>,
    pub(super) outside_ev: f64,
    pub(super) outside_ionisation_ev: f64,
    pub(super) outside_lattice_ev: f64,
}

impl GridAcc {
    pub(super) fn new(cells: usize) -> Self {
        Self {
            energy_ev: vec![0.0; cells],
            ionisation_ev: vec![0.0; cells],
            lattice_ev: vec![0.0; cells],
            outside_ev: 0.0,
            outside_ionisation_ev: 0.0,
            outside_lattice_ev: 0.0,
        }
    }

    pub(super) fn add(&mut self, cell: Option<usize>, e: f64, channel: Channel) {
        match cell {
            Some(i) => {
                self.energy_ev[i] += e;
                match channel {
                    Channel::Ionisation => self.ionisation_ev[i] += e,
                    Channel::Lattice => self.lattice_ev[i] += e,
                }
            }
            None => {
                self.outside_ev += e;
                match channel {
                    Channel::Ionisation => self.outside_ionisation_ev += e,
                    Channel::Lattice => self.outside_lattice_ev += e,
                }
            }
        }
    }

    pub(super) fn merge(&mut self, o: &GridAcc) {
        for (a, b) in self.energy_ev.iter_mut().zip(&o.energy_ev) {
            *a += b;
        }
        for (a, b) in self.ionisation_ev.iter_mut().zip(&o.ionisation_ev) {
            *a += b;
        }
        for (a, b) in self.lattice_ev.iter_mut().zip(&o.lattice_ev) {
            *a += b;
        }
        self.outside_ev += o.outside_ev;
        self.outside_ionisation_ev += o.outside_ionisation_ev;
        self.outside_lattice_ev += o.outside_lattice_ev;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct ClassAcc {
    pub(super) count: u64,
    pub(super) energy_ev: f64,
    pub(super) polar_hist: Histogram,
}

impl ClassAcc {
    fn new(c: &ElectronTallyConfig) -> Self {
        Self {
            count: 0,
            energy_ev: 0.0,
            polar_hist: Histogram::new(c.escape_polar),
        }
    }

    pub(super) fn merge(&mut self, o: &ClassAcc) {
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
pub(super) struct FaceAcc {
    pub(super) slow: ClassAcc,
    pub(super) fast: ClassAcc,
    pub(super) energy_hist: Histogram,
}

impl FaceAcc {
    pub(super) fn new(c: &ElectronTallyConfig) -> Self {
        Self {
            slow: ClassAcc::new(c),
            fast: ClassAcc::new(c),
            energy_hist: Histogram::new(c.escape_energy),
        }
    }

    pub(super) fn merge(&mut self, o: &FaceAcc) {
        self.slow.merge(&o.slow);
        self.fast.merge(&o.fast);
        self.energy_hist.merge(&o.energy_hist);
    }

    pub(super) fn report(&self, histories: u64) -> FaceEmission {
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

/// Which [`VoxelCell`] field a deposit goes to.
#[derive(Debug, Clone, Copy)]
pub(super) enum DepositKind {
    /// Left at an inelastic event: `electronic_ev`.
    Electronic,
    /// Given to the lattice (LO-phonon emission, polaron trapping):
    /// `lattice_ev`.
    Lattice,
    /// The energy of an electron that fell below the stopping threshold:
    /// `rest_ev`.
    Rest,
}

impl DepositKind {
    /// The ionisation/lattice channel this deposit is booked to. The
    /// sub-cutoff residual goes to ionisation by convention.
    pub(super) fn channel(self) -> Channel {
        match self {
            DepositKind::Electronic | DepositKind::Rest => Channel::Ionisation,
            DepositKind::Lattice => Channel::Lattice,
        }
    }
}

/// Sparse per-voxel energy deposition on a [`VoxelTallyGrid`], keyed on the
/// flat voxel index (`BTreeMap`, so merging and export follow one fixed
/// order).
#[derive(Debug, Clone, PartialEq)]
pub(super) struct VoxelAcc {
    pub(super) grid: VoxelTallyGrid,
    pub(super) cells: BTreeMap<u64, VoxelCell>,
    pub(super) outside: VoxelCell,
}

impl VoxelAcc {
    pub(super) fn new(grid: VoxelTallyGrid) -> Self {
        Self {
            grid,
            cells: BTreeMap::new(),
            outside: VoxelCell::default(),
        }
    }

    pub(super) fn add(&mut self, pos: [f64; 3], e: f64, kind: DepositKind) {
        let cell = match self.grid.locate(pos) {
            Some(ijk) => self.cells.entry(self.grid.flat_index(ijk)).or_default(),
            None => &mut self.outside,
        };
        match kind {
            DepositKind::Electronic => cell.electronic_ev += e,
            DepositKind::Lattice => cell.lattice_ev += e,
            DepositKind::Rest => cell.rest_ev += e,
        }
    }

    pub(super) fn merge(&mut self, o: &VoxelAcc) {
        assert!(
            self.grid == o.grid,
            "merging electron tallies built on different voxel grids"
        );
        self.outside.add(&o.outside);
        for (&flat, cell) in &o.cells {
            self.cells.entry(flat).or_default().add(cell);
        }
    }
}

pub(super) fn per(x: f64, n: u64) -> f64 {
    if n == 0 {
        0.0
    } else {
        x / n as f64
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
}
