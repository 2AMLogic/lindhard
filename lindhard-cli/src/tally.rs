//! The tally the CLI runs: the library's [`SummaryTally`] (counts, energy
//! budget, stopped-primary depth histogram), the library's [`IonTally`]
//! (range moments and fits, damage, escapes) and the final state of every
//! primary. Every event goes to all of them.

use lindhard::ion::bca::{
    BcaTally, ElectronicChannel, EnergyBudget, Face, LatticeDeposit, Particle, SummaryTally,
};
use lindhard::tally::{IonTally, VoxelTally};

/// How a primary history ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// Came to rest in the target.
    Stopped,
    /// Left through the front face.
    Backscattered,
    /// Left through the back face.
    Transmitted,
    /// Left through a lateral face (voxel grids only).
    Lateral,
}

impl Fate {
    /// Every fate, in the order of [`Fate::code`].
    pub const ALL: [Fate; 4] = [
        Fate::Stopped,
        Fate::Backscattered,
        Fate::Transmitted,
        Fate::Lateral,
    ];

    /// Index of this fate in [`Fate::ALL`].
    pub fn code(self) -> u8 {
        match self {
            Fate::Stopped => 0,
            Fate::Backscattered => 1,
            Fate::Transmitted => 2,
            Fate::Lateral => 3,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Fate::Stopped => "stopped",
            Fate::Backscattered => "backscattered",
            Fate::Transmitted => "transmitted",
            Fate::Lateral => "lateral",
        }
    }
}

/// Final state of one primary.
#[derive(Debug, Clone, PartialEq)]
pub struct FinalState {
    pub index: u64,
    pub fate: Fate,
    /// Position, m (on the face for escaped primaries).
    pub pos: [f64; 3],
    /// Kinetic energy, eV (outside the target for escaped primaries).
    pub energy_ev: f64,
    pub dir: [f64; 3],
    pub layer: usize,
}

/// Summary, ion tally and per-primary final states. `merge` appends in chunk order, so
/// `finals` is sorted by history index and identical at any thread count.
#[derive(Debug)]
pub struct CliTally {
    pub summary: SummaryTally,
    pub ion: IonTally,
    /// The per-element voxel tally, when `[tally.voxel]` asks for it.
    pub voxel: Option<VoxelTally>,
    pub per_ion: bool,
    pub finals: Vec<FinalState>,
    current: u64,
}

impl CliTally {
    /// An empty tally; clones of `ion` (an empty [`IonTally`]) are the
    /// per-chunk tallies.
    pub fn new(bin_width_m: f64, n_bins: usize, per_ion: bool, ion: IonTally) -> Self {
        Self {
            summary: SummaryTally::new(bin_width_m, n_bins),
            ion,
            voxel: None,
            per_ion,
            finals: Vec::new(),
            current: 0,
        }
    }

    /// Also feed every event to `voxel` (an empty, species-aware tally).
    pub fn with_voxel(mut self, voxel: Option<VoxelTally>) -> Self {
        self.voxel = voxel;
        self
    }

    fn record(&mut self, p: &Particle, fate: Fate) {
        if self.per_ion && p.is_primary() {
            self.finals.push(FinalState {
                index: self.current,
                fate,
                pos: p.pos,
                energy_ev: p.energy_ev,
                dir: p.dir,
                layer: p.layer,
            });
        }
    }
}

impl BcaTally for CliTally {
    fn begin_history(&mut self, index: u64) {
        self.current = index;
        self.summary.begin_history(index);
        self.ion.begin_history(index);
        if let Some(v) = &mut self.voxel {
            v.begin_history(index);
        }
    }

    fn electronic(
        &mut self,
        p: &Particle,
        from: [f64; 3],
        channel: ElectronicChannel,
        energy_ev: f64,
    ) {
        self.ion.electronic(p, from, channel, energy_ev);
        if let Some(v) = &mut self.voxel {
            v.electronic(p, from, channel, energy_ev);
        }
    }

    fn lattice(&mut self, at: [f64; 3], layer: usize, kind: LatticeDeposit, energy_ev: f64) {
        self.summary.lattice(at, layer, kind, energy_ev);
        self.ion.lattice(at, layer, kind, energy_ev);
        if let Some(v) = &mut self.voxel {
            v.lattice(at, layer, kind, energy_ev);
        }
    }

    fn recoil(&mut self, r: &Particle) {
        self.summary.recoil(r);
        self.ion.recoil(r);
        if let Some(v) = &mut self.voxel {
            v.recoil(r);
        }
    }

    fn stopped(&mut self, p: &Particle) {
        self.summary.stopped(p);
        self.ion.stopped(p);
        if let Some(v) = &mut self.voxel {
            v.stopped(p);
        }
        self.record(p, Fate::Stopped);
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        self.summary.escaped(p, face);
        self.ion.escaped(p, face);
        if let Some(v) = &mut self.voxel {
            v.escaped(p, face);
        }
        let fate = match face {
            Face::Front => Fate::Backscattered,
            Face::Back => Fate::Transmitted,
            Face::Side => Fate::Lateral,
        };
        self.record(p, fate);
    }

    fn end_history(&mut self, index: u64, budget: &EnergyBudget) {
        self.summary.end_history(index, budget);
        self.ion.end_history(index, budget);
        if let Some(v) = &mut self.voxel {
            v.end_history(index, budget);
        }
    }

    fn merge(&mut self, other: Self) {
        self.summary.merge(other.summary);
        self.ion.merge(other.ion);
        if let (Some(a), Some(b)) = (&mut self.voxel, other.voxel) {
            a.merge(b);
        }
        self.finals.extend(other.finals);
    }
}
