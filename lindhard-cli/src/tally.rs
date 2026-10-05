//! The tally the CLI runs today: the library's [`SummaryTally`] (counts,
//! energy budget, stopped-primary depth histogram) plus the final state of
//! every primary. The range-moment, damage and sputtering tallies of the
//! library are added alongside it as they land, each as a new output section.

use lindhard::ion::bca::{BcaTally, EnergyBudget, Face, Particle, SummaryTally};

/// How a primary history ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// Came to rest in the target.
    Stopped,
    /// Left through the front face.
    Backscattered,
    /// Left through the back face.
    Transmitted,
}

impl Fate {
    pub fn as_str(self) -> &'static str {
        match self {
            Fate::Stopped => "stopped",
            Fate::Backscattered => "backscattered",
            Fate::Transmitted => "transmitted",
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

/// Summary plus per-primary final states. `merge` appends in chunk order, so
/// `finals` is sorted by history index and identical at any thread count.
pub struct CliTally {
    pub summary: SummaryTally,
    pub per_ion: bool,
    pub finals: Vec<FinalState>,
    current: u64,
}

impl CliTally {
    pub fn new(bin_width_m: f64, n_bins: usize, per_ion: bool) -> Self {
        Self {
            summary: SummaryTally::new(bin_width_m, n_bins),
            per_ion,
            finals: Vec::new(),
            current: 0,
        }
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
    }

    fn recoil(&mut self, r: &Particle) {
        self.summary.recoil(r);
    }

    fn stopped(&mut self, p: &Particle) {
        self.summary.stopped(p);
        self.record(p, Fate::Stopped);
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        self.summary.escaped(p, face);
        let fate = match face {
            Face::Front => Fate::Backscattered,
            Face::Back => Fate::Transmitted,
        };
        self.record(p, fate);
    }

    fn end_history(&mut self, index: u64, budget: &EnergyBudget) {
        self.summary.end_history(index, budget);
    }

    fn merge(&mut self, other: Self) {
        self.summary.merge(other.summary);
        self.finals.extend(other.finals);
    }
}
