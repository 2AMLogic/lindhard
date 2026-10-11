//! The tally the CLI runs: the library's [`SummaryTally`] (counts, energy
//! budget, stopped-primary depth histogram), the library's [`IonTally`]
//! (range moments and fits, damage, escapes) and the final state of every
//! primary. Every event goes to all of them.

use lindhard::input::Resolved;
use lindhard::ion::bca::{BcaTally, EnergyBudget, Face, LatticeDeposit, Particle, SummaryTally};
use lindhard::tally::{Binning, IonTally, IonTallyConfig};

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
        self.ion.begin_history(index);
    }

    fn lattice(&mut self, at: [f64; 3], layer: usize, kind: LatticeDeposit, energy_ev: f64) {
        self.summary.lattice(at, layer, kind, energy_ev);
        self.ion.lattice(at, layer, kind, energy_ev);
    }

    fn recoil(&mut self, r: &Particle) {
        self.summary.recoil(r);
        self.ion.recoil(r);
    }

    fn stopped(&mut self, p: &Particle) {
        self.summary.stopped(p);
        self.ion.stopped(p);
        self.record(p, Fate::Stopped);
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        self.summary.escaped(p, face);
        self.ion.escaped(p, face);
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
    }

    fn merge(&mut self, other: Self) {
        self.summary.merge(other.summary);
        self.ion.merge(other.ion);
        self.finals.extend(other.finals);
    }
}

/// The [`IonTallyConfig`] for a resolved input. Depth uses the
/// `depth_bin_nm`/`depth_bins` grid from zero; lateral and radial use
/// `lateral_bin_nm`/`lateral_bins` as documented in `docs/cli.md`; escape
/// energies run from 0 to `escape_energy_max_ev` (default: the beam energy)
/// and polar angles from 0 to 90 degrees.
pub fn ion_tally_config(r: &Resolved) -> anyhow::Result<IonTallyConfig> {
    let t = &r.input.tally;
    let nm = 1e-9;
    let depth_w = t.depth_bin_nm * nm;
    let lat_half = t.lateral_bins as f64 * t.lateral_bin_nm * nm;
    let e_max = t.escape_energy_max_ev.unwrap_or(r.beam.energy_ev);
    let b = |lo: f64, hi: f64, n: usize| Binning::new(lo, hi, n);
    Ok(IonTallyConfig {
        depth: b(0.0, depth_w * t.depth_bins as f64, t.depth_bins)?,
        lateral: b(-lat_half, lat_half, 2 * t.lateral_bins)?,
        radial: b(0.0, lat_half, t.lateral_bins)?,
        escape_energy: b(0.0, e_max, t.escape_energy_bins)?,
        escape_polar: b(0.0, std::f64::consts::FRAC_PI_2, t.escape_polar_bins)?,
    })
}
