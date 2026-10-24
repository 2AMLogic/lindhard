//! Opt-in recorder of ion paths: the collision points of chosen histories.
//!
//! [`TrajectoryRecorder`] is a [`BcaTally`] like any other, so the engine pays
//! nothing for it unless a run uses it: no field on [`Particle`], no flag and
//! no branch in the transport loop. [`Recorded`] forwards every hook to an
//! ordinary tally as well, so the usual results come out of the same run and
//! are identical with recording on or off (the recorder draws no random
//! numbers and never touches the engine state).
//!
//! # Convention
//!
//! A history is a list of *tracks*: the primary (generation 0) and every
//! displaced atom (a recoil, generation = parent + 1). A track is an ordered
//! list of points, positions in metres (`x` the depth) and energies in eV:
//!
//! * a recoil starts with its state at creation;
//! * every collision partner adds the particle's **incoming** state at the
//!   collision (the [`BcaTally::partner`] convention). A crystal step with
//!   several simultaneous partners adds one point. The primary's own track
//!   therefore starts at its first collision, not at the surface;
//! * a track that stops or leaves the target ends with that state (for an
//!   escape, the state outside after refraction by the surface barrier), the
//!   same position the other tallies see. A recoil that is created below the
//!   recoil cutoff, or with cascades off, is stopped at creation and has that
//!   one point.
//!
//! Each non-primary track names its `parent` track in the same history. The
//! engine reports no track-start event; the recorder rebuilds the link by
//! mirroring the engine's last-in first-out queue of recoils: a recoil
//! reported and then stopped with the identical state was never queued, any
//! other is followed after its parent ends, the most recent first.
//!
//! # Bounds
//!
//! Memory is bounded by two caps. Only histories with index below
//! [`TrajectoryConfig::max_histories`] are recorded, so the lowest indices win
//! whatever the thread count; the others are counted in
//! [`TrajectoryRecorder::histories_skipped`]. A history keeps at most
//! [`TrajectoryConfig::max_points`] points over all its tracks; further points
//! (and the tracks of further recoils) are dropped and the history is flagged
//! [`History::truncated`]. Histories are merged in chunk order, so the output
//! is bit-identical at any thread count.

use crate::ion::bca::{BcaTally, ElectronicChannel, EnergyBudget, Face, LatticeDeposit, Particle};

/// Caps on what a [`TrajectoryRecorder`] keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrajectoryConfig {
    /// Record histories with index `< max_histories`.
    pub max_histories: u64,
    /// Keep at most this many points per history, over all its tracks.
    pub max_points: usize,
}

/// One point of a track.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrajectoryPoint {
    /// Position `[x, y, z]`, m.
    pub pos: [f64; 3],
    /// Kinetic energy, eV.
    pub energy_ev: f64,
}

/// How a track ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrackEnd {
    /// Came to rest in the target.
    Stopped,
    /// Left the target through a face.
    Escaped(Face),
}

/// The path of one particle.
#[derive(Debug, Clone, PartialEq)]
pub struct Track {
    /// Species index, as [`Particle::species`].
    pub species: usize,
    /// 0 for the primary, parent's generation + 1 for a recoil.
    pub generation: u32,
    /// Index in [`History::tracks`] of the track that displaced this atom;
    /// `None` for the primary.
    pub parent: Option<usize>,
    /// How the track ended; `None` if the end was dropped by the point cap.
    pub end: Option<TrackEnd>,
    /// The points, in order of occurrence.
    pub points: Vec<TrajectoryPoint>,
}

/// The recorded tracks of one primary history. `tracks[0]` is the primary;
/// every `parent` index is smaller than the track's own.
#[derive(Debug, Clone, PartialEq)]
pub struct History {
    /// Primary index.
    pub index: u64,
    /// Tracks in order of creation.
    pub tracks: Vec<Track>,
    /// True if the point cap dropped points or tracks of this history.
    pub truncated: bool,
}

/// A recoil reported but not yet followed, as the engine's queue holds it.
#[derive(Debug, Clone, Copy)]
struct Queued {
    track: Option<usize>,
    generation: u32,
    pos: [f64; 3],
    energy_ev: f64,
}

/// The track being followed. `track` is `None` once the point cap stopped
/// the recorder from keeping it (its events are still followed so the queue
/// mirror stays in step).
#[derive(Debug, Clone, Copy)]
struct Current {
    track: Option<usize>,
    generation: u32,
}

#[derive(Debug, Clone)]
struct Active {
    history: History,
    points: usize,
    queue: Vec<Queued>,
    cur: Current,
}

/// A [`BcaTally`] that records ion and recoil paths. See the module
/// documentation for the convention and the caps.
#[derive(Debug, Clone)]
pub struct TrajectoryRecorder {
    cfg: TrajectoryConfig,
    /// Recorded histories, by increasing index.
    pub histories: Vec<History>,
    /// Histories seen but not recorded because of `max_histories`.
    pub histories_skipped: u64,
    active: Option<Active>,
}

impl TrajectoryRecorder {
    /// Empty recorder with the given caps.
    pub fn new(cfg: TrajectoryConfig) -> Self {
        Self {
            cfg,
            histories: Vec::new(),
            histories_skipped: 0,
            active: None,
        }
    }

    /// The caps in force.
    pub fn config(&self) -> TrajectoryConfig {
        self.cfg
    }
}

impl Active {
    fn push_point(&mut self, max_points: usize, pos: [f64; 3], energy_ev: f64) {
        let Some(t) = self.cur.track else { return };
        let track = &mut self.history.tracks[t];
        let pt = TrajectoryPoint { pos, energy_ev };
        // The partners of one crystal step share the incoming state.
        if track.points.last() == Some(&pt) {
            return;
        }
        if self.points >= max_points {
            self.history.truncated = true;
            return;
        }
        self.points += 1;
        track.points.push(pt);
    }

    /// The current track ends; follow the most recently queued recoil next.
    fn finish(&mut self, max_points: usize, p: &Particle, end: TrackEnd) {
        if let Some(t) = self.cur.track {
            let before = self.history.tracks[t].points.len();
            self.push_point(max_points, p.pos, p.energy_ev);
            let track = &mut self.history.tracks[t];
            // A dropped end point leaves the end unknown.
            if track.points.len() > before || track.points.last().map(|l| l.pos) == Some(p.pos) {
                track.end = Some(end);
            }
        }
        self.cur = match self.queue.pop() {
            Some(q) => Current {
                track: q.track,
                generation: q.generation,
            },
            None => Current {
                track: None,
                generation: 0,
            },
        };
    }

    /// A recoil that is stopped right after being reported was never queued.
    fn take_unqueued(&mut self, p: &Particle) -> bool {
        match self.queue.last() {
            Some(q)
                if q.generation == p.generation
                    && q.generation == self.cur.generation + 1
                    && q.pos == p.pos
                    && q.energy_ev == p.energy_ev =>
            {
                let q = self.queue.pop().expect("checked non-empty");
                if let Some(t) = q.track {
                    self.history.tracks[t].end = Some(TrackEnd::Stopped);
                }
                true
            }
            _ => false,
        }
    }
}

impl BcaTally for TrajectoryRecorder {
    fn begin_history(&mut self, index: u64) {
        if index >= self.cfg.max_histories {
            self.histories_skipped += 1;
            self.active = None;
            return;
        }
        self.active = Some(Active {
            history: History {
                index,
                tracks: vec![Track {
                    species: 0,
                    generation: 0,
                    parent: None,
                    end: None,
                    points: Vec::new(),
                }],
                truncated: false,
            },
            points: 0,
            queue: Vec::new(),
            cur: Current {
                track: Some(0),
                generation: 0,
            },
        });
    }

    fn partner(&mut self, p: &Particle, _impact_parameter_m: f64, _partners: usize) {
        let max = self.cfg.max_points;
        if let Some(a) = &mut self.active {
            a.push_point(max, p.pos, p.energy_ev);
        }
    }

    fn recoil(&mut self, r: &Particle) {
        let max = self.cfg.max_points;
        let Some(a) = &mut self.active else { return };
        let mut track = None;
        if let Some(parent) = a.cur.track {
            if a.points < max {
                a.history.tracks.push(Track {
                    species: r.species,
                    generation: r.generation,
                    parent: Some(parent),
                    end: None,
                    points: Vec::new(),
                });
                track = Some(a.history.tracks.len() - 1);
                // The new track is current only for the creation point.
                let saved = a.cur;
                a.cur.track = track;
                a.push_point(max, r.pos, r.energy_ev);
                a.cur = saved;
            } else {
                a.history.truncated = true;
            }
        }
        a.queue.push(Queued {
            track,
            generation: r.generation,
            pos: r.pos,
            energy_ev: r.energy_ev,
        });
    }

    fn stopped(&mut self, p: &Particle) {
        let max = self.cfg.max_points;
        let Some(a) = &mut self.active else { return };
        if !a.take_unqueued(p) {
            a.finish(max, p, TrackEnd::Stopped);
        }
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        let max = self.cfg.max_points;
        if let Some(a) = &mut self.active {
            a.finish(max, p, TrackEnd::Escaped(face));
        }
    }

    fn end_history(&mut self, _index: u64, _budget: &EnergyBudget) {
        if let Some(a) = self.active.take() {
            self.histories.push(a.history);
        }
    }

    fn merge(&mut self, other: Self) {
        self.histories.extend(other.histories);
        self.histories_skipped += other.histories_skipped;
    }
}

/// An ordinary tally and a [`TrajectoryRecorder`] fed by the same events of
/// one run. `inner` sees exactly what it would see alone.
#[derive(Debug, Clone)]
pub struct Recorded<T> {
    /// The tally whose results the run is for.
    pub inner: T,
    /// The recorder.
    pub recorder: TrajectoryRecorder,
}

impl<T: BcaTally> BcaTally for Recorded<T> {
    fn begin_history(&mut self, index: u64) {
        self.inner.begin_history(index);
        self.recorder.begin_history(index);
    }
    fn electronic(&mut self, p: &Particle, from: [f64; 3], c: ElectronicChannel, e: f64) {
        self.inner.electronic(p, from, c, e);
    }
    fn partner(&mut self, p: &Particle, b: f64, n: usize) {
        self.inner.partner(p, b, n);
        self.recorder.partner(p, b, n);
    }
    fn lattice(&mut self, at: [f64; 3], layer: usize, kind: LatticeDeposit, e: f64) {
        self.inner.lattice(at, layer, kind, e);
    }
    fn recoil(&mut self, r: &Particle) {
        self.inner.recoil(r);
        self.recorder.recoil(r);
    }
    fn stopped(&mut self, p: &Particle) {
        self.inner.stopped(p);
        self.recorder.stopped(p);
    }
    fn escaped(&mut self, p: &Particle, face: Face) {
        self.inner.escaped(p, face);
        self.recorder.escaped(p, face);
    }
    fn end_history(&mut self, index: u64, budget: &EnergyBudget) {
        self.inner.end_history(index, budget);
        self.recorder.end_history(index, budget);
    }
    fn merge(&mut self, other: Self) {
        self.inner.merge(other.inner);
        self.recorder.merge(other.recorder);
    }
}
