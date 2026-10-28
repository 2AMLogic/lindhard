//! The trajectory recorder (`lindhard::tally::trajectory`): endpoints agree
//! with the other tallies, recording changes no result, the caps flag
//! truncation, and the output does not depend on the thread count.

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, EnergyBudget, Face, Particle, SummaryTally,
};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{History, Recorded, TrackEnd, TrajectoryConfig, TrajectoryRecorder};

fn table() -> &'static ScatteringTable {
    static T: OnceLock<ScatteringTable> = OnceLock::new();
    T.get_or_init(|| {
        ScatteringTable::build(
            &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e4,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 16,
            },
        )
    })
}

fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    m
}

fn beam(count: u64) -> Beam {
    Beam {
        ion: Ion::new(5).unwrap(),
        energy_ev: 2.0e3,
        polar_rad: 0.3,
        azimuth_rad: 0.0,
        count,
    }
}

fn config(follow_recoils: bool) -> BcaConfig {
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 91;
    cfg.chunk_size = 16;
    cfg.follow_recoils = follow_recoils;
    cfg
}

fn rec(max_histories: u64, max_points: usize) -> TrajectoryRecorder {
    TrajectoryRecorder::new(TrajectoryConfig {
        max_histories,
        max_points,
    })
}

fn run<T: BcaTally>(
    count: u64,
    follow: bool,
    threads: usize,
    make: impl Fn() -> T + Sync + Send,
) -> T {
    let ls = LindhardScharff::new();
    // A semi-infinite target: ions stop, and recoils and ions can sputter.
    let stack = Stack::new(vec![], Some(si())).unwrap();
    let bca = Bca::new(beam(count), &stack, config(follow), &ls, table()).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| bca.run(make).unwrap())
}

/// Records the stopping or escape state of each primary, for comparison.
#[derive(Default)]
struct Ends(Vec<(u64, TrackEnd, [f64; 3], f64)>);

struct EndsTally {
    ends: Ends,
    index: u64,
}

impl BcaTally for EndsTally {
    fn begin_history(&mut self, i: u64) {
        self.index = i;
    }
    fn stopped(&mut self, p: &Particle) {
        if p.is_primary() {
            self.ends
                .0
                .push((self.index, TrackEnd::Stopped, p.pos, p.energy_ev));
        }
    }
    fn escaped(&mut self, p: &Particle, f: Face) {
        if p.is_primary() {
            self.ends
                .0
                .push((self.index, TrackEnd::Escaped(f), p.pos, p.energy_ev));
        }
    }
    fn merge(&mut self, o: Self) {
        self.ends.0.extend(o.ends.0);
    }
}

fn pair(max_histories: u64, max_points: usize) -> Recorded<EndsTally> {
    Recorded {
        inner: EndsTally {
            ends: Ends::default(),
            index: 0,
        },
        recorder: rec(max_histories, max_points),
    }
}

fn structure_ok(h: &History) {
    assert_eq!(h.tracks[0].generation, 0);
    assert_eq!(h.tracks[0].parent, None);
    for (i, t) in h.tracks.iter().enumerate().skip(1) {
        let parent = t.parent.expect("recoil has a parent");
        assert!(parent < i);
        assert_eq!(t.generation, h.tracks[parent].generation + 1);
        if !h.truncated {
            // The recoil is born at a point the engine reported for its parent.
            let born = t.points[0].pos;
            let pts = &h.tracks[parent].points;
            assert!(pts.iter().any(|p| p.pos == born), "history {}", h.index);
        }
    }
}

#[test]
fn last_point_is_the_reported_end() {
    let n = 300;
    let r = run(n, true, 2, || pair(n, usize::MAX));
    assert_eq!(r.recorder.histories.len() as u64, n);
    assert_eq!(r.recorder.histories_skipped, 0);
    assert_eq!(r.inner.ends.0.len() as u64, n);
    let (mut stopped, mut escaped, mut cascades) = (0, 0, 0);
    for (h, e) in r.recorder.histories.iter().zip(&r.inner.ends.0) {
        structure_ok(h);
        assert!(!h.truncated);
        assert_eq!(h.index, e.0);
        let t = &h.tracks[0];
        assert_eq!(t.end, Some(e.1));
        let last = t.points.last().unwrap();
        assert_eq!(last.pos, e.2);
        assert_eq!(last.energy_ev, e.3);
        match e.1 {
            TrackEnd::Stopped => stopped += 1,
            TrackEnd::Escaped(_) => escaped += 1,
        }
        cascades += h.tracks.len() - 1;
        for t in &h.tracks {
            assert!(t.end.is_some() && !t.points.is_empty());
        }
    }
    assert!(
        stopped > 0 && cascades > 0,
        "{stopped} {escaped} {cascades}"
    );
}

#[test]
fn cascades_off_makes_single_point_recoils() {
    let r = run(100, false, 1, || pair(100, usize::MAX));
    let mut recoils = 0;
    for h in &r.recorder.histories {
        structure_ok(h);
        for t in h.tracks.iter().skip(1) {
            recoils += 1;
            assert_eq!(t.points.len(), 1);
            assert_eq!(t.end, Some(TrackEnd::Stopped));
            assert_eq!(t.parent, Some(0));
        }
    }
    assert!(recoils > 0);
}

#[test]
fn recording_changes_no_tally() {
    let n = 200;
    let plain = run(n, true, 2, || SummaryTally::new(1e-9, 50));
    let both = run(n, true, 2, || Recorded {
        inner: SummaryTally::new(1e-9, 50),
        recorder: rec(n, usize::MAX),
    });
    assert_eq!(plain, both.inner);
}

#[test]
fn output_is_identical_across_thread_counts() {
    let n = 200;
    let r1 = run(n, true, 1, || rec(n, 40));
    for threads in [2, 8] {
        let r = run(n, true, threads, || rec(n, 40));
        assert_eq!(r1.histories, r.histories, "{threads} threads");
    }
}

#[test]
fn history_cap_keeps_the_lowest_indices() {
    let r = run(100, true, 4, || rec(7, usize::MAX));
    let idx: Vec<u64> = r.histories.iter().map(|h| h.index).collect();
    assert_eq!(idx, (0..7).collect::<Vec<_>>());
    assert_eq!(r.histories_skipped, 93);
    let none = run(20, true, 1, || rec(0, usize::MAX));
    assert!(none.histories.is_empty());
    assert_eq!(none.histories_skipped, 20);
}

#[test]
fn point_cap_flags_truncation() {
    let full = run(50, true, 1, || rec(50, usize::MAX));
    let total = |h: &History| h.tracks.iter().map(|t| t.points.len()).sum::<usize>();
    let capped = run(50, true, 1, || rec(50, 1));
    for (c, f) in capped.histories.iter().zip(&full.histories) {
        assert!(total(c) <= 1);
        assert_eq!(c.truncated, total(f) > 1);
        structure_ok(c);
    }
    assert!(capped.histories.iter().any(|h| h.truncated));
    let mid = run(50, true, 1, || rec(50, 10));
    for (c, f) in mid.histories.iter().zip(&full.histories) {
        assert!(total(c) <= 10);
        assert_eq!(c.truncated, total(f) > 10);
        structure_ok(c);
    }
}

#[test]
fn empty_history_has_a_primary_track() {
    let mut r = rec(1, 1);
    r.begin_history(0);
    r.end_history(0, &EnergyBudget::default());
    assert_eq!(r.histories.len(), 1);
    assert_eq!(r.histories[0].tracks.len(), 1);
}
