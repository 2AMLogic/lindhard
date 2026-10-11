//! Electron transport loop on synthetic cross-section tables: free-path
//! distribution, CSDA path length, interface attenuation, determinism across
//! thread counts and run metadata. The tables are computed here from simple
//! formulas (constant mean free paths, isotropic angles), not measured or
//! tabulated data.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{
    ElectronState, ElectronTally, EscapeRule, Face, Fate, LayerTables, Primary, SummaryTally,
    Transport, TransportConfig, TransportError,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::rng;
use rand_core::Rng;

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn parts(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    let quantiles = GRID
        .iter()
        .map(|&e| {
            if rate == 0.0 {
                vec![]
            } else {
                PROB.iter().map(|&p| row(e, p)).collect()
            }
        })
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/electron_transport.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// Constant inelastic mean free path, every loss exactly `dw` eV.
fn inelastic_const(lambda_m: f64, dw: f64) -> CrossSectionTable {
    parts(SamplingAxis::InelasticEnergyLoss, 1.0 / lambda_m, |_, _| dw)
}

/// Inelastic loss uniform on `[0, frac * E]`.
fn inelastic_frac(lambda_m: f64, frac: f64) -> CrossSectionTable {
    parts(
        SamplingAxis::InelasticEnergyLoss,
        1.0 / lambda_m,
        move |e, p| frac * e * p,
    )
}

fn elastic_none() -> CrossSectionTable {
    parts(SamplingAxis::ElasticPolarAngle, 0.0, |_, _| 0.0)
}

/// Isotropic elastic scattering: `cos θ = 1 - 2p`.
fn elastic_isotropic(lambda_m: f64) -> CrossSectionTable {
    parts(SamplingAxis::ElasticPolarAngle, 1.0 / lambda_m, |_, p| {
        (1.0 - 2.0 * p).acos()
    })
}

fn pair(elastic: CrossSectionTable, inelastic: CrossSectionTable) -> LayerTables {
    LayerTables { elastic, inelastic }
}

fn uniform(r: &mut impl Rng) -> f64 {
    (r.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

// ---------------------------------------------------------------------------
// Free path
// ---------------------------------------------------------------------------

#[derive(Default)]
struct FirstFlight {
    lengths: Vec<f64>,
    seen_first: bool,
}

impl ElectronTally for FirstFlight {
    fn begin_history(&mut self, _i: u64, _s: &ElectronState) {
        self.seen_first = false;
    }
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, length_m: f64) {
        if !self.seen_first {
            self.seen_first = true;
            self.lengths.push(length_m);
        }
    }
    fn merge(&mut self, o: Self) {
        self.lengths.extend(o.lengths);
    }
}

#[test]
fn first_flight_distance_is_exponential_by_ks_at_1_percent() {
    let lambda = 3.0e-9;
    let stack = Stack::semi_infinite(material());
    let t = Transport::new(
        stack,
        vec![pair(elastic_none(), inelastic_const(lambda, 1.0))],
        // The first collision takes the electron below the cutoff, so a
        // history is one flight plus one collision.
        TransportConfig::new(99.5),
    )
    .unwrap();
    let n = 20_000u64;
    let run = t
        .run(
            0xF1EE,
            n,
            256,
            &Primary::normal(100.0),
            FirstFlight::default,
        )
        .unwrap();
    let mut s = run.tally.lengths;
    assert_eq!(s.len() as u64, n);
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let nf = n as f64;
    let d = s
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let cdf = 1.0 - (-x / lambda).exp();
            ((i as f64 + 1.0) / nf - cdf)
                .abs()
                .max((cdf - i as f64 / nf).abs())
        })
        .fold(0.0, f64::max);
    // Asymptotic KS critical value at alpha = 0.01.
    let crit = 1.6276 / nf.sqrt();
    assert!(d < crit, "KS statistic {d} >= critical {crit}");
}

// ---------------------------------------------------------------------------
// CSDA
// ---------------------------------------------------------------------------

#[derive(Default)]
struct PathLog {
    steps: Vec<f64>,
    collisions: u64,
    losses: Vec<f64>,
    final_energy: f64,
}

impl ElectronTally for PathLog {
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, length_m: f64) {
        self.steps.push(length_m);
    }
    fn inelastic(&mut self, _a: &ElectronState, w: f64) {
        self.collisions += 1;
        self.losses.push(w);
    }
    fn stopped(&mut self, at: &ElectronState) {
        self.final_energy = at.energy_ev;
    }
    fn merge(&mut self, _o: Self) {
        unreachable!("single history");
    }
}

/// Constant loss per collision, no elastic scattering, constant mean free
/// path: the path length to the cutoff is the sum of the drawn free paths, and
/// the number of them is the smallest `n` with `E0 - n dW < cutoff`. The drawn paths
/// are reproduced here from the documented draw order (path, then channel,
/// then `W`), independently of the engine's bookkeeping.
#[test]
fn csda_path_length_equals_sum_of_drawn_paths() {
    let (lambda, dw, e0, cutoff) = (2.5e-9, 1.0, 100.0, 10.0);
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![pair(elastic_none(), inelastic_const(lambda, dw))],
        TransportConfig::new(cutoff),
    )
    .unwrap();
    for index in 0..20u64 {
        let seed = 77;
        let mut log = PathLog::default();
        let mut r = rng::stream(seed, index);
        let fate = t
            .history(&mut log, &mut r, index, &Primary::normal(e0))
            .unwrap();
        assert_eq!(fate, Fate::Stopped);

        // Smallest n with E0 - n dW < cutoff.
        let n = ((e0 - cutoff) / dw).floor() as usize + 1;
        assert_eq!(log.collisions as usize, n);
        assert_eq!(log.steps.len(), n);
        assert!(log.losses.iter().all(|&w| w == dw));
        assert!((log.final_energy - (e0 - n as f64 * dw)).abs() < 1e-12);

        // Replay the stream: one path draw, one channel draw, one W draw each.
        let mut replay = rng::stream(seed, index);
        let mut analytic = 0.0;
        for _ in 0..n {
            analytic += -(1.0 - uniform(&mut replay)).ln() / (1.0 / lambda);
            uniform(&mut replay); // channel
            uniform(&mut replay); // W
        }
        let recorded: f64 = log.steps.iter().sum();
        assert!(
            (recorded - analytic).abs() <= 1e-9 * analytic,
            "history {index}: recorded {recorded} vs analytic {analytic}"
        );
    }
}

// ---------------------------------------------------------------------------
// Interface
// ---------------------------------------------------------------------------

#[derive(Default)]
struct Attenuation {
    collisions: u64,
    crossed_clean: u64,
    through_clean: u64,
    histories: u64,
}

impl ElectronTally for Attenuation {
    fn begin_history(&mut self, _i: u64, _s: &ElectronState) {
        self.collisions = 0;
        self.histories += 1;
    }
    fn elastic(&mut self, _a: &ElectronState, _t: f64) {
        self.collisions += 1;
    }
    fn inelastic(&mut self, _a: &ElectronState, _w: f64) {
        self.collisions += 1;
    }
    fn interface(&mut self, at: &ElectronState, from: usize, to: usize) {
        assert_eq!((from, to), (0, 1));
        assert_eq!(at.layer, 1);
        if self.collisions == 0 {
            self.crossed_clean += 1;
        }
    }
    fn escaped(&mut self, _at: &ElectronState, face: Face) {
        assert_eq!(face, Face::Back);
        if self.collisions == 0 {
            self.through_clean += 1;
        }
    }
    fn merge(&mut self, o: Self) {
        self.crossed_clean += o.crossed_clean;
        self.through_clean += o.through_clean;
        self.histories += o.histories;
    }
}

#[test]
fn two_layer_attenuation_matches_exp_minus_d_over_lambda() {
    let (d1, l1) = (1.0e-8, 2.0e-8);
    let (d2, l2) = (1.5e-8, 1.0e-8);
    let t = Transport::new(
        Stack::new(vec![(material(), d1), (material(), d2)], None).unwrap(),
        vec![
            pair(elastic_none(), inelastic_const(l1, 1.0)),
            pair(elastic_none(), inelastic_const(l2, 1.0)),
        ],
        // Any collision ends the history; only clean flights matter here.
        TransportConfig::new(99.5),
    )
    .unwrap();
    let n = 60_000u64;
    let a = t
        .run(5, n, 512, &Primary::normal(100.0), Attenuation::default)
        .unwrap()
        .tally;
    assert_eq!(a.histories, n);
    let nf = n as f64;
    let p1 = (-d1 / l1).exp();
    let sigma1 = (p1 * (1.0 - p1) / nf).sqrt();
    let got1 = a.crossed_clean as f64 / nf;
    assert!((got1 - p1).abs() < 4.0 * sigma1, "layer 1: {got1} vs {p1}");

    // Redraw at the interface: the survival through layer 2 depends on l2.
    let p2 = (-d2 / l2).exp();
    let m = a.crossed_clean as f64;
    let sigma2 = (p2 * (1.0 - p2) / m).sqrt();
    let got2 = a.through_clean as f64 / m;
    assert!((got2 - p2).abs() < 4.0 * sigma2, "layer 2: {got2} vs {p2}");
}

#[test]
fn front_only_rule_absorbs_at_the_back_face() {
    let mut cfg = TransportConfig::new(1.0);
    cfg.escape_rule = EscapeRule::FrontOnly;
    let t = Transport::new(
        Stack::new(vec![(material(), 1e-9)], None).unwrap(),
        vec![pair(elastic_none(), inelastic_const(1.0, 1e-3))],
        cfg,
    )
    .unwrap();
    let r = t
        .run(1, 100, 16, &Primary::normal(100.0), SummaryTally::default)
        .unwrap();
    assert_eq!(r.tally.absorbed, 100);
    assert_eq!(r.tally.escaped_back, 0);
    assert_eq!(r.metadata.escape_rule, EscapeRule::FrontOnly);
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// Records every event with the exact bits of its numbers.
#[derive(Default)]
struct EventLog {
    index: u64,
    events: Vec<(u64, u8, [u64; 5])>,
    summary: SummaryTally,
}

impl EventLog {
    fn push(&mut self, kind: u8, s: &ElectronState, extra: f64) {
        self.events.push((
            self.index,
            kind,
            [
                s.pos[0].to_bits(),
                s.pos[1].to_bits(),
                s.dir[2].to_bits(),
                s.energy_ev.to_bits(),
                extra.to_bits(),
            ],
        ));
    }
}

impl ElectronTally for EventLog {
    fn begin_history(&mut self, i: u64, _s: &ElectronState) {
        self.index = i;
    }
    fn step(&mut self, _f: [f64; 3], e: &ElectronState, l: f64) {
        self.push(0, e, l);
        self.summary.step(_f, e, l);
    }
    fn elastic(&mut self, a: &ElectronState, t: f64) {
        self.push(1, a, t);
        self.summary.elastic(a, t);
    }
    fn inelastic(&mut self, a: &ElectronState, w: f64) {
        self.push(2, a, w);
        self.summary.inelastic(a, w);
    }
    fn interface(&mut self, a: &ElectronState, f: usize, t: usize) {
        self.push(3, a, (f * 10 + t) as f64);
        self.summary.interface(a, f, t);
    }
    fn stopped(&mut self, a: &ElectronState) {
        self.push(4, a, 0.0);
        self.summary.stopped(a);
    }
    fn escaped(&mut self, a: &ElectronState, face: Face) {
        self.push(5, a, if face == Face::Front { 0.0 } else { 1.0 });
        self.summary.escaped(a, face);
    }
    fn end_history(&mut self, i: u64, f: Fate) {
        self.summary.end_history(i, f);
    }
    fn merge(&mut self, o: Self) {
        self.events.extend(o.events);
        self.summary.merge(o.summary);
    }
}

fn mixed_transport() -> Transport {
    let tables = vec![
        pair(elastic_isotropic(4e-9), inelastic_frac(6e-9, 0.2)),
        pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.3)),
        pair(elastic_isotropic(5e-9), inelastic_frac(9e-9, 0.1)),
    ];
    Transport::new(
        Stack::new(
            vec![(material(), 8e-9), (material(), 1.2e-8)],
            Some(material()),
        )
        .unwrap(),
        tables,
        TransportConfig::new(5.0),
    )
    .unwrap()
}

fn run_threads(threads: usize) -> EventLog {
    let t = mixed_transport();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let primary = Primary {
        energy_ev: 500.0,
        direction: [1.0, 0.3, -0.2],
    };
    pool.install(|| {
        t.run(0xBEEF, 600, 7, &primary, EventLog::default)
            .unwrap()
            .tally
    })
}

#[test]
fn event_streams_and_tallies_are_bit_identical_across_thread_counts() {
    let r1 = run_threads(1);
    assert!(r1.events.len() > 10_000, "run too small to mean anything");
    assert!(r1.summary.interface_crossings > 0);
    assert!(r1.summary.escaped_front > 0 && r1.summary.stopped > 0);
    for n in [2, 8] {
        let r = run_threads(n);
        assert_eq!(r1.events, r.events, "event stream differs on {n} threads");
        let (a, b) = (&r1.summary, &r.summary);
        assert_eq!(a.path_m.to_bits(), b.path_m.to_bits());
        assert_eq!(a.inelastic_loss_ev.to_bits(), b.inelastic_loss_ev.to_bits());
        assert_eq!(a.escaped_energy_ev.to_bits(), b.escaped_energy_ev.to_bits());
        assert_eq!(a.rest_energy_ev.to_bits(), b.rest_energy_ev.to_bits());
        assert_eq!(a, b);
    }
}

#[test]
fn energy_is_conserved_in_the_mixed_run() {
    let s = run_threads(2).summary;
    assert_eq!(s.histories, 600);
    let accounted = s.inelastic_loss_ev + s.escaped_energy_ev + s.rest_energy_ev;
    let incident = 600.0 * 500.0;
    assert!((accounted - incident).abs() < 1e-9 * incident);
}

// ---------------------------------------------------------------------------
// Metadata and validation
// ---------------------------------------------------------------------------

#[test]
fn cutoff_and_escape_rule_are_in_the_run_metadata() {
    let mut cfg = TransportConfig::new(12.5);
    cfg.escape_rule = EscapeRule::FrontOnly;
    cfg.max_events = 1234;
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![pair(elastic_none(), inelastic_const(1e-9, 1.0))],
        cfg,
    )
    .unwrap();
    let r = t
        .run(9, 10, 3, &Primary::normal(50.0), SummaryTally::default)
        .unwrap();
    let m = &r.metadata;
    assert_eq!(m.cutoff_ev, 12.5);
    assert_eq!(m.escape_rule, EscapeRule::FrontOnly);
    assert_eq!(m.max_events, 1234);
    assert_eq!((m.seed, m.n_histories, m.chunk_size), (9, 10, 3));
    assert_eq!(m.layers.len(), 1);
    assert_eq!(m.layers[0].back_m, None);
    assert_eq!(m.layers[0].inelastic_model, "synthetic");
    assert!(m.layers[0]
        .elastic_provenance
        .contains("tests/electron_transport.rs"));
    let json = serde_json::to_string(m).unwrap();
    assert!(json.contains("\"cutoff_ev\":12.5") && json.contains("\"front-only\""));
}

#[test]
fn rejects_inconsistent_inputs() {
    let stack = || Stack::semi_infinite(material());
    let ok = || pair(elastic_none(), inelastic_const(1e-9, 1.0));
    assert!(matches!(
        Transport::new(stack(), vec![], TransportConfig::new(1.0)),
        Err(TransportError::LayerCount { .. })
    ));
    let swapped = pair(inelastic_const(1e-9, 1.0), elastic_none());
    assert!(matches!(
        Transport::new(stack(), vec![swapped], TransportConfig::new(1.0)),
        Err(TransportError::WrongAxis { layer: 0, .. })
    ));
    for c in [0.0, -1.0, f64::NAN] {
        assert!(Transport::new(stack(), vec![ok()], TransportConfig::new(c)).is_err());
    }
    let t = Transport::new(stack(), vec![ok()], TransportConfig::new(10.0)).unwrap();
    assert!(t
        .run(1, 1, 1, &Primary::normal(10.0), SummaryTally::default)
        .is_err());
    let zero_dir = Primary {
        energy_ev: 100.0,
        direction: [0.0; 3],
    };
    assert!(t.run(1, 1, 1, &zero_dir, SummaryTally::default).is_err());
}

#[test]
fn zero_rate_heading_inward_is_trapped_and_outward_escapes() {
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![pair(
            elastic_none(),
            parts(SamplingAxis::InelasticEnergyLoss, 0.0, |_, _| 0.0),
        )],
        TransportConfig::new(1.0),
    )
    .unwrap();
    let r = t
        .run(1, 4, 2, &Primary::normal(100.0), SummaryTally::default)
        .unwrap();
    assert_eq!(r.tally.trapped, 4);
    let out = Primary {
        energy_ev: 100.0,
        direction: [-1.0, 0.0, 0.0],
    };
    let r = t.run(1, 4, 2, &out, SummaryTally::default).unwrap();
    assert_eq!(r.tally.escaped_front, 4);
}
