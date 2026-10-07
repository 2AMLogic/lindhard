//! Secondary electrons and surface barriers in the electron loop, on
//! synthetic cross-section tables and synthetic band parameters (made-up
//! numbers chosen to exercise the code, not material data): per-event energy
//! conservation, the secondary's direction model, the step transmission
//! against the analytic formula, reflection below the barrier, the vacuum
//! cutoff, run metadata and determinism across thread counts.

use lindhard::electron::boundary::{BandModel, BandStructure};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::secondary::{SecondaryEvent, SecondaryModel};
use lindhard::electron::transport::{
    Boundary, BoundaryModel, CutoffReference, ElectronState, ElectronTally, EscapeRule, Face, Fate,
    LayerTables, Primary, SummaryTally, Transport, TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];
const PROVENANCE: &str = "synthetic, tests/electron_secondaries.rs";

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
        provenance: PROVENANCE.into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// Inelastic loss uniform on `[0, frac * E]`.
fn inelastic_frac(lambda_m: f64, frac: f64) -> CrossSectionTable {
    parts(
        SamplingAxis::InelasticEnergyLoss,
        1.0 / lambda_m,
        move |e, p| frac * e * p,
    )
}

fn inelastic_none() -> CrossSectionTable {
    parts(SamplingAxis::InelasticEnergyLoss, 0.0, |_, _| 0.0)
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

/// Synthetic metal: `E_F = 5 eV`, `Φ = 4 eV`, so `U = 9 eV`.
fn metal() -> BandStructure {
    BandStructure::new(
        BandModel::Metal {
            fermi_ev: 5.0,
            work_function_ev: 4.0,
        },
        PROVENANCE,
    )
    .unwrap()
}

/// Synthetic insulator: `W_v = 6 eV`, `E_g = 3 eV`, `χ = 1 eV`, so
/// `E_F = 7.5 eV` and `U = 10 eV`.
fn insulator() -> BandStructure {
    BandStructure::new(
        BandModel::Insulator {
            valence_band_width_ev: 6.0,
            band_gap_ev: 3.0,
            affinity_ev: 1.0,
        },
        PROVENANCE,
    )
    .unwrap()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---------------------------------------------------------------------------
// Energy conservation per event
// ---------------------------------------------------------------------------

/// Checks every inelastic event against the band parameters of its layer.
struct EventCheck {
    bands: Vec<BandStructure>,
    /// Energy and direction just before the current collision.
    before: (f64, [f64; 3]),
    w: f64,
    events: u64,
    liberated: u64,
    subgap: u64,
    created: u64,
    worst: f64,
}

impl EventCheck {
    fn new(bands: Vec<BandStructure>) -> Self {
        Self {
            bands,
            before: (0.0, [0.0; 3]),
            w: 0.0,
            events: 0,
            liberated: 0,
            subgap: 0,
            created: 0,
            worst: 0.0,
        }
    }
}

impl ElectronTally for EventCheck {
    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _l: f64) {
        self.before = (end.energy_ev, end.dir);
    }
    fn inelastic(&mut self, after: &ElectronState, w: f64) {
        assert_eq!(
            after.energy_ev,
            self.before.0 - w,
            "primary loses exactly W"
        );
        self.w = w;
    }
    fn secondary(&mut self, p: &ElectronState, e: &SecondaryEvent, c: Option<&ElectronState>) {
        let band = &self.bands[p.layer];
        let ef = band.fermi_ev();
        self.events += 1;
        assert_eq!(e.loss_ev, self.w);
        // Primary loss = secondary energy + binding + energy left in the
        // solid, to a few ulps of the terms involved.
        let scale = ef.abs() + e.loss_ev.abs() + band.band_gap_ev().unwrap_or(0.0);
        let residual = (e.loss_ev - (e.secondary_ev + e.binding_ev + e.deposited_ev)).abs();
        assert!(
            residual <= 4.0 * f64::EPSILON * scale,
            "residual {residual} for {e:?}"
        );
        self.worst = self.worst.max(residual / scale);
        match band.band_gap_ev() {
            Some(gap) if e.loss_ev <= gap => {
                self.subgap += 1;
                assert!(!e.liberated && c.is_none());
                assert_eq!((e.secondary_ev, e.binding_ev), (0.0, 0.0));
                assert_eq!(e.deposited_ev, e.loss_ev);
            }
            gap => {
                assert!(e.liberated, "{e:?}");
                self.liberated += 1;
                let b = gap.unwrap_or(0.0);
                assert_eq!(e.binding_ev, b - ef);
                // Verduin Eq. 3.86: E_SE = E_F + W - B.
                let e_se = ef + e.loss_ev - b;
                match c {
                    Some(s) => {
                        self.created += 1;
                        assert_eq!(e.secondary_ev, e_se);
                        assert_eq!(s.energy_ev, e_se);
                        assert_eq!((s.pos, s.layer), (p.pos, p.layer));
                        assert!((dot(s.dir, s.dir) - 1.0).abs() < 1e-12);
                        assert_eq!(e.deposited_ev, 0.0);
                    }
                    None => {
                        assert_eq!(e.secondary_ev, 0.0);
                        assert_eq!(e.deposited_ev, e_se);
                    }
                }
            }
        }
    }
    fn merge(&mut self, o: Self) {
        self.events += o.events;
        self.liberated += o.liberated;
        self.subgap += o.subgap;
        self.created += o.created;
        self.worst = self.worst.max(o.worst);
    }
}

fn metal_on_insulator(config: TransportConfig) -> Transport {
    Transport::with_band_structures(
        Stack::new(vec![(material(), 6e-9)], Some(material())).unwrap(),
        vec![
            pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3)),
            pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.25)),
        ],
        vec![metal(), insulator()],
        config,
    )
    .unwrap()
}

/// Secondaries and the step barrier, with the cutoff 0.5 eV above the vacuum
/// level (thresholds 9.5 eV in the metal, 10.5 eV in the insulator).
fn full_physics() -> TransportConfig {
    let mut cfg = TransportConfig::new(0.5);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    cfg
}

#[test]
fn secondary_energy_is_conserved_per_event_to_machine_precision() {
    let cfg = full_physics();
    let t = metal_on_insulator(cfg);
    let bands = vec![metal(), insulator()];
    let r = t
        .run(11, 400, 16, &Primary::normal(400.0), || {
            EventCheck::new(bands.clone())
        })
        .unwrap()
        .tally;
    assert!(r.events > 10_000, "too few events: {}", r.events);
    assert!(r.subgap > 100 && r.liberated > 1_000 && r.created > 1_000);
    assert!(
        r.created < r.liberated,
        "some secondaries fall below the cutoff"
    );
}

#[test]
fn run_energy_balances_with_secondaries_and_barriers() {
    let cfg = full_physics();
    let t = metal_on_insulator(cfg);
    let (n, e0) = (500u64, 300.0);
    let s = t
        .run(3, n, 32, &Primary::normal(e0), SummaryTally::default)
        .unwrap()
        .tally;
    assert_eq!(s.histories, n);
    assert_eq!(s.event_capped, 0);
    assert!(s.secondaries > 1_000 && s.secondaries_escaped > 0 && s.reflections > 0);
    let source = n as f64 * e0 + s.secondary_energy_ev;
    let sink = s.inelastic_loss_ev + s.escaped_energy_ev + s.rest_energy_ev + s.barrier_ev;
    assert!(
        (source - sink).abs() < 1e-9 * source,
        "source {source} vs sink {sink}"
    );
    // The per-event identity, summed.
    let events = s.secondary_energy_ev + s.binding_ev + s.deposited_ev;
    assert!((events - s.inelastic_loss_ev).abs() < 1e-9 * s.inelastic_loss_ev);
}

// ---------------------------------------------------------------------------
// Direction model
// ---------------------------------------------------------------------------

/// Checks the secondary's angle to the primary (Verduin Eqs. 3.100, 3.105,
/// 3.106) and the primary's deflection (Eq. 3.111).
struct AngleCheck {
    band: BandStructure,
    before: (f64, [f64; 3]),
    after_dir: [f64; 3],
    checked: u64,
}

impl ElectronTally for AngleCheck {
    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _l: f64) {
        self.before = (end.energy_ev, end.dir);
    }
    fn inelastic(&mut self, after: &ElectronState, _w: f64) {
        self.after_dir = after.dir;
    }
    fn secondary(&mut self, _p: &ElectronState, e: &SecondaryEvent, c: Option<&ElectronState>) {
        let Some(s) = c else { return };
        let (e_before, d) = self.before;
        let b = self.band.band_gap_ev().unwrap_or(0.0);
        let ratio = (e.loss_ev + b) / (e_before - self.band.fermi_ev() + 2.0 * b);
        if !(ratio > 0.0 && ratio < 1.0) {
            return;
        }
        let cos_b = ratio.sqrt();
        assert!(
            (dot(d, s.dir) - cos_b).abs() < 1e-12,
            "cos beta {} vs {cos_b}",
            dot(d, s.dir)
        );
        // Momentum conservation: the new primary direction is along
        // d - cos(beta) s.
        let p = [
            d[0] - cos_b * s.dir[0],
            d[1] - cos_b * s.dir[1],
            d[2] - cos_b * s.dir[2],
        ];
        let np = dot(p, p).sqrt();
        if np > 1e-6 {
            let q = self.after_dir;
            assert!((dot(q, p) / np - 1.0).abs() < 1e-10, "primary not along p");
        }
        self.checked += 1;
    }
    fn merge(&mut self, o: Self) {
        self.checked += o.checked;
    }
}

#[test]
fn secondary_direction_follows_the_binary_collision_angle() {
    for band in [metal(), insulator()] {
        let mut cfg = TransportConfig::new(10.0);
        // The instantaneous momentum is off so the angle is exact; it would
        // vanish anyway for the metal (B = 0).
        cfg.secondaries = SecondaryModel::KieftBosch {
            instantaneous_momentum: false,
            momentum_conservation: true,
        };
        let t = Transport::with_band_structures(
            Stack::semi_infinite(material()),
            vec![pair(elastic_isotropic(5e-9), inelastic_frac(4e-9, 0.3))],
            vec![band.clone()],
            cfg,
        )
        .unwrap();
        let r = t
            .run(21, 200, 8, &Primary::normal(500.0), || AngleCheck {
                band: band.clone(),
                before: (0.0, [0.0; 3]),
                after_dir: [0.0; 3],
                checked: 0,
            })
            .unwrap()
            .tally;
        assert!(r.checked > 1_000, "{:?}: only {} checked", band, r.checked);
    }
}

#[test]
fn instantaneous_momentum_keeps_directions_unit_and_spreads_the_angle() {
    // Insulator (B = E_g > 0): the added momentum moves the secondary off
    // the binary-collision cone.
    #[derive(Default)]
    struct Spread {
        before: [f64; 3],
        off_cone: u64,
        n: u64,
    }
    impl ElectronTally for Spread {
        fn step(&mut self, _f: [f64; 3], end: &ElectronState, _l: f64) {
            self.before = end.dir;
        }
        fn secondary(&mut self, p: &ElectronState, _e: &SecondaryEvent, c: Option<&ElectronState>) {
            assert!((dot(p.dir, p.dir) - 1.0).abs() < 1e-12);
            if let Some(s) = c {
                assert!((dot(s.dir, s.dir) - 1.0).abs() < 1e-12);
                self.n += 1;
                if dot(self.before, s.dir) < 0.0 {
                    self.off_cone += 1;
                }
            }
        }
        fn merge(&mut self, o: Self) {
            self.off_cone += o.off_cone;
            self.n += o.n;
        }
    }
    // Above the insulator's Fermi energy, 7.5 eV.
    let mut cfg = TransportConfig::new(8.0);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    let t = Transport::with_band_structures(
        Stack::semi_infinite(material()),
        vec![pair(elastic_none(), inelastic_frac(4e-9, 0.3))],
        vec![insulator()],
        cfg,
    )
    .unwrap();
    let r = t
        .run(5, 300, 10, &Primary::normal(300.0), Spread::default)
        .unwrap()
        .tally;
    // Without the added momentum cos(beta) >= 0 always; with it some
    // secondaries go backwards.
    assert!(r.n > 1_000 && r.off_cone > 0, "{} of {}", r.off_cone, r.n);
}

// ---------------------------------------------------------------------------
// Step barrier
// ---------------------------------------------------------------------------

/// Verduin Eq. 3.144, written with the wave numbers on the two sides,
/// `T = 4 k k' / (k + k')²`, `k ∝ sqrt(E_n)`, `k' ∝ sqrt(E_n + ΔU)`.
fn analytic_t(normal_ev: f64, du: f64) -> f64 {
    if normal_ev + du <= 0.0 {
        return 0.0;
    }
    let (k, k2) = (normal_ev.sqrt(), (normal_ev + du).sqrt());
    4.0 * k * k2 / ((k + k2) * (k + k2))
}

#[test]
fn entry_transmission_versus_normal_energy_matches_the_step_formula() {
    // A collisionless slab whose back face absorbs: a primary either gets in
    // (and is absorbed at the back) or is reflected at the front.
    let mut cfg = TransportConfig::new(0.01);
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.escape_rule = EscapeRule::FrontOnly;
    let t = Transport::with_band_structures(
        Stack::new(vec![(material(), 1e-8)], None).unwrap(),
        vec![pair(elastic_none(), inelastic_none())],
        vec![metal()],
        cfg,
    )
    .unwrap();
    let u = metal().inner_potential_ev();
    let n = 40_000u64;
    let cos = 0.6f64;
    for normal in [0.05, 0.3, 1.0, 4.0, 20.0] {
        let e = normal / (cos * cos);
        let primary = Primary {
            energy_ev: e,
            direction: [cos, (1.0 - cos * cos).sqrt(), 0.0],
        };
        let s = t
            .run(99, n, 500, &primary, SummaryTally::default)
            .unwrap()
            .tally;
        assert_eq!(s.absorbed + s.escaped_front, n);
        assert_eq!(s.reflections, s.escaped_front);
        let p = analytic_t(normal, u);
        let got = s.absorbed as f64 / n as f64;
        let sigma = (p * (1.0 - p) / n as f64).sqrt();
        assert!(
            (got - p).abs() < 4.0 * sigma,
            "E_n = {normal}: transmitted {got}, formula {p}"
        );
    }
}

/// Records, at every face, whether the electron had the normal energy to get
/// over and what happened.
#[derive(Default)]
struct FaceLog {
    inner: f64,
    at_face: Option<(f64, f64)>,
    below_reflected: u64,
    below_transmitted: u64,
    above_transmitted: u64,
}

impl ElectronTally for FaceLog {
    fn begin_history(&mut self, _i: u64, _s: &ElectronState) {
        // The entry from vacuum has no flight before it.
        self.at_face = None;
    }
    fn step(&mut self, _f: [f64; 3], end: &ElectronState, _l: f64) {
        self.at_face = Some((end.energy_ev, end.dir[0]));
    }
    fn barrier(&mut self, at: &ElectronState, b: Boundary, du: f64) {
        if let (Boundary::Surface(_), Some((e, mu))) = (b, self.at_face.take()) {
            assert_eq!(du, -self.inner);
            assert!((at.energy_ev - (e - self.inner)).abs() < 1e-12);
            if e * mu * mu <= self.inner {
                self.below_transmitted += 1;
            } else {
                self.above_transmitted += 1;
            }
        }
    }
    fn reflected(&mut self, at: &ElectronState, b: Boundary) {
        if let (Boundary::Surface(_), Some((e, mu))) = (b, self.at_face.take()) {
            assert_eq!(at.dir[0], -mu);
            if e * mu * mu <= self.inner {
                self.below_reflected += 1;
            }
        }
    }
    fn merge(&mut self, o: Self) {
        self.below_reflected += o.below_reflected;
        self.below_transmitted += o.below_transmitted;
        self.above_transmitted += o.above_transmitted;
    }
}

#[test]
fn an_electron_below_the_barrier_is_always_reflected() {
    // Isotropic elastic scattering, no loss: inside the electron has
    // E0 + U, and every direction is visited, so many arrive at a face with
    // E cos^2 below U.
    let mut cfg = TransportConfig::new(0.1);
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.escape_rule = EscapeRule::FrontOnly;
    let t = Transport::with_band_structures(
        Stack::new(vec![(material(), 4e-9)], None).unwrap(),
        vec![pair(elastic_isotropic(1e-9), inelastic_none())],
        vec![metal()],
        cfg,
    )
    .unwrap();
    let inner = metal().inner_potential_ev();
    let r = t
        .run(8, 2_000, 50, &Primary::normal(3.0), || FaceLog {
            inner,
            ..FaceLog::default()
        })
        .unwrap()
        .tally;
    assert!(r.below_reflected > 1_000, "{}", r.below_reflected);
    assert!(r.above_transmitted > 100, "{}", r.above_transmitted);
    assert_eq!(r.below_transmitted, 0);
}

#[test]
fn interface_step_changes_energy_by_the_inner_potential_difference() {
    #[derive(Default)]
    struct Steps {
        seen: u64,
    }
    impl ElectronTally for Steps {
        fn barrier(&mut self, at: &ElectronState, b: Boundary, du: f64) {
            if let Boundary::Interface { from, to } = b {
                assert_eq!(at.layer, to);
                // Layer 0 is the metal (U = 9), layer 1 the insulator (10).
                let expect = if (from, to) == (0, 1) { 1.0 } else { -1.0 };
                assert_eq!(du, expect);
                self.seen += 1;
            }
        }
        fn merge(&mut self, o: Self) {
            self.seen += o.seen;
        }
    }
    let mut cfg = TransportConfig::new(2.0);
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    let t = metal_on_insulator(cfg);
    let r = t
        .run(4, 300, 20, &Primary::normal(200.0), Steps::default)
        .unwrap()
        .tally;
    assert!(r.seen > 100);
}

// ---------------------------------------------------------------------------
// Cutoff reference
// ---------------------------------------------------------------------------

#[test]
fn vacuum_cutoff_stops_electrons_below_the_inner_potential_plus_cutoff() {
    #[derive(Default)]
    struct Stops {
        min_created: f64,
        max_stopped: f64,
        stopped: u64,
    }
    impl ElectronTally for Stops {
        fn secondary(
            &mut self,
            _p: &ElectronState,
            _e: &SecondaryEvent,
            c: Option<&ElectronState>,
        ) {
            if let Some(s) = c {
                self.min_created = self.min_created.min(s.energy_ev);
            }
        }
        fn stopped(&mut self, at: &ElectronState) {
            self.max_stopped = self.max_stopped.max(at.energy_ev);
            self.stopped += 1;
        }
        fn merge(&mut self, o: Self) {
            self.min_created = self.min_created.min(o.min_created);
            self.max_stopped = self.max_stopped.max(o.max_stopped);
            self.stopped += o.stopped;
        }
    }
    let cutoff = 1.5;
    let mut cfg = TransportConfig::new(cutoff);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    let t = Transport::with_band_structures(
        Stack::semi_infinite(material()),
        vec![pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3))],
        vec![metal()],
        cfg,
    )
    .unwrap();
    let threshold = metal().inner_potential_ev() + cutoff;
    let r = t
        .run(6, 200, 10, &Primary::normal(150.0), || Stops {
            min_created: f64::INFINITY,
            ..Stops::default()
        })
        .unwrap()
        .tally;
    assert!(r.stopped > 100);
    assert!(
        r.max_stopped < threshold,
        "{} >= {threshold}",
        r.max_stopped
    );
    assert!(
        r.min_created >= threshold,
        "{} < {threshold}",
        r.min_created
    );
}

// ---------------------------------------------------------------------------
// Defaults, validation and metadata
// ---------------------------------------------------------------------------

#[test]
fn band_parameters_change_nothing_while_the_models_are_off() {
    let cfg = TransportConfig::new(5.0);
    assert_eq!(cfg.secondaries, SecondaryModel::Off);
    assert_eq!(cfg.boundary, BoundaryModel::Transparent);
    assert_eq!(cfg.cutoff_reference, CutoffReference::BandBottom);
    let stack = || Stack::new(vec![(material(), 6e-9)], Some(material())).unwrap();
    let tables = || {
        vec![
            pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3)),
            pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.25)),
        ]
    };
    let plain = Transport::new(stack(), tables(), cfg).unwrap();
    let banded =
        Transport::with_band_structures(stack(), tables(), vec![metal(), insulator()], cfg)
            .unwrap();
    let p = Primary::normal(250.0);
    let a = plain.run(1, 300, 9, &p, SummaryTally::default).unwrap();
    let b = banded.run(1, 300, 9, &p, SummaryTally::default).unwrap();
    assert_eq!(a.tally, b.tally);
    assert_eq!(a.tally.secondaries + a.tally.reflections, 0);
    assert_eq!(a.tally.barrier_ev, 0.0);
}

#[test]
fn models_that_need_band_parameters_refuse_to_run_without_them() {
    let stack = || Stack::semi_infinite(material());
    let tables = || vec![pair(elastic_none(), inelastic_frac(1e-9, 0.1))];
    let mut se = TransportConfig::new(1.0);
    se.secondaries = SecondaryModel::KIEFT_BOSCH;
    assert!(Transport::new(stack(), tables(), se).is_err());
    let mut step = TransportConfig::new(1.0);
    step.boundary = BoundaryModel::STEP_BARRIER;
    assert!(Transport::new(stack(), tables(), step).is_err());
    assert!(Transport::with_band_structures(stack(), tables(), vec![], step).is_err());
    assert!(
        Transport::with_band_structures(stack(), tables(), vec![metal(), metal()], step).is_err()
    );
    // Under the secondary model the stopping threshold must clear the Fermi
    // energy (5 eV for the metal), measured either way.
    let mut low = se;
    low.cutoff_ev = 5.0;
    assert!(Transport::with_band_structures(stack(), tables(), vec![metal()], low).is_err());
    low.cutoff_ev = 5.5;
    assert!(Transport::with_band_structures(stack(), tables(), vec![metal()], low).is_ok());
    low.cutoff_ev = 0.1;
    low.cutoff_reference = CutoffReference::VacuumLevel;
    assert!(Transport::with_band_structures(stack(), tables(), vec![metal()], low).is_ok());
    let t = Transport::with_band_structures(stack(), tables(), vec![metal()], step).unwrap();
    // From vacuum the primary must head into the target.
    let outward = Primary {
        energy_ev: 100.0,
        direction: [-1.0, 0.0, 0.0],
    };
    assert!(t.run(1, 1, 1, &outward, SummaryTally::default).is_err());
    // With the vacuum cutoff the primary must clear U + cutoff once inside;
    // without the barrier its energy is already the inside one.
    let mut vac = TransportConfig::new(1.0);
    vac.cutoff_reference = CutoffReference::VacuumLevel;
    let t = Transport::with_band_structures(stack(), tables(), vec![metal()], vac).unwrap();
    assert!(t
        .run(1, 1, 1, &Primary::normal(9.5), SummaryTally::default)
        .is_err());
    assert!(t
        .run(1, 1, 1, &Primary::normal(10.5), SummaryTally::default)
        .is_ok());
}

#[test]
fn models_and_band_parameters_are_in_the_run_metadata() {
    let t = metal_on_insulator(full_physics());
    let m = t.metadata(1, 2, 3, &Primary::normal(100.0));
    assert_eq!(m.secondaries, SecondaryModel::KIEFT_BOSCH);
    assert_eq!(m.boundary, BoundaryModel::STEP_BARRIER);
    assert_eq!(m.layers[1].band_structure.as_ref(), Some(&insulator()));
    let json = serde_json::to_string(&m).unwrap();
    for needle in [
        "\"cutoff_reference\":\"vacuum-level\"",
        "\"model\":\"kieft-bosch\"",
        "\"instantaneous_momentum\":true",
        "\"model\":\"step-barrier\"",
        "\"kind\":\"insulator\"",
        "\"band_gap_ev\":3.0",
        PROVENANCE,
    ] {
        assert!(json.contains(needle), "{needle} missing from {json}");
    }
}

#[test]
fn free_electron_metal_uses_the_layer_density() {
    let m = material();
    let b = BandStructure::free_electron_metal(&m, 2.0, 4.0, PROVENANCE).unwrap();
    let n = 2.0 * m.atom_number_density();
    let ef = lindhard::electron::boundary::free_electron_fermi_energy_ev(n);
    assert_eq!(b.fermi_ev(), ef);
    assert_eq!(b.inner_potential_ev(), ef + 4.0);
    assert!(BandStructure::free_electron_metal(&m, 0.0, 4.0, PROVENANCE).is_err());
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// Every event with the exact bits of its numbers.
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
                s.pos[2].to_bits(),
                s.dir[1].to_bits(),
                s.energy_ev.to_bits(),
                extra.to_bits(),
            ],
        ));
    }
}

impl ElectronTally for EventLog {
    fn begin_history(&mut self, i: u64, s: &ElectronState) {
        self.index = i;
        self.push(9, s, 0.0);
    }
    fn step(&mut self, f: [f64; 3], e: &ElectronState, l: f64) {
        self.push(0, e, l);
        self.summary.step(f, e, l);
    }
    fn elastic(&mut self, a: &ElectronState, t: f64) {
        self.push(1, a, t);
        self.summary.elastic(a, t);
    }
    fn inelastic(&mut self, a: &ElectronState, w: f64) {
        self.push(2, a, w);
        self.summary.inelastic(a, w);
    }
    fn secondary(&mut self, p: &ElectronState, e: &SecondaryEvent, c: Option<&ElectronState>) {
        self.push(3, p, e.secondary_ev + 2.0 * e.deposited_ev);
        if let Some(s) = c {
            self.push(4, s, 0.0);
        }
        self.summary.secondary(p, e, c);
    }
    fn interface(&mut self, a: &ElectronState, f: usize, t: usize) {
        self.push(5, a, (f * 10 + t) as f64);
        self.summary.interface(a, f, t);
    }
    fn barrier(&mut self, a: &ElectronState, b: Boundary, du: f64) {
        self.push(6, a, du);
        self.summary.barrier(a, b, du);
    }
    fn reflected(&mut self, a: &ElectronState, b: Boundary) {
        self.push(7, a, 0.0);
        self.summary.reflected(a, b);
    }
    fn begin_secondary(&mut self, s: &ElectronState, g: u32) {
        self.push(8, s, g as f64);
    }
    fn end_secondary(&mut self, f: Fate) {
        self.summary.end_secondary(f);
    }
    fn stopped(&mut self, a: &ElectronState) {
        self.push(10, a, 0.0);
        self.summary.stopped(a);
    }
    fn escaped(&mut self, a: &ElectronState, face: Face) {
        self.push(11, a, if face == Face::Front { 0.0 } else { 1.0 });
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

fn run_threads(threads: usize) -> EventLog {
    let cfg = full_physics();
    let t = Transport::with_band_structures(
        Stack::new(
            vec![(material(), 5e-9), (material(), 8e-9)],
            Some(material()),
        )
        .unwrap(),
        vec![
            pair(elastic_isotropic(3e-9), inelastic_frac(5e-9, 0.25)),
            pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.3)),
            pair(elastic_isotropic(4e-9), inelastic_frac(6e-9, 0.2)),
        ],
        vec![metal(), insulator(), metal()],
        cfg,
    )
    .unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let primary = Primary {
        energy_ev: 400.0,
        direction: [1.0, 0.4, -0.1],
    };
    pool.install(|| {
        t.run(0x5EC0, 300, 7, &primary, EventLog::default)
            .unwrap()
            .tally
    })
}

#[test]
fn secondaries_and_barriers_are_bit_identical_across_thread_counts() {
    let r1 = run_threads(1);
    let s = &r1.summary;
    assert!(r1.events.len() > 50_000, "run too small to mean anything");
    assert!(s.secondaries > 1_000 && s.secondaries_escaped > 0);
    assert!(s.reflections > 0 && s.interface_crossings > 0 && s.escaped_front > 0);
    for n in [2, 8] {
        let r = run_threads(n);
        assert_eq!(r1.events, r.events, "event stream differs on {n} threads");
        let (a, b) = (&r1.summary, &r.summary);
        for (x, y) in [
            (a.path_m, b.path_m),
            (a.inelastic_loss_ev, b.inelastic_loss_ev),
            (a.escaped_energy_ev, b.escaped_energy_ev),
            (a.rest_energy_ev, b.rest_energy_ev),
            (a.secondary_energy_ev, b.secondary_energy_ev),
            (a.barrier_ev, b.barrier_ev),
        ] {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        assert_eq!(a, b);
    }
}
