//! Secondary electrons and surface barriers in the electron loop, on
//! synthetic cross-section tables and synthetic band parameters (made-up
//! numbers chosen to exercise the code, not material data): per-event energy
//! conservation, the secondary's direction model, the step transmission
//! against the analytic formula, reflection below the barrier, the vacuum
//! cutoff, run metadata, determinism across thread counts, and the energy
//! balance of `tally::electron` with both models on.

use lindhard::electron::boundary::{BandModel, BandStructure};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::secondary::{SecondaryEvent, SecondaryModel};
use lindhard::electron::transport::{
    Boundary, BoundaryModel, CutoffReference, ElectronState, ElectronTally, EscapeRule, Face, Fate,
    LayerTables, Primary, SummaryTally, Transport, TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{
    Binning, CartesianGrid, CylindricalGrid, ElectronReport, ElectronTallyConfig, FullElectronTally,
};

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
        let scale = ef.abs()
            + e.loss_ev.abs()
            + band.band_gap_ev().unwrap_or(0.0)
            + band.valence_binding_ev().unwrap_or(0.0);
        let residual = (e.loss_ev - (e.secondary_ev + e.binding_ev + e.deposited_ev)).abs();
        assert!(
            residual <= 4.0 * f64::EPSILON * scale,
            "residual {residual} for {e:?}"
        );
        self.worst = self.worst.max(residual / scale);
        // The binding of the liberated electron: an insulator's gap, or a
        // metal's valence binding (the Fermi level, 0, by default).
        let b = band
            .band_gap_ev()
            .or_else(|| band.valence_binding_ev())
            .unwrap_or(0.0);
        let below = match band.band_gap_ev() {
            Some(gap) => e.loss_ev <= gap,
            None => band.valence_binding_ev().is_some_and(|vb| e.loss_ev <= vb),
        };
        if below {
            self.subgap += 1;
            assert!(!e.liberated && c.is_none());
            assert_eq!((e.secondary_ev, e.binding_ev), (0.0, 0.0));
            assert_eq!(e.deposited_ev, e.loss_ev);
        } else {
            assert!(e.liberated, "{e:?}");
            self.liberated += 1;
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
// Energy reference of the inelastic table (#173)
// ---------------------------------------------------------------------------

/// Inelastic loss `W = frac * E` at every probability (a deterministic loss).
fn inelastic_fixed_frac(lambda_m: f64, frac: f64) -> CrossSectionTable {
    parts(
        SamplingAxis::InelasticEnergyLoss,
        1.0 / lambda_m,
        move |e, _p| frac * e,
    )
}

/// Every inelastic event as `(energy before, loss, energy after)`.
#[derive(Default)]
struct Losses {
    before: f64,
    events: Vec<(f64, f64, f64)>,
    secondaries: Vec<(f64, f64)>,
}

impl ElectronTally for Losses {
    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _l: f64) {
        self.before = end.energy_ev;
    }
    fn inelastic(&mut self, after: &ElectronState, w: f64) {
        self.events.push((self.before, w, after.energy_ev));
    }
    fn secondary(&mut self, _p: &ElectronState, e: &SecondaryEvent, c: Option<&ElectronState>) {
        if c.is_some() {
            self.secondaries.push((self.before, e.secondary_ev));
        }
    }
    fn merge(&mut self, mut o: Self) {
        self.events.append(&mut o.events);
        self.secondaries.append(&mut o.secondaries);
    }
}

fn losses_in_metal(frac: f64, cutoff: f64) -> Losses {
    let mut cfg = TransportConfig::new(cutoff);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    let t = Transport::with_band_structures(
        Stack::semi_infinite(material()),
        vec![pair(
            elastic_isotropic(3e-9),
            inelastic_fixed_frac(4e-9, frac),
        )],
        vec![metal()],
        cfg,
    )
    .unwrap();
    t.run(8, 50, 10, &Primary::normal(300.0), Losses::default)
        .unwrap()
        .tally
}

/// The transport reads the inelastic table at the electron's kinetic energy
/// measured from the band bottom (the vacuum energy plus `U` inside the
/// target), not from the Fermi level or the vacuum level: a table whose loss
/// is `E/2` at every row gives `W = E_before / 2` for every event. The
/// threshold, `U + 2 = 11` eV, keeps every event above `2 E_F = 10` eV, where
/// `E/2 < E - E_F` and the clamp below does not act.
#[test]
fn inelastic_table_is_read_at_the_band_bottom_energy() {
    let r = losses_in_metal(0.5, 2.0);
    assert!(r.events.len() > 500, "{}", r.events.len());
    // The primary's first event is at its entry energy, 300 eV + U.
    let u = metal().inner_potential_ev();
    assert!(r.events.iter().any(|&(e, _, _)| e == 300.0 + u));
    for &(e, w, after) in &r.events {
        assert!((w - 0.5 * e).abs() <= 1e-12 * e, "W {w} at E {e}");
        assert_eq!(after, e - w);
    }
}

/// With the Kieft-Bosch secondary model the transport clamps the sampled
/// loss to `E - E_F` (the primary cannot end below the Fermi level). A table
/// whose loss is the whole energy `E` at every row is therefore always
/// clamped: the primary is left at `E_F` and the secondary, `E_F + W`, gets
/// the primary's whole energy. This is the mismatch measured in #173: the
/// default tables are built with the model's Fermi energy 0, so their losses
/// reach `E`, beyond the `E - E_F` the transport allows (module docs of
/// `electron::transport`, "Energy reference of the inelastic table").
#[test]
fn losses_beyond_the_fermi_level_are_clamped_to_it() {
    let ef = metal().fermi_ev();
    let r = losses_in_metal(1.0, 2.0);
    assert!(r.events.len() > 100, "{}", r.events.len());
    for &(e, w, after) in &r.events {
        assert!((w - (e - ef)).abs() <= 1e-12 * e, "W {w} at E {e}");
        assert!((after - ef).abs() <= 1e-12 * e, "after {after}");
    }
    assert!(!r.secondaries.is_empty());
    for &(e, s) in &r.secondaries {
        assert!((s - e).abs() <= 1e-12 * e, "secondary {s} from E {e}");
    }
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
    run_threads_with(threads, metal())
}

/// [`run_threads`] with `metal` as the band of both metal layers.
fn run_threads_with(threads: usize, metal: BandStructure) -> EventLog {
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
        vec![metal.clone(), insulator(), metal],
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

// ---------------------------------------------------------------------------
// A metal's valence binding (BandStructure::with_valence_binding_ev)
// ---------------------------------------------------------------------------

/// The synthetic metal with its valence electron bound 7 eV below the Fermi
/// level (a made-up value).
fn bound_metal() -> BandStructure {
    metal().with_valence_binding_ev(7.0).unwrap()
}

#[test]
fn valence_binding_conserves_energy_and_frees_nothing_at_or_below_it() {
    let cfg = full_physics();
    let t = Transport::with_band_structures(
        Stack::new(vec![(material(), 6e-9)], Some(material())).unwrap(),
        vec![
            pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3)),
            pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.25)),
        ],
        vec![bound_metal(), bound_metal()],
        cfg,
    )
    .unwrap();
    let bands = vec![bound_metal(), bound_metal()];
    let r = t
        .run(11, 400, 16, &Primary::normal(400.0), || {
            EventCheck::new(bands.clone())
        })
        .unwrap()
        .tally;
    // EventCheck asserts, per event, E_SE = E_F + W - 7 eV for W > 7 eV and
    // nothing liberated for W <= 7 eV.
    assert!(r.events > 1_000, "too few events: {}", r.events);
    assert!(r.subgap > 100 && r.liberated > 100 && r.created > 100);
}

#[test]
fn zero_valence_binding_is_bit_identical_to_the_default() {
    let default = run_threads(1);
    let zero = run_threads_with(1, metal().with_valence_binding_ev(0.0).unwrap());
    assert_eq!(default.events, zero.events);
    assert_eq!(default.summary, zero.summary);
    // And a binding does change the run.
    let bound = run_threads_with(1, bound_metal());
    assert_ne!(default.events, bound.events);
}

#[test]
fn valence_binding_is_bit_identical_across_thread_counts() {
    let r1 = run_threads_with(1, bound_metal());
    assert!(r1.summary.secondaries > 100);
    for n in [2, 8] {
        let r = run_threads_with(n, bound_metal());
        assert_eq!(r1.events, r.events, "event stream differs on {n} threads");
        assert_eq!(r1.summary, r.summary);
    }
}

#[test]
fn valence_binding_is_refused_for_an_insulator_and_out_of_range() {
    assert!(insulator().with_valence_binding_ev(1.0).is_err());
    for b in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(metal().with_valence_binding_ev(b).is_err(), "{b}");
    }
    assert_eq!(metal().valence_binding_ev(), None);
    assert_eq!(bound_metal().valence_binding_ev(), Some(7.0));
}

// ---------------------------------------------------------------------------
// The electron tally with secondaries and barriers
// ---------------------------------------------------------------------------

const NM: f64 = 1e-9;

fn tally_config(e0: f64) -> ElectronTallyConfig {
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 60).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    c.cartesian = Some(CartesianGrid {
        x: Binning::new(0.0, 30.0 * NM, 6).unwrap(),
        y: Binning::new(-15.0 * NM, 15.0 * NM, 5).unwrap(),
        z: Binning::new(-15.0 * NM, 15.0 * NM, 5).unwrap(),
    });
    c.cylindrical = Some(CylindricalGrid {
        r: Binning::new(0.0, 15.0 * NM, 5).unwrap(),
        depth: Binning::new(0.0, 30.0 * NM, 6).unwrap(),
    });
    c
}

/// The run the review measured: 500 primaries at 300 eV, seed 3.
fn tally_run(t: &Transport) -> ElectronReport {
    let proto = FullElectronTally::new(t, tally_config(300.0)).unwrap();
    t.run(3, 500, 32, &Primary::normal(300.0), || proto.clone())
        .unwrap()
        .tally
        .report()
}

fn near(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
}

/// The balance `incident + fermi_sea = deposited + escaped + trapped +
/// barrier` and the sums every report must satisfy.
fn check_balance(r: &ElectronReport) {
    let b = &r.budget;
    assert!(
        b.relative_imbalance < 1e-9,
        "energy balance off by {} (relative): {b:?}",
        b.relative_imbalance
    );
    assert!(near(
        b.incident_ev + b.fermi_sea_ev,
        b.deposited_ev + b.escaped_ev + b.trapped_ev + b.barrier_ev,
        1e-9
    ));
    assert!(b.fermi_sea_ev >= 0.0);
    assert!(near(b.deposited_ev, b.inelastic_ev + b.residual_ev, 1e-12));
    let d = &r.deposition;
    assert!(near(d.per_layer_ev.iter().sum(), b.deposited_ev, 1e-9));
    for (cells, outside) in [
        d.cartesian.as_ref().map(|g| (&g.energy_ev, g.outside_ev)),
        d.cylindrical.as_ref().map(|g| (&g.energy_ev, g.outside_ev)),
    ]
    .into_iter()
    .flatten()
    {
        assert!(cells.iter().all(|&e| e >= 0.0));
        assert!(near(
            cells.iter().sum::<f64>() + outside,
            b.deposited_ev,
            1e-9
        ));
    }
    if let Some(g) = &r.generation_volume {
        assert!(near(g.energy_ev, b.deposited_ev, 1e-9));
    }
    assert!(near(r.front.energy_ev, b.escaped_front_ev, 1e-12));
    assert!(near(r.back.energy_ev, b.escaped_back_ev, 1e-12));
    let f = &r.fates;
    assert_eq!(
        f.stopped + f.escaped_front + f.escaped_back + f.absorbed + f.trapped + f.event_capped,
        r.histories
    );
}

#[test]
fn full_tally_balances_with_secondaries() {
    // Band-bottom cutoff 8 eV, above both Fermi energies (5 and 7.5 eV).
    let mut cfg = TransportConfig::new(8.0);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    let r = tally_run(&metal_on_insulator(cfg));
    check_balance(&r);
    let b = &r.budget;
    // Both layers liberate conduction-band electrons (B - E_F < 0).
    assert!(b.fermi_sea_ev > 0.0);
    assert_eq!(b.barrier_ev, 0.0);
    assert_eq!(b.trapped_ev, 0.0);
    // Secondaries end in the tally: more stops than primaries, and slow
    // electrons out of the front face.
    assert!(r.stopping_points.stopped > r.fates.stopped);
    // The primary-only stopping points are exactly the stopped histories.
    assert_eq!(r.stopping_points.primaries.stopped, r.fates.stopped);
    assert!(r.front.slow.count > 0);
    assert_eq!(r.metadata.stopping_threshold_ev, vec![8.0, 8.0]);
}

#[test]
fn full_tally_balances_with_the_barrier_from_the_band_bottom() {
    let mut cfg = TransportConfig::new(8.0);
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    let r = tally_run(&metal_on_insulator(cfg));
    check_balance(&r);
    let b = &r.budget;
    assert_eq!(b.fermi_sea_ev, 0.0);
    // Every primary that got in paid -U = -9 eV, every escape got +U back.
    assert!(b.barrier_ev != 0.0);
    assert_eq!(b.trapped_ev, 0.0);
    assert_eq!(r.stopping_points.stopped, r.fates.stopped);
}

#[test]
fn full_tally_balances_with_the_barrier_from_the_vacuum_level() {
    let mut cfg = TransportConfig::new(0.5);
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    let r = tally_run(&metal_on_insulator(cfg));
    check_balance(&r);
    let b = &r.budget;
    assert!(b.barrier_ev != 0.0);
    // Stops between the cutoff and U + cutoff are stops, not trapped
    // electrons: no layer here has a zero rate.
    assert_eq!(b.trapped_ev, 0.0);
    assert!(r.fates.stopped > 0 && b.residual_ev > 0.0);
    assert_eq!(r.stopping_points.stopped, r.fates.stopped);
    assert_eq!(r.metadata.stopping_threshold_ev, vec![9.5, 10.5]);
}

#[test]
fn full_tally_balances_with_secondaries_and_barriers() {
    let r = tally_run(&metal_on_insulator(full_physics()));
    check_balance(&r);
    let b = &r.budget;
    assert!(b.fermi_sea_ev > 0.0 && b.barrier_ev != 0.0);
    assert_eq!(b.trapped_ev, 0.0);
    assert!(r.stopping_points.stopped > r.fates.stopped);
    assert!(r.front.slow.count > 0 && r.front.fast.count > 0);
}

#[test]
fn full_tally_balances_when_primaries_and_secondaries_hit_the_event_cap() {
    let mut cfg = full_physics();
    cfg.max_events = 6;
    cfg.escape_rule = EscapeRule::FrontOnly;
    let t = Transport::with_band_structures(
        Stack::new(vec![(material(), 6e-9), (material(), 4e-9)], None).unwrap(),
        vec![
            pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3)),
            pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.25)),
        ],
        vec![metal(), insulator()],
        cfg,
    )
    .unwrap();
    let r = tally_run(&t);
    check_balance(&r);
    assert!(r.fates.event_capped > 0 && r.budget.event_cap_ev > 0.0);
    assert!(r.budget.absorbed_ev > 0.0);
}

/// A transport whose small event cap cuts off primaries and secondaries.
fn capped_transport() -> Transport {
    let mut cfg = full_physics();
    cfg.max_events = 6;
    Transport::with_band_structures(
        Stack::new(vec![(material(), 6e-9), (material(), 4e-9)], None).unwrap(),
        vec![
            pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3)),
            pair(elastic_isotropic(2e-9), inelastic_frac(3e-9, 0.25)),
        ],
        vec![metal(), insulator()],
        cfg,
    )
    .unwrap()
}

fn capped_run<T: ElectronTally + Clone + Sync>(threads: usize, proto: T) -> T {
    let t = capped_transport();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| {
            t.run(0xCA9, 400, 16, &Primary::normal(300.0), || proto.clone())
                .unwrap()
                .tally
        })
}

#[test]
fn event_cap_counts_merge_identically_across_threads_and_match_the_summary() {
    let t = capped_transport();
    let proto = FullElectronTally::new(&t, tally_config(300.0)).unwrap();
    let r1 = capped_run(1, proto.clone()).report();
    check_balance(&r1);
    let c = r1.event_caps;
    assert!(c.secondary_tracks > 0, "{c:?}");
    assert!(r1.fates.event_capped > 0);
    // Every history with a capped primary is affected; no history is counted
    // twice.
    assert!(c.affected_histories >= r1.fates.event_capped);
    assert!(c.affected_histories <= r1.histories);
    assert!(c.affected_histories <= r1.fates.event_capped + c.secondary_tracks);
    assert!(r1.budget.event_cap_ev > 0.0);
    let r2 = capped_run(2, proto).report();
    assert_eq!(r1, r2, "report differs on 2 threads");
    // The lightweight tally sees the same secondary caps.
    let s1 = capped_run(1, SummaryTally::default());
    let s2 = capped_run(2, SummaryTally::default());
    assert_eq!(s1, s2);
    assert_eq!(s1.secondaries_event_capped, c.secondary_tracks);
    assert_eq!(s1.event_capped, r1.fates.event_capped);
}

#[test]
fn full_tally_balances_with_a_binding_below_the_band_bottom() {
    // Synthetic insulator with W_v = 1, E_g = 4, χ = 1 eV: E_F = 3 eV, so
    // B - E_F = +1 eV and the liberated electron brings nothing; that 1 eV of
    // each event stays in the solid.
    let deep = BandStructure::new(
        BandModel::Insulator {
            valence_band_width_ev: 1.0,
            band_gap_ev: 4.0,
            affinity_ev: 1.0,
        },
        PROVENANCE,
    )
    .unwrap();
    let mut cfg = TransportConfig::new(4.0);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    let t = Transport::with_band_structures(
        Stack::semi_infinite(material()),
        vec![pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3))],
        vec![deep],
        cfg,
    )
    .unwrap();
    let r = tally_run(&t);
    check_balance(&r);
    assert_eq!(r.budget.fermi_sea_ev, 0.0);
    assert!(r.stopping_points.stopped > r.fates.stopped);
}

#[test]
fn full_tally_balances_with_a_metal_valence_binding() {
    // The synthetic metal with B = 7 eV below its 5 eV Fermi level: each
    // liberated electron's initial state is 2 eV below the band bottom, and
    // a loss W <= 7 eV stays in the solid.
    let mut cfg = TransportConfig::new(8.0);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    let t = Transport::with_band_structures(
        Stack::semi_infinite(material()),
        vec![pair(elastic_isotropic(3e-9), inelastic_frac(4e-9, 0.3))],
        vec![bound_metal()],
        cfg,
    )
    .unwrap();
    let r = tally_run(&t);
    check_balance(&r);
    assert_eq!(r.budget.fermi_sea_ev, 0.0);
    assert!(r.stopping_points.stopped > r.fates.stopped);
}

fn at(x_nm: f64, energy_ev: f64, layer: usize) -> ElectronState {
    ElectronState {
        pos: [x_nm * NM, 0.0, 0.0],
        dir: [1.0, 0.0, 0.0],
        energy_ev,
        layer,
    }
}

/// A conduction electron liberated in the metal (`E_F = 5`, `B = 0`).
fn metal_event(w: f64) -> SecondaryEvent {
    SecondaryEvent {
        loss_ev: w,
        liberated: true,
        secondary_ev: 5.0 + w,
        binding_ev: -5.0,
        deposited_ev: 0.0,
    }
}

#[test]
fn full_tally_keeps_a_capped_primary_apart_from_its_secondaries() {
    let t = metal_on_insulator(full_physics());
    let mut tally = FullElectronTally::new(&t, tally_config(300.0)).unwrap();
    // The primary enters (+9 eV), loses 100 eV to a 105 eV secondary and is
    // cut off by the event cap at 209 eV; the secondary then escapes.
    let mut vac = at(0.0, 300.0, 0);
    tally.begin_history(0, &vac);
    vac.energy_ev = 309.0;
    tally.barrier(&vac, Boundary::Surface(Face::Front), 9.0);
    tally.step([0.0; 3], &at(1.0, 309.0, 0), 1.0 * NM);
    tally.inelastic(&at(1.0, 209.0, 0), 100.0);
    let s = at(1.0, 105.0, 0);
    tally.secondary(&at(1.0, 209.0, 0), &metal_event(100.0), Some(&s));
    tally.begin_secondary(&s, 1);
    let mut out = at(0.0, 105.0, 0);
    out.dir = [-1.0, 0.0, 0.0];
    tally.step(s.pos, &out, 1.0 * NM);
    out.energy_ev = 96.0;
    tally.barrier(&out, Boundary::Surface(Face::Front), -9.0);
    tally.escaped(&out, Face::Front);
    tally.end_secondary(Fate::Escaped(Face::Front));
    tally.end_history(0, Fate::EventCap);
    let r = tally.report();
    check_balance(&r);
    assert_eq!(r.budget.event_cap_ev, 209.0);
    assert_eq!(r.budget.escaped_front_ev, 96.0);
    assert_eq!(r.budget.fermi_sea_ev, 5.0);
    assert_eq!(r.budget.barrier_ev, 0.0);
    assert_eq!(r.budget.deposited_ev, 0.0);
    assert_eq!(r.fates.event_capped, 1);
    // The primary was capped, its secondary escaped.
    assert_eq!(r.event_caps.secondary_tracks, 0);
    assert_eq!(r.event_caps.affected_histories, 1);
}

#[test]
fn full_tally_counts_capped_secondaries_and_vacuum_level_stops() {
    let t = metal_on_insulator(full_physics());
    let mut tally = FullElectronTally::new(&t, tally_config(300.0)).unwrap();
    // The primary enters, loses 300 eV and stops at 9 eV, below the metal's
    // threshold U + cutoff = 9.5 eV (but above the 0.5 eV cutoff itself);
    // its 305 eV secondary is cut off by the event cap.
    let mut vac = at(0.0, 300.0, 0);
    tally.begin_history(0, &vac);
    vac.energy_ev = 309.0;
    tally.barrier(&vac, Boundary::Surface(Face::Front), 9.0);
    tally.step([0.0; 3], &at(2.0, 309.0, 0), 2.0 * NM);
    tally.inelastic(&at(2.0, 9.0, 0), 300.0);
    let s = at(2.0, 305.0, 0);
    tally.secondary(&at(2.0, 9.0, 0), &metal_event(300.0), Some(&s));
    tally.stopped(&at(2.0, 9.0, 0));
    tally.begin_secondary(&s, 1);
    tally.step(s.pos, &at(3.0, 305.0, 0), 1.0 * NM);
    tally.end_secondary(Fate::EventCap);
    tally.end_history(0, Fate::Stopped);
    let r = tally.report();
    check_balance(&r);
    assert_eq!(r.budget.residual_ev, 9.0);
    assert_eq!(r.budget.no_interaction_ev, 0.0);
    assert_eq!(r.budget.event_cap_ev, 305.0);
    assert_eq!(r.budget.barrier_ev, -9.0);
    assert_eq!(r.stopping_points.stopped, 1);
    assert_eq!(r.fates.stopped, 1);
    // A completed primary with a capped secondary: no primary cap, one
    // capped track, one affected history.
    assert_eq!(r.fates.event_capped, 0);
    assert_eq!(r.event_caps.secondary_tracks, 1);
    assert_eq!(r.event_caps.affected_histories, 1);
}

/// A secondary of `energy_ev` liberated at `x_nm` that is cut off by the
/// event cap after one step.
fn capped_secondary(tally: &mut FullElectronTally, x_nm: f64, energy_ev: f64, generation: u32) {
    let s = at(x_nm, energy_ev, 0);
    tally.begin_secondary(&s, generation);
    tally.step(s.pos, &at(x_nm + 0.5, energy_ev, 0), 0.5 * NM);
    tally.end_secondary(Fate::EventCap);
}

#[test]
fn full_tally_counts_each_capped_track_but_each_history_once() {
    let t = metal_on_insulator(full_physics());
    let mut tally = FullElectronTally::new(&t, tally_config(300.0)).unwrap();

    // History 0: the primary stops at 9 eV after three
    // 100 eV losses (all to metal conduction electrons); two of its three
    // secondaries are cut off by the cap.
    let mut vac = at(0.0, 300.0, 0);
    tally.begin_history(0, &vac);
    vac.energy_ev = 309.0;
    tally.barrier(&vac, Boundary::Surface(Face::Front), 9.0);
    tally.step([0.0; 3], &at(1.0, 309.0, 0), 1.0 * NM);
    let mut secs = Vec::new();
    for (e_after, w) in [(209.0, 100.0), (109.0, 100.0), (9.0, 100.0)] {
        tally.inelastic(&at(1.0, e_after, 0), w);
        let s = at(1.0, 5.0 + w, 0);
        tally.secondary(&at(1.0, e_after, 0), &metal_event(w), Some(&s));
        secs.push(s);
    }
    tally.stopped(&at(1.0, 9.0, 0));
    capped_secondary(&mut tally, 1.0, secs[0].energy_ev, 1);
    capped_secondary(&mut tally, 1.0, secs[1].energy_ev, 1);
    // The third finds no interaction and is trapped with its energy.
    tally.begin_secondary(&secs[2], 1);
    tally.stopped(&secs[2]);
    tally.end_secondary(Fate::Trapped);
    tally.end_history(0, Fate::Stopped);
    let r = tally.report();
    assert_eq!(r.fates.event_capped, 0);
    assert_eq!(r.event_caps.secondary_tracks, 2);
    assert_eq!(r.event_caps.affected_histories, 1);
    assert_eq!(r.budget.event_cap_ev, 210.0);

    // History 1: the primary is capped at 209 eV and both of its
    // secondaries are capped too: still one affected history.
    let mut vac = at(0.0, 300.0, 0);
    tally.begin_history(1, &vac);
    vac.energy_ev = 309.0;
    tally.barrier(&vac, Boundary::Surface(Face::Front), 9.0);
    tally.step([0.0; 3], &at(1.0, 309.0, 0), 1.0 * NM);
    tally.inelastic(&at(1.0, 259.0, 0), 50.0);
    tally.secondary(
        &at(1.0, 259.0, 0),
        &metal_event(50.0),
        Some(&at(1.0, 55.0, 0)),
    );
    tally.inelastic(&at(1.0, 209.0, 0), 50.0);
    tally.secondary(
        &at(1.0, 209.0, 0),
        &metal_event(50.0),
        Some(&at(1.0, 55.0, 0)),
    );
    capped_secondary(&mut tally, 1.0, 55.0, 1);
    capped_secondary(&mut tally, 1.0, 55.0, 1);
    tally.end_history(1, Fate::EventCap);
    let r = tally.report();
    assert_eq!(r.fates.event_capped, 1);
    assert_eq!(r.event_caps.secondary_tracks, 4);
    assert_eq!(r.event_caps.affected_histories, 2);
    // The primary's own last state (209 eV) is still the one counted for it.
    assert_eq!(r.budget.event_cap_ev, 210.0 + 110.0 + 209.0);

    // History 2: the primary escapes without secondaries: not affected.
    let vac = at(0.0, 300.0, 0);
    tally.begin_history(2, &vac);
    let mut out = at(0.0, 300.0, 0);
    out.dir = [-1.0, 0.0, 0.0];
    tally.reflected(&out, Boundary::Surface(Face::Front));
    tally.escaped(&out, Face::Front);
    tally.end_history(2, Fate::Escaped(Face::Front));
    let r = tally.report();
    check_balance(&r);
    assert_eq!(r.histories, 3);
    assert_eq!(r.fates.event_capped, 1);
    assert_eq!(r.event_caps.secondary_tracks, 4);
    assert_eq!(r.event_caps.affected_histories, 2);
}

fn tally_threaded(threads: usize) -> ElectronReport {
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
        full_physics(),
    )
    .unwrap();
    let proto = FullElectronTally::new(&t, tally_config(400.0)).unwrap();
    let primary = Primary {
        energy_ev: 400.0,
        direction: [1.0, 0.4, -0.1],
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| {
            t.run(0x5EC0, 600, 7, &primary, || proto.clone())
                .unwrap()
                .tally
                .report()
        })
}

#[test]
fn full_tally_report_is_bit_identical_for_1_2_and_8_threads_with_both_models() {
    let r1 = tally_threaded(1);
    check_balance(&r1);
    assert!(r1.budget.fermi_sea_ev > 0.0 && r1.budget.barrier_ev != 0.0);
    assert!(r1.front.slow.count > 100);
    let j1 = serde_json::to_string(&r1).unwrap();
    for n in [2, 8] {
        let r = tally_threaded(n);
        assert_eq!(r1, r, "report differs on {n} threads");
        assert_eq!(j1, serde_json::to_string(&r).unwrap());
    }
}

/// Deposits are measured from the band bottom, so with secondaries in a metal
/// they include the Fermi-sea energy of liberated conduction electrons and can
/// exceed the energy the beam imparted. The net energy imparted is
/// `incident - escaped = deposited + trapped + barrier - fermi_sea -
/// phonon_absorbed` (module docs of `lindhard::tally::electron`, "Deposits and
/// the energy imparted"). The run is the one the review of #116 measured.
#[test]
fn net_energy_imparted_subtracts_the_fermi_sea_from_the_deposits() {
    let r = tally_run(&metal_on_insulator(full_physics()));
    check_balance(&r);
    let b = &r.budget;
    let imparted = b.incident_ev - b.escaped_ev;
    let net = b.deposited_ev + b.trapped_ev + b.barrier_ev - b.fermi_sea_ev - b.phonon_absorbed_ev;
    // Same tolerance as the balance, relative to the energy entering it.
    let scale = b.incident_ev + b.fermi_sea_ev + b.phonon_absorbed_ev;
    assert!(
        (imparted - net).abs() <= 1e-9 * scale,
        "imparted {imparted} eV, net from deposits {net} eV: {b:?}"
    );
    // The documented caveat, pinned by a real case: the deposits exceed what
    // the beam gave the solid.
    assert!(b.fermi_sea_ev > 0.0);
    assert!(
        b.deposited_ev > imparted,
        "deposited {} eV, imparted {imparted} eV",
        b.deposited_ev
    );
}
