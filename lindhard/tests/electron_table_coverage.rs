//! Table-coverage diagnostics of the electron transport (#253): how many
//! rate evaluations read each cross-section table inside its energy grid, and
//! how many used the constant continuation below or above it. Exact counts on
//! histories whose energies are fixed by construction (lower and upper
//! overflow, exact endpoints, wholly in-grid transport, elastic and inelastic
//! grids that differ, a two-layer crossing, secondaries), the identity with
//! the number of flights in a random run, determinism across thread counts,
//! the claim that counting changes nothing, and reading reports written before
//! the field existed. The tables and band parameters are made up for the
//! tests (constant or proportional losses, isotropic angles), not material
//! data.

use lindhard::electron::boundary::{BandModel, BandStructure};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::secondary::SecondaryModel;
use lindhard::electron::transport::{
    BoundaryModel, CutoffReference, ElectronState, ElectronTally, Fate, GridCoverage, LayerTables,
    Primary, SummaryTally, TableChannel, Transport, TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{
    Binning, ElectronReport, ElectronTallyConfig, FullElectronTally, LayerTableCoverage,
    TableCoverageCounts, TableCoverageTally,
};

const PROB: [f64; 3] = [0.0, 0.5, 1.0];
const NM: f64 = 1e-9;
const PROVENANCE: &str = "synthetic, tests/electron_table_coverage.rs";

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

/// A table on `grid` with the same rate at every energy and quantile rows
/// from `row(E, p)` (none where the rate is zero).
fn table(
    axis: SamplingAxis,
    grid: &[f64],
    rate: f64,
    row: impl Fn(f64, f64) -> f64,
) -> CrossSectionTable {
    let quantiles = grid
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
        energy_ev: grid.to_vec(),
        inverse_mfp_per_m: vec![rate; grid.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

fn elastic_none(grid: &[f64]) -> CrossSectionTable {
    table(SamplingAxis::ElasticPolarAngle, grid, 0.0, |_, _| 0.0)
}

fn inelastic_none(grid: &[f64]) -> CrossSectionTable {
    table(SamplingAxis::InelasticEnergyLoss, grid, 0.0, |_, _| 0.0)
}

fn elastic_isotropic(grid: &[f64], lambda_m: f64) -> CrossSectionTable {
    table(
        SamplingAxis::ElasticPolarAngle,
        grid,
        1.0 / lambda_m,
        |_, p| (1.0 - 2.0 * p).acos(),
    )
}

/// Loss `dw` at every energy and probability: the same inside and beyond the
/// grid, so the continuation does not change the history.
fn inelastic_const(grid: &[f64], lambda_m: f64, dw: f64) -> CrossSectionTable {
    table(
        SamplingAxis::InelasticEnergyLoss,
        grid,
        1.0 / lambda_m,
        move |_, _| dw,
    )
}

/// Loss `frac * E` at every probability.
fn inelastic_frac(grid: &[f64], lambda_m: f64, frac: f64) -> CrossSectionTable {
    table(
        SamplingAxis::InelasticEnergyLoss,
        grid,
        1.0 / lambda_m,
        move |e, _| frac * e,
    )
}

fn pair(elastic: CrossSectionTable, inelastic: CrossSectionTable) -> LayerTables {
    LayerTables { elastic, inelastic }
}

fn finite(thicknesses_nm: &[f64]) -> Stack {
    Stack::new(
        thicknesses_nm
            .iter()
            .map(|&t| (material(), t * NM))
            .collect(),
        None,
    )
    .unwrap()
}

/// Synthetic metal: `E_F = 5 eV`, `Φ = 4 eV`.
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

fn counts(min: f64, max: f64, below: u64, within: u64, above: u64) -> TableCoverageCounts {
    TableCoverageCounts {
        energy_min_ev: min,
        energy_max_ev: max,
        below,
        within,
        above,
    }
}

fn coverage_of(t: &Transport, n: u64, primary: &Primary) -> Vec<LayerTableCoverage> {
    let proto = TableCoverageTally::new(t);
    t.run(7, n, 4, primary, || proto.clone())
        .unwrap()
        .tally
        .layers()
        .to_vec()
}

// ---------------------------------------------------------------------------
// Exact counts
// ---------------------------------------------------------------------------

#[test]
fn grid_coverage_includes_both_endpoints() {
    let t = elastic_none(&[10.0, 100.0, 1_000.0]);
    assert_eq!(GridCoverage::of(&t, 9.999), GridCoverage::Below);
    assert_eq!(GridCoverage::of(&t, 10.0), GridCoverage::Within);
    assert_eq!(GridCoverage::of(&t, 55.0), GridCoverage::Within);
    assert_eq!(GridCoverage::of(&t, 1_000.0), GridCoverage::Within);
    assert_eq!(GridCoverage::of(&t, 1_000.001), GridCoverage::Above);
}

/// Zero rates in a finite layer: one pass of the loop, which evaluates both
/// (zero) rates once and flies straight to the back face. So every history is
/// exactly one evaluation per channel, a boundary-limited flight with no
/// collision. The primary's energy against the grid `[10, 1000]` eV decides
/// the class.
#[test]
fn zero_rate_flight_is_one_evaluation_classed_by_energy() {
    let grid = [10.0, 100.0, 1_000.0];
    let t = Transport::new(
        finite(&[5.0]),
        vec![pair(elastic_none(&grid), inelastic_none(&grid))],
        TransportConfig::new(1.0),
    )
    .unwrap();
    let n = 9;
    for (e0, below, within, above) in [
        (5.0, n, 0, 0),
        (10.0, 0, n, 0),
        (55.0, 0, n, 0),
        (1_000.0, 0, n, 0),
        (2_000.0, 0, 0, n),
    ] {
        let c = coverage_of(&t, n, &Primary::normal(e0));
        let expect = counts(10.0, 1_000.0, below, within, above);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].layer, 0);
        assert_eq!(c[0].elastic, expect, "elastic at {e0} eV");
        assert_eq!(c[0].inelastic, expect, "inelastic at {e0} eV");
    }
}

/// A semi-infinite substrate, no elastic scattering, a constant loss of
/// 50 eV, no secondaries: the primary moves straight in and is evaluated at
/// exactly 500, 450, ..., 50 eV (ten collisions; after the last it has 0 eV,
/// below the 25 eV cutoff). The inelastic grid `[100, 400]` sees 50 below,
/// 100 to 400 within (both endpoints hit exactly) and 450, 500 above; the
/// elastic grid `[50, 500]` sees all ten within, at both endpoints.
#[test]
fn collisions_across_the_grid_ends_with_exact_endpoints() {
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![pair(
            elastic_none(&[50.0, 200.0, 500.0]),
            inelastic_const(&[100.0, 250.0, 400.0], 2.0 * NM, 50.0),
        )],
        TransportConfig::new(25.0),
    )
    .unwrap();
    let n = 6;
    let p = t
        .run(3, n, 4, &Primary::normal(500.0), || Probe::new(&t))
        .unwrap()
        .tally;
    assert_eq!(p.inelastic, 10 * n);
    assert_eq!(p.stopped, n);
    let c = p.cov.layers();
    assert_eq!(c[0].elastic, counts(50.0, 500.0, 0, 10 * n, 0));
    assert_eq!(c[0].inelastic, counts(100.0, 400.0, n, 7 * n, 2 * n));
}

/// Two finite layers with zero rates and different grids per layer and per
/// channel: the primary is evaluated once in layer 0, crosses the interface,
/// is evaluated once in layer 1 and leaves through the back face. At 300 eV:
/// layer 0 elastic `[10, 1000]` within, inelastic `[500, 5000]` below; layer
/// 1 elastic `[20, 200]` above, inelastic `[300, 400]` within (exact lower
/// endpoint).
#[test]
fn multilayer_crossing_counts_each_layer_and_channel_on_its_own_grid() {
    let t = Transport::new(
        finite(&[3.0, 4.0]),
        vec![
            pair(
                elastic_none(&[10.0, 1_000.0]),
                inelastic_none(&[500.0, 5_000.0]),
            ),
            pair(
                elastic_none(&[20.0, 200.0]),
                inelastic_none(&[300.0, 400.0]),
            ),
        ],
        TransportConfig::new(1.0),
    )
    .unwrap();
    let n = 5;
    let c = coverage_of(&t, n, &Primary::normal(300.0));
    assert_eq!(c.len(), 2);
    assert_eq!(c[0].elastic, counts(10.0, 1_000.0, 0, n, 0));
    assert_eq!(c[0].inelastic, counts(500.0, 5_000.0, n, 0, 0));
    assert_eq!(c[1].layer, 1);
    assert_eq!(c[1].elastic, counts(20.0, 200.0, 0, 0, n));
    assert_eq!(c[1].inelastic, counts(300.0, 400.0, 0, n, 0));
}

/// Which electron each evaluation belongs to.
#[derive(Default)]
struct BySource {
    in_secondary: bool,
    primary: Vec<(TableChannel, GridCoverage)>,
    secondary: Vec<(TableChannel, GridCoverage)>,
}

impl ElectronTally for BySource {
    fn begin_history(&mut self, _i: u64, _s: &ElectronState) {
        self.in_secondary = false;
    }
    fn begin_secondary(&mut self, _s: &ElectronState, _g: u32) {
        self.in_secondary = true;
    }
    fn table_lookup(&mut self, _at: &ElectronState, ch: TableChannel, c: GridCoverage) {
        if self.in_secondary {
            self.secondary.push((ch, c));
        } else {
            self.primary.push((ch, c));
        }
    }
    fn merge(&mut self, mut o: Self) {
        self.primary.append(&mut o.primary);
        self.secondary.append(&mut o.secondary);
    }
}

fn tally_config(e0: f64) -> ElectronTallyConfig {
    ElectronTallyConfig::new(
        Binning::new(0.0, e0, 10).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    )
}

/// Secondaries on a fixed energy ladder. Metal with `E_F = 5 eV`, cutoff
/// 45 eV from the band bottom, loss `W = E/2` (the inelastic grid covers
/// every energy, so this holds), no elastic scattering, no deflection of the
/// primary and no instantaneous momentum, in a semi-infinite substrate. A
/// secondary gets `E_F + W` and starts at a positive angle to `+x`, so the
/// primary and its secondaries move deeper and never escape. Energies at
/// each evaluation (up to rounding of the interpolated loss):
///
/// - primary: 200, 100, 50 (then 25, stopped; its secondaries 105, 55, 30,
///   the last below the cutoff and not created);
/// - secondary 55: 55 (then 27.5, stopped; its secondary 32.5 not created);
/// - secondary 105: 105, 52.5 (then 26.25; secondaries 57.5, created, and
///   31.25, not);
/// - secondary 57.5: 57.5, once whether it collides (then 28.75, stopped,
///   and 33.75 not created) or leaves through the front face first.
///
/// So 3 primary and 4 secondary evaluations per history. Against the
/// elastic grid `[60, 150]`: 50, 55, 52.5 and 57.5 below, 100 and 105
/// within, 200 above.
#[test]
fn secondaries_are_counted_with_the_primary() {
    let mut cfg = TransportConfig::new(45.0);
    cfg.secondaries = SecondaryModel::KieftBosch {
        instantaneous_momentum: false,
        momentum_conservation: false,
    };
    let t = Transport::with_band_structures(
        Stack::semi_infinite(material()),
        vec![pair(
            elastic_none(&[60.0, 150.0]),
            inelastic_frac(&[1.0, 10.0, 100.0, 1_000.0], 3.0 * NM, 0.5),
        )],
        vec![metal()],
        cfg,
    )
    .unwrap();
    let n = 8;
    let p = Primary::normal(200.0);

    let by = t.run(5, n, 3, &p, BySource::default).unwrap().tally;
    let class = |v: &[(TableChannel, GridCoverage)], ch: TableChannel, c: GridCoverage| {
        v.iter().filter(|&&x| x == (ch, c)).count() as u64
    };
    use GridCoverage::{Above, Below, Within};
    use TableChannel::{Elastic, Inelastic};
    assert_eq!(by.primary.len() as u64, 2 * 3 * n);
    assert_eq!(by.secondary.len() as u64, 2 * 4 * n);
    assert_eq!(class(&by.primary, Elastic, Below), n);
    assert_eq!(class(&by.primary, Elastic, Within), n);
    assert_eq!(class(&by.primary, Elastic, Above), n);
    assert_eq!(class(&by.primary, Inelastic, Within), 3 * n);
    assert_eq!(class(&by.secondary, Elastic, Below), 3 * n);
    assert_eq!(class(&by.secondary, Elastic, Within), n);
    assert_eq!(class(&by.secondary, Inelastic, Within), 4 * n);

    // The full tally reports the same counts, primaries and secondaries
    // together.
    let proto = FullElectronTally::new(&t, tally_config(200.0)).unwrap();
    let r = t.run(5, n, 3, &p, || proto.clone()).unwrap().tally.report();
    assert_eq!(r.table_coverage.len(), 1);
    assert_eq!(
        r.table_coverage[0].elastic,
        counts(60.0, 150.0, 4 * n, 2 * n, n)
    );
    assert_eq!(
        r.table_coverage[0].inelastic,
        counts(1.0, 1_000.0, 0, 7 * n, 0)
    );
}

// ---------------------------------------------------------------------------
// A random run: identity with the flights, determinism, no side effects
// ---------------------------------------------------------------------------

/// Two finite layers, elastic and inelastic scattering, secondaries and the
/// step barrier, with grids narrower than the energies reached, so that all
/// three classes occur.
fn mixed() -> Transport {
    let mut cfg = TransportConfig::new(2.0);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    Transport::with_band_structures(
        finite(&[4.0, 20.0]),
        vec![
            pair(
                elastic_isotropic(&[30.0, 100.0, 300.0], 2.0 * NM),
                inelastic_frac(&[20.0, 200.0], 4.0 * NM, 0.3),
            ),
            pair(
                elastic_isotropic(&[15.0, 150.0], 3.0 * NM),
                inelastic_frac(&[40.0, 400.0], 3.0 * NM, 0.3),
            ),
        ],
        vec![metal(), metal()],
        cfg,
    )
    .unwrap()
}

/// Flights, collisions and evaluations, with the coverage counts.
#[derive(Clone)]
struct Probe {
    /// Flights started, per layer (each ends in a `step`).
    flights: Vec<u64>,
    /// Evaluations per layer, `[elastic, inelastic]`.
    lookups: Vec<[u64; 2]>,
    collisions: u64,
    inelastic: u64,
    stopped: u64,
    secondaries: u64,
    reflections: u64,
    interfaces: u64,
    trapped: u64,
    cov: TableCoverageTally,
}

impl Probe {
    fn new(t: &Transport) -> Self {
        let n = t.stack().layers().len();
        Self {
            flights: vec![0; n],
            lookups: vec![[0; 2]; n],
            collisions: 0,
            inelastic: 0,
            stopped: 0,
            secondaries: 0,
            reflections: 0,
            interfaces: 0,
            trapped: 0,
            cov: TableCoverageTally::new(t),
        }
    }
}

impl ElectronTally for Probe {
    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _l: f64) {
        // The layer the flight was in: a flight that reaches a layer face
        // ends on it, still in its own layer (the crossing comes after).
        self.flights[end.layer] += 1;
    }
    fn table_lookup(&mut self, at: &ElectronState, ch: TableChannel, c: GridCoverage) {
        self.lookups[at.layer][ch as usize] += 1;
        self.cov.table_lookup(at, ch, c);
    }
    fn elastic(&mut self, _a: &ElectronState, _th: f64) {
        self.collisions += 1;
    }
    fn inelastic(&mut self, _a: &ElectronState, _w: f64) {
        self.collisions += 1;
        self.inelastic += 1;
    }
    fn reflected(&mut self, _a: &ElectronState, _b: lindhard::electron::transport::Boundary) {
        self.reflections += 1;
    }
    fn interface(&mut self, _a: &ElectronState, _f: usize, _t: usize) {
        self.interfaces += 1;
    }
    fn begin_secondary(&mut self, _s: &ElectronState, _g: u32) {
        self.secondaries += 1;
    }
    fn end_secondary(&mut self, f: Fate) {
        self.trapped += u64::from(f == Fate::Trapped);
    }
    fn end_history(&mut self, _i: u64, f: Fate) {
        self.trapped += u64::from(f == Fate::Trapped);
        self.stopped += u64::from(f == Fate::Stopped);
    }
    fn merge(&mut self, o: Self) {
        for (a, b) in self.flights.iter_mut().zip(&o.flights) {
            *a += b;
        }
        for (a, b) in self.lookups.iter_mut().zip(&o.lookups) {
            a[0] += b[0];
            a[1] += b[1];
        }
        self.collisions += o.collisions;
        self.inelastic += o.inelastic;
        self.stopped += o.stopped;
        self.secondaries += o.secondaries;
        self.reflections += o.reflections;
        self.interfaces += o.interfaces;
        self.trapped += o.trapped;
        self.cov.merge_from(&o.cov);
    }
}

/// Every pass of the loop evaluates both tables of the layer once, and a
/// pass ends in a flight (to a collision or a face) or, with zero total
/// rate and no face ahead, in [`Fate::Trapped`] (impossible here: every
/// layer scatters). So the evaluations of each channel equal the flights
/// started in that layer: more than the collisions, by the flights cut short
/// at a face (each crossing or reflection is followed by a new evaluation).
#[test]
fn evaluations_equal_flights_not_collisions() {
    let t = mixed();
    let p = t
        .run(11, 200, 16, &Primary::normal(300.0), || Probe::new(&t))
        .unwrap()
        .tally;
    assert_eq!(p.trapped, 0);
    let mut total = 0;
    for (l, f) in p.flights.iter().enumerate() {
        assert_eq!(p.lookups[l], [*f, *f], "layer {l}");
        let c = &p.cov.layers()[l];
        assert_eq!(c.elastic.total(), *f);
        assert_eq!(c.inelastic.total(), *f);
        total += f;
    }
    assert!(
        total > p.collisions,
        "{total} evaluations, {} collisions",
        p.collisions
    );
    assert!(p.secondaries > 0 && p.reflections > 0 && p.interfaces > 0);
}

fn mixed_report(threads: usize) -> ElectronReport {
    let t = mixed();
    let proto = FullElectronTally::new(&t, tally_config(300.0)).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| {
            t.run(0x253, 400, 16, &Primary::normal(300.0), || proto.clone())
                .unwrap()
                .tally
                .report()
        })
}

#[test]
fn coverage_is_identical_for_1_2_and_8_threads() {
    let r1 = mixed_report(1);
    let c = &r1.table_coverage;
    assert_eq!(c.len(), 2);
    for l in c {
        assert_eq!(l.elastic.total(), l.inelastic.total());
        assert!(l.elastic.total() > 0);
    }
    // All three classes occur somewhere.
    let all = c.iter().flat_map(|l| [l.elastic, l.inelastic]);
    let (b, w, a) = all.fold((0, 0, 0), |(b, w, a), x| {
        (b + x.below, w + x.within, a + x.above)
    });
    assert!(b > 0 && w > 0 && a > 0, "below {b}, within {w}, above {a}");
    assert_eq!(c[0].elastic.energy_min_ev, 30.0);
    assert_eq!(c[0].elastic.energy_max_ev, 300.0);
    assert_eq!(c[1].inelastic.energy_min_ev, 40.0);
    assert_eq!(c[1].inelastic.energy_max_ev, 400.0);
    let j1 = serde_json::to_string(&r1).unwrap();
    for n in [2, 8] {
        let r = mixed_report(n);
        assert_eq!(r1, r, "report differs on {n} threads");
        assert_eq!(j1, serde_json::to_string(&r).unwrap());
    }
}

/// A tally that forwards every hook to a [`SummaryTally`] and, optionally,
/// counts table coverage too.
#[derive(Clone)]
struct Wrapped {
    summary: SummaryTally,
    cov: Option<TableCoverageTally>,
    draws_seen: Vec<u64>,
}

impl ElectronTally for Wrapped {
    fn step(&mut self, from: [f64; 3], end: &ElectronState, l: f64) {
        self.summary.step(from, end, l);
        self.draws_seen
            .push(end.pos[0].to_bits() ^ end.energy_ev.to_bits());
    }
    fn table_lookup(&mut self, at: &ElectronState, ch: TableChannel, c: GridCoverage) {
        if let Some(cov) = &mut self.cov {
            cov.table_lookup(at, ch, c);
        }
    }
    fn elastic(&mut self, a: &ElectronState, th: f64) {
        self.summary.elastic(a, th);
    }
    fn inelastic(&mut self, a: &ElectronState, w: f64) {
        self.summary.inelastic(a, w);
    }
    fn secondary(
        &mut self,
        p: &ElectronState,
        e: &lindhard::electron::secondary::SecondaryEvent,
        c: Option<&ElectronState>,
    ) {
        self.summary.secondary(p, e, c);
    }
    fn reflected(&mut self, a: &ElectronState, b: lindhard::electron::transport::Boundary) {
        self.summary.reflected(a, b);
    }
    fn end_history(&mut self, i: u64, f: Fate) {
        self.summary.end_history(i, f);
    }
    fn merge(&mut self, mut o: Self) {
        self.summary.merge(o.summary);
        if let (Some(a), Some(b)) = (&mut self.cov, &o.cov) {
            a.merge_from(b);
        }
        self.draws_seen.append(&mut o.draws_seen);
    }
}

/// Counting coverage draws no random number and changes no state: a run
/// that counts gives bit-identical flights and summary to one that ignores
/// the hook. (That the physical results of the transport are those from
/// before the hook existed is pinned by `tests/golden.rs`, whose electron
/// batch is unchanged.)
#[test]
fn counting_changes_nothing() {
    let t = mixed();
    let p = Primary::normal(300.0);
    let run = |cov: Option<TableCoverageTally>| {
        let proto = Wrapped {
            summary: SummaryTally::default(),
            cov,
            draws_seen: Vec::new(),
        };
        t.run(19, 150, 8, &p, || proto.clone()).unwrap().tally
    };
    let off = run(None);
    let on = run(Some(TableCoverageTally::new(&t)));
    assert_eq!(off.summary, on.summary);
    assert_eq!(off.draws_seen, on.draws_seen);
    assert!(on.cov.unwrap().layers()[0].elastic.total() > 0);
}

// ---------------------------------------------------------------------------
// Serialization
// ---------------------------------------------------------------------------

#[test]
fn report_serializes_coverage_and_reads_older_reports_without_it() {
    let r = mixed_report(1);
    let mut v = serde_json::to_value(&r).unwrap();
    let c = &v["table_coverage"];
    assert_eq!(c.as_array().unwrap().len(), 2);
    for key in ["energy_min_ev", "energy_max_ev", "below", "within", "above"] {
        assert!(c[1]["inelastic"].get(key).is_some(), "{key}");
    }
    assert_eq!(c[1]["layer"], 1);
    assert_eq!(
        c[0]["elastic"]["below"].as_u64().unwrap(),
        r.table_coverage[0].elastic.below
    );
    let back: ElectronReport = serde_json::from_value(v.clone()).unwrap();
    assert_eq!(back.table_coverage, r.table_coverage);

    v.as_object_mut().unwrap().remove("table_coverage");
    let old: ElectronReport = serde_json::from_value(v).unwrap();
    assert!(old.table_coverage.is_empty());
    assert_eq!(old.histories, r.histories);
}

#[test]
#[should_panic(expected = "different tables")]
fn merging_coverage_for_different_tables_panics() {
    let a = TableCoverageTally::new(&mixed());
    let other = Transport::new(
        finite(&[5.0]),
        vec![pair(elastic_none(&[1.0, 2.0]), inelastic_none(&[1.0, 2.0]))],
        TransportConfig::new(0.5),
    )
    .unwrap();
    let mut b = TableCoverageTally::new(&other);
    b.merge_from(&a);
}
