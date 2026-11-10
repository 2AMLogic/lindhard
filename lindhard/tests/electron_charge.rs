//! The charge tally of `tally::electron`: every electron is counted exactly
//! once, so `incident + created = stopped + escaped` in integers, in every
//! ending a history can have, with and without the secondary model, in every
//! grid, on any thread count. Synthetic cross-section tables and band
//! parameters (made-up numbers chosen to exercise the code, not material
//! data).

use lindhard::electron::boundary::{BandModel, BandStructure};
use lindhard::electron::data::SamplingAxis;
use lindhard::electron::secondary::SecondaryModel;
use lindhard::electron::transport::{
    BoundaryModel, CutoffReference, EscapeRule, LayerTables, Primary, Transport, TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::tally::{
    Binning, CartesianGrid, ChargeCells, CylindricalGrid, ElectronReport, ElectronTallyConfig,
    FullElectronTally,
};

mod common;
use common::{bare_si, synthetic_table};

const NM: f64 = 1e-9;
const PROVENANCE: &str = "synthetic, tests/electron_charge.rs";

/// Isotropic elastic scattering and inelastic losses uniform on `[0, frac E]`.
fn scattering(lambda_el: f64, lambda_inel: f64, frac: f64) -> LayerTables {
    LayerTables {
        elastic: synthetic_table(SamplingAxis::ElasticPolarAngle, 1.0 / lambda_el, |_, p| {
            (1.0 - 2.0 * p).acos()
        }),
        inelastic: synthetic_table(
            SamplingAxis::InelasticEnergyLoss,
            1.0 / lambda_inel,
            move |e, p| frac * e * p,
        ),
    }
}

fn empty() -> LayerTables {
    LayerTables {
        elastic: synthetic_table(SamplingAxis::ElasticPolarAngle, 0.0, |_, _| 0.0),
        inelastic: synthetic_table(SamplingAxis::InelasticEnergyLoss, 0.0, |_, _| 0.0),
    }
}

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

fn config(e0: f64) -> ElectronTallyConfig {
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 50).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    c.cartesian = Some(CartesianGrid {
        x: Binning::new(0.0, 8.0 * NM, 8).unwrap(),
        y: Binning::new(-4.0 * NM, 4.0 * NM, 4).unwrap(),
        z: Binning::new(-4.0 * NM, 4.0 * NM, 4).unwrap(),
    });
    c.cylindrical = Some(CylindricalGrid {
        r: Binning::new(0.0, 6.0 * NM, 6).unwrap(),
        depth: Binning::new(0.0, 8.0 * NM, 8).unwrap(),
    });
    c
}

fn run(t: &Transport, n: u64, e0: f64, seed: u64) -> ElectronReport {
    let proto = FullElectronTally::new(t, config(e0)).unwrap();
    t.run(seed, n, 64, &Primary::normal(e0), || proto.clone())
        .unwrap()
        .tally
        .report()
}

fn sum_cells(c: &ChargeCells) -> (u64, u64) {
    let mut s = 0;
    let mut k = 0;
    for x in c.counts.iter().chain(std::iter::once(&c.outside)) {
        s += x.stopped();
        k += x.created();
    }
    (s, k)
}

/// The invariants of every report: the integer balance, the identity
/// `stopped - created = incident - escaped`, and that the layer, Cartesian and
/// cylindrical counts each add up to the totals.
fn check(r: &ElectronReport) {
    let c = &r.charge;
    assert_eq!(c.incident, r.histories);
    let escaped: u64 = [&c.escaped_front, &c.escaped_back]
        .iter()
        .map(|e| e.primary + e.secondary)
        .sum();
    assert_eq!(c.imbalance, 0, "charge not conserved: {c:?}");
    assert_eq!(
        c.total.net(),
        c.incident as i64 - escaped as i64,
        "net deposited charge differs from incident minus escaped"
    );
    // Escapes agree with the yields.
    assert_eq!(
        c.escaped_front.primary + c.escaped_front.secondary,
        r.front.count
    );
    assert_eq!(
        c.escaped_back.primary + c.escaped_back.secondary,
        r.back.count
    );
    let want = (c.total.stopped(), c.total.created());
    assert_eq!(sum_cells(&c.per_layer), want);
    assert_eq!(sum_cells(&c.cartesian.as_ref().unwrap().cells), want);
    assert_eq!(sum_cells(&c.cylindrical.as_ref().unwrap().cells), want);
    // A history's primary is the only primary: at most one primary ends in
    // the target per history.
    assert!(c.total.stopped_primary <= r.histories);
    // The squared-net sums bound the sums: sum d^2 >= |sum d| per cell.
    for g in [
        &c.per_layer,
        &c.cartesian.as_ref().unwrap().cells,
        &c.cylindrical.as_ref().unwrap().cells,
    ] {
        for (n, q) in g.counts.iter().zip(&g.net_sum_sq) {
            assert!(*q >= n.net().unsigned_abs());
        }
    }
}

#[test]
fn closes_without_secondaries_in_a_multilayer_with_both_faces() {
    let t = Transport::new(
        Stack::new(vec![(bare_si(), 3.0 * NM), (bare_si(), 6.0 * NM)], None).unwrap(),
        vec![
            scattering(2.0 * NM, 3.0 * NM, 0.3),
            scattering(2.0 * NM, 4.0 * NM, 0.3),
        ],
        TransportConfig::new(10.0),
    )
    .unwrap();
    let r = run(&t, 3_000, 500.0, 11);
    check(&r);
    let c = &r.charge;
    assert_eq!(c.total.created(), 0);
    assert_eq!(c.total.stopped_secondary, 0);
    assert_eq!(c.total.stopped_primary, r.fates.stopped);
    assert!(r.fates.escaped_back > 0 && r.fates.escaped_front > 0);
    assert_eq!(c.per_layer.counts.len(), 2);
    assert!(c.per_layer.counts.iter().all(|x| x.stopped() > 0));
}

#[test]
fn closes_for_each_terminal_fate_without_secondaries() {
    // Absorbed at the back face.
    let mut cfg = TransportConfig::new(10.0);
    cfg.escape_rule = EscapeRule::FrontOnly;
    let t = Transport::new(
        Stack::new(vec![(bare_si(), 6.0 * NM)], None).unwrap(),
        vec![scattering(3.0 * NM, 6.0 * NM, 0.2)],
        cfg,
    )
    .unwrap();
    let r = run(&t, 2_000, 500.0, 12);
    check(&r);
    assert!(r.fates.absorbed > 0);
    assert_eq!(
        r.charge.total.stopped_primary,
        r.histories - r.fates.escaped_front
    );

    // Trapped in an interaction-free substrate.
    let t = Transport::new(
        Stack::new(vec![(bare_si(), 5.0 * NM)], Some(bare_si())).unwrap(),
        vec![scattering(3.0 * NM, 6.0 * NM, 0.2), empty()],
        TransportConfig::new(10.0),
    )
    .unwrap();
    let r = run(&t, 2_000, 500.0, 13);
    check(&r);
    assert!(r.fates.trapped > 0);

    // Cut off by the collision cap.
    let mut cfg = TransportConfig::new(10.0);
    cfg.max_events = 4;
    let t = Transport::new(
        Stack::semi_infinite(bare_si()),
        vec![scattering(2.0 * NM, 3.0 * NM, 0.05)],
        cfg,
    )
    .unwrap();
    let r = run(&t, 1_000, 500.0, 14);
    check(&r);
    assert!(r.fates.event_capped > 0);
    assert_eq!(
        r.charge.total.stopped_primary,
        r.histories - r.fates.escaped_front
    );
}

fn with_secondaries(cutoff_ev: f64, max_events: u64) -> Transport {
    let mut cfg = TransportConfig::new(cutoff_ev);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    cfg.max_events = max_events;
    Transport::with_band_structures(
        Stack::new(vec![(bare_si(), 6.0 * NM)], Some(bare_si())).unwrap(),
        vec![
            scattering(3.0 * NM, 4.0 * NM, 0.3),
            scattering(2.0 * NM, 3.0 * NM, 0.25),
        ],
        vec![metal(), insulator()],
        cfg,
    )
    .unwrap()
}

#[test]
fn closes_with_secondaries_a_barrier_and_reflection() {
    let r = run(&with_secondaries(0.5, 100_000), 400, 300.0, 21);
    check(&r);
    let c = &r.charge;
    assert!(c.total.created() > 1_000);
    assert!(c.total.created_by_secondary > 0 && c.total.created_by_primary > 0);
    assert!(c.total.stopped_secondary > 0);
    assert!(c.escaped_front.secondary > 0, "no secondary escaped");
}

#[test]
fn closes_with_secondaries_cut_off_by_the_event_cap() {
    let r = run(&with_secondaries(0.5, 6), 400, 300.0, 22);
    check(&r);
    assert!(r.event_caps.secondary_tracks > 0 || r.fates.event_capped > 0);
}

#[test]
fn closes_at_every_cutoff() {
    // The cutoff sets which liberated electrons are followed and which are
    // dropped (and counted) where they are made; the balance holds at every
    // cutoff.
    let mut followed = Vec::new();
    for cutoff in [0.5, 3.0, 8.0] {
        let r = run(&with_secondaries(cutoff, 100_000), 300, 200.0, 23);
        check(&r);
        followed.push(r.charge.total.created() as f64 / r.histories as f64);
    }
    assert!(followed.iter().all(|&x| x > 0.0));
}

#[test]
fn report_is_bit_identical_on_1_2_and_4_threads_with_secondaries() {
    let t = with_secondaries(0.5, 100_000);
    let go = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| run(&t, 300, 200.0, 0xC4A9))
    };
    let one = go(1);
    check(&one);
    for n in [2, 4] {
        let r = go(n);
        assert_eq!(one.charge, r.charge, "charge differs on {n} threads");
    }
}

#[test]
fn reports_written_before_the_charge_tally_read_back() {
    let t = with_secondaries(0.5, 100_000);
    let r = run(&t, 50, 100.0, 24);
    let mut v = serde_json::to_value(&r).unwrap();
    v.as_object_mut().unwrap().remove("charge");
    let back: ElectronReport = serde_json::from_value(v).unwrap();
    assert_eq!(back.charge.incident, 0);
    assert!(back.charge.per_layer.counts.is_empty());
}
