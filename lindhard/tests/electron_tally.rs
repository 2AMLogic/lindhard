//! `tally::electron` on the electron transport loop: energy balance in every
//! ending a history can have, grid and layer sums, the SE/BSE split and its
//! metadata, polar angles, generation-volume moments, and bit-identical
//! reports across thread counts. The cross-section tables are computed here
//! from simple formulas (constant mean free paths, isotropic angles), not
//! measured or tabulated data.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{
    ElectronState, ElectronTally, EscapeRule, Face, Fate, LayerTables, Primary, Transport,
    TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{
    Binning, CartesianGrid, CylindricalGrid, ElectronReport, ElectronTallyConfig,
    ElectronTallyError, FullElectronTally, SE_BSE_SPLIT_EV,
};

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];
const NM: f64 = 1e-9;

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn table(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
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
        provenance: "computed in tests/electron_tally.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// Isotropic elastic scattering (`cos θ = 1 - 2p`) and inelastic losses
/// uniform on `[0, frac E]`.
fn scattering(lambda_el: f64, lambda_inel: f64, frac: f64) -> LayerTables {
    LayerTables {
        elastic: table(SamplingAxis::ElasticPolarAngle, 1.0 / lambda_el, |_, p| {
            (1.0 - 2.0 * p).acos()
        }),
        inelastic: table(
            SamplingAxis::InelasticEnergyLoss,
            1.0 / lambda_inel,
            move |e, p| frac * e * p,
        ),
    }
}

/// No interaction at all.
fn empty() -> LayerTables {
    LayerTables {
        elastic: table(SamplingAxis::ElasticPolarAngle, 0.0, |_, _| 0.0),
        inelastic: table(SamplingAxis::InelasticEnergyLoss, 0.0, |_, _| 0.0),
    }
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

fn config(e0: f64) -> ElectronTallyConfig {
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 50).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    c.cartesian = Some(CartesianGrid {
        x: Binning::new(0.0, 40.0 * NM, 8).unwrap(),
        y: Binning::new(-20.0 * NM, 20.0 * NM, 6).unwrap(),
        z: Binning::new(-20.0 * NM, 20.0 * NM, 6).unwrap(),
    });
    c.cylindrical = Some(CylindricalGrid {
        r: Binning::new(0.0, 20.0 * NM, 10).unwrap(),
        depth: Binning::new(0.0, 40.0 * NM, 8).unwrap(),
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

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
}

/// The invariants every report must satisfy.
fn check_report(r: &ElectronReport) {
    let b = &r.budget;
    assert!(
        b.relative_imbalance < 1e-9,
        "energy balance off by {} (relative)",
        b.relative_imbalance
    );
    assert!(close(
        b.deposited_ev + b.escaped_ev + b.trapped_ev,
        b.incident_ev,
        1e-9
    ));
    let f = &r.fates;
    assert_eq!(
        f.stopped + f.escaped_front + f.escaped_back + f.absorbed + f.trapped + f.event_capped,
        r.histories
    );
    let d = &r.deposition;
    let layers: f64 = d.per_layer_ev.iter().sum();
    assert!(close(layers, b.deposited_ev, 1e-9));
    let cart = d.cartesian.as_ref().unwrap();
    assert!(close(
        cart.energy_ev.iter().sum::<f64>() + cart.outside_ev,
        b.deposited_ev,
        1e-9
    ));
    let cyl = d.cylindrical.as_ref().unwrap();
    assert!(close(
        cyl.energy_ev.iter().sum::<f64>() + cyl.outside_ev,
        b.deposited_ev,
        1e-9
    ));
    if let Some(g) = &r.generation_volume {
        assert!(close(g.energy_ev, b.deposited_ev, 1e-9));
    }
    for face in [&r.front, &r.back] {
        assert_eq!(face.count, face.slow.count + face.fast.count);
        assert_eq!(face.energy_histogram.total(), face.count);
        assert_eq!(face.slow.polar_histogram.total(), face.slow.count);
        assert_eq!(face.fast.polar_histogram.total(), face.fast.count);
    }
    assert!(close(r.front.energy_ev, b.escaped_front_ev, 1e-12));
    assert!(close(r.back.energy_ev, b.escaped_back_ev, 1e-12));
    assert!(close(
        r.yields.total_sigma,
        r.yields.backscatter_eta + r.yields.secondary_delta,
        1e-12
    ));
}

// ---------------------------------------------------------------------------
// Energy balance, one test per way a history can end
// ---------------------------------------------------------------------------

#[test]
fn balance_closes_with_stopping_and_backscatter_in_a_substrate() {
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![scattering(2.0 * NM, 3.0 * NM, 0.3)],
        TransportConfig::new(10.0),
    )
    .unwrap();
    let r = run(&t, 2_000, 1_000.0, 1);
    check_report(&r);
    assert!(r.fates.stopped > 0 && r.fates.escaped_front > 0);
    assert!(r.budget.residual_ev > 0.0 && r.budget.escaped_front_ev > 0.0);
    assert!(r.front.slow.count > 0 && r.front.fast.count > 0);
    assert_eq!(r.back.count, 0);
    assert_eq!(r.budget.trapped_ev, 0.0);
    assert_eq!(r.stopping_points.stopped, r.fates.stopped);
    // Without secondaries every stopping point is a primary's.
    assert_eq!(r.stopping_points.primaries.stopped, r.fates.stopped);
    assert_eq!(r.stopping_points.primaries.depth, r.stopping_points.depth);
    let g = r.generation_volume.unwrap();
    assert!(g.mean_m[0] > 0.0 && g.std_dev_m[0] > 0.0 && g.rms_radius_m > 0.0);
}

#[test]
fn balance_closes_with_transmission_through_a_finite_stack() {
    let t = Transport::new(
        finite(&[5.0, 10.0]),
        vec![
            scattering(4.0 * NM, 6.0 * NM, 0.2),
            scattering(3.0 * NM, 5.0 * NM, 0.2),
        ],
        TransportConfig::new(10.0),
    )
    .unwrap();
    let r = run(&t, 2_000, 1_000.0, 2);
    check_report(&r);
    assert!(r.fates.escaped_back > 0 && r.budget.escaped_back_ev > 0.0);
    assert!(r.yields.transmitted_fast > 0.0);
    assert!(r.deposition.per_layer_ev.iter().all(|&e| e > 0.0));
}

#[test]
fn balance_closes_with_back_face_absorption() {
    let mut cfg = TransportConfig::new(10.0);
    cfg.escape_rule = EscapeRule::FrontOnly;
    let t = Transport::new(
        finite(&[10.0]),
        vec![scattering(4.0 * NM, 6.0 * NM, 0.2)],
        cfg,
    )
    .unwrap();
    let r = run(&t, 2_000, 1_000.0, 3);
    check_report(&r);
    assert!(r.fates.absorbed > 0 && r.budget.absorbed_ev > 0.0);
    assert_eq!(r.back.count, 0);
}

#[test]
fn balance_closes_with_electrons_trapped_in_an_empty_substrate() {
    // Electrons entering the interaction-free substrate heading inward have
    // nothing to hit and no face to reach.
    let t = Transport::new(
        Stack::new(vec![(material(), 5.0 * NM)], Some(material())).unwrap(),
        vec![scattering(3.0 * NM, 6.0 * NM, 0.2), empty()],
        TransportConfig::new(10.0),
    )
    .unwrap();
    let r = run(&t, 2_000, 1_000.0, 4);
    check_report(&r);
    assert!(r.fates.trapped > 0 && r.budget.no_interaction_ev > 0.0);
    assert_eq!(r.deposition.per_layer_ev[1], 0.0);
}

#[test]
fn balance_closes_with_histories_cut_off_by_the_event_cap() {
    let mut cfg = TransportConfig::new(10.0);
    cfg.max_events = 5;
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![scattering(2.0 * NM, 3.0 * NM, 0.05)],
        cfg,
    )
    .unwrap();
    let r = run(&t, 1_000, 1_000.0, 5);
    check_report(&r);
    assert!(r.fates.event_capped > 0 && r.budget.event_cap_ev > 0.0);
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

fn threaded(threads: usize) -> ElectronReport {
    let t = Transport::new(
        finite(&[5.0, 30.0]),
        vec![
            scattering(3.0 * NM, 5.0 * NM, 0.3),
            scattering(2.0 * NM, 4.0 * NM, 0.3),
        ],
        TransportConfig::new(5.0),
    )
    .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| run(&t, 3_000, 2_000.0, 0xE1EC))
}

#[test]
fn report_is_bit_identical_for_1_2_and_8_threads() {
    let r1 = threaded(1);
    check_report(&r1);
    let j1 = serde_json::to_string(&r1).unwrap();
    for n in [2, 8] {
        let r = threaded(n);
        assert_eq!(r1, r, "report differs on {n} threads");
        assert_eq!(j1, serde_json::to_string(&r).unwrap());
    }
    // The serialized report reads back. serde_json's default float parser
    // may be one ulp off, so floats are compared to rounding only.
    let back: ElectronReport = serde_json::from_str(&j1).unwrap();
    assert_eq!(back.histories, r1.histories);
    assert_eq!(back.fates, r1.fates);
    assert_eq!(back.metadata.config, r1.metadata.config);
    assert_eq!(back.front.energy_histogram, r1.front.energy_histogram);
    assert!(close(
        back.budget.deposited_ev,
        r1.budget.deposited_ev,
        1e-15
    ));
    assert!(close(
        back.yields.backscatter_eta,
        r1.yields.backscatter_eta,
        1e-15
    ));
}

// ---------------------------------------------------------------------------
// SE/BSE split
// ---------------------------------------------------------------------------

fn transport() -> Transport {
    Transport::new(
        Stack::semi_infinite(material()),
        vec![scattering(2.0 * NM, 3.0 * NM, 0.3)],
        TransportConfig::new(10.0),
    )
    .unwrap()
}

fn state(pos: [f64; 3], dir: [f64; 3], energy_ev: f64) -> ElectronState {
    ElectronState {
        pos,
        dir,
        energy_ev,
        layer: 0,
    }
}

#[test]
fn split_defaults_to_50_ev_and_is_recorded() {
    let tally = FullElectronTally::new(&transport(), config(1_000.0)).unwrap();
    let m = tally.report().metadata;
    assert_eq!(SE_BSE_SPLIT_EV, 50.0);
    assert_eq!(m.se_bse_split_ev, 50.0);
    assert_eq!(m.config.se_bse_split_ev, 50.0);
    assert!(m.se_bse_split_source.contains("10.1038/s41598-022-20466-3"));
    assert!(m.se_bse_split_rule.contains("E < split"));
    assert_eq!(m.cutoff_ev, 10.0);
    assert_eq!(m.layers, 1);
}

#[test]
fn split_classes_below_as_slow_and_at_or_above_as_fast() {
    let mut tally = FullElectronTally::new(&transport(), config(1_000.0)).unwrap();
    let out = [-1.0, 0.0, 0.0];
    for (i, e) in [49.999, 50.0, 50.001, 3.0, 900.0].into_iter().enumerate() {
        tally.begin_history(i as u64, &state([0.0; 3], [1.0, 0.0, 0.0], e));
        tally.escaped(&state([0.0; 3], out, e), Face::Front);
        tally.end_history(i as u64, Fate::Escaped(Face::Front));
    }
    let r = tally.report();
    assert_eq!(r.front.slow.count, 2);
    assert_eq!(r.front.fast.count, 3);
    assert_eq!(r.yields.secondary_delta, 2.0 / 5.0);
    assert_eq!(r.yields.backscatter_eta, 3.0 / 5.0);
    assert!(close(
        r.front.slow.mean_energy_ev,
        (49.999 + 3.0) / 2.0,
        1e-15
    ));
    check_report(&r);
}

#[test]
fn split_is_configurable() {
    let t = transport();
    let reports: Vec<ElectronReport> = [SE_BSE_SPLIT_EV, 0.0, 1e9]
        .into_iter()
        .map(|split| {
            let mut c = config(1_000.0);
            c.se_bse_split_ev = split;
            let proto = FullElectronTally::new(&t, c).unwrap();
            t.run(9, 1_000, 64, &Primary::normal(1_000.0), || proto.clone())
                .unwrap()
                .tally
                .report()
        })
        .collect();
    let sigma = reports[0].yields.total_sigma;
    assert!(sigma > 0.0);
    for r in &reports {
        assert_eq!(r.yields.total_sigma, sigma);
        check_report(r);
    }
    assert_eq!(reports[1].metadata.se_bse_split_ev, 0.0);
    assert_eq!(reports[1].front.slow.count, 0);
    assert_eq!(reports[2].metadata.se_bse_split_ev, 1e9);
    assert_eq!(reports[2].front.fast.count, 0);
}

// ---------------------------------------------------------------------------
// Hand-driven events: polar angles, grids, moments, trapped classification
// ---------------------------------------------------------------------------

#[test]
fn polar_angle_is_measured_from_the_outward_normal() {
    let mut c = config(1_000.0);
    c.escape_polar = Binning::new(0.0, std::f64::consts::FRAC_PI_2, 90).unwrap();
    let mut tally = FullElectronTally::new(&transport(), c).unwrap();
    let a = 60f64.to_radians();
    // Front face: outward normal is -x. Back face: +x.
    tally.escaped(
        &state([0.0; 3], [-a.cos(), a.sin(), 0.0], 500.0),
        Face::Front,
    );
    tally.escaped(&state([0.0; 3], [-1.0, 0.0, 0.0], 500.0), Face::Front);
    tally.escaped(
        &state([1.0, 0.0, 0.0], [a.cos(), 0.0, -a.sin()], 5.0),
        Face::Back,
    );
    let r = tally.report();
    let fast = &r.front.fast.polar_histogram.counts;
    assert_eq!((fast[0], fast[60]), (1, 1));
    assert_eq!(r.back.slow.polar_histogram.counts[60], 1);
}

#[test]
fn deposits_land_in_the_right_cells_with_exact_moments() {
    let mut tally = FullElectronTally::new(&transport(), config(1_000.0)).unwrap();
    let dir = [1.0, 0.0, 0.0];
    // Two inelastic deposits and one sub-cutoff stop; no other events.
    tally.begin_history(0, &state([0.0; 3], dir, 100.0));
    let a = state([2.0 * NM, NM, -NM], dir, 70.0);
    tally.inelastic(&a, 30.0);
    let b = state([12.0 * NM, -6.0 * NM, 8.0 * NM], dir, 8.0);
    tally.inelastic(&b, 62.0);
    tally.stopped(&b);
    tally.end_history(0, Fate::Stopped);
    let r = tally.report();
    check_report(&r);
    assert_eq!(r.budget.inelastic_ev, 92.0);
    assert_eq!(r.budget.residual_ev, 8.0);
    assert_eq!(r.budget.deposited_ev, 100.0);
    let cart = r.deposition.cartesian.as_ref().unwrap();
    // x bins of 5 nm, y and z bins of 20/3 nm from -20 nm.
    assert_eq!(cart.at(0, 3, 2), 30.0);
    assert_eq!(cart.at(2, 2, 4), 70.0);
    assert_eq!(cart.outside_ev, 0.0);
    let cyl = r.deposition.cylindrical.as_ref().unwrap();
    // r bins of 2 nm: |a| = 2^(1/2) nm (bin 0), |b| = 10 nm (bin 5).
    assert_eq!(cyl.at(0, 0), 30.0);
    assert_eq!(cyl.at(5, 2), 70.0);
    let g = r.generation_volume.unwrap();
    let mean_x = (30.0 * 2.0 + 70.0 * 12.0) / 100.0 * NM;
    assert!(close(g.mean_m[0], mean_x, 1e-12));
    let var_x = (30.0 * (2.0 * NM - mean_x).powi(2) + 70.0 * (12.0 * NM - mean_x).powi(2)) / 100.0;
    assert!(close(g.std_dev_m[0], var_x.sqrt(), 1e-12));
    let rms = ((30.0 * 2.0 + 70.0 * 100.0) / 100.0f64).sqrt() * NM;
    assert!(close(g.rms_radius_m, rms, 1e-12));
    assert_eq!(r.stopping_points.stopped, 1);
}

#[test]
fn stopped_at_or_above_the_cutoff_counts_as_trapped() {
    let mut tally = FullElectronTally::new(&transport(), config(1_000.0)).unwrap();
    tally.begin_history(0, &state([0.0; 3], [1.0, 0.0, 0.0], 400.0));
    tally.stopped(&state([3.0 * NM, 0.0, 0.0], [1.0, 0.0, 0.0], 400.0));
    tally.end_history(0, Fate::Trapped);
    let r = tally.report();
    check_report(&r);
    assert_eq!(r.budget.no_interaction_ev, 400.0);
    assert_eq!(r.budget.deposited_ev, 0.0);
    assert_eq!(r.stopping_points.stopped, 0);
    assert!(r.generation_volume.is_none());
}

#[test]
fn event_cap_counts_the_last_seen_energy_as_trapped() {
    let mut tally = FullElectronTally::new(&transport(), config(1_000.0)).unwrap();
    let dir = [1.0, 0.0, 0.0];
    tally.begin_history(0, &state([0.0; 3], dir, 400.0));
    tally.step([0.0; 3], &state([1.0 * NM, 0.0, 0.0], dir, 400.0), 1.0 * NM);
    tally.inelastic(&state([1.0 * NM, 0.0, 0.0], dir, 250.0), 150.0);
    tally.end_history(0, Fate::EventCap);
    let r = tally.report();
    check_report(&r);
    assert_eq!(r.budget.event_cap_ev, 250.0);
    assert_eq!(r.fates.event_capped, 1);
}

// ---------------------------------------------------------------------------
// Configuration errors and merge guard
// ---------------------------------------------------------------------------

#[test]
fn bad_configurations_are_rejected() {
    let t = transport();
    let mut c = config(1_000.0);
    c.se_bse_split_ev = -1.0;
    assert_eq!(
        FullElectronTally::new(&t, c).unwrap_err(),
        ElectronTallyError::Split(-1.0)
    );
    c.se_bse_split_ev = f64::NAN;
    assert!(matches!(
        FullElectronTally::new(&t, c),
        Err(ElectronTallyError::Split(_))
    ));
    let mut c = config(1_000.0);
    c.cylindrical = Some(CylindricalGrid {
        r: Binning::new(-1.0, 1.0, 2).unwrap(),
        depth: Binning::new(0.0, 1.0, 2).unwrap(),
    });
    assert_eq!(
        FullElectronTally::new(&t, c).unwrap_err(),
        ElectronTallyError::NegativeRadius(-1.0)
    );
    let mut c = config(1_000.0);
    c.escape_polar = Binning::new(0.0, 4.0, 2).unwrap();
    assert!(matches!(
        FullElectronTally::new(&t, c),
        Err(ElectronTallyError::PolarRange { .. })
    ));
    let mut c = config(1_000.0);
    let huge = Binning::new(0.0, 1.0, usize::MAX / 2).unwrap();
    c.cartesian = Some(CartesianGrid {
        x: huge,
        y: huge,
        z: huge,
    });
    assert_eq!(
        FullElectronTally::new(&t, c).unwrap_err(),
        ElectronTallyError::GridTooLarge("Cartesian")
    );
}

#[test]
#[should_panic(expected = "different runs")]
fn merging_tallies_with_different_splits_panics() {
    let t = transport();
    let mut a = FullElectronTally::new(&t, config(1_000.0)).unwrap();
    let mut c = config(1_000.0);
    c.se_bse_split_ev = 20.0;
    let b = FullElectronTally::new(&t, c).unwrap();
    a.merge(b);
}
