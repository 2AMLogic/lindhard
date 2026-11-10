//! Lateral implant / damage profiles, snapshots, the film cut-through run and
//! the published-profile comparison of the 2-D dynamic run.

use lindhard::geometry::{Boundary, Stack};
use lindhard::ion::bca::{BcaConfig, Beam};
use lindhard::ion::dynamic::{
    compare_profile, BeamProfile, CompositionGrid, CutThrough, Dynamic2dConfig, DynamicRun2d,
    ProfileSnapshot, Relaxation, StepPolicy,
};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::tally::voxel::VoxelTallyGrid;

const NM: f64 = 1e-9;
const SI: u8 = 14;
const AS: u8 = 33;

mod common;
use common::si;
use common::zbl_si_table as table;

fn column(n_slabs: usize, slab_nm: f64, substrate: bool) -> CompositionGrid {
    let layers = (0..n_slabs).map(|_| (si(), slab_nm * NM)).collect();
    let stack = Stack::new(layers, substrate.then(si)).unwrap();
    let rel = Relaxation::ideal_mixing_elemental(&[SI, AS]).unwrap();
    let mut g = CompositionGrid::from_stack(&stack, rel).unwrap();
    g.seed_energies(AS, Some(15.0), None, Some(2.0));
    g
}

fn grid(n_cols: usize, h_nm: f64) -> VoxelTallyGrid {
    // 2 nm in depth over 40 nm, one voxel per column laterally.
    VoxelTallyGrid::new([0.0, 0.0, 0.0], [2.0 * NM, h_nm * NM, 1.0], [20, n_cols, 1]).unwrap()
}

struct Setup {
    cols: Vec<CompositionGrid>,
    h_nm: f64,
    count: u64,
    fluence: f64,
    step: u64,
    profile: BeamProfile,
    lateral: Boundary,
}

fn run<R>(
    s: Setup,
    snaps: &[f64],
    f: impl FnOnce(&mut DynamicRun2d<'_>) -> R,
) -> (R, DynamicRun2d<'static>) {
    // The stopping model and table must outlive the run.
    let ls: &'static LindhardScharff = Box::leak(Box::new(LindhardScharff::new()));
    let tab: &'static _ = Box::leak(Box::new(table()));
    let n = s.cols.len();
    let mut c = BcaConfig::new(5.0, 2.0);
    c.seed = 4;
    let mut r = DynamicRun2d::new(
        s.cols,
        s.h_nm * NM,
        Beam::normal(Ion::new(AS).unwrap(), 1000.0, s.count),
        c,
        ls,
        tab,
        Dynamic2dConfig {
            fluence_m2: s.fluence,
            policy: StepPolicy::Adaptive {
                max_ions_per_step: s.step,
                min_ions_per_step: 2,
                max_change: 0.3,
            },
            profile: s.profile,
            lateral: s.lateral,
        },
    )
    .unwrap()
    .with_profiles(grid(n, s.h_nm), snaps)
    .unwrap();
    let out = f(&mut r);
    (out, r)
}

fn uniform(count: u64, fluence: f64) -> Setup {
    Setup {
        cols: vec![column(4, 3.0, true); 4],
        h_nm: 5.0,
        count,
        fluence,
        step: 100,
        profile: BeamProfile::Uniform,
        lateral: Boundary::Periodic,
    }
}

#[test]
fn profile_totals_match_the_event_counts() {
    let (_, r) = run(uniform(300, 4e19), &[], |r| r.run_all(|_, _| {}).unwrap());
    let snap = r.profile_snapshot().unwrap();
    let t = snap.grand_total();
    assert_eq!(snap.beam_z, AS);
    // Uniform weights: every beam atom at rest (slab or substrate) is in the
    // profile, each standing for fluence/count atoms/m² over the width h.
    let y = r.cumulative();
    let w = 4e19 / 300.0 * 5.0 * NM;
    let rest = (y.primaries_in_slabs + y.primaries_in_substrate) as f64;
    assert!((t.implanted_per_m - rest * w).abs() < 1e-9 * rest * w);
    // Same displacement definition as the engine count (recoils created).
    assert!((t.displacements_per_m - y.displaced as f64 * w).abs() < 1e-9 * y.displaced as f64 * w);
    assert!(t.substrate_displacements_per_m <= t.displacements_per_m);
    // The beam stops near the surface: the depth profile peaks above 10 nm.
    let d = snap.depth_profile();
    let peak = (0..d.len())
        .max_by(|&a, &b| d[a].implanted_per_m.total_cmp(&d[b].implanted_per_m))
        .unwrap();
    assert!(peak < 5, "peak voxel {peak}");
}

#[test]
fn snapshots_land_on_the_requested_fluences() {
    let (_, r) = run(uniform(300, 4e19), &[1e19, 2.5e19, 4e19], |r| {
        r.run_all(|_, _| {}).unwrap()
    });
    let s = r.snapshots();
    assert_eq!(s.len(), 3);
    for (snap, want) in s.iter().zip([1e19, 2.5e19, 4e19]) {
        assert!((snap.peak_fluence_m2 - want).abs() < 0.01 * want);
    }
    // Cumulative: implanted atoms only grow.
    let tot: Vec<f64> = s.iter().map(|x| x.grand_total().implanted_per_m).collect();
    assert!(tot[0] < tot[1] && tot[1] < tot[2], "{tot:?}");
    // The last snapshot is the final state.
    assert_eq!(&s[2], &r.profile_snapshot().unwrap());
}

#[test]
fn snapshot_fluence_out_of_range_is_rejected() {
    let ls = LindhardScharff::new();
    let tab = table();
    let mk = |f: f64| {
        DynamicRun2d::new(
            vec![column(2, 3.0, true)],
            5.0 * NM,
            Beam::normal(Ion::new(AS).unwrap(), 1000.0, 10),
            BcaConfig::new(5.0, 2.0),
            &ls,
            tab,
            Dynamic2dConfig {
                fluence_m2: 1e19,
                policy: StepPolicy::Fixed { ions_per_step: 5 },
                profile: BeamProfile::Uniform,
                lateral: Boundary::Periodic,
            },
        )
        .unwrap()
        .with_profiles(grid(1, 5.0), &[f])
        .is_err()
    };
    assert!(mk(0.0) && mk(2e19) && mk(f64::NAN) && !mk(5e18));
}

fn thin_film_beam() -> Setup {
    Setup {
        cols: vec![column(2, 2.0, true); 9],
        h_nm: 4.0,
        count: 400,
        fluence: 3e21,
        step: 40,
        profile: BeamProfile::Gaussian {
            center_m: 18.0 * NM,
            fwhm_m: 12.0 * NM,
        },
        lateral: Boundary::Vacuum,
    }
}

fn csv_of(s: &ProfileSnapshot) -> String {
    s.to_csv()
}

#[test]
fn cut_through_stops_at_the_substrate_and_reports_the_damage_below() {
    let (ct, r): (CutThrough, _) = run(thin_film_beam(), &[], |r| {
        r.run_until_cut_through(|_, _| {}).unwrap()
    });
    assert!(ct.reached, "no cut-through: {ct:?}");
    let c = ct.column.unwrap();
    assert!((3..=5).contains(&c), "first bare column {c}");
    assert!(r.columns()[c].n_slabs() == 0);
    assert!(ct.fluence_before_m2 < ct.fluence_m2 && ct.fluence_m2 < 3e21);
    assert!(ct.ions_done < 400);
    let snap = ct.profiles.unwrap();
    let t = snap.grand_total();
    assert!(t.substrate_displacements_per_m > 0.0, "{t:?}");
    assert!(t.substrate_displacements_per_m < t.displacements_per_m);
    assert!(csv_of(&snap).contains("substrate_displacements_per_m3"));
}

#[test]
fn a_beam_that_never_cuts_through_reports_so() {
    let mut s = thin_film_beam();
    s.fluence = 3e19;
    let (ct, r) = run(s, &[], |r| r.run_until_cut_through(|_, _| {}).unwrap());
    assert!(!ct.reached && ct.column.is_none());
    assert_eq!(ct.ions_done, r.ions_total_per_column());
    assert!((ct.fluence_m2 - 3e19).abs() < 1e-6 * 3e19);
}

#[test]
fn cut_through_needs_a_substrate() {
    let mut s = thin_film_beam();
    s.cols = vec![column(2, 2.0, false); 9];
    let ls = LindhardScharff::new();
    let tab = table();
    let mut r = DynamicRun2d::new(
        s.cols,
        4.0 * NM,
        Beam::normal(Ion::new(AS).unwrap(), 1000.0, 10),
        BcaConfig::new(5.0, 2.0),
        &ls,
        tab,
        Dynamic2dConfig {
            fluence_m2: 1e19,
            policy: StepPolicy::Fixed { ions_per_step: 5 },
            profile: BeamProfile::Uniform,
            lateral: Boundary::Periodic,
        },
    )
    .unwrap();
    assert!(r.run_until_cut_through(|_, _| {}).is_err());
}

#[test]
fn profiles_are_bit_identical_at_1_and_8_threads() {
    let go = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let (ct, _) = run(thin_film_beam(), &[1e20], |r| {
                    r.run_until_cut_through(|_, _| {}).unwrap()
                });
                (ct.fluence_m2.to_bits(), ct.profiles.unwrap().to_csv())
            })
    };
    assert_eq!(go(1), go(8));
}

#[test]
fn comparison_metric_and_tolerance() {
    let y = [0.0, 10e-9, 20e-9, 30e-9];
    let d = [0.0, 10e-9, 10e-9, 0.0];
    // Points on the simulated curve: exact.
    let c = compare_profile(&y, &d, &[(5e-9, 5e-9), (15e-9, 10e-9)], 0.1).unwrap();
    assert!(c.rms_m < 1e-18 && c.within_tolerance);
    // A uniform 2 nm offset against 8 nm published depth: 25 % relative RMS.
    let c = compare_profile(&y, &d, &[(10e-9, 8e-9), (20e-9, 8e-9)], 0.1).unwrap();
    assert!((c.relative_rms - 0.25).abs() < 1e-12 && !c.within_tolerance);
    assert!((c.max_abs_m - 2e-9).abs() < 1e-18);
    // Out of range, bad input.
    assert!(compare_profile(&y, &d, &[(40e-9, 1e-9)], 0.1).is_err());
    assert!(compare_profile(&y, &d, &[], 0.1).is_err());
    assert!(compare_profile(&y, &d, &[(1e-9, 1e-9)], 0.0).is_err());
}

/// One small run per beam species (He, Ne, Ga, Xe) into a Si column: the report
/// schema, species tag and totals hold for each. Gas beams are given the Si
/// number density (fixed-density relaxation) since they have no solid volume
/// here; this checks bookkeeping, not a physical agreement.
#[test]
fn each_beam_species_reports_a_consistent_profile() {
    for (z, name) in [(2u8, "He"), (10, "Ne"), (31, "Ga"), (54, "Xe")] {
        let stack = Stack::new(vec![(si(), 3.0 * NM); 4], Some(si())).unwrap();
        let rel = Relaxation::fixed_number_density(4.996e28).unwrap();
        let mut g = CompositionGrid::from_stack(&stack, rel).unwrap();
        g.seed_energies(z, Some(15.0), None, Some(2.0));
        let ls = LindhardScharff::new();
        let tab = table();
        let mut c = BcaConfig::new(5.0, 2.0);
        c.seed = 4;
        let count = 200u64;
        let fluence = 2e19;
        let mut r = DynamicRun2d::new(
            vec![g; 3],
            5.0 * NM,
            Beam::normal(Ion::new(z).unwrap(), 1000.0, count),
            c,
            &ls,
            tab,
            Dynamic2dConfig {
                fluence_m2: fluence,
                policy: StepPolicy::Fixed { ions_per_step: 100 },
                profile: BeamProfile::Uniform,
                lateral: Boundary::Periodic,
            },
        )
        .unwrap()
        .with_profiles(grid(3, 5.0), &[fluence])
        .unwrap();
        r.run_all(|_, _| {}).unwrap();
        let snap = r.profile_snapshot().unwrap();
        assert_eq!(snap.beam_z, z, "{name}");
        let t = snap.grand_total();
        let y = r.cumulative();
        let w = fluence / count as f64 * 5.0 * NM;
        let rest = (y.primaries_in_slabs + y.primaries_in_substrate) as f64;
        assert!(rest > 0.0, "{name}: no ion came to rest");
        assert!(
            (t.implanted_per_m - rest * w).abs() < 1e-9 * rest * w,
            "{name}"
        );
        assert!(
            (t.displacements_per_m - y.displaced as f64 * w).abs()
                < 1e-9 * (y.displaced as f64 * w).max(1.0),
            "{name}"
        );
        assert!(
            snap.to_csv().contains("substrate_displacements_per_m3"),
            "{name}"
        );
        assert_eq!(r.snapshots().len(), 1, "{name}");
    }
}
