//! The 2-D dynamic run (`ion::dynamic::DynamicRun2d`): the uniform-beam limit
//! against the 1-D run, finite-beam surface profile, atom conservation with
//! redeposition, thread-count determinism and edge cases.

use lindhard::geometry::Boundary;
use lindhard::ion::bca::Beam;
use lindhard::ion::dynamic::{
    BeamProfile, CompositionGrid, Dynamic2dConfig, DynamicConfig, DynamicRun, DynamicRun2d,
    Relaxation, Step2d, StepPolicy,
};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;

const NM: f64 = 1e-9;
const SI: u8 = 14;
const AS: u8 = 33;

mod common;
use common::bca_config as config;
use common::with_threads;
use common::zbl_si_table as table;

fn column(
    n_slabs: usize,
    slab_nm: f64,
    substrate: bool,
    relaxation: Relaxation,
) -> CompositionGrid {
    common::implanted_grid(n_slabs, slab_nm, substrate, relaxation, AS)
}

fn ideal() -> Relaxation {
    Relaxation::ideal_mixing_elemental(&[SI, AS]).unwrap()
}

fn beam(count: u64, polar_rad: f64) -> Beam {
    Beam {
        polar_rad,
        ..Beam::normal(Ion::new(AS).unwrap(), 1000.0, count)
    }
}

struct Out2d {
    columns: Vec<CompositionGrid>,
    surface: Vec<f64>,
    steps: Vec<Step2d>,
}

#[allow(clippy::too_many_arguments)]
fn run2d(
    columns: Vec<CompositionGrid>,
    h_nm: f64,
    count: u64,
    fluence: f64,
    ions_per_step: u64,
    profile: BeamProfile,
    lateral: Boundary,
    polar_rad: f64,
    seed: u64,
) -> Out2d {
    let ls = LindhardScharff::new();
    let mut run = DynamicRun2d::new(
        columns,
        h_nm * NM,
        beam(count, polar_rad),
        config(seed),
        &ls,
        table(),
        Dynamic2dConfig {
            fluence_m2: fluence,
            policy: StepPolicy::Adaptive {
                max_ions_per_step: ions_per_step,
                min_ions_per_step: 2,
                max_change: 0.3,
            },
            profile,
            lateral,
        },
    )
    .unwrap();
    let mut steps = Vec::new();
    run.run_all(|s, _| steps.push(s.clone())).unwrap();
    assert_eq!(run.ions_done(), count);
    Out2d {
        columns: run.columns().to_vec(),
        surface: run.surface_m().to_vec(),
        steps,
    }
}

fn fingerprint(o: &Out2d) -> Vec<u64> {
    let mut v: Vec<u64> = o.surface.iter().map(|r| r.to_bits()).collect();
    for g in &o.columns {
        for k in 0..g.n_slabs() {
            for (z, a) in g.inventory(k).unwrap() {
                v.extend([u64::from(z), a.to_bits()]);
            }
        }
        v.extend(g.thicknesses_m().iter().map(|t| t.to_bits()));
    }
    for s in &o.steps {
        v.extend([
            s.max_change.to_bits(),
            s.sputtered_atoms_per_m.to_bits(),
            s.redeposited_atoms_per_m.to_bits(),
            s.yields.sputtered_total(),
        ]);
    }
    v
}

#[test]
fn uniform_beam_reproduces_the_1d_run() {
    const N_COLS: usize = 4;
    const PER_COL: u64 = 500;
    let fluence = 4e19;
    let ls = LindhardScharff::new();
    // 1-D reference with the same total number of primaries.
    let mut one = DynamicRun::new(
        column(6, 3.0, true, ideal()),
        beam(N_COLS as u64 * PER_COL, 0.0),
        config(3),
        &ls,
        table(),
        DynamicConfig {
            fluence_m2: fluence,
            policy: StepPolicy::Fixed { ions_per_step: 400 },
            erosion: true,
        },
    )
    .unwrap();
    one.run_all(|_, _| {}).unwrap();
    let r1 = one.recession_m();
    let as1 = one.grid().total_inventory(AS);
    assert!(r1 > 0.2 * NM, "the reference must erode measurably: {r1}");

    let o = run2d(
        vec![column(6, 3.0, true, ideal()); N_COLS],
        6.0,
        PER_COL,
        fluence,
        100,
        BeamProfile::Uniform,
        Boundary::Periodic,
        0.0,
        3,
    );
    let r2 = o.surface.iter().sum::<f64>() / N_COLS as f64;
    let as2 = o.columns.iter().map(|g| g.total_inventory(AS)).sum::<f64>() / N_COLS as f64;
    // Recession and implanted dose agree within the counting noise
    // (about 1/sqrt(N) of a few thousand events each, here with margin).
    assert!(
        (r2 - r1).abs() < 0.15 * r1,
        "recession 2-D {r2:e} vs 1-D {r1:e}"
    );
    assert!(
        (as2 - as1).abs() < 0.08 * as1,
        "retained As 2-D {as2:e} vs 1-D {as1:e}"
    );
    // Composition-depth: As per slab depth (slabs of 3 nm), averaged over the
    // columns, against the 1-D profile. Compare the retained fractions in the
    // first two slabs, where most of the dose sits.
    let depth = |slabs: Vec<Vec<(u8, f64)>>| -> f64 {
        slabs
            .iter()
            .take(2)
            .flatten()
            .filter(|(z, _)| *z == AS)
            .map(|(_, a)| a)
            .sum()
    };
    let d1 = depth(
        (0..one.grid().n_slabs())
            .map(|k| one.grid().inventory(k).unwrap())
            .collect(),
    );
    let d2 = o
        .columns
        .iter()
        .map(|g| depth((0..g.n_slabs()).map(|k| g.inventory(k).unwrap()).collect()))
        .sum::<f64>()
        / N_COLS as f64;
    assert!(
        (d2 - d1).abs() < 0.12 * d1,
        "As in the top slabs {d2:e} vs {d1:e}"
    );
    // The surface stays flat to within the noise.
    let spread = o.surface.iter().map(|r| (r - r2).abs()).fold(0.0, f64::max);
    assert!(spread < 0.5 * r2.abs(), "spread {spread:e} vs mean {r2:e}");
}

#[test]
fn finite_beam_digs_a_profile() {
    let n = 12;
    let o = run2d(
        vec![column(6, 3.0, true, ideal()); n],
        5.0,
        300,
        1.6e20,
        100,
        BeamProfile::Gaussian {
            center_m: 30.0 * NM,
            fwhm_m: 16.0 * NM,
        },
        Boundary::Vacuum,
        0.0,
        5,
    );
    // Deepest at the beam centre, shallow in the wings.
    let (jmax, rmax) = o
        .surface
        .iter()
        .enumerate()
        .fold((0, f64::MIN), |a, (j, &r)| if r > a.1 { (j, r) } else { a });
    assert!((5..=6).contains(&jmax), "deepest column {jmax}");
    assert!(rmax > 0.5 * NM, "centre recession {rmax:e}");
    assert!(o.surface[0] < 0.25 * rmax && o.surface[n - 1] < 0.25 * rmax);
    // The surface profile is reported against fluence: the centre deepens
    // monotonically (step records carry it).
    let centre: Vec<f64> = o.steps.iter().map(|s| s.surface_m[5]).collect();
    assert!(
        centre.windows(2).all(|w| w[1] >= w[0] - 1e-12),
        "{centre:?}"
    );
    // Local incidence angle: zero on the flat first step, non-zero on the
    // flanks later.
    assert!(o.steps[0].incidence_rad.iter().all(|&a| a.abs() < 1e-12));
    let last = o.steps.last().unwrap();
    assert!(
        last.incidence_rad.iter().any(|&a| a > 0.05),
        "{:?}",
        last.incidence_rad
    );
    let f: f64 = o
        .steps
        .iter()
        .map(|s| s.redeposited_atoms_per_m)
        .sum::<f64>()
        / o.steps.iter().map(|s| s.sputtered_atoms_per_m).sum::<f64>();
    assert!((0.0..=1.0).contains(&f));
}

#[test]
fn atoms_are_conserved_with_redeposition() {
    // A film without a substrate: nothing can leave but through the front.
    let n = 10;
    let o = run2d(
        vec![column(8, 3.0, false, ideal()); n],
        5.0,
        200,
        6e20,
        100,
        BeamProfile::Gaussian {
            center_m: 25.0 * NM,
            fwhm_m: 8.0 * NM,
        },
        Boundary::Periodic,
        0.0,
        9,
    );
    let initial: f64 = {
        let g = column(8, 3.0, false, ideal());
        g.total_inventory(SI) * 5.0 * NM * n as f64
    };
    let mut total = initial;
    let mut redeposited = 0.0;
    for s in &o.steps {
        assert_eq!(s.yields.transmitted, 0);
        assert_eq!(s.yields.recoils_transmitted, 0);
        assert_eq!(s.lateral_escapes, 0);
        assert_eq!(s.clamped, 0);
        total += s.implanted_atoms_per_m - s.sputtered_atoms_per_m + s.redeposited_atoms_per_m;
        redeposited += s.redeposited_atoms_per_m;
    }
    assert!(redeposited > 0.0, "the crater must redeposit something");
    let now: f64 = o
        .columns
        .iter()
        .map(|g| (g.total_inventory(SI) + g.total_inventory(AS)) * 5.0 * NM)
        .sum();
    assert!(
        (now - total).abs() < 1e-9 * initial,
        "inventory {now:e} vs bookkeeping {total:e}"
    );
}

#[test]
fn bit_identical_at_1_2_and_8_threads() {
    let go = || {
        run2d(
            vec![column(5, 3.0, true, ideal()); 8],
            5.0,
            120,
            6e19,
            40,
            BeamProfile::Gaussian {
                center_m: 20.0 * NM,
                fwhm_m: 14.0 * NM,
            },
            Boundary::Periodic,
            0.3,
            21,
        )
    };
    let reference = fingerprint(&with_threads(1, go));
    for n in [2, 8] {
        assert_eq!(reference, fingerprint(&with_threads(n, go)), "{n} threads");
    }
}

#[test]
fn fixed_number_density_convention_is_supported() {
    let rel = Relaxation::fixed_number_density(5.0e28).unwrap();
    let o = run2d(
        vec![column(5, 3.0, true, rel); 4],
        6.0,
        100,
        4e19,
        50,
        BeamProfile::Uniform,
        Boundary::Periodic,
        0.0,
        4,
    );
    assert!(o.surface.iter().all(|r| r.is_finite() && *r > 0.0));
}

#[test]
fn hole_reaching_the_substrate_and_wide_beams_and_grazing_incidence() {
    // A thin film that the beam cuts through: the column becomes the bare
    // substrate and recession stops there.
    let o = run2d(
        vec![column(2, 1.5, true, ideal()); 5],
        4.0,
        150,
        6e20,
        50,
        BeamProfile::Gaussian {
            center_m: 10.0 * NM,
            fwhm_m: 6.0 * NM,
        },
        Boundary::Vacuum,
        0.0,
        2,
    );
    assert_eq!(
        o.columns[2].n_slabs(),
        0,
        "the hole must reach the substrate"
    );
    assert!(o.surface.iter().all(|r| r.is_finite()));
    // Beam much wider than the target: all columns illuminated, near uniform.
    let ls = LindhardScharff::new();
    let run = DynamicRun2d::new(
        vec![column(3, 3.0, true, ideal()); 4],
        5.0 * NM,
        beam(10, 0.0),
        config(1),
        &ls,
        table(),
        Dynamic2dConfig {
            fluence_m2: 1e19,
            policy: StepPolicy::Fixed { ions_per_step: 10 },
            profile: BeamProfile::Gaussian {
                center_m: 10.0 * NM,
                fwhm_m: 1.0e-3,
            },
            lateral: Boundary::Vacuum,
        },
    )
    .unwrap();
    assert!(run.profile_weights().iter().all(|&w| w > 0.999));
    // Grazing incidence on a structured surface stays finite.
    let o = run2d(
        vec![column(4, 3.0, true, ideal()); 6],
        5.0,
        80,
        8e19,
        40,
        BeamProfile::Gaussian {
            center_m: 15.0 * NM,
            fwhm_m: 10.0 * NM,
        },
        Boundary::Periodic,
        1.45,
        6,
    );
    assert!(o.surface.iter().all(|r| r.is_finite()));
    assert!(o.steps.iter().all(|s| s.max_change.is_finite()));
}
