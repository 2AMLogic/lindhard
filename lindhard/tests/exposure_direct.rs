//! Direct Monte Carlo exposure (`lindhard::exposure::direct_dose`) against
//! the Gaussian-sum superposition: agreement on a laterally uniform stack,
//! flagged differences over a buried lateral step, shaped spots, bit-identical
//! results at any thread count and the edge cases.
//!
//! The cross-section tables are computed here from simple formulas (constant
//! mean free paths, isotropic angles), not measured or tabulated data. The
//! film and the light substrate have no elastic scattering, so a primary
//! deposits its energy at the lateral position it starts from and the
//! point-spread function of the uniform stack is exactly the beam profile
//! (the superposition is then an exact model, with no fit error to allow
//! for); the heavy substrate scatters strongly and returns energy to the
//! film, which the superposition cannot know about.

use lindhard::electron::data::SamplingAxis;
use lindhard::electron::transport::{LayerTables, Primary, Transport, TransportConfig};
use lindhard::exposure::{
    compare_with_convolution, comparison_csv, direct_dose, DirectConfig, Exposure, Grid, Shot,
    SpotShape,
};
use lindhard::geometry::{Boundary, Stack, VoxelGrid};
use lindhard::material::Material;
use lindhard::tally::{
    Binning, ElectronTallyConfig, FullElectronTally, GaussianPsf, LogRadialBinning, PsfConfig,
};

mod common;
use common::{synthetic_scattering, synthetic_table, with_threads};

const NM: f64 = 1e-9;
const E0: f64 = 500.0;
const FILM_NM: f64 = 20.0;
/// Tolerance of the agreement and difference tests, in combined standard
/// errors.
const TOL_SIGMAS: f64 = 5.0;
/// Relative uncertainty assigned to the superposition: the pencil-beam run
/// that sets `E` has 60 000 histories (about 0.3 % statistical error on `E`),
/// so 2 % covers it generously.
const CONV_REL_ERR: f64 = 0.02;
/// Points below this fraction of the grid's peak dose are not judged: the run
/// (10 000 histories per replica) does not resolve them.
const FLOOR: f64 = 0.01;
const PSF_HISTORIES: u64 = 60_000;
/// Global beam blur `b` (in `exp(-r²/b²)`): wide against the 20 nm cell and
/// against the sub-sampling pitch of the cell average.
const BLUR_NM: f64 = 15.0;

fn material(z: u8) -> Material {
    Material::from_atom_fractions(&[(z, 1.0)], None).unwrap()
}

/// Film (0) and light substrate (1): no elastic scattering (a zero rate);
/// heavy substrate (2): strong elastic scattering. Inelastic losses uniform
/// on `[0, 0.1 E]`.
fn tables() -> Vec<LayerTables> {
    let straight = || LayerTables {
        elastic: synthetic_table(SamplingAxis::ElasticPolarAngle, 0.0, |_, p| p),
        inelastic: synthetic_table(
            SamplingAxis::InelasticEnergyLoss,
            1.0 / (15.0 * NM),
            |e, p| 0.1 * e * p,
        ),
    };
    vec![
        straight(),
        straight(),
        synthetic_scattering(2.0 * NM, 15.0 * NM, 0.1),
    ]
}

fn config() -> TransportConfig {
    TransportConfig::new(10.0)
}

/// Film on a light substrate as a laterally infinite stack: the reference
/// for the pencil-beam PSF.
fn light_stack() -> Transport {
    let stack = Stack::new(
        vec![(material(14), FILM_NM * NM), (material(14), 80.0 * NM)],
        None,
    )
    .unwrap();
    Transport::new(stack, tables()[..2].to_vec(), config()).unwrap()
}

/// Voxel target: 20 voxels of 5 nm in depth (film 20 nm, substrate 80 nm),
/// `ny` columns of 10 nm along `y` (offset 0 at the middle), periodic
/// laterally; the substrate is heavy where `heavy(j)`.
fn voxel_target(ny: usize, heavy: impl Fn(usize) -> bool) -> Transport<VoxelGrid> {
    let (nx, nz) = (20, 2);
    let mut cells = Vec::new();
    for _k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                cells.push(if i < 4 {
                    0
                } else if heavy(j) {
                    2
                } else {
                    1
                });
            }
        }
    }
    let g = VoxelGrid::new(
        vec![material(14), material(14), material(29)],
        [nx, ny, nz],
        [5.0 * NM, 10.0 * NM, 1.0e-6],
        cells,
        [Boundary::Vacuum, Boundary::Periodic, Boundary::Periodic],
    )
    .unwrap();
    Transport::on_geometry(g, tables(), config()).unwrap()
}

/// The PSF of the uniform stack: the energy one primary leaves in the film
/// (a pencil-beam run on the laterally infinite stack, slab total over
/// histories) in a delta-like double Gaussian (widths far below every
/// length used here, negligible second term). The beam blur is the
/// exposure's.
fn psf() -> GaussianPsf {
    let t = light_stack();
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, E0, 1).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 1).unwrap(),
    );
    c.psf = Some(
        PsfConfig::new(
            LogRadialBinning::new(1.0, 2.0, 1).unwrap(),
            0.0,
            FILM_NM * NM,
        )
        .unwrap(),
    );
    let proto = FullElectronTally::new(&t, c).unwrap();
    let profile = t
        .run(0xF17, PSF_HISTORIES, 64, &Primary::normal(E0), || {
            proto.clone()
        })
        .unwrap()
        .tally
        .report()
        .deposition
        .psf
        .unwrap();
    let e = profile.total_ev / profile.histories as f64;
    GaussianPsf::double(e, 0.3 * NM, 0.6 * NM, 1e-9).unwrap()
}

fn direct_config(n: u64) -> DirectConfig {
    DirectConfig {
        histories_per_shot: n,
        replicas: 8,
        chunk_size: 64,
        seed: 20_261_110,
        depth_lo_m: 0.0,
        depth_hi_m: FILM_NM * NM,
    }
}

/// 36 × 3 cells of 20 nm, x from -380 nm, y centred on 0.
fn grid() -> Grid {
    Grid {
        x0_m: -380.0 * NM,
        y0_m: -20.0 * NM,
        dx_m: 20.0 * NM,
        dy_m: 20.0 * NM,
        nx: 36,
        ny: 3,
    }
}

fn shots() -> Vec<Shot> {
    vec![
        Shot::point(-200.0 * NM, 0.0, 1.0),
        Shot::point(300.0 * NM, 0.0, 2.0),
    ]
}

fn blurred(psf: GaussianPsf) -> Exposure {
    Exposure::new(psf, shots(), BLUR_NM * NM).unwrap()
}

fn run_direct(t: &Transport<VoxelGrid>, n: u64) -> lindhard::exposure::DirectDose {
    direct_dose(
        t,
        &Primary::normal(E0),
        &shots(),
        BLUR_NM * NM,
        &grid(),
        &direct_config(n),
    )
    .unwrap()
}

#[test]
fn uniform_stack_small_pattern_direct_agrees_with_convolution() {
    let t = voxel_target(120, |_| false);
    let d = run_direct(&t, 80_000);
    let cmp =
        compare_with_convolution(&d, &blurred(psf()), CONV_REL_ERR, TOL_SIGMAS, FLOOR).unwrap();
    println!("{}", comparison_csv(&cmp));
    let flagged: Vec<_> = cmp.iter().filter(|c| c.flagged).collect();
    assert!(
        flagged.is_empty(),
        "flagged on a uniform stack: {flagged:?}"
    );
    // Not vacuous: the doses are substantial and carry a statistical error.
    let peak = cmp
        .iter()
        .map(|c| c.convolution_ev_per_m2)
        .fold(0.0, f64::max);
    let busy = cmp
        .iter()
        .filter(|c| c.convolution_ev_per_m2 > FLOOR * peak)
        .count();
    assert!(busy >= 6, "only {busy} points carry dose");
    assert!(d.std_err_ev_per_m2.iter().any(|&e| e > 0.0));
}

#[test]
fn buried_lateral_step_is_flagged_where_the_results_differ() {
    // Heavy substrate for y < 0 (offset from the middle), light beyond.
    let t = voxel_target(120, |j| j < 60);
    let d = run_direct(&t, 80_000);
    let cmp =
        compare_with_convolution(&d, &blurred(psf()), CONV_REL_ERR, TOL_SIGMAS, FLOOR).unwrap();
    println!("{}", comparison_csv(&cmp));
    let near =
        |x_nm: f64, c: &&lindhard::exposure::PointComparison| (c.x_m / NM - x_nm).abs() <= 40.0;
    // Over the heavy substrate the backscatter adds dose the uniform
    // convolution lacks: flagged, with the direct dose the larger.
    let heavy: Vec<_> = cmp.iter().filter(|c| near(-200.0, c)).collect();
    assert!(heavy.iter().any(|c| c.flagged && c.deviation_sigmas > 0.0));
    // Over the light region, far from the step, the two agree.
    let light: Vec<_> = cmp.iter().filter(|c| near(300.0, c)).collect();
    assert!(light.iter().all(|c| !c.flagged), "{light:?}");
}

#[test]
fn step_directly_under_a_shot_gives_finite_flagged_difference() {
    let t = voxel_target(120, |j| j < 60);
    let on_edge = [Shot::point(0.0, 0.0, 1.0)];
    let g = Grid {
        x0_m: -100.0 * NM,
        y0_m: 0.0,
        dx_m: 20.0 * NM,
        dy_m: 20.0 * NM,
        nx: 11,
        ny: 1,
    };
    let d = direct_dose(
        &t,
        &Primary::normal(E0),
        &on_edge,
        BLUR_NM * NM,
        &g,
        &direct_config(80_000),
    )
    .unwrap();
    let conv = Exposure::new(psf(), on_edge.to_vec(), BLUR_NM * NM).unwrap();
    let cmp = compare_with_convolution(&d, &conv, CONV_REL_ERR, TOL_SIGMAS, FLOOR).unwrap();
    assert!(cmp.iter().all(|c| c.direct_ev_per_m2.is_finite()));
    assert!(cmp.iter().any(|c| c.flagged && c.deviation_sigmas > 0.0));
}

#[test]
fn rectangle_spot_fills_its_footprint_and_nothing_else() {
    // No blur and no lateral scatter: a uniform rectangle of 80 x 40 nm
    // deposits E per primary spread over its area, in cells of 20 nm.
    let t = voxel_target(120, |_| false);
    let mut rect = Shot::point(0.0, 0.0, 3.0);
    rect.shape = SpotShape::Rect {
        width_m: 80.0 * NM,
        height_m: 40.0 * NM,
    };
    let g = Grid {
        x0_m: -90.0 * NM,
        y0_m: -70.0 * NM,
        dx_m: 20.0 * NM,
        dy_m: 20.0 * NM,
        nx: 10,
        ny: 8,
    };
    let d = direct_dose(
        &t,
        &Primary::normal(E0),
        &[rect],
        0.0,
        &g,
        &direct_config(160_000),
    )
    .unwrap();
    let e = psf().total_ev;
    let want = 3.0 * e / (80.0 * NM * 40.0 * NM);
    for j in 0..g.ny {
        for i in 0..g.nx {
            let [x, y] = g.point(i, j);
            let k = j * g.nx + i;
            let inside = x.abs() < 40.0 * NM && y.abs() < 20.0 * NM;
            if inside {
                let tol = TOL_SIGMAS * d.std_err_ev_per_m2[k] + CONV_REL_ERR * want;
                assert!(
                    (d.dose_ev_per_m2[k] - want).abs() <= tol,
                    "cell ({i},{j}): {} vs {want}",
                    d.dose_ev_per_m2[k]
                );
            } else {
                assert_eq!(d.dose_ev_per_m2[k], 0.0, "cell ({i},{j}) outside the spot");
            }
        }
    }
}

#[test]
fn result_is_bit_identical_for_1_2_and_4_threads() {
    let t = voxel_target(120, |j| j < 60);
    let run = |threads| {
        with_threads(threads, || {
            direct_dose(
                &t,
                &Primary::normal(E0),
                &shots(),
                0.0,
                &grid(),
                &direct_config(4_000),
            )
            .unwrap()
        })
    };
    let r1 = run(1);
    assert!(r1.dose_ev_per_m2.iter().any(|&v| v > 0.0));
    assert_eq!(r1, run(2));
    assert_eq!(r1, run(4));
}

#[test]
fn shot_outside_the_grid_and_zero_weight_shot_deposit_nothing_on_it() {
    let t = voxel_target(120, |_| false);
    let g = Grid {
        x0_m: -40.0 * NM,
        y0_m: 0.0,
        dx_m: 20.0 * NM,
        dy_m: 20.0 * NM,
        nx: 5,
        ny: 1,
    };
    let cfg = direct_config(800);
    let run = |shots: &[Shot]| direct_dose(&t, &Primary::normal(E0), shots, 0.0, &g, &cfg).unwrap();
    let far = run(&[Shot::point(400.0 * NM, 0.0, 1.0)]);
    assert!(far.dose_ev_per_m2.iter().all(|&v| v == 0.0));
    let zero = run(&[Shot::point(0.0, 0.0, 0.0)]);
    assert!(zero.dose_ev_per_m2.iter().all(|&v| v == 0.0));
    // A zero-weight shot adds nothing to a real one.
    let base = run(&[Shot::point(0.0, 0.0, 1.0)]);
    assert!(base.dose_ev_per_m2.iter().any(|&v| v > 0.0));
    let with_zero = run(&[
        Shot::point(0.0, 0.0, 1.0),
        Shot::point(100.0 * NM, 0.0, 0.0),
    ]);
    assert_eq!(base.dose_ev_per_m2, with_zero.dose_ev_per_m2);
}

#[test]
fn shaped_spot_with_blur_and_bad_settings_are_rejected() {
    let t = voxel_target(120, |_| false);
    let g = grid();
    let ok = direct_config(800);
    let p = Primary::normal(E0);
    let mut disc = Shot::point(0.0, 0.0, 1.0);
    disc.shape = SpotShape::Disc {
        radius_m: 10.0 * NM,
    };
    assert!(direct_dose(&t, &p, &[disc], 1.0 * NM, &g, &ok).is_err());
    let one = [Shot::point(0.0, 0.0, 1.0)];
    for bad in [
        DirectConfig { replicas: 1, ..ok },
        DirectConfig {
            histories_per_shot: 4,
            ..ok
        },
        DirectConfig {
            chunk_size: 0,
            ..ok
        },
    ] {
        assert!(direct_dose(&t, &p, &one, 0.0, &g, &bad).is_err());
    }
}
