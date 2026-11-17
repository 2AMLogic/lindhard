//! Position-dependent PSF: the energy left in a thin film as the beam's
//! landing point is scanned across a buried edge. The cross-section tables
//! are synthetic (constant mean free paths, isotropic angles), not measured
//! or tabulated data, so the step is a physical-consistency check and not a
//! comparison with a published measurement (see `docs/validation.md`).

use lindhard::electron::transport::{BeamSpread, LayerTables, Primary, Transport, TransportConfig};
use lindhard::geometry::{Boundary, VoxelGrid};
use lindhard::material::Material;
use lindhard::tally::{scan_landing_points, PsfScanConfig, PsfScanReport};

mod common;
use common::synthetic_scattering;
use common::with_threads;

const NM: f64 = 1e-9;
const E0: f64 = 500.0;
const FILM_NM: f64 = 20.0;

fn material(z: u8) -> Material {
    Material::from_atom_fractions(&[(z, 1.0)], None).unwrap()
}

/// Film (0), light substrate (1: few elastic events) and heavy substrate
/// (2: many elastic events per inelastic one, hence strong backscatter).
fn tables() -> Vec<LayerTables> {
    vec![
        synthetic_scattering(20.0 * NM, 15.0 * NM, 0.1),
        synthetic_scattering(40.0 * NM, 15.0 * NM, 0.1),
        synthetic_scattering(2.0 * NM, 15.0 * NM, 0.1),
    ]
}

/// A 20 nm film (5 nm voxels along `x`) on a substrate that is heavy for
/// `y < edge` and light beyond it. 40 columns of 10 nm along `y`.
fn grid(heavy: impl Fn(usize) -> bool, lateral: Boundary) -> VoxelGrid {
    let (nx, ny, nz) = (10, 40, 2);
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
    VoxelGrid::new(
        vec![material(14), material(14), material(29)],
        [nx, ny, nz],
        [5.0 * NM, 10.0 * NM, 1.0e-6],
        cells,
        [lateral, lateral, Boundary::Periodic],
    )
    .unwrap()
}

fn transport(g: VoxelGrid) -> Transport<VoxelGrid> {
    Transport::on_geometry(g, tables(), TransportConfig::new(10.0)).unwrap()
}

fn scan_config(ys_nm: &[f64], n: u64) -> PsfScanConfig {
    PsfScanConfig {
        depth_lo_m: 0.0,
        depth_hi_m: FILM_NM * NM,
        // Landing points are offsets from the centre of the front face
        // (y = 200 nm); the arguments are absolute y positions.
        landing_m: ys_nm.iter().map(|y| [(y - 200.0) * NM, 0.0]).collect(),
        histories_per_point: n,
        chunk_size: 64,
        seed: 20_261_110,
    }
}

fn run(t: &Transport<VoxelGrid>, ys_nm: &[f64], n: u64) -> PsfScanReport {
    scan_landing_points(
        t,
        &Primary::normal(E0),
        &BeamSpread::default(),
        &scan_config(ys_nm, n),
        "synthetic tables, test",
    )
    .unwrap()
}

#[test]
fn uniform_target_gives_position_independent_slab_energy() {
    let t = transport(grid(|_| false, Boundary::Periodic));
    let ys = [60.0, 120.0, 180.0, 240.0, 300.0];
    let r = run(&t, &ys, 4_000);
    let m = r
        .points
        .iter()
        .map(|p| p.energy_per_electron_ev)
        .sum::<f64>()
        / ys.len() as f64;
    for p in &r.points {
        assert!(p.std_err_ev > 0.0);
        let d = (p.energy_per_electron_ev - m).abs();
        assert!(d <= 4.5 * p.std_err_ev, "{:?} vs mean {m}", p);
    }
}

#[test]
fn buried_heavy_edge_shows_a_monotone_step_of_the_expected_sign() {
    // Heavy substrate for y < 200 nm, light beyond.
    let t = transport(grid(|j| j < 20, Boundary::Vacuum));
    let ys = [60.0, 110.0, 160.0, 190.0, 200.0, 210.0, 240.0, 290.0, 340.0];
    let r = run(&t, &ys, 8_000);
    let e: Vec<f64> = r.points.iter().map(|p| p.energy_per_electron_ev).collect();
    let se: Vec<f64> = r.points.iter().map(|p| p.std_err_ev).collect();
    let heavy = e[0];
    let light = *e.last().unwrap();
    // Sign and size: backscatter from the heavy substrate returns energy
    // through the film.
    assert!(
        heavy > 1.1 * light,
        "heavy {heavy} eV vs light {light} eV per electron"
    );
    // Plateaus on both sides, a monotone step between them (each pair of
    // neighbours may fluctuate by their combined error).
    for w in 0..e.len() - 1 {
        let tol = 3.0 * (se[w] * se[w] + se[w + 1] * se[w + 1]).sqrt();
        assert!(e[w + 1] <= e[w] + tol, "not monotone at {w}: {e:?}");
    }
    assert!((e[1] - heavy).abs() <= 4.5 * (se[0] * se[0] + se[1] * se[1]).sqrt());
    let n = e.len();
    assert!(
        (e[n - 2] - light).abs() <= 4.5 * (se[n - 1] * se[n - 1] + se[n - 2] * se[n - 2]).sqrt()
    );
    // The edge point is between the plateaus.
    let mid = e[4];
    assert!(
        mid < heavy && mid > light,
        "edge {mid}, plateaus {heavy}/{light}"
    );
    println!(
        "step: heavy {heavy:.4} +- {:.4} eV, light {light:.4} +- {:.4} eV, ratio {:.3}",
        se[0],
        se[n - 1],
        heavy / light
    );
}

#[test]
fn report_is_bit_identical_across_thread_counts_and_point_order_seeds() {
    let t = transport(grid(|j| j < 20, Boundary::Vacuum));
    let ys = [150.0, 200.0, 250.0];
    let a = with_threads(1, || run(&t, &ys, 1_500));
    let b = with_threads(2, || run(&t, &ys, 1_500));
    assert_eq!(a, b);
    assert_eq!(a.to_csv(), b.to_csv());
    // A point's result does not depend on the other points when its index
    // and the base seed are the same.
    let c = run(&t, &ys[..1], 1_500);
    assert_eq!(a.points[0], c.points[0]);
}

#[test]
fn csv_and_json_have_one_row_per_landing_point() {
    let t = transport(grid(|j| j < 20, Boundary::Vacuum));
    let r = run(&t, &[100.0, 300.0], 200);
    let csv = r.to_csv();
    assert_eq!(csv.lines().count(), 3);
    let json = serde_json::to_string(&r).unwrap();
    let back: PsfScanReport = serde_json::from_str(&json).unwrap();
    // serde_json's default float parsing is not exactly round-trip (1 ulp).
    assert_eq!(back.points.len(), 2);
    assert_eq!(back.source, r.source);
    for (a, b) in back.points.iter().zip(&r.points) {
        assert_eq!((a.seed, a.histories), (b.seed, b.histories));
        for (x, y) in [
            (a.energy_per_electron_ev, b.energy_per_electron_ev),
            (a.std_err_ev, b.std_err_ev),
            (a.landing_m[0], b.landing_m[0]),
        ] {
            assert!((x - y).abs() <= 1e-12 * y.abs());
        }
    }
}

#[test]
fn invalid_requests_are_rejected() {
    let t = transport(grid(|_| false, Boundary::Periodic));
    let go = |c: PsfScanConfig| {
        scan_landing_points(&t, &Primary::normal(E0), &BeamSpread::default(), &c, "x")
    };
    let mut c = scan_config(&[], 10);
    assert!(go(c.clone()).is_err());
    c = scan_config(&[100.0], 0);
    assert!(go(c).is_err());
    c = scan_config(&[100.0], 10);
    c.landing_m[0][0] = f64::NAN;
    assert!(go(c).is_err());
    c = scan_config(&[100.0], 10);
    c.depth_hi_m = c.depth_lo_m;
    assert!(go(c).is_err());
}
