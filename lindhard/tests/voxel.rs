//! The BCA engine on a voxel grid (`lindhard::geometry::VoxelGrid`): a grid
//! that encodes a 1D stack must reproduce the stack's profiles, periodic
//! lateral boundaries must be invisible, vacuum side faces must conserve
//! energy, and results must not depend on the thread count.

use std::sync::OnceLock;

use lindhard::geometry::{Boundary, Geometry, Stack, VoxelGrid};
use lindhard::ion::bca::{Bca, BcaConfig, BcaError, Beam, SummaryTally};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;

const NM: f64 = 1e-9;

fn table() -> &'static ScatteringTable {
    static T: OnceLock<ScatteringTable> = OnceLock::new();
    T.get_or_init(|| {
        ScatteringTable::build(
            &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e4,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 16,
            },
        )
    })
}

/// Test parameters (not data): `E_d` = 15 eV, `E_s` = 2 eV where unset.
fn ready(mut m: Material) -> Material {
    for (z, _) in m.atom_fractions() {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}

fn si() -> Material {
    ready(Material::from_atom_fractions(&[(14, 1.0)], None).unwrap())
}

fn sio2() -> Material {
    ready(Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap())
}

fn beam(count: u64) -> Beam {
    Beam {
        ion: Ion::new(5).unwrap(),
        energy_ev: 2.0e3,
        polar_rad: 0.3,
        azimuth_rad: 0.0,
        count,
    }
}

fn summary(bca: &Bca, bins: usize) -> SummaryTally {
    bca.run(|| SummaryTally::new(0.5 * NM, bins)).unwrap()
}

/// A column grid: `nx` voxels of `h` along depth with `material_at(i)`,
/// `ny x nz` transverse voxels of width `w`.
fn column(
    materials: Vec<Material>,
    nx: usize,
    h: f64,
    transverse: ([usize; 2], f64),
    lateral: Boundary,
    material_at: impl Fn(usize) -> u32,
) -> VoxelGrid {
    let ([ny, nz], w) = transverse;
    let mut cells = Vec::with_capacity(nx * ny * nz);
    for _ in 0..ny * nz {
        cells.extend((0..nx).map(&material_at));
    }
    VoxelGrid::new(
        materials,
        [nx, ny, nz],
        [h, w, w],
        cells,
        [Boundary::Vacuum, lateral, lateral],
    )
    .unwrap()
}

fn ks(a: &[u64], b: &[u64]) -> f64 {
    let (na, nb) = (a.iter().sum::<u64>() as f64, b.iter().sum::<u64>() as f64);
    let (mut ca, mut cb, mut d) = (0.0, 0.0, 0.0f64);
    for (x, y) in a.iter().zip(b) {
        ca += *x as f64 / na;
        cb += *y as f64 / nb;
        d = d.max((ca - cb).abs());
    }
    d
}

/// Same-sample-size agreement of two summaries (KS at alpha = 0.001, means
/// within 3 standard errors).
fn assert_same_profile(a: &SummaryTally, b: &SummaryTally, what: &str) {
    let n = a.primaries_stopped.min(b.primaries_stopped) as f64;
    assert!(n > 4000.0, "{what}: too few stopped ({n})");
    let d = ks(&a.depth_hist, &b.depth_hist);
    let d_crit = 1.95 * (2.0 / n).sqrt();
    assert!(d < d_crit, "{what}: KS {d} >= {d_crit}");
    let se = (a.depth_std().powi(2) / a.primaries_stopped as f64
        + b.depth_std().powi(2) / b.primaries_stopped as f64)
        .sqrt();
    assert!(
        (a.mean_depth() - b.mean_depth()).abs() < 3.0 * se,
        "{what}: {} vs {} (se {se})",
        a.mean_depth(),
        b.mean_depth()
    );
}

#[test]
fn voxel_encoded_stack_reproduces_the_layer_engine() {
    let ls = LindhardScharff::new();
    // 4 nm SiO2 on 26 nm Si, finite: back face at 30 nm. Voxels of 1 nm
    // along depth; the transverse axes are periodic.
    let stack = Stack::new(vec![(sio2(), 4.0 * NM), (si(), 26.0 * NM)], None).unwrap();
    let grid = column(
        vec![sio2(), si()],
        30,
        NM,
        ([2, 2], 20.0 * NM),
        Boundary::Periodic,
        |i| u32::from(i >= 4),
    );
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.follow_recoils = false;
    cfg.seed = 21;
    let b = beam(8000);
    let a = summary(&Bca::new(b, &stack, cfg, &ls, table()).unwrap(), 70);
    let v = summary(&Bca::new(b, &grid, cfg, &ls, table()).unwrap(), 70);
    assert_same_profile(&a, &v, "stack vs voxel");
    assert_eq!(a.histories, v.histories);
    // Same escape channels (backscatter within 4 sigma of a binomial).
    let p = a.backscattered as f64 / a.histories as f64;
    let sigma = (p * (1.0 - p) * a.histories as f64 * 2.0).sqrt().max(1.0);
    assert!(
        (a.backscattered as f64 - v.backscattered as f64).abs() < 4.0 * sigma,
        "{} vs {}",
        a.backscattered,
        v.backscattered
    );
    assert_eq!(v.lateral + v.recoils_lateral, 0);
    assert_eq!(v.budget.lateral, 0.0);

    // Negative control: a grid with a 20 % denser lower material must differ.
    let mut dense = si();
    dense.set_mass_density(1.2 * dense.mass_density()).unwrap();
    let other = column(
        vec![sio2(), dense],
        30,
        NM,
        ([2, 2], 20.0 * NM),
        Boundary::Periodic,
        |i| u32::from(i >= 4),
    );
    let c = summary(&Bca::new(b, &other, cfg, &ls, table()).unwrap(), 70);
    let d = ks(&a.depth_hist, &c.depth_hist);
    let n = a.primaries_stopped.min(c.primaries_stopped) as f64;
    assert!(d > 1.95 * (2.0 / n).sqrt(), "control not detected, KS {d}");
}

#[test]
fn periodic_lateral_wrapping_is_invisible() {
    // A 2 nm transverse period (the beam crosses it many times) against a
    // 100 um one (it never does): same depth profile, same energy budget
    // fractions.
    let ls = LindhardScharff::new();
    let mk = |w: f64| column(vec![si()], 40, NM, ([1, 1], w), Boundary::Periodic, |_| 0);
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.follow_recoils = false;
    cfg.seed = 5;
    let b = beam(8000);
    let (narrow, wide) = (mk(2.0 * NM), mk(1e-4));
    let a = summary(&Bca::new(b, &narrow, cfg, &ls, table()).unwrap(), 80);
    let w = summary(&Bca::new(b, &wide, cfg, &ls, table()).unwrap(), 80);
    assert_same_profile(&a, &w, "narrow vs wide period");
    assert_eq!(a.lateral, 0);
}

#[test]
fn vacuum_side_faces_conserve_energy_and_count_escapes() {
    let ls = LindhardScharff::new();
    // A 6 x 6 x 6 nm vacuum-bounded block of two materials: every face is
    // reachable, cascades on.
    let n = 6;
    let mut cells = Vec::new();
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                cells.push(u32::from((i + j + k) % 2 == 0 && i >= 2));
            }
        }
    }
    let grid = VoxelGrid::new(
        vec![si(), sio2()],
        [n; 3],
        [NM; 3],
        cells,
        [Boundary::Vacuum; 3],
    )
    .unwrap();
    let mut cfg = BcaConfig::new(5.0, 1.0);
    cfg.seed = 3;
    let mut b = beam(600);
    b.ion = Ion::new(14).unwrap();
    b.polar_rad = 0.9;
    let bca = Bca::new(b, &grid, cfg, &ls, table()).unwrap();
    let t = summary(&bca, 8);
    assert!(
        t.max_relative_residual < 1e-9,
        "{}",
        t.max_relative_residual
    );
    assert!(t.lateral + t.recoils_lateral > 0, "no lateral escapes");
    assert!(t.budget.lateral > 0.0);
    assert!(t.transmitted + t.recoils_transmitted > 0);
    assert_eq!(
        t.primaries_stopped + t.backscattered + t.transmitted + t.lateral,
        t.histories
    );
}

#[test]
fn entry_point_is_explicit() {
    let ls = LindhardScharff::new();
    let grid = column(
        vec![si()],
        10,
        NM,
        ([4, 4], 5.0 * NM),
        Boundary::Vacuum,
        |_| 0,
    );
    let cfg = BcaConfig::new(5.0, 2.0);
    let bca = Bca::new(beam(10), &grid, cfg, &ls, table()).unwrap();
    assert_eq!(grid.entry_point(), [0.0, 10.0 * NM, 10.0 * NM]);
    assert!(matches!(
        bca.with_entry_point([-NM, 0.0, 0.0]),
        Err(BcaError::InvalidBeam(_))
    ));
    let bca = Bca::new(beam(10), &grid, cfg, &ls, table()).unwrap();
    // The corner voxel instead of the centre.
    let bca = bca.with_entry_point([0.0, 0.5 * NM, 0.5 * NM]).unwrap();
    let t = summary(&bca, 8);
    assert_eq!(t.histories, 10);
}

fn bits(t: &SummaryTally) -> Vec<u64> {
    let b = &t.budget;
    let mut v = vec![
        t.histories,
        t.primaries_stopped,
        t.backscattered,
        t.transmitted,
        t.lateral,
        t.sputtered,
        t.recoils_transmitted,
        t.recoils_lateral,
        t.recoils,
        t.depth_sum.to_bits(),
        t.depth_sq_sum.to_bits(),
    ];
    for x in [
        b.incident,
        b.electronic_nonlocal,
        b.electronic_local,
        b.lattice,
        b.surface_barrier,
        b.backscattered,
        b.sputtered,
        b.transmitted,
        b.lateral,
        b.rest,
    ] {
        v.push(x.to_bits());
    }
    v.extend(&t.depth_hist);
    v
}

#[test]
fn voxel_runs_are_bit_identical_across_thread_counts() {
    let ls = LindhardScharff::new();
    let grid = column(
        vec![sio2(), si()],
        12,
        NM,
        ([3, 3], 3.0 * NM),
        Boundary::Vacuum,
        |i| u32::from(i >= 3),
    );
    let mut cfg = BcaConfig::new(5.0, 1.0);
    cfg.seed = 0xC0FFEE;
    cfg.chunk_size = 16;
    cfg.weak_collisions = 2;
    let mut b = beam(300);
    b.ion = Ion::new(18).unwrap();
    let bca = Bca::new(b, &grid, cfg, &ls, table()).unwrap();
    let run = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| summary(&bca, 12))
    };
    let r1 = run(1);
    assert!(r1.recoils > 0 && r1.lateral + r1.recoils_lateral > 0);
    for n in [2, 8] {
        assert_eq!(bits(&r1), bits(&run(n)), "differs on {n} threads");
    }
}
