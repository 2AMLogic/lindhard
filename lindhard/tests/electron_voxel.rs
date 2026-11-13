//! Electron transport on a voxel target: a laterally uniform voxel grid
//! against the layered stack it encodes, the per-voxel energy deposition
//! tally and its export, lateral escapes, and bit-identical results across
//! thread counts. The cross-section tables are computed here from simple
//! formulas (constant mean free paths, isotropic angles), not measured or
//! tabulated data.

use lindhard::electron::data::SamplingAxis;
use lindhard::electron::transport::{
    BoundaryModel, LayerTables, Primary, Transport, TransportConfig,
};
use lindhard::geometry::{Boundary, Stack, VoxelGrid};
use lindhard::material::Material;
use lindhard::tally::{
    Binning, ElectronReport, ElectronTallyConfig, FullElectronTally, VoxelReport, VoxelTallyGrid,
};

mod common;
use common::synthetic_table;
use common::with_threads;

const NM: f64 = 1e-9;
const E0: f64 = 500.0;

fn material(z: u8) -> Material {
    Material::from_atom_fractions(&[(z, 1.0)], None).unwrap()
}

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

/// Two materials with different scattering.
fn two_tables() -> Vec<LayerTables> {
    vec![
        scattering(12.0 * NM, 25.0 * NM, 0.1),
        scattering(9.0 * NM, 18.0 * NM, 0.1),
    ]
}

fn config() -> TransportConfig {
    TransportConfig::new(10.0)
}

/// 20 nm of material 0 on 30 nm of material 1, as a stack.
fn layered() -> Stack {
    Stack::new(
        vec![(material(14), 20.0 * NM), (material(29), 30.0 * NM)],
        None,
    )
    .unwrap()
}

/// The same layers as 10 voxels of 5 nm along `x`, laterally uniform and
/// infinite (periodic `y` and `z`).
fn uniform_voxels() -> VoxelGrid {
    let (nx, ny, nz) = (10, 2, 2);
    let mut cells = Vec::new();
    for _k in 0..nz {
        for _j in 0..ny {
            for i in 0..nx {
                cells.push(u32::from(i >= 4));
            }
        }
    }
    VoxelGrid::new(
        vec![material(14), material(29)],
        [nx, ny, nz],
        [5.0 * NM, 1.0e-6, 1.0e-6],
        cells,
        [Boundary::Vacuum, Boundary::Periodic, Boundary::Periodic],
    )
    .unwrap()
}

fn tally_config() -> ElectronTallyConfig {
    ElectronTallyConfig::new(
        Binning::new(0.0, E0, 25).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 6).unwrap(),
    )
}

/// Ten depth slabs of 5 nm and one huge lateral voxel: the laterally
/// integrated depth profile, whatever the target.
fn depth_grid() -> VoxelTallyGrid {
    VoxelTallyGrid::new([0.0, -0.5, -0.5], [5.0 * NM, 1.0, 1.0], [10, 1, 1]).unwrap()
}

fn tally<G: lindhard::geometry::Geometry>(
    t: &Transport<G>,
    grid: VoxelTallyGrid,
) -> FullElectronTally {
    FullElectronTally::new(t, tally_config())
        .unwrap()
        .with_voxel_grid(grid)
}

fn run_report<G: lindhard::geometry::Geometry>(
    t: &Transport<G>,
    n: u64,
    seed: u64,
) -> (ElectronReport, VoxelReport) {
    let proto = tally(t, depth_grid());
    let tally = t
        .run(seed, n, 64, &Primary::normal(E0), || proto.clone())
        .unwrap()
        .tally;
    (tally.report(), tally.voxel_report().unwrap())
}

/// Two counts from independent runs agree within `sigmas` Poisson standard
/// deviations of their difference.
fn counts_agree(name: &str, a: u64, b: u64, sigmas: f64) {
    let tol = sigmas * ((a + b) as f64).sqrt() + 1.0;
    assert!(
        (a as f64 - b as f64).abs() <= tol,
        "{name}: {a} vs {b} (tolerance {tol})"
    );
}

#[test]
fn uniform_voxel_target_matches_the_layered_stack_within_statistics() {
    let stack = Transport::new(layered(), two_tables(), config()).unwrap();
    let voxel = Transport::on_geometry(uniform_voxels(), two_tables(), config()).unwrap();
    let n = 20_000;
    // Independent seeds: the agreement is statistical, not an artefact of
    // identical random streams.
    let (rs, vs) = run_report(&stack, n, 0xA11CE);
    let (rv, vv) = run_report(&voxel, n, 0xB0B);

    // Yields and fates.
    counts_agree("stopped", rs.fates.stopped, rv.fates.stopped, 5.0);
    counts_agree("front", rs.front.count, rv.front.count, 5.0);
    counts_agree("back", rs.back.count, rv.back.count, 5.0);
    counts_agree(
        "fate front",
        rs.fates.escaped_front,
        rv.fates.escaped_front,
        5.0,
    );
    counts_agree(
        "fate back",
        rs.fates.escaped_back,
        rv.fates.escaped_back,
        5.0,
    );
    assert_eq!(rv.fates.escaped_side, 0, "periodic sides cannot leak");
    assert!(rv.side.is_none());
    assert!(rs.front.count > 100 && rs.back.count > 100);

    // Escape energy spectra, bin by bin.
    for (name, a, b) in [
        ("front", &rs.front, &rv.front),
        ("back", &rs.back, &rv.back),
    ] {
        for (i, (&x, &y)) in a
            .energy_histogram
            .counts
            .iter()
            .zip(&b.energy_histogram.counts)
            .enumerate()
        {
            counts_agree(&format!("{name} spectrum bin {i}"), x, y, 5.0);
        }
    }

    // Depth profile of the deposition: each slab's share of the total.
    let (ps, pv) = (vs.profile(0), vv.profile(0));
    let (ts, tv) = (
        vs.grand_totals().deposited_ev(),
        vv.grand_totals().deposited_ev(),
    );
    assert!(ts > 0.0 && tv > 0.0);
    for (i, (a, b)) in ps.iter().zip(&pv).enumerate() {
        let (fa, fb) = (a.deposited_ev() / ts, b.deposited_ev() / tv);
        assert!(
            (fa - fb).abs() <= 0.01,
            "slab {i}: deposited fraction {fa} vs {fb}"
        );
    }
    // Per-material totals are the per-layer totals of the stack.
    for m in 0..2 {
        let (a, b) = (rs.deposition.per_layer_ev[m], rv.deposition.per_layer_ev[m]);
        assert!(
            (a - b).abs() <= 0.02 * a.max(b),
            "material {m}: deposited {a} vs {b}"
        );
    }
}

#[test]
fn same_seed_gives_the_same_histories_on_stack_and_uniform_voxels() {
    // A voxel face inside a homogeneous block is not an event, so the random
    // stream is consumed exactly as in the layered run: the two runs follow
    // the same histories (up to rounding of the crossing points).
    let stack = Transport::new(layered(), two_tables(), config()).unwrap();
    let voxel = Transport::on_geometry(uniform_voxels(), two_tables(), config()).unwrap();
    let (rs, vs) = run_report(&stack, 3_000, 7);
    let (rv, vv) = run_report(&voxel, 3_000, 7);
    assert_eq!(rs.fates, rv.fates);
    assert_eq!(rs.front.count, rv.front.count);
    assert_eq!(rs.back.count, rv.back.count);
    assert!(
        (rs.budget.deposited_ev - rv.budget.deposited_ev).abs() <= 1e-9 * rs.budget.deposited_ev
    );
    assert_eq!(vs.voxels.len(), vv.voxels.len());
}

#[test]
fn voxel_deposition_sums_to_the_energy_budget() {
    let voxel = Transport::on_geometry(uniform_voxels(), two_tables(), config()).unwrap();
    let (r, v) = run_report(&voxel, 5_000, 11);
    assert!(r.budget.relative_imbalance < 1e-9);
    let total = v.grand_totals();
    let dep = r.budget.deposited_ev;
    assert!((total.deposited_ev() - dep).abs() <= 1e-9 * dep);
    assert!((total.electronic_ev - r.budget.inelastic_ev).abs() <= 1e-9 * dep);
    assert!((total.rest_ev - r.budget.residual_ev).abs() <= 1e-9 * dep);
    assert_eq!(v.histories, 5_000);
    // The shared export formats round-trip the electron deposition.
    let csv = v.to_csv();
    assert_eq!(VoxelReport::read_csv(&csv).unwrap(), v);
    let prov = lindhard::tally::ExportProvenance::new("electron", 11, &"0".repeat(64));
    let (back, p) = VoxelReport::read_export_csv(&v.to_export_csv(&prov)).unwrap();
    assert_eq!((back, p), (v.clone(), prov));
    let vtk = v.to_vtk();
    assert_eq!(VoxelReport::read_vtk(&vtk).unwrap(), v);
}

/// A 3 x 3 x 3 block target of two materials, 15 nm voxels, vacuum on every
/// face, so electrons leave through the sides as well.
fn small_block() -> VoxelGrid {
    let cells = (0..27u32)
        .map(|c| {
            let (i, j, k) = (c % 3, (c / 3) % 3, c / 9);
            (i + j + k) % 2
        })
        .collect();
    VoxelGrid::new(
        vec![material(14), material(29)],
        [3, 3, 3],
        [15.0 * NM; 3],
        cells,
        [Boundary::Vacuum; 3],
    )
    .unwrap()
}

#[test]
fn electrons_leave_a_finite_block_through_its_lateral_faces() {
    let t = Transport::on_geometry(small_block(), two_tables(), config()).unwrap();
    let grid = VoxelTallyGrid::from_voxel_grid(t.geometry());
    let proto = tally(&t, grid);
    let run = t
        .run(5, 4_000, 64, &Primary::normal(E0), || proto.clone())
        .unwrap();
    let r = run.tally.report();
    assert!(r.fates.escaped_side > 0 && r.fates.escaped_front > 0);
    let side = r.side.expect("side escapes are reported");
    assert!(side.count > 0 && r.budget.escaped_side_ev > 0.0);
    assert!(r.budget.relative_imbalance < 1e-9);
    assert_eq!(
        r.budget.escaped_ev,
        r.budget.escaped_front_ev + r.budget.escaped_back_ev + r.budget.escaped_side_ev
    );
    // Nothing is deposited outside the block, and the voxels hold it all.
    let v = run.tally.voxel_report().unwrap();
    assert_eq!(v.outside.deposited_ev(), 0.0);
    let dep = r.budget.deposited_ev;
    assert!((v.grand_totals().deposited_ev() - dep).abs() <= 1e-9 * dep);
}

fn run_threads(threads: usize) -> FullElectronTally {
    let t = Transport::on_geometry(small_block(), two_tables(), config()).unwrap();
    let grid = VoxelTallyGrid::from_voxel_grid(t.geometry());
    let proto = tally(&t, grid);
    with_threads(threads, || {
        t.run(0xDE7, 3_000, 32, &Primary::normal(E0), || proto.clone())
            .unwrap()
            .tally
    })
}

#[test]
fn voxel_results_are_bit_identical_across_thread_counts() {
    let one = run_threads(1);
    let two = run_threads(2);
    assert_eq!(one, two);
    assert_eq!(one.report(), two.report());
    assert_eq!(one.voxel_report(), two.voxel_report());
    assert_eq!(
        one.voxel_report().unwrap().to_csv(),
        two.voxel_report().unwrap().to_csv()
    );
}

#[test]
fn the_step_barrier_is_refused_off_the_stack() {
    let mut cfg = config();
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    assert!(Transport::on_geometry(uniform_voxels(), two_tables(), cfg).is_err());
}
