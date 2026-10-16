//! The sparse voxel tally (`lindhard::tally::voxel`): conservation against the
//! engine's 1D tallies on a laterally uniform target, bit-identical results at
//! any thread count, sparse storage, segment splitting, and the CSV and VTK
//! exporters (including the checked-in example files).

use std::path::PathBuf;
use std::sync::OnceLock;

use lindhard::geometry::{Boundary, Stack, VoxelGrid};
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, ElectronicChannel, EnergyBudget, Particle, SummaryTally,
};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::Binning;
use lindhard::tally::{
    IonTally, IonTallyConfig, VoxelCell, VoxelRecord, VoxelReport, VoxelTally, VoxelTallyError,
    VoxelTallyGrid,
};

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

/// Test parameters (not data): `E_d` = 15 eV, `E_s` = 2 eV.
fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    m
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

/// `nx` voxels of 1 nm along depth in front of vacuum, `2 x 2` transverse
/// voxels of 20 nm with periodic sides: laterally uniform.
fn column(nx: usize) -> VoxelGrid {
    VoxelGrid::new(
        vec![si()],
        [nx, 2, 2],
        [NM, 20.0 * NM, 20.0 * NM],
        vec![0; nx * 4],
        [Boundary::Vacuum, Boundary::Periodic, Boundary::Periodic],
    )
    .unwrap()
}

/// A voxel tally and a 1D summary fed by the same events of one run.
#[derive(Clone)]
struct Both {
    voxel: VoxelTally,
    summary: SummaryTally,
}

impl BcaTally for Both {
    fn begin_history(&mut self, i: u64) {
        self.voxel.begin_history(i);
        self.summary.begin_history(i);
    }
    fn electronic(&mut self, p: &Particle, from: [f64; 3], c: ElectronicChannel, e: f64) {
        self.voxel.electronic(p, from, c, e);
        self.summary.electronic(p, from, c, e);
    }
    fn lattice(&mut self, at: [f64; 3], l: usize, k: lindhard::ion::bca::LatticeDeposit, e: f64) {
        self.voxel.lattice(at, l, k, e);
        self.summary.lattice(at, l, k, e);
    }
    fn recoil(&mut self, r: &Particle) {
        self.voxel.recoil(r);
        self.summary.recoil(r);
    }
    fn stopped(&mut self, p: &Particle) {
        self.voxel.stopped(p);
        self.summary.stopped(p);
    }
    fn escaped(&mut self, p: &Particle, f: lindhard::ion::bca::Face) {
        self.voxel.escaped(p, f);
        self.summary.escaped(p, f);
    }
    fn end_history(&mut self, i: u64, b: &EnergyBudget) {
        self.voxel.end_history(i, b);
        self.summary.end_history(i, b);
    }
    fn merge(&mut self, o: Self) {
        self.voxel.merge(o.voxel);
        self.summary.merge(o.summary);
    }
}

fn run_both(grid: &VoxelGrid, count: u64, chunk: u64) -> Both {
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 33;
    cfg.chunk_size = chunk;
    let bca = Bca::new(beam(count), grid, cfg, &ls, table()).unwrap();
    let tally_grid = VoxelTallyGrid::from_voxel_grid(grid);
    bca.run(|| Both {
        voxel: VoxelTally::new(tally_grid),
        // Depth bins of one voxel; the last bin collects overflow.
        summary: SummaryTally::new(NM, grid.dims()[0]),
    })
    .unwrap()
}

fn rel_close(a: f64, b: f64, what: &str) {
    let scale = a.abs().max(b.abs()).max(1e-300);
    assert!((a - b).abs() <= 1e-9 * scale, "{what}: {a} vs {b}");
}

#[test]
fn totals_equal_the_one_dimensional_tallies() {
    let grid = column(40);
    let both = run_both(&grid, 3000, 64);
    let r = both.voxel.report();
    let s = &both.summary;
    let b = &s.budget;
    assert_eq!(r.histories, s.histories);
    assert!(s.primaries_stopped > 1000 && s.recoils > 1000);

    // Nothing leaves the grid on the periodic sides, and the beam enters at
    // x = 0: every event of a stopped particle is inside the grid, except
    // transmission through the back face, which deposits only inside.
    assert_eq!(r.outside.implanted, 0);
    let t = r.totals();
    assert_eq!(t.implanted, s.primaries_stopped);
    assert_eq!(t.displacements, s.recoils);
    rel_close(
        t.electronic_ev,
        b.electronic_nonlocal + b.electronic_local,
        "electronic",
    );
    rel_close(t.lattice_ev, b.lattice, "lattice");
    rel_close(t.rest_ev, b.rest, "rest");
    rel_close(
        r.grand_totals().deposited_ev(),
        b.electronic_nonlocal + b.electronic_local + b.lattice + b.rest,
        "deposited",
    );

    // The laterally integrated depth profile is the 1D depth histogram bin
    // for bin (voxel depth = histogram bin; the grid's last voxel and the
    // histogram's overflow bin both take everything deeper).
    let profile = r.profile(0);
    let implanted: Vec<u64> = profile.iter().map(|c| c.implanted).collect();
    assert_eq!(implanted, s.depth_hist);
    // The profile is genuinely non-trivial and the tally is sparse.
    assert!(implanted.iter().filter(|&&n| n > 0).count() > 5);
    assert!(r.voxels.len() <= grid.dims().iter().product::<usize>());
}

#[test]
fn matches_the_ion_tally_of_the_layer_engine() {
    // The same laterally uniform target as a stack of one semi-infinite
    // layer: the voxel tally totals agree with `IonTally`, within the
    // statistics of independent runs (the engines draw equivalent but not
    // identical histories).
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 44;
    let n = 6000;
    let stack = Stack::new(vec![], Some(si())).unwrap();
    let ion_cfg = IonTallyConfig {
        depth: Binning::new(0.0, 200.0 * NM, 40).unwrap(),
        lateral: Binning::new(-100.0 * NM, 100.0 * NM, 20).unwrap(),
        radial: Binning::new(0.0, 100.0 * NM, 20).unwrap(),
        escape_energy: Binning::new(0.0, 2.0e3, 20).unwrap(),
        escape_polar: Binning::new(0.0, 1.6, 8).unwrap(),
    };
    let bca = Bca::new(beam(n), &stack, cfg, &ls, table()).unwrap();
    let proto = IonTally::new(&stack, &bca.species_z(), ion_cfg).unwrap();
    let ion = bca.run(|| proto.clone()).unwrap().report(false);

    // Deep voxels (2 um) so nothing is transmitted.
    let grid = column(2000);
    let both = run_both(&grid, n, 64);
    let t = both.voxel.report().totals();
    let stopped = ion.range.stopped as f64;
    assert!((t.implanted as f64 - stopped).abs() < 0.02 * stopped);
    let disp = ion.damage.cascade.displacements as f64;
    assert!((t.displacements as f64 - disp).abs() < 0.05 * disp);
}

#[test]
fn results_are_bit_identical_across_thread_counts() {
    let grid = column(12);
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| run_both(&grid, 400, 16).voxel.report())
    };
    let r1 = run(1);
    assert!(r1.voxels.len() > 10 && r1.totals().displacements > 0);
    for n in [2, 8] {
        let rn = run(n);
        assert_eq!(r1.voxels.len(), rn.voxels.len());
        for (a, b) in r1.voxels.iter().zip(&rn.voxels) {
            assert_eq!(a.index, b.index);
            assert_eq!(a.cell.implanted, b.cell.implanted);
            assert_eq!(a.cell.displacements, b.cell.displacements);
            for (x, y) in [
                (a.cell.electronic_ev, b.cell.electronic_ev),
                (a.cell.lattice_ev, b.cell.lattice_ev),
                (a.cell.rest_ev, b.cell.rest_ev),
            ] {
                assert_eq!(x.to_bits(), y.to_bits(), "differs on {n} threads");
            }
        }
        assert_eq!(r1.to_csv(), rn.to_csv());
    }
}

fn particle(pos: [f64; 3], energy_ev: f64, generation: u32) -> Particle {
    Particle {
        species: 0,
        z: 14,
        mass_amu: 28.0,
        energy_ev,
        pos,
        dir: [1.0, 0.0, 0.0],
        layer: 0,
        origin_layer: 0,
        generation,
    }
}

#[test]
fn empty_voxels_cost_no_memory() {
    // 10^12 voxels: a dense array would need terabytes.
    let grid = VoxelTallyGrid::new([0.0; 3], [NM; 3], [10_000; 3]).unwrap();
    assert_eq!(grid.n_voxels(), 1_000_000_000_000);
    let mut t = VoxelTally::new(grid);
    t.stopped(&particle([5.5 * NM, 1.5 * NM, 9_999.5 * NM], 3.0, 0));
    t.stopped(&particle([5.5 * NM, 1.5 * NM, 9_999.5 * NM], 3.0, 0));
    assert_eq!(t.occupied(), 1);
    let r = t.report();
    assert_eq!(r.voxels.len(), 1);
    assert_eq!(r.voxels[0].index, [5, 1, 9_999]);
    assert_eq!(r.voxels[0].cell.implanted, 2);
    assert_eq!(
        grid.unflatten(grid.flat_index([5, 1, 9_999])),
        [5, 1, 9_999]
    );
}

#[test]
fn flight_segments_split_between_voxels_and_conserve_energy() {
    let grid = VoxelTallyGrid::new([0.0; 3], [NM; 3], [4, 4, 1]).unwrap();
    let mut t = VoxelTally::new(grid);
    // Along x through voxels 0..3 at 0.5 nm: 0.25, 1, 1, 0.5 nm.
    let p = particle([2.5 * NM, 0.5 * NM, 0.5 * NM], 100.0, 0);
    t.electronic(
        &p,
        [0.25 * NM, 0.5 * NM, 0.5 * NM],
        ElectronicChannel::NonLocal,
        9.0,
    );
    let r = t.report();
    let by_x: Vec<(usize, f64)> = r
        .voxels
        .iter()
        .map(|v| (v.index[0], v.cell.electronic_ev))
        .collect();
    assert_eq!(by_x.len(), 3);
    // 2.25 nm total: 0.75/2.25, 1/2.25 and 0.5/2.25 of 9 eV.
    for (i, want) in [(0, 3.0), (1, 4.0), (2, 2.0)] {
        let got = by_x.iter().find(|(x, _)| *x == i).unwrap().1;
        assert!((got - want).abs() < 1e-12, "voxel {i}: {got} vs {want}");
    }

    // A diagonal segment through a corner, and one that leaves the grid: the
    // part beyond the grid goes to `outside`, and the sum is conserved.
    let mut t = VoxelTally::new(grid);
    let p = particle([5.0 * NM, 5.0 * NM, 0.5 * NM], 100.0, 0);
    t.electronic(&p, [-NM, -NM, 0.5 * NM], ElectronicChannel::NonLocal, 12.0);
    let r = t.report();
    rel_close(r.grand_totals().electronic_ev, 12.0, "conserved");
    // 6 nm diagonal, of which the middle 4 nm (x = 0..4) is inside.
    rel_close(r.totals().electronic_ev, 8.0, "inside");
    rel_close(r.outside.electronic_ev, 4.0, "outside");
    assert_eq!(r.voxels.len(), 4);
    assert!(r.voxels.iter().all(|v| v.index[0] == v.index[1]));
}

#[test]
fn point_events_go_to_the_voxel_that_owns_them() {
    let grid = VoxelTallyGrid::new([0.0; 3], [NM; 3], [2, 2, 2]).unwrap();
    // Half-open voxels: an interior face belongs to the higher index; the
    // far end of the grid, a negative coordinate and NaN are outside.
    assert_eq!(grid.locate([NM, 0.0, 0.0]), Some([1, 0, 0]));
    assert_eq!(grid.locate([2.0 * NM, 0.0, 0.0]), None);
    assert_eq!(grid.locate([-1e-30, 0.0, 0.0]), None);
    assert_eq!(grid.locate([f64::NAN, 0.0, 0.0]), None);
    let mut t = VoxelTally::new(grid);
    t.recoil(&particle([1.5 * NM, 0.5 * NM, 0.5 * NM], 20.0, 1));
    t.recoil(&particle([3.0 * NM, 0.5 * NM, 0.5 * NM], 20.0, 1));
    let r = t.report();
    assert_eq!(r.totals().displacements, 1);
    assert_eq!(r.outside.displacements, 1);
}

#[test]
fn invalid_grids_are_rejected() {
    assert_eq!(
        VoxelTallyGrid::new([0.0; 3], [NM; 3], [1, 0, 1]),
        Err(VoxelTallyError::ZeroDimension { axis: 1 })
    );
    assert!(matches!(
        VoxelTallyGrid::new([0.0; 3], [NM, -NM, NM], [1; 3]),
        Err(VoxelTallyError::InvalidSpacing { axis: 1, .. })
    ));
    assert!(matches!(
        VoxelTallyGrid::new([0.0, f64::NAN, 0.0], [NM; 3], [1; 3]),
        Err(VoxelTallyError::InvalidOrigin { axis: 1, .. })
    ));
    assert!(matches!(
        VoxelTallyGrid::new([0.0; 3], [NM; 3], [usize::MAX; 3]),
        Err(VoxelTallyError::TooManyVoxels { .. })
    ));
}

#[test]
#[should_panic(expected = "different grids")]
fn merging_tallies_of_different_grids_panics() {
    let a = VoxelTallyGrid::new([0.0; 3], [NM; 3], [2; 3]).unwrap();
    let b = VoxelTallyGrid::new([0.0; 3], [2.0 * NM; 3], [2; 3]).unwrap();
    VoxelTally::new(a).merge(VoxelTally::new(b));
}

/// A fixed, engine-independent report: the checked-in example files are its
/// exports, so they do not change with the physics.
fn fixture() -> VoxelReport {
    let grid = VoxelTallyGrid::new([0.0; 3], [5.0 * NM; 3], [4, 2, 2]).unwrap();
    let rec = |index: [usize; 3], implanted, displacements, e, l, r| VoxelRecord {
        index,
        cell: VoxelCell {
            implanted,
            displacements,
            electronic_ev: e,
            lattice_ev: l,
            rest_ev: r,
        },
    };
    VoxelReport {
        grid,
        histories: 100,
        voxels: vec![
            rec([0, 0, 0], 0, 12, 1800.5, 96.25, 0.0),
            rec([1, 0, 0], 3, 41, 5210.125, 310.5, 4.75),
            rec([2, 0, 0], 9, 17, 2410.0, 120.0, 21.5),
            rec([1, 1, 0], 1, 5, 640.25, 33.0, 2.5),
            rec([2, 1, 1], 0, 2, 0.1, 0.2, 0.0),
        ],
        outside: VoxelCell {
            implanted: 0,
            displacements: 1,
            electronic_ev: 3.5,
            lattice_ev: 0.0,
            rest_ev: 0.0,
        },
    }
}

fn example(name: &str) -> String {
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "..",
        "examples",
        "voxel-tally",
        name,
    ]
    .iter()
    .collect();
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

#[test]
fn the_checked_in_examples_are_the_exports_and_read_back() {
    let r = fixture();
    // The files are the exports of the fixture (to regenerate them, write
    // `fixture().to_csv()` and `.to_vtk()` over the files).
    assert_eq!(example("example.csv"), r.to_csv());
    assert_eq!(example("example.vtk"), r.to_vtk());

    // Read back the checked-in files themselves.
    assert_eq!(VoxelReport::read_csv(&example("example.csv")).unwrap(), r);
    let from_vtk = VoxelReport::read_vtk(&example("example.vtk")).unwrap();
    assert_eq!(from_vtk.voxels, r.voxels);
    assert_eq!(from_vtk.grid, r.grid);
    assert_eq!(from_vtk.histories, r.histories);
}

#[test]
fn engine_reports_round_trip_through_csv_and_vtk() {
    let report = run_both(&column(12), 300, 16).voxel.report();
    assert!(report.voxels.len() > 10);
    // Floats survive bit for bit.
    assert_eq!(VoxelReport::read_csv(&report.to_csv()).unwrap(), report);
    let vtk = VoxelReport::read_vtk(&report.to_vtk()).unwrap();
    assert_eq!(vtk.voxels, report.voxels);
    assert_eq!(vtk.grid, report.grid);
    assert_eq!(vtk.histories, report.histories);
}

#[test]
fn vtk_file_has_the_legacy_ascii_layout() {
    let vtk = fixture().to_vtk();
    let lines: Vec<&str> = vtk.lines().collect();
    assert_eq!(lines[0], "# vtk DataFile Version 3.0");
    assert!(lines[1].len() <= 255);
    assert_eq!(&lines[2..4], ["ASCII", "DATASET UNSTRUCTURED_GRID"]);
    assert_eq!(lines[4], "POINTS 40 double");
    assert!(vtk.contains("CELLS 5 45\n8 0 1 2 3 4 5 6 7\n"));
    assert!(vtk.contains("CELL_TYPES 5\n11\n"));
    assert!(vtk
        .contains("CELL_DATA 5\nSCALARS implanted unsigned_long 1\nLOOKUP_TABLE default\n0\n3\n"));
    // First cell: the voxel [0, 0, 0], corners in VTK_VOXEL order.
    assert_eq!(lines[5], "0e0 0e0 0e0");
    assert_eq!(lines[6], "5e-9 0e0 0e0");
    assert_eq!(lines[7], "0e0 5e-9 0e0");
    assert_eq!(lines[9], "0e0 0e0 5e-9");
}

#[test]
fn malformed_files_are_rejected_with_a_line() {
    let csv = fixture().to_csv();
    assert!(VoxelReport::read_csv("").is_err());
    // Out-of-order rows.
    let swapped = csv.replacen("1,0,0,", "9,0,0,", 1);
    assert!(matches!(
        VoxelReport::read_csv(&swapped),
        Err(VoxelTallyError::Parse { .. })
    ));
    let truncated: String = csv.lines().take(9).collect::<Vec<_>>().join("\n") + "\n1,2";
    assert!(matches!(
        VoxelReport::read_csv(&truncated),
        Err(VoxelTallyError::Parse { line: 10, .. })
    ));
    let vtk = fixture().to_vtk();
    assert!(VoxelReport::read_vtk("not vtk").is_err());
    assert!(VoxelReport::read_vtk(&vtk[..vtk.len() / 2]).is_err());
}
