//! The fluence stepping loop (`ion::dynamic::DynamicRun`): low-fluence limit,
//! convergence under step refinement, index continuity, determinism at any
//! thread count, adaptive reject/retry and conservation.

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, BcaTally, Beam, SummaryTally};
use lindhard::ion::dynamic::{
    CompositionGrid, DynamicConfig, DynamicRun, InventoryTally, Relaxation, StepPolicy, StepRecord,
};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;

const NM: f64 = 1e-9;
const SI: u8 = 14;
const AS: u8 = 33;

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

/// Test parameters (not data): E_d = 15 eV for every element, E_s = 2 eV where
/// the element table has none.
fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(SI, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(SI, 15.0).unwrap();
    m
}

/// `n_slabs` slabs of `slab_nm` of Si, with a Si substrate if asked.
fn grid(n_slabs: usize, slab_nm: f64, substrate: bool) -> CompositionGrid {
    let layers = (0..n_slabs).map(|_| (si(), slab_nm * NM)).collect();
    let stack = Stack::new(layers, substrate.then(si)).unwrap();
    let mut g = CompositionGrid::from_stack(
        &stack,
        Relaxation::ideal_mixing_elemental(&[SI, AS]).unwrap(),
    )
    .unwrap();
    // The beam species enters the slabs by implantation.
    g.seed_energies(AS, Some(15.0), None, Some(2.0));
    g
}

fn beam(count: u64) -> Beam {
    Beam::normal(Ion::new(AS).unwrap(), 1000.0, count)
}

fn config(seed: u64) -> BcaConfig {
    let mut c = BcaConfig::new(5.0, 2.0);
    c.seed = seed;
    c
}

fn with_threads<R: Send>(n: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
        .install(f)
}

struct Outcome {
    grid: CompositionGrid,
    records: Vec<StepRecord>,
}

fn drive(
    g: CompositionGrid,
    count: u64,
    fluence_m2: f64,
    policy: StepPolicy,
    seed: u64,
) -> Outcome {
    let ls = LindhardScharff::new();
    let mut run = DynamicRun::new(
        g,
        beam(count),
        config(seed),
        &ls,
        table(),
        DynamicConfig { fluence_m2, policy },
    )
    .unwrap();
    let mut records = Vec::new();
    run.run_all(|r, _| records.push(r.clone())).unwrap();
    assert_eq!(run.ions_done(), count);
    Outcome {
        grid: run.grid().clone(),
        records,
    }
}

fn fixed(n: u64) -> StepPolicy {
    StepPolicy::Fixed { ions_per_step: n }
}

fn total_thickness(g: &CompositionGrid) -> f64 {
    g.thicknesses_m().iter().sum()
}

/// A bit-exact fingerprint of a finished run.
fn fingerprint(o: &Outcome) -> Vec<u64> {
    let mut v = Vec::new();
    for r in &o.records {
        v.extend([
            r.first_index,
            r.ions,
            u64::from(r.attempts),
            r.max_change.to_bits(),
        ]);
        v.push(r.yields.sputtered_total());
        v.push(r.yields.backscattered);
        v.push(r.removed_slabs.len() as u64);
    }
    v.extend(o.grid.thicknesses_m().iter().map(|t| t.to_bits()));
    for i in 0..o.grid.n_slabs() {
        for (z, a) in o.grid.inventory(i).unwrap() {
            v.push(u64::from(z));
            v.push(a.to_bits());
        }
    }
    v
}

#[test]
fn low_fluence_single_step_matches_the_static_engine() {
    let n = 400;
    let g = grid(4, 5.0, true);
    let stack = g.to_stack().unwrap();
    let ls = LindhardScharff::new();
    let bca = Bca::new(beam(n), &stack, config(7), &ls, table()).unwrap();
    let st = bca.run(|| SummaryTally::new(NM, 40)).unwrap();

    let o = drive(g, n, 1e14, fixed(n), 7);
    assert_eq!(o.records.len(), 1);
    let y = &o.records[0].yields;
    assert_eq!(o.records[0].first_index, 0);
    assert_eq!(y.histories, st.histories);
    assert_eq!(y.displaced, st.recoils);
    assert_eq!(
        y.primaries_in_slabs + y.primaries_in_substrate,
        st.primaries_stopped
    );
    assert_eq!(y.backscattered, st.backscattered);
    assert_eq!(y.transmitted, st.transmitted);
    assert_eq!(y.sputtered_total(), st.sputtered);
    // The target is practically unchanged at this fluence.
    let t0 = total_thickness(&grid(4, 5.0, true));
    assert!((total_thickness(&o.grid) / t0 - 1.0).abs() < 1e-3);
}

#[test]
fn step_ranges_continue_the_global_stream() {
    // Running 0..n in one go or as consecutive ranges draws the same histories.
    let g = grid(3, 6.0, true);
    let stack = g.to_stack().unwrap();
    let ls = LindhardScharff::new();
    let bca = Bca::new(beam(300), &stack, config(3), &ls, table()).unwrap();
    let whole = bca.run_range(0, 300, || InventoryTally::new(3)).unwrap();
    let mut parts = bca.run_range(0, 120, || InventoryTally::new(3)).unwrap();
    parts.merge(bca.run_range(120, 180, || InventoryTally::new(3)).unwrap());
    assert_eq!(whole, parts);
    // ...and a range that does not start at 0 is not the replay of 0..n.
    let tail = bca.run_range(120, 180, || InventoryTally::new(3)).unwrap();
    let replay = bca.run_range(0, 180, || InventoryTally::new(3)).unwrap();
    assert_ne!(tail, replay);
}

#[test]
fn bit_identical_at_1_2_and_8_threads() {
    let policies = [
        fixed(50),
        StepPolicy::Adaptive {
            max_ions_per_step: 200,
            min_ions_per_step: 5,
            max_change: 0.05,
        },
    ];
    for policy in policies {
        let reference = with_threads(1, || drive(grid(4, 4.0, true), 600, 4e19, policy, 11));
        let fp = fingerprint(&reference);
        for n in [2, 8] {
            let o = with_threads(n, || drive(grid(4, 4.0, true), 600, 4e19, policy, 11));
            assert_eq!(fp, fingerprint(&o), "{policy:?} differs on {n} threads");
        }
    }
}

#[test]
fn adaptive_retry_consumes_no_indices() {
    let o = drive(
        grid(4, 4.0, true),
        800,
        8e19,
        StepPolicy::Adaptive {
            max_ions_per_step: 400,
            min_ions_per_step: 2,
            max_change: 0.02,
        },
        5,
    );
    assert!(
        o.records.iter().any(|r| r.attempts > 1),
        "the bound was meant to reject at least one step: {:?}",
        o.records
            .iter()
            .map(|r| (r.ions, r.attempts))
            .collect::<Vec<_>>()
    );
    let mut next = 0;
    for r in &o.records {
        assert_eq!(r.first_index, next, "accepted steps must tile 0..count");
        next += r.ions;
        assert!(r.max_change <= 0.02 || r.ions == 2, "{r:?}");
    }
    assert_eq!(next, 800);
}

#[test]
fn inventory_is_conserved_per_the_event_conventions() {
    // No substrate and a target thick enough to stop everything: the net
    // change of atoms is the retained primaries minus the escaped recoils.
    let g0 = grid(4, 5.0, false);
    let atoms = |g: &CompositionGrid| g.total_inventory(SI) + g.total_inventory(AS);
    let before = atoms(&g0);
    let fluence = 2e19;
    let o = drive(g0, 2000, fluence, fixed(100), 2);
    let per_ion = fluence / 2000.0;
    let (mut kept, mut lost) = (0u64, 0u64);
    for r in &o.records {
        assert_eq!(r.clamped, 0);
        kept += r.yields.primaries_in_slabs;
        lost += r.yields.sputtered_total() + r.yields.recoils_transmitted;
    }
    let expected = (kept as f64 - lost as f64) * per_ion;
    let got = atoms(&o.grid) - before;
    assert!(
        (got - expected).abs() <= 1e-9 * before,
        "got {got}, expected {expected}"
    );
    assert!(o.grid.total_inventory(AS) > 0.0 && lost > 0);
}

#[test]
fn a_slab_that_empties_is_removed_and_the_run_continues() {
    // A film eroded through its front slabs: indices are slab indices of the
    // grid before each step, so the run must keep mapping layers to slabs
    // after removals.
    let g0 = grid(4, 1.0, false);
    let o = drive(
        g0,
        3000,
        4e20,
        StepPolicy::Adaptive {
            max_ions_per_step: 300,
            min_ions_per_step: 1,
            max_change: 0.25,
        },
        2,
    );
    assert!(
        o.records.iter().any(|r| !r.removed_slabs.is_empty()),
        "the 4 nm film was meant to lose a slab"
    );
    assert!(o.grid.n_slabs() < 4);
    // The grid still converts to a consistent stack of the remaining slabs.
    assert_eq!(o.grid.to_stack().unwrap().layers().len(), o.grid.n_slabs());
    let mut next = 0;
    for r in &o.records {
        assert_eq!(r.first_index, next);
        next += r.ions;
    }
    assert_eq!(next, 3000);
}

#[test]
fn halving_the_step_converges() {
    let count = 6000;
    let fluence = 1.5e19;
    let t = |k: u64| {
        let o = drive(grid(5, 3.0, true), count, fluence, fixed(count / k), 21);
        (total_thickness(&o.grid), o.grid.total_inventory(SI))
    };
    let reference = t(60);
    let d = |k| (t(k).0 - reference.0).abs();
    let (d2, d4, d8) = (d(2), d(4), d(8));
    eprintln!("thickness differences from the 60-step run: {d2:e} {d4:e} {d8:e}");
    assert!(d4 < d2 && 2.0 * d8 < d2, "{d2:e} {d4:e} {d8:e}");
}
