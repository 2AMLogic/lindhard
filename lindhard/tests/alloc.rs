//! Steady-state collision loop: no heap allocation.
//!
//! A counting global allocator (per thread, so the test harness does not add
//! noise) counts allocations while histories run in reused
//! [`HistoryBuffers`]. The histories are first run once to let the buffers
//! grow to the largest cascade they meet; the same histories (same streams)
//! are then run again and must allocate nothing, whatever the number of
//! collisions.

// The one place unsafe is needed: a counting `GlobalAlloc` (workspace lint is `deny`).
#![allow(unsafe_code)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, ElectronicLoss, EnergyBudget, HistoryBuffers, LatticeDeposit,
    MeanFreePath, Particle,
};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::rng::stream;

thread_local! {
    static ALLOCS: Cell<u64> = const { Cell::new(0) };
}

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let _ = ALLOCS.try_with(|c| c.set(c.get() + 1));
        System.alloc(l)
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l)
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        let _ = ALLOCS.try_with(|c| c.set(c.get() + 1));
        System.realloc(p, l, n)
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        let _ = ALLOCS.try_with(|c| c.set(c.get() + 1));
        System.alloc_zeroed(l)
    }
}

#[global_allocator]
static A: Counting = Counting;

fn allocs() -> u64 {
    ALLOCS.with(Cell::get)
}

/// Counts events; holds no heap data.
#[derive(Default)]
struct Counts {
    collisions: u64,
    recoils: u64,
    energy: f64,
}

impl BcaTally for Counts {
    fn lattice(&mut self, _: [f64; 3], _: usize, _: LatticeDeposit, e: f64) {
        self.collisions += 1;
        self.energy += e;
    }
    fn recoil(&mut self, _: &Particle) {
        self.recoils += 1;
    }
    fn end_history(&mut self, _: u64, b: &EnergyBudget) {
        self.energy += b.rest;
    }
    fn merge(&mut self, o: Self) {
        self.collisions += o.collisions;
        self.recoils += o.recoils;
        self.energy += o.energy;
    }
}

fn material(z: u8) -> Material {
    let mut m = Material::from_atom_fractions(&[(z, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(z, 15.0).unwrap();
    if m.surface_binding_energy_ev(z).is_err() {
        m.set_surface_binding_energy_ev(z, 2.0).unwrap();
    }
    m
}

fn table() -> ScatteringTable {
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
}

fn check(z_ion: u8, energy_ev: f64, z_target: u8, edit: impl Fn(&mut BcaConfig)) {
    const N: u64 = 150;
    let stack = Stack::semi_infinite(material(z_target));
    let table = table();
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 7;
    edit(&mut cfg);
    let beam = Beam {
        ion: Ion::new(z_ion).unwrap(),
        energy_ev,
        polar_rad: 7f64.to_radians(),
        azimuth_rad: 0.0,
        count: N,
    };
    let bca = Bca::new(beam, &stack, cfg, &ls, &table).unwrap();
    let mut buffers = HistoryBuffers::new();
    let run = |buffers: &mut HistoryBuffers| {
        let mut t = Counts::default();
        for i in 0..N {
            let mut rng = stream(cfg.seed, i);
            bca.history_in(buffers, &mut t, &mut rng, i).unwrap();
        }
        t
    };
    let warm = run(&mut buffers);
    let before = allocs();
    let hot = run(&mut buffers);
    let n = allocs() - before;
    assert!(hot.collisions > 10_000, "too few collisions to mean much");
    assert_eq!(hot.collisions, warm.collisions);
    assert_eq!(hot.energy.to_bits(), warm.energy.to_bits());
    assert_eq!(
        n, 0,
        "{n} allocations in {} collisions ({} recoils)",
        hot.collisions, hot.recoils
    );
}

#[test]
fn constant_free_path_cascade_allocates_nothing() {
    check(18, 1.0e3, 29, |_| {});
}

#[test]
fn energy_dependent_path_equipartition_allocates_nothing() {
    check(5, 5.0e3, 14, |c| {
        c.mean_free_path = MeanFreePath::EnergyDependent {
            min_cm_angle_rad: 0.01,
        };
        c.electronic = ElectronicLoss::EquipartitionLsOr;
    });
}

#[test]
fn weak_collisions_allocate_nothing() {
    check(18, 1.0e3, 29, |c| c.weak_collisions = 3);
}

/// The counter itself works: `Bca::history` makes fresh buffers each call, so
/// it must be seen allocating.
#[test]
fn counter_sees_allocations_of_fresh_buffers() {
    let stack = Stack::semi_infinite(material(29));
    let table = table();
    let ls = LindhardScharff::new();
    let cfg = BcaConfig::new(5.0, 2.0);
    let beam = Beam::normal(Ion::new(18).unwrap(), 1.0e3, 10);
    let bca = Bca::new(beam, &stack, cfg, &ls, &table).unwrap();
    let before = allocs();
    for i in 0..10 {
        let mut t = Counts::default();
        bca.history(&mut t, &mut stream(1, i), i).unwrap();
    }
    assert!(allocs() - before >= 10);
}
