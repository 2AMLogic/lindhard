//! Composition grid and volume relaxation (`ion::dynamic`).

use lindhard::geometry::Stack;
use lindhard::ion::dynamic::{
    atomic_volume_from_density, CompositionGrid, DynamicError, InventoryDelta, Relaxation,
};
use lindhard::material::Material;

const SI: u8 = 14;
const O: u8 = 8;
const AR: u8 = 18;

fn si() -> Material {
    Material::from_atom_fractions(&[(SI, 1.0)], None).unwrap()
}

fn d(slab: usize, z: u8, delta: f64) -> InventoryDelta {
    InventoryDelta {
        slab,
        z,
        delta_atoms_m2: delta,
    }
}

fn rel_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs())
}

fn ideal_si_o() -> Relaxation {
    // Explicit volumes (not tied to any table): Si from its solid, O given.
    Relaxation::ideal_mixing(&[(SI, atomic_volume_from_density(SI).unwrap()), (O, 2.0e-29)])
        .unwrap()
}

fn two_slab_grid(r: Relaxation) -> CompositionGrid {
    let s = Stack::new(vec![(si(), 10e-9), (si(), 20e-9)], Some(si())).unwrap();
    CompositionGrid::from_stack(&s, r).unwrap()
}

#[test]
fn initial_inventory_is_density_times_thickness() {
    let g = two_slab_grid(ideal_si_o());
    let n = si().number_density_of(SI).unwrap();
    let inv = g.inventory(0).unwrap();
    assert_eq!(inv.len(), 1);
    assert!(rel_eq(inv[0].1, n * 10e-9, 1e-12));
    assert!(rel_eq(g.total_inventory(SI), n * 30e-9, 1e-12));
}

#[test]
fn zero_step_is_identity_when_consistent() {
    let mut g = two_slab_grid(ideal_si_o());
    let before_t = g.thicknesses_m();
    let before_inv = g.inventory(0).unwrap();
    let out = g.apply(&[]).unwrap();
    assert!(out.removed_slabs.is_empty());
    for (a, b) in before_t.iter().zip(g.thicknesses_m()) {
        assert!(rel_eq(*a, b, 1e-12));
    }
    assert_eq!(g.inventory(0).unwrap(), before_inv);
    // Zero deltas explicitly: same.
    g.apply(&[d(0, SI, 0.0)]).unwrap();
    assert!(rel_eq(before_t[0], g.thicknesses_m()[0], 1e-12));
}

#[test]
fn relaxation_keeps_pure_element_at_its_density() {
    let rho = si().mass_density();
    for r in [
        ideal_si_o(),
        Relaxation::fixed_number_density(si().atom_number_density()).unwrap(),
    ] {
        let mut g = two_slab_grid(r);
        // Add 10 % more Si to slab 0: thickness scales, density stays.
        let a0 = g.inventory(0).unwrap()[0].1;
        g.apply(&[d(0, SI, 0.1 * a0)]).unwrap();
        let stack = g.to_stack().unwrap();
        let l0 = &stack.layers()[0];
        assert!(rel_eq(l0.thickness_m(), 11e-9, 1e-9));
        assert!(rel_eq(l0.material().mass_density(), rho, 1e-9));
        let n = l0.material().number_density_of(SI).unwrap();
        assert!(rel_eq(n * l0.thickness_m(), 1.1 * a0, 1e-12));
    }
}

#[test]
fn binary_mixture_hand_computed() {
    let mut g = two_slab_grid(ideal_si_o());
    let vsi = atomic_volume_from_density(SI).unwrap();
    let a_si = g.inventory(0).unwrap()[0].1;
    let a_o = 0.5 * a_si;
    g.apply(&[d(0, O, a_o)]).unwrap();
    let stack = g.to_stack().unwrap();
    let l = &stack.layers()[0];
    let t = a_si * vsi + a_o * 2.0e-29;
    assert!(rel_eq(l.thickness_m(), t, 1e-12));
    let m = l.material();
    let fr = m.atom_fractions();
    assert_eq!(fr[0].0, O);
    assert!(rel_eq(fr[0].1, a_o / (a_o + a_si), 1e-12));
    assert!(rel_eq(fr[1].1, a_si / (a_o + a_si), 1e-12));
    // n_i * thickness recovers the inventories.
    assert!(rel_eq(m.number_density_of(SI).unwrap() * t, a_si, 1e-12));
    assert!(rel_eq(m.number_density_of(O).unwrap() * t, a_o, 1e-12));
    // Mass density: mass per area over thickness.
    let na = lindhard::constants::AVOGADRO;
    let mass = (a_si * 28.085 + a_o * 15.999) * 1e-3 / na;
    assert!(rel_eq(m.mass_density(), mass / t, 1e-3));
}

#[test]
fn fixed_number_density_option() {
    let n_mix = 6.0e28;
    let mut g = two_slab_grid(Relaxation::fixed_number_density(n_mix).unwrap());
    let a = g.inventory(0).unwrap()[0].1;
    g.apply(&[d(0, AR, 0.25 * a)]).unwrap();
    let t = 1.25 * a / n_mix;
    assert!(rel_eq(g.thicknesses_m()[0], t, 1e-12));
    let stack = g.to_stack().unwrap();
    assert!(rel_eq(
        stack.layers()[0].material().atom_number_density(),
        n_mix,
        1e-9
    ));
}

#[test]
fn pure_implant_adds_exactly_the_dose() {
    let mut g = two_slab_grid(Relaxation::fixed_number_density(5.0e28).unwrap());
    let dose = 3.0e19;
    g.apply(&[d(1, AR, dose)]).unwrap();
    assert_eq!(g.total_inventory(AR), dose);
    g.apply(&[d(1, AR, dose)]).unwrap();
    assert_eq!(g.total_inventory(AR), 2.0 * dose);
}

#[test]
fn sputter_and_target_transmission_remove_exactly() {
    let mut g = two_slab_grid(Relaxation::fixed_number_density(5.0e28).unwrap());
    let before = g.total_inventory(SI);
    let loss = 1.0e18;
    g.apply(&[d(0, SI, -loss)]).unwrap();
    assert_eq!(g.total_inventory(SI), before - loss);
}

#[test]
fn relocation_preserves_global_total() {
    let mut g = two_slab_grid(Relaxation::fixed_number_density(5.0e28).unwrap());
    // Beam Si implanted into slab 1 as a retained primary, and one target Si
    // recoil relocated from slab 0 to slab 1: roles are not conflated.
    let before = g.total_inventory(SI);
    let a0 = g.inventory(0).unwrap()[0].1;
    let a1 = g.inventory(1).unwrap()[0].1;
    let moved = 1.0e17;
    let beam = 2.0e17;
    g.apply(&[d(0, SI, -moved), d(1, SI, moved), d(1, SI, beam)])
        .unwrap();
    assert_eq!(g.inventory(0).unwrap()[0].1, a0 - moved);
    assert_eq!(g.inventory(1).unwrap()[0].1, a1 + moved + beam);
    assert!(rel_eq(g.total_inventory(SI), before + beam, 1e-14));
}

#[test]
fn invalid_updates_are_atomic() {
    let mut g = two_slab_grid(ideal_si_o());
    let snapshot = g.clone();
    let a = g.inventory(0).unwrap()[0].1;

    // Valid first delta, bad second: nothing applied.
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let e = g.apply(&[d(0, SI, 1.0e10), d(1, SI, bad)]).unwrap_err();
        assert!(matches!(e, DynamicError::NonFiniteDelta { .. }));
        assert_eq!(g, snapshot);
    }
    let e = g.apply(&[d(0, SI, 1.0e10), d(1, SI, -1e30)]).unwrap_err();
    assert!(matches!(e, DynamicError::NegativeInventory { slab: 1, .. }));
    assert_eq!(g, snapshot);
    let e = g.apply(&[d(0, O, -1.0)]).unwrap_err();
    assert!(matches!(e, DynamicError::NegativeInventory { .. }));
    assert_eq!(g, snapshot);
    // Missing volume (Ar has none in this convention).
    let e = g.apply(&[d(0, SI, 1.0), d(1, AR, 1.0e18)]).unwrap_err();
    assert_eq!(e, DynamicError::MissingVolume(AR));
    assert_eq!(g, snapshot);
    let e = g.apply(&[d(5, SI, 1.0)]).unwrap_err();
    assert!(matches!(e, DynamicError::SlabOutOfRange { .. }));
    let e = g.apply(&[d(0, 200, 1.0)]).unwrap_err();
    assert_eq!(e, DynamicError::UnknownAtomicNumber(200));
    // Overflow.
    let e = g
        .apply(&[d(0, SI, f64::MAX), d(0, SI, f64::MAX)])
        .unwrap_err();
    assert!(matches!(e, DynamicError::Overflow { .. }));
    assert_eq!(g, snapshot);
    let _ = a;
}

#[test]
fn bad_relaxation_parameters_rejected() {
    for v in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            Relaxation::ideal_mixing(&[(SI, v)]),
            Err(DynamicError::InvalidParameter(_))
        ));
        assert!(matches!(
            Relaxation::fixed_number_density(v),
            Err(DynamicError::InvalidParameter(_))
        ));
    }
    // A gas has no elemental solid density to derive a volume from.
    assert_eq!(
        atomic_volume_from_density(AR),
        Err(DynamicError::NoElementalDensity(AR))
    );
}

#[test]
fn depleted_slabs_are_removed() {
    let mut g = two_slab_grid(Relaxation::fixed_number_density(5.0e28).unwrap());
    let a0 = g.inventory(0).unwrap()[0].1;
    let out = g.apply(&[d(0, SI, -a0)]).unwrap();
    assert_eq!(out.removed_slabs, vec![0]);
    assert_eq!(g.n_slabs(), 1);
    let stack = g.to_stack().unwrap();
    assert_eq!(stack.layers()[0].front_m(), 0.0);

    // All slabs depleted: error without substrate, atomic.
    let s = Stack::new(vec![(si(), 10e-9)], None).unwrap();
    let mut f =
        CompositionGrid::from_stack(&s, Relaxation::fixed_number_density(5.0e28).unwrap()).unwrap();
    let snap = f.clone();
    let a = f.inventory(0).unwrap()[0].1;
    assert_eq!(f.apply(&[d(0, SI, -a)]), Err(DynamicError::EmptyTarget));
    assert_eq!(f, snap);

    // With a substrate the grid becomes the bare substrate.
    let mut h = two_slab_grid(Relaxation::fixed_number_density(5.0e28).unwrap());
    let a0 = h.inventory(0).unwrap()[0].1;
    let a1 = h.inventory(1).unwrap()[0].1;
    let out = h.apply(&[d(0, SI, -a0), d(1, SI, -a1)]).unwrap();
    assert_eq!(out.removed_slabs, vec![0, 1]);
    let st = h.to_stack().unwrap();
    assert_eq!(st.layers().len(), 1);
    assert!(st.layers()[0].thickness_m().is_infinite());
}

#[test]
fn rebuilt_stack_is_contiguous_with_substrate_and_energies() {
    let mut m = si();
    m.set_displacement_energy_ev(SI, 31.5).unwrap();
    m.set_surface_binding_energy_ev(SI, 4.2).unwrap();
    m.set_lattice_binding_energy_ev(SI, 1.5).unwrap();
    let s = Stack::new(vec![(m, 10e-9), (si(), 20e-9)], Some(si())).unwrap();
    let mut g =
        CompositionGrid::from_stack(&s, Relaxation::fixed_number_density(5.0e28).unwrap()).unwrap();
    let a = g.inventory(0).unwrap()[0].1;
    g.apply(&[d(0, SI, 0.3 * a), d(0, AR, 1.0e18)]).unwrap();
    let st = g.to_stack().unwrap();
    assert_eq!(st.layers().len(), 3);
    let mut x = 0.0;
    for l in &st.layers()[..2] {
        assert_eq!(l.front_m(), x);
        assert!(l.thickness_m() > 0.0 && l.thickness_m().is_finite());
        x = l.back_m();
    }
    assert!(st.layers()[2].back_m().is_infinite());
    let m0 = st.layers()[0].material();
    assert_eq!(m0.displacement_energy_ev(SI).unwrap(), 31.5);
    assert_eq!(m0.surface_binding_energy_ev(SI).unwrap(), 4.2);
    assert_eq!(m0.lattice_binding_energy_ev(SI).unwrap(), 1.5);
    // The substrate is untouched by updates.
    assert_eq!(st.layers()[2].material(), &si());
}

#[test]
fn updates_are_deterministic() {
    let run = || {
        let mut g = two_slab_grid(ideal_si_o());
        let a = g.inventory(0).unwrap()[0].1;
        g.apply(&[d(0, O, 0.3 * a), d(1, O, 0.1 * a), d(0, SI, -0.05 * a)])
            .unwrap();
        g
    };
    assert_eq!(run(), run());
}
