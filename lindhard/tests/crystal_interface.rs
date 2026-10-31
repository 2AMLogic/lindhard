//! The amorphous-to-crystal interface of a mixed target (issue #514, step 23
//! of the crystalline-implant epic): an ion that crosses an amorphous
//! overlayer into the crystal continues exactly as one started at the crystal
//! surface, and a crystal-in-stack run reports recoils, lattice deposits and
//! escapes with the layer they belong to. See `docs/screen-oxide.md`.

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, CrystalTarget, ElectronicChannel, EnergyBudget, Face,
    LatticeDeposit, Particle, Thermal,
};
use lindhard::ion::crystal::debye::THETA_D_SI;
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;

const NM: f64 = 1e-9;
const SEED: u64 = 7;
/// Leading events of each history that are compared. Collisions in a crystal
/// amplify rounding differences of the two geometries exponentially, so only
/// the start of a history can be compared at a tolerance that still means
/// "identical"; the entry state is the part the interface decides.
const PREFIX: usize = 6;

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

fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

fn sio2() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
    for z in [14, 8] {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}

/// An overlayer so dilute that no collision or measurable electronic loss
/// happens across it: it only supplies an interface.
fn transparent() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], Some(1e-15)).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

fn crystal_si(tilt_deg: f64) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(&lat, [1, 0, 0], [0, 1, 0], tilt_deg.to_radians(), 0.0, 0.0).unwrap();
    CrystalTarget::new(lat, o).with_thermal(Thermal::new(300.0, THETA_D_SI))
}

fn bca<'a>(
    stack: &'a Stack,
    region: usize,
    ion_z: u8,
    energy_ev: f64,
    tilt_deg: f64,
    count: u64,
    ls: &'a LindhardScharff,
) -> Bca<'a> {
    let b = Beam {
        ion: Ion::new(ion_z).unwrap(),
        energy_ev,
        polar_rad: tilt_deg.to_radians(),
        azimuth_rad: 0.0,
        count,
    };
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = SEED;
    Bca::new(b, stack, cfg, ls, table())
        .unwrap()
        .with_crystal(crystal_si(tilt_deg), &[region])
        .unwrap()
}

/// What happened to the primary and its recoils, in event order, with depths
/// measured from `origin` (the crystal surface).
///
/// `log` entries are `(history, kind, position - origin, energy, layer)`;
/// kind 0 = primary stop, 1 = recoil, 2 = primary electronic loss.
#[derive(Default)]
struct Events {
    origin: f64,
    log: Vec<(u64, u8, [f64; 3], f64, usize)>,
    history: u64,
}

impl BcaTally for Events {
    fn begin_history(&mut self, index: u64) {
        self.history = index;
    }
    fn electronic(&mut self, p: &Particle, _from: [f64; 3], _c: ElectronicChannel, e: f64) {
        if p.is_primary() && e > 1e-3 {
            let mut at = p.pos;
            at[0] -= self.origin;
            self.log.push((self.history, 2, at, p.energy_ev, p.layer));
        }
    }
    fn recoil(&mut self, r: &Particle) {
        let mut at = r.pos;
        at[0] -= self.origin;
        self.log.push((self.history, 1, at, r.energy_ev, r.layer));
    }
    fn stopped(&mut self, p: &Particle) {
        if p.is_primary() {
            let mut at = p.pos;
            at[0] -= self.origin;
            self.log.push((self.history, 0, at, p.energy_ev, p.layer));
        }
    }
    fn merge(&mut self, mut o: Self) {
        self.log.append(&mut o.log);
    }
}

/// An ion that crosses a transparent overlayer arrives at the crystal with
/// the beam's state and then behaves exactly as an ion started at the bare
/// crystal's surface: same recoils at the same sites, same stopping point.
///
/// The overlayer is a whole number of lattice constants thick and the beam is
/// along the surface normal ([100]), so the lattice seen from the interface is
/// the lattice seen from the surface, and the per-history lattice shift (drawn
/// from a stream that does not depend on the transport) is the same too.
#[test]
fn crossing_the_interface_equals_starting_at_it() {
    let a = Lattice::silicon().orthogonal_cell().0[0];
    let ls = LindhardScharff::new();
    let (z, e, n) = (5u8, 3.0e3, 60);
    let bare_stack = Stack::semi_infinite(si());
    let bare = bca(&bare_stack, 0, z, e, 0.0, n, &ls);
    let want = bare.run(Events::default).unwrap();
    for cells in [1.0, 4.0] {
        let t = cells * a;
        let stack = Stack::new(vec![(transparent(), t)], Some(si())).unwrap();
        let b = bca(&stack, 1, z, e, 0.0, n, &ls);
        let got = b
            .run(|| Events {
                origin: t,
                ..Events::default()
            })
            .unwrap();
        let mut compared = 0;
        for h in 0..n {
            let first = |l: &Events| -> Vec<_> {
                l.log
                    .iter()
                    .filter(|ev| ev.0 == h)
                    .take(PREFIX)
                    .copied()
                    .collect()
            };
            let (g, w) = (first(&got), first(&want));
            assert_eq!(g.len(), w.len(), "{cells} cells, history {h}: event count");
            for (g, w) in g.iter().zip(&w) {
                assert_eq!(g.1, w.1, "{cells} cells, history {h}: event kind");
                // Layer: the crystal is layer 1 behind the overlayer, 0 bare.
                assert_eq!(g.4, w.4 + 1, "{cells} cells, history {h}: layer");
                for k in 0..3 {
                    assert!(
                        (g.2[k] - w.2[k]).abs() < 1e-6 * NM,
                        "{cells} cells, history {h}, axis {k}: {} vs {}",
                        g.2[k],
                        w.2[k]
                    );
                }
                assert!(
                    (g.3 - w.3).abs() <= 1e-9 * w.3.max(1.0),
                    "{cells} cells, history {h}: energy {} vs {}",
                    g.3,
                    w.3
                );
                compared += 1;
            }
        }
        assert!(
            compared >= n as usize * PREFIX / 2,
            "{cells} cells: only {compared} events"
        );
    }
}

/// Bookkeeping by layer in a crystal-in-stack run: recoils, lattice deposits
/// and escaped atoms carry the layer of the place they belong to, for the
/// amorphous oxide and for the crystal substrate alike.
#[derive(Default)]
struct Layers {
    thickness: f64,
    recoil_by_layer: [u64; 2],
    lattice_by_layer: [u64; 2],
    sputtered_by_origin: [u64; 2],
    bad: Vec<String>,
}

impl Layers {
    fn expect(&mut self, what: &str, x: f64, layer: usize) {
        // Exactly on the interface the layer is either side's; elsewhere the
        // depth decides.
        if (x - self.thickness).abs() < 1e-13 {
            return;
        }
        let want = usize::from(x > self.thickness);
        if layer != want {
            self.bad
                .push(format!("{what} at x = {x} has layer {layer}, want {want}"));
        }
    }
}

impl BcaTally for Layers {
    fn recoil(&mut self, r: &Particle) {
        self.expect("recoil", r.pos[0], r.layer);
        if r.origin_layer != r.layer {
            self.bad.push("recoil origin_layer != layer".into());
        }
        self.recoil_by_layer[r.layer.min(1)] += 1;
    }
    fn lattice(&mut self, at: [f64; 3], layer: usize, _kind: LatticeDeposit, _e: f64) {
        self.expect("lattice deposit", at[0], layer);
        self.lattice_by_layer[layer.min(1)] += 1;
    }
    fn electronic(&mut self, _p: &Particle, _f: [f64; 3], _c: ElectronicChannel, _e: f64) {}
    fn escaped(&mut self, p: &Particle, face: Face) {
        if !p.is_primary() && matches!(face, Face::Front) {
            self.sputtered_by_origin[p.origin_layer.min(1)] += 1;
        }
    }
    fn end_history(&mut self, _i: u64, _b: &EnergyBudget) {}
    fn merge(&mut self, o: Self) {
        for k in 0..2 {
            self.recoil_by_layer[k] += o.recoil_by_layer[k];
            self.lattice_by_layer[k] += o.lattice_by_layer[k];
            self.sputtered_by_origin[k] += o.sputtered_by_origin[k];
        }
        self.bad.extend(o.bad);
    }
}

#[test]
fn recoils_and_deposits_carry_their_layer_in_a_crystal_stack() {
    let ls = LindhardScharff::new();
    let t = 3.0 * NM;
    let stack = Stack::new(vec![(sio2(), t)], Some(si())).unwrap();
    let b = bca(&stack, 1, 33, 3.0e4, 45.0, 200, &ls);
    let tally = b
        .run(|| Layers {
            thickness: t,
            ..Layers::default()
        })
        .unwrap();
    assert!(
        tally.bad.is_empty(),
        "{:?}",
        &tally.bad[..tally.bad.len().min(5)]
    );
    // Both the amorphous oxide and the crystal produced recoils and deposits.
    assert!(tally.recoil_by_layer[0] > 0, "{:?}", tally.recoil_by_layer);
    assert!(tally.recoil_by_layer[1] > 0, "{:?}", tally.recoil_by_layer);
    assert!(tally.lattice_by_layer[0] > 0 && tally.lattice_by_layer[1] > 0);
    // Atoms leaving the front face started in the oxide or came back from the
    // crystal; the origin is a real layer either way.
    println!("sputtered by origin layer: {:?}", tally.sputtered_by_origin);
    assert!(tally.sputtered_by_origin[0] > 0);
}
