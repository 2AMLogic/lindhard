//! Amorphous runs must stay bit-identical to the engine before crystal
//! regions existed (issue #180). The digests below were produced by the
//! engine as of `4538a48` (before the crystal flight model) and are a SHA-256
//! of the complete serialised `IonReport` (budget, range histograms and
//! moments, damage, escapes), so any change to any tallied bit fails here.
//!
//! The two cases mirror `examples/b_5keV_si.toml` (B 5 keV into Si at 7
//! degrees, Lindhard-Scharff, constant free path) and a variant that turns on
//! every optional branch of the flight loop (local electronic loss, weak
//! collisions, a layered target with an interface).

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam, ElectronicLoss};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonTally, IonTallyConfig};
use sha2::{Digest, Sha256};

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

fn elemental(z: u8) -> Material {
    let mut m = Material::from_atom_fractions(&[(z, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(z, 15.0).unwrap();
    if m.surface_binding_energy_ev(z).is_err() {
        m.set_surface_binding_energy_ev(z, 2.0).unwrap();
    }
    m
}

fn digest(stack: &Stack, z_ion: u8, energy_ev: f64, cfg: BcaConfig, n: u64) -> String {
    let ls = LindhardScharff::new();
    let beam = Beam {
        ion: Ion::new(z_ion).unwrap(),
        energy_ev,
        polar_rad: 7f64.to_radians(),
        azimuth_rad: 0.0,
        count: n,
    };
    let bca = Bca::new(beam, stack, cfg, &ls, table()).unwrap();
    let species = bca.species_z();
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 200.0 * NM, 200).unwrap(),
        lateral: Binning::new(-100.0 * NM, 100.0 * NM, 100).unwrap(),
        radial: Binning::new(0.0, 100.0 * NM, 50).unwrap(),
        escape_energy: Binning::new(0.0, 5.0e4, 50).unwrap(),
        escape_polar: Binning::new(0.0, std::f64::consts::FRAC_PI_2, 18).unwrap(),
    };
    let tally = bca
        .run(|| IonTally::new(stack, &species, tc).unwrap())
        .unwrap();
    let json = serde_json::to_string(&tally.report(false)).unwrap();
    Sha256::digest(json.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

#[test]
fn b_5kev_si_example_is_unchanged() {
    let stack = Stack::semi_infinite(elemental(14));
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 1;
    let d = digest(&stack, 5, 5.0e3, cfg, 3000);
    assert_eq!(
        d, "358ec00db8250c22f99a4a0c200442dafed8457ffa47d8eb5738f983f6964db0",
        "amorphous B 5 keV Si report changed"
    );
}

#[test]
fn all_optional_branches_are_unchanged() {
    let stack = Stack::new(vec![(elemental(14), 30.0 * NM)], Some(elemental(29))).unwrap();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 7;
    cfg.electronic = ElectronicLoss::EquipartitionLsOr;
    cfg.weak_collisions = 2;
    let d = digest(&stack, 18, 2.0e3, cfg, 1500);
    assert_eq!(
        d, "ad542f70448a84fb590cab530d66202274ff33f216f4c3423a7b9506683eec23",
        "amorphous multi-branch report changed"
    );
}
