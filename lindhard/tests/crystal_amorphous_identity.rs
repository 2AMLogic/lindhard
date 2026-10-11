//! Amorphous transport must be unchanged by the crystal flight model
//! (issue #180). Two kinds of check:
//!
//! * **Portable, in-process A/B** (every platform): the complete serialised
//!   `IonReport` of a plain `Bca` is compared with the same run after
//!   `with_crystal` attaches a crystal to a region the particles never reach.
//!   Both reports come from the same process, so they must be bit-identical
//!   everywhere. This tests the property itself: attaching a crystal does not
//!   touch transport, or the random streams, in amorphous regions.
//! * **Fixed digests** (x86_64 Linux only): SHA-256 of the report produced
//!   by the engine as of `4538a48` (before the crystal flight model), so any
//!   change to any tallied bit fails here. They were recorded on x86_64
//!   Linux; libm `exp`/`ln`/`pow`/trigonometry differ in the last bit on
//!   other targets (macOS, aarch64 Linux), so the serialised floats, and the
//!   digests, are platform dependent. The gate keeps the before/after
//!   evidence for the refactor on the platform where it was measured.
//!
//! The cases mirror `examples/b_5keV_si.toml` (B 5 keV into Si at 7 degrees,
//! Lindhard-Scharff, constant free path) and a variant that turns on the
//! optional branches of the flight loop (local electronic loss, weak
//! collisions in the digest case, a layered target with an interface).

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam, CrystalTarget, ElectronicLoss};
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonTally, IonTallyConfig};

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

/// Serialised `IonReport` of a 7-degree run; with `crystal`, a silicon
/// crystal is attached to the given regions.
fn report_json(
    stack: &Stack,
    z_ion: u8,
    energy_ev: f64,
    cfg: BcaConfig,
    n: u64,
    crystal: Option<&[usize]>,
) -> String {
    let ls = LindhardScharff::new();
    let beam = Beam {
        ion: Ion::new(z_ion).unwrap(),
        energy_ev,
        polar_rad: 7f64.to_radians(),
        azimuth_rad: 0.0,
        count: n,
    };
    let mut bca = Bca::new(beam, stack, cfg, &ls, table()).unwrap();
    if let Some(regions) = crystal {
        bca = bca.with_crystal(crystal_si(), regions).unwrap();
    }
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
    serde_json::to_string(&tally.report(false)).unwrap()
}

/// Silicon cut on (100), oriented along the beam's 7-degree tilt.
fn crystal_si() -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(&lat, [1, 0, 0], [0, 1, 0], 7f64.to_radians(), 0.0, 0.0).unwrap();
    CrystalTarget::new(lat, o)
}

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
fn digest(stack: &Stack, z_ion: u8, energy_ev: f64, cfg: BcaConfig, n: u64) -> String {
    use sha2::{Digest, Sha256};
    let json = report_json(stack, z_ion, energy_ev, cfg, n, None);
    Sha256::digest(json.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

#[test]
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
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
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
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

/// B 5 keV (range of order 20 nm) into 1 um of amorphous Si on a Si crystal
/// substrate that no particle reaches.
#[test]
fn unreachable_crystal_leaves_b_5kev_si_identical() {
    let stack = Stack::new(vec![(elemental(14), 1000.0 * NM)], Some(elemental(14))).unwrap();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 1;
    let plain = report_json(&stack, 5, 5.0e3, cfg, 3000, None);
    let with = report_json(&stack, 5, 5.0e3, cfg, 3000, Some(&[1]));
    assert!(
        plain == with,
        "attaching an unreachable crystal changed the amorphous report"
    );
}

/// Ar 2 keV through 30 nm Si into 1 um Cu (an interface, local electronic
/// loss) with a Si crystal substrate behind it that no particle reaches.
/// Weak collisions are amorphous-only and cannot be combined with a crystal.
#[test]
fn unreachable_crystal_leaves_layered_run_identical() {
    let stack = Stack::new(
        vec![(elemental(14), 30.0 * NM), (elemental(29), 1000.0 * NM)],
        Some(elemental(14)),
    )
    .unwrap();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 7;
    cfg.electronic = ElectronicLoss::EquipartitionLsOr;
    let plain = report_json(&stack, 18, 2.0e3, cfg, 1500, None);
    let with = report_json(&stack, 18, 2.0e3, cfg, 1500, Some(&[2]));
    assert!(
        plain == with,
        "attaching an unreachable crystal changed the amorphous report"
    );
}
