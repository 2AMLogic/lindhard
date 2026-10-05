//! Shared fixtures for the benchmarks: the three representative problems and
//! a scattering table. Everything is deterministic (fixed seeds, fixed ion
//! counts) so only the timing varies between runs.

#![allow(dead_code)]

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{BcaConfig, Beam};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonTallyConfig};

const NM: f64 = 1e-9;

/// One ZBL table for every benchmark (the angle depends on the screening
/// function only). Same grid as `tests/bca.rs`.
pub fn table() -> &'static ScatteringTable {
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

/// Elemental target with illustrative model parameters: `E_d` = 15 eV, and
/// `E_s` = 2 eV where the element table has no default. Benchmark inputs,
/// not recommended values.
pub fn elemental(z: u8) -> Material {
    let mut m = Material::from_atom_fractions(&[(z, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(z, 15.0).unwrap();
    if m.surface_binding_energy_ev(z).is_err() {
        m.set_surface_binding_energy_ev(z, 2.0).unwrap();
    }
    m
}

/// A representative problem: ion, energy (eV), substrate Z and a label.
pub struct Problem {
    pub label: &'static str,
    pub z_ion: u8,
    pub energy_ev: f64,
    pub z_target: u8,
}

pub const PROBLEMS: [Problem; 3] = [
    Problem {
        label: "B_5keV_Si",
        z_ion: 5,
        energy_ev: 5.0e3,
        z_target: 14,
    },
    Problem {
        label: "As_50keV_Si",
        z_ion: 33,
        energy_ev: 5.0e4,
        z_target: 14,
    },
    Problem {
        label: "Ar_1keV_Cu",
        z_ion: 18,
        energy_ev: 1.0e3,
        z_target: 29,
    },
];

pub fn stack(p: &Problem) -> Stack {
    Stack::semi_infinite(elemental(p.z_target))
}

pub fn beam(p: &Problem, count: u64) -> Beam {
    Beam {
        ion: Ion::new(p.z_ion).unwrap(),
        energy_ev: p.energy_ev,
        polar_rad: 7f64.to_radians(),
        azimuth_rad: 0.0,
        count,
    }
}

pub fn config(seed: u64) -> BcaConfig {
    let mut c = BcaConfig::new(5.0, 2.0);
    c.seed = seed;
    c
}

pub fn tally_config() -> IonTallyConfig {
    IonTallyConfig {
        depth: Binning::new(0.0, 200.0 * NM, 200).unwrap(),
        lateral: Binning::new(-100.0 * NM, 100.0 * NM, 100).unwrap(),
        radial: Binning::new(0.0, 100.0 * NM, 50).unwrap(),
        escape_energy: Binning::new(0.0, 5.0e4, 50).unwrap(),
        escape_polar: Binning::new(0.0, std::f64::consts::FRAC_PI_2, 18).unwrap(),
    }
}
