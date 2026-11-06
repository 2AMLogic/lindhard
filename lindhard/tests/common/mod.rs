//! Helpers shared by the integration tests. Cargo does not build this
//! directory as its own test target; each test binary pulls it in with
//! `mod common;`.

#![allow(dead_code)]

use std::sync::OnceLock;

use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::material::Material;

/// One Si-Si ZBL table shared by every test in a binary (the angle depends
/// on the screening function only, so one table serves every pair). The grid
/// is deliberately coarser than `TableSpec::default()` for test speed.
pub(crate) fn zbl_si_table() -> &'static ScatteringTable {
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

// The Si and SiO2 targets below differ only in how the surface binding energy
// is set, which moves sputtering and any binding-dependent tally. Each test
// keeps the variant it always used so no recorded value changes. The 15 eV
// displacement energy and the 2 eV surface binding are test parameters, not
// data.

/// Si with `E_d` = 15 eV and no surface binding energy set.
pub(crate) fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

/// Si with `E_d` = 15 eV and `E_s` = 2 eV set unconditionally.
pub(crate) fn si_with_surface_binding() -> Material {
    let mut m = si();
    m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    m
}

/// Si with `E_d` = 15 eV and `E_s` = 2 eV only where it is not already set.
pub(crate) fn si_surface_binding_if_unset() -> Material {
    let mut m = si();
    if m.surface_binding_energy_ev(14).is_err() {
        m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    }
    m
}

/// Amorphous SiO2 at 2200 kg/m3 with `E_d` = 15 eV for both elements and
/// `E_s` = 2 eV where it is not already set.
pub(crate) fn sio2() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
    for z in [14, 8] {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}
