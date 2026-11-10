//! Helpers shared by the integration tests. Cargo does not build this
//! directory as its own test target; each test binary pulls it in with
//! `mod common;`.

#![allow(dead_code)]

use std::sync::OnceLock;

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
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

/// Bare Si (Z = 14) with no displacement or surface binding energy set, as
/// the electron tests use it: they never read either.
pub(crate) fn bare_si() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

/// Energy grid (eV) of the synthetic electron tables.
pub(crate) const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];

/// Cumulative probabilities of the synthetic table rows.
pub(crate) const PROB: [f64; 3] = [0.0, 0.5, 1.0];

/// A synthetic cross-section table on the energy grid `grid`: constant
/// inverse mean free path `rate` (1/m) and `row(E, p)` as the `p`-quantile at
/// energy `E` over [`PROB`]. A zero `rate` gives empty rows (no interaction).
/// The numbers are made up to exercise the code, not material data.
pub(crate) fn synthetic_table_on(
    grid: &[f64],
    axis: SamplingAxis,
    rate: f64,
    row: impl Fn(f64, f64) -> f64,
) -> CrossSectionTable {
    let quantiles = grid
        .iter()
        .map(|&e| {
            if rate == 0.0 {
                vec![]
            } else {
                PROB.iter().map(|&p| row(e, p)).collect()
            }
        })
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "synthetic test table".into(),
        axis,
        energy_ev: grid.to_vec(),
        inverse_mfp_per_m: vec![rate; grid.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// [`synthetic_table_on`] on the standard [`GRID`].
pub(crate) fn synthetic_table(
    axis: SamplingAxis,
    rate: f64,
    row: impl Fn(f64, f64) -> f64,
) -> CrossSectionTable {
    synthetic_table_on(&GRID, axis, rate, row)
}
