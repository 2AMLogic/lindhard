//! Helpers shared by the integration tests. Cargo does not build this
//! directory as its own test target; each test binary pulls it in with
//! `mod common;`.

#![allow(dead_code)]

use std::sync::OnceLock;

use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};

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
