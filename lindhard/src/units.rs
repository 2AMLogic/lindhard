//! Unit conventions.
//!
//! All physical quantities are SI internally (metres, kilograms, seconds,
//! joules, kelvin). Two places deliberately use customary atomic-physics units
//! at the API boundary, and say so in their names: energies carry an `_ev`
//! suffix (electronvolts) and table densities are quoted in g/cm³ in the
//! source tables. The factors below convert between them; they are exact
//! because the SI 2019 redefinition fixes the elementary charge and the gram is
//! a defined fraction of the kilogram.

use crate::constants::ELEMENTARY_CHARGE;

/// Joules per electronvolt. One eV is the elementary charge times one volt
/// (CODATA 2022 value of `e`, exact; see [`ELEMENTARY_CHARGE`]).
pub const J_PER_EV: f64 = ELEMENTARY_CHARGE;

/// Kilograms per cubic metre for a density of 1 g/cm³ (exact: 1e-3 kg / 1e-6 m³).
pub const KG_M3_PER_G_CM3: f64 = 1.0e3;

/// Metres per ångström (exact by definition of the prefix).
pub const M_PER_ANGSTROM: f64 = 1.0e-10;

/// Convert electronvolts to joules.
#[inline]
pub fn ev_to_j(e_ev: f64) -> f64 {
    e_ev * J_PER_EV
}

/// Convert joules to electronvolts.
#[inline]
pub fn j_to_ev(e_j: f64) -> f64 {
    e_j / J_PER_EV
}

/// Convert g/cm³ to kg/m³.
#[inline]
pub fn g_cm3_to_kg_m3(rho: f64) -> f64 {
    rho * KG_M3_PER_G_CM3
}

/// Convert kg/m³ to g/cm³.
#[inline]
pub fn kg_m3_to_g_cm3(rho: f64) -> f64 {
    rho / KG_M3_PER_G_CM3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips() {
        assert!((j_to_ev(ev_to_j(25.0)) - 25.0).abs() < 1e-12);
        assert!((kg_m3_to_g_cm3(g_cm3_to_kg_m3(2.33)) - 2.33).abs() < 1e-12);
    }
}
