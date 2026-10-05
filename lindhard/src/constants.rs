//! Physical constants, SI units.
//!
//! Source for every value: CODATA 2022 recommended values, E. Tiesinga,
//! P. J. Mohr, D. B. Newell, B. N. Taylor, "CODATA recommended values of the
//! fundamental physical constants: 2022", Rev. Mod. Phys. 97, 025002 (2025)
//! (also NIST SP 961 / the NIST CODATA web listing). Values marked "exact" are
//! fixed by the 2019 SI redefinition; the rest carry their standard
//! uncertainty in the doc comment.

/// Elementary charge `e`, C. CODATA 2022: 1.602 176 634e-19 (exact).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Planck constant `h`, J s. CODATA 2022: 6.626 070 15e-34 (exact).
pub const PLANCK: f64 = 6.626_070_15e-34;

/// Reduced Planck constant `h / 2π`, J s. Derived from the exact [`PLANCK`].
pub const HBAR: f64 = PLANCK / std::f64::consts::TAU;

/// Speed of light in vacuum `c`, m/s. CODATA 2022: 299 792 458 (exact).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Boltzmann constant `k_B`, J/K. CODATA 2022: 1.380 649e-23 (exact).
pub const BOLTZMANN: f64 = 1.380_649e-23;

/// Avogadro constant `N_A`, 1/mol. CODATA 2022: 6.022 140 76e23 (exact).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Vacuum electric permittivity `ε₀`, F/m. CODATA 2022: 8.854 187 8188(14)e-12.
pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_818_8e-12;

/// Fine-structure constant `α`, dimensionless. CODATA 2022: 7.297 352 5643(11)e-3.
pub const FINE_STRUCTURE: f64 = 7.297_352_564_3e-3;

/// Electron mass `m_e`, kg. CODATA 2022: 9.109 383 7139(28)e-31.
pub const ELECTRON_MASS: f64 = 9.109_383_713_9e-31;

/// Proton mass `m_p`, kg. CODATA 2022: 1.672 621 925 95(52)e-27.
pub const PROTON_MASS: f64 = 1.672_621_925_95e-27;

/// Atomic mass constant `m_u` (unified atomic mass unit), kg.
/// CODATA 2022: 1.660 539 068 92(52)e-27.
pub const ATOMIC_MASS_UNIT: f64 = 1.660_539_068_92e-27;

/// Bohr radius `a₀`, m. CODATA 2022: 5.291 772 105 44(82)e-11.
pub const BOHR_RADIUS: f64 = 5.291_772_105_44e-11;

/// Hartree energy `E_h`, J. CODATA 2022: 4.359 744 722 2060(48)e-18.
pub const HARTREE_ENERGY: f64 = 4.359_744_722_206e-18;

/// Electron rest energy `m_e c²`, J. CODATA 2022: 8.187 105 7880(26)e-14
/// (equivalently 0.510 998 950 69(16) MeV).
pub const ELECTRON_REST_ENERGY: f64 = 8.187_105_788_0e-14;

/// Coulomb constant `e² / (4π ε₀)`, J m. Derived from the constants above
/// (about 2.307e-28 J m, i.e. 14.3996 eV Å); this is the combination the
/// screened-Coulomb interatomic potentials in the ion engine are written with.
pub const COULOMB_E2: f64 =
    ELEMENTARY_CHARGE * ELEMENTARY_CHARGE / (4.0 * std::f64::consts::PI * VACUUM_PERMITTIVITY);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coulomb_constant_in_ev_angstrom() {
        // e^2/(4 pi eps0) = 14.3996 eV Angstrom (derived value, 4-digit check).
        let v = COULOMB_E2 / ELEMENTARY_CHARGE / 1.0e-10;
        assert!((v - 14.3996).abs() < 1e-3, "{v}");
    }

    #[test]
    fn hartree_is_alpha_squared_electron_rest_energy() {
        // E_h = alpha^2 m_e c^2 (identity, holds to CODATA rounding).
        let rel =
            (HARTREE_ENERGY / (FINE_STRUCTURE * FINE_STRUCTURE * ELECTRON_REST_ENERGY) - 1.0).abs();
        assert!(rel < 1e-9, "{rel}");
    }
}
