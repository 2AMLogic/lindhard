//! Physical constants, SI units.
//!
//! Source for every value: CODATA 2022 recommended values, P. J. Mohr,
//! D. B. Newell, B. N. Taylor, E. Tiesinga, "CODATA recommended values of the
//! fundamental physical constants: 2022", Rev. Mod. Phys. 97, 025002 (2025),
//! doi:10.1103/RevModPhys.97.025002. Values marked "exact" are fixed by the
//! 2019 SI redefinition; the rest carry their standard uncertainty in the doc
//! comment.
//!
//! Verified digit by digit on 2026-10-06 (#15) against NIST's listing of the
//! 2022 CODATA adjustment, <https://physics.nist.gov/cuu/Constants/Table/allascii.txt>
//! (rows "elementary charge", "Planck constant", "speed of light in vacuum",
//! "Boltzmann constant", "Avogadro constant", "vacuum electric permittivity",
//! "fine-structure constant", "electron mass", "proton mass", "atomic mass
//! constant", "Bohr radius", "Hartree energy", "electron mass energy
//! equivalent"). Every stored value equals the listed one. [`HBAR`] and
//! [`COULOMB_E2`] are not stored values but expressions of the stored ones;
//! the tests check them, and four stored values, against other rows of the
//! same listing.

/// Elementary charge `e`, C. CODATA 2022: 1.602 176 634e-19 (exact).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Planck constant `h`, J s. CODATA 2022: 6.626 070 15e-34 (exact).
pub const PLANCK: f64 = 6.626_070_15e-34;

/// Reduced Planck constant `h / 2π`, J s. Derived from the exact [`PLANCK`];
/// CODATA 2022 lists 1.054 571 817...e-34 (exact).
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

/// Coulomb constant `e² / (4π ε₀)`, J m. Derived from [`ELEMENTARY_CHARGE`]
/// (C) and [`VACUUM_PERMITTIVITY`] (F/m = C²/(J m)), so the unit is J m
/// (about 2.307e-28 J m, i.e. 14.3996 eV Å). CODATA lists no row for this
/// combination; its relative uncertainty is that of ε₀, 1.6e-10. This is the
/// combination the screened-Coulomb interatomic potentials in the ion engine
/// are written with.
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

    /// Rows of the NIST listing of the 2022 CODATA adjustment
    /// (<https://physics.nist.gov/cuu/Constants/Table/allascii.txt>, read
    /// 2026-10-06) that are not stored here, recomputed from the stored values.
    /// Where both sides are printed to 11 or more significant digits the
    /// tolerance only covers the rounding of those digits (a few 1e-12), so a
    /// wrong digit in a stored input fails; the identities (`a0`, `alpha`,
    /// `m_e c^2`) hold to the CODATA rounding of their inputs, 1e-9.
    #[test]
    fn derived_values_match_codata_2022_listing() {
        let rel = |a: f64, b: f64| (a / b - 1.0).abs();
        // "reduced Planck constant": 1.054 571 817... e-34 J s (exact, truncated).
        assert!(rel(HBAR, 1.054_571_817e-34) < 1e-9);
        // "reduced Planck constant in eV s": 6.582 119 569... e-16 (exact, truncated).
        assert!(rel(HBAR / ELEMENTARY_CHARGE, 6.582_119_569e-16) < 1e-9);
        // "electron mass energy equivalent in MeV": 0.510 998 950 69(16).
        let mec2_mev = ELECTRON_REST_ENERGY / ELEMENTARY_CHARGE / 1.0e6;
        assert!(rel(mec2_mev, 0.510_998_950_69) < 3e-11, "{mec2_mev}");
        // "electron mass energy equivalent" = m_e c^2: 8.187 105 7880(26) e-14 J.
        let mec2 = ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
        assert!(rel(mec2, ELECTRON_REST_ENERGY) < 1e-9, "{mec2}");
        // "Hartree energy in eV": 27.211 386 245 981(30).
        let eh_ev = HARTREE_ENERGY / ELEMENTARY_CHARGE;
        assert!(rel(eh_ev, 27.211_386_245_981) < 1e-12, "{eh_ev}");
        // "molar mass constant" M_u = N_A m_u: 1.000 000 001 05(31) e-3 kg/mol.
        let mu = AVOGADRO * ATOMIC_MASS_UNIT;
        assert!(rel(mu, 1.000_000_001_05e-3) < 2e-11, "{mu}");
        // Bohr radius a0 = hbar / (alpha m_e c), to CODATA rounding.
        let a0 = HBAR / (FINE_STRUCTURE * ELECTRON_MASS * SPEED_OF_LIGHT);
        assert!(rel(a0, BOHR_RADIUS) < 1e-9, "{a0}");
        // alpha = e^2 / (4 pi eps0 hbar c), i.e. COULOMB_E2 = alpha hbar c.
        assert!(rel(COULOMB_E2, FINE_STRUCTURE * HBAR * SPEED_OF_LIGHT) < 1e-9);
    }
}
