//! Element table, Z = 1..=92.
//!
//! The table is `const` data with no I/O, so table builders and a WASM front
//! end can use it directly. Provenance is also recorded in
//! `docs/data-provenance.md`.
//!
//! # Sources
//!
//! * **Atomic weight** (`atomic_weight`, g/mol, i.e. relative atomic mass):
//!   the CIAAW / IUPAC standard atomic weights (2021 table; J. Meija et al.,
//!   "Atomic weights of the elements 2013", Pure Appl. Chem. 88, 265 (2016), and
//!   the CIAAW 2021 revisions at <https://www.ciaaw.org>). For elements whose
//!   standard weight is an interval (H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl),
//!   the **conventional** single value published by CIAAW is used. Elements
//!   with no standard atomic weight (Tc, Pm, Po, At, Rn, Fr, Ra, Ac) carry the
//!   mass number of their longest-lived isotope instead and have
//!   `weight_is_mass_number == true`; that is a convention, not a measurement
//!   of a natural mixture.
//! * **Density** (`density_g_cm3`): density of the solid at about 20 °C
//!   (25 °C or the stated allotrope where the handbook says so) as listed in
//!   W. M. Haynes (ed.), *CRC Handbook of Chemistry and Physics*, "Physical
//!   Constants of the Elements" (the common room-temperature allotrope:
//!   graphite for C, white P, alpha-S, white Sn). It is `None` where no solid
//!   density at ambient conditions exists (gases: H, He, N, O, F, Ne, Cl, Ar,
//!   Kr, Xe, Rn; liquids: Br, Hg; no data: At, Fr). Densities are a default for a pure
//!   element only; compounds need an explicit density override.
//! * **Sublimation (cohesive) enthalpy** (`default_surface_binding_ev`): where
//!   a measured value is tabulated, the cohesive energy per atom (eV) of
//!   C. Kittel, *Introduction to Solid State Physics*, 8th ed. (Wiley, 2005),
//!   chapter 3 (cohesive energies of the elements), which is the enthalpy of
//!   atomization at 0 K. Using it as the surface binding energy `E_s` is the
//!   usual binary-collision convention (Sigmund, Phys. Rev. 184, 383 (1969)),
//!   **not** a measured surface barrier. Elements without an entry are `None`:
//!   the user must set `E_s` explicitly rather than have one invented.
//! * **Displacement energy** (`default_displacement_ev`): threshold
//!   displacement energies as recommended in ASTM E521 ("Standard Practice for
//!   Investigating the Effects of Neutron Radiation Damage Using Charged-Particle
//!   Irradiation"), the convention underlying the NRT damage model
//!   (Norgett, Robinson, Torrens, Nucl. Eng. Des. 33, 50 (1975)). These are
//!   conventions for damage accounting, not measured properties of a specific
//!   crystal direction. `None` where no convention value is recorded.

/// One element of the periodic table with its default data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Element {
    /// Atomic number.
    pub z: u8,
    /// Element symbol, e.g. `"Si"`.
    pub symbol: &'static str,
    /// Standard (conventional) atomic weight in g/mol; see the module docs.
    pub atomic_weight: f64,
    /// True if `atomic_weight` is the mass number of the longest-lived isotope
    /// because IUPAC publishes no standard atomic weight for the element.
    pub weight_is_mass_number: bool,
    /// Default solid density in g/cm³, if one is defined; see the module docs.
    pub density_g_cm3: Option<f64>,
    /// Default surface binding energy `E_s` in eV (cohesive-energy convention).
    pub default_surface_binding_ev: Option<f64>,
    /// Default displacement energy `E_d` in eV (ASTM E521 convention).
    pub default_displacement_ev: Option<f64>,
}

impl Element {
    /// Atomic mass in kg: `atomic_weight` times the atomic mass constant.
    ///
    /// (The molar-mass constant differs from 1 g/mol only at the 1e-10 level;
    /// the standard-atomic-weight uncertainty is far larger.)
    pub fn mass_kg(&self) -> f64 {
        self.atomic_weight * crate::constants::ATOMIC_MASS_UNIT
    }

    /// Default density in kg/m³, if defined.
    pub fn density_kg_m3(&self) -> Option<f64> {
        self.density_g_cm3.map(crate::units::g_cm3_to_kg_m3)
    }
}

/// Number of elements in the table.
pub const NUM_ELEMENTS: usize = 92;

/// Look an element up by atomic number (1..=92).
pub fn element(z: u8) -> Option<&'static Element> {
    if z == 0 {
        return None;
    }
    ELEMENTS.get(usize::from(z) - 1)
}

/// Look an element up by case-sensitive symbol (`"Si"`, not `"si"`).
pub fn element_by_symbol(symbol: &str) -> Option<&'static Element> {
    ELEMENTS.iter().find(|e| e.symbol == symbol)
}

const fn el(
    z: u8,
    symbol: &'static str,
    atomic_weight: f64,
    weight_is_mass_number: bool,
    density_g_cm3: Option<f64>,
    default_surface_binding_ev: Option<f64>,
    default_displacement_ev: Option<f64>,
) -> Element {
    Element {
        z,
        symbol,
        atomic_weight,
        weight_is_mass_number,
        density_g_cm3,
        default_surface_binding_ev,
        default_displacement_ev,
    }
}

/// The table, indexed by `Z - 1`.
// Sn's cohesive energy (3.14 eV) is data, not pi.
#[allow(clippy::approx_constant)]
pub static ELEMENTS: [Element; NUM_ELEMENTS] = [
    el(1, "H", 1.008, false, None, None, None),
    el(2, "He", 4.002602, false, None, None, None),
    el(3, "Li", 6.94, false, Some(0.534), Some(1.63), None),
    el(4, "Be", 9.0121831, false, Some(1.85), Some(3.32), None),
    el(5, "B", 10.81, false, Some(2.34), Some(5.81), None),
    el(6, "C", 12.011, false, Some(2.267), Some(7.37), None),
    el(7, "N", 14.007, false, None, None, None),
    el(8, "O", 15.999, false, None, None, None),
    el(9, "F", 18.998403162, false, None, None, None),
    el(10, "Ne", 20.1797, false, None, None, None),
    el(11, "Na", 22.98976928, false, Some(0.97), Some(1.11), None),
    el(12, "Mg", 24.305, false, Some(1.74), Some(1.51), None),
    el(
        13,
        "Al",
        26.9815384,
        false,
        Some(2.70),
        Some(3.39),
        Some(25.0),
    ),
    el(14, "Si", 28.085, false, Some(2.329), Some(4.63), None),
    el(15, "P", 30.973761998, false, Some(1.82), None, None),
    el(16, "S", 32.06, false, Some(2.07), None, None),
    el(17, "Cl", 35.45, false, None, None, None),
    el(18, "Ar", 39.95, false, None, None, None),
    el(19, "K", 39.0983, false, Some(0.862), Some(0.93), None),
    el(20, "Ca", 40.078, false, Some(1.55), Some(1.84), None),
    el(21, "Sc", 44.955908, false, Some(2.989), None, None),
    el(22, "Ti", 47.867, false, Some(4.54), Some(4.85), Some(30.0)),
    el(23, "V", 50.9415, false, Some(6.11), Some(5.31), Some(40.0)),
    el(24, "Cr", 51.9961, false, Some(7.15), Some(4.1), Some(40.0)),
    el(25, "Mn", 54.938043, false, Some(7.44), None, None),
    el(26, "Fe", 55.845, false, Some(7.874), Some(4.28), Some(40.0)),
    el(
        27,
        "Co",
        58.933194,
        false,
        Some(8.90),
        Some(4.39),
        Some(40.0),
    ),
    el(
        28,
        "Ni",
        58.6934,
        false,
        Some(8.908),
        Some(4.44),
        Some(40.0),
    ),
    el(29, "Cu", 63.546, false, Some(8.96), Some(3.49), Some(30.0)),
    el(30, "Zn", 65.38, false, Some(7.134), Some(1.35), None),
    el(31, "Ga", 69.723, false, Some(5.904), Some(2.81), None),
    el(32, "Ge", 72.630, false, Some(5.323), Some(3.85), None),
    el(33, "As", 74.921595, false, Some(5.727), Some(2.96), None),
    el(34, "Se", 78.971, false, Some(4.81), None, None),
    el(35, "Br", 79.904, false, None, None, None),
    el(36, "Kr", 83.798, false, None, None, None),
    el(37, "Rb", 85.4678, false, Some(1.532), None, None),
    el(38, "Sr", 87.62, false, Some(2.64), None, None),
    el(39, "Y", 88.90584, false, Some(4.469), None, None),
    el(40, "Zr", 91.224, false, Some(6.506), Some(6.25), Some(40.0)),
    el(
        41,
        "Nb",
        92.90637,
        false,
        Some(8.57),
        Some(7.57),
        Some(60.0),
    ),
    el(42, "Mo", 95.95, false, Some(10.28), Some(6.82), Some(60.0)),
    el(43, "Tc", 98.0, false, Some(11.0), None, None),
    el(44, "Ru", 101.07, false, Some(12.45), None, None),
    el(45, "Rh", 102.90549, false, Some(12.41), None, None),
    el(46, "Pd", 106.42, false, Some(12.023), Some(3.89), None),
    el(
        47,
        "Ag",
        107.8682,
        false,
        Some(10.49),
        Some(2.95),
        Some(30.0),
    ),
    el(48, "Cd", 112.414, false, Some(8.65), Some(1.16), None),
    el(49, "In", 114.818, false, Some(7.31), Some(2.52), None),
    el(50, "Sn", 118.710, false, Some(7.287), Some(3.14), None),
    el(51, "Sb", 121.760, false, Some(6.685), Some(2.75), None),
    el(52, "Te", 127.60, false, Some(6.232), None, None),
    el(53, "I", 126.90447, false, Some(4.93), None, None),
    el(54, "Xe", 131.293, false, None, None, None),
    el(55, "Cs", 132.90545196, false, Some(1.93), None, None),
    el(56, "Ba", 137.327, false, Some(3.51), None, None),
    el(57, "La", 138.90547, false, Some(6.162), None, None),
    el(58, "Ce", 140.116, false, Some(6.770), None, None),
    el(59, "Pr", 140.90766, false, Some(6.77), None, None),
    el(60, "Nd", 144.242, false, Some(7.01), None, None),
    el(61, "Pm", 145.0, false, Some(7.26), None, None),
    el(62, "Sm", 150.36, false, Some(7.52), None, None),
    el(63, "Eu", 151.964, false, Some(5.264), None, None),
    el(64, "Gd", 157.25, false, Some(7.90), None, None),
    el(65, "Tb", 158.925354, false, Some(8.23), None, None),
    el(66, "Dy", 162.500, false, Some(8.54), None, None),
    el(67, "Ho", 164.930329, false, Some(8.79), None, None),
    el(68, "Er", 167.259, false, Some(9.066), None, None),
    el(69, "Tm", 168.934219, false, Some(9.32), None, None),
    el(70, "Yb", 173.045, false, Some(6.90), None, None),
    el(71, "Lu", 174.9668, false, Some(9.841), None, None),
    el(72, "Hf", 178.486, false, Some(13.31), None, None),
    el(
        73,
        "Ta",
        180.94788,
        false,
        Some(16.65),
        Some(8.1),
        Some(90.0),
    ),
    el(74, "W", 183.84, false, Some(19.25), Some(8.9), Some(90.0)),
    el(75, "Re", 186.207, false, Some(21.02), None, None),
    el(76, "Os", 190.23, false, Some(22.59), None, None),
    el(77, "Ir", 192.217, false, Some(22.56), None, None),
    el(78, "Pt", 195.084, false, Some(21.45), Some(5.84), None),
    el(
        79,
        "Au",
        196.966570,
        false,
        Some(19.30),
        Some(3.81),
        Some(30.0),
    ),
    el(80, "Hg", 200.592, false, None, None, None),
    el(81, "Tl", 204.38, false, Some(11.85), None, None),
    el(82, "Pb", 207.2, false, Some(11.34), Some(2.03), Some(25.0)),
    el(83, "Bi", 208.98040, false, Some(9.78), None, None),
    el(84, "Po", 209.0, false, Some(9.196), None, None),
    el(85, "At", 210.0, false, None, None, None),
    el(86, "Rn", 222.0, false, None, None, None),
    el(87, "Fr", 223.0, false, None, None, None),
    el(88, "Ra", 226.0, false, Some(5.0), None, None),
    el(89, "Ac", 227.0, false, Some(10.07), None, None),
    el(90, "Th", 232.0377, false, Some(11.72), None, None),
    el(91, "Pa", 231.03588, false, Some(15.37), None, None),
    el(92, "U", 238.02891, false, Some(19.1), None, None),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_ordered_and_complete() {
        for (i, e) in ELEMENTS.iter().enumerate() {
            assert_eq!(usize::from(e.z), i + 1);
            assert!(e.atomic_weight > 0.9);
            if let Some(d) = e.density_g_cm3 {
                assert!(d > 0.0 && d < 30.0);
            }
        }
        assert_eq!(element(14).unwrap().symbol, "Si");
        assert_eq!(element(92).unwrap().symbol, "U");
        assert!(element(0).is_none() && element(93).is_none());
        assert_eq!(element_by_symbol("Ga").unwrap().z, 31);
        assert!(element_by_symbol("Xx").is_none());
    }

    #[test]
    fn symbols_are_unique() {
        for (i, a) in ELEMENTS.iter().enumerate() {
            for b in &ELEMENTS[i + 1..] {
                assert_ne!(a.symbol, b.symbol);
            }
        }
    }
}
