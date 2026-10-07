//! Element table, Z = 1..=92.
//!
//! The table is `const` data with no I/O, so table builders and a WASM front
//! end can use it directly. Provenance is also recorded in
//! `docs/data-provenance.md`.
//!
//! # Sources
//!
//! * **Atomic weight** (`atomic_weight`, g/mol, i.e. relative atomic mass):
//!   the CIAAW / IUPAC standard atomic weights as tabulated at
//!   <https://www.ciaaw.org/atomic-weights.htm> ("Standard Atomic Weights
//!   2024": the 2021 report, T. Prohaska et al., "Standard atomic weights of
//!   the elements 2021", Pure Appl. Chem. 94, 573 (2022),
//!   doi:10.1515/pac-2019-0603, plus the CIAAW 2024 revisions of Gd, Lu and
//!   Zr). For the 14 elements whose standard weight is an interval (H, Li, B,
//!   C, N, O, Mg, Si, S, Cl, Ar, Br, Tl, Pb) the single value is the
//!   abridged value of <https://www.ciaaw.org/abridged-atomic-weights.htm>,
//!   which the 2021 report (Table 1, column 7) says replaces, and equals, the
//!   former conventional atomic weights; for H that value is printed 1.0080,
//!   the former conventional 1.008. Elements with no standard atomic weight
//!   (Tc, Pm, Po, At, Rn, Fr, Ra, Ac) carry the mass number of their
//!   longest-lived isotope as listed by CIAAW at
//!   <https://www.ciaaw.org/radioactive-elements.htm> (from NUBASE2020,
//!   F. G. Kondev et al., Chin. Phys. C 45, 030001 (2021)) and have
//!   `weight_is_mass_number == true`; that is a convention, not a measurement
//!   of a natural mixture. For Tc that page lists 97 and 98 as equally
//!   long-lived; 98 is used. All 92 entries checked against these pages on
//!   2026-10-06 (#15).
//! * **Density** (`density_g_cm3`): density of the solid, g/cm³, as printed
//!   in A. Thompson et al., *X-Ray Data Booklet*, LBNL/PUB-490 Rev. 3
//!   (Lawrence Berkeley National Laboratory, October 2009), section 5.2,
//!   Table 5-2 "Properties of the elements", pp. 5-5 to 5-8
//!   (<https://xdb.lbl.gov/xdb-new.pdf>, PDF pages 153-156). Section 5.2 says
//!   the data were taken mostly from D. R. Lide (ed.), *CRC Handbook of
//!   Chemistry and Physics*, 80th ed. (CRC Press, 1999), and that densities
//!   are specific gravities at 20 °C unless a superscript gives another
//!   temperature (25 °C for Si, Sc, Ni, Zn, Ge, Y, La, Ce, Nd, Pm, Sm, Eu,
//!   Gd, Dy, Ho, Er, Tm, Lu; V 18.7 °C, Ga 29.6 °C, Ir 17 °C). Tc and Pa are
//!   marked "calculated" there, U "~18.95". The table names no allotrope
//!   except graphite for C. All populated entries checked on 2026-10-07
//!   (#15; the PDF's SHA-256 is recorded in `docs/data-provenance.md`), and
//!   37 values set to the printed ones. Four entries the booklet does not pin
//!   were checked on 2026-10-07 (#15) against other open primary
//!   compilations, all U.S. government reports:
//!   - Cr 7.20 (was 7.15, which lay **outside** the booklet's printed range
//!     7.18-7.20): the upper end of that range, and the density "calculated
//!     from the NBS lattice constant", 7.200 at 25 °C, in H. E. Swanson,
//!     R. K. Fuyat and G. M. Ugrinic, *Standard X-ray Diffraction Powder
//!     Patterns*, NBS Circular 539, Vol. V (1955), p. 20
//!     (<https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular539v5.pdf>).
//!     This is an X-ray (crystal) density.
//!   - Au 19.30: the booklet prints "~19.3"; NBS Circular 539, Vol. I
//!     (Swanson and Tatge, 1953), p. 33, prints 19.302 at 25 °C from its
//!     lattice constant
//!     (<https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular539v1.pdf>).
//!   - Ra 5.0: no density in the booklet; "a specific gravity of
//!     approximately 5.0" in H. W. Kirby and M. L. Salutsky, *The
//!     Radiochemistry of Radium*, NAS-NS 3057 (1964), p. 3, citing the
//!     International Critical Tables (<https://www.osti.gov/servlets/purl/4560824>).
//!   - Ac 10.07: no density in the booklet; the crystal density of fcc Ac
//!     metal (a = 5.311 Å), Table 4, p. 11, of J. D. Farr, A. L. Giorgi,
//!     M. G. Bowman and R. K. Money, "The crystal structure of actinium metal
//!     and actinium hydride", Los Alamos report LA-1545 (1953)
//!     (<https://www.osti.gov/servlets/purl/4397640>), published as
//!     J. Inorg. Nucl. Chem. 18, 42 (1961).
//!
//!   Still **not verified** to their stored digits: C 2.267 (printed as the
//!   range 1.9-2.3) and Mn 7.44 (the upper end of the printed range
//!   7.21-7.44). For Mn, the X-ray density of alpha-Mn, "(calculated) 7.475
//!   g/cm³" at 25 °C (a = 8.9121 Å, Z = 58), is printed in M. C. Morris
//!   et al., *Standard X-ray Diffraction Powder Patterns*, NBS Monograph 25,
//!   Section 17 (1980), p. 50
//!   (<https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nbsmonograph25-17.pdf>,
//!   checked 2026-10-07, #15). It lies above the booklet's range, so unlike
//!   Cr it does not pin a value inside that range; 7.44 is kept, unverified.
//!   Si keeps 2.329: the booklet prints 2.33
//!   at 25 °C, and 2.329 is the crystal density `M(Si) / V_m(Si)` from the
//!   CODATA 2022 "molar volume of silicon" (1.205 883 199e-5 m³/mol,
//!   <https://physics.nist.gov/cuu/Constants/Table/allascii.txt>) for any
//!   weight in the standard interval [28.084, 28.086]. The specific gravities
//!   are used as g/cm³; the booklet does not state the reference water
//!   state. Later CRC editions may differ; none was seen. It is `None` where
//!   no solid density at ambient conditions exists (gases: H, He, N, O, F,
//!   Ne, Cl, Ar, Kr, Xe, Rn; liquids: Br, Hg; no data: At, Fr). Densities are
//!   a default for a pure element only; compounds need an explicit density
//!   override.
//! * **Sublimation (cohesive) enthalpy** (`default_surface_binding_ev`): where
//!   a measured value is tabulated, the cohesive energy per atom (eV) of
//!   C. Kittel, *Introduction to Solid State Physics*, 8th ed. (Wiley, 2005),
//!   chapter 3 (cohesive energies of the elements), which is the enthalpy of
//!   atomization at 0 K. Using it as the surface binding energy `E_s` is the
//!   usual binary-collision convention (Sigmund, Phys. Rev. 184, 383 (1969)),
//!   **not** a measured surface barrier. Elements without an entry are `None`:
//!   the user must set `E_s` explicitly rather than have one invented.
//!   **Not verified against Kittel** (the book could not be reached; the
//!   publisher's free excerpts of the 8th edition, its contents and index,
//!   place the cohesive-energy table on p. 50 of chapter 3 but do not
//!   include that page). Checked
//!   instead against L. Brewer, "The cohesive energies of the elements",
//!   report LBL-3720 (1975), Table I, energy of atomization to the gaseous
//!   ground state at 0 K in kcal/gram-atom (<https://www.osti.gov/servlets/purl/7187973>):
//!   32 of the 35 populated values agree with it to the printed digit; Li
//!   and Mg differ by one in the last digit and B is 5.81 here against
//!   5.77 eV there (see `docs/data-provenance.md`).
//! * **Displacement energy** (`default_displacement_ev`): threshold
//!   displacement energies as recommended in ASTM E521 ("Standard Practice for
//!   Investigating the Effects of Neutron Radiation Damage Using Charged-Particle
//!   Irradiation"), the convention underlying the NRT damage model
//!   (Norgett, Robinson, Torrens, Nucl. Eng. Des. 33, 50 (1975)). These are
//!   conventions for damage accounting, not measured properties of a specific
//!   crystal direction. `None` where no convention value is recorded.
//!   **Not verified**: the standard could not be reached (see
//!   `docs/data-provenance.md`). The values are also not those of the
//!   Greenwood-Smither compilation used by SPECTER and NJOY (they differ for
//!   Al, Ti, Cu, Nb and Ag), which is not E521 either.

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
    /// (The molar-mass constant differs from 1 g/mol by about 1e-9, CODATA
    /// 2022 `M_u` = 1.000 000 001 05e-3 kg/mol; the standard-atomic-weight
    /// uncertainty is far larger.)
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
    el(4, "Be", 9.0121831, false, Some(1.848), Some(3.32), None),
    el(5, "B", 10.81, false, Some(2.34), Some(5.81), None),
    el(6, "C", 12.011, false, Some(2.267), Some(7.37), None),
    el(7, "N", 14.007, false, None, None, None),
    el(8, "O", 15.999, false, None, None, None),
    el(9, "F", 18.998403162, false, None, None, None),
    el(10, "Ne", 20.1797, false, None, None, None),
    el(11, "Na", 22.98976928, false, Some(0.971), Some(1.11), None),
    el(12, "Mg", 24.305, false, Some(1.738), Some(1.51), None),
    el(
        13,
        "Al",
        26.9815384,
        false,
        Some(2.6989),
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
    el(21, "Sc", 44.955907, false, Some(2.989), None, None),
    el(22, "Ti", 47.867, false, Some(4.54), Some(4.85), Some(30.0)),
    el(23, "V", 50.9415, false, Some(6.11), Some(5.31), Some(40.0)),
    el(24, "Cr", 51.9961, false, Some(7.20), Some(4.1), Some(40.0)),
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
        Some(8.902),
        Some(4.44),
        Some(40.0),
    ),
    el(29, "Cu", 63.546, false, Some(8.96), Some(3.49), Some(30.0)),
    el(30, "Zn", 65.38, false, Some(7.133), Some(1.35), None),
    el(31, "Ga", 69.723, false, Some(5.904), Some(2.81), None),
    el(32, "Ge", 72.630, false, Some(5.323), Some(3.85), None),
    el(33, "As", 74.921595, false, Some(5.73), Some(2.96), None),
    el(34, "Se", 78.971, false, Some(4.79), None, None),
    el(35, "Br", 79.904, false, None, None, None),
    el(36, "Kr", 83.798, false, None, None, None),
    el(37, "Rb", 85.4678, false, Some(1.532), None, None),
    el(38, "Sr", 87.62, false, Some(2.54), None, None),
    el(39, "Y", 88.905838, false, Some(4.469), None, None),
    el(40, "Zr", 91.222, false, Some(6.506), Some(6.25), Some(40.0)),
    el(
        41,
        "Nb",
        92.90637,
        false,
        Some(8.57),
        Some(7.57),
        Some(60.0),
    ),
    el(42, "Mo", 95.95, false, Some(10.22), Some(6.82), Some(60.0)),
    el(43, "Tc", 98.0, true, Some(11.50), None, None),
    el(44, "Ru", 101.07, false, Some(12.41), None, None),
    el(45, "Rh", 102.90549, false, Some(12.41), None, None),
    el(46, "Pd", 106.42, false, Some(12.02), Some(3.89), None),
    el(
        47,
        "Ag",
        107.8682,
        false,
        Some(10.50),
        Some(2.95),
        Some(30.0),
    ),
    el(48, "Cd", 112.414, false, Some(8.65), Some(1.16), None),
    el(49, "In", 114.818, false, Some(7.31), Some(2.52), None),
    el(50, "Sn", 118.710, false, Some(7.31), Some(3.14), None),
    el(51, "Sb", 121.760, false, Some(6.691), Some(2.75), None),
    el(52, "Te", 127.60, false, Some(6.24), None, None),
    el(53, "I", 126.90447, false, Some(4.93), None, None),
    el(54, "Xe", 131.293, false, None, None, None),
    el(55, "Cs", 132.90545196, false, Some(1.873), None, None),
    el(56, "Ba", 137.327, false, Some(3.5), None, None),
    el(57, "La", 138.90547, false, Some(6.145), None, None),
    el(58, "Ce", 140.116, false, Some(6.770), None, None),
    el(59, "Pr", 140.90766, false, Some(6.773), None, None),
    el(60, "Nd", 144.242, false, Some(7.008), None, None),
    el(61, "Pm", 145.0, true, Some(7.264), None, None),
    el(62, "Sm", 150.36, false, Some(7.52), None, None),
    el(63, "Eu", 151.964, false, Some(5.244), None, None),
    el(64, "Gd", 157.249, false, Some(7.901), None, None),
    el(65, "Tb", 158.925354, false, Some(8.23), None, None),
    el(66, "Dy", 162.500, false, Some(8.551), None, None),
    el(67, "Ho", 164.930329, false, Some(8.795), None, None),
    el(68, "Er", 167.259, false, Some(9.066), None, None),
    el(69, "Tm", 168.934219, false, Some(9.321), None, None),
    el(70, "Yb", 173.045, false, Some(6.966), None, None),
    el(71, "Lu", 174.96669, false, Some(9.841), None, None),
    el(72, "Hf", 178.486, false, Some(13.31), None, None),
    el(
        73,
        "Ta",
        180.94788,
        false,
        Some(16.654),
        Some(8.1),
        Some(90.0),
    ),
    el(74, "W", 183.84, false, Some(19.3), Some(8.9), Some(90.0)),
    el(75, "Re", 186.207, false, Some(21.02), None, None),
    el(76, "Os", 190.23, false, Some(22.57), None, None),
    el(77, "Ir", 192.217, false, Some(22.42), None, None),
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
    el(82, "Pb", 207.2, false, Some(11.35), Some(2.03), Some(25.0)),
    el(83, "Bi", 208.98040, false, Some(9.747), None, None),
    el(84, "Po", 209.0, true, Some(9.32), None, None),
    el(85, "At", 210.0, true, None, None, None),
    el(86, "Rn", 222.0, true, None, None, None),
    el(87, "Fr", 223.0, true, None, None, None),
    el(88, "Ra", 226.0, true, Some(5.0), None, None),
    el(89, "Ac", 227.0, true, Some(10.07), None, None),
    el(90, "Th", 232.0377, false, Some(11.72), None, None),
    el(91, "Pa", 231.03588, false, Some(15.37), None, None),
    el(92, "U", 238.02891, false, Some(18.95), None, None),
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
    fn mass_number_flag_is_exactly_the_elements_without_standard_weight() {
        // Tc, Pm, Po, At, Rn, Fr, Ra, Ac have no CIAAW standard atomic weight.
        let flagged: Vec<u8> = ELEMENTS
            .iter()
            .filter(|e| e.weight_is_mass_number)
            .map(|e| e.z)
            .collect();
        assert_eq!(flagged, [43, 61, 84, 85, 86, 87, 88, 89]);
        // A mass number is an integer.
        for e in ELEMENTS.iter().filter(|e| e.weight_is_mass_number) {
            assert_eq!(e.atomic_weight, e.atomic_weight.trunc());
        }
    }

    #[test]
    fn symbols_are_unique() {
        for (i, a) in ELEMENTS.iter().enumerate() {
            for b in &ELEMENTS[i + 1..] {
                assert_ne!(a.symbol, b.symbol);
            }
        }
    }

    fn by(symbol: &str) -> &'static Element {
        element_by_symbol(symbol).unwrap()
    }

    /// Standard atomic weights as printed by CIAAW, "Standard Atomic Weights
    /// 2024", <https://www.ciaaw.org/atomic-weights.htm> (read 2026-10-06),
    /// to the digits printed there. For interval elements the value is the
    /// conventional (abridged) one of
    /// <https://www.ciaaw.org/abridged-atomic-weights.htm>. That value need
    /// not lie inside the interval (O: 15.999 against [15.999 03, 15.999 77]);
    /// the 2021 report (Prohaska et al. 2022, text above Table 1) defines its
    /// printed uncertainty as the smallest symmetric one covering the
    /// interval, so the check is `value +- uncertainty` contains the interval.
    #[test]
    fn atomic_weights_match_ciaaw_spot_values() {
        // Single-valued standard weights, value(uncertainty) as printed.
        for (sym, w) in [
            ("P", 30.973_761_998), // 30.973 761 998(5)
            ("Ga", 69.723),        // 69.723(1)
            ("As", 74.921_595),    // 74.921 595(6)
            ("Au", 196.966_570),   // 196.966 570(4)
            // Corrected in #15 (were 44.955908, 88.90584, 91.224, 157.25, 174.9668).
            ("Sc", 44.955_907), // 44.955 907(4)
            ("Y", 88.905_838),  // 88.905 838(2)
            ("Zr", 91.222),     // 91.222(3), CIAAW 2024 revision
            ("Gd", 157.249),    // 157.249(2), CIAAW 2024 revision
            ("Lu", 174.966_69), // 174.966 69(5), CIAAW 2024 revision
        ] {
            assert_eq!(by(sym).atomic_weight, w, "{sym}");
            assert!(!by(sym).weight_is_mass_number, "{sym}");
        }
        // Interval elements: (conventional value +- its abridged uncertainty,
        // standard-weight interval).
        for (sym, w, u, lo, hi) in [
            ("Si", 28.085, 0.001, 28.084, 28.086),
            ("O", 15.999, 0.001, 15.999_03, 15.999_77),
            ("B", 10.81, 0.02, 10.806, 10.821),
            ("Ar", 39.95, 0.16, 39.792, 39.963),
            ("Pb", 207.2, 1.1, 206.14, 207.94),
        ] {
            let e = by(sym);
            assert_eq!(e.atomic_weight, w, "{sym}");
            assert!(w - u <= lo && hi <= w + u, "{sym}");
            assert!(!e.weight_is_mass_number, "{sym}");
        }
    }

    /// Mass numbers of the longest-lived isotopes, CIAAW "Radioactive
    /// elements", <https://www.ciaaw.org/radioactive-elements.htm> (NUBASE2020;
    /// read 2026-10-06). Tc is listed with 97 and 98 as equally long-lived.
    #[test]
    fn mass_number_fallbacks_match_ciaaw() {
        for (sym, a) in [
            ("Tc", 98.0),
            ("Pm", 145.0),
            ("Po", 209.0),
            ("At", 210.0),
            ("Rn", 222.0),
            ("Fr", 223.0),
            ("Ra", 226.0),
            ("Ac", 227.0),
        ] {
            assert_eq!(by(sym).atomic_weight, a, "{sym}");
            assert!(by(sym).weight_is_mass_number, "{sym}");
        }
    }

    /// Cross-check of `E_s` against L. Brewer, LBL-3720 (1975), Table I,
    /// atomization energy to the gaseous ground state at 0 K, kcal/gram-atom
    /// (<https://www.osti.gov/servlets/purl/7187973>). This is not the cited
    /// Kittel edition, which could not be reached; see the module docs.
    /// 1 kcal/mol = 4184 J / (e N_A) eV per atom (thermochemical calorie).
    #[test]
    fn surface_binding_cross_checks_against_brewer_lbl_3720() {
        use crate::constants::{AVOGADRO, ELEMENTARY_CHARGE};
        let kcal_ev = 4184.0 / (ELEMENTARY_CHARGE * AVOGADRO);
        // Agree to the printed digit (half a unit in the last place, 0.005 eV).
        for (sym, kcal) in [("Si", 106.7), ("Ga", 64.8), ("As", 68.2), ("Au", 87.96)] {
            let es = by(sym).default_surface_binding_ev.unwrap();
            assert!(
                (es - kcal * kcal_ev).abs() <= 0.005,
                "{sym}: {es} vs {}",
                kcal * kcal_ev
            );
        }
        // B (beta): Brewer 133 +- 3 kcal = 5.77 +- 0.13 eV; the table's 5.81 is
        // inside that uncertainty but not equal to it. Unresolved until the
        // Kittel edition is seen.
        let b = by("B").default_surface_binding_ev.unwrap();
        assert!((b - 133.0 * kcal_ev).abs() <= 3.0 * kcal_ev, "{b}");
        // O is a gas; P has no entry. Both stay unset by design.
        assert_eq!(by("O").default_surface_binding_ev, None);
        assert_eq!(by("P").default_surface_binding_ev, None);
    }

    /// Densities as printed in the X-Ray Data Booklet, LBNL/PUB-490 Rev. 3
    /// (2009), Table 5-2, pp. 5-5 to 5-8 (<https://xdb.lbl.gov/xdb-new.pdf>,
    /// read 2026-10-07; from the CRC Handbook, 80th ed.), g/cm³, to the
    /// digits printed there. Superscripts give the temperature where it is
    /// not 20 °C.
    #[test]
    fn densities_match_x_ray_data_booklet_spot_values() {
        for (sym, rho) in [
            ("B", 2.34),
            ("P", 1.82),
            ("Ga", 5.904), // 5.904^29.6
            ("Cu", 8.96),
            ("Sc", 2.989), // 2.989^25; reviewer-flagged, unchanged
            // Corrected in #15 (were 5.727, 10.49, 1.93, 22.59, 19.1).
            ("As", 5.73),
            ("Ag", 10.50),
            ("Cs", 1.873),
            ("Os", 22.57),
            ("U", 18.95), // printed "~18.95"
        ] {
            assert_eq!(by(sym).density_g_cm3, Some(rho), "{sym}");
        }
        // Au is printed only as "~19.3"; the stored 19.30 agrees to that
        // precision, its second decimal is not verified.
        let au = by("Au").density_g_cm3.unwrap();
        assert!((au - 19.3).abs() < 0.05, "{au}");
    }

    /// Densities the X-Ray Data Booklet does not pin, checked against other
    /// open primary compilations (read 2026-10-07; see the module docs and
    /// `docs/data-provenance.md`), g/cm³.
    #[test]
    fn densities_match_other_primary_compilations() {
        // NBS Circular 539 Vol. I (1953), p. 33: Au 19.302 at 25 °C from the
        // lattice constant; the booklet prints "~19.3". Agrees to the stored
        // digits.
        let au = by("Au").density_g_cm3.unwrap();
        assert!((au - 19.302).abs() < 0.005, "{au}");
        // NBS Circular 539 Vol. V (1955), p. 20: Cr 7.200 at 25 °C from the
        // lattice constant, the upper end of the booklet's 7.18-7.20.
        // Corrected in #15 (was 7.15, outside that range).
        let cr = by("Cr").density_g_cm3.unwrap();
        assert_eq!(cr, 7.20);
        assert!((7.18..=7.20).contains(&cr));
        // Kirby and Salutsky, NAS-NS 3057 (1964), p. 3: Ra "approximately
        // 5.0". Reviewer-flagged; unchanged.
        assert_eq!(by("Ra").density_g_cm3, Some(5.0));
        // Farr et al., LA-1545 (1953), Table 4: Ac 10.07 (fcc, a = 5.311 Å).
        assert_eq!(by("Ac").density_g_cm3, Some(10.07));
    }

    /// The Ac entry is the crystal density of its fcc cell as LA-1545 prints
    /// it: 4 atoms of mass number 227 in a cube of edge 5.311 Å, rounded to
    /// the printed 10.07 g/cm³ (the report's a carries +- 0.010 Å).
    #[test]
    fn actinium_density_matches_its_lattice_constant() {
        use crate::constants::AVOGADRO;
        let ac = by("Ac");
        let a_cm = 5.311e-8;
        let rho = 4.0 * ac.atomic_weight / (AVOGADRO * a_cm * a_cm * a_cm);
        assert!((rho - ac.density_g_cm3.unwrap()).abs() < 0.01, "{rho}");
    }

    /// Mn is not pinned by any source read so far. The X-Ray Data Booklet
    /// prints the range 7.21-7.44; NBS Monograph 25, Section 17 (1980), p. 50,
    /// prints the alpha-Mn X-ray density "(calculated) 7.475 g/cm³" at 25 °C
    /// from a = 8.9121 Å and Z = 58 atoms per cell (volume printed as
    /// 707.85 Å³). That page is reproduced here from its printed lattice
    /// constant with CODATA 2022 N_A and the CIAAW weight; it lies above the
    /// booklet's range, so the stored upper end, 7.44, is kept and stays
    /// unverified.
    #[test]
    fn manganese_density_is_bounded_but_not_pinned() {
        use crate::constants::AVOGADRO;
        let mn = by("Mn");
        let a_cm: f64 = 8.9121e-8;
        assert!((a_cm * a_cm * a_cm * 1e24 - 707.85).abs() < 0.005);
        let rho_x = 58.0 * mn.atomic_weight / (AVOGADRO * a_cm * a_cm * a_cm);
        assert!((rho_x - 7.475).abs() < 0.0005, "{rho_x}");
        let stored = mn.density_g_cm3.unwrap();
        assert_eq!(stored, 7.44);
        assert!((7.21..=7.44).contains(&stored));
        assert!(rho_x > 7.44);
    }

    /// Si: the CODATA 2022 "molar volume of silicon", V_m(Si) =
    /// 1.205 883 199e-5 m³/mol
    /// (<https://physics.nist.gov/cuu/Constants/Table/allascii.txt>), which
    /// equals N_A a³ / 8 for the listed lattice parameter a = 5.431 020 511e-10
    /// m, gives the crystal density M(Si) / V_m(Si). For every weight in the
    /// CIAAW standard interval [28.084, 28.086] it rounds to the stored 2.329
    /// g/cm³; the X-Ray Data Booklet prints 2.33 at 25 °C.
    #[test]
    fn silicon_density_matches_codata_molar_volume() {
        use crate::constants::AVOGADRO;
        let a = 5.431_020_511e-10; // m
        let v_m = 1.205_883_199e-5; // m^3/mol
        assert!((AVOGADRO * a * a * a / 8.0 / v_m - 1.0).abs() < 1e-9);
        let si = by("Si").density_g_cm3.unwrap();
        assert_eq!(si, 2.329);
        for m_g in [28.084, 28.085, 28.086] {
            // g/mol / (m^3/mol) = g/m^3; 1e-6 converts to g/cm^3.
            let rho = m_g / v_m * 1e-6;
            assert!((rho - si).abs() < 0.0005, "{m_g}: {rho}");
        }
        assert!((si - 2.33).abs() < 0.005);
    }

    /// Defaults that are unset by design stay unset. `E_d` values are not
    /// pinned: ASTM E521 could not be reached, so no populated value is
    /// source-verified.
    #[test]
    fn unset_defaults_are_preserved() {
        // Gases have no solid density.
        for sym in ["H", "N", "O", "F", "Cl", "Ar", "Rn"] {
            assert_eq!(by(sym).density_g_cm3, None, "{sym}");
        }
        // No E_d convention value is recorded for the semiconductor and dopant
        // elements; the user must set one.
        for sym in ["Si", "O", "Ga", "As", "B", "P"] {
            assert_eq!(by(sym).default_displacement_ev, None, "{sym}");
        }
    }
}
