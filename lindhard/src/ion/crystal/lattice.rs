//! Crystal lattices: Bravais vectors, a basis of sites with species, and
//! lattice parameters stated with their temperature.
//!
//! # Cubic structures
//!
//! Both cubic structures are a face-centred cubic Bravais lattice with a
//! two-site basis. The primitive vectors and the basis are taken as printed in
//! M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. Hart and
//! S. Curtarolo, "The AFLOW Library of Crystallographic Prototypes: Part 1",
//! Comput. Mater. Sci. 136, S1 (2017), doi:10.1016/j.commatsci.2017.01.017
//! (read as arXiv:1607.02532v1):
//!
//! * face-centred cubic primitive vectors `a1 = (0, a/2, a/2)`,
//!   `a2 = (a/2, 0, a/2)`, `a3 = (a/2, a/2, 0)` (both entries below);
//! * **diamond** (Strukturbericht A4, prototype `A_cF8_227_a`, space group
//!   Fd-3m, No. 227), p. 616: two sites of one species, Wyckoff 8a, at
//!   lattice coordinates `(1/8, 1/8, 1/8)` and `(7/8, 7/8, 7/8)`; the page
//!   lists Si, Ge and Sn among the elements with this structure;
//! * **zincblende** (B3, prototype `AB_cF8_216_c_a`, space group F-43m,
//!   No. 216), p. 543: Zn at `(0, 0, 0)` (Wyckoff 4a) and S at
//!   `(1/4, 1/4, 1/4)` (4c).
//!
//! For a zincblende compound the first species named in
//! [`Lattice::zincblende`] takes the Zn (4a) site and the second the S (4c)
//! site. Swapping them is the same crystal inverted through the origin; the
//! choice only matters for the polarity of `{111}` faces, which this step
//! does not use.
//!
//! The lattice coordinates `u` of a site give its Cartesian position
//! `u1 a1 + u2 a2 + u3 a3` (Mehl et al., eq. (8)). The primitive cell volume
//! is `V = a1 . (a2 x a3)` (eq. (6)), here `a^3 / 4`, so the atom number
//! density is `(number of basis sites) / V = 8 / a^3` for both structures.
//!
//! # Hexagonal structures
//!
//! Wurtzite and the SiC polytypes 4H and 6H share the space group P6_3mc
//! (No. 186) and the hexagonal primitive vectors of Mehl et al. (2017),
//! pp. 423, 425 and 427:
//!
//! ```text
//! a1 = (a/2, -sqrt(3) a/2, 0),  a2 = (a/2, sqrt(3) a/2, 0),  a3 = (0, 0, c),
//! ```
//!
//! so the crystal-frame `z` axis is the `c` axis `[0001]` and `x` is along
//! `a1 + a2`. The sites are on the Wyckoff positions 2a, `(0, 0, z)` and
//! `(0, 0, 1/2 + z)`, and 2b, `(1/3, 2/3, z)` and `(2/3, 1/3, 1/2 + z)`, in
//! lattice coordinates of `(a1, a2, a3)`:
//!
//! * **wurtzite** (B4, `AB_hP4_186_b_b`), p. 425: Zn on 2b with `z = 0` and S
//!   on 2b with `z = u`. Here the first species named in
//!   [`Lattice::wurtzite`] takes the Zn site and the second the S site.
//!   `u` is the anion-cation bond length along `[0001]` in units of `c`, the
//!   definition of F. Bernardini, V. Fiorentini and D. Vanderbilt, Phys. Rev.
//!   B 56, R10024 (1997), doi:10.1103/PhysRevB.56.R10024, p. 1 (read as
//!   arXiv:cond-mat/9705105v1);
//! * **4H** (B5, `AB_hP8_186_ab_ab`), p. 423, stacking ABAC: Si and C on 2a
//!   and 2b, four parameters `z1..z4`;
//! * **6H** (B6, `AB_hP12_186_ab_a2b`), p. 427, stacking ABCACB: Si and C on
//!   2a and on two 2b sets, six parameters `z1..z6`.
//!
//! The pages name A, B and C the three close-packed column positions of the
//! layers ("one possible stacking (ABAC) for tetrahedral structures. Compare
//! this to Zincblende (ABCABC), Wurtzite (ABABAB), 6H (ABCACB)", p. 423). In
//! lattice coordinates they are A = `(0, 0)` (the 2a column), B =
//! `(1/3, 2/3)` and C = `(2/3, 1/3)` (the two columns of 2b); this matches the
//! printed 4H and 6H sites. [`Lattice::polytype`] builds the ideal crystal of
//! any such sequence; see its docs.
//!
//! The primitive cell volume is `a1 . (a2 x a3) = (sqrt(3)/2) a^2 c`, and the
//! hexagonal cell is primitive, so it is also the conventional cell.
//!
//! # Miller indices
//!
//! Miller indices `(hkl)` and direction indices `[uvw]` refer to the
//! **conventional** cell `(c1, c2, c3)`: the cube `(a x, a y, a z)` for the
//! cubic structures (not the primitive vectors above), and `(a1, a2, a3)` for
//! the hexagonal ones. The reciprocal vectors of a cell are
//! `b_i = 2 pi (c_j x c_k) / V` (Mehl et al., eqs. (9)-(10)). The planes
//! `(hkl)` are orthogonal to the reciprocal lattice vector
//! `g = h b1 + k b2 + l b3` (the definition of Miller indices via the
//! reciprocal lattice, as given in the Wikipedia article "Miller index",
//! revision 1378193236 of 2026-10-03, which cites N. W. Ashcroft and
//! N. D. Mermin, *Solid State Physics* (1976); the book itself was not
//! opened). A direction `[uvw]` is `u c1 + v c2 + w c3`, and by
//! `c_i . b_j = 2 pi delta_ij` (eq. (9)) it lies in the plane `(hkl)` exactly
//! when `h u + k v + l w = 0` (the Weiss zone law, derived here from eq. (9)).
//! For the cubic cell `g` is parallel to `[hkl]`.
//!
//! # Miller-Bravais (four-index) indices
//!
//! For the hexagonal structures, planes and directions may also be given with
//! four indices. The convention, stated once here:
//!
//! * the three basal axes are `a1`, `a2` and `a3' = -(a1 + a2)` (the
//!   Wikipedia "Miller index" article, revision 1378193236, section
//!   "Hexagonal and rhombohedral structures": the `[100]`, `[010]` and
//!   `[-1-10]` directions are equivalent, and the third index is the inverse
//!   intercept on the `[-1-10]` axis), and the fourth axis is `c = a3`;
//! * a plane `(hkil)` must satisfy `h + k + i = 0`, and `h`, `k`, `l` are its
//!   three-index Miller indices (same article: "h, k and l are identical to
//!   the corresponding Miller indices, and i is a redundant index"); its
//!   normal is `g = h b1 + k b2 + l b3` ([`plane_from_hkil`]);
//! * a direction `[uvtw]` is the vector `u a1 + v a2 + t a3' + w c`, with the
//!   constraint `u + v + t = 0` that makes the four indices unique. Putting
//!   `a3' = -(a1 + a2)` gives the three-index direction
//!   `[U V W] = [u - t, v - t, w]` ([`direction_from_uvtw`]; derived here).
//!   The same substitution turns the four-index zone law
//!   `h u + k v + i t + l w = 0` into the three-index one, since
//!   `i = -(h + k)`.
//!
//! So `[0001]` is the `c` axis, `(0, 0, 1)` in the crystal frame, and
//! `[11-20]` is `3 (a1 + a2) = 3 a x`, along `+x`. The article notes that ad
//! hoc four-index schemes for lattice vectors exist (in the electron
//! microscopy literature) that do not follow this rule; this crate uses only
//! the one above. The original statement of the four-index direction symbol
//! (F. C. Frank, Acta Cryst. 18, 862 (1965), doi:10.1107/S0365110X65002116)
//! was not opened.
//!
//! # Lattice parameters
//!
//! A lattice parameter is a measured fact, so every preset cites the record it
//! was read from and states the temperature it refers to where the record
//! does; see [`LatticeConstant`], [`HexagonalConstants`] and the
//! `docs/data-provenance.md` rows. The lattice is rigid: no thermal expansion
//! is applied, and the temperature is carried as metadata of the cited value.

use std::f64::consts::PI;

use super::{add, cross, dot, scale, unit, CrystalError};
use crate::elements::element;
use crate::units::M_PER_ANGSTROM;

/// A cited lattice constant: value, the temperature it refers to, and where
/// it was read.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatticeConstant {
    /// Material the value belongs to (display only).
    pub material: &'static str,
    /// Cubic lattice constant (conventional cell edge), metres.
    pub a_m: f64,
    /// Temperature the value refers to, kelvin.
    pub temperature_k: f64,
    /// Short citation of the record the value was read from.
    pub source: &'static str,
}

/// Cited hexagonal cell parameters: `a`, `c`, the temperature they refer to
/// (`None` when the record read does not state it) and where they were read.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HexagonalConstants {
    /// Material the values belong to (display only).
    pub material: &'static str,
    /// Basal lattice parameter `a`, metres.
    pub a_m: f64,
    /// Lattice parameter `c` (the `[0001]` period), metres.
    pub c_m: f64,
    /// Temperature the values refer to, kelvin, if the record states it.
    pub temperature_k: Option<f64>,
    /// Short citation of the record the values were read from.
    pub source: &'static str,
}

/// Silicon: `a = 5.431 020 511(89) x 10^-10 m`, in vacuum at 22.5 °C
/// (295.65 K).
///
/// CODATA 2022: P. J. Mohr, D. B. Newell, B. N. Taylor and E. Tiesinga,
/// "CODATA Recommended Values of the Fundamental Physical Constants: 2022",
/// Rev. Mod. Phys. 97, 025002 (2025), read as arXiv:2409.03787v1,
/// Table XXXIV "Values of some x-ray-related quantities", p. 56: "lattice
/// parameter of Si (in vacuum, 22.5 °C)", with the footnote that it is the
/// unit cell edge of an ideal single crystal of naturally occurring Si free
/// of impurities and imperfections. The same value is in the NIST CODATA
/// table (<https://physics.nist.gov/cuu/Constants/Table/allascii.txt>).
pub const SI_LATTICE_CONSTANT: LatticeConstant = LatticeConstant {
    material: "Si",
    a_m: 5.431_020_511e-10,
    temperature_k: 295.65,
    source: "CODATA 2022 (Mohr et al., Rev. Mod. Phys. 97, 025002), Table XXXIV, \
             in vacuum at 22.5 °C",
};

/// Germanium: `a = 5.6576 Å` at 25 °C (298.15 K).
///
/// H. E. Swanson and E. Tatge, *Standard X-ray Diffraction Powder Patterns*,
/// NBS Circular 539, Vol. I (1953), section 2.6 "Germanium (Cubic)",
/// pp. 18-19: the NBS lattice constant ("Swanson and Tatge, 5.6576") in the
/// table of unit cells "in angstrom units at 25 °C"
/// (<https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular539v1.pdf>).
pub const GE_LATTICE_CONSTANT: LatticeConstant = LatticeConstant {
    material: "Ge",
    a_m: 5.6576 * M_PER_ANGSTROM,
    temperature_k: 298.15,
    source: "NBS Circular 539, Vol. I (1953), pp. 18-19, at 25 °C",
};

/// Gallium arsenide: `a = 5.652 Å` at 25 °C (298.15 K).
///
/// H. E. Swanson, M. C. Morris, E. H. Evans and L. Ulmer, *Standard X-ray
/// Diffraction Powder Patterns*, NBS Monograph 25, Section 3 (1964),
/// "Gallium Arsenide, GaAs (cubic)", p. 33: "1963 National Bureau of
/// Standards at 25 °C, 5.652" Å, the average of the last five lines; the
/// same page states the zinc sulfide structure, space group F-43m, and
/// "The density of gallium arsenide calculated from the NBS lattice constant
/// is 5.321 g/cm³ at 25 °C"
/// (<https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nbsmonograph25-3.pdf>).
pub const GAAS_LATTICE_CONSTANT: LatticeConstant = LatticeConstant {
    material: "GaAs",
    a_m: 5.652 * M_PER_ANGSTROM,
    temperature_k: 298.15,
    source: "NBS Monograph 25, Section 3 (1964), p. 33, at 25 °C",
};

/// The GaAs density printed with [`GAAS_LATTICE_CONSTANT`] on the same NBS
/// page, g/cm³ at 25 °C. Used to cross-check the lattice number density
/// against the [`crate::material`] module.
pub const GAAS_DENSITY_G_CM3: f64 = 5.321;

/// Cubic (3C) silicon carbide: `a = 4.3596 Å` at 297 K.
///
/// Read from the Ioffe Institute "New Semiconductor Materials" archive,
/// "Basic Parameters of Silicon Carbide (SiC)"
/// (<http://www.ioffe.ru/SVA/NSM/Semicond/SiC/basic.html>, read 2026-10-07):
/// "Lattice constant, 3C-SiC, a = 4.3596 A, 297 K, Debye-Scherrer", citing
/// A. Taylor and R. M. Jones (1960). This is a secondary compilation; the
/// primary paper was **not** opened, see `docs/data-provenance.md`.
pub const SIC_3C_LATTICE_CONSTANT: LatticeConstant = LatticeConstant {
    material: "3C-SiC",
    a_m: 4.3596 * M_PER_ANGSTROM,
    temperature_k: 297.0,
    source: "Ioffe NSM archive (citing Taylor & Jones 1960), at 297 K; \
             primary not opened",
};

/// Wurtzite gallium nitride: `a = 3.189 Å`, `c = 5.178 Å` at 300 K.
///
/// Read from the Ioffe Institute "New Semiconductor Materials" archive,
/// "Basic Parameters of Gallium Nitride (GaN)", section "Basic Parameters for
/// Wurtzite crystal structure"
/// (<http://www.ioffe.ru/SVA/NSM/Semicond/GaN/basic.html>, read 2026-10-08):
/// "Lattice constant, a: 3.189 A, 300 K" and "Lattice constant, c:
/// 5.178 A, 300 K", both credited to Qian et al. (1996) (W. Qian,
/// M. Skowronski and G. R. Rohrer, in *III-Nitride, SiC, and Diamond
/// Materials for Electronic Devices*, Materials Research Society Symposium
/// Proceedings 423 (1996), 475-486, per the archive's reference list). Both values are taken from that one row so
/// that `a` and `c` belong to the same measurement; the page also lists
/// `c = 5.186 A` at 300 K from a handbook (Bougrov et al. 2001). This is a
/// secondary compilation and the primary was **not** opened; the archive's
/// GaN density and atom count do not match these parameters, see
/// `docs/data-provenance.md` (open questions).
pub const GAN_LATTICE_CONSTANTS: HexagonalConstants = HexagonalConstants {
    material: "GaN (wurtzite)",
    a_m: 3.189 * M_PER_ANGSTROM,
    c_m: 5.178 * M_PER_ANGSTROM,
    temperature_k: Some(300.0),
    source: "Ioffe NSM archive (citing Qian et al. 1996), at 300 K; primary not opened",
};

/// The wurtzite GaN density printed on the same Ioffe page as
/// [`GAN_LATTICE_CONSTANTS`], "Density 6.15 g cm-3, 300 K" (no reference
/// given). It is **not** consistent with the archive's own `a` and `c`
/// (about 0.85 % higher than `4 M(GaN) / 2 / (N_A V)`); the test pins the
/// difference and the gap is recorded in `docs/data-provenance.md`.
pub const GAN_DENSITY_G_CM3: f64 = 6.15;

/// 4H silicon carbide: `a = 3.08051 Å`, `c = 10.0848 Å`; temperature not
/// stated in the record read.
///
/// A. Bauer, P. Reischauer, J. Kräusslich, N. Schell, W. Matz and K. Goetz,
/// "Structure refinement of the silicon carbide polytypes 4H and 6H:
/// unambiguous determination of the refinement parameters", Acta Cryst. A 57,
/// 60 (2001), doi:10.1107/S0108767300012915, as reproduced in Mehl et al.
/// (2017), prototype `AB_hP8_186_ab_ab` ("Moissanite-4H SiC (B5)"),
/// pp. 423-424 and its CIF on p. 733: `_cell_length_a 3.08051`,
/// `_cell_length_c 10.08480`. The primary paper is closed access and was
/// **not** opened, so its measurement temperature is not known here.
pub const SIC_4H_LATTICE_CONSTANTS: HexagonalConstants = HexagonalConstants {
    material: "4H-SiC",
    a_m: 3.08051 * M_PER_ANGSTROM,
    c_m: 10.0848 * M_PER_ANGSTROM,
    temperature_k: None,
    source: "Bauer et al., Acta Cryst. A 57, 60 (2001), via AFLOW AB_hP8_186_ab_ab \
             (Mehl et al. 2017, pp. 423, 733); temperature not stated",
};

/// The Wyckoff parameters `z1..z4` of 4H-SiC as printed with
/// [`SIC_4H_LATTICE_CONSTANTS`] (Mehl et al. 2017, p. 733,
/// `_aflow_params_values`): Si I (2a) `0.0`, C I (2a) `0.18784`, Si II (2b)
/// `0.24982`, C II (2b) `0.43671`, from Bauer et al. (2001).
pub const SIC_4H_Z: [f64; 4] = [0.0, 0.18784, 0.24982, 0.43671];

/// 6H silicon carbide: `a = 3.08129 Å`, `c = 15.11976 Å`; temperature not
/// stated in the record read.
///
/// Bauer et al. (2001) (see [`SIC_4H_LATTICE_CONSTANTS`]) as reproduced in
/// Mehl et al. (2017), prototype `AB_hP12_186_ab_a2b` ("Moissanite-6H SiC
/// (B6)"), pp. 427-428 and its CIF on p. 734: `_cell_length_a 3.08129`,
/// `_cell_length_c 15.11976`. Primary not opened.
pub const SIC_6H_LATTICE_CONSTANTS: HexagonalConstants = HexagonalConstants {
    material: "6H-SiC",
    a_m: 3.08129 * M_PER_ANGSTROM,
    c_m: 15.11976 * M_PER_ANGSTROM,
    temperature_k: None,
    source: "Bauer et al., Acta Cryst. A 57, 60 (2001), via AFLOW AB_hP12_186_ab_a2b \
             (Mehl et al. 2017, pp. 427, 734); temperature not stated",
};

/// The Wyckoff parameters `z1..z6` of 6H-SiC as printed with
/// [`SIC_6H_LATTICE_CONSTANTS`] (Mehl et al. 2017, p. 734): Si I (2a) `0.0`,
/// C I (2a) `0.1254`, Si II (2b) `0.16675`, C II (2b) `0.29215`, Si III (2b)
/// `0.8335`, C III (2b) `-0.0415`, from Bauer et al. (2001).
pub const SIC_6H_Z: [f64; 6] = [0.0, 0.1254, 0.16675, 0.29215, 0.8335, -0.0415];

/// 4H-SiC density, "3.211 g cm-3, 300 K", credited to Gomes de Mesquita
/// (1967) by the Ioffe NSM archive SiC basic-parameters page (read
/// 2026-10-08; primary not opened). Used in the number-density test.
pub const SIC_4H_DENSITY_G_CM3: f64 = 3.211;

/// 6H-SiC density, "3.21 g cm-3, 300 K", credited to Harris et al. (1995b)
/// by the same Ioffe page (primary not opened). Used in the number-density
/// test.
pub const SIC_6H_DENSITY_G_CM3: f64 = 3.21;

/// The wurtzite `u` at which the bond along `[0001]` and the three other
/// bonds of each atom have the same length: `u = 1/4 + a^2 / (3 c^2)`.
///
/// Derived here from the wurtzite sites above: a cation on column B at
/// height 0 has its anion neighbour at `u c` straight above it, and three
/// anions of the C layer at height `(u - 1/2) c`, displaced in the plane by
/// the B-C column spacing `a / sqrt(3)`. Equal lengths,
/// `(u c)^2 = a^2/3 + (1/2 - u)^2 c^2`, give the formula. At the ideal ratio
/// `c/a = sqrt(8/3)` it is `3/8`.
///
/// This is a **geometric model value**, not a measurement. It is what
/// [`Lattice::gallium_nitride`] uses until a measured `u` is sourced.
pub fn equal_bond_wurtzite_u(a_m: f64, c_m: f64) -> f64 {
    0.25 + a_m * a_m / (3.0 * c_m * c_m)
}

/// Which structure a [`Lattice`] has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Structure {
    /// Diamond (A4, Fd-3m): one species on both fcc sublattices.
    Diamond,
    /// Zincblende (B3, F-43m): two species, one per fcc sublattice.
    Zincblende,
    /// Wurtzite (B4, P6_3mc): the two-layer hexagonal stacking, 2H.
    Wurtzite,
    /// A hexagonal tetrahedral polytype with `layers` cation-anion bilayers
    /// per `c` (4 for 4H, 6 for 6H).
    Polytype {
        /// Bilayers per `c` period.
        layers: u8,
    },
}

impl Structure {
    /// Whether the structure has the hexagonal cell (`a1`, `a2`, `c`).
    pub fn is_hexagonal(self) -> bool {
        matches!(self, Self::Wurtzite | Self::Polytype { .. })
    }
}

/// One basis site: species and lattice (fractional) coordinates with respect
/// to the primitive vectors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Site {
    /// Atomic number of the atom on this site.
    pub z: u8,
    /// Lattice coordinates `(u1, u2, u3)`: the position is
    /// `u1 a1 + u2 a2 + u3 a3`.
    pub fractional: [f64; 3],
}

/// A crystal: a Bravais lattice with a basis. Cubic (face-centred cubic
/// Bravais lattice) or hexagonal; see the module docs.
#[derive(Debug, Clone, PartialEq)]
pub struct Lattice {
    structure: Structure,
    a_m: f64,
    c_m: f64,
    temperature_k: Option<f64>,
    primitive: [[f64; 3]; 3],
    basis: Vec<Site>,
}

fn check_z(z: u8) -> Result<u8, CrystalError> {
    element(z)
        .map(|_| z)
        .ok_or(CrystalError::UnknownAtomicNumber(z))
}

fn check_length(x: f64) -> Result<f64, CrystalError> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(CrystalError::InvalidLatticeConstant(x))
    }
}

fn check_temperature(t: Option<f64>) -> Result<Option<f64>, CrystalError> {
    match t {
        Some(t) if !(t.is_finite() && t >= 0.0) => Err(CrystalError::InvalidTemperature(t)),
        _ => Ok(t),
    }
}

/// Reciprocal vectors `b_i = 2 pi (c_j x c_k) / V` of the cell `c`
/// (Mehl et al. 2017, eq. (10)).
pub fn reciprocal_vectors(c: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let v = dot(c[0], cross(c[1], c[2]));
    let f = 2.0 * PI / v;
    [
        scale(cross(c[1], c[2]), f),
        scale(cross(c[2], c[0]), f),
        scale(cross(c[0], c[1]), f),
    ]
}

/// The hexagonal primitive vectors `a1 = (a/2, -sqrt(3) a/2, 0)`,
/// `a2 = (a/2, sqrt(3) a/2, 0)`, `a3 = (0, 0, c)` (Mehl et al. 2017,
/// pp. 423, 425, 427).
pub fn hexagonal_vectors(a_m: f64, c_m: f64) -> [[f64; 3]; 3] {
    let h = 0.5 * a_m;
    let s = 0.5 * 3.0f64.sqrt() * a_m;
    [[h, -s, 0.0], [h, s, 0.0], [0.0, 0.0, c_m]]
}

/// Three-index Miller indices `(hkl)` of the Miller-Bravais plane `(hkil)`.
///
/// # Errors
/// [`CrystalError::InvalidMillerBravais`] unless `h + k + i = 0`.
pub fn plane_from_hkil(hkil: [i32; 4]) -> Result<[i32; 3], CrystalError> {
    let [h, k, i, l] = hkil;
    if i64::from(h) + i64::from(k) + i64::from(i) != 0 {
        return Err(CrystalError::InvalidMillerBravais {
            indices: hkil,
            why: "a plane (hkil) needs h + k + i = 0",
        });
    }
    Ok([h, k, l])
}

/// Three-index direction `[U V W] = [u - t, v - t, w]` of the Miller-Bravais
/// direction `[uvtw]` (see the module docs for the convention).
///
/// # Errors
/// [`CrystalError::InvalidMillerBravais`] unless `u + v + t = 0`, or if an
/// index overflows `i32`.
pub fn direction_from_uvtw(uvtw: [i32; 4]) -> Result<[i32; 3], CrystalError> {
    let [u, v, t, w] = uvtw;
    if i64::from(u) + i64::from(v) + i64::from(t) != 0 {
        return Err(CrystalError::InvalidMillerBravais {
            indices: uvtw,
            why: "a direction [uvtw] needs u + v + t = 0",
        });
    }
    match (u.checked_sub(t), v.checked_sub(t)) {
        (Some(big_u), Some(big_v)) => Ok([big_u, big_v, w]),
        _ => Err(CrystalError::InvalidMillerBravais {
            indices: uvtw,
            why: "index overflow",
        }),
    }
}

fn combine(c: [[f64; 3]; 3], n: [f64; 3]) -> [f64; 3] {
    add(add(scale(c[0], n[0]), scale(c[1], n[1])), scale(c[2], n[2]))
}

fn as_f64(n: [i32; 3]) -> [f64; 3] {
    [f64::from(n[0]), f64::from(n[1]), f64::from(n[2])]
}

/// Wrap a fraction into `[0, 1)`. Values within `1e-12` below 1 (rounding of
/// `1/3 + 2/3` and the like) become 0, so that every site has a single
/// representative.
fn wrap_unit(x: f64) -> f64 {
    let w = x - x.floor();
    if w >= 1.0 - 1e-12 {
        0.0
    } else {
        w
    }
}

/// Lattice coordinates `(f1, f2)` of the close-packed columns A, B, C.
fn column(letter: char) -> Option<[f64; 2]> {
    match letter {
        'A' => Some([0.0, 0.0]),
        'B' => Some([1.0 / 3.0, 2.0 / 3.0]),
        'C' => Some([2.0 / 3.0, 1.0 / 3.0]),
        _ => None,
    }
}

/// Wyckoff 2a of P6_3mc: `(0, 0, z)`, `(0, 0, 1/2 + z)`.
fn wyckoff_2a(z_atom: u8, z: f64) -> [(u8, [f64; 3]); 2] {
    [(z_atom, [0.0, 0.0, z]), (z_atom, [0.0, 0.0, 0.5 + z])]
}

/// Wyckoff 2b of P6_3mc: `(1/3, 2/3, z)`, `(2/3, 1/3, 1/2 + z)`.
fn wyckoff_2b(z_atom: u8, z: f64) -> [(u8, [f64; 3]); 2] {
    [
        (z_atom, [1.0 / 3.0, 2.0 / 3.0, z]),
        (z_atom, [2.0 / 3.0, 1.0 / 3.0, 0.5 + z]),
    ]
}

fn check_u(u: f64) -> Result<f64, CrystalError> {
    if u.is_finite() && u > 0.0 && u < 0.5 {
        Ok(u)
    } else {
        Err(CrystalError::InvalidInternalParameter(u))
    }
}

impl Lattice {
    fn build(
        structure: Structure,
        sites: [(u8, [f64; 3]); 2],
        a_m: f64,
        temperature_k: f64,
    ) -> Result<Self, CrystalError> {
        check_length(a_m)?;
        check_temperature(Some(temperature_k))?;
        let h = 0.5 * a_m;
        // Face-centred cubic primitive vectors, Mehl et al. (2017),
        // pp. 543 and 616.
        let primitive = [[0.0, h, h], [h, 0.0, h], [h, h, 0.0]];
        let basis = sites
            .iter()
            .map(|&(z, fractional)| {
                Ok(Site {
                    z: check_z(z)?,
                    fractional,
                })
            })
            .collect::<Result<Vec<_>, CrystalError>>()?;
        Ok(Self {
            structure,
            a_m,
            c_m: a_m,
            temperature_k: Some(temperature_k),
            primitive,
            basis,
        })
    }

    /// A hexagonal lattice on [`hexagonal_vectors`] with the given sites
    /// (lattice coordinates, wrapped into `[0, 1)`).
    fn build_hexagonal(
        structure: Structure,
        sites: &[(u8, [f64; 3])],
        a_m: f64,
        c_m: f64,
        temperature_k: Option<f64>,
    ) -> Result<Self, CrystalError> {
        check_length(a_m)?;
        check_length(c_m)?;
        check_temperature(temperature_k)?;
        let basis = sites
            .iter()
            .map(|&(z, f)| {
                Ok(Site {
                    z: check_z(z)?,
                    fractional: f.map(wrap_unit),
                })
            })
            .collect::<Result<Vec<_>, CrystalError>>()?;
        Ok(Self {
            structure,
            a_m,
            c_m,
            temperature_k,
            primitive: hexagonal_vectors(a_m, c_m),
            basis,
        })
    }

    /// Diamond structure of element `z` with lattice constant `a_m` (metres)
    /// stated at `temperature_k`. Basis as printed by Mehl et al. (2017),
    /// p. 616: `(1/8, 1/8, 1/8)` and `(7/8, 7/8, 7/8)`.
    pub fn diamond(z: u8, a_m: f64, temperature_k: f64) -> Result<Self, CrystalError> {
        Self::build(
            Structure::Diamond,
            [(z, [0.125; 3]), (z, [0.875; 3])],
            a_m,
            temperature_k,
        )
    }

    /// Zincblende structure with `z_a` on the Zn (4a) site `(0, 0, 0)` and
    /// `z_b` on the S (4c) site `(1/4, 1/4, 1/4)` (Mehl et al. 2017, p. 543).
    pub fn zincblende(
        z_a: u8,
        z_b: u8,
        a_m: f64,
        temperature_k: f64,
    ) -> Result<Self, CrystalError> {
        Self::build(
            Structure::Zincblende,
            [(z_a, [0.0; 3]), (z_b, [0.25; 3])],
            a_m,
            temperature_k,
        )
    }

    /// Wurtzite with `z_a` on the Zn sites and `z_b` on the S sites (Mehl et
    /// al. 2017, p. 425, Wyckoff 2b): `z_a` at `(1/3, 2/3, 0)` and
    /// `(2/3, 1/3, 1/2)`, `z_b` at `(1/3, 2/3, u)` and `(2/3, 1/3, 1/2 + u)`,
    /// in that order. `a_m`, `c_m` in metres; `temperature_k` is the
    /// temperature the parameters refer to, if known.
    ///
    /// # Errors
    /// Invalid lengths or temperature, an unknown atomic number, or `u`
    /// outside `(0, 1/2)` ([`CrystalError::InvalidInternalParameter`]).
    pub fn wurtzite(
        z_a: u8,
        z_b: u8,
        a_m: f64,
        c_m: f64,
        u: f64,
        temperature_k: Option<f64>,
    ) -> Result<Self, CrystalError> {
        let u = check_u(u)?;
        let [a1, a2] = wyckoff_2b(z_a, 0.0);
        let [b1, b2] = wyckoff_2b(z_b, u);
        Self::build_hexagonal(
            Structure::Wurtzite,
            &[a1, a2, b1, b2],
            a_m,
            c_m,
            temperature_k,
        )
    }

    /// The ideal tetrahedral polytype with the close-packed stacking
    /// `stacking` (letters `A`, `B`, `C`, one per bilayer, e.g. `"ABCB"` for
    /// 4H or `"ABCACB"` for 6H), cation `z_a` and anion `z_b`.
    ///
    /// With `n` letters, bilayer `k` (from 0) has its cation at the column of
    /// letter `k` and height `k c / n`, and its anion straight above it at
    /// height `k c / n + 2 u c / n`. So `u` is the bond length along
    /// `[0001]` in units of the two-bilayer height `2 c / n`: the wurtzite `u`
    /// for `n = 2`, and `3/8` for the ideal tetrahedron at any `n` (the bond
    /// is three quarters of the bilayer spacing, as in the ideal wurtzite).
    /// This generalisation is ours; the measured 4H and 6H presets use the
    /// printed sites instead ([`Self::silicon_carbide_4h`],
    /// [`Self::silicon_carbide_6h`]). Sites are listed per bilayer, cation
    /// then anion, in stacking order.
    ///
    /// The Ramsdell number is `n` (Mehl et al. 2017, p. 423: "The 4H refers to
    /// the fact that there are 4 CSi dimers in a hexagonal unit cell"). The
    /// stacking `"BC"` with the same `u` gives the sites of
    /// [`Self::wurtzite`] (listed in a different order).
    ///
    /// # Errors
    /// [`CrystalError::InvalidStacking`] for a sequence shorter than 2 or
    /// longer than 255 letters, a letter other than `A`, `B`, `C`, or two
    /// equal neighbours (cyclically: the last letter is followed by the
    /// first); [`CrystalError::InvalidInternalParameter`] for `u` outside
    /// `(0, 1/2)`; and the errors of [`Self::wurtzite`].
    pub fn polytype(
        z_a: u8,
        z_b: u8,
        stacking: &str,
        a_m: f64,
        c_m: f64,
        u: f64,
        temperature_k: Option<f64>,
    ) -> Result<Self, CrystalError> {
        let bad = |why: &'static str| CrystalError::InvalidStacking {
            stacking: stacking.to_string(),
            why,
        };
        let letters: Vec<char> = stacking.chars().collect();
        let n = letters.len();
        if n < 2 {
            return Err(bad("needs at least two layers"));
        }
        let layers = u8::try_from(n).map_err(|_| bad("at most 255 layers"))?;
        let columns = letters
            .iter()
            .map(|&l| column(l).ok_or_else(|| bad("letters must be A, B or C")))
            .collect::<Result<Vec<_>, _>>()?;
        if (0..n).any(|k| letters[k] == letters[(k + 1) % n]) {
            return Err(bad("adjacent layers (cyclically) must differ"));
        }
        let u = check_u(u)?;
        let nf = n as f64;
        let mut sites = Vec::with_capacity(2 * n);
        for (k, xy) in columns.iter().enumerate() {
            let z = k as f64 / nf;
            sites.push((z_a, [xy[0], xy[1], z]));
            sites.push((z_b, [xy[0], xy[1], z + 2.0 * u / nf]));
        }
        Self::build_hexagonal(
            Structure::Polytype { layers },
            &sites,
            a_m,
            c_m,
            temperature_k,
        )
    }

    /// Diamond-structure Si with [`SI_LATTICE_CONSTANT`].
    pub fn silicon() -> Self {
        let c = SI_LATTICE_CONSTANT;
        Self::diamond(14, c.a_m, c.temperature_k).expect("valid preset")
    }

    /// Diamond-structure Ge with [`GE_LATTICE_CONSTANT`].
    pub fn germanium() -> Self {
        let c = GE_LATTICE_CONSTANT;
        Self::diamond(32, c.a_m, c.temperature_k).expect("valid preset")
    }

    /// Zincblende GaAs (Ga on 4a, As on 4c) with [`GAAS_LATTICE_CONSTANT`].
    pub fn gallium_arsenide() -> Self {
        let c = GAAS_LATTICE_CONSTANT;
        Self::zincblende(31, 33, c.a_m, c.temperature_k).expect("valid preset")
    }

    /// Zincblende 3C-SiC (Si on 4a, C on 4c) with
    /// [`SIC_3C_LATTICE_CONSTANT`].
    pub fn silicon_carbide_3c() -> Self {
        let c = SIC_3C_LATTICE_CONSTANT;
        Self::zincblende(14, 6, c.a_m, c.temperature_k).expect("valid preset")
    }

    /// Wurtzite GaN (Ga on the Zn sites, N on the S sites) with
    /// [`GAN_LATTICE_CONSTANTS`].
    ///
    /// **Placeholder `u`.** No measured `u` for GaN could be opened (see
    /// `docs/data-provenance.md`, open questions), so this preset uses the
    /// geometric [`equal_bond_wurtzite_u`] of the cited `a` and `c`, not a
    /// measurement. Build with [`Self::wurtzite`] to supply a measured `u`.
    pub fn gallium_nitride() -> Self {
        let c = GAN_LATTICE_CONSTANTS;
        let u = equal_bond_wurtzite_u(c.a_m, c.c_m);
        Self::wurtzite(31, 7, c.a_m, c.c_m, u, c.temperature_k).expect("valid preset")
    }

    /// 4H-SiC with [`SIC_4H_LATTICE_CONSTANTS`] and the measured sites
    /// [`SIC_4H_Z`] of Bauer et al. (2001) as printed by Mehl et al. (2017),
    /// p. 423: Si I on 2a (`z1`), C I on 2a (`z2`), Si II on 2b (`z3`),
    /// C II on 2b (`z4`), in that order. The stacking is ABAC, which is ABCB
    /// relabelled.
    pub fn silicon_carbide_4h() -> Self {
        let c = SIC_4H_LATTICE_CONSTANTS;
        let z = SIC_4H_Z;
        let s: Vec<(u8, [f64; 3])> = [
            wyckoff_2a(14, z[0]),
            wyckoff_2a(6, z[1]),
            wyckoff_2b(14, z[2]),
            wyckoff_2b(6, z[3]),
        ]
        .concat();
        Self::build_hexagonal(
            Structure::Polytype { layers: 4 },
            &s,
            c.a_m,
            c.c_m,
            c.temperature_k,
        )
        .expect("valid preset")
    }

    /// 6H-SiC with [`SIC_6H_LATTICE_CONSTANTS`] and the measured sites
    /// [`SIC_6H_Z`] of Bauer et al. (2001) as printed by Mehl et al. (2017),
    /// p. 427: Si I, C I on 2a (`z1`, `z2`); Si II, C II, Si III, C III on 2b
    /// (`z3..z6`), in that order. The stacking is ABCACB.
    pub fn silicon_carbide_6h() -> Self {
        let c = SIC_6H_LATTICE_CONSTANTS;
        let z = SIC_6H_Z;
        let s: Vec<(u8, [f64; 3])> = [
            wyckoff_2a(14, z[0]),
            wyckoff_2a(6, z[1]),
            wyckoff_2b(14, z[2]),
            wyckoff_2b(6, z[3]),
            wyckoff_2b(14, z[4]),
            wyckoff_2b(6, z[5]),
        ]
        .concat();
        Self::build_hexagonal(
            Structure::Polytype { layers: 6 },
            &s,
            c.a_m,
            c.c_m,
            c.temperature_k,
        )
        .expect("valid preset")
    }

    /// The structure type.
    pub fn structure(&self) -> Structure {
        self.structure
    }

    /// Whether the lattice is hexagonal (wurtzite or a polytype).
    pub fn is_hexagonal(&self) -> bool {
        self.structure.is_hexagonal()
    }

    /// Lattice constant `a`, metres: the cubic cell edge, or the hexagonal
    /// basal parameter.
    pub fn lattice_constant(&self) -> f64 {
        self.a_m
    }

    /// Lattice constant `c`, metres: the hexagonal `[0001]` period, or the
    /// cubic cell edge (equal to [`Self::lattice_constant`]).
    pub fn lattice_constant_c(&self) -> f64 {
        self.c_m
    }

    /// Temperature the lattice parameters refer to, kelvin, or `None` if the
    /// cited record does not state it.
    pub fn temperature_k(&self) -> Option<f64> {
        self.temperature_k
    }

    /// The primitive (Bravais) vectors `a1, a2, a3`, metres, in the crystal
    /// frame: face-centred cubic, or hexagonal ([`hexagonal_vectors`]).
    pub fn primitive_vectors(&self) -> [[f64; 3]; 3] {
        self.primitive
    }

    /// The conventional cell vectors, metres: the cube `a x, a y, a z`, or
    /// for a hexagonal lattice the primitive `a1, a2, c`. Miller and
    /// direction indices refer to this cell.
    pub fn conventional_vectors(&self) -> [[f64; 3]; 3] {
        if self.is_hexagonal() {
            return self.primitive;
        }
        let a = self.a_m;
        [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]
    }

    /// The basis sites, in the order the constructor lists them.
    pub fn basis(&self) -> &[Site] {
        &self.basis
    }

    /// Cartesian positions (metres, crystal frame) of the basis sites in the
    /// primitive cell, with their atomic numbers.
    pub fn basis_positions(&self) -> Vec<(u8, [f64; 3])> {
        self.basis
            .iter()
            .map(|s| (s.z, combine(self.primitive, s.fractional)))
            .collect()
    }

    /// The atoms of one conventional cell, as fractional coordinates of that
    /// cell wrapped into `[0, 1)`.
    ///
    /// Cubic: each basis site plus the three face-centring translations
    /// `(0, 1/2, 1/2)`, `(1/2, 0, 1/2)`, `(1/2, 1/2, 0)`; eight atoms for both
    /// structures, in basis order then translation order. Hexagonal: the
    /// basis itself (the hexagonal cell is primitive), in basis order.
    pub fn conventional_cell_sites(&self) -> Vec<(u8, [f64; 3])> {
        if self.is_hexagonal() {
            return self.basis.iter().map(|s| (s.z, s.fractional)).collect();
        }
        const T: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ];
        let mut out = Vec::with_capacity(4 * self.basis.len());
        for s in &self.basis {
            // Cartesian position in units of a.
            let p = scale(combine(self.primitive, s.fractional), 1.0 / self.a_m);
            for t in T {
                let q = add(p, t);
                out.push((s.z, q.map(|x| x - x.floor())));
            }
        }
        out
    }

    /// A rectangular cell aligned with the crystal axes that tiles the
    /// crystal: its edge lengths (metres, along `x`, `y`, `z`) and its atoms
    /// as fractions of those edges, wrapped into `[0, 1)`. The lattice
    /// neighbour search ([`super::search`]) walks these cells.
    ///
    /// Cubic: the conventional cube, edges `(a, a, a)`, and
    /// [`Self::conventional_cell_sites`]. Hexagonal: the orthohexagonal cell
    /// spanned by `a1 + a2 = (a, 0, 0)`, `a2 - a1 = (0, sqrt(3) a, 0)` and
    /// `a3 = (0, 0, c)`, which are lattice vectors spanning a sublattice of
    /// index 2 (determinant of `[[1, 1, 0], [-1, 1, 0], [0, 0, 1]]`), so the
    /// cell holds every basis site twice: at its own position and shifted by
    /// `a2`, which is not in the sublattice. A site at lattice coordinates
    /// `(f1, f2, f3)` has the fractions `((f1 + f2)/2, (f2 - f1)/2, f3)`
    /// (derived here from the vectors above). Basis order, then the shift.
    pub fn orthogonal_cell(&self) -> ([f64; 3], Vec<(u8, [f64; 3])>) {
        if !self.is_hexagonal() {
            return ([self.a_m; 3], self.conventional_cell_sites());
        }
        let edges = [self.a_m, 3.0f64.sqrt() * self.a_m, self.c_m];
        let mut out = Vec::with_capacity(2 * self.basis.len());
        for s in &self.basis {
            let [f1, f2, f3] = s.fractional;
            for shift in [0.0, 1.0] {
                let g2 = f2 + shift;
                out.push((
                    s.z,
                    [
                        wrap_unit(0.5 * (f1 + g2)),
                        wrap_unit(0.5 * (g2 - f1)),
                        wrap_unit(f3),
                    ],
                ));
            }
        }
        (edges, out)
    }

    /// Primitive cell volume `a1 . (a2 x a3)`, m³ (Mehl et al. 2017, eq. (6)):
    /// `a³/4` (cubic) or `(sqrt(3)/2) a² c` (hexagonal).
    pub fn primitive_cell_volume(&self) -> f64 {
        let p = self.primitive;
        dot(p[0], cross(p[1], p[2]))
    }

    /// Atom number density, atoms per m³: basis sites per primitive cell
    /// volume.
    pub fn atom_number_density(&self) -> f64 {
        self.basis.len() as f64 / self.primitive_cell_volume()
    }

    /// Number density of element `z`, atoms per m³ (zero if absent).
    pub fn number_density_of(&self, z: u8) -> f64 {
        let n = self.basis.iter().filter(|s| s.z == z).count();
        n as f64 / self.primitive_cell_volume()
    }

    /// Unit normal (crystal frame) of the lattice planes `(hkl)`, indexed on
    /// the conventional cell: `g = h b1 + k b2 + l b3` normalised. For a
    /// cubic cell this is `[hkl] / |[hkl]|`.
    pub fn plane_normal(&self, hkl: [i32; 3]) -> Result<[f64; 3], CrystalError> {
        if hkl == [0; 3] {
            return Err(CrystalError::ZeroIndex(hkl));
        }
        let b = reciprocal_vectors(self.conventional_vectors());
        Ok(unit(combine(b, as_f64(hkl))))
    }

    /// Unit vector (crystal frame) along the lattice direction `[uvw]`,
    /// indexed on the conventional cell: `u c1 + v c2 + w c3` normalised.
    pub fn direction(&self, uvw: [i32; 3]) -> Result<[f64; 3], CrystalError> {
        if uvw == [0; 3] {
            return Err(CrystalError::ZeroIndex(uvw));
        }
        Ok(unit(combine(self.conventional_vectors(), as_f64(uvw))))
    }

    /// Unit normal of the Miller-Bravais planes `(hkil)` of a hexagonal
    /// lattice: [`Self::plane_normal`] of [`plane_from_hkil`].
    ///
    /// # Errors
    /// [`CrystalError::NotHexagonal`] for a cubic lattice, and the errors of
    /// the two functions named.
    pub fn plane_normal_hkil(&self, hkil: [i32; 4]) -> Result<[f64; 3], CrystalError> {
        if !self.is_hexagonal() {
            return Err(CrystalError::NotHexagonal);
        }
        self.plane_normal(plane_from_hkil(hkil)?)
    }

    /// Unit vector along the Miller-Bravais direction `[uvtw]` of a hexagonal
    /// lattice: [`Self::direction`] of [`direction_from_uvtw`].
    ///
    /// # Errors
    /// [`CrystalError::NotHexagonal`] for a cubic lattice, and the errors of
    /// the two functions named.
    pub fn direction_uvtw(&self, uvtw: [i32; 4]) -> Result<[f64; 3], CrystalError> {
        if !self.is_hexagonal() {
            return Err(CrystalError::NotHexagonal);
        }
        self.direction(direction_from_uvtw(uvtw)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_carry_their_cited_constants() {
        let si = Lattice::silicon();
        assert_eq!(si.lattice_constant(), 5.431_020_511e-10);
        assert_eq!(si.lattice_constant_c(), 5.431_020_511e-10);
        assert_eq!(si.temperature_k(), Some(295.65));
        assert_eq!(si.structure(), Structure::Diamond);
        assert!(!si.is_hexagonal());
        let gaas = Lattice::gallium_arsenide();
        assert_eq!(gaas.structure(), Structure::Zincblende);
        assert_eq!(gaas.basis()[0].z, 31);
        assert_eq!(gaas.basis()[1].z, 33);
        assert_eq!(Lattice::germanium().basis()[0].z, 32);
        let sic = Lattice::silicon_carbide_3c();
        assert_eq!((sic.basis()[0].z, sic.basis()[1].z), (14, 6));
        assert_eq!(sic.temperature_k(), Some(297.0));

        let gan = Lattice::gallium_nitride();
        assert_eq!(gan.structure(), Structure::Wurtzite);
        assert_eq!(gan.lattice_constant(), 3.189 * M_PER_ANGSTROM);
        assert_eq!(gan.lattice_constant_c(), 5.178 * M_PER_ANGSTROM);
        assert_eq!(gan.temperature_k(), Some(300.0));
        let z: Vec<u8> = gan.basis().iter().map(|s| s.z).collect();
        assert_eq!(z, [31, 31, 7, 7]);
        let h4 = Lattice::silicon_carbide_4h();
        assert_eq!(h4.structure(), Structure::Polytype { layers: 4 });
        assert_eq!(h4.basis().len(), 8);
        assert_eq!(h4.temperature_k(), None);
        let h6 = Lattice::silicon_carbide_6h();
        assert_eq!(h6.structure(), Structure::Polytype { layers: 6 });
        assert_eq!(h6.basis().len(), 12);
        assert!(h6.is_hexagonal());
    }

    #[test]
    fn primitive_volume_is_a_quarter_cube() {
        let si = Lattice::silicon();
        let a = si.lattice_constant();
        assert!((si.primitive_cell_volume() / (a * a * a) - 0.25).abs() < 1e-14);
        assert!((si.atom_number_density() * a * a * a - 8.0).abs() < 1e-12);
        let gaas = Lattice::gallium_arsenide();
        let a = gaas.lattice_constant();
        assert!((gaas.number_density_of(31) * a * a * a - 4.0).abs() < 1e-12);
        assert!((gaas.number_density_of(33) * a * a * a - 4.0).abs() < 1e-12);
        assert_eq!(gaas.number_density_of(14), 0.0);
    }

    #[test]
    fn hexagonal_volume_and_density() {
        for l in [
            Lattice::gallium_nitride(),
            Lattice::silicon_carbide_4h(),
            Lattice::silicon_carbide_6h(),
        ] {
            let (a, c) = (l.lattice_constant(), l.lattice_constant_c());
            let v = 0.5 * 3.0f64.sqrt() * a * a * c;
            assert!((l.primitive_cell_volume() / v - 1.0).abs() < 1e-14);
            let n = l.basis().len() as f64;
            assert!((l.atom_number_density() * v / n - 1.0).abs() < 1e-14);
        }
    }

    #[test]
    fn diamond_basis_bond_is_a_quarter_body_diagonal() {
        // The two sites differ by (3/4)(1,1,1) a, which is -(1/4)(1,1,1) a
        // plus the lattice vector (1,1,1) a: the nearest-neighbour bond is
        // a sqrt(3) / 4 along <111>.
        let si = Lattice::silicon();
        let p = si.basis_positions();
        let a = si.lattice_constant();
        let d = [
            p[1].1[0] - p[0].1[0],
            p[1].1[1] - p[0].1[1],
            p[1].1[2] - p[0].1[2],
        ];
        for x in d {
            assert!((x / a - 0.75).abs() < 1e-14);
        }
    }

    #[test]
    fn conventional_cell_has_eight_distinct_atoms() {
        for l in [Lattice::silicon(), Lattice::gallium_arsenide()] {
            let s = l.conventional_cell_sites();
            assert_eq!(s.len(), 8);
            for i in 0..8 {
                for j in 0..i {
                    let d: f64 = (0..3).map(|k| (s[i].1[k] - s[j].1[k]).abs()).sum();
                    assert!(d > 0.1, "{:?} and {:?} coincide", s[i], s[j]);
                }
            }
        }
        // Zincblende: Ga on the fcc sites of the cube, As shifted by 1/4.
        let s = Lattice::gallium_arsenide().conventional_cell_sites();
        assert_eq!(s[0], (31, [0.0, 0.0, 0.0]));
        assert_eq!(s[1], (31, [0.0, 0.5, 0.5]));
        assert_eq!(s[4], (33, [0.25, 0.25, 0.25]));
        assert_eq!(s[7], (33, [0.75, 0.75, 0.25]));
    }

    #[test]
    fn cubic_orthogonal_cell_is_the_conventional_cube() {
        let l = Lattice::gallium_arsenide();
        let (e, s) = l.orthogonal_cell();
        assert_eq!(e, [l.lattice_constant(); 3]);
        assert_eq!(s, l.conventional_cell_sites());
    }

    #[test]
    fn reciprocal_vectors_are_dual() {
        let l = Lattice::germanium();
        let h = Lattice::silicon_carbide_4h();
        for c in [
            l.primitive_vectors(),
            l.conventional_vectors(),
            h.conventional_vectors(),
        ] {
            let b = reciprocal_vectors(c);
            for (i, ci) in c.iter().enumerate() {
                for (j, bj) in b.iter().enumerate() {
                    let want = if i == j { 2.0 * PI } else { 0.0 };
                    assert!((dot(*ci, *bj) - want).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn miller_normals_and_directions() {
        let si = Lattice::silicon();
        let n = si.plane_normal([1, 1, 0]).unwrap();
        let s = 0.5f64.sqrt();
        assert!((n[0] - s).abs() < 1e-15 && (n[1] - s).abs() < 1e-15 && n[2].abs() < 1e-15);
        let d = si.direction([1, 1, 1]).unwrap();
        let t = (1.0f64 / 3.0).sqrt();
        assert!(d.iter().all(|x| (x - t).abs() < 1e-15));
        assert_eq!(
            si.plane_normal([0, 0, 0]),
            Err(CrystalError::ZeroIndex([0; 3]))
        );
        assert_eq!(
            si.direction([0, 0, 0]),
            Err(CrystalError::ZeroIndex([0; 3]))
        );
    }

    #[test]
    fn miller_bravais_conversion() {
        assert_eq!(plane_from_hkil([1, 1, -2, 0]), Ok([1, 1, 0]));
        assert_eq!(plane_from_hkil([0, 0, 0, 1]), Ok([0, 0, 1]));
        assert_eq!(direction_from_uvtw([1, 1, -2, 0]), Ok([3, 3, 0]));
        assert_eq!(direction_from_uvtw([2, -1, -1, 0]), Ok([3, 0, 0]));
        assert_eq!(direction_from_uvtw([0, 0, 0, 1]), Ok([0, 0, 1]));
        assert!(matches!(
            plane_from_hkil([1, 1, 1, 0]),
            Err(CrystalError::InvalidMillerBravais { .. })
        ));
        assert!(matches!(
            direction_from_uvtw([1, 0, 0, 0]),
            Err(CrystalError::InvalidMillerBravais { .. })
        ));
        assert!(matches!(
            direction_from_uvtw([i32::MAX, 0, -i32::MAX, 0]),
            Err(CrystalError::InvalidMillerBravais { .. })
        ));
        assert_eq!(
            Lattice::silicon().direction_uvtw([0, 0, 0, 1]),
            Err(CrystalError::NotHexagonal)
        );
        assert_eq!(
            Lattice::silicon().plane_normal_hkil([0, 0, 0, 1]),
            Err(CrystalError::NotHexagonal)
        );
    }

    #[test]
    fn bad_inputs_are_rejected() {
        assert_eq!(
            Lattice::diamond(14, -1.0, 300.0),
            Err(CrystalError::InvalidLatticeConstant(-1.0))
        );
        assert!(matches!(
            Lattice::diamond(14, f64::NAN, 300.0),
            Err(CrystalError::InvalidLatticeConstant(_))
        ));
        assert_eq!(
            Lattice::diamond(14, 5e-10, -1.0),
            Err(CrystalError::InvalidTemperature(-1.0))
        );
        assert_eq!(
            Lattice::zincblende(31, 0, 5e-10, 300.0),
            Err(CrystalError::UnknownAtomicNumber(0))
        );
        assert_eq!(
            Lattice::wurtzite(31, 7, 3e-10, 0.0, 0.375, None),
            Err(CrystalError::InvalidLatticeConstant(0.0))
        );
        assert_eq!(
            Lattice::wurtzite(31, 7, 3e-10, 5e-10, 0.5, None),
            Err(CrystalError::InvalidInternalParameter(0.5))
        );
        assert_eq!(
            Lattice::wurtzite(31, 7, 3e-10, 5e-10, 0.375, Some(-2.0)),
            Err(CrystalError::InvalidTemperature(-2.0))
        );
        for bad in ["A", "AAB", "ABA", "ABD", ""] {
            assert!(
                matches!(
                    Lattice::polytype(14, 6, bad, 3e-10, 1e-9, 0.375, None),
                    Err(CrystalError::InvalidStacking { .. })
                ),
                "{bad}"
            );
        }
    }
}
