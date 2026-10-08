//! Cubic lattices: Bravais vectors, a basis of sites with species, and a
//! lattice constant stated with its temperature.
//!
//! # Structures
//!
//! Both structures are a face-centred cubic Bravais lattice with a two-site
//! basis. The primitive vectors and the basis are taken as printed in
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
//! # Miller indices
//!
//! For cubic crystals, Miller indices `(hkl)` and direction indices `[uvw]`
//! refer to the **conventional** cubic cell (edge `a`, axes along the crystal
//! frame axes), not to the primitive vectors above. The reciprocal vectors of
//! a cell `(c1, c2, c3)` are `b_i = 2 pi (c_j x c_k) / V` (Mehl et al.,
//! eqs. (9)-(10)). The planes `(hkl)` are orthogonal to the reciprocal
//! lattice vector `g = h b1 + k b2 + l b3` (the definition of Miller indices
//! via the reciprocal lattice, as given in the Wikipedia article "Miller
//! index", revision 1378193236 of 2026-10-03, which cites N. W. Ashcroft and
//! N. D. Mermin, *Solid State Physics* (1976); the book itself was not
//! opened). A direction `[uvw]` is `u c1 + v c2 + w c3`, and by
//! `c_i . b_j = 2 pi delta_ij` (eq. (9)) it lies in the plane `(hkl)` exactly
//! when `h u + k v + l w = 0` (the Weiss zone law, derived here from eq. (9)).
//! For the cubic cell `g` is parallel to `[hkl]`.
//!
//! # Lattice constants
//!
//! A lattice constant is a measured fact, so every preset cites the record it
//! was read from and states the temperature it refers to; see
//! [`LatticeConstant`] and the `docs/data-provenance.md` row. The lattice is
//! rigid: no thermal expansion is applied, and the temperature is carried as
//! metadata of the cited value.

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

/// Which cubic structure a [`Lattice`] has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Structure {
    /// Diamond (A4, Fd-3m): one species on both fcc sublattices.
    Diamond,
    /// Zincblende (B3, F-43m): two species, one per fcc sublattice.
    Zincblende,
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

/// A cubic crystal: a face-centred cubic Bravais lattice with a basis.
#[derive(Debug, Clone, PartialEq)]
pub struct Lattice {
    structure: Structure,
    a_m: f64,
    temperature_k: f64,
    primitive: [[f64; 3]; 3],
    basis: Vec<Site>,
}

fn check_z(z: u8) -> Result<u8, CrystalError> {
    element(z)
        .map(|_| z)
        .ok_or(CrystalError::UnknownAtomicNumber(z))
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

fn combine(c: [[f64; 3]; 3], n: [f64; 3]) -> [f64; 3] {
    add(add(scale(c[0], n[0]), scale(c[1], n[1])), scale(c[2], n[2]))
}

fn as_f64(n: [i32; 3]) -> [f64; 3] {
    [f64::from(n[0]), f64::from(n[1]), f64::from(n[2])]
}

impl Lattice {
    fn build(
        structure: Structure,
        sites: [(u8, [f64; 3]); 2],
        a_m: f64,
        temperature_k: f64,
    ) -> Result<Self, CrystalError> {
        if !(a_m.is_finite() && a_m > 0.0) {
            return Err(CrystalError::InvalidLatticeConstant(a_m));
        }
        if !(temperature_k.is_finite() && temperature_k >= 0.0) {
            return Err(CrystalError::InvalidTemperature(temperature_k));
        }
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
            temperature_k,
            primitive,
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

    /// The structure type.
    pub fn structure(&self) -> Structure {
        self.structure
    }

    /// Cubic lattice constant (conventional cell edge), metres.
    pub fn lattice_constant(&self) -> f64 {
        self.a_m
    }

    /// Temperature the lattice constant refers to, kelvin.
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// The face-centred cubic primitive (Bravais) vectors `a1, a2, a3`,
    /// metres, in the crystal frame.
    pub fn primitive_vectors(&self) -> [[f64; 3]; 3] {
        self.primitive
    }

    /// The conventional cubic cell vectors `a x, a y, a z`, metres. Miller
    /// and direction indices refer to this cell.
    pub fn conventional_vectors(&self) -> [[f64; 3]; 3] {
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

    /// The atoms of one conventional cubic cell, as fractional coordinates of
    /// that cell wrapped into `[0, 1)`: each basis site plus the three
    /// face-centring translations `(0, 1/2, 1/2)`, `(1/2, 0, 1/2)`,
    /// `(1/2, 1/2, 0)`. Eight atoms for both structures, in basis order then
    /// translation order.
    pub fn conventional_cell_sites(&self) -> Vec<(u8, [f64; 3])> {
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

    /// Primitive cell volume `a1 . (a2 x a3)`, m³ (Mehl et al. 2017, eq. (6)).
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_carry_their_cited_constants() {
        let si = Lattice::silicon();
        assert_eq!(si.lattice_constant(), 5.431_020_511e-10);
        assert_eq!(si.temperature_k(), 295.65);
        assert_eq!(si.structure(), Structure::Diamond);
        let gaas = Lattice::gallium_arsenide();
        assert_eq!(gaas.structure(), Structure::Zincblende);
        assert_eq!(gaas.basis()[0].z, 31);
        assert_eq!(gaas.basis()[1].z, 33);
        assert_eq!(Lattice::germanium().basis()[0].z, 32);
        let sic = Lattice::silicon_carbide_3c();
        assert_eq!((sic.basis()[0].z, sic.basis()[1].z), (14, 6));
        assert_eq!(sic.temperature_k(), 297.0);
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
    fn reciprocal_vectors_are_dual() {
        let l = Lattice::germanium();
        for c in [l.primitive_vectors(), l.conventional_vectors()] {
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
    }
}
