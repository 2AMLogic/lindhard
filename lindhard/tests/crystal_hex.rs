//! Acceptance tests for the hexagonal lattices (`lindhard::ion::crystal`,
//! issue #179): wurtzite GaN and the 4H/6H-SiC polytypes, four-index
//! (Miller-Bravais) orientation input, number densities against the
//! `material` module, the stacking period, and the lattice search through a
//! hexagonal crystal.

use lindhard::ion::crystal::lattice::{
    equal_bond_wurtzite_u, GAN_DENSITY_G_CM3, MAX_POLYTYPE_LAYERS, SIC_4H_DENSITY_G_CM3,
    SIC_6H_DENSITY_G_CM3,
};
use lindhard::ion::crystal::search::{compare, path_metrics, site_position_in};
use lindhard::ion::crystal::{Candidate, CrystalError, Lattice, LatticeSearch, Orientation};
use lindhard::material::Material;
use lindhard::rng::{stream, ParticleRng};
use lindhard::units::g_cm3_to_kg_m3;
use rand_core::Rng;

fn assert_close(got: [f64; 3], want: [f64; 3], tol: f64) {
    for i in 0..3 {
        assert!(
            (got[i] - want[i]).abs() <= tol,
            "component {i}: got {got:?}, want {want:?}"
        );
    }
}

fn hexagonal_presets() -> [Lattice; 3] {
    [
        Lattice::gallium_nitride(),
        Lattice::silicon_carbide_4h(),
        Lattice::silicon_carbide_6h(),
    ]
}

/// The `<0001>` and `<11-20>` channel directions in the crystal frame.
///
/// With `a1 = (a/2, -sqrt(3) a/2, 0)`, `a2 = (a/2, sqrt(3) a/2, 0)`,
/// `a3 = (0, 0, c)` (Mehl et al. 2017) and `[uvtw] -> [u - t, v - t, w]`:
///
/// ```text
/// [0001]  -> [001] = c z                            -> (0, 0, 1)
/// [11-20] -> [330] = 3 (a1 + a2) = 3 a x            -> (1, 0, 0)
/// [-2110] -> [-300] = -3 a1                         -> (-1/2, sqrt(3)/2, 0)
/// [1-100] -> [1-10] = a1 - a2 = (0, -sqrt(3) a, 0)  -> (0, -1, 0)
/// ```
///
/// The plane normals `(0001)` and `(11-20)` are parallel to the directions
/// of the same indices (`g = b1 + b2` is along `x`), and `(1-100)` is along
/// `-y`. Both directions must be independent of `a`, `c` and the basis.
#[test]
fn channel_directions_to_1e_12() {
    let s3 = 3.0f64.sqrt();
    for l in hexagonal_presets() {
        assert_close(
            l.direction_uvtw([0, 0, 0, 1]).unwrap(),
            [0.0, 0.0, 1.0],
            1e-12,
        );
        assert_close(
            l.direction_uvtw([1, 1, -2, 0]).unwrap(),
            [1.0, 0.0, 0.0],
            1e-12,
        );
        assert_close(
            l.direction_uvtw([-2, 1, 1, 0]).unwrap(),
            [-0.5, 0.5 * s3, 0.0],
            1e-12,
        );
        assert_close(
            l.direction_uvtw([1, -1, 0, 0]).unwrap(),
            [0.0, -1.0, 0.0],
            1e-12,
        );
        assert_close(
            l.plane_normal_hkil([0, 0, 0, 1]).unwrap(),
            [0.0, 0.0, 1.0],
            1e-12,
        );
        assert_close(
            l.plane_normal_hkil([1, 1, -2, 0]).unwrap(),
            [1.0, 0.0, 0.0],
            1e-12,
        );
        assert_close(
            l.plane_normal_hkil([1, -1, 0, 0]).unwrap(),
            [0.0, -1.0, 0.0],
            1e-12,
        );
    }
}

/// A c-plane wafer, `(0001)` with reference `[11-20]`, sends a zero-tilt beam
/// down `[0001]`, whatever the twist and wafer rotation; an a-plane wafer,
/// `(11-20)` with reference `[0001]`, down `[11-20]`.
///
/// Tilted case, derived from the orientation formula
/// `d = cos θ n + sin θ cos(φ - ω) r + sin θ sin(φ - ω) t`: on the c-plane,
/// `n = (0, 0, 1)`, `r = (1, 0, 0)` and `t = n x r = (0, 1, 0)`, so tilt 7°
/// and twist 30° (`ω = 0`) give
/// `d = (sin 7° cos 30°, sin 7° sin 30°, cos 7°)`.
#[test]
fn miller_bravais_wafer_cuts() {
    let gan = Lattice::gallium_nitride();
    for (twist, rot) in [(0.0, 0.0), (0.4, -1.1), (2.0, 3.0)] {
        let c = Orientation::new_miller_bravais(&gan, [0, 0, 0, 1], [1, 1, -2, 0], 0.0, twist, rot)
            .unwrap();
        assert_close(c.beam_crystal(), [0.0, 0.0, 1.0], 1e-12);
        assert_eq!(c.normal_hkl(), [0, 0, 1]);
        assert_eq!(c.reference_uvw(), [3, 3, 0]);
        let a = Orientation::new_miller_bravais(&gan, [1, 1, -2, 0], [0, 0, 0, 1], 0.0, twist, rot)
            .unwrap();
        assert_close(a.beam_crystal(), [1.0, 0.0, 0.0], 1e-12);
    }
    let (th, ph) = (7f64.to_radians(), 30f64.to_radians());
    let o = Orientation::new_miller_bravais(
        &Lattice::silicon_carbide_4h(),
        [0, 0, 0, 1],
        [1, 1, -2, 0],
        th,
        ph,
        0.0,
    )
    .unwrap();
    assert_close(
        o.beam_crystal(),
        [th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()],
        1e-12,
    );
    // A reference outside the surface plane: [1-100] . (1-100) zone sum 2.
    assert!(matches!(
        Orientation::new_miller_bravais(&gan, [1, -1, 0, 0], [1, -1, 0, 0], 0.0, 0.0, 0.0),
        Err(CrystalError::ReferenceNotInPlane { .. })
    ));
    assert_eq!(
        Orientation::new_miller_bravais(
            &Lattice::silicon(),
            [0, 0, 0, 1],
            [1, 1, -2, 0],
            0.0,
            0.0,
            0.0
        ),
        Err(CrystalError::NotHexagonal)
    );
    assert!(matches!(
        Orientation::new_miller_bravais(&gan, [1, 1, 1, 0], [0, 0, 0, 1], 0.0, 0.0, 0.0),
        Err(CrystalError::InvalidMillerBravais { .. })
    ));
}

/// Relative difference between the lattice atom number density and the
/// `material` module's `ρ N_A / M̄` at the cited density `rho_g_cm3`, for the
/// whole crystal and per species.
fn density_differences(l: &Lattice, rho_g_cm3: f64, z: [u8; 2]) -> [f64; 3] {
    let m =
        Material::from_atom_fractions(&[(z[0], 1.0), (z[1], 1.0)], Some(g_cm3_to_kg_m3(rho_g_cm3)))
            .unwrap();
    [
        l.atom_number_density() / m.atom_number_density() - 1.0,
        l.number_density_of(z[0]) / m.number_density_of(z[0]).unwrap() - 1.0,
        l.number_density_of(z[1]) / m.number_density_of(z[1]).unwrap() - 1.0,
    ]
}

/// Atom number density against the `material` module: the cases that agree.
///
/// * 6H-SiC: the Ioffe archive density 3.21 g/cm³ (300 K) is printed to
///   three figures, half a unit in the last place is 1.56e-3 relative; the
///   cell parameters (`a = 3.08129 Å`, `c = 15.11976 Å`) add under 2e-6 and
///   the Si and C standard-weight intervals about 5e-5. Tolerance 1.7e-3.
///   Found: +1.05e-3.
///
/// Issue #179 asks for agreement within the cited precision. Only 6H meets
/// it with the sources that could be opened; 4H and GaN are in
/// [`number_density_discrepancies_are_recorded_gaps`], which does not claim it.
#[test]
fn number_density_agrees_with_the_material_module_6h() {
    for rel in density_differences(
        &Lattice::silicon_carbide_6h(),
        SIC_6H_DENSITY_G_CM3,
        [14, 6],
    ) {
        assert!(rel.abs() < 1.7e-3, "6H relative difference {rel:e}");
    }
}

/// **Not an agreement test.** Regression pins on two known, unresolved
/// disagreements between sources, so that a change to a cell or a density is
/// noticed. The acceptance criterion of #179 (density within the cited
/// precision) is **not met** for these presets, and the test must not be read
/// as meeting it. Neither constant was adjusted to hide the offset; both are
/// listed under Open questions in `docs/data-provenance.md`.
///
/// * 4H-SiC: against the archive's 3.211 g/cm³ (half a unit: 1.6e-4) the
///   lattice is +7.52e-4 denser, outside the printed precision. The two
///   numbers are different measurements (Bauer et al. 2001 cell, Gomes de
///   Mesquita 1967 density, per the archive) and the cell's temperature is
///   not stated.
/// * GaN: the archive's 6.15 g/cm³ is 0.85 % above the X-ray density of the
///   archive's own `a` and `c` (lattice/material - 1 = -8.52e-3), far outside
///   its printed precision (8e-4).
#[test]
fn number_density_discrepancies_are_recorded_gaps() {
    for rel in density_differences(
        &Lattice::silicon_carbide_4h(),
        SIC_4H_DENSITY_G_CM3,
        [14, 6],
    ) {
        assert!((rel - 7.52e-4).abs() < 5e-6, "4H pinned difference {rel:e}");
    }
    for rel in density_differences(&Lattice::gallium_nitride(), GAN_DENSITY_G_CM3, [31, 7]) {
        assert!(
            (rel + 8.52e-3).abs() < 5e-6,
            "GaN pinned difference {rel:e}"
        );
    }
}

/// The close-packed column (A, B or C) of lattice coordinates `(f1, f2)`.
fn column_letter(f: [f64; 3]) -> char {
    let near = |x: f64, y: f64| {
        let d = |p: f64, q: f64| {
            let r = (p - q).rem_euclid(1.0);
            r.min(1.0 - r)
        };
        d(f[0], x) < 1e-9 && d(f[1], y) < 1e-9
    };
    if near(0.0, 0.0) {
        'A'
    } else if near(1.0 / 3.0, 2.0 / 3.0) {
        'B'
    } else if near(2.0 / 3.0, 1.0 / 3.0) {
        'C'
    } else {
        panic!("{f:?} is not on a close-packed column")
    }
}

/// Count the layers of species `z` in one `c` period: the distinct heights
/// (to `1e-3 c`) of its sites, in order, with their column letters.
fn layers(l: &Lattice, z: u8) -> Vec<(f64, char)> {
    let mut v: Vec<(f64, char)> = l
        .basis()
        .iter()
        .filter(|s| s.z == z)
        .map(|s| (s.fractional[2], column_letter(s.fractional)))
        .collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    v.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-3);
    v
}

/// Relabel a stacking so that it starts `A`, `B` (the letters are only
/// defined up to a permutation and a cyclic shift of the origin).
fn canonical(seq: &[char]) -> String {
    let (x, y) = (seq[0], seq[1]);
    let w = ['A', 'B', 'C']
        .into_iter()
        .find(|&c| c != x && c != y)
        .unwrap();
    seq.iter()
        .map(|&c| {
            if c == x {
                'A'
            } else if c == y {
                'B'
            } else {
                assert_eq!(c, w);
                'C'
            }
        })
        .collect()
}

/// The stacking period, by counting layers per `c`: 2 for wurtzite, 4 for
/// 4H and 6 for 6H, for both species, evenly spaced by `c/n` (to the
/// `1e-3 c` the measured sites deviate from the ideal), with each anion
/// layer `3/4` of a layer spacing above its cation layer (to `2e-3 c`), and
/// the sequences ABAB, ABAC (ABCB relabelled) and ABCACB as printed by Mehl
/// et al. (2017), pp. 423-427.
#[test]
fn stacking_period_from_layers_per_c() {
    let cases = [
        (Lattice::gallium_nitride(), [31u8, 7], 2usize, "AB"),
        (Lattice::silicon_carbide_4h(), [14, 6], 4, "ABAC"),
        (Lattice::silicon_carbide_6h(), [14, 6], 6, "ABCACB"),
    ];
    for (l, [cation, anion], n, seq) in cases {
        let c_layers = layers(&l, cation);
        let a_layers = layers(&l, anion);
        assert_eq!(c_layers.len(), n, "{:?}", l.structure());
        assert_eq!(a_layers.len(), n, "{:?}", l.structure());
        let nf = n as f64;
        let z0 = c_layers[0].0;
        for (k, &(z, _)) in c_layers.iter().enumerate() {
            assert!((z - z0 - k as f64 / nf).abs() < 1e-3, "{k}: {z}");
        }
        // Every anion layer sits on its cation's column, 3/4 of a spacing up.
        for &(zc, col) in &c_layers {
            let above = a_layers
                .iter()
                .find(|&&(za, _)| (za - zc).rem_euclid(1.0) < 1.0 / nf)
                .unwrap();
            assert_eq!(above.1, col);
            let dz = (above.0 - zc).rem_euclid(1.0);
            assert!((dz * nf - 0.75).abs() < 2e-3 * nf, "bond {dz}");
        }
        let letters: Vec<char> = c_layers.iter().map(|x| x.1).collect();
        assert_eq!(canonical(&letters), seq);
    }
    // The ideal constructor gives the same counts.
    for (stacking, n) in [("AB", 2), ("ABCB", 4), ("ABCACB", 6), ("ABCBCACAB", 9)] {
        let l = Lattice::polytype(14, 6, stacking, 3e-10, 2.5e-10 * n as f64, 0.375, None).unwrap();
        assert_eq!(layers(&l, 14).len(), n);
        assert_eq!(layers(&l, 6).len(), n);
        assert_eq!(l.basis().len(), 2 * n);
    }
}

fn same_sites(x: &Lattice, y: &Lattice, tol_z: f64) {
    assert_eq!(x.basis().len(), y.basis().len());
    for s in x.basis() {
        let hit = y.basis().iter().any(|t| {
            let d = |i: usize| {
                let r = (s.fractional[i] - t.fractional[i]).rem_euclid(1.0);
                r.min(1.0 - r)
            };
            t.z == s.z && d(0) < 1e-12 && d(1) < 1e-12 && d(2) < tol_z
        });
        assert!(hit, "{s:?} has no partner");
    }
}

/// The measured 4H and 6H sites lie within `1e-3 c` of the ideal polytype
/// with the same stacking and `u = 3/8`, and the wurtzite constructor is the
/// two-layer polytype `"BC"` exactly.
#[test]
fn measured_polytypes_are_near_ideal() {
    let h4 = Lattice::silicon_carbide_4h();
    let ideal4 = Lattice::polytype(
        14,
        6,
        "ABAC",
        h4.lattice_constant(),
        h4.lattice_constant_c(),
        0.375,
        None,
    )
    .unwrap();
    same_sites(&h4, &ideal4, 1e-3);
    let h6 = Lattice::silicon_carbide_6h();
    let ideal6 = Lattice::polytype(
        14,
        6,
        "ABCACB",
        h6.lattice_constant(),
        h6.lattice_constant_c(),
        0.375,
        None,
    )
    .unwrap();
    same_sites(&h6, &ideal6, 1e-3);
    let w = Lattice::wurtzite(31, 7, 3.2e-10, 5.2e-10, 0.377, Some(300.0)).unwrap();
    let p = Lattice::polytype(31, 7, "BC", 3.2e-10, 5.2e-10, 0.377, Some(300.0)).unwrap();
    same_sites(&w, &p, 1e-15);
}

/// The placeholder GaN `u` is the equal-bond value: every Ga has four N
/// neighbours at one distance (to 1e-12 relative), the shortest Ga-N
/// distances in the crystal. At the ideal `c/a = sqrt(8/3)` it is 3/8.
#[test]
fn gan_placeholder_u_gives_four_equal_bonds() {
    assert!((equal_bond_wurtzite_u(1.0, (8.0f64 / 3.0).sqrt()) - 0.375).abs() < 1e-15);
    let gan = Lattice::gallium_nitride();
    let p = gan.primitive_vectors();
    let pos = gan.basis_positions();
    let ga = pos[0].1;
    let mut d: Vec<f64> = Vec::new();
    for i in -2i32..=2 {
        for j in -2i32..=2 {
            for k in -2i32..=2 {
                for &(z, r) in &pos {
                    if z != 7 {
                        continue;
                    }
                    let t = [0, 1, 2].map(|m| {
                        r[m] + f64::from(i) * p[0][m]
                            + f64::from(j) * p[1][m]
                            + f64::from(k) * p[2][m]
                    });
                    d.push(((0..3).map(|m| (t[m] - ga[m]).powi(2)).sum::<f64>()).sqrt());
                }
            }
        }
    }
    d.sort_by(f64::total_cmp);
    for x in &d[1..4] {
        assert!((x / d[0] - 1.0).abs() < 1e-12, "{:?}", &d[..5]);
    }
    assert!(d[4] / d[0] > 1.5, "{:?}", &d[..5]);
}

fn uniform(rng: &mut ParticleRng) -> f64 {
    (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64
}

/// The orthohexagonal cell tiles the crystal: in a block of 4 x 4 x 4 cells
/// it holds exactly the sites of the hexagonal lattice `r + i a1 + j a2 +
/// k c` that fall in the block, each once (to 1e-12 a).
#[test]
fn orthogonal_cell_tiles_the_hexagonal_crystal() {
    for l in hexagonal_presets() {
        let (e, sites) = l.orthogonal_cell();
        assert_eq!(sites.len(), 2 * l.basis().len());
        let mut from_box: Vec<(u8, [f64; 3])> = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    for &(z, f) in &sites {
                        from_box.push((z, site_position_in(e, [i, j, k], f)));
                    }
                }
            }
        }
        let p = l.primitive_vectors();
        let tol = 1e-12 * l.lattice_constant();
        let inside = |r: [f64; 3]| (0..3).all(|m| r[m] > -tol && r[m] < 4.0 * e[m] - tol);
        let mut from_hex = 0usize;
        for i in -8i32..=12 {
            for j in -8i32..=12 {
                for k in -1i32..=5 {
                    for &(z, r) in &l.basis_positions() {
                        let t = [0, 1, 2].map(|m| {
                            r[m] + f64::from(i) * p[0][m]
                                + f64::from(j) * p[1][m]
                                + f64::from(k) * p[2][m]
                        });
                        if !inside(t) {
                            continue;
                        }
                        from_hex += 1;
                        let hits = from_box
                            .iter()
                            .filter(|&&(zb, b)| {
                                zb == z && (0..3).all(|m| (b[m] - t[m]).abs() < tol)
                            })
                            .count();
                        assert_eq!(hits, 1, "{t:?}");
                    }
                }
            }
        }
        assert_eq!(from_hex, from_box.len());
    }
}

/// O(N) reference over a block of orthogonal cells, through the same
/// `site_position_in` and `path_metrics` as the search.
fn brute(
    l: &Lattice,
    hi: [i64; 3],
    origin: [f64; 3],
    d: [f64; 3],
    p_max: f64,
    length: f64,
) -> Vec<Candidate> {
    let (e, sites) = l.orthogonal_cell();
    let mut out = Vec::new();
    for i in 0..hi[0] {
        for j in 0..hi[1] {
            for k in 0..hi[2] {
                for (bi, &(z, f)) in sites.iter().enumerate() {
                    let cell = [i, j, k];
                    let pos = site_position_in(e, cell, f);
                    let (s, p2) = path_metrics(origin, d, pos);
                    if s >= 0.0 && s <= length && p2 <= p_max * p_max {
                        out.push(Candidate {
                            z,
                            cell,
                            basis: bi as u8,
                            position: pos,
                            s,
                            p: p2.sqrt(),
                        });
                    }
                }
            }
        }
    }
    out.sort_by(compare);
    out
}

/// The cell-walk search agrees bit for bit with the brute-force scan on
/// random segments through GaN, 4H- and 6H-SiC, and finds both species.
#[test]
fn search_matches_brute_force_hexagonal() {
    for (seed, l) in hexagonal_presets().into_iter().enumerate() {
        let search = LatticeSearch::new(&l);
        let e = search.cell_edges();
        assert_eq!(search.lattice_constant(), l.lattice_constant());
        // About 8 a along every axis.
        let n_cells = [8i64, 5, (8.0 * e[0] / e[2]).ceil().max(3.0) as i64];
        let mut total = 0;
        for idx in 0..2_000u64 {
            let mut rng = stream(179 + seed as u64, idx);
            let origin = [0, 1, 2].map(|k| e[k] * n_cells[k] as f64 * uniform(&mut rng));
            let dir = [0, 1, 2].map(|_| 2.0 * uniform(&mut rng) - 1.0);
            let norm = dir.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-3 {
                continue;
            }
            let d = dir.map(|x| x / norm);
            let a = e[0];
            let p_max = a * (0.02 + 1.5 * uniform(&mut rng));
            let length = a * (0.01 + 6.0 * uniform(&mut rng).powi(2));
            let want = brute(&l, n_cells, origin, d, p_max, length);
            let mut got = search.search(origin, dir, p_max, length).unwrap();
            got.retain(|c| (0..3).all(|k| c.cell[k] >= 0 && c.cell[k] < n_cells[k]));
            assert_eq!(got, want, "ray {idx}");
            total += want.len();
        }
        assert!(total > 2_000, "{total}");
        let a = e[0];
        let v = search
            .search([0.3 * a, 0.2 * a, 0.1 * a], [1.0, 2.0, 3.0], a, 5.0 * a)
            .unwrap();
        let zs: Vec<u8> = l.basis().iter().map(|s| s.z).collect();
        for z in zs {
            assert!(v.iter().any(|c| c.z == z), "{z}");
        }
    }
}

/// Down the `[0001]` axis through a Ga column, the search meets one Ga and
/// one N per half period: 2 bilayers per `c`, at heights `0, u c, c/2,
/// (1/2 + u) c` in the column of `(1/3, 2/3)` and `(2/3, 1/3)` alternately.
/// Along an empty column (A, `(0, 0)`) it meets nothing at a small radius.
#[test]
fn search_along_c_counts_the_layers() {
    let gan = Lattice::gallium_nitride();
    let search = LatticeSearch::new(&gan);
    let (a, c) = (gan.lattice_constant(), gan.lattice_constant_c());
    let col_b = gan.basis_positions()[0].1;
    let v = search
        .search(
            [col_b[0], col_b[1], -0.01 * c],
            [0.0, 0.0, 1.0],
            0.05 * a,
            3.0 * c,
        )
        .unwrap();
    // The B column holds Ga at k c and N at (k + u) c; C-column atoms are
    // a / sqrt(3) away. Three periods: 3 Ga and 3 N.
    let ga = v.iter().filter(|x| x.z == 31).count();
    let n = v.iter().filter(|x| x.z == 7).count();
    assert_eq!((ga, n), (3, 3), "{v:?}");
    let empty = search
        .search([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.05 * a, 3.0 * c)
        .unwrap();
    assert!(empty.is_empty(), "{empty:?}");
}

/// `Candidate::basis` is a `u8`, so a polytype must not have more than
/// `MAX_POLYTYPE_LAYERS` (64) letters: its orthohexagonal cell has `4n`
/// sites. The longest accepted stacking keeps every site index distinct in a
/// search; one letter pair more is rejected rather than wrapped.
#[test]
fn long_polytype_basis_indices_do_not_alias() {
    let n = MAX_POLYTYPE_LAYERS;
    let make = |stacking: &str| Lattice::polytype(14, 6, stacking, 3.0e-10, 1.0e-9, 0.375, None);
    assert!(matches!(
        make(&"AB".repeat(n / 2 + 1)),
        Err(CrystalError::InvalidStacking { .. })
    ));
    let l = make(&"AB".repeat(n / 2)).unwrap();
    assert_eq!(l.orthogonal_cell().1.len(), 4 * n);
    let (a, c) = (l.lattice_constant(), l.lattice_constant_c());
    let v = LatticeSearch::new(&l)
        .search([0.0, 0.0, -0.01 * c], [0.0, 0.0, 1.0], 3.0 * a, 1.02 * c)
        .unwrap();
    assert!(v.iter().any(|x| x.basis >= 128), "high indices reached");
    let mut keys: Vec<_> = v.iter().map(|x| (x.cell, x.basis)).collect();
    let total = keys.len();
    keys.sort_unstable();
    keys.dedup();
    assert_eq!(keys.len(), total, "a (cell, basis) pair repeats");
}
