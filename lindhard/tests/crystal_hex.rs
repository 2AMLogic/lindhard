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
use lindhard::rng::{stream, uniform};
use lindhard::units::g_cm3_to_kg_m3;

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

/// One row of the density comparison: what each source printed, so that a
/// failure message (or `--nocapture`) states the evidence.
struct DensityCase {
    name: &'static str,
    lattice: Lattice,
    rho_g_cm3: f64,
    /// Elements of the compound.
    z: [u8; 2],
    /// Printed density and its temperature, as read from the source.
    density_source: &'static str,
    /// Cell source and its temperature (or that none is stated).
    cell_source: &'static str,
    /// Half a unit in the last printed place of the density, relative.
    printed_precision: f64,
}

fn density_cases() -> [DensityCase; 3] {
    [
        DensityCase {
            name: "4H-SiC",
            lattice: Lattice::silicon_carbide_4h(),
            rho_g_cm3: SIC_4H_DENSITY_G_CM3,
            z: [14, 6],
            density_source: "Ioffe NSM 3.211 g/cm3, 300 K, credited to Gomes de Mesquita 1967",
            cell_source: "Bauer 2001 via Mehl 2017: a=3.08051, c=10.0848 A, temperature not stated",
            printed_precision: 0.0005 / 3.211,
        },
        DensityCase {
            name: "6H-SiC",
            lattice: Lattice::silicon_carbide_6h(),
            rho_g_cm3: SIC_6H_DENSITY_G_CM3,
            z: [14, 6],
            density_source: "Ioffe NSM 3.21 g/cm3, 300 K, credited to Harris 1995b",
            cell_source:
                "Bauer 2001 via Mehl 2017: a=3.08129, c=15.11976 A, temperature not stated",
            printed_precision: 0.005 / 3.21,
        },
        DensityCase {
            name: "GaN",
            lattice: Lattice::gallium_nitride(),
            rho_g_cm3: GAN_DENSITY_G_CM3,
            z: [31, 7],
            density_source: "Juza & Hahn 1938 via COD 1010168: measured 6.1 g/cm3, \
                             temperature not stated",
            cell_source: "Melnik 1997 Table 1, annealed bulk: a=3.1890 A (+-0.00005), \
                          c=5.18500 A (+-0.00002), 300 K",
            printed_precision: 0.05 / 6.1,
        },
    ]
}

/// Prints, for each material, the source values, temperatures, printed
/// precision, the computed densities with units, and the relative difference
/// of the lattice against the `material` module (whole crystal and per
/// species). The computed densities are the lattice and `material` atom number
/// densities in atoms/m³ (whole crystal, then per species), and the cited mass
/// density next to the mass density the lattice implies, ρ_cited (1 + Δ), in
/// g/cm³; both sides share one mean atomic weight, so the ratio of number
/// densities is the ratio of mass densities. Returns the whole-crystal
/// difference.
fn report(c: &DensityCase) -> f64 {
    let d = density_differences(&c.lattice, c.rho_g_cm3, c.z);
    let m = Material::from_atom_fractions(
        &[(c.z[0], 1.0), (c.z[1], 1.0)],
        Some(g_cm3_to_kg_m3(c.rho_g_cm3)),
    )
    .unwrap();
    println!(
        "{}: density [{}]; cell [{}]; printed precision (half last digit) {:.2e}; \
         atom density lattice {:.6e} /m3 (Z={} {:.6e}, Z={} {:.6e}), \
         material {:.6e} /m3 (Z={} {:.6e}, Z={} {:.6e}); \
         mass density cited {:.5} g/cm3, implied by lattice {:.5} g/cm3; \
         lattice/material - 1 = {:+.3e} (species {:+.3e}, {:+.3e})",
        c.name,
        c.density_source,
        c.cell_source,
        c.printed_precision,
        c.lattice.atom_number_density(),
        c.z[0],
        c.lattice.number_density_of(c.z[0]),
        c.z[1],
        c.lattice.number_density_of(c.z[1]),
        m.atom_number_density(),
        c.z[0],
        m.number_density_of(c.z[0]).unwrap(),
        c.z[1],
        m.number_density_of(c.z[1]).unwrap(),
        c.rho_g_cm3,
        c.rho_g_cm3 * (1.0 + d[0]),
        d[0],
        d[1],
        d[2]
    );
    d[0]
}

/// Every material states its evidence; the species differences equal the
/// whole-crystal one because both sides share the stoichiometry.
#[test]
fn number_density_report_states_the_sources() {
    for c in density_cases() {
        let d = density_differences(&c.lattice, c.rho_g_cm3, c.z);
        // Same stoichiometry on both sides, so species and total agree.
        assert!(
            (d[0] - d[1]).abs() < 1e-12 && (d[0] - d[2]).abs() < 1e-12,
            "{}",
            c.name
        );
        report(&c);
    }
}

/// Atom number density against the `material` module: the cases that agree.
///
/// * 6H-SiC: the Ioffe archive density 3.21 g/cm³ (300 K) is printed to
///   three figures, half a unit in the last place is 1.56e-3 relative; the
///   cell parameters (`a = 3.08129 Å`, `c = 15.11976 Å`) add under 2e-6 and
///   the Si and C standard-weight intervals about 5e-5. Tolerance 1.7e-3.
///   Found: +1.05e-3.
/// * GaN: the measured density 6.1 g/cm³ of Juza and Hahn (1938, via COD
///   1010168; temperature not stated) is printed to two figures, half a unit
///   in the last place is 8.20e-3 relative. The cell of Melnik et al. (1997)
///   at 300 K adds 3.5e-5 (its stated accuracy, 0.00005 Å on `a` twice and
///   0.00002 Å on `c`) and the Ga and N standard-weight intervals under 1e-4.
///   Tolerance 8.3e-3. Found: −1.74e-3. The density is the only open
///   measured GaN bulk density found, and its two figures make this a coarse
///   check; see `docs/data-provenance.md`.
///
/// Issue #179 asks for agreement within the cited precision. 6H and GaN meet
/// it with the sources that could be opened; 4H meets it against the Ioffe
/// archive value only in the sense of
/// [`number_density_4h_agrees_within_measurement_spread`] (provisional, spread
/// as uncertainty), and the printed-precision miss is pinned in
/// [`number_density_discrepancies_are_recorded_gaps`].
#[test]
fn number_density_agrees_with_the_material_module() {
    for rel in density_differences(
        &Lattice::silicon_carbide_6h(),
        SIC_6H_DENSITY_G_CM3,
        [14, 6],
    ) {
        assert!(rel.abs() < 1.7e-3, "6H relative difference {rel:e}");
    }
    for rel in density_differences(&Lattice::gallium_nitride(), GAN_DENSITY_G_CM3, [31, 7]) {
        assert!(rel.abs() < 8.3e-3, "GaN relative difference {rel:e}");
    }
}

/// **Not an agreement test.** Regression pin on a known, unresolved
/// disagreement between sources, so that a change to the cell or the density
/// is noticed. Against the Ioffe value the printed-precision criterion is
/// **not met** for 4H-SiC, and this test must not be read as meeting it; the
/// ruled, spread-based agreement is
/// [`number_density_4h_agrees_within_measurement_spread`]. Neither constant was adjusted to hide the offset; the gap is
/// listed under Open questions in `docs/data-provenance.md`.
///
/// * 4H-SiC: against the archive's 3.211 g/cm³ (half a unit: 1.6e-4) the
///   lattice is +7.52e-4 denser, outside the printed precision. The two
///   numbers are different measurements (Bauer et al. 2001 cell, Gomes de
///   Mesquita 1967 density, per the archive) and the cell's temperature is
///   not stated.
///
/// The open, primary measured 4H densities found on 2026-10-11 (third pass,
/// see Open questions in `docs/data-provenance.md`) do not close the gap
/// either; they are pinned in
/// [`opened_measured_sic_densities_are_recorded_gaps`].
#[test]
fn number_density_discrepancies_are_recorded_gaps() {
    for rel in density_differences(
        &Lattice::silicon_carbide_4h(),
        SIC_4H_DENSITY_G_CM3,
        [14, 6],
    ) {
        assert!((rel - 7.52e-4).abs() < 5e-6, "4H pinned difference {rel:e}");
    }
}

/// The three opened primary measured 4H-SiC bulk densities, g/cm³: Lundqvist
/// 1948 (3.201 ± 0.005), ORNL/TM-2023/3059 wafer A (3.17) and wafer B (3.21).
/// See [`opened_measured_sic_densities_are_recorded_gaps`] for the sources.
const MEASURED_4H_DENSITIES_G_CM3: [f64; 3] = [3.201, 3.17, 3.21];

/// 4H-SiC density agreement, with the measurement spread as the stated
/// uncertainty (operator ruling on #179, 2026-10-28, option 3). **Provisional**,
/// not a validation: the three measurements differ in material quality and
/// method, and their spread is taken as the uncertainty by ruling, not shown to
/// be one.
///
/// Derivation, fixed before the cell was compared:
///
/// * reference density ρ̄ = mean of the three measured values
///   = 9.581 / 3 = 3.19367 g/cm³;
/// * uncertainty = half the range of the three values
///   = (3.21 − 3.17) / 2 = 0.020 g/cm³;
/// * relative tolerance = 0.020 / ρ̄ = 6.26e-3.
///
/// The Bauer 4H cell implies 3.21342 g/cm³, +6.18e-3 from ρ̄, inside the
/// tolerance by 7.8e-5 (1.2 %). The margin is thin and depends on the choice of
/// reference: against the midrange 3.19 the cell would sit at +7.3e-3, outside
/// the same half-range. The reference (mean) and the half-range were chosen as
/// the plain central value and spread of three numbers; the tolerance is
/// computed here from the three values, not set from the result. The 6H-SiC and
/// GaN tolerances and all lattice constants are unchanged.
#[test]
fn number_density_4h_agrees_within_measurement_spread() {
    let v = MEASURED_4H_DENSITIES_G_CM3;
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let min = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let half_range = (max - min) / 2.0;
    let tolerance = half_range / mean;
    let lattice = Lattice::silicon_carbide_4h();
    let d = density_differences(&lattice, mean, [14, 6]);
    println!(
        "4H-SiC: measured densities {v:?} g/cm3 (Lundqvist 1948 3.201 +-0.005; ORNL/TM-2023/3059 \
         wafers A 3.17, B 3.21; temperatures not stated); mean {mean:.5}, half-range \
         {half_range:.4} g/cm3, relative tolerance {tolerance:.3e}; cell Bauer 2001 via Mehl 2017 \
         implies {:.5} g/cm3; lattice/material - 1 = {:+.3e} (species {:+.3e}, {:+.3e}); provisional",
        mean * (1.0 + d[0]),
        d[0],
        d[1],
        d[2]
    );
    assert!(
        (d[0] - d[1]).abs() < 1e-12 && (d[0] - d[2]).abs() < 1e-12,
        "species differences equal the total"
    );
    assert!(
        (tolerance - 6.26e-3).abs() < 1e-5,
        "tolerance {tolerance:e}"
    );
    for rel in d {
        assert!(
            rel.abs() < tolerance,
            "4H relative difference {rel:e} vs {tolerance:e}"
        );
    }
}

/// **Not an agreement test.** The Bauer et al. (2001) cells against measured
/// SiC bulk densities from sources that were opened, pinned so that the
/// evidence behind the 4H gap is reproducible. None of these values is used by
/// a preset, and no tolerance or constant was changed because of them.
///
/// * D. Lundqvist, Acta Chem. Scand. 2, 177 (1948),
///   doi:10.3891/acta.chem.scand.02-0177, p. 182: type III (4H) on "very
///   impure" material (Al 0.30 %, Fe 0.07 %) "which contained only type III
///   and a little of type I", `σ_obs. = 3.201 + 0.005` (printed with `+`),
///   temperature not stated. The 4H cell is +3.88e-3 denser, 2.5 times the
///   stated 0.005 (1.56e-3 relative); the paper itself says the specific
///   weights of the black type-III-rich carbides "show slightly greater
///   deviations owing to their content of foreign matter". Same page, type II
///   (6H) on "very pure" material (Al and Fe below 0.01 %): `σ_obs. = 3.209 +
///   0.001`; the 6H cell is +1.36e-3 denser, 4.4 times 0.001 (3.1e-4).
/// * P. L. Mulligan et al., ORNL/TM-2023/3059 (2023), doi:10.2172/2001384,
///   section 4.1 and Table 4, p. 7: densities "measured via direct measurement
///   of the parallelepiped specimen dimensions and massing", no temperature
///   and no uncertainty stated, for two 4H wafers: material A (N-type) 3.17
///   g/cm³ and material B (semi-insulating) 3.21 g/cm³. The 4H cell is
///   +1.37e-2 from A and +1.06e-3 from B. B alone would sit within half a unit
///   of its last digit (1.56e-3), but the two 4H wafers of the same method
///   differ by 1.3e-2, about eight times that, so the printed digit does not
///   measure the method's precision and B is not adopted (choosing it over A
///   would be a selection).
#[test]
fn opened_measured_sic_densities_are_recorded_gaps() {
    let h4: fn() -> Lattice = Lattice::silicon_carbide_4h;
    let h6: fn() -> Lattice = Lattice::silicon_carbide_6h;
    let cases = [
        ("4H Lundqvist 1948", h4, 3.201, 3.879e-3),
        ("4H ORNL 2023 A", h4, 3.17, 1.3696e-2),
        ("4H ORNL 2023 B", h4, 3.21, 1.064e-3),
        ("6H Lundqvist 1948", h6, 3.209, 1.362e-3),
        // The archive value of the test above, for scale.
        ("4H Ioffe 3.211", h4, 3.211, 7.52e-4),
    ];
    for (name, lattice, rho, pinned) in cases {
        for rel in density_differences(&lattice(), rho, [14, 6]) {
            assert!(
                (rel - pinned).abs() < 5e-6,
                "{name}: pinned difference {pinned:e}, found {rel:e}"
            );
        }
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

/// A hexagonal crystal runs in the BCA engine (`Bca::with_crystal`): the
/// thermal amplitudes are indexed by the sites of the orthohexagonal search
/// cell (twice the basis), the default `p_max` is the nearest-neighbour
/// distance of the hexagonal crystal, and the lattice-constant temperature
/// of a source that states none is `None` in the run metadata. A plumbing
/// check only: no hexagonal channeling result is validated here.
#[test]
fn hexagonal_crystal_runs_in_the_bca_engine() {
    use lindhard::geometry::Stack;
    use lindhard::ion::bca::{Bca, BcaConfig, Beam, CrystalTarget, SummaryTally, Thermal};
    use lindhard::ion::potential::{Potential, Screening};
    use lindhard::ion::scattering::{ScatteringTable, TableSpec};
    use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
    use lindhard::ion::stopping::Ion;

    let lat = Lattice::silicon_carbide_6h();
    let (a, c) = (lat.lattice_constant(), lat.lattice_constant_c());
    let rho = g_cm3_to_kg_m3(SIC_6H_DENSITY_G_CM3);
    let mut m = Material::from_atom_fractions(&[(14, 1.0), (6, 1.0)], Some(rho)).unwrap();
    for z in [14, 6] {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 5.0).unwrap();
        }
    }
    let st = Stack::new(vec![(m, 1e-6)], None).unwrap();
    let table = ScatteringTable::build(
        &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
        &TableSpec {
            eps_min: 1e-6,
            eps_max: 1e4,
            beta_min: 1e-5,
            beta_max: 1e2,
            per_decade: 8,
        },
    );
    let ls = LindhardScharff::new();
    let o =
        Orientation::new_miller_bravais(&lat, [0, 0, 0, 1], [1, 1, -2, 0], 0.0, 0.0, 0.0).unwrap();
    // 300 K and an arbitrary 1000 K Debye temperature: test inputs, not
    // SiC properties.
    let target = CrystalTarget::new(lat, o).with_thermal(Thermal::new(300.0, 1000.0));
    let beam = Beam {
        ion: Ion::new(7).unwrap(),
        energy_ev: 5e3,
        polar_rad: 7f64.to_radians(),
        azimuth_rad: 0.0,
        count: 20,
    };
    let bca = Bca::new(beam, &st, BcaConfig::new(5.0, 2.0), &ls, &table)
        .unwrap()
        .with_crystal(target, &[0])
        .unwrap();

    let meta = &bca.crystal_metadata()[0];
    assert_eq!(meta.lattice_constant_m, a);
    assert_eq!(meta.lattice_constant_temperature_k, None);
    // The Si-C bond of the near-ideal 6H crystal: a sqrt(3/8) in the basal
    // network, about 3c/24 along c (both within 1 %).
    let bond = a * (3.0f64 / 8.0).sqrt();
    assert!(
        (meta.p_max_m / bond - 1.0).abs() < 1e-2,
        "p_max {} m, bond {bond} m (c/8 = {})",
        meta.p_max_m,
        c / 8.0
    );

    let t = bca.run(|| SummaryTally::new(1e-9, 1000)).unwrap();
    assert_eq!(t.histories, 20);
}
