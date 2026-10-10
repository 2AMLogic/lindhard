//! Tests for the hydrogenic inner-shell ELF builder
//! (`lindhard::electron::inelastic::shell_elf`).
//!
//! The binding energies and occupancies are the committed EADL2017 table
//! (`SubshellBindingTable::eadl2017()`); the atom densities and the valence
//! ELF (a Drude-Lorentz plasmon) are **synthetic**, round numbers chosen to
//! exercise the code, not material data.

use lindhard::constants::{
    BOHR_RADIUS, COULOMB_E2, ELECTRON_MASS, ELEMENTARY_CHARGE, FINE_STRUCTURE, HARTREE_ENERGY,
    HBAR, SPEED_OF_LIGHT,
};
use lindhard::electron::data::{OpticalElf, Subshell, SubshellBindingTable};
use lindhard::electron::inelastic::{
    hydrogenic_2p_oscillator_strength_density_per_ev,
    hydrogenic_2s_oscillator_strength_density_per_ev,
    hydrogenic_k_oscillator_strength_density_per_ev, hydrogenic_shell_elf, hydrogenic_shell_elfs,
    DrudeLorentz, DrudeLorentzOscillator, ShellElfGrid, ShellResolvedChannels, SinglePolePenn,
    SumRuleReport,
};
use rayon::prelude::*;
use std::f64::consts::PI;

/// A synthetic atom density, m⁻³ (a round number, not a material's).
const ATOMS_PER_M3: f64 = 5.0e28;
const MATERIAL: &str = "synthetic target (EADL edges, made-up density)";

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

/// The Stobbe K-shell cross section per electron of a hydrogen-like ion of
/// charge `z` at photon energy `eps` times its ionization energy, m²
/// (Mikhailov, Nefiodov and Plunien, Phys. Lett. A 368, 391 (2007),
/// arXiv:0704.2180, eq. (1)), written independently of the library.
fn stobbe_sigma_m2(z: f64, eps: f64) -> f64 {
    let xi = 1.0 / (eps - 1.0).sqrt();
    let acot = (1.0 / xi).atan();
    FINE_STRUCTURE * BOHR_RADIUS * BOHR_RADIUS * 512.0 * PI * PI / (3.0 * z * z)
        * (-4.0 * xi * acot).exp()
        / ((1.0 + xi.powi(-2)).powi(4) * (1.0 - (-2.0 * PI * xi).exp()))
}

/// Salvat, arXiv:2605.22442, eq. (6.2), `df/dW = m_e c σ / (2π² e² ħ)`
/// (Gaussian `e²`, i.e. `e²/(4π ε₀)` in SI), in eV⁻¹.
fn salvat_df_dw_per_ev(sigma_m2: f64) -> f64 {
    let per_j = ELECTRON_MASS * SPEED_OF_LIGHT * sigma_m2 / (2.0 * PI * PI * COULOMB_E2 * HBAR);
    per_j * ELEMENTARY_CHARGE
}

/// The library's closed form against eq. (6.2) applied to eq. (1) with CODATA
/// constants, for hydrogen-like ions (`I = Z² E_h / 2`).
#[test]
fn oscillator_strength_density_is_eq_6_2_of_the_stobbe_cross_section() {
    for z in [1.0, 6.0, 29.0] {
        let i_ev = z * z * HARTREE_ENERGY / ELEMENTARY_CHARGE / 2.0;
        for eps in [1.000_001, 1.01, 1.5, 2.0, 10.0, 137.0, 1.0e4] {
            let want = salvat_df_dw_per_ev(stobbe_sigma_m2(z, eps));
            let got = hydrogenic_k_oscillator_strength_density_per_ev(i_ev, eps * i_ev);
            assert!(
                rel(got, want) < 1e-9,
                "Z = {z}, eps = {eps}: {got:e} vs {want:e}"
            );
        }
        // Eq. (3): the threshold value α a0² 2⁹π² / (3 e⁴ Z²).
        let thr = FINE_STRUCTURE * BOHR_RADIUS * BOHR_RADIUS * 512.0 * PI * PI
            / (3.0 * 1f64.exp().powi(4) * z * z);
        let got = hydrogenic_k_oscillator_strength_density_per_ev(i_ev, i_ev);
        assert!(rel(got, salvat_df_dw_per_ev(thr)) < 1e-9, "Z = {z}");
        // Eq. (2): the high-energy limit α a0² 2⁸π / (3 Z²) ε^(-7/2).
        let eps: f64 = 1.0e8;
        let asym = FINE_STRUCTURE * BOHR_RADIUS * BOHR_RADIUS * 256.0 * PI / (3.0 * z * z)
            * eps.powf(-3.5);
        let got = hydrogenic_k_oscillator_strength_density_per_ev(i_ev, eps * i_ev);
        assert!(rel(got, salvat_df_dw_per_ev(asym)) < 1e-3, "Z = {z}");
    }
}

/// `∫_B^∞ df/dW dW` per electron, by Simpson's rule in `u = ln(W/B - 1)`
/// on a fine grid (independent of the library's sampling grid).
fn continuum_strength_per_electron(b: f64, top_factor: f64) -> f64 {
    let (lo, hi) = ((1e-14f64).ln(), (top_factor - 1.0).ln());
    let n = 200_000;
    let f = |u: f64| {
        let x = u.exp();
        hydrogenic_k_oscillator_strength_density_per_ev(b, b * (1.0 + x)) * b * x
    };
    let h = (hi - lo) / n as f64;
    let mut s = 0.0;
    for k in 0..n {
        let u0 = lo + h * k as f64;
        s += h / 6.0 * (f(u0) + 4.0 * f(u0 + 0.5 * h) + f(u0 + h));
    }
    s
}

/// The Bethe-sum check: the `N_eff` of the built ELF (Shinotsuka et al.
/// 2017 eq. (18), `SumRuleReport`) is the shell occupancy times the
/// oscillator strength per electron the formula implies above the edge.
#[test]
fn shell_elf_f_sum_is_the_occupancy_times_the_continuum_oscillator_strength() {
    let table = SubshellBindingTable::eadl2017();
    // Per electron, from the edge to infinity: independent of B.
    let s_inf = continuum_strength_per_electron(538.0, 1.0e12);
    let s_inf_h = continuum_strength_per_electron(13.6, 1.0e12);
    assert!(rel(s_inf, s_inf_h) < 1e-12);
    // The value the formula gives (a computed number, not a reference one).
    assert!((s_inf - 0.435).abs() < 5e-4, "{s_inf}");
    let grid = ShellElfGrid::default();
    for (z, occupancy) in [(1u8, 1.0), (6, 2.0), (8, 2.0), (29, 2.0)] {
        let (shell, elf) =
            hydrogenic_shell_elf(MATERIAL, &table, z, Subshell::K, ATOMS_PER_M3, grid).unwrap();
        let b = shell.binding_energy_ev;
        let s_grid = continuum_strength_per_electron(b, grid.max_energy_factor);
        // Truncation at 10^3 B.
        assert!(rel(s_grid, s_inf) < 1e-6, "Z = {z}");
        let r = SumRuleReport::with_target_density(&elf, ATOMS_PER_M3).unwrap();
        let n_eff = r.effective_electrons().unwrap();
        let want = occupancy * s_grid;
        // The linear interpolation between the log-spaced samples.
        assert!(rel(n_eff, want) < 1.5e-4, "Z = {z}: {n_eff} vs {want}");
    }
    // A finer grid converges on the formula.
    let fine = ShellElfGrid {
        points_per_decade: 1000,
        max_energy_factor: 1.0e3,
    };
    let (_, elf) =
        hydrogenic_shell_elf(MATERIAL, &table, 8, Subshell::K, ATOMS_PER_M3, fine).unwrap();
    let n_eff = SumRuleReport::with_target_density(&elf, ATOMS_PER_M3)
        .unwrap()
        .effective_electrons()
        .unwrap();
    assert!(rel(n_eff, 2.0 * s_inf) < 1.5e-5, "{n_eff}");
    // N_eff does not depend on the density (the ELF scales with it).
    let (_, elf2) =
        hydrogenic_shell_elf(MATERIAL, &table, 8, Subshell::K, 3.0 * ATOMS_PER_M3, fine).unwrap();
    let n2 = SumRuleReport::with_target_density(&elf2, 3.0 * ATOMS_PER_M3)
        .unwrap()
        .effective_electrons()
        .unwrap();
    assert!(rel(n2, n_eff) < 1e-12);
}

#[test]
fn shell_elf_starts_at_the_edge_and_carries_its_provenance() {
    let table = SubshellBindingTable::eadl2017();
    let shells = hydrogenic_shell_elfs(
        MATERIAL,
        &table,
        &[
            (8, Subshell::K, ATOMS_PER_M3),
            (1, Subshell::K, 2.0 * ATOMS_PER_M3),
        ],
        ShellElfGrid::default(),
    )
    .unwrap();
    assert_eq!(shells.len(), 2);
    let (o, elf) = &shells[0];
    assert_eq!((o.z, o.subshell), (8, Subshell::K));
    assert_eq!(o.binding_energy_ev, 538.0);
    assert_eq!(elf.energy_range_ev().0, 538.0);
    assert!(elf.elf(537.9).is_err()); // nothing below the edge
    assert!(elf.elf(538.0).unwrap() > 0.0); // the step of dV2022 eq. (2)
    assert_eq!(elf.material(), MATERIAL);
    for needle in [
        "Stobbe",
        "0704.2180",
        "2605.22442",
        "6121",
        "EADL2017",
        "Z = 8 K",
    ] {
        assert!(elf.provenance().contains(needle), "{needle}");
    }
    assert_eq!(shells[1].0.binding_energy_ev, 13.6);
    // The edge value: (π/2) (ħΩ_a)² q (2⁷/(3B)) e⁻⁴ / B.
    let hw2 = HBAR * HBAR * ATOMS_PER_M3 * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE
        / (lindhard::constants::VACUUM_PERMITTIVITY * ELECTRON_MASS)
        / (ELEMENTARY_CHARGE * ELEMENTARY_CHARGE);
    let want = PI / 2.0 * hw2 * 2.0 * 128.0 / (3.0 * 538.0) * (-4.0f64).exp() / 538.0;
    assert!(rel(elf.elf(538.0).unwrap(), want) < 1e-12);
}

#[test]
fn builder_rejects_what_it_cannot_build() {
    let table = SubshellBindingTable::eadl2017();
    let g = ShellElfGrid::default();
    // The M shell and above have no formula here (Cu has M1 ... M5).
    for label in ["M1", "M2", "M3", "M4", "M5", "N1"] {
        let s = Subshell::from_label(label).unwrap();
        let e = hydrogenic_shell_elf(MATERIAL, &table, 29, s, ATOMS_PER_M3, g).unwrap_err();
        assert!(e.to_string().contains("K, L1, L2 and L3"), "{label}: {e}");
    }
    let m1 = Subshell::from_label("M1").unwrap();
    // Hydrogen has no L shell in the table.
    let l1 = Subshell::from_label("L1").unwrap();
    assert!(hydrogenic_shell_elf(MATERIAL, &table, 1, l1, ATOMS_PER_M3, g).is_err());
    assert!(hydrogenic_shell_elf(MATERIAL, &table, 93, Subshell::K, ATOMS_PER_M3, g).is_err());
    for n in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(hydrogenic_shell_elf(MATERIAL, &table, 8, Subshell::K, n, g).is_err());
    }
    let bad = ShellElfGrid {
        points_per_decade: 0,
        max_energy_factor: 10.0,
    };
    assert!(hydrogenic_shell_elf(MATERIAL, &table, 8, Subshell::K, ATOMS_PER_M3, bad).is_err());
    assert!(hydrogenic_shell_elf("", &table, 8, Subshell::K, ATOMS_PER_M3, g).is_err());
    // One bad shell fails the whole list.
    assert!(hydrogenic_shell_elfs(
        MATERIAL,
        &table,
        &[(29, Subshell::K, ATOMS_PER_M3), (29, m1, ATOMS_PER_M3)],
        g
    )
    .is_err());
}

/// A per-electron oscillator-strength density `df/dW(B, W)`, eV⁻¹.
type Density = fn(f64, f64) -> f64;

/// `∫_B^(top B) df/dW dW` per electron for any of the hydrogenic formulas,
/// by Simpson's rule in `u = ln(W/B - 1)` on a fine grid (independent of
/// the library's sampling grid).
fn continuum_strength_of(density: Density, b: f64, top_factor: f64) -> f64 {
    let (lo, hi) = ((1e-14f64).ln(), (top_factor - 1.0).ln());
    let n = 200_000;
    let f = |u: f64| {
        let x = u.exp();
        density(b, b * (1.0 + x)) * b * x
    };
    let h = (hi - lo) / n as f64;
    let mut s = 0.0;
    for k in 0..n {
        let u0 = lo + h * k as f64;
        s += h / 6.0 * (f(u0) + 4.0 * f(u0 + 0.5 * h) + f(u0 + h));
    }
    s
}

/// Bound-free Gaunt factors of W. J. Karzas and R. Latter, Astrophys. J.
/// Suppl. 6, 167 (1961), Table 1 (p. 178, "n = 1, n = 2"), transcribed
/// 2026-10-09 from the page image of the ADS scan
/// (<https://articles.adsabs.harvard.edu/pdf/1961ApJS....6..167K>): the
/// electron energy `E / Z² Ry` and the columns 1s, 2s and 2p, four printed
/// digits each. The rows below `E = 10⁻⁴`, which repeat the last one, are
/// left out.
const KARZAS_LATTER_TABLE_1: [(f64, f64, f64, f64); 39] = [
    (0.1000e13, 0.6928e-5, 0.2771e-4, 0.6928e-17),
    (0.1111e12, 0.2078e-4, 0.8314e-4, 0.1871e-15),
    (0.1000e11, 0.6928e-4, 0.2771e-3, 0.6928e-14),
    (0.1111e10, 0.2078e-3, 0.8313e-3, 0.1870e-12),
    (0.1000e9, 0.6926e-3, 0.2770e-2, 0.6926e-11),
    (0.1111e8, 0.2076e-2, 0.8306e-2, 0.1869e-9),
    (0.1000e7, 0.6906e-2, 0.2763e-1, 0.6906e-8),
    (0.1111e6, 0.2059e-1, 0.8236e-1, 0.1853e-6),
    (0.4000e5, 0.3410e-1, 0.1364e0, 0.8525e-6),
    (0.2041e5, 0.4745e-1, 0.1898e0, 0.2325e-5),
    (0.1000e5, 0.6715e-1, 0.2686e0, 0.6714e-5),
    (0.4444e4, 0.9917e-1, 0.3966e0, 0.2231e-4),
    (0.2500e4, 0.1302e0, 0.5207e0, 0.5206e-4),
    (0.1111e4, 0.1894e0, 0.7572e0, 0.1703e-3),
    (0.4000e3, 0.2971e0, 0.1187e1, 0.7411e-3),
    (0.2041e3, 0.3918e0, 0.1563e1, 0.1912e-2),
    (0.1000e3, 0.5129e0, 0.2042e1, 0.5087e-2),
    (0.4444e2, 0.6687e0, 0.2645e1, 0.1477e-1),
    (0.2500e2, 0.7800e0, 0.3059e1, 0.3019e-1),
    (0.1600e2, 0.8585e0, 0.3331e1, 0.5100e-1),
    (0.1111e2, 0.9129e0, 0.3497e1, 0.7642e-1),
    (0.6250e1, 0.9729e0, 0.3612e1, 0.1373e0),
    (0.4000e1, 0.9939e0, 0.3553e1, 0.2055e0),
    (0.2778e1, 0.9948e0, 0.3408e1, 0.2752e0),
    (0.2041e1, 0.9856e0, 0.3225e1, 0.3423e0),
    (0.1562e1, 0.9721e0, 0.3031e1, 0.4045e0),
    (0.1235e1, 0.9571e0, 0.2842e1, 0.4607e0),
    (0.1000e1, 0.9423e0, 0.2664e1, 0.5106e0),
    (0.6944e0, 0.9157e0, 0.2354e1, 0.5925e0),
    (0.4444e0, 0.8850e0, 0.2000e1, 0.6784e0),
    (0.2500e0, 0.8531e0, 0.1626e1, 0.7587e0),
    (0.1111e0, 0.8246e0, 0.1279e1, 0.8192e0),
    (0.4000e-1, 0.8076e0, 0.1067e1, 0.8459e0),
    (0.2041e-1, 0.8026e0, 0.1003e1, 0.8518e0),
    (0.1000e-1, 0.7999e0, 0.9686e0, 0.8545e0),
    (0.4444e-2, 0.7985e0, 0.9498e0, 0.8558e0),
    (0.2500e-2, 0.7980e0, 0.9431e0, 0.8562e0),
    (0.1111e-2, 0.7976e0, 0.9384e0, 0.8565e0),
    (0.1000e-3, 0.7973e0, 0.9349e0, 0.8567e0),
];

/// The Gaunt factor the library's `df/dW` implies: Karzas and Latter 1961
/// eq. (40), `g = σ/σ^K`, with Kramers' cross section of eq. (39),
/// `σ^K = (2⁴/(3√3)) (e²/(m c ν)) (1/n) (ρ²/(1+ρ²))²`, and Salvat's eq.
/// (6.2), which gives `σ = (π e²/(m c ν)) W df/dW`; `(1+ρ²)/ρ² = W/I` and
/// `W/I = 1 + n² E/(Z² Ry)`. So `g = (3√3 π n/16) (W/I)² W df/dW`.
fn gaunt_factor(density: Density, n: f64, electron_energy_ry: f64) -> f64 {
    let b = 250.0;
    let eps = 1.0 + n * n * electron_energy_ry;
    let w = eps * b;
    3.0 * 3f64.sqrt() * PI * n / 16.0 * eps * eps * w * density(b, w)
}

/// The closed forms reproduce every 1s, 2s and 2p entry of Karzas and
/// Latter's Table 1, from the edge to `E = 10¹² Z² Ry`, to the table's four
/// printed digits. The tolerance 1e-3 is half a unit of the fourth digit (at
/// most 5e-4 of the entry) plus the rounding of the printed energy to four
/// digits (at most 4.5e-4 of `E`; `|d ln g / d ln E|` is at most 3/2, the
/// high-energy slope of the 2p column).
#[test]
fn densities_reproduce_the_karzas_latter_table_1_gaunt_factors() {
    let mut worst = 0.0f64;
    for (e, g1s, g2s, g2p) in KARZAS_LATTER_TABLE_1 {
        let cases: [(Density, f64, f64, &str); 3] = [
            (
                hydrogenic_k_oscillator_strength_density_per_ev,
                1.0,
                g1s,
                "1s",
            ),
            (
                hydrogenic_2s_oscillator_strength_density_per_ev,
                2.0,
                g2s,
                "2s",
            ),
            (
                hydrogenic_2p_oscillator_strength_density_per_ev,
                2.0,
                g2p,
                "2p",
            ),
        ];
        for (density, n, want, name) in cases {
            let got = gaunt_factor(density, n, e);
            worst = worst.max(rel(got, want));
            assert!(
                rel(got, want) < 1e-3,
                "{name} at E = {e}: {got:e} vs {want:e}"
            );
        }
    }
    // The agreement is in fact twice as good as the bound.
    assert!(worst < 5e-4, "{worst}");
    // The edge values of the table (its rows from E = 1.111e-5 down): .7973,
    // .9346 and .8567, against the closed-form limits (2⁷/3) e⁻⁴,
    // (2¹²/3) e⁻⁸ and (11 · 2¹⁰/9) e⁻⁸ of W df/dW at W = I.
    let k = 3.0 * 3f64.sqrt() * PI / 16.0;
    for (got, want) in [
        (k * 128.0 / 3.0 * (-4.0f64).exp(), 0.7973),
        (2.0 * k * 4096.0 / 3.0 * (-8.0f64).exp(), 0.9346),
        (2.0 * k * 11264.0 / 9.0 * (-8.0f64).exp(), 0.8567),
    ] {
        assert!((got - want).abs() < 5.1e-5, "{got} vs {want}");
    }
    for (density, n, want) in [
        (
            hydrogenic_2s_oscillator_strength_density_per_ev as Density,
            2.0,
            0.9346,
        ),
        (
            hydrogenic_2p_oscillator_strength_density_per_ev,
            2.0,
            0.8567,
        ),
    ] {
        assert!((gaunt_factor(density, n, 0.0) - want).abs() < 5.1e-5);
    }
}

/// The Bethe-sum check for the L subshells: the `N_eff` of each built ELF is
/// the subshell occupancy times the continuum oscillator strength per
/// electron its formula implies, and the 2p strength is shared between L2
/// and L3 as their occupancies.
#[test]
fn l_shell_elf_f_sums_are_the_occupancies_times_the_continuum_oscillator_strengths() {
    let table = SubshellBindingTable::eadl2017();
    let [l1, l2, l3] = ["L1", "L2", "L3"].map(|l| Subshell::from_label(l).unwrap());
    let s2s: Density = hydrogenic_2s_oscillator_strength_density_per_ev;
    let s2p: Density = hydrogenic_2p_oscillator_strength_density_per_ev;
    // Per electron, from the edge to infinity: independent of B.
    let s2s_inf = continuum_strength_of(s2s, 149.7, 1.0e12);
    let s2p_inf = continuum_strength_of(s2p, 99.2, 1.0e12);
    assert!(rel(continuum_strength_of(s2s, 3.4, 1.0e12), s2s_inf) < 1e-12);
    assert!(rel(continuum_strength_of(s2p, 3.4, 1.0e12), s2p_inf) < 1e-12);
    // The values the formulas give (computed numbers, not reference ones).
    assert!((s2s_inf - 0.351).abs() < 5e-4, "{s2s_inf}");
    assert!((s2p_inf - 0.191).abs() < 5e-4, "{s2p_inf}");
    let grid = ShellElfGrid::default();
    let n_eff = |z: u8, s: Subshell, grid: ShellElfGrid| {
        let (shell, elf) = hydrogenic_shell_elf(MATERIAL, &table, z, s, ATOMS_PER_M3, grid)
            .unwrap_or_else(|e| panic!("Z = {z} {}: {e}", s.label()));
        assert_eq!((shell.z, shell.subshell), (z, s));
        let b = shell.binding_energy_ev;
        assert_eq!(elf.energy_range_ev().0, b);
        assert!(elf.elf(b * (1.0 - 1e-6)).is_err()); // nothing below the edge
        assert!(elf.elf_values().iter().all(|&x| x > 0.0));
        let r = SumRuleReport::with_target_density(&elf, ATOMS_PER_M3).unwrap();
        (r.effective_electrons().unwrap(), b)
    };
    let occupancy = |z: u8, s: Subshell| table.atom(z).unwrap().shell(s).unwrap().occupancy();
    // Si and Cu have closed L shells; B and N have a partly filled 2p.
    for z in [5u8, 7, 14, 29] {
        let (q1, q2) = (occupancy(z, l1), occupancy(z, l2));
        let q3 = table
            .atom(z)
            .unwrap()
            .shell(l3)
            .map_or(0.0, |s| s.occupancy());
        if z >= 14 {
            assert_eq!((q1, q2, q3), (2.0, 2.0, 4.0), "Z = {z}");
        }
        let mut p_total = 0.0;
        for (s, q, density, s_inf) in [
            (l1, q1, s2s, s2s_inf),
            (l2, q2, s2p, s2p_inf),
            (l3, q3, s2p, s2p_inf),
        ] {
            if q == 0.0 {
                continue;
            }
            let (got, b) = n_eff(z, s, grid);
            let s_grid = continuum_strength_of(density, b, grid.max_energy_factor);
            // Truncation at 10^3 B.
            assert!(rel(s_grid, s_inf) < 1e-6, "Z = {z} {}", s.label());
            // The linear interpolation between the log-spaced samples.
            let want = q * s_grid;
            assert!(
                rel(got, want) < 1.5e-4,
                "Z = {z} {}: {got} vs {want}",
                s.label()
            );
            if s != l1 {
                p_total += got;
            }
        }
        // L2 + L3 is one 2p shell of q(L2) + q(L3) electrons ...
        assert!(rel(p_total, (q2 + q3) * s2p_inf) < 1.5e-4, "Z = {z}");
        // ... shared as the occupancies (the binding energies differ, the
        // dimensionless integral does not).
        if q3 > 0.0 {
            let ratio = n_eff(z, l2, grid).0 / n_eff(z, l3, grid).0;
            assert!(rel(ratio, q2 / q3) < 1e-9, "Z = {z}: {ratio}");
        }
    }
    // A finer grid converges on the formulas.
    let fine = ShellElfGrid {
        points_per_decade: 1000,
        max_energy_factor: 1.0e3,
    };
    assert!(rel(n_eff(14, l1, fine).0, 2.0 * s2s_inf) < 1.5e-5);
    assert!(rel(n_eff(14, l3, fine).0, 4.0 * s2p_inf) < 1.5e-5);
}

#[test]
fn l_shell_elfs_start_at_their_edges_and_carry_their_provenance() {
    let table = SubshellBindingTable::eadl2017();
    let [l1, l2, l3] = ["L1", "L2", "L3"].map(|l| Subshell::from_label(l).unwrap());
    let shells = hydrogenic_shell_elfs(
        MATERIAL,
        &table,
        &[
            (14, Subshell::K, ATOMS_PER_M3),
            (14, l1, ATOMS_PER_M3),
            (14, l2, ATOMS_PER_M3),
            (14, l3, ATOMS_PER_M3),
        ],
        ShellElfGrid::default(),
    )
    .unwrap();
    assert_eq!(shells.len(), 4);
    // The K shell keeps its own citation.
    assert!(shells[0].1.provenance().contains("Stobbe"));
    assert!(!shells[0].1.provenance().contains("Karzas"));
    let hw2 = HBAR * HBAR * ATOMS_PER_M3 * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE
        / (lindhard::constants::VACUUM_PERMITTIVITY * ELECTRON_MASS)
        / (ELEMENTARY_CHARGE * ELEMENTARY_CHARGE);
    // Edge values of W df/dW per electron: (2¹²/3) e⁻⁸ (2s), (11 · 2¹⁰/9) e⁻⁸ (2p).
    let e8 = (-8.0f64).exp();
    for ((shell, elf), (s, edge)) in shells[1..].iter().zip([
        (l1, 4096.0 / 3.0 * e8),
        (l2, 11264.0 / 9.0 * e8),
        (l3, 11264.0 / 9.0 * e8),
    ]) {
        assert_eq!((shell.z, shell.subshell), (14, s));
        let binding = table.atom(14).unwrap().shell(s).unwrap();
        let b = binding.binding_energy_ev();
        assert_eq!(shell.binding_energy_ev, b);
        assert_eq!(elf.energy_range_ev().0, b);
        let want = PI / 2.0 * hw2 * binding.occupancy() * edge / (b * b);
        assert!(rel(elf.elf(b).unwrap(), want) < 1e-12, "{}", s.label());
        for needle in [
            "Karzas and Latter",
            "(36), (37)",
            "2605.22442",
            "6121",
            "EADL2017",
            &format!("Z = 14 {}", s.label()),
        ] {
            assert!(elf.provenance().contains(needle), "{needle}");
        }
        assert!(!elf.provenance().contains("Stobbe"));
    }
    // The edges are ordered K > L1 > L2 >= L3.
    let b: Vec<f64> = shells.iter().map(|s| s.0.binding_energy_ev).collect();
    assert!(b[0] > b[1] && b[1] > b[2] && b[2] >= b[3], "{b:?}");
}

/// The L-subshell tables are bit-identical whatever the number of threads.
#[test]
fn built_l_shells_are_bit_identical_across_thread_counts() {
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let table = SubshellBindingTable::eadl2017();
            let jobs: Vec<(u8, Subshell)> = (13u8..=30)
                .flat_map(|z| ["L1", "L2", "L3"].map(|l| (z, Subshell::from_label(l).unwrap())))
                .collect();
            jobs.into_par_iter()
                .map(|(z, s)| {
                    let grid = ShellElfGrid::default();
                    hydrogenic_shell_elf(MATERIAL, &table, z, s, ATOMS_PER_M3, grid)
                        .unwrap()
                        .1
                        .elf_values()
                        .to_vec()
                })
                .collect::<Vec<Vec<f64>>>()
        })
    };
    let one = run(1);
    assert_eq!(one.len(), 54);
    for t in [2, 8] {
        let other = run(t);
        assert_eq!(one.len(), other.len());
        for (a, b) in one.iter().zip(&other) {
            assert!(a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits()));
        }
    }
}

fn drude() -> DrudeLorentz {
    DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)]).unwrap()
}

/// The O K shell (538 eV) to 30 B and a valence Drude ELF sampled on one
/// merged grid: Drude knots below the edge, a knot `B (1 - 1e-9)`, then the
/// shell's own knots. Returns (valence, shell, sum) on that grid.
fn valence_shell_and_sum() -> (
    OpticalElf,
    (lindhard::electron::inelastic::InnerShell, OpticalElf),
    OpticalElf,
) {
    let table = SubshellBindingTable::eadl2017();
    let grid = ShellElfGrid {
        points_per_decade: 100,
        max_energy_factor: 30.0,
    };
    let (shell, selph) =
        hydrogenic_shell_elf(MATERIAL, &table, 8, Subshell::K, ATOMS_PER_M3, grid).unwrap();
    let b = shell.binding_energy_ev;
    let n = 200;
    let (lo, hi): (f64, f64) = (0.05, b * (1.0 - 1e-6));
    let mut knots: Vec<f64> = (0..n)
        .map(|i| lo * (hi / lo).powf(i as f64 / (n - 1) as f64))
        .collect();
    knots.push(b * (1.0 - 1e-9));
    knots.extend_from_slice(selph.energy_ev());
    let d = drude();
    let v: Vec<f64> = knots.iter().map(|&e| d.elf(e)).collect();
    let valence =
        OpticalElf::new("synthetic Drude", "synthetic", knots.clone(), v.clone()).unwrap();
    let sum: Vec<f64> = knots
        .iter()
        .zip(&v)
        .map(|(&e, &x)| x + if e >= b { selph.elf(e).unwrap() } else { 0.0 })
        .collect();
    let sum = OpticalElf::new("synthetic sum", "synthetic", knots, sum).unwrap();
    (valence, (shell, selph), sum)
}

/// Without exchange the channels (valence + built K shell) add up to the
/// model of the summed ELF (the SPA is linear in the ELF); the only
/// difference is the 1e-9 B ramp the summed table needs at the edge.
#[test]
fn built_shells_feed_the_channels_and_add_up_to_the_summed_model() {
    let (valence, shell, sum) = valence_shell_and_sum();
    let c = ShellResolvedChannels::new(SinglePolePenn::new(valence), vec![shell]).unwrap();
    let total = SinglePolePenn::new(sum);
    let mut opened = false;
    for e in [300.0, 2000.0, 6000.0, 12_000.0] {
        let r = c.inverse_imfps(e).unwrap();
        opened |= r.shells_per_m[0] > 0.0;
        let want = total.imfp_and_stopping(e).unwrap().inverse_imfp_per_m;
        assert!(rel(r.total_per_m(), want) < 1e-6, "{e} eV");
        for loss in [37.0, 600.0, 1500.0] {
            if loss >= e {
                continue;
            }
            let d = c.diimfps_per_m_ev(e, loss).unwrap();
            let want = total.diimfp_per_m_ev(e, loss).unwrap();
            assert!(
                (d.total_per_m_ev() - want).abs() <= 1e-6 * want.abs().max(1e-300),
                "{e} eV, {loss} eV"
            );
        }
    }
    assert!(opened, "the K channel never opened");
    // Below the edge the shell's DIIMFP is zero.
    let d = c.diimfps_per_m_ev(6000.0, 500.0).unwrap();
    assert_eq!(d.shells_per_m_ev[0], 0.0);
    assert!(d.valence_per_m_ev > 0.0);
}

/// Bit-identical results whatever the number of threads.
#[test]
fn built_shells_are_bit_identical_across_thread_counts() {
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let table = SubshellBindingTable::eadl2017();
            let elfs: Vec<Vec<f64>> = (1u8..=30)
                .into_par_iter()
                .map(|z| {
                    hydrogenic_shell_elf(
                        MATERIAL,
                        &table,
                        z,
                        Subshell::K,
                        ATOMS_PER_M3,
                        ShellElfGrid::default(),
                    )
                    .unwrap()
                    .1
                    .elf_values()
                    .to_vec()
                })
                .collect();
            let (valence, shell, _) = valence_shell_and_sum();
            let c = ShellResolvedChannels::new(SinglePolePenn::new(valence), vec![shell]).unwrap();
            let imfps: Vec<f64> = [2000.0, 6000.0, 12_000.0]
                .par_iter()
                .map(|&e| c.inverse_imfps(e).unwrap().shells_per_m[0])
                .collect();
            (elfs, imfps)
        })
    };
    let one = run(1);
    for t in [2, 8] {
        let other = run(t);
        assert_eq!(one.0.len(), other.0.len());
        for (a, b) in one.0.iter().zip(&other.0) {
            assert!(a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits()));
        }
        assert!(one
            .1
            .iter()
            .zip(&other.1)
            .all(|(x, y)| x.to_bits() == y.to_bits()));
    }
}
