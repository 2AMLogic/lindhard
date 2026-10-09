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
    let l1 = Subshell::from_label("L1").unwrap();
    assert!(hydrogenic_shell_elf(MATERIAL, &table, 14, l1, ATOMS_PER_M3, g).is_err());
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
        &[(8, Subshell::K, ATOMS_PER_M3), (8, l1, ATOMS_PER_M3)],
        g
    )
    .is_err());
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
