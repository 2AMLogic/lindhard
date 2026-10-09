//! Inner-shell optical ELFs built from the hydrogenic K-shell
//! photoionization formula, for [`ShellResolvedChannels`].
//!
//! # Sources (all opened 2026-10-08)
//!
//! * **The shell ELF.** P. de Vera, S. Taioli, P. E. Trevisanutto, M. Dapor,
//!   I. Abril, S. Simonucci and R. Garcia-Molina, Int. J. Mol. Sci. 23, 6121
//!   (2022), doi:10.3390/ijms23116121 (open access, PMC9181504; "dV2022", as
//!   in [`super::inner_shell`]), eq. (2): the inner-shell part of the ELF of a
//!   target with `N` molecules per unit volume, `ν_j` atoms of element `j`
//!   per molecule, is
//!
//!   ```text
//!   Im[-1/ε(k, E)]_inner = (2π² ħ² e² N / E) Σ_j ν_j Σ_nℓ df_nℓ^j(k, E)/dE Θ(E - B_nℓ^j),
//!   ```
//!
//!   written in atomic units (the electron mass is set to 1). The prefactor
//!   with the mass restored is fixed by the f-sum rule of H. Shinotsuka et
//!   al., Surf. Interface Anal. 49, 238 (2017), doi:10.1002/sia.6123 (open
//!   copy PMC5524379), eq. (18), `N_eff = (2 / (π (ħΩ_a)²)) ∫ W ELF dW` with
//!   `(ħΩ_a)² = ħ² n_a e² / (ε₀ m_e)` (SI, the form [`super::sum_rules`]
//!   uses): requiring `N_eff = ∫ df` for one shell of atoms at density `n_a`
//!   gives, in the optical limit `k -> 0`,
//!
//!   ```text
//!   ELF_nℓ(W) = (π/2) (ħΩ_a)² (1/W) df_nℓ/dW,   W >= B_nℓ,   zero below,
//!   ```
//!
//!   which is dV2022 eq. (2) for one shell (`(π/2) 4π = 2π²` with the
//!   Gaussian `Ω² = 4π n e²/m`). The shell's oscillator strength is then
//!   exactly the shell's `N_eff` (tested).
//!
//! * **The oscillator-strength density.** F. Salvat, "Inelastic collisions of
//!   fast charged particles with atoms. Relativistic plane-wave Born
//!   approximation", arXiv:2605.22442 (2026), eq. (6.2): the optical
//!   oscillator strength is proportional to the dipole photoabsorption cross
//!   section, `df/dW = (m_e c / (2π² e² ħ)) σ_ph(W)` (Gaussian `e²`); and eq.
//!   (6.1a): a subshell with `q_a` of its `2|κ_a|` places filled contributes
//!   `q_a / 2|κ_a|` times the closed-subshell value, i.e. the oscillator
//!   strength scales with the occupancy.
//!
//! * **The hydrogenic K-shell cross section.** M. Stobbe, Ann. Phys. (Leipzig)
//!   7, 661 (1930) (not read), as printed by A. I. Mikhailov, A. V. Nefiodov
//!   and G. Plunien, Phys. Lett. A 368, 391 (2007),
//!   doi:10.1016/j.physleta.2007.04.027 (open copy arXiv:0704.2180), eq. (1):
//!   for one electron in the 1s state of a hydrogen-like ion of charge `Z`,
//!
//!   ```text
//!   σ_K = α a₀² (2⁹ π² / (3 Z²)) exp(-4ξ cot⁻¹ ξ) / ((1 + ξ⁻²)⁴ [1 - exp(-2πξ)]),
//!   ```
//!
//!   with `ξ = 1/sqrt(ε_γ - 1)`, `ε_γ = W/I` the photon energy over the
//!   ionization energy `I = Z² E_h / 2`, valid for `1 <= ε_γ << 2 (αZ)⁻²`
//!   (dipole approximation, nonrelativistic). Its eq. (3) is the threshold
//!   value (`ξ -> ∞`, `exp(-4)`), its eq. (2) the high-energy limit
//!   `α a₀² (2⁸ π / (3 Z²)) ε_γ^(-7/2)`, and its eq. (6) says that for two
//!   independent K electrons the cross section is `2 σ_K`.
//!
//! # The formula used
//!
//! Inserting eq. (1) into eq. (6.2) with `e² = α ħ c`, `a₀ = ħ / (m_e c α)`
//! and `I = Z² E_h / 2`, `E_h = ħ² / (m_e a₀²)` (our algebra; the test
//! `oscillator_strength_density_is_eq_6_2_of_the_stobbe_cross_section`
//! redoes it numerically with CODATA constants), the oscillator-strength
//! density per K electron is
//!
//! ```text
//! df_K/dW = (2⁷ / (3 I)) exp(-4 atan(t) / t) / ((1 + t²)⁴ [1 - exp(-2π/t)]),   t = 1/ξ = sqrt(W/I - 1),
//! ```
//!
//! equal to `(2⁷ / (3 I)) e⁻⁴` at the edge. Integrated from the edge to
//! infinity it gives about 0.435 per electron (computed here, by quadrature
//! of the formula; the rest of the hydrogenic 1s oscillator strength goes to
//! discrete excitations below the edge, which dV2022 eq. (2) leaves out with
//! its step `Θ(E - B)`).
//!
//! **Our choices, not from the sources:**
//!
//! * The ionization energy `I` of the formula is set to the shell's binding
//!   energy `B` from the [`SubshellBindingTable`] (for example EADL2017),
//!   so the edge of the ELF is the tabulated edge. For a hydrogen-like ion
//!   this is the formula itself; for a many-electron atom it amounts to a
//!   hydrogenic shell with the effective charge `Z_eff = sqrt(2B/E_h)` that
//!   reproduces `B`. No screening rule is applied beyond that.
//! * The oscillator strength is the per-electron value times the subshell
//!   occupancy from the same table (Salvat eq. (6.1a), Mikhailov et al. eq.
//!   (6)).
//! * The ELF is sampled on a log-spaced grid from `B` to `B` times
//!   [`ShellElfGrid::max_energy_factor`] and interpolated linearly between
//!   the samples, like every [`OpticalElf`]. The formula is nonrelativistic;
//!   the caller keeps the top of the grid well below `2 B (αZ)⁻²`.
//!
//! # Scope
//!
//! Only the K shell is built: the hydrogenic L- and M-subshell formulas
//! have not been found in an openly readable source, so any other subshell
//! is an error. The valence remainder (total ELF minus the inner shells) is
//! not built here.
//!
//! Nothing here is random: the tables are bit-identical whatever the number
//! of threads.
//!
//! [`ShellResolvedChannels`]: super::ShellResolvedChannels

use super::inner_shell::InnerShell;
use crate::constants::{ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY};
use crate::electron::data::{ElectronDataError, OpticalElf, Subshell, SubshellBindingTable};

type Result<T> = std::result::Result<T, ElectronDataError>;

fn invalid(what: &'static str, reason: String) -> ElectronDataError {
    ElectronDataError::Invalid { what, reason }
}

/// The citation recorded in the provenance of every built ELF.
const SOURCE: &str = "hydrogenic K-shell optical ELF: Stobbe (1930) cross section as printed \
     in Mikhailov, Nefiodov and Plunien, Phys. Lett. A 368, 391 (2007), arXiv:0704.2180, eq. (1); \
     df/dW from Salvat, arXiv:2605.22442 (2026), eqs. (6.1a), (6.2); ELF from de Vera et al., \
     Int. J. Mol. Sci. 23, 6121 (2022), eq. (2), prefactor from Shinotsuka et al., Surf. Interface \
     Anal. 49, 238 (2017), eq. (18); ionization energy set to the tabulated binding energy";

/// The hydrogenic K-shell optical oscillator-strength density per electron,
/// `df/dW` in eV⁻¹, at energy loss `energy_ev` for a shell of binding energy
/// `binding_energy_ev` (module docs, "The formula used"). Zero below the
/// edge; at the edge, its limit `(2⁷/(3B)) e⁻⁴`.
pub fn hydrogenic_k_oscillator_strength_density_per_ev(
    binding_energy_ev: f64,
    energy_ev: f64,
) -> f64 {
    // Also zero for a NaN loss.
    if energy_ev.is_nan() || energy_ev < binding_energy_ev {
        return 0.0;
    }
    let prefactor = 128.0 / (3.0 * binding_energy_ev);
    let eps = energy_ev / binding_energy_ev;
    let t2 = eps - 1.0;
    if t2 <= 0.0 {
        return prefactor * (-4.0f64).exp();
    }
    let t = t2.sqrt();
    // exp(-4 ξ cot⁻¹ξ) with ξ = 1/t: ξ cot⁻¹ ξ = atan(t)/t.
    let num = (-4.0 * t.atan() / t).exp();
    // (1 + ξ⁻²)⁴ = (1 + t²)⁴ = ε⁴; 1 - exp(-2πξ) = -expm1(-2π/t).
    let den = eps.powi(4) * -(-std::f64::consts::TAU / t).exp_m1();
    prefactor * num / den
}

/// The sampling grid of a built shell ELF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellElfGrid {
    /// Log-spaced samples per decade of energy (at least 2).
    pub points_per_decade: usize,
    /// The top of the grid over the binding energy, `W_max / B` (above 1).
    pub max_energy_factor: f64,
}

impl Default for ShellElfGrid {
    /// 300 points per decade up to `10³ B`. The oscillator strength above
    /// `10³ B` is below 1e-6 of the shell's, and the linear interpolation
    /// raises the f-sum by about 1e-4 (both tested; the error falls as the
    /// square of the spacing).
    fn default() -> Self {
        Self {
            points_per_decade: 300,
            max_energy_factor: 1.0e3,
        }
    }
}

impl ShellElfGrid {
    fn energies_ev(&self, binding_energy_ev: f64) -> Result<Vec<f64>> {
        let f = self.max_energy_factor;
        if self.points_per_decade < 2 || !(f.is_finite() && f > 1.0) {
            return Err(invalid(
                "shell ELF grid",
                format!(
                    "need at least 2 points per decade and a finite top factor above 1, got {} \
                     and {f}",
                    self.points_per_decade
                ),
            ));
        }
        let decades = f.log10();
        let n = ((decades * self.points_per_decade as f64).ceil() as usize).max(1);
        let mut e: Vec<f64> = (0..=n)
            .map(|i| binding_energy_ev * 10f64.powf(decades * i as f64 / n as f64))
            .collect();
        // The first sample is the edge itself, bit for bit.
        e[0] = binding_energy_ev;
        Ok(e)
    }
}

/// `(ħΩ_a)² = ħ² n_a e² / (ε₀ m_e)` in eV² for `n_a` per m³ (Shinotsuka et
/// al. 2017, eq. (18)).
fn plasma_energy_sq_ev2(density_per_m3: f64) -> f64 {
    let j2 = HBAR * HBAR * density_per_m3 * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE
        / (VACUUM_PERMITTIVITY * ELECTRON_MASS);
    j2 / (ELEMENTARY_CHARGE * ELEMENTARY_CHARGE)
}

/// The hydrogenic optical ELF of one inner shell (module docs): the shell
/// `subshell` of element `z`, with binding energy and occupancy from
/// `table`, for `atoms_per_m3` atoms of that element per m³ in the target
/// `material`. Returns the [`InnerShell`] and its [`OpticalElf`], ready for
/// [`super::ShellResolvedChannels::new`]: the ELF starts at the binding
/// energy and is zero below it.
///
/// Errors if the subshell is not the K shell (no other hydrogenic formula is
/// implemented), if `(z, subshell)` is not in `table`, if `atoms_per_m3` is
/// not finite and positive, or if the grid is invalid.
pub fn hydrogenic_shell_elf(
    material: &str,
    table: &SubshellBindingTable,
    z: u8,
    subshell: Subshell,
    atoms_per_m3: f64,
    grid: ShellElfGrid,
) -> Result<(InnerShell, OpticalElf)> {
    if subshell != Subshell::K {
        return Err(invalid(
            "hydrogenic shell ELF",
            format!(
                "Z = {z} {}: only the K shell is implemented (no hydrogenic formula for other \
                 subshells from an open source yet)",
                subshell.label()
            ),
        ));
    }
    if !(atoms_per_m3.is_finite() && atoms_per_m3 > 0.0) {
        return Err(invalid(
            "hydrogenic shell ELF",
            format!("atom density must be finite and positive, got {atoms_per_m3} m^-3"),
        ));
    }
    let binding = table
        .atom(z)
        .and_then(|a| a.shell(subshell))
        .ok_or_else(|| {
            invalid(
                "hydrogenic shell ELF",
                format!("Z = {z} {} is not in the binding table", subshell.label()),
            )
        })?;
    let b = binding.binding_energy_ev();
    let occupancy = binding.occupancy();
    let shell = InnerShell {
        z,
        subshell,
        binding_energy_ev: b,
    };
    let energy_ev = grid.energies_ev(b)?;
    let scale = std::f64::consts::FRAC_PI_2 * plasma_energy_sq_ev2(atoms_per_m3) * occupancy;
    let elf = energy_ev
        .iter()
        .map(|&w| scale * hydrogenic_k_oscillator_strength_density_per_ev(b, w) / w)
        .collect();
    let provenance = format!(
        "{SOURCE}; Z = {z} {}, B = {b} eV and occupancy {occupancy} from: {}; {atoms_per_m3:e} \
         atoms/m^3; {} points per decade up to {} B",
        subshell.label(),
        table.provenance(),
        grid.points_per_decade,
        grid.max_energy_factor,
    );
    let elf = OpticalElf::new(material, provenance, energy_ev, elf)?;
    Ok((shell, elf))
}

/// [`hydrogenic_shell_elf`] for several shells of a target `material`, each
/// given as `(Z, subshell, atoms of that element per m³)` and kept in the
/// given order. For a compound with `N` molecules per m³ and `ν_j` atoms of
/// element `j` per molecule, the density of element `j` is `ν_j N` (dV2022
/// eq. (2)).
pub fn hydrogenic_shell_elfs(
    material: &str,
    table: &SubshellBindingTable,
    shells: &[(u8, Subshell, f64)],
    grid: ShellElfGrid,
) -> Result<Vec<(InnerShell, OpticalElf)>> {
    shells
        .iter()
        .map(|&(z, s, n)| hydrogenic_shell_elf(material, table, z, s, n, grid))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn density_is_zero_below_the_edge_and_continuous_above_it() {
        let b = 284.0;
        assert_eq!(
            hydrogenic_k_oscillator_strength_density_per_ev(b, 283.9),
            0.0
        );
        assert_eq!(
            hydrogenic_k_oscillator_strength_density_per_ev(b, f64::NAN),
            0.0
        );
        let at = hydrogenic_k_oscillator_strength_density_per_ev(b, b);
        let just_above = hydrogenic_k_oscillator_strength_density_per_ev(b, b * (1.0 + 1e-10));
        assert!((just_above / at - 1.0).abs() < 1e-4, "{at} {just_above}");
        // Monotonically decreasing above the edge.
        let mut prev = at;
        for k in 1..200 {
            let w = b * (1.0 + 0.05 * k as f64);
            let x = hydrogenic_k_oscillator_strength_density_per_ev(b, w);
            assert!(x < prev && x > 0.0, "{w}");
            prev = x;
        }
    }

    #[test]
    fn grid_starts_at_the_edge() {
        let e = ShellElfGrid::default().energies_ev(1838.9).unwrap();
        assert_eq!(e[0], 1838.9);
        assert_eq!(e.len(), 901);
        assert!((e[900] / 1.8389e6 - 1.0).abs() < 1e-12);
        assert!(e.windows(2).all(|w| w[1] > w[0]));
        for bad in [
            ShellElfGrid {
                points_per_decade: 1,
                max_energy_factor: 10.0,
            },
            ShellElfGrid {
                points_per_decade: 10,
                max_energy_factor: 1.0,
            },
            ShellElfGrid {
                points_per_decade: 10,
                max_energy_factor: f64::INFINITY,
            },
        ] {
            assert!(bad.energies_ev(100.0).is_err());
        }
    }
}
