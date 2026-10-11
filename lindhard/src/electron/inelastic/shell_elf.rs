//! Inner-shell optical ELFs built from the hydrogenic K-shell (1s) and
//! L-subshell (2s, 2p) photoionization formulas, for
//! [`ShellResolvedChannels`].
//!
//! # Sources (opened 2026-10-08; Karzas and Latter opened 2026-10-09)
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
//! * **The hydrogenic L-subshell cross sections.** W. J. Karzas and R.
//!   Latter, "Electron Radiative Transitions in a Coulomb Field", Astrophys.
//!   J. Suppl. 6, 167 (1961) ("KL1961"; open scan at
//!   <https://articles.adsabs.harvard.edu/pdf/1961ApJS....6..167K>, read
//!   from the page images, the text layer of the scan being unreliable in
//!   the equations). For a dipole transition from the bound state `(n, ℓ)`
//!   of a hydrogen-like ion to the continuum, averaged over the magnetic
//!   substates (eq. (27), one electron), the cross section is the sum, eq.
//!   (38), of the `ℓ -> ℓ-1` part, eq. (36),
//!
//!   ```text
//!   σ_(nℓ -> E,ℓ-1) = (π e² / (m c ν)) (2^(4ℓ) / 3)
//!       ℓ² (n+ℓ)! {[1² + η²][2² + η²] ... [(ℓ-1)² + η²]} / ((2ℓ+1)! (2ℓ-1)! (n-ℓ-1)!)
//!       × (exp[-4η cot⁻¹ρ] / (1 - e^(-2πη))) (ρ^(2ℓ+2) / (1 + ρ²)^(2n-2))
//!       × [G_ℓ(ℓ+1-n; η; ρ) - (1 + ρ²)⁻² G_ℓ(ℓ-1-n; η; ρ)]²
//!   ```
//!
//!   (the expression in braces is 1 for `ℓ = 1`), and the `ℓ -> ℓ+1` part,
//!   eq. (37),
//!
//!   ```text
//!   σ_(nℓ -> E,ℓ+1) = (π e² / (m c ν)) (2^(4ℓ+6) / 3)
//!       (ℓ+1)² (n+ℓ)! (1² + η²)(2² + η²) ... [(ℓ+1)² + η²]
//!         / ((2ℓ+1) (2ℓ+1)! (2ℓ+2)! (n-ℓ-1)! [(ℓ+1)² + η²]²)
//!       × (exp[-4η cot⁻¹ρ] / (1 - e^(-2πη))) (ρ^(2ℓ+4) η² / (1 + ρ²)^(2n))
//!       × [(ℓ+1-n) G_(ℓ+1)(ℓ+1-n; η; ρ) + ((ℓ+1+n) / (1 + ρ²)) G_(ℓ+1)(ℓ-n; η; ρ)]²,
//!   ```
//!
//!   with `η = (Z² Ry / E)^(1/2)` for a free electron of energy `E` (eq.
//!   (33)), `ρ = η/n`, and the real polynomials of its Appendix C,
//!   `G_ℓ(-m; η; ρ) = Σ_(s=0)^(2m) b_s ρ^s` (eq. (C.5)) with, eq. (C.8),
//!
//!   ```text
//!   b_0 = 1,   b_1 = 2mη/ℓ,
//!   b_s = -(1 / (s (s + 2ℓ - 1))) [4η (s-1-m) b_(s-1) + (2m+2-s)(2m+2ℓ+1-s) b_(s-2)].
//!   ```
//!
//!   Kramers' cross section is its eq. (39),
//!   `σ^K = (2⁴ / (3√3)) (e² / (m c ν)) (1/n) (ρ² / (1 + ρ²))²`, and the
//!   Gaunt factor its eq. (40), `g = σ / σ^K`, tabulated in its Table 1 for
//!   `n = 1 ... 6` and every `ℓ` against `E / Z² Ry`.
//!
//! # The formulas used
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
//! **The L subshells.** The photon energy is `W = hν = E + I` with
//! `I = Z² Ry / n²`, so `W/I = 1 + 1/ρ²` and `ρ = 1/t` with the same
//! `t = sqrt(W/I - 1)`. Salvat's eq. (6.2) turns the prefactor
//! `π e² / (m c ν)` of KL1961 eqs. (36) and (37) into `1/W`: `df/dW` is the
//! rest of the right-hand side over `W`. For `n = 2` (`η = 2ρ`) the
//! polynomials of eq. (C.8) are (our algebra)
//!
//! ```text
//! G_1(0) = G_2(0) = 1,
//! G_1(-1; η; ρ) = 1 + 2ηρ - ρ²                              = 1 + 3ρ²,
//! G_1(-2; η; ρ) = 1 + 4ηρ + ((8η² - 10)/3) ρ² - 4ηρ³ + ρ⁴ = (1 + ρ²)(1 + (11/3) ρ²),
//! G_2(-1; η; ρ) = 1 + ηρ - ρ²                               = 1 + ρ²,
//! ```
//!
//! the bracket of eq. (37) is `2 (1 + 4ρ²)` for 2s and `4` for 2p, and that
//! of eq. (36) is `-(8/3) ρ² / (1 + ρ²)` for 2p. Per electron, with
//! `ε = W/I = 1 + t²`,
//!
//! ```text
//! df_2s/dW = (2¹⁰ / (3 I)) ((ε + 3) / ε⁵)    exp(-8 atan(t) / t) / [1 - exp(-4π/t)],
//! df_2p/dW = (2¹⁰ / (9 I)) ((3ε + 8) / ε⁶)   exp(-8 atan(t) / t) / [1 - exp(-4π/t)],
//! ```
//!
//! equal to `(2¹² / (3 I)) e⁻⁸` and `(11 · 2¹⁰ / (9 I)) e⁻⁸` at the edge.
//! Of the 2p value, `(2¹³ / (27 I)) (ε + 3) / ε⁶` times the same last factor
//! is the transition to the d continuum, eq. (37), and
//! `(2¹⁰ / (27 I)) / ε⁵` times it the transition to the s continuum, eq.
//! (36). The same
//! reduction for `n = 1` gives back the K-shell formula above, an
//! independent check of the transcription. The unit tests of this module
//! evaluate eqs. (36), (37) and (C.8) as printed and compare them with these
//! closed forms; `lindhard/tests/penn_shell_elf.rs` compares the closed
//! forms with the 2s and 2p Gaunt factors of KL1961 Table 1. Integrated from
//! the edge to infinity they give about 0.351 (2s) and 0.191 (2p) per
//! electron (computed here, by quadrature of the formulas).
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
//!   (6)). L1 is the 2s formula. L2 (2p1/2) and L3 (2p3/2) are both the
//!   nonrelativistic 2p formula, each with its own tabulated binding energy
//!   and occupancy, so the 2p oscillator strength is shared between them in
//!   the ratio of their occupancies.
//! * The ELF is sampled on a log-spaced grid from `B` to `B` times
//!   [`ShellElfGrid::max_energy_factor`] and interpolated linearly between
//!   the samples, like every [`OpticalElf`]. The formula is nonrelativistic;
//!   the caller keeps the top of the grid well below `2 B (αZ)⁻²`.
//!
//! # Scope
//!
//! Only the K, L1, L2 and L3 subshells are built; any other subshell (M and
//! above) is an error. KL1961 eqs. (36) and (37) hold for every `(n, ℓ)`,
//! but their reduction for `n >= 3` has not been done or checked here. The
//! valence remainder (total ELF minus the inner shells) is not built here.
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

/// The citation recorded in the provenance of every built K-shell ELF.
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

/// The citation recorded in the provenance of every built L-subshell ELF.
const SOURCE_L: &str = "hydrogenic L-subshell optical ELF: 2s and 2p cross sections from \
     Karzas and Latter, Astrophys. J. Suppl. 6, 167 (1961), eqs. (36), (37), (38) and (C.8), \
     reduced for n = 2; df/dW from Salvat, arXiv:2605.22442 (2026), eqs. (6.1a), (6.2); ELF from \
     de Vera et al., Int. J. Mol. Sci. 23, 6121 (2022), eq. (2), prefactor from Shinotsuka et al., \
     Surf. Interface Anal. 49, 238 (2017), eq. (18); ionization energy set to the tabulated \
     binding energy";

/// The factor `exp(-8 atan(t)/t) / [1 - exp(-4π/t)]` common to the two
/// `n = 2` formulas (`exp[-4η cot⁻¹ρ] / (1 - e^(-2πη))` of Karzas and Latter
/// 1961, eqs. (36) and (37), with `η = 2ρ = 2/t`), for `t² = W/B - 1`; its
/// limit `e⁻⁸` at the edge.
fn n2_coulomb_factor(t2: f64) -> f64 {
    if t2 <= 0.0 {
        return (-8.0f64).exp();
    }
    let t = t2.sqrt();
    (-8.0 * t.atan() / t).exp() / -(-2.0 * std::f64::consts::TAU / t).exp_m1()
}

/// The hydrogenic 2s (L1) optical oscillator-strength density per electron,
/// `df/dW` in eV⁻¹, at energy loss `energy_ev` for a subshell of binding
/// energy `binding_energy_ev`: Karzas and Latter, Astrophys. J. Suppl. 6,
/// 167 (1961), eq. (37) for `n = 2`, `ℓ = 0`, converted by Salvat,
/// arXiv:2605.22442, eq. (6.2) (module docs, "The formulas used"). Zero
/// below the edge; at the edge, its limit `(2¹²/(3B)) e⁻⁸`.
pub fn hydrogenic_2s_oscillator_strength_density_per_ev(
    binding_energy_ev: f64,
    energy_ev: f64,
) -> f64 {
    // Also zero for a NaN loss.
    if energy_ev.is_nan() || energy_ev < binding_energy_ev {
        return 0.0;
    }
    let eps = energy_ev / binding_energy_ev;
    1024.0 / (3.0 * binding_energy_ev) * (eps + 3.0) / eps.powi(5) * n2_coulomb_factor(eps - 1.0)
}

/// The hydrogenic 2p (L2, L3) optical oscillator-strength density per
/// electron, `df/dW` in eV⁻¹, at energy loss `energy_ev` for a subshell of
/// binding energy `binding_energy_ev`: Karzas and Latter, Astrophys. J.
/// Suppl. 6, 167 (1961), eqs. (36), (37) and (38) for `n = 2`, `ℓ = 1`
/// (the s and d continua together), converted by Salvat, arXiv:2605.22442,
/// eq. (6.2) (module docs, "The formulas used"). Zero below the edge; at
/// the edge, its limit `(11 · 2¹⁰/(9B)) e⁻⁸`.
pub fn hydrogenic_2p_oscillator_strength_density_per_ev(
    binding_energy_ev: f64,
    energy_ev: f64,
) -> f64 {
    // Also zero for a NaN loss.
    if energy_ev.is_nan() || energy_ev < binding_energy_ev {
        return 0.0;
    }
    let eps = energy_ev / binding_energy_ev;
    1024.0 / (9.0 * binding_energy_ev) * (3.0 * eps + 8.0) / eps.powi(6)
        * n2_coulomb_factor(eps - 1.0)
}

/// The hydrogenic orbital a subshell is built from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Orbital {
    /// K.
    S1,
    /// L1.
    S2,
    /// L2 and L3.
    P2,
}

impl Orbital {
    /// K, L1, L2 and L3 (ENDF designators 1 to 4); nothing else.
    fn of(subshell: Subshell) -> Option<Self> {
        match subshell.designator() {
            1 => Some(Self::S1),
            2 => Some(Self::S2),
            3 | 4 => Some(Self::P2),
            _ => None,
        }
    }

    fn source(self) -> &'static str {
        match self {
            Self::S1 => SOURCE,
            Self::S2 | Self::P2 => SOURCE_L,
        }
    }

    /// `df/dW` per electron, eV⁻¹.
    fn oscillator_strength_density_per_ev(self, binding_energy_ev: f64, energy_ev: f64) -> f64 {
        match self {
            Self::S1 => {
                hydrogenic_k_oscillator_strength_density_per_ev(binding_energy_ev, energy_ev)
            }
            Self::S2 => {
                hydrogenic_2s_oscillator_strength_density_per_ev(binding_energy_ev, energy_ev)
            }
            Self::P2 => {
                hydrogenic_2p_oscillator_strength_density_per_ev(binding_energy_ev, energy_ev)
            }
        }
    }
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
    /// `10³ B` is below 1e-6 of the shell's (K, L1, L2 and L3), and the
    /// linear interpolation raises the f-sum by about 1e-4 (both tested; the
    /// error falls as the square of the spacing).
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
/// K is built from the 1s formula, L1 from the 2s formula, and L2 and L3
/// each from the 2p formula with its own binding energy and occupancy.
///
/// Errors if the subshell is not K, L1, L2 or L3 (no hydrogenic formula is
/// implemented for the M shell and above), if `(z, subshell)` is not in
/// `table`, if `atoms_per_m3` is
/// not finite and positive, or if the grid is invalid.
pub fn hydrogenic_shell_elf(
    material: &str,
    table: &SubshellBindingTable,
    z: u8,
    subshell: Subshell,
    atoms_per_m3: f64,
    grid: ShellElfGrid,
) -> Result<(InnerShell, OpticalElf)> {
    let Some(orbital) = Orbital::of(subshell) else {
        return Err(invalid(
            "hydrogenic shell ELF",
            format!(
                "Z = {z} {}: only the K, L1, L2 and L3 subshells are implemented (the \
                 hydrogenic formulas for the M shell and above are not reduced or checked here)",
                subshell.label()
            ),
        ));
    };
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
        .map(|&w| scale * orbital.oscillator_strength_density_per_ev(b, w) / w)
        .collect();
    let provenance = format!(
        "{}; Z = {z} {}, B = {b} eV and occupancy {occupancy} from: {}; {atoms_per_m3:e} \
         atoms/m^3; {} points per decade up to {} B",
        orbital.source(),
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

    /// `G_ℓ(-m; η; ρ)`, Karzas and Latter 1961 eqs. (C.5) and (C.8), by the
    /// recurrence as printed.
    fn kl_g(l: u32, m: u32, eta: f64, rho: f64) -> f64 {
        if m == 0 {
            return 1.0;
        }
        let (l, m) = (f64::from(l), f64::from(m));
        let mut b = vec![1.0, 2.0 * m * eta / l];
        for s in 2..=(2.0 * m) as usize {
            let sf = s as f64;
            let x = 4.0 * eta * (sf - 1.0 - m) * b[s - 1]
                + (2.0 * m + 2.0 - sf) * (2.0 * m + 2.0 * l + 1.0 - sf) * b[s - 2];
            b.push(-x / (sf * (sf + 2.0 * l - 1.0)));
        }
        b.iter().rev().fold(0.0, |acc, &c| acc * rho + c)
    }

    fn factorial(k: u32) -> f64 {
        (1..=k).map(f64::from).product()
    }

    /// Karzas and Latter 1961 eqs. (36) and (37) as printed, without the
    /// prefactor `π e² / (m c ν)`: the `ℓ -> ℓ-1` and `ℓ -> ℓ+1` parts.
    fn kl_eqs_36_37(n: u32, l: u32, eta: f64) -> (f64, f64) {
        let nf = f64::from(n);
        let lf = f64::from(l);
        let rho = eta / nf;
        let r2 = 1.0 + rho * rho;
        // cot⁻¹ρ = atan(1/ρ).
        let coulomb =
            (-4.0 * eta * (1.0 / rho).atan()).exp() / (1.0 - (-std::f64::consts::TAU * eta).exp());
        let product =
            |top: u32| -> f64 { (1..=top).map(|j| f64::from(j * j) + eta * eta).product() };
        let down = if l == 0 {
            0.0
        } else {
            let bracket = kl_g(l, n - l - 1, eta, rho) - kl_g(l, n - l + 1, eta, rho) / (r2 * r2);
            2f64.powi(4 * l as i32) / 3.0 * lf * lf * factorial(n + l) * product(l - 1)
                / (factorial(2 * l + 1) * factorial(2 * l - 1) * factorial(n - l - 1))
                * coulomb
                * rho.powi(2 * l as i32 + 2)
                / r2.powi(2 * n as i32 - 2)
                * bracket
                * bracket
        };
        let bracket = (lf + 1.0 - nf) * kl_g(l + 1, n - l - 1, eta, rho)
            + (lf + 1.0 + nf) / r2 * kl_g(l + 1, n - l, eta, rho);
        let up = 2f64.powi(4 * l as i32 + 6) / 3.0
            * (lf + 1.0).powi(2)
            * factorial(n + l)
            * product(l + 1)
            / ((2.0 * lf + 1.0)
                * factorial(2 * l + 1)
                * factorial(2 * l + 2)
                * factorial(n - l - 1)
                * ((lf + 1.0).powi(2) + eta * eta).powi(2))
            * coulomb
            * rho.powi(2 * l as i32 + 4)
            * eta
            * eta
            / r2.powi(2 * n as i32)
            * bracket
            * bracket;
        (down, up)
    }

    /// The closed forms are eqs. (36) + (37) over `W` (Salvat eq. (6.2)
    /// cancels the prefactor against `1/W`), with `W/B = 1 + 1/ρ²`. For
    /// `n = 1` this is the Stobbe formula, from a different source.
    #[test]
    fn closed_forms_are_karzas_latter_eqs_36_and_37() {
        let b = 99.2;
        type Density = fn(f64, f64) -> f64;
        let cases: [(u32, u32, Density); 3] = [
            (1, 0, hydrogenic_k_oscillator_strength_density_per_ev),
            (2, 0, hydrogenic_2s_oscillator_strength_density_per_ev),
            (2, 1, hydrogenic_2p_oscillator_strength_density_per_ev),
        ];
        for (n, l, closed) in cases {
            for rho in [30.0, 3.0, 1.0, 0.5, 0.1, 0.03] {
                let w = b * (1.0 + 1.0 / (rho * rho));
                let (down, up) = kl_eqs_36_37(n, l, rho * f64::from(n));
                let want = (down + up) / w;
                let got = closed(b, w);
                // The bracket of eq. (36) for 2p cancels to O(ρ²).
                assert!(
                    (got / want - 1.0).abs() < 1e-9,
                    "n = {n}, l = {l}, rho = {rho}: {got:e} vs {want:e}"
                );
            }
        }
        // The two parts of 2p: (2¹⁰/(27B)) / ε⁵ to the s continuum, eq. (36),
        // and (2¹³/(27B)) (ε + 3)/ε⁶ to the d continuum, eq. (37), each times
        // the Coulomb factor.
        let eps: f64 = 5.0;
        let (down, up) = kl_eqs_36_37(2, 1, 2.0 / (eps - 1.0).sqrt());
        let c = n2_coulomb_factor(eps - 1.0);
        let want_s = 1024.0 / (27.0 * b) / eps.powi(5) * c;
        let want_d = 8192.0 / (27.0 * b) * (eps + 3.0) / eps.powi(6) * c;
        assert!((down / (eps * b) / want_s - 1.0).abs() < 1e-12);
        assert!((up / (eps * b) / want_d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn l_densities_are_zero_below_the_edge_and_continuous_above_it() {
        let b = 99.2;
        type Density = fn(f64, f64) -> f64;
        let cases: [(Density, f64); 2] = [
            (
                hydrogenic_2s_oscillator_strength_density_per_ev,
                4096.0 / 3.0,
            ),
            (
                hydrogenic_2p_oscillator_strength_density_per_ev,
                11264.0 / 9.0,
            ),
        ];
        for (f, edge) in cases {
            assert_eq!(f(b, 99.1), 0.0);
            assert_eq!(f(b, f64::NAN), 0.0);
            let at = f(b, b);
            assert!((at / (edge / b * (-8.0f64).exp()) - 1.0).abs() < 1e-14);
            let just_above = f(b, b * (1.0 + 1e-10));
            assert!((just_above / at - 1.0).abs() < 1e-4, "{at} {just_above}");
            let mut prev = at;
            for k in 1..200 {
                let w = b * (1.0 + 0.05 * f64::from(k));
                let x = f(b, w);
                assert!(x < prev && x > 0.0, "{w}");
                prev = x;
            }
        }
    }

    #[test]
    fn orbitals_are_k_and_the_three_l_subshells() {
        let of = |label: &str| Orbital::of(Subshell::from_label(label).unwrap());
        assert_eq!(of("K"), Some(Orbital::S1));
        assert_eq!(of("L1"), Some(Orbital::S2));
        assert_eq!(of("L2"), Some(Orbital::P2));
        assert_eq!(of("L3"), Some(Orbital::P2));
        assert_eq!(of("M1"), None);
        assert_eq!(of("M5"), None);
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
