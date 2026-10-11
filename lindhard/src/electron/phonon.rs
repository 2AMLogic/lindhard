//! Insulator energy-loss channels for low-energy electrons: Fröhlich
//! longitudinal-optical (LO) phonon scattering (emission and absorption) and
//! polaron trapping.
//!
//! Both are extra channels of the event loop in [`crate::electron::transport`],
//! each with its own inverse mean free path next to the elastic and inelastic
//! tables. They are **opt-in per layer** ([`InsulatorChannels`]) and off by
//! default: they describe an electron in a polar insulator, and a layer that
//! is a metal (or any layer the caller does not opt in) has neither. The choice
//! is recorded in the run metadata.
//!
//! # Sources
//!
//! The model is Fröhlich's electron–LO-phonon coupling (H. Fröhlich,
//! "Electrons in lattice fields", Adv. Phys. 3, 325 (1954),
//! doi:10.1080/00018735400101213) as applied to low-energy electron transport
//! by J. Llacer and E. L. Garwin, J. Appl. Phys. 40, 2766 and 2776 (1969)
//! (doi:10.1063/1.1658076 for part II), and the polaron-trapping rate of
//! J. P. Ganachaud and A. Mokrani, Surf. Sci. 334, 329 (1995),
//! doi:10.1016/0039-6028(95)00474-2. None of these three papers was open to
//! us (no open-access copy; the publisher pages answer a bot challenge, see
//! `docs/data-provenance.md`). The equations implemented here are taken from
//! two open-access papers that restate them, and the two agree:
//!
//! - **Ding et al. (2021):** Z. J. Ding, C. Li, B. Da, J. Liu, "Charging effect induced by
//!   electron beam irradiation: a review", Sci. Technol. Adv. Mater. 22, 932
//!   (2021), doi:10.1080/14686996.2021.1976597 (CC BY 4.0), Section 2.3:
//!   Eq. (14) (LO-phonon creation frequency), Eq. (15) (annihilation
//!   frequency), the Bose occupation in the text after Eq. (14), and Eq. (16)
//!   (angular distribution). The review credits these to H.-J. Fitting and
//!   J.-U. Friemann, Phys. Status Solidi A 69, 349 (1982), building on
//!   Llacer and Garwin.
//! - **Taioli & Dapor (2025):** S. Taioli and M. Dapor, "Advancements in Secondary and
//!   Backscattered Electron Energy Spectra and Yields Analysis: from Theory
//!   to Applications", arXiv:2404.07521v5 (2025): Eq. (33) (phonon-emission
//!   inverse mean free path), Eq. (34) (Bose occupation) and Eq. (37)
//!   (polaron inverse mean free path `C exp(-γE)`), both credited there to
//!   Ganachaud and Mokrani (1995).
//!
//! # Fröhlich LO-phonon scattering
//!
//! With `x = ħω/E`, `n` the Bose occupation of the mode at temperature `T`
//! (Ding et al. (2021) text after Eq. (14); Taioli & Dapor Eq. (34)),
//!
//! ```text
//! n(T) = 1 / (exp(ħω / k_B T) - 1)
//! ```
//!
//! the inverse mean free paths for emitting (`+`, needs `E > ħω`) and
//! absorbing (`-`) one phonon of energy `ħω` are
//!
//! ```text
//! λ₊⁻¹(E) = (n + 1) / a₀ · (1/ε∞ - 1/ε₀) · x/2 · ln[(1 + √(1 - x)) / (1 - √(1 - x))]
//! λ₋⁻¹(E) =  n      / a₀ · (1/ε∞ - 1/ε₀) · x/2 · ln[(√(1 + x) + 1) / (√(1 + x) - 1)]
//! ```
//!
//! `λ₊⁻¹` is Taioli & Dapor Eq. (33) (where the coupling is written
//! `(ε₀ - ε∞)/(ε₀ ε∞)`, the same quantity). Both are also Ding et al. (2021) Eqs. (14)
//! and (15), which give the rates per unit time `f±`; dividing by the speed
//! `v = √(2E/m*)` gives the forms above with `a₀ = 4πε₀ħ²/(m e²)` when the
//! effective mass `m*` is the free-electron mass, which is the choice made
//! here (as in Taioli & Dapor Eq. (33)). `ε₀` and `ε∞` are the static and
//! high-frequency relative permittivities of the material, `a₀` the Bohr
//! radius.
//!
//! **Detailed balance.** With `E' = E - ħω` the two logarithms coincide:
//! both equal `ln[(√E + √E')² / ħω]`. Hence
//! `E λ₊⁻¹(E) = ((n + 1)/n) E' λ₋⁻¹(E')`: the rate per unit time times the
//! density of states (`∝ √E`) of a transition `E → E'` and of its reverse
//! `E' → E` stand in the ratio `(n + 1)/n = exp(ħω/k_B T)`. This follows from
//! the two formulas above by algebra; the tests check it, analytically and by
//! Monte Carlo.
//!
//! **Angle.** The polar deflection `θ` in either process follows Ding et al. (2021)
//! Eq. (16), with `E` and `E'` the energies before and after the event
//! (`E' = E ∓ ħω`):
//!
//! ```text
//! dw/dΩ = (EE')^(1/2) / [π (E + E' - 2(EE')^(1/2) cos θ) ln R],
//! R = (E + E' + 2(EE')^(1/2)) / (E + E' - 2(EE')^(1/2))
//! ```
//!
//! It is normalised over the sphere. Its cumulative distribution in
//! `μ = cos θ` from `-1` is `F(μ) = ln[(a + b)/(a - bμ)] / ln R`
//! (`a = E + E'`, `b = 2(EE')^(1/2)`), which inverts to
//! `μ = (a - (a + b) R^(-u)) / b` for `u` uniform on `[0, 1)`. That inversion
//! is our own algebra on Eq. (16). The azimuth is uniform.
//!
//! # Polaron trapping
//!
//! The trapping inverse mean free path is (Taioli & Dapor Eq. (37), from Ganachaud and
//! Mokrani 1995)
//!
//! ```text
//! λ_pol⁻¹(E) = C exp(-γ E)
//! ```
//!
//! with `C` (here in m⁻¹) and `γ` (eV⁻¹) material constants. A trapped
//! electron deposits its whole remaining energy at the trapping point and its
//! history ends.
//!
//! # Material parameters
//!
//! Every parameter is a validated input with a provenance string, recorded in
//! the run metadata. The only committed material values are the SiO₂ ones of
//! [`FrohlichPhonon::sio2`] (Ding et al. (2021) Section 2.3). No `C`, `γ` is committed for
//! any material: the source that tabulates them (Ganachaud and Mokrani 1995)
//! could not be opened, so [`PolaronTrapping::new`] takes them from the
//! caller. The gaps are listed in `docs/data-provenance.md`.

use serde::Serialize;

use crate::constants::{BOHR_RADIUS, BOLTZMANN, ELEMENTARY_CHARGE};

/// Errors from building an insulator channel.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum PhononError {
    /// A parameter is out of range.
    #[error("invalid {what}: {why}")]
    Invalid {
        /// The parameter.
        what: &'static str,
        /// Why it was rejected.
        why: String,
    },
}

fn invalid<T>(what: &'static str, why: impl Into<String>) -> Result<T, PhononError> {
    Err(PhononError::Invalid {
        what,
        why: why.into(),
    })
}

fn check_provenance(p: &str) -> Result<(), PhononError> {
    if p.trim().is_empty() {
        return invalid("provenance", "must name the source of the parameters");
    }
    Ok(())
}

/// Boltzmann constant in eV/K (`k_B / e`, both exact in CODATA 2022, see
/// [`crate::constants`]).
const BOLTZMANN_EV_PER_K: f64 = BOLTZMANN / ELEMENTARY_CHARGE;

/// The two LO-phonon modes of SiO₂ quoted in Ding et al. (2021) Section 2.3 (Ding, Li,
/// Da, Liu, Sci. Technol. Adv. Mater. 22, 932 (2021),
/// doi:10.1080/14686996.2021.1976597), who cite J.-C. Kuhr and H.-J. Fitting,
/// J. Electron Spectrosc. Relat. Phenom. 105, 257 (1999) for them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Sio2LoMode {
    /// `ħω = 63 meV`.
    Mev63,
    /// `ħω = 153 meV`.
    Mev153,
}

/// One Fröhlich LO-phonon mode of a polar material at a temperature: the
/// parameters of the emission and absorption channels. See the module docs
/// for the equations and their sources.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FrohlichPhonon {
    hbar_omega_ev: f64,
    eps_static: f64,
    eps_high_frequency: f64,
    temperature_k: f64,
    provenance: String,
}

impl FrohlichPhonon {
    /// Check and build. `hbar_omega_ev` is the LO-phonon energy `ħω` (eV,
    /// finite and positive), `eps_static` and `eps_high_frequency` the static
    /// and high-frequency relative permittivities `ε₀ > ε∞ > 0` (finite), and
    /// `temperature_k` the lattice temperature (K, finite, `>= 0`; at `0` there
    /// is no absorption). `provenance` names the source of the values and must
    /// not be blank.
    pub fn new(
        hbar_omega_ev: f64,
        eps_static: f64,
        eps_high_frequency: f64,
        temperature_k: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, PhononError> {
        if !(hbar_omega_ev.is_finite() && hbar_omega_ev > 0.0) {
            return invalid("phonon energy", "must be finite and positive");
        }
        if !(eps_high_frequency.is_finite() && eps_high_frequency > 0.0) {
            return invalid("high-frequency permittivity", "must be finite and positive");
        }
        if !(eps_static.is_finite() && eps_static > eps_high_frequency) {
            return invalid(
                "static permittivity",
                format!(
                    "{eps_static} must be finite and exceed the high-frequency value {eps_high_frequency}"
                ),
            );
        }
        if !(temperature_k.is_finite() && temperature_k >= 0.0) {
            return invalid("temperature", "must be finite and non-negative");
        }
        let provenance = provenance.into();
        check_provenance(&provenance)?;
        Ok(Self {
            hbar_omega_ev,
            eps_static,
            eps_high_frequency,
            temperature_k,
            provenance,
        })
    }

    /// One LO mode of SiO₂ at `temperature_k`, with `ε₀ = 3.9` and
    /// `ε∞ = 2.25` and the mode energies of [`Sio2LoMode`], all as quoted in
    /// Ding, Li, Da, Liu, Sci. Technol. Adv. Mater. 22, 932 (2021),
    /// doi:10.1080/14686996.2021.1976597, Section 2.3 (the text after
    /// Eq. (14) for `ε₀`, `ε∞`; the text before it for the modes). The review
    /// treats the two modes one at a time; this model has one mode per layer,
    /// and how to combine both modes is not modelled here.
    pub fn sio2(mode: Sio2LoMode, temperature_k: f64) -> Result<Self, PhononError> {
        let hbar_omega_ev = match mode {
            Sio2LoMode::Mev63 => 0.063,
            Sio2LoMode::Mev153 => 0.153,
        };
        Self::new(
            hbar_omega_ev,
            3.9,
            2.25,
            temperature_k,
            "SiO2: eps_static 3.9, eps_inf 2.25, LO modes 63 and 153 meV as quoted in \
             Ding, Li, Da, Liu, Sci. Technol. Adv. Mater. 22, 932 (2021), \
             doi:10.1080/14686996.2021.1976597, Sec. 2.3 (modes credited there to \
             Kuhr & Fitting, J. Electron Spectrosc. Relat. Phenom. 105, 257 (1999))",
        )
    }

    /// LO-phonon energy `ħω`, eV.
    pub fn hbar_omega_ev(&self) -> f64 {
        self.hbar_omega_ev
    }

    /// Static relative permittivity `ε₀`.
    pub fn eps_static(&self) -> f64 {
        self.eps_static
    }

    /// High-frequency relative permittivity `ε∞`.
    pub fn eps_high_frequency(&self) -> f64 {
        self.eps_high_frequency
    }

    /// Lattice temperature, K.
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Source of the parameter values.
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Bose occupation `n(T) = 1/(exp(ħω/k_B T) - 1)` of the mode (Ding et al. (2021)
    /// text after Eq. (14); Taioli & Dapor Eq. (34)); `0` at `T = 0`.
    pub fn occupation(&self) -> f64 {
        if self.temperature_k == 0.0 {
            return 0.0;
        }
        1.0 / (self.hbar_omega_ev / (BOLTZMANN_EV_PER_K * self.temperature_k)).exp_m1()
    }

    /// `(1/ε∞ - 1/ε₀) / (2 a₀)`, m⁻¹: the prefactor shared by both rates.
    fn coupling_per_m(&self) -> f64 {
        (1.0 / self.eps_high_frequency - 1.0 / self.eps_static) / (2.0 * BOHR_RADIUS)
    }

    /// Inverse mean free path for emitting one phonon at energy `e_ev`, m⁻¹
    /// (Taioli & Dapor Eq. (33); Ding et al. (2021) Eq. (14) divided by the speed). Zero for
    /// `e_ev <= ħω`. The logarithm is evaluated as
    /// `ln[(1 + √(1 - x))² / x]`, the same quantity without the cancellation
    /// in `1 - √(1 - x)` at small `x`.
    pub fn emission_inverse_mfp_per_m(&self, e_ev: f64) -> f64 {
        if !e_ev.is_finite() || e_ev <= self.hbar_omega_ev {
            return 0.0;
        }
        let x = self.hbar_omega_ev / e_ev;
        let r = (1.0 - x).sqrt();
        (self.occupation() + 1.0) * self.coupling_per_m() * x * ((1.0 + r).powi(2) / x).ln()
    }

    /// Inverse mean free path for absorbing one phonon at energy `e_ev`, m⁻¹
    /// (Ding et al. (2021) Eq. (15) divided by the speed). Zero at `T = 0` or for
    /// `e_ev <= 0`. The logarithm is evaluated as `ln[(√(1 + x) + 1)² / x]`.
    pub fn absorption_inverse_mfp_per_m(&self, e_ev: f64) -> f64 {
        let n = self.occupation();
        if n == 0.0 || !e_ev.is_finite() || e_ev <= 0.0 {
            return 0.0;
        }
        let x = self.hbar_omega_ev / e_ev;
        let q = (1.0 + x).sqrt();
        n * self.coupling_per_m() * x * ((q + 1.0).powi(2) / x).ln()
    }
}

/// `cos θ` of a phonon event that takes the electron from `e_ev` to
/// `e_after_ev` (both positive and different), at cumulative probability `u`
/// in `[0, 1)`: the inverse of the distribution of Ding et al. (2021)
/// Eq. (16), see the module docs. `ln R` is evaluated as
/// `2 ln[(√E + √E')² / |E - E'|]` to avoid the cancellation in
/// `E + E' - 2√(EE')`.
pub fn sample_cos_theta(e_ev: f64, e_after_ev: f64, u: f64) -> f64 {
    let (s, s2) = (e_ev.sqrt(), e_after_ev.sqrt());
    let a = e_ev + e_after_ev;
    let b = 2.0 * s * s2;
    let apb = (s + s2) * (s + s2);
    let ln_r = 2.0 * (apb / (e_ev - e_after_ev).abs()).ln();
    ((a - apb * (-u * ln_r).exp()) / b).clamp(-1.0, 1.0)
}

/// The cumulative distribution `F(μ)` of `μ = cos θ` under Ding et al. (2021) Eq. (16)
/// for a phonon event `E -> E'` (the inverse of [`sample_cos_theta`]).
pub fn cos_theta_cdf(e_ev: f64, e_after_ev: f64, mu: f64) -> f64 {
    let (s, s2) = (e_ev.sqrt(), e_after_ev.sqrt());
    let b = 2.0 * s * s2;
    let apb = (s + s2) * (s + s2);
    let ln_r = 2.0 * (apb / (e_ev - e_after_ev).abs()).ln();
    let mu = mu.clamp(-1.0, 1.0);
    // a - b·mu = (√E - √E')² + b(1 - mu), without the cancellation at mu = 1.
    let denom = (s - s2) * (s - s2) + b * (1.0 - mu);
    ((apb / denom).ln() / ln_r).clamp(0.0, 1.0)
}

/// Polaron trapping with inverse mean free path `C exp(-γE)` (Taioli & Dapor Eq. (37),
/// Taioli & Dapor, arXiv:2404.07521v5, from Ganachaud & Mokrani, Surf. Sci.
/// 334, 329 (1995)).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PolaronTrapping {
    c_per_m: f64,
    gamma_per_ev: f64,
    provenance: String,
}

impl PolaronTrapping {
    /// Check and build. `c_per_m` is `C` in m⁻¹ (finite, positive),
    /// `gamma_per_ev` is `γ` in eV⁻¹ (finite, `>= 0`). No material values are
    /// built in (see the module docs): the caller supplies them with a
    /// non-blank `provenance`.
    pub fn new(
        c_per_m: f64,
        gamma_per_ev: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, PhononError> {
        if !(c_per_m.is_finite() && c_per_m > 0.0) {
            return invalid("polaron C", "must be finite and positive");
        }
        if !(gamma_per_ev.is_finite() && gamma_per_ev >= 0.0) {
            return invalid("polaron gamma", "must be finite and non-negative");
        }
        let provenance = provenance.into();
        check_provenance(&provenance)?;
        Ok(Self {
            c_per_m,
            gamma_per_ev,
            provenance,
        })
    }

    /// `C`, m⁻¹.
    pub fn c_per_m(&self) -> f64 {
        self.c_per_m
    }

    /// `γ`, eV⁻¹.
    pub fn gamma_per_ev(&self) -> f64 {
        self.gamma_per_ev
    }

    /// Source of the parameter values.
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Trapping inverse mean free path `C exp(-γE)` at `e_ev`, m⁻¹.
    pub fn inverse_mfp_per_m(&self, e_ev: f64) -> f64 {
        self.c_per_m * (-self.gamma_per_ev * e_ev).exp()
    }
}

/// The insulator channels switched on in one layer. The default has none,
/// which is what a metal (or any layer not opted in) gets.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct InsulatorChannels {
    /// Fröhlich LO-phonon emission and absorption, if on.
    pub phonon: Option<FrohlichPhonon>,
    /// Polaron trapping, if on.
    pub polaron: Option<PolaronTrapping>,
}

impl InsulatorChannels {
    /// No insulator channel (the default, and the only choice for a metal).
    pub fn none() -> Self {
        Self::default()
    }

    /// Whether any channel is on.
    pub fn any(&self) -> bool {
        self.phonon.is_some() || self.polaron.is_some()
    }

    /// The three rates at `e_ev`, m⁻¹: (phonon emission, phonon absorption,
    /// polaron trapping); zero for a channel that is off.
    pub fn rates(&self, e_ev: f64) -> (f64, f64, f64) {
        let (em, ab) = match &self.phonon {
            Some(p) => (
                p.emission_inverse_mfp_per_m(e_ev),
                p.absorption_inverse_mfp_per_m(e_ev),
            ),
            None => (0.0, 0.0),
        };
        let tr = self
            .polaron
            .as_ref()
            .map_or(0.0, |p| p.inverse_mfp_per_m(e_ev));
        (em, ab, tr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mode() -> FrohlichPhonon {
        // Synthetic parameters, not material data.
        FrohlichPhonon::new(0.1, 4.0, 2.0, 600.0, "synthetic").unwrap()
    }

    #[test]
    fn stable_logs_equal_the_textbook_forms() {
        let p = mode();
        let pref = (1.0 / 2.0 - 1.0 / 4.0) / (2.0 * BOHR_RADIUS);
        let n = p.occupation();
        for e in [0.15, 1.0, 30.0] {
            let x: f64 = 0.1 / e;
            let em =
                (n + 1.0) * pref * x * ((1.0 + (1.0 - x).sqrt()) / (1.0 - (1.0 - x).sqrt())).ln();
            let ab = n * pref * x * (((1.0 + x).sqrt() + 1.0) / ((1.0 + x).sqrt() - 1.0)).ln();
            assert!((p.emission_inverse_mfp_per_m(e) / em - 1.0).abs() < 1e-12);
            assert!((p.absorption_inverse_mfp_per_m(e) / ab - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn detailed_balance_holds_exactly_in_the_formulas() {
        let p = mode();
        let n = p.occupation();
        for e in [0.11, 0.5, 2.0, 40.0] {
            let ep = e - 0.1;
            let ratio =
                e * p.emission_inverse_mfp_per_m(e) / (ep * p.absorption_inverse_mfp_per_m(ep));
            assert!(
                (ratio / ((n + 1.0) / n) - 1.0).abs() < 1e-12,
                "{e}: {ratio}"
            );
        }
    }

    #[test]
    fn emission_needs_energy_above_the_phonon_and_absorption_needs_heat() {
        let p = mode();
        assert_eq!(p.emission_inverse_mfp_per_m(0.1), 0.0);
        assert_eq!(p.emission_inverse_mfp_per_m(0.05), 0.0);
        let cold = FrohlichPhonon::new(0.1, 4.0, 2.0, 0.0, "synthetic").unwrap();
        assert_eq!(cold.occupation(), 0.0);
        assert_eq!(cold.absorption_inverse_mfp_per_m(1.0), 0.0);
        assert!(cold.emission_inverse_mfp_per_m(1.0) > 0.0);
    }

    #[test]
    fn angular_inverse_cdf_round_trips_and_is_normalised() {
        for (e, ep) in [(1.0, 0.9), (0.2, 0.1), (5.0, 5.1), (50.0, 49.95)] {
            assert!(sample_cos_theta(e, ep, 0.0) <= -1.0 + 1e-9);
            assert!((cos_theta_cdf(e, ep, 1.0) - 1.0).abs() < 1e-12);
            assert!(cos_theta_cdf(e, ep, -1.0).abs() < 1e-12);
            for u in [0.01, 0.3, 0.5, 0.9, 0.999] {
                let mu = sample_cos_theta(e, ep, u);
                assert!((cos_theta_cdf(e, ep, mu) - u).abs() < 1e-9, "{e} {ep} {u}");
            }
            // Midpoint-rule integral of Eq. (16) over the sphere (skipped for
            // a forward spike narrower than the quadrature step).
            if (e - ep).abs() / e < 0.01 {
                continue;
            }
            let (a, b) = (e + ep, 2.0 * f64::sqrt(e * ep));
            let ln_r = ((a + b) / (a - b)).ln();
            let k = 200_000;
            let mut sum = 0.0;
            for i in 0..k {
                let mu = -1.0 + (i as f64 + 0.5) * 2.0 / k as f64;
                sum += 2.0 * std::f64::consts::PI * f64::sqrt(e * ep)
                    / (std::f64::consts::PI * (a - b * mu) * ln_r)
                    * (2.0 / k as f64);
            }
            assert!((sum - 1.0).abs() < 1e-3, "{e} {ep}: {sum}");
        }
    }

    #[test]
    fn validation_rejects_bad_parameters() {
        assert!(FrohlichPhonon::new(0.0, 4.0, 2.0, 300.0, "s").is_err());
        assert!(FrohlichPhonon::new(0.1, 2.0, 2.0, 300.0, "s").is_err());
        assert!(FrohlichPhonon::new(0.1, 4.0, 0.0, 300.0, "s").is_err());
        assert!(FrohlichPhonon::new(0.1, 4.0, 2.0, -1.0, "s").is_err());
        assert!(FrohlichPhonon::new(0.1, 4.0, 2.0, 300.0, " ").is_err());
        assert!(PolaronTrapping::new(0.0, 1.0, "s").is_err());
        assert!(PolaronTrapping::new(1e9, -1.0, "s").is_err());
        assert!(PolaronTrapping::new(1e9, 1.0, "").is_err());
        assert!(PolaronTrapping::new(1e9, 0.0, "s").is_ok());
    }

    #[test]
    fn sio2_preset_carries_the_quoted_values() {
        let p = FrohlichPhonon::sio2(Sio2LoMode::Mev153, 300.0).unwrap();
        assert_eq!(p.hbar_omega_ev(), 0.153);
        assert_eq!(p.eps_static(), 3.9);
        assert_eq!(p.eps_high_frequency(), 2.25);
        assert!(p.provenance().contains("10.1080/14686996.2021.1976597"));
        let q = FrohlichPhonon::sio2(Sio2LoMode::Mev63, 300.0).unwrap();
        assert_eq!(q.hbar_omega_ev(), 0.063);
    }
}
