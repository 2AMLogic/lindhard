//! Choosing the ELF-extension model per material.
//!
//! [`PennAlgorithm`] names the three models (the single-pole approximation,
//! [`super::penn::SinglePolePenn`], the full Penn algorithm,
//! [`super::full_penn::FullPenn`], and the Mermin-ELF oscillator fit,
//! [`super::mermin::MerminPenn`]); [`PennInelastic`] holds one of them built
//! from one material's optical ELF and gives all the same interface. The
//! choice is made per material (one [`PennInelastic`] per
//! [`OpticalElf`]) and is recorded in the metadata of whatever is derived
//! from it: [`PennInelastic::model_identity`] is the string for the `model`
//! field of a [`crate::electron::data::CrossSectionTable`], and it includes
//! the algorithm's label, so a table cannot be mistaken for one of the other
//! algorithm.

use super::full_penn::FullPenn;
use super::mermin::MerminPenn;
use super::mermin_fit::MerminFitOptions;
use super::penn::{InelasticPoint, SinglePolePenn};
use crate::electron::data::{ElectronDataError, OpticalElf};

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Which ELF-extension model to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PennAlgorithm {
    /// The single-pole approximation (Shinotsuka et al. 2017 eqs. (7)-(13)).
    SinglePole,
    /// The full Penn algorithm over Lindhard dielectric functions (eq. (4)).
    Full,
    /// A fit of Mermin oscillators to the optical ELF (MELF-GOS; dV2022
    /// eqs. (3)-(5), see [`super::mermin`]), extended to finite momentum by
    /// the Mermin dielectric function.
    Mermin,
}

impl PennAlgorithm {
    /// The stable label, `"penn-single-pole"`, `"penn-full"` or
    /// `"mermin-melf"`.
    pub fn label(self) -> &'static str {
        match self {
            Self::SinglePole => "penn-single-pole",
            Self::Full => "penn-full",
            Self::Mermin => "mermin-melf",
        }
    }

    /// The algorithm with this label.
    pub fn from_label(label: &str) -> Option<Self> {
        [Self::SinglePole, Self::Full, Self::Mermin]
            .into_iter()
            .find(|a| a.label() == label)
    }
}

impl std::str::FromStr for PennAlgorithm {
    type Err = ElectronDataError;

    fn from_str(s: &str) -> Result<Self> {
        Self::from_label(s).ok_or_else(|| ElectronDataError::Invalid {
            what: "Penn algorithm",
            reason: format!("unknown label {s:?}; expected \"penn-single-pole\", \"penn-full\" or \"mermin-melf\""),
        })
    }
}

/// One of the three models, built from one material's optical ELF.
#[derive(Debug, Clone)]
pub enum PennInelastic {
    /// The single-pole model.
    SinglePole(SinglePolePenn),
    /// The full model.
    Full(FullPenn),
    /// The Mermin-ELF model (a fit of oscillators).
    Mermin(MerminPenn),
}

impl PennInelastic {
    /// The model of `algorithm` for `elf` (Fermi energy zero, default
    /// tolerance). For [`PennAlgorithm::Mermin`] the oscillators are fitted
    /// with [`MerminFitOptions::default`]; use [`Self::try_new`] to handle a
    /// failed fit, or [`Self::mermin`] to choose the fit options.
    ///
    /// # Panics
    ///
    /// Only for [`PennAlgorithm::Mermin`], if the default fit fails (an ELF
    /// table too short for three oscillators, or an ELF that is zero).
    pub fn new(algorithm: PennAlgorithm, elf: OpticalElf) -> Self {
        Self::try_new(algorithm, elf).expect("the default Mermin fit failed; use try_new")
    }

    /// As [`Self::new`], returning the error of a failed Mermin fit.
    pub fn try_new(algorithm: PennAlgorithm, elf: OpticalElf) -> Result<Self> {
        Ok(match algorithm {
            PennAlgorithm::SinglePole => Self::SinglePole(SinglePolePenn::new(elf)),
            PennAlgorithm::Full => Self::Full(FullPenn::new(elf)),
            PennAlgorithm::Mermin => Self::mermin(elf, &MerminFitOptions::default())?,
        })
    }

    /// The Mermin-ELF model with oscillators fitted to `elf` with `options`.
    pub fn mermin(elf: OpticalElf, options: &MerminFitOptions) -> Result<Self> {
        Ok(Self::Mermin(MerminPenn::fit(elf, options)?))
    }

    /// Which algorithm this is.
    pub fn algorithm(&self) -> PennAlgorithm {
        match self {
            Self::SinglePole(_) => PennAlgorithm::SinglePole,
            Self::Full(_) => PennAlgorithm::Full,
            Self::Mermin(_) => PennAlgorithm::Mermin,
        }
    }

    /// The optical ELF the model is built from.
    pub fn optical_elf(&self) -> &OpticalElf {
        match self {
            Self::SinglePole(m) => m.optical_elf(),
            Self::Full(m) => m.optical_elf(),
            Self::Mermin(m) => m.optical_elf(),
        }
    }

    /// Set the Fermi energy, eV (see the two models).
    pub fn with_fermi_energy_ev(self, fermi_ev: f64) -> Result<Self> {
        Ok(match self {
            Self::SinglePole(m) => Self::SinglePole(m.with_fermi_energy_ev(fermi_ev)?),
            Self::Full(m) => Self::Full(m.with_fermi_energy_ev(fermi_ev)?),
            Self::Mermin(m) => Self::Mermin(m.with_fermi_energy_ev(fermi_ev)?),
        })
    }

    /// The Fermi energy, eV (see the two models).
    pub fn fermi_energy_ev(&self) -> f64 {
        match self {
            Self::SinglePole(m) => m.fermi_energy_ev(),
            Self::Full(m) => m.fermi_energy_ev(),
            Self::Mermin(m) => m.fermi_energy_ev(),
        }
    }

    /// The identity string to record as the `model` of derived tables:
    /// the algorithm label, the nonrelativistic kinematics, the Fermi energy
    /// and the integration tolerance, e.g.
    /// `penn-full; nonrelativistic; E_F = 0 eV; rel tol 1e-4`. For the Mermin
    /// model the fitted oscillators follow the label (their parameters
    /// decide the table), e.g. `mermin-melf (3 oscillators (A=...))`.
    pub fn model_identity(&self) -> String {
        let (tol, fermi) = match self {
            Self::SinglePole(m) => (m.relative_tolerance(), m.fermi_energy_ev()),
            Self::Full(m) => (m.relative_tolerance(), m.fermi_energy_ev()),
            Self::Mermin(m) => (m.relative_tolerance(), m.fermi_energy_ev()),
        };
        let label = match self {
            Self::Mermin(m) => format!(
                "{} ({})",
                self.algorithm().label(),
                m.fit_report().summary()
            ),
            _ => self.algorithm().label().to_string(),
        };
        format!("{label}; nonrelativistic; E_F = {fermi} eV; rel tol {tol:e}")
    }

    /// The DIIMFP, m⁻¹ eV⁻¹ (see the two models).
    pub fn diimfp_per_m_ev(&self, energy_ev: f64, loss_ev: f64) -> Result<f64> {
        match self {
            Self::SinglePole(m) => m.diimfp_per_m_ev(energy_ev, loss_ev),
            Self::Full(m) => m.diimfp_per_m_ev(energy_ev, loss_ev),
            Self::Mermin(m) => m.diimfp_per_m_ev(energy_ev, loss_ev),
        }
    }

    /// The inverse IMFP and the stopping power at `energy_ev`.
    pub fn imfp_and_stopping(&self, energy_ev: f64) -> Result<InelasticPoint> {
        match self {
            Self::SinglePole(m) => m.imfp_and_stopping(energy_ev),
            Self::Full(m) => m.imfp_and_stopping(energy_ev),
            Self::Mermin(m) => m.imfp_and_stopping(energy_ev),
        }
    }

    /// IMFP and stopping power on a grid of energies.
    pub fn tabulate(&self, energies_ev: &[f64]) -> Result<Vec<InelasticPoint>> {
        match self {
            Self::SinglePole(m) => m.tabulate(energies_ev),
            Self::Full(m) => m.tabulate(energies_ev),
            Self::Mermin(m) => m.tabulate(energies_ev),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_round_trip_and_unknown_labels_are_rejected() {
        for a in [
            PennAlgorithm::SinglePole,
            PennAlgorithm::Full,
            PennAlgorithm::Mermin,
        ] {
            assert_eq!(PennAlgorithm::from_label(a.label()), Some(a));
            assert_eq!(a.label().parse::<PennAlgorithm>().unwrap(), a);
        }
        assert!(PennAlgorithm::from_label("penn").is_none());
        assert!("penn".parse::<PennAlgorithm>().is_err());
    }
}
