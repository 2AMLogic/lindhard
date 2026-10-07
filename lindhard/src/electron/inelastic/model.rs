//! Choosing the Penn algorithm per material.
//!
//! [`PennAlgorithm`] names the two algorithms (the single-pole approximation,
//! [`super::penn::SinglePolePenn`], and the full algorithm,
//! [`super::full_penn::FullPenn`]); [`PennInelastic`] holds one of them built
//! from one material's optical ELF and gives both the same interface. The
//! choice is made per material (one [`PennInelastic`] per
//! [`OpticalElf`]) and is recorded in the metadata of whatever is derived
//! from it: [`PennInelastic::model_identity`] is the string for the `model`
//! field of a [`crate::electron::data::CrossSectionTable`], and it includes
//! the algorithm's label, so a table cannot be mistaken for one of the other
//! algorithm.

use super::full_penn::FullPenn;
use super::penn::{InelasticPoint, SinglePolePenn};
use crate::electron::data::{ElectronDataError, OpticalElf};

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Which Penn algorithm to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PennAlgorithm {
    /// The single-pole approximation (Shinotsuka et al. 2017 eqs. (7)-(13)).
    SinglePole,
    /// The full Penn algorithm over Lindhard dielectric functions (eq. (4)).
    Full,
}

impl PennAlgorithm {
    /// The stable label, `"penn-single-pole"` or `"penn-full"`.
    pub fn label(self) -> &'static str {
        match self {
            Self::SinglePole => "penn-single-pole",
            Self::Full => "penn-full",
        }
    }

    /// The algorithm with this label.
    pub fn from_label(label: &str) -> Option<Self> {
        [Self::SinglePole, Self::Full]
            .into_iter()
            .find(|a| a.label() == label)
    }
}

impl std::str::FromStr for PennAlgorithm {
    type Err = ElectronDataError;

    fn from_str(s: &str) -> Result<Self> {
        Self::from_label(s).ok_or_else(|| ElectronDataError::Invalid {
            what: "Penn algorithm",
            reason: format!("unknown label {s:?}; expected \"penn-single-pole\" or \"penn-full\""),
        })
    }
}

/// One of the two Penn models, built from one material's optical ELF.
#[derive(Debug, Clone)]
pub enum PennInelastic {
    /// The single-pole model.
    SinglePole(SinglePolePenn),
    /// The full model.
    Full(FullPenn),
}

impl PennInelastic {
    /// The model of `algorithm` for `elf` (Fermi energy zero, default
    /// tolerance).
    pub fn new(algorithm: PennAlgorithm, elf: OpticalElf) -> Self {
        match algorithm {
            PennAlgorithm::SinglePole => Self::SinglePole(SinglePolePenn::new(elf)),
            PennAlgorithm::Full => Self::Full(FullPenn::new(elf)),
        }
    }

    /// Which algorithm this is.
    pub fn algorithm(&self) -> PennAlgorithm {
        match self {
            Self::SinglePole(_) => PennAlgorithm::SinglePole,
            Self::Full(_) => PennAlgorithm::Full,
        }
    }

    /// The optical ELF the model is built from.
    pub fn optical_elf(&self) -> &OpticalElf {
        match self {
            Self::SinglePole(m) => m.optical_elf(),
            Self::Full(m) => m.optical_elf(),
        }
    }

    /// Set the Fermi energy, eV (see the two models).
    pub fn with_fermi_energy_ev(self, fermi_ev: f64) -> Result<Self> {
        Ok(match self {
            Self::SinglePole(m) => Self::SinglePole(m.with_fermi_energy_ev(fermi_ev)?),
            Self::Full(m) => Self::Full(m.with_fermi_energy_ev(fermi_ev)?),
        })
    }

    /// The identity string to record as the `model` of derived tables:
    /// the algorithm label, the nonrelativistic kinematics, the Fermi energy
    /// and the integration tolerance, e.g.
    /// `penn-full; nonrelativistic; E_F = 0 eV; rel tol 1e-4`.
    pub fn model_identity(&self) -> String {
        let (tol, fermi) = match self {
            Self::SinglePole(m) => (m.relative_tolerance(), m.fermi_energy_ev()),
            Self::Full(m) => (m.relative_tolerance(), m.fermi_energy_ev()),
        };
        format!(
            "{}; nonrelativistic; E_F = {fermi} eV; rel tol {tol:e}",
            self.algorithm().label()
        )
    }

    /// The DIIMFP, m⁻¹ eV⁻¹ (see the two models).
    pub fn diimfp_per_m_ev(&self, energy_ev: f64, loss_ev: f64) -> Result<f64> {
        match self {
            Self::SinglePole(m) => m.diimfp_per_m_ev(energy_ev, loss_ev),
            Self::Full(m) => m.diimfp_per_m_ev(energy_ev, loss_ev),
        }
    }

    /// The inverse IMFP and the stopping power at `energy_ev`.
    pub fn imfp_and_stopping(&self, energy_ev: f64) -> Result<InelasticPoint> {
        match self {
            Self::SinglePole(m) => m.imfp_and_stopping(energy_ev),
            Self::Full(m) => m.imfp_and_stopping(energy_ev),
        }
    }

    /// IMFP and stopping power on a grid of energies.
    pub fn tabulate(&self, energies_ev: &[f64]) -> Result<Vec<InelasticPoint>> {
        match self {
            Self::SinglePole(m) => m.tabulate(energies_ev),
            Self::Full(m) => m.tabulate(energies_ev),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_round_trip_and_unknown_labels_are_rejected() {
        for a in [PennAlgorithm::SinglePole, PennAlgorithm::Full] {
            assert_eq!(PennAlgorithm::from_label(a.label()), Some(a));
            assert_eq!(a.label().parse::<PennAlgorithm>().unwrap(), a);
        }
        assert!(PennAlgorithm::from_label("penn").is_none());
        assert!("penn".parse::<PennAlgorithm>().is_err());
    }
}
