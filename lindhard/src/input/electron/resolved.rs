//! The resolved electron run: [`ResolvedElectron`] and its parts.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::path::PathBuf;

use serde::Serialize;

use crate::electron::boundary::BandStructure;
use crate::electron::data::OpticalElf;
use crate::electron::elastic::corrections::{
    Corrections, CorrelationPolarization, PolarizationCutoff, CORRELATION_POLARIZATION_MODEL,
    EXCHANGE_MODEL,
};
use crate::electron::inelastic::{MerminFitOptions, PennAlgorithm};
use crate::electron::phonon::InsulatorChannels;
use crate::electron::secondary::SecondaryModel;
use crate::electron::transport::{BoundaryModel, Primary, TransportConfig};
use crate::geometry::Stack;
use crate::input::{ModelInfo, ResolvedLayer};
use crate::material::Material;
use crate::tally::electron::SE_BSE_SPLIT_SOURCE;
use crate::tally::{ElectronTallyConfig, PsfModel, PsfNormalization};

use super::schema::*;

/// A data file a run read, identified for provenance.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DataFile {
    /// The path as written in the input.
    pub path: String,
    /// Where it was read from (absolute where possible).
    pub resolved_path: PathBuf,
    /// SHA-256 of the file's bytes, lowercase hex.
    pub sha256: String,
}

/// The elastic choices after validation.
#[derive(Debug, Clone, PartialEq)]
pub struct ElasticChoice {
    /// The elastic model.
    pub model: ElasticModelChoice,
    /// The atomic potential.
    pub potential: PotentialChoice,
    /// Furness-McCarthy exchange on or off.
    pub exchange: bool,
    /// Per-element correlation-polarization inputs (empty when off).
    pub correlation_polarization: BTreeMap<u8, CorrelationPolarization>,
}

impl ElasticChoice {
    /// The corrections for element `z`.
    pub fn corrections(&self, z: u8) -> Corrections {
        Corrections {
            exchange: self.exchange,
            correlation_polarization: self.correlation_polarization.get(&z).cloned(),
        }
    }

    /// Whether any correction is on.
    pub fn any_correction(&self) -> bool {
        self.exchange || !self.correlation_polarization.is_empty()
    }

    /// One line naming the corrections and their inputs, appended to the
    /// potential description of every element's table (empty when none is
    /// on).
    pub fn corrections_description(&self) -> String {
        if !self.any_correction() {
            return String::new();
        }
        let mut parts = Vec::new();
        if self.exchange {
            parts.push(format!("exchange: {EXCHANGE_MODEL}"));
        }
        if !self.correlation_polarization.is_empty() {
            let each: Vec<String> = self
                .correlation_polarization
                .iter()
                .map(|(z, c)| {
                    let sym = crate::elements::element(*z).map_or("?", |e| e.symbol);
                    let cut = match c.cutoff {
                        PolarizationCutoff::Seltzer => "Seltzer b_pol^2".to_string(),
                        PolarizationCutoff::BPolSquared(b) => format!("b_pol^2 = {b}"),
                    };
                    format!(
                        "{sym}: alpha_p = {} bohr^3 ({}), {cut}, outer radius {} bohr",
                        c.polarizability, c.polarizability_source, c.outer_radius
                    )
                })
                .collect();
            parts.push(format!(
                "correlation-polarization: {CORRELATION_POLARIZATION_MODEL}; {}",
                each.join("; ")
            ));
        }
        format!("; corrections: {}", parts.join("; "))
    }
}

/// The electron data of one material after loading and validation.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedElectronMaterial {
    /// The name the target uses.
    pub name: String,
    /// The material (composition and density).
    pub material: Material,
    /// The optical ELF, validated (with its provenance).
    pub optical_elf: OpticalElf,
    /// The ELF file.
    pub optical_elf_file: DataFile,
    /// Band parameters, if given.
    pub band: Option<BandStructure>,
    /// Insulator channels (none by default).
    pub channels: InsulatorChannels,
    /// Options of the Mermin oscillator fit of the `mermin-melf` table:
    /// [`MerminFitOptions::default`] unless `mermin_fit` is given.
    pub mermin_fit: MerminFitOptions,
}

/// A validated electron input, in the engine's types.
#[derive(Debug, Clone)]
pub struct ResolvedElectron {
    /// The input as echoed: defaults filled in, `run.threads` removed.
    pub input: ElectronInput,
    /// The target.
    pub stack: Stack,
    /// Layers in stack order, with where each material came from.
    pub layers: Vec<ResolvedLayer>,
    /// For each layer, the index of its material in `materials`.
    pub layer_material: Vec<usize>,
    /// The materials the target uses, by name, with their electron data.
    pub materials: Vec<ResolvedElectronMaterial>,
    /// Transport configuration.
    pub config: TransportConfig,
    /// The primary electron.
    pub primary: Primary,
    /// Elastic choices.
    pub elastic: ElasticChoice,
    /// Inelastic model.
    pub inelastic: PennAlgorithm,
    /// Inelastic Fermi energy, eV, of the tables of materials without a band
    /// (0 if any material has one).
    pub inelastic_fermi_ev: f64,
    /// Energy grid of every table, eV.
    pub table_energy_ev: Vec<f64>,
    /// Tally configuration (SI units).
    pub tally: ElectronTallyConfig,
    /// The PSF forms to fit to the profile of `tally.psf`, in input order
    /// (duplicates removed).
    pub psf_fits: Vec<PsfModel>,
    /// How the PSF fits set their energy scale.
    pub psf_normalization: PsfNormalization,
    /// Non-fatal advice.
    pub warnings: Vec<String>,
}

impl ResolvedElectron {
    /// Every model in use, with its published source, for the output
    /// metadata. Data provenance (the ELF files, band, phonon and polaron
    /// parameters, and the tables built from them) is recorded separately.
    pub fn models(&self) -> Vec<ModelInfo> {
        let mut m = vec![
            ModelInfo {
                role: "electron transport",
                name: "event-by-event",
                citation: Cow::Borrowed(
                    "Kieft and Bosch, J. Phys. D: Appl. Phys. 41, 215310 (2008), \
                     doi:10.1088/0022-3727/41/21/215310",
                ),
            },
            ModelInfo {
                role: "elastic scattering",
                name: "mott-partial-waves",
                citation: Cow::Borrowed(
                    "Mott cross section from radial-Dirac phase shifts (electron::elastic), \
                     independent-atom additivity (electron::elastic::table)",
                ),
            },
        ];
        m.push(match self.elastic.potential {
            PotentialChoice::ThomasFermiYukawa => ModelInfo {
                role: "elastic potential",
                name: "thomas-fermi-yukawa",
                citation: Cow::Borrowed(
                    "STAND-IN Yukawa potential with the Thomas-Fermi screening length \
                     a = 0.8853 a0 Z^(-1/3) (Firsov 1958; Lindhard, Scharff and Schiott 1963); \
                     not a DHFS atomic potential",
                ),
            },
            PotentialChoice::SalvatDhfs => ModelInfo {
                role: "elastic potential",
                name: "salvat-dhfs",
                citation: Cow::Borrowed(
                    "Salvat, Martinez, Mayol and Parellada, Phys. Rev. A 36, 467 (1987), \
                     doi:10.1103/PhysRevA.36.467",
                ),
            },
        });
        if self.elastic.exchange {
            m.push(ModelInfo {
                role: "elastic exchange correction",
                name: "furness-mccarthy",
                citation: Cow::Borrowed(EXCHANGE_MODEL),
            });
        }
        if !self.elastic.correlation_polarization.is_empty() {
            m.push(ModelInfo {
                role: "elastic correlation-polarization correction",
                name: "salvat-2003",
                citation: Cow::Borrowed(CORRELATION_POLARIZATION_MODEL),
            });
        }
        m.push(ModelInfo {
            role: "inelastic scattering",
            name: self.inelastic.label(),
            citation: Cow::Borrowed(match self.inelastic {
                PennAlgorithm::SinglePole => {
                    "single-pole Penn algorithm, Penn, Phys. Rev. B 35, 482 (1987), as written \
                     in Shinotsuka et al., Surf. Interface Anal. 49, 238 (2017), \
                     doi:10.1002/sia.6123, eqs. (7)-(13)"
                }
                PennAlgorithm::Full => {
                    "full Penn algorithm, Penn, Phys. Rev. B 35, 482 (1987), as written in \
                     Shinotsuka et al., Surf. Interface Anal. 49, 238 (2017), \
                     doi:10.1002/sia.6123, eqs. (4)-(6)"
                }
                PennAlgorithm::Mermin => {
                    "Mermin-ELF (MELF-GOS) oscillator fit, Mermin, Phys. Rev. B 1, 2362 (1970), \
                     as written in de Vera et al., Int. J. Mol. Sci. 23, 6121 (2022), \
                     doi:10.3390/ijms23116121, eqs. (3)-(5)"
                }
            }),
        });
        if self.config.secondaries != SecondaryModel::Off {
            m.push(ModelInfo {
                role: "secondary electrons",
                name: "kieft-bosch",
                citation: Cow::Borrowed(
                    "Kieft and Bosch, J. Phys. D 41, 215310 (2008), as written out in Verduin, \
                     PhD thesis, TU Delft (2017), doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c, \
                     Eqs. 3.86, 3.100-3.111; ported from Nebula (BSD-3-Clause)",
                ),
            });
        }
        if self.config.boundary != BoundaryModel::Transparent {
            m.push(ModelInfo {
                role: "surface and interface barrier",
                name: "step-barrier",
                citation: Cow::Borrowed(
                    "inner-potential step, Verduin, PhD thesis, TU Delft (2017), Eqs. 3.136, \
                     3.139, 3.145; ported from Nebula (BSD-3-Clause)",
                ),
            });
        }
        if self.materials.iter().any(|x| x.channels.phonon.is_some()) {
            m.push(ModelInfo {
                role: "LO-phonon scattering",
                name: "frohlich",
                citation: Cow::Borrowed(
                    "Frohlich, Adv. Phys. 3, 325 (1954), as written in Ding, Li, Da, Liu, Sci. \
                     Technol. Adv. Mater. 22, 932 (2021), Eqs. (14)-(16), and Taioli and Dapor, \
                     arXiv:2404.07521v5, Eqs. (33)-(34)",
                ),
            });
        }
        if self.materials.iter().any(|x| x.channels.polaron.is_some()) {
            m.push(ModelInfo {
                role: "polaron trapping",
                name: "ganachaud-mokrani",
                citation: Cow::Borrowed(
                    "C exp(-gamma E), Ganachaud and Mokrani, Surf. Sci. 334, 329 (1995), as \
                     written in Taioli and Dapor, arXiv:2404.07521v5, Eq. (37)",
                ),
            });
        }
        m.push(ModelInfo {
            role: "SE/BSE split",
            name: "energy-split",
            citation: Cow::Borrowed(SE_BSE_SPLIT_SOURCE),
        });
        m
    }
}
