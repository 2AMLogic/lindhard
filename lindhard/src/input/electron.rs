//! The electron run description: one TOML document whose `[electron]` table
//! describes the beam, the transport choices, the elastic and inelastic
//! models, the per-material electron data and the tallies, next to the same
//! `[materials]` and `[target]` tables as an ion run.
//!
//! This module only **exposes** what the library already does
//! ([`crate::electron`], [`crate::tally::electron`]); it adds no physics. The
//! schema is documented, with the rules for extending it, in `docs/cli.md`
//! (section "Electron runs"). In short:
//!
//! ```toml
//! [electron.beam]
//! energy_ev = 10000.0
//!
//! [electron.transport]
//! cutoff_ev = 50.0
//!
//! [electron.elastic]
//! potential = "thomas-fermi-yukawa"
//!
//! [electron.materials.Si]
//! optical_elf = "si_elf.toml"
//!
//! [target]
//! substrate = "Si"
//!
//! [run]
//! histories = 1000
//! seed = 1
//! ```
//!
//! [`ElectronInput::resolve_in`] validates the description, reads the data
//! files it names and checks their provenance with the loaders of
//! [`crate::electron::data`] (an optical ELF without a provenance is refused;
//! so are band, phonon and polaron parameters without one), and returns the
//! engine's types ([`ResolvedElectron`]). Building the cross-section tables
//! and running the transport is the caller's (`lindhard-cli`).
//!
//! Like the ion input, every table rejects unknown keys and the echoed input
//! ([`ResolvedElectron::input`]) has every default filled in, so a result can be
//! reproduced from its own header.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use super::{invalid, resolve_material, InputError, ModelInfo, ResolvedLayer, TargetSpec};
use crate::electron::boundary::{BandModel, BandStructure};
use crate::electron::data::OpticalElf;
use crate::electron::elastic::corrections::{
    Corrections, CorrelationPolarization, PolarizationCutoff, CORRELATION_POLARIZATION_MODEL,
    EXCHANGE_MODEL,
};
use crate::electron::elastic::table::{
    DEFAULT_MAX_ENERGY_EV, DEFAULT_MIN_ENERGY_EV, DEFAULT_POINTS_PER_DECADE,
};
use crate::electron::inelastic::PennAlgorithm;
use crate::electron::phonon::{FrohlichPhonon, InsulatorChannels, PolaronTrapping, Sio2LoMode};
use crate::electron::secondary::SecondaryModel;
use crate::electron::transport::{
    BoundaryModel, CutoffReference, EscapeRule, Primary, TransportConfig,
};
use crate::elements::element_by_symbol;
use crate::geometry::Stack;
use crate::material::{Material, MaterialSpec};
use crate::tally::electron::SE_BSE_SPLIT_SOURCE;
use crate::tally::{
    Binning, CartesianGrid, CylindricalGrid, ElectronTallyConfig, LogRadialBinning, PsfConfig,
    PsfModel, PsfNormalization, SE_BSE_SPLIT_EV,
};

const NM: f64 = 1e-9;

fn finite_pos(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

/// The whole electron run description.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronInput {
    /// Everything electron-specific.
    pub electron: ElectronSpec,
    /// Named materials (as in an ion run).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, MaterialSpec>,
    /// The layered target (as in an ion run).
    pub target: TargetSpec,
    /// History count, seed and threads.
    pub run: ElectronRunSpec,
}

/// `[electron]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronSpec {
    /// The primary electrons.
    pub beam: ElectronBeamSpec,
    /// Cutoff, faces, secondaries.
    pub transport: TransportSpec,
    /// Elastic model.
    pub elastic: ElasticSpec,
    /// Inelastic model.
    #[serde(default)]
    pub inelastic: InelasticSpec,
    /// Energy grid of the cross-section tables.
    #[serde(default)]
    pub tables: TableGridSpec,
    /// Electron data of each material the target uses, by the name the
    /// target gives it (a `[materials]` key or an element symbol).
    pub materials: BTreeMap<String, ElectronMaterialSpec>,
    /// What the run records.
    #[serde(default)]
    pub tally: ElectronTallySpec,
}

/// `[electron.beam]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronBeamSpec {
    /// Kinetic energy of the primaries, eV: the vacuum energy with
    /// `boundary = "step-barrier"`, the energy inside the first layer
    /// otherwise.
    pub energy_ev: f64,
    /// Polar angle from the surface normal, degrees, in `[0, 90)`.
    #[serde(default)]
    pub tilt_deg: f64,
    /// Azimuth of the incidence plane, degrees.
    #[serde(default)]
    pub azimuth_deg: f64,
}

/// Where the cutoff is measured from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum CutoffReferenceChoice {
    /// [`CutoffReference::BandBottom`].
    #[default]
    BandBottom,
    /// [`CutoffReference::VacuumLevel`].
    VacuumLevel,
}

/// Escape rule at the faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum EscapeRuleChoice {
    /// [`EscapeRule::BothFaces`].
    #[default]
    BothFaces,
    /// [`EscapeRule::FrontOnly`].
    FrontOnly,
}

/// Secondary-electron model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SecondaryChoice {
    /// [`SecondaryModel::Off`].
    #[default]
    Off,
    /// [`SecondaryModel::KieftBosch`].
    KieftBosch,
}

/// Face model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryChoice {
    /// [`BoundaryModel::Transparent`].
    #[default]
    Transparent,
    /// [`BoundaryModel::StepBarrier`].
    StepBarrier,
}

fn default_max_events() -> u64 {
    TransportConfig::new(1.0).max_events
}

/// `[electron.transport]`: the [`TransportConfig`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TransportSpec {
    /// Stopping cutoff, eV.
    pub cutoff_ev: f64,
    /// Where the cutoff is measured from.
    #[serde(default)]
    pub cutoff_reference: CutoffReferenceChoice,
    /// What the faces of the target do.
    #[serde(default)]
    pub escape_rule: EscapeRuleChoice,
    /// Collision and reflection cap per electron.
    #[serde(default = "default_max_events")]
    pub max_events: u64,
    /// Secondary-electron generation.
    #[serde(default)]
    pub secondaries: SecondaryChoice,
    /// Kieft-Bosch option: instantaneous momentum of the bound electron.
    /// Default `true` with `"kieft-bosch"`; not allowed with `"off"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub instantaneous_momentum: Option<bool>,
    /// Kieft-Bosch option: deflect the primary to conserve momentum.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub momentum_conservation: Option<bool>,
    /// Face model.
    #[serde(default)]
    pub boundary: BoundaryChoice,
    /// Step-barrier option: quantum-mechanical transmission. Default `true`
    /// with `"step-barrier"`; not allowed with `"transparent"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub quantum_transmission: Option<bool>,
    /// Step-barrier option: refraction on transmission.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub refraction: Option<bool>,
}

/// Elastic model. Only the Mott cross section from the radial-Dirac partial
/// waves exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ElasticModelChoice {
    /// [`crate::electron::elastic`] and [`crate::electron::elastic::table`].
    #[default]
    Mott,
}

impl ElasticModelChoice {
    /// The stable label (the input spelling).
    pub fn label(self) -> &'static str {
        match self {
            Self::Mott => "mott",
        }
    }
}

/// Atomic potential of the elastic model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum PotentialChoice {
    /// The Thomas-Fermi-length Yukawa **stand-in**
    /// ([`crate::electron::elastic::table::ThomasFermiYukawa`]).
    ThomasFermiYukawa,
    /// Salvat et al. (1987) DHFS
    /// ([`crate::electron::elastic::table::SalvatDhfsTable`]), Table I
    /// coefficients for Z = 1..92.
    SalvatDhfs,
}

impl PotentialChoice {
    /// The stable label.
    pub fn label(self) -> &'static str {
        match self {
            Self::ThomasFermiYukawa => "thomas-fermi-yukawa",
            Self::SalvatDhfs => "salvat-dhfs",
        }
    }
}

/// `[electron.elastic]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElasticSpec {
    /// The elastic model.
    #[serde(default)]
    pub model: ElasticModelChoice,
    /// The atomic potential (required: the only one that runs is a
    /// stand-in, and the input must say so).
    pub potential: PotentialChoice,
    /// Furness-McCarthy exchange correction.
    #[serde(default)]
    pub exchange: bool,
    /// Correlation-polarization correction, absent for off.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub correlation_polarization: Option<CorrelationPolarizationSpec>,
}

/// One element's dipole polarizability and its source.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolarizabilitySpec {
    /// Static dipole polarizability, bohr³.
    pub bohr3: f64,
    /// Citation of the value (required, non-blank).
    pub source: String,
}

fn default_outer_radius() -> f64 {
    crate::electron::elastic::corrections::DEFAULT_OUTER_RADIUS
}

/// `[electron.elastic.correlation_polarization]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CorrelationPolarizationSpec {
    /// `b_pol²` of the polarization cutoff; absent: Seltzer's
    /// `(E - 50 eV)/(16 eV)`, which needs every table energy above 50 eV.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub b_pol_squared: Option<f64>,
    /// Radius beyond which the polarization tail is dropped, bohr.
    #[serde(default = "default_outer_radius")]
    pub outer_radius_bohr: f64,
    /// Polarizability of every element of the target, by symbol.
    pub polarizability: BTreeMap<String, PolarizabilitySpec>,
}

/// `[electron.inelastic]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InelasticSpec {
    /// `"penn-single-pole"`, `"penn-full"` or `"mermin-melf"`.
    #[serde(default = "default_inelastic_model")]
    pub model: String,
    /// Fermi energy of the model, eV, for the inelastic tables of materials
    /// without a band. Must be 0 if a material has a band: the table of a
    /// material with a band takes the band's (#241).
    #[serde(default)]
    pub fermi_energy_ev: f64,
}

fn default_inelastic_model() -> String {
    PennAlgorithm::SinglePole.label().to_string()
}

impl Default for InelasticSpec {
    fn default() -> Self {
        Self {
            model: default_inelastic_model(),
            fermi_energy_ev: 0.0,
        }
    }
}

fn default_min_energy() -> f64 {
    DEFAULT_MIN_ENERGY_EV
}

fn default_points_per_decade() -> f64 {
    DEFAULT_POINTS_PER_DECADE
}

/// `[electron.tables]`: the energy grid both tables of every layer share.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TableGridSpec {
    /// Lowest grid energy, eV.
    #[serde(default = "default_min_energy")]
    pub min_energy_ev: f64,
    /// Highest grid energy, eV. Default: the beam energy (plus the largest
    /// inner potential with the step barrier).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_energy_ev: Option<f64>,
    /// Minimum points per decade.
    #[serde(default = "default_points_per_decade")]
    pub points_per_decade: f64,
}

impl Default for TableGridSpec {
    fn default() -> Self {
        Self {
            min_energy_ev: default_min_energy(),
            max_energy_ev: None,
            points_per_decade: default_points_per_decade(),
        }
    }
}

/// Band parameters of a material ([`BandStructure`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum BandSpec {
    /// [`BandModel::Metal`].
    Metal {
        /// Fermi energy above the band bottom, eV.
        fermi_ev: f64,
        /// Work function, eV.
        work_function_ev: f64,
        /// Source of the values.
        provenance: String,
    },
    /// [`BandModel::Insulator`].
    Insulator {
        /// Valence band width, eV.
        valence_band_width_ev: f64,
        /// Band gap, eV.
        band_gap_ev: f64,
        /// Electron affinity, eV.
        affinity_ev: f64,
        /// Source of the values.
        provenance: String,
    },
    /// [`BandStructure::free_electron_metal`]: the Fermi energy from the
    /// free-electron relation and the material's atom density.
    FreeElectronMetal {
        /// Valence electrons per atom (a model choice).
        valence_electrons_per_atom: f64,
        /// Work function, eV.
        work_function_ev: f64,
        /// Source of the work function and of the valence count.
        provenance: String,
    },
}

/// A Fröhlich LO-phonon channel: explicit parameters, or a preset the library
/// carries with its citation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhononSpec {
    /// `"sio2-63mev"` or `"sio2-153mev"` ([`FrohlichPhonon::sio2`]); then
    /// only `temperature_k` may be given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub preset: Option<String>,
    /// LO-phonon energy, eV.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hbar_omega_ev: Option<f64>,
    /// Static relative permittivity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_static: Option<f64>,
    /// High-frequency relative permittivity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_high_frequency: Option<f64>,
    /// Lattice temperature, K.
    pub temperature_k: f64,
    /// Source of the values (required without a preset).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<String>,
}

/// A polaron-trapping channel ([`PolaronTrapping`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolaronSpec {
    /// `C`, nm⁻¹.
    pub c_per_nm: f64,
    /// `γ`, eV⁻¹.
    pub gamma_per_ev: f64,
    /// Source of the values.
    pub provenance: String,
}

/// `[electron.materials.<name>]`: the electron data of one material.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronMaterialSpec {
    /// Path of the optical ELF file ([`OpticalElf`] TOML form), relative to
    /// the input file's directory.
    pub optical_elf: String,
    /// Band parameters; required with secondaries, the step barrier or the
    /// vacuum-level cutoff.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub band: Option<BandSpec>,
    /// Fröhlich LO-phonon channel (polar insulators only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phonon: Option<PhononSpec>,
    /// Polaron trapping (polar insulators only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub polaron: Option<PolaronSpec>,
}

/// A uniform binning in nm.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BinsNm {
    /// Lower edge, nm.
    pub lo_nm: f64,
    /// Upper edge, nm.
    pub hi_nm: f64,
    /// Number of bins.
    pub bins: usize,
}

/// `[electron.tally.cartesian]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CartesianSpec {
    /// Depth bins.
    pub x: BinsNm,
    /// Lateral `y` bins.
    pub y: BinsNm,
    /// Lateral `z` bins.
    pub z: BinsNm,
}

/// `[electron.tally.cylindrical]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CylindricalSpec {
    /// Radial bins (from `r >= 0`).
    pub r: BinsNm,
    /// Depth bins.
    pub depth: BinsNm,
}

/// A PSF form to fit (`fits` of `[electron.tally.psf]`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum PsfFitSpec {
    /// Double Gaussian ([`PsfModel::DoubleGaussian`]).
    Double,
    /// Triple Gaussian ([`PsfModel::TripleGaussian`]).
    Triple,
}

impl From<PsfFitSpec> for PsfModel {
    fn from(f: PsfFitSpec) -> Self {
        match f {
            PsfFitSpec::Double => PsfModel::DoubleGaussian,
            PsfFitSpec::Triple => PsfModel::TripleGaussian,
        }
    }
}

/// How the energy scale of a PSF fit is set (`normalization`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum PsfNormalizationSpec {
    /// Fixed to the slab energy ([`PsfNormalization::SlabTotal`]).
    #[default]
    SlabTotal,
    /// A free parameter ([`PsfNormalization::Free`]).
    Free,
}

impl From<PsfNormalizationSpec> for PsfNormalization {
    fn from(n: PsfNormalizationSpec) -> Self {
        match n {
            PsfNormalizationSpec::SlabTotal => PsfNormalization::SlabTotal,
            PsfNormalizationSpec::Free => PsfNormalization::Free,
        }
    }
}

fn default_psf_fits() -> Vec<PsfFitSpec> {
    vec![PsfFitSpec::Double, PsfFitSpec::Triple]
}

/// `[electron.tally.psf]`: the radial profile of the energy deposited in a
/// depth slab, on log radial bins, and the PSF fits to it
/// ([`crate::tally::psf`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PsfSpec {
    /// Lower depth of the slab, nm.
    pub depth_lo_nm: f64,
    /// Upper depth of the slab (excluded), nm.
    pub depth_hi_nm: f64,
    /// Outer edge of the central disc and inner edge of the first log bin, nm.
    pub r_min_nm: f64,
    /// Outer edge of the last bin, nm.
    pub r_max_nm: f64,
    /// Log bins between `r_min_nm` and `r_max_nm`.
    pub bins: usize,
    /// Forms to fit (`"double"`, `"triple"`); may be empty.
    #[serde(default = "default_psf_fits")]
    pub fits: Vec<PsfFitSpec>,
    /// How the energy scale of the fits is set.
    #[serde(default)]
    pub normalization: PsfNormalizationSpec,
}

fn default_split() -> f64 {
    SE_BSE_SPLIT_EV
}
fn default_energy_bins() -> usize {
    100
}
fn default_polar_max() -> f64 {
    90.0
}
fn default_polar_bins() -> usize {
    18
}

/// `[electron.tally]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronTallySpec {
    /// SE/BSE split, eV.
    #[serde(default = "default_split")]
    pub se_bse_split_ev: f64,
    /// Upper edge of the escape-energy spectra, eV (from 0). Default: the
    /// beam energy.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub escape_energy_max_ev: Option<f64>,
    /// Escape-energy bins.
    #[serde(default = "default_energy_bins")]
    pub escape_energy_bins: usize,
    /// Upper edge of the escape polar-angle spectra, degrees from the
    /// outward normal (from 0, at most 180).
    #[serde(default = "default_polar_max")]
    pub escape_polar_max_deg: f64,
    /// Escape polar-angle bins.
    #[serde(default = "default_polar_bins")]
    pub escape_polar_bins: usize,
    /// Cartesian deposition grid.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cartesian: Option<CartesianSpec>,
    /// Cylindrical deposition grid.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cylindrical: Option<CylindricalSpec>,
    /// Radial profile of a depth slab and its PSF fits.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub psf: Option<PsfSpec>,
}

impl Default for ElectronTallySpec {
    fn default() -> Self {
        Self {
            se_bse_split_ev: default_split(),
            escape_energy_max_ev: None,
            escape_energy_bins: default_energy_bins(),
            escape_polar_max_deg: default_polar_max(),
            escape_polar_bins: default_polar_bins(),
            cartesian: None,
            cylindrical: None,
            psf: None,
        }
    }
}

/// `[run]` of an electron run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronRunSpec {
    /// Number of primary electrons (histories).
    pub histories: u64,
    /// Run seed.
    pub seed: u64,
    /// Worker threads (not echoed; never changes results).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub threads: Option<usize>,
}

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

fn read_data_file(
    base_dir: &Path,
    field: &str,
    path: &str,
) -> Result<(String, DataFile), InputError> {
    let joined = base_dir.join(path);
    let bytes = std::fs::read(&joined)
        .map_err(|e| invalid(field, format!("cannot read {}: {e}", joined.display())))?;
    let text = String::from_utf8(bytes.clone()).map_err(|_| {
        invalid(
            field,
            format!("{} is not valid UTF-8 text", joined.display()),
        )
    })?;
    let sha256 = Sha256::digest(&bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    Ok((
        text,
        DataFile {
            path: path.to_string(),
            resolved_path: std::fs::canonicalize(&joined).unwrap_or(joined),
            sha256,
        },
    ))
}

fn binning_nm(field: &str, b: &BinsNm) -> Result<Binning, InputError> {
    Binning::new(b.lo_nm * NM, b.hi_nm * NM, b.bins)
        .map_err(|e| invalid(field, format!("{e} (edges in nm)")))
}

/// The smallest stopping threshold over the layers (the transport's rule:
/// the cutoff, or the inner potential plus the cutoff from the vacuum level).
fn lowest_threshold(
    materials: &[ResolvedElectronMaterial],
    layer_material: &[usize],
    reference: CutoffReference,
    cutoff_ev: f64,
) -> f64 {
    layer_material
        .iter()
        .map(|&i| match reference {
            CutoffReference::BandBottom => cutoff_ev,
            CutoffReference::VacuumLevel => {
                materials[i]
                    .band
                    .as_ref()
                    .map_or(0.0, BandStructure::inner_potential_ev)
                    + cutoff_ev
            }
        })
        .fold(f64::INFINITY, f64::min)
}

impl ElectronInput {
    /// Whether a TOML document is an electron input: it has a top-level
    /// `electron` table.
    pub fn is_electron_toml(text: &str) -> bool {
        text.parse::<toml::Table>()
            .map(|t| t.contains_key("electron"))
            .unwrap_or(false)
    }

    /// Parse a TOML document (schema errors name the key and line).
    pub fn from_toml_str(text: &str) -> Result<Self, InputError> {
        if text
            .parse::<toml::Table>()
            .is_ok_and(|t| t.contains_key("crystal"))
        {
            return Err(invalid(
                "crystal",
                "not supported with an [electron] run (crystal transport is for ion runs)",
            ));
        }
        toml::from_str(text).map_err(|e| InputError::Parse(e.to_string()))
    }

    /// The input as a TOML document in the same schema.
    pub fn to_toml_string(&self) -> Result<String, InputError> {
        toml::to_string(self).map_err(|e| InputError::Parse(e.to_string()))
    }

    /// Validate everything with relative data paths taken from the current
    /// directory.
    pub fn resolve(&self) -> Result<ResolvedElectron, InputError> {
        self.resolve_in(Path::new("."))
    }

    /// Validate everything, read the data files (paths relative to
    /// `base_dir`, the input file's directory) and build the engine's
    /// inputs. Errors name the offending key.
    pub fn resolve_in(&self, base_dir: &Path) -> Result<ResolvedElectron, InputError> {
        let mut warnings = Vec::new();
        let e = &self.electron;
        let mut echo = self.clone();
        echo.run.threads = None;

        // Run.
        if self.run.histories == 0 {
            return Err(invalid("run.histories", "must be at least 1"));
        }
        if self.run.threads == Some(0) {
            return Err(invalid(
                "run.threads",
                "must be at least 1 (omit for all cores)",
            ));
        }

        // Beam.
        let b = &e.beam;
        if !finite_pos(b.energy_ev) {
            return Err(invalid(
                "electron.beam.energy_ev",
                "must be finite and positive",
            ));
        }
        if !(b.tilt_deg.is_finite() && (0.0..90.0).contains(&b.tilt_deg)) {
            return Err(invalid(
                "electron.beam.tilt_deg",
                "must be in [0, 90) degrees",
            ));
        }
        if !b.azimuth_deg.is_finite() {
            return Err(invalid("electron.beam.azimuth_deg", "must be finite"));
        }
        let (s, c) = b.tilt_deg.to_radians().sin_cos();
        let (sa, ca) = b.azimuth_deg.to_radians().sin_cos();
        let primary = Primary {
            energy_ev: b.energy_ev,
            direction: [c, s * ca, s * sa],
        };

        // Transport configuration.
        let t = &e.transport;
        if !finite_pos(t.cutoff_ev) {
            return Err(invalid(
                "electron.transport.cutoff_ev",
                "must be finite and positive",
            ));
        }
        if t.max_events == 0 {
            return Err(invalid(
                "electron.transport.max_events",
                "must be at least 1",
            ));
        }
        let secondaries = match t.secondaries {
            SecondaryChoice::Off => {
                for (k, v) in [
                    ("instantaneous_momentum", t.instantaneous_momentum),
                    ("momentum_conservation", t.momentum_conservation),
                ] {
                    if v.is_some() {
                        return Err(invalid(
                            format!("electron.transport.{k}"),
                            "only applies with secondaries = \"kieft-bosch\"",
                        ));
                    }
                }
                SecondaryModel::Off
            }
            SecondaryChoice::KieftBosch => {
                let im = t.instantaneous_momentum.unwrap_or(true);
                let mc = t.momentum_conservation.unwrap_or(true);
                echo.electron.transport.instantaneous_momentum = Some(im);
                echo.electron.transport.momentum_conservation = Some(mc);
                SecondaryModel::KieftBosch {
                    instantaneous_momentum: im,
                    momentum_conservation: mc,
                }
            }
        };
        let boundary = match t.boundary {
            BoundaryChoice::Transparent => {
                for (k, v) in [
                    ("quantum_transmission", t.quantum_transmission),
                    ("refraction", t.refraction),
                ] {
                    if v.is_some() {
                        return Err(invalid(
                            format!("electron.transport.{k}"),
                            "only applies with boundary = \"step-barrier\"",
                        ));
                    }
                }
                BoundaryModel::Transparent
            }
            BoundaryChoice::StepBarrier => {
                let qt = t.quantum_transmission.unwrap_or(true);
                let rf = t.refraction.unwrap_or(true);
                echo.electron.transport.quantum_transmission = Some(qt);
                echo.electron.transport.refraction = Some(rf);
                BoundaryModel::StepBarrier {
                    quantum_transmission: qt,
                    refraction: rf,
                }
            }
        };
        let cutoff_reference = match t.cutoff_reference {
            CutoffReferenceChoice::BandBottom => CutoffReference::BandBottom,
            CutoffReferenceChoice::VacuumLevel => CutoffReference::VacuumLevel,
        };
        let config = TransportConfig {
            cutoff_ev: t.cutoff_ev,
            escape_rule: match t.escape_rule {
                EscapeRuleChoice::BothFaces => EscapeRule::BothFaces,
                EscapeRuleChoice::FrontOnly => EscapeRule::FrontOnly,
            },
            max_events: t.max_events,
            secondaries,
            boundary,
            cutoff_reference,
        };
        let needs_bands = secondaries != SecondaryModel::Off
            || boundary != BoundaryModel::Transparent
            || cutoff_reference == CutoffReference::VacuumLevel;

        // Materials (all of them, so a broken unused entry is still reported).
        for (name, spec) in &self.materials {
            Material::try_from(spec.clone())
                .map_err(|err| invalid(format!("materials.{name}"), err.to_string()))?;
        }

        // Target.
        let tg = &self.target;
        if tg.layers.is_empty() && tg.substrate.is_none() {
            return Err(invalid(
                "target",
                "needs at least one [[target.layers]] entry or a substrate",
            ));
        }
        let mut layers = Vec::new();
        let mut fields = Vec::new();
        let mut thick = Vec::new();
        for (i, l) in tg.layers.iter().enumerate() {
            let field = format!("target.layers[{i}]");
            if !finite_pos(l.thickness_nm) {
                return Err(invalid(
                    format!("{field}.thickness_nm"),
                    format!("{} nm must be finite and positive", l.thickness_nm),
                ));
            }
            let f = format!("{field}.material");
            layers.push(resolve_material(&self.materials, &f, &l.material)?);
            fields.push(f);
            thick.push(l.thickness_nm * NM);
        }
        if let Some(sub) = &tg.substrate {
            layers.push(resolve_material(&self.materials, "target.substrate", sub)?);
            fields.push("target.substrate".to_string());
        }
        let substrate = tg.substrate.as_ref().map(|_| {
            layers
                .last()
                .expect("substrate pushed last")
                .material
                .clone()
        });
        let finite: Vec<(Material, f64)> = layers
            .iter()
            .zip(&thick)
            .map(|(l, &t)| (l.material.clone(), t))
            .collect();
        let stack =
            Stack::new(finite, substrate).map_err(|err| invalid("target", err.to_string()))?;

        // Electron data, one entry per material name the target uses.
        let mut materials: Vec<ResolvedElectronMaterial> = Vec::new();
        let mut layer_material = Vec::with_capacity(layers.len());
        for (l, field) in layers.iter().zip(&fields) {
            if l.source == "inline" {
                return Err(invalid(
                    field,
                    "an electron run needs a named material (a [materials] key or an \
                     element symbol), so that [electron.materials.<name>] can give its \
                     electron data",
                ));
            }
            if let Some(i) = materials.iter().position(|m| m.name == l.source) {
                layer_material.push(i);
                continue;
            }
            let name = &l.source;
            let mfield = format!("electron.materials.{name}");
            let spec = e.materials.get(name).ok_or_else(|| {
                invalid(
                    field,
                    format!(
                        "material {name:?} has no electron data; add [{mfield}] with at \
                         least optical_elf"
                    ),
                )
            })?;
            materials.push(self.resolve_electron_material(
                base_dir,
                &mfield,
                name,
                &l.material,
                spec,
                needs_bands,
            )?);
            layer_material.push(materials.len() - 1);
        }
        for name in e.materials.keys() {
            if !materials.iter().any(|m| &m.name == name) {
                warnings.push(format!(
                    "electron.materials.{name}: not used by any target layer"
                ));
            }
        }

        // Elastic.
        let el = &e.elastic;
        let mut cp = BTreeMap::new();
        if let Some(c) = &el.correlation_polarization {
            let f = "electron.elastic.correlation_polarization";
            let cutoff = match c.b_pol_squared {
                None => PolarizationCutoff::Seltzer,
                Some(b) if b.is_finite() && b >= 0.0 => PolarizationCutoff::BPolSquared(b),
                Some(b) => {
                    return Err(invalid(
                        format!("{f}.b_pol_squared"),
                        format!("{b} must be finite and non-negative"),
                    ))
                }
            };
            if !finite_pos(c.outer_radius_bohr) {
                return Err(invalid(
                    format!("{f}.outer_radius_bohr"),
                    "must be finite and positive",
                ));
            }
            let mut zs: Vec<u8> = Vec::new();
            for l in &layers {
                for comp in l.material.components() {
                    if comp.atom_fraction() > 0.0 && !zs.contains(&comp.z()) {
                        zs.push(comp.z());
                    }
                }
            }
            for (sym, p) in &c.polarizability {
                let pf = format!("{f}.polarizability.{sym}");
                let elm = element_by_symbol(sym)
                    .ok_or_else(|| invalid(&pf, format!("unknown element symbol {sym:?}")))?;
                if !zs.contains(&elm.z) {
                    return Err(invalid(&pf, format!("{sym} is not in any target layer")));
                }
                if !finite_pos(p.bohr3) {
                    return Err(invalid(
                        format!("{pf}.bohr3"),
                        "must be finite and positive",
                    ));
                }
                if p.source.trim().is_empty() {
                    return Err(invalid(
                        format!("{pf}.source"),
                        "a polarizability needs a source (a citation of the value)",
                    ));
                }
                cp.insert(
                    elm.z,
                    CorrelationPolarization {
                        polarizability: p.bohr3,
                        polarizability_source: p.source.trim().to_string(),
                        cutoff,
                        outer_radius: c.outer_radius_bohr,
                    },
                );
            }
            for z in zs {
                if !cp.contains_key(&z) {
                    let sym = crate::elements::element(z).expect("valid").symbol;
                    return Err(invalid(
                        format!("{f}.polarizability"),
                        format!("no polarizability for {sym}; give every target element one"),
                    ));
                }
            }
        }
        let elastic = ElasticChoice {
            model: el.model,
            potential: el.potential,
            exchange: el.exchange,
            correlation_polarization: cp,
        };

        // Inelastic.
        let inelastic: PennAlgorithm = e.inelastic.model.parse().map_err(
            |err: crate::electron::data::ElectronDataError| {
                invalid("electron.inelastic.model", err.to_string())
            },
        )?;
        if !(e.inelastic.fermi_energy_ev.is_finite() && e.inelastic.fermi_energy_ev >= 0.0) {
            return Err(invalid(
                "electron.inelastic.fermi_energy_ev",
                "must be finite and non-negative",
            ));
        }
        // The inelastic table of a material with a band takes its Fermi
        // energy from the band (#241); a second one would count it twice.
        if e.inelastic.fermi_energy_ev != 0.0 {
            if let Some(m) = materials.iter().find(|m| m.band.is_some()) {
                return Err(invalid(
                    "electron.inelastic.fermi_energy_ev",
                    format!(
                        "must be 0 when a material has a band (material {} does): the \
                         inelastic table of a material with a band takes its Fermi energy \
                         from the band",
                        m.name
                    ),
                ));
            }
        }

        // Table grid.
        let g = &e.tables;
        let max_inner = materials
            .iter()
            .filter_map(|m| m.band.as_ref().map(BandStructure::inner_potential_ev))
            .fold(0.0, f64::max);
        let top = match boundary {
            BoundaryModel::StepBarrier { .. } => b.energy_ev + max_inner,
            BoundaryModel::Transparent => b.energy_ev,
        };
        let max_energy = g.max_energy_ev.unwrap_or(top);
        echo.electron.tables.max_energy_ev = Some(max_energy);
        if !(finite_pos(g.min_energy_ev) && max_energy.is_finite() && max_energy > g.min_energy_ev)
        {
            return Err(invalid(
                "electron.tables",
                format!(
                    "need 0 < min_energy_ev < max_energy_ev, got {} and {max_energy} eV",
                    g.min_energy_ev
                ),
            ));
        }
        let table_energy_ev = crate::electron::elastic::table::log_energy_grid(
            g.min_energy_ev,
            max_energy,
            g.points_per_decade,
        )
        .map_err(|err| invalid("electron.tables", err.to_string()))?;
        if max_energy < top {
            warnings.push(format!(
                "electron.tables.max_energy_ev: {max_energy} eV is below the largest energy \
                 an electron can have ({top} eV); rates above it are held at the last row \
                 (results.table_coverage counts how often)"
            ));
        }
        // The lowest energy an electron is followed at: the cutoff, or
        // `U + cutoff` from the vacuum level (the transport's thresholds).
        let lowest = lowest_threshold(&materials, &layer_material, cutoff_reference, t.cutoff_ev);
        if g.min_energy_ev > lowest {
            warnings.push(format!(
                "electron.tables.min_energy_ev: {} eV is above the lowest stopping threshold \
                 ({lowest} eV); rates below it are held at the first row \
                 (results.table_coverage counts how often)",
                g.min_energy_ev
            ));
        }
        if max_energy > DEFAULT_MAX_ENERGY_EV {
            warnings.push(format!(
                "electron.tables.max_energy_ev: {max_energy} eV is above the {DEFAULT_MAX_ENERGY_EV} eV \
                 the elastic tables are tested to; the inelastic models are nonrelativistic"
            ));
        }

        // Tally.
        let ts = &e.tally;
        let f = "electron.tally";
        if !(ts.se_bse_split_ev.is_finite() && ts.se_bse_split_ev >= 0.0) {
            return Err(invalid(
                format!("{f}.se_bse_split_ev"),
                "must be finite and non-negative",
            ));
        }
        let emax = ts.escape_energy_max_ev.unwrap_or(b.energy_ev);
        echo.electron.tally.escape_energy_max_ev = Some(emax);
        let escape_energy = Binning::new(0.0, emax, ts.escape_energy_bins).map_err(|err| {
            invalid(
                format!("{f}.escape_energy_max_ev / escape_energy_bins"),
                err.to_string(),
            )
        })?;
        if !(finite_pos(ts.escape_polar_max_deg) && ts.escape_polar_max_deg <= 180.0) {
            return Err(invalid(
                format!("{f}.escape_polar_max_deg"),
                "must be in (0, 180] degrees",
            ));
        }
        let escape_polar = Binning::new(
            0.0,
            ts.escape_polar_max_deg
                .to_radians()
                .min(std::f64::consts::PI),
            ts.escape_polar_bins,
        )
        .map_err(|err| invalid(format!("{f}.escape_polar_bins"), err.to_string()))?;
        let mut tally = ElectronTallyConfig::new(escape_energy, escape_polar);
        tally.se_bse_split_ev = ts.se_bse_split_ev;
        if let Some(c) = &ts.cartesian {
            tally.cartesian = Some(CartesianGrid {
                x: binning_nm(&format!("{f}.cartesian.x"), &c.x)?,
                y: binning_nm(&format!("{f}.cartesian.y"), &c.y)?,
                z: binning_nm(&format!("{f}.cartesian.z"), &c.z)?,
            });
        }
        if let Some(c) = &ts.cylindrical {
            if c.r.lo_nm < 0.0 {
                return Err(invalid(format!("{f}.cylindrical.r.lo_nm"), "must be >= 0"));
            }
            tally.cylindrical = Some(CylindricalGrid {
                r: binning_nm(&format!("{f}.cylindrical.r"), &c.r)?,
                depth: binning_nm(&format!("{f}.cylindrical.depth"), &c.depth)?,
            });
        }
        let mut psf_fits: Vec<PsfModel> = Vec::new();
        let mut psf_normalization = PsfNormalization::default();
        if let Some(p) = &ts.psf {
            let nm = |key: &str, v: f64| -> Result<f64, InputError> {
                if v.is_finite() {
                    Ok(v * NM)
                } else {
                    Err(invalid(format!("{f}.psf.{key}"), "must be finite"))
                }
            };
            let depth_lo = nm("depth_lo_nm", p.depth_lo_nm)?;
            let depth_hi = nm("depth_hi_nm", p.depth_hi_nm)?;
            if depth_hi <= depth_lo {
                return Err(invalid(
                    format!("{f}.psf.depth_hi_nm"),
                    "must be greater than depth_lo_nm",
                ));
            }
            let r_min = nm("r_min_nm", p.r_min_nm)?;
            let r_max = nm("r_max_nm", p.r_max_nm)?;
            if r_min <= 0.0 {
                return Err(invalid(format!("{f}.psf.r_min_nm"), "must be positive"));
            }
            if r_max <= r_min {
                return Err(invalid(
                    format!("{f}.psf.r_max_nm"),
                    "must be greater than r_min_nm",
                ));
            }
            if p.bins == 0 {
                return Err(invalid(format!("{f}.psf.bins"), "must be at least 1"));
            }
            let radial = LogRadialBinning::new(r_min, r_max, p.bins)
                .map_err(|e| invalid(format!("{f}.psf.bins"), e.to_string()))?;
            tally.psf = Some(
                PsfConfig::new(radial, depth_lo, depth_hi)
                    .map_err(|e| invalid(format!("{f}.psf"), e.to_string()))?,
            );
            for m in p.fits.iter().copied().map(PsfModel::from) {
                if !psf_fits.contains(&m) {
                    psf_fits.push(m);
                }
            }
            psf_normalization = p.normalization.into();
        }

        Ok(ResolvedElectron {
            input: echo,
            stack,
            layers,
            layer_material,
            materials,
            config,
            primary,
            elastic,
            inelastic,
            inelastic_fermi_ev: e.inelastic.fermi_energy_ev,
            table_energy_ev,
            tally,
            psf_fits,
            psf_normalization,
            warnings,
        })
    }

    fn resolve_electron_material(
        &self,
        base_dir: &Path,
        field: &str,
        name: &str,
        material: &Material,
        spec: &ElectronMaterialSpec,
        needs_bands: bool,
    ) -> Result<ResolvedElectronMaterial, InputError> {
        let ef = format!("{field}.optical_elf");
        let (text, file) = read_data_file(base_dir, &ef, &spec.optical_elf)?;
        // The loader refuses a missing or blank provenance.
        let optical_elf = OpticalElf::from_toml_str(&text)
            .map_err(|err| invalid(&ef, format!("{}: {err}", spec.optical_elf)))?;

        let band = match &spec.band {
            None if needs_bands => {
                return Err(invalid(
                    format!("{field}.band"),
                    "required with secondaries, the step barrier or the vacuum-level \
                     cutoff (they need the material's Fermi energy and inner potential)",
                ))
            }
            None => None,
            Some(bs) => {
                let r = match bs {
                    BandSpec::Metal {
                        fermi_ev,
                        work_function_ev,
                        provenance,
                    } => BandStructure::new(
                        BandModel::Metal {
                            fermi_ev: *fermi_ev,
                            work_function_ev: *work_function_ev,
                        },
                        provenance.trim(),
                    ),
                    BandSpec::Insulator {
                        valence_band_width_ev,
                        band_gap_ev,
                        affinity_ev,
                        provenance,
                    } => BandStructure::new(
                        BandModel::Insulator {
                            valence_band_width_ev: *valence_band_width_ev,
                            band_gap_ev: *band_gap_ev,
                            affinity_ev: *affinity_ev,
                        },
                        provenance.trim(),
                    ),
                    BandSpec::FreeElectronMetal {
                        valence_electrons_per_atom,
                        work_function_ev,
                        provenance,
                    } => BandStructure::free_electron_metal(
                        material,
                        *valence_electrons_per_atom,
                        *work_function_ev,
                        provenance.trim(),
                    ),
                };
                Some(r.map_err(|err| invalid(format!("{field}.band"), err.to_string()))?)
            }
        };

        let mut channels = InsulatorChannels::none();
        if let Some(p) = &spec.phonon {
            let pf = format!("{field}.phonon");
            let ph = match p.preset.as_deref() {
                Some(preset) => {
                    if p.hbar_omega_ev.is_some()
                        || p.eps_static.is_some()
                        || p.eps_high_frequency.is_some()
                        || p.provenance.is_some()
                    {
                        return Err(invalid(
                            &pf,
                            "a preset fixes the parameters and their source; give only \
                             preset and temperature_k",
                        ));
                    }
                    let mode = match preset {
                        "sio2-63mev" => Sio2LoMode::Mev63,
                        "sio2-153mev" => Sio2LoMode::Mev153,
                        other => {
                            return Err(invalid(
                                format!("{pf}.preset"),
                                format!(
                                    "unknown preset {other:?}; expected \"sio2-63mev\" or \
                                     \"sio2-153mev\""
                                ),
                            ))
                        }
                    };
                    FrohlichPhonon::sio2(mode, p.temperature_k)
                }
                None => {
                    let need = |k: &str, v: Option<f64>| {
                        v.ok_or_else(|| invalid(format!("{pf}.{k}"), "required without a preset"))
                    };
                    FrohlichPhonon::new(
                        need("hbar_omega_ev", p.hbar_omega_ev)?,
                        need("eps_static", p.eps_static)?,
                        need("eps_high_frequency", p.eps_high_frequency)?,
                        p.temperature_k,
                        p.provenance.as_deref().unwrap_or("").trim(),
                    )
                }
            };
            channels.phonon = Some(ph.map_err(|err| invalid(&pf, err.to_string()))?);
        }
        if let Some(p) = &spec.polaron {
            channels.polaron = Some(
                PolaronTrapping::new(p.c_per_nm / NM, p.gamma_per_ev, p.provenance.trim())
                    .map_err(|err| invalid(format!("{field}.polaron"), err.to_string()))?,
            );
        }

        Ok(ResolvedElectronMaterial {
            name: name.to_string(),
            material: material.clone(),
            optical_elf,
            optical_elf_file: file,
            band,
            channels,
        })
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    const ELF: &str = "material = \"synthetic\"\nprovenance = \"synthetic test ELF\"\n\
                       energy_ev = [1.0, 10.0, 100.0]\nelf = [0.1, 1.0, 0.1]\n";

    fn dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "lindhard-electron-input-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    const GOOD: &str = r#"
[electron.beam]
energy_ev = 1000.0

[electron.transport]
cutoff_ev = 20.0

[electron.elastic]
potential = "thomas-fermi-yukawa"

[electron.materials.Si]
optical_elf = "elf.toml"

[target]
substrate = "Si"

[run]
histories = 10
seed = 1
"#;

    #[test]
    fn detects_electron_documents() {
        assert!(ElectronInput::is_electron_toml(GOOD));
        assert!(!ElectronInput::is_electron_toml("[beam]\nion = \"B\"\n"));
    }

    #[test]
    fn resolves_and_echoes_defaults() {
        let d = dir("good");
        std::fs::write(d.join("elf.toml"), ELF).unwrap();
        let i = ElectronInput::from_toml_str(GOOD).unwrap();
        let r = i.resolve_in(&d).unwrap();
        assert_eq!(r.layers.len(), 1);
        assert_eq!(
            r.materials[0].optical_elf.provenance(),
            "synthetic test ELF"
        );
        assert_eq!(r.materials[0].optical_elf_file.sha256.len(), 64);
        assert_eq!(r.input.electron.tables.max_energy_ev, Some(1000.0));
        assert_eq!(r.input.electron.tally.escape_energy_max_ev, Some(1000.0));
        assert_eq!(r.inelastic, PennAlgorithm::SinglePole);
        assert_eq!(*r.table_energy_ev.last().unwrap(), 1000.0);
        // The echo round-trips through TOML.
        let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
        assert_eq!(back, r.input);
    }

    #[test]
    fn data_without_provenance_is_refused() {
        let d = dir("noprov");
        std::fs::write(d.join("elf.toml"), ELF.replace("synthetic test ELF", " ")).unwrap();
        let e = ElectronInput::from_toml_str(GOOD)
            .unwrap()
            .resolve_in(&d)
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("electron.materials.Si.optical_elf") && e.contains("provenance"),
            "{e}"
        );
        // Band parameters without a provenance are refused too.
        std::fs::write(d.join("elf.toml"), ELF).unwrap();
        let text = GOOD.replace(
            "optical_elf = \"elf.toml\"",
            "optical_elf = \"elf.toml\"\nband = { kind = \"metal\", fermi_ev = 5.0, \
             work_function_ev = 4.0, provenance = \"\" }",
        );
        let e = ElectronInput::from_toml_str(&text)
            .unwrap()
            .resolve_in(&d)
            .unwrap_err()
            .to_string();
        assert!(e.contains("electron.materials.Si.band"), "{e}");
    }

    #[test]
    fn errors_name_the_field() {
        let d = dir("errors");
        std::fs::write(d.join("elf.toml"), ELF).unwrap();
        let bad = |text: &str| {
            ElectronInput::from_toml_str(text)
                .and_then(|i| i.resolve_in(&d))
                .unwrap_err()
                .to_string()
        };
        let e = bad(&GOOD.replace("cutoff_ev = 20.0", "cutoff_ev = -1.0"));
        assert!(e.contains("electron.transport.cutoff_ev"), "{e}");
        let e = bad(&GOOD.replace("[electron.materials.Si]", "[electron.materials.Ge]"));
        assert!(
            e.contains("target.substrate") && e.contains("electron data"),
            "{e}"
        );
        let e = bad(&GOOD.replace(
            "cutoff_ev = 20.0",
            "cutoff_ev = 20.0\nsecondaries = \"kieft-bosch\"",
        ));
        assert!(e.contains("electron.materials.Si.band"), "{e}");
        let e = bad(&GOOD.replace("potential = \"thomas-fermi-yukawa\"", "potential = \"x\""));
        assert!(e.contains("thomas-fermi-yukawa"), "{e}");
        let e = bad(&GOOD.replace(
            "[electron.elastic]",
            "[electron.inelastic]\nmodel = \"penn\"\n[electron.elastic]",
        ));
        assert!(e.contains("electron.inelastic.model"), "{e}");
        let e = bad(&GOOD.replace("energy_ev = 1000.0", "energy_kev = 1.0"));
        assert!(e.contains("energy_kev"), "{e}");
        let e = bad(&GOOD.replace(
            "potential = \"thomas-fermi-yukawa\"",
            "potential = \"thomas-fermi-yukawa\"\n[electron.elastic.correlation_polarization.polarizability]\n",
        ));
        assert!(e.contains("no polarizability for Si"), "{e}");
    }

    /// `fermi_energy_ev` sets the inelastic tables of materials without a
    /// band; with a band the table takes the band's Fermi energy, so a
    /// nonzero value is refused instead of counted twice (#241).
    #[test]
    fn inelastic_fermi_energy_is_refused_with_a_band() {
        let d = dir("fermi-band");
        std::fs::write(d.join("elf.toml"), ELF).unwrap();
        let fermi = |text: &str, ev: &str| {
            text.replace(
                "[electron.elastic]",
                &format!("[electron.inelastic]\nfermi_energy_ev = {ev}\n[electron.elastic]"),
            )
        };
        let r = ElectronInput::from_toml_str(&fermi(GOOD, "1.5"))
            .unwrap()
            .resolve_in(&d)
            .unwrap();
        assert_eq!(r.inelastic_fermi_ev, 1.5);
        let banded = GOOD.replace(
            "optical_elf = \"elf.toml\"",
            "optical_elf = \"elf.toml\"\nband = { kind = \"metal\", fermi_ev = 5.0, \
             work_function_ev = 4.0, provenance = \"synthetic test band\" }",
        );
        let r = ElectronInput::from_toml_str(&fermi(&banded, "0.0"))
            .unwrap()
            .resolve_in(&d)
            .unwrap();
        assert_eq!(r.inelastic_fermi_ev, 0.0);
        let e = ElectronInput::from_toml_str(&fermi(&banded, "1.5"))
            .unwrap()
            .resolve_in(&d)
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("electron.inelastic.fermi_energy_ev") && e.contains("material Si"),
            "{e}"
        );
    }

    #[test]
    fn psf_section_defaults_echo_and_errors_name_fields() {
        let d = dir("psf");
        std::fs::write(d.join("elf.toml"), ELF).unwrap();
        let psf = "\n[electron.tally.psf]\ndepth_lo_nm = 0.0\ndepth_hi_nm = 5.0\nr_min_nm = 0.5\nr_max_nm = 50.0\nbins = 8\n";
        let text = format!("{GOOD}{psf}");
        let r = ElectronInput::from_toml_str(&text)
            .unwrap()
            .resolve_in(&d)
            .unwrap();
        let c = r.tally.psf.unwrap();
        assert_eq!(c.depth_hi_m, 5e-9);
        assert_eq!(c.radial.bins, 8);
        assert_eq!(
            r.psf_fits,
            vec![PsfModel::DoubleGaussian, PsfModel::TripleGaussian]
        );
        assert_eq!(r.psf_normalization, PsfNormalization::SlabTotal);
        let echo = r.input.electron.tally.psf.clone().unwrap();
        assert_eq!(echo.fits, default_psf_fits());
        assert_eq!(echo.normalization, PsfNormalizationSpec::SlabTotal);
        let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
        assert_eq!(back, r.input);
        // No PSF: nothing resolved.
        let r0 = ElectronInput::from_toml_str(GOOD)
            .unwrap()
            .resolve_in(&d)
            .unwrap();
        assert!(r0.tally.psf.is_none() && r0.psf_fits.is_empty());
        // Empty and explicit fit lists, free normalisation.
        let r = ElectronInput::from_toml_str(&text.replace(
            "bins = 8",
            "bins = 8\nfits = [\"triple\", \"triple\"]\nnormalization = \"free\"",
        ))
        .unwrap()
        .resolve_in(&d)
        .unwrap();
        assert_eq!(r.psf_fits, vec![PsfModel::TripleGaussian]);
        assert_eq!(r.psf_normalization, PsfNormalization::Free);
        let r = ElectronInput::from_toml_str(&text.replace("bins = 8", "bins = 8\nfits = []"))
            .unwrap()
            .resolve_in(&d)
            .unwrap();
        assert!(r.psf_fits.is_empty());
        // Errors.
        let bad = |from: &str, to: &str| {
            ElectronInput::from_toml_str(&text.replace(from, to))
                .and_then(|i| i.resolve_in(&d))
                .unwrap_err()
                .to_string()
        };
        let e = bad("depth_hi_nm = 5.0", "depth_hi_nm = 0.0");
        assert!(e.contains("electron.tally.psf.depth_hi_nm"), "{e}");
        let e = bad("r_max_nm = 50.0", "r_max_nm = 0.5");
        assert!(e.contains("electron.tally.psf.r_max_nm"), "{e}");
        let e = bad("r_min_nm = 0.5", "r_min_nm = 0.0");
        assert!(e.contains("electron.tally.psf.r_min_nm"), "{e}");
        let e = bad("bins = 8", "bins = 0");
        assert!(e.contains("electron.tally.psf.bins"), "{e}");
        let e = bad("bins = 8", "bins = 8\nfits = [\"quadruple\"]");
        assert!(e.contains("quadruple") && e.contains("fits"), "{e}");
        let e = bad("bins = 8", "bins = 8\nnormalization = \"x\"");
        assert!(e.contains("normalization"), "{e}");
        let e = bad("bins = 8", "bins = 8\nbogus = 1");
        assert!(e.contains("bogus"), "{e}");
    }

    #[test]
    fn every_option_resolves_and_is_recorded() {
        let d = dir("rich");
        std::fs::write(d.join("elf.toml"), ELF).unwrap();
        let text = r#"
[electron.beam]
energy_ev = 2000.0
tilt_deg = 30.0
azimuth_deg = 10.0

[electron.transport]
cutoff_ev = 2.0
cutoff_reference = "vacuum-level"
escape_rule = "front-only"
max_events = 1000
secondaries = "kieft-bosch"
instantaneous_momentum = false
boundary = "step-barrier"
refraction = false

[electron.elastic]
potential = "thomas-fermi-yukawa"
exchange = true
[electron.elastic.correlation_polarization]
b_pol_squared = 1.0
polarizability.Si = { bohr3 = 40.0, source = "synthetic test value" }
polarizability.O = { bohr3 = 5.0, source = "synthetic test value" }

[electron.inelastic]
model = "mermin-melf"

[electron.materials.Ox]
optical_elf = "elf.toml"
band = { kind = "insulator", valence_band_width_ev = 6.0, band_gap_ev = 3.0, affinity_ev = 1.0, provenance = "synthetic test band" }
phonon = { preset = "sio2-63mev", temperature_k = 300.0 }
polaron = { c_per_nm = 0.1, gamma_per_ev = 0.5, provenance = "synthetic test polaron" }

[electron.materials.Si]
optical_elf = "elf.toml"
band = { kind = "free-electron-metal", valence_electrons_per_atom = 4.0, work_function_ev = 4.0, provenance = "synthetic test band" }

[electron.tally]
cartesian = { x = { lo_nm = 0.0, hi_nm = 10.0, bins = 2 }, y = { lo_nm = -5.0, hi_nm = 5.0, bins = 2 }, z = { lo_nm = -5.0, hi_nm = 5.0, bins = 2 } }
cylindrical = { r = { lo_nm = 0.0, hi_nm = 5.0, bins = 2 }, depth = { lo_nm = 0.0, hi_nm = 10.0, bins = 2 } }

[materials.Ox]
density_g_cm3 = 2.2
elements = [{ symbol = "Si", atom_fraction = 1.0 }, { symbol = "O", atom_fraction = 2.0 }]

[[target.layers]]
material = "Ox"
thickness_nm = 5.0

[target]
substrate = "Si"

[run]
histories = 10
seed = 3
"#;
        let r = ElectronInput::from_toml_str(text)
            .unwrap()
            .resolve_in(&d)
            .unwrap();
        assert_eq!(r.layer_material, vec![0, 1]);
        assert_eq!(r.materials[0].name, "Ox");
        assert!(r.materials[0].channels.phonon.is_some());
        assert!(r.materials[0].channels.polaron.is_some());
        assert!(r.materials[1].band.as_ref().unwrap().fermi_ev() > 0.0);
        assert_eq!(r.inelastic, PennAlgorithm::Mermin);
        assert_eq!(
            r.config.secondaries,
            SecondaryModel::KieftBosch {
                instantaneous_momentum: false,
                momentum_conservation: true
            }
        );
        assert_eq!(
            r.config.boundary,
            BoundaryModel::StepBarrier {
                quantum_transmission: true,
                refraction: false
            }
        );
        assert!(r.tally.cartesian.is_some() && r.tally.cylindrical.is_some());
        assert_eq!(r.elastic.correlation_polarization.len(), 2);
        let desc = r.elastic.corrections_description();
        assert!(desc.contains("exchange") && desc.contains("synthetic test value"));
        let names: Vec<_> = r.models().iter().map(|m| m.name).collect();
        for n in [
            "furness-mccarthy",
            "salvat-2003",
            "mermin-melf",
            "frohlich",
            "ganachaud-mokrani",
        ] {
            assert!(names.contains(&n), "{n}: {names:?}");
        }
        // The echo has the defaulted options filled in and round-trips.
        let t = &r.input.electron.transport;
        assert_eq!(t.momentum_conservation, Some(true));
        assert_eq!(t.quantum_transmission, Some(true));
        let max = r.input.electron.tables.max_energy_ev.unwrap();
        assert!(
            max > 2000.0,
            "beam energy plus the largest inner potential: {max}"
        );
        let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
        assert_eq!(back, r.input);
        // A preset fixes its own source.
        let e = ElectronInput::from_toml_str(&text.replace(
            "preset = \"sio2-63mev\", temperature_k = 300.0",
            "preset = \"sio2-63mev\", temperature_k = 300.0, provenance = \"x\"",
        ))
        .unwrap()
        .resolve_in(&d)
        .unwrap_err()
        .to_string();
        assert!(e.contains("electron.materials.Ox.phonon"), "{e}");
    }
}
