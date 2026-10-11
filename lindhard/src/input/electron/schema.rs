//! The serde schema of the electron input: the `*Spec` and `*Choice` types
//! and their defaults.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::electron::elastic::table::{DEFAULT_MIN_ENERGY_EV, DEFAULT_POINTS_PER_DECADE};
use crate::electron::inelastic::{FitWeighting, MerminFitOptions, PennAlgorithm};
use crate::electron::transport::TransportConfig;
use crate::input::TargetSpec;
use crate::material::MaterialSpec;
use crate::tally::{PsfModel, PsfNormalization, SE_BSE_SPLIT_EV};

// Named only by intra-doc links in this module.
#[cfg(doc)]
use crate::electron::{
    boundary::{BandModel, BandStructure},
    data::OpticalElf,
    phonon::{FrohlichPhonon, PolaronTrapping},
    secondary::SecondaryModel,
    transport::{BoundaryModel, CutoffReference, EscapeRule},
};

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
        /// Binding energy below the Fermi level of the electron a valence
        /// event liberates, eV ([`BandStructure::with_valence_binding_ev`]);
        /// absent: the Fermi level, the default.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        valence_binding_ev: Option<f64>,
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
        /// As for [`BandSpec::Metal`].
        #[serde(default, skip_serializing_if = "Option::is_none")]
        valence_binding_ev: Option<f64>,
        /// Source of the work function and of the valence count (and of
        /// `valence_binding_ev`, if given).
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

fn default_mermin_oscillators() -> usize {
    MerminFitOptions::default().n_oscillators
}

/// The `floor` of the default relative weighting of
/// [`MerminFitOptions::default`].
pub(super) fn default_relative_floor() -> f64 {
    match MerminFitOptions::default().weighting {
        FitWeighting::Relative { floor } => floor,
        FitWeighting::Uniform => unreachable!("the default fit weighting is relative"),
    }
}

/// Weighting of the residuals of the Mermin oscillator fit
/// ([`FitWeighting`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MerminWeightingChoice {
    /// [`FitWeighting::Relative`] (the default).
    #[default]
    Relative,
    /// [`FitWeighting::Uniform`].
    Uniform,
}

/// `[electron.materials.<name>.mermin_fit]`: the options of the oscillator
/// fit of the `mermin-melf` model to this material's optical ELF
/// ([`MerminFitOptions`],
/// [`crate::electron::inelastic::fit_mermin_oscillators`]). Absent: exactly
/// [`MerminFitOptions::default`]. Only `n_oscillators` and `weighting` are
/// exposed; the starting point and the iteration cap keep their defaults.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MerminFitSpec {
    /// Number of oscillators, `1..=MAX_OSCILLATORS` (16); default 3.
    #[serde(default = "default_mermin_oscillators")]
    pub oscillators: usize,
    /// `"relative"` (default) or `"uniform"`.
    #[serde(default)]
    pub weighting: MerminWeightingChoice,
    /// `floor` of the relative weighting, in `(0, 1]`; default 0.01. Not
    /// allowed with `"uniform"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relative_floor: Option<f64>,
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
    /// Options of the Mermin oscillator fit (`mermin-melf` only); absent:
    /// [`MerminFitOptions::default`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mermin_fit: Option<MerminFitSpec>,
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

pub(super) fn default_psf_fits() -> Vec<PsfFitSpec> {
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
