//! The run description: one TOML document describing the beam, the layered
//! target and its materials, the physics choices and the run size.
//!
//! The same types are serialized back into the output metadata
//! (`summary.json`), so a result carries the input that produced it. Every
//! table uses `deny_unknown_fields`, so a misspelt key is an error naming the
//! key rather than a silently ignored setting. Defaults are filled in on
//! parsing, so the echoed input shows every choice explicitly.
//!
//! The schema is documented with an example in `docs/cli.md`. In short:
//!
//! ```toml
//! [beam]
//! ion = "B"              # element symbol; mass_amu optional
//! energy_ev = 5000.0
//! tilt_deg = 7.0          # polar angle from the surface normal
//!
//! [target]
//! substrate = "Si"        # a name from [materials], or an element symbol
//!
//! [physics]
//! primary_cutoff_ev = 5.0
//! recoil_cutoff_ev = 2.0
//! [physics.energies.Si]
//! e_d_ev = 15.0
//!
//! [run]
//! ions = 1000
//! seed = 1
//! ```
//!
//! [`Input::resolve`] validates the whole description and turns it into the
//! engine's types ([`Resolved`]); [`Resolved::models`] lists every model in
//! use with its citation, for the output metadata.
//!
//! An electron run is described by a different document, with an
//! `[electron]` table and the same `[materials]` and `[target]`; its schema is
//! [`electron`].

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::path::{Path, PathBuf};

use serde::de::{self, MapAccess, Visitor};
use serde::{Deserialize, Deserializer, Serialize};
use sha2::{Digest, Sha256};

use crate::elements::element_by_symbol;
use crate::geometry::Stack;
use crate::ion::bca::{BcaConfig, Beam, ElectronicLoss, MeanFreePath, MAX_WEAK_COLLISIONS};
use crate::ion::potential::{Screening, ScreeningLength};
use crate::ion::scattering::TableSpec;
use crate::ion::stopping::bethe::BetheBloch;
use crate::ion::stopping::lindhard_scharff::LindhardScharff;
use crate::ion::stopping::table::{StoppingTable, TableOverride};
use crate::ion::stopping::{ElectronicStopping, Ion};
use crate::material::{EnergyKind, Material, MaterialSpec};

pub mod electron;

/// Errors from reading or validating an [`Input`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum InputError {
    /// The document is not valid TOML or does not match the schema (unknown
    /// key, wrong type, missing required key). The message from the parser
    /// names the key and the line.
    #[error("{0}")]
    Parse(String),
    /// A value is out of range or inconsistent. `field` is the dotted path of
    /// the offending key, e.g. `target.layers[1].thickness_nm`.
    #[error("{field}: {message}")]
    Invalid {
        /// Dotted path of the key.
        field: String,
        /// What is wrong with it.
        message: String,
    },
}

fn invalid(field: impl Into<String>, message: impl Into<String>) -> InputError {
    InputError::Invalid {
        field: field.into(),
        message: message.into(),
    }
}

/// The whole run description.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    /// The incident beam.
    pub beam: BeamSpec,
    /// Named materials, referenced by name from the target layers.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, MaterialSpec>,
    /// The layered target.
    pub target: TargetSpec,
    /// Physics model choices and cutoffs.
    pub physics: PhysicsSpec,
    /// User-supplied electronic stopping tables (optional).
    #[serde(default, skip_serializing_if = "StoppingSpec::is_empty")]
    pub stopping: StoppingSpec,
    /// Run size, seed and threads.
    pub run: RunSpec,
    /// What the run records.
    #[serde(default)]
    pub tally: TallySpec,
    /// Fluence-dependent target (optional). Absent: the static run.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dynamic: Option<DynamicSpec>,
}

/// `[beam]`: species, energy and direction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BeamSpec {
    /// Element symbol of the projectile, e.g. `"B"`.
    pub ion: String,
    /// Projectile mass, u. Default: the standard atomic weight.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_amu: Option<f64>,
    /// Incident kinetic energy, eV.
    pub energy_ev: f64,
    /// Polar angle of incidence from the surface normal, degrees, in `[0, 90)`.
    #[serde(default)]
    pub tilt_deg: f64,
    /// Azimuth of the incidence plane about the normal, degrees.
    #[serde(default)]
    pub azimuth_deg: f64,
}

/// A material given by name or inline.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(untagged)]
pub enum MaterialRef {
    /// A key of `[materials]`, or else an element symbol (a pure element at
    /// its tabulated density).
    Name(String),
    /// An inline material table.
    Inline(MaterialSpec),
}

impl<'de> Deserialize<'de> for MaterialRef {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = MaterialRef;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a material name (string) or an inline material table")
            }
            fn visit_str<E: de::Error>(self, s: &str) -> Result<MaterialRef, E> {
                Ok(MaterialRef::Name(s.to_owned()))
            }
            fn visit_map<A: MapAccess<'de>>(self, m: A) -> Result<MaterialRef, A::Error> {
                MaterialSpec::deserialize(de::value::MapAccessDeserializer::new(m))
                    .map(MaterialRef::Inline)
            }
        }
        d.deserialize_any(V)
    }
}

/// One finite layer, front to back.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LayerSpec {
    /// The layer material.
    pub material: MaterialRef,
    /// Thickness, nm.
    pub thickness_nm: f64,
}

/// `[target]`: finite layers front to back, then an optional semi-infinite
/// substrate. At least one of the two must be given.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TargetSpec {
    /// Finite layers, front first.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub layers: Vec<LayerSpec>,
    /// Semi-infinite substrate behind the layers. Without one the target has
    /// a back face and particles can be transmitted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub substrate: Option<MaterialRef>,
}

/// Screening function choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum PotentialChoice {
    /// [`Screening::ZblUniversal`].
    #[default]
    Zbl,
    /// [`Screening::KrC`].
    KrC,
    /// [`Screening::Moliere`].
    Moliere,
    /// [`Screening::LenzJensen`].
    LenzJensen,
}

impl PotentialChoice {
    /// The engine's screening function.
    pub fn screening(self) -> Screening {
        match self {
            Self::Zbl => Screening::ZblUniversal,
            Self::KrC => Screening::KrC,
            Self::Moliere => Screening::Moliere,
            Self::LenzJensen => Screening::LenzJensen,
        }
    }
}

/// Screening length choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum LengthChoice {
    /// [`ScreeningLength::Universal`].
    Universal,
    /// [`ScreeningLength::Firsov`].
    Firsov,
    /// [`ScreeningLength::Lindhard`].
    Lindhard,
}

impl LengthChoice {
    fn from_length(l: ScreeningLength) -> Self {
        match l {
            ScreeningLength::Universal => Self::Universal,
            ScreeningLength::Firsov => Self::Firsov,
            ScreeningLength::Lindhard => Self::Lindhard,
        }
    }

    /// The engine's screening length.
    pub fn length(self) -> ScreeningLength {
        match self {
            Self::Universal => ScreeningLength::Universal,
            Self::Firsov => ScreeningLength::Firsov,
            Self::Lindhard => ScreeningLength::Lindhard,
        }
    }
}

/// `[stopping]`: user-supplied electronic stopping tables.
///
/// ```toml
/// [stopping]
/// tables = ["tables/b_in_si.toml"]
/// ```
///
/// Each entry is a file in the [`StoppingTable`] format (one ion and target
/// element, `eV 1e-15 cm^2` per atom, mandatory `provenance`). Relative paths
/// resolve against the directory of the input file. A table replaces the
/// `[physics] stopping` model for the (ion, target element) pair it declares;
/// every other pair uses the `[physics] stopping` model. Queries outside a
/// table's energy range are errors, never extrapolated.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StoppingSpec {
    /// Paths of table files, one per (ion, target element) pair.
    #[serde(default)]
    pub tables: Vec<String>,
}

impl StoppingSpec {
    /// No tables declared.
    pub fn is_empty(&self) -> bool {
        self.tables.is_empty()
    }
}

/// A user stopping table after loading, with what a run record needs to
/// identify the exact file.
#[derive(Debug, Clone)]
pub struct LoadedStoppingTable {
    /// The path as written in the input (`stopping.tables[i]`).
    pub path: String,
    /// The path the file was read from (the input directory joined with
    /// `path`, made absolute when possible).
    pub resolved_path: PathBuf,
    /// SHA-256 of the file's bytes, lowercase hex.
    pub sha256: String,
    /// The validated table.
    pub table: StoppingTable,
}

/// Electronic stopping choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum StoppingChoice {
    /// Lindhard-Scharff, all loss nonlocal along the free flight.
    #[default]
    LindhardScharff,
    /// Bethe-Bloch (bare charge, Bloch-rule mean excitation energy), all loss
    /// nonlocal. Only meaningful well above the Lindhard-Scharff regime.
    BetheBloch,
    /// Equipartition of Lindhard-Scharff: half nonlocal, half local at each
    /// collision by Oen-Robinson ([`ElectronicLoss::EquipartitionLsOr`]).
    EquipartitionLsOr,
}

/// Free-path convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum FreePathChoice {
    /// [`MeanFreePath::Constant`].
    #[default]
    Constant,
    /// [`MeanFreePath::EnergyDependent`]; needs `min_cm_angle_deg`.
    EnergyDependent,
}

/// Per-element energy overrides, eV.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergyOverride {
    /// Displacement energy `E_d`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_d_ev: Option<f64>,
    /// Lattice binding energy `E_b`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_b_ev: Option<f64>,
    /// Surface binding energy `E_s`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_s_ev: Option<f64>,
}

fn default_true() -> bool {
    true
}

/// `[physics]`: model choices, cutoffs and energy overrides.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhysicsSpec {
    /// Screening function. Default `"zbl"`.
    #[serde(default)]
    pub potential: PotentialChoice,
    /// Screening length. Default: the length conventionally paired with the
    /// screening function ([`Screening::default_length`]); filled in on
    /// resolution and echoed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub screening_length: Option<LengthChoice>,
    /// Electronic stopping. Default `"lindhard-scharff"`.
    #[serde(default)]
    pub stopping: StoppingChoice,
    /// Free-path convention. Default `"constant"`.
    #[serde(default)]
    pub free_path: FreePathChoice,
    /// Smallest centre-of-mass deflection treated as a collision, degrees;
    /// required with `free_path = "energy-dependent"`, rejected otherwise.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_cm_angle_deg: Option<f64>,
    /// Weak collisions per collision step, `0..=3`
    /// ([`crate::ion::bca::BcaConfig::weak_collisions`]); constant free path
    /// only. Default 0.
    #[serde(default, skip_serializing_if = "is_zero_u8")]
    pub weak_collisions: u8,
    /// The primary stops below this energy, eV. Required.
    pub primary_cutoff_ev: f64,
    /// Recoils stop below this energy, eV. Required. Keep it below the
    /// smallest `E_s` when sputtering matters.
    pub recoil_cutoff_ev: f64,
    /// Follow displaced atoms as full cascades. Default `true`.
    #[serde(default = "default_true")]
    pub follow_recoils: bool,
    /// Planar surface barrier for the beam species, eV. Default 0.
    #[serde(default)]
    pub primary_surface_binding_ev: f64,
    /// `E_d`/`E_b`/`E_s` overrides by element symbol, applied to that element
    /// in every layer after the material's own values (so they win).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub energies: BTreeMap<String, EnergyOverride>,
    /// Opt-in phenomenological calibration: `"none"` (default) or the name of
    /// a versioned factor set ([`TuningSet`]). Not a published model choice.
    /// See [`TuningSet`] for what it scales and where it applies.
    #[serde(default = "default_tuning", skip_serializing_if = "is_no_tuning")]
    pub tuning: String,
}

/// The `tuning` value that leaves the physics untouched.
pub const NO_TUNING: &str = "none";

fn default_tuning() -> String {
    NO_TUNING.to_string()
}

fn is_no_tuning(s: &String) -> bool {
    s == NO_TUNING
}

/// A named, immutable, versioned set of phenomenological surface-binding
/// energy multipliers (the pilot of the opt-in "experiment-tuned" mode).
///
/// This is calibration, not a published model choice: the factors stand for
/// the uncertain planar-barrier convention (`E_s` = cohesive energy). Each
/// record names its measured source data and fit recipe in `provenance`
/// (and in `docs/data-provenance.md`). Rules, all enforced on resolution:
///
/// - Only `[physics] tuning = "<name>"` selects a set; omission and `"none"`
///   change nothing, and the echoed input is then identical to an untuned one.
/// - A factor multiplies the *resolved* `E_s` of its element (the explicit
///   `[physics.energies]` / material value if present, else the elemental
///   default), exactly once per layer component, after overrides. The global
///   element table and the collision algorithm are never touched.
/// - Pilot scope: ion runs on single-element layers only, with a beam species
///   the set was fitted for ([`TuningSet::ions`]) and a target element the set
///   lists. Compounds, a `[dynamic]` target, other beams, unlisted elements and
///   electron input are rejected. A beam energy outside the fitted range or a
///   tilted beam is allowed but warned about (an extrapolation).
/// - Both the original and the effective energies are reported in
///   `summary.json` (`physics.tuning`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TuningSet {
    /// Name used in `[physics] tuning`. Names are never reused: a new version
    /// ships under a new name, so a name always selects the same factors.
    pub name: &'static str,
    /// Version of this set; a changed factor is a new version, never an edit.
    pub version: u32,
    /// Beam species (element symbols) the factors were fitted for.
    pub ions: &'static [&'static str],
    /// Beam energies (eV, inclusive) the fit data covered, at normal incidence.
    pub energy_range_ev: (f64, f64),
    /// `E_s` multiplier by element symbol; a target element not listed is
    /// rejected.
    pub e_s_factors: &'static [(&'static str, f64)],
    /// Source measurements and fit procedure, in one line.
    pub provenance: &'static str,
}

/// `es-sputter-ar-v1`: one `E_s` multiplier per element for Ar sputtering of
/// Si, Cu, Ag and Au, fitted by `validation/experiments/run.py --fit-tuning`
/// to measured yields only (the sets in `validation/data/sputtering/`, each
/// cited to its original publication in `docs/data-provenance.md`), with a
/// held-out evaluation in `validation/experiments/tuning_results.json`.
///
/// Phenomenological calibration, not a published value: the factors absorb
/// whatever the untuned matched problem misses (not shown to be `E_s` itself)
/// and were fitted under those matched settings (ZBL, Lindhard-Scharff,
/// constant free path, `E_d` = untuned `E_s`, `E_b` = 0, cutoffs 2 / 1 eV,
/// no weak collisions); their transfer to other settings is untested.
pub const ES_SPUTTER_AR_V1: TuningSet = TuningSet {
    name: "es-sputter-ar-v1",
    version: 1,
    ions: &["Ar"],
    energy_range_ev: (196.0, 10020.0),
    e_s_factors: &[("Si", 0.65), ("Cu", 0.55), ("Ag", 0.30), ("Au", 0.45)],
    provenance: "E_s multipliers fitted to measured Ar sputter yields (validation/data/sputtering, \
                 training sets only; held-out sets scored separately) by validation/experiments/run.py \
                 --fit-tuning; record: validation/experiments/tuning_results.json, docs/data-provenance.md \
                 (Tuning factor sets); issue #80",
};

/// The sets that ship with the engine. Each has a reproducible fit record and
/// a held-out evaluation (`docs/data-provenance.md`, "Tuning factor sets").
pub const TUNING_SETS: &[TuningSet] = &[ES_SPUTTER_AR_V1];

/// One component's tuning, as applied.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TunedEnergy {
    /// Index into [`Resolved::layers`].
    pub layer: usize,
    /// Element symbol.
    pub element: String,
    /// Resolved `E_s` before the factor, eV.
    pub e_s_original_ev: f64,
    /// The multiplier (1 when the set has no entry for the element).
    pub factor: f64,
    /// `E_s` used by the run, eV.
    pub e_s_effective_ev: f64,
}

/// What a tuning set did to a run; present in [`Resolved`] only when enabled.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TuningReport {
    /// Set name.
    pub set: String,
    /// Set version.
    pub version: u32,
    /// The set's provenance line.
    pub provenance: String,
    /// Always `"surface-binding-energy"` in the pilot.
    pub quantity: &'static str,
    /// Per layer component.
    pub components: Vec<TunedEnergy>,
}

/// `[run]`: history count, seed and threads.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunSpec {
    /// Number of primary ions (histories).
    pub ions: u64,
    /// Run seed. Results depend only on the seed and the input, never on the
    /// thread count.
    pub seed: u64,
    /// Worker threads. Default: all available cores. Does not change results,
    /// so it is not part of the echoed input.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub threads: Option<usize>,
}

/// Volume-relaxation convention of a dynamic run
/// ([`crate::ion::dynamic::Relaxation`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum RelaxationChoice {
    /// Ideal mixing of atomic volumes; every element at its elemental solid
    /// volume unless `atomic_volume_nm3` gives one.
    #[default]
    IdealMixing,
    /// One total atom number density for every slab; needs
    /// `number_density_cm3`.
    FixedNumberDensity,
}

fn default_min_ions() -> u64 {
    1
}

/// `[dynamic]`: a fluence-dependent target. The run's `ions` are delivered in
/// steps; after each step the target composition is updated from the
/// transport events (module docs of [`crate::ion::dynamic`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DynamicSpec {
    /// Total fluence of the whole run, ions/cm². The `run.ions` histories
    /// represent it, each standing for `fluence_cm2 / ions` ions/cm².
    pub fluence_cm2: f64,
    /// Ions per step; with `max_change` set, the largest (and first) step.
    pub ions_per_step: u64,
    /// Adaptive steps: the largest relative composition change per step (see
    /// [`crate::ion::dynamic::DynamicRun`]). Absent: fixed steps.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_change: Option<f64>,
    /// Adaptive steps: the smallest step. Default 1.
    #[serde(default = "default_min_ions")]
    pub min_ions_per_step: u64,
    /// Split every finite layer into slabs at most this thick, nm. Default
    /// (absent): each layer is one slab.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub slab_nm: Option<f64>,
    /// Volume relaxation convention. Default `"ideal-mixing"`.
    #[serde(default)]
    pub relaxation: RelaxationChoice,
    /// Total atom number density, atoms/cm³ (`"fixed-number-density"` only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub number_density_cm3: Option<f64>,
    /// Atomic volume per element symbol, nm³/atom (`"ideal-mixing"` only);
    /// overrides the elemental solid volume. Required for an element with no
    /// tabulated solid density (a gas, for example).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub atomic_volume_nm3: BTreeMap<String, f64>,
    /// Energies for elements that only enter the target during the run (the
    /// beam species, for example); an element already in a layer keeps the
    /// energies of that layer. Falls back to `[physics.energies]`, then to
    /// the element defaults.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub energies: BTreeMap<String, EnergyOverride>,
    /// Sputter erosion and surface recession: sputtered atoms are removed from
    /// the front of the target and the surface follows (see
    /// [`crate::ion::dynamic::DynamicRun`], "Erosion"). Default `false`, which
    /// keeps the surface fixed at `x = 0`.
    #[serde(default, skip_serializing_if = "is_false")]
    pub erosion: bool,
}

fn is_false(v: &bool) -> bool {
    !*v
}

fn is_zero_u8(v: &u8) -> bool {
    *v == 0
}

fn default_bin_nm() -> f64 {
    1.0
}
fn default_bins() -> usize {
    1000
}
fn default_lateral_bins() -> usize {
    100
}
fn default_escape_energy_bins() -> usize {
    100
}
fn default_escape_polar_bins() -> usize {
    30
}

/// `[tally]`: what the run records.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TallySpec {
    /// Depth-bin width of the stopped-primary profile, nm. Default 1.
    #[serde(default = "default_bin_nm")]
    pub depth_bin_nm: f64,
    /// Number of depth bins; the last also collects everything deeper.
    /// Default 1000.
    #[serde(default = "default_bins")]
    pub depth_bins: usize,
    /// Write the final state of every primary ion. Default `true`.
    #[serde(default = "default_true")]
    pub per_ion: bool,
    /// Bin width of the lateral (`y`, `z`) and radial profiles of stopped
    /// primaries, nm. Default 1.
    #[serde(default = "default_bin_nm")]
    pub lateral_bin_nm: f64,
    /// Lateral bins per side of the beam axis: `y` and `z` are binned over
    /// `[-lateral_bins * lateral_bin_nm, +lateral_bins * lateral_bin_nm)` and
    /// the radial distance over `[0, lateral_bins * lateral_bin_nm)`.
    /// Default 100.
    #[serde(default = "default_lateral_bins")]
    pub lateral_bins: usize,
    /// Upper edge of the escape-energy spectra, eV. Default (absent): the
    /// beam energy.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub escape_energy_max_ev: Option<f64>,
    /// Number of escape-energy bins over `[0, escape_energy_max_ev)`.
    /// Default 100.
    #[serde(default = "default_escape_energy_bins")]
    pub escape_energy_bins: usize,
    /// Number of escape polar-angle bins over `[0, 90)` degrees from the
    /// outward surface normal. Default 30.
    #[serde(default = "default_escape_polar_bins")]
    pub escape_polar_bins: usize,
    /// Also fit a dual-Pearson profile to the depth histogram. Default
    /// `false`.
    #[serde(default)]
    pub dual_pearson: bool,
}

impl Default for TallySpec {
    fn default() -> Self {
        Self {
            depth_bin_nm: default_bin_nm(),
            depth_bins: default_bins(),
            per_ion: true,
            lateral_bin_nm: default_bin_nm(),
            lateral_bins: default_lateral_bins(),
            escape_energy_max_ev: None,
            escape_energy_bins: default_escape_energy_bins(),
            escape_polar_bins: default_escape_polar_bins(),
            dual_pearson: false,
        }
    }
}

/// One model in use, for the output metadata.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct ModelInfo {
    /// What the model is for, e.g. `"screening function"`.
    pub role: &'static str,
    /// Stable model name.
    pub name: &'static str,
    /// Published source (for a user table: its path and provenance).
    pub citation: Cow<'static, str>,
}

/// One target layer after resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedLayer {
    /// How the input referred to the material (`[materials]` key, element
    /// symbol, or `"inline"`).
    pub source: String,
    /// The material with every energy set.
    pub material: Material,
}

/// A validated input, in the engine's types.
#[derive(Debug, Clone)]
pub struct Resolved {
    /// The input with defaults filled in (including the screening length).
    pub input: Input,
    /// The beam.
    pub beam: Beam,
    /// The target.
    pub stack: Stack,
    /// Layers in stack order, with where each material came from.
    pub layers: Vec<ResolvedLayer>,
    /// Engine configuration (including the seed).
    pub config: BcaConfig,
    /// Screening function.
    pub screening: Screening,
    /// Screening length.
    pub screening_length: ScreeningLength,
    /// Scattering-table grid.
    pub table_spec: TableSpec,
    /// User stopping tables, in input order (empty without `[stopping]`).
    pub stopping_tables: Vec<LoadedStoppingTable>,
    /// Non-fatal advice (e.g. beam energy outside a model's validity range).
    pub warnings: Vec<String>,
    /// The tuning applied, if `[physics] tuning` named a set.
    pub tuning: Option<TuningReport>,
}

/// Scattering-table grid used for every run: `per_decade` points per decade
/// of reduced energy and impact parameter. Angles outside the energy range
/// fall back to direct quadrature in the engine; impact parameters above
/// `beta_max` give no deflection (the constant-free-path `p_max` is about 20
/// screening lengths or less for solid densities). The measured interpolation
/// error is reported in the output metadata.
pub const TABLE_SPEC: TableSpec = TableSpec {
    eps_min: 1e-6,
    eps_max: 1e4,
    beta_min: 1e-5,
    beta_max: 1e2,
    per_decade: 32,
};

/// Resolve a [`MaterialRef`] against the `[materials]` table: a key of
/// `materials`, else an element symbol (the pure element at its tabulated
/// density), or an inline table. `field` names the key in errors. Shared by
/// the ion input and the electron input ([`electron`]).
pub(crate) fn resolve_material(
    materials: &BTreeMap<String, MaterialSpec>,
    field: &str,
    r: &MaterialRef,
) -> Result<ResolvedLayer, InputError> {
    match r {
        MaterialRef::Inline(spec) => Material::try_from(spec.clone())
            .map(|material| ResolvedLayer {
                source: "inline".into(),
                material,
            })
            .map_err(|e| invalid(field, e.to_string())),
        MaterialRef::Name(name) => {
            if let Some(spec) = materials.get(name) {
                let mut material = Material::try_from(spec.clone())
                    .map_err(|e| invalid(format!("materials.{name}"), e.to_string()))?;
                if material.name().is_none() {
                    material = material.with_name(name.clone());
                }
                return Ok(ResolvedLayer {
                    source: name.clone(),
                    material,
                });
            }
            let el = element_by_symbol(name).ok_or_else(|| {
                invalid(
                    field,
                    format!(
                        "unknown material {name:?}: not a key of [materials] \
                         and not an element symbol"
                    ),
                )
            })?;
            let material = Material::from_atom_fractions(&[(el.z, 1.0)], None)
                .map_err(|e| {
                    invalid(
                        field,
                        format!("element {name}: {e}; define it in [materials] with a density"),
                    )
                })?
                .with_name(name.clone());
            Ok(ResolvedLayer {
                source: name.clone(),
                material,
            })
        }
    }
}

fn finite_pos(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

impl Input {
    /// Parse a TOML document. Schema errors (unknown or missing keys, wrong
    /// types) are reported with the key and line; values are checked by
    /// [`Input::resolve`].
    pub fn from_toml_str(text: &str) -> Result<Self, InputError> {
        toml::from_str(text).map_err(|e| InputError::Parse(e.to_string()))
    }

    /// The input as a TOML document in the same schema [`Input::from_toml_str`]
    /// reads: parsing the text gives back an equal [`Input`].
    pub fn to_toml_string(&self) -> Result<String, InputError> {
        toml::to_string(self).map_err(|e| InputError::Parse(e.to_string()))
    }

    /// The input as echoed into output metadata: defaults filled in, and
    /// `run.threads` removed (it does not affect results).
    pub fn echo(&self) -> Input {
        let mut e = self.clone();
        e.run.threads = None;
        if e.physics.screening_length.is_none() {
            e.physics.screening_length = Some(LengthChoice::from_length(
                e.physics.potential.screening().default_length(),
            ));
        }
        e
    }

    /// Validate `[dynamic]`.
    fn check_dynamic(&self, d: &DynamicSpec) -> Result<(), InputError> {
        if !finite_pos(d.fluence_cm2) {
            return Err(invalid(
                "dynamic.fluence_cm2",
                "must be finite and positive",
            ));
        }
        if d.ions_per_step == 0 {
            return Err(invalid("dynamic.ions_per_step", "must be at least 1"));
        }
        match d.max_change {
            Some(c) if !finite_pos(c) => {
                return Err(invalid("dynamic.max_change", "must be finite and positive"))
            }
            Some(_) => {
                if d.min_ions_per_step == 0 || d.min_ions_per_step > d.ions_per_step {
                    return Err(invalid(
                        "dynamic.min_ions_per_step",
                        "must be between 1 and ions_per_step",
                    ));
                }
            }
            None => {
                if d.min_ions_per_step != default_min_ions() {
                    return Err(invalid(
                        "dynamic.min_ions_per_step",
                        "only used with max_change (adaptive steps)",
                    ));
                }
            }
        }
        if let Some(t) = d.slab_nm {
            if !finite_pos(t) {
                return Err(invalid("dynamic.slab_nm", "must be finite and positive"));
            }
        }
        match (d.relaxation, d.number_density_cm3) {
            (RelaxationChoice::FixedNumberDensity, None) => {
                return Err(invalid(
                    "dynamic.number_density_cm3",
                    "required with relaxation = \"fixed-number-density\"",
                ))
            }
            (RelaxationChoice::FixedNumberDensity, Some(n)) if !finite_pos(n) => {
                return Err(invalid(
                    "dynamic.number_density_cm3",
                    "must be finite and positive",
                ))
            }
            (RelaxationChoice::IdealMixing, Some(_)) => {
                return Err(invalid(
                    "dynamic.number_density_cm3",
                    "only used with relaxation = \"fixed-number-density\"",
                ))
            }
            _ => {}
        }
        if d.relaxation == RelaxationChoice::FixedNumberDensity && !d.atomic_volume_nm3.is_empty() {
            return Err(invalid(
                "dynamic.atomic_volume_nm3",
                "only used with relaxation = \"ideal-mixing\"",
            ));
        }
        for (sym, v) in &d.atomic_volume_nm3 {
            let field = format!("dynamic.atomic_volume_nm3.{sym}");
            if element_by_symbol(sym).is_none() {
                return Err(invalid(&field, format!("unknown element symbol {sym:?}")));
            }
            if !finite_pos(*v) {
                return Err(invalid(&field, "must be finite and positive"));
            }
        }
        for (sym, o) in &d.energies {
            let field = format!("dynamic.energies.{sym}");
            if element_by_symbol(sym).is_none() {
                return Err(invalid(&field, format!("unknown element symbol {sym:?}")));
            }
            for (key, v) in [
                ("e_d_ev", o.e_d_ev),
                ("e_b_ev", o.e_b_ev),
                ("e_s_ev", o.e_s_ev),
            ] {
                if v.is_some_and(|v| !(v.is_finite() && v >= 0.0)) {
                    return Err(invalid(
                        format!("{field}.{key}"),
                        "must be finite and non-negative",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Read and validate every `[stopping]` table file.
    fn load_stopping_tables(
        &self,
        base_dir: &Path,
    ) -> Result<Vec<LoadedStoppingTable>, InputError> {
        let mut loaded: Vec<LoadedStoppingTable> = Vec::new();
        for (i, path) in self.stopping.tables.iter().enumerate() {
            let field = format!("stopping.tables[{i}]");
            let joined = base_dir.join(path);
            let bytes = std::fs::read(&joined)
                .map_err(|e| invalid(&field, format!("cannot read {}: {e}", joined.display())))?;
            let text = std::str::from_utf8(&bytes).map_err(|_| {
                invalid(
                    &field,
                    format!("{} is not valid UTF-8 text", joined.display()),
                )
            })?;
            let table = StoppingTable::from_toml_str(text)
                .map_err(|e| invalid(&field, format!("{path}: {e}")))?;
            if let Some(prev) = loaded
                .iter()
                .position(|l| l.table.covers_pair(table.ion_z(), table.target_z()))
            {
                return Err(invalid(
                    &field,
                    format!(
                        "{path} declares the same ion/target pair (Z1={}, Z2={}) as \
                         stopping.tables[{prev}]; give each pair one table",
                        table.ion_z(),
                        table.target_z()
                    ),
                ));
            }
            let sha256 = Sha256::digest(&bytes)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            loaded.push(LoadedStoppingTable {
                path: path.clone(),
                resolved_path: std::fs::canonicalize(&joined).unwrap_or(joined),
                sha256,
                table,
            });
        }
        Ok(loaded)
    }

    /// Checks of the loaded tables against the rest of the run.
    fn check_stopping_tables(
        &self,
        tables: &[LoadedStoppingTable],
        ion: &Ion,
        layers: &[ResolvedLayer],
        warnings: &mut Vec<String>,
    ) -> Result<(), InputError> {
        if tables.is_empty() {
            return Ok(());
        }
        if self.physics.stopping == StoppingChoice::EquipartitionLsOr {
            return Err(invalid(
                "stopping.tables",
                "not usable with physics.stopping = \"equipartition-ls-or\", which carries \
                 its own Lindhard-Scharff/Oen-Robinson loss and would ignore the tables; \
                 use \"lindhard-scharff\" (or another model) as the fallback",
            ));
        }
        let e0 = self.beam.energy_ev;
        let cutoff = self.physics.primary_cutoff_ev;
        let recoil_cutoff = self.physics.recoil_cutoff_ev;
        let follow = self.physics.follow_recoils;
        let tol = crate::ion::stopping::table::MASS_RELATIVE_TOLERANCE;
        let target_elements: std::collections::BTreeSet<u8> = layers
            .iter()
            .flat_map(|ly| ly.material.components().iter().map(|c| c.z()))
            .collect();
        let symbol = |z: u8| crate::elements::element(z).map_or("?", |e| e.symbol);
        for (i, l) in tables.iter().enumerate() {
            let field = format!("stopping.tables[{i}]");
            let t = &l.table;
            let in_target = target_elements.contains(&t.target_z());
            // A table serves the beam ion of its element and, when recoils
            // are followed, every recoil of its element (`TableOverride`
            // keys on Z1, Z2 only).
            let serves_beam = t.ion_z() == ion.z();
            let serves_recoils = follow && target_elements.contains(&t.ion_z());
            if !in_target || !(serves_beam || serves_recoils) {
                let recoils = if follow {
                    ""
                } else {
                    " and physics.follow_recoils = false"
                };
                warnings.push(format!(
                    "{field}: the table ({}->{}, Z1={}, Z2={}) is for a pair that does not \
                     occur in this run (beam {}, target elements {}{recoils}); it is unused",
                    symbol(t.ion_z()),
                    symbol(t.target_z()),
                    t.ion_z(),
                    t.target_z(),
                    self.beam.ion,
                    target_elements
                        .iter()
                        .map(|&z| symbol(z))
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
                continue;
            }
            let (lo, hi) = t.energy_range_ev();
            if serves_recoils {
                // Recoils are built with the standard atomic weight
                // (`Ion::new`), so a table for another mass can never serve
                // them, and the run would stop at the first such recoil.
                let recoil = Ion::new(t.ion_z()).map_err(|e| invalid(&field, e.to_string()))?;
                if (recoil.mass_amu() - t.ion_mass_amu()).abs() > tol * t.ion_mass_amu() {
                    return Err(invalid(
                        &field,
                        format!(
                            "{}: table is for an ion of mass {} u, but with \
                             physics.follow_recoils = true it also serves the {} recoils, \
                             which have the standard atomic weight {} u; a table serves one \
                             mass per pair (use the standard weight in the table and the beam, \
                             drop the table, or set physics.follow_recoils = false)",
                            l.path,
                            t.ion_mass_amu(),
                            symbol(t.ion_z()),
                            recoil.mass_amu()
                        ),
                    ));
                }
                // Every recoil is followed down to recoil_cutoff_ev.
                if recoil_cutoff < lo {
                    return Err(invalid(
                        &field,
                        format!(
                            "{}: the table starts at {lo} eV, above physics.recoil_cutoff_ev = \
                             {recoil_cutoff} eV; with physics.follow_recoils = true it serves \
                             the {} recoils, which slow down to the recoil cutoff, so the run \
                             would stop outside the table range (extend the table down to the \
                             cutoff or raise physics.recoil_cutoff_ev)",
                            l.path,
                            symbol(t.ion_z()),
                        ),
                    ));
                }
                // The most a beam collision can hand a recoil of this element.
                let (m1, m2) = (ion.mass_amu(), recoil.mass_amu());
                let e_max = 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2)) * e0;
                if hi < e_max {
                    warnings.push(format!(
                        "{field}: the table ends at {hi} eV, below the largest energy a \
                         {} recoil can receive from the beam ({e_max:.4e} eV); the run fails \
                         if a recoil starts above the table range",
                        symbol(t.ion_z())
                    ));
                }
            }
            if !serves_beam {
                continue;
            }
            // The beam ion's pair: mass and the energy range it travels.
            if (ion.mass_amu() - t.ion_mass_amu()).abs() > tol * t.ion_mass_amu() {
                return Err(invalid(
                    &field,
                    format!(
                        "{}: table is for an ion of mass {} u but the beam ion has mass {} u \
                         (set `ion_mass_amu` in the table, or beam.mass_amu)",
                        l.path,
                        t.ion_mass_amu(),
                        ion.mass_amu()
                    ),
                ));
            }
            if e0 < lo || e0 > hi {
                return Err(invalid(
                    &field,
                    format!(
                        "{}: beam energy {e0} eV is outside the table range [{lo}, {hi}] eV",
                        l.path
                    ),
                ));
            }
            if cutoff < lo {
                warnings.push(format!(
                    "{field}: the table starts at {lo} eV, above physics.primary_cutoff_ev = \
                     {cutoff} eV; the run fails if a projectile slows below the table range"
                ));
            }
        }
        Ok(())
    }

    fn material(&self, field: &str, r: &MaterialRef) -> Result<ResolvedLayer, InputError> {
        resolve_material(&self.materials, field, r)
    }

    /// Validate everything and build the engine's inputs. Errors name the
    /// offending key.
    pub fn resolve(&self) -> Result<Resolved, InputError> {
        self.resolve_in(Path::new("."))
    }

    /// Like [`Input::resolve`], with relative `[stopping]` table paths taken
    /// relative to `base_dir` (the directory of the input file).
    pub fn resolve_in(&self, base_dir: &Path) -> Result<Resolved, InputError> {
        self.resolve_with_sets(base_dir, TUNING_SETS)
    }

    /// [`Input::resolve_in`] against an explicit registry of tuning sets
    /// (tests use synthetic sets; production uses [`TUNING_SETS`]).
    pub fn resolve_with_sets(
        &self,
        base_dir: &Path,
        sets: &[TuningSet],
    ) -> Result<Resolved, InputError> {
        let mut warnings = Vec::new();

        // Beam.
        let b = &self.beam;
        let el = element_by_symbol(&b.ion).ok_or_else(|| {
            invalid(
                "beam.ion",
                format!(
                    "unknown element symbol {:?} (symbols are case-sensitive)",
                    b.ion
                ),
            )
        })?;
        let ion = match b.mass_amu {
            None => Ion::new(el.z),
            Some(m) => Ion::with_mass(el.z, m),
        }
        .map_err(|e| invalid("beam.mass_amu", e.to_string()))?;
        if !finite_pos(b.energy_ev) {
            return Err(invalid("beam.energy_ev", "must be finite and positive"));
        }
        if !(b.tilt_deg.is_finite() && (0.0..90.0).contains(&b.tilt_deg)) {
            return Err(invalid("beam.tilt_deg", "must be in [0, 90) degrees"));
        }
        if !b.azimuth_deg.is_finite() {
            return Err(invalid("beam.azimuth_deg", "must be finite"));
        }
        if self.run.ions == 0 {
            return Err(invalid("run.ions", "must be at least 1"));
        }
        if self.run.threads == Some(0) {
            return Err(invalid(
                "run.threads",
                "must be at least 1 (omit for all cores)",
            ));
        }
        let beam = Beam {
            ion,
            energy_ev: b.energy_ev,
            polar_rad: b.tilt_deg.to_radians(),
            azimuth_rad: b.azimuth_deg.to_radians(),
            count: self.run.ions,
        };

        // Materials (all of them, so a broken unused entry is still reported).
        for (name, spec) in &self.materials {
            Material::try_from(spec.clone())
                .map_err(|e| invalid(format!("materials.{name}"), e.to_string()))?;
        }

        // Target.
        let t = &self.target;
        if t.layers.is_empty() && t.substrate.is_none() {
            return Err(invalid(
                "target",
                "needs at least one [[target.layers]] entry or a substrate",
            ));
        }
        let mut layers = Vec::with_capacity(t.layers.len() + 1);
        let mut fields = Vec::with_capacity(t.layers.len() + 1);
        let mut thick = Vec::with_capacity(t.layers.len());
        for (i, l) in t.layers.iter().enumerate() {
            let field = format!("target.layers[{i}]");
            if !finite_pos(l.thickness_nm) {
                return Err(invalid(
                    format!("{field}.thickness_nm"),
                    format!("{} nm must be finite and positive", l.thickness_nm),
                ));
            }
            layers.push(self.material(&format!("{field}.material"), &l.material)?);
            fields.push(format!("{field}.material"));
            thick.push(l.thickness_nm * 1e-9);
        }
        if let Some(s) = &t.substrate {
            layers.push(self.material("target.substrate", s)?);
            fields.push("target.substrate".to_string());
        }

        if let Some(d) = &self.dynamic {
            self.check_dynamic(d)?;
        }

        // Energy overrides.
        let p = &self.physics;
        for (sym, o) in &p.energies {
            let field = format!("physics.energies.{sym}");
            let el = element_by_symbol(sym)
                .ok_or_else(|| invalid(&field, format!("unknown element symbol {sym:?}")))?;
            let mut used = false;
            for l in &mut layers {
                if l.material.components().iter().all(|c| c.z() != el.z) {
                    continue;
                }
                used = true;
                let m = &mut l.material;
                let set = |r: Result<(), _>, key: &str| {
                    r.map_err(|e: crate::material::MaterialError| {
                        invalid(format!("{field}.{key}"), e.to_string())
                    })
                };
                if let Some(v) = o.e_d_ev {
                    set(m.set_displacement_energy_ev(el.z, v), "e_d_ev")?;
                }
                if let Some(v) = o.e_b_ev {
                    set(m.set_lattice_binding_energy_ev(el.z, v), "e_b_ev")?;
                }
                if let Some(v) = o.e_s_ev {
                    set(m.set_surface_binding_energy_ev(el.z, v), "e_s_ev")?;
                }
            }
            if !used {
                return Err(invalid(&field, format!("{sym} is not in any target layer")));
            }
        }
        for (l, field) in layers.iter().zip(&fields) {
            if let Some(&(z, kind)) = l.material.unset_energies().first() {
                let sym = crate::elements::element(z).expect("validated").symbol;
                let key = match kind {
                    EnergyKind::Displacement => "e_d_ev",
                    EnergyKind::LatticeBinding => "e_b_ev",
                    EnergyKind::SurfaceBinding => "e_s_ev",
                };
                return Err(invalid(
                    field.as_str(),
                    format!(
                        "{kind} of {sym} has no default; set `{key}` for {sym} in the \
                         material, or in [physics.energies.{sym}]"
                    ),
                ));
            }
        }

        let tuning = self.apply_tuning(sets, &mut layers, &mut warnings)?;

        let mut finite = Vec::with_capacity(thick.len());
        let mut it = layers.iter();
        for &th in &thick {
            finite.push((it.next().expect("one per layer").material.clone(), th));
        }
        let substrate = it.next().map(|l| l.material.clone());
        let stack = Stack::new(finite, substrate).map_err(|e| invalid("target", e.to_string()))?;

        // Physics.
        let mean_free_path = match (p.free_path, p.min_cm_angle_deg) {
            (FreePathChoice::Constant, None) => MeanFreePath::Constant,
            (FreePathChoice::Constant, Some(_)) => {
                return Err(invalid(
                    "physics.min_cm_angle_deg",
                    "only used with free_path = \"energy-dependent\"",
                ))
            }
            (FreePathChoice::EnergyDependent, None) => {
                return Err(invalid(
                    "physics.min_cm_angle_deg",
                    "required with free_path = \"energy-dependent\"",
                ))
            }
            (FreePathChoice::EnergyDependent, Some(a)) => {
                if !(a > 0.0 && a < 180.0) {
                    return Err(invalid(
                        "physics.min_cm_angle_deg",
                        format!("{a} must be in (0, 180) degrees"),
                    ));
                }
                MeanFreePath::EnergyDependent {
                    min_cm_angle_rad: a * PI / 180.0,
                }
            }
        };
        for (key, v) in [
            ("primary_cutoff_ev", p.primary_cutoff_ev),
            ("recoil_cutoff_ev", p.recoil_cutoff_ev),
        ] {
            if !finite_pos(v) {
                return Err(invalid(
                    format!("physics.{key}"),
                    format!("{v} must be finite and positive"),
                ));
            }
        }
        if !(p.primary_surface_binding_ev.is_finite() && p.primary_surface_binding_ev >= 0.0) {
            return Err(invalid(
                "physics.primary_surface_binding_ev",
                "must be finite and non-negative",
            ));
        }
        if p.primary_cutoff_ev >= b.energy_ev {
            return Err(invalid(
                "physics.primary_cutoff_ev",
                format!(
                    "{} eV is not below the beam energy {} eV",
                    p.primary_cutoff_ev, b.energy_ev
                ),
            ));
        }
        if p.follow_recoils {
            let min_es = layers
                .iter()
                .flat_map(|l| {
                    l.material
                        .components()
                        .iter()
                        .map(|c| l.material.surface_binding_energy_ev(c.z()).expect("set"))
                })
                .fold(f64::INFINITY, f64::min);
            if p.recoil_cutoff_ev > min_es {
                warnings.push(format!(
                    "physics.recoil_cutoff_ev = {} eV is above the smallest surface binding \
                     energy ({min_es} eV); sputtering will be underestimated",
                    p.recoil_cutoff_ev
                ));
            }
        }
        if p.weak_collisions > MAX_WEAK_COLLISIONS {
            return Err(invalid(
                "physics.weak_collisions",
                format!(
                    "{} must be at most {MAX_WEAK_COLLISIONS}",
                    p.weak_collisions
                ),
            ));
        }
        if p.weak_collisions > 0 && mean_free_path != MeanFreePath::Constant {
            return Err(invalid(
                "physics.weak_collisions",
                "only used with free_path = \"constant\"",
            ));
        }
        let mut config = BcaConfig::new(p.primary_cutoff_ev, p.recoil_cutoff_ev);
        config.mean_free_path = mean_free_path;
        config.weak_collisions = p.weak_collisions;
        config.electronic = match p.stopping {
            StoppingChoice::EquipartitionLsOr => ElectronicLoss::EquipartitionLsOr,
            _ => ElectronicLoss::NonLocal,
        };
        config.follow_recoils = p.follow_recoils;
        config.primary_surface_binding_ev = p.primary_surface_binding_ev;
        config.seed = self.run.seed;

        // Tally.
        if !finite_pos(self.tally.depth_bin_nm) {
            return Err(invalid("tally.depth_bin_nm", "must be finite and positive"));
        }
        if self.tally.depth_bins == 0 {
            return Err(invalid("tally.depth_bins", "must be at least 1"));
        }
        if !finite_pos(self.tally.lateral_bin_nm) {
            return Err(invalid(
                "tally.lateral_bin_nm",
                "must be finite and positive",
            ));
        }
        for (key, n) in [
            ("lateral_bins", self.tally.lateral_bins),
            ("escape_energy_bins", self.tally.escape_energy_bins),
            ("escape_polar_bins", self.tally.escape_polar_bins),
        ] {
            if n == 0 {
                return Err(invalid(format!("tally.{key}"), "must be at least 1"));
            }
        }
        if let Some(e) = self.tally.escape_energy_max_ev {
            if !finite_pos(e) {
                return Err(invalid(
                    "tally.escape_energy_max_ev",
                    "must be finite and positive",
                ));
            }
        }

        let screening = p.potential.screening();
        let screening_length = p
            .screening_length
            .map_or(screening.default_length(), LengthChoice::length);
        let stopping_tables = self.load_stopping_tables(base_dir)?;
        let fallback = stopping_model(p.stopping);
        let validity = fallback.validity(&ion);
        // The advisory range of the fallback model only matters for the pairs
        // it actually serves; a beam pair with a table has its range checked
        // as an error above.
        let beam_pair_uncovered = layers.iter().any(|l| {
            l.material.components().iter().any(|c| {
                !stopping_tables
                    .iter()
                    .any(|t| t.table.covers_pair(ion.z(), c.z()))
            })
        });
        if beam_pair_uncovered && !validity.contains(b.energy_ev) {
            warnings.push(format!(
                "beam energy {} eV is outside the advisory validity range of the {} \
                 stopping model for {} ([{}, {}] eV)",
                b.energy_ev,
                fallback.name(),
                b.ion,
                validity.min_energy_ev,
                validity.max_energy_ev
            ));
        }
        self.check_stopping_tables(&stopping_tables, &ion, &layers, &mut warnings)?;

        Ok(Resolved {
            input: self.echo(),
            beam,
            stack,
            layers,
            config,
            screening,
            screening_length,
            table_spec: TABLE_SPEC,
            stopping_tables,
            warnings,
            tuning,
        })
    }

    /// Apply `[physics] tuning` to the resolved layers (see [`TuningSet`]).
    fn apply_tuning(
        &self,
        sets: &[TuningSet],
        layers: &mut [ResolvedLayer],
        warnings: &mut Vec<String>,
    ) -> Result<Option<TuningReport>, InputError> {
        let name = self.physics.tuning.as_str();
        if name == NO_TUNING {
            return Ok(None);
        }
        let f = "physics.tuning";
        let set = sets.iter().find(|s| s.name == name).ok_or_else(|| {
            let known: Vec<&str> = sets.iter().map(|s| s.name).collect();
            invalid(
                f,
                format!("unknown tuning set {name:?}; use \"none\" or one of {known:?}"),
            )
        })?;
        if self.dynamic.is_some() {
            return Err(invalid(
                f,
                "tuning is not supported with a [dynamic] target (pilot: static single-element layers)",
            ));
        }
        let (e_lo, e_hi) = set.energy_range_ev;
        if !(e_lo.is_finite() && e_hi.is_finite() && 0.0 < e_lo && e_lo <= e_hi) {
            return Err(invalid(
                f,
                format!(
                    "set {name:?}: energy range {e_lo}..{e_hi} eV is not a finite positive range"
                ),
            ));
        }
        if !set.ions.contains(&self.beam.ion.as_str()) {
            return Err(invalid(
                f,
                format!(
                    "set {name:?} was fitted for beam ion(s) {:?}, not {:?}; use \"none\"",
                    set.ions, self.beam.ion
                ),
            ));
        }
        let e = self.beam.energy_ev;
        if !(e_lo..=e_hi).contains(&e) {
            warnings.push(format!(
                "tuning set {name:?} was fitted at {e_lo}..{e_hi} eV; {e} eV is an extrapolation"
            ));
        }
        if self.beam.tilt_deg != 0.0 {
            warnings.push(format!(
                "tuning set {name:?} was fitted at normal incidence; tilt {} deg is an extrapolation",
                self.beam.tilt_deg
            ));
        }
        let mut factors = BTreeMap::new();
        for &(sym, k) in set.e_s_factors {
            let el = element_by_symbol(sym)
                .ok_or_else(|| invalid(f, format!("set {name:?} names unknown element {sym:?}")))?;
            if !(k.is_finite() && k > 0.0) {
                return Err(invalid(
                    f,
                    format!("set {name:?}: factor for {sym} is {k}; must be finite and positive"),
                ));
            }
            if factors.insert(el.z, k).is_some() {
                return Err(invalid(f, format!("set {name:?} lists {sym} twice")));
            }
        }
        for (i, l) in layers.iter().enumerate() {
            if l.material.components().len() != 1 {
                return Err(invalid(
                    f,
                    format!(
                        "tuning supports single-element layers only; layer {i} ({}) is a compound",
                        l.source
                    ),
                ));
            }
        }
        let mut components = Vec::with_capacity(layers.len());
        for (i, l) in layers.iter_mut().enumerate() {
            let z = l.material.components()[0].z();
            let sym = crate::elements::element(z).expect("validated").symbol;
            let original = l.material.surface_binding_energy_ev(z).expect("set");
            let factor = *factors.get(&z).ok_or_else(|| {
                invalid(
                    f,
                    format!(
                        "set {name:?} has no factor for {sym} (layer {i}); it supports {:?}",
                        set.e_s_factors.iter().map(|p| p.0).collect::<Vec<_>>()
                    ),
                )
            })?;
            let effective = original * factor;
            l.material
                .set_surface_binding_energy_ev(z, effective)
                .map_err(|e| invalid(f, format!("layer {i}: effective E_s of {sym}: {e}")))?;
            components.push(TunedEnergy {
                layer: i,
                element: sym.to_string(),
                e_s_original_ev: original,
                factor,
                e_s_effective_ev: effective,
            });
        }
        if components.iter().all(|c| c.factor == 1.0) {
            warnings.push(format!(
                "tuning set {name:?} has factor 1 for every element in the target; nothing changed"
            ));
        }
        Ok(Some(TuningReport {
            set: set.name.to_string(),
            version: set.version,
            provenance: set.provenance.to_string(),
            quantity: "surface-binding-energy",
            components,
        }))
    }
}

/// The electronic stopping model for a choice. For
/// [`StoppingChoice::EquipartitionLsOr`] the engine does not use this model
/// (it carries its own Lindhard-Scharff/Oen-Robinson mix); Lindhard-Scharff
/// is returned so the validity range can still be checked.
pub fn stopping_model(choice: StoppingChoice) -> Box<dyn ElectronicStopping + Send + Sync> {
    match choice {
        StoppingChoice::LindhardScharff | StoppingChoice::EquipartitionLsOr => {
            Box::new(LindhardScharff::new())
        }
        StoppingChoice::BetheBloch => Box::new(BetheBloch::new()),
    }
}

impl Resolved {
    /// The electronic stopping model of the run: the `[physics] stopping`
    /// choice, with any `[stopping]` tables layered over it for the pairs
    /// they declare. This is what to pass to the engine.
    pub fn stopping_model(&self) -> Box<dyn ElectronicStopping + Send + Sync> {
        let base = stopping_model(self.input.physics.stopping);
        if self.stopping_tables.is_empty() {
            base
        } else {
            Box::new(TableOverride::new(
                self.stopping_tables
                    .iter()
                    .map(|l| l.table.clone())
                    .collect(),
                base,
            ))
        }
    }

    /// Every model in use, with its published source, in a fixed order.
    pub fn models(&self) -> Vec<ModelInfo> {
        let mut v = vec![ModelInfo {
            role: "transport",
            name: "amorphous-bca",
            citation: Cow::Borrowed(
                "J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980) 257; \
                       M. T. Robinson and I. M. Torrens, Phys. Rev. B 9 (1974) 5008; \
                       W. Eckstein, Computer Simulation of Ion-Solid Interactions (Springer, 1991)",
            ),
        }];
        v.push(match self.screening {
            Screening::ZblUniversal => ModelInfo {
                role: "screening function",
                name: "zbl-universal",
                citation: Cow::Borrowed(
                    "J. F. Ziegler, J. P. Biersack, U. Littmark, The Stopping and Range \
                           of Ions in Solids (Pergamon, 1985), ch. 2",
                ),
            },
            Screening::KrC => ModelInfo {
                role: "screening function",
                name: "kr-c",
                citation: Cow::Borrowed(
                    "W. D. Wilson, L. G. Haggmark, J. P. Biersack, Phys. Rev. B 15 (1977) 2458",
                ),
            },
            Screening::Moliere => ModelInfo {
                role: "screening function",
                name: "moliere",
                citation: Cow::Borrowed("G. Moliere, Z. Naturforsch. A 2 (1947) 133"),
            },
            Screening::LenzJensen => ModelInfo {
                role: "screening function",
                name: "lenz-jensen",
                citation: Cow::Borrowed(
                    "W. Lenz, Z. Phys. 77 (1932) 713; H. Jensen, Z. Phys. 77 (1932) 722",
                ),
            },
        });
        v.push(match self.screening_length {
            ScreeningLength::Universal => ModelInfo {
                role: "screening length",
                name: "universal",
                citation: Cow::Borrowed(
                    "J. F. Ziegler, J. P. Biersack, U. Littmark, The Stopping and Range \
                           of Ions in Solids (Pergamon, 1985), ch. 2",
                ),
            },
            ScreeningLength::Firsov => ModelInfo {
                role: "screening length",
                name: "firsov",
                citation: Cow::Borrowed("O. B. Firsov, Sov. Phys. JETP 6 (1958) 534"),
            },
            ScreeningLength::Lindhard => ModelInfo {
                role: "screening length",
                name: "lindhard",
                citation: Cow::Borrowed(
                    "J. Lindhard, M. Scharff, H. E. Schiott, Mat. Fys. Medd. Dan. Vid. \
                           Selsk. 33 (14) (1963)",
                ),
            },
        });
        v.push(ModelInfo {
            role: "scattering angle",
            name: "gauss-mehler-quadrature-table",
            citation: Cow::Borrowed(
                "scattering integral by Gauss-Mehler quadrature, tabulated in \
                       (reduced energy, reduced impact parameter); see lindhard::ion::scattering",
            ),
        });
        v.extend(match self.input.physics.stopping {
            StoppingChoice::LindhardScharff => vec![ModelInfo {
                role: "electronic stopping",
                name: "lindhard-scharff",
                citation: Cow::Borrowed("J. Lindhard and M. Scharff, Phys. Rev. 124 (1961) 128"),
            }],
            StoppingChoice::BetheBloch => vec![ModelInfo {
                role: "electronic stopping",
                name: "bethe-bloch",
                citation: Cow::Borrowed(
                    "H. Bethe, Ann. Phys. 5 (1930) 325; F. Bloch, Ann. Phys. 408 (1933) 285; \
                           mean excitation energy by the Bloch rule I = 10 eV Z2",
                ),
            }],
            StoppingChoice::EquipartitionLsOr => vec![
                ModelInfo {
                    role: "electronic stopping",
                    name: "lindhard-scharff",
                    citation: Cow::Borrowed(
                        "J. Lindhard and M. Scharff, Phys. Rev. 124 (1961) 128",
                    ),
                },
                ModelInfo {
                    role: "electronic loss partition",
                    name: "equipartition-ls-or",
                    citation: Cow::Borrowed(
                        "half nonlocal Lindhard-Scharff, half local Oen-Robinson: \
                               O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132 (1976) 647",
                    ),
                },
            ],
        });
        for l in &self.stopping_tables {
            v.push(ModelInfo {
                role: "electronic stopping (user table)",
                name: "user-table",
                citation: Cow::Owned(format!("{}: {}", l.path, l.table.provenance())),
            });
        }
        v.push(ModelInfo {
            role: "compound stopping",
            name: "bragg-additivity",
            citation: Cow::Borrowed(
                "W. H. Bragg and R. Kleeman, Phil. Mag. 10 (1905) 318; no compound correction",
            ),
        });
        v.push(match self.config.mean_free_path {
            MeanFreePath::Constant => ModelInfo {
                role: "free path",
                name: "constant",
                citation: Cow::Borrowed(
                    "J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980) 257",
                ),
            },
            MeanFreePath::EnergyDependent { .. } => ModelInfo {
                role: "free path",
                name: "energy-dependent",
                citation: Cow::Borrowed(
                    "W. Eckstein, Computer Simulation of Ion-Solid Interactions \
                           (Springer, 1991)",
                ),
            },
        });
        v.push(ModelInfo {
            role: "displacement criterion",
            name: "e_d-e_b",
            citation: Cow::Borrowed(
                "J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980) 257; \
                       W. Eckstein (1991)",
            ),
        });
        v.push(ModelInfo {
            role: "surface barrier",
            name: "planar",
            citation: Cow::Borrowed(
                "W. Eckstein, Computer Simulation of Ion-Solid Interactions (Springer, 1991)",
            ),
        });
        v.push(ModelInfo {
            role: "random numbers",
            name: "chacha8-per-history-stream",
            citation: Cow::Borrowed(
                "D. J. Bernstein, ChaCha, a variant of Salsa20 (2008); \
                       J. K. Salmon et al., Proc. SC'11 (2011)",
            ),
        });
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const B_SI: &str = r#"
[beam]
ion = "B"
energy_ev = 5000.0
tilt_deg = 7.0

[target]
substrate = "Si"

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 10
seed = 1
"#;

    fn err(text: &str) -> InputError {
        Input::from_toml_str(text)
            .and_then(|i| i.resolve().map(|_| ()))
            .unwrap_err()
    }

    fn field(e: &InputError) -> &str {
        match e {
            InputError::Invalid { field, .. } => field,
            InputError::Parse(_) => panic!("expected a value error, got {e}"),
        }
    }

    #[test]
    fn minimal_input_resolves_with_defaults() {
        let i = Input::from_toml_str(B_SI).unwrap();
        let r = i.resolve().unwrap();
        assert_eq!(r.beam.ion.z(), 5);
        assert!((r.beam.polar_rad - 7f64.to_radians()).abs() < 1e-15);
        assert_eq!(r.stack.layers().len(), 1);
        assert_eq!(r.config.mean_free_path, MeanFreePath::Constant);
        assert_eq!(r.config.electronic, ElectronicLoss::NonLocal);
        assert_eq!(r.config.seed, 1);
        assert_eq!(r.screening, Screening::ZblUniversal);
        assert_eq!(r.screening_length, ScreeningLength::Universal);
        assert_eq!(
            r.layers[0]
                .material
                .displacement_energy_ev(14)
                .unwrap()
                .to_bits(),
            15f64.to_bits()
        );
        assert_eq!(
            r.input.physics.screening_length,
            Some(LengthChoice::Universal)
        );
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    }

    #[test]
    fn echo_round_trips_and_drops_threads() {
        let mut i = Input::from_toml_str(B_SI).unwrap();
        i.run.threads = Some(3);
        let e = i.echo();
        assert_eq!(e.run.threads, None);
        let text = toml::to_string(&e).unwrap();
        assert_eq!(Input::from_toml_str(&text).unwrap(), e);
    }

    #[test]
    fn layers_by_name_and_inline() {
        let text = r#"
[beam]
ion = "As"
energy_ev = 5.0e4

[materials.SiO2]
density_g_cm3 = 2.2
elements = [
  { symbol = "Si", atom_fraction = 1.0 },
  { symbol = "O", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },
]

[target]
substrate = "Si"

[[target.layers]]
material = "SiO2"
thickness_nm = 10.0

[[target.layers]]
material = { elements = [ { symbol = "Si", atom_fraction = 1.0 } ] }
thickness_nm = 5.0

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 1.0
free_path = "energy-dependent"
min_cm_angle_deg = 0.5
stopping = "equipartition-ls-or"
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 10
seed = 2
"#;
        let r = Input::from_toml_str(text).unwrap().resolve().unwrap();
        assert_eq!(r.stack.layers().len(), 3);
        assert_eq!(r.layers[0].source, "SiO2");
        assert_eq!(r.layers[1].source, "inline");
        assert_eq!(r.layers[2].source, "Si");
        // The override reached Si in every layer.
        for l in &r.layers {
            assert_eq!(l.material.displacement_energy_ev(14).unwrap(), 15.0);
        }
        assert_eq!(r.config.electronic, ElectronicLoss::EquipartitionLsOr);
        assert!(matches!(
            r.config.mean_free_path,
            MeanFreePath::EnergyDependent { .. }
        ));
        assert!((r.stack.layers()[1].back_m() - 15e-9).abs() < 1e-20);
        // Echo of a mixed input round-trips through TOML.
        let back = Input::from_toml_str(&toml::to_string(&r.input).unwrap()).unwrap();
        assert_eq!(back, r.input);
    }

    #[test]
    fn weak_collisions_reach_the_engine_and_echo() {
        let r = Input::from_toml_str(&B_SI.replace(
            "primary_cutoff_ev = 5.0",
            "primary_cutoff_ev = 5.0\nweak_collisions = 3",
        ))
        .unwrap()
        .resolve()
        .unwrap();
        assert_eq!(r.config.weak_collisions, 3);
        let echo = toml::to_string(&r.input).unwrap();
        assert!(echo.contains("weak_collisions = 3"), "{echo}");
        // The default is 0 and is not echoed, so older echoes are unchanged.
        let r0 = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
        assert_eq!(r0.config.weak_collisions, 0);
        assert!(!toml::to_string(&r0.input)
            .unwrap()
            .contains("weak_collisions"));
    }

    #[test]
    fn unknown_key_is_named() {
        let e = err(&B_SI.replace("tilt_deg", "tilt_degrees"));
        assert!(matches!(e, InputError::Parse(_)));
        assert!(e.to_string().contains("tilt_degrees"), "{e}");
    }

    #[test]
    fn missing_displacement_energy_is_named() {
        let e = err(&B_SI.replace("[physics.energies.Si]\ne_d_ev = 15.0\n", ""));
        assert_eq!(field(&e), "target.substrate");
        assert!(e.to_string().contains("e_d_ev"), "{e}");
        assert!(e.to_string().contains("physics.energies.Si"), "{e}");
    }

    #[test]
    fn negative_thickness_is_named() {
        let text = B_SI.replace(
            "[target]\nsubstrate = \"Si\"\n",
            "[[target.layers]]\nmaterial = \"Si\"\nthickness_nm = -3.0\n",
        );
        let e = err(&text);
        assert_eq!(field(&e), "target.layers[0].thickness_nm");
    }

    #[test]
    fn value_errors_name_their_field() {
        let cases = [
            (B_SI.replace("ion = \"B\"", "ion = \"b\""), "beam.ion"),
            (B_SI.replace("5000.0", "-1.0"), "beam.energy_ev"),
            (B_SI.replace("7.0", "90.0"), "beam.tilt_deg"),
            (B_SI.replace("ions = 10", "ions = 0"), "run.ions"),
            (
                B_SI.replace("substrate = \"Si\"", "substrate = \"Unobtainium\""),
                "target.substrate",
            ),
            (
                B_SI.replace("substrate = \"Si\"", "substrate = \"O\""),
                "target.substrate",
            ),
            (
                B_SI.replace("recoil_cutoff_ev = 2.0", "recoil_cutoff_ev = 0.0"),
                "physics.recoil_cutoff_ev",
            ),
            (
                B_SI.replace(
                    "primary_cutoff_ev = 5.0",
                    "primary_cutoff_ev = 5.0\nmin_cm_angle_deg = 1.0",
                ),
                "physics.min_cm_angle_deg",
            ),
            (
                B_SI.replace(
                    "primary_cutoff_ev = 5.0",
                    "primary_cutoff_ev = 5.0\nfree_path = \"energy-dependent\"",
                ),
                "physics.min_cm_angle_deg",
            ),
            (
                B_SI.replace(
                    "primary_cutoff_ev = 5.0",
                    "primary_cutoff_ev = 5.0\nweak_collisions = 4",
                ),
                "physics.weak_collisions",
            ),
            (
                B_SI.replace(
                    "primary_cutoff_ev = 5.0",
                    "primary_cutoff_ev = 5.0\nweak_collisions = 2\nfree_path = \"energy-dependent\"\nmin_cm_angle_deg = 1.0",
                ),
                "physics.weak_collisions",
            ),
            (
                B_SI.replace("[physics.energies.Si]", "[physics.energies.Ge]"),
                "physics.energies.Ge",
            ),
            (
                B_SI.replace("e_d_ev = 15.0", "e_d_ev = -15.0"),
                "physics.energies.Si.e_d_ev",
            ),
            (
                format!("{B_SI}\n[tally]\ndepth_bin_nm = 0.0\n"),
                "tally.depth_bin_nm",
            ),
            (
                format!("{B_SI}\n[tally]\nlateral_bin_nm = -1.0\n"),
                "tally.lateral_bin_nm",
            ),
            (
                format!("{B_SI}\n[tally]\nlateral_bins = 0\n"),
                "tally.lateral_bins",
            ),
            (
                format!("{B_SI}\n[tally]\nescape_energy_bins = 0\n"),
                "tally.escape_energy_bins",
            ),
            (
                format!("{B_SI}\n[tally]\nescape_polar_bins = 0\n"),
                "tally.escape_polar_bins",
            ),
            (
                format!("{B_SI}\n[tally]\nescape_energy_max_ev = 0.0\n"),
                "tally.escape_energy_max_ev",
            ),
        ];
        for (text, want) in cases {
            let e = err(&text);
            assert_eq!(field(&e), want, "{e}");
        }
    }

    #[test]
    fn unknown_model_name_is_a_parse_error() {
        let e = err(&B_SI.replace(
            "primary_cutoff_ev = 5.0",
            "primary_cutoff_ev = 5.0\npotential = \"zbl2\"",
        ));
        assert!(matches!(e, InputError::Parse(_)));
        assert!(e.to_string().contains("zbl2"), "{e}");
    }

    #[test]
    fn warns_outside_stopping_validity() {
        let text = B_SI.replace(
            "primary_cutoff_ev = 5.0",
            "primary_cutoff_ev = 5.0\nstopping = \"bethe-bloch\"",
        );
        let r = Input::from_toml_str(&text).unwrap().resolve().unwrap();
        assert!(
            r.warnings.iter().any(|w| w.contains("bethe-bloch")),
            "{:?}",
            r.warnings
        );
    }

    #[test]
    fn models_list_every_choice() {
        let r = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
        let names: Vec<_> = r.models().iter().map(|m| m.name).collect();
        for n in ["zbl-universal", "universal", "lindhard-scharff", "constant"] {
            assert!(names.contains(&n), "{names:?}");
        }
        assert!(r.models().iter().all(|m| !m.citation.is_empty()));
    }

    // Tuning plumbing. The sets below are synthetic fixtures, not fitted
    // coefficients.
    const FIX: TuningSet = TuningSet {
        name: "fixture",
        version: 3,
        ions: &["B"],
        energy_range_ev: (1000.0, 10000.0),
        e_s_factors: &[("Si", 1.25), ("Ag", 0.5)],
        provenance: "synthetic test fixture",
    };

    fn with_tuning(text: &str, name: &str) -> Input {
        let mut i = Input::from_toml_str(text).unwrap();
        i.physics.tuning = name.to_string();
        i
    }

    fn resolve_fix(i: &Input, sets: &[TuningSet]) -> Result<Resolved, InputError> {
        i.resolve_with_sets(Path::new("."), sets)
    }

    fn es(r: &Resolved, layer: usize, z: u8) -> f64 {
        r.layers[layer]
            .material
            .surface_binding_energy_ev(z)
            .unwrap()
    }

    #[test]
    fn tuning_none_and_omitted_are_identical() {
        let a = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
        let b = with_tuning(B_SI, "none").resolve().unwrap();
        assert!(a.tuning.is_none() && b.tuning.is_none());
        assert_eq!(es(&a, 0, 14), es(&b, 0, 14));
        assert_eq!(
            toml::to_string(&a.input).unwrap(),
            toml::to_string(&b.input).unwrap()
        );
        assert!(!toml::to_string(&a.input).unwrap().contains("tuning"));
    }

    #[test]
    fn tuning_multiplies_resolved_default_once() {
        let base = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
        let r = resolve_fix(&with_tuning(B_SI, "fixture"), &[FIX]).unwrap();
        let t = r.tuning.as_ref().unwrap();
        assert_eq!((t.set.as_str(), t.version), ("fixture", 3));
        assert_eq!(t.components.len(), 1);
        let c = &t.components[0];
        assert_eq!(c.e_s_original_ev, es(&base, 0, 14));
        assert_eq!(c.factor, 1.25);
        assert_eq!(c.e_s_effective_ev, es(&r, 0, 14));
        assert_eq!(c.e_s_effective_ev, c.e_s_original_ev * 1.25);
        assert_eq!(
            toml::to_string(&r.input)
                .unwrap()
                .matches("fixture")
                .count(),
            1
        );
    }

    #[test]
    fn tuning_acts_on_explicit_override_and_each_layer() {
        let text = B_SI.replace("e_d_ev = 15.0", "e_d_ev = 15.0\ne_s_ev = 4.0").replace(
            "[target]\nsubstrate = \"Si\"\n",
            "[target]\nsubstrate = \"Si\"\n[[target.layers]]\nmaterial = \"Si\"\nthickness_nm = 5.0\n",
        );
        let r = resolve_fix(&with_tuning(&text, "fixture"), &[FIX]).unwrap();
        let t = r.tuning.unwrap();
        assert_eq!(t.components.len(), 2);
        for c in &t.components {
            assert_eq!(c.e_s_original_ev, 4.0);
            assert_eq!(c.e_s_effective_ev, 5.0);
        }
    }

    #[test]
    fn tuning_rejections() {
        let fe = |i: &Input, sets: &[TuningSet]| resolve_fix(i, sets).unwrap_err();
        let e = fe(&with_tuning(B_SI, "nope"), &[FIX]);
        assert_eq!(field(&e), "physics.tuning");
        assert!(e.to_string().contains("unknown tuning set"), "{e}");
        // A fixture name is unknown to the shipped registry.
        assert_eq!(
            field(&err(&B_SI.replace("[run]", "[run]\n").replace(
                "recoil_cutoff_ev = 2.0",
                "recoil_cutoff_ev = 2.0\ntuning = \"fixture\""
            ))),
            "physics.tuning"
        );
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let factors: &'static [(&str, f64)] = Box::leak(Box::new([("Si", bad)]));
            let s = TuningSet {
                e_s_factors: factors,
                ..FIX
            };
            assert_eq!(
                field(&fe(&with_tuning(B_SI, "fixture"), &[s])),
                "physics.tuning"
            );
        }
        let unknown_el = TuningSet {
            e_s_factors: &[("Xx", 1.0)],
            ..FIX
        };
        assert!(fe(&with_tuning(B_SI, "fixture"), &[unknown_el])
            .to_string()
            .contains("Xx"));
        let dup = TuningSet {
            e_s_factors: &[("Si", 1.0), ("Si", 2.0)],
            ..FIX
        };
        assert!(fe(&with_tuning(B_SI, "fixture"), &[dup])
            .to_string()
            .contains("twice"));
        // Compound layer.
        let compound = B_SI.replace(
            "[target]\nsubstrate = \"Si\"\n",
            "[materials.SiO2]\ndensity_g_cm3 = 2.2\nelements = [\n  { symbol = \"Si\", atom_fraction = 1.0 },\n  { symbol = \"O\", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },\n]\n[target]\nsubstrate = \"SiO2\"\n",
        );
        let e = fe(&with_tuning(&compound, "fixture"), &[FIX]);
        assert!(e.to_string().contains("single-element"), "{e}");
        // Unset energy stays an error, not a default.
        let none_set = TuningSet {
            e_s_factors: &[("Si", 1.0)],
            ..FIX
        };
        let r = resolve_fix(&with_tuning(B_SI, "fixture"), &[none_set]).unwrap();
        assert!(r.warnings.iter().any(|w| w.contains("nothing changed")));
        // A beam the set was not fitted for.
        let e = fe(
            &with_tuning(&B_SI.replace("ion = \"B\"", "ion = \"P\""), "fixture"),
            &[FIX],
        );
        assert_eq!(field(&e), "physics.tuning");
        assert!(e.to_string().contains("beam ion"), "{e}");
        // A target element the set does not list.
        let ge = B_SI
            .replace("substrate = \"Si\"", "substrate = \"Ge\"")
            .replace("[physics.energies.Si]", "[physics.energies.Ge]");
        let e = fe(&with_tuning(&ge, "fixture"), &[FIX]);
        assert!(e.to_string().contains("no factor for Ge"), "{e}");
        // A malformed energy range.
        for range in [(0.0, 1.0), (2.0, 1.0), (1.0, f64::INFINITY)] {
            let s = TuningSet {
                energy_range_ev: range,
                ..FIX
            };
            assert_eq!(
                field(&fe(&with_tuning(B_SI, "fixture"), &[s])),
                "physics.tuning"
            );
        }
    }

    #[test]
    fn tuning_warns_outside_the_fitted_domain() {
        // B_SI: 5 keV (inside 1..10 keV) at 7 deg tilt.
        let r = resolve_fix(&with_tuning(B_SI, "fixture"), &[FIX]).unwrap();
        assert!(r.warnings.iter().any(|w| w.contains("normal incidence")));
        assert!(!r.warnings.iter().any(|w| w.contains("fitted at 1000")));
        let low = B_SI
            .replace("energy_ev = 5000.0", "energy_ev = 500.0")
            .replace("tilt_deg = 7.0", "tilt_deg = 0.0");
        let r = resolve_fix(&with_tuning(&low, "fixture"), &[FIX]).unwrap();
        assert!(r.warnings.iter().any(|w| w.contains("extrapolation")));
        assert!(!r.warnings.iter().any(|w| w.contains("normal incidence")));
        assert_eq!(r.tuning.unwrap().components[0].factor, 1.25);
    }

    #[test]
    fn shipped_sets_are_well_formed() {
        let mut names = std::collections::BTreeSet::new();
        for s in TUNING_SETS {
            assert!(names.insert(s.name), "duplicate set name {}", s.name);
            assert_ne!(s.name, NO_TUNING);
            assert!(!s.ions.is_empty() && !s.e_s_factors.is_empty());
            assert!(!s.provenance.is_empty());
            for &(sym, k) in s.e_s_factors {
                assert!(element_by_symbol(sym).is_some(), "{sym}");
                assert!(k.is_finite() && k > 0.0, "{sym}: {k}");
            }
            // Every listed element resolves for every listed ion, once.
            for ion in s.ions {
                for &(sym, k) in s.e_s_factors {
                    let text = format!(
                        "[beam]\nion = \"{ion}\"\nenergy_ev = {}\n[target]\nsubstrate = \"{sym}\"\n\
                         [physics]\nprimary_cutoff_ev = 2.0\nrecoil_cutoff_ev = 1.0\ntuning = \"{}\"\n\
                         [physics.energies.{sym}]\ne_d_ev = 10.0\n\
                         [run]\nions = 1\nseed = 1\n",
                        s.energy_range_ev.0, s.name
                    );
                    let r = Input::from_toml_str(&text).unwrap().resolve().unwrap();
                    let c = &r.tuning.as_ref().unwrap().components[0];
                    assert_eq!(c.factor, k);
                    assert_eq!(c.e_s_effective_ev, c.e_s_original_ev * k);
                    assert_eq!(
                        es(&r, 0, element_by_symbol(sym).unwrap().z),
                        c.e_s_effective_ev
                    );
                }
            }
        }
    }
}
