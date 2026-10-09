//! Schema types: the serde structs and enums of the run description.
//!
//! Moved out of the parent module unchanged; every item is re-exported there.

use crate::ion::potential::{Screening, ScreeningLength};
use crate::ion::stopping::table::StoppingTable;
use crate::material::MaterialSpec;
use serde::de::{self, MapAccess, Visitor};
use serde::{Deserialize, Deserializer, Serialize};
use std::collections::BTreeMap;
use std::path::PathBuf;

// Only for intra-doc links.
#[cfg(doc)]
use crate::ion::bca::{ElectronicLoss, MeanFreePath};
#[cfg(doc)]
use super::Resolved;

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
    pub(super) fn from_length(l: ScreeningLength) -> Self {
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

pub(super) fn default_min_ions() -> u64 {
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
