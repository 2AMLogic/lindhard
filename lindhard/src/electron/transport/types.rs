//! The value types of the electron transport: the run configuration and its
//! error, the state of an electron and how its history ended, and the
//! metadata a run records. The engine that uses them is
//! [`Transport`](super::Transport); the hooks it reports through are in
//! [`super::tally`].

use serde::Serialize;

use crate::electron::beam::BeamSpread;
use crate::electron::boundary::BandStructure;
use crate::electron::data::{CrossSectionTable, SamplingAxis};
use crate::electron::phonon::{FrohlichPhonon, PolaronTrapping};
use crate::electron::secondary::SecondaryModel;
use crate::geometry::Face;

/// Errors from building a [`Transport`](super::Transport) or starting a run.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum TransportError {
    /// The number of table pairs differs from the number of layers.
    #[error("the stack has {layers} layers but {tables} table pairs were given")]
    LayerCount {
        /// Layers in the stack.
        layers: usize,
        /// Table pairs given.
        tables: usize,
    },
    /// A table has the wrong sampling axis for its slot.
    #[error("layer {layer}: the {slot} table samples {found:?}, expected {expected:?}")]
    WrongAxis {
        /// Layer index.
        layer: usize,
        /// `"elastic"` or `"inelastic"`.
        slot: &'static str,
        /// Axis found.
        found: SamplingAxis,
        /// Axis expected.
        expected: SamplingAxis,
    },
    /// A configuration or primary value is out of range.
    #[error("invalid {what}: {why}")]
    Invalid {
        /// The quantity.
        what: &'static str,
        /// Why it was rejected.
        why: String,
    },
}

/// An [`TransportError::Invalid`] for `what`, as an `Err`. Shared with
/// [`crate::electron::beam`], whose validators report the same error.
pub(crate) fn invalid<T>(what: &'static str, why: impl Into<String>) -> Result<T, TransportError> {
    Err(TransportError::Invalid {
        what,
        why: why.into(),
    })
}

/// What happens when an electron reaches a face of the target. Whether it can
/// get out at all is decided by the [`BoundaryModel`]; an electron that does
/// leave never returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum EscapeRule {
    /// Both the front face (`x = 0`) and the back face of a finite stack
    /// release the electron ([`Fate::Escaped`] with [`Face::Front`] or
    /// [`Face::Back`]).
    BothFaces,
    /// The front face releases the electron; the back face of a finite stack
    /// absorbs it ([`Fate::Absorbed`]).
    FrontOnly,
}

/// A face an electron met: between two layers, or a face of the target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// The face between layer `from` (where the electron came from) and
    /// layer `to`.
    Interface {
        /// Layer the electron was in.
        from: usize,
        /// Layer on the other side.
        to: usize,
    },
    /// A face of the target, with vacuum on the other side.
    Surface(Face),
}

/// What a face does to an electron.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(tag = "model", rename_all = "kebab-case")]
pub enum BoundaryModel {
    /// Faces change nothing (no inner potential, no reflection).
    Transparent,
    /// Each face is a step in the inner potential
    /// ([`crate::electron::boundary`]); every layer needs a
    /// [`BandStructure`].
    StepBarrier {
        /// Transmit with the quantum-mechanical probability `T` (Verduin
        /// Eq. 3.145); if false, every electron that can get over does.
        quantum_transmission: bool,
        /// Refract on transmission (Verduin Eq. 3.139); if false, keep the
        /// direction.
        refraction: bool,
    },
}

impl BoundaryModel {
    /// The step barrier with quantum transmission and refraction, Nebula's
    /// defaults (`boundary_intersect` template arguments at the commit named
    /// in [`crate::electron::boundary`]).
    pub const STEP_BARRIER: Self = Self::StepBarrier {
        quantum_transmission: true,
        refraction: true,
    };
}

/// Where the energy cutoff is measured from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum CutoffReference {
    /// From the band bottom of the electron's layer: an electron stops below
    /// `cutoff`.
    BandBottom,
    /// From the vacuum level: an electron stops below `U + cutoff`, `U` the
    /// inner potential of its layer (zero without a [`BandStructure`]), so
    /// electrons that could never leave are not followed. Nebula's rule.
    VacuumLevel,
}

/// How a history (or one electron of it) ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// The energy fell below the stopping threshold inside the target.
    Stopped,
    /// The electron left the target through a face.
    Escaped(Face),
    /// The electron reached the back face under [`EscapeRule::FrontOnly`].
    Absorbed,
    /// The electron has no interaction available and no face to reach (zero
    /// total rate heading parallel to the faces, or away from every face).
    Trapped,
    /// The electron hit [`TransportConfig::max_events`] and was cut off.
    EventCap,
    /// The electron was trapped as a polaron (an opt-in insulator channel,
    /// [`crate::electron::phonon::PolaronTrapping`]) and deposited its
    /// remaining energy there.
    PolaronTrapped,
}

/// Which way a Fröhlich LO-phonon event went.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhononEvent {
    /// The electron emitted a phonon and lost `ħω`.
    Emission,
    /// The electron absorbed a phonon and gained `ħω`.
    Absorption,
}

/// One of the two cross-section tables of a layer ([`LayerTables`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum TableChannel {
    /// The elastic table.
    Elastic,
    /// The inelastic table.
    Inelastic,
}

/// Where an energy lies against the energy grid of one table, as reported to
/// [`ElectronTally::table_lookup`](super::tally::ElectronTally::table_lookup).
///
/// The grid is the table's own [`CrossSectionTable::energy_ev`], from its
/// first value `E_min` to its last `E_max`. Inside `[E_min, E_max]` the rate
/// and the inverse CDF are interpolated between rows (an energy equal to an
/// endpoint reads that row exactly); outside it they are the endpoint row's,
/// held constant ([the module docs](super#model), "Free flight"). This enum says which of the
/// two happened; it measures how often that continuation is used, not
/// whether it is accurate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum GridCoverage {
    /// `E < E_min`: the first row's values were used.
    Below,
    /// `E_min <= E <= E_max` (both endpoints included).
    Within,
    /// `E > E_max`: the last row's values were used.
    Above,
}

impl GridCoverage {
    /// Where `energy_ev` lies against the grid of `table`.
    pub fn of(table: &CrossSectionTable, energy_ev: f64) -> Self {
        let g = table.energy_ev();
        if energy_ev < g[0] {
            Self::Below
        } else if energy_ev > g[g.len() - 1] {
            Self::Above
        } else {
            Self::Within
        }
    }
}

/// Run configuration, recorded in the [`RunMetadata`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TransportConfig {
    /// An electron stops when its energy falls below this (measured as
    /// [`TransportConfig::cutoff_reference`] says), eV.
    pub cutoff_ev: f64,
    /// What a face of the target does to an electron reaching it.
    pub escape_rule: EscapeRule,
    /// Safety cap on collisions and face reflections per electron (the
    /// primary and each secondary count separately). An electron that
    /// reaches it ends with [`Fate::EventCap`].
    pub max_events: u64,
    /// Secondary-electron generation at inelastic events.
    pub secondaries: SecondaryModel,
    /// What the faces do.
    pub boundary: BoundaryModel,
    /// Where the cutoff is measured from.
    pub cutoff_reference: CutoffReference,
}

impl TransportConfig {
    /// A configuration with the given cutoff (from the band bottom),
    /// [`EscapeRule::BothFaces`], a cap of 10 million events per electron, no
    /// secondaries and transparent faces.
    pub fn new(cutoff_ev: f64) -> Self {
        Self {
            cutoff_ev,
            escape_rule: EscapeRule::BothFaces,
            max_events: 10_000_000,
            secondaries: SecondaryModel::Off,
            boundary: BoundaryModel::Transparent,
            cutoff_reference: CutoffReference::BandBottom,
        }
    }

    pub(super) fn needs_bands(&self) -> bool {
        self.secondaries != SecondaryModel::Off || self.boundary != BoundaryModel::Transparent
    }
}

/// The tables of one layer.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerTables {
    /// Elastic table ([`SamplingAxis::ElasticPolarAngle`]).
    pub elastic: CrossSectionTable,
    /// Inelastic table ([`SamplingAxis::InelasticEnergyLoss`]).
    pub inelastic: CrossSectionTable,
}

/// A primary electron: it starts on the front face (`x = 0`, `y = z = 0`) in
/// the first layer, or, with [`BoundaryModel::StepBarrier`], just outside it
/// in vacuum (then `energy_ev` is the vacuum energy and the direction must
/// point into the target).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Primary {
    /// Kinetic energy, eV; must exceed the cutoff and, once inside, the
    /// stopping threshold of the first layer.
    pub energy_ev: f64,
    /// Direction of flight (any nonzero finite vector; normalised on use).
    /// Normal incidence is `[1, 0, 0]`; a negative `x` component starts the
    /// electron moving out of the target.
    pub direction: [f64; 3],
}

impl Primary {
    /// Normal incidence at `energy_ev`.
    pub fn normal(energy_ev: f64) -> Self {
        Self {
            energy_ev,
            direction: [1.0, 0.0, 0.0],
        }
    }
}

/// The state of an electron after an event. Positions are `[x, y, z]` in
/// metres, `x` the depth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectronState {
    /// Position, m.
    pub pos: [f64; 3],
    /// Unit direction of flight.
    pub dir: [f64; 3],
    /// Kinetic energy, eV.
    pub energy_ev: f64,
    /// Index of the layer the electron is in (or was last in, after an
    /// escape).
    pub layer: usize,
}

/// Identity of the tables used in one layer, for [`RunMetadata`].
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LayerMetadata {
    /// Layer index (0 = front).
    pub index: usize,
    /// Depth of the front face, m; `None` for a material of a geometry other
    /// than a [`Stack`](crate::geometry::Stack) (a voxel grid has no layer
    /// extent).
    pub front_m: Option<f64>,
    /// Depth of the back face, m; `None` for a semi-infinite substrate (or
    /// for a geometry other than a [`Stack`](crate::geometry::Stack)).
    pub back_m: Option<f64>,
    /// Elastic table model identity.
    pub elastic_model: String,
    /// Elastic table provenance.
    pub elastic_provenance: String,
    /// Inelastic table model identity.
    pub inelastic_model: String,
    /// Inelastic table provenance.
    pub inelastic_provenance: String,
    /// Band parameters of the layer, with their provenance, if given.
    pub band_structure: Option<BandStructure>,
    /// Fröhlich LO-phonon channel of this layer (its parameters and their
    /// provenance), or `None` when off (the default, and the choice for a
    /// metal).
    pub phonon: Option<FrohlichPhonon>,
    /// Polaron-trapping channel of this layer, or `None` when off.
    pub polaron: Option<PolaronTrapping>,
    /// Inner-shell channels of this layer, in sampling order
    /// ([`Transport::with_inner_shells`](super::Transport::with_inner_shells));
    /// empty (and not serialized) when none.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub inner_shells: Vec<InnerShellMetadata>,
}

/// Identity of one inner-shell channel table of a layer, for
/// [`LayerMetadata`].
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct InnerShellMetadata {
    /// Atomic number.
    pub z: u8,
    /// Subshell label (ENDF, e.g. `"K"`).
    pub subshell: String,
    /// Binding energy, eV.
    pub binding_energy_ev: f64,
    /// Where the binding energy came from.
    pub binding_provenance: String,
    /// The channel table's model identity.
    pub model: String,
    /// The channel table's provenance.
    pub provenance: String,
}

/// What a run was configured with, for output metadata.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RunMetadata {
    /// Energy cutoff, eV.
    pub cutoff_ev: f64,
    /// Where the cutoff is measured from.
    pub cutoff_reference: CutoffReference,
    /// Escape rule at the target faces.
    pub escape_rule: EscapeRule,
    /// Collision and reflection cap per electron.
    pub max_events: u64,
    /// Secondary-electron model.
    pub secondaries: SecondaryModel,
    /// Face model.
    pub boundary: BoundaryModel,
    /// Run seed.
    pub seed: u64,
    /// Primaries run.
    pub n_histories: u64,
    /// Chunk size of the parallel driver.
    pub chunk_size: u64,
    /// The primary electron.
    pub primary: Primary,
    /// Landing profile and divergence; left out for a point beam.
    #[serde(skip_serializing_if = "BeamSpread::is_point")]
    pub beam: BeamSpread,
    /// Per-layer table identities.
    pub layers: Vec<LayerMetadata>,
}

/// The result of [`Transport::run`](super::Transport::run).
#[derive(Debug, Clone, PartialEq)]
pub struct TransportRun<T> {
    /// The merged tally.
    pub tally: T,
    /// The run configuration.
    pub metadata: RunMetadata,
}
