//! Event-by-event electron transport in a layered stack, on precomputed
//! cross-section tables, with secondary electrons and surface barriers.
//!
//! This is the low-energy electron Monte Carlo loop of Kieft and Bosch,
//! J. Phys. D: Appl. Phys. 41, 215310 (2008), doi:10.1088/0022-3727/41/21/215310
//! (the scheme on which Nebula is built): an electron at `(r, direction, E)`
//! in a layer flies a free path drawn from the total inverse mean free path
//! (the sum of the elastic and inelastic [`CrossSectionTable`] rates of the
//! layer at `E`), then picks the elastic or the inelastic channel with
//! probability proportional to the two rates, samples the polar angle `θ` or
//! the energy loss `W` from that channel's inverse CDF, and updates its state.
//! An electron ends when `E` falls below the stopping threshold or it leaves
//! the target. The loop itself is written from the paper; the secondary
//! electron kinematics ([`crate::electron::secondary`]) and the potential step
//! ([`crate::electron::boundary`]) are ported from Nebula with attribution.
//!
//! # Model
//!
//! - **Free flight.** `s = -ln(1 - u) / (λ_el⁻¹ + λ_inel⁻¹)`, `u` uniform in
//!   `[0, 1)`. Rates are interpolated linearly in energy on each table's own
//!   grid and held constant beyond its ends.
//! - **Layer boundaries are crossed exactly.** If the drawn path is longer
//!   than the distance to the layer face along the flight direction, the
//!   electron moves to the face (its depth is set to the face value exactly),
//!   meets the face (below) and **redraws** the path with the rates of the
//!   layer it is then in. This is exact because the free path is memoryless.
//! - **Elastic collision.** `θ` from the table, azimuth uniform on `[0, 2π)`;
//!   the energy is unchanged (no recoil energy).
//! - **Inelastic collision.** `W` from the table (clamped to `E`); the energy
//!   becomes `E - W`. With [`SecondaryModel::Off`] the flight direction is
//!   unchanged and nothing else happens. With
//!   [`SecondaryModel::KieftBosch`] `W` is further clamped to `E - E_F` (the
//!   primary cannot end below the Fermi level, as in Nebula's loss table),
//!   and the event may liberate a secondary electron (energy `E_F + W - B`,
//!   direction by the Ivanchenko method) and deflect the primary; see
//!   [`crate::electron::secondary`]. The stopping threshold of every layer
//!   must then exceed its Fermi energy, or secondaries would multiply without
//!   end.
//! - **Insulator channels (opt-in per layer).** A layer given
//!   [`InsulatorChannels`] through [`Transport::with_insulator_channels`]
//!   adds up to three more rates to the total: Fröhlich LO-phonon emission
//!   and absorption, and polaron trapping (models and sources in
//!   [`crate::electron::phonon`]). A phonon event changes the energy by `∓ħω`
//!   and deflects the electron by `θ` drawn from the Fröhlich angular
//!   distribution, azimuth uniform. A trapped electron deposits all of its
//!   remaining energy where it is trapped and its history ends with
//!   [`Fate::PolaronTrapped`]. Layers default to no insulator channel, which
//!   is the right choice for a metal; the per-layer choice is in
//!   [`LayerMetadata`].
//! - **Channel choice.** With the rates in the fixed order elastic,
//!   inelastic, phonon emission, phonon absorption, polaron trapping, the
//!   channel is the first whose cumulative rate exceeds `u · total`; channels
//!   with zero rate are never chosen (a rounding overshoot falls to the last
//!   channel with a positive rate).
//! - **Faces.** With [`BoundaryModel::Transparent`] a face changes nothing:
//!   the electron enters the next layer or leaves the target at its energy
//!   and direction. With [`BoundaryModel::StepBarrier`] every face is a step
//!   in the inner potential (vacuum outside the target has `U = 0`), crossed
//!   with the quantum-mechanical transmission probability and refraction of
//!   [`crate::electron::boundary`], or reflected specularly; an electron
//!   below the barrier is always reflected. Energies inside a layer are then
//!   measured from that layer's band bottom and an escaped electron carries
//!   its vacuum energy `E - U`. The back face under [`EscapeRule::FrontOnly`]
//!   absorbs without a barrier.
//! - **Inverse CDFs between grid energies** are interpolated linearly in
//!   energy between the two bracketing rows, at the same cumulative
//!   probability (a zero-rate bracketing row is ignored).
//!
//! # Electrons of one history
//!
//! A history starts with its primary. With [`BoundaryModel::StepBarrier`] the
//! primary is incident from vacuum on the front face, its [`Primary`] energy
//! is a vacuum energy, and it first meets the front-face step (it may be
//! reflected there, ending the history as [`Fate::Escaped`] through
//! [`Face::Front`]). Secondaries are pushed on a last-in, first-out stack and
//! each is transported in full, like the primary, once the electron before it
//! has ended; secondaries of secondaries are pushed on the same stack. The
//! history ends when the stack is empty.
//!
//! An electron stops when its energy falls below the **stopping threshold**
//! of its layer: the cutoff, measured from the band bottom
//! ([`CutoffReference::BandBottom`]) or from the vacuum level
//! ([`CutoffReference::VacuumLevel`], threshold `U + cutoff`, Nebula's rule
//! in `source/drivers/cpu/cpu_driver.inl` at the commit named in
//! [`crate::electron::boundary`]). A secondary below the threshold of its
//! layer is not created; its energy is recorded as deposited.
//!
//! # Energy reference of the inelastic table
//!
//! The transport looks up both tables of a layer at the electron's kinetic
//! energy `E` as it measures it: from the band bottom of the layer (with the
//! step barrier, the vacuum energy plus the inner potential `U`). An
//! inelastic table built by [`crate::electron::inelastic::table`] has the
//! model's own axis instead: `T`, the energy above the *model's* Fermi
//! energy, with kinematics on `T' = T + E_F(model)` and losses up to `T`
//! (Shinotsuka et al. 2017, eqs. (2)-(3), as in
//! [`crate::electron::inelastic::penn`]). A table does not know the band.
//!
//! - With the model's Fermi energy at its default, 0, a row read at the
//!   band-bottom energy `E` has the right kinematics (`T' = E`) but losses
//!   up to `E` instead of `E - E_F`. Under [`SecondaryModel::KieftBosch`] the
//!   clamp above takes every loss beyond `E - E_F` to `E - E_F`: the primary
//!   is left at the Fermi level and the secondary takes its whole energy.
//!   This is how `lindhard run` uses the tables
//!   (`[electron.inelastic] fermi_energy_ev = 0`), and the tests
//!   `inelastic_table_is_read_at_the_band_bottom_energy` and
//!   `losses_beyond_the_fermi_level_are_clamped_to_it`
//!   (`tests/electron_secondaries.rs`) pin it.
//! - Setting the model's Fermi energy to the band's while still reading rows
//!   at `E` counts `E_F` twice (`T' = E + E_F`).
//! - The consistent convention, that of the cstool table compiler which builds
//!   Nebula's tables (rows on the band-bottom axis `K`, kinematics on `K`,
//!   losses below `K - E_F`; `compile_full_imfp_icdf` in
//!   `cstool/dielectric_function/compile.py`), would be rows at
//!   `T = E - E_F` with the model's Fermi energy set to the band's. It is not
//!   what `lindhard run` builds.
//!
//! **Measured effect (#173).** With the single-pole Penn default, the Al
//! input of `validation/experiments/se_yield.py` (`E_F` = 11.66 eV) sends
//! 78 to 82 % of the inelastic events of electrons 5 to 20 eV above the Fermi
//! level into the clamp (Al 400 eV, 200 histories). Rebuilding the table in
//! the cstool convention removes every clamped event but moves δ only within
//! the run's noise: Al 6.89 to 6.68 (400 eV), 7.44 to 7.96 (800 eV); Au 3.43
//! to 3.63 and Cu 1.67 to 1.80 (800 eV). Setting
//! `fermi_energy_ev` to the band value raises Al δ (6.89 to 7.63 at 400 eV,
//! 7.44 to 9.01 at 800 eV). The reference mismatch is therefore real but is
//! not the cause of the δ overestimate of #173; that is the single-pole
//! model's inelastic mean free path at low energy
//! ([`crate::electron::inelastic::penn`], "Low energies"). The commands are
//! in `lindhard-cli/examples/inelastic_low_energy.rs`.
//!
//! # Random draw order
//!
//! The order of draws from the per-history stream is part of the contract (the
//! tests replay it). With the step barrier, the primary first draws one
//! uniform for its entry. Then, per flight with a nonzero total rate: one
//! uniform for the path. If a collision happens, one for the channel, then
//! `θ` and the azimuth (elastic), `W` (inelastic), or `cos θ` and the azimuth
//! (phonon emission or absorption) or nothing (polaron trapping). Under the
//! Kieft-Bosch secondary model an inelastic event that liberates an electron
//! then draws the secondary's azimuth and, with the instantaneous momentum on,
//! two more uniforms. With no insulator channel switched on, the draws are
//! exactly those of a run without them. A path drawn for a flight that ends at a layer face is discarded;
//! with the step barrier, the face draws one uniform if the electron has the
//! normal energy to get over (and none if it is surely reflected). A flight
//! with zero total rate draws nothing. Secondaries continue on the same
//! stream, in stack order, after the electron before them ends.
//!
//! # Reproducibility
//!
//! [`Transport::run`] goes through [`crate::rng::run_particles`]: history `i`
//! uses the stream `(seed, i)`, and chunk tallies are merged in chunk order,
//! so a tally whose [`ElectronTally::merge`] is order-faithful (a sum, or a
//! concatenation) is bit-identical at any thread count.
//!
//! # Hooks
//!
//! The [`ElectronTally`] trait mirrors [`crate::ion::bca::BcaTally`]: every
//! event hook has a no-op default. Per-electron hooks (`step`, `elastic`,
//! `inelastic`, `secondary`, `interface`, `barrier`, `reflected`, `stopped`,
//! `escaped`, `absorbed`, `phonon`, `polaron_trapped`) fire for the primary and
//! for every secondary;
//! [`ElectronTally::begin_secondary`] and [`ElectronTally::end_secondary`]
//! bracket each secondary, and [`ElectronTally::end_history`] comes last, with
//! the primary's fate. The insulator channels report through
//! [`ElectronTally::phonon`] and [`ElectronTally::polaron_trapped`].

use serde::Serialize;

use rand_core::Rng;

use crate::electron::boundary::{cross_step, BandStructure, StepOutcome};
use crate::electron::data::{CrossSectionTable, SamplingAxis};
use crate::electron::phonon::{
    sample_cos_theta, FrohlichPhonon, InsulatorChannels, PolaronTrapping,
};
use crate::electron::secondary::{kieft_bosch, SecondaryEvent, SecondaryModel};
use crate::geometry::Stack;
use crate::rng::run_particles;

/// Errors from building a [`Transport`] or starting a run.
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

fn invalid<T>(what: &'static str, why: impl Into<String>) -> Result<T, TransportError> {
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

/// Which face of the target an electron left through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Face {
    /// The front surface, `x = 0`.
    Front,
    /// The back face of a finite stack.
    Back,
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

    fn needs_bands(&self) -> bool {
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

/// Event hooks called by [`Transport`]. All default to no-ops except `merge`.
/// The `ElectronState` passed is the state **after** the event.
pub trait ElectronTally: Send {
    /// A primary history starts. With [`BoundaryModel::StepBarrier`], `start`
    /// is the primary in vacuum, before it meets the front face.
    fn begin_history(&mut self, _index: u64, _start: &ElectronState) {}

    /// A flight segment of `length_m` ended at `end.pos` (at a collision, a
    /// layer face or the target face). It started at `from` with energy
    /// `end.energy_ev` (the energy does not change in flight).
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, _length_m: f64) {}

    /// An elastic collision deflected the electron by `theta` rad.
    fn elastic(&mut self, _after: &ElectronState, _theta: f64) {}

    /// An inelastic collision took `w_ev` from the electron. `after` carries
    /// the new direction if the secondary model deflected it.
    fn inelastic(&mut self, _after: &ElectronState, _w_ev: f64) {}

    /// Called right after [`ElectronTally::inelastic`] for every inelastic
    /// event when a [`SecondaryModel`] is on: the event's energy bookkeeping,
    /// and the secondary that was created (`None` if none was). The
    /// secondary is transported later, between
    /// [`ElectronTally::begin_secondary`] and
    /// [`ElectronTally::end_secondary`].
    fn secondary(
        &mut self,
        _primary: &ElectronState,
        _event: &SecondaryEvent,
        _created: Option<&ElectronState>,
    ) {
    }

    /// A Fröhlich LO-phonon event changed the energy by `-ħω` (emission) or
    /// `+ħω` (absorption), `hbar_omega_ev` given, and deflected the electron by
    /// `theta` rad. `after` is the state after the event.
    fn phonon(
        &mut self,
        _after: &ElectronState,
        _kind: PhononEvent,
        _hbar_omega_ev: f64,
        _theta: f64,
    ) {
    }

    /// The electron was trapped as a polaron at `at` and deposits its whole
    /// remaining energy `at.energy_ev` there. The history (or the secondary)
    /// then ends with [`Fate::PolaronTrapped`] (no [`ElectronTally::stopped`]
    /// call).
    fn polaron_trapped(&mut self, _at: &ElectronState) {}

    /// The electron crossed from layer `from_layer` into `to_layer`; `at` is
    /// its state on the face, inside the new layer (after refraction, with
    /// [`BoundaryModel::StepBarrier`]).
    fn interface(&mut self, _at: &ElectronState, _from_layer: usize, _to_layer: usize) {}

    /// With [`BoundaryModel::StepBarrier`]: the electron got through the face
    /// `boundary`, and its kinetic energy changed by `delta_u_ev` (the step
    /// `U' - U`). `at` is its state just past the face (in vacuum, for a
    /// surface). Comes before [`ElectronTally::interface`] or
    /// [`ElectronTally::escaped`].
    fn barrier(&mut self, _at: &ElectronState, _boundary: Boundary, _delta_u_ev: f64) {}

    /// With [`BoundaryModel::StepBarrier`]: the electron was reflected at the
    /// face `boundary`; `at` is its state after the reflection. At the
    /// primary's entry this is followed by [`ElectronTally::escaped`].
    fn reflected(&mut self, _at: &ElectronState, _boundary: Boundary) {}

    /// A secondary of the given `generation` (1 for a secondary of the
    /// primary, 2 for one of a secondary, ...) starts at `start`.
    fn begin_secondary(&mut self, _start: &ElectronState, _generation: u32) {}

    /// The current secondary ended with `fate`.
    fn end_secondary(&mut self, _fate: Fate) {}

    /// The electron fell below the stopping threshold (or was trapped) at
    /// `at`. A stop has `at.energy_ev` below the threshold of `at.layer`, a
    /// trapped electron at least that energy
    /// ([`Transport::stopping_thresholds_ev`]).
    fn stopped(&mut self, _at: &ElectronState) {}

    /// The electron left the target through `face`; `at.energy_ev` is the
    /// energy it leaves with (its vacuum energy, with the step barrier).
    fn escaped(&mut self, _at: &ElectronState, _face: Face) {}

    /// The electron was absorbed at the back face ([`EscapeRule::FrontOnly`]).
    fn absorbed(&mut self, _at: &ElectronState) {}

    /// The history ended (its last secondary included); `fate` is the
    /// primary's.
    fn end_history(&mut self, _index: u64, _fate: Fate) {}

    /// Fold `other` (a later chunk) into `self`.
    fn merge(&mut self, other: Self)
    where
        Self: Sized;
}

/// A minimal summary tally: counts, path length and energy bookkeeping.
///
/// The fate counts (`stopped` to `event_capped`) are per history, by the
/// primary's fate; the event counts and energies cover every electron. With
/// every electron ending stopped, trapped or escaped, the energies balance:
/// incident energy plus `secondary_energy_ev` equals `inelastic_loss_ev +
/// escaped_energy_ev + rest_energy_ev + barrier_ev`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SummaryTally {
    /// Histories run.
    pub histories: u64,
    /// Histories whose primary stopped below the threshold.
    pub stopped: u64,
    /// Primaries that left through the front face.
    pub escaped_front: u64,
    /// Primaries that left through the back face.
    pub escaped_back: u64,
    /// Primaries absorbed at the back face.
    pub absorbed: u64,
    /// Trapped primaries.
    pub trapped: u64,
    /// Primaries cut off by the event cap.
    pub event_capped: u64,
    /// Primaries trapped as polarons.
    pub polaron_trapped: u64,
    /// Secondaries created (and transported).
    pub secondaries: u64,
    /// Secondaries that left the target through either face.
    pub secondaries_escaped: u64,
    /// Face reflections (step barrier).
    pub reflections: u64,
    /// Elastic collisions.
    pub elastic_events: u64,
    /// Inelastic collisions.
    pub inelastic_events: u64,
    /// Layer-face crossings between layers.
    pub interface_crossings: u64,
    /// LO-phonon emissions.
    pub phonon_emissions: u64,
    /// LO-phonon absorptions.
    pub phonon_absorptions: u64,
    /// Total flight path, m.
    pub path_m: f64,
    /// Total energy lost in inelastic collisions, eV.
    pub inelastic_loss_ev: f64,
    /// Energy carried out of the target by escaped electrons, eV.
    pub escaped_energy_ev: f64,
    /// Energy of stopped electrons when they stopped, eV.
    pub rest_energy_ev: f64,
    /// Energy given to the lattice by phonon emission, eV.
    pub phonon_emitted_ev: f64,
    /// Energy taken from the lattice by phonon absorption, eV.
    pub phonon_absorbed_ev: f64,
    /// Energy deposited by electrons trapped as polarons, eV.
    pub polaron_deposited_ev: f64,
    /// Kinetic energy of the secondaries created, eV.
    pub secondary_energy_ev: f64,
    /// Sum of [`SecondaryEvent::binding_ev`], eV.
    pub binding_ev: f64,
    /// Sum of [`SecondaryEvent::deposited_ev`], eV.
    pub deposited_ev: f64,
    /// Kinetic energy taken by potential steps, the sum of `-ΔU` over face
    /// transmissions, eV.
    pub barrier_ev: f64,
}

impl ElectronTally for SummaryTally {
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, length_m: f64) {
        self.path_m += length_m;
    }
    fn elastic(&mut self, _after: &ElectronState, _theta: f64) {
        self.elastic_events += 1;
    }
    fn inelastic(&mut self, _after: &ElectronState, w_ev: f64) {
        self.inelastic_events += 1;
        self.inelastic_loss_ev += w_ev;
    }
    fn secondary(&mut self, _p: &ElectronState, e: &SecondaryEvent, c: Option<&ElectronState>) {
        if c.is_some() {
            self.secondaries += 1;
            self.secondary_energy_ev += e.secondary_ev;
        }
        self.binding_ev += e.binding_ev;
        self.deposited_ev += e.deposited_ev;
    }
    fn phonon(
        &mut self,
        _after: &ElectronState,
        kind: PhononEvent,
        hbar_omega_ev: f64,
        _theta: f64,
    ) {
        match kind {
            PhononEvent::Emission => {
                self.phonon_emissions += 1;
                self.phonon_emitted_ev += hbar_omega_ev;
            }
            PhononEvent::Absorption => {
                self.phonon_absorptions += 1;
                self.phonon_absorbed_ev += hbar_omega_ev;
            }
        }
    }
    fn polaron_trapped(&mut self, at: &ElectronState) {
        self.polaron_deposited_ev += at.energy_ev;
    }
    fn interface(&mut self, _at: &ElectronState, _from: usize, _to: usize) {
        self.interface_crossings += 1;
    }
    fn barrier(&mut self, _at: &ElectronState, _b: Boundary, delta_u_ev: f64) {
        self.barrier_ev -= delta_u_ev;
    }
    fn reflected(&mut self, _at: &ElectronState, _b: Boundary) {
        self.reflections += 1;
    }
    fn end_secondary(&mut self, fate: Fate) {
        if matches!(fate, Fate::Escaped(_)) {
            self.secondaries_escaped += 1;
        }
    }
    fn stopped(&mut self, at: &ElectronState) {
        self.rest_energy_ev += at.energy_ev;
    }
    fn escaped(&mut self, at: &ElectronState, _face: Face) {
        self.escaped_energy_ev += at.energy_ev;
    }
    fn end_history(&mut self, _index: u64, fate: Fate) {
        self.histories += 1;
        match fate {
            Fate::Stopped => self.stopped += 1,
            Fate::Escaped(Face::Front) => self.escaped_front += 1,
            Fate::Escaped(Face::Back) => self.escaped_back += 1,
            Fate::Absorbed => self.absorbed += 1,
            Fate::Trapped => self.trapped += 1,
            Fate::EventCap => self.event_capped += 1,
            Fate::PolaronTrapped => self.polaron_trapped += 1,
        }
    }
    fn merge(&mut self, o: Self) {
        self.histories += o.histories;
        self.stopped += o.stopped;
        self.escaped_front += o.escaped_front;
        self.escaped_back += o.escaped_back;
        self.absorbed += o.absorbed;
        self.trapped += o.trapped;
        self.event_capped += o.event_capped;
        self.polaron_trapped += o.polaron_trapped;
        self.secondaries += o.secondaries;
        self.secondaries_escaped += o.secondaries_escaped;
        self.reflections += o.reflections;
        self.elastic_events += o.elastic_events;
        self.inelastic_events += o.inelastic_events;
        self.interface_crossings += o.interface_crossings;
        self.phonon_emissions += o.phonon_emissions;
        self.phonon_absorptions += o.phonon_absorptions;
        self.path_m += o.path_m;
        self.inelastic_loss_ev += o.inelastic_loss_ev;
        self.escaped_energy_ev += o.escaped_energy_ev;
        self.rest_energy_ev += o.rest_energy_ev;
        self.phonon_emitted_ev += o.phonon_emitted_ev;
        self.phonon_absorbed_ev += o.phonon_absorbed_ev;
        self.polaron_deposited_ev += o.polaron_deposited_ev;
        self.secondary_energy_ev += o.secondary_energy_ev;
        self.binding_ev += o.binding_ev;
        self.deposited_ev += o.deposited_ev;
        self.barrier_ev += o.barrier_ev;
    }
}

/// Identity of the tables used in one layer, for [`RunMetadata`].
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LayerMetadata {
    /// Layer index (0 = front).
    pub index: usize,
    /// Depth of the front face, m.
    pub front_m: f64,
    /// Depth of the back face, m; `None` for a semi-infinite substrate.
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
    /// Per-layer table identities.
    pub layers: Vec<LayerMetadata>,
}

/// The result of [`Transport::run`].
#[derive(Debug, Clone, PartialEq)]
pub struct TransportRun<T> {
    /// The merged tally.
    pub tally: T,
    /// The run configuration.
    pub metadata: RunMetadata,
}

/// The electron transport engine: a stack, the tables of each layer, their
/// band parameters (optional) and a configuration.
#[derive(Debug, Clone)]
pub struct Transport {
    stack: Stack,
    tables: Vec<LayerTables>,
    channels: Vec<InsulatorChannels>,
    bands: Option<Vec<BandStructure>>,
    /// Inner potential per layer, eV (zero without band parameters).
    inner: Vec<f64>,
    /// Stopping threshold per layer, eV.
    threshold: Vec<f64>,
    config: TransportConfig,
}

/// An electron waiting on the history's stack, with its generation.
type Pending = (ElectronState, u32);

impl Transport {
    /// Check and assemble, without band parameters. `tables[i]` belongs to
    /// `stack.layers()[i]`. Fails if the configuration needs band parameters
    /// (a [`SecondaryModel`] other than `Off`, or the step barrier); use
    /// [`Transport::with_band_structures`] then.
    pub fn new(
        stack: Stack,
        tables: Vec<LayerTables>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        Self::build(stack, tables, None, config)
    }

    /// Check and assemble with band parameters: `bands[i]` belongs to
    /// `stack.layers()[i]`. They supply the Fermi energy and band gap of the
    /// secondary model and the inner potential of the step barrier and of
    /// [`CutoffReference::VacuumLevel`].
    pub fn with_band_structures(
        stack: Stack,
        tables: Vec<LayerTables>,
        bands: Vec<BandStructure>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        Self::build(stack, tables, Some(bands), config)
    }

    fn build(
        stack: Stack,
        tables: Vec<LayerTables>,
        bands: Option<Vec<BandStructure>>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        if tables.len() != stack.layers().len() {
            return Err(TransportError::LayerCount {
                layers: stack.layers().len(),
                tables: tables.len(),
            });
        }
        for (layer, t) in tables.iter().enumerate() {
            for (slot, table, expected) in [
                ("elastic", &t.elastic, SamplingAxis::ElasticPolarAngle),
                ("inelastic", &t.inelastic, SamplingAxis::InelasticEnergyLoss),
            ] {
                if table.axis() != expected {
                    return Err(TransportError::WrongAxis {
                        layer,
                        slot,
                        found: table.axis(),
                        expected,
                    });
                }
            }
        }
        if !(config.cutoff_ev.is_finite() && config.cutoff_ev > 0.0) {
            return invalid("cutoff", "must be finite and positive");
        }
        if config.max_events == 0 {
            return invalid("max_events", "must be at least 1");
        }
        match &bands {
            Some(b) if b.len() != stack.layers().len() => {
                return invalid(
                    "band structures",
                    format!(
                        "the stack has {} layers but {} band structures were given",
                        stack.layers().len(),
                        b.len()
                    ),
                );
            }
            None if config.needs_bands() => {
                return invalid(
                    "band structures",
                    "secondary generation and the step barrier need one per layer \
                     (Transport::with_band_structures)",
                );
            }
            _ => {}
        }
        let inner: Vec<f64> = match &bands {
            Some(b) => b.iter().map(BandStructure::inner_potential_ev).collect(),
            None => vec![0.0; tables.len()],
        };
        let threshold: Vec<f64> = inner
            .iter()
            .map(|&u| match config.cutoff_reference {
                CutoffReference::BandBottom => config.cutoff_ev,
                CutoffReference::VacuumLevel => u + config.cutoff_ev,
            })
            .collect();
        if config.secondaries != SecondaryModel::Off {
            // Every secondary starts at E_F + W, so a threshold at or below
            // the Fermi level would let them multiply without end; with it
            // above, the energy above the Fermi level, which inelastic events
            // conserve or reduce, bounds the number of electrons followed.
            let b = bands.as_ref().expect("checked above");
            for (layer, (band, &t)) in b.iter().zip(&threshold).enumerate() {
                if t <= band.fermi_ev() {
                    return invalid(
                        "cutoff",
                        format!(
                            "layer {layer}: the stopping threshold {t} eV must exceed the \
                             Fermi energy {} eV under the secondary model",
                            band.fermi_ev()
                        ),
                    );
                }
            }
        }
        let channels = vec![InsulatorChannels::none(); tables.len()];
        Ok(Self {
            stack,
            tables,
            channels,
            bands,
            inner,
            threshold,
            config,
        })
    }

    /// Switch on the insulator channels `channels` in layer `layer` (they
    /// replace that layer's previous choice). Every layer starts with
    /// [`InsulatorChannels::none`]; only opt a layer in when it is a polar
    /// insulator (never a metal). The choice is recorded in the
    /// [`LayerMetadata`] of the run.
    pub fn with_insulator_channels(
        mut self,
        layer: usize,
        channels: InsulatorChannels,
    ) -> Result<Self, TransportError> {
        if layer >= self.channels.len() {
            return invalid(
                "layer",
                format!(
                    "{layer} is out of range (the stack has {} layers)",
                    self.channels.len()
                ),
            );
        }
        self.channels[layer] = channels;
        Ok(self)
    }

    /// The insulator channels of each layer.
    pub fn insulator_channels(&self) -> &[InsulatorChannels] {
        &self.channels
    }

    /// The configuration.
    pub fn config(&self) -> &TransportConfig {
        &self.config
    }

    /// The target.
    pub fn stack(&self) -> &Stack {
        &self.stack
    }

    /// The band parameters per layer, if given.
    pub fn band_structures(&self) -> Option<&[BandStructure]> {
        self.bands.as_deref()
    }

    /// The stopping threshold of each layer, eV: the cutoff under
    /// [`CutoffReference::BandBottom`], `U + cutoff` under
    /// [`CutoffReference::VacuumLevel`] (module docs). An electron reported
    /// by [`ElectronTally::stopped`] with less energy than its layer's
    /// threshold fell below it ([`Fate::Stopped`]); one with at least that
    /// energy is [`Fate::Trapped`].
    pub fn stopping_thresholds_ev(&self) -> &[f64] {
        &self.threshold
    }

    fn step_barrier(&self) -> Option<(bool, bool)> {
        match self.config.boundary {
            BoundaryModel::Transparent => None,
            BoundaryModel::StepBarrier {
                quantum_transmission,
                refraction,
            } => Some((quantum_transmission, refraction)),
        }
    }

    fn check_primary(&self, p: &Primary) -> Result<[f64; 3], TransportError> {
        if !(p.energy_ev.is_finite() && p.energy_ev > self.config.cutoff_ev) {
            return invalid(
                "primary energy",
                format!("{} eV must exceed the cutoff", p.energy_ev),
            );
        }
        let n = (p.direction[0].powi(2) + p.direction[1].powi(2) + p.direction[2].powi(2)).sqrt();
        if !(n.is_finite() && n > 0.0) {
            return invalid("primary direction", "must be finite and nonzero");
        }
        let dir = [p.direction[0] / n, p.direction[1] / n, p.direction[2] / n];
        let inside = if self.step_barrier().is_some() {
            if dir[0] <= 0.0 {
                return invalid(
                    "primary direction",
                    "must point into the target (positive x) with the step barrier",
                );
            }
            p.energy_ev + self.inner[0]
        } else {
            p.energy_ev
        };
        if inside <= self.threshold[0] {
            return invalid(
                "primary energy",
                format!(
                    "{inside} eV inside the first layer must exceed its stopping threshold {} eV",
                    self.threshold[0]
                ),
            );
        }
        Ok(dir)
    }

    /// The metadata of a run with these arguments.
    pub fn metadata(
        &self,
        seed: u64,
        n_histories: u64,
        chunk_size: u64,
        primary: &Primary,
    ) -> RunMetadata {
        RunMetadata {
            cutoff_ev: self.config.cutoff_ev,
            cutoff_reference: self.config.cutoff_reference,
            escape_rule: self.config.escape_rule,
            max_events: self.config.max_events,
            secondaries: self.config.secondaries,
            boundary: self.config.boundary,
            seed,
            n_histories,
            chunk_size,
            primary: *primary,
            layers: self
                .stack
                .layers()
                .iter()
                .zip(&self.tables)
                .zip(&self.channels)
                .enumerate()
                .map(|(index, ((l, t), c))| LayerMetadata {
                    index,
                    front_m: l.front_m(),
                    back_m: l.back_m().is_finite().then_some(l.back_m()),
                    elastic_model: t.elastic.model().to_string(),
                    elastic_provenance: t.elastic.provenance().to_string(),
                    inelastic_model: t.inelastic.model().to_string(),
                    inelastic_provenance: t.inelastic.provenance().to_string(),
                    band_structure: self.bands.as_ref().map(|b| b[index].clone()),
                    phonon: c.phonon.clone(),
                    polaron: c.polaron.clone(),
                })
                .collect(),
        }
    }

    /// Run `n_histories` primaries in parallel (on the current rayon pool),
    /// bit-identically at any thread count. `new_tally` makes an empty tally;
    /// one is made per chunk of `chunk_size` histories.
    pub fn run<T, N>(
        &self,
        seed: u64,
        n_histories: u64,
        chunk_size: u64,
        primary: &Primary,
        new_tally: N,
    ) -> Result<TransportRun<T>, TransportError>
    where
        T: ElectronTally,
        N: Fn() -> T + Sync,
    {
        self.check_primary(primary)?;
        let tally = run_particles(
            seed,
            n_histories,
            chunk_size,
            new_tally,
            |t, rng, i| {
                self.history(t, rng, i, primary)
                    .expect("primary validated before the run");
            },
            |total, part| total.merge(part),
        );
        Ok(TransportRun {
            tally,
            metadata: self.metadata(seed, n_histories, chunk_size, primary),
        })
    }

    /// Simulate one primary, and every secondary it leads to, with the given
    /// random stream. Returns the primary's fate.
    pub fn history<T, R>(
        &self,
        tally: &mut T,
        rng: &mut R,
        index: u64,
        primary: &Primary,
    ) -> Result<Fate, TransportError>
    where
        T: ElectronTally,
        R: Rng,
    {
        let dir = self.check_primary(primary)?;
        let mut st = ElectronState {
            pos: [0.0; 3],
            dir,
            energy_ev: primary.energy_ev,
            layer: 0,
        };
        tally.begin_history(index, &st);
        let mut pending: Vec<Pending> = Vec::new();
        let fate = match self.enter(tally, rng, &mut st) {
            Some(f) => f,
            None => self.follow(tally, rng, &mut st, 0, &mut pending),
        };
        while let Some((mut s, generation)) = pending.pop() {
            tally.begin_secondary(&s, generation);
            let f = self.follow(tally, rng, &mut s, generation, &mut pending);
            tally.end_secondary(f);
        }
        tally.end_history(index, fate);
        Ok(fate)
    }

    /// With the step barrier, take the primary through the front face from
    /// vacuum. Returns its fate if it was reflected.
    fn enter<T: ElectronTally, R: Rng>(
        &self,
        tally: &mut T,
        rng: &mut R,
        st: &mut ElectronState,
    ) -> Option<Fate> {
        let (quantum, refraction) = self.step_barrier()?;
        let du = self.inner[0];
        let boundary = Boundary::Surface(Face::Front);
        // The normal energy is positive and the step is up, so the electron
        // can always get over: one draw.
        let u = uniform(rng);
        match cross_step(st.dir, st.energy_ev, du, quantum, refraction, u) {
            StepOutcome::Transmitted { dir, energy_ev } => {
                st.dir = dir;
                st.energy_ev = energy_ev;
                tally.barrier(st, boundary, du);
                None
            }
            StepOutcome::Reflected { dir } => {
                st.dir = dir;
                tally.reflected(st, boundary);
                tally.escaped(st, Face::Front);
                Some(Fate::Escaped(Face::Front))
            }
        }
    }

    fn follow<T: ElectronTally, R: Rng>(
        &self,
        tally: &mut T,
        rng: &mut R,
        st: &mut ElectronState,
        generation: u32,
        pending: &mut Vec<Pending>,
    ) -> Fate {
        let layers = self.stack.layers();
        let mut events = 0u64;
        loop {
            if events >= self.config.max_events {
                return Fate::EventCap;
            }
            let tabs = &self.tables[st.layer];
            let chans = &self.channels[st.layer];
            let (el, el_at) = rate(&tabs.elastic, st.energy_ev);
            let (inel, inel_at) = rate(&tabs.inelastic, st.energy_ev);
            let (em, ab, tr) = chans.rates(st.energy_ev);
            let total = el + inel + em + ab + tr;
            let layer = &layers[st.layer];
            let mu = st.dir[0];
            let to_face = if mu > 0.0 {
                ((layer.back_m() - st.pos[0]) / mu).max(0.0)
            } else if mu < 0.0 {
                ((layer.front_m() - st.pos[0]) / mu).max(0.0)
            } else {
                f64::INFINITY
            };
            let free = if total > 0.0 {
                -(1.0 - uniform(rng)).ln() / total
            } else {
                f64::INFINITY
            };
            if free.is_infinite() && to_face.is_infinite() {
                tally.stopped(st);
                return Fate::Trapped;
            }
            let from = st.pos;
            if free < to_face {
                for k in 0..3 {
                    st.pos[k] += free * st.dir[k];
                }
                tally.step(from, st, free);
                events += 1;
                match choose(uniform(rng) * total, [el, inel, em, ab, tr]) {
                    0 => {
                        let u = uniform(rng);
                        let theta =
                            sample(&tabs.elastic, el_at, u).clamp(0.0, std::f64::consts::PI);
                        let phi = std::f64::consts::TAU * uniform(rng);
                        st.dir = deflect(st.dir, theta, phi);
                        tally.elastic(st, theta);
                    }
                    1 => {
                        let u = uniform(rng);
                        let w = sample(&tabs.inelastic, inel_at, u).clamp(0.0, st.energy_ev);
                        let threshold = self.threshold[st.layer];
                        match self.config.secondaries {
                            SecondaryModel::Off => {
                                st.energy_ev -= w;
                                tally.inelastic(st, w);
                            }
                            SecondaryModel::KieftBosch {
                                instantaneous_momentum,
                                momentum_conservation,
                            } => {
                                let band =
                                    &self.bands.as_ref().expect("checked in build")[st.layer];
                                // The primary cannot end below the Fermi level:
                                // Nebula clamps its loss table to `K - E_F`
                                // (`kieft_inelastic::create`, commit named in
                                // `electron::boundary`). The threshold exceeds
                                // E_F, so the bound is positive.
                                let w = w.min(st.energy_ev - band.fermi_ev());
                                let out = kieft_bosch(
                                    instantaneous_momentum,
                                    momentum_conservation,
                                    band,
                                    st.dir,
                                    st.energy_ev,
                                    w,
                                    threshold,
                                    rng,
                                );
                                st.energy_ev -= w;
                                st.dir = out.primary_dir;
                                tally.inelastic(st, w);
                                let created = out.secondary.map(|(dir, energy_ev)| ElectronState {
                                    pos: st.pos,
                                    dir,
                                    energy_ev,
                                    layer: st.layer,
                                });
                                tally.secondary(st, &out.event, created.as_ref());
                                if let Some(s) = created {
                                    pending.push((s, generation + 1));
                                }
                            }
                        }
                        if st.energy_ev < threshold {
                            tally.stopped(st);
                            return Fate::Stopped;
                        }
                    }
                    k @ (2 | 3) => {
                        let ph = chans
                            .phonon
                            .as_ref()
                            .expect("a positive phonon rate implies the channel is on");
                        let hw = ph.hbar_omega_ev();
                        let (kind, after) = if k == 2 {
                            (PhononEvent::Emission, st.energy_ev - hw)
                        } else {
                            (PhononEvent::Absorption, st.energy_ev + hw)
                        };
                        let mu = sample_cos_theta(st.energy_ev, after, uniform(rng));
                        let theta = mu.acos();
                        let phi = std::f64::consts::TAU * uniform(rng);
                        st.dir = deflect(st.dir, theta, phi);
                        st.energy_ev = after;
                        tally.phonon(st, kind, hw, theta);
                        if st.energy_ev < self.threshold[st.layer] {
                            tally.stopped(st);
                            return Fate::Stopped;
                        }
                    }
                    _ => {
                        tally.polaron_trapped(st);
                        return Fate::PolaronTrapped;
                    }
                }
                continue;
            }
            // Reach the face exactly; the path is redrawn afterwards.
            let face_x = if mu > 0.0 {
                layer.back_m()
            } else {
                layer.front_m()
            };
            for k in 0..3 {
                st.pos[k] += to_face * st.dir[k];
            }
            st.pos[0] = face_x;
            tally.step(from, st, to_face);
            let last = layers.len() - 1;
            let (to, boundary) = if mu > 0.0 {
                if st.layer < last {
                    let to = st.layer + 1;
                    (Some(to), Boundary::Interface { from: st.layer, to })
                } else if self.config.escape_rule == EscapeRule::BothFaces {
                    (None, Boundary::Surface(Face::Back))
                } else {
                    tally.absorbed(st);
                    return Fate::Absorbed;
                }
            } else if st.layer > 0 {
                let to = st.layer - 1;
                (Some(to), Boundary::Interface { from: st.layer, to })
            } else {
                (None, Boundary::Surface(Face::Front))
            };
            let from_layer = st.layer;
            if let Some((quantum, refraction)) = self.step_barrier() {
                let du = to.map_or(0.0, |t| self.inner[t]) - self.inner[from_layer];
                // A draw only if the electron has the normal energy to get
                // over (cross_step ignores it otherwise).
                let can = st.energy_ev * mu * mu + du > 0.0;
                let u = if can { uniform(rng) } else { 1.0 };
                match cross_step(st.dir, st.energy_ev, du, quantum, refraction, u) {
                    StepOutcome::Transmitted { dir, energy_ev } => {
                        st.dir = dir;
                        st.energy_ev = energy_ev;
                        if let Some(t) = to {
                            st.layer = t;
                        }
                        tally.barrier(st, boundary, du);
                    }
                    StepOutcome::Reflected { dir } => {
                        st.dir = dir;
                        events += 1;
                        tally.reflected(st, boundary);
                        continue;
                    }
                }
            } else if let Some(t) = to {
                st.layer = t;
            }
            match boundary {
                Boundary::Interface { from, to } => {
                    tally.interface(st, from, to);
                    if st.energy_ev < self.threshold[to] {
                        tally.stopped(st);
                        return Fate::Stopped;
                    }
                }
                Boundary::Surface(face) => {
                    tally.escaped(st, face);
                    return Fate::Escaped(face);
                }
            }
        }
    }
}

/// The channel for `x = u · total`: the first index whose cumulative rate
/// exceeds `x`, skipping zero rates; if rounding leaves `x` at or past the
/// sum, the last channel with a positive rate. With only the first two rates
/// positive this is `x < el` → elastic, else inelastic, as before the
/// insulator channels existed. Only called when the total is positive.
fn choose(x: f64, rates: [f64; 5]) -> usize {
    let mut cum = 0.0;
    let mut last = 0;
    for (k, &r) in rates.iter().enumerate() {
        if r > 0.0 {
            cum += r;
            last = k;
            if x < cum {
                return k;
            }
        }
    }
    last
}

/// Uniform on `[0, 1)` with 53 random bits.
pub(crate) fn uniform<R: Rng>(rng: &mut R) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Where an energy sits on a table's grid: lower row and weight of the upper.
#[derive(Clone, Copy)]
struct GridPos {
    i: usize,
    f: f64,
}

/// Inverse mean free path at `e` (linear in energy, constant beyond the ends)
/// and the grid position.
fn rate(t: &CrossSectionTable, e: f64) -> (f64, GridPos) {
    let g = t.energy_ev();
    let r = t.inverse_mfp_per_m();
    let n = g.len();
    let pos = if e <= g[0] {
        GridPos { i: 0, f: 0.0 }
    } else if e >= g[n - 1] {
        GridPos { i: n - 2, f: 1.0 }
    } else {
        let i = g.partition_point(|&x| x <= e) - 1;
        GridPos {
            i,
            f: (e - g[i]) / (g[i + 1] - g[i]),
        }
    };
    let v = if pos.f == 0.0 {
        r[pos.i]
    } else if pos.f == 1.0 {
        r[pos.i + 1]
    } else {
        r[pos.i] + pos.f * (r[pos.i + 1] - r[pos.i])
    };
    (v, pos)
}

/// The inverse CDF at probability `u`, interpolated between the bracketing
/// rows. Only called where the rate is positive, so one of the rows exists.
fn sample(t: &CrossSectionTable, at: GridPos, u: f64) -> f64 {
    let row = |i: usize| t.inverse_cdf(i, u).ok();
    match (at.f, row(at.i), row(at.i + 1)) {
        (f, Some(a), Some(b)) if f > 0.0 && f < 1.0 => a + f * (b - a),
        (f, Some(a), _) if f < 1.0 => a,
        (_, _, Some(b)) => b,
        (_, Some(a), None) => a,
        (_, None, None) => unreachable!("a positive rate implies a stored row"),
    }
}

/// Rotate the unit vector `d` by polar angle `theta` about itself and azimuth
/// `phi` about `d` (measured from an arbitrary but deterministic reference
/// axis), returning a unit vector.
fn deflect(d: [f64; 3], theta: f64, phi: f64) -> [f64; 3] {
    let (st, ct) = theta.sin_cos();
    deflect_cs(d, ct, st, phi)
}

/// [`deflect`] with the polar angle given by its cosine `ct` and sine `st`.
pub(crate) fn deflect_cs(d: [f64; 3], ct: f64, st: f64, phi: f64) -> [f64; 3] {
    let helper = if d[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let e1 = normalize(cross(helper, d));
    let e2 = cross(d, e1);
    let (sp, cp) = phi.sin_cos();
    normalize([
        ct * d[0] + st * (cp * e1[0] + sp * e2[0]),
        ct * d[1] + st * (cp * e1[1] + sp * e2[1]),
        ct * d[2] + st * (cp * e1[2] + sp * e2[2]),
    ])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deflect_keeps_unit_length_and_angle() {
        for d in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            normalize([1.0, 2.0, -3.0]),
        ] {
            for (theta, phi) in [(0.0, 0.0), (0.7, 1.9), (3.0, 5.0)] {
                let n = deflect(d, theta, phi);
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                assert!((len - 1.0).abs() < 1e-12);
                let dot = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
                assert!((dot - theta.cos()).abs() < 1e-12);
            }
        }
    }
}
