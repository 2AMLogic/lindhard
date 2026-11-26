//! What the electron transport reports, and a small built-in summary tally.
//!
//! The engine ([`Transport`](super::Transport)) calls the [`ElectronTally`]
//! hooks as events happen. Every hook has a no-op default except `merge`, so
//! a tally implements only what it needs; the full tally of a run is
//! [`crate::tally::electron`]. Tallies are merged in chunk order by
//! [`crate::rng::run_particles`], so a tally whose `merge` is order-faithful
//! is bit-identical at any thread count.

use super::types::{Boundary, ElectronState, Fate, GridCoverage, PhononEvent, TableChannel};
use crate::electron::inelastic::InnerShell;
use crate::electron::secondary::SecondaryEvent;
use crate::geometry::Face;

/// Event hooks called by [`Transport`](super::Transport). All default to no-ops except `merge`.
/// The `ElectronState` passed is the state **after** the event.
pub trait ElectronTally: Send {
    /// A primary history starts. With
    /// [`BoundaryModel::StepBarrier`](super::types::BoundaryModel::StepBarrier),
    /// `start` is the primary in vacuum, before it meets the front face.
    fn begin_history(&mut self, _index: u64, _start: &ElectronState) {}

    /// A flight segment of `length_m` ended at `end.pos` (at a collision, a
    /// layer face or the target face). It started at `from` with energy
    /// `end.energy_ev` (the energy does not change in flight).
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, _length_m: f64) {}

    /// The rate of `channel`'s table in layer `at.layer` was evaluated at
    /// `at.energy_ev`, which lies `coverage` against that table's grid.
    ///
    /// Called twice (elastic, then inelastic) at the start of every pass of
    /// the transport loop that evaluates the rates, before the free path is
    /// drawn: for every flight, whether it ends at a collision, at a layer
    /// face (the path is then discarded and the rates are evaluated again in
    /// the next pass) or nowhere (zero total rate, [`Fate::Trapped`]); after
    /// a face reflection; and for a channel whose rate is zero. It is not
    /// called for the primary's entry through the front face (no table is
    /// read), nor once an electron has reached the event cap, nor for the
    /// tables of inner-shell channels. The inverse CDF
    /// sampled at a collision uses the grid position of the same evaluation
    /// and is not reported again. So these are **rate evaluations**, not
    /// collisions: their number exceeds the number of collisions by the
    /// number of boundary-limited and rate-less flights.
    fn table_lookup(
        &mut self,
        _at: &ElectronState,
        _channel: TableChannel,
        _coverage: GridCoverage,
    ) {
    }

    /// An elastic collision deflected the electron by `theta` rad.
    fn elastic(&mut self, _after: &ElectronState, _theta: f64) {}

    /// An inelastic collision took `w_ev` from the electron. `after` carries
    /// the new direction if the secondary model deflected it.
    fn inelastic(&mut self, _after: &ElectronState, _w_ev: f64) {}

    /// Called right after [`ElectronTally::inelastic`] for every inelastic
    /// event drawn from an inner-shell channel
    /// ([`Transport::with_inner_shells`](super::Transport::with_inner_shells)),
    /// with or without a secondary model: the shell (its element, subshell
    /// and binding energy) and the loss `w_ev` (the same as passed to
    /// [`ElectronTally::inelastic`]). Not called for valence events. Comes
    /// before [`ElectronTally::secondary`].
    fn inner_shell(&mut self, _after: &ElectronState, _shell: &InnerShell, _w_ev: f64) {}

    /// Called right after [`ElectronTally::inelastic`] (and
    /// [`ElectronTally::inner_shell`], for a shell event) for every inelastic
    /// event when a [`SecondaryModel`](crate::electron::secondary::SecondaryModel)
    /// is on: the event's energy bookkeeping, and the secondary that was
    /// created (`None` if none was). The secondary is transported later,
    /// between [`ElectronTally::begin_secondary`] and
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
    /// [`BoundaryModel::StepBarrier`](super::types::BoundaryModel::StepBarrier)).
    fn interface(&mut self, _at: &ElectronState, _from_layer: usize, _to_layer: usize) {}

    /// With [`BoundaryModel::StepBarrier`](super::types::BoundaryModel::StepBarrier):
    /// the electron got through the face `boundary`, and its kinetic energy
    /// changed by `delta_u_ev` (the step `U' - U`). `at` is its state just
    /// past the face (in vacuum, for a surface). Comes before
    /// [`ElectronTally::interface`] or [`ElectronTally::escaped`].
    fn barrier(&mut self, _at: &ElectronState, _boundary: Boundary, _delta_u_ev: f64) {}

    /// With [`BoundaryModel::StepBarrier`](super::types::BoundaryModel::StepBarrier):
    /// the electron was reflected at the face `boundary`; `at` is its state
    /// after the reflection. At the primary's entry this is followed by
    /// [`ElectronTally::escaped`].
    fn reflected(&mut self, _at: &ElectronState, _boundary: Boundary) {}

    /// A secondary of the given `generation` (1 for a secondary of the
    /// primary, 2 for one of a secondary, ...) starts at `start`.
    fn begin_secondary(&mut self, _start: &ElectronState, _generation: u32) {}

    /// The current secondary ended with `fate`.
    fn end_secondary(&mut self, _fate: Fate) {}

    /// The electron fell below the stopping threshold (or was trapped) at
    /// `at`. A stop has `at.energy_ev` below the threshold of `at.layer`, a
    /// trapped electron at least that energy
    /// ([`Transport::stopping_thresholds_ev`](super::Transport::stopping_thresholds_ev)).
    fn stopped(&mut self, _at: &ElectronState) {}

    /// The electron left the target through `face`; `at.energy_ev` is the
    /// energy it leaves with (its vacuum energy, with the step barrier).
    fn escaped(&mut self, _at: &ElectronState, _face: Face) {}

    /// The electron was absorbed at the back face
    /// ([`EscapeRule::FrontOnly`](super::types::EscapeRule::FrontOnly)).
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
    /// Primaries that left through a lateral face of a voxel grid.
    pub escaped_side: u64,
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
    /// Secondaries cut off by the event cap ([`Fate::EventCap`]), one per
    /// capped secondary track. A numerical-truncation diagnostic: zero
    /// `event_capped` primaries does not mean no secondary was capped.
    pub secondaries_event_capped: u64,
    /// Face reflections (step barrier).
    pub reflections: u64,
    /// Elastic collisions.
    pub elastic_events: u64,
    /// Inelastic collisions (valence and inner-shell).
    pub inelastic_events: u64,
    /// Inelastic collisions in an inner-shell channel (part of
    /// `inelastic_events`).
    pub inner_shell_events: u64,
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
    /// Energy lost in inner-shell collisions, eV (part of
    /// `inelastic_loss_ev`).
    pub inner_shell_loss_ev: f64,
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
    fn inner_shell(&mut self, _after: &ElectronState, _shell: &InnerShell, w_ev: f64) {
        self.inner_shell_events += 1;
        self.inner_shell_loss_ev += w_ev;
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
        match fate {
            Fate::Escaped(_) => self.secondaries_escaped += 1,
            Fate::EventCap => self.secondaries_event_capped += 1,
            _ => {}
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
            Fate::Escaped(Face::Side) => self.escaped_side += 1,
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
        self.escaped_side += o.escaped_side;
        self.absorbed += o.absorbed;
        self.trapped += o.trapped;
        self.event_capped += o.event_capped;
        self.polaron_trapped += o.polaron_trapped;
        self.secondaries += o.secondaries;
        self.secondaries_escaped += o.secondaries_escaped;
        self.secondaries_event_capped += o.secondaries_event_capped;
        self.reflections += o.reflections;
        self.elastic_events += o.elastic_events;
        self.inelastic_events += o.inelastic_events;
        self.inner_shell_events += o.inner_shell_events;
        self.interface_crossings += o.interface_crossings;
        self.phonon_emissions += o.phonon_emissions;
        self.phonon_absorptions += o.phonon_absorptions;
        self.path_m += o.path_m;
        self.inelastic_loss_ev += o.inelastic_loss_ev;
        self.inner_shell_loss_ev += o.inner_shell_loss_ev;
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
