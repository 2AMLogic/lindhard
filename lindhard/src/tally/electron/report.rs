//! Plain-data serde report types of the electron tally: everything a run
//! measured, as serializable structs with no behaviour beyond accessors.
//! The accumulators that fill them live in [`super::acc`] (energy) and
//! [`super::charge`] (charge).

use serde::{Deserialize, Serialize};

use super::charge::ChargeReport;
use super::{CartesianGrid, CylindricalGrid, ElectronTallyConfig};
use crate::tally::dose::DoseReport;
use crate::tally::exposure::ExposureReport;
use crate::tally::hist::Histogram;
use crate::tally::moments::MomentSummary;
use crate::tally::psf::RadialProfile;
use crate::tally::table_coverage::LayerTableCoverage;

/// Everything a [`FullElectronTally`](super::FullElectronTally) measured, as
/// plain data.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronReport {
    /// Primary histories.
    pub histories: u64,
    /// The tally settings and conventions.
    pub metadata: ElectronTallyMetadata,
    /// How the primary histories ended.
    pub fates: FateCounts,
    /// Electrons cut off by the collision cap: capped secondary tracks and
    /// histories with any capped electron. A numerical diagnostic; all zero
    /// in reports written before it existed.
    #[serde(default)]
    pub event_caps: EventCapCounts,
    /// Summed energy balance of all histories, eV.
    pub budget: ElectronEnergyBudget,
    /// Emission yields per primary.
    pub yields: Yields,
    /// Electrons leaving the front face (`x = 0`).
    pub front: FaceEmission,
    /// Electrons leaving the back face of a finite stack.
    pub back: FaceEmission,
    /// Electrons leaving through a lateral face of a voxel target; `None`
    /// when none did (always, for a stack).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub side: Option<FaceEmission>,
    /// Where the energy was deposited.
    pub deposition: DepositionReport,
    /// Dose per layer, depth-dose and interface transmission, if configured
    /// ([`ElectronTallyConfig::dose`]).
    #[serde(default)]
    pub dose: Option<DoseReport>,
    /// Energy-weighted moments of the deposition position; `None` if nothing
    /// was deposited.
    pub generation_volume: Option<GenerationVolume>,
    /// Where electrons stopped and where holes were left, and the balance
    /// against the incident electrons. Empty in reports written before it
    /// existed.
    #[serde(default)]
    pub charge: ChargeReport,
    /// Where electrons fell below the stopping threshold.
    pub stopping_points: StoppingPoints,
    /// Per layer, the rate evaluations of the elastic and inelastic tables
    /// inside and beyond each table's energy grid, with the grid bounds
    /// (primaries and secondaries together). These count evaluations, not
    /// collisions; see [`table_coverage`](crate::tally::table_coverage).
    /// Empty in reports written before it existed.
    #[serde(default)]
    pub table_coverage: Vec<LayerTableCoverage>,
}

/// The settings and conventions a report was made with.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronTallyMetadata {
    /// SE/BSE split, eV.
    pub se_bse_split_ev: f64,
    /// How an escaping electron is classed against the split.
    pub se_bse_split_rule: String,
    /// Where the split convention comes from.
    pub se_bse_split_source: String,
    /// Transport cutoff, eV, as configured (measured from the band bottom or
    /// the vacuum level, see the run's
    /// [`RunMetadata`](crate::electron::transport::RunMetadata)).
    pub cutoff_ev: f64,
    /// Stopping threshold per layer, eV: a `stopped` electron below its
    /// layer's threshold is a deposit, at or above it a trapped electron.
    pub stopping_threshold_ev: Vec<f64>,
    /// Layers in the stack.
    pub layers: usize,
    /// The full tally configuration (grids and spectrum binnings).
    pub config: ElectronTallyConfig,
}

/// How the primary histories ended (one per history, from its
/// [`Fate`](crate::electron::transport::Fate)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct FateCounts {
    /// Fell below the stopping threshold in the target.
    pub stopped: u64,
    /// Left through the front face.
    pub escaped_front: u64,
    /// Left through the back face.
    pub escaped_back: u64,
    /// Left through a lateral face of a voxel target.
    #[serde(default)]
    pub escaped_side: u64,
    /// Absorbed at the back face.
    pub absorbed: u64,
    /// No interaction available and no face to reach.
    pub trapped: u64,
    /// Cut off by the collision cap.
    pub event_capped: u64,
    /// Trapped as a polaron (an opt-in insulator channel).
    #[serde(default)]
    pub polaron_trapped: u64,
}

impl FateCounts {
    pub(super) fn merge(&mut self, o: &FateCounts) {
        self.stopped += o.stopped;
        self.escaped_front += o.escaped_front;
        self.escaped_back += o.escaped_back;
        self.escaped_side += o.escaped_side;
        self.absorbed += o.absorbed;
        self.trapped += o.trapped;
        self.event_capped += o.event_capped;
        self.polaron_trapped += o.polaron_trapped;
    }
}

/// Electrons cut off by the collision cap
/// ([`Fate::EventCap`](crate::electron::transport::Fate::EventCap)), a
/// numerical truncation diagnostic rather than a physical fate. The
/// primaries' own caps are [`FateCounts::event_capped`] (one fate per
/// primary); these counts add the secondaries. Both are plain counts
/// (tracks, histories), not energies; the energy the capped electrons still
/// carried is [`ElectronEnergyBudget::event_cap_ev`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct EventCapCounts {
    /// Secondary tracks cut off by the cap, one per capped secondary.
    #[serde(default)]
    pub secondary_tracks: u64,
    /// Primary histories in which the primary or at least one secondary was
    /// cut off by the cap, each counted once.
    #[serde(default)]
    pub affected_histories: u64,
}

impl EventCapCounts {
    pub(super) fn merge(&mut self, o: &EventCapCounts) {
        self.secondary_tracks += o.secondary_tracks;
        self.affected_histories += o.affected_histories;
    }
}

/// The energy balance, eV, summed over all histories and every electron of
/// them (see [the module docs](super#the-energy-balance)):
/// `incident + fermi_sea + phonon_absorbed = deposited + escaped + trapped +
/// barrier`, with
///
/// - `deposited = inelastic + residual + phonon_emitted + polaron`: energy
///   left in the solid at inelastic events, plus the energy electrons had
///   left when they fell below the stopping threshold, plus the `ħω` given to
///   the lattice at LO-phonon emissions, plus the whole remaining energy of
///   electrons trapped as polarons;
/// - `phonon_absorbed`: the `ħω` electrons took from the lattice at LO-phonon
///   absorptions (a source; zero without the insulator channels);
/// - `escaped = escaped_front + escaped_back + escaped_side`;
/// - `trapped = no_interaction + absorbed + event_cap`: energy still carried
///   by electrons that ended in the target without falling below the
///   threshold (no interaction available, absorbed at the back face, or cut
///   off by the collision cap);
/// - `fermi_sea`: the kinetic energy liberated secondaries already had (a
///   source; zero without a secondary model);
/// - `barrier`: the kinetic energy taken by the potential steps at the faces,
///   the sum of `-ΔU` (signed; zero with transparent faces).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElectronEnergyBudget {
    /// Energy of the primaries.
    pub incident_ev: f64,
    /// Deposited in the target, measured from the band bottom. With a
    /// secondary model in a layer with a Fermi energy this includes the
    /// Fermi-sea energy of liberated conduction electrons (`fermi_sea_ev`),
    /// so it can exceed the energy the beam imparted; the net energy imparted
    /// is `incident_ev - escaped_ev = deposited_ev + trapped_ev + barrier_ev -
    /// fermi_sea_ev - phonon_absorbed_ev` (see
    /// [the module docs](super#deposits-and-the-energy-imparted)).
    pub deposited_ev: f64,
    /// Carried out of the target.
    pub escaped_ev: f64,
    /// Carried by electrons that ended in the target at or above the
    /// stopping threshold.
    pub trapped_ev: f64,
    /// Energy left in the solid at inelastic events (part of
    /// `deposited_ev`): the loss `W` without a secondary model; with one,
    /// [`SecondaryEvent::deposited_ev`](crate::electron::secondary::SecondaryEvent::deposited_ev)
    /// plus any positive
    /// [`SecondaryEvent::binding_ev`](crate::electron::secondary::SecondaryEvent::binding_ev).
    pub inelastic_ev: f64,
    /// Remaining energy of electrons that fell below the stopping threshold
    /// (part of `deposited_ev`).
    pub residual_ev: f64,
    /// Out of the front face (part of `escaped_ev`).
    pub escaped_front_ev: f64,
    /// Out of the back face (part of `escaped_ev`).
    pub escaped_back_ev: f64,
    /// Out of a lateral face of a voxel target (part of `escaped_ev`); zero
    /// for a stack.
    #[serde(default)]
    pub escaped_side_ev: f64,
    /// Electrons with no interaction and no face to reach (part of
    /// `trapped_ev`).
    pub no_interaction_ev: f64,
    /// Electrons absorbed at the back face (part of `trapped_ev`).
    pub absorbed_ev: f64,
    /// Electrons cut off by the collision cap (part of `trapped_ev`).
    pub event_cap_ev: f64,
    /// Energy the liberated secondaries brought with them, the sum of
    /// `-binding_ev` over events where it is negative (a source).
    pub fermi_sea_ev: f64,
    /// Kinetic energy taken by the potential steps, the sum of `-ΔU` over
    /// face transmissions (signed).
    pub barrier_ev: f64,
    /// `ħω` given to the lattice by LO-phonon emission (part of
    /// `deposited_ev`).
    #[serde(default)]
    pub phonon_emitted_ev: f64,
    /// `ħω` taken from the lattice by LO-phonon absorption (a source).
    #[serde(default)]
    pub phonon_absorbed_ev: f64,
    /// Remaining energy of electrons trapped as polarons (part of
    /// `deposited_ev`).
    #[serde(default)]
    pub polaron_ev: f64,
    /// Ionisation channel of `deposited_ev`: `inelastic_ev + residual_ev`
    /// (the residual is placed here by convention, see
    /// [the channel split](super#the-ionisation-and-lattice-channels)).
    #[serde(default)]
    pub ionisation_ev: f64,
    /// Lattice channel of `deposited_ev`: `phonon_emitted_ev + polaron_ev`.
    #[serde(default)]
    pub lattice_ev: f64,
    /// `|deposited + escaped + trapped + barrier - incident - fermi_sea| /
    /// (incident + fermi_sea)` (0 with no incident energy).
    pub relative_imbalance: f64,
}

/// Emission yields per primary.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Yields {
    /// Backscatter yield `η`: fast electrons out of the front face per
    /// primary.
    pub backscatter_eta: f64,
    /// Secondary yield `δ`: slow electrons out of the front face per primary.
    pub secondary_delta: f64,
    /// Total yield `σ = η + δ`.
    pub total_sigma: f64,
    /// Fast electrons out of the back face per primary.
    pub transmitted_fast: f64,
    /// Slow electrons out of the back face per primary.
    pub transmitted_slow: f64,
}

/// Electrons leaving through one face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FaceEmission {
    /// Electrons.
    pub count: u64,
    /// Electrons per primary.
    pub per_primary: f64,
    /// Sum of their energies, eV.
    pub energy_ev: f64,
    /// Energy spectrum, eV, both classes.
    pub energy_histogram: Histogram,
    /// Slow electrons (`E < split`, the secondary class).
    pub slow: EmissionClass,
    /// Fast electrons (`E >= split`, the backscattered class).
    pub fast: EmissionClass,
}

/// One energy class of the electrons leaving through a face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EmissionClass {
    /// Electrons.
    pub count: u64,
    /// Electrons per primary.
    pub per_primary: f64,
    /// Sum of their energies, eV.
    pub energy_ev: f64,
    /// Mean energy, eV (0 if none).
    pub mean_energy_ev: f64,
    /// Polar angle from the outward surface normal, rad.
    pub polar_histogram: Histogram,
}

/// Deposited energy by layer and on the configured grids, eV.
///
/// Deposits are measured from the band bottom. With a secondary model in a
/// layer with a Fermi energy they include the Fermi-sea energy of liberated
/// conduction electrons, so they sum to more than the energy the beam
/// imparted; see [`ElectronEnergyBudget::deposited_ev`] and
/// [the module docs](super#deposits-and-the-energy-imparted).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DepositionReport {
    /// By the layer the deposit was made in.
    pub per_layer_ev: Vec<f64>,
    /// The ionisation part of `per_layer_ev`, eV (electronic excitation,
    /// including the residual of sub-cutoff electrons; see
    /// [the channel split](super#the-ionisation-and-lattice-channels)). Empty
    /// in reports written before it existed.
    #[serde(default)]
    pub per_layer_ionisation_ev: Vec<f64>,
    /// The lattice part of `per_layer_ev`, eV (LO-phonon emission and polaron
    /// traps); `ionisation + lattice = per_layer_ev` per layer.
    #[serde(default)]
    pub per_layer_lattice_ev: Vec<f64>,
    /// On the Cartesian grid, if configured.
    pub cartesian: Option<CartesianDeposition>,
    /// On the cylindrical grid, if configured.
    pub cylindrical: Option<CylindricalDeposition>,
    /// Radial profile in the PSF depth slab, with per-history errors, if
    /// configured ([`ElectronTallyConfig::psf`]); fit it with
    /// [`fit_psf`](crate::tally::psf::fit_psf).
    #[serde(default)]
    pub psf: Option<RadialProfile>,
    /// Exposure statistics, if configured ([`ElectronTallyConfig::exposure`]).
    #[serde(default)]
    pub exposure: Option<ExposureReport>,
}

/// Deposited energy on a [`CartesianGrid`].
///
/// Measured from the band bottom: with a secondary model in a layer with a
/// Fermi energy the cells include the Fermi-sea energy of liberated
/// conduction electrons, so they sum to more than the energy the beam
/// imparted (see [`ElectronEnergyBudget::deposited_ev`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CartesianDeposition {
    /// The grid.
    pub grid: CartesianGrid,
    /// Energy per voxel, eV, at [`CartesianGrid::index`].
    pub energy_ev: Vec<f64>,
    /// Energy deposited outside the grid, eV.
    pub outside_ev: f64,
    /// Ionisation part of `energy_ev` per voxel, eV; `ionisation + lattice =
    /// energy` per voxel. Empty in reports written before it existed.
    #[serde(default)]
    pub ionisation_ev: Vec<f64>,
    /// Lattice part of `energy_ev` per voxel, eV.
    #[serde(default)]
    pub lattice_ev: Vec<f64>,
    /// Ionisation part of `outside_ev`, eV.
    #[serde(default)]
    pub outside_ionisation_ev: f64,
    /// Lattice part of `outside_ev`, eV.
    #[serde(default)]
    pub outside_lattice_ev: f64,
}

impl CartesianDeposition {
    /// Energy in voxel `(ix, iy, iz)`, eV.
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.energy_ev[self.grid.index(ix, iy, iz)]
    }

    /// Energy per volume in each voxel, eV/m³.
    pub fn density_ev_per_m3(&self) -> Vec<f64> {
        let v = self.grid.voxel_volume_m3();
        self.energy_ev.iter().map(|e| e / v).collect()
    }
}

/// Deposited energy on a [`CylindricalGrid`].
///
/// Measured from the band bottom: with a secondary model in a layer with a
/// Fermi energy the cells include the Fermi-sea energy of liberated
/// conduction electrons, so they sum to more than the energy the beam
/// imparted (see [`ElectronEnergyBudget::deposited_ev`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CylindricalDeposition {
    /// The grid.
    pub grid: CylindricalGrid,
    /// Energy per cell, eV, at [`CylindricalGrid::index`].
    pub energy_ev: Vec<f64>,
    /// Energy deposited outside the grid, eV.
    pub outside_ev: f64,
    /// Ionisation part of `energy_ev` per cell, eV; `ionisation + lattice =
    /// energy` per cell. Empty in reports written before it existed.
    #[serde(default)]
    pub ionisation_ev: Vec<f64>,
    /// Lattice part of `energy_ev` per cell, eV.
    #[serde(default)]
    pub lattice_ev: Vec<f64>,
    /// Ionisation part of `outside_ev`, eV.
    #[serde(default)]
    pub outside_ionisation_ev: f64,
    /// Lattice part of `outside_ev`, eV.
    #[serde(default)]
    pub outside_lattice_ev: f64,
}

impl CylindricalDeposition {
    /// Energy in cell `(ir, ix)`, eV.
    pub fn at(&self, ir: usize, ix: usize) -> f64 {
        self.energy_ev[self.grid.index(ir, ix)]
    }

    /// Energy per volume in each cell, eV/m³.
    pub fn density_ev_per_m3(&self) -> Vec<f64> {
        self.energy_ev
            .iter()
            .enumerate()
            .map(|(i, e)| e / self.grid.cell_volume_m3(i / self.grid.depth.bins))
            .collect()
    }
}

/// Energy-weighted moments of the deposition position (weights: the energy
/// of each deposit). Lengths in metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GenerationVolume {
    /// Total weight: the deposited energy, eV.
    pub energy_ev: f64,
    /// Weighted mean of `[x, y, z]`.
    pub mean_m: [f64; 3],
    /// Weighted standard deviation of `[x, y, z]`.
    pub std_dev_m: [f64; 3],
    /// Weighted mean distance from the beam axis.
    pub mean_radius_m: f64,
    /// Weighted standard deviation of the distance from the beam axis.
    pub radius_std_dev_m: f64,
    /// Weighted RMS distance from the beam axis, `<y² + z²>^(1/2)`.
    pub rms_radius_m: f64,
}

/// Where electrons fell below the stopping threshold (unweighted, one sample
/// each).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoppingPoints {
    /// Electrons that fell below the stopping threshold.
    pub stopped: u64,
    /// Depth moments, m; `None` with fewer than two.
    pub depth: Option<MomentSummary>,
    /// Moments of the distance from the beam axis, m.
    pub radial: Option<MomentSummary>,
    /// The primaries among `stopped`: those that fell below the threshold
    /// before any secondary of their history began (the transport follows a
    /// history's primary to its end first), i.e. histories whose fate is
    /// [`Fate::Stopped`](crate::electron::transport::Fate::Stopped). Absent
    /// from reports written before it existed.
    #[serde(default)]
    pub primaries: PrimaryStoppingPoints,
}

/// Where the primaries fell below the stopping threshold (unweighted, one
/// sample per stopped primary): the penetration depth of stopped primaries,
/// as compared with other codes in `docs/validation.md`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PrimaryStoppingPoints {
    /// Primaries that fell below the stopping threshold.
    pub stopped: u64,
    /// Depth moments, m; `None` with fewer than two.
    pub depth: Option<MomentSummary>,
    /// Moments of the distance from the beam axis, m.
    pub radial: Option<MomentSummary>,
}
