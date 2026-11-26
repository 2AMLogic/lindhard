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
//! # Targets other than a stack
//!
//! [`Transport<G>`] is generic over [`Geometry`]; [`Transport::new`] and
//! [`Transport::with_band_structures`] build the layered `Transport<Stack>`,
//! [`Transport::on_geometry`] and
//! [`Transport::on_geometry_with_band_structures`] any other geometry (for
//! example a [`crate::geometry::VoxelGrid`]). Tables, band parameters,
//! channels and thresholds are indexed by **material**
//! ([`Geometry::material_index`]), and [`ElectronState::layer`] is the
//! electron's material index (the layer index, for a stack). The flight uses
//! [`Geometry::flight`]: a face between two regions of one material is not an
//! event and the free path is not redrawn there. A primary starts at
//! [`Geometry::entry_point`]. An electron leaving through a lateral face of a
//! voxel grid ends with [`Fate::Escaped`]`(`[`Face::Side`]`)`, whatever the
//! [`EscapeRule`]. The [`Boundary::Interface`] indices are region indices.
//! The step barrier is applied to the `x` component of the direction, the
//! face normal of a stack only, so the `on_geometry` constructors accept
//! [`BoundaryModel::Transparent`] only. For a stack the results are those of
//! before the generalisation, bit for bit.
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
//! - **Inner-shell channels (opt-in per layer).** A layer given
//!   [`ShellChannelTable`]s through [`Transport::with_inner_shells`] adds one
//!   rate per shell `j`, read from that shell's table like the others; the
//!   layer's inelastic table is then the valence channel. The construction
//!   is **channel-first**: the inelastic rate is the sum of the valence and
//!   shell rates, a channel is chosen with probability proportional to its
//!   rate (below), and the loss `ω` is drawn from the chosen channel's own
//!   inverse CDF. (The conditional choice of
//!   [`crate::electron::inelastic::ShellResolvedChannels::sample_channel`]
//!   belongs with a loss drawn from the total DIIMFP and is not used.) A
//!   shell channel is **closed** (rate zero) when the electron cannot pay its
//!   binding energy `B_j`: when `E <= B_j` without a secondary model and
//!   `E - E_F <= B_j` with one. An open channel's `ω` is kept in
//!   `[B_j, E]` (`[B_j, E - E_F]` with a secondary model), so the primary
//!   never ends below the Fermi level and the shell's secondary energy is
//!   never negative. With [`SecondaryModel::Off`] the event only takes `ω`;
//!   with [`SecondaryModel::KieftBosch`] it liberates an electron of energy
//!   `E_F + ω - B_j` and leaves `B_j - E_F` in the solid
//!   ([`crate::electron::secondary`], "Inner-shell events"). Layers default
//!   to no shell channel, which is exactly the valence-only transport.
//! - **Channel choice.** With the rates in the fixed order elastic,
//!   inelastic (valence), phonon emission, phonon absorption, polaron
//!   trapping, then the inner shells in the order given, the channel is the
//!   first whose cumulative rate exceeds `u · total`; channels with zero
//!   rate are never chosen (a rounding overshoot falls to the last channel
//!   with a positive rate). A layer without shells sums and compares exactly
//!   the five rates of before.
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
//! step barrier, the vacuum energy plus the inner potential `U`). The models
//! of [`crate::electron::inelastic`] have their own axis instead: `T`, the
//! energy above the *model's* Fermi energy, with kinematics on
//! `T' = T + E_F(model)` and losses up to `T` (Shinotsuka et al. 2017, eqs.
//! (2)-(3), as in [`crate::electron::inelastic::penn`]). The table has to
//! bridge the two, and it does so when it is built
//! ([`crate::electron::inelastic::table`], "Energy axis"):
//!
//! - **Band-bottom axis** ([`EnergyAxis::BandBottom`], with the model's Fermi
//!   energy set to the band's): the row at `E` is the model's at
//!   `T = E - E_F`, so its kinematics are on `T' = E` and its losses stop at
//!   `E - E_F`. This is the convention of the cstool table compiler which
//!   builds Nebula's tables (rows at the kinetic energy `K`, kinematics on
//!   `K`, losses below `K - F`, the "Fermi correction" of
//!   `compile_full_imfp_icdf` in `cstool/dielectric_function/compile.py`,
//!   Nebula-simulator/cstool commit `0c739eb`), and `F` is what cstool
//!   passes there, [`BandStructure::min_excitation_ev`]: the Fermi energy of
//!   a metal, the conduction-band bottom of an insulator. It is how
//!   `lindhard run` builds the table of every material with a band. The
//!   losses of such a table are below the clamp of
//!   [`SecondaryModel::KieftBosch`] above, which therefore does not act
//!   (test `band_bottom_table_never_reaches_the_fermi_clamp`,
//!   `tests/electron_secondaries.rs`). One exception remains, for either
//!   axis: between the last empty row of a table and its first row with
//!   losses the rate is interpolated from zero and the loss is drawn from
//!   that first row, so it can exceed what an electron below that row may
//!   lose; there the clamps (`W <= E`, and `W <= E - E_F` under Kieft-Bosch)
//!   are still the bound.
//! - **Model axis** ([`EnergyAxis::ModelFermiLevel`], the builder's default,
//!   and what `lindhard run` builds for a material without a band, with
//!   `[electron.inelastic] fermi_energy_ev`): the row read at `E` is the
//!   model's at `T = E`. With the model's Fermi energy 0 the kinematics are
//!   right (`T' = E`) but the losses reach `E`; in a layer with a band and
//!   the Kieft-Bosch model the clamp takes every loss beyond `E - E_F` to
//!   `E - E_F`, leaving the primary at the Fermi level and giving the
//!   secondary its whole energy. Setting the model's Fermi energy to the
//!   band's on this axis counts `E_F` twice (`T' = E + E_F`). The tests
//!   `inelastic_table_is_read_at_the_band_bottom_energy` and
//!   `losses_beyond_the_fermi_level_are_clamped_to_it` pin how the transport
//!   reads and clamps such a table.
//!
//! Inner-shell channel tables
//! ([`crate::electron::inelastic::table::build_shell_channel_tables`]) are
//! built on the axis of the options they are given, like the valence table,
//! and are read the same way, at the band-bottom energy `E`. `lindhard run`
//! does not build them yet (no input, #156); a library caller picks the axis.
//! On top of the table, the transport applies the shell's own limits on the
//! band-bottom axis: the channel is closed unless the largest allowed loss
//! (`E`, or `E - E_F` under the secondary model) exceeds `B_j`, and `ω` is
//! clamped into `[B_j`, that loss`]`. Near an edge, a table on the model
//! axis (losses up to `E` above the model's Fermi energy) has losses above
//! `E - E_F` that this clamp moves down to `E - E_F`, as for valence losses;
//! a band-bottom table's losses are already below `E - E_F`, apart from the
//! interpolation exception above.
//!
//! **Measured effect (#173).** Until #241 `lindhard run` built every table
//! on the model axis with Fermi energy 0. With the single-pole Penn default,
//! the Al input of `validation/experiments/se_yield.py` (`E_F` = 11.66 eV)
//! then sent 78 to 82 % of the inelastic events of electrons 5 to 20 eV above
//! the Fermi level into the clamp (Al 400 eV, 200 histories). A diagnostic
//! table in the cstool convention (the model's rows at `T = E - E_F`,
//! relabelled `E`, so on a grid shifted by `E_F`) removed every clamped
//! event and changed δ by a few per cent, against a factor of about 5 to
//! explain (1000 histories, seed 1): single pole Al 7.61 to 8.17 and Au 3.23
//! to 3.37 (800 eV); Mermin Al 1.34 to 1.26 (400 eV). With 200 histories:
//! single pole Cu 1.67 to 1.80, Mermin Au 2.62 to 2.58 and Cu 2.02 (800 eV).
//! The gap between the single pole and Mermin was the same in both. The axis
//! was therefore a real mismatch but not the cause of the δ overestimate of
//! #173; that is the single-pole model's inelastic mean free path at low
//! energy ([`crate::electron::inelastic::penn`], "Low energies").
//!
//! With the band-bottom table `lindhard run` builds since #241 (on the
//! input's own grid), the same Al input at 400 eV (200 histories, seed 1, so
//! about ±0.1 to 0.2 in δ) has no clamped event in any energy bin, and δ is
//! 6.44 against 6.89 with the table built as before. The commands are in
//! `lindhard-cli/examples/inelastic_low_energy.rs` (modes `clamp` and
//! `legacy`).
//!
//! [`CrossSectionTable`]: crate::electron::data::CrossSectionTable
//! [`EnergyAxis::BandBottom`]: crate::electron::inelastic::table::EnergyAxis::BandBottom
//! [`EnergyAxis::ModelFermiLevel`]: crate::electron::inelastic::table::EnergyAxis::ModelFermiLevel
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
//! two more uniforms. An inner-shell event draws one uniform for `ω`, then,
//! under the Kieft-Bosch model, the secondary's azimuth and, with the
//! instantaneous momentum on, two more uniforms (a shell event always
//! liberates an electron). With no insulator or inner-shell channel switched
//! on, the draws are exactly those of a run without them. A path drawn for a flight that ends at a layer face is discarded;
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
//! [`ElectronTally::phonon`] and [`ElectronTally::polaron_trapped`]. An
//! inner-shell event calls [`ElectronTally::inelastic`] like a valence one,
//! then [`ElectronTally::inner_shell`] with the shell, then (with a
//! secondary model) [`ElectronTally::secondary`].
//!
//! [`ElectronTally::table_lookup`] is a diagnostic hook, also with a no-op
//! default: it reports, per channel, whether each rate evaluation read the
//! table inside its energy grid or used the constant continuation beyond an
//! end ([`GridCoverage`]). It draws nothing and changes nothing, so a tally
//! that uses it sees the same histories as one that does not.
//! [`crate::tally::TableCoverageTally`] counts these evaluations. The rates
//! of inner-shell channel tables are not reported to it.

use rand_core::Rng;

use crate::electron::boundary::{cross_step, BandStructure, StepOutcome};
use crate::electron::data::{SamplingAxis, ShellChannelTable};
use crate::electron::inelastic::InnerShell;
use crate::electron::phonon::{sample_cos_theta, InsulatorChannels};
use crate::electron::secondary::{kieft_bosch, kieft_bosch_shell, Outcome, SecondaryModel};
use crate::geometry::{ExitOutcome, Flight, Geometry, Stack};
use crate::rng::{run_particles, uniform};

mod sampling;
mod tally;
mod types;

/// The beam of a run: the lateral profile on the front face and the
/// angular spread about [`Primary::direction`]. Both lived here until
/// they were split into [`crate::electron::beam`]; the re-export keeps
/// `electron::transport::BeamProfile` and `BeamSpread` resolving.
pub use crate::electron::beam::{BeamProfile, BeamSpread};
pub use crate::geometry::Face;
pub use tally::{ElectronTally, SummaryTally};
pub use types::{
    Boundary, BoundaryModel, CutoffReference, ElectronState, EscapeRule, Fate, GridCoverage,
    InnerShellMetadata, LayerMetadata, LayerTables, PhononEvent, Primary, RunMetadata,
    TableChannel, TransportConfig, TransportError, TransportRun,
};

pub(crate) use sampling::deflect_cs;
pub(crate) use types::invalid;

use sampling::{choose, deflect, rate, sample, GridPos};

/// The electron transport engine: a stack, the tables of each layer, their
/// band parameters (optional) and a configuration.
#[derive(Debug, Clone)]
pub struct Transport<G = Stack> {
    geometry: G,
    /// `(front, back)` depth of each material's layer, for the metadata of a
    /// stack run (`None` for other geometries).
    extents: Vec<(Option<f64>, Option<f64>)>,
    tables: Vec<LayerTables>,
    channels: Vec<InsulatorChannels>,
    /// Inner-shell channel tables per layer (empty by default).
    shells: Vec<Vec<ShellChannelTable>>,
    bands: Option<Vec<BandStructure>>,
    /// Inner potential per layer, eV (zero without band parameters).
    inner: Vec<f64>,
    /// Stopping threshold per layer, eV.
    threshold: Vec<f64>,
    config: TransportConfig,
}

/// An electron waiting on the history's stack, with its region and generation.
type Pending = (ElectronState, usize, u32);

impl Transport<Stack> {
    /// Check and assemble, without band parameters. `tables[i]` belongs to
    /// `stack.layers()[i]`. Fails if the configuration needs band parameters
    /// (a [`SecondaryModel`] other than `Off`, or the step barrier); use
    /// [`Transport::with_band_structures`] then.
    pub fn new(
        stack: Stack,
        tables: Vec<LayerTables>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        let extents = stack_extents(&stack);
        Self::build(stack, extents, tables, None, config)
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
        let extents = stack_extents(&stack);
        Self::build(stack, extents, tables, Some(bands), config)
    }

    /// The target.
    pub fn stack(&self) -> &Stack {
        &self.geometry
    }
}

/// `(front, back)` depth of each layer of `stack`, for [`LayerMetadata`].
fn stack_extents(stack: &Stack) -> Vec<(Option<f64>, Option<f64>)> {
    stack
        .layers()
        .iter()
        .map(|l| {
            (
                Some(l.front_m()),
                l.back_m().is_finite().then_some(l.back_m()),
            )
        })
        .collect()
}

impl<G: Geometry> Transport<G> {
    /// [`Transport::new`] on any [`Geometry`] (a [`crate::geometry::VoxelGrid`],
    /// a mesh, a CSG solid): `tables[m]` belongs to material `m` of the
    /// geometry ([`Geometry::material`]), whatever the regions that hold it.
    /// Everything per material (tables, band parameters, insulator and
    /// inner-shell channels, stopping thresholds) is indexed by material; the
    /// electron's [`ElectronState::layer`] is then its material index.
    ///
    /// The primary starts at [`Geometry::entry_point`].
    pub fn on_geometry(
        geometry: G,
        tables: Vec<LayerTables>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        Self::on_geometry_checked(geometry, tables, None, config)
    }

    /// [`Transport::with_band_structures`] on any [`Geometry`]; `bands[m]`
    /// belongs to material `m`.
    pub fn on_geometry_with_band_structures(
        geometry: G,
        tables: Vec<LayerTables>,
        bands: Vec<BandStructure>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        Self::on_geometry_checked(geometry, tables, Some(bands), config)
    }

    fn on_geometry_checked(
        geometry: G,
        tables: Vec<LayerTables>,
        bands: Option<Vec<BandStructure>>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        if config.boundary != BoundaryModel::Transparent {
            // The step is applied to the normal component along x
            // (`cross_step`), which is the face normal of a stack only.
            return invalid(
                "boundary",
                "the step barrier is implemented for the layered Stack only; \
                 use BoundaryModel::Transparent on other geometries",
            );
        }
        let extents = vec![(None, None); geometry.n_materials()];
        Self::build(geometry, extents, tables, bands, config)
    }

    /// The target geometry.
    pub fn geometry(&self) -> &G {
        &self.geometry
    }

    /// Number of materials, which is the number of table pairs.
    pub fn n_materials(&self) -> usize {
        self.tables.len()
    }

    fn build(
        geometry: G,
        extents: Vec<(Option<f64>, Option<f64>)>,
        tables: Vec<LayerTables>,
        bands: Option<Vec<BandStructure>>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        if tables.len() != geometry.n_materials() {
            return Err(TransportError::LayerCount {
                layers: geometry.n_materials(),
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
            Some(b) if b.len() != geometry.n_materials() => {
                return invalid(
                    "band structures",
                    format!(
                        "the stack has {} layers but {} band structures were given",
                        geometry.n_materials(),
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
        let shells = vec![Vec::new(); tables.len()];
        Ok(Self {
            geometry,
            extents,
            tables,
            channels,
            shells,
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

    /// Add the inner-shell ionisation channels `shells` to layer `layer`
    /// (they replace that layer's previous shells; an empty list removes
    /// them). Every layer starts with none, which is the valence-only
    /// transport of before. The layer's inelastic table is then its
    /// **valence** channel: it must not already contain these shells (build
    /// both with
    /// [`crate::electron::inelastic::table::build_shell_channel_tables`]).
    /// The channels are sampled in the given order (module docs, "Inner-shell
    /// channels") and recorded in the [`LayerMetadata`] of the run.
    ///
    /// Fails, leaving nothing changed, if `layer` is out of range, if a
    /// shell's table names a different material from the layer's inelastic
    /// table (tables built for another target), if a `(Z, subshell)` repeats,
    /// or, with band parameters, if a binding energy does not exceed the
    /// layer's Fermi energy (an inner shell lies below the Fermi level).
    pub fn with_inner_shells(
        mut self,
        layer: usize,
        shells: Vec<ShellChannelTable>,
    ) -> Result<Self, TransportError> {
        if layer >= self.shells.len() {
            return invalid(
                "layer",
                format!(
                    "{layer} is out of range (the stack has {} layers)",
                    self.shells.len()
                ),
            );
        }
        let material = self.tables[layer].inelastic.material();
        for (i, s) in shells.iter().enumerate() {
            let name = format!(
                "layer {layer}: inner shell Z = {} {}",
                s.z(),
                s.subshell().label()
            );
            if s.table().material() != material {
                return invalid(
                    "inner-shell channel",
                    format!(
                        "{name}: its table is for {:?}, the layer's inelastic table for {material:?}",
                        s.table().material()
                    ),
                );
            }
            if shells[..i]
                .iter()
                .any(|o| o.z() == s.z() && o.subshell() == s.subshell())
            {
                return invalid("inner-shell channel", format!("{name} listed twice"));
            }
            if let Some(bands) = &self.bands {
                let ef = bands[layer].fermi_ev();
                if s.binding_energy_ev() <= ef {
                    return invalid(
                        "inner-shell channel",
                        format!(
                            "{name}: binding energy {} eV does not exceed the Fermi energy \
                             {ef} eV",
                            s.binding_energy_ev()
                        ),
                    );
                }
            }
        }
        self.shells[layer] = shells;
        Ok(self)
    }

    /// The inner-shell channel tables of each layer.
    pub fn inner_shells(&self) -> &[Vec<ShellChannelTable>] {
        &self.shells
    }

    /// The configuration.
    pub fn config(&self) -> &TransportConfig {
        &self.config
    }

    /// The tables of each layer.
    pub fn layer_tables(&self) -> &[LayerTables] {
        &self.tables
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

    /// The largest inelastic loss an electron in state `st` may take: `E`
    /// without a secondary model, `E - E_F` (the Fermi-level clamp) with one.
    fn max_loss_ev(&self, st: &ElectronState) -> f64 {
        match (&self.config.secondaries, &self.bands) {
            (SecondaryModel::KieftBosch(_), Some(b)) => st.energy_ev - b[st.layer].fermi_ev(),
            _ => st.energy_ev,
        }
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

    /// The region and material a primary starts in: at the geometry's
    /// [`Geometry::entry_point`].
    fn start(&self) -> Result<(usize, usize), TransportError> {
        match self.geometry.locate(self.geometry.entry_point()) {
            Some(region) => Ok((region, self.geometry.material_index(region))),
            None => invalid("entry point", "lies outside the target"),
        }
    }

    fn check_primary(&self, p: &Primary) -> Result<[f64; 3], TransportError> {
        let (_, m0) = self.start()?;
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
            p.energy_ev + self.inner[m0]
        } else {
            p.energy_ev
        };
        if inside <= self.threshold[m0] {
            return invalid(
                "primary energy",
                format!(
                    "{inside} eV inside the first layer must exceed its stopping threshold {} eV",
                    self.threshold[m0]
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
            beam: BeamSpread::default(),
            layers: self
                .extents
                .iter()
                .zip(&self.tables)
                .zip(&self.channels)
                .enumerate()
                .map(|(index, ((&(front_m, back_m), t), c))| LayerMetadata {
                    index,
                    front_m,
                    back_m,
                    elastic_model: t.elastic.model().to_string(),
                    elastic_provenance: t.elastic.provenance().to_string(),
                    inelastic_model: t.inelastic.model().to_string(),
                    inelastic_provenance: t.inelastic.provenance().to_string(),
                    band_structure: self.bands.as_ref().map(|b| b[index].clone()),
                    phonon: c.phonon.clone(),
                    polaron: c.polaron.clone(),
                    inner_shells: self.shells[index]
                        .iter()
                        .map(|s| InnerShellMetadata {
                            z: s.z(),
                            subshell: s.subshell().label().to_string(),
                            binding_energy_ev: s.binding_energy_ev(),
                            binding_provenance: s.binding_provenance().to_string(),
                            model: s.table().model().to_string(),
                            provenance: s.table().provenance().to_string(),
                        })
                        .collect(),
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
        self.run_beam(
            seed,
            n_histories,
            chunk_size,
            primary,
            &BeamSpread::default(),
            new_tally,
        )
    }

    /// [`Self::run`] with a finite beam: each history starts at a position
    /// and direction drawn from `beam` ([`BeamSpread`]).
    pub fn run_beam<T, N>(
        &self,
        seed: u64,
        n_histories: u64,
        chunk_size: u64,
        primary: &Primary,
        beam: &BeamSpread,
        new_tally: N,
    ) -> Result<TransportRun<T>, TransportError>
    where
        T: ElectronTally,
        N: Fn() -> T + Sync,
    {
        self.check_primary(primary)?;
        beam.validate()?;
        let tally = run_particles(
            seed,
            n_histories,
            chunk_size,
            new_tally,
            |t, rng, i| {
                self.history_beam(t, rng, i, primary, beam)
                    .expect("primary validated before the run");
            },
            |total, part| total.merge(part),
        );
        let mut metadata = self.metadata(seed, n_histories, chunk_size, primary);
        metadata.beam = *beam;
        Ok(TransportRun { tally, metadata })
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
        self.history_beam(tally, rng, index, primary, &BeamSpread::default())
    }

    /// [`Self::history`] for a finite beam: the start position, then the
    /// direction, are drawn from `rng` before the transport (no draws for the
    /// default [`BeamSpread`]).
    pub fn history_beam<T, R>(
        &self,
        tally: &mut T,
        rng: &mut R,
        index: u64,
        primary: &Primary,
        beam: &BeamSpread,
    ) -> Result<Fate, TransportError>
    where
        T: ElectronTally,
        R: Rng,
    {
        let dir = self.check_primary(primary)?;
        let [y, z] = beam.sample_position(rng);
        let dir = beam.sample_direction(dir, rng)?;
        // The beam's lateral landing point is an offset from the
        // geometry's entry point, and the start region is the one
        // containing the landing point itself.
        let entry = self.geometry.entry_point();
        let pos = [entry[0], entry[1] + y, entry[2] + z];
        let Some(region) = self.geometry.locate(pos) else {
            return invalid("beam landing point", "lies outside the target");
        };
        let material = self.geometry.material_index(region);
        let mut st = ElectronState {
            pos,
            dir,
            energy_ev: primary.energy_ev,
            layer: material,
        };
        tally.begin_history(index, &st);
        let mut pending: Vec<Pending> = Vec::new();
        let fate = match self.enter(tally, rng, &mut st) {
            Some(f) => f,
            None => self.follow(tally, rng, &mut st, region, 0, &mut pending),
        };
        while let Some((mut s, region, generation)) = pending.pop() {
            tally.begin_secondary(&s, generation);
            let f = self.follow(tally, rng, &mut s, region, generation, &mut pending);
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
        let du = self.inner[st.layer];
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
        mut region: usize,
        generation: u32,
        pending: &mut Vec<Pending>,
    ) -> Fate {
        let mut events = 0u64;
        // Rates and grid positions of the layer's inner-shell channels (an
        // empty Vec does not allocate, so a valence-only run pays nothing).
        let mut shell_rates: Vec<(f64, GridPos)> = Vec::new();
        loop {
            if events >= self.config.max_events {
                return Fate::EventCap;
            }
            let tabs = &self.tables[st.layer];
            let chans = &self.channels[st.layer];
            let (el, el_at) = rate(&tabs.elastic, st.energy_ev);
            let (inel, inel_at) = rate(&tabs.inelastic, st.energy_ev);
            // Diagnostics only: no draw, no change of state.
            tally.table_lookup(
                st,
                TableChannel::Elastic,
                GridCoverage::of(&tabs.elastic, st.energy_ev),
            );
            tally.table_lookup(
                st,
                TableChannel::Inelastic,
                GridCoverage::of(&tabs.inelastic, st.energy_ev),
            );
            let (em, ab, tr) = chans.rates(st.energy_ev);
            let mut total = el + inel + em + ab + tr;
            // Inner-shell channels, after the others (module docs). A
            // channel is closed when the primary cannot pay its binding
            // energy: the largest loss allowed is `E` without a secondary
            // model and `E - E_F` with one (the Fermi-level clamp below).
            let shells = &self.shells[st.layer];
            shell_rates.clear();
            if !shells.is_empty() {
                let max_loss = self.max_loss_ev(st);
                for s in shells {
                    let (r, at) = rate(s.table(), st.energy_ev);
                    let r = if max_loss > s.binding_energy_ev() {
                        r
                    } else {
                        0.0
                    };
                    total += r;
                    shell_rates.push((r, at));
                }
            }
            let mu = st.dir[0];
            let free = if total > 0.0 {
                -(1.0 - uniform(rng)).ln() / total
            } else {
                f64::INFINITY
            };
            // The first event of the flight within the free path (a face
            // reached exactly at the free path counts as reached). With an
            // unbounded free path the limit is the largest finite length, so
            // a face at infinity (the back of a substrate) is no event.
            let flight = self
                .geometry
                .flight(region, st.pos, st.dir, free.min(f64::MAX));
            let exit = match flight {
                Flight::Event(e) => Some(e),
                Flight::Clear { .. } if free.is_infinite() => {
                    tally.stopped(st);
                    return Fate::Trapped;
                }
                Flight::Clear { region: end } => {
                    region = end;
                    None
                }
            };
            let from = st.pos;
            let Some(exit) = exit else {
                for k in 0..3 {
                    st.pos[k] += free * st.dir[k];
                }
                tally.step(from, st, free);
                events += 1;
                let x = uniform(rng) * total;
                let base = [el, inel, em, ab, tr];
                match choose(x, base.into_iter().chain(shell_rates.iter().map(|r| r.0))) {
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
                            SecondaryModel::KieftBosch(opts) => {
                                let band =
                                    &self.bands.as_ref().expect("checked in build")[st.layer];
                                // The primary cannot end below the Fermi level:
                                // Nebula clamps its loss table to `K - E_F`
                                // (`kieft_inelastic::create`, commit named in
                                // `electron::boundary`). The threshold exceeds
                                // E_F, so the bound is positive.
                                let w = w.min(st.energy_ev - band.fermi_ev());
                                let out = kieft_bosch(opts, band, st, w, threshold, rng);
                                st.energy_ev -= w;
                                st.dir = out.primary_dir;
                                tally.inelastic(st, w);
                                spawn_secondary(tally, st, &out, region, generation, pending);
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
                    4 => {
                        tally.polaron_trapped(st);
                        return Fate::PolaronTrapped;
                    }
                    k => {
                        let j = k - 5;
                        let sh = &shells[j];
                        let at = shell_rates[j].1;
                        let b = sh.binding_energy_ev();
                        let shell = InnerShell::from(sh);
                        let threshold = self.threshold[st.layer];
                        // The channel is open, so `b < max_loss`, and the
                        // loss lies in `[b, max_loss]` (module docs).
                        let max_loss = self.max_loss_ev(st);
                        let u = uniform(rng);
                        let w = sample(sh.table(), at, u).max(b).min(max_loss);
                        match self.config.secondaries {
                            SecondaryModel::Off => {
                                st.energy_ev -= w;
                                tally.inelastic(st, w);
                                tally.inner_shell(st, &shell, w);
                            }
                            SecondaryModel::KieftBosch(opts) => {
                                let band =
                                    &self.bands.as_ref().expect("checked in build")[st.layer];
                                let out = kieft_bosch_shell(opts, band, b, st, w, threshold, rng);
                                st.energy_ev -= w;
                                st.dir = out.primary_dir;
                                tally.inelastic(st, w);
                                tally.inner_shell(st, &shell, w);
                                spawn_secondary(tally, st, &out, region, generation, pending);
                            }
                        }
                        if st.energy_ev < threshold {
                            tally.stopped(st);
                            return Fate::Stopped;
                        }
                    }
                }
                continue;
            };
            // Reach the face exactly; the path is redrawn afterwards.
            let to_face = exit.distance.max(0.0);
            st.pos = exit.at;
            tally.step(from, st, to_face);
            let (to, boundary) = match exit.outcome {
                ExitOutcome::Enter { region: to, .. } => {
                    (Some(to), Boundary::Interface { from: region, to })
                }
                ExitOutcome::Escape { face, .. } => match face {
                    Face::Front => (None, Boundary::Surface(Face::Front)),
                    Face::Side => (None, Boundary::Surface(Face::Side)),
                    Face::Back if self.config.escape_rule == EscapeRule::BothFaces => {
                        (None, Boundary::Surface(Face::Back))
                    }
                    Face::Back => {
                        tally.absorbed(st);
                        return Fate::Absorbed;
                    }
                },
            };
            let from_material = st.layer;
            if let Some((quantum, refraction)) = self.step_barrier() {
                let du = to.map_or(0.0, |t| self.inner[self.geometry.material_index(t)])
                    - self.inner[from_material];
                // A draw only if the electron has the normal energy to get
                // over (cross_step ignores it otherwise).
                let can = st.energy_ev * mu * mu + du > 0.0;
                let u = if can { uniform(rng) } else { 1.0 };
                match cross_step(st.dir, st.energy_ev, du, quantum, refraction, u) {
                    StepOutcome::Transmitted { dir, energy_ev } => {
                        st.dir = dir;
                        st.energy_ev = energy_ev;
                        if let Some(t) = to {
                            region = t;
                            st.layer = self.geometry.material_index(t);
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
                region = t;
                st.layer = self.geometry.material_index(t);
            }
            // A periodic crossing re-enters at the opposite face; the
            // segment above ended on the face it left through.
            if let ExitOutcome::Enter { pos, .. } = exit.outcome {
                st.pos = pos;
            }
            match boundary {
                Boundary::Interface { from, to } => {
                    tally.interface(st, from, to);
                    if st.energy_ev < self.threshold[st.layer] {
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

/// Hand the secondary of one inelastic event to the tally and, if it is to be
/// followed, onto the queue: the tail shared by the valence and inner-shell
/// [`SecondaryModel::KieftBosch`] arms of `follow`.
///
/// `st` is the primary **after** its loss and deflection, as every tally hook
/// takes it; the secondary starts at the primary's position and in its layer
/// (no delocalisation, [`crate::electron::secondary`] module docs). One
/// `pending.push` for both channels keeps the queue in a single order — the
/// order the thread-count determinism invariant depends on.
fn spawn_secondary<T: ElectronTally>(
    tally: &mut T,
    st: &ElectronState,
    out: &Outcome,
    region: usize,
    generation: u32,
    pending: &mut Vec<Pending>,
) {
    let created = out.secondary.map(|(dir, energy_ev)| ElectronState {
        pos: st.pos,
        dir,
        energy_ev,
        layer: st.layer,
    });
    tally.secondary(st, &out.event, created.as_ref());
    if let Some(s) = created {
        pending.push((s, region, generation + 1));
    }
}
