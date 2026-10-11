//! Amorphous binary-collision-approximation (BCA) transport in a layered or
//! voxelised target, with full recoil cascades.
//!
//! # Model
//!
//! The target is any [`Geometry`]: a [`Stack`](crate::geometry::Stack) of
//! layers or a [`VoxelGrid`](crate::geometry::VoxelGrid), whose regions are
//! homogeneous, structureless (amorphous). A particle alternates straight
//! free flights and binary elastic collisions with target atoms, losing
//! energy continuously to electrons along each flight. This is the
//! amorphous-target BCA of J. P. Biersack and L. G. Haggmark, Nucl. Instrum.
//! Methods 174 (1980) 257 (TRIM), within the general BCA framework of M. T.
//! Robinson and I. M. Torrens, Phys. Rev. B 9 (1974) 5008 and the treatment
//! in W. Eckstein, *Computer Simulation of Ion-Solid Interactions*
//! (Springer, 1991). Everything here is implemented from those
//! publications; no code from the Tier B/C programs named in
//! `CONTRIBUTING.md` was consulted.
//!
//! One history:
//!
//! 1. The primary enters at the geometry's entry point (the origin of the
//!    front face of a stack; the centre of the front face of a voxel grid
//!    unless [`Bca::with_entry_point`] chooses another) with the beam's
//!    direction.
//! 2. **Free flight** of length `tau * lambda` ([`MeanFreePath`]), with the
//!    electronic loss `N S_e(E) s` taken at the energy at the start of the
//!    segment (Bragg additivity over the layer composition,
//!    [`crate::ion::stopping::bragg`]).
//! 3. **Collision** with a partner drawn by stoichiometry, at impact parameter
//!    `p = p_max sqrt(R)` (uniform over the disc of radius `p_max`) and a
//!    uniform azimuth. The centre-of-mass angle comes from the
//!    [`ScatteringTable`]; lab angles and the transfer `T` from
//!    [`kinematics`].
//! 4. **Recoil**: if `T > E_d` (the partner's displacement energy) the atom is
//!    displaced with energy `T - E_b` (the lattice binding `E_b` stays in the
//!    lattice). Otherwise `T` stays in the lattice at the site. This is the
//!    displacement criterion of Biersack and Haggmark (1980) and Eckstein
//!    (1991). Displaced atoms are followed in turn as full cascades
//!    ([`BcaConfig::follow_recoils`]), from an explicit stack (no recursion),
//!    and draw from the same per-history random stream, so the stream is
//!    keyed on the primary index only and results do not depend on the thread
//!    count.
//! 5. A particle stops when its energy falls below its **cutoff**
//!    ([`BcaConfig::primary_cutoff_ev`], [`BcaConfig::recoil_cutoff_ev`]),
//!    or leaves through the front or back face if it can overcome the planar
//!    surface barrier ([`kinematics::refract_out`]); otherwise it is reflected
//!    specularly back into the target.
//!
//! # Free-path convention and layer boundaries
//!
//! [`Particle::layer`] is the particle's *region* in the geometry: the layer
//! index in a stack, the flat voxel index in a voxel grid. Material data
//! (scattering, stopping, energies) is cached per material, shared by all
//! voxels of that material. The surface barrier of an escaping atom uses the
//! `E_s` of the material it leaves from and the actual outward normal of the
//! face (see [`kinematics::refract_out_normal`]); a particle that cannot
//! escape is reflected specularly about that normal. Escapes through the
//! lateral faces of a voxel grid are reported as [`Face::Side`] and counted
//! in [`EnergyBudget::lateral`]. [`crate::tally::ion::IonTally`] and the CLI
//! tallies describe stacks only; use [`SummaryTally`] or a custom
//! [`BcaTally`] for voxel runs.
//!
//! A flight is drawn as a dimensionless number of mean free paths `tau` and
//! converted to a length with the local `lambda`. When a flight reaches a
//! change of material (or a periodic wrap of a voxel grid) it is truncated
//! exactly there, the electronic loss for the truncated length is applied,
//! the material is switched, and the flight continues in the new region
//! with the **unused** part `tau - s/lambda` (the flight is not redrawn). A
//! voxel face between voxels of one material is not an event; the region a
//! flight ends in is still tracked ([`Geometry::flight`]), so the collision,
//! its tallies and its recoils carry the region of the collision site. The
//! collision happens in the material where the flight ends, with a partner
//! and `p_max` from that material. For exponential paths this is exact by
//! memorylessness; for the constant path it means a layer split of one
//! material changes nothing but floating-point rounding,
//! which the `split_layer_is_equivalent` test checks.
//!
//! # Weak collisions (optional)
//!
//! With the constant free path every flight ends in one collision with
//! `p <= p_max`, and the nuclear loss of collisions with `p > p_max` is
//! dropped, while the electronic loss of the flight is charged in full. At
//! cascade-tail energies (a few to tens of eV) the dropped part is large: for
//! Cu on Cu with the ZBL potential the disc `p <= p_max` carries 45 % of the
//! nuclear stopping cross section at 5 eV and 54 % at 10 eV, so the
//! electronic share of a cascade comes out far too high (issue #64; the level-1
//! check `damage.cascade_electronic_share` measures it).
//!
//! [`BcaConfig::weak_collisions`] `= K` adds the "weak" collisions of
//! W. Moller and W. Eckstein, *TRIDYN - Binary collision simulation of atomic
//! collisions and dynamic composition changes in solids*, report IPP 9/64,
//! Max-Planck-Institut fur Plasmaphysik, Garching (1988) (the long write-up
//! of Comput. Phys. Commun. 51 (1988) 355), read at
//! <https://pure.mpg.de/rest/items/item_2131703/component/file_2131702/content>.
//! Read only the report body, PDF pp. 1-47: never Appendix 1 (the TRIDYN
//! program listing, PDF pp. 48-86), which is Tier C (see `CONTRIBUTING.md`).
//! Quoted:
//!
//! * p. 14: "'weak' collisions might occur with more distant atoms which
//!   might contribute to energy loss and angular deflection. The present
//!   version allows up to three additional weak collisions with impact
//!   parameters larger than `p_max`, each of them representing one additional
//!   atomic volume", with eq. (9) replaced by eq. (26),
//!   `p_k^weak = p_max sqrt(k + r_p)`, `k = 1, 2, 3` (`r_p` uniform in
//!   `[0, 1)`): the `k`-th partner is uniform over the annulus between
//!   `p_max sqrt(k)` and `p_max sqrt(k + 1)`, of area `pi p_max^2`. Fig. 3
//!   (p. 11) draws the first-order partner at its own azimuth.
//! * p. 26: "the weak collision loop is entered which finally defines new
//!   directions after each of the simultaneous collisions. (Actually, the
//!   last passage of the weak collision loop represents the hard collision.)
//!   A primary recoil may be generated and stored for each hard collision."
//! * p. 35: the program sums the elastic energy transfers of the weak and
//!   hard collisions, and their local inelastic (Oen-Robinson) losses.
//! * Figs. 7 and 8, pp. 29-30 (flow charts of the projectile and cascade
//!   loops): inside the weak-collision loop, after the partner's species,
//!   azimuth and impact parameter are drawn, a "Target beyond Surf.?" test
//!   skips the passage when the partner would lie outside the target; and
//!   "Elastic Energy Loss" is applied once, after the loop.
//!
//! So each collision step here is: the `K` weak collisions, `k = 1..K`, each
//! with a partner drawn by stoichiometry, `p = p_max sqrt(k + R)` and a
//! uniform azimuth, then the hard collision with `p = p_max sqrt(R)`. Each
//! deflects the particle and takes its transfer `T`; under
//! [`ElectronicLoss::EquipartitionLsOr`] each also takes its Oen-Robinson
//! local loss, so the local half is sampled out to `p_max sqrt(K + 1)`. A weak
//! collision never makes a recoil: its `T` stays in the lattice at the site
//! ([`LatticeDeposit::Weak`], counted in [`EnergyBudget::lattice`]). Random
//! draws per weak collision: partner, `R`, azimuth, in that order and before
//! the hard collision's, so `K = 0` draws exactly what the engine drew before
//! the option existed and reproduces it bit for bit.
//!
//! **Surface test.** As in Figs. 7 and 8, a weak collision whose partner would
//! lie in front of the front surface (in vacuum) is skipped; its draws are
//! still made. The partner sits at distance `p` from the path on the side
//! opposite to the deflection (where a hard collision sends its recoil). In
//! TRIDYN the test also covers the hard collision (the loop's last passage);
//! here it does not, because that would change the engine without weak
//! collisions, which this option leaves untouched. Whether the hard
//! collision should be tested too is left as a separate model decision. The
//! back face of a finite target is not tested: the flow charts test "the
//! surface" only.
//!
//! **Deviation from TRIDYN: the energy of each collision.** The flow charts
//! apply the elastic loss once, after the loop, and p. 35 sums the
//! transfers, so in TRIDYN every collision of a step is evaluated at the
//! energy the step starts with. Here each is evaluated at the energy left
//! after the previous one (and in the direction after it), so no transfer
//! can exceed the particle's energy. Evaluated at the starting energy, the
//! summed transfers of a slow atom often exceed its energy, and the report
//! shows no rule for that case: for Ar 1 keV on Cu with `E_d = E_s` and
//! `K = 3`, 18 % of the binary collisions would need a clamp, about 200 eV per
//! ion in all, and the sputter yield would depend on how the clamp is
//! applied. Measured on that problem (20 000 ions, seed 1): 0.882 as here;
//! 0.875 with only the weak collisions at the starting energy (each transfer
//! clamped to the energy left, the hard collision at the energy left); 0.976
//! with every collision at the starting energy and each transfer clamped in
//! loop order. Since the clamp is not in the source, the sequential choice
//! is kept.
//!
//! With `K` weak collisions the nuclear loss per path is the stopping cross
//! section integrated to `p_max sqrt(K + 1)` instead of `p_max` (93 % of the
//! full value at 5 eV for Cu on Cu with `K = 3`). It is available with
//! [`MeanFreePath::Constant`] only (`p_max` is the constant-path radius);
//! [`MeanFreePath::EnergyDependent`] rejects it, since that convention drops
//! small-angle collisions by design.
//!
//! **Default 0**, because there is no single published convention: Biersack
//! and Haggmark (1980), whose constant free path this engine follows, use one
//! collision per flight; TRIDYN allows up to three (the report read does not
//! give a default); SDTrimSP, a TRIDYN descendant, defaults to two for
//! projectiles and two for recoils (`iwc`, `iwcr`: A. Mutzke et al.,
//! *SDTrimSP Version 6.00*, report IPP 2019-02 (2019), Table 13, p. 73,
//! "number of ring cylinders for weak simultaneous collisions"); and the
//! RustBCA manual documents `weak_collision_order` with a default of 0 (its
//! GitHub wiki, "Standalone Code: Input File"; manual only, Tier B).
//!
//! The option is not a free improvement, which is the other reason it is off
//! by default. It brings the electronic share of a cascade in line with the
//! Lindhard partition integral equation solved for the engine's own inputs
//! (level-1 rows `damage.cascade_electronic_share.*`: for Cu on Cu at 1 keV
//! about 0.79 without weak collisions and 0.35 with `K = 3`, against 0.38).
//! But a weak transfer stays in the lattice, and at a few eV it is a large
//! part of the energy (Cu on Cu at 5 eV, ZBL: on average 46, 25 and 13 % in
//! the first three annuli), so slow atoms near the surface stop sooner: the
//! Ar 1 keV on Cu sputter yield with `E_d = E_s` falls from 1.93 to 0.88 with
//! `K = 3` (0.66 without the surface test; `docs/validation.md`, level 2).
//! Use it for the energy partition (damage energy, electronic share); its
//! effect on yields is large and has not been checked against measured
//! yields.
//!
//! # Energy bookkeeping
//!
//! Every energy change is subtracted from the particle and added to exactly
//! one field of an [`EnergyBudget`], so each history conserves energy up to
//! rounding (`EnergyBudget::residual`). Electronic losses larger than the
//! particle's energy are clamped to it.
//!
//! Vacancy, interstitial and replacement counts, and the NRT damage estimate,
//! are kept by the tally [`crate::tally::IonTally`] from the events reported
//! here ([`crate::ion::damage`] has the models).
//!
//! # Crystal regions
//!
//! Everything above is the amorphous partner model. [`Bca::with_crystal`]
//! switches chosen regions to the crystal flight model (partners from an
//! explicit lattice, in [`crystal`]); the other regions, and engines that do
//! not call it, are untouched and bit-identical to the amorphous engine.
//!
//! # Not modelled here (see the issue tracker)
//!
//! Damage accumulation in crystal regions (the lattice there is perfect,
//! static or thermally vibrating; the local Oen-Robinson loss of
//! [`ElectronicLoss::EquipartitionLsOr`] is taken at every lattice partner,
//! see [`crystal`]), target composition changes with fluence, and refraction
//! of the incident beam at the entrance surface (negligible at keV energies).

pub mod crystal;
pub mod kinematics;
pub mod tally;

pub use crystal::{CrystalMetadata, CrystalTarget, Thermal, ThermalMetadata};
pub use tally::{BcaTally, ElectronicChannel, EnergyBudget, Face, LatticeDeposit, SummaryTally};

use std::collections::HashMap;
use std::f64::consts::PI;

use rand_core::Rng;

use crate::geometry::{ExitOutcome, Flight, Geometry};
use crate::ion::crystal::Divergence;
use crate::ion::potential::{Potential, Screening};
use crate::ion::scattering::{closest_approach, theta_quadrature, ScatteringTable};
use crate::ion::stopping::bragg::{bragg_cross_section_per_atom, NoCorrection};
use crate::ion::stopping::mix::EquipartitionMix;
use crate::ion::stopping::{ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::material::EnergyKind;
use crate::rng::{run_particles_range, ParticleRng};
use crate::units::J_PER_EV;

/// Largest [`BcaConfig::weak_collisions`]: "up to three additional weak
/// collisions" (Moller and Eckstein, IPP 9/64 (1988), p. 14).
pub const MAX_WEAK_COLLISIONS: u8 = 3;

/// The incident beam.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beam {
    /// Projectile species (atomic number and mass).
    pub ion: Ion,
    /// Incident kinetic energy, eV.
    pub energy_ev: f64,
    /// Polar angle of incidence from the surface normal, radians, in `[0, pi/2)`.
    pub polar_rad: f64,
    /// Azimuth of incidence about the normal, radians.
    pub azimuth_rad: f64,
    /// Number of primary histories.
    pub count: u64,
}

impl Beam {
    /// Normal incidence.
    pub fn normal(ion: Ion, energy_ev: f64, count: u64) -> Self {
        Self {
            ion,
            energy_ev,
            polar_rad: 0.0,
            azimuth_rad: 0.0,
            count,
        }
    }

    /// Unit direction `[x, y, z]` of incidence.
    pub fn direction(&self) -> [f64; 3] {
        let (s, c) = self.polar_rad.sin_cos();
        let (sa, ca) = self.azimuth_rad.sin_cos();
        [c, s * ca, s * sa]
    }
}

/// Word position (in the 32-bit words of the ChaCha stream) of the segment
/// of a history's random stream that beam-divergence draws come from:
/// `2^65`. The transport draws start at word 0 and never reach it, and it is
/// disjoint from the thermal-displacement segment (`2^66`) and the
/// crystal-shift segment (`2^67`) (each is far shorter than the gap to the
/// next), so enabling divergence leaves every other draw of the history
/// unchanged.
pub const DIVERGENCE_STREAM_WORD: u128 = 1u128 << 65;

/// Most deflected directions drawn for one primary before
/// [`Bca::with_divergence`]'s inward conditioning gives up with
/// [`BcaError::BeamDivergence`]. For the small-angle spreads the input
/// front end accepts, one draw is inward with probability about one half or
/// more, so exhausting 1000 attempts has probability below `2^-1000`; the
/// bound only guarantees that a pathological setting cannot hang a run.
pub const MAX_DIVERGENCE_ATTEMPTS: u32 = 1000;

/// Run metadata of a beam divergence ([`Bca::divergence_metadata`]): the
/// distribution, its width and the incidence policy, so a result states the
/// spread it was run with.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct DivergenceMetadata {
    /// `"gaussian"` or `"uniform-cone"`.
    pub model: &'static str,
    /// Which angle `width_rad` is: `"sigma_per_plane"` (Gaussian standard
    /// deviation of each of two orthogonal plane angles) or
    /// `"cone_half_angle"` (uniform in solid angle).
    pub width_kind: &'static str,
    /// The width, radians.
    pub width_rad: f64,
    /// The width, degrees (the unit the input uses).
    pub width_deg: f64,
    /// Incidence policy: `"inward-conditioned"`, i.e. the distribution
    /// conditioned on directions into the target by rejection sampling.
    pub incidence: &'static str,
    /// Rejection-sampling attempt bound per primary.
    pub max_attempts: u32,
    /// Word position of the random-stream segment the draws come from.
    pub stream_word: &'static str,
}

impl DivergenceMetadata {
    /// Metadata of `divergence`; `None` for [`Divergence::None`].
    pub fn new(divergence: &Divergence) -> Option<Self> {
        let (model, width_kind, width_rad) = match *divergence {
            Divergence::None => return None,
            Divergence::Gaussian { sigma_rad } => ("gaussian", "sigma_per_plane", sigma_rad),
            Divergence::UniformCone { half_angle_rad } => {
                ("uniform-cone", "cone_half_angle", half_angle_rad)
            }
        };
        Some(Self {
            model,
            width_kind,
            width_rad,
            width_deg: width_rad.to_degrees(),
            incidence: "inward-conditioned",
            max_attempts: MAX_DIVERGENCE_ATTEMPTS,
            stream_word: "2^65",
        })
    }
}

/// Why a single history failed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum HistoryError {
    /// The stopping model failed.
    #[error(transparent)]
    Stopping(#[from] StoppingError),
    /// No inward direction was drawn from the beam divergence within
    /// [`MAX_DIVERGENCE_ATTEMPTS`] attempts.
    #[error("no inward beam direction in {attempts} divergence draws")]
    BeamDivergence {
        /// Attempts made.
        attempts: u32,
    },
}

/// Draw `divergence` about `central` until the direction points into the
/// target (`dir[0] > 0`), at most [`MAX_DIVERGENCE_ATTEMPTS`] times. The
/// result is the divergence law conditioned on inward incidence.
fn sample_inward<R: Rng + ?Sized>(
    divergence: &Divergence,
    central: [f64; 3],
    rng: &mut R,
) -> Result<[f64; 3], HistoryError> {
    for _ in 0..MAX_DIVERGENCE_ATTEMPTS {
        let d = divergence.sample_direction(central, rng);
        if d[0] > 0.0 {
            return Ok(d);
        }
    }
    Err(HistoryError::BeamDivergence {
        attempts: MAX_DIVERGENCE_ATTEMPTS,
    })
}

/// How free-flight lengths and the impact-parameter limit are chosen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MeanFreePath {
    /// Constant flight length `l = N^(-1/3)` (the mean interatomic distance)
    /// and impact parameters uniform over a disc of radius
    /// `p_max = (pi N^(2/3))^(-1/2)`, so that `N pi p_max^2 l = 1`: each flight
    /// sweeps exactly one atom's worth of target. This is the amorphous-target
    /// convention of Biersack and Haggmark, Nucl. Instrum. Methods 174 (1980)
    /// 257. To avoid every primary making its first collision at the same depth
    /// `l cos(theta_in)`, the primary's first flight is `R l` with `R` uniform
    /// in `[0, 1)`; recoils start at an atom site and fly a full `l`.
    /// Collisions with `p > p_max` are dropped unless
    /// [`BcaConfig::weak_collisions`] adds them (module docs, "Weak
    /// collisions").
    Constant,
    /// Energy-dependent free path with Poisson-distributed flight lengths.
    /// Collisions deflecting by less than `min_cm_angle_rad` (centre-of-mass)
    /// are neglected: for each element `i`, `p_i(E)` is the impact parameter at
    /// which the centre-of-mass angle equals that minimum, capped at the
    /// constant-convention `p_max` (so the free path is never shorter than
    /// `N^(-1/3)`). The mean free path is `lambda = 1 / (N pi Sum_i x_i p_i^2)`,
    /// the flight length is exponential with that mean, the partner is drawn
    /// with probability `x_i p_i^2 / Sum_j x_j p_j^2`, and `p = p_i sqrt(R)`.
    /// This is the standard cross-section cut-off treatment of a Poisson
    /// collision process (Eckstein 1991). It reduces to `Constant` (with
    /// exponential instead of fixed flight lengths) when every `p_i` hits the
    /// cap. The nuclear energy loss of the neglected small-angle collisions is
    /// dropped, so choose the minimum angle small; and since `p_i` never
    /// exceeds the constant-convention `p_max`, collisions beyond it are
    /// dropped too, as in `Constant` without weak collisions (which this
    /// variant does not support).
    EnergyDependent {
        /// Smallest centre-of-mass deflection treated as a collision, radians,
        /// in `(0, pi)`.
        min_cm_angle_rad: f64,
    },
}

/// Electronic energy-loss treatment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElectronicLoss {
    /// All electronic loss is continuous (nonlocal) along the free flight,
    /// from the stopping model passed to [`Bca::new`].
    NonLocal,
    /// Equipartition of Lindhard-Scharff stopping: half nonlocal along the
    /// flight, half local at each collision by the Oen-Robinson
    /// impact-parameter-dependent loss evaluated at the distance of closest
    /// approach (O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132
    /// (1976) 647; [`EquipartitionMix`]). The stopping model passed to
    /// [`Bca::new`] is not used in this mode. The local part is only sampled
    /// at the collisions, out to `p_max` (to `p_max sqrt(K + 1)` with
    /// [`BcaConfig::weak_collisions`] `= K`, each weak collision taking its
    /// own local loss), so where that radius is not large compared with
    /// `a / 0.3` (the decay length of the Oen-Robinson loss) the total
    /// electronic stopping is somewhat below the Lindhard-Scharff value.
    EquipartitionLsOr,
}

/// Run configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BcaConfig {
    /// Free-path convention. Default [`MeanFreePath::Constant`].
    pub mean_free_path: MeanFreePath,
    /// Electronic-loss treatment. Default [`ElectronicLoss::NonLocal`].
    pub electronic: ElectronicLoss,
    /// Number of simultaneous weak collisions per collision step, `0..=3`
    /// ([`MeanFreePath::Constant`] only; see the module docs, "Weak
    /// collisions"). Default 0, the single-collision convention of Biersack
    /// and Haggmark (1980).
    pub weak_collisions: u8,
    /// Follow displaced atoms as full cascades. If `false`, a displaced atom
    /// stops where it was created and its energy counts as rest energy (an
    /// ion-only run). Default `true`.
    pub follow_recoils: bool,
    /// The primary stops when its energy falls below this, eV.
    pub primary_cutoff_ev: f64,
    /// Recoils stop when their energy falls below this, eV. A recoil cutoff
    /// above the surface binding energies suppresses sputtering, so keep it
    /// below the smallest `E_s` when sputtering matters (Eckstein 1991).
    pub recoil_cutoff_ev: f64,
    /// Planar surface barrier for the primary's species at both faces, eV.
    /// Default 0 (backscattered ions are not bound). Target atoms use the
    /// `E_s` of the material at the face (see [`Bca::new`]).
    pub primary_surface_binding_ev: f64,
    /// Run seed for [`crate::rng::stream`].
    pub seed: u64,
    /// Histories per work chunk for [`crate::rng::run_particles`]; fixes the summation
    /// order, so keep it independent of the thread count. Default 64.
    pub chunk_size: u64,
}

impl BcaConfig {
    /// Defaults with the two cutoffs, which have no defensible default and
    /// must be chosen for the problem.
    pub fn new(primary_cutoff_ev: f64, recoil_cutoff_ev: f64) -> Self {
        Self {
            mean_free_path: MeanFreePath::Constant,
            electronic: ElectronicLoss::NonLocal,
            weak_collisions: 0,
            follow_recoils: true,
            primary_cutoff_ev,
            recoil_cutoff_ev,
            primary_surface_binding_ev: 0.0,
            seed: 0,
            chunk_size: 64,
        }
    }
}

/// Errors from setting up or running the engine.
#[derive(Debug, thiserror::Error)]
pub enum BcaError {
    /// A beam parameter is out of range.
    #[error("invalid beam: {0}")]
    InvalidBeam(String),
    /// A configuration parameter is out of range.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),
    /// A material has an energy with no value; set it before running.
    #[error("layer {layer}: {kind} for element Z={z} is not set")]
    EnergyNotSet {
        /// Material index (the layer index for a stack).
        layer: usize,
        /// Atomic number.
        z: u8,
        /// Which energy.
        kind: EnergyKind,
    },
    /// The stopping model failed during a history (for example a user table
    /// queried outside its range). `index` is the lowest failing history.
    #[error("history {index}: {source}")]
    Stopping {
        /// Primary index.
        index: u64,
        /// The model error.
        source: StoppingError,
    },
    /// Beam divergence produced no inward direction within the attempt
    /// bound ([`MAX_DIVERGENCE_ATTEMPTS`]). `index` is the lowest failing
    /// history.
    #[error("history {index}: no inward beam direction in {attempts} divergence draws")]
    BeamDivergence {
        /// Primary index.
        index: u64,
        /// Attempts made.
        attempts: u32,
    },
}

/// A moving or stopped particle.
#[derive(Debug, Clone, PartialEq)]
pub struct Particle {
    /// Species index: 0 is the beam species, then the target elements in
    /// order of first appearance in the stack.
    pub species: usize,
    /// Atomic number.
    pub z: u8,
    /// Mass, u.
    pub mass_amu: f64,
    /// Kinetic energy, eV.
    pub energy_ev: f64,
    /// Position `[x, y, z]`, m.
    pub pos: [f64; 3],
    /// Unit direction.
    pub dir: [f64; 3],
    /// Index of the region the particle is in: the layer of a stack, the flat
    /// voxel index of a voxel grid (see [`Geometry`]).
    pub layer: usize,
    /// Region the particle started in: the entry region of a primary, the
    /// region of the displacing collision for a recoil. Unlike `layer` it does
    /// not change as the particle moves, so a tally that sees the particle
    /// leave can tell where the atom came from.
    pub origin_layer: usize,
    /// 0 for the primary, parent's generation + 1 for recoils.
    pub generation: u32,
}

impl Particle {
    /// True for the beam particle of a history.
    pub fn is_primary(&self) -> bool {
        self.generation == 0
    }
}

#[derive(Debug, Clone)]
struct SpeciesData {
    ion: Ion,
    cutoff_ev: f64,
    /// `E_s` used when the exit material does not contain the element: that
    /// of the first material that does (the primary's own binding for the
    /// beam).
    barrier_default_ev: f64,
}

#[derive(Debug, Clone, Copy)]
struct PairData {
    /// Screening length, m.
    a: f64,
    /// Reduced energy per eV of lab energy.
    eps_per_ev: f64,
    gamma: f64,
    /// M1 / M2.
    mu: f64,
}

#[derive(Debug, Clone)]
struct ElemData {
    z: u8,
    species: usize,
    fraction: f64,
    e_d_ev: f64,
    e_b_ev: f64,
}

#[derive(Debug, Clone)]
struct LayerData {
    /// Atom density, m^-3.
    n: f64,
    /// Constant flight length N^(-1/3), m.
    ell: f64,
    /// Constant-convention p_max, m.
    p_const: f64,
    elems: Vec<ElemData>,
    /// Per projectile species: `c` with `Σ_j x_j S_j(E) = c sqrt(E)` over the
    /// layer's elements for the configured nonlocal model, if every element's
    /// stopping has that form (see
    /// [`ElectronicStopping::sqrt_energy_coefficient`]). The fractional
    /// powers in the coefficient depend only on the ion and target, so they
    /// are evaluated here once instead of at every flight.
    loss_coef: Vec<Option<f64>>,
}

/// `ln beta_max` on a grid uniform in `ln eps`, for the energy-dependent
/// free path.
#[derive(Debug, Clone)]
struct BetaMaxGrid {
    ln_eps0: f64,
    d: f64,
    ln_beta: Vec<f64>,
}

impl BetaMaxGrid {
    const PER_DECADE: f64 = 16.0;

    fn build(table: &ScatteringTable, theta_min: f64) -> Self {
        let spec = table.spec();
        let d = std::f64::consts::LN_10 / Self::PER_DECADE;
        let ln_eps0 = spec.eps_min.ln();
        let n = ((spec.eps_max.ln() - ln_eps0) / d).ceil() as usize + 1;
        let (lb0, lb1) = (spec.beta_min.ln(), spec.beta_max.ln());
        let ln_beta = (0..n)
            .map(|i| {
                let eps = (ln_eps0 + i as f64 * d).exp().min(spec.eps_max);
                let th = |lb: f64| table.theta(eps, lb.exp()).unwrap_or(0.0);
                if th(lb1) >= theta_min {
                    return lb1;
                }
                // theta decreases with beta; bisect in ln beta.
                let (mut lo, mut hi) = (lb0, lb1);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if th(mid) >= theta_min {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                0.5 * (lo + hi)
            })
            .collect();
        Self {
            ln_eps0,
            d,
            ln_beta,
        }
    }

    /// `beta_max(eps)`, linear in `ln eps`, clamped at the grid ends.
    fn beta(&self, eps: f64) -> f64 {
        let f = ((eps.ln() - self.ln_eps0) / self.d).max(0.0);
        let last = self.ln_beta.len() - 1;
        let i = (f.floor() as usize).min(last.saturating_sub(1));
        let t = (f - i as f64).min(1.0);
        let j = (i + 1).min(last);
        ((1.0 - t) * self.ln_beta[i] + t * self.ln_beta[j]).exp()
    }
}

/// What [`Bca::binary`] reports: the energy transfer `T`, the direction of
/// the projectile before the collision, and the recoil's lab deflection as
/// `(sin, cos)`.
type BinaryOutcome = (f64, [f64; 3], (f64, f64));

/// Per-flight scratch space for the per-element `p_max` and partner weights.
#[derive(Default)]
struct Scratch {
    p_max: Vec<f64>,
    cum: Vec<f64>,
    /// Crystal flight model: the lattice sites found by the last search, the
    /// sites hit by the previous collision, and the partners of the current
    /// collision step.
    cands: Vec<crystal::Candidate>,
    last: Vec<crystal::SiteId>,
    targets: Vec<crystal::Partner>,
    /// Thermal crystals: the static site of each entry of `cands` (whose
    /// positions are displaced), a sort buffer, and the displacements of the
    /// sites met on the current flight line (`displaced_line`: crystal and
    /// direction bits), with the buffer that replaces them after a search.
    homes: Vec<[f64; 3]>,
    encounter: Vec<(crystal::Candidate, [f64; 3])>,
    displaced: HashMap<crystal::SiteId, [f64; 3]>,
    displaced_next: HashMap<crystal::SiteId, [f64; 3]>,
    displaced_line: Option<(usize, [u64; 3])>,
}

impl Scratch {
    /// Forget the displaced sites: the next crystal search starts a new
    /// encounter.
    fn new_line(&mut self) {
        self.displaced.clear();
        self.displaced_line = None;
    }
}

/// Reusable working memory for [`Bca::history_in`]: the stack of particles
/// still to be followed and the per-flight scratch. A history that runs in
/// buffers which have already seen a cascade as large allocates nothing, so a
/// worker thread allocates only while its buffers grow to the largest cascade
/// it meets, not per history and not per collision.
#[derive(Default)]
pub struct HistoryBuffers {
    pending: Vec<Particle>,
    scratch: Scratch,
    /// Crystal regions: this history's random lattice translation per
    /// crystal (empty when there are none).
    shifts: Vec<[f64; 3]>,
    /// Thermal crystals: this history's stream of site displacements
    /// (`None` when no crystal vibrates).
    thermal_rng: Option<ParticleRng>,
}

// Manual: the per-history scratch space is internal and has no `Debug`.
impl std::fmt::Debug for HistoryBuffers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HistoryBuffers")
            .field("pending", &self.pending.len())
            .field("shifts", &self.shifts.len())
            .finish_non_exhaustive()
    }
}

impl HistoryBuffers {
    /// Empty buffers.
    pub fn new() -> Self {
        Self::default()
    }
}

/// The nonlocal half of [`EquipartitionMix`] as an [`ElectronicStopping`]
/// model, so the Bragg sum can be reused.
struct MixNonLocal(EquipartitionMix);

impl ElectronicStopping for MixNonLocal {
    fn name(&self) -> &'static str {
        "equipartition-ls-or-nonlocal"
    }
    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        self.0.nonlocal_stopping(ion, target_z, energy_ev)
    }
    fn validity(&self, ion: &Ion) -> ValidityRange {
        self.0.validity(ion)
    }
    fn sqrt_energy_coefficient(&self, ion: &Ion, target_z: u8) -> Option<f64> {
        self.0.nonlocal_sqrt_energy_coefficient(ion, target_z)
    }
}

/// The BCA engine for one beam, target and configuration.
///
/// Holds borrowed, read-only inputs and precomputed per-layer and per-pair
/// data, so it is `Sync` and shared by all worker threads.
pub struct Bca<'a> {
    beam: Beam,
    geometry: &'a dyn Geometry,
    entry_pos: [f64; 3],
    entry_region: usize,
    config: BcaConfig,
    stopping: &'a (dyn ElectronicStopping + Sync),
    mix: MixNonLocal,
    table: &'a ScatteringTable,
    screening: Screening,
    species: Vec<SpeciesData>,
    /// Per-material data (indexed by material index, not region).
    layers: Vec<LayerData>,
    /// Surface barrier `E_s` by `material * n_species + species`.
    barriers: Vec<f64>,
    /// `pairs[projectile species * n_species + target species]`.
    pairs: Vec<PairData>,
    beta_max: Option<BetaMaxGrid>,
    /// Crystal partner models ([`Bca::with_crystal`]); empty for an
    /// amorphous target.
    crystals: Vec<crystal::CrystalData>,
    /// Index into `crystals` by region; empty when there are no crystals.
    region_crystal: Vec<Option<usize>>,
    /// Beam divergence about the nominal direction ([`Bca::with_divergence`]);
    /// [`Divergence::None`] (the default) draws nothing.
    divergence: Divergence,
}

fn finite_nonneg(v: f64) -> bool {
    v.is_finite() && v >= 0.0
}

// Manual: the geometry and stopping model are trait objects without `Debug`,
// and the precomputed tables are large.
impl std::fmt::Debug for Bca<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Bca")
            .field("beam", &self.beam)
            .field("entry_pos", &self.entry_pos)
            .field("entry_region", &self.entry_region)
            .field("config", &self.config)
            .field("screening", &self.screening)
            .field("crystals", &self.crystals.len())
            .finish_non_exhaustive()
    }
}

impl<'a> Bca<'a> {
    /// Validate the inputs and precompute per-layer and per-pair data.
    ///
    /// * Every layer material must have `E_d`, `E_b` and `E_s` set for every
    ///   element (checked with [`crate::material::Material::unset_energies`]),
    ///   so a run never fails half way for a missing parameter.
    /// * The scattering table fixes the screening function and length
    ///   convention; the same table serves every projectile/target pair,
    ///   since the angle depends on the reduced variables only.
    /// * Target atoms use, at each face, the `E_s` of their element in the
    ///   material they leave from; an element absent from that material uses
    ///   its `E_s` in the first material that contains it.
    /// * `geometry` is a [`Stack`](crate::geometry::Stack) or a
    ///   [`VoxelGrid`](crate::geometry::VoxelGrid) (or any other
    ///   [`Geometry`]); its default entry point must lie in the target
    ///   ([`Bca::with_entry_point`] changes it).
    pub fn new(
        beam: Beam,
        geometry: &'a dyn Geometry,
        config: BcaConfig,
        stopping: &'a (dyn ElectronicStopping + Sync),
        table: &'a ScatteringTable,
    ) -> Result<Self, BcaError> {
        if !(beam.energy_ev.is_finite() && beam.energy_ev > 0.0) {
            return Err(BcaError::InvalidBeam(format!(
                "energy {} eV must be finite and positive",
                beam.energy_ev
            )));
        }
        if !(beam.polar_rad >= 0.0 && beam.polar_rad < 0.5 * PI) {
            return Err(BcaError::InvalidBeam(format!(
                "polar angle {} rad must be in [0, pi/2)",
                beam.polar_rad
            )));
        }
        if !beam.azimuth_rad.is_finite() {
            return Err(BcaError::InvalidBeam("azimuth must be finite".into()));
        }
        for (name, v) in [
            ("primary_cutoff_ev", config.primary_cutoff_ev),
            ("recoil_cutoff_ev", config.recoil_cutoff_ev),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(BcaError::InvalidConfig(format!(
                    "{name} = {v} must be finite and positive"
                )));
            }
        }
        if !finite_nonneg(config.primary_surface_binding_ev) {
            return Err(BcaError::InvalidConfig(format!(
                "primary_surface_binding_ev = {} must be finite and non-negative",
                config.primary_surface_binding_ev
            )));
        }
        if let MeanFreePath::EnergyDependent { min_cm_angle_rad } = config.mean_free_path {
            if !(min_cm_angle_rad > 0.0 && min_cm_angle_rad < PI) {
                return Err(BcaError::InvalidConfig(format!(
                    "min_cm_angle_rad = {min_cm_angle_rad} must be in (0, pi)"
                )));
            }
            if config.weak_collisions > 0 {
                return Err(BcaError::InvalidConfig(
                    "weak_collisions requires the constant free path".into(),
                ));
            }
        }
        if config.weak_collisions > MAX_WEAK_COLLISIONS {
            return Err(BcaError::InvalidConfig(format!(
                "weak_collisions = {} must be at most {MAX_WEAK_COLLISIONS}",
                config.weak_collisions
            )));
        }
        let n_materials = geometry.n_materials();
        for layer in 0..n_materials {
            if let Some(&(z, kind)) = geometry.material(layer).unset_energies().first() {
                return Err(BcaError::EnergyNotSet { layer, z, kind });
            }
        }
        let entry_pos = geometry.entry_point();
        let entry_region = geometry.locate(entry_pos).ok_or_else(|| {
            BcaError::InvalidBeam(format!("entry point {entry_pos:?} is outside the target"))
        })?;

        // Species: the beam, then target elements in order of appearance.
        let mut species = vec![SpeciesData {
            ion: beam.ion,
            cutoff_ev: config.primary_cutoff_ev,
            barrier_default_ev: config.primary_surface_binding_ev,
        }];
        let mut species_z: Vec<u8> = vec![0];
        let mut layers = Vec::with_capacity(n_materials);
        for mi in 0..n_materials {
            let m = geometry.material(mi);
            let n = m.atom_number_density();
            let mut elems = Vec::with_capacity(m.components().len());
            for c in m.components() {
                let z = c.z();
                let s = match species_z.iter().skip(1).position(|&q| q == z) {
                    Some(i) => i + 1,
                    None => {
                        let e_s = m.surface_binding_energy_ev(z).expect("validated");
                        species.push(SpeciesData {
                            ion: Ion::new(z).expect("element in table"),
                            cutoff_ev: config.recoil_cutoff_ev,
                            barrier_default_ev: e_s,
                        });
                        species_z.push(z);
                        species.len() - 1
                    }
                };
                elems.push(ElemData {
                    z,
                    species: s,
                    fraction: c.atom_fraction(),
                    e_d_ev: m.displacement_energy_ev(z).expect("validated"),
                    e_b_ev: m.lattice_binding_energy_ev(z).expect("validated"),
                });
            }
            layers.push(LayerData {
                n,
                ell: n.powf(-1.0 / 3.0),
                p_const: 1.0 / (PI * n.powf(2.0 / 3.0)).sqrt(),
                elems,
                loss_coef: Vec::new(),
            });
        }
        let mix = MixNonLocal(EquipartitionMix::new());
        {
            let model: &dyn ElectronicStopping = match config.electronic {
                ElectronicLoss::NonLocal => stopping,
                ElectronicLoss::EquipartitionLsOr => &mix,
            };
            for lay in &mut layers {
                lay.loss_coef = species
                    .iter()
                    .map(|sp| {
                        // Same order as the Bragg sum.
                        let mut sum = 0.0;
                        for e in &lay.elems {
                            sum += e.fraction * model.sqrt_energy_coefficient(&sp.ion, e.z)?;
                        }
                        Some(sum)
                    })
                    .collect();
            }
        }
        // Surface barriers by exit material; the beam always uses its own.
        let mut barriers = Vec::with_capacity(n_materials * species.len());
        for mi in 0..n_materials {
            let m = geometry.material(mi);
            for (si, &z) in species_z.iter().enumerate() {
                barriers.push(if si == 0 {
                    species[0].barrier_default_ev
                } else {
                    m.surface_binding_energy_ev(z)
                        .unwrap_or(species[si].barrier_default_ev)
                });
            }
        }

        let tp = table.potential();
        let ns = species.len();
        let mut pairs = Vec::with_capacity(ns * ns);
        for p in &species {
            for t in &species {
                let (z1, z2) = (f64::from(p.ion.z()), f64::from(t.ion.z()));
                let (m1, m2) = (p.ion.mass_amu(), t.ion.mass_amu());
                let pot = Potential::new(tp.screening, z1, z2).with_length(tp.length);
                pairs.push(PairData {
                    a: pot.screening_length(),
                    eps_per_ev: pot.reduced_energy(Potential::cm_energy(J_PER_EV, m1, m2)),
                    gamma: kinematics::gamma(m1, m2),
                    mu: m1 / m2,
                });
            }
        }
        let beta_max = match config.mean_free_path {
            MeanFreePath::Constant => None,
            MeanFreePath::EnergyDependent { min_cm_angle_rad } => {
                Some(BetaMaxGrid::build(table, min_cm_angle_rad))
            }
        };
        Ok(Self {
            beam,
            geometry,
            entry_pos,
            entry_region,
            config,
            stopping,
            mix,
            table,
            screening: tp.screening,
            species,
            layers,
            barriers,
            pairs,
            beta_max,
            crystals: Vec::new(),
            region_crystal: Vec::new(),
            divergence: Divergence::None,
        })
    }

    /// Give the beam an angular spread about its nominal direction (opt-in;
    /// the default is [`Divergence::None`], which changes nothing).
    ///
    /// Each primary's initial direction is the nominal
    /// [`Beam::direction`] deflected by [`Divergence::sample_direction`] and
    /// *conditioned on pointing into the target* (positive depth component):
    /// a deflected direction with `dir[0] <= 0` is rejected and redrawn, up
    /// to [`MAX_DIVERGENCE_ATTEMPTS`] attempts, after which the history
    /// fails with [`BcaError::BeamDivergence`]. The draws come from a copy of
    /// the history's stream positioned at word `2^65`
    /// ([`DIVERGENCE_STREAM_WORD`]), a segment disjoint from the transport
    /// draws, the thermal segment (`2^66`) and the crystal-shift segment
    /// (`2^67`), so the spread changes no other draw of the history. The
    /// spread is a property of the beam: it applies to the primary from the
    /// entry point, including while it crosses an amorphous layer; recoils
    /// keep their collision-generated directions. The nominal orientation
    /// convention (and the crystal orientation built from it) is unchanged.
    ///
    /// Errors with [`BcaError::InvalidBeam`] when the parameters fail
    /// [`Divergence::validate`].
    pub fn with_divergence(mut self, divergence: Divergence) -> Result<Self, BcaError> {
        divergence
            .validate()
            .map_err(|e| BcaError::InvalidBeam(format!("divergence: {e}")))?;
        self.divergence = divergence;
        Ok(self)
    }

    /// Metadata of the beam divergence; `None` without one.
    pub fn divergence_metadata(&self) -> Option<DivergenceMetadata> {
        DivergenceMetadata::new(&self.divergence)
    }

    /// The beam divergence ([`Divergence::None`] unless set by
    /// [`Bca::with_divergence`]).
    pub fn divergence(&self) -> Divergence {
        self.divergence
    }

    /// The initial direction of primary `index` of this engine's run: what
    /// [`Bca::run`] gives that history (the stream is
    /// [`crate::rng::stream`]`(seed, index)`). The nominal direction without
    /// divergence; otherwise the inward-conditioned divergent sample.
    pub fn primary_direction(&self, index: u64) -> Result<[f64; 3], BcaError> {
        self.initial_direction(&crate::rng::stream(self.config.seed, index))
            .map_err(|e| match e {
                HistoryError::BeamDivergence { attempts } => {
                    BcaError::BeamDivergence { index, attempts }
                }
                HistoryError::Stopping(source) => BcaError::Stopping { index, source },
            })
    }

    /// The direction of primary history `rng` starts with: the nominal beam
    /// direction, or the inward-conditioned divergent sample (see
    /// [`Bca::with_divergence`]).
    fn initial_direction(&self, rng: &ParticleRng) -> Result<[f64; 3], HistoryError> {
        let central = self.beam.direction();
        if matches!(self.divergence, Divergence::None) {
            return Ok(central);
        }
        let mut brng = rng.clone();
        brng.set_word_pos(DIVERGENCE_STREAM_WORD);
        sample_inward(&self.divergence, central, &mut brng)
    }

    /// Start primaries at `pos` instead of the geometry's default entry point.
    /// `pos` must lie in the target (for a voxel grid this selects the
    /// incident voxel explicitly; on a shared face the voxel with the higher
    /// index owns the point).
    pub fn with_entry_point(mut self, pos: [f64; 3]) -> Result<Self, BcaError> {
        self.entry_region = self.geometry.locate(pos).ok_or_else(|| {
            BcaError::InvalidBeam(format!("entry point {pos:?} is outside the target"))
        })?;
        self.entry_pos = pos;
        Ok(self)
    }

    /// The beam.
    pub fn beam(&self) -> &Beam {
        &self.beam
    }

    /// The configuration.
    pub fn config(&self) -> &BcaConfig {
        &self.config
    }

    /// Atomic numbers of the particle species, by species index (index 0 is
    /// the beam).
    pub fn species_z(&self) -> Vec<u8> {
        self.species.iter().map(|s| s.ion.z()).collect()
    }

    /// Run all `beam.count` histories in parallel on the current rayon pool
    /// and return the merged tally. Bit-identical at any thread count when the
    /// tally's `merge` is order-independent of scheduling (it is called in
    /// chunk order). If a history fails, the error of the lowest failing
    /// index is returned.
    pub fn run<T, N>(&self, new_tally: N) -> Result<T, BcaError>
    where
        T: BcaTally,
        N: Fn() -> T + Sync,
    {
        self.run_range(0, self.beam.count, new_tally)
    }

    /// [`Bca::run`] for the global primary indices
    /// `first_index..first_index + count` (instead of `0..beam.count`):
    /// primary `first_index + k` uses the stream
    /// `stream(seed, first_index + k)` and reports that index to the tally
    /// hooks. A run cut into consecutive ranges draws exactly the histories of
    /// one longer run; this is what lets the fluence steps of a dynamic run
    /// continue the global stream. `run` is `run_range(0, beam.count, ..)`.
    pub fn run_range<T, N>(&self, first_index: u64, count: u64, new_tally: N) -> Result<T, BcaError>
    where
        T: BcaTally,
        N: Fn() -> T + Sync,
    {
        struct Acc<T> {
            tally: T,
            err: Option<(u64, HistoryError)>,
            buffers: HistoryBuffers,
        }
        let acc = run_particles_range(
            self.config.seed,
            first_index,
            count,
            self.config.chunk_size,
            || Acc {
                tally: new_tally(),
                err: None,
                buffers: HistoryBuffers::new(),
            },
            |acc: &mut Acc<T>, rng, i| {
                if acc.err.is_some() {
                    return;
                }
                if let Err(e) = self.history_in(&mut acc.buffers, &mut acc.tally, rng, i) {
                    acc.err = Some((i, e));
                }
            },
            |a, b| {
                if a.err.is_none() {
                    a.err = b.err;
                }
                a.tally.merge(b.tally);
            },
        );
        match acc.err {
            Some((index, HistoryError::Stopping(source))) => {
                Err(BcaError::Stopping { index, source })
            }
            Some((index, HistoryError::BeamDivergence { attempts })) => {
                Err(BcaError::BeamDivergence { index, attempts })
            }
            None => Ok(acc.tally),
        }
    }

    /// Simulate primary `index` and its cascade with the given stream (use
    /// [`crate::rng::stream`]`(seed, index)` for results consistent with
    /// [`Bca::run`]). Returns the history's energy budget.
    pub fn history<T: BcaTally>(
        &self,
        tally: &mut T,
        rng: &mut ParticleRng,
        index: u64,
    ) -> Result<EnergyBudget, HistoryError> {
        self.history_in(&mut HistoryBuffers::new(), tally, rng, index)
    }

    /// [`Bca::history`] in caller-owned working memory. Identical results;
    /// reusing `buffers` across histories removes the per-history
    /// allocations ([`Bca::run`] does this per worker chunk).
    pub fn history_in<T: BcaTally>(
        &self,
        buffers: &mut HistoryBuffers,
        tally: &mut T,
        rng: &mut ParticleRng,
        index: u64,
    ) -> Result<EnergyBudget, HistoryError> {
        tally.begin_history(index);
        let mut budget = EnergyBudget {
            incident: self.beam.energy_ev,
            ..EnergyBudget::default()
        };
        let HistoryBuffers {
            pending,
            scratch,
            shifts,
            thermal_rng,
        } = buffers;
        // Each crystal's lattice is translated by a random vector for every
        // history, uniformly over the rectangular cell that tiles it
        // (`Lattice::orthogonal_cell`), so the beam samples all positions of
        // the unit cell. The shifts come from a copy of this history's stream
        // positioned at word 2^67, half the ChaCha period: a segment the
        // transport draws below can never reach. So the transport stream is the same with or
        // without crystals, and a particle that never enters a crystal
        // region sees exactly the draws of the amorphous engine
        // (`tests/crystal_amorphous_identity.rs`).
        shifts.clear();
        if !self.crystals.is_empty() {
            let mut srng = rng.clone();
            srng.set_word_pos(1u128 << 67);
            for cr in &self.crystals {
                let e = cr.cell_edges();
                shifts.push([
                    e[0] * Self::uniform(&mut srng),
                    e[1] * Self::uniform(&mut srng),
                    e[2] * Self::uniform(&mut srng),
                ]);
            }
        }
        // Thermal displacements of lattice sites come from a third copy of
        // the stream, at word 2^66: like the shifts, a segment neither the
        // transport draws nor the shifts reach, so turning vibration on
        // changes no other draw of the history.
        *thermal_rng = if self.crystals.iter().any(|c| c.is_thermal()) {
            let mut trng = rng.clone();
            trng.set_word_pos(1u128 << 66);
            Some(trng)
        } else {
            None
        };
        pending.clear();
        let dir = self.initial_direction(rng)?;
        pending.push(Particle {
            species: 0,
            z: self.beam.ion.z(),
            mass_amu: self.beam.ion.mass_amu(),
            energy_ev: self.beam.energy_ev,
            pos: self.entry_pos,
            dir,
            layer: self.entry_region,
            origin_layer: self.entry_region,
            generation: 0,
        });
        while let Some(p) = pending.pop() {
            self.transport(
                p,
                rng,
                pending,
                &mut budget,
                tally,
                scratch,
                shifts,
                thermal_rng,
            )?;
        }
        tally.end_history(index, &budget);
        Ok(budget)
    }

    /// Per-material data of the material filling `region`.
    fn lay(&self, region: usize) -> &LayerData {
        &self.layers[self.geometry.material_index(region)]
    }

    fn uniform(rng: &mut ParticleRng) -> f64 {
        (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Number of mean free paths for a new flight.
    fn draw_tau(&self, rng: &mut ParticleRng, first_primary_flight: bool) -> f64 {
        match self.config.mean_free_path {
            MeanFreePath::Constant if first_primary_flight => Self::uniform(rng),
            MeanFreePath::Constant => 1.0,
            MeanFreePath::EnergyDependent { .. } => -(1.0 - Self::uniform(rng)).ln(),
        }
    }

    /// Mean free path at the particle's energy in its layer; fills
    /// `scratch` with the per-element `p_max` and cumulative partner weights.
    fn flight(&self, p: &Particle, lay: &LayerData, scratch: &mut Scratch) -> f64 {
        scratch.p_max.clear();
        scratch.cum.clear();
        let mut acc = 0.0;
        match &self.beta_max {
            None => {
                for e in &lay.elems {
                    acc += e.fraction;
                    scratch.p_max.push(lay.p_const);
                    scratch.cum.push(acc);
                }
                lay.ell
            }
            Some(grid) => {
                let ns = self.species.len();
                for e in &lay.elems {
                    let pair = &self.pairs[p.species * ns + e.species];
                    let pm = (grid.beta(p.energy_ev * pair.eps_per_ev) * pair.a).min(lay.p_const);
                    acc += e.fraction * pm * pm;
                    scratch.p_max.push(pm);
                    scratch.cum.push(acc);
                }
                1.0 / (lay.n * PI * acc)
            }
        }
    }

    fn electronic_nonlocal<T: BcaTally>(
        &self,
        p: &mut Particle,
        from: [f64; 3],
        len: f64,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> Result<(), StoppingError> {
        if len <= 0.0 {
            return Ok(());
        }
        let mi = self.geometry.material_index(p.layer);
        let lay = &self.layers[mi];
        let s = match lay.loss_coef[p.species] {
            // `c sqrt(E)`; the `E > 0` test keeps a non-physical energy on the
            // general path, which reports it.
            Some(c) if p.energy_ev > 0.0 && p.energy_ev.is_finite() => c * p.energy_ev.sqrt(),
            _ => {
                let model: &dyn ElectronicStopping = match self.config.electronic {
                    ElectronicLoss::NonLocal => self.stopping,
                    ElectronicLoss::EquipartitionLsOr => &self.mix,
                };
                let material = self.geometry.material(mi);
                let ion = &self.species[p.species].ion;
                bragg_cross_section_per_atom(model, &NoCorrection, ion, material, p.energy_ev)?
            }
        };
        let de = (lay.n * s * len / J_PER_EV).min(p.energy_ev);
        if de > 0.0 {
            p.energy_ev -= de;
            budget.electronic_nonlocal += de;
            tally.electronic(p, from, ElectronicChannel::NonLocal, de);
        }
        Ok(())
    }

    fn advance(p: &mut Particle, len: f64) {
        for k in 0..3 {
            p.pos[k] += len * p.dir[k];
        }
    }

    /// Follow one particle until it stops or escapes; push its recoils.
    #[allow(clippy::too_many_arguments)]
    fn transport<T: BcaTally>(
        &self,
        mut p: Particle,
        rng: &mut ParticleRng,
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
        scratch: &mut Scratch,
        shifts: &[[f64; 3]],
        thermal_rng: &mut Option<ParticleRng>,
    ) -> Result<(), StoppingError> {
        let cutoff = self.species[p.species].cutoff_ev;
        let mut tau = self.draw_tau(rng, p.is_primary());
        scratch.last.clear();
        scratch.new_line();
        loop {
            if p.energy_ev < cutoff {
                budget.rest += p.energy_ev;
                tally.stopped(&p);
                return Ok(());
            }
            if let Some(ci) = self.crystal_at(p.layer) {
                if self.crystal_step(
                    ci,
                    shifts[ci],
                    &mut p,
                    pending,
                    budget,
                    tally,
                    scratch,
                    thermal_rng,
                )? {
                    return Ok(());
                }
                if self.crystal_at(p.layer).is_none() {
                    // Left the crystal for an amorphous region: a fresh flight.
                    tau = self.draw_tau(rng, false);
                    scratch.last.clear();
                    scratch.new_line();
                }
                continue;
            }
            let lay = self.lay(p.layer);
            let lambda = self.flight(&p, lay, scratch);
            let s = tau * lambda;

            let ex = match self.geometry.flight(p.layer, p.pos, p.dir, s) {
                Flight::Event(ex) => ex,
                Flight::Clear { region } => {
                    // The whole flight is in one material, but it may have
                    // crossed same-material faces (voxels of one block):
                    // take the region it ends in before colliding, so the
                    // next flight, the tallies and any recoil start from the
                    // right one.
                    let from = p.pos;
                    Self::advance(&mut p, s);
                    p.layer = region;
                    self.electronic_nonlocal(&mut p, from, s, budget, tally)?;
                    if p.energy_ev >= cutoff {
                        self.collide(&mut p, rng, scratch, pending, budget, tally)?;
                        tau = self.draw_tau(rng, false);
                    }
                    continue;
                }
            };
            // Truncate at the event; keep the unused part of the flight.
            let d_b = ex.distance;
            let from = p.pos;
            p.pos = ex.at;
            self.electronic_nonlocal(&mut p, from, d_b, budget, tally)?;
            tau = (tau - d_b / lambda).max(0.0);
            if p.energy_ev < cutoff {
                continue;
            }
            if self.apply_event(&mut p, ex.outcome, budget, tally) {
                return Ok(());
            }
        }
    }

    /// Apply the outcome of a geometry event to `p`, which is already at the
    /// event point with the electronic loss of the segment taken: enter the
    /// next region, or reach a surface and escape (returns `true`, the
    /// particle is finished) or reflect specularly (returns `false`).
    fn apply_event<T: BcaTally>(
        &self,
        p: &mut Particle,
        outcome: ExitOutcome,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> bool {
        match outcome {
            ExitOutcome::Enter { region, pos } => {
                p.layer = region;
                p.pos = pos;
            }
            ExitOutcome::Escape { face, normal } => {
                let mi = self.geometry.material_index(p.layer);
                let e_s = self.barriers[mi * self.species.len() + p.species];
                match kinematics::refract_out_normal(p.energy_ev, p.dir, normal, e_s) {
                    Some((e_out, dir_out)) => {
                        budget.surface_barrier += p.energy_ev - e_out;
                        p.energy_ev = e_out;
                        p.dir = dir_out;
                        match (face, p.is_primary()) {
                            (Face::Front, true) => budget.backscattered += e_out,
                            (Face::Front, false) => budget.sputtered += e_out,
                            (Face::Back, _) => budget.transmitted += e_out,
                            (Face::Side, _) => budget.lateral += e_out,
                        }
                        tally.escaped(p, face);
                        return true;
                    }
                    None => {
                        // Specular reflection about the face normal.
                        let dn = p.dir[0] * normal[0] + p.dir[1] * normal[1] + p.dir[2] * normal[2];
                        for (d, n) in p.dir.iter_mut().zip(normal) {
                            *d -= 2.0 * dn * n;
                        }
                    }
                }
            }
        }
        false
    }

    /// `tan(theta / 2)` for the centre-of-mass angle at `(eps, beta)`.
    fn half_angle_tan(&self, eps: f64, beta: f64) -> f64 {
        // Outside the tabulated energy range fall back to direct quadrature
        // (deterministic, just slower).
        self.table
            .half_angle_tan(eps, beta)
            .unwrap_or_else(|| (0.5 * theta_quadrature(self.screening, eps, beta)).tan())
    }

    /// Partner element index, by stoichiometry (weighted by `p_max^2` for
    /// the energy-dependent free path).
    fn draw_partner(rng: &mut ParticleRng, scratch: &Scratch) -> usize {
        let total = *scratch.cum.last().expect("layer has elements");
        let u = Self::uniform(rng) * total;
        scratch
            .cum
            .iter()
            .position(|&c| u < c)
            .unwrap_or(scratch.cum.len() - 1)
    }

    /// One binary collision of `p`, at its current energy and direction,
    /// with element `j` of its layer at impact parameter `b` and azimuth
    /// azimuth `(sin, cos)` = `azimuth`. Deflects `p`, removes the transfer `T`
    /// and, under [`ElectronicLoss::EquipartitionLsOr`], the local electronic
    /// loss (clamped to the energy left), and returns `(T, incoming
    /// direction, recoil lab angle as (sin, cos))`.
    ///
    /// The angles are carried as sines and cosines: with `t = tan(theta/2)`
    /// from the table, `sin(theta/2) = t / sqrt(1 + t^2)` and `cos(theta/2) =
    /// 1 / sqrt(1 + t^2)` give the energy transfer, the projectile deflection
    /// ([`kinematics::lab_projectile_sc`]) and the recoil deflection
    /// `phi = (pi - theta) / 2`, `(sin, cos) = (cos(theta/2), sin(theta/2))`
    /// with one square root and no `atan`, `atan2`, `sin` or `cos`.
    #[allow(clippy::too_many_arguments)]
    fn binary<T: BcaTally>(
        &self,
        p: &mut Particle,
        j: usize,
        b: f64,
        azimuth: (f64, f64),
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> Result<BinaryOutcome, StoppingError> {
        tally.partner(p, b, 1);
        let elem = &self.lay(p.layer).elems[j];
        let pair = self.pairs[p.species * self.species.len() + elem.species];
        let e0 = p.energy_ev;
        let eps = e0 * pair.eps_per_ev;
        let beta = b / pair.a;
        let tan_half = self.half_angle_tan(eps, beta);
        let c = 1.0 / (1.0 + tan_half * tan_half).sqrt();
        let s = tan_half * c;
        let t = (pair.gamma * e0 * s * s).min(e0);
        let psi = kinematics::lab_projectile_sc(s, c, pair.mu);
        let incoming = p.dir;
        p.dir = kinematics::rotate_sc(incoming, psi, azimuth);
        p.energy_ev -= t;

        if self.config.electronic == ElectronicLoss::EquipartitionLsOr {
            let r_min = closest_approach(self.screening, eps, beta) * pair.a;
            let ion = &self.species[p.species].ion;
            let q = (self.mix.0.local_loss(ion, elem.z, e0, r_min)? / J_PER_EV).min(p.energy_ev);
            if q > 0.0 {
                p.energy_ev -= q;
                budget.electronic_local += q;
                tally.electronic(p, p.pos, ElectronicChannel::Local, q);
            }
        }
        Ok((t, incoming, (c, s)))
    }

    /// Whether the partner of a weak collision of `p` at impact parameter
    /// `b` and azimuth `azimuth` lies in front of the target's front surface,
    /// in vacuum. The partner sits at distance `b` from the path on the side
    /// opposite to the deflection, where a hard collision sends its recoil:
    /// at `b rotate(dir, pi/2, azimuth + pi)` from the collision site.
    fn partner_beyond_surface(&self, p: &Particle, b: f64, azimuth: (f64, f64)) -> bool {
        // sin and cos of `azimuth + pi` are the negatives.
        let offset = kinematics::rotate_sc(p.dir, (1.0, 0.0), (-azimuth.0, -azimuth.1));
        self.geometry.in_vacuum([
            p.pos[0] + b * offset[0],
            p.pos[1] + b * offset[1],
            p.pos[2] + b * offset[2],
        ])
    }

    /// One collision step: the `weak_collisions` simultaneous weak
    /// collisions, then the hard collision (see the module docs).
    fn collide<T: BcaTally>(
        &self,
        p: &mut Particle,
        rng: &mut ParticleRng,
        scratch: &Scratch,
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> Result<(), StoppingError> {
        // Weak collisions (Moller and Eckstein, IPP 9/64 (1988), p. 14,
        // eq. (26)): the k-th partner sits in the k-th annulus,
        // p = p_max sqrt(k + R), with its own species and azimuth. They come
        // before the hard collision, the direction is updated after each, and
        // only the hard collision makes a recoil (p. 26); a weak transfer
        // stays in the lattice at the site. A weak partner that would lie
        // beyond the front surface is skipped, after its draws are made (the
        // "Target beyond Surf.?" test inside the weak-collision loop of the
        // flow charts, Figs. 7 and 8, pp. 29-30). Each collision is evaluated
        // at the energy left after the previous one, a deliberate deviation
        // from TRIDYN (module docs).
        for k in 1..=self.config.weak_collisions {
            let j = Self::draw_partner(rng, scratch);
            let b = scratch.p_max[j] * (f64::from(k) + Self::uniform(rng)).sqrt();
            let azimuth = (2.0 * PI * Self::uniform(rng)).sin_cos();
            if self.partner_beyond_surface(p, b, azimuth) {
                continue;
            }
            let (t, _, _) = self.binary(p, j, b, azimuth, budget, tally)?;
            if t > 0.0 {
                budget.lattice += t;
                tally.lattice(p.pos, p.layer, LatticeDeposit::Weak, t);
            }
        }

        // The hard collision: impact parameter uniform over the disc of
        // radius p_max, p = p_max sqrt(R). Its partner is not tested against
        // the surface, unlike in TRIDYN: that would change the engine
        // without weak collisions (module docs).
        let j = Self::draw_partner(rng, scratch);
        let b = scratch.p_max[j] * Self::uniform(rng).sqrt();
        let azimuth = (2.0 * PI * Self::uniform(rng)).sin_cos();
        let (t, incoming, phi) = self.binary(p, j, b, azimuth, budget, tally)?;

        self.emit_recoil(
            p.generation,
            p.pos,
            p.layer,
            j,
            t,
            || kinematics::rotate_sc(incoming, phi, (-azimuth.0, -azimuth.1)),
            pending,
            budget,
            tally,
        );
        Ok(())
    }

    /// Dispose of the energy `t` transferred to element `j` of the material
    /// of `layer` at `pos`: a displaced atom (if `t` exceeds both `E_d` and
    /// `E_b`; it starts with `t - E_b` along `recoil_dir()` and is followed
    /// or stopped as configured), otherwise heat in the lattice. The same
    /// criterion for the amorphous and the crystal flight model.
    #[allow(clippy::too_many_arguments)]
    fn emit_recoil<T: BcaTally>(
        &self,
        parent_generation: u32,
        pos: [f64; 3],
        layer: usize,
        j: usize,
        t: f64,
        recoil_dir: impl FnOnce() -> [f64; 3],
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) {
        let elem = &self.lay(layer).elems[j];
        if t > elem.e_d_ev && t > elem.e_b_ev {
            let e_r = t - elem.e_b_ev;
            if elem.e_b_ev > 0.0 {
                budget.lattice += elem.e_b_ev;
                tally.lattice(pos, layer, LatticeDeposit::Binding, elem.e_b_ev);
            }
            let sp = &self.species[elem.species];
            let r = Particle {
                species: elem.species,
                z: elem.z,
                mass_amu: sp.ion.mass_amu(),
                energy_ev: e_r,
                pos,
                dir: recoil_dir(),
                layer,
                origin_layer: layer,
                generation: parent_generation + 1,
            };
            tally.recoil(&r);
            if self.config.follow_recoils && e_r >= sp.cutoff_ev {
                pending.push(r);
            } else {
                budget.rest += e_r;
                tally.stopped(&r);
            }
        } else if t > 0.0 {
            budget.lattice += t;
            tally.lattice(pos, layer, LatticeDeposit::Subthreshold, t);
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn inward_sampling_is_bounded_and_cannot_hang() {
        use crate::ion::crystal::Divergence;
        // A zero-width spread about an outward direction never yields an
        // inward sample: the bound turns the loop into an error.
        let mut rng = crate::rng::stream(1, 0);
        let r = sample_inward(
            &Divergence::Gaussian { sigma_rad: 0.0 },
            [-1.0, 0.0, 0.0],
            &mut rng,
        );
        assert_eq!(
            r,
            Err(HistoryError::BeamDivergence {
                attempts: MAX_DIVERGENCE_ATTEMPTS
            })
        );
    }

    use super::*;
    use crate::geometry::Stack;
    use crate::ion::scattering::TableSpec;
    use crate::material::Material;

    fn si() -> Material {
        let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        m.set_displacement_energy_ev(14, 15.0).unwrap();
        m
    }

    fn small_table() -> ScatteringTable {
        let pot = Potential::new(Screening::ZblUniversal, 14.0, 14.0);
        ScatteringTable::build(
            &pot,
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e3,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 8,
            },
        )
    }

    #[test]
    fn beam_direction_is_unit_and_inward() {
        let b = Beam {
            ion: Ion::new(5).unwrap(),
            energy_ev: 1e3,
            polar_rad: 0.7,
            azimuth_rad: 2.0,
            count: 1,
        };
        let d = b.direction();
        assert!((d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - 1.0).abs() < 1e-15);
        assert!((d[0] - 0.7f64.cos()).abs() < 1e-15);
    }

    #[test]
    fn constant_convention_sweeps_one_atom_per_flight() {
        let stack = Stack::semi_infinite(si());
        let table = small_table();
        let ls = crate::ion::stopping::lindhard_scharff::LindhardScharff::new();
        let bca = Bca::new(
            Beam::normal(Ion::new(5).unwrap(), 1e3, 1),
            &stack,
            BcaConfig::new(5.0, 2.0),
            &ls,
            &table,
        )
        .unwrap();
        let l = &bca.layers[0];
        assert!((l.n * PI * l.p_const * l.p_const * l.ell - 1.0).abs() < 1e-14);
    }

    #[test]
    fn beta_max_grid_inverts_the_table() {
        let table = small_table();
        let th = 0.02;
        let g = BetaMaxGrid::build(&table, th);
        for &eps in &[1e-3, 0.1, 10.0, 300.0] {
            let b = g.beta(eps);
            let got = table.theta(eps, b).unwrap();
            assert!((got / th - 1.0).abs() < 0.05, "eps={eps}: {got}");
        }
    }

    #[test]
    fn weak_partner_surface_test() {
        let stack = Stack::semi_infinite(si());
        let table = small_table();
        let ls = crate::ion::stopping::lindhard_scharff::LindhardScharff::new();
        let bca = Bca::new(
            Beam::normal(Ion::new(14).unwrap(), 1e3, 1),
            &stack,
            BcaConfig::new(5.0, 2.0),
            &ls,
            &table,
        )
        .unwrap();
        let b = 2e-10;
        let at = |x: f64, dir: [f64; 3]| Particle {
            species: 0,
            z: 14,
            mass_amu: 28.0,
            energy_ev: 10.0,
            pos: [x, 0.0, 0.0],
            dir,
            layer: 0,
            origin_layer: 0,
            generation: 0,
        };
        let n = 3600;
        let azimuths = (0..n).map(|i| 2.0 * PI * (i as f64 + 0.5) / n as f64);
        let beyond = |p: &Particle| {
            azimuths
                .clone()
                .filter(|&az| bca.partner_beyond_surface(p, b, az.sin_cos()))
                .count() as f64
                / n as f64
        };
        // Moving parallel to the surface: on the surface half of the
        // partners are in vacuum, at depth b/2 a third (offset_x < -1/2),
        // at depth b none. Moving along the normal, none.
        let par = [0.0, 1.0, 0.0];
        assert!((beyond(&at(0.0, par)) - 0.5).abs() < 1e-3);
        assert!((beyond(&at(0.5 * b, par)) - 1.0 / 3.0).abs() < 1e-3);
        assert_eq!(beyond(&at(1.0001 * b, par)), 0.0);
        assert_eq!(beyond(&at(0.0, [1.0, 0.0, 0.0])), 0.0);
        // The partner is on the side opposite to the deflection.
        let p = at(0.0, par);
        for az in azimuths {
            let deflected = kinematics::rotate(p.dir, 0.3, az);
            if deflected[0].abs() > 1e-6 {
                assert_eq!(
                    bca.partner_beyond_surface(&p, b, az.sin_cos()),
                    deflected[0] > 0.0
                );
            }
        }
    }

    #[test]
    fn rejects_bad_setup() {
        let table = small_table();
        let ls = crate::ion::stopping::lindhard_scharff::LindhardScharff::new();
        let si = Stack::semi_infinite(si());
        let ion = Ion::new(5).unwrap();
        let ok = BcaConfig::new(5.0, 2.0);
        let mk = |beam: Beam, cfg: BcaConfig| Bca::new(beam, &si, cfg, &ls, &table).map(|_| ());
        assert!(matches!(
            mk(Beam::normal(ion, 0.0, 1), ok),
            Err(BcaError::InvalidBeam(_))
        ));
        let mut b = Beam::normal(ion, 1e3, 1);
        b.polar_rad = 0.5 * PI;
        assert!(matches!(mk(b, ok), Err(BcaError::InvalidBeam(_))));
        assert!(matches!(
            mk(Beam::normal(ion, 1e3, 1), BcaConfig::new(0.0, 2.0)),
            Err(BcaError::InvalidConfig(_))
        ));
        let mut c = ok;
        c.mean_free_path = MeanFreePath::EnergyDependent {
            min_cm_angle_rad: 0.0,
        };
        assert!(matches!(
            mk(Beam::normal(ion, 1e3, 1), c),
            Err(BcaError::InvalidConfig(_))
        ));
        assert!(mk(Beam::normal(ion, 1e3, 1), ok).is_ok());
        let mut c = ok;
        c.weak_collisions = MAX_WEAK_COLLISIONS;
        assert!(mk(Beam::normal(ion, 1e3, 1), c).is_ok());
        c.weak_collisions = MAX_WEAK_COLLISIONS + 1;
        assert!(matches!(
            mk(Beam::normal(ion, 1e3, 1), c),
            Err(BcaError::InvalidConfig(_))
        ));
        c.weak_collisions = 1;
        c.mean_free_path = MeanFreePath::EnergyDependent {
            min_cm_angle_rad: 0.01,
        };
        assert!(matches!(
            mk(Beam::normal(ion, 1e3, 1), c),
            Err(BcaError::InvalidConfig(_))
        ));
    }
}
