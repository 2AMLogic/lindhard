//! Amorphous binary-collision-approximation (BCA) transport in a 1D layered
//! target, with full recoil cascades.
//!
//! # Model
//!
//! The target is a [`Stack`] of homogeneous, structureless (amorphous)
//! layers. A particle alternates straight free flights and binary elastic
//! collisions with target atoms, losing energy continuously to electrons along
//! each flight. This is the amorphous-target BCA of J. P. Biersack and L. G.
//! Haggmark, Nucl. Instrum. Methods 174 (1980) 257 (TRIM), within the general
//! BCA framework of M. T. Robinson and I. M. Torrens, Phys. Rev. B 9 (1974)
//! 5008 and the treatment in W. Eckstein, *Computer Simulation of Ion-Solid
//! Interactions* (Springer, 1991). Everything here is implemented from those
//! publications; no code from the Tier B/C programs named in
//! `CONTRIBUTING.md` was consulted.
//!
//! One history:
//!
//! 1. The primary enters at the origin of the front face with the beam's
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
//! A flight is drawn as a dimensionless number of mean free paths `tau` and
//! converted to a length with the local `lambda`. When a flight reaches an
//! interface it is truncated exactly there, the electronic loss for the
//! truncated length is applied, the material is switched, and the flight
//! continues in the new layer with the **unused** part `tau - s/lambda` (the
//! flight is not redrawn). The collision then happens in the layer where the
//! flight ends, with a partner and `p_max` from that layer. For exponential
//! paths this is exact by memorylessness; for the constant path it means a
//! layer split of one material changes nothing but floating-point rounding,
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
//! <https://pure.mpg.de/rest/items/item_2131703/component/file_2131702/content>:
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
//!
//! So each collision step here is: the `K` weak collisions, `k = 1..K`, each
//! with a partner drawn by stoichiometry, `p = p_max sqrt(k + R)` and a
//! uniform azimuth, then the hard collision with `p = p_max sqrt(R)`. Each
//! deflects the particle and takes its transfer `T`; under
//! [`ElectronicLoss::EquipartitionLsOr`] each also takes its Oen-Robinson
//! local loss, so the local half is sampled out to `p_max sqrt(K + 1)`. A weak
//! collision never makes a recoil: its `T` stays in the lattice at the site
//! ([`LatticeDeposit::Weak`], counted in [`EnergyBudget::lattice`]). The
//! report does not say whether the simultaneous collisions use the energy the
//! step began with or the energy left after the previous one; here each uses
//! the energy left (and the direction after the previous one), so no transfer
//! can exceed the particle's energy. Random draws per weak collision: partner,
//! `R`, azimuth, in that order and before the hard collision's, so `K = 0`
//! draws exactly what the engine drew before the option existed and
//! reproduces it bit for bit.
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
//! Ar 1 keV on Cu sputter yield with `E_d = E_s` falls from 1.93 to 0.66 with
//! `K = 3` (`docs/validation.md`, level 2). Use it for the energy partition
//! (damage energy, electronic share); its effect on yields is large and has
//! not been checked against measured yields.
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
//! # Not modelled here (see the issue tracker)
//!
//! Crystal structure, target composition changes with fluence, and refraction
//! of the incident beam at the entrance surface (negligible at keV energies).

pub mod kinematics;
pub mod tally;

pub use tally::{BcaTally, ElectronicChannel, EnergyBudget, Face, LatticeDeposit, SummaryTally};

use std::f64::consts::PI;

use rand_core::Rng;

use crate::geometry::Stack;
use crate::ion::potential::{Potential, Screening};
use crate::ion::scattering::{closest_approach, theta_quadrature, ScatteringTable};
use crate::ion::stopping::bragg::{bragg_cross_section_per_atom, NoCorrection};
use crate::ion::stopping::mix::EquipartitionMix;
use crate::ion::stopping::{ElectronicStopping, Ion, StoppingError, ValidityRange};
use crate::material::EnergyKind;
use crate::rng::{run_particles, ParticleRng};
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
    /// Histories per work chunk for [`run_particles`]; fixes the summation
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
    /// A layer material has an energy with no value; set it before running.
    #[error("layer {layer}: {kind} for element Z={z} is not set")]
    EnergyNotSet {
        /// Layer index.
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
    /// Index of the layer the particle is in.
    pub layer: usize,
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
    barrier_front_ev: f64,
    barrier_back_ev: f64,
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

/// Per-history scratch space for the energy-dependent free path.
#[derive(Default)]
struct Scratch {
    p_max: Vec<f64>,
    cum: Vec<f64>,
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
}

/// The BCA engine for one beam, target and configuration.
///
/// Holds borrowed, read-only inputs and precomputed per-layer and per-pair
/// data, so it is `Sync` and shared by all worker threads.
pub struct Bca<'a> {
    beam: Beam,
    stack: &'a Stack,
    config: BcaConfig,
    stopping: &'a (dyn ElectronicStopping + Sync),
    mix: MixNonLocal,
    table: &'a ScatteringTable,
    screening: Screening,
    species: Vec<SpeciesData>,
    layers: Vec<LayerData>,
    /// `pairs[projectile species * n_species + target species]`.
    pairs: Vec<PairData>,
    beta_max: Option<BetaMaxGrid>,
}

fn finite_nonneg(v: f64) -> bool {
    v.is_finite() && v >= 0.0
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
    ///   material at that face; an element absent from that material uses its
    ///   `E_s` in the first layer that contains it.
    pub fn new(
        beam: Beam,
        stack: &'a Stack,
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
        for (layer, l) in stack.layers().iter().enumerate() {
            if let Some(&(z, kind)) = l.material().unset_energies().first() {
                return Err(BcaError::EnergyNotSet { layer, z, kind });
            }
        }

        // Species: the beam, then target elements in order of appearance.
        let mut species = vec![SpeciesData {
            ion: beam.ion,
            cutoff_ev: config.primary_cutoff_ev,
            barrier_front_ev: config.primary_surface_binding_ev,
            barrier_back_ev: config.primary_surface_binding_ev,
        }];
        let mut species_z: Vec<u8> = vec![0];
        let mut layers = Vec::with_capacity(stack.layers().len());
        for l in stack.layers() {
            let m = l.material();
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
                            barrier_front_ev: e_s,
                            barrier_back_ev: e_s,
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
            });
        }
        // Face barriers from the materials at the faces.
        let front = stack.layers().first().expect("non-empty").material();
        let back = stack.layers().last().expect("non-empty").material();
        for (s, &z) in species.iter_mut().zip(&species_z).skip(1) {
            if let Ok(e) = front.surface_binding_energy_ev(z) {
                s.barrier_front_ev = e;
            }
            if let Ok(e) = back.surface_binding_energy_ev(z) {
                s.barrier_back_ev = e;
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
            stack,
            config,
            stopping,
            mix: MixNonLocal(EquipartitionMix::new()),
            table,
            screening: tp.screening,
            species,
            layers,
            pairs,
            beta_max,
        })
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
        struct Acc<T> {
            tally: T,
            err: Option<(u64, StoppingError)>,
        }
        let acc = run_particles(
            self.config.seed,
            self.beam.count,
            self.config.chunk_size,
            || Acc {
                tally: new_tally(),
                err: None,
            },
            |acc: &mut Acc<T>, rng, i| {
                if acc.err.is_some() {
                    return;
                }
                if let Err(e) = self.history(&mut acc.tally, rng, i) {
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
            Some((index, source)) => Err(BcaError::Stopping { index, source }),
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
    ) -> Result<EnergyBudget, StoppingError> {
        tally.begin_history(index);
        let mut budget = EnergyBudget {
            incident: self.beam.energy_ev,
            ..EnergyBudget::default()
        };
        let mut pending = vec![Particle {
            species: 0,
            z: self.beam.ion.z(),
            mass_amu: self.beam.ion.mass_amu(),
            energy_ev: self.beam.energy_ev,
            pos: [0.0; 3],
            dir: self.beam.direction(),
            layer: 0,
            generation: 0,
        }];
        let mut scratch = Scratch::default();
        while let Some(p) = pending.pop() {
            self.transport(p, rng, &mut pending, &mut budget, tally, &mut scratch)?;
        }
        tally.end_history(index, &budget);
        Ok(budget)
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
        let model: &dyn ElectronicStopping = match self.config.electronic {
            ElectronicLoss::NonLocal => self.stopping,
            ElectronicLoss::EquipartitionLsOr => &self.mix,
        };
        let material = self.stack.layers()[p.layer].material();
        let ion = &self.species[p.species].ion;
        let s = bragg_cross_section_per_atom(model, &NoCorrection, ion, material, p.energy_ev)?;
        let de = (self.layers[p.layer].n * s * len / J_PER_EV).min(p.energy_ev);
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
    fn transport<T: BcaTally>(
        &self,
        mut p: Particle,
        rng: &mut ParticleRng,
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
        scratch: &mut Scratch,
    ) -> Result<(), StoppingError> {
        let cutoff = self.species[p.species].cutoff_ev;
        let mut tau = self.draw_tau(rng, p.is_primary());
        loop {
            if p.energy_ev < cutoff {
                budget.rest += p.energy_ev;
                tally.stopped(&p);
                return Ok(());
            }
            let lay = &self.layers[p.layer];
            let lambda = self.flight(&p, lay, scratch);
            let s = tau * lambda;
            let geo = &self.stack.layers()[p.layer];
            let (d_b, boundary) = if p.dir[0] > 0.0 {
                ((geo.back_m() - p.pos[0]) / p.dir[0], geo.back_m())
            } else if p.dir[0] < 0.0 {
                ((geo.front_m() - p.pos[0]) / p.dir[0], geo.front_m())
            } else {
                (f64::INFINITY, f64::NAN)
            };

            if d_b <= s {
                // Truncate at the interface; keep the unused part of the flight.
                let from = p.pos;
                Self::advance(&mut p, d_b);
                p.pos[0] = boundary;
                self.electronic_nonlocal(&mut p, from, d_b, budget, tally)?;
                tau = (tau - d_b / lambda).max(0.0);
                if p.energy_ev < cutoff {
                    continue;
                }
                let outward = p.dir[0] > 0.0;
                let n_layers = self.layers.len();
                let face = match (outward, p.layer) {
                    (false, 0) => Some(Face::Front),
                    (true, l) if l + 1 == n_layers => Some(Face::Back),
                    _ => None,
                };
                match face {
                    None if outward => p.layer += 1,
                    None => p.layer -= 1,
                    Some(face) => {
                        let sp = &self.species[p.species];
                        let e_s = match face {
                            Face::Front => sp.barrier_front_ev,
                            Face::Back => sp.barrier_back_ev,
                        };
                        match kinematics::refract_out(p.energy_ev, p.dir, e_s) {
                            Some((e_out, dir_out)) => {
                                budget.surface_barrier += p.energy_ev - e_out;
                                p.energy_ev = e_out;
                                p.dir = dir_out;
                                match (face, p.is_primary()) {
                                    (Face::Front, true) => budget.backscattered += e_out,
                                    (Face::Front, false) => budget.sputtered += e_out,
                                    (Face::Back, _) => budget.transmitted += e_out,
                                }
                                tally.escaped(&p, face);
                                return Ok(());
                            }
                            None => p.dir[0] = -p.dir[0],
                        }
                    }
                }
                continue;
            }

            let from = p.pos;
            Self::advance(&mut p, s);
            self.electronic_nonlocal(&mut p, from, s, budget, tally)?;
            if p.energy_ev >= cutoff {
                self.collide(&mut p, rng, scratch, pending, budget, tally)?;
                tau = self.draw_tau(rng, false);
            }
        }
    }

    fn theta(&self, eps: f64, beta: f64) -> f64 {
        // Outside the tabulated energy range fall back to direct quadrature
        // (deterministic, just slower).
        self.table
            .theta(eps, beta)
            .unwrap_or_else(|| theta_quadrature(self.screening, eps, beta))
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
    /// `azimuth`. Deflects `p`, removes the transfer `T` and, under
    /// [`ElectronicLoss::EquipartitionLsOr`], the local electronic loss (clamped
    /// to the energy left), and returns `(T, incoming direction, recoil lab
    /// angle)`.
    #[allow(clippy::too_many_arguments)]
    fn binary<T: BcaTally>(
        &self,
        p: &mut Particle,
        j: usize,
        b: f64,
        azimuth: f64,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> Result<(f64, [f64; 3], f64), StoppingError> {
        let elem = &self.layers[p.layer].elems[j];
        let pair = self.pairs[p.species * self.species.len() + elem.species];
        let e0 = p.energy_ev;
        let eps = e0 * pair.eps_per_ev;
        let beta = b / pair.a;
        let theta = self.theta(eps, beta);
        let half = (0.5 * theta).sin();
        let t = (pair.gamma * e0 * half * half).min(e0);
        let (psi, phi) = kinematics::lab_angles(theta, pair.mu);
        let incoming = p.dir;
        p.dir = kinematics::rotate(incoming, psi, azimuth);
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
        Ok((t, incoming, phi))
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
        // stays in the lattice at the site. Each collision is evaluated at
        // the energy left after the previous one (see the module docs).
        for k in 1..=self.config.weak_collisions {
            let j = Self::draw_partner(rng, scratch);
            let b = scratch.p_max[j] * (f64::from(k) + Self::uniform(rng)).sqrt();
            let azimuth = 2.0 * PI * Self::uniform(rng);
            let (t, _, _) = self.binary(p, j, b, azimuth, budget, tally)?;
            if t > 0.0 {
                budget.lattice += t;
                tally.lattice(p.pos, p.layer, LatticeDeposit::Weak, t);
            }
        }

        // The hard collision: impact parameter uniform over the disc of
        // radius p_max, p = p_max sqrt(R).
        let j = Self::draw_partner(rng, scratch);
        let b = scratch.p_max[j] * Self::uniform(rng).sqrt();
        let azimuth = 2.0 * PI * Self::uniform(rng);
        let (t, incoming, phi) = self.binary(p, j, b, azimuth, budget, tally)?;

        let elem = &self.layers[p.layer].elems[j];
        if t > elem.e_d_ev && t > elem.e_b_ev {
            let e_r = t - elem.e_b_ev;
            if elem.e_b_ev > 0.0 {
                budget.lattice += elem.e_b_ev;
                tally.lattice(p.pos, p.layer, LatticeDeposit::Binding, elem.e_b_ev);
            }
            let sp = &self.species[elem.species];
            let r = Particle {
                species: elem.species,
                z: elem.z,
                mass_amu: sp.ion.mass_amu(),
                energy_ev: e_r,
                pos: p.pos,
                dir: kinematics::rotate(incoming, phi, azimuth + PI),
                layer: p.layer,
                generation: p.generation + 1,
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
            tally.lattice(p.pos, p.layer, LatticeDeposit::Subthreshold, t);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
