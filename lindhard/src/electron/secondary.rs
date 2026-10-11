//! Secondary-electron generation at inelastic events.
//!
//! The model is Kieft and Bosch's, J. Phys. D: Appl. Phys. 41, 215310 (2008),
//! doi:10.1088/0022-3727/41/21/215310, as written out by T. Verduin,
//! *Quantum Noise Effects in e-Beam Lithography and Metrology*, PhD thesis,
//! Delft University of Technology (2017),
//! doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c ("Verduin" below,
//! with the thesis's printed page and equation numbers). The kinematics are
//! ported from Nebula, `source/physics/kieft/inelastic.h` at commit
//! `a50a8e83207980ed221ac83cfe874751d39b86a9`
//! (<https://github.com/Nebula-simulator/nebula>), BSD-3-Clause,
//! Copyright (c) 2020, Nebula-simulator; the notice is in
//! `THIRD_PARTY_LICENSES.md`.
//!
//! # Energy
//!
//! An inelastic event takes `W` from the primary (from the inelastic
//! [`CrossSectionTable`](crate::electron::data::CrossSectionTable)) and
//! lifts an electron of binding energy `B` (measured below the Fermi level)
//! into the conduction band with kinetic energy
//! `E_SE = E_F + W - B` (Verduin Eq. 3.86, p. 78), measured, like every
//! kinetic energy inside a material, from the band bottom
//! ([`crate::electron::boundary`]). For a valence event (no shell channel
//! supplies `B`):
//!
//! - in a metal, `B = 0`: the electron comes from the Fermi level (Nebula's
//!   handling when no binding energy applies), and every event with `W > 0`
//!   makes one. A band built with
//!   [`BandStructure::with_valence_binding_ev`] sets `B` to that value
//!   instead (the rule of Azzolini et al., arXiv:1809.00859, p. 6; cited
//!   there), and then only an event with `W > B` makes one, a loss
//!   `W <= B` staying in the solid as below; `B = 0` is the default
//!   unchanged;
//! - in a semiconductor or insulator, `B = E_g` when `W > E_g`, an excitation
//!   across the gap (Verduin p. 78; Nebula sets the binding energy to the
//!   band gap); an event with `W <= E_g` makes no secondary and its energy
//!   stays in the solid (Verduin p. 78 and section 3.4 assign it to optical
//!   phonons).
//!
//! The bookkeeping of one event is a [`SecondaryEvent`], which satisfies
//! `loss = secondary + binding + deposited` with `binding = B - E_F`, the
//! energy that the liberated electron's initial state lies below the band
//! bottom (negative for a conduction electron, which already had `E_F`).
//!
//! # Inner-shell events
//!
//! A loss `ω` drawn from the channel of inner shell `j`
//! ([`crate::electron::transport::Transport::with_inner_shells`]) uses the
//! same Eq. 3.86 with `B = B_j`, the shell's binding energy:
//! `E_SE = E_F + ω - B_j`. This is the shell-resolved form of de Vera et al.,
//! Int. J. Mol. Sci. 23, 6121 (2022), eqs. (32) and (35) (the reference of
//! [`crate::electron::inelastic::inner_shell`]): a loss `ω = W + B_j` in shell
//! `j` gives the emitted electron `W = ω - B_j`, which here is its kinetic
//! energy above the Fermi level, so `E_SE = E_F + W` on the band-bottom axis.
//!
//! **Binding reference.** Eq. 3.86 measures `B` below the Fermi level, and so
//! does this module for a shell: `B_j` is used as given, as a binding energy
//! below `E_F`. The subshell energies of a
//! [`crate::electron::data::SubshellBindingTable`] (EADL) are free-atom
//! values; they are not converted to a solid-state, Fermi-referenced edge
//! (no work-function or chemical shift is applied). The event's
//! [`SecondaryEvent::binding_ev`] is `B_j - E_F`, positive for every shell
//! that lies below the Fermi level, and the tallies keep it in the solid at
//! the event (`tally::electron`): atomic relaxation (Auger electrons,
//! fluorescence photons) is not transported.
//!
//! **Limits.** The transport draws `ω` within `[B_j, E - E_F]` (`E` the
//! primary's band-bottom energy), so `E_SE` lies in `[E_F, E - B_j]` and is
//! never negative; a shell channel is closed when `E - E_F <= B_j`. Every
//! shell event liberates an electron (no gap test); the electron is followed
//! only if `E_SE` reaches the stopping threshold, else its energy is
//! deposited, as for a valence event.
//!
//! **Direction.** The Ivanchenko transformation below with `B = B_j`
//! (Verduin Eqs. 3.105-3.111 are written for a bound electron of any
//! binding).
//!
//! # Direction
//!
//! The method of Ivanchenko, as used by Kieft and Bosch (Verduin pp. 82-85):
//! the binding is transformed away, `E -> E - E_F + 2B` (Eq. 3.105) and
//! `ΔE -> W + B` (Eq. 3.106), and the free-electron binary-collision angle of
//! the secondary, `cos β = sqrt(ΔE / E)` (Eq. 3.100, non-relativistic), gives
//! its direction about the primary's with a uniform azimuth (Eq. 3.107).
//! Optionally a randomly oriented instantaneous momentum of relative size
//! `sqrt(B / ΔE)` is added (Eq. 3.108), and the primary is deflected to
//! conserve momentum, `p̂'_i ∝ p̂_i - cos β p̂''_t` (Eq. 3.111). Both options
//! are on in Nebula's defaults and in [`SecondaryModel::KIEFT_BOSCH`]. The
//! secondary starts at the primary's position (no delocalisation, Verduin
//! p. 85).

use rand_core::Rng;
use serde::Serialize;

use crate::electron::boundary::BandStructure;
use crate::electron::transport::{deflect_cs, normalize, uniform};

/// Whether and how inelastic events make secondary electrons.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(tag = "model", rename_all = "kebab-case")]
pub enum SecondaryModel {
    /// No secondaries; the inelastic event only takes `W` from the primary
    /// and keeps its direction.
    Off,
    /// Kieft and Bosch (module docs).
    KieftBosch {
        /// Add the instantaneous momentum of the bound electron to the
        /// secondary's direction (Verduin Eq. 3.108).
        instantaneous_momentum: bool,
        /// Deflect the primary to conserve momentum (Verduin Eq. 3.111).
        momentum_conservation: bool,
    },
}

impl SecondaryModel {
    /// Kieft and Bosch with both options on, Nebula's defaults
    /// (`kieft_inelastic` template arguments at the commit in the module
    /// docs).
    pub const KIEFT_BOSCH: Self = Self::KieftBosch {
        instantaneous_momentum: true,
        momentum_conservation: true,
    };
}

/// The energy bookkeeping of one inelastic event under a secondary model, all
/// in eV. `loss_ev = secondary_ev + binding_ev + deposited_ev` up to
/// rounding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondaryEvent {
    /// `W`, the energy the primary lost.
    pub loss_ev: f64,
    /// Whether an electron was lifted into the conduction band (`false` for a
    /// sub-gap loss in an insulator or a zero loss).
    pub liberated: bool,
    /// Kinetic energy of the secondary that was created and will be
    /// transported, `E_F + W - B`; zero if none was (nothing liberated, or
    /// the liberated electron was below the stopping threshold, in which case
    /// its energy is in `deposited_ev`).
    pub secondary_ev: f64,
    /// `B - E_F` if an electron was liberated, else zero: how far its initial
    /// state lay below the band bottom.
    pub binding_ev: f64,
    /// Energy left in the solid at the event: `W` for a sub-gap loss, the
    /// liberated electron's energy if it was not followed, else zero.
    pub deposited_ev: f64,
}

/// What [`kieft_bosch`] decided for one event.
pub(crate) struct Outcome {
    pub event: SecondaryEvent,
    /// The new primary direction.
    pub primary_dir: [f64; 3],
    /// The secondary's direction and energy, if one is to be followed.
    pub secondary: Option<([f64; 3], f64)>,
}

/// Apply the Kieft-Bosch secondary model to one valence inelastic event,
/// ported from Nebula's `kieft_inelastic::execute` (file and commit in the
/// module docs).
///
/// `dir` is the primary's unit direction and `e_ev` its energy before the
/// loss `w_ev`; `threshold_ev` is the stopping threshold of the layer (a
/// secondary below it is not followed). The binding energy comes from the
/// band (module docs, "Energy"). Draws, in order, only when an electron is
/// liberated: the azimuth of the secondary, then (with
/// `instantaneous_momentum`) `U1` and `U2` of Verduin Eq. 3.108.
#[allow(clippy::too_many_arguments)]
pub(crate) fn kieft_bosch<R: Rng>(
    instantaneous_momentum: bool,
    momentum_conservation: bool,
    band: &BandStructure,
    dir: [f64; 3],
    e_ev: f64,
    w_ev: f64,
    threshold_ev: f64,
    rng: &mut R,
) -> Outcome {
    let binding = match band.band_gap_ev() {
        None => {
            // The Fermi level (0) unless the band sets a valence binding.
            let b = band.valence_binding_ev().unwrap_or(0.0);
            (w_ev > b).then_some(b)
        }
        Some(gap) if w_ev > gap => Some(gap),
        Some(_) => None,
    };
    let Some(b) = binding else {
        return Outcome {
            event: SecondaryEvent {
                loss_ev: w_ev,
                liberated: false,
                secondary_ev: 0.0,
                binding_ev: 0.0,
                deposited_ev: w_ev,
            },
            primary_dir: dir,
            secondary: None,
        };
    };
    liberate(
        instantaneous_momentum,
        momentum_conservation,
        band.fermi_ev(),
        b,
        dir,
        e_ev,
        w_ev,
        threshold_ev,
        rng,
    )
}

/// Apply the Kieft-Bosch secondary model to one **inner-shell** ionisation:
/// a loss `w_ev` in the channel of a shell with binding energy
/// `shell_binding_ev` (module docs, "Inner-shell events").
///
/// The caller guarantees `shell_binding_ev <= w_ev <= e_ev - E_F` (the
/// channel's kinematic limits on the band-bottom axis), so the secondary
/// energy `E_F + w - B` is at least `E_F` and at most `e_ev - B`. The event
/// always liberates an electron; it is followed only if its energy reaches
/// `threshold_ev`. Draws as [`kieft_bosch`] does for a liberated electron.
#[allow(clippy::too_many_arguments)]
pub(crate) fn kieft_bosch_shell<R: Rng>(
    instantaneous_momentum: bool,
    momentum_conservation: bool,
    band: &BandStructure,
    shell_binding_ev: f64,
    dir: [f64; 3],
    e_ev: f64,
    w_ev: f64,
    threshold_ev: f64,
    rng: &mut R,
) -> Outcome {
    debug_assert!(w_ev >= shell_binding_ev && shell_binding_ev > 0.0);
    liberate(
        instantaneous_momentum,
        momentum_conservation,
        band.fermi_ev(),
        shell_binding_ev,
        dir,
        e_ev,
        w_ev,
        threshold_ev,
        rng,
    )
}

/// The event of a loss `w_ev` that lifts an electron of binding energy `b`
/// (below the Fermi level `fermi`) into the conduction band: Verduin
/// Eqs. 3.86 and 3.100-3.111 as ported from Nebula (module docs).
#[allow(clippy::too_many_arguments)]
fn liberate<R: Rng>(
    instantaneous_momentum: bool,
    momentum_conservation: bool,
    fermi: f64,
    b: f64,
    dir: [f64; 3],
    e_ev: f64,
    w_ev: f64,
    threshold_ev: f64,
    rng: &mut R,
) -> Outcome {
    // Verduin Eqs. 3.105, 3.106 and 3.100. `dk > 0` because `w > 0` and
    // `b >= 0`; the ratio is clamped to [0, 1] (a loss near the full energy,
    // or an energy below the Fermi level, would otherwise leave it).
    let k = e_ev - fermi + 2.0 * b;
    let dk = b + w_ev;
    let cos2 = if k > 0.0 {
        (dk / k).clamp(0.0, 1.0)
    } else {
        1.0
    };
    let cos_b = cos2.sqrt();
    let sin_b = (1.0 - cos2).sqrt();

    // Eq. 3.107: about the primary, uniform azimuth.
    let phi = std::f64::consts::TAU * uniform(rng);
    let mut sdir = deflect_cs(dir, cos_b, sin_b, phi);
    if instantaneous_momentum {
        // Eq. 3.108: cos(theta) = 2 U1 - 1, azimuth 2 pi U2.
        let ct = 2.0 * uniform(rng) - 1.0;
        let st = (1.0 - ct * ct).max(0.0).sqrt();
        let (sp, cp) = (std::f64::consts::TAU * uniform(rng)).sin_cos();
        let a = (b / dk).sqrt();
        sdir = normalize([
            sdir[0] + a * st * cp,
            sdir[1] + a * st * sp,
            sdir[2] + a * ct,
        ]);
    }

    let primary_dir = if momentum_conservation {
        // Eq. 3.111, up to normalisation.
        let p = [
            dir[0] - cos_b * sdir[0],
            dir[1] - cos_b * sdir[1],
            dir[2] - cos_b * sdir[2],
        ];
        let n2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        if n2 > 1e-24 && n2.is_finite() {
            normalize(p)
        } else {
            // A head-on transfer of all the momentum leaves no direction;
            // keep the old one.
            dir
        }
    } else {
        dir
    };

    // Eq. 3.86.
    let e_se = fermi + w_ev - b;
    let binding_ev = b - fermi;
    let followed = e_se >= threshold_ev;
    Outcome {
        event: SecondaryEvent {
            loss_ev: w_ev,
            liberated: true,
            secondary_ev: if followed { e_se } else { 0.0 },
            binding_ev,
            deposited_ev: if followed { 0.0 } else { e_se },
        },
        primary_dir,
        secondary: followed.then_some((sdir, e_se)),
    }
}
