//! Inner-shell ionization channels from shell-resolved loss functions.
//!
//! # Source
//!
//! P. de Vera, S. Taioli, P. E. Trevisanutto, M. Dapor, I. Abril, S.
//! Simonucci and R. Garcia-Molina, "Energy Deposition around Swift Carbon-Ion
//! Tracks in Liquid Water", Int. J. Mol. Sci. 23, 6121 (2022),
//! doi:10.3390/ijms23116121 (open access, PMC9181504; "dV2022", as in
//! [`super::penn`]):
//!
//! * eq. (1) splits the ELF into an outer-shell and an inner-shell part, and
//!   eq. (2) writes the inner part as a sum over elements `j` (weighted by the
//!   atoms per molecule `ν_j`) and subshells `nℓ` of atomic generalised
//!   oscillator strengths, each multiplied by the step `Θ(E - B_nℓ^j)`: an
//!   inner shell contributes only to losses above its binding energy;
//! * eqs. (32) and (35) write the secondary-electron spectrum and the
//!   ionization cross section as a sum of one outer-shell (valence) term and
//!   one term per inner shell `j`, each built from **that shell's own loss
//!   function** `Im[-1/ε(k, W + B_j)]_j`, where `W` is the kinetic energy of
//!   the emitted electron. So a loss `ω = W + B_j` in shell `j` hands the
//!   secondary `W = ω - B_j` ([`ShellResolvedChannels::secondary_energy_ev`]).
//!
//! # The model here
//!
//! [`ShellResolvedChannels`] takes what those equations take: a valence
//! optical ELF and **one optical ELF per inner shell**. The shell ELFs are
//! either supplied by the caller or built from atomic oscillator strengths
//! as in dV2022 eq. (2) by [`super::shell_elf`] (the hydrogenic K-shell
//! formula, with binding energies and occupancies from a
//! [`SubshellBindingTable`]; other subshells are not built there). A single
//! total ELF is never split: no partition rule of ours is applied. Each shell's ELF must
//! start at or above its binding energy (the step of eq. (2)). Each channel's
//! DIIMFP is the single-pole Penn DIIMFP ([`SinglePolePenn`]) of its own ELF;
//! with the exchange correction applying, shell `j` uses its binding energy
//! in the Born-Ochkur denominator `T' - W = T' - ω + B_j` and is limited to
//! `ω <= (T' + B_j)/2`, and the valence channel has `B = 0` (dV2022 eqs. (32),
//! (34), (35); see the `super::penn` module docs for the terms of eq. (32)
//! not followed: the Coulomb-field correction of the direct term, and the
//! energy-dependent valence binding energy `B(T)`).
//!
//! **Our choices, not from dV2022:** the momentum dependence of every
//! channel's loss function is the Penn SPA of S2017 applied to that
//! channel's optical ELF (dV2022 uses Mermin functions for the outer shells
//! and hydrogenic GOS for the inner ones); and the channel models share the
//! valence model's Fermi energy, tolerance and exchange setting. The SPA is
//! linear in the optical ELF, so without exchange the channel inverse mean
//! free paths add up to that of one model built from the sum of the ELFs
//! (when the sum is representable on one grid; tested).
//!
//! The binding energies come from a [`SubshellBindingTable`]
//! ([`InnerShell::from_table`]) or are given directly.

use super::penn::SinglePolePenn;
use crate::electron::data::{ElectronDataError, OpticalElf, Subshell, SubshellBindingTable};

type Result<T> = std::result::Result<T, ElectronDataError>;

fn invalid(what: &'static str, reason: String) -> ElectronDataError {
    ElectronDataError::Invalid { what, reason }
}

/// An inner-shell channel: one subshell of one element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InnerShell {
    /// Atomic number.
    pub z: u8,
    /// The subshell.
    pub subshell: Subshell,
    /// Binding energy `B`, eV.
    pub binding_energy_ev: f64,
}

impl InnerShell {
    /// The shell `subshell` of element `z`, with its binding energy from
    /// `table`.
    pub fn from_table(table: &SubshellBindingTable, z: u8, subshell: Subshell) -> Result<Self> {
        let binding_energy_ev = table
            .atom(z)
            .and_then(|a| a.binding_energy_ev(subshell))
            .ok_or_else(|| {
                invalid(
                    "inner shell",
                    format!("Z = {z} {} is not in the binding table", subshell.label()),
                )
            })?;
        Ok(Self {
            z,
            subshell,
            binding_energy_ev,
        })
    }
}

/// A loss channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    /// Valence excitation and ionization.
    Valence,
    /// Ionization of inner shell `i` (an index into
    /// [`ShellResolvedChannels::inner_shells`]).
    InnerShell(usize),
}

/// Inverse mean free path by channel, m⁻¹.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelInverseImfp {
    /// Valence channel.
    pub valence_per_m: f64,
    /// Inner-shell channels, in the order of
    /// [`ShellResolvedChannels::inner_shells`].
    pub shells_per_m: Vec<f64>,
}

impl ChannelInverseImfp {
    /// The sum over channels.
    pub fn total_per_m(&self) -> f64 {
        self.valence_per_m + self.shells_per_m.iter().sum::<f64>()
    }
}

/// DIIMFP by channel at one loss, m⁻¹ eV⁻¹.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelDiimfp {
    /// Valence channel.
    pub valence_per_m_ev: f64,
    /// Inner-shell channels, in the order of
    /// [`ShellResolvedChannels::inner_shells`].
    pub shells_per_m_ev: Vec<f64>,
}

impl ChannelDiimfp {
    /// The sum over channels.
    pub fn total_per_m_ev(&self) -> f64 {
        self.valence_per_m_ev + self.shells_per_m_ev.iter().sum::<f64>()
    }
}

/// A valence channel and inner-shell channels, each with its own optical ELF
/// (module docs).
#[derive(Debug, Clone)]
pub struct ShellResolvedChannels {
    valence: SinglePolePenn,
    shells: Vec<InnerShell>,
    models: Vec<SinglePolePenn>,
}

impl ShellResolvedChannels {
    /// Channels from a `valence` model and one optical ELF per inner shell,
    /// kept in the given order. The shell models take the valence model's
    /// Fermi energy, tolerance and exchange setting. Each binding energy must
    /// be finite and positive, each shell's ELF must start at or above its
    /// binding energy, and a `(Z, subshell)` may not repeat.
    pub fn new(valence: SinglePolePenn, shells: Vec<(InnerShell, OpticalElf)>) -> Result<Self> {
        let mut inner = Vec::with_capacity(shells.len());
        let mut models = Vec::with_capacity(shells.len());
        for (shell, elf) in shells {
            let b = shell.binding_energy_ev;
            let name = format!("Z = {} {}", shell.z, shell.subshell.label());
            if !(b.is_finite() && b > 0.0) {
                return Err(invalid(
                    "inner shell",
                    format!("{name}: binding energy must be finite and positive, got {b} eV"),
                ));
            }
            let first = elf.energy_range_ev().0;
            if first < b {
                return Err(invalid(
                    "inner-shell ELF",
                    format!(
                        "{name}: the ELF starts at {first} eV, below the binding energy {b} eV \
                         (a shell contributes only above its edge)"
                    ),
                ));
            }
            if inner
                .iter()
                .any(|s: &InnerShell| s.z == shell.z && s.subshell == shell.subshell)
            {
                return Err(invalid("inner shell", format!("{name} listed twice")));
            }
            models.push(valence.same_settings_for(elf));
            inner.push(shell);
        }
        Ok(Self {
            valence,
            shells: inner,
            models,
        })
    }

    /// The valence model.
    pub fn valence(&self) -> &SinglePolePenn {
        &self.valence
    }

    /// The inner-shell channels, in the order given.
    pub fn inner_shells(&self) -> &[InnerShell] {
        &self.shells
    }

    /// The model of inner shell `i` (its own ELF, the valence settings).
    pub fn shell_model(&self, i: usize) -> &SinglePolePenn {
        &self.models[i]
    }

    /// The channel-resolved inverse IMFP at kinetic energy `energy_ev` above
    /// the Fermi level: the valence model's `λ⁻¹`, and for shell `i`
    /// [`SinglePolePenn::imfp_and_stopping_with_binding`] of its model with
    /// its binding energy (module docs).
    pub fn inverse_imfps(&self, energy_ev: f64) -> Result<ChannelInverseImfp> {
        let valence_per_m = self
            .valence
            .imfp_and_stopping(energy_ev)?
            .inverse_imfp_per_m;
        let shells_per_m = self
            .shells
            .iter()
            .zip(&self.models)
            .map(|(s, m)| {
                m.imfp_and_stopping_with_binding(energy_ev, s.binding_energy_ev)
                    .map(|p| p.inverse_imfp_per_m)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(ChannelInverseImfp {
            valence_per_m,
            shells_per_m,
        })
    }

    /// The channel-resolved DIIMFP at kinetic energy `energy_ev` and loss
    /// `loss_ev` ([`SinglePolePenn::diimfp_with_binding_per_m_ev`] per
    /// channel).
    pub fn diimfps_per_m_ev(&self, energy_ev: f64, loss_ev: f64) -> Result<ChannelDiimfp> {
        let valence_per_m_ev = self.valence.diimfp_per_m_ev(energy_ev, loss_ev)?;
        let shells_per_m_ev = self
            .shells
            .iter()
            .zip(&self.models)
            .map(|(s, m)| m.diimfp_with_binding_per_m_ev(energy_ev, loss_ev, s.binding_energy_ev))
            .collect::<Result<Vec<_>>>()?;
        Ok(ChannelDiimfp {
            valence_per_m_ev,
            shells_per_m_ev,
        })
    }

    /// The channel of a loss `loss_ev` at kinetic energy `energy_ev`, for a
    /// uniform variate `u` in `[0, 1)`: each channel with probability
    /// proportional to its DIIMFP there, the inner shells first in order,
    /// then valence. `None` if every channel's DIIMFP is zero there.
    /// Deterministic in its arguments.
    pub fn sample_channel(&self, energy_ev: f64, loss_ev: f64, u: f64) -> Result<Option<Channel>> {
        let d = self.diimfps_per_m_ev(energy_ev, loss_ev)?;
        let total = d.total_per_m_ev();
        if total <= 0.0 {
            return Ok(None);
        }
        let mut acc = 0.0;
        for (i, p) in d.shells_per_m_ev.iter().enumerate() {
            acc += p / total;
            if *p > 0.0 && u < acc {
                return Ok(Some(Channel::InnerShell(i)));
            }
        }
        Ok(Some(Channel::Valence))
    }

    /// The kinetic energy `W = ω - B` handed to the secondary electron by an
    /// inner-shell loss `loss_ev`, eV (dV2022 eqs. (32), (34)); `None` for the
    /// valence channel (its secondary energy depends on the band structure,
    /// see `electron::secondary`) and below the edge.
    pub fn secondary_energy_ev(&self, channel: Channel, loss_ev: f64) -> Option<f64> {
        match channel {
            Channel::Valence => None,
            Channel::InnerShell(i) => {
                let w = loss_ev - self.shells.get(i)?.binding_energy_ev;
                (w >= 0.0).then_some(w)
            }
        }
    }
}
