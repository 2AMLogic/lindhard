//! Attribution of energy losses to inner-shell ionization above subshell
//! edges.
//!
//! # Source and what is ours
//!
//! The structure follows the open account in M. A. Quinto *et al.*, Int. J.
//! Mol. Sci. 23, 6121 (2022), doi:10.3390/ijms23116121 (PMC9181504), Section
//! 3, eqs. (32) and (35): the secondary-electron spectrum of a material is the
//! sum of an outer-shell (valence) term, built from the loss function at
//! energy `W + B(T)`, and one term per inner shell `j`, built from that
//! shell's loss function at `W + B_j`, where `B_j` is the shell's binding
//! energy and `W` the kinetic energy of the emitted electron. So a loss
//! `ω = W + B_j` attributed to shell `j` hands the secondary the energy
//! `W = ω - B_j` ([`ChannelPartition::secondary_energy_ev`]). Shell `j`
//! contributes only for `ω >= B_j`.
//!
//! **The partition of one optical ELF between channels is our own model
//! assumption, not a published rule.** The cited account has a separate loss
//! function per shell (from atomic or TDDFT data), which a single optical ELF
//! does not provide. Here, at the loss `ω`, the optical ELF is divided among
//! all subshells with `B <= ω` in proportion to `x_a n_s`, the number of
//! atoms `x_a` of the element per formula unit times the subshell's ground
//! state occupancy `n_s` (from [`SubshellBindingTable`]): shell `s` takes the
//! fraction `x_a n_s / Σ_{B_t <= ω} x_b n_t`. Subshells with `B` below a
//! user-given valence cutoff stay in the valence channel; the others are the
//! inner-shell channels. This counting rule makes each channel's share
//! proportional to the electrons that can take part, nothing more; it
//! ignores that atomic oscillator strength is not proportional to occupancy
//! near an edge. The table must therefore list the valence subshells too
//! (EADL does), or the valence channel is starved above the first inner
//! edge. The partition is applied to the energy loss `ω` of the DIIMFP, not to
//! the Penn plasma frequency `ω₀(q, ω)`.
//!
//! Channel inverse mean free paths are the integrals of
//! `f_channel(ω) p(T, ω)` over the allowed losses, where `p` is the DIIMFP of
//! [`SinglePolePenn`], with the optional exchange factor `1 + F` and the
//! indistinguishability limit `ω <= (T' + B)/2` (valence: `B = 0`) of the
//! source, eq. (35). They add up to the total by construction.

use super::penn::{hartree_ev, SinglePolePenn};
use super::quadrature::{integrate_segments, GaussLegendre};
use crate::constants::BOHR_RADIUS;
use crate::electron::data::{ElectronDataError, Subshell, SubshellBindingTable};

type Result<T> = std::result::Result<T, ElectronDataError>;

/// An inner-shell channel: one subshell of one element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InnerShell {
    /// Atomic number.
    pub z: u8,
    /// The subshell.
    pub subshell: Subshell,
    /// Binding energy `E_B`, eV.
    pub binding_energy_ev: f64,
    /// Partition weight `x_a n_s` (atoms per formula unit times occupancy).
    pub weight: f64,
}

/// A loss channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    /// Valence excitation and ionization (and every subshell below the
    /// valence cutoff).
    Valence,
    /// Ionization of inner shell `i` (an index into
    /// [`ChannelPartition::inner_shells`]).
    InnerShell(usize),
}

/// Inverse mean free path by channel, m⁻¹.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelInverseImfp {
    /// Valence channel.
    pub valence_per_m: f64,
    /// Inner-shell channels, in the order of
    /// [`ChannelPartition::inner_shells`].
    pub shells_per_m: Vec<f64>,
}

impl ChannelInverseImfp {
    /// The sum over channels.
    pub fn total_per_m(&self) -> f64 {
        self.valence_per_m + self.shells_per_m.iter().sum::<f64>()
    }
}

/// The partition of an optical ELF between the valence channel and the
/// inner-shell channels of a target (module docs).
#[derive(Debug, Clone)]
pub struct ChannelPartition {
    valence_cutoff_ev: f64,
    inner: Vec<InnerShell>,
    /// Every subshell `(E_B, weight)`, ascending in `E_B`.
    all: Vec<(f64, f64)>,
    /// Prefix sums of the weight over `all`, and over its inner members.
    cum_all: Vec<f64>,
    cum_inner: Vec<f64>,
    /// Per inner shell, its position in `all`.
    inner_pos: Vec<usize>,
}

impl ChannelPartition {
    /// The partition for a target with `composition` = `(Z, atoms per formula
    /// unit)`, with binding energies from `table`; subshells with
    /// `E_B < valence_cutoff_ev` are valence. Every `Z` must be tabulated and
    /// the amounts finite and positive; `Z` may not repeat.
    pub fn new(
        table: &SubshellBindingTable,
        composition: &[(u8, f64)],
        valence_cutoff_ev: f64,
    ) -> Result<Self> {
        let bad = |what: &'static str, reason: String| ElectronDataError::Invalid { what, reason };
        if !(valence_cutoff_ev.is_finite() && valence_cutoff_ev >= 0.0) {
            return Err(bad(
                "valence cutoff",
                format!("must be finite and non-negative, got {valence_cutoff_ev} eV"),
            ));
        }
        if composition.is_empty() {
            return Err(bad("composition", "is empty".into()));
        }
        // (E_B, z, subshell, weight, is_inner)
        let mut rows: Vec<(f64, u8, Subshell, f64, bool)> = Vec::new();
        for (k, &(z, amount)) in composition.iter().enumerate() {
            if !(amount.is_finite() && amount > 0.0) {
                return Err(bad(
                    "composition",
                    format!("Z = {z}: amount must be finite and positive, got {amount}"),
                ));
            }
            if composition[..k].iter().any(|&(zz, _)| zz == z) {
                return Err(bad("composition", format!("Z = {z} listed twice")));
            }
            let atom = table.atom(z).ok_or_else(|| {
                bad(
                    "composition",
                    format!("Z = {z} is not in the binding table"),
                )
            })?;
            for s in atom.shells() {
                let eb = s.binding_energy_ev();
                rows.push((
                    eb,
                    z,
                    s.subshell(),
                    amount * s.occupancy(),
                    eb >= valence_cutoff_ev,
                ));
            }
        }
        rows.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut inner = Vec::new();
        let mut inner_pos = Vec::new();
        let mut all = Vec::new();
        let (mut cum_all, mut cum_inner) = (Vec::new(), Vec::new());
        let (mut ta, mut ti) = (0.0, 0.0);
        for (i, &(eb, z, subshell, weight, is_inner)) in rows.iter().enumerate() {
            ta += weight;
            if is_inner {
                ti += weight;
                inner.push(InnerShell {
                    z,
                    subshell,
                    binding_energy_ev: eb,
                    weight,
                });
                inner_pos.push(i);
            }
            all.push((eb, weight));
            cum_all.push(ta);
            cum_inner.push(ti);
        }
        Ok(Self {
            valence_cutoff_ev,
            inner,
            all,
            cum_all,
            cum_inner,
            inner_pos,
        })
    }

    /// The valence cutoff, eV.
    pub fn valence_cutoff_ev(&self) -> f64 {
        self.valence_cutoff_ev
    }

    /// The inner-shell channels, ascending in binding energy (ties by `Z`,
    /// then subshell).
    pub fn inner_shells(&self) -> &[InnerShell] {
        &self.inner
    }

    /// Number of subshells with `E_B <= loss_ev`.
    fn active(&self, loss_ev: f64) -> usize {
        self.all.partition_point(|&(eb, _)| eb <= loss_ev)
    }

    /// The share of the ELF at loss `loss_ev` that goes to inner shell `i`
    /// (zero below its edge).
    pub fn shell_fraction(&self, loss_ev: f64, i: usize) -> f64 {
        let n = self.active(loss_ev);
        if n == 0 || self.inner_pos[i] >= n {
            return 0.0;
        }
        self.inner[i].weight / self.cum_all[n - 1]
    }

    /// The share that stays in the valence channel (one if no subshell is
    /// open yet).
    pub fn valence_fraction(&self, loss_ev: f64) -> f64 {
        let n = self.active(loss_ev);
        if n == 0 {
            return 1.0;
        }
        let total = self.cum_all[n - 1];
        ((total - self.cum_inner[n - 1]) / total).max(0.0)
    }

    /// The channel for a uniform variate `u` in `[0, 1)`: the inner shells
    /// open at `loss_ev` in order, then the valence channel takes the rest.
    /// Deterministic in `(loss_ev, u)`.
    pub fn sample_channel(&self, loss_ev: f64, u: f64) -> Channel {
        let mut acc = 0.0;
        for i in 0..self.inner.len() {
            acc += self.shell_fraction(loss_ev, i);
            if u < acc {
                return Channel::InnerShell(i);
            }
        }
        Channel::Valence
    }

    /// The kinetic energy `W = ω - E_B` handed to the secondary electron by
    /// an inner-shell loss `loss_ev`, eV; `None` for the valence channel
    /// (its secondary energy depends on the band structure, see
    /// `electron::secondary`) and below the edge.
    pub fn secondary_energy_ev(&self, channel: Channel, loss_ev: f64) -> Option<f64> {
        match channel {
            Channel::Valence => None,
            Channel::InnerShell(i) => {
                let w = loss_ev - self.inner[i].binding_energy_ev;
                (w >= 0.0).then_some(w)
            }
        }
    }

    /// The channel-resolved inverse IMFP of `model` at kinetic energy
    /// `energy_ev` above the Fermi level (module docs). With the exchange
    /// correction enabled and applicable at this energy the valence channel
    /// is limited to `ω <= T'/2` and shell `i` to `ω <= (T' + E_B)/2`; without
    /// it every channel reaches `ω = T`, and the total is the model's
    /// `λ⁻¹` ([`SinglePolePenn::imfp_and_stopping`]) up to the integration
    /// tolerance.
    pub fn inverse_imfps(
        &self,
        model: &SinglePolePenn,
        energy_ev: f64,
    ) -> Result<ChannelInverseImfp> {
        if !(energy_ev.is_finite() && energy_ev > 0.0) {
            return Err(ElectronDataError::Invalid {
                what: "electron energy",
                reason: format!("must be finite and positive, got {energy_ev} eV"),
            });
        }
        let h = hartree_ev();
        let t_au = energy_ev / h;
        let tp = energy_ev + model.fermi_energy_ev();
        let exchange = model.exchange_applies_at(energy_ev);
        let gl = GaussLegendre::new(5);
        let tol = if exchange {
            (model.relative_tolerance() * 100.0).min(1e-2)
        } else {
            model.relative_tolerance()
        };
        let knots = model.optical_elf().energy_ev();
        let edges: Vec<f64> = self.all.iter().map(|&(eb, _)| eb).collect();

        let integrate = |lo: f64, hi: f64, frac: &dyn Fn(f64) -> f64| -> f64 {
            if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
                return 0.0;
            }
            let mut breaks: Vec<f64> = knots
                .iter()
                .chain(edges.iter())
                .copied()
                .filter(|&x| x > lo && x < hi)
                .collect();
            breaks.sort_by(f64::total_cmp);
            breaks.dedup();
            breaks.insert(0, lo);
            breaks.push(hi);
            integrate_segments(
                &gl,
                &mut |w: f64| {
                    let f = frac(w);
                    if f == 0.0 {
                        return [0.0];
                    }
                    // per bohr per hartree -> per m per eV; dω is in eV.
                    let p = model.diimfp_core_au(t_au, w / h, exchange) / (BOHR_RADIUS * h);
                    [f * p]
                },
                &breaks,
                tol,
            )[0]
        };

        let top = |b: f64| {
            if exchange {
                (0.5 * (tp + b)).min(energy_ev)
            } else {
                energy_ev
            }
        };
        let valence_per_m = integrate(0.0, top(0.0), &|w| self.valence_fraction(w));
        let shells_per_m = self
            .inner
            .iter()
            .enumerate()
            .map(|(i, s)| {
                integrate(s.binding_energy_ev, top(s.binding_energy_ev), &|w| {
                    self.shell_fraction(w, i)
                })
            })
            .collect();
        Ok(ChannelInverseImfp {
            valence_per_m,
            shells_per_m,
        })
    }
}
