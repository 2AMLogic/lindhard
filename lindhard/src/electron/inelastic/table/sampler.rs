//! The momentum-transfer sampler of the inelastic table (parent module
//! docs).

use super::super::penn::{hartree_ev, SinglePolePenn};
use super::density::{adaptive_density, LinearDensity};
use super::{check_no_exchange, invalid, InelasticTableError, DEFAULT_Q_TOLERANCE};
use crate::constants::BOHR_RADIUS;

/// Samples the momentum transfer `q` of an inelastic collision with incident
/// kinetic energy `E` and loss `W` (module docs). Built once per `(E, W)`.
#[derive(Debug, Clone, PartialEq)]
pub struct MomentumTransferSampler {
    density: LinearDensity,
    /// Hartree-unit wavenumbers `k` and `k'`, and `T' = E + E_F`.
    k: f64,
    kp: f64,
    tp: f64,
}

impl MomentumTransferSampler {
    /// Build the sampler with the default tolerance [`DEFAULT_Q_TOLERANCE`].
    pub fn new(
        penn: &SinglePolePenn,
        energy_ev: f64,
        loss_ev: f64,
    ) -> Result<Self, InelasticTableError> {
        Self::with_tolerance(penn, energy_ev, loss_ev, DEFAULT_Q_TOLERANCE)
    }

    /// Build the sampler with density-model tolerance `tol` in
    /// `[1e-9, 1e-2]`. Fails if `0 < W <= E` does not hold or the slice has
    /// no support (`W` below the ELF table, or no `q` allowed).
    pub fn with_tolerance(
        penn: &SinglePolePenn,
        energy_ev: f64,
        loss_ev: f64,
        tol: f64,
    ) -> Result<Self, InelasticTableError> {
        check_no_exchange(penn)?;
        if !(energy_ev.is_finite() && energy_ev > 0.0) {
            return Err(invalid("electron energy", format!("{energy_ev} eV")));
        }
        if !(loss_ev.is_finite() && loss_ev > 0.0 && loss_ev <= energy_ev) {
            return Err(invalid(
                "energy loss",
                format!("{loss_ev} eV (need 0 < W <= E = {energy_ev} eV)"),
            ));
        }
        if !(tol.is_finite() && (1.0e-9..=1.0e-2).contains(&tol)) {
            return Err(invalid(
                "q tolerance",
                format!("{tol} (need 1e-9 ..= 1e-2)"),
            ));
        }
        let h = hartree_ev();
        let (t, w) = (energy_ev / h, loss_ev / h);
        let breaks = penn.q_slice_breaks_au(t, w).ok_or_else(|| {
            invalid(
                "momentum-transfer slice",
                format!("no allowed q at E = {energy_ev} eV, W = {loss_ev} eV"),
            )
        })?;
        let f = |u: f64| penn.q_integrand_au(w, u);
        // The node bound is an error for the loss rows of a table (#345); a
        // momentum-transfer slice keeps the density it reached, as before.
        let (density, _) = adaptive_density(&f, &breaks, 8, tol, false, &[], 0.0);
        if !(density.total() > 0.0 && density.total().is_finite()) {
            return Err(invalid(
                "momentum-transfer slice",
                format!("zero loss function at E = {energy_ev} eV, W = {loss_ev} eV"),
            ));
        }
        let tp = t + penn.fermi_energy_ev() / h;
        Ok(Self {
            density,
            k: (2.0 * tp).sqrt(),
            kp: (2.0 * (tp - w)).max(0.0).sqrt(),
            tp,
        })
    }

    /// The momentum transfer (as a wavenumber, m⁻¹) at cumulative
    /// probability `u` in `[0, 1]` (clamped).
    pub fn sample(&self, u: f64) -> f64 {
        self.density.quantile(u.clamp(0.0, 1.0)).exp() / BOHR_RADIUS
    }

    /// The support `(q_min, q_max)` of the sampler, m⁻¹: the kinematic
    /// limits `q±` cut to the tabulated ELF range.
    pub fn q_range_per_m(&self) -> (f64, f64) {
        let n = self.density.x.len();
        (
            self.density.x[0].exp() / BOHR_RADIUS,
            self.density.x[n - 1].exp() / BOHR_RADIUS,
        )
    }

    /// The kinematic limits `q- = k - k'` and `q+ = k + k'`, m⁻¹.
    pub fn kinematic_limits_per_m(&self) -> (f64, f64) {
        (
            (self.k - self.kp) / BOHR_RADIUS,
            (self.k + self.kp) / BOHR_RADIUS,
        )
    }

    /// The cosine of the primary's deflection angle for momentum transfer
    /// `q` (m⁻¹), `(k² + k'² - q²)/(2 k k')`, clamped to `[-1, 1]`.
    pub fn cos_theta(&self, q_per_m: f64) -> f64 {
        let q = q_per_m * BOHR_RADIUS;
        let c = (self.k * self.k + self.kp * self.kp - q * q) / (2.0 * self.k * self.kp);
        c.clamp(-1.0, 1.0)
    }

    /// The DIIMFP `p(E, W)` (m⁻¹ eV⁻¹) of the sampler's own density model:
    /// the integral of the density over `ln q`, divided by `π T'`. It agrees
    /// with [`SinglePolePenn::diimfp_per_m_ev`] to the model tolerance.
    pub fn diimfp_per_m_ev(&self) -> f64 {
        let h = hartree_ev();
        self.density.total() / (std::f64::consts::PI * self.tp) / (BOHR_RADIUS * h)
    }
}
