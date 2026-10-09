//! Optional corrections to the static potential: local exchange and
//! correlation-polarization (issue #91, step 3 of the Mott elastic model).
//!
//! Both are **off by default**. They need the atomic electron density
//! `rho(r)`, which a [`ScreenedPotential`] does not provide, so the density is a
//! separate, explicit input ([`ElectronDensity`]). For the analytic potentials
//! of this module the density is the one that generates the potential through
//! Poisson's equation (Salvat, Martinez, Mayol & Parellada, Phys. Rev. A 36, 467
//! (1987), doi:10.1103/PhysRevA.36.467, Eqs. (2) and (12)); it is never
//! inferred from an arbitrary user potential.
//!
//! # Sources
//!
//! The model is the one of F. Salvat, Phys. Rev. A 68, 012708 (2003),
//! doi:10.1103/PhysRevA.68.012708, Sec. II (read in the open copy of the
//! University of Barcelona repository, diposit.ub.edu). That paper is the
//! published account of the optical-model potential that the ELSEPA paper
//! (Salvat, Jablonski & Powell, Comput. Phys. Commun. 165, 157 (2005)) uses;
//! the ELSEPA paper itself is closed access and was not opened, and no ELSEPA
//! code was read (Tier C).
//!
//! * **Exchange**, Furness & McCarthy, J. Phys. B 6, 2280 (1973),
//!   doi:10.1088/0022-3700/6/11/021, as written in Salvat (2003) Eq. (2) (and
//!   identically in Jablonski, Salvat & Powell, J. Phys. Chem. Ref. Data 33,
//!   409 (2004), Eq. (15)), in hartree atomic units:
//!
//!   ```text
//!   V_ex(r) = 1/2 [E - V_st(r)] - 1/2 { [E - V_st(r)]^2 + 4 pi rho(r) }^(1/2)
//!   ```
//!
//!   with `E` the kinetic energy of the projectile and `V_st` the static
//!   potential energy. **The equation number in the 1973 paper itself is not
//!   cited**: that paper is closed access (IOP subscription; Unpaywall lists
//!   only a metadata record), so only the two secondary accounts above were
//!   read. We evaluate the algebraically identical form
//!   `-2 pi rho / (D + sqrt(D^2 + 4 pi rho))`, `D = E - V_st >= 0`, which avoids
//!   the cancellation of the printed form when `4 pi rho << D^2` (high energy,
//!   outer atom). `V_ex <= 0` everywhere and `V_ex = 0` where `rho = 0`.
//!   Salvat (2003) quotes Bransden et al. (1976) that it is accurate above
//!   about 1 hartree.
//! * **Correlation-polarization**, Salvat (2003) Eqs. (3)-(7) and (9):
//!   - long-range polarization (Buckingham) `V_pol = -alpha_p/(2 (r^2 + d^2)^2)`,
//!     Eq. (3), with `d^4 = (1/2) alpha_p Z^(-1/3) b_pol^2` (bohr units),
//!     Eq. (4), and Seltzer's empirical `b_pol^2 = (E - 50 eV)/(16 eV)`,
//!     Eq. (5) (proposed for `E >= 100 eV`);
//!   - the local-density correlation potential of Perdew & Zunger, Phys. Rev.
//!     B 23, 5048 (1981), in the parametrization printed as Salvat (2003)
//!     Eqs. (6)-(7), with `r_s = (3/(4 pi rho))^(1/3)`:
//!     `V_co = 0.0311 ln r_s - 0.0584 + 0.00133 r_s ln r_s - 0.0084 r_s`
//!     for `r_s < 1` and
//!     `V_co = -beta_0 (1 + 7/6 beta_1 r_s^(1/2) + 4/3 beta_2 r_s)/(1 + beta_1 r_s^(1/2) + beta_2 r_s)^2`
//!     for `r_s >= 1`, `beta_0 = 0.1423`, `beta_1 = 1.0529`, `beta_2 = 0.3334`.
//!     **Sign of Eq. (7a).** Salvat prints Eq. (7a) with an overall minus
//!     sign in front of the bracket. With that sign the two branches jump from
//!     `+0.0668` to `-0.0668` hartree at `r_s = 1` and the high-density
//!     correlation potential would be repulsive. We use the sign that makes
//!     the branches continuous (tested), which is also the sign of the same
//!     Perdew-Zunger expression `A ln r_s + (B - A/3) + (2/3) C r_s ln r_s +
//!     (1/3)(2D - C) r_s` with `A = 0.0311`, `B = -0.048`, `C = 0.002`,
//!     `D = -0.0116` as reproduced (citing Perdew & Zunger) on the Delta
//!     Science Institute DFT course page
//!     <https://www.dsedu.org/courses/dft/xc>, read 2026-10-07. The 1981 paper
//!     itself was not opened (APS bot challenge; Unpaywall's link is the
//!     publisher PDF);
//!   - the join, Eq. (9): `V_cp = max(V_co, V_pol)` for `r < r_cp` and
//!     `V_cp = V_pol` for `r >= r_cp`, where `r_cp` is the outer radius at
//!     which `V_co` and `V_pol` cross.
//!
//! # Numerical choices (ours, not from the papers)
//!
//! * **Finding `r_cp`.** `V_co - V_pol` is sampled on a geometric grid (ratio
//!   1.01) from [`JOIN_SCAN_FLOOR`] to the outer radius and the outermost sign
//!   change is refined by bisection. Two crossings closer than one grid ratio
//!   would be missed. If `V_co > V_pol` on the whole grid the crossing lies
//!   below the floor (`V_co` diverges only logarithmically at the nucleus);
//!   then Eq. (9) is `V_pol` on the whole grid and [`r_cp`] is reported as
//!   `None`. If `V_co < V_pol` at the outer radius, the outer crossing lies
//!   beyond it and construction fails.
//! * **Outer radius.** The `r^-4` tail is integrated out to
//!   [`CorrelationPolarization::outer_radius`] and set to zero beyond it (a
//!   breakpoint and the matching radius). The neglected tail is not covered by
//!   the `|r V| < threshold` rule of [`ScreenedPotential::matching_radius`],
//!   which would put the matching radius at thousands of bohr; its effect is
//!   measured by the convergence tests instead (doubling the outer radius).
//!   The matching radius of the corrected potential is the larger of the outer
//!   radius and the default scan applied to `V_st + V_ex` (both `<= 0`, so no
//!   cancellation can stop the scan early).
//! * **Start radius.** With the polarization tail, partial waves of `l` up to
//!   about `k` times the outer radius have their turning points where the
//!   potential is not negligible. The step-1 start rule is not accurate there
//!   (measured: errors up to 0.25 rad at `l ~ 1250`, Cu, 10 keV, 100 bohr), so
//!   [`CorrectedPotential`] reports [`ScreenedPotential::long_range`] and the
//!   solver uses a WKB start criterion instead. The uncorrected path keeps
//!   the step-1 rule and its bits.
//!
//! [`r_cp`]: CorrelationPolarizationInfo::join_radius_bohr

use super::{
    invalid, solve, ElasticError, ElasticResult, ElasticSolver, SalvatDhfs, ScreenedPotential,
    SolverOptions, Yukawa, HARTREE_EV,
};
use std::f64::consts::PI;

/// Model identifier recorded in the metadata when exchange is on.
pub const EXCHANGE_MODEL: &str =
    "Furness-McCarthy local exchange, J. Phys. B 6, 2280 (1973), as Salvat, PRA 68, 012708 (2003) Eq. (2)";

/// Model identifier recorded in the metadata when correlation-polarization is on.
pub const CORRELATION_POLARIZATION_MODEL: &str =
    "Salvat, PRA 68, 012708 (2003) Eqs. (3)-(7), (9): \
     Perdew-Zunger LDA correlation joined to a Buckingham polarization tail";

/// Lower end of the grid on which the join radius `r_cp` is searched, bohr.
pub const JOIN_SCAN_FLOOR: f64 = 1.0e-8;

/// Default outer radius of the polarization tail, bohr (see the module docs).
pub const DEFAULT_OUTER_RADIUS: f64 = 50.0;

/// Atomic electron density, the input the corrections need beyond the
/// potential.
pub trait ElectronDensity: Sync {
    /// Electron number density at radius `r` bohr, electrons per bohr^3
    /// (`>= 0`).
    fn density(&self, r: f64) -> f64;

    /// Where the density comes from (recorded in the result metadata).
    fn density_source(&self) -> String;
}

/// Poisson density of the analytic screening function, Salvat et al. (1987)
/// Eq. (12): `rho(r) = (Z/(4 pi r)) sum_i A_i alpha_i^2 exp(-alpha_i r)`,
/// summed over the terms in use (two for the asterisked, `A_3 = 0` rows of
/// Table I).
impl ElectronDensity for SalvatDhfs {
    fn density(&self, r: f64) -> f64 {
        let s: f64 = self
            .amplitudes()
            .iter()
            .zip(self.alphas())
            .map(|(a, al)| a * al * al * (-al * r).exp())
            .sum();
        self.z * s / (4.0 * PI * r)
    }
    fn density_source(&self) -> String {
        format!(
            "Poisson density of the Salvat et al. (1987) analytic screening function, \
             PRA 36, 467 Eq. (12), Z={}, A={:?}, alpha={:?} 1/bohr",
            self.z,
            self.amplitudes(),
            self.alphas()
        )
    }
}

/// Poisson density of `V = -Z exp(-r/a)/r` (the one-term case of Salvat et
/// al. (1987) Eq. (12)): `rho = Z exp(-r/a)/(4 pi a^2 r)`. A fixture density.
impl ElectronDensity for Yukawa {
    fn density(&self, r: f64) -> f64 {
        let a = self.screening_length;
        self.z * (-r / a).exp() / (4.0 * PI * a * a * r)
    }
    fn density_source(&self) -> String {
        format!(
            "Poisson density of a Yukawa potential (fixture), Z={}, a={} bohr",
            self.z, self.screening_length
        )
    }
}

/// How the cutoff `d` of the polarization potential is set (Salvat 2003
/// Eqs. (4)-(5)).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PolarizationCutoff {
    /// Seltzer's recipe `b_pol^2 = (E - 50 eV)/(16 eV)`, Eq. (5); needs
    /// `E > 50 eV` (proposed for `E >= 100 eV`).
    Seltzer,
    /// A caller-chosen `b_pol^2 >= 0` (dimensionless) in Eq. (4).
    BPolSquared(f64),
}

/// Inputs of the correlation-polarization correction.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationPolarization {
    /// Static dipole polarizability `alpha_p` of the target atom, bohr^3
    /// (atomic units of polarizability).
    pub polarizability: f64,
    /// Citation of `polarizability` (recorded in the metadata; must not be
    /// empty).
    pub polarizability_source: String,
    /// Cutoff parameter of the Buckingham potential.
    pub cutoff: PolarizationCutoff,
    /// Radius beyond which the polarization tail is dropped, bohr (a numerical
    /// setting, see the module docs).
    pub outer_radius: f64,
}

impl CorrelationPolarization {
    /// Seltzer's cutoff and [`DEFAULT_OUTER_RADIUS`].
    pub fn new(polarizability: f64, polarizability_source: impl Into<String>) -> Self {
        Self {
            polarizability,
            polarizability_source: polarizability_source.into(),
            cutoff: PolarizationCutoff::Seltzer,
            outer_radius: DEFAULT_OUTER_RADIUS,
        }
    }
}

/// Which corrections to add to the static potential. The default is none.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Corrections {
    /// Furness-McCarthy exchange on/off.
    pub exchange: bool,
    /// Correlation-polarization, `None` for off.
    pub correlation_polarization: Option<CorrelationPolarization>,
}

impl Corrections {
    /// Both corrections off.
    pub fn none() -> Self {
        Self::default()
    }

    /// True when neither correction is on.
    pub fn is_none(&self) -> bool {
        !self.exchange && self.correlation_polarization.is_none()
    }
}

/// Metadata of the exchange correction.
#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeInfo {
    /// Model identifier ([`EXCHANGE_MODEL`]).
    pub model: &'static str,
}

/// Metadata of the correlation-polarization correction.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationPolarizationInfo {
    /// Model identifier ([`CORRELATION_POLARIZATION_MODEL`]).
    pub model: &'static str,
    /// Dipole polarizability used, bohr^3.
    pub polarizability_bohr3: f64,
    /// Its citation, as supplied.
    pub polarizability_source: String,
    /// `b_pol^2` used in Eq. (4), dimensionless.
    pub b_pol_squared: f64,
    /// Cutoff `d` of Eq. (3), bohr.
    pub cutoff_d_bohr: f64,
    /// Outer crossing `r_cp` of Eq. (9), bohr; `None` when it lies below
    /// [`JOIN_SCAN_FLOOR`] (then `V_cp = V_pol` on the whole grid).
    pub join_radius_bohr: Option<f64>,
    /// Radius beyond which the tail is dropped, bohr.
    pub outer_radius_bohr: f64,
}

/// Which corrections a result includes. The default reports both off.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CorrectionMetadata {
    /// `Some` when exchange is on.
    pub exchange: Option<ExchangeInfo>,
    /// `Some` when correlation-polarization is on.
    pub correlation_polarization: Option<CorrelationPolarizationInfo>,
    /// [`ElectronDensity::density_source`] of the density used, when any
    /// correction is on.
    pub density_source: Option<String>,
}

#[derive(Debug, Clone, Copy)]
struct CpState {
    alpha: f64,
    d2: f64,
    r_cp: f64,
    r_out: f64,
}

/// The static potential plus the selected corrections, **at one kinetic
/// energy** (the exchange and the polarization cutoff depend on it). Build the
/// solver with [`CorrectedPotential::solver`] so the solver energy cannot
/// differ from the one the potential was built for.
pub struct CorrectedPotential<'a> {
    stat: &'a dyn ScreenedPotential,
    dens: &'a dyn ElectronDensity,
    energy_ev: f64,
    e: f64,
    exchange: bool,
    cp: Option<CpState>,
    meta: CorrectionMetadata,
}

/// Perdew-Zunger correlation potential, Salvat (2003) Eqs. (6)-(7) (sign of
/// (7a) as discussed in the module docs), hartree, for density `rho`
/// (bohr^-3). Zero for `rho <= 0` (the `r_s -> infinity` limit).
pub fn correlation_potential(rho: f64) -> f64 {
    if rho.is_nan() || rho <= 0.0 {
        return 0.0;
    }
    let rs = (3.0 / (4.0 * PI * rho)).cbrt();
    if rs < 1.0 {
        let l = rs.ln();
        0.0311 * l - 0.0584 + 0.00133 * rs * l - 0.0084 * rs
    } else {
        let (b0, b1, b2) = (0.1423, 1.0529, 0.3334);
        let sq = rs.sqrt();
        let den = 1.0 + b1 * sq + b2 * rs;
        -b0 * (1.0 + 7.0 / 6.0 * b1 * sq + 4.0 / 3.0 * b2 * rs) / (den * den)
    }
}

/// Furness-McCarthy exchange potential, Salvat (2003) Eq. (2), hartree, for
/// local kinetic-energy argument `d = E - V_st` (hartree) and density `rho`
/// (bohr^-3), in the cancellation-free form of the module docs.
pub fn exchange_potential(d: f64, rho: f64) -> f64 {
    let rho = rho.max(0.0);
    let s = (d * d + 4.0 * PI * rho).sqrt();
    if d >= 0.0 {
        if rho == 0.0 {
            return 0.0;
        }
        -2.0 * PI * rho / (d + s)
    } else {
        0.5 * (d - s)
    }
}

impl<'a> CorrectedPotential<'a> {
    /// Build the corrected potential of `stat` with density `dens` at kinetic
    /// energy `energy_ev` (eV, > 0).
    pub fn new(
        stat: &'a dyn ScreenedPotential,
        dens: &'a dyn ElectronDensity,
        energy_ev: f64,
        corrections: &Corrections,
    ) -> Result<Self, ElasticError> {
        if !energy_ev.is_finite() || energy_ev <= 0.0 {
            return Err(invalid(
                "energy",
                format!("{energy_ev} eV (need finite, > 0)"),
            ));
        }
        let e = energy_ev / HARTREE_EV;
        let mut meta = CorrectionMetadata::default();
        if corrections.exchange {
            meta.exchange = Some(ExchangeInfo {
                model: EXCHANGE_MODEL,
            });
        }
        let mut this = Self {
            stat,
            dens,
            energy_ev,
            e,
            exchange: corrections.exchange,
            cp: None,
            meta,
        };
        if let Some(cp) = &corrections.correlation_polarization {
            let (state, info) = this.build_cp(cp)?;
            this.cp = Some(state);
            this.meta.correlation_polarization = Some(info);
        }
        if !corrections.is_none() {
            this.meta.density_source = Some(dens.density_source());
        }
        Ok(this)
    }

    fn build_cp(
        &self,
        cp: &CorrelationPolarization,
    ) -> Result<(CpState, CorrelationPolarizationInfo), ElasticError> {
        let alpha = cp.polarizability;
        if !alpha.is_finite() || alpha <= 0.0 {
            return Err(invalid(
                "dipole polarizability",
                format!("{alpha} bohr^3 (need finite, > 0)"),
            ));
        }
        if cp.polarizability_source.trim().is_empty() {
            return Err(invalid(
                "polarizability source",
                "a citation for the polarizability is required",
            ));
        }
        if !cp.outer_radius.is_finite() || cp.outer_radius <= 0.0 {
            return Err(invalid(
                "outer radius",
                format!("{} bohr (need finite, > 0)", cp.outer_radius),
            ));
        }
        let b2 = match cp.cutoff {
            PolarizationCutoff::Seltzer => {
                if self.energy_ev <= 50.0 {
                    return Err(invalid(
                        "polarization cutoff",
                        format!(
                            "Seltzer's b_pol^2 = (E - 50 eV)/16 eV needs E > 50 eV, got {} eV",
                            self.energy_ev
                        ),
                    ));
                }
                (self.energy_ev - 50.0) / 16.0
            }
            PolarizationCutoff::BPolSquared(b2) => {
                if !b2.is_finite() || b2 < 0.0 {
                    return Err(invalid("b_pol^2", format!("{b2} (need finite, >= 0)")));
                }
                b2
            }
        };
        let z = self.stat.nuclear_charge();
        if !(z.is_finite() && z > 0.0) {
            return Err(invalid(
                "nuclear charge",
                "the polarization cutoff of Eq. (4) needs the atomic number Z > 0",
            ));
        }
        let d4 = 0.5 * alpha * z.powf(-1.0 / 3.0) * b2;
        let d2 = d4.sqrt();
        let r_out = cp.outer_radius;
        let vpol = |r: f64| -alpha / (2.0 * (r * r + d2).powi(2));
        let f = |r: f64| correlation_potential(self.dens.density(r)) - vpol(r);
        if f(r_out) < 0.0 {
            return Err(invalid(
                "outer radius",
                format!(
                    "{r_out} bohr: the correlation potential is still below the polarization \
                     tail there, so the outer crossing r_cp of Eq. (9) lies beyond it"
                ),
            ));
        }
        // Outermost sign change of f on a geometric grid.
        let mut grid = Vec::new();
        let mut r = JOIN_SCAN_FLOOR;
        while r < r_out {
            grid.push(r);
            r *= 1.01;
        }
        grid.push(r_out);
        let vals: Vec<f64> = grid.iter().map(|&r| f(r)).collect();
        let mut r_cp = None;
        for i in (1..grid.len()).rev() {
            if (vals[i - 1] < 0.0) != (vals[i] < 0.0) {
                let (mut lo, mut hi) = (grid[i - 1], grid[i]);
                let neg_lo = vals[i - 1] < 0.0;
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if (f(mid) < 0.0) == neg_lo {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                    if hi - lo <= 1e-15 * hi {
                        break;
                    }
                }
                r_cp = Some(0.5 * (lo + hi));
                break;
            }
        }
        let state = CpState {
            alpha,
            d2,
            r_cp: r_cp.unwrap_or(JOIN_SCAN_FLOOR),
            r_out,
        };
        let info = CorrelationPolarizationInfo {
            model: CORRELATION_POLARIZATION_MODEL,
            polarizability_bohr3: alpha,
            polarizability_source: cp.polarizability_source.clone(),
            b_pol_squared: b2,
            cutoff_d_bohr: d2.sqrt(),
            join_radius_bohr: r_cp,
            outer_radius_bohr: r_out,
        };
        Ok((state, info))
    }

    /// Kinetic energy the potential was built for, eV.
    pub fn energy_ev(&self) -> f64 {
        self.energy_ev
    }

    /// What this potential includes.
    pub fn metadata(&self) -> &CorrectionMetadata {
        &self.meta
    }

    /// Exchange potential energy at `r` (0 when exchange is off), hartree.
    pub fn exchange_energy(&self, r: f64) -> f64 {
        if !self.exchange {
            return 0.0;
        }
        exchange_potential(self.e - self.stat.energy(r), self.dens.density(r))
    }

    /// Polarization potential `V_pol(r)`, Eq. (3), without the join or the
    /// outer cut (0 when the correction is off), hartree.
    pub fn polarization_energy(&self, r: f64) -> f64 {
        match self.cp {
            Some(s) => -s.alpha / (2.0 * (r * r + s.d2).powi(2)),
            None => 0.0,
        }
    }

    /// Correlation-polarization potential `V_cp(r)`, Eq. (9), with the outer
    /// cut (0 when off), hartree.
    pub fn correlation_polarization_energy(&self, r: f64) -> f64 {
        let Some(s) = self.cp else {
            return 0.0;
        };
        if r >= s.r_out {
            return 0.0;
        }
        let vpol = -s.alpha / (2.0 * (r * r + s.d2).powi(2));
        if r < s.r_cp {
            vpol.max(correlation_potential(self.dens.density(r)))
        } else {
            vpol
        }
    }

    /// A solver at the energy this potential was built for.
    pub fn solver(&self, opts: SolverOptions) -> Result<ElasticSolver<'_>, ElasticError> {
        ElasticSolver::new(self, self.energy_ev, opts)
    }
}

impl ScreenedPotential for CorrectedPotential<'_> {
    fn energy(&self, r: f64) -> f64 {
        let v = self.stat.energy(r);
        if !self.exchange && self.cp.is_none() {
            return v;
        }
        let mut total = v;
        if self.exchange {
            total += exchange_potential(self.e - v, self.dens.density(r));
        }
        if self.cp.is_some() {
            total += self.correlation_polarization_energy(r);
        }
        total
    }
    fn nuclear_charge(&self) -> f64 {
        self.stat.nuclear_charge()
    }
    fn breakpoints(&self) -> Vec<f64> {
        let mut b = self.stat.breakpoints();
        if let Some(s) = self.cp {
            b.push(s.r_out);
        }
        b
    }
    fn length_scale(&self) -> f64 {
        self.stat.length_scale()
    }
    fn long_range(&self) -> bool {
        self.cp.is_some() || self.stat.long_range()
    }
    fn bound_energy_ev(&self) -> Option<f64> {
        Some(self.energy_ev)
    }
    fn correction_metadata(&self) -> CorrectionMetadata {
        self.meta.clone()
    }
    fn matching_radius(&self, threshold: f64) -> f64 {
        if !self.exchange && self.cp.is_none() {
            return self.stat.matching_radius(threshold);
        }
        // The default rule on V_st + V_ex (both <= 0), beyond every breakpoint
        // including the outer radius of the polarization tail.
        let mut r = self.length_scale();
        let floor = self.breakpoints().into_iter().fold(0.0_f64, f64::max);
        let short = |r: f64| self.stat.energy(r) + self.exchange_energy(r);
        while (r < floor || (r * short(r)).abs() >= threshold) && r < 1.0e7 {
            r *= 1.05;
        }
        r
    }
}

/// [`solve`] with optional corrections. With both off this **is** [`solve`]
/// on `pot` (bit-identical results, metadata reporting both off); otherwise
/// the corrected potential is built at `energy_ev` and solved at the same
/// energy, and the result's [`ElasticResult::corrections`] records the
/// choices.
pub fn solve_corrected(
    pot: &dyn ScreenedPotential,
    density: &dyn ElectronDensity,
    energy_ev: f64,
    corrections: &Corrections,
    thetas: &[f64],
    opts: SolverOptions,
) -> Result<ElasticResult, ElasticError> {
    if corrections.is_none() {
        return solve(pot, energy_ev, thetas, opts);
    }
    let corrected = CorrectedPotential::new(pot, density, energy_ev, corrections)?;
    solve(&corrected, energy_ev, thetas, opts)
}
