//! Elastic scattering of electrons by a static central potential: the radial
//! Dirac equation, partial-wave phase shifts and the Mott cross sections.
//!
//! This is step 1 of the Mott elastic model (issue #17): phase shifts, the
//! differential cross section `DCS(theta)`, the total elastic cross section
//! `sigma_el`, the first transport cross section `sigma_tr1` and the Sherman
//! function, all at **one** kinetic energy and for a potential given as a
//! function. Optional Furness-McCarthy exchange and correlation-polarization
//! corrections (issue #91), off by default, are in [`corrections`]
//! ([`solve_corrected`]). Energy-grid tables for a material are in [`table`]
//! (issue #90); condensed phases are a follow-up issue.
//!
//! # Units
//!
//! Hartree atomic units inside this module: lengths in bohr, energies in
//! hartree, the electron mass is 1 and the speed of light is `c = 1/alpha`
//! (`alpha` from [`crate::constants::FINE_STRUCTURE`]). Cross sections come out
//! in bohr^2 (`DCS` in bohr^2/sr); [`BOHR2_TO_M2`] converts. The public energy
//! argument is the electron **kinetic** energy in eV.
//!
//! # Method
//!
//! With `P(r) = r g(r)` and `Q(r) = r f(r)` the large and small radial
//! components, the radial Dirac equations for a central potential energy
//! `V(r)` and Dirac quantum number `kappa` are (kinetic energy `E`, `c = 1/alpha`)
//!
//! ```text
//! dP/dr = -kappa/r P + (E - V + 2c^2)/c Q
//! dQ/dr = -(E - V)/c P + kappa/r Q
//! ```
//!
//! (the form used in the ELSEPA paper, Salvat, Jablonski & Powell, Comput.
//! Phys. Commun. 165, 157 (2005); here re-derived from the Dirac equation and
//! checked against the closed-form free solution below, not copied from any
//! code). `kappa = -(l+1)` is the spin-up channel `j = l + 1/2` and `kappa = l`
//! the spin-down channel `j = l - 1/2` (`l >= 1`).
//!
//! * **Free solutions.** For `V = 0` the regular and irregular solutions are
//!   `r (u_l(kr), sigma s u_lbar(kr))`, with `u` a spherical Bessel `j` or
//!   `y`, `k = sqrt(E (E + 2c^2))/c`, `s = c k/(E + 2c^2)`, `sigma = -1` and
//!   `lbar = l + 1` for `kappa < 0`, `sigma = +1` and `lbar = l - 1` for
//!   `kappa > 0`. (Verified by substitution in the equations above and by the
//!   tests.)
//! * **Integration.** The regular solution is integrated outward from a small
//!   radius `r_s`. Inside the centrifugal barrier the irregular solution
//!   decays outward, so an imprecise starting ratio `Q/P` is harmless as long
//!   as `r_s` is well below the classical turning point: `r_s` is placed by a
//!   WKB criterion on the bisected turning point, so that the irregular
//!   admixture is suppressed by `exp(-60)` there (see `wkb_start`), but
//!   never below a floor of `1e-10/(1+Z)` bohr, where the leading-order
//!   point-Coulomb or regular power-series start is accurate. (Until #131 the
//!   static potentials used `r_s = r_t exp(-60/|kappa|)`, which assumes
//!   power-law growth; measured against the WKB start for Au, Salvat DHFS
//!   screening, its phase-shift error was below `2e-12` rad at 10 and 20 keV,
//!   `6.3e-10` rad at 30 keV, `1.6e-7` rad at 50 keV and `2.5e-4` rad at
//!   100 keV (`|kappa|` 779, `sigma_tr1` off by `2.4e-5` relative).) The
//!   integrator
//!   is the Gragg-Bulirsch-Stoer modified-midpoint method with polynomial
//!   extrapolation in `h^2` (Stoer & Bulirsch, *Introduction to Numerical
//!   Analysis*), substep sequence 2, 4, 6, 8 (order 8), with a step set from
//!   the local eigenvalue of the system,
//!   `h = step_scale / (sqrt|q| + (1 + (2 kappa^2)^(1/3))/r)`
//!   (the cube-root term resolves the Airy region at a turning point), capped
//!   at a quarter of `step_scale` times the potential's length scale,
//!   `q = (E-V)(E-V+2c^2)/c^2 - kappa^2/r^2`. The amplitude is rescaled when
//!   it grows past `1e100`; only the ratio `Q/P` matters.
//! * **Matching.** At the matching radius `R`, where `|r V(r)|` has fallen
//!   below [`SolverOptions::matching_threshold`] (a potential that has
//!   a finite range ends there exactly), `tan(delta) = (Q J_p - P J_q) /
//!   (Q Y_p - P Y_q)` with `J`, `Y` the regular and irregular free solutions.
//!   The potentials here are screened, so no Coulomb tail is matched. The
//!   neglected tail changes `delta` by roughly the threshold divided by `k`.
//! * **Partial-wave sum.** `f(theta)` and `g(theta)` are the Mott direct and
//!   spin-flip amplitudes, `f = (1/2ik) sum_l [(l+1)(e^{2 i d_-} - 1) +
//!   l (e^{2 i d_+} - 1)] P_l(cos theta)`, `g = (1/2ik) sum_l
//!   [e^{2 i d_-} - e^{2 i d_+}] P_l^1(cos theta)` with `P_l^1 = sin theta
//!   dP_l/d cos theta`, and `DCS = |f|^2 + |g|^2`. These are the standard
//!   textbook (Mott 1929; Motz, Olsen & Koch 1964) expressions written from
//!   memory of the form, **not** checked digit for digit against a paper in
//!   this work; the tests that pin them are the optical theorem, the exact
//!   `1 - beta^2 sin^2(theta/2)` first-Born limit (which fixes the relative
//!   size of `g`) and the square-well closed form. The Sherman function is
//!   `S = i (f g* - f* g) / (|f|^2 + |g|^2)`; its **sign convention is not
//!   verified** against a published definition.
//! * **Cutoff in `l`.** Phase shifts are computed for `l = 0, 1, 2, ...` and
//!   the series is cut after `consecutive` (5) successive `l` for which both
//!   `|delta_-|` and `|delta_+|` are below `phase_tolerance` (1e-8 rad);
//!   [`PartialWaves::l_max`] is the highest `l` computed. When `x = kR` is
//!   below `l` by enough that the free Bessel functions leave the range of
//!   `f64`, the wave never reaches the potential and its shift is set to
//!   zero. At high energy `l_max` grows into the thousands; this module
//!   **integrates every partial wave directly** (cost linear in `l_max` times
//!   `k R`). A WKB / asymptotic phase for large `l`, as in the literature on
//!   which the issue draws (Salvat et al. 2005), is **not implemented**: that
//!   paper was not opened for this work (closed access) and the method is
//!   not reproduced from memory. See follow-up work.
//! * **Cross sections.** `sigma_el = (4 pi/k^2) sum_l [(l+1) sin^2 d_- +
//!   l sin^2 d_+]` (equal to `(4 pi/k) Im f(0)`, the optical theorem). The
//!   `DCS` times `(1 - cos theta)` is a polynomial of degree at most
//!   `2 l_max + 1` in `cos theta`, so `sigma_tr1` is integrated *exactly*
//!   with an `(l_max + 2)`-point Gauss-Legendre rule.
//!
//! # Potentials
//!
//! [`Yukawa`] and [`SquareWell`] are analytically tractable. [`SalvatDhfs`] is
//! the analytic Dirac-Hartree-Fock-Slater screening function of Salvat,
//! Martinez, Mayol & Parellada, Phys. Rev. A 36, 467 (1987),
//! doi:10.1103/PhysRevA.36.467, with its Table I coefficients for Z = 1..92
//! (transcribed from the open repository copy of the paper, #130), see the
//! type and `docs/data-provenance.md`.
//!
//! # Determinism
//!
//! Each phase shift is an independent, deterministic computation. Blocks of
//! partial waves are evaluated in parallel with rayon and collected in order,
//! so the result is bit-identical on any thread count (tested).

use crate::constants::{BOHR_RADIUS, ELEMENTARY_CHARGE, FINE_STRUCTURE, HARTREE_ENERGY};
use rayon::prelude::*;

pub mod corrections;
pub use corrections::{
    solve_corrected, CorrectedPotential, CorrectionMetadata, Corrections, CorrelationPolarization,
    ElectronDensity, PolarizationCutoff,
};

/// Square metres per bohr squared (`a0^2`, [`BOHR_RADIUS`] squared).
pub const BOHR2_TO_M2: f64 = BOHR_RADIUS * BOHR_RADIUS;

/// Hartree energy in eV (CODATA 2022 `E_h` over `e`).
pub const HARTREE_EV: f64 = HARTREE_ENERGY / ELEMENTARY_CHARGE;

/// Errors from the elastic solver.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ElasticError {
    /// An input violates an invariant.
    #[error("invalid {what}: {reason}")]
    Invalid {
        /// Which quantity.
        what: &'static str,
        /// What is wrong with it.
        reason: String,
    },
    /// No screening coefficients for this atomic number: Table I of Salvat
    /// et al. (1987) covers Z = 1..92 only.
    #[error(
        "Salvat et al. (1987) DHFS screening coefficients for Z={0} are not available: \
         Table I covers Z = 1..92"
    )]
    ScreeningCoefficientsUnavailable(u32),
    /// The integration produced a non-finite value.
    #[error("numerical failure integrating the radial equations for kappa = {0}")]
    NumericalFailure(i32),
    /// The adaptive partial-wave sum did not converge by `max_l`.
    #[error("partial-wave series not converged by l = {0}")]
    NotConverged(usize),
}

fn invalid(what: &'static str, reason: impl Into<String>) -> ElasticError {
    ElasticError::Invalid {
        what,
        reason: reason.into(),
    }
}

/// A static, spherically symmetric potential energy of the electron.
pub trait ScreenedPotential: Sync {
    /// Potential energy `V(r)` of the electron, hartree, at radius `r` bohr.
    /// Attractive potentials are negative.
    fn energy(&self, r: f64) -> f64;

    /// Coefficient `Z` of the point-Coulomb singularity `V -> -Z/r` at the
    /// origin (0 for a regular potential). Selects the series used to start
    /// the integration.
    fn nuclear_charge(&self) -> f64 {
        0.0
    }

    /// Radii at which `V` is discontinuous; the integrator restarts there.
    fn breakpoints(&self) -> Vec<f64> {
        Vec::new()
    }

    /// A length over which `V` changes appreciably, bohr. It caps the step.
    fn length_scale(&self) -> f64;

    /// The kinetic energy (eV) this potential was built for, when it depends
    /// on one (such as [`CorrectedPotential`], whose exchange and polarization
    /// cutoff use the energy). [`ElasticSolver::new`] rejects any other solver
    /// energy. The default, `None`, accepts every energy.
    fn bound_energy_ev(&self) -> Option<f64> {
        None
    }

    /// Which corrections the potential includes; copied into
    /// [`ElasticResult::corrections`] by [`solve`]. The default reports both
    /// off.
    fn correction_metadata(&self) -> CorrectionMetadata {
        CorrectionMetadata::default()
    }

    /// Radius beyond which the potential is neglected: the smallest radius on
    /// a geometric grid (ratio 1.05, from `length_scale`) with
    /// `|r V(r)| < threshold` that also lies beyond every breakpoint.
    fn matching_radius(&self, threshold: f64) -> f64 {
        let mut r = self.length_scale();
        let floor = self.breakpoints().into_iter().fold(0.0_f64, f64::max);
        while (r < floor || (r * self.energy(r)).abs() >= threshold) && r < 1.0e7 {
            r *= 1.05;
        }
        r
    }
}

/// Exponentially screened Coulomb potential `V(r) = -Z exp(-r/a)/r`
/// (hartree, bohr). Its first-Born cross sections are analytic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Yukawa {
    z: f64,
    screening_length: f64,
}

impl Yukawa {
    /// `z` is the (positive for attraction) charge, `screening_length` the
    /// screening length `a` in bohr.
    pub fn new(z: f64, screening_length: f64) -> Result<Self, ElasticError> {
        if !z.is_finite() || z < 0.0 {
            return Err(invalid("Yukawa charge", format!("{z} (need finite, >= 0)")));
        }
        if !screening_length.is_finite() || screening_length <= 0.0 {
            return Err(invalid(
                "Yukawa screening length",
                format!("{screening_length} (need finite, > 0)"),
            ));
        }
        Ok(Self {
            z,
            screening_length,
        })
    }

    /// Charge `Z`.
    pub fn z(&self) -> f64 {
        self.z
    }

    /// Screening length, bohr.
    pub fn screening_length(&self) -> f64 {
        self.screening_length
    }
}

impl ScreenedPotential for Yukawa {
    fn energy(&self, r: f64) -> f64 {
        -self.z * (-r / self.screening_length).exp() / r
    }
    fn nuclear_charge(&self) -> f64 {
        self.z
    }
    fn length_scale(&self) -> f64 {
        self.screening_length
    }
}

/// Spherical square well `V = -depth` for `r < radius`, 0 beyond (hartree,
/// bohr). A negative `depth` is a barrier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SquareWell {
    depth: f64,
    radius: f64,
}

impl SquareWell {
    /// Build a well; `radius` must be positive and `depth` finite.
    pub fn new(depth: f64, radius: f64) -> Result<Self, ElasticError> {
        if !depth.is_finite() {
            return Err(invalid("square-well depth", format!("{depth}")));
        }
        if !radius.is_finite() || radius <= 0.0 {
            return Err(invalid(
                "square-well radius",
                format!("{radius} (need finite, > 0)"),
            ));
        }
        Ok(Self { depth, radius })
    }

    /// Depth `V0` (the potential energy inside is `-V0`), hartree.
    pub fn depth(&self) -> f64 {
        self.depth
    }

    /// Well radius, bohr.
    pub fn radius(&self) -> f64 {
        self.radius
    }
}

impl ScreenedPotential for SquareWell {
    fn energy(&self, r: f64) -> f64 {
        if r < self.radius {
            -self.depth
        } else {
            0.0
        }
    }
    fn breakpoints(&self) -> Vec<f64> {
        vec![self.radius]
    }
    fn length_scale(&self) -> f64 {
        self.radius
    }
    fn matching_radius(&self, _threshold: f64) -> f64 {
        self.radius
    }
}

/// The three-term analytic screening function of Salvat, Martinez, Mayol &
/// Parellada, Phys. Rev. A 36, 467 (1987), doi:10.1103/PhysRevA.36.467:
/// `V(r) = -(Z/r) sum_i A_i exp(-alpha_i r)` with `sum_i A_i = 1`
/// (neutral atom), `r` in bohr.
///
/// The form is the paper's Eq. (11) (screening function) inserted in its
/// Eq. (1); its Eq. (12) gives the electron density used by
/// [`ElectronDensity`]. Both were checked (issue #91) in the open copy of the
/// paper in the University of Barcelona repository (diposit.ub.edu; the
/// publisher's copy is closed). [`SalvatDhfs::for_element`] returns the
/// paper's Table I coefficients for Z = 1..92 (pp. 470-471, transcribed in
/// #130, see `salvat_table.rs` and `docs/data-provenance.md`);
/// [`SalvatDhfs::from_coefficients`] takes caller-supplied ones.
///
/// Rows marked with an asterisk in Table I (H..P, Ca, Sc, Se, Br, Kr, Xe)
/// have `A_3 = 0` and no `alpha_3` (p. 471); they are held as **two-term**
/// potentials, so no `alpha_3` value enters [`ScreenedPotential::energy`],
/// [`ElectronDensity::density`] or [`ScreenedPotential::length_scale`].
///
/// Some amplitudes are negative (H, He: `A_1`; S, Cl, Ar: `A_2`). For H
/// (`A_1 = -184.39`, `A_2 = 185.39`, nearly equal `alpha`s) the two terms
/// cancel to about 1 part in 370 near the nucleus, so `phi(r)` there keeps
/// roughly 13 significant digits in `f64` rather than 16.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SalvatDhfs {
    z: f64,
    a: [f64; 3],
    alpha: [f64; 3],
    /// Number of terms in use (2 or 3); entries past it are unused.
    terms: usize,
}

impl SalvatDhfs {
    /// Build from caller-supplied coefficients (`a` dimensionless summing to
    /// 1 within 1e-6, `alpha` in 1/bohr and positive).
    pub fn from_coefficients(z: u32, a: [f64; 3], alpha: [f64; 3]) -> Result<Self, ElasticError> {
        if !(1..=92).contains(&z) {
            return Err(invalid("atomic number", format!("{z} (need 1..=92)")));
        }
        if a.iter().chain(alpha.iter()).any(|v| !v.is_finite()) {
            return Err(invalid("screening coefficients", "not finite"));
        }
        if (a.iter().sum::<f64>() - 1.0).abs() > 1e-6 {
            return Err(invalid(
                "screening coefficients",
                "the amplitudes A_i must sum to 1 (neutral atom)",
            ));
        }
        if alpha.iter().any(|&v| v <= 0.0) {
            return Err(invalid("screening coefficients", "alpha_i must be > 0"));
        }
        Ok(Self {
            z: f64::from(z),
            a,
            alpha,
            terms: 3,
        })
    }

    /// Build a two-term potential (`A_3 = 0`, the form of the asterisked rows
    /// of Salvat et al. (1987) Table I): `a` sums to 1 within 1e-6, `alpha`
    /// in 1/bohr and positive.
    pub fn from_two_terms(z: u32, a: [f64; 2], alpha: [f64; 2]) -> Result<Self, ElasticError> {
        // Validate through the three-term path with a zero third term; the
        // third alpha only satisfies the check and is never used.
        let mut p = Self::from_coefficients(z, [a[0], a[1], 0.0], [alpha[0], alpha[1], 1.0])?;
        p.alpha[2] = 0.0;
        p.terms = 2;
        Ok(p)
    }

    /// The Salvat et al. (1987) Table I potential of element `z`
    /// (pp. 470-471, `A_3 = 1 - A_1 - A_2`; two-term for the asterisked
    /// rows). [`ElasticError::ScreeningCoefficientsUnavailable`] outside
    /// Z = 1..92.
    pub fn for_element(z: u32) -> Result<Self, ElasticError> {
        let row = z
            .checked_sub(1)
            .and_then(|i| salvat_table::TABLE_I.get(i as usize))
            .ok_or(ElasticError::ScreeningCoefficientsUnavailable(z))?;
        debug_assert_eq!(row.z, z);
        match row.alpha3 {
            Some(alpha3) => Self::from_coefficients(
                z,
                [row.a1, row.a2, 1.0 - row.a1 - row.a2],
                [row.alpha1, row.alpha2, alpha3],
            ),
            None => Self::from_two_terms(z, [row.a1, row.a2], [row.alpha1, row.alpha2]),
        }
    }

    /// Amplitudes `A_i` of the terms in use (2 or 3).
    pub fn amplitudes(&self) -> &[f64] {
        &self.a[..self.terms]
    }

    /// Screening constants `alpha_i` (1/bohr) of the terms in use (2 or 3).
    pub fn alphas(&self) -> &[f64] {
        &self.alpha[..self.terms]
    }
}

impl ScreenedPotential for SalvatDhfs {
    fn energy(&self, r: f64) -> f64 {
        let phi: f64 = self
            .amplitudes()
            .iter()
            .zip(self.alphas())
            .map(|(a, al)| a * (-al * r).exp())
            .sum();
        -self.z * phi / r
    }
    fn nuclear_charge(&self) -> f64 {
        self.z
    }
    fn length_scale(&self) -> f64 {
        1.0 / self.alphas().iter().copied().fold(0.0, f64::max)
    }
}

/// Numerical settings of the solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverOptions {
    /// Matching radius criterion: `|r V(r)| < threshold` (hartree bohr).
    pub matching_threshold: f64,
    /// A partial wave is "small" when `|delta| < phase_tolerance` (rad).
    pub phase_tolerance: f64,
    /// Number of consecutive small partial waves that ends the series.
    pub consecutive: usize,
    /// Radians of local phase per integrator step (smaller is more accurate).
    pub step_scale: f64,
    /// Safety cap on `l`.
    pub max_l: usize,
    /// Speed of light, atomic units. Defaults to `1/alpha`; the
    /// non-relativistic limit is `c -> infinity`.
    pub light_speed: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            matching_threshold: 1.0e-10,
            phase_tolerance: 1.0e-8,
            consecutive: 5,
            step_scale: 0.25,
            max_l: 20_000,
            light_speed: 1.0 / FINE_STRUCTURE,
        }
    }
}

// ---------------------------------------------------------------------------
// Spherical Bessel functions

/// Spherical Bessel `j_n(x)` and `y_n(x)` for `n = 0..=n_max`. `None` when a
/// value leaves a safe range of `f64` (wave deep in the forbidden region).
fn spherical_bessel(n_max: usize, x: f64) -> Option<(Vec<f64>, Vec<f64>)> {
    if x.is_nan() || x <= 0.0 || !x.is_finite() {
        return None;
    }
    let (s, c) = x.sin_cos();
    let mut y = vec![0.0; n_max + 2];
    y[0] = -c / x;
    y[1] = -c / (x * x) - s / x;
    for n in 1..=n_max {
        y[n + 1] = (2 * n + 1) as f64 / x * y[n] - y[n - 1];
        if !y[n + 1].is_finite() || y[n + 1].abs() > 1.0e280 {
            return None;
        }
    }
    y.truncate(n_max + 1);

    // j_0, j_1 anchors (series for small x, where the closed forms cancel).
    let (j0, j1) = if x < 0.1 {
        let x2 = x * x;
        (
            1.0 - x2 / 6.0 * (1.0 - x2 / 20.0 * (1.0 - x2 / 42.0 * (1.0 - x2 / 72.0))),
            x / 3.0 * (1.0 - x2 / 10.0 * (1.0 - x2 / 28.0 * (1.0 - x2 / 54.0))),
        )
    } else {
        (s / x, s / (x * x) - c / x)
    };

    // Miller downward recurrence from far above both n_max and x.
    let start = n_max + 30 + (2.0 * x) as usize + (40.0 * (n_max as f64)).sqrt() as usize;
    let mut j = vec![0.0; start + 2];
    j[start] = 1.0e-30;
    for n in (1..=start).rev() {
        j[n - 1] = (2 * n + 1) as f64 / x * j[n] - j[n + 1];
        if j[n - 1].abs() > 1.0e250 {
            for v in &mut j[(n - 1)..] {
                *v *= 1.0e-250;
            }
        }
    }
    let norm = if j0.abs() >= j1.abs() {
        j0 / j[0]
    } else {
        j1 / j[1]
    };
    let jn: Vec<f64> = j[..=n_max].iter().map(|v| v * norm).collect();
    if jn.iter().any(|v| !v.is_finite()) {
        return None;
    }
    Some((jn, y))
}

// ---------------------------------------------------------------------------
// Gauss-Legendre nodes

/// Gauss-Legendre nodes and weights on `[-1, 1]` (Newton iteration on the
/// three-term recurrence), nodes ascending.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    let m = n.div_ceil(2);
    for i in 0..m {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p1, mut p2) = (1.0, 0.0);
            for l in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * l + 1) as f64 * z * p2 - l as f64 * p3) / (l + 1) as f64;
            }
            dp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1.0e-15 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        let wi = 2.0 / ((1.0 - z * z) * dp * dp);
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
    (x, w)
}

// ---------------------------------------------------------------------------
// Radial integration

type Pair = [f64; 2];

/// One Gragg-Bulirsch-Stoer step of size `h` for `y' = f(r, y)`.
fn bs_step(f: &impl Fn(f64, Pair) -> Pair, r: f64, y: Pair, h: f64) -> Pair {
    const SEQ: [usize; 4] = [2, 4, 6, 8];
    let f0 = f(r, y);
    let mut tab = [[[0.0_f64; 2]; 4]; 4];
    for (i, &n) in SEQ.iter().enumerate() {
        let hs = h / n as f64;
        let mut zp = y;
        let mut zc = [y[0] + hs * f0[0], y[1] + hs * f0[1]];
        for m in 1..n {
            let fz = f(r + m as f64 * hs, zc);
            let zn = [zp[0] + 2.0 * hs * fz[0], zp[1] + 2.0 * hs * fz[1]];
            zp = zc;
            zc = zn;
        }
        let fz = f(r + h, zc);
        tab[i][0] = [
            0.5 * (zc[0] + zp[0] + hs * fz[0]),
            0.5 * (zc[1] + zp[1] + hs * fz[1]),
        ];
        for j in 1..=i {
            let ratio = (SEQ[i] as f64 / SEQ[i - j] as f64).powi(2);
            let (cur, prev) = (tab[i][j - 1], tab[i - 1][j - 1]);
            tab[i][j] = [
                cur[0] + (cur[0] - prev[0]) / (ratio - 1.0),
                cur[1] + (cur[1] - prev[1]) / (ratio - 1.0),
            ];
        }
    }
    tab[3][3]
}

/// WKB amplitude exponent required between the start radius and the turning
/// point: the irregular solution admixed by an imprecise start is suppressed
/// by `exp(-2 W)`, here `exp(-60)`.
const WKB_START_EXPONENT: f64 = 30.0;

/// Start radius of the outward integration, for every potential (#91, #131).
///
/// In a classically forbidden region (`q < 0`) the WKB approximation gives the
/// two independent solutions as `|q|^(-1/4) exp(+-\int sqrt(-q) dr)` (Bender &
/// Orszag, *Advanced Mathematical Methods for Scientists and Engineers*,
/// McGraw-Hill 1978, ch. 10), so the solution decaying outward (the irregular
/// one, admixed by an imprecise starting ratio `Q/P`) loses `exp(-2 W)`
/// relative to the regular one between the start and the turning point,
/// `W = \int_{r_s}^{r_t} sqrt(-q) dr`.
///
/// The step-1 rule `r_t exp(-60/|kappa|)` assumed the power-law growth of a
/// free wave deep inside the barrier; for `|kappa|` of several hundred it lands
/// within a few per cent of the turning point, and the 5 % grid on which `r_t`
/// is found can even put it past the turning point, so the start error is not
/// suppressed (measured: phase errors up to 0.25 rad for `l ~ 1250` with a
/// polarization tail, #91; `2.5e-4` rad at `|kappa| = 779` for static Au at
/// 100 keV, #131). Here the turning point is refined by bisection and the start
/// is moved inward until `W >= WKB_START_EXPONENT` (midpoint rule on a
/// geometric grid of ratio 0.99), or to `r_floor`.
fn wkb_start(q_of: &impl Fn(f64) -> f64, r_t: f64, r_floor: f64, r_match: f64) -> f64 {
    // Refine the first sign change of q below the grid value r_t.
    let mut hi = r_t;
    let mut lo = (r_t / 1.05).max(r_floor);
    if q_of(lo) < 0.0 && q_of(hi) > 0.0 {
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if q_of(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
    } else {
        lo = r_t;
    }
    let mut r = lo;
    let mut w = 0.0;
    while w < WKB_START_EXPONENT && r > r_floor {
        let next = (r * 0.99).max(r_floor);
        let q = q_of(0.5 * (r + next));
        if q < 0.0 {
            w += (-q).sqrt() * (r - next);
        }
        r = next;
    }
    r.min(0.5 * r_match)
}

/// Where [`ElasticSolver::phase_shift`] starts the outward integration. Only
/// [`StartRule::Wkb`] is used outside the tests; the others are kept so the
/// tests can measure it against them (#131).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StartRule {
    /// `wkb_start` (every potential since #131).
    Wkb,
    /// The step-1 rule `r_t exp(-60/|kappa|)`, used by potentials without a
    /// polarization tail before #131.
    #[cfg_attr(not(test), allow(dead_code))]
    Legacy,
    /// The floor radius `1e-10/(1+Z)` bohr: the most expensive start and the
    /// least sensitive to the starting ratio, a reference for the others.
    #[cfg_attr(not(test), allow(dead_code))]
    Floor,
}

/// Mott partial-wave data at one energy: phase shifts and the amplitudes built
/// from them.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialWaves {
    k: f64,
    /// `delta_-[l]` is the phase shift of `kappa = -(l+1)`.
    delta_minus: Vec<f64>,
    /// `delta_+[l]` is the phase shift of `kappa = l` (`delta_+[0] = 0`: no such channel).
    delta_plus: Vec<f64>,
}

impl PartialWaves {
    /// Build from phase shifts for `l = 0..=l_max`. `delta_plus[0]` is ignored
    /// (set to 0). The two vectors must have equal, nonzero length.
    pub fn from_phase_shifts(
        k: f64,
        delta_minus: Vec<f64>,
        mut delta_plus: Vec<f64>,
    ) -> Result<Self, ElasticError> {
        if !k.is_finite() || k <= 0.0 {
            return Err(invalid("wave number", format!("{k}")));
        }
        if delta_minus.is_empty() || delta_minus.len() != delta_plus.len() {
            return Err(invalid("phase-shift arrays", "need equal, nonzero lengths"));
        }
        if delta_minus
            .iter()
            .chain(&delta_plus)
            .any(|v| !v.is_finite())
        {
            return Err(invalid("phase shifts", "not finite"));
        }
        delta_plus[0] = 0.0;
        Ok(Self {
            k,
            delta_minus,
            delta_plus,
        })
    }

    /// Wave number `k`, 1/bohr.
    pub fn wave_number(&self) -> f64 {
        self.k
    }

    /// Highest `l` for which phase shifts are held.
    pub fn l_max(&self) -> usize {
        self.delta_minus.len() - 1
    }

    /// Phase shift `delta_kappa` (rad); `kappa` is `-(l+1)` or `l >= 1`.
    pub fn phase_shift(&self, kappa: i32) -> Option<f64> {
        if kappa < 0 {
            self.delta_minus.get((-kappa - 1) as usize).copied()
        } else if kappa >= 1 {
            self.delta_plus.get(kappa as usize).copied()
        } else {
            None
        }
    }

    /// Direct and spin-flip amplitudes `(f, g)` at `cos(theta) = x`, each as
    /// `[re, im]`, bohr.
    pub fn amplitudes(&self, x: f64) -> (Pair, Pair) {
        let sin_t = (1.0 - x * x).max(0.0).sqrt();
        let (mut p_prev, mut p) = (0.0_f64, 1.0_f64);
        let (mut q_prev, mut q) = (0.0_f64, 0.0_f64);
        let (mut f, mut g) = ([0.0; 2], [0.0; 2]);
        let pref = 1.0 / (2.0 * self.k);
        for l in 0..=self.l_max() {
            let lf = l as f64;
            let (dm, dp) = (self.delta_minus[l], self.delta_plus[l]);
            // (e^{2 i d} - 1) = 2 i sin d e^{i d}; divide by 2 i k:  sin d e^{i d}/k.
            let (sm, cm) = dm.sin_cos();
            let (sp, cp) = dp.sin_cos();
            // fcoef = [(l+1) sin dm e^{i dm} + l sin dp e^{i dp}]/k
            let fr = ((lf + 1.0) * sm * cm + lf * sp * cp) / self.k;
            let fi = ((lf + 1.0) * sm * sm + lf * sp * sp) / self.k;
            f[0] += fr * p;
            f[1] += fi * p;
            // gcoef = (e^{2 i dm} - e^{2 i dp})/(2 i k)
            let gr = pref * ((2.0 * dm).sin() - (2.0 * dp).sin());
            let gi = -pref * ((2.0 * dm).cos() - (2.0 * dp).cos());
            g[0] += gr * q;
            g[1] += gi * q;
            // advance Legendre recurrences to l + 1
            let p_next = ((2.0 * lf + 1.0) * x * p - lf * p_prev) / (lf + 1.0);
            p_prev = p;
            p = p_next;
            let q_next = if l == 0 {
                sin_t
            } else {
                ((2.0 * lf + 1.0) * x * q - (lf + 1.0) * q_prev) / lf
            };
            q_prev = q;
            q = q_next;
        }
        (f, g)
    }

    /// Differential cross section `|f|^2 + |g|^2` at polar angle `theta`
    /// (rad), bohr^2/sr.
    pub fn dcs(&self, theta: f64) -> f64 {
        let (f, g) = self.amplitudes(theta.cos());
        f[0] * f[0] + f[1] * f[1] + g[0] * g[0] + g[1] * g[1]
    }

    /// Sherman function `S = i (f g* - f* g)/(|f|^2 + |g|^2)` at `theta`.
    /// The sign convention is not verified against a published definition.
    pub fn sherman(&self, theta: f64) -> f64 {
        let (f, g) = self.amplitudes(theta.cos());
        let den = f[0] * f[0] + f[1] * f[1] + g[0] * g[0] + g[1] * g[1];
        if den == 0.0 {
            return 0.0;
        }
        // i (f g* - f* g) = -2 Im(f g*)
        -2.0 * (f[1] * g[0] - f[0] * g[1]) / den
    }

    /// Total elastic cross section from the phase shifts, bohr^2.
    pub fn sigma_el(&self) -> f64 {
        let mut s = 0.0;
        for l in 0..=self.l_max() {
            let lf = l as f64;
            s += (lf + 1.0) * self.delta_minus[l].sin().powi(2)
                + lf * self.delta_plus[l].sin().powi(2);
        }
        4.0 * std::f64::consts::PI * s / (self.k * self.k)
    }

    /// `(4 pi/k) Im f(0)`: the optical-theorem value of `sigma_el`, bohr^2.
    pub fn sigma_optical(&self) -> f64 {
        let (f, _) = self.amplitudes(1.0);
        4.0 * std::f64::consts::PI / self.k * f[1]
    }

    /// `2 pi \int DCS(x) w(x) dx` by an exact Gauss-Legendre rule: `w` must be
    /// a polynomial of degree <= 1.
    fn integrate_dcs(&self, w: impl Fn(f64) -> f64) -> f64 {
        let (xs, ws) = gauss_legendre(self.l_max() + 2);
        let mut s = 0.0;
        for (x, wt) in xs.iter().zip(&ws) {
            let (f, g) = self.amplitudes(*x);
            let d = f[0] * f[0] + f[1] * f[1] + g[0] * g[0] + g[1] * g[1];
            s += wt * d * w(*x);
        }
        2.0 * std::f64::consts::PI * s
    }

    /// `sigma_el` by integrating the DCS (a consistency check of
    /// [`Self::sigma_el`]), bohr^2.
    pub fn sigma_el_quadrature(&self) -> f64 {
        self.integrate_dcs(|_| 1.0)
    }

    /// First transport cross section `\int (1 - cos theta) dsigma`, bohr^2.
    pub fn sigma_tr1(&self) -> f64 {
        self.integrate_dcs(|x| 1.0 - x)
    }
}

/// Result of [`solve`]: everything the issue asks for at one energy.
#[derive(Debug, Clone, PartialEq)]
pub struct ElasticResult {
    /// Kinetic energy, eV.
    pub energy_ev: f64,
    /// Highest partial wave computed.
    pub l_max: usize,
    /// The phase shifts behind the numbers below.
    pub partial_waves: PartialWaves,
    /// `DCS` at each requested angle, bohr^2/sr.
    pub dcs: Vec<f64>,
    /// Sherman function at each requested angle.
    pub sherman: Vec<f64>,
    /// Total elastic cross section, bohr^2.
    pub sigma_el: f64,
    /// First transport cross section, bohr^2.
    pub sigma_tr1: f64,
    /// Which exchange / correlation-polarization corrections the potential
    /// included ([`solve`] reports both off).
    pub corrections: CorrectionMetadata,
}

/// Radial Dirac solver at one kinetic energy.
pub struct ElasticSolver<'a> {
    pot: &'a dyn ScreenedPotential,
    e: f64,
    k: f64,
    c: f64,
    opts: SolverOptions,
    r_match: f64,
    breaks: Vec<f64>,
    /// Start rule. Always [`StartRule::Wkb`] outside the tests, which switch
    /// it to measure the other rules against it (#131).
    start_rule: StartRule,
}

impl<'a> ElasticSolver<'a> {
    /// Prepare a solver for kinetic energy `energy_ev` (eV, > 0).
    pub fn new(
        pot: &'a dyn ScreenedPotential,
        energy_ev: f64,
        opts: SolverOptions,
    ) -> Result<Self, ElasticError> {
        if !energy_ev.is_finite() || energy_ev <= 0.0 {
            return Err(invalid(
                "energy",
                format!("{energy_ev} eV (need finite, > 0)"),
            ));
        }
        if let Some(bound) = pot.bound_energy_ev() {
            if bound != energy_ev {
                return Err(invalid(
                    "energy",
                    format!(
                        "{energy_ev} eV does not match the potential's bound energy {bound} eV"
                    ),
                ));
            }
        }
        if opts.step_scale.is_nan() || opts.step_scale <= 0.0 || opts.step_scale > 4.0 {
            return Err(invalid(
                "step_scale",
                format!("{} (need 0 < s <= 4)", opts.step_scale),
            ));
        }
        if !opts.light_speed.is_finite()
            || opts.light_speed <= 0.0
            || opts.matching_threshold.is_nan()
            || opts.matching_threshold <= 0.0
            || opts.phase_tolerance.is_nan()
            || opts.phase_tolerance <= 0.0
            || opts.consecutive == 0
        {
            return Err(invalid("solver options", "non-positive setting"));
        }
        let e = energy_ev / HARTREE_EV;
        let c = opts.light_speed;
        let k = (e * (e + 2.0 * c * c)).sqrt() / c;
        let r_match = pot.matching_radius(opts.matching_threshold);
        let mut breaks: Vec<f64> = pot
            .breakpoints()
            .into_iter()
            .filter(|b| *b > 0.0 && *b < r_match)
            .collect();
        breaks.sort_by(f64::total_cmp);
        Ok(Self {
            pot,
            e,
            k,
            c,
            opts,
            r_match,
            breaks,
            start_rule: StartRule::Wkb,
        })
    }

    /// Wave number, 1/bohr.
    pub fn wave_number(&self) -> f64 {
        self.k
    }

    /// Matching radius, bohr.
    pub fn matching_radius(&self) -> f64 {
        self.r_match
    }

    /// Phase shift (rad, reduced to `(-pi/2, pi/2]`) of Dirac channel `kappa`
    /// (`kappa <= -1` or `kappa >= 1`).
    pub fn phase_shift(&self, kappa: i32) -> Result<f64, ElasticError> {
        if kappa == 0 {
            return Err(invalid("kappa", "0 is not a Dirac quantum number"));
        }
        let (e, c) = (self.e, self.c);
        let kap = f64::from(kappa);
        let l = if kappa < 0 { -kappa - 1 } else { kappa } as usize;
        let lbar = if kappa < 0 { l + 1 } else { l - 1 };
        let sigma = if kappa < 0 { -1.0 } else { 1.0 };
        let s = c * self.k / (e + 2.0 * c * c);

        // Free Bessel functions at the matching radius.
        let Some((jb, yb)) = spherical_bessel(l.max(lbar), self.k * self.r_match) else {
            return Ok(0.0);
        };
        let (jp, jq) = (jb[l], sigma * s * jb[lbar]);
        let (yp, yq) = (yb[l], sigma * s * yb[lbar]);

        let z = self.pot.nuclear_charge();
        let pot = self.pot;
        let q_of = |r: f64| {
            let d = e - pot.energy(r);
            d * (d + 2.0 * c * c) / (c * c) - kap * kap / (r * r)
        };

        // Turning point: first radius where the local wave is oscillatory.
        let r_floor = 1.0e-10 / (1.0 + z);
        let mut r_t = self.r_match;
        let mut r = r_floor;
        while r < self.r_match {
            if q_of(r) > 0.0 {
                r_t = r;
                break;
            }
            r *= 1.05;
        }
        let r_s = match self.start_rule {
            StartRule::Wkb => wkb_start(&q_of, r_t, r_floor, self.r_match),
            StartRule::Legacy => (r_t * (-60.0 / kap.abs()).exp())
                .max(r_floor)
                .min(0.5 * self.r_match),
            StartRule::Floor => r_floor,
        };

        // Starting values.
        let y0: Pair = if z > 0.0 {
            let alpha = c.recip();
            let gamma = (kap * kap - (alpha * z).powi(2)).sqrt();
            if kappa < 0 {
                [1.0, -(z / c) / (gamma - kap)]
            } else {
                [z / (c * (gamma + kap)), 1.0]
            }
        } else {
            let v = pot.energy(r_s);
            if kappa < 0 {
                [1.0, -(e - v) / (c * (2.0 * l as f64 + 3.0)) * r_s]
            } else {
                [(e - v + 2.0 * c * c) / (c * (2.0 * kap + 1.0)) * r_s, 1.0]
            }
        };

        let [p, q] = self.integrate(kap, r_s, y0)?;
        let num = q * jp - p * jq;
        let den = q * yp - p * yq;
        if !(num.is_finite() && den.is_finite()) {
            return Err(ElasticError::NumericalFailure(kappa));
        }
        Ok((num / den).atan())
    }

    fn integrate(&self, kap: f64, r_s: f64, y0: Pair) -> Result<Pair, ElasticError> {
        let (e, c) = (self.e, self.c);
        let pot = self.pot;
        let scale = pot.length_scale();
        let mut y = y0;
        let mut r;
        let mut edges: Vec<f64> = self.breaks.iter().copied().filter(|b| *b > r_s).collect();
        edges.push(self.r_match);
        let mut a = r_s;
        for &b in &edges {
            let lo = a + 1.0e-12 * b;
            let hi = b - 1.0e-12 * b;
            let f = |rr: f64, yy: Pair| -> Pair {
                let v = pot.energy(rr.clamp(lo, hi.max(lo)));
                let (d1, d2) = ((e - v + 2.0 * c * c) / c, (e - v) / c);
                [
                    -kap / rr * yy[0] + d1 * yy[1],
                    -d2 * yy[0] + kap / rr * yy[1],
                ]
            };
            r = a;
            while r < b {
                let v = pot.energy(r.clamp(lo, hi.max(lo)));
                let d = e - v;
                let q = d * (d + 2.0 * c * c) / (c * c) - kap * kap / (r * r);
                // the cube-root term resolves the Airy-like region at a turning point
                let rate = q.abs().sqrt() + (1.0 + (2.0 * kap * kap).cbrt()) / r;
                let mut h = (self.opts.step_scale / rate).min(0.25 * self.opts.step_scale * scale);
                if r + 1.25 * h >= b {
                    h = b - r;
                }
                y = bs_step(&f, r, y, h);
                r = if h == b - r { b } else { r + h };
                let m = y[0].abs().max(y[1].abs());
                if !m.is_finite() {
                    return Err(ElasticError::NumericalFailure(kap as i32));
                }
                if m > 1.0e100 {
                    y = [y[0] * 1.0e-100, y[1] * 1.0e-100];
                }
            }
            a = b;
        }
        Ok(y)
    }

    /// Phase shifts for `l = 0..=l_max` (both spin channels), in parallel and
    /// deterministically.
    pub fn phase_shifts_up_to(&self, l_max: usize) -> Result<PartialWaves, ElasticError> {
        let ls: Vec<usize> = (0..=l_max).collect();
        let pairs: Result<Vec<(f64, f64)>, ElasticError> = ls
            .par_iter()
            .map(|&l| {
                let dm = self.phase_shift(-(l as i32) - 1)?;
                let dp = if l == 0 {
                    0.0
                } else {
                    self.phase_shift(l as i32)?
                };
                Ok((dm, dp))
            })
            .collect();
        let pairs = pairs?;
        PartialWaves::from_phase_shifts(
            self.k,
            pairs.iter().map(|p| p.0).collect(),
            pairs.iter().map(|p| p.1).collect(),
        )
    }

    /// Phase shifts with the adaptive cutoff of the module documentation.
    pub fn partial_waves(&self) -> Result<PartialWaves, ElasticError> {
        const BLOCK: usize = 32;
        let tol = self.opts.phase_tolerance;
        let (mut dm, mut dp) = (Vec::new(), Vec::new());
        let mut run = 0usize;
        let mut start = 0usize;
        while start <= self.opts.max_l {
            let end = (start + BLOCK).min(self.opts.max_l + 1);
            let ls: Vec<usize> = (start..end).collect();
            let block: Result<Vec<(f64, f64)>, ElasticError> = ls
                .par_iter()
                .map(|&l| {
                    let m = self.phase_shift(-(l as i32) - 1)?;
                    let p = if l == 0 {
                        0.0
                    } else {
                        self.phase_shift(l as i32)?
                    };
                    Ok((m, p))
                })
                .collect();
            for (m, p) in block? {
                dm.push(m);
                dp.push(p);
                if m.abs() < tol && p.abs() < tol {
                    run += 1;
                } else {
                    run = 0;
                }
                if run >= self.opts.consecutive {
                    return PartialWaves::from_phase_shifts(self.k, dm, dp);
                }
            }
            start = end;
        }
        Err(ElasticError::NotConverged(self.opts.max_l))
    }
}

/// Solve at one energy and evaluate the observables on a caller-supplied
/// angle grid (`thetas`, rad).
pub fn solve(
    pot: &dyn ScreenedPotential,
    energy_ev: f64,
    thetas: &[f64],
    opts: SolverOptions,
) -> Result<ElasticResult, ElasticError> {
    if thetas
        .iter()
        .any(|t| !t.is_finite() || *t < 0.0 || *t > std::f64::consts::PI)
    {
        return Err(invalid("theta grid", "angles must lie in [0, pi]"));
    }
    let solver = ElasticSolver::new(pot, energy_ev, opts)?;
    let pw = solver.partial_waves()?;
    Ok(ElasticResult {
        energy_ev,
        l_max: pw.l_max(),
        dcs: thetas.iter().map(|t| pw.dcs(*t)).collect(),
        sherman: thetas.iter().map(|t| pw.sherman(*t)).collect(),
        sigma_el: pw.sigma_el(),
        sigma_tr1: pw.sigma_tr1(),
        partial_waves: pw,
        corrections: pot.correction_metadata(),
    })
}

mod salvat_table;
pub mod table;

#[cfg(test)]
mod tests;
