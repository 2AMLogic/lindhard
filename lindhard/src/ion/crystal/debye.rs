//! Debye model of thermal vibration: RMS displacement and its sampling.
//!
//! # Model
//!
//! A monatomic harmonic crystal with the Debye density of states, cut off at
//! the Debye frequency `omega_D = k_B Theta_D / hbar`, has a mean-square
//! displacement of an atom along one Cartesian axis
//!
//! ```text
//! u1^2 = (3 hbar^2 / (m k_B Theta_D)) * [ (T/Theta_D) D1(Theta_D/T) + 1/4 ]
//! D1(x) = (1/x) * integral_0^x  t / (exp(t) - 1)  dt
//! ```
//!
//! where `m` is the atomic mass and `D1` is the Debye function of order 1
//! (called `Phi` in the channeling literature). The `1/4` is the zero-point
//! term; the `D1` term is the thermal part.
//!
//! **Source.** The expression is the Debye-model evaluation of the harmonic
//! mean-square displacement, `2W(q) = q^2 <u_x^2>` with
//! `<u_x^2> = (hbar / (6 m)) * int g(w) coth(hbar w / 2 k_B T) / w dw`
//! (w is angular frequency; g(w) is the per-atom density of states normalised
//! to three modes, `int_0^{w_D} g dw = 3`),
//! as in N. W. Ashcroft and N. D. Mermin, *Solid State Physics* (Holt,
//! Rinehart and Winston, 1976), Appendix L. We could not open the printed
//! book; the equation was read in the form given in the "Mean squared
//! displacement" section of the Wikipedia article "Debye function"
//! (<https://en.wikipedia.org/wiki/Debye_function>, retrieved 2026-10-07),
//! which cites that appendix and writes
//! `2W(q) = (3/2) (hbar^2 q^2 / (m hbar omega_D)) [2 (k_B T/hbar omega_D)
//! D1(hbar omega_D / k_B T) + 1/2]`; dividing by `q^2` gives the form above.
//! The same quantity in the channeling literature is D. S. Gemmell, Rev. Mod.
//! Phys. 46, 129 (1974), doi:10.1103/RevModPhys.46.129; we could not open that
//! paper either, so its equation number is not cited here and its text was
//! not consulted. As an independent check of the transcription, the unit
//! tests integrate the `coth` form directly with the Debye density of states
//! `g(w) = 9 w^2 / w_D^3` (so `int_0^{w_D} g dw = 3`) and compare it with the closed form.
//!
//! **Limits.** `T -> 0`: `u1^2 -> 3 hbar^2 / (4 m k_B Theta_D)` (zero-point).
//! `T >> Theta_D`: `D1(x) -> 1 - x/4 + x^2/36 - ...`, so
//! `u1^2 -> 3 hbar^2 T / (m k_B Theta_D^2) = 3 k_B T / (m omega_D^2)` (the
//! classical equipartition value). The `-x/4` term of `D1` cancels the
//! zero-point `1/4` exactly, so the relative correction is `+x^2/36` with
//! `x = Theta_D/T`, i.e. about 0.03 % at `T = 10 Theta_D`.
//!
//! # One-dimensional versus three-dimensional amplitude
//!
//! [`ThermalVibration::rms_1d`] is `u1 = sqrt(u1^2)`, the RMS displacement
//! along *one* Cartesian axis (the quantity named `u_1` in the channeling
//! literature, which sets the Gaussian width of the transverse displacement
//! distribution). [`ThermalVibration::rms_3d`] is the RMS length of the full
//! displacement vector, `u = sqrt(3) u1` for an isotropic (cubic) crystal.
//! Sampling draws each Cartesian component independently from a Gaussian of
//! standard deviation `u1`.
//!
//! # Debye temperature
//!
//! A Debye temperature is not unique: it depends on the property it was
//! fitted to (low-temperature specific heat, elastic constants, or the
//! temperature dependence of X-ray or neutron Bragg intensities), and the
//! values differ by tens of kelvin for the same crystal. The constants
//! [`THETA_D_SI`] and friends are the room-temperature "Debye temperature"
//! entries of the Ioffe Institute *New Semiconductor Materials* handbook;
//! that source does not state the method for Si, Ge and GaAs, and for 3C-SiC
//! it gives Goldberg et al. (2001) as its reference (also without a method).
//! We therefore do **not** claim a specific-heat or X-ray determination for
//! them. They are the only openable, citable values we found, and they are
//! defaults only: the model takes `Theta_D` as an argument so a value
//! matched to a given measurement can be used instead. A lattice-vibration
//! amplitude fitted from X-ray Debye-Waller factors would be the more direct
//! choice for channeling work; no openable source for such values was found
//! (see `docs/data-provenance.md`, "Debye temperatures").
//!
//! For a compound (GaAs, SiC) one `Theta_D` is applied to every species with
//! that species' own mass, which is the usual single-parameter approximation
//! and not a property of the Debye model itself.
//!
//! # Sampling and determinism
//!
//! [`ThermalVibration::sample_displacement`] draws from the caller's
//! per-particle stream ([`crate::rng::ParticleRng`]) and always consumes
//! exactly six `u64` words (three Box-Muller pairs, of which only the cosine
//! member is used), so a history's stream position does not depend on the
//! sampled values and results do not depend on the thread count.
//! Box-Muller: G. E. P. Box and M. E. Muller, Ann. Math. Statist. 29, 610
//! (1958).

use crate::constants::{ATOMIC_MASS_UNIT, BOLTZMANN, HBAR};
use crate::rng::ParticleRng;
use rand_core::Rng;

/// Debye temperature of Si in kelvin: 640 K, "Debye temperature" row of
/// <https://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html> (retrieved
/// 2026-10-07). Determination method not stated by the source; see the module
/// docs.
pub const THETA_D_SI: f64 = 640.0;

/// Debye temperature of Ge in kelvin: 374 K, "Debye temperature" row of
/// <https://www.ioffe.ru/SVA/NSM/Semicond/Ge/basic.html> (retrieved
/// 2026-10-07). Method not stated by the source.
pub const THETA_D_GE: f64 = 374.0;

/// Debye temperature of GaAs in kelvin: 360 K, "Debye temperature" row of
/// <https://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html> (retrieved
/// 2026-10-07). Method not stated by the source.
pub const THETA_D_GAAS: f64 = 360.0;

/// Debye temperature of 3C-SiC in kelvin: 1200 K, "Debye temperature" row of
/// <https://www.ioffe.ru/SVA/NSM/Semicond/SiC/basic.html> (retrieved
/// 2026-10-07), which credits Goldberg et al. (2001). Method not stated.
pub const THETA_D_3C_SIC: f64 = 1200.0;

/// Errors from building a [`ThermalVibration`].
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum DebyeError {
    /// The Debye temperature must be finite and positive.
    #[error("Debye temperature must be finite and > 0 K, got {0}")]
    BadDebyeTemperature(f64),
    /// The atomic mass must be finite and positive.
    #[error("atomic mass must be finite and > 0 u, got {0}")]
    BadMass(f64),
    /// The temperature must be finite and non-negative.
    #[error("temperature must be finite and >= 0 K, got {0}")]
    BadTemperature(f64),
}

/// Debye function of order 1, `D1(x) = (1/x) int_0^x t/(e^t - 1) dt`, for
/// `x > 0`. `D1(0) = 1`, `D1(x) -> pi^2/(6x)` as `x -> infinity`.
///
/// For `x <= 3` the integrand's Taylor series `t/(e^t-1) = sum b_n t^n`
/// (radius of convergence `2 pi`) is integrated term by term; the `b_n` are
/// generated from `(t/(e^t-1)) * ((e^t-1)/t) = 1`. For larger `x` the integral
/// is `pi^2/6 - int_x^inf`, and the tail is the convergent series
/// `sum_k e^{-kx} (x/k + 1/k^2)`.
pub fn debye_d1(x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    if x <= 3.0 {
        const N: usize = 80;
        // a_k = 1/(k+1)!, so that sum a_k t^k = (e^t - 1)/t.
        let mut a = [0.0f64; N + 1];
        a[0] = 1.0;
        for k in 1..=N {
            a[k] = a[k - 1] / (k as f64 + 1.0);
        }
        let mut b = [0.0f64; N + 1];
        b[0] = 1.0;
        for n in 1..=N {
            let mut s = 0.0;
            for k in 1..=n {
                s += a[k] * b[n - k];
            }
            b[n] = -s;
        }
        // Horner in x for sum b_n x^n / (n+1).
        let mut acc = 0.0;
        for n in (0..=N).rev() {
            acc = acc * x + b[n] / (n as f64 + 1.0);
        }
        acc
    } else {
        let mut tail = 0.0;
        let mut k = 1.0f64;
        loop {
            let term = (-k * x).exp() * (x / k + 1.0 / (k * k));
            tail += term;
            if term < 1e-18 * tail || k > 2000.0 {
                break;
            }
            k += 1.0;
        }
        (std::f64::consts::PI.powi(2) / 6.0 - tail) / x
    }
}

/// Debye-model vibration of one atomic species at one temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalVibration {
    theta_d_k: f64,
    mass_amu: f64,
    temperature_k: f64,
    u1_sq: f64,
}

impl ThermalVibration {
    /// Build the model for Debye temperature `theta_d_k` (K), atomic mass
    /// `mass_amu` (g/mol, numerically the mass in unified atomic mass units)
    /// and target temperature `temperature_k` (K, may be 0).
    pub fn new(theta_d_k: f64, mass_amu: f64, temperature_k: f64) -> Result<Self, DebyeError> {
        if !(theta_d_k.is_finite() && theta_d_k > 0.0) {
            return Err(DebyeError::BadDebyeTemperature(theta_d_k));
        }
        if !(mass_amu.is_finite() && mass_amu > 0.0) {
            return Err(DebyeError::BadMass(mass_amu));
        }
        if !(temperature_k.is_finite() && temperature_k >= 0.0) {
            return Err(DebyeError::BadTemperature(temperature_k));
        }
        Ok(Self {
            theta_d_k,
            mass_amu,
            temperature_k,
            u1_sq: mean_square_1d(theta_d_k, mass_amu, temperature_k),
        })
    }

    /// Debye temperature in kelvin.
    pub fn debye_temperature(&self) -> f64 {
        self.theta_d_k
    }

    /// Target temperature in kelvin.
    pub fn temperature(&self) -> f64 {
        self.temperature_k
    }

    /// Atomic mass in u.
    pub fn mass_amu(&self) -> f64 {
        self.mass_amu
    }

    /// Mean-square displacement along one Cartesian axis, `u1^2`, in m^2.
    pub fn mean_square_1d(&self) -> f64 {
        self.u1_sq
    }

    /// One-dimensional RMS displacement `u1`, in metres.
    pub fn rms_1d(&self) -> f64 {
        self.u1_sq.sqrt()
    }

    /// Three-dimensional RMS displacement `u = sqrt(3) u1`, in metres.
    pub fn rms_3d(&self) -> f64 {
        (3.0 * self.u1_sq).sqrt()
    }

    /// Zero-point mean-square displacement `3 hbar^2 / (4 m k_B Theta_D)`
    /// (the `T -> 0` limit of [`Self::mean_square_1d`]), in m^2.
    pub fn zero_point_mean_square_1d(&self) -> f64 {
        mean_square_1d(self.theta_d_k, self.mass_amu, 0.0)
    }

    /// Classical (equipartition) mean-square displacement
    /// `3 hbar^2 T / (m k_B Theta_D^2)` (the `T >> Theta_D` limit), in m^2.
    pub fn classical_mean_square_1d(&self) -> f64 {
        let m = self.mass_amu * ATOMIC_MASS_UNIT;
        3.0 * HBAR * HBAR * self.temperature_k / (m * BOLTZMANN * self.theta_d_k * self.theta_d_k)
    }

    /// Draw an independent displacement of one atom from its lattice site:
    /// three Cartesian components, each Gaussian with standard deviation
    /// [`Self::rms_1d`], in metres. Consumes exactly six `u64` words from
    /// `rng`.
    pub fn sample_displacement(&self, rng: &mut ParticleRng) -> [f64; 3] {
        let s = self.rms_1d();
        let mut out = [0.0; 3];
        for c in &mut out {
            // u1 in (0, 1] so that ln is finite; u2 in [0, 1).
            let u1 = 1.0 - unit(rng);
            let u2 = unit(rng);
            *c = s * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
        }
        out
    }
}

/// `u1^2` in m^2 for Debye temperature (K), mass (u), temperature (K).
fn mean_square_1d(theta_d_k: f64, mass_amu: f64, temperature_k: f64) -> f64 {
    let m = mass_amu * ATOMIC_MASS_UNIT;
    let pref = 3.0 * HBAR * HBAR / (m * BOLTZMANN * theta_d_k);
    let thermal = if temperature_k > 0.0 {
        let x = theta_d_k / temperature_k;
        debye_d1(x) / x
    } else {
        0.0
    };
    pref * (thermal + 0.25)
}

#[inline]
fn unit(rng: &mut ParticleRng) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::stream;

    /// Composite Simpson reference, independent of the series and tail
    /// expansions used in `debye_d1`.
    fn simpson<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
        let h = (b - a) / n as f64;
        let mut s = f(a) + f(b);
        for i in 1..n {
            s += f(a + i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
        }
        s * h / 3.0
    }

    fn bose_integrand(t: f64) -> f64 {
        if t < 1e-8 {
            1.0 - 0.5 * t
        } else {
            t / t.exp_m1()
        }
    }

    fn d1_quadrature(x: f64) -> f64 {
        simpson(bose_integrand, 0.0, x, 20_000) / x
    }

    #[test]
    fn debye_function_matches_quadrature() {
        for &x in &[1e-6, 0.01, 0.3, 1.0, 2.9, 3.0, 3.1, 5.0, 10.0, 25.0, 40.0] {
            let want = d1_quadrature(x);
            let got = debye_d1(x);
            assert!((got - want).abs() <= 1e-12 * want, "x={x}: {got} vs {want}");
        }
        assert_eq!(debye_d1(0.0), 1.0);
    }

    #[test]
    fn debye_function_known_series_terms() {
        // D1(x) = 1 - x/4 + x^2/36 - x^4/3600 + ...
        let x = 0.1f64;
        let s = 1.0 - x / 4.0 + x * x / 36.0 - x.powi(4) / 3600.0;
        assert!((debye_d1(x) - s).abs() < 1e-9);
    }

    /// The formula as written in the module docs, evaluated independently.
    fn formula(theta: f64, t: f64, mass_amu: f64) -> f64 {
        let m = mass_amu * ATOMIC_MASS_UNIT;
        let x = theta / t;
        3.0 * HBAR * HBAR / (m * BOLTZMANN * theta) * (d1_quadrature(x) / x + 0.25)
    }

    #[test]
    fn matches_cited_formula() {
        for &(theta, t, mass) in &[
            (640.0, 300.0, 28.085),
            (374.0, 77.0, 72.63),
            (360.0, 600.0, 69.723),
            (1200.0, 1000.0, 12.011),
            (500.0, 5.0, 40.0),
        ] {
            let v = ThermalVibration::new(theta, mass, t).unwrap();
            let want = formula(theta, t, mass);
            assert!(
                (v.mean_square_1d() - want).abs() <= 1e-12 * want,
                "{theta} {t}: {} vs {want}",
                v.mean_square_1d()
            );
            assert!((v.rms_3d() - 3f64.sqrt() * want.sqrt()).abs() <= 1e-12 * want.sqrt() * 2.0);
        }
    }

    /// Direct integration of the harmonic coth form with the Debye density of
    /// states g(w) = 9 w^2 / w_D^3 (w angular frequency, normalised so int_0^{w_D} g dw = 3
    /// modes per atom): u1^2 = (hbar/(6 m)) int g(w) coth(hbar w/2kT) / w dw.
    #[test]
    fn matches_coth_integral_with_debye_dos() {
        let (theta, mass) = (640.0, 28.085);
        let m = mass * ATOMIC_MASS_UNIT;
        for &t in &[20.0, 100.0, 300.0, 900.0] {
            // Variable y = hbar w / (k T) in (0, x]; w = y kT/hbar, w_D = x kT/hbar.
            let x = theta / t;
            let kt = BOLTZMANN * t;
            // g normalised to 3 modes: int_0^{wD} 9 w^2/wD^3 dw = 3.
            // u1^2 = (hbar/(6 m)) * int g(w) coth(y/2)/w dw
            let integrand = |y: f64| {
                if y < 1e-9 {
                    // y coth(y/2) -> 2 as y -> 0.
                    return 9.0 * 2.0 / (x * x * x) / (kt / HBAR);
                }
                // g(w) dw / w = 9 w dw / wD^3 = 9 y dy / (x^3 kT/hbar), using:
                // w dw / wD^3 = y dy / (x^3 (kT/hbar)).
                9.0 * y / (x * x * x) / (kt / HBAR) * (y / 2.0).tanh().recip()
            };
            let val = HBAR / (2.0 * m) / 3.0 * simpson(integrand, 0.0, x, 200_000);
            let got = ThermalVibration::new(theta, mass, t)
                .unwrap()
                .mean_square_1d();
            assert!((got - val).abs() <= 1e-9 * val, "T={t}: {got} vs {val}");
        }
    }

    #[test]
    fn zero_temperature_is_zero_point() {
        for &(theta, mass) in &[(640.0, 28.085), (374.0, 72.63), (1200.0, 12.011)] {
            let v = ThermalVibration::new(theta, mass, 0.0).unwrap();
            let m = mass * ATOMIC_MASS_UNIT;
            let zp = 3.0 * HBAR * HBAR / (4.0 * m * BOLTZMANN * theta);
            assert!((v.mean_square_1d() - zp).abs() <= 1e-12 * zp);
            // And the T -> 0+ limit of the finite-T evaluation.
            let tiny = ThermalVibration::new(theta, mass, theta * 1e-8).unwrap();
            assert!((tiny.mean_square_1d() - zp).abs() <= 1e-12 * zp);
        }
    }

    #[test]
    fn classical_limit_at_ten_theta() {
        for &(theta, mass) in &[(640.0, 28.085), (374.0, 72.63), (360.0, 69.723)] {
            let v = ThermalVibration::new(theta, mass, 10.0 * theta).unwrap();
            let cl = v.classical_mean_square_1d();
            let rel = (v.mean_square_1d() - cl).abs() / cl;
            assert!(rel < 0.01, "{rel}");
        }
    }

    #[test]
    fn grows_with_temperature_and_rejects_bad_input() {
        let a = ThermalVibration::new(640.0, 28.085, 100.0).unwrap();
        let b = ThermalVibration::new(640.0, 28.085, 300.0).unwrap();
        assert!(b.rms_1d() > a.rms_1d());
        assert!(ThermalVibration::new(0.0, 28.0, 300.0).is_err());
        assert!(ThermalVibration::new(f64::NAN, 28.0, 300.0).is_err());
        assert!(ThermalVibration::new(640.0, -1.0, 300.0).is_err());
        assert!(ThermalVibration::new(640.0, 28.0, -1.0).is_err());
        assert!(ThermalVibration::new(640.0, 28.0, f64::INFINITY).is_err());
    }

    #[test]
    fn sampled_rms_matches_model() {
        let v = ThermalVibration::new(THETA_D_SI, 28.085, 300.0).unwrap();
        let mut rng = stream(7, 0);
        let n = 1_000_000;
        let (mut s1, mut s3) = (0.0, 0.0);
        for _ in 0..n {
            let d = v.sample_displacement(&mut rng);
            s1 += d[0] * d[0];
            s3 += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        }
        let rms1 = (s1 / n as f64).sqrt();
        let rms3 = (s3 / n as f64).sqrt();
        assert!((rms1 / v.rms_1d() - 1.0).abs() < 0.01, "{rms1}");
        assert!((rms3 / v.rms_3d() - 1.0).abs() < 0.01, "{rms3}");
    }

    #[test]
    fn sampling_consumes_fixed_stream_length() {
        let v = ThermalVibration::new(THETA_D_SI, 28.085, 300.0).unwrap();
        let mut a = stream(3, 9);
        let mut b = stream(3, 9);
        let _ = v.sample_displacement(&mut a);
        for _ in 0..6 {
            b.next_u64();
        }
        assert_eq!(a.next_u64(), b.next_u64());
    }
}
