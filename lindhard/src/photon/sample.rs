//! Samplers for the angular distributions of incoherent and coherent
//! scattering, and the rotation of a direction by a scattering angle.

use crate::constants::{ELECTRON_REST_ENERGY, ELEMENTARY_CHARGE, PLANCK, SPEED_OF_LIGHT};
use crate::rng::uniform;
use crate::vec3::{normalize, V3};
use rand_core::Rng;

use super::data::Curve;

/// `h c`, eV·Å (CODATA 2022 constants from [`crate::constants`]).
pub(crate) const HC_EV_ANGSTROM: f64 = PLANCK * SPEED_OF_LIGHT / ELEMENTARY_CHARGE * 1.0e10;

/// Electron rest energy, eV.
pub(crate) const MC2_EV: f64 = ELECTRON_REST_ENERGY / ELEMENTARY_CHARGE;

/// Sample the incoherent scattering of a photon of `e_ev` on an atom of atomic
/// number `z` with scattering function `s` (EPDL MF=27 MT=504): returns
/// `(cos θ, E')`.
///
/// The angular density is `dσ/dΩ ∝ KN(E, θ) S(x, Z)` with `x = sin(θ/2)/λ`
/// in 1/Å, the form of the EPDL incoherent cross section. The Klein-Nishina
/// factor is sampled by H. Kahn's rejection method (H. Kahn, "Applications of
/// Monte Carlo", RAND Corporation report RM-1237-AEC (1956)): with
/// `α = E/mc²` and `ξ = E/E'`, take `ξ = 1 + 2α r₁` and accept with
/// probability `4 (1/ξ - 1/ξ²)` (chosen with probability `(1 + 2α)/(9 + 2α)`),
/// or `ξ = (1 + 2α)/(1 + 2α r₁)` and accept with `(μ² + 1/ξ)/2`, `μ =
/// 1 - (ξ - 1)/α`. (The mixture reproduces `ε + 1/ε - 1 + μ²`, `ε = 1/ξ`, the
/// Klein-Nishina density in `ε`; checked by the unit tests against a
/// numerical integral.) The result is then accepted with probability
/// `S(x, Z)/Z` (`S ≤ Z`), which weights it by the binding correction.
pub(crate) fn sample_incoherent<R: Rng + ?Sized>(
    rng: &mut R,
    e_ev: f64,
    z: u8,
    s: &Curve,
) -> (f64, f64) {
    let alpha = e_ev / MC2_EV;
    let p_first = (1.0 + 2.0 * alpha) / (9.0 + 2.0 * alpha);
    let zf = f64::from(z);
    loop {
        let (r1, r2, r3) = (uniform(rng), uniform(rng), uniform(rng));
        let xi = if r1 <= p_first {
            let xi = 1.0 + 2.0 * alpha * r2;
            if r3 > 4.0 * (1.0 / xi - 1.0 / (xi * xi)) {
                continue;
            }
            xi
        } else {
            let xi = (1.0 + 2.0 * alpha) / (1.0 + 2.0 * alpha * r2);
            let mu = 1.0 - (xi - 1.0) / alpha;
            if r3 > 0.5 * (mu * mu + 1.0 / xi) {
                continue;
            }
            xi
        };
        let cos_t = (1.0 - (xi - 1.0) / alpha).clamp(-1.0, 1.0);
        let x = e_ev / HC_EV_ANGSTROM * ((1.0 - cos_t) * 0.5).sqrt();
        if uniform(rng) * zf <= s.eval(x.max(f64::MIN_POSITIVE)) {
            return (cos_t, e_ev / xi);
        }
    }
}

/// Inverse-CDF sampler of the coherent (Rayleigh) scattering angle of one
/// element.
///
/// The density is `dσ/dΩ ∝ (1 + μ²)/2 · F(x, Z)²` (Thomson times the form
/// factor; the anomalous scattering factors of EPDL MF=27 MT=505/506 are not
/// used), `μ = cos θ`, `x = x_m sin(θ/2)` with `x_m = E/hc`. Since `dμ ∝
/// d(x²)`, the angle is drawn by sampling `t = x²` from `F²` on `[0, x_m²]`,
/// then accepting with `(1 + μ²)/2` (acceptance at least 1/2). The integral
/// of `F²` over `t` is exact: between table points `F` is a power law in `x`
/// (the table's log-log interpolation), so `F²` is a power law in `t`, and
/// below the first point `F` is held constant.
#[derive(Debug, Clone)]
pub(crate) struct CoherentSampler {
    /// `t = x²` at the table points (increasing, `t[0] > 0`).
    t: Vec<f64>,
    /// `F²` at the table points.
    f2: Vec<f64>,
    /// Power-law exponent `s` of `F²(t)` on `[t[i], t[i+1]]`.
    s: Vec<f64>,
    /// `∫₀^{t[i]} F² dt`.
    cum: Vec<f64>,
}

/// `∫ F² dt` from `t0` to `t` for `F² = f0 (t/t0)^s`.
fn power_integral(f0: f64, t0: f64, s: f64, t: f64) -> f64 {
    let u = t / t0;
    if (s + 1.0).abs() < 1e-9 {
        f0 * t0 * u.ln()
    } else {
        f0 * t0 * (u.powf(s + 1.0) - 1.0) / (s + 1.0)
    }
}

impl CoherentSampler {
    /// Build from the form-factor curve `F(x)` (1/Å).
    pub(crate) fn new(form_factor: &Curve) -> Self {
        let x = form_factor.x();
        let y = form_factor.y();
        let n = x.len();
        let t: Vec<f64> = x.iter().map(|v| v * v).collect();
        let f2: Vec<f64> = y.iter().map(|v| v * v).collect();
        let mut s = Vec::with_capacity(n - 1);
        let mut cum = Vec::with_capacity(n);
        cum.push(f2[0] * t[0]); // F constant on [0, t0]
        for i in 0..n - 1 {
            let si = if t[i + 1] > t[i] {
                (f2[i + 1] / f2[i]).ln() / (t[i + 1] / t[i]).ln()
            } else {
                0.0
            };
            s.push(si);
            let seg = if t[i + 1] > t[i] {
                power_integral(f2[i], t[i], si, t[i + 1])
            } else {
                0.0
            };
            cum.push(cum[i] + seg);
        }
        Self { t, f2, s, cum }
    }

    /// `∫₀^t F² dt'` for `t` within the table.
    fn integral(&self, t: f64) -> f64 {
        let n = self.t.len();
        if t <= self.t[0] {
            return self.f2[0] * t;
        }
        let i = (self.t.partition_point(|&v| v <= t) - 1).min(n - 2);
        self.cum[i] + power_integral(self.f2[i], self.t[i], self.s[i], t)
    }

    /// The `t` with `∫₀^t F² = target` (inverse of [`Self::integral`]).
    fn invert(&self, target: f64) -> f64 {
        if target <= self.cum[0] {
            return target / self.f2[0];
        }
        let n = self.t.len();
        let i = (self.cum.partition_point(|&v| v <= target) - 1).min(n - 2);
        let rem = target - self.cum[i];
        let (f0, t0, s) = (self.f2[i], self.t[i], self.s[i]);
        if (s + 1.0).abs() < 1e-9 {
            t0 * (rem / (f0 * t0)).exp()
        } else {
            t0 * (1.0 + rem * (s + 1.0) / (f0 * t0))
                .max(0.0)
                .powf(1.0 / (s + 1.0))
        }
    }

    /// Sample `cos θ` for a photon of `e_ev`.
    pub(crate) fn sample_cos<R: Rng + ?Sized>(&self, rng: &mut R, e_ev: f64) -> f64 {
        let t_m = (e_ev / HC_EV_ANGSTROM).powi(2);
        // `E <= 1 MeV` keeps `t_m` inside the table (`x_m = 80.7 < 100`); for
        // safety, clamp to the last point.
        let t_m = t_m.min(self.t[self.t.len() - 1]);
        let total = self.integral(t_m);
        loop {
            let t = self.invert(uniform(rng) * total).min(t_m);
            let mu = 1.0 - 2.0 * t / t_m;
            if uniform(rng) <= 0.5 * (1.0 + mu * mu) {
                return mu.clamp(-1.0, 1.0);
            }
        }
    }
}

/// Rotate the unit vector `dir` by the polar angle with cosine `cos_t` and the
/// azimuth `phi` about it.
pub(crate) fn rotate(dir: V3, cos_t: f64, phi: f64) -> V3 {
    let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
    let (sp, cp) = phi.sin_cos();
    let [u, v, w] = dir;
    let out = if w.abs() > 0.999_999 {
        [sin_t * cp, sin_t * sp, cos_t * w.signum()]
    } else {
        let a = (1.0 - w * w).sqrt();
        [
            u * cos_t + sin_t * (u * w * cp - v * sp) / a,
            v * cos_t + sin_t * (v * w * cp + u * sp) / a,
            w * cos_t - sin_t * cp * a,
        ]
    };
    normalize(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::photon::PhotonTable;
    use crate::rng::stream;

    /// Mean `E'/E` of Klein-Nishina scattering by direct quadrature of
    /// `ε + 1/ε - 1 + μ²` over `μ` weighted by `ε²`.
    fn kn_mean_epsilon(alpha: f64) -> f64 {
        let n = 200_000;
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..n {
            let mu = -1.0 + 2.0 * (f64::from(i) + 0.5) / f64::from(n);
            let eps = 1.0 / (1.0 + alpha * (1.0 - mu));
            let w = eps * eps * (eps + 1.0 / eps - (1.0 - mu * mu));
            num += w * eps;
            den += w;
        }
        num / den
    }

    #[test]
    fn kahn_sampler_reproduces_the_klein_nishina_mean_energy() {
        // A scattering function equal to Z everywhere switches binding off.
        let flat = Curve::new(vec![1e-3, 100.0], vec![1.0, 1.0]).unwrap();
        for e in [3.0e4, 1.0e5, 6.62e5] {
            let alpha = e / MC2_EV;
            let mut rng = stream(7, 0);
            let n = 200_000;
            let mean: f64 = (0..n)
                .map(|_| sample_incoherent(&mut rng, e, 1, &flat).1 / e)
                .sum::<f64>()
                / f64::from(n);
            let want = kn_mean_epsilon(alpha);
            assert!((mean - want).abs() < 4e-3, "E = {e}: {mean} vs {want}");
        }
    }

    #[test]
    fn binding_suppresses_forward_incoherent_scattering() {
        let t = PhotonTable::epdl2017();
        let pb = t.atom(82).unwrap();
        let mut rng = stream(3, 1);
        let n = 20_000;
        let mean_cos = |s: &Curve, z, rng: &mut _| {
            (0..n)
                .map(|_| sample_incoherent(rng, 2.0e4, z, s).0)
                .sum::<f64>()
                / f64::from(n)
        };
        let flat = Curve::new(vec![1e-3, 100.0], vec![82.0, 82.0]).unwrap();
        let free = mean_cos(&flat, 82, &mut rng);
        let bound = mean_cos(&pb.scattering_function, 82, &mut rng);
        assert!(bound < free - 0.05, "free {free}, bound {bound}");
    }

    #[test]
    fn coherent_sampler_matches_the_numerical_distribution() {
        let t = PhotonTable::epdl2017();
        let a = t.atom(13).unwrap();
        let e = 6.0e4;
        let sampler = CoherentSampler::new(&a.form_factor);
        let x_m = e / HC_EV_ANGSTROM;
        // Numerical mean of cos θ with density (1 + μ²) F²(x_m sqrt((1-μ)/2)).
        let n = 400_000;
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..n {
            let mu = -1.0 + 2.0 * (f64::from(i) + 0.5) / f64::from(n);
            let x = x_m * ((1.0 - mu) * 0.5).sqrt();
            let w = (1.0 + mu * mu) * a.form_factor.eval(x).powi(2);
            num += w * mu;
            den += w;
        }
        let want = num / den;
        let mut rng = stream(11, 2);
        let m = 200_000;
        let got: f64 = (0..m).map(|_| sampler.sample_cos(&mut rng, e)).sum::<f64>() / f64::from(m);
        assert!((got - want).abs() < 3e-3, "{got} vs {want}");
        assert!(want > 0.5, "coherent scattering is forward peaked");
    }

    #[test]
    fn coherent_integral_and_inverse_agree() {
        let t = PhotonTable::epdl2017();
        let s = CoherentSampler::new(&t.atom(29).unwrap().form_factor);
        for &tt in &[1e-6, 3e-4, 0.02, 1.0, 40.0, 5000.0] {
            let a = s.integral(tt);
            assert!((s.invert(a) / tt - 1.0).abs() < 1e-9, "t = {tt}");
        }
    }

    #[test]
    fn rotation_makes_the_requested_angle_and_stays_unit() {
        for dir in [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            normalize([0.3, -0.5, 0.8]),
        ] {
            for (c, phi) in [(0.3, 0.7), (-0.9, 4.0), (1.0, 1.0), (-1.0, 2.0)] {
                let r = rotate(dir, c, phi);
                let n = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
                assert!((n - 1.0).abs() < 1e-12);
                let d = r[0] * dir[0] + r[1] * dir[1] + r[2] * dir[2];
                assert!((d - c).abs() < 1e-9, "{d} vs {c}");
            }
        }
    }
}
