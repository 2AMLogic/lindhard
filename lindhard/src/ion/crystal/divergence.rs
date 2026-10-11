//! Beam divergence: the spread of incidence directions about the nominal
//! beam direction, sampled on the per-particle random stream.
//!
//! A divergence model gives each particle a deflection from the central
//! direction: a polar angle `ϑ` from it and an azimuth `ψ` about it. The
//! azimuth is uniform on `[0, 2π)` for both models, and the polar angle has:
//!
//! * [`Divergence::Gaussian`] with standard deviation `σ`: the two angular
//!   deviations `(ϑ cos ψ, ϑ sin ψ)` in orthogonal planes through the central
//!   direction are independent normal variables of mean 0 and standard
//!   deviation `σ`. Changing the bivariate normal density
//!   `exp(-(ξ² + η²) / 2σ²) / (2π σ²) dξ dη` to polar coordinates gives
//!   `(1 / 2π) dψ · (ϑ / σ²) exp(-ϑ² / 2σ²) dϑ`: `ψ` uniform and independent
//!   of `ϑ`, and `ϑ` Rayleigh distributed with CDF
//!   `F(ϑ) = 1 - exp(-ϑ² / 2σ²)`. It is drawn by inversion,
//!   `ϑ = σ sqrt(-2 ln(1 - U))`. The angular deviations are treated as plane
//!   angles, which is the small-angle reading of a Gaussian divergence; a
//!   draw with `ϑ > π` has probability `exp(-π² / 2σ²)` and is still
//!   applied as a rotation by `ϑ`.
//! * [`Divergence::UniformCone`] with half-angle `α`: directions uniform in
//!   solid angle inside the cone `ϑ <= α`. The solid-angle element is
//!   `sin ϑ dϑ dψ`, so `F(ϑ) = (1 - cos ϑ) / (1 - cos α)`. Using
//!   `1 - cos x = 2 sin²(x / 2)`, inversion gives
//!   `ϑ = 2 asin(sqrt(U) sin(α / 2))`, which keeps full relative precision
//!   for small angles.
//!
//! # Random numbers
//!
//! Draws come from the caller's generator, which for a transport run is the
//! particle's own stream [`crate::rng::stream`]`(seed, index)`, so the
//! result depends only on `(seed, index)` and not on the thread count. Each
//! sample takes exactly two 64-bit draws, in this order: `U` for the polar
//! angle, then `V` for the azimuth `ψ = 2π V`. Uniforms are the top 53 bits
//! of a draw scaled to `[0, 1)`, as in [`crate::ion::bca`].
//! [`Divergence::None`] takes no draws, so a beam without divergence leaves
//! the stream untouched.
//!
//! The deflected direction is the central direction rotated by `ϑ` at
//! azimuth `ψ` with [`crate::ion::bca::kinematics::rotate`]; its azimuth
//! zero is a fixed, documented frame, which is immaterial because `ψ` is
//! uniform.

use std::f64::consts::PI;

use rand_core::Rng;

use super::CrystalError;
use crate::ion::bca::kinematics::rotate;

/// Beam divergence model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Divergence {
    /// No divergence: every particle gets the central direction.
    None,
    /// Gaussian divergence: independent normal angular deviations with
    /// standard deviation `sigma_rad` in two orthogonal planes.
    Gaussian {
        /// Standard deviation per plane, radians (finite, `>= 0`).
        sigma_rad: f64,
    },
    /// Directions uniform in solid angle within a cone.
    UniformCone {
        /// Cone half-angle, radians, in `[0, π]`.
        half_angle_rad: f64,
    },
}

fn uniform<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

impl Divergence {
    /// Check the parameters.
    pub fn validate(&self) -> Result<(), CrystalError> {
        match *self {
            Self::None => Ok(()),
            Self::Gaussian { sigma_rad: s } if s.is_finite() && s >= 0.0 => Ok(()),
            Self::Gaussian { sigma_rad } => Err(CrystalError::InvalidAngle {
                name: "divergence sigma",
                value: sigma_rad,
                why: "must be finite and non-negative",
            }),
            Self::UniformCone { half_angle_rad: a } if (0.0..=PI).contains(&a) => Ok(()),
            Self::UniformCone { half_angle_rad } => Err(CrystalError::InvalidAngle {
                name: "divergence half-angle",
                value: half_angle_rad,
                why: "must be in [0, pi]",
            }),
        }
    }

    /// Target CDF of the polar deflection `ϑ` (see the module docs).
    /// [`Divergence::None`] is a point mass at 0.
    pub fn polar_cdf(&self, theta: f64) -> f64 {
        if theta < 0.0 {
            return 0.0;
        }
        match *self {
            Self::None => 1.0,
            Self::Gaussian { sigma_rad: s } => {
                if s == 0.0 {
                    1.0
                } else {
                    -(-theta * theta / (2.0 * s * s)).exp_m1()
                }
            }
            Self::UniformCone { half_angle_rad: a } => {
                if theta >= a {
                    1.0
                } else {
                    let h = (0.5 * theta).sin() / (0.5 * a).sin();
                    h * h
                }
            }
        }
    }

    /// Draw a deflection `(ϑ, ψ)`: polar angle from the central direction
    /// and azimuth about it, radians. Takes two draws (none for
    /// [`Divergence::None`]); the parameters must be valid
    /// ([`Self::validate`]).
    pub fn sample_deflection<R: Rng + ?Sized>(&self, rng: &mut R) -> (f64, f64) {
        let polar = match *self {
            Self::None => return (0.0, 0.0),
            Self::Gaussian { sigma_rad } => {
                let u = uniform(rng);
                // 1 - u is in (0, 1], so the logarithm is finite.
                sigma_rad * (-2.0 * (-u).ln_1p()).sqrt()
            }
            Self::UniformCone { half_angle_rad } => {
                let u = uniform(rng);
                2.0 * (u.sqrt() * (0.5 * half_angle_rad).sin()).asin()
            }
        };
        let azimuth = 2.0 * PI * uniform(rng);
        (polar, azimuth)
    }

    /// Deflect the unit vector `central` by a sampled divergence. For
    /// [`Divergence::None`] this returns `central` unchanged and draws
    /// nothing.
    pub fn sample_direction<R: Rng + ?Sized>(&self, central: [f64; 3], rng: &mut R) -> [f64; 3] {
        match self {
            Self::None => central,
            _ => {
                let (polar, azimuth) = self.sample_deflection(rng);
                rotate(central, polar, azimuth)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::stream;

    #[test]
    fn none_draws_nothing() {
        let mut a = stream(3, 7);
        let d = Divergence::None.sample_direction([1.0, 0.0, 0.0], &mut a);
        assert_eq!(d, [1.0, 0.0, 0.0]);
        assert_eq!(a.next_u64(), stream(3, 7).next_u64());
    }

    #[test]
    fn each_sample_takes_two_draws() {
        for div in [
            Divergence::Gaussian { sigma_rad: 0.01 },
            Divergence::UniformCone {
                half_angle_rad: 0.02,
            },
        ] {
            let mut a = stream(3, 7);
            div.sample_deflection(&mut a);
            let mut b = stream(3, 7);
            b.next_u64();
            b.next_u64();
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn cone_stays_inside_and_zero_width_is_exact() {
        let div = Divergence::UniformCone {
            half_angle_rad: 0.05,
        };
        let mut rng = stream(1, 0);
        for _ in 0..1000 {
            let (p, a) = div.sample_deflection(&mut rng);
            assert!((0.0..=0.05).contains(&p));
            assert!((0.0..2.0 * PI).contains(&a));
        }
        let (p, _) = Divergence::Gaussian { sigma_rad: 0.0 }.sample_deflection(&mut rng);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn cdfs_are_normalised() {
        let g = Divergence::Gaussian { sigma_rad: 0.01 };
        assert_eq!(g.polar_cdf(0.0), 0.0);
        assert!((g.polar_cdf(1.0) - 1.0).abs() < 1e-15);
        let c = Divergence::UniformCone {
            half_angle_rad: 0.3,
        };
        assert_eq!(c.polar_cdf(0.0), 0.0);
        assert_eq!(c.polar_cdf(0.3), 1.0);
        // Half the solid angle: 1 - cos(t) = (1 - cos(0.3)) / 2.
        let t = (1.0 - 0.5 * (1.0 - 0.3f64.cos())).acos();
        assert!((c.polar_cdf(t) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn validation() {
        assert!(Divergence::None.validate().is_ok());
        assert!(Divergence::Gaussian { sigma_rad: 0.0 }.validate().is_ok());
        assert!(Divergence::Gaussian { sigma_rad: -1.0 }.validate().is_err());
        assert!(Divergence::Gaussian {
            sigma_rad: f64::NAN
        }
        .validate()
        .is_err());
        assert!(Divergence::UniformCone { half_angle_rad: PI }
            .validate()
            .is_ok());
        assert!(Divergence::UniformCone {
            half_angle_rad: 4.0
        }
        .validate()
        .is_err());
    }
}
