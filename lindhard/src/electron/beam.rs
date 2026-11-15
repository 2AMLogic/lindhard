//! Where a run's primaries start: the lateral beam profile on the front face
//! and the angular spread about the nominal direction.
//!
//! [`BeamSpread`] is the beam of an electron run — a landing point on the
//! front face, a [`BeamProfile`] about it, and a [`Divergence`] about the
//! nominal [`Primary::direction`](crate::electron::transport::Primary::direction).
//! Every profile is a plain geometric law with no tabulated data, and the
//! default is a point beam without divergence, which draws no random numbers.
//! The draw order per history (position, then direction, both from that
//! history's own stream) is what makes a run identical at any thread count;
//! see [`BeamSpread`] for it.
//!
//! The same beam laws appear twice in the input formats, in metres here and
//! in nanometres in
//! [`BeamProfileSpec`](crate::input::electron::BeamProfileSpec), under two
//! historical tag spellings; `beam_serde_module` is the one deserializer
//! both use. [`crate::electron::transport`] re-exports [`BeamProfile`] and
//! [`BeamSpread`], which lived there until this module was split out.

use serde::Serialize;

use rand_core::Rng;

use crate::electron::transport::{invalid, TransportError};
use crate::ion::crystal::Divergence;
use crate::rng::uniform;

/// Lateral profile of the beam on the front face, about its landing point.
/// Lengths are metres.
///
/// Each profile is a plain geometric law (no tabulated data): a Gaussian spot
/// has independent normal `y` and `z` offsets with standard deviation
/// `sigma_m` each (sampled by the Box-Muller transform, two draws); a disc is
/// uniform in area (`r = R sqrt(U)`, azimuth uniform; two draws); a rectangle
/// is uniform in `y` and `z` over its full widths (two draws). A point takes
/// no draws.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[serde(tag = "model", rename_all = "kebab-case")]
pub enum BeamProfile {
    /// Every primary lands exactly at the landing point.
    Point,
    /// Gaussian spot, standard deviation per axis.
    Gaussian {
        /// Standard deviation of each of `y` and `z`, m (finite, `> 0`).
        sigma_m: f64,
    },
    /// Uniform disc.
    Disc {
        /// Disc radius, m (finite, `> 0`).
        radius_m: f64,
    },
    /// Uniform rectangle centred on the landing point.
    Rectangle {
        /// Full width along `y`, m (finite, `> 0`).
        width_y_m: f64,
        /// Full width along `z`, m (finite, `> 0`).
        width_z_m: f64,
    },
}

/// Serde of an exposure beam under the two tag spellings, for one beam type.
///
/// The report's [`BeamProfile`] (metres) and the input's
/// [`BeamProfileSpec`](crate::input::electron::BeamProfileSpec) (nanometres)
/// carry the same beam laws and the same compatibility rule, so they share
/// one implementation: the input is read under either the
/// `[electron.beam.profile]` spelling (`model`, `disc`) or the historical
/// one (`kind`, `uniform-disc`), but not both tags together, and the output
/// is always the historical spelling, so reports and echoed inputs are
/// unchanged. A missing size field, or one belonging to another profile, is
/// an error.
///
/// Serde field names have to be literal identifiers, which is why this is a
/// macro and not a generic function. It expands to a private module with the
/// `serialize` and `deserialize` a `#[serde(with = "…")]` wants, from the
/// module name, the beam enum (in scope where the macro is invoked) and its
/// four size fields:
///
/// ```ignore
/// beam_serde_module!(beam_serde, BeamProfile, sigma_m, radius_m, width_y_m, width_z_m);
/// ```
macro_rules! beam_serde_module {
    ($name:ident, $beam:ident, $sigma:ident, $radius:ident, $width_y:ident, $width_z:ident) => {
        mod $name {
            use super::$beam;
            use serde::{Deserialize, Deserializer, Serialize, Serializer};

            #[derive(Serialize)]
            #[serde(tag = "kind", rename_all = "kebab-case")]
            enum Out {
                Point,
                Gaussian { $sigma: f64 },
                UniformDisc { $radius: f64 },
                Rectangle { $width_y: f64, $width_z: f64 },
            }

            #[derive(Deserialize)]
            #[serde(deny_unknown_fields)]
            struct Raw {
                kind: Option<String>,
                model: Option<String>,
                $sigma: Option<f64>,
                $radius: Option<f64>,
                $width_y: Option<f64>,
                $width_z: Option<f64>,
            }

            pub(super) fn serialize<S: Serializer>(b: &$beam, s: S) -> Result<S::Ok, S::Error> {
                match *b {
                    $beam::Point => Out::Point,
                    $beam::Gaussian { $sigma } => Out::Gaussian { $sigma },
                    $beam::Disc { $radius } => Out::UniformDisc { $radius },
                    $beam::Rectangle { $width_y, $width_z } => {
                        Out::Rectangle { $width_y, $width_z }
                    }
                }
                .serialize(s)
            }

            pub(super) fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<$beam, D::Error> {
                use serde::de::Error;
                let r = Raw::deserialize(d)?;
                let tag = match (&r.kind, &r.model) {
                    (Some(_), Some(_)) => {
                        return Err(D::Error::custom("give `model` or `kind`, not both"))
                    }
                    (Some(t), None) | (None, Some(t)) => t.as_str(),
                    (None, None) => {
                        return Err(D::Error::custom("missing `model` (or legacy `kind`)"))
                    }
                };
                let need = |v: Option<f64>, name: &str| {
                    v.ok_or_else(|| D::Error::custom(format!("`{tag}` beam needs `{name}`")))
                };
                let (used, beam) = match tag {
                    "point" => (0, $beam::Point),
                    "gaussian" => (
                        1,
                        $beam::Gaussian {
                            $sigma: need(r.$sigma, stringify!($sigma))?,
                        },
                    ),
                    "disc" | "uniform-disc" => (
                        2,
                        $beam::Disc {
                            $radius: need(r.$radius, stringify!($radius))?,
                        },
                    ),
                    "rectangle" => (
                        3,
                        $beam::Rectangle {
                            $width_y: need(r.$width_y, stringify!($width_y))?,
                            $width_z: need(r.$width_z, stringify!($width_z))?,
                        },
                    ),
                    other => return Err(D::Error::custom(format!("unknown beam model `{other}`"))),
                };
                let extra = [
                    r.$sigma.is_some() && used != 1,
                    r.$radius.is_some() && used != 2,
                    (r.$width_y.is_some() || r.$width_z.is_some()) && used != 3,
                ];
                if extra.iter().any(|&e| e) {
                    return Err(D::Error::custom(format!("`{tag}` beam: unexpected field")));
                }
                Ok(beam)
            }
        }
    };
}
pub(crate) use beam_serde_module;

/// Where a run's primaries land and how they diverge: a finite lateral
/// profile about a landing point on the front face, and an angular spread
/// about the nominal
/// [`Primary::direction`](crate::electron::transport::Primary::direction).
/// The default is a point beam without divergence, which draws no random
/// numbers, so a run with it is bit-identical to one made before this type
/// existed.
///
/// Per history the position is sampled first, then the direction, both from
/// that history's own stream, so a run is identical at any thread count. The
/// divergence is the ion beam's
/// ([`Divergence`]: Gaussian per-plane
/// angles, or uniform in solid angle within a cone). Where the nominal
/// direction points into the target (`x > 0`), deflections that would send
/// the primary out of it are rejected and redrawn (the inward conditioning of
/// the ion side), so the draw count is then variable but still a function of
/// the history's stream alone.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct BeamSpread {
    /// Landing point `(y, z)` on the front face, m.
    pub landing_m: [f64; 2],
    /// Lateral profile about the landing point.
    pub profile: BeamProfile,
    /// Angular spread about the nominal direction.
    pub divergence: Divergence,
}

impl Default for BeamSpread {
    fn default() -> Self {
        Self {
            landing_m: [0.0; 2],
            profile: BeamProfile::Point,
            divergence: Divergence::None,
        }
    }
}

/// Whether a beam size is finite and strictly positive (shared by the
/// transport and exposure beam validators).
pub(crate) fn finite_positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

/// Lateral offset `(y, z)` of a Gaussian spot from two uniform draws
/// `u1, u2` in `[0, 1)`: the Box-Muller transform (Box & Muller, Ann. Math.
/// Stat. 29, 610 (1958)), `r = sigma sqrt(-2 ln(1 - u1))`, azimuth
/// `phi = 2 pi u2`, offset `(r cos phi, r sin phi)`. `1 - u1` is in `(0, 1]`,
/// so the radius is finite; `ln(1 - u1)` is evaluated as `ln_1p(-u1)`.
pub(crate) fn gaussian_offset(sigma_m: f64, u1: f64, u2: f64) -> [f64; 2] {
    let r = sigma_m * (-2.0 * (-u1).ln_1p()).sqrt();
    let (s, c) = (std::f64::consts::TAU * u2).sin_cos();
    [r * c, r * s]
}

/// Lateral offset `(y, z)` uniform in area on a disc of radius `radius_m`
/// from two uniform draws `u1, u2` in `[0, 1)`: `r = R sqrt(u1)`, azimuth
/// `phi = 2 pi u2`, offset `(r cos phi, r sin phi)`.
pub(crate) fn disc_offset(radius_m: f64, u1: f64, u2: f64) -> [f64; 2] {
    let r = radius_m * u1.sqrt();
    let (s, c) = (std::f64::consts::TAU * u2).sin_cos();
    [r * c, r * s]
}

/// Redraw cap of the inward conditioning of [`BeamSpread`].
const MAX_DIVERGENCE_REDRAWS: u32 = 10_000;

impl BeamSpread {
    /// Whether this is the default point beam (it is then left out of the
    /// run metadata).
    pub fn is_point(&self) -> bool {
        *self == Self::default()
    }

    /// Check the parameters.
    pub fn validate(&self) -> Result<(), TransportError> {
        if !self.landing_m.iter().all(|v| v.is_finite()) {
            return invalid("beam landing point", "must be finite");
        }
        let ok = finite_positive;
        let good = match self.profile {
            BeamProfile::Point => true,
            BeamProfile::Gaussian { sigma_m } => ok(sigma_m),
            BeamProfile::Disc { radius_m } => ok(radius_m),
            BeamProfile::Rectangle {
                width_y_m,
                width_z_m,
            } => ok(width_y_m) && ok(width_z_m),
        };
        if !good {
            return invalid("beam profile", "sizes must be finite and positive");
        }
        self.divergence
            .validate()
            .or_else(|e| invalid("beam divergence", e.to_string()))
    }

    /// Draw the start `(y, z)`, m: the landing point plus a profile offset.
    pub fn sample_position<R: Rng + ?Sized>(&self, rng: &mut R) -> [f64; 2] {
        let [y0, z0] = self.landing_m;
        let (dy, dz) = match self.profile {
            BeamProfile::Point => (0.0, 0.0),
            BeamProfile::Gaussian { sigma_m } => {
                let (u1, u2) = (uniform(rng), uniform(rng));
                let [dy, dz] = gaussian_offset(sigma_m, u1, u2);
                (dy, dz)
            }
            BeamProfile::Disc { radius_m } => {
                let (u1, u2) = (uniform(rng), uniform(rng));
                let [dy, dz] = disc_offset(radius_m, u1, u2);
                (dy, dz)
            }
            BeamProfile::Rectangle {
                width_y_m,
                width_z_m,
            } => (
                width_y_m * (uniform(rng) - 0.5),
                width_z_m * (uniform(rng) - 0.5),
            ),
        };
        [y0 + dy, z0 + dz]
    }

    /// Draw the direction about the unit `central` direction.
    pub(crate) fn sample_direction<R: Rng + ?Sized>(
        &self,
        central: [f64; 3],
        rng: &mut R,
    ) -> Result<[f64; 3], TransportError> {
        if self.divergence == Divergence::None {
            return Ok(central);
        }
        for _ in 0..MAX_DIVERGENCE_REDRAWS {
            let d = self.divergence.sample_direction(central, rng);
            if central[0] <= 0.0 || d[0] > 0.0 {
                return Ok(d);
            }
        }
        invalid(
            "beam divergence",
            "no inward direction was drawn; the nominal direction is too close to the surface",
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn beam_offset_helpers_match_the_inline_formulas() {
        let tau = std::f64::consts::TAU;
        for &(u1, u2) in &[(0.0, 0.0), (0.25, 0.75), (0.5, 0.1), (0.999_999_999, 0.9)] {
            // Disc: exact.
            let r = 3.0e-6 * f64::sqrt(u1);
            let want = [r * (tau * u2).cos(), r * (tau * u2).sin()];
            assert_eq!(disc_offset(3.0e-6, u1, u2), want);
            // Gaussian: the old exposure form `ln(1 - u1)`, to rounding.
            let r = 2.0e-6 * (-2.0 * (1.0 - u1).ln()).sqrt();
            let want = [r * (tau * u2).cos(), r * (tau * u2).sin()];
            let got = gaussian_offset(2.0e-6, u1, u2);
            for k in 0..2 {
                assert!((got[k] - want[k]).abs() <= 1e-15 * want[k].abs().max(1e-9));
            }
            assert!(got.iter().all(|v| v.is_finite()));
        }
        assert_eq!(gaussian_offset(1.0, 0.0, 0.3), [0.0, 0.0]);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(!finite_positive(bad));
        }
        assert!(finite_positive(1e-12));
    }
}
