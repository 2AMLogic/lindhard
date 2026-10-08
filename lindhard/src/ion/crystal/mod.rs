//! Crystalline targets: the cubic lattice model and the geometry that sets
//! the ion direction in the crystal frame (step 19a of the M2 crystal plan).
//!
//! This module is the data model. The binary-collision engine
//! ([`crate::ion::bca`]) uses it for the regions named in
//! [`Bca::with_crystal`](crate::ion::bca::Bca::with_crystal) (the crystal
//! flight model, see [`crate::ion::bca::crystal`]); everywhere else it is
//! unchanged and amorphous runs are bit-identical. It provides three pieces:
//!
//! * [`Lattice`] ([`lattice`]): Bravais (primitive) vectors, a basis of sites
//!   with their species, and a lattice constant stated with its temperature.
//!   Constructors exist for the **diamond** (Si, Ge) and **zincblende**
//!   (GaAs, 3C-SiC) structures, and four presets carry cited lattice
//!   constants ([`LatticeConstant`]). Miller-index plane normals and
//!   lattice directions are computed from the reciprocal lattice.
//! * [`Orientation`] ([`orientation`]): the wafer cut (surface normal and an
//!   in-plane reference direction, both as indices), and the beam tilt, twist
//!   and wafer rotation. It gives the lab-to-crystal rotation and the beam
//!   direction in the crystal frame. The conventions are written out in the
//!   [`orientation`] module docs and in `docs/crystal-orientation.md`, with a
//!   figure.
//! * [`Divergence`] ([`divergence`]): beam divergence sampled on the
//!   per-particle random stream ([`crate::rng::stream`]): a Gaussian
//!   (independent normal deviations in two orthogonal angles) or a uniform
//!   cone.
//!
//! # Frames
//!
//! * **Crystal frame**: Cartesian axes along the cube edges of the
//!   conventional cell, so the direction `[100]` is `(1, 0, 0)`. Lengths are
//!   metres ([`crate::units`]).
//! * **Lab frame**: the frame of [`crate::geometry`]: `x` is depth (the inward
//!   surface normal), `y` and `z` are lateral. A beam direction with polar
//!   angle `θ` from `+x` and azimuth `φ` from `+y` towards `+z` is
//!   `(cos θ, sin θ cos φ, sin θ sin φ)`, the same as
//!   [`crate::ion::bca::Beam::direction`].
//!
//! Thermal vibration of lattice atoms is the Debye model in [`debye`]
//! (displacement sampling only; it is not yet used by the engine, which
//! treats the lattice as static).
//!
//! [`search`] finds the lattice sites within an impact parameter of a path
//! segment by walking unit cells, ordered along the path (geometry only; the
//! crystal flight model of the engine calls it for every segment).
//!
//! Hexagonal lattices (wurtzite GaN, 4H/6H-SiC) are a later step and are not
//! here.

pub mod debye;
pub mod divergence;
pub mod lattice;
pub mod orientation;
pub mod search;

pub use divergence::Divergence;
pub use lattice::{Lattice, LatticeConstant, Site, Structure};
pub use orientation::Orientation;
pub use search::{Candidate, LatticeSearch};

/// Errors from building lattices, orientations and divergence models.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CrystalError {
    /// The lattice constant was zero, negative or not finite.
    #[error("invalid lattice constant {0} m: must be finite and positive")]
    InvalidLatticeConstant(f64),
    /// The temperature was negative or not finite.
    #[error("invalid temperature {0} K: must be finite and non-negative")]
    InvalidTemperature(f64),
    /// Atomic number outside the element table.
    #[error("unknown atomic number {0}")]
    UnknownAtomicNumber(u8),
    /// A Miller index or direction index triple was all zero.
    #[error("index triple {0:?} is all zero")]
    ZeroIndex([i32; 3]),
    /// The in-plane reference direction does not lie in the surface plane
    /// (the zone law `h u + k v + l w = 0` fails).
    #[error(
        "reference direction {reference:?} does not lie in the ({normal:?}) plane: \
         h u + k v + l w = {dot}, must be 0"
    )]
    ReferenceNotInPlane {
        /// Miller indices of the surface plane.
        normal: [i32; 3],
        /// Indices of the reference direction.
        reference: [i32; 3],
        /// The offending zone-law sum.
        dot: i64,
    },
    /// A lattice-search argument was invalid.
    #[error("invalid search {name}: {why}")]
    InvalidSearch {
        /// Which argument.
        name: &'static str,
        /// The rule it broke.
        why: &'static str,
    },
    /// An angle was out of range or not finite.
    #[error("invalid {name} {value} rad: {why}")]
    InvalidAngle {
        /// Which angle.
        name: &'static str,
        /// Offending value.
        value: f64,
        /// The rule it broke.
        why: &'static str,
    },
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

pub(crate) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub(crate) fn unit(v: [f64; 3]) -> [f64; 3] {
    scale(v, 1.0 / dot(v, v).sqrt())
}
