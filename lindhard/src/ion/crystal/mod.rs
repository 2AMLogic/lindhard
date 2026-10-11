//! Crystalline targets: the cubic and hexagonal lattice models and the
//! geometry that sets the ion direction in the crystal frame (steps 19a and
//! 19b of the M2 crystal plan).
//!
//! This module is the data model. The binary-collision engine
//! ([`crate::ion::bca`]) uses it for the regions named in
//! [`Bca::with_crystal`](crate::ion::bca::Bca::with_crystal) (the crystal
//! flight model, see [`crate::ion::bca::crystal`]); everywhere else it is
//! unchanged and amorphous runs are bit-identical. It provides three pieces:
//!
//! * [`Lattice`] ([`lattice`]): Bravais (primitive) vectors, a basis of sites
//!   with their species, and lattice parameters stated with their
//!   temperature. Constructors exist for the cubic **diamond** (Si, Ge) and
//!   **zincblende** (GaAs, 3C-SiC) structures and the hexagonal **wurtzite**
//!   and **polytype** (stacking sequence, e.g. 4H, 6H) structures. Presets
//!   carry cited lattice parameters ([`LatticeConstant`],
//!   [`HexagonalConstants`]): Si, Ge, GaAs, 3C-SiC, GaN, 4H-SiC and 6H-SiC.
//!   Miller-index plane normals and lattice directions are computed from the
//!   reciprocal lattice; hexagonal lattices also take four-index
//!   (Miller-Bravais) planes `(hkil)` and directions `[uvtw]`, with the
//!   convention in the [`lattice`] module docs.
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
//!   conventional cell, so the direction `[100]` is `(1, 0, 0)`. For a
//!   hexagonal lattice, `z` is the `c` axis `[0001]` and `x` is along
//!   `[11-20]` (`a1 + a2`). Lengths are metres ([`crate::units`]).
//! * **Lab frame**: the frame of [`crate::geometry`]: `x` is depth (the inward
//!   surface normal), `y` and `z` are lateral. A beam direction with polar
//!   angle `θ` from `+x` and azimuth `φ` from `+y` towards `+z` is
//!   `(cos θ, sin θ cos φ, sin θ sin φ)`, the same as
//!   [`crate::ion::bca::Beam::direction`].
//!
//! Thermal vibration of lattice atoms is the Debye model in [`debye`]
//! (amplitude and displacement sampling). The crystal flight model of the
//! engine uses it when a crystal is given a temperature
//! ([`CrystalTarget::thermal`](crate::ion::bca::CrystalTarget::thermal)).
//!
//! [`search`] finds the lattice sites within an impact parameter of a path
//! segment by walking rectangular cells ([`Lattice::orthogonal_cell`]: the
//! cube, or the orthohexagonal cell), ordered along the path (geometry only;
//! the crystal flight model of the engine calls it for every segment).

pub mod debye;
pub mod divergence;
pub mod lattice;
pub mod orientation;
pub mod search;

pub use divergence::Divergence;
pub use lattice::{HexagonalConstants, Lattice, LatticeConstant, Site, Structure};
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
    /// A four-index (Miller-Bravais) plane or direction broke its constraint.
    #[error("invalid Miller-Bravais indices {indices:?}: {why}")]
    InvalidMillerBravais {
        /// The four indices as given.
        indices: [i32; 4],
        /// The rule they broke.
        why: &'static str,
    },
    /// Four-index input was used with a lattice that is not hexagonal.
    #[error("Miller-Bravais (four-index) input needs a hexagonal lattice")]
    NotHexagonal,
    /// A wurtzite or polytype internal parameter `u` was outside `(0, 1/2)`
    /// or not finite.
    #[error("invalid internal parameter u = {0}: must be finite and in (0, 1/2)")]
    InvalidInternalParameter(f64),
    /// A polytype stacking sequence was not valid.
    #[error("invalid stacking sequence {stacking:?}: {why}")]
    InvalidStacking {
        /// The sequence as given.
        stacking: String,
        /// The rule it broke.
        why: &'static str,
    },
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
