//! Target geometry: the engine-facing [`Geometry`] trait, the 1D layered
//! [`Stack`], the 3D [`VoxelGrid`], the triangle-mesh [`MeshGeometry`] and the
//! constructive-solid-geometry [`CsgGeometry`].
//!
//! # Coordinates
//!
//! `x` is depth, measured from the front surface (`x = 0`) into the target;
//! `y` and `z` are lateral. Lengths are metres (SI, see [`crate::units`]).
//! A [`Stack`] is an ordered list of layers, each homogeneous and laterally
//! infinite (the planar slab model of amorphous-target binary-collision codes,
//! e.g. J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980)
//! 257). The last layer may be semi-infinite (a substrate); otherwise the
//! target has a back face at its total thickness.
//!
//! # The `Geometry` contract
//!
//! A geometry is a set of numbered *regions* (layers of a stack, voxels of a
//! grid), each filled by one of a numbered list of *materials*. A transport
//! engine keeps a particle's current region, asks [`Geometry::exit`] how far
//! a straight flight can go before something changes, truncates the free path
//! there, and carries the unused optical depth into the next region (the
//! same contract as the layer interface of the M0 engine). "Something
//! changes" means the material changes, the particle leaves the target, or
//! it wraps through a periodic boundary; crossing into a region of the same
//! material is *not* an event, so a voxel boundary inside a homogeneous block
//! costs nothing. Per-material data (scattering, stopping) is cached by
//! material index, never per region.
//!
//! # Voxel grids
//!
//! [`VoxelGrid`] is a regular box `[0, Lx] x [0, Ly] x [0, Lz]`, with
//! `L = n h` for `n` voxels of spacing `h` on each axis, and the origin at the
//! corner of the front face (`x = 0` is the front surface, so a grid can
//! encode a [`Stack`] with the layers along `x`). Voxel `(i, j, k)` has the
//! flat region index `i + nx (j + ny k)` (x fastest). It owns the half-open
//! box `[i h, (i + 1) h)` on each axis, so a point exactly on an interior
//! face belongs to the voxel with the higher index, as a depth exactly on a
//! layer interface belongs to the deeper layer of a stack.
//!
//! Each axis is either [`Boundary::Vacuum`] (a particle reaching either end
//! of the axis leaves the target) or [`Boundary::Periodic`] (it re-enters at
//! the opposite face with the same direction; the travelled length and so the
//! energy loss are those of the unwrapped path). Faces on the `x` axis are
//! reported as [`Face::Front`] (`x = 0`) and [`Face::Back`] (`x = Lx`); faces
//! on `y` and `z` are [`Face::Side`].
//!
//! The exact ray traversal is the 3D digital differential analyser of
//! J. Amanatides and A. Woo, "A fast voxel traversal algorithm for ray
//! tracing", Eurographics '87 (1987) 3-10: the distance to the next face of
//! each axis is `(face - position) / direction`, the axis with the smallest
//! distance is crossed, and the face distances are recomputed from the face
//! coordinates (not accumulated), so rounding does not drift along a long
//! ray. Axes with a zero direction component are never crossed. When several
//! axes reach a face at exactly the same distance (an edge or corner), they
//! are crossed together, so the ray passes directly to the diagonal voxel
//! without visiting the voxels it only touches. If a tied crossing leaves the
//! target through several vacuum faces, the face of the lowest axis (`x`,
//! then `y`, then `z`) is reported.

//! # Triangle meshes
//!
//! [`MeshGeometry`] (module [`mesh`]) is a set of closed STL/OBJ triangle
//! solids, one region each, with a BVH for ray queries. Its overlap, surface
//! and tolerance rules, the `Face` mapping of its escapes, the loaders and the
//! watertightness rules are in the [`mesh`] module docs.

//! # CSG
//!
//! [`CsgGeometry`] (module [`csg`]) is a set of [`Csg`] solids built from box,
//! cylinder and half-space [`Primitive`]s by union, intersection and
//! difference, one region per top-level solid, classified along each flight by
//! ray casting. It follows the mesh ownership and tolerance rules, so a box
//! gives the same events as a CSG solid, a mesh or a voxel grid; the rules,
//! the handling of unbounded half-spaces and the algorithms are in the
//! [`csg`] module docs.

mod bvh;
pub mod csg;
pub mod mesh;
mod stack;
#[cfg(test)]
mod test_util;
mod voxel;

pub use csg::{Csg, CsgGeometry, Primitive};
pub use mesh::{MeshGeometry, TriMesh};
pub use stack::{Layer, Stack};
pub use voxel::{Boundary, VoxelGrid};

use crate::material::Material;

/// Errors from building a [`Stack`], a [`VoxelGrid`], a [`MeshGeometry`] or a
/// [`CsgGeometry`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum GeometryError {
    /// No layers and no substrate.
    #[error("a stack needs at least one layer or a substrate")]
    Empty,
    /// A layer thickness was zero, negative or not finite.
    #[error("layer {index} has invalid thickness {thickness_m} m: must be finite and positive")]
    InvalidThickness {
        /// Layer index (0 = front).
        index: usize,
        /// Offending value.
        thickness_m: f64,
    },
    /// A voxel grid needs at least one material.
    #[error("a voxel grid needs at least one material")]
    NoMaterials,
    /// A voxel grid dimension was zero.
    #[error("voxel grid dimension on axis {axis} is zero")]
    ZeroDimension {
        /// Axis (0 = x, 1 = y, 2 = z).
        axis: usize,
    },
    /// A voxel spacing was zero, negative, not finite, or so large that the
    /// grid extent overflows.
    #[error("voxel spacing {spacing_m} m on axis {axis} is invalid: must be finite and positive, with a finite extent")]
    InvalidSpacing {
        /// Axis (0 = x, 1 = y, 2 = z).
        axis: usize,
        /// Offending value.
        spacing_m: f64,
    },
    /// The number of voxels does not fit in `usize`.
    #[error("voxel grid dimensions {dims:?} overflow the voxel count")]
    TooManyVoxels {
        /// Requested dimensions.
        dims: [usize; 3],
    },
    /// The material-index array does not have one entry per voxel.
    #[error("voxel grid needs {expected} material indices, got {actual}")]
    CellCountMismatch {
        /// `nx ny nz`.
        expected: usize,
        /// Length supplied.
        actual: usize,
    },
    /// A voxel names a material that does not exist.
    #[error("voxel {cell} has material index {index}, but there are only {n_materials} materials")]
    MaterialIndexOutOfRange {
        /// Flat voxel index.
        cell: usize,
        /// Offending material index.
        index: u32,
        /// Number of materials supplied.
        n_materials: usize,
    },
    /// A mesh or CSG geometry needs at least one solid.
    #[error("a mesh or CSG geometry needs at least one solid")]
    NoSolids,
    /// A solid names a material that does not exist.
    #[error(
        "solid {solid} has material index {index}, but there are only {n_materials} materials"
    )]
    SolidMaterialOutOfRange {
        /// Solid index (listing order).
        solid: usize,
        /// Offending material index.
        index: u32,
        /// Number of materials supplied.
        n_materials: usize,
    },
    /// A mesh file could not be parsed. `line` is 1-based (0 for a binary file).
    #[error("{format} parse error (line {line}): {message}")]
    MeshParse {
        /// `"STL"` or `"OBJ"`.
        format: &'static str,
        /// Line number, or 0 when not applicable.
        line: usize,
        /// What was wrong.
        message: String,
    },
    /// The file-unit scale of a mesh loader was not finite and positive.
    #[error("mesh unit scale {unit_m} m is invalid: must be finite and positive")]
    MeshInvalidScale {
        /// Offending value.
        unit_m: f64,
    },
    /// A mesh has no triangles.
    #[error("the mesh has no triangles")]
    MeshEmpty,
    /// A mesh vertex coordinate is NaN or infinite.
    #[error("triangle {triangle} has a non-finite coordinate")]
    MeshNonFinite {
        /// Triangle index (file order).
        triangle: usize,
    },
    /// A triangle has zero (or negligible) area.
    #[error("triangle {triangle} is degenerate (zero area): vertices {vertices:?}")]
    MeshDegenerateTriangle {
        /// Triangle index (file order).
        triangle: usize,
        /// Its three vertex positions, m.
        vertices: [[f64; 3]; 3],
    },
    /// An edge belongs to only one triangle: the mesh is not watertight.
    #[error("mesh is not watertight: edge {edge:?} has no neighbouring triangle traversing it in the opposite direction")]
    MeshEdgeNotShared {
        /// The two vertex positions, m.
        edge: [[f64; 3]; 2],
    },
    /// Two triangles traverse an edge in the same direction (inconsistent
    /// orientation, or more than two triangles on the edge).
    #[error(
        "mesh orientation is inconsistent: edge {edge:?} is traversed twice in the same direction"
    )]
    MeshInconsistentOrientation {
        /// The two vertex positions, m.
        edge: [[f64; 3]; 2],
    },
    /// The signed volume is not positive: the normals point inward.
    #[error("mesh has signed volume {volume_m3} m^3: triangles must be counter-clockwise seen from outside")]
    MeshNotOutward {
        /// The signed volume, m^3.
        volume_m3: f64,
    },
    /// A CSG solid has an invalid primitive (non-finite or non-positive
    /// size, zero axis or normal), an empty union or intersection, or an
    /// empty bounding box.
    #[error("CSG solid {solid} is invalid: {reason}")]
    CsgInvalid {
        /// Solid index (listing order).
        solid: usize,
        /// What was wrong.
        reason: &'static str,
    },
    /// A top-level CSG solid is unbounded (a half-space not bounded by the
    /// combination it is in; see the [`csg`] module docs).
    #[error("CSG solid {solid} is unbounded: half-spaces must be bounded by an intersection or difference")]
    CsgUnbounded {
        /// Solid index (listing order).
        solid: usize,
    },
}

/// Which face of the target a particle left through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Face {
    /// The front surface, `x = 0` (backscattering, sputtering).
    Front,
    /// The back face, `x = Lx` of a finite stack or voxel grid (transmission).
    Back,
    /// A lateral face of a voxel grid, `y = 0`, `y = Ly`, `z = 0` or
    /// `z = Lz` (only reachable with a [`Boundary::Vacuum`] lateral axis).
    Side,
}

/// What happens at the end of a flight segment returned by
/// [`Geometry::exit`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExitOutcome {
    /// The particle continues in `region` at `pos` (the crossing point,
    /// wrapped into the domain if a periodic boundary was crossed).
    Enter {
        /// The region entered.
        region: usize,
        /// Position after the crossing, m.
        pos: [f64; 3],
    },
    /// The particle reaches the target surface.
    Escape {
        /// Which face.
        face: Face,
        /// Outward unit normal of the face.
        normal: [f64; 3],
    },
}

/// The first event along a straight flight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Exit {
    /// Path length from the start of the flight to the event, m (the flight is
    /// truncated here).
    pub distance: f64,
    /// The point where the event happens, on the crossed face (before any
    /// periodic wrap), m.
    pub at: [f64; 3],
    /// What the event is.
    pub outcome: ExitOutcome,
}

/// The outcome of a flight query, [`Geometry::flight`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Flight {
    /// An event happens within the flight limit (see [`Geometry::exit`]).
    Event(Exit),
    /// No event within the limit: the whole flight stays in one material,
    /// and ends in `region`. This may differ from the starting region when
    /// the flight crossed faces between regions of the same material (voxels
    /// of one block). A face reached exactly at the limit counts as crossed,
    /// so `region` is the region the particle is entering there, the same
    /// directional ownership as an [`ExitOutcome::Enter`].
    Clear {
        /// The region containing the end of the flight.
        region: usize,
    },
}

/// The engine-facing description of a target. See the module docs for the
/// contract. Implemented by [`Stack`], [`VoxelGrid`], [`MeshGeometry`] and
/// [`CsgGeometry`].
pub trait Geometry: Sync {
    /// Number of distinct materials.
    fn n_materials(&self) -> usize;

    /// Material `index` (`0..n_materials()`).
    fn material(&self, index: usize) -> &Material;

    /// Number of regions.
    fn n_regions(&self) -> usize;

    /// Index of the material filling `region`.
    fn material_index(&self, region: usize) -> usize;

    /// The region containing `pos`, or `None` outside the target (or for a
    /// non-finite position). Ownership of shared faces: see the module docs.
    fn locate(&self, pos: [f64; 3]) -> Option<usize>;

    /// Where a beam enters by default: the origin for a [`Stack`], the centre
    /// of the front face for a [`VoxelGrid`].
    fn entry_point(&self) -> [f64; 3];

    /// The first event along the straight flight from `pos` (inside `region`,
    /// or on one of its faces) in the unit direction `dir`, if it happens at
    /// a path length of at most `limit`; otherwise `None`. The event is a
    /// change of material, a periodic wrap, or leaving the target.
    ///
    /// `dir` must be a finite unit vector; a zero or non-finite direction
    /// returns `None`.
    fn exit(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Option<Exit>;

    /// [`Geometry::exit`], also reporting the region the flight ends in when
    /// there is no event, so a transport engine can keep its particle's
    /// region consistent with its position after crossing same-material
    /// faces. The default suits a geometry whose regions are each a whole
    /// material run along every flight (a [`Stack`], where leaving a layer is
    /// always an event): without an event the region is unchanged.
    fn flight(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Flight {
        match self.exit(region, pos, dir, limit) {
            Some(e) => Flight::Event(e),
            None => Flight::Clear { region },
        }
    }

    /// Whether `pos` lies outside the target in vacuum, used to skip weak-
    /// collision partners that would lie beyond a surface. For a [`Stack`]
    /// this is in front of the front face only (as in the M0 engine).
    fn in_vacuum(&self, pos: [f64; 3]) -> bool;
}
