//! Target geometry: the engine-facing [`Geometry`] trait, the 1D layered
//! [`Stack`], the 3D [`VoxelGrid`] and the triangle-mesh [`MeshGeometry`].
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

mod bvh;
pub mod mesh;

pub use mesh::{MeshGeometry, TriMesh};

use crate::material::Material;

/// Errors from building a [`Stack`], a [`VoxelGrid`] or a [`MeshGeometry`].
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
    /// A mesh geometry needs at least one solid.
    #[error("a mesh geometry needs at least one solid")]
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
/// contract. Implemented by [`Stack`] and [`VoxelGrid`].
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

/// One homogeneous layer.
#[derive(Debug, Clone, PartialEq)]
pub struct Layer {
    material: Material,
    /// Depth of the front face, m.
    front_m: f64,
    /// Depth of the back face, m (`INFINITY` for a substrate).
    back_m: f64,
}

impl Layer {
    /// The layer material.
    pub fn material(&self) -> &Material {
        &self.material
    }

    /// Depth of the front face, m.
    pub fn front_m(&self) -> f64 {
        self.front_m
    }

    /// Depth of the back face, m (`f64::INFINITY` for a semi-infinite substrate).
    pub fn back_m(&self) -> f64 {
        self.back_m
    }

    /// Thickness, m (`f64::INFINITY` for a semi-infinite substrate).
    pub fn thickness_m(&self) -> f64 {
        self.back_m - self.front_m
    }
}

/// A 1D layered target: finite layers front to back, optionally ending in a
/// semi-infinite substrate.
#[derive(Debug, Clone, PartialEq)]
pub struct Stack {
    layers: Vec<Layer>,
}

impl Stack {
    /// Build from `(material, thickness in m)` pairs, front first, and an
    /// optional semi-infinite substrate behind them.
    pub fn new(
        layers: Vec<(Material, f64)>,
        substrate: Option<Material>,
    ) -> Result<Self, GeometryError> {
        if layers.is_empty() && substrate.is_none() {
            return Err(GeometryError::Empty);
        }
        let mut out = Vec::with_capacity(layers.len() + 1);
        let mut x = 0.0f64;
        for (index, (material, thickness_m)) in layers.into_iter().enumerate() {
            if !(thickness_m.is_finite() && thickness_m > 0.0) {
                return Err(GeometryError::InvalidThickness { index, thickness_m });
            }
            let back = x + thickness_m;
            out.push(Layer {
                material,
                front_m: x,
                back_m: back,
            });
            x = back;
        }
        if let Some(material) = substrate {
            out.push(Layer {
                material,
                front_m: x,
                back_m: f64::INFINITY,
            });
        }
        Ok(Self { layers: out })
    }

    /// A single semi-infinite material.
    pub fn semi_infinite(material: Material) -> Self {
        Self {
            layers: vec![Layer {
                material,
                front_m: 0.0,
                back_m: f64::INFINITY,
            }],
        }
    }

    /// Layers front to back.
    pub fn layers(&self) -> &[Layer] {
        &self.layers
    }

    /// Depth of the back face, m, or `None` if the last layer is semi-infinite.
    pub fn back_face_m(&self) -> Option<f64> {
        let b = self.layers.last().expect("non-empty").back_m;
        b.is_finite().then_some(b)
    }

    /// Index of the layer containing depth `x`, or `None` outside the target.
    /// A depth exactly on an interface belongs to the deeper layer.
    pub fn layer_at(&self, x: f64) -> Option<usize> {
        if x.is_nan() || x < 0.0 {
            return None;
        }
        self.layers.iter().position(|l| x < l.back_m)
    }
}

impl Geometry for Stack {
    fn n_materials(&self) -> usize {
        self.layers.len()
    }

    fn material(&self, index: usize) -> &Material {
        &self.layers[index].material
    }

    fn n_regions(&self) -> usize {
        self.layers.len()
    }

    /// Each layer is its own material entry, so layer data is cached per
    /// layer.
    fn material_index(&self, region: usize) -> usize {
        region
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        self.layer_at(pos[0])
    }

    fn entry_point(&self) -> [f64; 3] {
        [0.0; 3]
    }

    fn exit(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Option<Exit> {
        let layer = &self.layers[region];
        let (d, boundary) = if dir[0] > 0.0 {
            ((layer.back_m - pos[0]) / dir[0], layer.back_m)
        } else if dir[0] < 0.0 {
            ((layer.front_m - pos[0]) / dir[0], layer.front_m)
        } else {
            return None;
        };
        if d.is_nan() || d > limit {
            return None;
        }
        let at = [boundary, pos[1] + d * dir[1], pos[2] + d * dir[2]];
        let outward = dir[0] > 0.0;
        let outcome = if !outward && region == 0 {
            ExitOutcome::Escape {
                face: Face::Front,
                normal: [-1.0, 0.0, 0.0],
            }
        } else if outward && region + 1 == self.layers.len() {
            ExitOutcome::Escape {
                face: Face::Back,
                normal: [1.0, 0.0, 0.0],
            }
        } else {
            ExitOutcome::Enter {
                region: if outward { region + 1 } else { region - 1 },
                pos: at,
            }
        };
        Some(Exit {
            distance: d,
            at,
            outcome,
        })
    }

    fn in_vacuum(&self, pos: [f64; 3]) -> bool {
        pos[0] < 0.0
    }
}

/// What a voxel grid does to a particle that reaches the ends of one axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// The particle leaves the target (vacuum beyond the face).
    Vacuum,
    /// The particle re-enters at the opposite face, unchanged in direction
    /// and energy.
    Periodic,
}

/// A regular 3D grid of voxels, each filled by one material. See the module
/// docs for the origin, flattening order, face ownership, boundary
/// behaviour and the traversal algorithm.
///
/// ```
/// use lindhard::geometry::{Boundary, Geometry, VoxelGrid};
/// use lindhard::material::Material;
///
/// let si = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
/// let c = Material::from_atom_fractions(&[(6, 1.0)], None).unwrap();
/// // 2 x 1 x 1 voxels of 10 nm: silicon in front of carbon.
/// let grid = VoxelGrid::new(
///     vec![si, c],
///     [2, 1, 1],
///     [1e-8; 3],
///     vec![0, 1],
///     [Boundary::Vacuum, Boundary::Periodic, Boundary::Periodic],
/// )
/// .unwrap();
/// assert_eq!(grid.locate([1.5e-8, 0.5e-8, 0.5e-8]), Some(1));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelGrid {
    materials: Vec<Material>,
    dims: [usize; 3],
    spacing_m: [f64; 3],
    extent_m: [f64; 3],
    boundaries: [Boundary; 3],
    cells: Vec<u32>,
}

impl VoxelGrid {
    /// Build a grid of `dims` voxels of edge lengths `spacing_m`, where voxel
    /// `i + nx (j + ny k)` has material `cells[...]` (an index into
    /// `materials`), with `boundaries` per axis `[x, y, z]`.
    ///
    /// Every dimension must be nonzero, every spacing finite and positive
    /// with a finite extent `n h`, the voxel count must not overflow, and
    /// every material index must exist.
    pub fn new(
        materials: Vec<Material>,
        dims: [usize; 3],
        spacing_m: [f64; 3],
        cells: Vec<u32>,
        boundaries: [Boundary; 3],
    ) -> Result<Self, GeometryError> {
        if materials.is_empty() {
            return Err(GeometryError::NoMaterials);
        }
        let mut extent_m = [0.0; 3];
        for axis in 0..3 {
            if dims[axis] == 0 {
                return Err(GeometryError::ZeroDimension { axis });
            }
            let h = spacing_m[axis];
            let l = dims[axis] as f64 * h;
            if !(h.is_finite() && h > 0.0 && l.is_finite()) {
                return Err(GeometryError::InvalidSpacing { axis, spacing_m: h });
            }
            extent_m[axis] = l;
        }
        let expected = dims[0]
            .checked_mul(dims[1])
            .and_then(|v| v.checked_mul(dims[2]))
            .ok_or(GeometryError::TooManyVoxels { dims })?;
        if cells.len() != expected {
            return Err(GeometryError::CellCountMismatch {
                expected,
                actual: cells.len(),
            });
        }
        if let Some((cell, &index)) = cells
            .iter()
            .enumerate()
            .find(|&(_, &m)| m as usize >= materials.len())
        {
            return Err(GeometryError::MaterialIndexOutOfRange {
                cell,
                index,
                n_materials: materials.len(),
            });
        }
        Ok(Self {
            materials,
            dims,
            spacing_m,
            extent_m,
            boundaries,
            cells,
        })
    }

    /// Voxels per axis `[nx, ny, nz]`.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Voxel edge lengths `[hx, hy, hz]`, m.
    pub fn spacing_m(&self) -> [f64; 3] {
        self.spacing_m
    }

    /// Box size `[Lx, Ly, Lz]` = `n h`, m.
    pub fn extent_m(&self) -> [f64; 3] {
        self.extent_m
    }

    /// Boundary behaviour of axis `axis` (0 = x, 1 = y, 2 = z).
    pub fn boundary(&self, axis: usize) -> Boundary {
        self.boundaries[axis]
    }

    /// Flat region index of voxel `(i, j, k)`, or `None` out of range.
    pub fn flat_index(&self, ijk: [usize; 3]) -> Option<usize> {
        (ijk[0] < self.dims[0] && ijk[1] < self.dims[1] && ijk[2] < self.dims[2])
            .then(|| ijk[0] + self.dims[0] * (ijk[1] + self.dims[1] * ijk[2]))
    }

    /// Voxel `(i, j, k)` of a flat region index.
    pub fn voxel_of(&self, region: usize) -> [usize; 3] {
        let [nx, ny, _] = self.dims;
        [region % nx, (region / nx) % ny, region / (nx * ny)]
    }

    /// Material index of voxel `(i, j, k)`, or `None` out of range.
    pub fn cell_material(&self, ijk: [usize; 3]) -> Option<usize> {
        self.flat_index(ijk).map(|r| self.cells[r] as usize)
    }

    fn flat_unchecked(&self, ijk: [usize; 3]) -> usize {
        ijk[0] + self.dims[0] * (ijk[1] + self.dims[1] * ijk[2])
    }
}

impl Geometry for VoxelGrid {
    fn n_materials(&self) -> usize {
        self.materials.len()
    }

    fn material(&self, index: usize) -> &Material {
        &self.materials[index]
    }

    fn n_regions(&self) -> usize {
        self.cells.len()
    }

    fn material_index(&self, region: usize) -> usize {
        self.cells[region] as usize
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        let mut ijk = [0usize; 3];
        for a in 0..3 {
            let x = pos[a];
            let l = self.extent_m[a];
            if !(x >= 0.0 && x <= l) {
                return None;
            }
            ijk[a] = if x == l {
                // The upper face belongs to the periodic image (voxel 0).
                if self.boundaries[a] == Boundary::Vacuum {
                    return None;
                }
                0
            } else {
                // `x / h` can round up to `n` just below the face.
                ((x / self.spacing_m[a]) as usize).min(self.dims[a] - 1)
            };
        }
        Some(self.flat_unchecked(ijk))
    }

    fn entry_point(&self) -> [f64; 3] {
        [0.0, 0.5 * self.extent_m[1], 0.5 * self.extent_m[2]]
    }

    fn exit(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Option<Exit> {
        match self.flight(region, pos, dir, limit) {
            Flight::Event(e) => Some(e),
            Flight::Clear { .. } => None,
        }
    }

    /// The traversal behind [`Geometry::exit`]. Without an event, the region
    /// is the voxel the DDA has reached at the limit: same-material faces up
    /// to and including the limit have been crossed, so a flight ending
    /// exactly on a face is in the voxel beyond it along `dir`. A periodic
    /// wrap is always an event, so a clear flight never wraps.
    fn flight(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Flight {
        if !dir.iter().all(|d| d.is_finite()) || dir == [0.0; 3] {
            return Flight::Clear { region };
        }
        let home = self.cells[region];
        let mut cell = self.voxel_of(region);
        let mut step = [0i64; 3];
        for a in 0..3 {
            step[a] = if dir[a] > 0.0 {
                1
            } else if dir[a] < 0.0 {
                -1
            } else {
                0
            };
        }
        // Face crossed next on each axis, as a coordinate: the upper face of
        // the voxel when moving up, the lower face when moving down.
        let face_coord = |a: usize, c: usize| -> f64 {
            let k = if step[a] > 0 { c + 1 } else { c };
            k as f64 * self.spacing_m[a]
        };
        let t_of = |a: usize, c: usize| -> f64 {
            if step[a] == 0 {
                f64::INFINITY
            } else {
                // A start slightly past the face (rounding) crosses at once.
                ((face_coord(a, c) - pos[a]) / dir[a]).max(0.0)
            }
        };
        let mut t_next = [t_of(0, cell[0]), t_of(1, cell[1]), t_of(2, cell[2])];
        loop {
            let t = t_next[0].min(t_next[1]).min(t_next[2]);
            if t.is_nan() || t > limit {
                return Flight::Clear {
                    region: self.flat_unchecked(cell),
                };
            }
            let mut at = [
                pos[0] + t * dir[0],
                pos[1] + t * dir[1],
                pos[2] + t * dir[2],
            ];
            let mut wrapped = false;
            let mut escape: Option<(usize, bool)> = None;
            let mut new_pos = at;
            for a in 0..3 {
                // Tied axes are crossed together (exact equality).
                if step[a] == 0 || t_next[a] != t {
                    continue;
                }
                at[a] = face_coord(a, cell[a]);
                new_pos[a] = at[a];
                let n = self.dims[a];
                if step[a] > 0 {
                    if cell[a] + 1 == n {
                        match self.boundaries[a] {
                            Boundary::Vacuum => {
                                escape.get_or_insert((a, true));
                            }
                            Boundary::Periodic => {
                                wrapped = true;
                                cell[a] = 0;
                                new_pos[a] = 0.0;
                            }
                        }
                    } else {
                        cell[a] += 1;
                    }
                } else if cell[a] == 0 {
                    match self.boundaries[a] {
                        Boundary::Vacuum => {
                            escape.get_or_insert((a, false));
                        }
                        Boundary::Periodic => {
                            wrapped = true;
                            cell[a] = n - 1;
                            new_pos[a] = self.extent_m[a];
                        }
                    }
                } else {
                    cell[a] -= 1;
                }
                t_next[a] = t_of(a, cell[a]);
            }
            if let Some((axis, positive)) = escape {
                let face = match (axis, positive) {
                    (0, false) => Face::Front,
                    (0, true) => Face::Back,
                    _ => Face::Side,
                };
                let mut normal = [0.0; 3];
                normal[axis] = if positive { 1.0 } else { -1.0 };
                return Flight::Event(Exit {
                    distance: t,
                    at,
                    outcome: ExitOutcome::Escape { face, normal },
                });
            }
            let next = self.flat_unchecked(cell);
            if wrapped || self.cells[next] != home {
                // A periodic image is the same voxel set: ownership puts the
                // wrapped coordinate in the voxel just entered.
                return Flight::Event(Exit {
                    distance: t,
                    at,
                    outcome: ExitOutcome::Enter {
                        region: next,
                        pos: new_pos,
                    },
                });
            }
        }
    }

    fn in_vacuum(&self, pos: [f64; 3]) -> bool {
        (0..3).any(|a| {
            self.boundaries[a] == Boundary::Vacuum && (pos[a] < 0.0 || pos[a] > self.extent_m[a])
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn si() -> Material {
        Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
    }

    #[test]
    fn boundaries_accumulate() {
        let s = Stack::new(vec![(si(), 1e-8), (si(), 2e-8)], Some(si())).unwrap();
        let b: Vec<_> = s
            .layers()
            .iter()
            .map(|l| (l.front_m(), l.back_m()))
            .collect();
        assert_eq!(b[0], (0.0, 1e-8));
        assert_eq!(b[1], (1e-8, 1e-8 + 2e-8));
        assert_eq!(b[2].1, f64::INFINITY);
        assert_eq!(s.back_face_m(), None);
        assert_eq!(s.layer_at(0.0), Some(0));
        assert_eq!(s.layer_at(1e-8), Some(1));
        assert_eq!(s.layer_at(1.0), Some(2));
        assert_eq!(s.layer_at(-1e-12), None);
        let f = Stack::new(vec![(si(), 1e-8)], None).unwrap();
        assert_eq!(f.back_face_m(), Some(1e-8));
        assert_eq!(f.layer_at(1e-8), None);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(Stack::new(vec![], None), Err(GeometryError::Empty));
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                Stack::new(vec![(si(), t)], None),
                Err(GeometryError::InvalidThickness { index: 0, .. })
            ));
        }
    }

    fn c() -> Material {
        Material::from_atom_fractions(&[(6, 1.0)], None).unwrap()
    }

    const V: Boundary = Boundary::Vacuum;
    const P: Boundary = Boundary::Periodic;

    /// `dims` voxels of unit spacing, material `f(i, j, k)` (0 or 1).
    fn grid(
        dims: [usize; 3],
        b: [Boundary; 3],
        f: impl Fn(usize, usize, usize) -> u32,
    ) -> VoxelGrid {
        let mut cells = Vec::new();
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    cells.push(f(i, j, k));
                }
            }
        }
        VoxelGrid::new(vec![si(), c()], dims, [1.0; 3], cells, b).unwrap()
    }

    fn unit(v: [f64; 3]) -> [f64; 3] {
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / n, v[1] / n, v[2] / n]
    }

    fn enter(e: &Exit) -> (usize, [f64; 3]) {
        match e.outcome {
            ExitOutcome::Enter { region, pos } => (region, pos),
            ExitOutcome::Escape { .. } => panic!("expected Enter, got {e:?}"),
        }
    }

    #[test]
    fn voxel_validation() {
        let m = || vec![si()];
        let ok = |dims, sp: [f64; 3], cells: Vec<u32>| VoxelGrid::new(m(), dims, sp, cells, [V; 3]);
        assert!(ok([1, 1, 1], [1.0; 3], vec![0]).is_ok());
        assert_eq!(
            VoxelGrid::new(vec![], [1; 3], [1.0; 3], vec![0], [V; 3]),
            Err(GeometryError::NoMaterials)
        );
        assert_eq!(
            ok([1, 0, 1], [1.0; 3], vec![]),
            Err(GeometryError::ZeroDimension { axis: 1 })
        );
        for h in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            let r = ok([2, 1, 1], [1.0, 1.0, h], vec![0; 2]);
            let r = if h == f64::MAX {
                ok([2, 1, 1], [h, 1.0, 1.0], vec![0; 2])
            } else {
                r
            };
            assert!(
                matches!(r, Err(GeometryError::InvalidSpacing { .. })),
                "{h}: {r:?}"
            );
        }
        assert!(matches!(
            ok([usize::MAX, 2, 1], [1.0; 3], vec![]),
            Err(GeometryError::TooManyVoxels { .. })
        ));
        assert_eq!(
            ok([2, 2, 1], [1.0; 3], vec![0; 3]),
            Err(GeometryError::CellCountMismatch {
                expected: 4,
                actual: 3
            })
        );
        assert!(matches!(
            ok([2, 1, 1], [1.0; 3], vec![0, 1]),
            Err(GeometryError::MaterialIndexOutOfRange {
                cell: 1,
                index: 1,
                n_materials: 1
            })
        ));
    }

    #[test]
    fn flattening_and_locate() {
        let g = grid([3, 2, 2], [V, P, V], |i, j, k| {
            u32::from(i == 2 && j == 1 && k == 1)
        });
        // x fastest.
        assert_eq!(g.flat_index([2, 1, 1]), Some(11));
        assert_eq!(g.voxel_of(11), [2, 1, 1]);
        assert_eq!(g.flat_index([3, 0, 0]), None);
        assert_eq!(g.cell_material([2, 1, 1]), Some(1));
        assert_eq!(g.material_index(11), 1);
        assert_eq!(g.n_regions(), 12);
        assert_eq!(g.locate([2.5, 1.5, 1.5]), Some(11));
        // Interior faces belong to the higher voxel.
        assert_eq!(g.locate([1.0, 0.0, 0.0]), Some(1));
        // Upper face: vacuum is outside, a periodic axis wraps to voxel 0.
        assert_eq!(g.locate([3.0, 0.5, 0.5]), None);
        assert_eq!(g.locate([0.5, 2.0, 0.5]), Some(0));
        assert_eq!(g.locate([-1e-12, 0.5, 0.5]), None);
        assert_eq!(g.locate([0.5, 0.5, f64::NAN]), None);
        // Just below the face never rounds into a nonexistent voxel.
        assert_eq!(
            g.locate([f64::from_bits(3.0f64.to_bits() - 1), 0.5, 0.5]),
            Some(2)
        );
        assert_eq!(g.entry_point(), [0.0, 1.0, 1.0]);
        assert!(g.in_vacuum([-0.1, 0.0, 0.0]));
        assert!(g.in_vacuum([0.0, 0.0, 2.1]));
        assert!(!g.in_vacuum([0.0, -5.0, 0.0])); // periodic axis
    }

    #[test]
    fn oblique_ray_stops_at_a_material_change_and_skips_equal_voxels() {
        // Material 1 for x >= 2; y faces are crossed on the way but are inside
        // one material, so they are not events.
        let g = grid([4, 3, 1], [V, V, P], |i, _, _| u32::from(i >= 2));
        let d = unit([2.0, 0.5, 0.0]);
        let start = [0.5, 0.2, 0.5];
        let e = g.exit(0, start, d, 100.0).unwrap();
        let t = 1.5 / d[0];
        assert!((e.distance - t).abs() < 1e-14 * t);
        assert_eq!(e.at[0], 2.0);
        assert!((e.at[1] - (0.2 + t * d[1])).abs() < 1e-14);
        let (r, pos) = enter(&e);
        assert_eq!(r, g.flat_index([2, 0, 0]).unwrap());
        assert_eq!(pos, e.at);
        // Not reached within a shorter limit; reached at exactly the limit.
        assert!(g.exit(0, start, d, 0.99 * t).is_none());
        assert!(g.exit(0, start, d, e.distance).is_some());
    }

    #[test]
    fn oblique_ray_changes_material_across_a_lateral_face() {
        let g = grid([3, 3, 1], [V, V, P], |_, j, _| u32::from(j >= 1));
        let d = unit([1.0, 2.0, 0.0]);
        let e = g.exit(0, [0.5, 0.5, 0.5], d, 10.0).unwrap();
        assert!((e.distance - 0.5 / d[1]).abs() < 1e-14);
        assert_eq!(e.at[1], 1.0);
        assert!((e.at[0] - (0.5 + e.distance * d[0])).abs() < 1e-14);
        assert_eq!(enter(&e).0, g.flat_index([0, 1, 0]).unwrap());
    }

    #[test]
    fn negative_direction_and_starting_on_a_face() {
        let g = grid([4, 1, 1], [V, P, P], |i, _, _| u32::from(i >= 2));
        // From voxel 3 going back: the face at x = 2.
        let e = g.exit(3, [3.5, 0.5, 0.5], [-1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 1.5);
        assert_eq!(enter(&e), (1, [2.0, 0.5, 0.5]));
        // Starting on the face x = 2 (owned by voxel 2): moving down crosses
        // at once, moving up does not see that face.
        assert_eq!(g.locate([2.0, 0.5, 0.5]), Some(2));
        let e = g.exit(2, [2.0, 0.5, 0.5], [-1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 0.0);
        assert_eq!(enter(&e).0, 1);
        let e = g.exit(2, [2.0, 0.5, 0.5], [1.0, 0.0, 0.0], 10.0).unwrap();
        // Same material up to the back face: vacuum escape at x = 4.
        assert_eq!(e.distance, 2.0);
        assert_eq!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Back,
                normal: [1.0, 0.0, 0.0]
            }
        );
    }

    #[test]
    fn edge_and_corner_ties_cross_together() {
        // Only the diagonal voxel (1,1,0) differs: a 45 degree ray from the
        // centre of (0,0,0) hits the corner at (1,1) and enters it directly.
        let g = grid([3, 3, 1], [V, V, P], |i, j, _| u32::from(i == 1 && j == 1));
        let d = unit([1.0, 1.0, 0.0]);
        let e = g.exit(0, [0.5, 0.5, 0.5], d, 10.0).unwrap();
        assert!((e.distance - 0.5 * 2f64.sqrt()).abs() < 1e-15);
        assert_eq!(e.at[..2], [1.0, 1.0]);
        assert_eq!(enter(&e).0, g.flat_index([1, 1, 0]).unwrap());

        // Only the two voxels the ray merely touches differ: it passes the
        // corner without any event and leaves through the far corner, which
        // is reported as the back face (lowest axis first).
        let g = grid([3, 3, 1], [V, V, P], |i, j, _| u32::from(i != j));
        let e = g.exit(0, [0.5, 0.5, 0.5], d, 10.0).unwrap();
        assert!((e.distance - 2.5 * 2f64.sqrt()).abs() < 1e-14);
        assert_eq!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Back,
                normal: [1.0, 0.0, 0.0]
            }
        );
        assert_eq!(e.at[..2], [3.0, 3.0]);
    }

    #[test]
    fn side_faces_are_reported_as_side() {
        let g = grid([2, 2, 2], [P, V, V], |_, _, _| 0);
        let e = g.exit(0, [0.5, 0.5, 0.5], [0.0, -1.0, 0.0], 10.0).unwrap();
        assert_eq!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Side,
                normal: [0.0, -1.0, 0.0]
            }
        );
        assert_eq!(e.distance, 0.5);
        let e = g.exit(0, [0.5, 0.5, 0.5], [0.0, 0.0, 1.0], 10.0).unwrap();
        assert_eq!(e.distance, 1.5);
        assert!(matches!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Side,
                normal: [0.0, 0.0, 1.0]
            }
        ));
    }

    #[test]
    fn periodic_wrap_is_an_event_that_keeps_the_path_length() {
        let g = grid([2, 2, 1], [P, P, P], |_, _, _| 0);
        let e = g.exit(1, [1.5, 0.5, 0.5], [1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 0.5);
        assert_eq!(e.at, [2.0, 0.5, 0.5]);
        assert_eq!(enter(&e), (0, [0.0, 0.5, 0.5]));
        let e = g.exit(0, [0.5, 0.5, 0.5], [-1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 0.5);
        assert_eq!(enter(&e), (1, [2.0, 0.5, 0.5]));
        // Oblique, wrapping on y after crossing (same-material) x faces.
        let d = unit([1.0, 3.0, 0.0]);
        let e = g.exit(0, [0.1, 0.5, 0.5], d, 10.0).unwrap();
        assert!((e.distance - 1.5 / d[1]).abs() < 1e-14);
        let (r, pos) = enter(&e);
        assert_eq!(pos[1], 0.0);
        assert!((pos[0] - (0.1 + e.distance * d[0])).abs() < 1e-14);
        assert_eq!(r, g.flat_index([0, 0, 0]).unwrap());
        // A single-voxel periodic axis wraps into itself.
        let g1 = VoxelGrid::new(vec![si()], [1, 1, 1], [2.0; 3], vec![0], [P; 3]).unwrap();
        let e = g1.exit(0, [0.5, 1.0, 1.0], [-1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 0.5);
        assert_eq!(enter(&e), (0, [2.0, 1.0, 1.0]));
    }

    #[test]
    fn clear_flight_reports_the_region_it_ends_in() {
        // Regression: [2, 2, 1] voxels, material 1 only in voxel (1, 1, 0).
        let g = grid([2, 2, 1], [V; 3], |i, j, _| u32::from(i == 1 && j == 1));
        let (x, y) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        // One unit along +x crosses only a same-material face: no event, but
        // the flight ends in voxel (1, 0, 0).
        let f = g.flight(0, [0.5; 3], x, 1.0);
        assert_eq!(f, Flight::Clear { region: 1 });
        assert!(g.exit(0, [0.5; 3], x, 1.0).is_none());
        // Turning to +y from there meets the material interface at y = 1.
        let e = g.exit(1, [1.5, 0.5, 0.5], y, 1.0).unwrap();
        assert_eq!(e.distance, 0.5);
        assert_eq!(enter(&e), (3, [1.5, 1.0, 0.5]));
        // The stale starting region would miss it (x-cell 0 has no y change).
        assert!(g.exit(0, [1.5, 0.5, 0.5], y, 1.0).is_none());

        // Directional ownership: a flight ending exactly on a face is in the
        // voxel beyond it along the direction, as for an `Enter`.
        assert_eq!(g.flight(0, [0.5; 3], x, 0.5), Flight::Clear { region: 1 });
        let back = [-1.0, 0.0, 0.0];
        assert_eq!(
            g.flight(1, [1.5, 0.5, 0.5], back, 0.5),
            Flight::Clear { region: 0 }
        );
        // Short of the face: unchanged. Degenerate direction: unchanged.
        assert_eq!(g.flight(0, [0.5; 3], x, 0.25), Flight::Clear { region: 0 });
        assert_eq!(
            g.flight(1, [1.5; 3], [0.0; 3], 9.0),
            Flight::Clear { region: 1 }
        );

        // Periodic: after a wrap the particle sits on the upper face in the
        // last voxel (`locate` would give voxel 0); a clear flight from there
        // keeps that ownership, and reaching the wrap is always an event.
        let p = grid([2, 1, 1], [P, P, P], |_, _, _| 0);
        assert_eq!(p.locate([2.0, 0.5, 0.5]), Some(0));
        assert_eq!(
            p.flight(1, [2.0, 0.5, 0.5], back, 0.25),
            Flight::Clear { region: 1 }
        );
        assert!(matches!(
            p.flight(1, [1.5, 0.5, 0.5], x, 0.5),
            Flight::Event(Exit {
                outcome: ExitOutcome::Enter { region: 0, .. },
                ..
            })
        ));

        // A stack keeps its region: leaving a layer is always an event.
        let s = Stack::new(vec![(si(), 1.0), (c(), 2.0)], None).unwrap();
        assert_eq!(s.flight(0, [0.5; 3], x, 0.25), Flight::Clear { region: 0 });
        assert!(matches!(s.flight(0, [0.5; 3], x, 0.5), Flight::Event(_)));
    }

    #[test]
    fn degenerate_directions_have_no_events() {
        let g = grid([2, 1, 1], [V, V, V], |_, _, _| 0);
        for d in [[0.0; 3], [f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]] {
            assert!(g.exit(0, [0.5; 3], d, f64::INFINITY).is_none());
        }
    }

    #[test]
    fn long_rays_do_not_drift() {
        // 1000 voxels along x with a 1e-3 tilt: the crossing of the last
        // x-face is computed from the face coordinate, not accumulated.
        let n = 1000;
        let g = VoxelGrid::new(
            vec![si()],
            [n, 1, 1],
            [0.1, 1e3, 1e3],
            vec![0; n],
            [V, P, P],
        )
        .unwrap();
        let d = unit([1.0, 1e-3, 0.0]);
        let e = g.exit(0, [0.05, 0.05, 0.05], d, 1e9).unwrap();
        assert_eq!(e.at[0], 100.0);
        assert!((e.distance - (100.0 - 0.05) / d[0]).abs() < 1e-12);
    }

    #[test]
    fn stack_exit_matches_the_layer_contract() {
        let s = Stack::new(vec![(si(), 1.0), (c(), 2.0)], None).unwrap();
        assert_eq!(s.n_materials(), 2);
        assert_eq!(s.locate([0.0; 3]), Some(0));
        assert_eq!(s.entry_point(), [0.0; 3]);
        assert!(s.in_vacuum([-1e-9, 0.0, 0.0]) && !s.in_vacuum([4.0, 0.0, 0.0]));
        let d = unit([1.0, 1.0, 0.0]);
        let e = s.exit(0, [0.25, 0.0, 0.0], d, 10.0).unwrap();
        assert_eq!(e.at[0], 1.0);
        assert_eq!(enter(&e).0, 1);
        let e = s.exit(1, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 2.0);
        assert!(matches!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Back,
                ..
            }
        ));
        let e = s.exit(0, [0.5, 0.0, 0.0], [-1.0, 0.0, 0.0], 10.0).unwrap();
        assert!(matches!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Front,
                normal: [-1.0, 0.0, 0.0]
            }
        ));
        // Parallel to the interfaces, or too short: no event.
        assert!(s.exit(0, [0.5; 3], [0.0, 1.0, 0.0], 1e9).is_none());
        assert!(s.exit(0, [0.5; 3], [1.0, 0.0, 0.0], 0.1).is_none());
        // Substrate: no back face.
        let sub = Stack::semi_infinite(si());
        assert!(sub.exit(0, [0.5; 3], [1.0, 0.0, 0.0], 1e9).is_none());
    }
}
