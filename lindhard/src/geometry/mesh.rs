//! Triangle-mesh solids: loaders (STL, OBJ), the watertightness check, and the
//! [`MeshGeometry`] target built on a BVH.
//!
//! # Solids and regions
//!
//! A [`TriMesh`] is one closed, outward-oriented triangle surface (possibly
//! several shells, for example an outer surface and an inward-facing cavity).
//! A [`MeshGeometry`] is an ordered list of solids, each tagged with an index
//! into a shared material list. **One region per solid**, in the order listed.
//! Space outside every solid is vacuum.
//!
//! # Ownership rules
//!
//! These are the overlap, boundary and tolerance rules of the mesh target.
//! Other solid representations of the same target (such as a
//! constructive-solid-geometry one) are meant to follow them.
//!
//! 1. **Overlap.** Where solids overlap, the first listed solid owns the
//!    space.
//! 2. **Surface points, directional.** A flight is owned, at each point, by the
//!    solid it is heading into. Concretely, the state of a solid at the start
//!    of a flight is its state at `pos + tol * dir`, so a start on a surface
//!    belongs to the solid on the side the flight enters, as the faces of a
//!    [`VoxelGrid`](super::VoxelGrid) do. [`Geometry::locate`] has no
//!    direction; it uses the fixed offset `tol * (1, 1, 1) / sqrt(3)`, so
//!    lower faces (and edges, vertices) of an axis-aligned solid belong to it
//!    and upper ones do not, exactly as voxel faces do.
//! 3. **No re-hit.** Crossings at `t <= tol` along the flight are ignored, so a
//!    particle that has just crossed a surface and sits on it (to rounding)
//!    does not cross it again.
//! 4. **Edges and vertices.** A flight through an edge or vertex shared by
//!    several triangles crosses the surface once: the ray-triangle test
//!    resolves exact ties by a fixed rule that is consistent between
//!    neighbours (see the `bvh` module source), so no crossing is doubled or
//!    lost. Crossings of different solids within `tol` of each other are one
//!    event, so a flight through a surface shared by two solids passes
//!    directly from one to the other and never through vacuum between them.
//!
//! The tolerance is `tol = 1e-10 x` the diagonal of the bounding box of all
//! solids (see [`MeshGeometry::tolerance_m`]); it is far above rounding
//! (about `1e-16` of a coordinate) and below any feature worth resolving.
//! Two crossings of one solid closer than `tol` cancel.
//!
//! # Events
//!
//! [`Geometry::exit`] reports `Enter { region, pos }` when the owner changes
//! to a solid of a *different material* (a change between solids of the same
//! material is not an event, as crossing a voxel face inside one block is
//! not), and `Escape` when the flight leaves every solid. [`Geometry::flight`]
//! is overridden to report the solid a clear flight ends in. The `Face` of an
//! escape comes from the outward unit normal `n` of the triangle crossed:
//! [`Face::Front`] if `n_x` is negative and the largest component in
//! magnitude, [`Face::Back`] if positive and largest, [`Face::Side`]
//! otherwise (for flat axis-aligned faces this is the voxel convention; on a
//! curved surface the target is split by the dominant axis).
//!
//! The start state of every other solid costs one unbounded ray query, so a
//! flight in a scene of `n` solids costs `O(n log m)`. A single solid takes
//! the fast path.
//!
//! # Point in solid
//!
//! The parity of the number of crossings of a ray from the point to infinity
//! (J. Foley et al., *Computer Graphics: Principles and Practice*, 2nd ed.,
//! Addison-Wesley 1990, section 15.10, ray-casting inside tests). The ray
//! direction is fixed, and exact edge or vertex hits are resolved by the same
//! tie-break as above, so the result is deterministic.
//!
//! # File formats and units
//!
//! Lengths are metres at the API. A file's coordinates are multiplied by the
//! `unit_m` argument of the loaders (`1e-9` for a file in nm, `1.0` for
//! metres); neither format records units.
//!
//! * **STL** (3D Systems, *StereoLithography Interface Specification*,
//!   1989): binary (80-byte header, `u32` count `n`, then `n` records of 12
//!   little-endian `f32` -- normal and three vertices -- and a 2-byte
//!   attribute) and ASCII (`solid` ... `facet normal` ... `vertex x y z` ...).
//!   Binary is recognised by the length `84 + 50 n`, not by the header text
//!   (binary files often begin with `solid`). The stored normals are ignored:
//!   orientation comes from the vertex order (counter-clockwise seen from
//!   outside), and is checked.
//! * **OBJ** (Wavefront; the file format specification of Alias|Wavefront,
//!   1992): `v x y z [w]` and `f` with `i`, `i/t`, `i//n` or `i/t/n`
//!   vertices, 1-based or negative (relative to the vertices read so far).
//!   Polygons are fan-triangulated from the first vertex. All faces of a file
//!   form one mesh; other statements (`vt`, `vn`, `o`, `g`, `usemtl`, `s`,
//!   ...) are skipped.
//!
//! # Watertightness
//!
//! Vertices are welded by exact coordinate equality (the formats repeat
//! them). Then every undirected edge must belong to exactly two triangles that
//! traverse it in opposite directions (a closed, consistently oriented
//! 2-manifold), no triangle may be degenerate (zero area or a repeated
//! vertex, or a sliver thinner than `1e-12` of its longest edge), and the
//! signed volume must be positive, which means the normals point outward.
//! Violations return a [`GeometryError`] naming a bad edge by the positions of
//! its two vertices.

use std::collections::BTreeMap;

use super::bvh::{Bvh, Ray, V3};
use super::{Exit, ExitOutcome, Face, Flight, Geometry, GeometryError};
use crate::material::Material;

/// Relative size of the geometric tolerance: `tol = TOL_REL x` the scene
/// diagonal.
const TOL_REL: f64 = 1e-10;

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// A closed, watertight, outward-oriented triangle mesh with welded vertices.
#[derive(Debug, Clone, PartialEq)]
pub struct TriMesh {
    vertices: Vec<V3>,
    triangles: Vec<[u32; 3]>,
    volume_m3: f64,
    lo: V3,
    hi: V3,
}

impl TriMesh {
    /// Build from a triangle soup (counter-clockwise seen from outside),
    /// coordinates in metres. Welds vertices and runs the watertightness
    /// checks of the module docs.
    pub fn from_triangles(soup: &[[V3; 3]]) -> Result<Self, GeometryError> {
        if soup.is_empty() {
            return Err(GeometryError::MeshEmpty);
        }
        let mut index: BTreeMap<[u64; 3], u32> = BTreeMap::new();
        let mut vertices: Vec<V3> = Vec::new();
        let mut triangles = Vec::with_capacity(soup.len());
        for (ti, tri) in soup.iter().enumerate() {
            let mut t = [0u32; 3];
            for (k, v) in tri.iter().enumerate() {
                if !v.iter().all(|c| c.is_finite()) {
                    return Err(GeometryError::MeshNonFinite { triangle: ti });
                }
                // `+ 0.0` turns -0.0 into 0.0 so both weld.
                let key = [
                    (v[0] + 0.0).to_bits(),
                    (v[1] + 0.0).to_bits(),
                    (v[2] + 0.0).to_bits(),
                ];
                t[k] = *index.entry(key).or_insert_with(|| {
                    vertices.push([v[0] + 0.0, v[1] + 0.0, v[2] + 0.0]);
                    (vertices.len() - 1) as u32
                });
            }
            triangles.push(t);
        }
        let pos = |i: u32| vertices[i as usize];

        for (ti, t) in triangles.iter().enumerate() {
            let [a, b, c] = [pos(t[0]), pos(t[1]), pos(t[2])];
            let n = norm(cross(sub(b, a), sub(c, a)));
            let longest = norm(sub(b, a)).max(norm(sub(c, b))).max(norm(sub(a, c)));
            if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] || n <= 1e-12 * longest * longest {
                return Err(GeometryError::MeshDegenerateTriangle {
                    triangle: ti,
                    vertices: [a, b, c],
                });
            }
        }

        // Directed edges: each at most once, and each with its reverse.
        let mut directed: BTreeMap<(u32, u32), ()> = BTreeMap::new();
        for t in &triangles {
            for k in 0..3 {
                let e = (t[k], t[(k + 1) % 3]);
                if directed.insert(e, ()).is_some() {
                    return Err(GeometryError::MeshInconsistentOrientation {
                        edge: [pos(e.0), pos(e.1)],
                    });
                }
            }
        }
        for &(a, b) in directed.keys() {
            if !directed.contains_key(&(b, a)) {
                return Err(GeometryError::MeshEdgeNotShared {
                    edge: [pos(a), pos(b)],
                });
            }
        }

        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for v in &vertices {
            for k in 0..3 {
                lo[k] = lo[k].min(v[k]);
                hi[k] = hi[k].max(v[k]);
            }
        }
        // Signed volume as a sum of tetrahedra from the box centre.
        let r = [
            0.5 * (lo[0] + hi[0]),
            0.5 * (lo[1] + hi[1]),
            0.5 * (lo[2] + hi[2]),
        ];
        let volume_m3 = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = [sub(pos(t[0]), r), sub(pos(t[1]), r), sub(pos(t[2]), r)];
                dot(a, cross(b, c)) / 6.0
            })
            .sum::<f64>();
        if volume_m3.is_nan() || volume_m3 <= 0.0 {
            return Err(GeometryError::MeshNotOutward { volume_m3 });
        }
        Ok(Self {
            vertices,
            triangles,
            volume_m3,
            lo,
            hi,
        })
    }

    /// Parse an STL file (ASCII or binary), scaling coordinates by `unit_m`
    /// (metres per file unit).
    pub fn from_stl(bytes: &[u8], unit_m: f64) -> Result<Self, GeometryError> {
        check_unit(unit_m)?;
        let perr = |line: usize, message: String| GeometryError::MeshParse {
            format: "STL",
            line,
            message,
        };
        let is_binary = bytes.len() >= 84 && {
            let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as u64;
            bytes.len() as u64 == 84 + 50 * n
        };
        let mut soup: Vec<[V3; 3]> = Vec::new();
        if is_binary {
            for rec in bytes[84..].chunks_exact(50) {
                let f = |o: usize| {
                    f64::from(f32::from_le_bytes([
                        rec[o],
                        rec[o + 1],
                        rec[o + 2],
                        rec[o + 3],
                    ])) * unit_m
                };
                let v = |k: usize| [f(12 + 12 * k), f(16 + 12 * k), f(20 + 12 * k)];
                soup.push([v(0), v(1), v(2)]);
            }
        } else {
            let text = std::str::from_utf8(bytes).map_err(|_| {
                perr(
                    0,
                    "neither a binary STL of the declared length nor UTF-8 text".into(),
                )
            })?;
            if !text.trim_start().starts_with("solid") {
                return Err(perr(1, "not an STL file: expected a `solid` header".into()));
            }
            let mut cur: Vec<V3> = Vec::new();
            for (ln, line) in text.lines().enumerate() {
                let mut it = line.split_whitespace();
                if it.next() != Some("vertex") {
                    continue;
                }
                let mut c = [0.0; 3];
                for slot in &mut c {
                    *slot = it
                        .next()
                        .and_then(|s| s.parse::<f64>().ok())
                        .ok_or_else(|| perr(ln + 1, "malformed `vertex` line".into()))?
                        * unit_m;
                }
                cur.push(c);
                if cur.len() == 3 {
                    soup.push([cur[0], cur[1], cur[2]]);
                    cur.clear();
                }
            }
            if !cur.is_empty() {
                return Err(perr(
                    text.lines().count(),
                    "`vertex` lines are not a multiple of three".into(),
                ));
            }
        }
        Self::from_triangles(&soup)
    }

    /// Parse a Wavefront OBJ file, scaling coordinates by `unit_m` (metres per
    /// file unit). Faces with more than three vertices are fan-triangulated.
    pub fn from_obj(text: &str, unit_m: f64) -> Result<Self, GeometryError> {
        check_unit(unit_m)?;
        let perr = |line: usize, message: String| GeometryError::MeshParse {
            format: "OBJ",
            line,
            message,
        };
        let mut verts: Vec<V3> = Vec::new();
        let mut soup: Vec<[V3; 3]> = Vec::new();
        for (ln, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap_or("");
            let mut it = line.split_whitespace();
            match it.next() {
                Some("v") => {
                    let mut c = [0.0; 3];
                    for slot in &mut c {
                        *slot = it
                            .next()
                            .and_then(|s| s.parse::<f64>().ok())
                            .ok_or_else(|| perr(ln + 1, "malformed `v` line".into()))?
                            * unit_m;
                    }
                    verts.push(c);
                }
                Some("f") => {
                    let mut poly: Vec<V3> = Vec::new();
                    for tok in it {
                        let first = tok.split('/').next().unwrap_or("");
                        let i: i64 = first
                            .parse()
                            .map_err(|_| perr(ln + 1, format!("bad face index `{tok}`")))?;
                        let k = if i > 0 {
                            i - 1
                        } else if i < 0 {
                            verts.len() as i64 + i
                        } else {
                            return Err(perr(ln + 1, "face index 0 is invalid".into()));
                        };
                        if k < 0 || k as usize >= verts.len() {
                            return Err(perr(
                                ln + 1,
                                format!("face index {i} refers to a vertex not yet defined"),
                            ));
                        }
                        poly.push(verts[k as usize]);
                    }
                    if poly.len() < 3 {
                        return Err(perr(ln + 1, "a face needs at least three vertices".into()));
                    }
                    for k in 1..poly.len() - 1 {
                        soup.push([poly[0], poly[k], poly[k + 1]]);
                    }
                }
                _ => {}
            }
        }
        Self::from_triangles(&soup)
    }

    /// Welded vertex positions, m.
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }

    /// Triangles as indices into [`TriMesh::vertices`], counter-clockwise
    /// seen from outside.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// Enclosed volume, m^3 (positive for an outward-oriented mesh).
    pub fn volume_m3(&self) -> f64 {
        self.volume_m3
    }

    /// Axis-aligned bounds `(min, max)`, m.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        (self.lo, self.hi)
    }

    fn soup(&self) -> Vec<[V3; 3]> {
        self.triangles
            .iter()
            .map(|t| t.map(|i| self.vertices[i as usize]))
            .collect()
    }
}

fn check_unit(unit_m: f64) -> Result<(), GeometryError> {
    if unit_m.is_finite() && unit_m > 0.0 {
        Ok(())
    } else {
        Err(GeometryError::MeshInvalidScale { unit_m })
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Solid {
    mesh: TriMesh,
    material: u32,
    bvh: Bvh,
}

/// A target of closed triangle-mesh solids, each tagged with a material, with
/// a BVH per solid. See the module docs for the ownership rules.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshGeometry {
    materials: Vec<Material>,
    solids: Vec<Solid>,
    tol: f64,
    lo: V3,
    hi: V3,
}

impl MeshGeometry {
    /// Build from the material list and `(mesh, material index)` solids, in
    /// priority order (the first listed owns overlaps). Meshes are in metres.
    pub fn new(
        materials: Vec<Material>,
        solids: Vec<(TriMesh, u32)>,
    ) -> Result<Self, GeometryError> {
        if materials.is_empty() {
            return Err(GeometryError::NoMaterials);
        }
        if solids.is_empty() {
            return Err(GeometryError::NoSolids);
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (i, (mesh, m)) in solids.iter().enumerate() {
            if *m as usize >= materials.len() {
                return Err(GeometryError::SolidMaterialOutOfRange {
                    solid: i,
                    index: *m,
                    n_materials: materials.len(),
                });
            }
            for k in 0..3 {
                lo[k] = lo[k].min(mesh.lo[k]);
                hi[k] = hi[k].max(mesh.hi[k]);
            }
        }
        let tol = TOL_REL * norm(sub(hi, lo));
        let solids = solids
            .into_iter()
            .map(|(mesh, material)| {
                let bvh = Bvh::build(&mesh.soup(), tol);
                Solid {
                    mesh,
                    material,
                    bvh,
                }
            })
            .collect();
        Ok(Self {
            materials,
            solids,
            tol,
            lo,
            hi,
        })
    }

    /// The solids' meshes, in region order.
    pub fn mesh(&self, region: usize) -> &TriMesh {
        &self.solids[region].mesh
    }

    /// The geometric tolerance of the ownership rules, m.
    pub fn tolerance_m(&self) -> f64 {
        self.tol
    }

    /// Axis-aligned bounds of all solids, `(min, max)`, m.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        (self.lo, self.hi)
    }

    /// Whether the point `p` is inside solid `i`, by the parity of the
    /// crossings of a fixed ray (module docs).
    fn inside_at(&self, i: usize, p: V3, dir: V3) -> bool {
        let mut n = 0usize;
        let ray = Ray::new(p, dir);
        self.solids[i]
            .bvh
            .for_each_hit(&ray, 0.0, f64::INFINITY, |_, _| n += 1);
        n % 2 == 1
    }

    /// [`Geometry::locate`] with the offset of rule 2 applied with the given
    /// sign on each axis (all +1 for `locate`).
    fn locate_shifted(&self, pos: [f64; 3], sign: [f64; 3]) -> Option<usize> {
        if !pos.iter().all(|c| c.is_finite()) {
            return None;
        }
        // The fixed offset of rule 2 and a fixed, off-axis ray direction.
        let s = self.tol / 3f64.sqrt();
        let p = [
            pos[0] + sign[0] * s,
            pos[1] + sign[1] * s,
            pos[2] + sign[2] * s,
        ];
        let dir = {
            let d = [1.0, 2.0, 3.0];
            let l = norm(d);
            [d[0] / l, d[1] / l, d[2] / l]
        };
        (0..self.solids.len()).find(|&i| {
            let (lo, hi) = self.solids[i].bvh.bounds().expect("non-empty mesh");
            (0..3).all(|k| p[k] >= lo[k] && p[k] <= hi[k]) && self.inside_at(i, p, dir)
        })
    }

    /// Outward unit normal of triangle `tri` (leaf order) of solid `i`.
    fn normal(&self, i: usize, tri: u32) -> V3 {
        let [a, b, c] = *self.solids[i].bvh.triangle(tri);
        let n = cross(sub(b, a), sub(c, a));
        let l = norm(n);
        [n[0] / l, n[1] / l, n[2] / l]
    }
}

/// The face class of an outward normal (module docs).
fn face_of(n: V3) -> Face {
    let (ax, ay, az) = (n[0].abs(), n[1].abs(), n[2].abs());
    if ax > ay && ax > az {
        if n[0] < 0.0 {
            Face::Front
        } else {
            Face::Back
        }
    } else {
        Face::Side
    }
}

impl Geometry for MeshGeometry {
    fn n_materials(&self) -> usize {
        self.materials.len()
    }

    fn material(&self, index: usize) -> &Material {
        &self.materials[index]
    }

    fn n_regions(&self) -> usize {
        self.solids.len()
    }

    fn material_index(&self, region: usize) -> usize {
        self.solids[region].material as usize
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        self.locate_shifted(pos, [1.0; 3])
    }

    fn entry_point(&self) -> [f64; 3] {
        [
            self.lo[0],
            0.5 * (self.lo[1] + self.hi[1]),
            0.5 * (self.lo[2] + self.hi[2]),
        ]
    }

    fn exit(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Option<Exit> {
        match self.flight(region, pos, dir, limit) {
            Flight::Event(e) => Some(e),
            Flight::Clear { .. } => None,
        }
    }

    fn flight(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Flight {
        let clear = Flight::Clear { region };
        if !dir.iter().all(|d| d.is_finite())
            || dir == [0.0; 3]
            || !pos.iter().all(|c| c.is_finite())
            || limit.is_nan()
        {
            return clear;
        }
        let tol = self.tol;
        let ray = Ray::new(pos, dir);
        let n = self.solids.len();

        // Crossings (t, solid, triangle) in (tol, limit + tol], and the
        // state of every solid at t = tol. The flight's own region is taken
        // to be inside.
        let mut hits: Vec<(f64, u32, u32)> = Vec::new();
        let mut inside = vec![false; n];
        inside[region] = true;
        for (j, solid) in self.solids.iter().enumerate() {
            let own = j == region;
            let mut count = 0usize;
            // Other solids need the whole ray for the parity of their start
            // state; the flight's own solid only the flight.
            let tmax = if n == 1 || own {
                limit + tol
            } else {
                f64::INFINITY
            };
            solid.bvh.for_each_hit(&ray, tol, tmax, |t, tri| {
                count += 1;
                if t <= limit + tol {
                    hits.push((t, j as u32, tri));
                }
            });
            if !own {
                inside[j] = count % 2 == 1;
            }
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));

        let owner_of = |inside: &[bool]| inside.iter().position(|&b| b);
        let mut current = region;
        // Start-state mismatch (a surface start heading into a solid that
        // has priority or is not `region`): an event at distance 0.
        let mut pending_start = true;
        let mut i = 0usize;
        loop {
            let (t0, cluster_end) = if pending_start {
                (0.0, i)
            } else {
                if i >= hits.len() || hits[i].0 > limit {
                    return Flight::Clear { region: current };
                }
                let t0 = hits[i].0;
                let mut e = i;
                while e < hits.len() && hits[e].0 - t0 <= tol {
                    inside[hits[e].1 as usize] ^= true;
                    e += 1;
                }
                (t0, e)
            };
            let cluster = &hits[i..cluster_end];
            pending_start = false;
            i = cluster_end;
            let at = [
                pos[0] + t0 * dir[0],
                pos[1] + t0 * dir[1],
                pos[2] + t0 * dir[2],
            ];
            match owner_of(&inside) {
                Some(k) if k == current => {}
                Some(k) => {
                    if self.solids[k].material == self.solids[current].material {
                        current = k;
                    } else {
                        return Flight::Event(Exit {
                            distance: t0,
                            at,
                            outcome: ExitOutcome::Enter { region: k, pos: at },
                        });
                    }
                }
                None => {
                    // The triangle of the solid being left, else any.
                    let (_, sol, tri) = cluster
                        .iter()
                        .copied()
                        .find(|h| h.1 as usize == current)
                        .or_else(|| cluster.first().copied())
                        .unwrap_or((t0, current as u32, u32::MAX));
                    let normal = if tri == u32::MAX {
                        [-dir[0], -dir[1], -dir[2]]
                    } else {
                        self.normal(sol as usize, tri)
                    };
                    return Flight::Event(Exit {
                        distance: t0,
                        at,
                        outcome: ExitOutcome::Escape {
                            face: face_of(normal),
                            normal,
                        },
                    });
                }
            }
        }
    }

    /// A point on the surface of a solid is not in vacuum, whichever of the
    /// eight diagonal offsets (the one of [`Geometry::locate`] and its
    /// mirror images) would put it outside.
    fn in_vacuum(&self, pos: [f64; 3]) -> bool {
        (0..8).all(|m| {
            let sign = |k: usize| if m >> k & 1 == 0 { 1.0 } else { -1.0 };
            self.locate_shifted(pos, [sign(0), sign(1), sign(2)])
                .is_none()
        })
    }
}
