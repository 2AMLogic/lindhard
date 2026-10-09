//! Constructive solid geometry: box, cylinder and half-space primitives,
//! combined by union, intersection and difference, and the [`CsgGeometry`]
//! target built from them.
//!
//! # Solids and regions
//!
//! A [`Csg`] is a tree whose leaves are [`Primitive`]s and whose inner nodes
//! are [`Csg::Union`], [`Csg::Intersection`] and [`Csg::Difference`]; the
//! operators nest freely. A [`CsgGeometry`] is an ordered list of top-level
//! solids, each tagged with an index into a shared material list. **One region
//! per top-level solid**, in the order listed. Space outside every solid is
//! vacuum ([`Geometry::locate`] returns `None`).
//!
//! # Ownership rules
//!
//! These are the rules of the triangle-mesh target, adopted verbatim (see the
//! [`mesh`](super::mesh) module docs), so that a CSG solid, the same solid as
//! a mesh and the same solid as a [`VoxelGrid`](super::VoxelGrid) agree on
//! every shared face.
//!
//! 1. **Overlap.** Where solids overlap, the first listed solid owns the
//!    space, so a mask listed before a substrate cuts into it. No overlap
//!    search is needed.
//! 2. **Surface points, directional.** A flight is owned, at each point, by
//!    the solid it is heading into: the state of a solid at the start of a
//!    flight is its state at `pos + tol * dir`. [`Geometry::locate`] has no
//!    direction; it uses the fixed offset `tol * (1, 1, 1) / sqrt(3)`, so the
//!    lower faces (and edges, vertices) of an axis-aligned box belong to it and
//!    the upper ones do not, exactly as voxel faces do. Where that offset does
//!    not move a point off a surface (a plane containing the diagonal
//!    `(1, 1, 1)`), and for a flight lying exactly in a surface, the point
//!    belongs to the side that `(1, 1, 1)` points into, and for a plane
//!    containing the diagonal to the side its normal's first non-zero
//!    component points away from. Each plane thus gives its points to exactly
//!    one side.
//! 3. **No re-hit.** Crossings at `t <= tol` along the flight are ignored, so
//!    a particle that has just crossed a surface and sits on it (to rounding)
//!    does not cross it again.
//! 4. **Coincident and tangent surfaces.** Along a flight every primitive and
//!    every boolean node is a sorted list of inside spans (below). At each
//!    node, spans no longer than `tol` are removed and gaps no longer than
//!    `tol` are closed. So a ray that grazes a cylinder or clips a box corner
//!    within `tol` does not enter it; a hole cut flush with a face leaves no
//!    sliver; and two boxes united along a shared face are one span with no
//!    crossing at the face. Crossings of different top-level solids within
//!    `tol` of each other are one event, so a flight through a surface shared
//!    by two solids passes directly from one to the other and never through
//!    vacuum between them. This is the CSG form of the mesh rule that two
//!    crossings of one solid closer than `tol` cancel.
//!
//! The tolerance is `tol = 1e-10 x` the diagonal of the bounding box of all
//! solids, the same constant as the mesh target (see
//! [`CsgGeometry::tolerance_m`]).
//!
//! # Bounded solids
//!
//! A half-space is unbounded, so it may only appear where a combination
//! bounds it: every top-level solid must have a finite axis-aligned bounding
//! box, otherwise [`CsgGeometry::new`] returns
//! [`GeometryError::CsgUnbounded`]. Bounds are taken conservatively: a box's
//! own; a cylinder's exact box (its end discs, each of which extends
//! `r sqrt(1 - a_k^2)` along axis `k` for a unit axis `a`); a half-space whose
//! normal is a coordinate axis is bounded on that one side, any other is
//! unbounded; a union takes the hull of its operands, an intersection their
//! overlap, a difference the bounds of what it cuts from. So a box cut by an
//! oblique half-space, or six axis-aligned half-spaces intersected, are
//! bounded, but a lone half-space or a union containing one is not. A
//! top-level bounding box that is empty on some axis is rejected as
//! [`GeometryError::CsgInvalid`]. The scene bounding box (for `tol` and
//! [`Geometry::entry_point`]) is the hull of the top-level boxes.
//!
//! # Events
//!
//! As for meshes: [`Geometry::exit`] reports `Enter { region, pos }` when the
//! owner changes to a solid of a *different material* (a change between
//! solids of the same material is not an event), and `Escape` when the flight
//! leaves every solid. [`Geometry::flight`] is overridden to report the solid
//! a clear flight ends in. The `Face` of an escape comes from the outward
//! unit normal `n` of the surface crossed (for the subtracted operand of a
//! difference, its inward normal): [`Face::Front`](super::Face::Front) if
//! `n_x` is negative and the largest component in magnitude,
//! [`Face::Back`](super::Face::Back) if positive and largest,
//! [`Face::Side`](super::Face::Side) otherwise.
//!
//! # Ray classification
//!
//! The ray-casting method of S. D. Roth, "Ray casting for modeling solids",
//! Computer Graphics and Image Processing 18 (1982) 109-144: each primitive
//! is classified along the (whole, unbounded) line of the flight into the
//! parameter spans where the line is inside it, and the span lists are
//! combined up the tree by the set operation of each node. Every span end
//! carries the outward normal of the surface it lies on; a difference
//! negates the normals of the subtracted operand. Here union and difference
//! are built from intersection and complement (`A u B` is the complement of
//! `A' n B'`, `A - B` is `A n B'`), and the rule 4 clean-up is applied at
//! every node. Each primitive is convex, so its classification is a single
//! span:
//!
//! * **Box**, axis aligned: the intersection of the three slab intervals of
//!   T. L. Kay and J. T. Kajiya, "Ray tracing complex scenes", Computer
//!   Graphics (Proc. SIGGRAPH '86) 20 (4) (1986) 269-278.
//! * **Cylinder**, finite and capped: the slab interval of its two end planes
//!   (as for the box) intersected with the interval inside the infinite
//!   cylinder, whose ends are the roots of `|w + t v|^2 = r^2` (`w`, `v` the
//!   components of the offset from the base and of the direction normal to
//!   the axis). The roots use the cancellation-free form of the quadratic
//!   formula, `q = -(b + sgn(b) sqrt(b^2 - a c))`, `t = q / a` and `t = c / q`
//!   (W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
//!   *Numerical Recipes*, 3rd ed., Cambridge University Press 2007, section
//!   5.6). A tangent line (zero discriminant) does not enter.
//! * **Half-space** `{p : n . p <= offset}`: the half-line on the inside of
//!   the plane crossing.
//!
//! [`Geometry::locate`] classifies the offset point directly: a union holds it
//! if any operand does, an intersection if all do, a difference if the first
//! operand does and the second does not.
//!
//! A flight classifies every primitive of every solid once, so its cost is
//! linear in the total number of primitives. Nothing is cached between
//! queries and no query depends on thread count or evaluation order.

use super::mesh::{face_of, TOL_REL};
use super::{Exit, ExitOutcome, Flight, Geometry, GeometryError};
use crate::material::Material;

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn neg(a: V3) -> V3 {
    [-a[0], -a[1], -a[2]]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

fn axis(k: usize) -> V3 {
    let mut e = [0.0; 3];
    e[k] = 1.0;
    e
}

/// Whether a point exactly on a surface with outward normal `n` belongs to
/// the solid (rule 2 of the module docs): yes if the diagonal `(1, 1, 1)`
/// points inward; for a surface containing the diagonal, yes if the first
/// non-zero component of `n` is negative. `owns(n)` and `owns(-n)` always
/// differ for a non-zero `n` (negation is exact, so the sum flips sign
/// exactly).
fn owns(n: V3) -> bool {
    let e = n[0] + n[1] + n[2];
    if e != 0.0 {
        return e < 0.0;
    }
    n.iter().find(|c| **c != 0.0).is_some_and(|c| *c < 0.0)
}

/// Whether the coordinate `x` lies in the slab `[lo, hi]` measured along the
/// unit vector `u` (outward normals `-u` at `lo`, `u` at `hi`), with rule-2
/// ownership of the two planes.
fn in_slab(x: f64, lo: f64, hi: f64, u: V3) -> bool {
    (x > lo || (x == lo && owns(neg(u)))) && (x < hi || (x == hi && owns(u)))
}

/// One axis-aligned primitive solid. Lengths in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Primitive {
    /// The axis-aligned box `[lo, hi]` (each `lo[k] < hi[k]`).
    Cuboid {
        /// Lower corner, m.
        lo: [f64; 3],
        /// Upper corner, m.
        hi: [f64; 3],
    },
    /// A finite cylinder with flat end caps: the points `base + s a + w` with
    /// `0 <= s <= height` and `w` normal to `a`, `|w| <= radius`, where `a` is
    /// `axis` normalised.
    Cylinder {
        /// Centre of the first end cap, m.
        base: [f64; 3],
        /// Direction of the axis (any non-zero length; normalised on build).
        axis: [f64; 3],
        /// Radius, m.
        radius: f64,
        /// Length along the axis, m.
        height: f64,
    },
    /// The half-space `{p : normal . p <= offset}`, with `normal` its outward
    /// normal. Scaling `normal` and `offset` together describes the same
    /// half-space (normalised on build). Unbounded: see the module docs.
    HalfSpace {
        /// Outward normal (any non-zero length).
        normal: [f64; 3],
        /// Offset, in units of `|normal|` m.
        offset: f64,
    },
}

/// A constructive-solid-geometry tree. See the module docs.
#[derive(Debug, Clone, PartialEq)]
pub enum Csg {
    /// A single primitive.
    Primitive(Primitive),
    /// Points inside any operand (at least one operand).
    Union(Vec<Csg>),
    /// Points inside every operand (at least one operand).
    Intersection(Vec<Csg>),
    /// Points inside the first operand and not inside the second.
    Difference(Box<Csg>, Box<Csg>),
}

impl Csg {
    /// The axis-aligned box `[lo, hi]`, m.
    pub fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> Self {
        Self::Primitive(Primitive::Cuboid { lo, hi })
    }

    /// A capped cylinder from the centre of its first end cap along `axis`,
    /// m.
    pub fn cylinder(base: [f64; 3], axis: [f64; 3], radius: f64, height: f64) -> Self {
        Self::Primitive(Primitive::Cylinder {
            base,
            axis,
            radius,
            height,
        })
    }

    /// The half-space `{p : normal . p <= offset}`.
    pub fn half_space(normal: [f64; 3], offset: f64) -> Self {
        Self::Primitive(Primitive::HalfSpace { normal, offset })
    }

    /// `self` united with `other`.
    pub fn union(self, other: Csg) -> Self {
        Self::Union(vec![self, other])
    }

    /// `self` intersected with `other`.
    pub fn intersect(self, other: Csg) -> Self {
        Self::Intersection(vec![self, other])
    }

    /// `self` with `other` removed.
    pub fn minus(self, other: Csg) -> Self {
        Self::Difference(Box::new(self), Box::new(other))
    }

    /// Check every primitive and normalise cylinder axes and half-space
    /// normals in place.
    fn validate(&mut self) -> Result<(), &'static str> {
        let finite = |v: &V3| v.iter().all(|c| c.is_finite());
        match self {
            Self::Primitive(Primitive::Cuboid { lo, hi }) => {
                if !finite(lo) || !finite(hi) {
                    return Err("a box corner is not finite");
                }
                if (0..3).any(|k| lo[k] >= hi[k]) {
                    return Err("a box needs lo < hi on every axis");
                }
            }
            Self::Primitive(Primitive::Cylinder {
                base,
                axis,
                radius,
                height,
            }) => {
                if !finite(base) {
                    return Err("a cylinder base is not finite");
                }
                let l = norm(*axis);
                if !finite(axis) || !(l > 0.0 && l.is_finite()) {
                    return Err("a cylinder axis must be finite and non-zero");
                }
                if !(radius.is_finite() && *radius > 0.0) {
                    return Err("a cylinder radius must be finite and positive");
                }
                if !(height.is_finite() && *height > 0.0) {
                    return Err("a cylinder height must be finite and positive");
                }
                *axis = scale(*axis, 1.0 / l);
            }
            Self::Primitive(Primitive::HalfSpace { normal, offset }) => {
                let l = norm(*normal);
                if !finite(normal) || !(l > 0.0 && l.is_finite()) {
                    return Err("a half-space normal must be finite and non-zero");
                }
                if !offset.is_finite() {
                    return Err("a half-space offset is not finite");
                }
                // Divide, not multiply by 1/l: x / sqrt(x*x) is exactly ±1 for an axis normal.
                *normal = [normal[0] / l, normal[1] / l, normal[2] / l];
                *offset /= l;
            }
            Self::Union(v) | Self::Intersection(v) => {
                if v.is_empty() {
                    return Err("a union or intersection needs at least one operand");
                }
                for c in v {
                    c.validate()?;
                }
            }
            Self::Difference(a, b) => {
                a.validate()?;
                b.validate()?;
            }
        }
        Ok(())
    }

    /// Conservative axis-aligned bounds `(lo, hi)` (module docs), possibly
    /// infinite. Requires a validated tree.
    fn bounds(&self) -> (V3, V3) {
        match self {
            Self::Primitive(Primitive::Cuboid { lo, hi }) => (*lo, *hi),
            Self::Primitive(Primitive::Cylinder {
                base,
                axis,
                radius,
                height,
            }) => {
                let top = [
                    base[0] + height * axis[0],
                    base[1] + height * axis[1],
                    base[2] + height * axis[2],
                ];
                let mut lo = [0.0; 3];
                let mut hi = [0.0; 3];
                for k in 0..3 {
                    let e = radius * (1.0 - axis[k] * axis[k]).max(0.0).sqrt();
                    lo[k] = base[k].min(top[k]) - e;
                    hi[k] = base[k].max(top[k]) + e;
                }
                (lo, hi)
            }
            Self::Primitive(Primitive::HalfSpace { normal, offset }) => {
                let mut lo = [f64::NEG_INFINITY; 3];
                let mut hi = [f64::INFINITY; 3];
                for k in 0..3 {
                    if normal[k] == 1.0 {
                        hi[k] = *offset;
                    } else if normal[k] == -1.0 {
                        lo[k] = -*offset;
                    }
                }
                (lo, hi)
            }
            Self::Union(v) => v.iter().map(Csg::bounds).fold(
                ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
                |(lo, hi), (l, h)| {
                    (
                        [lo[0].min(l[0]), lo[1].min(l[1]), lo[2].min(l[2])],
                        [hi[0].max(h[0]), hi[1].max(h[1]), hi[2].max(h[2])],
                    )
                },
            ),
            Self::Intersection(v) => v.iter().map(Csg::bounds).fold(
                ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]),
                |(lo, hi), (l, h)| {
                    (
                        [lo[0].max(l[0]), lo[1].max(l[1]), lo[2].max(l[2])],
                        [hi[0].min(h[0]), hi[1].min(h[1]), hi[2].min(h[2])],
                    )
                },
            ),
            Self::Difference(a, _) => a.bounds(),
        }
    }

    /// Whether `p` is inside, with rule-2 ownership of points exactly on a
    /// surface. Requires a validated tree.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        match self {
            Self::Primitive(prim) => prim.contains(p),
            Self::Union(v) => v.iter().any(|c| c.contains(p)),
            Self::Intersection(v) => v.iter().all(|c| c.contains(p)),
            Self::Difference(a, b) => a.contains(p) && !b.contains(p),
        }
    }

    /// The inside spans of the line `o + t d` (`t` over all reals), sorted,
    /// with the rule-4 clean-up at tolerance `tol`. Requires a validated tree.
    fn spans(&self, o: V3, d: V3, tol: f64) -> Vec<Span> {
        let out = match self {
            Self::Primitive(prim) => prim.span(o, d).into_iter().collect(),
            Self::Union(v) => {
                let mut it = v.iter().map(|c| c.spans(o, d, tol));
                let first = it.next().unwrap_or_default();
                it.fold(first, |a, b| union(&a, &b))
            }
            Self::Intersection(v) => {
                let mut it = v.iter().map(|c| c.spans(o, d, tol));
                let first = it.next().unwrap_or_default();
                it.fold(first, |a, b| intersect(&a, &b))
            }
            Self::Difference(a, b) => {
                intersect(&a.spans(o, d, tol), &complement(&b.spans(o, d, tol)))
            }
        };
        clean(out, tol)
    }
}

impl Primitive {
    fn contains(&self, p: V3) -> bool {
        match *self {
            Self::Cuboid { lo, hi } => (0..3).all(|k| in_slab(p[k], lo[k], hi[k], axis(k))),
            Self::Cylinder {
                base,
                axis: a,
                radius,
                height,
            } => {
                let rel = sub(p, base);
                let z = dot(rel, a);
                if !in_slab(z, 0.0, height, a) {
                    return false;
                }
                let w = sub(rel, scale(a, z));
                let c = dot(w, w) - radius * radius;
                c < 0.0 || (c == 0.0 && owns(w))
            }
            Self::HalfSpace { normal, offset } => {
                let s = dot(normal, p) - offset;
                s < 0.0 || (s == 0.0 && owns(normal))
            }
        }
    }

    /// The single inside span of the line `o + t d`, if any.
    fn span(&self, o: V3, d: V3) -> Option<Span> {
        match *self {
            Self::Cuboid { lo, hi } => (0..3).try_fold(Span::ALL, |s, k| {
                clip(s, slab(o[k], d[k], lo[k], hi[k], axis(k))?)
            }),
            Self::Cylinder {
                base,
                axis: a,
                radius,
                height,
            } => {
                let rel = sub(o, base);
                let (z0, dz) = (dot(rel, a), dot(d, a));
                let caps = slab(z0, dz, 0.0, height, a)?;
                let w = sub(rel, scale(a, z0));
                let v = sub(d, scale(a, dz));
                let qa = dot(v, v);
                let qb = dot(w, v);
                let qc = dot(w, w) - radius * radius;
                let side = if qa == 0.0 {
                    // Parallel to the axis: inside the tube or not, all along.
                    if qc < 0.0 || (qc == 0.0 && owns(w)) {
                        Span::ALL
                    } else {
                        return None;
                    }
                } else {
                    let disc = qb * qb - qa * qc;
                    if disc.is_nan() || disc <= 0.0 {
                        return None;
                    }
                    // Numerical Recipes (3rd ed.) section 5.6.
                    let q = -(qb + disc.sqrt().copysign(qb));
                    let (r1, r2) = (q / qa, qc / q);
                    let (t0, t1) = if r1 <= r2 { (r1, r2) } else { (r2, r1) };
                    let radial = |t: f64| {
                        let n = [w[0] + t * v[0], w[1] + t * v[1], w[2] + t * v[2]];
                        scale(n, 1.0 / norm(n))
                    };
                    Span {
                        t0,
                        n0: radial(t0),
                        t1,
                        n1: radial(t1),
                    }
                };
                clip(caps, side)
            }
            Self::HalfSpace { normal, offset } => {
                let s0 = dot(normal, o) - offset;
                let ds = dot(normal, d);
                if ds == 0.0 {
                    (s0 < 0.0 || (s0 == 0.0 && owns(normal))).then_some(Span::ALL)
                } else {
                    let t = -s0 / ds;
                    Some(if ds > 0.0 {
                        Span {
                            t0: f64::NEG_INFINITY,
                            n0: [0.0; 3],
                            t1: t,
                            n1: normal,
                        }
                    } else {
                        Span {
                            t0: t,
                            n0: normal,
                            t1: f64::INFINITY,
                            n1: [0.0; 3],
                        }
                    })
                }
            }
        }
    }
}

/// An inside span `[t0, t1]` of a line, with the outward unit normals of the
/// surfaces at its ends (unused at an infinite end).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Span {
    t0: f64,
    n0: V3,
    t1: f64,
    n1: V3,
}

impl Span {
    const ALL: Span = Span {
        t0: f64::NEG_INFINITY,
        n0: [0.0; 3],
        t1: f64::INFINITY,
        n1: [0.0; 3],
    };
}

/// The span of the slab `[lo, hi]` along the unit vector `u`, for a line at
/// coordinate `x0` moving at rate `dx` (Kay and Kajiya 1986).
fn slab(x0: f64, dx: f64, lo: f64, hi: f64, u: V3) -> Option<Span> {
    if dx == 0.0 {
        return in_slab(x0, lo, hi, u).then_some(Span::ALL);
    }
    let (ta, tb) = ((lo - x0) / dx, (hi - x0) / dx);
    Some(if dx > 0.0 {
        Span {
            t0: ta,
            n0: neg(u),
            t1: tb,
            n1: u,
        }
    } else {
        Span {
            t0: tb,
            n0: u,
            t1: ta,
            n1: neg(u),
        }
    })
}

/// The overlap of two spans, if not empty. Ties in either end go to `a`.
fn clip(a: Span, b: Span) -> Option<Span> {
    let (t0, n0) = if b.t0 > a.t0 {
        (b.t0, b.n0)
    } else {
        (a.t0, a.n0)
    };
    let (t1, n1) = if b.t1 < a.t1 {
        (b.t1, b.n1)
    } else {
        (a.t1, a.n1)
    };
    (t0 < t1).then_some(Span { t0, n0, t1, n1 })
}

/// The intersection of two sorted, disjoint span lists.
fn intersect(a: &[Span], b: &[Span]) -> Vec<Span> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if let Some(s) = clip(a[i], b[j]) {
            out.push(s);
        }
        if a[i].t1 < b[j].t1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}

/// The complement of a sorted, disjoint span list; the normals are negated so
/// they stay outward.
fn complement(a: &[Span]) -> Vec<Span> {
    let mut out = Vec::with_capacity(a.len() + 1);
    let (mut t, mut n) = (f64::NEG_INFINITY, [0.0; 3]);
    for s in a {
        if s.t0 > t {
            out.push(Span {
                t0: t,
                n0: n,
                t1: s.t0,
                n1: neg(s.n0),
            });
        }
        (t, n) = (s.t1, neg(s.n1));
    }
    if t < f64::INFINITY {
        out.push(Span {
            t0: t,
            n0: n,
            t1: f64::INFINITY,
            n1: [0.0; 3],
        });
    }
    out
}

/// The union of two sorted, disjoint span lists (De Morgan).
fn union(a: &[Span], b: &[Span]) -> Vec<Span> {
    complement(&intersect(&complement(a), &complement(b)))
}

/// Rule 4: drop spans no longer than `tol`, then close gaps no longer than
/// `tol`.
fn clean(spans: Vec<Span>, tol: f64) -> Vec<Span> {
    let mut out: Vec<Span> = Vec::with_capacity(spans.len());
    for s in spans.into_iter().filter(|s| s.t1 - s.t0 > tol) {
        match out.last_mut() {
            Some(last) if s.t0 - last.t1 <= tol => {
                last.t1 = s.t1;
                last.n1 = s.n1;
            }
            _ => out.push(s),
        }
    }
    out
}

/// A target of CSG solids, each tagged with a material. See the module docs
/// for the ownership rules.
#[derive(Debug, Clone, PartialEq)]
pub struct CsgGeometry {
    materials: Vec<Material>,
    solids: Vec<(Csg, u32)>,
    tol: f64,
    lo: V3,
    hi: V3,
}

impl CsgGeometry {
    /// Build from the material list and `(solid, material index)` pairs, in
    /// priority order (the first listed owns overlaps). Lengths in metres.
    /// Cylinder axes and half-space normals are normalised.
    pub fn new(materials: Vec<Material>, solids: Vec<(Csg, u32)>) -> Result<Self, GeometryError> {
        if materials.is_empty() {
            return Err(GeometryError::NoMaterials);
        }
        if solids.is_empty() {
            return Err(GeometryError::NoSolids);
        }
        let mut solids = solids;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (i, (csg, m)) in solids.iter_mut().enumerate() {
            if *m as usize >= materials.len() {
                return Err(GeometryError::SolidMaterialOutOfRange {
                    solid: i,
                    index: *m,
                    n_materials: materials.len(),
                });
            }
            csg.validate()
                .map_err(|reason| GeometryError::CsgInvalid { solid: i, reason })?;
            let (l, h) = csg.bounds();
            if !l.iter().chain(&h).all(|c| c.is_finite()) {
                return Err(GeometryError::CsgUnbounded { solid: i });
            }
            if (0..3).any(|k| l[k] > h[k]) {
                return Err(GeometryError::CsgInvalid {
                    solid: i,
                    reason: "the solid's bounding box is empty",
                });
            }
            for k in 0..3 {
                lo[k] = lo[k].min(l[k]);
                hi[k] = hi[k].max(h[k]);
            }
        }
        let tol = TOL_REL * norm(sub(hi, lo));
        Ok(Self {
            materials,
            solids,
            tol,
            lo,
            hi,
        })
    }

    /// The solid of `region` (validated: unit cylinder axes and half-space
    /// normals).
    pub fn solid(&self, region: usize) -> &Csg {
        &self.solids[region].0
    }

    /// The geometric tolerance of the ownership rules, m.
    pub fn tolerance_m(&self) -> f64 {
        self.tol
    }

    /// Axis-aligned bounds of all solids, `(min, max)`, m.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        (self.lo, self.hi)
    }

    /// [`Geometry::locate`] with the offset of rule 2 applied with the given
    /// sign on each axis (all +1 for `locate`).
    fn locate_shifted(&self, pos: V3, sign: V3) -> Option<usize> {
        if !pos.iter().all(|c| c.is_finite()) {
            return None;
        }
        let s = self.tol / 3f64.sqrt();
        let p = [
            pos[0] + sign[0] * s,
            pos[1] + sign[1] * s,
            pos[2] + sign[2] * s,
        ];
        self.solids.iter().position(|(c, _)| c.contains(p))
    }
}

impl Geometry for CsgGeometry {
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
        self.solids[region].1 as usize
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
        let n = self.solids.len();

        // Span ends (t, solid, entering, normal) in (tol, limit + tol], and
        // the state of every solid at t = tol (rules 2 and 3). The flight's
        // own region is taken to be inside, as for meshes.
        let mut events: Vec<(f64, u32, bool, V3)> = Vec::new();
        let mut inside = vec![false; n];
        for (j, (csg, _)) in self.solids.iter().enumerate() {
            for s in csg.spans(pos, dir, tol) {
                if s.t0 <= tol && s.t1 > tol {
                    inside[j] = true;
                }
                for (t, entering, normal) in [(s.t0, true, s.n0), (s.t1, false, s.n1)] {
                    if t > tol && t <= limit + tol {
                        events.push((t, j as u32, entering, normal));
                    }
                }
            }
        }
        inside[region] = true;
        events.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));

        let owner_of = |inside: &[bool]| inside.iter().position(|&b| b);
        let mut current = region;
        // Start-state mismatch (a surface start heading into a solid that has
        // priority or is not `region`): an event at distance 0.
        let mut pending_start = true;
        let mut i = 0usize;
        loop {
            let (t0, cluster_end) = if pending_start {
                (0.0, i)
            } else {
                if i >= events.len() || events[i].0 > limit {
                    return Flight::Clear { region: current };
                }
                let t0 = events[i].0;
                let mut e = i;
                while e < events.len() && events[e].0 - t0 <= tol {
                    inside[events[e].1 as usize] = events[e].2;
                    e += 1;
                }
                (t0, e)
            };
            let cluster = &events[i..cluster_end];
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
                    if self.solids[k].1 == self.solids[current].1 {
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
                    // The surface of the solid being left, else any.
                    let normal = cluster
                        .iter()
                        .find(|e| e.1 as usize == current && !e.2)
                        .or_else(|| cluster.first())
                        .map_or(neg(dir), |e| e.3);
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
    /// eight diagonal offsets (the one of [`Geometry::locate`] and its mirror
    /// images) would put it outside.
    fn in_vacuum(&self, pos: [f64; 3]) -> bool {
        (0..8).all(|m| {
            let sign = |k: usize| if m >> k & 1 == 0 { 1.0 } else { -1.0 };
            self.locate_shifted(pos, [sign(0), sign(1), sign(2)])
                .is_none()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_util::{c, enter, si, unit};
    use super::*;

    const TOL: f64 = 1e-12;

    fn spans(c: &Csg, o: V3, d: V3) -> Vec<(f64, f64)> {
        let mut c = c.clone();
        c.validate().unwrap();
        c.spans(o, d, TOL).iter().map(|s| (s.t0, s.t1)).collect()
    }

    fn close(a: f64, b: f64) -> bool {
        a == b || (a - b).abs() <= 1e-12 * (1.0 + a.abs().max(b.abs()))
    }

    fn same(got: &[(f64, f64)], want: &[(f64, f64)]) -> bool {
        got.len() == want.len()
            && got
                .iter()
                .zip(want)
                .all(|(g, w)| close(g.0, w.0) && close(g.1, w.1))
    }

    const INF: f64 = f64::INFINITY;

    #[test]
    fn box_span_and_containment() {
        let b = Csg::cuboid([0.0; 3], [2.0, 1.0, 1.0]);
        let s = spans(&b, [-1.0, 0.5, 0.5], [1.0, 0.0, 0.0]);
        assert!(same(&s, &[(1.0, 3.0)]), "{s:?}");
        let mut b2 = b.clone();
        b2.validate().unwrap();
        let sp = b2.spans([-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], TOL);
        assert_eq!(sp[0].n0, [-1.0, 0.0, 0.0]);
        assert_eq!(sp[0].n1, [1.0, 0.0, 0.0]);
        // Oblique.
        let d = unit([1.0, 1.0, 0.0]);
        let s = spans(&b, [0.0, 0.0, 0.5], d);
        assert!(same(&s, &[(0.0, 2f64.sqrt())]), "{s:?}");
        // Misses, and parallel outside a slab.
        assert!(spans(&b, [-1.0, 2.0, 0.5], [1.0, 0.0, 0.0]).is_empty());
        // Parallel inside: unbounded in that direction only via other axes.
        let s = spans(&b, [1.0, 0.5, 0.5], [0.0, 0.0, 1.0]);
        assert!(same(&s, &[(-0.5, 0.5)]), "{s:?}");
        // A ray lying in the lower face is inside, in the upper face outside.
        assert_eq!(spans(&b, [-1.0, 0.0, 0.5], [1.0, 0.0, 0.0]).len(), 1);
        assert!(spans(&b, [-1.0, 1.0, 0.5], [1.0, 0.0, 0.0]).is_empty());
        // Containment: half-open, lower faces own.
        assert!(b.contains([0.0, 0.0, 0.0]));
        assert!(b.contains([1.0, 0.5, 0.5]));
        assert!(!b.contains([2.0, 0.5, 0.5]));
        assert!(!b.contains([1.0, 1.0, 0.5]));
    }

    #[test]
    fn cylinder_span_caps_and_tangent() {
        // Axis z, radius 1, from z = 0 to 2.
        let cyl = Csg::cylinder([0.0; 3], [0.0, 0.0, 5.0], 1.0, 2.0);
        // Across the axis.
        let s = spans(&cyl, [-3.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        assert!(same(&s, &[(2.0, 4.0)]), "{s:?}");
        // Off-centre chord 2 sqrt(1 - b^2).
        let b = 0.6;
        let s = spans(&cyl, [-3.0, b, 1.0], [1.0, 0.0, 0.0]);
        let h = (1.0 - b * b).sqrt();
        assert!(same(&s, &[(3.0 - h, 3.0 + h)]), "{s:?}");
        // Along the axis: through both caps.
        let s = spans(&cyl, [0.5, 0.0, -1.0], [0.0, 0.0, 1.0]);
        assert!(same(&s, &[(1.0, 3.0)]), "{s:?}");
        let mut v = cyl.clone();
        v.validate().unwrap();
        let sp = v.spans([0.5, 0.0, -1.0], [0.0, 0.0, 1.0], TOL);
        assert_eq!((sp[0].n0, sp[0].n1), ([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]));
        // Oblique, through a cap and the side: enters the bottom cap at
        // (0, 0, 0), leaves the side at x = 1, z = 1.
        let d = unit([1.0, 0.0, 1.0]);
        let s = spans(&cyl, [-1.0, 0.0, -1.0], d);
        assert!(same(&s, &[(2f64.sqrt(), 2.0 * 2f64.sqrt())]), "{s:?}");
        // Tangent: no span.
        assert!(spans(&cyl, [-3.0, 1.0, 1.0], [1.0, 0.0, 0.0]).is_empty());
        // Parallel to the axis outside: none.
        assert!(spans(&cyl, [1.5, 0.0, -1.0], [0.0, 0.0, 1.0]).is_empty());
        // Containment (on the validated, unit-axis tree) and a tilted axis.
        let mut cyl = cyl;
        cyl.validate().unwrap();
        assert!(cyl.contains([0.0, 0.0, 0.0]));
        assert!(cyl.contains([0.5, 0.5, 1.0]));
        assert!(!cyl.contains([0.0, 0.0, 2.0]));
        assert!(!cyl.contains([0.8, 0.8, 1.0]));
        let tilted = Csg::cylinder([0.0; 3], [1.0, 1.0, 0.0], 0.5, 3.0);
        let s = spans(&tilted, [0.0, 0.0, 0.0], unit([1.0, 1.0, 0.0]));
        assert!(same(&s, &[(0.0, 3.0)]), "{s:?}");
    }

    #[test]
    fn half_space_span_and_parallel() {
        // x <= 2.
        let h = Csg::half_space([2.0, 0.0, 0.0], 4.0);
        assert!(same(&spans(&h, [0.0; 3], [1.0, 0.0, 0.0]), &[(-INF, 2.0)]));
        assert!(same(&spans(&h, [0.0; 3], [-1.0, 0.0, 0.0]), &[(-2.0, INF)]));
        // Parallel to the plane: all or nothing.
        assert!(same(
            &spans(&h, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            &[(-INF, INF)]
        ));
        assert!(spans(&h, [3.0, 0.0, 0.0], [0.0, 1.0, 0.0]).is_empty());
        // In the plane: its outward normal +x points along the diagonal, so
        // the plane is outside.
        assert!(spans(&h, [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]).is_empty());
        assert!(!h.contains([2.0, 0.0, 0.0]));
        assert!(h.contains([1.999, 5.0, 5.0]));
        // A plane containing the diagonal: exactly one side owns it.
        let a = Csg::half_space([1.0, -1.0, 0.0], 0.0);
        let b = Csg::half_space([-1.0, 1.0, 0.0], 0.0);
        let p = [1.0, 1.0, 7.0];
        assert_ne!(a.contains(p), b.contains(p));
    }

    #[test]
    fn boolean_spans() {
        let x = [1.0, 0.0, 0.0];
        let o = [-10.0, 0.5, 0.5];
        let a = Csg::cuboid([0.0; 3], [2.0, 1.0, 1.0]);
        let b = Csg::cuboid([1.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        let far = Csg::cuboid([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]);
        // Overlapping union: one span; disjoint union: two.
        let s = spans(&a.clone().union(b.clone()), o, x);
        assert!(same(&s, &[(10.0, 13.0)]), "{s:?}");
        let s = spans(&a.clone().union(far.clone()), o, x);
        assert!(same(&s, &[(10.0, 12.0), (15.0, 16.0)]), "{s:?}");
        // Abutting: merged (rule 4).
        let ab = Csg::cuboid([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        let s = spans(&a.clone().union(ab), o, x);
        assert!(same(&s, &[(10.0, 13.0)]), "{s:?}");
        // Intersection and difference.
        let s = spans(&a.clone().intersect(b.clone()), o, x);
        assert!(same(&s, &[(11.0, 12.0)]), "{s:?}");
        let s = spans(&a.clone().minus(b.clone()), o, x);
        assert!(same(&s, &[(10.0, 11.0)]), "{s:?}");
        let s = spans(&b.clone().minus(a.clone()), o, x);
        assert!(same(&s, &[(12.0, 13.0)]), "{s:?}");
        // A hole in the middle splits the span; its normals are inward to
        // the hole.
        let hole = Csg::cuboid([0.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mut d = a.clone().minus(hole);
        d.validate().unwrap();
        let sp = d.spans(o, x, TOL);
        assert_eq!(sp.len(), 2);
        assert!(close(sp[0].t1, 10.5) && close(sp[1].t0, 11.0));
        assert_eq!(sp[0].n1, [1.0, 0.0, 0.0]);
        assert_eq!(sp[1].n0, [-1.0, 0.0, 0.0]);
        // Nested: (A u B) - C with C in the middle of the union.
        let cc = Csg::cuboid([1.5, 0.0, 0.0], [2.5, 1.0, 1.0]);
        let s = spans(&a.clone().union(b.clone()).minus(cc), o, x);
        assert!(same(&s, &[(10.0, 11.5), (12.5, 13.0)]), "{s:?}");
        // A flush cut leaves no sliver (rule 4).
        let flush = Csg::cuboid([1.0 - 1e-14, 0.0, 0.0], [5.0, 1.0, 1.0]);
        let s = spans(&b.clone().minus(flush), o, x);
        assert!(s.is_empty(), "{s:?}");
        // Six axis-aligned half-spaces equal the box.
        let six = Csg::Intersection(vec![
            Csg::half_space([-1.0, 0.0, 0.0], 0.0),
            Csg::half_space([1.0, 0.0, 0.0], 2.0),
            Csg::half_space([0.0, -1.0, 0.0], 0.0),
            Csg::half_space([0.0, 1.0, 0.0], 1.0),
            Csg::half_space([0.0, 0.0, -1.0], 0.0),
            Csg::half_space([0.0, 0.0, 1.0], 1.0),
        ]);
        let d = unit([0.3, 0.2, 0.1]);
        assert!(same(
            &spans(&six, [0.1, 0.1, 0.1], d),
            &spans(&a, [0.1, 0.1, 0.1], d)
        ));
        let mut v = six.clone();
        v.validate().unwrap();
        assert_eq!(v.bounds(), ([0.0; 3], [2.0, 1.0, 1.0]));
    }

    #[test]
    fn non_unit_axis_normals_stay_axis_aligned() {
        // 49 * (1/49) != 1 in f64, so normalising must divide, not scale.
        for k in [49.0, 98.0, 103.0, 107.0, 161.0, 187.0] {
            let six = Csg::Intersection(vec![
                Csg::half_space([-k, 0.0, 0.0], 0.0),
                Csg::half_space([k, 0.0, 0.0], 2.0 * k),
                Csg::half_space([0.0, -k, 0.0], 0.0),
                Csg::half_space([0.0, k, 0.0], k),
                Csg::half_space([0.0, 0.0, -k], 0.0),
                Csg::half_space([0.0, 0.0, k], k),
            ]);
            let mut v = six.clone();
            v.validate().unwrap();
            assert_eq!(v.bounds(), ([0.0; 3], [2.0, 1.0, 1.0]), "k = {k}");
        }
    }

    #[test]
    fn validation_and_bounds() {
        let ok = |c: Csg| CsgGeometry::new(vec![si()], vec![(c, 0)]);
        assert!(matches!(
            ok(Csg::cuboid([0.0; 3], [0.0, 1.0, 1.0])),
            Err(GeometryError::CsgInvalid { solid: 0, .. })
        ));
        assert!(matches!(
            ok(Csg::cylinder([0.0; 3], [0.0; 3], 1.0, 1.0)),
            Err(GeometryError::CsgInvalid { .. })
        ));
        assert!(matches!(
            ok(Csg::cylinder([0.0; 3], [0.0, 0.0, 1.0], -1.0, 1.0)),
            Err(GeometryError::CsgInvalid { .. })
        ));
        assert!(matches!(
            ok(Csg::Union(vec![])),
            Err(GeometryError::CsgInvalid { .. })
        ));
        assert_eq!(
            ok(Csg::half_space([1.0, 0.0, 0.0], 1.0)),
            Err(GeometryError::CsgUnbounded { solid: 0 })
        );
        assert_eq!(
            ok(Csg::cuboid([0.0; 3], [1.0; 3]).union(Csg::half_space([1.0, 1.0, 0.0], 1.0))),
            Err(GeometryError::CsgUnbounded { solid: 0 })
        );
        assert!(matches!(
            ok(Csg::cuboid([0.0; 3], [1.0; 3]).intersect(Csg::cuboid([2.0; 3], [3.0; 3]))),
            Err(GeometryError::CsgInvalid { .. })
        ));
        // A box cut by an oblique plane is bounded by the box.
        let g = ok(Csg::cuboid([0.0; 3], [1.0; 3]).minus(Csg::half_space([1.0, 1.0, 1.0], 0.5)))
            .unwrap();
        assert_eq!(g.bounds(), ([0.0; 3], [1.0; 3]));
        // A tilted cylinder's exact box.
        let s = 0.5f64.sqrt();
        let g = ok(Csg::cylinder([0.0; 3], [1.0, 1.0, 0.0], 1.0, 2.0)).unwrap();
        let (lo, hi) = g.bounds();
        let want_lo = [-s, -s, -1.0];
        let want_hi = [2.0 * s + s, 2.0 * s + s, 1.0];
        for k in 0..3 {
            assert!(
                close(lo[k], want_lo[k]) && close(hi[k], want_hi[k]),
                "{lo:?} {hi:?}"
            );
        }
        // Validation normalises the half-space.
        let g =
            ok(Csg::cuboid([0.0; 3], [4.0; 3]).intersect(Csg::half_space([2.0, 0.0, 0.0], 6.0)))
                .unwrap();
        assert_eq!(g.bounds(), ([0.0; 3], [3.0, 4.0, 4.0]));
        assert!(CsgGeometry::new(vec![], vec![(Csg::cuboid([0.0; 3], [1.0; 3]), 0)]).is_err());
        assert_eq!(
            CsgGeometry::new(vec![si()], vec![]),
            Err(GeometryError::NoSolids)
        );
        assert_eq!(
            CsgGeometry::new(vec![si()], vec![(Csg::cuboid([0.0; 3], [1.0; 3]), 1)]),
            Err(GeometryError::SolidMaterialOutOfRange {
                solid: 0,
                index: 1,
                n_materials: 1
            })
        );
    }

    #[test]
    fn mask_listed_first_cuts_into_substrate() {
        // A mask block over part of a substrate, listed first.
        let g = CsgGeometry::new(
            vec![si(), c()],
            vec![
                (Csg::cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 1),
                (Csg::cuboid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]), 0),
            ],
        )
        .unwrap();
        assert_eq!(g.locate([0.5, 0.5, 0.5]), Some(0));
        assert_eq!(g.locate([2.0, 0.5, 0.5]), Some(1));
        let e = g.exit(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
        assert!(close(e.distance, 0.5));
        assert_eq!(enter(&e).0, 1);
    }
}
