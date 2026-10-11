//! A bounding-volume hierarchy over triangles and an exact ray-triangle test.
//!
//! # Construction
//!
//! Binned surface-area-heuristic build: the SAH is from J. D. MacDonald and
//! K. S. Booth, "Heuristics for ray tracing using space subdivision", The
//! Visual Computer 6 (1990) 153-166, with the binned evaluation of I. Wald,
//! "On fast construction of SAH-based bounding volume hierarchies", Proc. IEEE
//! Symposium on Interactive Ray Tracing (2007) 33-40. Each node tries 16 bins
//! on each axis of the triangle centroids and takes the cheapest plane (ties go
//! to the lowest axis, then the lowest bin). A node of at most
//! [`MAX_LEAF`] triangles becomes a leaf when that is no dearer than the
//! best split. Below depth [`SAH_DEPTH`] (and whenever the centroids
//! coincide) the node is split at the median instead, which bounds the tree
//! depth by `SAH_DEPTH + log2(n)`.
//!
//! The build is single threaded and uses only ordered containers and stable
//! partitions, so the tree depends on the input triangles alone, never on a
//! thread count or a hash seed.
//!
//! # Layout and traversal
//!
//! A flat depth-first node array: the first child of an interior node is the
//! next node, the second is stored explicitly; triangles are copied into leaf
//! order. Traversal keeps an explicit stack and visits every node whose box
//! the ray segment overlaps (a hit list, not only the nearest hit, because
//! the mesh queries need every crossing). Node boxes are padded by an absolute
//! amount, so rounding in the slab test never drops a triangle that touches a
//! box face.
//!
//! # Ray-triangle test
//!
//! The watertight algorithm of S. Woop, C. Benthin and I. Wald, "Watertight
//! ray/triangle intersection", J. Comput. Graph. Tech. 2 (2013) 65-82: the
//! triangle is translated to the ray origin and sheared so the ray runs along
//! `+z`, and the three signed edge functions are evaluated from the sheared
//! vertices. An edge shared by two triangles is evaluated from the same two
//! vertices (swapped), so the two results are exact negatives and a ray
//! cannot pass between them. That paper accepts an edge function of exactly
//! zero from both neighbours (a ray through an edge may hit twice). Here a
//! zero is resolved by the top-left tie-break of rasterisation (see
//! [`edge_owns`]), which is the sign a fixed infinitesimal shift of the ray
//! would give. With that shift every ray is generic, so a ray through a shared
//! edge or vertex crosses a closed surface exactly as often as a ray beside
//! it. The treatment is exact in the sense that the decisions are consistent
//! between neighbouring triangles; it assumes the edge functions of
//! neighbours are computed from identical vertex coordinates, which welded
//! meshes guarantee.

pub(super) type V3 = [f64; 3];

/// Largest number of triangles in a leaf.
const MAX_LEAF: usize = 4;
/// Bins per axis in the SAH sweep.
const BINS: usize = 16;
/// Depth below which the SAH is used; deeper nodes split at the median.
const SAH_DEPTH: usize = 32;
/// Traversal stack size: depth at most `SAH_DEPTH + 32` and one pending
/// sibling per level.
const STACK: usize = 96;

#[derive(Debug, Clone, PartialEq)]
struct Node {
    min: V3,
    max: V3,
    /// Leaf: index of the first triangle. Interior: index of the second child.
    first: u32,
    /// Triangles in a leaf; 0 for an interior node.
    count: u32,
}

/// A BVH over a fixed triangle list.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Bvh {
    nodes: Vec<Node>,
    /// Triangles in leaf order.
    tris: Vec<[V3; 3]>,
}

struct Prim {
    min: V3,
    max: V3,
    centroid: V3,
}

fn union(a: (V3, V3), b: (V3, V3)) -> (V3, V3) {
    let mut r = a;
    for k in 0..3 {
        r.0[k] = r.0[k].min(b.0[k]);
        r.1[k] = r.1[k].max(b.1[k]);
    }
    r
}

const EMPTY: (V3, V3) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);

fn area(b: (V3, V3)) -> f64 {
    let d = [b.1[0] - b.0[0], b.1[1] - b.0[1], b.1[2] - b.0[2]];
    if d[0] < 0.0 {
        return 0.0;
    }
    2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
}

impl Bvh {
    /// Build over `tris`. `pad` (>= 0, metres) inflates every node box.
    pub(super) fn build(tris: &[[V3; 3]], pad: f64) -> Self {
        let prims: Vec<Prim> = tris
            .iter()
            .map(|t| {
                let mut min = t[0];
                let mut max = t[0];
                for v in &t[1..] {
                    for k in 0..3 {
                        min[k] = min[k].min(v[k]);
                        max[k] = max[k].max(v[k]);
                    }
                }
                let centroid = [
                    0.5 * (min[0] + max[0]),
                    0.5 * (min[1] + max[1]),
                    0.5 * (min[2] + max[2]),
                ];
                Prim { min, max, centroid }
            })
            .collect();
        let mut ids: Vec<u32> = (0..tris.len() as u32).collect();
        let mut nodes = Vec::with_capacity(2 * tris.len().max(1));
        if !tris.is_empty() {
            build_rec(&mut nodes, &prims, &mut ids, 0, 0, pad);
        }
        let tris = ids.iter().map(|&i| tris[i as usize]).collect();
        Self { nodes, tris }
    }

    /// The bounds of everything (padded), or `None` for an empty tree.
    pub(super) fn bounds(&self) -> Option<(V3, V3)> {
        self.nodes.first().map(|n| (n.min, n.max))
    }

    /// Triangle `i` (leaf order, as passed to the visitor).
    pub(super) fn triangle(&self, i: u32) -> &[V3; 3] {
        &self.tris[i as usize]
    }

    /// Call `f(t, triangle)` for every triangle crossed by `ray` at a
    /// parameter `t` with `tmin < t <= tmax`, in traversal order.
    pub(super) fn for_each_hit(
        &self,
        ray: &Ray,
        tmin: f64,
        tmax: f64,
        mut f: impl FnMut(f64, u32),
    ) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = [0u32; STACK];
        let mut sp = 0usize;
        let mut cur = 0usize;
        loop {
            let node = &self.nodes[cur];
            if ray.hits_box(node, tmin, tmax) {
                if node.count == 0 {
                    stack[sp] = node.first;
                    sp += 1;
                    cur += 1;
                    continue;
                }
                let first = node.first as usize;
                for i in first..first + node.count as usize {
                    if let Some(t) = ray.intersect(&self.tris[i]) {
                        if t > tmin && t <= tmax {
                            f(t, i as u32);
                        }
                    }
                }
            }
            if sp == 0 {
                return;
            }
            sp -= 1;
            cur = stack[sp] as usize;
        }
    }
}

fn build_rec(
    nodes: &mut Vec<Node>,
    prims: &[Prim],
    ids: &mut [u32],
    first: usize,
    depth: usize,
    pad: f64,
) {
    let idx = nodes.len();
    let mut b = EMPTY;
    let mut cb = EMPTY;
    for &i in ids.iter() {
        let p = &prims[i as usize];
        b = union(b, (p.min, p.max));
        cb = union(cb, (p.centroid, p.centroid));
    }
    nodes.push(Node {
        min: [b.0[0] - pad, b.0[1] - pad, b.0[2] - pad],
        max: [b.1[0] + pad, b.1[1] + pad, b.1[2] + pad],
        first: first as u32,
        count: ids.len() as u32,
    });
    let n = ids.len();
    if n <= 1 {
        return;
    }
    let mut mid = None;
    if depth < SAH_DEPTH {
        if let Some((cost, axis, bin)) = best_sah(prims, ids, b, cb) {
            // Leaf cost n, traversal cost 1 (relative units).
            if n <= MAX_LEAF && n as f64 <= cost {
                return;
            }
            mid = Some(partition(prims, ids, axis, bin, cb));
        }
    }
    let mid = match mid {
        Some(m) if m > 0 && m < n => m,
        _ => {
            if n <= MAX_LEAF && depth < SAH_DEPTH {
                return;
            }
            // Median split along the longest centroid axis.
            let ext = [cb.1[0] - cb.0[0], cb.1[1] - cb.0[1], cb.1[2] - cb.0[2]];
            let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                0
            } else if ext[1] >= ext[2] {
                1
            } else {
                2
            };
            ids.sort_by(|&a, &c| {
                prims[a as usize].centroid[axis]
                    .total_cmp(&prims[c as usize].centroid[axis])
                    .then(a.cmp(&c))
            });
            n / 2
        }
    };
    let (l, r) = ids.split_at_mut(mid);
    build_rec(nodes, prims, l, first, depth + 1, pad);
    let second = nodes.len();
    nodes[idx].first = second as u32;
    nodes[idx].count = 0;
    build_rec(nodes, prims, r, first + mid, depth + 1, pad);
}

fn bin_of(c: f64, lo: f64, ext: f64) -> usize {
    (((c - lo) / ext * BINS as f64) as usize).min(BINS - 1)
}

/// The cheapest binned split as `(cost, axis, last bin of the left side)`.
fn best_sah(prims: &[Prim], ids: &[u32], b: (V3, V3), cb: (V3, V3)) -> Option<(f64, usize, usize)> {
    let parent = area(b);
    let mut best: Option<(f64, usize, usize)> = None;
    for axis in 0..3 {
        let ext = cb.1[axis] - cb.0[axis];
        if ext.is_nan() || ext <= 0.0 {
            continue;
        }
        let mut count = [0usize; BINS];
        let mut bbox = [EMPTY; BINS];
        for &i in ids {
            let p = &prims[i as usize];
            let k = bin_of(p.centroid[axis], cb.0[axis], ext);
            count[k] += 1;
            bbox[k] = union(bbox[k], (p.min, p.max));
        }
        // Suffix sweep: right-hand side of each plane.
        let mut right = [(0usize, 0.0f64); BINS];
        let mut acc = EMPTY;
        let mut c = 0usize;
        for k in (1..BINS).rev() {
            acc = union(acc, bbox[k]);
            c += count[k];
            right[k - 1] = (c, area(acc));
        }
        let mut acc = EMPTY;
        let mut c = 0usize;
        for k in 0..BINS - 1 {
            acc = union(acc, bbox[k]);
            c += count[k];
            let (nr, ar) = right[k];
            if c == 0 || nr == 0 {
                continue;
            }
            let cost = 1.0
                + if parent > 0.0 {
                    (area(acc) * c as f64 + ar * nr as f64) / parent
                } else {
                    (c + nr) as f64
                };
            if best.is_none_or(|(bc, _, _)| cost < bc) {
                best = Some((cost, axis, k));
            }
        }
    }
    best
}

/// Stable partition of `ids` by bin; returns the size of the left part.
fn partition(prims: &[Prim], ids: &mut [u32], axis: usize, bin: usize, cb: (V3, V3)) -> usize {
    let ext = cb.1[axis] - cb.0[axis];
    let (mut left, mut right) = (Vec::new(), Vec::new());
    for &i in ids.iter() {
        if bin_of(prims[i as usize].centroid[axis], cb.0[axis], ext) <= bin {
            left.push(i);
        } else {
            right.push(i);
        }
    }
    let m = left.len();
    left.extend(right);
    ids.copy_from_slice(&left);
    m
}

/// A ray prepared for the watertight test.
pub(super) struct Ray {
    o: V3,
    inv: V3,
    kx: usize,
    ky: usize,
    kz: usize,
    sx: f64,
    sy: f64,
    sz: f64,
}

/// The top-left tie-break for an edge function that evaluates to exactly
/// zero: whether the (sheared, projected) edge with direction `(dx, dy)`
/// counts the ray as inside. It is `e(s) = e + dx sy - dy sx` for the ray
/// origin shifted by `s = (-eps, -eps^2)`, which for `e = 0` is positive
/// exactly when `dy > 0`, or `dy = 0` and `dx < 0`. Two triangles sharing an
/// edge traverse it in opposite directions, so exactly one of them owns it.
fn edge_owns(e: f64, dx: f64, dy: f64) -> bool {
    e > 0.0 || (e == 0.0 && (dy > 0.0 || (dy == 0.0 && dx < 0.0)))
}

impl Ray {
    /// `d` must be finite and non-zero (it need not be normalised; `t` is in
    /// units of `d`).
    pub(super) fn new(o: V3, d: V3) -> Self {
        let mut kz = 0;
        if d[1].abs() > d[kz].abs() {
            kz = 1;
        }
        if d[2].abs() > d[kz].abs() {
            kz = 2;
        }
        let mut kx = (kz + 1) % 3;
        let mut ky = (kx + 1) % 3;
        if d[kz] < 0.0 {
            std::mem::swap(&mut kx, &mut ky);
        }
        Self {
            o,
            inv: [1.0 / d[0], 1.0 / d[1], 1.0 / d[2]],
            kx,
            ky,
            kz,
            sx: d[kx] / d[kz],
            sy: d[ky] / d[kz],
            sz: 1.0 / d[kz],
        }
    }

    /// Slab test against a padded box; a NaN slab (origin on its plane with a
    /// zero direction component) is ignored, which only errs on the side of
    /// visiting the node.
    fn hits_box(&self, n: &Node, tmin: f64, tmax: f64) -> bool {
        let mut lo = tmin;
        let mut hi = tmax;
        for a in 0..3 {
            let t1 = (n.min[a] - self.o[a]) * self.inv[a];
            let t2 = (n.max[a] - self.o[a]) * self.inv[a];
            lo = lo.max(t1.min(t2));
            hi = hi.min(t1.max(t2));
        }
        lo <= hi
    }

    /// The ray parameter at which the ray crosses the triangle, for either
    /// facing and any sign of `t`, or `None` if it does not (see the module
    /// docs for the tie-break on edges and vertices).
    pub(super) fn intersect(&self, tri: &[V3; 3]) -> Option<f64> {
        let (kx, ky, kz) = (self.kx, self.ky, self.kz);
        let mut x = [0.0; 3];
        let mut y = [0.0; 3];
        let mut z = [0.0; 3];
        for i in 0..3 {
            let p = [
                tri[i][0] - self.o[0],
                tri[i][1] - self.o[1],
                tri[i][2] - self.o[2],
            ];
            x[i] = p[kx] - self.sx * p[kz];
            y[i] = p[ky] - self.sy * p[kz];
            z[i] = self.sz * p[kz];
        }
        // Edge functions: u opposite A (edge BC), v opposite B, w opposite C.
        let u = x[1] * y[2] - y[1] * x[2];
        let v = x[2] * y[0] - y[2] * x[0];
        let w = x[0] * y[1] - y[0] * x[1];
        let det = u + v + w;
        if det == 0.0 || det.is_nan() {
            return None;
        }
        // Edge directions in the projected plane.
        let d = [
            (x[2] - x[1], y[2] - y[1]),
            (x[0] - x[2], y[0] - y[2]),
            (x[1] - x[0], y[1] - y[0]),
        ];
        // A back-facing (clockwise) projection is the mirror image: negate
        // the edge functions and reverse the edges.
        let s = if det < 0.0 { -1.0 } else { 1.0 };
        for (e, (dx, dy)) in [u, v, w].into_iter().zip(d) {
            if !edge_owns(s * e, s * dx, s * dy) {
                return None;
            }
        }
        Some((u * z[0] + v * z[1] + w * z[2]) / det)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad_pair() -> [[V3; 3]; 2] {
        // A unit square in z = 0 split along its diagonal.
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [1.0, 1.0, 0.0];
        let d = [0.0, 1.0, 0.0];
        [[a, b, c], [a, c, d]]
    }

    #[test]
    fn diagonal_hit_once_either_facing() {
        let tris = quad_pair();
        for dz in [1.0, -1.0] {
            // On the diagonal x = y, exactly.
            let ray = Ray::new([0.25, 0.25, -dz], [0.0, 0.0, dz]);
            let hits: Vec<_> = tris.iter().filter_map(|t| ray.intersect(t)).collect();
            assert_eq!(hits.len(), 1, "{hits:?}");
            assert!((hits[0] - 1.0).abs() < 1e-15);
        }
    }

    #[test]
    fn vertex_hit_once_for_a_fan() {
        // Four triangles around the origin, ray through the common vertex.
        let o = [0.0; 3];
        let p = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let tris: Vec<[V3; 3]> = (0..4).map(|i| [o, p[i], p[(i + 1) % 4]]).collect();
        let ray = Ray::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        let n = tris.iter().filter(|t| ray.intersect(t).is_some()).count();
        assert_eq!(n, 1);
    }

    #[test]
    fn bvh_matches_brute_force() {
        // A strip of triangles; every ray via the BVH equals the linear scan.
        let mut tris = Vec::new();
        for i in 0..50 {
            let x = i as f64;
            tris.push([[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.5 + 0.01 * x]]);
        }
        let bvh = Bvh::build(&tris, 1e-12);
        for i in 0..50 {
            let ray = Ray::new([i as f64 + 0.2, 0.2, -1.0], [0.0, 0.1, 1.0]);
            let mut a = Vec::new();
            bvh.for_each_hit(&ray, f64::NEG_INFINITY, f64::INFINITY, |t, _| a.push(t));
            let mut b: Vec<f64> = tris.iter().filter_map(|t| ray.intersect(t)).collect();
            a.sort_by(f64::total_cmp);
            b.sort_by(f64::total_cmp);
            assert_eq!(a, b);
        }
        assert_eq!(bvh, Bvh::build(&tris, 1e-12));
    }
}
