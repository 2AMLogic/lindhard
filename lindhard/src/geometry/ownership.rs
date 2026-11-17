//! The ownership walk shared by the triangle-mesh and the CSG target.
//!
//! [`MeshGeometry`](super::MeshGeometry) and [`CsgGeometry`](super::CsgGeometry)
//! implement one set of overlap, boundary and tolerance rules -- stated in the
//! [`mesh`](super::mesh) module docs and adopted by the [`csg`](super::csg)
//! one -- so that a solid as a mesh, the same solid as a CSG tree and the same
//! solid as a [`VoxelGrid`](super::VoxelGrid) agree on every shared face. The
//! parts of that implementation which do not depend on how a solid is
//! represented live here, so that a change to a rule is made once:
//!
//! * [`TOL_REL`], the relative size of the geometric tolerance, and
//!   [`face_of`], the face class of an outward normal.
//! * [`Segment`]: the argument guard of
//!   [`Geometry::flight`](super::Geometry::flight), and the point at a path
//!   length along the flight.
//! * [`shifted`], the fixed offset of rule 2, and [`in_vacuum`], the eight
//!   mirrored offsets built on it.
//! * [`entry_point`], the centre of the front face of the scene box.
//! * [`walk`], the cluster loop that turns a sorted list of crossings into the
//!   first event of a flight, and [`exit_of`], which narrows its answer to
//!   [`Geometry::exit`](super::Geometry::exit).
//!
//! What each target still owns is what depends on its representation: how a
//! flight becomes a sorted list of crossings and a start state (BVH triangle
//! hits with a parity count for meshes, Roth span ends for CSG), and which
//! crossing of a cluster supplies the normal of an escape.
//!
//! Nothing here dispatches dynamically: [`walk`] is monomorphised per crossing
//! type and per closure, so a flight runs the same code as the per-target
//! copies this module replaces.

use super::{Exit, ExitOutcome, Face, Flight};
use crate::vec3::V3;

/// Relative size of the geometric tolerance: `tol = TOL_REL x` the scene
/// diagonal.
pub(super) const TOL_REL: f64 = 1e-10;

/// The face class of an outward normal: [`Face::Front`] if `n_x` is negative
/// and the largest component in magnitude, [`Face::Back`] if positive and
/// largest, [`Face::Side`] otherwise (see the [`mesh`](super::mesh) module
/// docs).
pub(super) fn face_of(n: V3) -> Face {
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

/// Where a beam enters by default: the centre of the front face of the scene
/// bounding box `(lo, hi)`.
pub(super) fn entry_point(lo: V3, hi: V3) -> V3 {
    [lo[0], 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])]
}

/// `pos` moved by the fixed offset of rule 2, `tol * (1, 1, 1) / sqrt(3)` with
/// the given sign on each axis (all `+1` for
/// [`Geometry::locate`](super::Geometry::locate)), or `None` for a non-finite
/// position.
pub(super) fn shifted(pos: V3, sign: V3, tol: f64) -> Option<V3> {
    if !pos.iter().all(|c| c.is_finite()) {
        return None;
    }
    let s = tol / 3f64.sqrt();
    Some([
        pos[0] + sign[0] * s,
        pos[1] + sign[1] * s,
        pos[2] + sign[2] * s,
    ])
}

/// [`Geometry::in_vacuum`](super::Geometry::in_vacuum) for a set of solids: a
/// point on the surface of one is not in vacuum, whichever of the eight
/// diagonal offsets (the one of [`Geometry::locate`](super::Geometry::locate)
/// and its mirror images) would put it outside. `locate_shifted` is the
/// target's own `locate` under the sign pattern it is given.
pub(super) fn in_vacuum(mut locate_shifted: impl FnMut(V3) -> Option<usize>) -> bool {
    (0..8).all(|m| {
        let sign = |k: usize| if m >> k & 1 == 0 { 1.0 } else { -1.0 };
        locate_shifted([sign(0), sign(1), sign(2)]).is_none()
    })
}

/// [`Geometry::exit`](super::Geometry::exit) from the answer of
/// [`Geometry::flight`](super::Geometry::flight): both targets report the
/// region of a clear flight, and `exit` drops it.
pub(super) fn exit_of(flight: Flight) -> Option<Exit> {
    match flight {
        Flight::Event(e) => Some(e),
        Flight::Clear { .. } => None,
    }
}

/// The straight flight a [`walk`] runs along: the arguments of
/// [`Geometry::flight`](super::Geometry::flight) that describe the ray, plus
/// the geometry's tolerance.
pub(super) struct Segment {
    /// Start of the flight, m.
    pub(super) pos: V3,
    /// Unit direction.
    pub(super) dir: V3,
    /// Largest path length the flight may reach, m.
    pub(super) limit: f64,
    /// The geometric tolerance of the ownership rules, m.
    pub(super) tol: f64,
}

impl Segment {
    /// The argument guard of [`Geometry::flight`](super::Geometry::flight):
    /// `None` for a zero or non-finite direction, a non-finite start, or a
    /// NaN limit, which every target answers with
    /// [`Flight::Clear`](super::Flight::Clear) and no event.
    pub(super) fn new(pos: V3, dir: V3, limit: f64, tol: f64) -> Option<Self> {
        if !dir.iter().all(|d| d.is_finite())
            || dir == [0.0; 3]
            || !pos.iter().all(|c| c.is_finite())
            || limit.is_nan()
        {
            return None;
        }
        Some(Self {
            pos,
            dir,
            limit,
            tol,
        })
    }

    /// The point at path length `t` along the flight, m.
    pub(super) fn at(&self, t: f64) -> V3 {
        [
            self.pos[0] + t * self.dir[0],
            self.pos[1] + t * self.dir[1],
            self.pos[2] + t * self.dir[2],
        ]
    }
}

/// One crossing of the surface of one solid along a flight, as [`walk`] reads
/// it. Implemented for each target's own crossing record, so that neither the
/// contents of that record nor the order the target sorts them in is fixed
/// here.
pub(super) trait Crossing {
    /// Path length from the start of the flight to the crossing, m.
    fn t(&self) -> f64;

    /// The solid whose surface is crossed.
    fn solid(&self) -> usize;

    /// The state of that solid just after the crossing: `Some(b)` to set it
    /// (a CSG span end knows which side it is), `None` to toggle it (a mesh
    /// crossing only knows that the surface was crossed, so the state is the
    /// running parity).
    fn state_after(&self) -> Option<bool>;
}

/// A mesh crossing: `(t, solid, triangle)`. The surface was crossed, so the
/// solid's state toggles.
impl Crossing for (f64, u32, u32) {
    fn t(&self) -> f64 {
        self.0
    }

    fn solid(&self) -> usize {
        self.1 as usize
    }

    fn state_after(&self) -> Option<bool> {
        None
    }
}

/// A CSG span end: `(t, solid, entering, normal)`. The end names the side, so
/// the solid's state is assigned.
impl Crossing for (f64, u32, bool, V3) {
    fn t(&self) -> f64 {
        self.0
    }

    fn solid(&self) -> usize {
        self.1 as usize
    }

    fn state_after(&self) -> Option<bool> {
        Some(self.2)
    }
}

/// The first event along `seg`, given every crossing of every solid within
/// `(tol, limit + tol]` sorted by path length and the state `inside` of every
/// solid at `t = tol` (the flight's own `region` counted as inside).
///
/// Crossings within `tol` of each other are one cluster and so one event
/// (rule 4). After each cluster the owner is the first solid still inside
/// (rule 1): unchanged, a same-material solid (not an event, the walk carries
/// on in the new region), a different material
/// ([`ExitOutcome::Enter`]), or none ([`ExitOutcome::Escape`], with the
/// normal `escape_normal` picks from the cluster). A mismatch at the start of
/// the flight -- a surface start heading into a solid that has priority or is
/// not `region` -- is an event at distance 0 (rule 2). With no event before
/// `limit` the flight is [`Flight::Clear`] in the region it ends in.
///
/// `inside` is left in an unspecified state.
pub(super) fn walk<C: Crossing>(
    seg: &Segment,
    region: usize,
    crossings: &[C],
    inside: &mut [bool],
    material_of: impl Fn(usize) -> u32,
    escape_normal: impl Fn(&[C], usize) -> V3,
) -> Flight {
    let owner_of = |inside: &[bool]| inside.iter().position(|&b| b);
    let mut current = region;
    let mut pending_start = true;
    let mut i = 0usize;
    loop {
        let (t0, cluster_end) = if pending_start {
            (0.0, i)
        } else {
            if i >= crossings.len() || crossings[i].t() > seg.limit {
                return Flight::Clear { region: current };
            }
            let t0 = crossings[i].t();
            let mut e = i;
            while e < crossings.len() && crossings[e].t() - t0 <= seg.tol {
                let j = crossings[e].solid();
                inside[j] = crossings[e].state_after().unwrap_or(!inside[j]);
                e += 1;
            }
            (t0, e)
        };
        let cluster = &crossings[i..cluster_end];
        pending_start = false;
        i = cluster_end;
        let at = seg.at(t0);
        match owner_of(inside) {
            Some(k) if k == current => {}
            Some(k) => {
                if material_of(k) == material_of(current) {
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
                let normal = escape_normal(cluster, current);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn face_of_takes_the_dominant_axis() {
        assert_eq!(face_of([-1.0, 0.0, 0.0]), Face::Front);
        assert_eq!(face_of([1.0, 0.0, 0.0]), Face::Back);
        assert_eq!(face_of([0.0, -1.0, 0.0]), Face::Side);
        assert_eq!(face_of([0.0, 0.0, 1.0]), Face::Side);
        // Ties are a side: no axis is strictly largest.
        let s = 1.0 / 2f64.sqrt();
        assert_eq!(face_of([-s, s, 0.0]), Face::Side);
        assert_eq!(face_of([-0.8, 0.6, 0.0]), Face::Front);
    }

    #[test]
    fn entry_point_is_the_centre_of_the_front_face() {
        assert_eq!(
            entry_point([0.0, -1.0, 2.0], [3.0, 1.0, 4.0]),
            [0.0, 0.0, 3.0]
        );
    }

    #[test]
    fn shifted_moves_along_the_diagonal_and_rejects_non_finite() {
        let s = 1.0 / 3f64.sqrt();
        assert_eq!(shifted([0.0; 3], [1.0; 3], 1.0), Some([s, s, s]));
        assert_eq!(shifted([0.0; 3], [1.0, -1.0, 1.0], 1.0), Some([s, -s, s]));
        assert_eq!(shifted([f64::NAN, 0.0, 0.0], [1.0; 3], 1.0), None);
        assert_eq!(shifted([f64::INFINITY, 0.0, 0.0], [1.0; 3], 1.0), None);
    }

    #[test]
    fn in_vacuum_tries_all_eight_offsets() {
        let mut seen: Vec<V3> = Vec::new();
        assert!(in_vacuum(|sign| {
            seen.push(sign);
            None
        }));
        assert_eq!(seen.len(), 8);
        for s in &seen {
            assert!(s.iter().all(|c| c.abs() == 1.0));
        }
        seen.sort_by(|a, b| a.partial_cmp(b).unwrap());
        seen.dedup();
        assert_eq!(seen.len(), 8);
        // One offset finding a solid is enough to say the point is not in
        // vacuum.
        assert!(!in_vacuum(|sign| if sign == [1.0, 1.0, -1.0] {
            Some(0)
        } else {
            None
        }));
    }

    #[test]
    fn segment_guards_its_arguments() {
        assert!(Segment::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 1e-9).is_some());
        assert!(Segment::new([0.0; 3], [1.0, 0.0, 0.0], f64::INFINITY, 1e-9).is_some());
        assert!(Segment::new([0.0; 3], [0.0; 3], 1.0, 1e-9).is_none());
        assert!(Segment::new([0.0; 3], [f64::NAN, 0.0, 0.0], 1.0, 1e-9).is_none());
        assert!(Segment::new([f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 1e-9).is_none());
        assert!(Segment::new([0.0; 3], [1.0, 0.0, 0.0], f64::NAN, 1e-9).is_none());
    }

    #[test]
    fn segment_at_walks_the_ray() {
        let seg = Segment::new([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], 10.0, 1e-9).unwrap();
        assert_eq!(seg.at(0.0), [1.0, 2.0, 3.0]);
        assert_eq!(seg.at(2.5), [1.0, 4.5, 3.0]);
    }

    /// `walk` over mesh-shaped crossings (toggling), one material per solid.
    fn mesh_walk(
        region: usize,
        inside: &mut [bool],
        hits: &[(f64, u32, u32)],
        limit: f64,
    ) -> Flight {
        let seg = Segment::new([0.0; 3], [1.0, 0.0, 0.0], limit, 0.1).unwrap();
        walk(
            &seg,
            region,
            hits,
            inside,
            |i| i as u32,
            |_, _| [-1.0, 0.0, 0.0],
        )
    }

    #[test]
    fn walk_reports_entering_and_escaping() {
        // One solid, the flight starting inside it and leaving at 2.
        let mut inside = [true];
        let f = mesh_walk(0, &mut inside, &[(2.0, 0, 7)], 10.0);
        match f {
            Flight::Event(e) => {
                assert_eq!(e.distance, 2.0);
                assert_eq!(e.at, [2.0, 0.0, 0.0]);
                assert!(matches!(e.outcome, ExitOutcome::Escape { .. }));
            }
            Flight::Clear { .. } => panic!("expected an escape, got {f:?}"),
        }
        // The flight is in solid 1 and enters solid 0 at 3; solid 0 has
        // priority (rule 1) and a different material, so that is an event.
        let mut inside = [false, true];
        let f = mesh_walk(1, &mut inside, &[(3.0, 0, 0)], 10.0);
        match f {
            Flight::Event(e) => {
                assert_eq!(e.distance, 3.0);
                assert_eq!(
                    e.outcome,
                    ExitOutcome::Enter {
                        region: 0,
                        pos: [3.0, 0.0, 0.0],
                    }
                );
            }
            Flight::Clear { .. } => panic!("expected an entry, got {f:?}"),
        }
        // The same crossing of a solid listed *after* the current one does
        // not change the owner: rule 1 keeps the first listed solid.
        let mut inside = [true, false];
        let f = mesh_walk(0, &mut inside, &[(3.0, 1, 0)], 10.0);
        assert_eq!(f, Flight::Clear { region: 0 });
    }

    #[test]
    fn walk_clusters_crossings_within_tol() {
        // Solid 0 is left and solid 1 entered 0.05 later, inside tol = 0.1:
        // one event, and no escape through the vacuum between them.
        let mut inside = [true, false];
        let f = mesh_walk(0, &mut inside, &[(1.0, 0, 0), (1.05, 1, 0)], 10.0);
        match f {
            Flight::Event(e) => {
                assert_eq!(e.distance, 1.0);
                assert!(matches!(e.outcome, ExitOutcome::Enter { region: 1, .. }));
            }
            Flight::Clear { .. } => panic!("expected one event, got {f:?}"),
        }
    }

    #[test]
    fn walk_is_clear_past_the_limit_and_carries_same_material_regions() {
        // The crossing at 5 is beyond the limit.
        let mut inside = [true];
        assert_eq!(
            mesh_walk(0, &mut inside, &[(5.0, 0, 0)], 2.0),
            Flight::Clear { region: 0 }
        );
        // Two solids of one material: entering the second is not an event,
        // and the clear flight ends in it.
        let seg = Segment::new([0.0; 3], [1.0, 0.0, 0.0], 10.0, 0.1).unwrap();
        let mut inside = [false, true];
        let f = walk(
            &seg,
            1,
            &[(1.0, 0, 0)],
            &mut inside,
            |_| 0,
            |_, _| [-1.0, 0.0, 0.0],
        );
        assert_eq!(f, Flight::Clear { region: 0 });
    }

    #[test]
    fn walk_reports_a_start_mismatch_at_distance_zero() {
        // The flight starts in region 1 but solid 0 has priority and is also
        // inside: an event at 0.
        let mut inside = [true, true];
        let f = mesh_walk(1, &mut inside, &[], 10.0);
        match f {
            Flight::Event(e) => {
                assert_eq!(e.distance, 0.0);
                assert!(matches!(e.outcome, ExitOutcome::Enter { region: 0, .. }));
            }
            Flight::Clear { .. } => panic!("expected a start event, got {f:?}"),
        }
    }

    #[test]
    fn walk_assigns_state_for_csg_span_ends() {
        // A CSG span end sets the state rather than toggling it, so two
        // "leaving" ends in a row leave the solid outside either way.
        let seg = Segment::new([0.0; 3], [1.0, 0.0, 0.0], 10.0, 1e-3).unwrap();
        let n = [0.0, 1.0, 0.0];
        let mut inside = [true];
        let f = walk(
            &seg,
            0,
            &[(1.0, 0, false, n), (2.0, 0, false, n)],
            &mut inside,
            |i| i as u32,
            |cluster: &[(f64, u32, bool, V3)], _| cluster[0].3,
        );
        match f {
            Flight::Event(e) => {
                assert_eq!(e.distance, 1.0);
                assert_eq!(
                    e.outcome,
                    ExitOutcome::Escape {
                        face: Face::Side,
                        normal: n,
                    }
                );
            }
            Flight::Clear { .. } => panic!("expected an escape, got {f:?}"),
        }
    }
}
