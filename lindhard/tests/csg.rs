//! Constructive-solid-geometry targets (`lindhard::geometry::CsgGeometry`):
//! analytic chords through a box, a cylinder and a box with a cylindrical
//! hole; agreement with the voxel grid and the triangle mesh on the same
//! boxes; one test per ownership rule; rejection of unbounded solids; and
//! determinism across threads.

use lindhard::geometry::{
    Boundary, Csg, CsgGeometry, Exit, ExitOutcome, Face, Flight, Geometry, GeometryError,
    MeshGeometry, TriMesh, VoxelGrid,
};
use lindhard::material::Material;

type V3 = [f64; 3];

fn si() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn carbon() -> Material {
    Material::from_atom_fractions(&[(6, 1.0)], None).unwrap()
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn dir(&mut self) -> V3 {
        loop {
            let d = [
                2.0 * self.next() - 1.0,
                2.0 * self.next() - 1.0,
                2.0 * self.next() - 1.0,
            ];
            let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if l > 0.1 && l < 1.0 {
                return [d[0] / l, d[1] / l, d[2] / l];
            }
        }
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * (1.0 + a.abs().max(b.abs()))
}

fn escape(e: &Exit) -> (Face, V3) {
    match e.outcome {
        ExitOutcome::Escape { face, normal } => (face, normal),
        ExitOutcome::Enter { .. } => panic!("expected Escape, got {e:?}"),
    }
}

fn enter(e: &Exit) -> (usize, V3) {
    match e.outcome {
        ExitOutcome::Enter { region, pos } => (region, pos),
        ExitOutcome::Escape { .. } => panic!("expected Enter, got {e:?}"),
    }
}

fn one(c: Csg) -> CsgGeometry {
    CsgGeometry::new(vec![si()], vec![(c, 0)]).unwrap()
}

/// Distance from `o` along `d` to the far side of `[lo, hi]` (slab method).
fn box_chord(o: V3, d: V3, lo: V3, hi: V3) -> f64 {
    let mut t = f64::INFINITY;
    for k in 0..3 {
        if d[k] > 0.0 {
            t = t.min((hi[k] - o[k]) / d[k]);
        } else if d[k] < 0.0 {
            t = t.min((lo[k] - o[k]) / d[k]);
        }
    }
    t
}

/// 12 counter-clockwise (outward) triangles of the box `[lo, hi]`.
fn box_soup(lo: V3, hi: V3) -> Vec<[V3; 3]> {
    let p = |x: usize, y: usize, z: usize| {
        [
            if x == 0 { lo[0] } else { hi[0] },
            if y == 0 { lo[1] } else { hi[1] },
            if z == 0 { lo[2] } else { hi[2] },
        ]
    };
    let quads = [
        [p(0, 0, 0), p(0, 0, 1), p(0, 1, 1), p(0, 1, 0)],
        [p(1, 0, 0), p(1, 1, 0), p(1, 1, 1), p(1, 0, 1)],
        [p(0, 0, 0), p(1, 0, 0), p(1, 0, 1), p(0, 0, 1)],
        [p(0, 1, 0), p(0, 1, 1), p(1, 1, 1), p(1, 1, 0)],
        [p(0, 0, 0), p(0, 1, 0), p(1, 1, 0), p(1, 0, 0)],
        [p(0, 0, 1), p(1, 0, 1), p(1, 1, 1), p(0, 1, 1)],
    ];
    quads
        .iter()
        .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
        .collect()
}

fn box_mesh(lo: V3, hi: V3) -> TriMesh {
    TriMesh::from_triangles(&box_soup(lo, hi)).unwrap()
}

/// The box `[lo, hi]` as the intersection of six half-spaces.
fn six_planes(lo: V3, hi: V3) -> Csg {
    let mut v = Vec::new();
    for k in 0..3 {
        let mut n = [0.0; 3];
        n[k] = -1.0;
        v.push(Csg::half_space(n, -lo[k]));
        n[k] = 1.0;
        v.push(Csg::half_space(n, hi[k]));
    }
    Csg::Intersection(v)
}

/// Follow a particle to the surface, collecting `(kind, cumulative distance,
/// position, material)` of every event.
fn trace(g: &impl Geometry, mut pos: V3, d: V3) -> Vec<(char, f64, V3, usize)> {
    let mut region = g.locate(pos).unwrap();
    let mut total = 0.0;
    let mut out = Vec::new();
    for _ in 0..16 {
        let e = g.exit(region, pos, d, 1e9).expect("a closed target");
        total += e.distance;
        match e.outcome {
            ExitOutcome::Enter { region: r, pos: p } => {
                out.push(('E', total, p, g.material_index(r)));
                region = r;
                pos = p;
            }
            ExitOutcome::Escape { .. } => {
                out.push(('X', total, e.at, g.material_index(region)));
                return out;
            }
        }
    }
    panic!("too many events");
}

fn same_trace(a: &[(char, f64, V3, usize)], b: &[(char, f64, V3, usize)]) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(x, y)| {
            (x.0, x.3) == (y.0, y.3) && close(x.1, y.1) && (0..3).all(|k| close(x.2[k], y.2[k]))
        })
}

// ---------------------------------------------------------------- chords ----

#[test]
fn box_gives_exact_chords() {
    let (lo, hi) = ([0.0, 0.0, 0.0], [3.0, 2.0, 1.0]);
    let g = one(Csg::cuboid(lo, hi));
    assert_eq!(g.n_regions(), 1);
    assert_eq!(g.n_materials(), 1);
    assert_eq!(g.material_index(0), 0);
    assert_eq!(g.entry_point(), [0.0, 1.0, 0.5]);
    let e = g.exit(0, [0.0, 1.0, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, 3.0);
    assert_eq!(escape(&e), (Face::Back, [1.0, 0.0, 0.0]));
    let e = g.exit(0, [1.0, 1.0, 0.5], [-1.0, 0.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, 1.0);
    assert_eq!(escape(&e), (Face::Front, [-1.0, 0.0, 0.0]));
    let e = g.exit(0, [1.0, 1.0, 0.25], [0.0, 1.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, 1.0);
    assert_eq!(escape(&e), (Face::Side, [0.0, 1.0, 0.0]));
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for _ in 0..2000 {
        let o = [3.0 * rng.next(), 2.0 * rng.next(), rng.next()];
        let d = rng.dir();
        let e = g.exit(0, o, d, 1e9).unwrap();
        assert!(close(e.distance, box_chord(o, d, lo, hi)), "{o:?} {d:?}");
    }
    // A limit short of the face: no event.
    assert!(g.exit(0, [1.0, 1.0, 0.5], [1.0, 0.0, 0.0], 1.9).is_none());
    assert!(g.exit(0, [1.0, 1.0, 0.5], [1.0, 0.0, 0.0], 2.0).is_some());
    assert_eq!(
        g.flight(0, [1.0, 1.0, 0.5], [1.0, 0.0, 0.0], 1.0),
        Flight::Clear { region: 0 }
    );
    // Degenerate directions are ignored.
    assert!(g.exit(0, [1.0, 1.0, 0.5], [0.0; 3], 1e9).is_none());
    assert!(g
        .exit(0, [1.0, 1.0, 0.5], [f64::NAN, 0.0, 0.0], 1e9)
        .is_none());
    // Locate and vacuum, voxel convention.
    assert_eq!(g.locate([0.0, 0.0, 0.0]), Some(0));
    assert_eq!(g.locate([3.0, 1.0, 0.5]), None);
    assert_eq!(g.locate([1.0, 1.0, f64::NAN]), None);
    assert!(g.in_vacuum([5.0, 0.0, 0.0]));
    assert!(!g.in_vacuum([3.0, 1.0, 0.5]));
}

/// Analytic exit distance from inside the cylinder of radius `r` about the
/// `x` axis, `0 <= x <= len`, and whether it is through a cap.
fn cylinder_exit(o: V3, d: V3, r: f64, len: f64) -> (f64, bool) {
    let a = d[1] * d[1] + d[2] * d[2];
    let b = o[1] * d[1] + o[2] * d[2];
    let c = o[1] * o[1] + o[2] * o[2] - r * r;
    let side = if a > 0.0 {
        (-b + (b * b - a * c).sqrt()) / a
    } else {
        f64::INFINITY
    };
    let cap = if d[0] > 0.0 {
        (len - o[0]) / d[0]
    } else if d[0] < 0.0 {
        -o[0] / d[0]
    } else {
        f64::INFINITY
    };
    if cap < side {
        (cap, true)
    } else {
        (side, false)
    }
}

#[test]
fn cylinder_gives_analytic_chords() {
    let (r, len) = (0.5, 3.0);
    let g = one(Csg::cylinder([0.0; 3], [1.0, 0.0, 0.0], r, len));
    assert_eq!(g.bounds(), ([0.0, -r, -r], [len, r, r]));
    // Along the axis: the height.
    let e = g.exit(0, [0.0, 0.1, 0.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, len);
    assert_eq!(escape(&e), (Face::Back, [1.0, 0.0, 0.0]));
    // Across the axis: the diameter, then an offset chord 2 sqrt(r^2 - b^2).
    let e = g.exit(0, [1.0, -r, 0.0], [0.0, 1.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0 * r));
    assert_eq!(escape(&e), (Face::Side, [0.0, 1.0, 0.0]));
    let b = 0.3;
    let h = (r * r - b * b).sqrt();
    let e = g.exit(0, [1.0, -h, b], [0.0, 1.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0 * h), "{e:?}");
    let (_, n) = escape(&e);
    assert!(close(n[1], h / r) && close(n[2], b / r));
    // Oblique, from inside.
    let mut rng = Rng(31337);
    for _ in 0..2000 {
        let (rho, phi) = (r * rng.next().sqrt(), std::f64::consts::TAU * rng.next());
        let o = [len * rng.next(), rho * phi.cos(), rho * phi.sin()];
        let d = rng.dir();
        let e = g.exit(0, o, d, 1e9).unwrap();
        let (t, cap) = cylinder_exit(o, d, r, len);
        assert!(close(e.distance, t), "{o:?} {d:?} {e:?} {t}");
        let (face, _) = escape(&e);
        if cap {
            assert_ne!(face, Face::Side);
        }
    }
    // A tangent line does not enter.
    let tangent = one(Csg::cuboid([-1.0; 3], [4.0, 1.0, 1.0]).minus(Csg::cylinder(
        [0.0; 3],
        [1.0, 0.0, 0.0],
        r,
        len,
    )));
    let e = tangent
        .exit(0, [1.0, -1.0, r], [0.0, 1.0, 0.0], 1e9)
        .unwrap();
    assert!(close(e.distance, 2.0), "{e:?}");
}

#[test]
fn box_minus_cylinder_hole() {
    // A 2 x 2 x 2 box with a hole of radius 0.5 drilled along z through
    // its centre.
    let r = 0.5;
    let drill = || Csg::cylinder([1.0, 1.0, -1.0], [0.0, 0.0, 1.0], r, 4.0);
    let block = || Csg::cuboid([0.0; 3], [2.0; 3]);
    let g = one(block().minus(drill()));
    // Bounded by the box it is cut from.
    assert_eq!(g.bounds(), ([0.0; 3], [2.0; 3]));
    for b in [0.0, 0.2, 0.45] {
        let h = (r * r - b * b).sqrt();
        let o = [0.0, 1.0 + b, 1.0];
        assert_eq!(g.locate(o), Some(0));
        // Into the hole: an escape at the hole wall, whose outward normal (of
        // the drilled solid) points into the hole.
        let e = g.exit(0, o, [1.0, 0.0, 0.0], 1e9).unwrap();
        assert!(close(e.distance, 1.0 - h), "{b} {e:?}");
        let (face, n) = escape(&e);
        assert!(close(n[0], h / r) && close(n[1], -b / r), "{n:?}");
        assert_eq!(face, if h > b { Face::Back } else { Face::Side });
        // The hole is vacuum.
        assert!(g.in_vacuum([1.0, 1.0 + b, 1.0]));
        assert_eq!(g.locate([1.0, 1.0 + b, 1.0]), None);
        // From the far wall to the back face.
        let far = [1.0 + h, 1.0 + b, 1.0];
        assert_eq!(g.locate(far), Some(0));
        let e = g.exit(0, far, [1.0, 0.0, 0.0], 1e9).unwrap();
        assert!(close(e.distance, 1.0 - h));
    }
    // Missing the hole: the full chord.
    let e = g.exit(0, [0.0, 1.6, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0));
    // Parallel to the hole, beside it.
    let e = g.exit(0, [0.2, 0.2, 0.0], [0.0, 0.0, 1.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0));
    // The same hole filled with a plug of another material listed first:
    // the plug owns the overlap, so the events are the hole walls.
    let plugged = CsgGeometry::new(vec![si(), carbon()], vec![(drill(), 1), (block(), 0)]).unwrap();
    let b = 0.2;
    let h = (r * r - b * b).sqrt();
    let t = trace(&plugged, [0.0, 1.0 + b, 1.0], [1.0, 0.0, 0.0]);
    let want = [('E', 1.0 - h, 1), ('E', 1.0 + h, 0), ('X', 2.0, 0)];
    assert_eq!(t.len(), want.len(), "{t:?}");
    for (x, w) in t.iter().zip(&want) {
        assert_eq!((x.0, x.3), (w.0, w.2));
        assert!(close(x.1, w.1), "{t:?}");
    }
}

// ----------------------------------------------------------- equivalence ----

#[test]
fn csg_box_matches_voxel_grid_and_mesh() {
    let h = 0.5;
    let extent = [4.0 * h, 2.0 * h, 2.0 * h];
    let voxel = VoxelGrid::new(
        vec![si()],
        [4, 2, 2],
        [h; 3],
        vec![0; 16],
        [Boundary::Vacuum; 3],
    )
    .unwrap();
    let mesh = MeshGeometry::new(vec![si()], vec![(box_mesh([0.0; 3], extent), 0)]).unwrap();
    let csg = one(Csg::cuboid([0.0; 3], extent));
    let planes = one(six_planes([0.0; 3], extent));
    assert_eq!(csg.tolerance_m(), mesh.tolerance_m());
    assert_eq!(planes.bounds(), csg.bounds());
    for g in [&csg as &dyn Geometry, &planes] {
        assert_eq!(g.entry_point(), voxel.entry_point());
        assert_eq!(g.entry_point(), mesh.entry_point());
    }
    let mut rng = Rng(12345);
    for _ in 0..3000 {
        let o = [
            extent[0] * rng.next(),
            extent[1] * rng.next(),
            extent[2] * rng.next(),
        ];
        let d = rng.dir();
        let ev = voxel.exit(voxel.locate(o).unwrap(), o, d, 1e9).unwrap();
        let em = mesh.exit(mesh.locate(o).unwrap(), o, d, 1e9).unwrap();
        for g in [&csg as &dyn Geometry, &planes] {
            let ec = g.exit(g.locate(o).unwrap(), o, d, 1e9).unwrap();
            for other in [&ev, &em] {
                assert!(close(ec.distance, other.distance), "{ec:?} {other:?}");
                for k in 0..3 {
                    assert!(close(ec.at[k], other.at[k]));
                }
                assert_eq!(escape(&ec), escape(other));
            }
        }
    }
    // Points on and beside faces agree in `locate` and `in_vacuum`.
    for x in [0.0, 0.5, 1.0, 2.0, 2.5, -0.5] {
        for y in [0.0, 0.5, 1.0, 1.5] {
            for z in [0.0, 0.25, 1.0] {
                let p = [x, y, z];
                for g in [&csg as &dyn Geometry, &planes] {
                    assert_eq!(g.locate(p).is_some(), voxel.locate(p).is_some(), "{p:?}");
                    assert_eq!(g.locate(p).is_some(), mesh.locate(p).is_some(), "{p:?}");
                    assert_eq!(g.in_vacuum(p), voxel.in_vacuum(p), "{p:?}");
                    assert_eq!(g.in_vacuum(p), mesh.in_vacuum(p), "{p:?}");
                }
            }
        }
    }
}

#[test]
fn two_material_csg_matches_voxel_grid_and_mesh() {
    let h = 0.5;
    let cells: Vec<u32> = (0..16).map(|c| u32::from(c % 4 >= 2)).collect();
    let voxel = VoxelGrid::new(
        vec![si(), carbon()],
        [4, 2, 2],
        [h; 3],
        cells,
        [Boundary::Vacuum; 3],
    )
    .unwrap();
    let (a, b) = (([0.0; 3], [1.0; 3]), ([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]));
    let mesh = MeshGeometry::new(
        vec![si(), carbon()],
        vec![(box_mesh(a.0, a.1), 0), (box_mesh(b.0, b.1), 1)],
    )
    .unwrap();
    let csg = CsgGeometry::new(
        vec![si(), carbon()],
        vec![(Csg::cuboid(a.0, a.1), 0), (Csg::cuboid(b.0, b.1), 1)],
    )
    .unwrap();
    let mut rng = Rng(777);
    for _ in 0..3000 {
        let o = [2.0 * rng.next(), rng.next(), rng.next()];
        let d = rng.dir();
        let c = trace(&csg, o, d);
        assert!(same_trace(&c, &trace(&voxel, o, d)), "{o:?} {d:?}");
        assert!(same_trace(&c, &trace(&mesh, o, d)), "{o:?} {d:?}");
    }
    // Start on the shared face, either way: the same events as the grid and
    // the mesh.
    let p = [1.0, 0.3, 0.3];
    let rc = csg.locate(p).unwrap();
    assert_eq!(csg.material_index(rc), 1);
    for dx in [1.0, -1.0] {
        let d = [dx, 0.0, 0.0];
        let ec = csg.exit(rc, p, d, 1e9).unwrap();
        let ev = voxel.exit(voxel.locate(p).unwrap(), p, d, 1e9).unwrap();
        let em = mesh.exit(mesh.locate(p).unwrap(), p, d, 1e9).unwrap();
        assert!(close(ec.distance, ev.distance) && close(ec.distance, em.distance));
        let mat = |g: &dyn Geometry, e: &Exit| match e.outcome {
            ExitOutcome::Enter { region, .. } => Some(g.material_index(region)),
            ExitOutcome::Escape { .. } => None,
        };
        assert_eq!(mat(&csg, &ec), mat(&voxel, &ev));
        assert_eq!(mat(&csg, &ec), mat(&mesh, &em));
    }
}

// ----------------------------------------------------------------- rules ----

fn two_boxes(first_a: bool) -> CsgGeometry {
    // A = [0,2]^3 material 0; B = [1,3] x [0,2] x [0,2] material 1.
    let a = (Csg::cuboid([0.0; 3], [2.0; 3]), 0);
    let b = (Csg::cuboid([1.0, 0.0, 0.0], [3.0, 2.0, 2.0]), 1);
    let solids = if first_a { vec![a, b] } else { vec![b, a] };
    CsgGeometry::new(vec![si(), carbon()], solids).unwrap()
}

#[test]
fn rule_overlap_first_listed_owns() {
    let p = [1.5, 1.0, 1.0];
    let g = two_boxes(true);
    assert_eq!(g.material_index(g.locate(p).unwrap()), 0);
    let g2 = two_boxes(false);
    assert_eq!(g2.material_index(g2.locate(p).unwrap()), 1);
    // A first: B's face at x = 1 is not an event; at x = 2 B takes over.
    let e = g.exit(0, [0.5, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.5));
    assert_eq!(enter(&e).0, 1);
    // B first: the event is at x = 1, and B's end at x = 3 is an escape.
    let ra = g2.locate([0.5, 1.0, 1.0]).unwrap();
    let e = g2.exit(ra, [0.5, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 0.5));
    assert_eq!(g2.material_index(enter(&e).0), 1);
    let rb = g2.locate([2.5, 1.0, 1.0]).unwrap();
    let e = g2.exit(rb, [2.5, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 0.5));
    assert_eq!(escape(&e).0, Face::Back);
}

fn abutting() -> CsgGeometry {
    CsgGeometry::new(
        vec![si(), carbon()],
        vec![
            (Csg::cuboid([0.0; 3], [1.0; 3]), 0),
            (Csg::cuboid([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]), 1),
        ],
    )
    .unwrap()
}

#[test]
fn rule_surface_belongs_to_the_solid_being_entered() {
    let g = abutting();
    let on_face = [1.0, 0.5, 0.5];
    // `locate` on the shared face: the upper solid (lower faces own).
    assert_eq!(g.locate(on_face), Some(1));
    // On A's lower face it is A; on A's upper (outer) face, vacuum.
    assert_eq!(g.locate([0.0, 0.5, 0.5]), Some(0));
    assert_eq!(g.locate([2.0, 0.5, 0.5]), None);
    // Heading +x: stays in B, leaves at x = 2.
    let e = g.exit(1, on_face, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e).0, Face::Back);
    // Heading -x from the same point the flight is entering A: an
    // immediate event, with no vacuum between the solids.
    let e = g.exit(1, on_face, [-1.0, 0.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, 0.0);
    assert_eq!(enter(&e).0, 0);
    // From inside A through the face: one Enter, at the face.
    let e = g.exit(0, [0.25, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 0.75));
    let (r, p) = enter(&e);
    assert_eq!(r, 1);
    assert!(close(p[0], 1.0));
}

#[test]
fn rule_surface_just_crossed_is_not_hit_again() {
    let g = abutting();
    let mut rng = Rng(99);
    for _ in 0..500 {
        let mut d = rng.dir();
        d[0] = d[0].abs() + 0.2;
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let d = d.map(|c| c / l);
        let o = [0.5, 0.25 + 0.5 * rng.next(), 0.25 + 0.5 * rng.next()];
        let e = g.exit(0, o, d, 1e9).unwrap();
        if let ExitOutcome::Enter { region, pos } = e.outcome {
            let e2 = g.exit(region, pos, d, 1e9).unwrap();
            assert!(e2.distance > 1e-3, "{e2:?}");
            assert!(escape(&e2).1[0] != -1.0);
        }
    }
    // Crossings within `tol` of the start are ignored; just beyond, they
    // count.
    let b = one(Csg::cuboid([0.0; 3], [2.0; 3]));
    let tol = b.tolerance_m();
    assert!(close(tol, 1e-10 * 12f64.sqrt()));
    let e = b.exit(0, [0.0, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, 2.0);
    let e = b
        .exit(0, [2.0 - 0.5 * tol, 1.0, 1.0], [-1.0, 0.0, 0.0], 1e9)
        .unwrap();
    assert!(close(e.distance, 2.0 - 0.5 * tol));
    let e = b
        .exit(0, [2.0 - 3.0 * tol, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9)
        .unwrap();
    assert!(e.distance > tol && e.distance < 4.0 * tol, "{e:?}");
}

#[test]
fn rule_shared_face_passes_directly_and_coincident_surfaces_cancel() {
    // Two solids sharing a face: straight through, at an edge too.
    let g = CsgGeometry::new(
        vec![si(), carbon()],
        vec![
            (Csg::cuboid([0.0; 3], [2.0; 3]), 0),
            (Csg::cuboid([2.0, 0.0, 0.0], [4.0, 2.0, 2.0]), 1),
        ],
    )
    .unwrap();
    let c = [1.0, 1.0, 1.0];
    let e = g.exit(0, c, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    let (r, p) = enter(&e);
    assert_eq!(r, 1);
    let e = g.exit(1, p, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0));
    assert_eq!(escape(&e).0, Face::Back);
    let s = 0.5f64.sqrt();
    let e = g.exit(0, [1.0, 1.0, 1.0], [s, s, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2f64.sqrt()), "{e:?}");
    // Two halves of one solid united along a face: no event inside, even
    // with the second half cut a hair short of the first.
    let u = one(Csg::cuboid([0.0; 3], [1.0; 3])
        .union(Csg::cuboid([1.0 + 1e-14, 0.0, 0.0], [2.0, 1.0, 1.0])));
    let e = u.exit(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.5));
    // A cut that stops within `tol` of the far face leaves no sliver: the
    // flight out of the neighbouring solid escapes instead of entering it.
    let f = CsgGeometry::new(
        vec![si(), carbon()],
        vec![
            (
                Csg::cuboid([0.0; 3], [2.0; 3])
                    .minus(Csg::cuboid([-1.0; 3], [2.0 - 1e-14, 3.0, 3.0])),
                1,
            ),
            (Csg::cuboid([2.0, 0.0, 0.0], [4.0, 2.0, 2.0]), 0),
        ],
    )
    .unwrap();
    let e = f.exit(1, [3.0, 1.0, 1.0], [-1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e).0, Face::Front);
}

#[test]
fn same_material_solids_are_not_an_event() {
    let g = CsgGeometry::new(
        vec![si()],
        vec![
            (Csg::cuboid([0.0; 3], [1.0; 3]), 0),
            (Csg::cuboid([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]), 0),
        ],
    )
    .unwrap();
    let e = g.exit(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.5));
    assert_eq!(escape(&e).0, Face::Back);
    assert_eq!(
        g.flight(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1.0),
        Flight::Clear { region: 1 }
    );
    // A vacuum gap between same-material solids is an escape.
    let gap = CsgGeometry::new(
        vec![si()],
        vec![
            (Csg::cuboid([0.0; 3], [1.0; 3]), 0),
            (Csg::cuboid([1.5, 0.0, 0.0], [2.5, 1.0, 1.0]), 0),
        ],
    )
    .unwrap();
    let e = gap.exit(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 0.5));
    escape(&e);
}

// ------------------------------------------------------------- rejection ----

#[test]
fn unbounded_and_invalid_solids_are_rejected() {
    // z <= 0.5 (normalised on build).
    let hs = || Csg::half_space([0.0, 0.0, 2.0], 1.0);
    assert_eq!(
        CsgGeometry::new(vec![si()], vec![(hs(), 0)]),
        Err(GeometryError::CsgUnbounded { solid: 0 })
    );
    // Second solid unbounded: named by its index.
    assert_eq!(
        CsgGeometry::new(
            vec![si()],
            vec![(Csg::cuboid([0.0; 3], [1.0; 3]), 0), (hs(), 0)]
        ),
        Err(GeometryError::CsgUnbounded { solid: 1 })
    );
    // Bounded by an intersection: accepted (a half box).
    let g = one(Csg::cuboid([0.0; 3], [1.0; 3]).intersect(hs()));
    assert_eq!(g.bounds(), ([0.0; 3], [1.0, 1.0, 0.5]));
    let g = one(Csg::cuboid([0.0, 0.0, -1.0], [1.0; 3]).intersect(hs()));
    let e = g.exit(0, [0.5, 0.5, -0.5], [0.0, 0.0, 1.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e).1, [0.0, 0.0, 1.0]);
    let err = CsgGeometry::new(
        vec![si()],
        vec![(Csg::cylinder([0.0; 3], [0.0, 0.0, 1.0], 1.0, 0.0), 0)],
    )
    .unwrap_err();
    assert!(matches!(err, GeometryError::CsgInvalid { solid: 0, .. }));
    assert!(err.to_string().contains("height"), "{err}");
}

// ----------------------------------------------------------- determinism ----

#[test]
fn queries_are_deterministic_across_threads() {
    let build = || {
        CsgGeometry::new(
            vec![si(), carbon()],
            vec![
                (
                    Csg::cylinder([0.5, 0.5, -1.0], [0.1, 0.2, 1.0], 0.2, 3.0),
                    1,
                ),
                (
                    Csg::cuboid([0.0; 3], [1.0; 3])
                        .union(Csg::cuboid([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]))
                        .intersect(Csg::half_space([1.0, 1.0, 1.0], 2.8))
                        .minus(Csg::cylinder([1.2, 0.5, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0)),
                    0,
                ),
            ],
        )
        .unwrap()
    };
    let g = build();
    assert_eq!(g, build());
    let mut rng = Rng(2024);
    let rays: Vec<(V3, V3)> = (0..4000)
        .map(|_| ([0.2 + 0.6 * rng.next(), rng.next(), rng.next()], rng.dir()))
        .filter(|(o, _)| g.locate(*o).is_some())
        .collect();
    assert!(rays.len() > 1000);
    let run = |g: &CsgGeometry, rays: &[(V3, V3)]| -> Vec<Option<Exit>> {
        rays.iter()
            .map(|&(o, d)| g.exit(g.locate(o).unwrap(), o, d, 1e9))
            .collect()
    };
    let serial = run(&g, &rays);
    for threads in [2usize, 8] {
        let chunk = rays.len().div_ceil(threads);
        let parts: Vec<Vec<Option<Exit>>> = std::thread::scope(|s| {
            let handles: Vec<_> = rays
                .chunks(chunk)
                .map(|c| {
                    let g = &g;
                    s.spawn(move || run(g, c))
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let merged: Vec<_> = parts.into_iter().flatten().collect();
        assert_eq!(merged.len(), serial.len());
        for (a, b) in merged.iter().zip(&serial) {
            let (a, b) = (a.unwrap(), b.unwrap());
            assert_eq!(a.distance.to_bits(), b.distance.to_bits());
            assert_eq!(a, b);
        }
    }
}
