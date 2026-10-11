//! Triangle-mesh solids (`lindhard::geometry::MeshGeometry`): exact chords on a
//! meshed box, tessellation-bounded chords on an icosphere, agreement with the
//! voxel grid, one test per ownership rule, rejection of bad meshes, format
//! equivalence and determinism. All meshes are generated here in code.

use lindhard::geometry::{
    Boundary, Exit, ExitOutcome, Face, Flight, Geometry, GeometryError, MeshGeometry, TriMesh,
    VoxelGrid,
};
use lindhard::material::Material;

type V3 = [f64; 3];

fn si() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn carbon() -> Material {
    Material::from_atom_fractions(&[(6, 1.0)], None).unwrap()
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
        [p(0, 0, 0), p(0, 0, 1), p(0, 1, 1), p(0, 1, 0)], // -x
        [p(1, 0, 0), p(1, 1, 0), p(1, 1, 1), p(1, 0, 1)], // +x
        [p(0, 0, 0), p(1, 0, 0), p(1, 0, 1), p(0, 0, 1)], // -y
        [p(0, 1, 0), p(0, 1, 1), p(1, 1, 1), p(1, 1, 0)], // +y
        [p(0, 0, 0), p(0, 1, 0), p(1, 1, 0), p(1, 0, 0)], // -z
        [p(0, 0, 1), p(1, 0, 1), p(1, 1, 1), p(0, 1, 1)], // +z
    ];
    quads
        .iter()
        .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
        .collect()
}

fn box_mesh(lo: V3, hi: V3) -> TriMesh {
    TriMesh::from_triangles(&box_soup(lo, hi)).unwrap()
}

fn one_box(lo: V3, hi: V3) -> MeshGeometry {
    MeshGeometry::new(vec![si()], vec![(box_mesh(lo, hi), 0)]).unwrap()
}

/// Icosphere of radius `r` about the origin, `levels` midpoint subdivisions.
fn icosphere_soup(r: f64, levels: usize) -> Vec<[V3; 3]> {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<V3> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[usize; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let unit = |p: V3| {
        let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        [r * p[0] / l, r * p[1] / l, r * p[2] / l]
    };
    for p in v.iter_mut() {
        *p = unit(*p);
    }
    for _ in 0..levels {
        let mut cache = std::collections::BTreeMap::new();
        let mut mid = |a: usize, b: usize, v: &mut Vec<V3>| -> usize {
            *cache.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let m = [
                    0.5 * (v[a][0] + v[b][0]),
                    0.5 * (v[a][1] + v[b][1]),
                    0.5 * (v[a][2] + v[b][2]),
                ];
                v.push(unit(m));
                v.len() - 1
            })
        };
        let mut g = Vec::new();
        for [a, b, c] in f {
            let (ab, bc, ca) = (mid(a, b, &mut v), mid(b, c, &mut v), mid(c, a, &mut v));
            g.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        f = g;
    }
    f.iter().map(|t| t.map(|i| v[i])).collect()
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

// ---------------------------------------------------------------- boxes ----

#[test]
fn meshed_box_gives_exact_chords() {
    let (lo, hi) = ([0.0, 0.0, 0.0], [3.0, 2.0, 1.0]);
    let g = one_box(lo, hi);
    assert_eq!(g.n_regions(), 1);
    assert_eq!(g.n_materials(), 1);
    assert_eq!(g.material_index(0), 0);
    assert_eq!(g.entry_point(), [0.0, 1.0, 0.5]);
    // Axis-aligned.
    let e = g.exit(0, [0.0, 1.0, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 3.0));
    assert_eq!(escape(&e), (Face::Back, [1.0, 0.0, 0.0]));
    let e = g.exit(0, [1.0, 1.0, 0.5], [-1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e), (Face::Front, [-1.0, 0.0, 0.0]));
    let e = g.exit(0, [1.0, 1.0, 0.25], [0.0, 1.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e), (Face::Side, [0.0, 1.0, 0.0]));
    // Oblique, against the slab formula.
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for _ in 0..2000 {
        let o = [3.0 * rng.next(), 2.0 * rng.next(), rng.next()];
        let d = rng.dir();
        let e = g.exit(0, o, d, 1e9).unwrap();
        assert!(close(e.distance, box_chord(o, d, lo, hi)), "{o:?} {d:?}");
        let (_, n) = escape(&e);
        assert_eq!(n.iter().map(|c| c.abs()).sum::<f64>(), 1.0);
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
}

#[test]
fn locate_and_vacuum() {
    let g = one_box([0.0; 3], [2.0; 3]);
    assert_eq!(g.locate([1.0, 1.0, 1.0]), Some(0));
    assert_eq!(g.locate([3.0, 1.0, 1.0]), None);
    assert_eq!(g.locate([-1e-3, 1.0, 1.0]), None);
    assert_eq!(g.locate([1.0, 1.0, f64::NAN]), None);
    // Voxel convention: lower faces belong to the solid, upper ones do not.
    assert_eq!(g.locate([0.0, 1.0, 1.0]), Some(0));
    assert_eq!(g.locate([0.0, 0.0, 0.0]), Some(0));
    assert_eq!(g.locate([2.0, 1.0, 1.0]), None);
    assert!(g.in_vacuum([5.0, 0.0, 0.0]));
    assert!(!g.in_vacuum([1.0, 1.0, 1.0]));
}

#[test]
fn mesh_box_matches_voxel_grid() {
    let h = 0.5;
    let dims = [4, 2, 2];
    let extent = [4.0 * h, 2.0 * h, 2.0 * h];
    let voxel =
        VoxelGrid::new(vec![si()], dims, [h; 3], vec![0; 16], [Boundary::Vacuum; 3]).unwrap();
    let mesh = one_box([0.0; 3], extent);
    assert_eq!(mesh.entry_point(), voxel.entry_point());
    let mut rng = Rng(12345);
    for _ in 0..3000 {
        let o = [
            extent[0] * rng.next(),
            extent[1] * rng.next(),
            extent[2] * rng.next(),
        ];
        let d = rng.dir();
        let (rv, rm) = (voxel.locate(o).unwrap(), mesh.locate(o).unwrap());
        let (ev, em) = (
            voxel.exit(rv, o, d, 1e9).unwrap(),
            mesh.exit(rm, o, d, 1e9).unwrap(),
        );
        assert!(close(ev.distance, em.distance));
        for k in 0..3 {
            assert!(close(ev.at[k], em.at[k]));
        }
        assert_eq!(escape(&ev), escape(&em));
    }
    // Points on and beside faces agree in `locate` (lower faces own, upper
    // faces do not) and in `in_vacuum`.
    for x in [0.0, 0.5, 1.0, 2.0, 2.5, -0.5] {
        for y in [0.0, 0.5, 1.0, 1.5] {
            for z in [0.0, 0.25, 1.0] {
                let p = [x, y, z];
                assert_eq!(voxel.locate(p).is_some(), mesh.locate(p).is_some(), "{p:?}");
                assert_eq!(voxel.in_vacuum(p), mesh.in_vacuum(p), "{p:?}");
            }
        }
    }
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

#[test]
fn two_material_mesh_matches_voxel_grid() {
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
    let mesh = MeshGeometry::new(
        vec![si(), carbon()],
        vec![
            (box_mesh([0.0; 3], [1.0, 1.0, 1.0]), 0),
            (box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]), 1),
        ],
    )
    .unwrap();
    let mut rng = Rng(777);
    for _ in 0..3000 {
        let o = [2.0 * rng.next(), rng.next(), rng.next()];
        let d = rng.dir();
        let (a, b) = (trace(&voxel, o, d), trace(&mesh, o, d));
        assert_eq!(a.len(), b.len(), "{o:?} {d:?}");
        for (x, y) in a.iter().zip(&b) {
            assert_eq!((x.0, x.3), (y.0, y.3));
            assert!(close(x.1, y.1), "{x:?} {y:?}");
            for k in 0..3 {
                assert!(close(x.2[k], y.2[k]));
            }
        }
    }
    // Start on the shared face, either way: the same events as the grid.
    for dx in [1.0, -1.0] {
        let p = [1.0, 0.3, 0.3];
        let (rv, rm) = (voxel.locate(p).unwrap(), mesh.locate(p).unwrap());
        assert_eq!(voxel.material_index(rv), mesh.material_index(rm));
        let d = [dx, 0.0, 0.0];
        let (ev, em) = (voxel.exit(rv, p, d, 1e9), mesh.exit(rm, p, d, 1e9));
        let (ev, em) = (ev.unwrap(), em.unwrap());
        assert!(close(ev.distance, em.distance), "{ev:?} {em:?}");
        let mat = |g: &dyn Geometry, e: &Exit| match e.outcome {
            ExitOutcome::Enter { region, .. } => Some(g.material_index(region)),
            ExitOutcome::Escape { .. } => None,
        };
        assert_eq!(mat(&voxel, &ev), mat(&mesh, &em));
    }
}

// ------------------------------------------------------------- sphere ----

#[test]
fn icosphere_chords_within_the_tessellation_bound() {
    let r = 1.0;
    let soup = icosphere_soup(r, 2);
    // Largest circumradius of a face: its plane is at sqrt(r^2 - rho^2).
    let rho_max = soup
        .iter()
        .map(|t| {
            let sub = |a: V3, b: V3| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let len = |a: V3| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
            let (ab, bc, ca) = (
                len(sub(t[1], t[0])),
                len(sub(t[2], t[1])),
                len(sub(t[0], t[2])),
            );
            let u = sub(t[1], t[0]);
            let v = sub(t[2], t[0]);
            let c = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            ab * bc * ca / (2.0 * len(c))
        })
        .fold(0.0, f64::max);
    let r_in = (r * r - rho_max * rho_max).sqrt();
    assert!(r_in < r && r_in > 0.9 * r);
    let mesh = TriMesh::from_triangles(&soup).unwrap();
    // Volume of an inscribed polyhedron: below the sphere, above the insphere.
    let sphere_v = 4.0 / 3.0 * std::f64::consts::PI;
    assert!(mesh.volume_m3() < sphere_v && mesh.volume_m3() > sphere_v * r_in.powi(3));
    let g = MeshGeometry::new(vec![si()], vec![(mesh, 0)]).unwrap();
    // True distance to a sphere of radius `s` about the origin.
    let to_sphere = |o: V3, d: V3, s: f64| {
        let b = o[0] * d[0] + o[1] * d[1] + o[2] * d[2];
        let c = o[0] * o[0] + o[1] * o[1] + o[2] * o[2] - s * s;
        -b + (b * b - c).sqrt()
    };
    let mut rng = Rng(4242);
    for _ in 0..2000 {
        let o0 = rng.dir();
        let o = o0.map(|c| c * 0.9 * r_in * rng.next());
        let d = rng.dir();
        let e = g.exit(0, o, d, 1e9).unwrap();
        let (lo, hi) = (to_sphere(o, d, r_in), to_sphere(o, d, r));
        assert!(
            e.distance >= lo * (1.0 - 1e-12) && e.distance <= hi * (1.0 + 1e-12),
            "{} not in [{lo}, {hi}]",
            e.distance
        );
    }
    // A chord through the centre, along the axes.
    for axis in 0..3 {
        let mut d = [0.0; 3];
        d[axis] = 1.0;
        let e = g.exit(0, [0.0; 3], d, 1e9).unwrap();
        assert!(e.distance >= r_in && e.distance <= r);
    }
}

// -------------------------------------------------------------- rules ----

fn two_boxes(first_a: bool) -> MeshGeometry {
    // A = [0,2]^3 material 0; B = [1,3] x [0,2] x [0,2] material 1.
    let a = (box_mesh([0.0; 3], [2.0; 3]), 0);
    let b = (box_mesh([1.0, 0.0, 0.0], [3.0, 2.0, 2.0]), 1);
    let solids = if first_a { vec![a, b] } else { vec![b, a] };
    MeshGeometry::new(vec![si(), carbon()], solids).unwrap()
}

#[test]
fn rule_overlap_first_listed_owns() {
    let p = [1.5, 1.0, 1.0];
    let g = two_boxes(true);
    assert_eq!(g.material_index(g.locate(p).unwrap()), 0);
    let g2 = two_boxes(false);
    assert_eq!(g2.material_index(g2.locate(p).unwrap()), 1);
    // A flight through A first: B's face at x = 1 is not an event (A owns the
    // overlap); at x = 2 A ends and B takes over.
    let e = g.exit(0, [0.5, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.5));
    assert_eq!(enter(&e).0, 1);
    // Listed the other way, B owns from x = 1: the event is at x = 1.
    let a_region = g2.locate([0.5, 1.0, 1.0]).unwrap();
    assert_eq!(g2.material_index(a_region), 0);
    let e = g2
        .exit(a_region, [0.5, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9)
        .unwrap();
    assert!(close(e.distance, 0.5));
    assert_eq!(g2.material_index(enter(&e).0), 1);
    // ... and B's end at x = 3 is an escape.
    let b_region = g2.locate([2.5, 1.0, 1.0]).unwrap();
    let e = g2
        .exit(b_region, [2.5, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9)
        .unwrap();
    assert!(close(e.distance, 0.5));
    assert_eq!(escape(&e).0, Face::Back);
}

#[test]
fn rule_surface_belongs_to_the_solid_being_entered() {
    // Abutting boxes: A = [0,1]^3, B = [1,2] x [0,1]^2.
    let g = MeshGeometry::new(
        vec![si(), carbon()],
        vec![
            (box_mesh([0.0; 3], [1.0; 3]), 0),
            (box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]), 1),
        ],
    )
    .unwrap();
    let on_face = [1.0, 0.5, 0.5];
    // Standing on the shared face, region B (what `locate` says).
    assert_eq!(g.locate(on_face), Some(1));
    // Heading +x: stays in B, leaves at x = 2.
    let e = g.exit(1, on_face, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e).0, Face::Back);
    // Heading -x from the same point the ray is entering A: an immediate
    // event, with no vacuum between the solids.
    let e = g.exit(1, on_face, [-1.0, 0.0, 0.0], 1e9).unwrap();
    assert_eq!(e.distance, 0.0);
    assert_eq!(enter(&e).0, 0);
    // From inside A through the face: one Enter, at the right place.
    let e = g.exit(0, [0.25, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 0.75));
    let (r, p) = enter(&e);
    assert_eq!(r, 1);
    assert!(close(p[0], 1.0));
}

#[test]
fn rule_surface_just_crossed_is_not_hit_again() {
    let g = MeshGeometry::new(
        vec![si(), carbon()],
        vec![
            (box_mesh([0.0; 3], [1.0; 3]), 0),
            (box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]), 1),
        ],
    )
    .unwrap();
    let mut rng = Rng(99);
    for _ in 0..500 {
        let mut d = rng.dir();
        d[0] = d[0].abs() + 0.2;
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let d = d.map(|c| c / l);
        let o = [0.5, 0.25 + 0.5 * rng.next(), 0.25 + 0.5 * rng.next()];
        let e = g.exit(0, o, d, 1e9).unwrap();
        if let ExitOutcome::Enter { region, pos } = e.outcome {
            // Re-entered from the crossing point (rounded, on the face): the
            // next event is the far side of B, never the face just crossed.
            let e2 = g.exit(region, pos, d, 1e9).unwrap();
            assert!(e2.distance > 1e-3, "{e2:?}");
            assert!(escape(&e2).1[0] != -1.0);
        }
    }
    // A start on the exit face of a single box heading outward-adjacent
    // inward ignores the face it stands on.
    let b = one_box([0.0; 3], [2.0; 3]);
    let e = b.exit(0, [0.0, 1.0, 1.0], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0));
    let e = b
        .exit(0, [2.0 - 1e-15, 1.0, 1.0], [-1.0, 0.0, 0.0], 1e9)
        .unwrap();
    assert!(close(e.distance, 2.0));
}

#[test]
fn rule_shared_edge_or_vertex_crosses_once() {
    let g = one_box([0.0; 3], [2.0; 3]);
    let c = [1.0, 1.0, 1.0];
    // Through the diagonal edge of the +x and -x face triangulations.
    let e = g.exit(0, c, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e), (Face::Back, [1.0, 0.0, 0.0]));
    let e = g.exit(0, c, [-1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(escape(&e).0, Face::Front);
    // Through a box edge (two faces meet) and a vertex (three faces).
    let s = 0.5f64.sqrt();
    let e = g.exit(0, c, [s, s, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2f64.sqrt()), "{e:?}");
    escape(&e);
    let t = 1.0 / 3f64.sqrt();
    let e = g.exit(0, c, [t, t, t], 1e9).unwrap();
    assert!(close(e.distance, 3f64.sqrt()), "{e:?}");
    escape(&e);
    // The same through a vertex of the back and front sides.
    let e = g.exit(0, c, [-t, -t, -t], 1e9).unwrap();
    assert!(close(e.distance, 3f64.sqrt()));
    assert_eq!(escape(&e).0, Face::Side);
    // An icosphere vertex: the ray passes through a vertex shared by five
    // triangles, and leaves once.
    let soup = icosphere_soup(1.0, 1);
    let v = soup[0][0];
    let sph = MeshGeometry::new(
        vec![si()],
        vec![(TriMesh::from_triangles(&soup).unwrap(), 0)],
    )
    .unwrap();
    let e = sph.exit(0, [0.0; 3], v, 1e9).unwrap();
    assert!(close(e.distance, 1.0), "{e:?}");
    // Two solids sharing a face: the crossing at a shared diagonal passes
    // straight through.
    let two = MeshGeometry::new(
        vec![si(), carbon()],
        vec![
            (box_mesh([0.0; 3], [2.0; 3]), 0),
            (box_mesh([2.0, 0.0, 0.0], [4.0, 2.0, 2.0]), 1),
        ],
    )
    .unwrap();
    let e = two.exit(0, c, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.0));
    assert_eq!(enter(&e).0, 1);
    let e = two.exit(1, enter(&e).1, [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 2.0));
    assert_eq!(escape(&e).0, Face::Back);
}

#[test]
fn same_material_solids_are_not_an_event() {
    let g = MeshGeometry::new(
        vec![si()],
        vec![
            (box_mesh([0.0; 3], [1.0; 3]), 0),
            (box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]), 0),
        ],
    )
    .unwrap();
    let e = g.exit(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 1.5));
    assert_eq!(escape(&e).0, Face::Back);
    // A clear flight reports the solid it ends in.
    assert_eq!(
        g.flight(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1.0),
        Flight::Clear { region: 1 }
    );
    // A vacuum gap between same-material solids is an escape.
    let gap = MeshGeometry::new(
        vec![si()],
        vec![
            (box_mesh([0.0; 3], [1.0; 3]), 0),
            (box_mesh([1.5, 0.0, 0.0], [2.5, 1.0, 1.0]), 0),
        ],
    )
    .unwrap();
    let e = gap.exit(0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1e9).unwrap();
    assert!(close(e.distance, 0.5));
    escape(&e);
}

// ----------------------------------------------------------- rejection ----

fn corners(e: &[[f64; 3]; 2], soup: &[[V3; 3]]) -> bool {
    let known = |p: &V3| soup.iter().flatten().any(|q| q == p);
    known(&e[0]) && known(&e[1]) && e[0] != e[1]
}

#[test]
fn open_mesh_is_rejected_naming_an_edge() {
    let mut soup = box_soup([0.0; 3], [1.0; 3]);
    soup.pop();
    match TriMesh::from_triangles(&soup) {
        Err(GeometryError::MeshEdgeNotShared { edge }) => {
            assert!(corners(&edge, &soup));
            let msg = GeometryError::MeshEdgeNotShared { edge }.to_string();
            assert!(msg.contains("edge"), "{msg}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn inconsistent_orientation_is_rejected() {
    let mut soup = box_soup([0.0; 3], [1.0; 3]);
    soup[3].swap(0, 1);
    match TriMesh::from_triangles(&soup) {
        Err(GeometryError::MeshInconsistentOrientation { edge }) => {
            assert!(corners(&edge, &soup));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn inside_out_mesh_is_rejected() {
    let soup: Vec<[V3; 3]> = box_soup([0.0; 3], [1.0; 3])
        .into_iter()
        .map(|[a, b, c]| [a, c, b])
        .collect();
    match TriMesh::from_triangles(&soup) {
        Err(GeometryError::MeshNotOutward { volume_m3 }) => assert!(volume_m3 < 0.0),
        other => panic!("{other:?}"),
    }
}

#[test]
fn degenerate_and_bad_input_is_rejected() {
    let mut soup = box_soup([0.0; 3], [1.0; 3]);
    let p = soup[0][0];
    soup.push([p, p, [5.0, 5.0, 5.0]]);
    assert!(matches!(
        TriMesh::from_triangles(&soup),
        Err(GeometryError::MeshDegenerateTriangle { triangle: 12, .. })
    ));
    // Collinear.
    let mut soup = box_soup([0.0; 3], [1.0; 3]);
    soup[5] = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
    assert!(matches!(
        TriMesh::from_triangles(&soup),
        Err(GeometryError::MeshDegenerateTriangle { triangle: 5, .. })
    ));
    let mut soup = box_soup([0.0; 3], [1.0; 3]);
    soup[2][1][1] = f64::NAN;
    assert!(matches!(
        TriMesh::from_triangles(&soup),
        Err(GeometryError::MeshNonFinite { triangle: 2 })
    ));
    assert_eq!(TriMesh::from_triangles(&[]), Err(GeometryError::MeshEmpty));
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            TriMesh::from_obj("v 0 0 0\n", bad),
            Err(GeometryError::MeshInvalidScale { .. })
        ));
    }
    // Geometry-level validation.
    let m = box_mesh([0.0; 3], [1.0; 3]);
    assert_eq!(
        MeshGeometry::new(vec![], vec![(m.clone(), 0)]),
        Err(GeometryError::NoMaterials)
    );
    assert_eq!(
        MeshGeometry::new(vec![si()], vec![]),
        Err(GeometryError::NoSolids)
    );
    assert_eq!(
        MeshGeometry::new(vec![si()], vec![(m.clone(), 0), (m, 3)]),
        Err(GeometryError::SolidMaterialOutOfRange {
            solid: 1,
            index: 3,
            n_materials: 1
        })
    );
}

#[test]
fn parse_errors_name_the_line() {
    let e = TriMesh::from_obj("v 0 0 0\nv 1 0 0\nf 1 2 9\n", 1.0).unwrap_err();
    assert!(
        matches!(
            e,
            GeometryError::MeshParse {
                format: "OBJ",
                line: 3,
                ..
            }
        ),
        "{e:?}"
    );
    let e = TriMesh::from_obj("v 0 0 0\nf 0 1 1\n", 1.0).unwrap_err();
    assert!(
        matches!(e, GeometryError::MeshParse { line: 2, .. }),
        "{e:?}"
    );
    let e = TriMesh::from_obj("v 0 0\n", 1.0).unwrap_err();
    assert!(
        matches!(e, GeometryError::MeshParse { line: 1, .. }),
        "{e:?}"
    );
    let e = TriMesh::from_obj("v 0 0 0\nv 1 0 0\nf 1 2\n", 1.0).unwrap_err();
    assert!(
        matches!(e, GeometryError::MeshParse { line: 3, .. }),
        "{e:?}"
    );
    let e = TriMesh::from_stl(b"hello", 1.0).unwrap_err();
    assert!(
        matches!(e, GeometryError::MeshParse { format: "STL", .. }),
        "{e:?}"
    );
    let e = TriMesh::from_stl(
        b"solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nendloop\nendfacet\nendsolid",
        1.0,
    )
    .unwrap_err();
    assert!(
        matches!(e, GeometryError::MeshParse { format: "STL", .. }),
        "{e:?}"
    );
}

// ------------------------------------------------------------- formats ----

fn stl_ascii(soup: &[[V3; 3]]) -> String {
    let mut s = String::from("solid test\n");
    for t in soup {
        s.push_str("  facet normal 0 0 0\n    outer loop\n");
        for v in t {
            s.push_str(&format!("      vertex {} {} {}\n", v[0], v[1], v[2]));
        }
        s.push_str("    endloop\n  endfacet\n");
    }
    s.push_str("endsolid test\n");
    s
}

fn stl_binary(soup: &[[V3; 3]]) -> Vec<u8> {
    // A header that starts with `solid`, as many binary writers do.
    let mut b = b"solid binary header".to_vec();
    b.resize(80, 0);
    b.extend((soup.len() as u32).to_le_bytes());
    for t in soup {
        b.extend([0u8; 12]);
        for v in t {
            for c in v {
                b.extend((*c as f32).to_le_bytes());
            }
        }
        b.extend([0u8; 2]);
    }
    b
}

fn obj_text(soup: &[[V3; 3]], style: usize) -> String {
    let mut verts: Vec<V3> = Vec::new();
    let mut faces = Vec::new();
    for t in soup {
        let idx: Vec<usize> = t
            .iter()
            .map(|v| {
                verts.iter().position(|w| w == v).unwrap_or_else(|| {
                    verts.push(*v);
                    verts.len() - 1
                }) + 1
            })
            .collect();
        faces.push(idx);
    }
    let mut s = String::from("# test\no cube\n");
    for v in &verts {
        s.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
    }
    s.push_str("vt 0 0\nvn 0 0 1\n");
    for f in faces {
        let tok = |i: usize| match style {
            0 => format!("{i}"),
            1 => format!("{i}/1/1"),
            2 => format!("{i}//1"),
            _ => format!("{}", i as i64 - verts.len() as i64 - 1), // negative
        };
        s.push_str(&format!(
            "f {} {} {}  # comment\n",
            tok(f[0]),
            tok(f[1]),
            tok(f[2])
        ));
    }
    s
}

#[test]
fn stl_ascii_binary_and_obj_agree() {
    let soup = box_soup([0.0, 0.0, 0.0], [2.0, 1.0, 0.5]);
    let reference = TriMesh::from_triangles(&soup).unwrap();
    assert_eq!(reference.vertices().len(), 8);
    assert_eq!(reference.triangles().len(), 12);
    assert!(close(reference.volume_m3(), 1.0));
    assert_eq!(reference.bounds(), ([0.0; 3], [2.0, 1.0, 0.5]));
    assert_eq!(
        TriMesh::from_stl(stl_ascii(&soup).as_bytes(), 1.0).unwrap(),
        reference
    );
    assert_eq!(
        TriMesh::from_stl(&stl_binary(&soup), 1.0).unwrap(),
        reference
    );
    for style in 0..4 {
        assert_eq!(
            TriMesh::from_obj(&obj_text(&soup, style), 1.0).unwrap(),
            reference,
            "style {style}"
        );
    }
}

#[test]
fn obj_polygons_and_units() {
    // Six quads, fan-triangulated, in nm.
    let mut s = String::new();
    for x in [0, 1] {
        for y in [0, 1] {
            for z in [0, 1] {
                s.push_str(&format!("v {x} {y} {z}\n"));
            }
        }
    }
    // Vertex order: index = 1 + 4x + 2y + z.
    for q in [
        [1, 2, 4, 3], // -x : (0,0,0) (0,0,1) (0,1,1) (0,1,0)
        [5, 7, 8, 6], // +x
        [1, 5, 6, 2], // -y
        [3, 4, 8, 7], // +y
        [1, 3, 7, 5], // -z
        [2, 6, 8, 4], // +z
    ] {
        s.push_str(&format!("f {} {} {} {}\n", q[0], q[1], q[2], q[3]));
    }
    let nm = 1e-9;
    let mesh = TriMesh::from_obj(&s, nm).unwrap();
    assert_eq!(mesh.triangles().len(), 12);
    assert!(close(mesh.volume_m3() / (nm * nm * nm), 1.0));
    assert_eq!(mesh.bounds().1, [nm; 3]);
    // Negative indices on quads.
    let mut t = s.clone();
    t.truncate(t.find("f ").unwrap());
    t.push_str("f 1 2 4 3\nf 5 7 8 6\nf 1 5 6 2\nf 3 4 8 7\nf 1 3 7 5\nf 2 6 8 4\n");
    assert_eq!(TriMesh::from_obj(&t, nm).unwrap(), mesh);
    let neg = t.replace("f 1 2 4 3", "f -8 -7 -5 -6");
    assert_eq!(TriMesh::from_obj(&neg, nm).unwrap(), mesh);
    // The same solid as a geometry gives the same chord whatever the format.
    let g = MeshGeometry::new(vec![si()], vec![(mesh, 0)]).unwrap();
    let e = g
        .exit(0, [0.0, 0.5 * nm, 0.5 * nm], [1.0, 0.0, 0.0], 1.0)
        .unwrap();
    assert!(close(e.distance / nm, 1.0));
}

// --------------------------------------------------------- determinism ----

#[test]
fn build_and_queries_are_deterministic() {
    let soup = icosphere_soup(1.0, 3);
    let build = || {
        MeshGeometry::new(
            vec![si(), carbon()],
            vec![
                (TriMesh::from_triangles(&soup).unwrap(), 0),
                (box_mesh([0.5, -0.5, -0.5], [2.0, 0.5, 0.5]), 1),
            ],
        )
        .unwrap()
    };
    let g = build();
    assert_eq!(g, build());

    let mut rng = Rng(2024);
    let rays: Vec<(V3, V3)> = (0..4000)
        .map(|_| {
            (
                [0.3 * rng.next(), 0.3 * rng.next(), 0.3 * rng.next()],
                rng.dir(),
            )
        })
        .collect();
    let run = |g: &MeshGeometry, rays: &[(V3, V3)]| -> Vec<Option<Exit>> {
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
        // Bit-identical, not merely close.
        assert_eq!(merged.len(), serial.len());
        for (a, b) in merged.iter().zip(&serial) {
            let (a, b) = (a.unwrap(), b.unwrap());
            assert_eq!(a.distance.to_bits(), b.distance.to_bits());
            assert_eq!(a, b);
        }
    }
}
