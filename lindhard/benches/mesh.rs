//! Ray queries on a box target: the triangle-mesh geometry (BVH) against the
//! equivalent voxel grid. The mesh box is 12 triangles; the icosphere is a
//! 5120-triangle surface where the BVH matters.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use std::hint::black_box;

use lindhard::geometry::{Boundary, Geometry, MeshGeometry, TriMesh, VoxelGrid};
use lindhard::material::Material;

fn box_soup(hi: [f64; 3]) -> Vec<[[f64; 3]; 3]> {
    let p = |x: f64, y: f64, z: f64| [x * hi[0], y * hi[1], z * hi[2]];
    let quads = [
        [p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)],
        [p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)],
        [p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)],
        [p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)],
        [p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)],
        [p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)],
    ];
    quads
        .iter()
        .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
        .collect()
}

/// Subdivided octahedron projected to a sphere of radius `r`.
fn sphere_soup(r: f64, levels: usize) -> Vec<[[f64; 3]; 3]> {
    let unit = |p: [f64; 3]| {
        let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        [r * p[0] / l, r * p[1] / l, r * p[2] / l]
    };
    let (px, nx) = ([1., 0., 0.], [-1., 0., 0.]);
    let (py, ny) = ([0., 1., 0.], [0., -1., 0.]);
    let (pz, nz) = ([0., 0., 1.], [0., 0., -1.]);
    let mut tris = vec![
        [px, py, pz],
        [py, nx, pz],
        [nx, ny, pz],
        [ny, px, pz],
        [py, px, nz],
        [nx, py, nz],
        [ny, nx, nz],
        [px, ny, nz],
    ];
    for _ in 0..levels {
        let mid = |a: [f64; 3], b: [f64; 3]| unit([a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
        tris = tris
            .iter()
            .flat_map(|&[a, b, c]| {
                let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
                [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
            })
            .collect();
    }
    tris.iter().map(|t| t.map(unit)).collect()
}

fn rays(n: usize, extent: [f64; 3]) -> Vec<([f64; 3], [f64; 3])> {
    let mut s = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let o = [
                extent[0] * (0.1 + 0.8 * next()),
                extent[1] * (0.1 + 0.8 * next()),
                extent[2] * (0.1 + 0.8 * next()),
            ];
            let d = [next() - 0.5, next() - 0.5, next() - 0.5];
            let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            (o, [d[0] / l, d[1] / l, d[2] / l])
        })
        .collect()
}

fn bench(c: &mut Criterion) {
    let si = || Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    let n = 4;
    let h = 1e-8;
    let extent = [n as f64 * h; 3];
    let voxel = VoxelGrid::new(
        vec![si()],
        [n; 3],
        [h; 3],
        vec![0; n * n * n],
        [Boundary::Vacuum; 3],
    )
    .unwrap();
    let mesh = MeshGeometry::new(
        vec![si()],
        vec![(TriMesh::from_triangles(&box_soup(extent)).unwrap(), 0)],
    )
    .unwrap();
    let sphere = MeshGeometry::new(
        vec![si()],
        vec![(
            TriMesh::from_triangles(&sphere_soup(0.5 * extent[0], 5)).unwrap(),
            0,
        )],
    )
    .unwrap();
    let r = rays(1024, extent);
    // The sphere is centred on the origin: start inside it.
    let centred: Vec<_> = r
        .iter()
        .map(|&(o, d)| {
            let c = 0.5 * extent[0];
            ([0.6 * (o[0] - c), 0.6 * (o[1] - c), 0.6 * (o[2] - c)], d)
        })
        .collect();

    let mut g = c.benchmark_group("mesh_exit");
    g.throughput(Throughput::Elements(r.len() as u64));
    for (name, geo, rs) in [
        ("voxel_box", &voxel as &dyn Geometry, &r),
        ("mesh_box", &mesh as &dyn Geometry, &r),
        ("mesh_sphere_5120", &sphere as &dyn Geometry, &centred),
    ] {
        let regions: Vec<usize> = rs.iter().map(|&(o, _)| geo.locate(o).unwrap()).collect();
        g.bench_function(name, |b| {
            b.iter(|| {
                let mut acc = 0.0;
                for (&(o, d), &reg) in rs.iter().zip(&regions) {
                    if let Some(e) = geo.exit(black_box(reg), black_box(o), black_box(d), 1.0) {
                        acc += e.distance;
                    }
                }
                black_box(acc)
            })
        });
    }
    g.finish();

    let mut g = c.benchmark_group("mesh_locate");
    g.throughput(Throughput::Elements(r.len() as u64));
    for (name, geo) in [
        ("voxel_box", &voxel as &dyn Geometry),
        ("mesh_box", &mesh as &dyn Geometry),
    ] {
        g.bench_function(name, |b| {
            b.iter(|| {
                let mut n = 0usize;
                for &(o, _) in &r {
                    n += usize::from(geo.locate(black_box(o)).is_some());
                }
                black_box(n)
            })
        });
    }
    g.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
