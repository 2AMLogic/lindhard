//! One (Z, E) radial-Dirac solve (all partial waves, `sigma_el` and
//! `sigma_tr1`) at 100 eV, 1 keV and 50 keV, and the full elastic table build
//! (`electron::elastic::table`, default 75-point 10 eV to 50 keV grid,
//! adaptive probability grid) for Si and Au. The potential is a Thomas-Fermi
//! length Yukawa stand-in (the Salvat et al. 1987 table is not in the tree,
//! see `docs/data-provenance.md`).

use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use lindhard::electron::elastic::table::{
    build_elastic_table, ElasticTableOptions, ThomasFermiYukawa,
};
use lindhard::electron::elastic::{ElasticSolver, SolverOptions, Yukawa};
use lindhard::material::Material;

fn bench_elastic(c: &mut Criterion) {
    let z = 29.0_f64;
    let pot = Yukawa::new(z, 0.8853 / z.cbrt()).expect("valid potential");
    let mut group = c.benchmark_group("elastic_solve");
    group.sample_size(10);
    for (name, e_ev) in [("100eV", 100.0), ("1keV", 1.0e3), ("50keV", 5.0e4)] {
        group.bench_function(name, |b| {
            b.iter(|| {
                let solver = ElasticSolver::new(&pot, black_box(e_ev), SolverOptions::default())
                    .expect("solver");
                let pw = solver.partial_waves().expect("converged");
                black_box((pw.sigma_el(), pw.sigma_tr1()))
            })
        });
    }
    group.finish();
}

fn bench_elastic_table(c: &mut Criterion) {
    let opts = ElasticTableOptions::default();
    let mut group = c.benchmark_group("elastic_table");
    group.sample_size(10);
    for (name, z) in [("Si", 14u8), ("Au", 79)] {
        let m = Material::from_atom_fractions(&[(z, 1.0)], None).expect("tabulated density");
        group.bench_function(name, |b| {
            b.iter(|| {
                black_box(
                    build_elastic_table(black_box(&m), &ThomasFermiYukawa, &opts).expect("table"),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_elastic, bench_elastic_table);
criterion_main!(benches);
