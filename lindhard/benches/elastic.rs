//! One (Z, E) radial-Dirac solve (all partial waves, `sigma_el` and
//! `sigma_tr1`) at 100 eV, 1 keV and 50 keV. The potential is a Thomas-Fermi
//! length Yukawa stand-in for Z = 29 (the Salvat et al. 1987 table is not in
//! the tree, see `docs/data-provenance.md`).

use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use lindhard::electron::elastic::{ElasticSolver, SolverOptions, Yukawa};

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

criterion_group!(benches, bench_elastic);
criterion_main!(benches);
