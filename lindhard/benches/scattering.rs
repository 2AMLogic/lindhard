//! Scattering-angle evaluation: table lookup vs direct Gauss-Mehler quadrature
//! vs the Biersack-Haggmark magic formula, on the same set of (eps, beta)
//! points for the ZBL universal potential.

use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{
    theta_magic, theta_quadrature, MagicConstants, ScatteringTable, TableSpec,
};

/// Deterministic log-uniform sample of (eps, beta) inside the table range.
fn sample_points(n: usize) -> Vec<(f64, f64)> {
    let mut st = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        st = st
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let eps = 1e-3 * 10f64.powf(5.0 * next());
            let beta = 1e-3 * 10f64.powf(4.5 * next());
            (eps, beta)
        })
        .collect()
}

fn bench_scattering(c: &mut Criterion) {
    let s = Screening::ZblUniversal;
    let pot = Potential::new(s, 14.0, 14.0);
    let spec = TableSpec {
        eps_min: 1e-3,
        eps_max: 1e2,
        beta_min: 1e-3,
        beta_max: 1e2,
        per_decade: 32,
    };
    let table = ScatteringTable::build(&pot, &spec);
    let k = MagicConstants::ZBL;
    let pts = sample_points(1024);

    let mut g = c.benchmark_group("theta_1024_points");
    g.bench_function("table_lookup", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(e, bb) in &pts {
                acc += table.theta(black_box(e), black_box(bb)).unwrap();
            }
            acc
        })
    });
    g.bench_function("magic_formula", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(e, bb) in &pts {
                acc += theta_magic(s, &k, black_box(e), black_box(bb));
            }
            acc
        })
    });
    g.bench_function("gauss_mehler_quadrature", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(e, bb) in &pts {
                acc += theta_quadrature(s, black_box(e), black_box(bb));
            }
            acc
        })
    });
    g.finish();
}

criterion_group!(benches, bench_scattering);
criterion_main!(benches);
