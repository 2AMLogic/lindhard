//! Electron engine throughput (#152): the matched problems of #150, electrons
//! at normal incidence into bulk Si and Cu at 1, 5 and 20 keV, with the full
//! physics lindhard runs there (Mott elastic, single-pole Penn inelastic,
//! Kieft-Bosch secondaries, step barrier, cutoff at the vacuum level; see
//! `electron_common`).
//!
//! * `electron_tables_<material>`: the cross-section table build, timed
//!   separately from the transport. The elastic and the inelastic table of
//!   the material on the grid of its 20 keV problem (10 eV up to the beam
//!   energy plus the inner potential, 20 points per decade), on one thread.
//! * `electron_transport_<problem>`: a fixed number of primaries through
//!   `Transport::run` with the full tally, on explicit rayon pools of 1
//!   thread, or of each count in `LINDHARD_BENCH_THREADS` (for example
//!   `1,2,4`). `thrpt` is electrons (primaries) per second. The tables are
//!   built once, outside the timing.
//!
//! The material data are SYNTHETIC unless `LINDHARD_BENCH_ELECTRON_INPUTS`
//! names the matched inputs written by `validation/oracles/bench_electron.py
//! --write-inputs` (`electron_common` module docs). The thread-scaling curve
//! with its determinism check is `examples/electron_scaling.rs`.

mod electron_common;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use std::hint::black_box;

use electron_common as ec;

fn pool(n: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .expect("thread pool")
}

fn thread_counts() -> Vec<usize> {
    std::env::var("LINDHARD_BENCH_THREADS")
        .ok()
        .map(|s| {
            s.split(',')
                .map(|x| x.trim().parse().expect("LINDHARD_BENCH_THREADS: 1,2,..."))
                .filter(|&n| n > 0)
                .collect()
        })
        .unwrap_or_else(|| vec![1])
}

fn bench_tables(c: &mut Criterion) {
    let one = pool(1);
    for p in ec::PROBLEMS.iter().filter(|p| p.energy_ev == 2.0e4) {
        let (r, _) = ec::problem_input(p, p.bench_histories).expect("electron input");
        let mut g = c.benchmark_group(format!("electron_tables_{}", p.symbol));
        g.sample_size(10);
        g.throughput(Throughput::Elements(r.table_energy_ev.len() as u64));
        for m in &r.materials {
            g.bench_function("elastic", |b| {
                b.iter(|| one.install(|| ec::elastic_table(black_box(&r), &m.material).unwrap()))
            });
            g.bench_function("inelastic", |b| {
                b.iter(|| one.install(|| ec::inelastic_table(black_box(&r), m).unwrap()))
            });
        }
        g.finish();
    }
}

fn bench_transport(c: &mut Criterion) {
    let counts = thread_counts();
    let pools: Vec<_> = counts.iter().map(|&n| pool(n)).collect();
    let one = pool(1);
    for p in &ec::PROBLEMS {
        let (r, _) = ec::problem_input(p, p.bench_histories).expect("electron input");
        let tables = one.install(|| ec::tables(&r)).expect("tables");
        let (t, proto) = ec::transport(&r, &tables).expect("transport");
        let n = p.bench_histories;
        let mut g = c.benchmark_group(format!("electron_transport_{}", p.id));
        g.sample_size(10);
        g.throughput(Throughput::Elements(n));
        for (&k, pool) in counts.iter().zip(&pools) {
            g.bench_function(format!("threads_{k}"), |b| {
                b.iter(|| pool.install(|| ec::run(&r, &t, &proto, black_box(n)).unwrap()))
            });
        }
        g.finish();
    }
}

criterion_group!(benches, bench_tables, bench_transport);
criterion_main!(benches);
