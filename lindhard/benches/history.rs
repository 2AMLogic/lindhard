//! End-to-end BCA benchmarks: single histories for three representative
//! problems, free-path conventions, full cascade vs ions only, tally
//! overhead, and ions/s throughput of 10^4 ions across thread counts.
//!
//! Throughput is reported by criterion as elements/s, i.e. ions/s.

mod common;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion, Throughput};
use std::hint::black_box;

use common::{beam, config, stack, table, tally_config, Problem, PROBLEMS};
use lindhard::ion::bca::{Bca, BcaConfig, MeanFreePath, SummaryTally};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::rng::stream;
use lindhard::tally::IonTally;

const SEED: u64 = 1;

type Tweak = fn(&mut BcaConfig);

/// Time `count` histories one after another on this thread, without the
/// parallel driver, so the per-history cost is not mixed with scheduling.
fn histories(bca: &Bca, count: u64) -> f64 {
    let mut t = SummaryTally::new(1e-9, 100);
    let mut acc = 0.0;
    for i in 0..count {
        let mut rng = stream(SEED, i);
        acc += bca.history(&mut t, &mut rng, i).unwrap().rest;
    }
    acc
}

fn bench_single_history(c: &mut Criterion) {
    let ls = LindhardScharff::new();
    let mut g = c.benchmark_group("history_100_ions_single_thread");
    g.throughput(Throughput::Elements(100));
    for p in &PROBLEMS {
        let st = stack(p);
        let bca = Bca::new(beam(p, 1), &st, config(SEED), &ls, table()).unwrap();
        g.bench_function(p.label, |b| b.iter(|| histories(&bca, black_box(100))));
    }
    g.finish();
}

/// Variants of the config on B 5 keV -> Si, 1000 ions each.
fn bench_variants(c: &mut Criterion) {
    let ls = LindhardScharff::new();
    let p = &PROBLEMS[0];
    let st = stack(p);
    let variants: [(&str, Tweak); 4] = [
        ("cascade_constant_path", |_| {}),
        ("ions_only_constant_path", |c| c.follow_recoils = false),
        ("cascade_energy_dependent_path", |c| {
            c.mean_free_path = MeanFreePath::EnergyDependent {
                min_cm_angle_rad: 0.01,
            }
        }),
        ("ions_only_energy_dependent_path", |c| {
            c.follow_recoils = false;
            c.mean_free_path = MeanFreePath::EnergyDependent {
                min_cm_angle_rad: 0.01,
            }
        }),
    ];
    let mut g = c.benchmark_group("cascade_and_free_path_B_5keV_Si_1000_ions");
    g.throughput(Throughput::Elements(1000));
    for (name, tweak) in variants {
        let mut cfg = config(SEED);
        tweak(&mut cfg);
        let bca = Bca::new(beam(p, 1), &st, cfg, &ls, table()).unwrap();
        g.bench_function(name, |b| b.iter(|| histories(&bca, black_box(1000))));
    }
    g.finish();
}

/// The cheap summary tally against the full `IonTally` (histograms, moments,
/// damage, escapes) on the same histories.
fn bench_tally_overhead(c: &mut Criterion) {
    let ls = LindhardScharff::new();
    let p = &PROBLEMS[0];
    let st = stack(p);
    let bca = Bca::new(beam(p, 1), &st, config(SEED), &ls, table()).unwrap();
    let species = bca.species_z();
    let mut g = c.benchmark_group("tally_overhead_B_5keV_Si_1000_ions");
    g.throughput(Throughput::Elements(1000));
    g.bench_function("summary_tally", |b| b.iter(|| histories(&bca, 1000)));
    g.bench_function("ion_tally", |b| {
        b.iter_batched(
            || IonTally::new(&st, &species, tally_config()).unwrap(),
            |mut t| {
                let mut acc = 0.0;
                for i in 0..1000 {
                    let mut rng = stream(SEED, i);
                    acc += bca.history(&mut t, &mut rng, i).unwrap().rest;
                }
                (t, acc)
            },
            BatchSize::LargeInput,
        )
    });
    g.finish();
}

/// Thread counts to report: 1, 2, 4, ... up to the machine's parallelism
/// (always including it).
fn thread_counts() -> Vec<usize> {
    let max = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut v: Vec<usize> = std::iter::successors(Some(1usize), |n| Some(n * 2))
        .take_while(|&n| n < max)
        .collect();
    v.push(max);
    v
}

fn run_pool(bca: &Bca, pool: &rayon::ThreadPool) -> u64 {
    pool.install(|| bca.run(|| SummaryTally::new(1e-9, 100)).unwrap().histories)
}

/// 10^4 ions per problem through `Bca::run` on explicit rayon pools (not an
/// environment variable), for thread scaling.
fn bench_throughput(c: &mut Criterion) {
    let ls = LindhardScharff::new();
    let counts = thread_counts();
    let pools: Vec<_> = counts
        .iter()
        .map(|&n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
        })
        .collect();
    for p in &PROBLEMS {
        throughput_one(c, p, &ls, &counts, &pools);
    }
}

fn throughput_one(
    c: &mut Criterion,
    p: &Problem,
    ls: &LindhardScharff,
    counts: &[usize],
    pools: &[rayon::ThreadPool],
) {
    const IONS: u64 = 10_000;
    let st = stack(p);
    let bca = Bca::new(beam(p, IONS), &st, config(SEED), ls, table()).unwrap();
    let mut g = c.benchmark_group(format!("throughput_{}_10k_ions", p.label));
    g.sample_size(10);
    g.throughput(Throughput::Elements(IONS));
    for (&n, pool) in counts.iter().zip(pools) {
        g.bench_function(format!("threads_{n}"), |b| b.iter(|| run_pool(&bca, pool)));
    }
    g.finish();
}

criterion_group!(
    benches,
    bench_single_history,
    bench_variants,
    bench_tally_overhead,
    bench_throughput
);
criterion_main!(benches);
