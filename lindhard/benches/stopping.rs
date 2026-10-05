//! Electronic stopping evaluation over an energy sweep: Lindhard-Scharff,
//! Oen-Robinson (impact-averaged and local), the equipartition split and
//! Bethe-Bloch.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use std::hint::black_box;

use lindhard::ion::stopping::bethe::BetheBloch;
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::mix::EquipartitionMix;
use lindhard::ion::stopping::oen_robinson::OenRobinson;
use lindhard::ion::stopping::{ElectronicStopping, Ion};

const N: usize = 256;

/// `N` log-spaced energies from `lo` to `hi` eV.
fn sweep(lo: f64, hi: f64) -> Vec<f64> {
    (0..N)
        .map(|i| lo * (hi / lo).powf(i as f64 / (N - 1) as f64))
        .collect()
}

fn eval<S: ElectronicStopping>(model: &S, ion: &Ion, z2: u8, es: &[f64]) -> f64 {
    let mut acc = 0.0;
    for &e in es {
        acc += model.stopping(ion, z2, black_box(e)).unwrap();
    }
    acc
}

fn bench_stopping(c: &mut Criterion) {
    let boron = Ion::new(5).unwrap();
    let proton = Ion::proton();
    let nuclear = sweep(1e1, 1e6); // eV: the keV-implant range
    let fast = sweep(1e5, 1e8); // eV: proton, Bethe regime

    let ls = LindhardScharff::new();
    let or = OenRobinson::new();
    let mix = EquipartitionMix::new();
    let bethe = BetheBloch::new();

    let mut g = c.benchmark_group("stopping_256_energies");
    g.throughput(Throughput::Elements(N as u64));
    g.bench_function("lindhard_scharff_B_Si", |b| {
        b.iter(|| eval(&ls, &boron, 14, &nuclear))
    });
    g.bench_function("oen_robinson_B_Si", |b| {
        b.iter(|| eval(&or, &boron, 14, &nuclear))
    });
    g.bench_function("equipartition_mix_B_Si", |b| {
        b.iter(|| eval(&mix, &boron, 14, &nuclear))
    });
    g.bench_function("bethe_bloch_p_Si", |b| {
        b.iter(|| eval(&bethe, &proton, 14, &fast))
    });
    g.bench_function("oen_robinson_local_loss_B_Si", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &e in &nuclear {
                acc += or
                    .local_loss(&boron, 14, black_box(e), black_box(1.0e-11))
                    .unwrap();
            }
            acc
        })
    });
    g.finish();
}

criterion_group!(benches, bench_stopping);
criterion_main!(benches);
