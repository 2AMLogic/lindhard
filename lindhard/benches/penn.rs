//! Penn models. Group 1: single-pole IMFP λ(E) and stopping power S(E) on a
//! 200-point log grid (10 eV to 50 keV) for one material. Group 2: the same
//! quantities on 8 energies for the single-pole and the full Penn algorithm.
//!
//! The material is the synthetic Drude plasmon of the tests (E_p = 20 eV,
//! γ = 5 eV; not physical data), tabulated at 1000 log-spaced energies from
//! 10 meV to 100 keV. The time is dominated by the outer integral over the
//! table, so it scales with the number of table knots below E.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use std::hint::black_box;

use lindhard::electron::inelastic::{
    DrudeLorentz, DrudeLorentzOscillator, FullPenn, SinglePolePenn,
};

const N: usize = 200;
/// Energies of the single-pole vs full comparison (the full model is three
/// nested adaptive integrals and costs about a second per energy).
const N_COMPARE: usize = 8;

fn bench_penn(c: &mut Criterion) {
    let elf = DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])
        .unwrap()
        .to_optical_elf("synthetic Drude plasmon", 1e-2, 1e5, 1000)
        .unwrap();
    let penn = SinglePolePenn::new(elf);
    let (lo, hi) = (10.0f64, 5.0e4f64);
    let grid: Vec<f64> = (0..N)
        .map(|i| lo * (hi / lo).powf(i as f64 / (N - 1) as f64))
        .collect();

    let mut g = c.benchmark_group("penn_spa_200_energies");
    g.sample_size(10);
    g.throughput(Throughput::Elements(N as u64));
    g.bench_function("imfp_and_stopping_drude_1000_knots", |b| {
        b.iter(|| penn.tabulate(black_box(&grid)).unwrap())
    });
    g.finish();

    // The same material and energies, single pole against the full Penn
    // algorithm (default tolerance 1e-4), 8 log-spaced energies, 10 eV to 50 keV.
    let few: Vec<f64> = (0..N_COMPARE)
        .map(|i| lo * (hi / lo).powf(i as f64 / (N_COMPARE - 1) as f64))
        .collect();
    let full = FullPenn::new(penn.optical_elf().clone());
    let mut g = c.benchmark_group("penn_single_pole_vs_full_8_energies");
    g.sample_size(10);
    g.throughput(Throughput::Elements(N_COMPARE as u64));
    g.bench_function("single_pole", |b| {
        b.iter(|| penn.tabulate(black_box(&few)).unwrap())
    });
    g.bench_function("full", |b| {
        b.iter(|| full.tabulate(black_box(&few)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, bench_penn);
criterion_main!(benches);
