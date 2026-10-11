//! Crystal against amorphous flight model: ions per second for B 5 keV into
//! silicon (issue #180), the baseline for the performance risk of the
//! crystalline-implant epic. Same ion, energy, stopping model, cutoffs and
//! full recoil cascades in all cases; only the collision-partner model
//! differs. Single-threaded histories (the parallel driver is the same for
//! both), 7 degrees tilt and 22 degrees twist (a "random" implant direction)
//! and exactly along <110> (a channeled beam, where the ions travel far and
//! the search runs the most segments).
//!
//! Criterion reports elements/s, i.e. ions/s. The ratio is read off the
//! group's three lines.

mod common;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use std::hint::black_box;

use common::{config, elemental, table};
use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, Beam, CrystalTarget, SummaryTally};
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::rng::stream;

const SEED: u64 = 1;
const IONS: u64 = 200;

fn beam(tilt_deg: f64, twist_deg: f64) -> Beam {
    Beam {
        ion: Ion::new(5).unwrap(),
        energy_ev: 5.0e3,
        polar_rad: tilt_deg.to_radians(),
        azimuth_rad: twist_deg.to_radians(),
        count: 1,
    }
}

fn crystal(tilt_deg: f64, twist_deg: f64) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(
        &lat,
        [1, 0, 0],
        [0, 1, 0],
        tilt_deg.to_radians(),
        twist_deg.to_radians(),
        0.0,
    )
    .unwrap();
    CrystalTarget::new(lat, o)
}

fn histories(bca: &Bca, count: u64) -> f64 {
    let mut t = SummaryTally::new(1e-9, 100);
    let mut acc = 0.0;
    for i in 0..count {
        let mut rng = stream(SEED, i);
        acc += bca.history(&mut t, &mut rng, i).unwrap().rest;
    }
    acc
}

fn bench(c: &mut Criterion) {
    let ls = LindhardScharff::new();
    let st = Stack::semi_infinite(elemental(14));
    let mut g = c.benchmark_group("crystal_vs_amorphous_B_5keV_Si_200_ions");
    g.throughput(Throughput::Elements(IONS));
    g.sample_size(20);
    let amorphous = Bca::new(beam(7.0, 22.0), &st, config(SEED), &ls, table()).unwrap();
    g.bench_function("amorphous_7_22", |b| {
        b.iter(|| histories(&amorphous, black_box(IONS)))
    });
    let random = Bca::new(beam(7.0, 22.0), &st, config(SEED), &ls, table())
        .unwrap()
        .with_crystal(crystal(7.0, 22.0), &[0])
        .unwrap();
    g.bench_function("crystal_7_22", |b| {
        b.iter(|| histories(&random, black_box(IONS)))
    });
    let channel = Bca::new(beam(45.0, 0.0), &st, config(SEED), &ls, table())
        .unwrap()
        .with_crystal(crystal(45.0, 0.0), &[0])
        .unwrap();
    g.bench_function("crystal_110_channel", |b| {
        b.iter(|| histories(&channel, black_box(IONS)))
    });
    g.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
