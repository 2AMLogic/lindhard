//! Sampled thermal displacements must be bit-identical at any thread count.

use lindhard::ion::crystal::debye::{ThermalVibration, THETA_D_GE};
use lindhard::rng::run_particles;

#[derive(Default, PartialEq, Debug)]
struct Tally {
    sum: [f64; 3],
    sum_sq: f64,
    first: Vec<[f64; 3]>,
}

fn run(threads: usize) -> Tally {
    let v = ThermalVibration::new(THETA_D_GE, 72.63, 300.0).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| {
        run_particles(
            0xD3B7E,
            20_000,
            64,
            Tally::default,
            |t, rng, i| {
                // Several atoms per history, as a lattice search would draw.
                for _ in 0..5 {
                    let d = v.sample_displacement(rng);
                    for (s, x) in t.sum.iter_mut().zip(d) {
                        *s += x;
                    }
                    t.sum_sq += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                    if i < 3 {
                        t.first.push(d);
                    }
                }
            },
            |tot, p| {
                for (s, x) in tot.sum.iter_mut().zip(p.sum) {
                    *s += x;
                }
                tot.sum_sq += p.sum_sq;
                tot.first.extend(p.first);
            },
        )
    })
}

#[test]
fn displacements_are_deterministic_across_thread_counts() {
    let one = run(1);
    assert_eq!(one, run(2));
    assert_eq!(one, run(7));
    assert!(one.sum_sq > 0.0);
}
