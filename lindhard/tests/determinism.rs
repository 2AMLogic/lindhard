//! A toy random-walk tally must be bit-identical at any thread count.

use lindhard::rng::run_particles;
use rand_core::Rng;

#[derive(Default, PartialEq, Debug)]
struct Tally {
    sum: f64,
    sum_sq: f64,
    hist: [u64; 16],
}

fn run(threads: usize) -> Tally {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| {
        run_particles(
            0xC0FFEE,
            10_000,
            64,
            Tally::default,
            |t, rng, _i| {
                // 100-step walk with irrational-ish step sizes so that
                // summation order matters in floating point.
                let mut x = 0.0f64;
                for _ in 0..100 {
                    let u = (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
                    x += (u - 0.5) * 0.1234567;
                }
                t.sum += x;
                t.sum_sq += x * x;
                t.hist[((x * 2.0 + 8.0).clamp(0.0, 15.0)) as usize] += 1;
            },
            |a, b| {
                a.sum += b.sum;
                a.sum_sq += b.sum_sq;
                for (h, g) in a.hist.iter_mut().zip(b.hist) {
                    *h += g;
                }
            },
        )
    })
}

#[test]
fn tally_is_bit_identical_across_thread_counts() {
    let r1 = run(1);
    for n in [2, 8] {
        let r = run(n);
        assert_eq!(
            r1.sum.to_bits(),
            r.sum.to_bits(),
            "sum differs on {n} threads"
        );
        assert_eq!(r1.sum_sq.to_bits(), r.sum_sq.to_bits());
        assert_eq!(r1.hist, r.hist);
    }
    assert_ne!(r1.sum, 0.0);
}
