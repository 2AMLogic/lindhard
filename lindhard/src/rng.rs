//! Deterministic per-particle random streams and a parallel driver.
//!
//! Every particle history draws from its own stream, keyed on
//! `(run seed, particle index)`. The generator is ChaCha8 (the `rand_chacha`
//! crate, MIT OR Apache-2.0): a counter-based design in which the stream
//! identifier is part of the cipher state, so streams for different indices
//! are distinct keystreams by construction, not segments of one sequence
//! spaced apart. Counter-based generators for parallel Monte Carlo are
//! described in J. K. Salmon, M. A. Moraes, R. O. Dror, D. E. Shaw, "Parallel
//! random numbers: as easy as 1, 2, 3", Proc. SC'11 (2011); ChaCha is
//! D. J. Bernstein, "ChaCha, a variant of Salsa20" (2008), here with 8 rounds.
//!
//! The parallel driver [`run_particles`] fixes the floating-point summation
//! order independent of the thread count: particles are split into chunks of a
//! caller-chosen, thread-independent size, each chunk gets its own tally, and
//! the chunk tallies are merged in chunk-index order. The determinism test is
//! `lindhard/tests/determinism.rs`. Chunks are evaluated in consecutive
//! batches of a few chunks per worker and each batch is merged before the next
//! starts, so the number of partial tallies held at once is bounded by the
//! pool size, not by the history count; batching never changes a chunk, a
//! stream or the merge order.

use rand_chacha::ChaCha8Rng;
use rand_core::{Rng, SeedableRng};
use rayon::prelude::*;

/// The generator type returned by [`stream`].
pub type ParticleRng = ChaCha8Rng;

/// Uniform draw on `[0, 1)` with 53 random bits: the top 53 bits of one
/// `next_u64()`, scaled by the exact power of two `2^-53`.
///
/// Every sampler in the crate maps raw words to uniform variates through this
/// one function, so the mapping from a stream to physics is fixed in one place
/// (the thread-count invariant of [`run_particles`] rests on it).
pub fn uniform<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Random stream for particle `index' of the run with the given `seed`.
///
/// The 64-bit seed forms the ChaCha key (little-endian, zero padded) and the
/// particle index is the 64-bit ChaCha stream identifier. The result depends
/// only on `(seed, index)`.
pub fn stream(seed: u64, index: u64) -> ParticleRng {
    let mut key = [0u8; 32];
    key[..8].copy_from_slice(&seed.to_le_bytes());
    let mut rng = ChaCha8Rng::from_seed(key);
    rng.set_stream(index);
    rng
}

/// Run `n_particles` histories in parallel and return the merged tally.
///
/// * `chunk_size` sets the unit of work and of summation order. It must not be
///   derived from the thread count if bit-identical results across thread
///   counts are wanted. A value of 0 is treated as 1.
/// * `new_tally` makes an empty tally; one is created per chunk.
/// * `history(&mut tally, &mut rng, index)` simulates particle `index` using
///   the stream `stream(seed, index)`.
/// * `merge(&mut total, chunk_tally)` folds a chunk tally into the total. It is
///   called sequentially in increasing chunk order, starting from an empty
///   tally, so the result is independent of scheduling and thread count.
///
/// Chunks run in parallel in consecutive batches of a fixed small multiple of
/// the worker count; each batch is merged, in chunk order, before the next
/// begins. At most one batch of chunk tallies (plus the total) is alive at a
/// time, independent of `n_particles`. The batch length affects scheduling
/// only, never the result.
///
/// Runs on the current rayon pool (use `ThreadPool::install` to choose it).
/// Particle indices are `0..n_particles`; [`run_particles_range`] continues a
/// global stream from another first index.
pub fn run_particles<T, N, H, M>(
    seed: u64,
    n_particles: u64,
    chunk_size: u64,
    new_tally: N,
    history: H,
    merge: M,
) -> T
where
    T: Send,
    N: Fn() -> T + Sync,
    H: Fn(&mut T, &mut ParticleRng, u64) + Sync,
    M: FnMut(&mut T, T),
{
    run_particles_range(seed, 0, n_particles, chunk_size, new_tally, history, merge)
}

/// [`run_particles`] for the global particle indices
/// `first_index..first_index + n_particles`: particle `first_index + k` uses
/// the stream `stream(seed, first_index + k)` and is passed that index.
///
/// This is how a run split into several parts (the fluence steps of a dynamic
/// run) continues one global stream instead of replaying indices from zero.
/// Chunks are counted from `first_index`, so with `first_index = 0` this is
/// exactly [`run_particles`]. Because every history depends only on
/// `(seed, index)`, splitting `0..n` into consecutive ranges draws the same
/// random numbers as one run of `n` particles.
pub fn run_particles_range<T, N, H, M>(
    seed: u64,
    first_index: u64,
    n_particles: u64,
    chunk_size: u64,
    new_tally: N,
    history: H,
    mut merge: M,
) -> T
where
    T: Send,
    N: Fn() -> T + Sync,
    H: Fn(&mut T, &mut ParticleRng, u64) + Sync,
    M: FnMut(&mut T, T),
{
    let batch_len = rayon::current_num_threads().saturating_mul(BATCH_CHUNKS_PER_WORKER);
    run_chunk_batches(
        seed,
        first_index,
        n_particles,
        chunk_size,
        batch_len as u64,
        &new_tally,
        &history,
        &mut merge,
    )
}

/// Chunks per rayon worker in one batch of [`run_particles_range`].
///
/// The driver evaluates at most `BATCH_CHUNKS_PER_WORKER` times the worker
/// count chunk tallies at once, so the number of live partial tallies is
/// bounded by the pool size rather than by the history count. The value only
/// shapes scheduling: chunk boundaries, the stream of every history and the
/// merge order are fixed by `chunk_size`, so results do not depend on it.
/// A few chunks per worker leave rayon room to balance unequal chunk costs
/// inside a batch while keeping the idle tail at each batch boundary short.
const BATCH_CHUNKS_PER_WORKER: usize = 4;

/// Shared body of [`run_particles_range`] with an explicit batch length (in
/// chunks; 0 is treated as 1).
///
/// Consecutive batches of `batch_len` chunk indices are evaluated in parallel;
/// each batch is collected in chunk-index order and merged into the running
/// total before the next batch starts. Every chunk is therefore merged in the
/// same increasing order, into the same initially empty total, as a
/// collect-all-then-fold driver would, and at most `batch_len` chunk tallies
/// (plus the total) are alive at any time.
#[allow(clippy::too_many_arguments)]
fn run_chunk_batches<T, N, H, M>(
    seed: u64,
    first_index: u64,
    n_particles: u64,
    chunk_size: u64,
    batch_len: u64,
    new_tally: &N,
    history: &H,
    merge: &mut M,
) -> T
where
    T: Send,
    N: Fn() -> T + Sync,
    H: Fn(&mut T, &mut ParticleRng, u64) + Sync,
    M: FnMut(&mut T, T),
{
    let chunk = chunk_size.max(1);
    let n_chunks = n_particles.div_ceil(chunk);
    let batch_len = batch_len.max(1);
    let mut total = new_tally();
    let mut batch_start = 0;
    while batch_start < n_chunks {
        let batch_end = batch_start.saturating_add(batch_len).min(n_chunks);
        let partials: Vec<T> = (batch_start..batch_end)
            .into_par_iter()
            .map(|c| {
                let mut tally = new_tally();
                let start = c * chunk;
                let end = (start + chunk).min(n_particles);
                for k in start..end {
                    let i = first_index + k;
                    let mut rng = stream(seed, i);
                    history(&mut tally, &mut rng, i);
                }
                tally
            })
            .collect(); // indexed collect preserves chunk order
        for p in partials {
            merge(&mut total, p);
        }
        batch_start = batch_end;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    #[test]
    fn uniform_is_the_top_53_bits_scaled_by_two_to_the_minus_53() {
        let (mut a, mut b) = (stream(7, 3), stream(7, 3));
        for _ in 0..1000 {
            let u = uniform(&mut a);
            assert_eq!(u, (b.next_u64() >> 11) as f64 / (1u64 << 53) as f64);
            assert!((0.0..1.0).contains(&u));
        }
    }

    /// Test-only copy of the driver before batching (#393): every chunk tally
    /// is collected, then all are folded in chunk order into an empty total.
    fn collect_then_fold<T, N, H, M>(
        seed: u64,
        first_index: u64,
        n_particles: u64,
        chunk_size: u64,
        new_tally: N,
        history: H,
        mut merge: M,
    ) -> T
    where
        T: Send,
        N: Fn() -> T + Sync,
        H: Fn(&mut T, &mut ParticleRng, u64) + Sync,
        M: FnMut(&mut T, T),
    {
        let chunk = chunk_size.max(1);
        let n_chunks = n_particles.div_ceil(chunk);
        let partials: Vec<T> = (0..n_chunks)
            .into_par_iter()
            .map(|c| {
                let mut tally = new_tally();
                let start = c * chunk;
                let end = (start + chunk).min(n_particles);
                for k in start..end {
                    let i = first_index + k;
                    let mut rng = stream(seed, i);
                    history(&mut tally, &mut rng, i);
                }
                tally
            })
            .collect();
        let mut total = new_tally();
        for p in partials {
            merge(&mut total, p);
        }
        total
    }

    /// A tally whose merge is not associative: a floating-point sum whose
    /// rounding depends on the merge tree, a running hash over the merge
    /// sequence, and the merged history indices in order.
    #[derive(Default, Debug, PartialEq)]
    struct Order {
        sum: f64,
        hash: u64,
        indices: Vec<u64>,
        merges: u64,
    }

    fn order_history(t: &mut Order, rng: &mut ParticleRng, i: u64) {
        let u = uniform(rng);
        // Spread magnitudes so the rounding of the sum depends on order.
        t.sum += (u - 0.5) * 10f64.powi((i % 17) as i32 - 8);
        t.indices.push(i);
    }

    fn order_merge(a: &mut Order, b: Order) {
        a.sum += b.sum;
        a.hash =
            a.hash.rotate_left(7).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ b.sum.to_bits() ^ b.hash;
        a.indices.extend(b.indices);
        a.merges += 1 + b.merges;
    }

    #[test]
    fn batched_driver_matches_collect_then_fold() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        pool.install(|| {
            // (first index, histories, chunk size); covers nonzero first
            // indices, incomplete final chunks, zero histories and chunk 0.
            let cases = [
                (0, 1000, 7),
                (12_345, 1000, 7),
                (3, 999, 10),
                (0, 0, 5),
                (40, 0, 0),
                (5, 37, 0),
                (0, 1, 64),
                (u64::MAX / 2, 257, 16),
            ];
            for (first, n, chunk) in cases {
                let want = collect_then_fold(
                    0xBADC0DE,
                    first,
                    n,
                    chunk,
                    Order::default,
                    order_history,
                    order_merge,
                );
                assert_eq!(want.indices, (first..first + n).collect::<Vec<_>>());
                // Batch lengths of one chunk, a few, more than all chunks,
                // the pool default, and 0 (treated as 1).
                for batch in [0, 1, 2, 3, 8, 1_000_000] {
                    let got = run_chunk_batches(
                        0xBADC0DE,
                        first,
                        n,
                        chunk,
                        batch,
                        &Order::default,
                        &order_history,
                        &mut order_merge,
                    );
                    assert_eq!(
                        got.sum.to_bits(),
                        want.sum.to_bits(),
                        "{first} {n} {chunk} {batch}"
                    );
                    assert_eq!(got, want, "{first} {n} {chunk} {batch}");
                }
                let got = run_particles_range(
                    0xBADC0DE,
                    first,
                    n,
                    chunk,
                    Order::default,
                    order_history,
                    order_merge,
                );
                assert_eq!(got, want, "{first} {n} {chunk} default batch");
            }
        });
    }

    /// Counts live instances and records the peak.
    struct Live {
        live: Arc<AtomicUsize>,
        count: u64,
    }

    impl Live {
        fn new(live: &Arc<AtomicUsize>, peak: &Arc<AtomicUsize>) -> Self {
            let now = live.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            Live {
                live: Arc::clone(live),
                count: 0,
            }
        }
    }

    impl Drop for Live {
        fn drop(&mut self) {
            self.live.fetch_sub(1, Ordering::SeqCst);
        }
    }

    /// Peak number of live tallies over one run and the merged history count.
    fn live_peak(threads: usize, n: u64, chunk: u64, batch: Option<u64>) -> (usize, u64) {
        let live = Arc::new(AtomicUsize::new(0));
        let peak = Arc::new(AtomicUsize::new(0));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let new = || Live::new(&live, &peak);
        let history = |t: &mut Live, _: &mut ParticleRng, _: u64| t.count += 1;
        let mut merge = |a: &mut Live, b: Live| a.count += b.count;
        let total = pool.install(|| match batch {
            Some(b) => run_chunk_batches(1, 0, n, chunk, b, &new, &history, &mut merge),
            None => run_particles_range(1, 0, n, chunk, new, history, merge),
        });
        let count = total.count;
        drop(total);
        assert_eq!(live.load(Ordering::SeqCst), 0, "a tally leaked");
        (peak.load(Ordering::SeqCst), count)
    }

    #[test]
    fn live_partial_tallies_stay_within_the_batch_bound() {
        // Explicit batch of 5 chunks: at most 5 partials plus the total, for
        // any number of chunks, while the old driver held every chunk.
        for n_chunks in [1u64, 5, 6, 50, 2_000] {
            let (peak, count) = live_peak(2, n_chunks * 3, 3, Some(5));
            assert_eq!(count, n_chunks * 3);
            assert!(peak <= 5 + 1, "{n_chunks} chunks: peak {peak}");
        }
        // Public driver on a 2-worker pool, 2_000 chunks = 250 batches.
        let bound = 2 * BATCH_CHUNKS_PER_WORKER + 1;
        for n_chunks in [10u64, 2_000] {
            let (peak, count) = live_peak(2, n_chunks * 4, 4, None);
            assert_eq!(count, n_chunks * 4);
            assert!(peak <= bound, "{n_chunks} chunks: peak {peak} > {bound}");
        }
        // The reference holds all of them: the bound is a real change.
        let live = Arc::new(AtomicUsize::new(0));
        let peak = Arc::new(AtomicUsize::new(0));
        let total = collect_then_fold(
            1,
            0,
            2_000,
            1,
            || Live::new(&live, &peak),
            |t, _, _| t.count += 1,
            |a, b| a.count += b.count,
        );
        assert_eq!(total.count, 2_000);
        assert!(peak.load(Ordering::SeqCst) > 2_000);
    }

    #[test]
    fn ranges_continue_the_global_stream() {
        let collect = |first: u64, n: u64| {
            run_particles_range(
                9,
                first,
                n,
                4,
                Vec::new,
                |v: &mut Vec<(u64, u64)>, rng, i| v.push((i, rng.next_u64())),
                |a, b| a.extend(b),
            )
        };
        let whole = collect(0, 10);
        let mut parts = collect(0, 3);
        parts.extend(collect(3, 7));
        assert_eq!(whole, parts);
        assert_eq!(whole[5], (5, stream(9, 5).next_u64()));
    }

    #[test]
    fn stream_is_reproducible_and_distinct() {
        let x = stream(1, 5).next_u64();
        assert_eq!(x, stream(1, 5).next_u64());
        assert_ne!(x, stream(1, 6).next_u64());
        assert_ne!(x, stream(2, 5).next_u64());
    }
}
