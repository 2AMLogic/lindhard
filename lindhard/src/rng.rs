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
//! `lindhard/tests/determinism.rs`.

use rand_chacha::ChaCha8Rng;
use rand_core::SeedableRng;
use rayon::prelude::*;

/// The generator type returned by [`stream`].
pub type ParticleRng = ChaCha8Rng;

/// Random stream for particle `index` of the run with the given `seed`.
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
        .collect(); // indexed collect preserves chunk order
    let mut total = new_tally();
    for p in partials {
        merge(&mut total, p);
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::Rng;

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
