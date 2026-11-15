//! The thread-pool plumbing every front end shares: reject
//! `threads == Some(0)` and build the [`rayon::ThreadPool`], in one place
//! instead of once per entry point.
//!
//! The canonical rejection of `run.threads = 0` lives in the library's input
//! resolvers (`lindhard::input::Input::resolve_in` and the electron
//! `resolve_run`), so every resolved input is already clean. These helpers
//! keep the direct entry points (`sim::simulate`,
//! `dynamic::simulate_dynamic`, `electron::simulate_electron` and the
//! exposure commands) safe when their `threads` argument never went through
//! a resolver, and give every front end the same pool construction. The
//! thread count never changes the results.

use anyhow::{bail, Context, Result};
use rayon::ThreadPool;

/// Build a rayon pool with `threads` workers (`None`: all available cores,
/// rayon's default). Rejects `Some(0)` the same way the input resolvers
/// reject `run.threads = 0`.
pub fn build(threads: Option<usize>) -> Result<ThreadPool> {
    if threads == Some(0) {
        bail!("threads must be at least 1");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.unwrap_or(0))
        .build()
        .context("building the thread pool")
}

/// Run `f` on a pool of `threads` workers; with `None`, run `f` directly on
/// the current pool (all available cores). The thread count never changes the
/// results.
pub fn with_threads<T: Send>(threads: Option<usize>, f: impl FnOnce() -> T + Send) -> Result<T> {
    match threads {
        None => Ok(f()),
        Some(_) => Ok(build(threads)?.install(f)),
    }
}
