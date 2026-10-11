//! Stamps `git describe --always --dirty` into the binary as
//! `LINDHARD_GIT_DESCRIBE`, or `"unknown"` outside a git checkout (for example
//! a crates.io build). Never fails the build.
//!
//! A non-empty `LINDHARD_GIT_DESCRIBE` in the build environment overrides the
//! `git` lookup. Release CI sets it so that a wheel built in a container
//! without git (manylinux) and the `lindhard` command built on the host stamp
//! the same string into `summary.json`.
//!
//! Invalidation. Cargo reruns this script only when a watched path changes,
//! so the watch list has to cover everything `git describe --always --dirty`
//! reads: the workspace's tracked compile inputs (an unstaged edit changes the
//! `-dirty` suffix without touching `.git/HEAD` or the index), the index, the
//! per-worktree `HEAD`, the ref `HEAD` names, `packed-refs`, and the common
//! `refs/heads` and `refs/tags` (shared by linked worktrees). All locations
//! come from `git rev-parse --git-path`, never from a guessed `.git/` layout.
//! Untracked files do not affect `--dirty` and are not watched.

use std::path::{Path, PathBuf};
use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8(out.stdout).ok()?.trim().to_string();
    (!s.is_empty()).then_some(s)
}

/// Emits `rerun-if-changed` for `path` if it exists. A missing path would make
/// Cargo rerun the script on every build, so absent optional inputs are skipped.
fn watch(path: &Path) {
    if path.exists() {
        println!("cargo:rerun-if-changed={}", path.display());
    }
}

/// A path printed by `git rev-parse --git-path`, which is relative to the
/// working directory (the package directory) unless Git made it absolute.
fn git_path(name: &str) -> Option<PathBuf> {
    git(&["rev-parse", "--git-path", name]).map(PathBuf::from)
}

fn watch_git_metadata() {
    for name in ["HEAD", "index", "packed-refs", "refs/heads", "refs/tags"] {
        if let Some(p) = git_path(name) {
            watch(&p);
        }
    }
    // The branch tip lives in its own ref, which does not touch HEAD.
    if let Some(r) = git(&["symbolic-ref", "-q", "HEAD"]) {
        if let Some(p) = git_path(&r) {
            watch(&p);
        }
    }
}

fn watch_tracked_inputs() {
    let Some(top) = git(&["rev-parse", "--show-toplevel"]) else {
        return;
    };
    let top = PathBuf::from(top);
    // Everything compiled into the binary: the workspace manifests and
    // lockfile, the toolchain pin, and the sources of the CLI and the library
    // it links (Cargo scans directories recursively).
    for rel in [
        "Cargo.toml",
        "Cargo.lock",
        "rust-toolchain.toml",
        "lindhard/Cargo.toml",
        "lindhard/build.rs",
        "lindhard/src",
        "lindhard-cli/Cargo.toml",
        "lindhard-cli/src",
    ] {
        watch(&top.join(rel));
    }
}

fn main() {
    println!("cargo:rerun-if-env-changed=LINDHARD_GIT_DESCRIBE");
    println!("cargo:rerun-if-changed=build.rs");
    let override_ = std::env::var("LINDHARD_GIT_DESCRIBE")
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());
    let describe = match override_ {
        // The override is exact and Git is not consulted: nothing else to watch.
        Some(s) => s,
        None => {
            watch_git_metadata();
            watch_tracked_inputs();
            git(&["describe", "--always", "--dirty"]).unwrap_or_else(|| "unknown".into())
        }
    };
    println!("cargo:rustc-env=LINDHARD_GIT_DESCRIBE={describe}");
}
