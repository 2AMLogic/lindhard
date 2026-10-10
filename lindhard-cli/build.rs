//! Stamps `git describe --always --dirty` into the binary as
//! `LINDHARD_GIT_DESCRIBE`, or `"unknown"` outside a git checkout (for example
//! a crates.io build). Never fails the build.
//!
//! A non-empty `LINDHARD_GIT_DESCRIBE` in the build environment overrides the
//! `git` lookup. Release CI sets it so that a wheel built in a container
//! without git (manylinux) and the `lindhard` command built on the host stamp
//! the same string into `summary.json`.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8(out.stdout).ok()?.trim().to_string();
    (!s.is_empty()).then_some(s)
}

fn main() {
    println!("cargo:rerun-if-env-changed=LINDHARD_GIT_DESCRIBE");
    let describe = std::env::var("LINDHARD_GIT_DESCRIBE")
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .or_else(|| git(&["describe", "--always", "--dirty"]))
        .unwrap_or_else(|| "unknown".into());
    println!("cargo:rustc-env=LINDHARD_GIT_DESCRIBE={describe}");
    // Re-run when HEAD moves or the index changes (commits, checkouts).
    if let Some(dir) = git(&["rev-parse", "--absolute-git-dir"]) {
        println!("cargo:rerun-if-changed={dir}/HEAD");
        println!("cargo:rerun-if-changed={dir}/index");
    }
    println!("cargo:rerun-if-changed=build.rs");
}
