# Installing

lindhard is pre-release and is not yet published as a package; build it
from source. It is pure Rust with no C or Fortran dependencies in the
default build.

## Prerequisites

- A stable Rust toolchain. The repository's `rust-toolchain.toml` selects
  the stable channel, so [rustup](https://rustup.rs) installs the right one
  on first use.
- git.

## Build

```sh
git clone https://github.com/2AMLogic/lindhard.git
cd lindhard
cargo build --release -p lindhard-cli
```

The binary is `target/release/lindhard`. Build in release mode for real
runs; a debug build is many times slower. You can run it in place, through
cargo, or install it into cargo's binary directory:

```sh
target/release/lindhard --version
cargo run --release -p lindhard-cli -- --version
cargo install --path lindhard-cli      # puts `lindhard` on your PATH
```

`--version` prints the crate version and the `git describe` of the source it
was built from, for example `lindhard 0.0.1 (c2da637)`. The same two values
are written into every `summary.json`, so a result can be traced to the
binary that produced it.

## Commands

```text
lindhard check input.toml                  # parse and validate, no transport
lindhard run input.toml --out dir/         # run, write dir/summary.json and CSVs
lindhard run input.toml --out dir/ --ions 200 --seed 7 --threads 4
lindhard --version
```

`--ions` and `--seed` override `run.ions` and `run.seed`, and the output
records the override. `--threads` overrides `run.threads` and never changes
the results. The [first run](first-run.md) walks through both commands.

## The library

The physics is in the `lindhard` library crate; the command line is a thin
front end on it. To browse the library API:

```sh
cargo doc -p lindhard --open
```
