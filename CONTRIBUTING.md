# Contributing

## The one hard rule: clean room

lindhard is released under MIT. So the code, comments, data tables and test
fixtures in this tree must be either our own work from published methods, or
ported from a permissively licensed source with its notice kept. We learn from
other codes, but **how** we may learn depends on each code's license.

### Tier A: permissive. Read, learn, and port with attribution

| Code | License | Useful for |
|---|---|---|
| OpenTRIM (ir2-lab) | MIT | Architecture, 3D grid tallies, damage bookkeeping, GUI/HDF5 output ideas |
| DISPLATH (Julia) | MIT | BCA on explicit atom positions (a route to crystalline targets) |
| Nebula + cstool (TU Delft) | BSD-3-Clause | Low-energy electron MC, dielectric-function cross sections, GPU layout |
| ESPNN | MIT | A stopping model trained on experimental data |
| PyEPICS | BSD-3-Clause | EEDL/EADL/EPDL data access |

Porting from a Tier A source is allowed. Keep the upstream copyright notice,
add the source to `THIRD_PARTY_LICENSES.md`, and put a comment at the port
site naming the upstream file and commit. **Exception:** OpenTRIM's electronic
stopping tables come from SRIM-2013, which is Tier C data. Do not port them.

### Tier B: copyleft (GPL, LGPL, AGPL). Papers, manuals and outputs only

RustBCA, iradina, IM3D, libdEdx, CATIMA, JIBAL, Corteo, PenRed, EGSnrc.

Code derived from these could never ship under MIT. You may use:

- **The papers and manuals.** Most of these codes are well documented in
  journals (RustBCA in JOSS, iradina in NIM B, CATIMA and libdEdx in their
  docs). Implement from the paper.
- **Their input/output formats**, as a convenience. A file format is not code.
- **The programs themselves as oracles.** Run them unmodified on your own
  machine and compare our results against theirs (see "Oracles" below).

You may not have their source open while writing code here, transcribe it from
memory, or copy their tables, constants files or test fixtures.

### Tier C: closed or non-commercial. Published literature only

SRIM/TRIM (and every table SRIM produces), SDTrimSP, TRIDYN/TRI3DYN, IMSIL,
Crystal-TRIM, MARLOWE, IIS, CASINO, PHITS, FLUKA, MCNP, ELSEPA, DPASS/PASS,
CasP.

Use only their published papers and manuals. Their output must not be
committed as reference data. In particular, **no SRIM stopping table, and
nothing interpolated or fitted from one, may enter this tree**, and that
includes tables passed along through other projects.

### Oracles

Comparing against third-party programs is encouraged, and the harness under
`validation/` drives them. Their raw outputs are regenerated locally into
`validation/oracle-runs/`, which is gitignored. Commit only our own **summary
metrics** of a comparison (for example "Rp within 3% of oracle X at version
Y"), never their tables or curves.

### Data

Every number that enters the tree as data needs a provenance entry in
[`docs/data-provenance.md`](docs/data-provenance.md). It must be either
computed from a published formula, or an experimental measurement cited to its
paper. "It came from a well-known table" is not provenance.

### CI enforces part of this

The `clean-room` job greps every file for copyleft license notices and SRIM
output headers, and fails the build if it finds one. It catches copied files,
not paraphrase. Paraphrase is on you.

## Disclosure

This repository is **public**. Everything in its history, its issues and its
PRs is published the moment it is pushed or posted. This project implements
*published* physics. If you think you have devised a method that is genuinely
new (one not found in the literature), stop and raise it with the operator
**before** committing it, filing an issue about it, or describing it in a PR.
Do not describe the operator's downstream applications in this repo either.

## Everything else

- Rust stable, `cargo fmt`, `cargo clippy -- -D warnings`, `cargo test` green, and
  `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` clean.
- `#![forbid(unsafe_code)]` in the library stays. SIMD goes through safe
  crates (`wide`), not intrinsics.
- No C/Fortran dependencies in the default build.
- Every dependency must be permissively licensed and come from crates.io.
  `cargo deny check licenses sources` enforces both in CI, and `deny.toml`
  holds the allowlist.
- **Reproducibility is a feature.** Every random draw comes from a counter-based
  stream keyed on (seed, particle index), so a run gives the same bits on 1
  thread or 64. A PR that breaks this needs a test showing why it must. The invariant is
  tested by [`lindhard/tests/determinism.rs`](lindhard/tests/determinism.rs)
  (bit-identical tallies on 1, 2 and 8 threads); the helpers are in
  `lindhard::rng`.
- Physics constants and model choices cite their source in a doc comment.
- User-visible changes (behaviour, input/output formats, CLI flags, public API,
  packaging) add a line under `## [Unreleased]` in [`CHANGELOG.md`](CHANGELOG.md)
  ([Keep a Changelog](https://keepachangelog.com/) style). Internal refactors,
  tests and docs-only changes need no entry.
- The minimum supported Rust version is `rust-version` in the root
  `Cargo.toml`, and the CI `msrv` job checks the workspace on exactly that
  toolchain; raise it deliberately, in its own `CHANGELOG.md` entry, never as
  a side effect of using a newer language or std feature.
- One PR per issue. PRs are reviewed by Loom's Judge and merged by Champion,
  using merge commits.
