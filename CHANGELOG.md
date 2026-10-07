# Changelog

All notable changes to `lindhard` and `lindhard-cli` are recorded here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the crates use [Semantic Versioning](https://semver.org/) (both crates share
one version).

## [Unreleased]

### Added

- `[physics] stopping = "none"`: no electronic stopping (nuclear loss only),
  for like-for-like comparisons with codes run with electronic stopping off
  and for nuclear-only studies. Echoed in `summary.json` and listed in
  `physics.models`; `[stopping]` tables still serve the pairs they declare.
  Library: `ion::stopping::none::NoStopping`. Level-2 problems
  `b_5keV_si_nuclear`, `as_50keV_si_nuclear` and `ar_1keV_cu_nuclear`.
- Packaging metadata for both crates: `include` allowlists, crate-local
  README, LICENSE and CHANGELOG in each archive, and a CI `package` job that
  runs `cargo package` and `cargo publish --dry-run`.
