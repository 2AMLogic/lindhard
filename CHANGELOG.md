# Changelog

All notable changes to `lindhard` and `lindhard-cli` are recorded here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the crates use [Semantic Versioning](https://semver.org/) (both crates share
one version).

## [Unreleased]

### Changed

- Collision hot path (about 2.2 to 2.7 times the ions/s, see
  `docs/benchmarks.md`): the scattering angle is carried as `tan(theta/2)` and
  sines and cosines instead of angles, the Lindhard-Scharff coefficient is
  computed once per ion and target instead of at every flight, and
  `Bca::run` reuses its particle stack and scratch buffers across histories.
  Results are statistically unchanged but not bit-identical to 0.0.1 (the
  floating-point operations differ); they remain bit-identical across thread
  counts.

### Added

- Optical ELF datasets for Al and Cu (`validation/data/optical/`, Hagemann,
  Gudat and Kunz 1975, read through the CC0 refractiveindex.info
  transcription and checked against the scanned tables of DESY report
  SR-74/7), with a sum-rule test using power-law quadrature between knots.
  Cu passes both 5 % checks; Al passes the perfect-screening rule and misses
  the f-sum rule at +7.5 % (the authors' own value is +4.5 %), pinned and
  explained in the test (#98). Si is a documented gap.
- `Bca::history_in` and `HistoryBuffers`: `Bca::history` with caller-owned
  working memory, which allocates nothing in steady state.
- `ScatteringTable::half_angle_tan`, `kinematics::rotate_sc` and
  `kinematics::lab_projectile_sc`.
- `ElectronicStopping::sqrt_energy_coefficient` (defaulted, so existing
  implementations are unaffected), implemented by `LindhardScharff`.
- Packaging metadata for both crates: `include` allowlists, crate-local
  README, LICENSE and CHANGELOG in each archive, and a CI `package` job that
  runs `cargo package` and `cargo publish --dry-run`.
