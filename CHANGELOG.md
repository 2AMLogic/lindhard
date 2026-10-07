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

- `geometry::Geometry`, the engine-facing target trait, implemented by `Stack`
  and the new `geometry::VoxelGrid` (regular 3D grid of material indices,
  periodic or vacuum boundaries per axis, exact 3D DDA traversal that
  truncates the free path at a material change; `Geometry::flight` also
  reports the region a flight without an event ends in, so the particle's
  region follows it across same-material voxel faces). `Bca::new` now takes any
  `&dyn Geometry` (a `&Stack` still works) and `Bca::with_entry_point` picks
  the incident point. `Face::Side` and `EnergyBudget::lateral` account for
  escapes through the lateral faces of a voxel grid; `Face` moved to
  `geometry` (still re-exported from `ion::bca`) and gained that variant;
  `kinematics::refract_out_normal` applies the surface barrier along an
  arbitrary face normal. `Particle::layer` is now the region index (layer or
  flat voxel index). Stack runs are unchanged.
- `Bca::history_in` and `HistoryBuffers`: `Bca::history` with caller-owned
  working memory, which allocates nothing in steady state.
- `ScatteringTable::half_angle_tan`, `kinematics::rotate_sc` and
  `kinematics::lab_projectile_sc`.
- `ElectronicStopping::sqrt_energy_coefficient` (defaulted, so existing
  implementations are unaffected), implemented by `LindhardScharff`.
- Packaging metadata for both crates: `include` allowlists, crate-local
  README, LICENSE and CHANGELOG in each archive, and a CI `package` job that
  runs `cargo package` and `cargo publish --dry-run`.
