# Architecture

This is the design intent, not a contract. Each section names the milestone
that delivers it ([`../WORK_PLAN.md`](../WORK_PLAN.md)).

## Shape

One library crate, `lindhard`, holds all the physics. Front ends sit on top of
it:

- `lindhard-cli`: TOML input, JSON/CSV output (M0); an `[electron]` input runs the electron engine (M1).
- `lindhard-py`: pyo3/maturin Python package, NumPy arrays out (bindings built; PyPI wheels in M3). It reuses the CLI's driver, tallies and writers (`lindhard-cli` is also a small library), so a Python run and a command run of the same input and seed agree bit for bit; it adds no physics and no second copy of the TOML schema.
- `lindhard-wasm`: wasm-bindgen build of the ion engine for in-browser range calculators (M3), with a static page published with the book. It runs the library's transport on one thread and does not use `lindhard-cli` (whose driver needs `std::time::Instant`); a test compares it to the native build for the same seed. Besides the one-shot calculator it has a steppable `Simulation` session for the in-browser simulator: it resolves a full TOML run description, advances it with `Bca::run_range` one engine chunk at a time (so the report and the recorded trajectories do not depend on the step sizes), and returns the `IonReport`, the trajectories and the element table as JSON.

Within the library, the shared core (materials, geometry, random numbers,
tallies, input/output) is separate from the particle engines (ion BCA,
electron MC). The engines share the core and nothing else.

## Modules

| Module | Responsibility | Milestone |
|---|---|---|
| `units`, `constants` | SI internally; CODATA constants with citations | M0 |
| `material` | Elements (Z, mass, density), compounds, mixtures; per-element displacement, lattice and surface binding energies, user-set with documented defaults | M0 |
| `rng` | Counter-based streams keyed on (seed, particle index) so results do not depend on the thread count | M0 |
| `geometry` | `Geometry` trait (regions, per-material data, exit event along a flight) with four implementations: the 1D layered `Stack`, the 3D `VoxelGrid` (exact DDA traversal), `MeshGeometry` (closed STL/OBJ solids, BVH; `geometry::mesh`) and `CsgGeometry` (box, cylinder and half-space primitives under union, intersection and difference; `geometry::csg`). Ownership and tolerance rules: module docs | M0+ |
| `ion::potential` | Screening functions: ZBL universal, Kr-C, Molière, Lenz-Jensen; screening lengths | M0 |
| `ion::scattering` | Scattering integral solved by Gauss–Mehler quadrature; precomputed (ε, b) tables; magic formula kept as a cross-check | M0 |
| `ion::stopping` | Electronic stopping: Lindhard-Scharff, Oen-Robinson, Bethe-Bloch with corrections, user tables with provenance; Bragg additivity plus optional compound corrections. Validity ranges: [`stopping-models.md`](stopping-models.md) | M0 |
| `ion::bca` | Amorphous BCA: free-flight path selection, full recoil cascades, cutoffs, sputtering and backscatter | M0 |
| `ion::damage` | NRT/Kinchin-Pease damage energy (Lindhard partition), reported beside the full-cascade defect counts of `tally::ion`, in different types and never conflated | M0 |
| `ion::crystal` | Crystal data model: `Lattice` (cubic diamond and zincblende, hexagonal wurtzite and polytypes), `Orientation` (tilt, twist, wafer rotation), `Divergence`, Debye thermal displacements and the unit-cell neighbour search. The BCA flight model that uses them is `ion::bca::crystal`. Conventions: [`crystal-orientation.md`](crystal-orientation.md) and the module docs | M2 |
| `ion::dynamic` | Target composition updated with fluence: `CompositionGrid` (per-slab areal inventories and volume relaxation, converted to a `geometry::Stack`), `InventoryTally` (atom counts per slab and element from transport events) and `DynamicRun` (deterministic fluence stepping with optional sputter erosion and a per-species retention limit). Conventions: module docs | M3 |
| `electron::data` | Validated inputs of the electron engine: optical ELF tables, subshell binding energies (from a user-supplied ENDF-6 File 28 copy) and versioned cross-section caches. Every loader requires a provenance | M1 |
| `electron::elastic` | Mott cross sections from our own partial-wave solution (not ELSEPA tables); `elastic::table` builds per-material `CrossSectionTable`s on a log energy grid from an injected atomic potential (Salvat et al. DHFS or a Thomas-Fermi stand-in). Condensed-phase potentials are deferred: [`muffin-tin-deferral.md`](muffin-tin-deferral.md) | M1 |
| `electron::inelastic` | Dielectric-function model (Lindhard / Mermin, Penn algorithm) built from optical data with provenance; inner-shell ionization channels from shell-resolved optical ELFs and an optional Born-Ochkur exchange correction | M1 |
| `electron::secondary`, `electron::boundary` | Secondary-electron generation at inelastic events (Kieft-Bosch energy, Ivanchenko direction) and the inner-potential step at surfaces and interfaces (transmission, refraction, reflection) on per-layer band parameters with provenance; ported from Nebula and cstool (BSD-3) | M1 |
| `electron::transport` | Event-by-event MC from about 10 eV to 50 keV in layered targets on `CrossSectionTable` inputs: exact layer crossings, secondaries, surface barriers, cutoffs, and opt-in per-layer insulator (`electron::phonon`) and inner-shell channels. Models and options: module docs | M1 |
| `tally` | Histograms, moments and Pearson fits; `tally::ion` (range, damage, sputter and backscatter; `IonReport`) and `tally::electron` (deposition grids, yields, escape spectra, energy balance; `ElectronReport`); `tally::table_coverage` (table-range diagnostic); `tally::psf` (pencil-beam PSF profile with Gaussian fits, JSON and CSV export); `tally::voxel` (sparse 3D `VoxelTally` of implanted atoms, displacements and deposited energy on a regular grid, with CSV and VTK export); `tally::trajectory` (opt-in `TrajectoryRecorder`: per-history collision points of the ion and its recoil cascade with parent links, capped in points per history and number of histories, truncation flagged; a plain `BcaTally`, so it costs nothing when unused and `Recorded<T>` pairs it with an ordinary tally without changing that tally's results; convention: module docs) | M0/M1 |

## Performance plan

- Ions and electrons are independent histories, so the work is embarrassingly
  parallel. Use `rayon` over chunks of particles, give each chunk its own
  tally, and merge them deterministically in chunk order, one bounded batch
  of chunks at a time so partial tallies do not grow with the history count.
- Hot paths use lookup tables: the scattering angle from (reduced energy,
  impact parameter) and stopping from (ion, element, energy), both built once
  per run and shared read-only.
- Structure-of-arrays batching plus `wide` SIMD once profiling says it is
  worth doing. Do not do it first.
- A GPU backend is a later, measured decision, not an assumption. The study
  and its current answer (defer) are in
  [`gpu-feasibility.md`](gpu-feasibility.md).
- Benchmarks (criterion) track ions per second against RustBCA and OpenTRIM on
  matched problems (see [`validation.md`](validation.md)).

## Input and output

- Input is TOML. One file describes the beam (species, energy, dose, angles),
  the target (layers, materials) and the physics choices. Every physics choice
  is spelled out in the output's metadata so a result can be reproduced from
  its own header.
- Output is JSON for summaries and CSV for profiles. HDF5 is optional behind a
  feature flag; it pulls in a C library, so it is never the default.
- The input schema, output layout and the rules for extending them are in
  [`cli.md`](cli.md).

## Reproducibility

- **Same platform, any thread count: bit-identical.** Particle `i` draws from
  its own stream, and tallies merge in a fixed order, so results do not depend
  on the number of threads (`lindhard/tests/determinism.rs`).
- **Across platforms: floating-point outputs agree to a relative tolerance,
  not bit for bit.** `f64::exp`, `ln`, `sin`, `cos`, `powf`, `tan` and `atan2`
  defer to the platform libm, which can differ in the last bits between glibc,
  macOS libm and aarch64 builds. Such a difference can in principle flip a
  comparison inside a history and so change an integer count; the engine does
  not promise otherwise. "Same platform" means the same binary, C library
  and libm code path: on x86-64, glibc picks its `exp`, `log` and `pow` code
  by CPU feature (FMA, AVX2), and #319 measured that the FMA and non-FMA
  paths build different cached electron tables from one binary. How a
  validation row is compared across hosts is in
  [`validation.md`](validation.md), "Reproducing the baseline across hosts".
- **What is checked.** `lindhard/tests/golden.rs` pins small fixed-seed runs
  (amorphous BCA batch, crystal BCA batch, electron batch) to values
  generated on x86-64 Linux. The amorphous and electron runs are checked on
  every platform of the CI OS matrix (x86-64 Linux, aarch64 Linux, aarch64
  macOS) with integer outputs exact and floats to a relative tolerance of
  1e-9. The crystal run has a pinned set per platform: one for Linux (both
  architectures) and a separately labelled one for aarch64 macOS (generated on
  the `macos-latest` runner), each checked exactly in the same sense. Any
  other platform is only checked statistically: conservation invariants exact,
  counts and sums to a 15% relative tolerance. The two sets differ (67721
  recoils on aarch64 macOS against 67786 on Linux) because the platform libm
  differs in the last bit for `tan`, `acos`, `sin`, `cos`, `exp`, `ln`,
  `powf` and others. The differences enter through the BCA path shared with
  the amorphous run, where they move floats in the last bits and leave the
  integer counts alone at this sample size; the crystal path's lattice-site
  and simultaneous-collision comparisons amplify them into different
  histories. No single comparison flips. The promise stays tolerance-only:
  a pure-Rust libm that would make the sets identical was considered and
  declined (#217). If the macOS runner image's libm changes, the macOS set
  may need regenerating. The test file says how to regenerate the values and
  which platform each set came from.

## Decisions

- **Muffin-tin potential deferred (#92, 2026-10-07).** `electron::elastic`
  models free-atom potentials only. The muffin-tin definition (truncation
  radius and offset) comes from a closed-access paper that nobody on the
  project has opened, and no open comparison giving the size of the effect was
  found, so the option is not built from memory. The note
  [`muffin-tin-deferral.md`](muffin-tin-deferral.md) states the gap, gives no
  number, and lists what lifts the deferral. A `MuffinTin` wrapper over
  `ScreenedPotential` is the intended shape.
