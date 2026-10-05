# Architecture

This is the design intent, not a contract. Each section names the milestone
that delivers it ([`../WORK_PLAN.md`](../WORK_PLAN.md)).

## Shape

One library crate, `lindhard`, holds all the physics. Front ends sit on top of
it:

- `lindhard-cli`: TOML input, JSON/CSV output (M0).
- `lindhard-py`: pyo3/maturin wheels on PyPI (M3).
- WASM build for in-browser range calculators (M3).

Within the library, the shared core (materials, geometry, random numbers,
tallies, input/output) is separate from the particle engines (ion BCA,
electron MC). The engines share the core and nothing else.

## Modules

| Module | Responsibility | Milestone |
|---|---|---|
| `units`, `constants` | SI internally; CODATA constants with citations | M0 |
| `material` | Elements (Z, mass, density), compounds, mixtures; per-element displacement, lattice and surface binding energies, user-set with documented defaults | M0 |
| `rng` | Counter-based streams keyed on (seed, particle index) so results do not depend on the thread count | M0 |
| `geometry` | 1D layered target (M0); 2D/3D voxel and triangle mesh (M3); stacks taken from layout cross-sections (M3) | M0+ |
| `ion::potential` | Screening functions: ZBL universal, Kr-C, Molière, Lenz-Jensen; screening lengths | M0 |
| `ion::scattering` | Scattering integral solved by Gauss–Mehler quadrature; precomputed (ε, b) tables; magic formula kept as a cross-check | M0 |
| `ion::stopping` | Electronic stopping: Lindhard-Scharff, Oen-Robinson, Bethe-Bloch with corrections, Biersack-Varelas interpolation, user tables with provenance; Bragg additivity plus optional compound corrections. Validity ranges: [`stopping-models.md`](stopping-models.md) | M0 |
| `ion::bca` | Amorphous BCA: free-flight path selection, full recoil cascades, cutoffs, sputtering and backscatter | M0 |
| `ion::damage` | NRT/Kinchin-Pease damage energy, alongside full-cascade vacancy and interstitial bookkeeping (the two are reported side by side, never conflated) | M0 |
| `ion::crystal` | Lattice-site targets, thermal vibration (Debye), tilt/twist/rotation, screen oxide, dynamic damage → dechanneling → amorphization | M2 |
| `ion::dynamic` | Target composition updated with fluence (sputter erosion, build-up of implanted atoms) | M3 |
| `electron::elastic` | Mott cross sections by our own partial-wave solution (not ELSEPA tables) | M1 |
| `electron::inelastic` | Dielectric-function model (Lindhard / Mermin, Penn algorithm) built from optical data with provenance; SE generation; interface refraction | M1 |
| `electron::transport` | Event-by-event MC from about 10 eV to 50 keV in layered and voxel targets | M1 |
| `tally` | Depth and lateral histograms, moments (Rp, ΔRp, γ, β), Pearson IV / dual-Pearson fits, 3D grids, energy-deposition maps, PSF extraction (double/triple Gaussian α/β/η) | M0/M1 |

## Performance plan

- Ions and electrons are independent histories, so the work is embarrassingly
  parallel. Use `rayon` over chunks of particles, give each worker its own
  tallies, and merge them deterministically at the end.
- Hot paths use lookup tables: the scattering angle from (reduced energy,
  impact parameter) and stopping from (ion, element, energy), both built once
  per run and shared read-only.
- Structure-of-arrays batching plus `wide` SIMD once profiling says it is
  worth doing. Do not do it first.
- A GPU backend is a later, measured decision, not an assumption.
- Benchmarks (criterion) track ions per second against RustBCA and OpenTRIM on
  matched problems (see [`validation.md`](validation.md)).

## Input and output

- Input is TOML. One file describes the beam (species, energy, dose, angles),
  the target (layers, materials) and the physics choices. Every physics choice
  is spelled out in the output's metadata so a result can be reproduced from
  its own header.
- Output is JSON for summaries and CSV for profiles. HDF5 is optional behind a
  feature flag; it pulls in a C library, so it is never the default.
