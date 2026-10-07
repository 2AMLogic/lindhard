# Architecture

This is the design intent, not a contract. Each section names the milestone
that delivers it ([`../WORK_PLAN.md`](../WORK_PLAN.md)).

## Shape

One library crate, `lindhard`, holds all the physics. Front ends sit on top of
it:

- `lindhard-cli`: TOML input, JSON/CSV output (M0); an `[electron]` input runs the electron engine (M1).
- `lindhard-py`: pyo3/maturin Python package, NumPy arrays out (bindings built; PyPI wheels in M3). It reuses the CLI's driver, tallies and writers (`lindhard-cli` is also a small library), so a Python run and a command run of the same input and seed agree bit for bit; it adds no physics and no second copy of the TOML schema.
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
| `geometry` | `Geometry` trait (regions, per-material data, exit-event along a flight); 1D layered `Stack` (M0); 3D `VoxelGrid` with exact DDA traversal, periodic or vacuum per axis; triangle mesh and stacks taken from layout cross-sections (M3) | M0+ |
| `ion::potential` | Screening functions: ZBL universal, Kr-C, Molière, Lenz-Jensen; screening lengths | M0 |
| `ion::scattering` | Scattering integral solved by Gauss–Mehler quadrature; precomputed (ε, b) tables; magic formula kept as a cross-check | M0 |
| `ion::stopping` | Electronic stopping: Lindhard-Scharff, Oen-Robinson, Bethe-Bloch with corrections, user tables with provenance; Bragg additivity plus optional compound corrections. Validity ranges: [`stopping-models.md`](stopping-models.md) | M0 |
| `ion::bca` | Amorphous BCA: free-flight path selection, full recoil cascades, cutoffs, sputtering and backscatter | M0 |
| `ion::damage` | NRT/Kinchin-Pease damage energy (Lindhard partition), alongside full-cascade vacancy, interstitial and replacement counts kept by `tally::ion` (the two are reported side by side, in different types, never conflated) | M0 |
| `ion::crystal` | Lattice-site targets, thermal vibration (Debye), tilt/twist/rotation, screen oxide, dynamic damage → dechanneling → amorphization | M2 |
| `ion::dynamic` | Target composition updated with fluence (sputter erosion, build-up of implanted atoms). `CompositionGrid` holds finite slabs with per-element areal inventories (atoms/m²), changed only by explicit deltas, then relaxes thicknesses (ideal mixing of atomic volumes, or a fixed mixture number density) and converts to a `geometry::Stack`; the front surface stays at x = 0 and the optional substrate is immutable. `InventoryTally` turns the transport events of a block of primaries into integer atom counts per slab and element (a recoil is subtracted where it is created and added where it stops), and `DynamicRun` is the fluence stepping loop: fixed or adaptive steps, each covering a range of the global particle indices so that a run continues one random stream, with deterministic merges at any thread count. Conventions: module docs | M3 |
| `electron::data` | Validated data the electron engine consumes: optical ELF tables, subshell binding energies (read from an ENDF-6 File 28 copy the user supplies), and versioned cross-section caches (inverse mean free path plus inverse CDFs of the elastic angle or inelastic energy loss). Every loader requires a provenance; serde reads run the same checks | M1 |
| `electron::elastic` | Mott cross sections by our own partial-wave solution (not ELSEPA tables); `elastic::table` runs it over a log energy grid for every element of a `Material` and combines elements by independent-atom additivity into a `CrossSectionTable` (inverse mean free path plus the inverse CDF of the polar angle). The atomic potential is injected (`PotentialSource`); the DHFS table is still a gap, so a Thomas-Fermi Yukawa stand-in is what runs today; the muffin-tin option for condensed targets is deferred, see [`muffin-tin-deferral.md`](muffin-tin-deferral.md) | M1 |
| `electron::inelastic` | Dielectric-function model (Lindhard / Mermin, Penn algorithm) built from optical data with provenance; inner-shell ionization channels from shell-resolved optical ELFs and an optional Born-Ochkur exchange correction | M1 |
| `electron::secondary`, `electron::boundary` | Secondary-electron generation at inelastic events (Kieft-Bosch energy and Ivanchenko direction model) and the inner-potential step at surfaces and interfaces (quantum transmission, refraction, reflection), on per-layer band parameters (Fermi energy, work function or electron affinity, band gap) that carry a provenance; ported from Nebula and cstool (BSD-3) | M1 |
| `electron::transport` | Event-by-event MC from about 10 eV to 50 keV in layered and voxel targets. Layered stacks on `CrossSectionTable` inputs are implemented (exact layer-face crossings with path redraw, configurable cutoff and escape rule recorded in run metadata, `ElectronTally` hooks mirroring `ion::bca`), secondaries on a per-history stack, surface and interface barriers, and a cutoff from the band bottom or the vacuum level. Opt-in per-layer insulator channels (`electron::phonon`) are implemented: Fröhlich LO-phonon emission and absorption with the Fröhlich angular distribution, and polaron trapping (`C exp(-γE)`, the electron deposits its energy and the history ends); off by default and for metals, recorded per layer in the run metadata | M1 |
| `tally` | Depth and lateral histograms, moments (Rp, ΔRp, γ, β), Pearson IV / dual-Pearson fits, damage and sputter/backscatter tallies (M0, `tally::ion` on the BCA hooks, plain-data `IonReport`); electron deposition on Cartesian and cylindrical r-z grids, η/δ yields split at a configurable energy (50 eV default), escape spectra, generation-volume moments and the energy balance (M1, `tally::electron` on the transport hooks, plain-data `ElectronReport`); PSF extraction (M1, `tally::psf`): the radial profile of a pencil beam in a depth slab on log radial bins with per-history errors, fed by `tally::electron`, and double/triple-Gaussian fits (α, β, η[, γ, ν]) with covariance, reduced χ² and residuals, exported as JSON and CSV | M0/M1 |

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
- The input schema, output layout and the rules for extending them are in
  [`cli.md`](cli.md).

## Decisions

- **Muffin-tin potential deferred (#92, 2026-10-07).** `electron::elastic`
  models free-atom potentials only. The muffin-tin definition (truncation
  radius and offset) comes from a closed-access paper that nobody on the
  project has opened, and no open comparison giving the size of the effect was
  found, so the option is not built from memory. The note
  [`muffin-tin-deferral.md`](muffin-tin-deferral.md) states the gap, gives no
  number, and lists what lifts the deferral. A `MuffinTin` wrapper over
  `ScreenedPotential` is the intended shape.
