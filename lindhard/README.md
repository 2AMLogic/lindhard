# lindhard (library)

The library crate. All the physics lives here; the
[`lindhard` command](https://github.com/2AMLogic/lindhard/blob/main/lindhard-cli/README.md) and the Python bindings sit on
top of it. Design: [`docs/architecture.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/architecture.md).

## Modules (`src/`)

| Module | What |
|---|---|
| `units`, `constants` | Unit conventions and physical constants (SI) |
| `elements` | Element table, Z = 1 to 92 |
| `material` | Materials by atom or mass fraction, density, and per-element `E_d`, `E_b`, `E_s` |
| `geometry` | Target geometry: the 1D layered stack and the 3D voxel grid |
| `rng` | Deterministic per-particle random streams and the parallel driver |
| `input` | The run description (beam, target, materials, physics, run size), parsed from TOML |
| `ion::potential` | Screened-Coulomb interatomic potentials |
| `ion::scattering` | Classical elastic scattering: quadrature, magic formula, lookup table |
| `ion::stopping` | Electronic stopping from published formulas: Lindhard-Scharff, Oen-Robinson and their mix, Bethe-Bloch, Bragg additivity, straggling, user tables ([`docs/stopping-models.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/stopping-models.md)) |
| `ion::bca` | Amorphous BCA transport in a 1D layered target, with full recoil cascades |
| `ion::damage` | Displacement damage: Lindhard partition, NRT, Kinchin-Pease |
| `ion::dynamic` | Dynamic composition: a finite-slab target whose composition and thickness change with fluence |
| `electron` | Electron transport (M1): interaction data, Mott elastic and dielectric-function inelastic cross sections, secondaries, surface barriers, phonon and polaron losses, and the event loop |
| `tally` | Histograms, mergeable moments, Pearson IV and dual-Pearson fits, the full ion tally (range, damage, escapes) and the electron tally (deposition grids, yields, spectra) |

## Tests and benches

- `tests/`: integration tests, including determinism across thread counts.
- `tests/validation/`: the level-1 validation harness, a `cargo test` target
  ([`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md)).
- `examples/`: small example programs (`cargo run -p lindhard --example <name>`).
- `benches/`: criterion benches, local only
  ([`docs/benchmarks.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/benchmarks.md)).
