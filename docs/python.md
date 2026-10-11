# Python bindings

The `lindhard` Python package (crate `lindhard-py`, built with
[maturin](https://www.maturin.rs) and [pyo3](https://pyo3.rs)) runs the same
engine as the command line and returns the tallies as NumPy arrays. Usage and
build steps: [`../lindhard-py/README.md`](../lindhard-py/README.md); type stubs
are in `lindhard-py/python/lindhard/`.

## Classes

| Class | What |
|---|---|
| `Element`, `Material` | Inline material: elements with fractions and optional `e_d_ev` / `e_b_ev` / `e_s_ev`; `density_g_cm3` (required for compounds) |
| `Layer`, `Target` | Finite layers front to back (a material name, element symbol or `Material`, and `thickness_nm`) and an optional semi-infinite substrate |
| `Beam` | `ion`, `energy_ev`, `mass_amu`, `tilt_deg`, `azimuth_deg`, and optionally `divergence_model` (`"gaussian"` or `"uniform-cone"`) with `divergence_deg` (the `[beam.divergence]` table) |
| `Physics` | The `[physics]` table: models (as the TOML names), cutoffs, energy overrides |
| `Tally` | The `[tally]` table (CLI defaults) |
| `Run` | Beam, target, physics, tally, `ions`, `seed`, `threads`, named `materials` and `stopping_tables`; `from_toml()`, `from_toml_file()`, `to_toml()`, `validate()`, `run()`, `run_dynamic()` |
| `RunResult` | `depth`, `lateral_y`, `lateral_z`, `radial`, `vacancies`, `interstitials`, `replacements` (each a `Histogram` with `edges`, `counts`, `density`, `underflow`, `overflow`), `escapes`, `ions` (final state of each primary), `summary()` (dict), `write(dir)` |
| `DynamicRunResult` | Result of `Run.run_dynamic()`: `species`, `steps` (per-step arrays), `composition` (per-step slab profiles), `summary()` (dict), `write(dir)` |

## One schema, one engine

The classes are value wrappers over `lindhard::input`'s structs, so the TOML
schema ([`cli.md`](cli.md)) has a single definition: `Run.to_toml()` output is
valid CLI input and `Run.from_toml()` reads CLI input. `Run.run()` calls the
same driver as `lindhard run` (`lindhard_cli::sim::simulate`) and the arrays
and `summary()` are computed with the same expressions as the CSV and
`summary.json` writers. For the same input and seed the result equals the
command's output bit for bit, at any thread count; `tests/test_cli_parity.py`
checks this against the built binary. `Run.run()` releases the GIL while it
transports.

A `[[crystal]]` section (`cli.md`) is part of that schema: `Run.from_toml()` keeps it,
`Run.to_toml()` writes it back, and `Run.run()` attaches the crystal in the
same driver, so `summary()` carries `physics.crystal` as the command's
`summary.json` does. The Python classes have no crystal constructor of their
own yet; a crystal run starts from TOML.

## Fluence-dependent runs

An input with a `[dynamic]` section (`cli.md`) is a fluence-stepped run. Read
it with `Run.from_toml()` and call `Run.run_dynamic(ions=, seed=, threads=)`,
which applies the same per-call overrides and warnings as `run()`, calls the
CLI's driver (`lindhard_cli::dynamic::simulate_dynamic`) with the GIL
released and returns an immutable `DynamicRunResult`. `Run.run()` is
unchanged and still rejects a `[dynamic]` input. `run_dynamic()` raises
`InputError` for an input without `[dynamic]`.

- `res.species` / `res.species_z`: the elements the run can contain (beam
  species and target elements, sorted by atomic number).
- `res.steps`: a dict of arrays, one entry per step with index 0 the initial
  target: `step`, `first_index`, `ions`, `ions_done`, `fluence_cm2`
  (ions/cm^2), `attempts`, `max_change`, `clamped`, `removed_slabs`,
  `n_slabs`, `surface_nm` (cumulative recession), `thickness_nm`, the
  cumulative `cum_*` event counts, `sputter_yield` and
  `cum_sputtered_by_species` (`(n_steps, len(species))`). These are the
  columns of `dynamic_steps.csv`.
- `res.composition`: a list parallel to `steps`; each entry has `front_nm`,
  `back_nm`, `thickness_nm` (`(n_slabs,)`) and `atoms_per_cm2`, `fraction`
  (`(n_slabs, len(species))`), the columns of `dynamic_composition.csv`.
  Erosion can change the slab count between steps, so the snapshots are
  separate arrays rather than one rectangular tensor.
- `res.summary()` is `dynamic_summary.json` as a dict; `res.write(dir)`
  writes it and the two CSVs with the CLI's policy (the previous summary is
  removed before any CSV changes, the new one is published last).

Arrays equal the CSV values bit for bit for the same input and seed, at any
thread count (`tests/test_dynamic.py`). A walk-through is in
[`../examples/python/dynamic_composition.py`](../examples/python/dynamic_composition.py).

## Errors

| Exception | Base classes | Raised for |
|---|---|---|
| `InputError` | `LindhardError`, `ValueError` | malformed TOML, unknown key, invalid value (message names the field), unknown choice name |
| `RunError` | `LindhardError`, `RuntimeError` | a run that was accepted but failed in transport (for example a stopping table queried outside its range) |

Resolution warnings (the CLI prints them on stderr) are issued as
`UserWarning`; `Run.validate()` returns them as a list.
