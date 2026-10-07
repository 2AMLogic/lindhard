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
| `Beam` | `ion`, `energy_ev`, `mass_amu`, `tilt_deg`, `azimuth_deg` |
| `Physics` | The `[physics]` table: models (as the TOML names), cutoffs, energy overrides |
| `Tally` | The `[tally]` table (CLI defaults) |
| `Run` | Beam, target, physics, tally, `ions`, `seed`, `threads`, named `materials` and `stopping_tables`; `from_toml()`, `from_toml_file()`, `to_toml()`, `validate()`, `run()` |
| `RunResult` | `depth`, `lateral_y`, `lateral_z`, `radial`, `vacancies`, `interstitials`, `replacements` (each a `Histogram` with `edges`, `counts`, `density`, `underflow`, `overflow`), `escapes`, `ions` (final state of each primary), `summary()` (dict), `write(dir)` |

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

## Errors

| Exception | Base classes | Raised for |
|---|---|---|
| `InputError` | `LindhardError`, `ValueError` | malformed TOML, unknown key, invalid value (message names the field), unknown choice name |
| `RunError` | `LindhardError`, `RuntimeError` | a run that was accepted but failed in transport (for example a stopping table queried outside its range) |

Resolution warnings (the CLI prints them on stderr) are issued as
`UserWarning`; `Run.validate()` returns them as a list.
