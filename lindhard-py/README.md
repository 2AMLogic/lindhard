# lindhard (Python)

Python bindings for [lindhard](https://github.com/2AMLogic/lindhard), built
with [maturin](https://www.maturin.rs) and [pyo3](https://pyo3.rs).

```python
import lindhard as lh

run = lh.Run(
    beam=lh.Beam("Ar", 1000.0),
    target=lh.Target(substrate="Cu"),
    physics=lh.Physics(primary_cutoff_ev=2.0, recoil_cutoff_ev=1.0, potential="kr-c"),
    ions=10_000,
    seed=1,
    tally=lh.Tally(depth_bin_nm=0.1, depth_bins=100),
)
result = run.run()                 # releases the GIL while transporting
result.depth.counts                # NumPy array: stopped primaries per depth bin
result.vacancies.edges             # NumPy array: bin edges, nm
result.summary()["results"]["yields"]   # dict, the CLI's summary.json
open("run.toml", "w").write(run.to_toml())   # the CLI's input schema
```

`Run.from_toml()` / `Run.to_toml()` read and write the same TOML as the
`lindhard` command (`docs/cli.md`), and a result equals the command's output
for the same input and seed, bit for bit, at any thread count.

Errors: `lindhard.InputError` (a `ValueError`) for an invalid configuration,
`lindhard.RunError` (a `RuntimeError`) for a failed run; both derive from
`lindhard.LindhardError`.

## Build and test

From `lindhard-py/`, in a virtual environment:

```sh
python -m venv .venv && . .venv/bin/activate
pip install maturin numpy pytest
maturin develop
pytest
```
