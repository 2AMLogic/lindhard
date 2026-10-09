# lindhard (Python)

Python bindings for [lindhard](https://github.com/2AMLogic/lindhard), built
with [maturin](https://www.maturin.rs) and [pyo3](https://pyo3.rs).

## Install

Once a release has been published to PyPI:

```sh
pip install lindhard
```

Each release ships one `abi3` wheel per platform (`cp39-abi3-*`), which works
on every CPython from 3.9 on, for these platforms: Linux x86-64 and aarch64
(manylinux), macOS x86-64 and arm64, and Windows x86-64. A source
distribution is also published; it needs a Rust toolchain to build.

Not yet covered, and tracked as follow-ups: musllinux (for example Alpine) and
Windows arm64. On those, build from the source distribution or from a checkout
(see below).

## Use

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

## Releases (maintainers)

`.github/workflows/wheels.yml` builds the wheels. A pull request touching the
bindings or the library builds and tests one Linux wheel. Pushing a `v*` tag
builds every wheel and the sdist, installs each wheel into a clean virtual
environment, runs this test suite against it on Python 3.9 and 3.13 (using the
`lindhard` command built on the same runner, via `$LINDHARD_BIN`), and then
publishes to PyPI. The tag must equal the workspace version in `Cargo.toml`
(tag `v0.0.1` for version `0.0.1`), or the run stops before publishing. Running
the workflow by hand builds and tests but does not publish.

Publishing uses PyPI trusted publishing (OpenID Connect); no API token is
stored. One-time setup by a project owner, which the workflow cannot do:

1. Create the `lindhard` project on PyPI (or add a pending publisher for it).
2. Add a trusted publisher to it: owner `2AMLogic`, repository `lindhard`,
   workflow `wheels.yml`, environment `pypi`.
3. In the GitHub repository settings, create an environment named `pypi` and
   protect it with required reviewers, so that each upload needs an approval.
