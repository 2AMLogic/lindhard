# lindhard

**A clean-room, MIT-licensed Rust engine for Monte Carlo transport of ions and
electrons in matter.** It computes implantation range profiles, damage,
sputtering and backscatter for ions, and energy deposition for low-energy
electrons. SRIM/TRIM covers the same problem class; lindhard uses none of its
code or data.

> **Status: early, pre-release (0.0.1).** Most of milestone M0, the amorphous
> ion core, has landed: materials and the element table, deterministic
> per-ion random streams and a parallel driver, screened potentials and the
> scattering integral, electronic stopping, the 1D layered BCA engine with
> full recoil cascades, range, damage and sputter/backscatter tallies, and the
> CLI. Validation levels 1 (analytic checks, in CI) and 2 (code-to-code
> oracles) are in place, and level 3 has its first experimental dataset (B in
> amorphous Si; P and As remain open under #51). Known deviations are
> reported, not hidden: with the default inputs the computed B ranges run long
> against that measurement, and the Ar → Cu sputter yield sits about 2x below
> RustBCA's. The investigation in #61 found that the default engine sends too
> much energy to electronic loss at low energy, because nuclear loss is cut off
> at p_max while electronic loss is not. #64 added TRIDYN-style weak collisions
> as an opt-in (`physics.weak_collisions`), which brings the cascade's
> electronic share into line with the LNST partition but lowers the sputter
> yield further, so the gap to RustBCA remains open under #61. See
> [`docs/validation.md`](docs/validation.md),
> [`WORK_PLAN.md`](WORK_PLAN.md) for milestones and
> [`docs/architecture.md`](docs/architecture.md) for the design. The electron
> engine (M1) and crystalline targets (M2) are not started.

## Quick start

```sh
cargo run -p lindhard-cli -- check examples/b_5keV_si.toml
cargo run --release -p lindhard-cli -- run examples/b_5keV_si.toml --out out/b_5keV_si
```

`run` writes `summary.json` and CSV profiles into the output directory. The
input schema and output layout are in [`docs/cli.md`](docs/cli.md); more
inputs are in [`examples/`](examples/README.md).

## Why another code

Some open codes already handle parts of this problem well. The gaps are in
license and scope:

- **SRIM/TRIM** is closed source, Windows-only and non-commercial. It treats
  every target as amorphous and as fixed, so the target never changes with
  dose.
- **Crystalline implant simulation** is not available in any open, maintained
  code. That covers channeling under tilt and twist, screen oxides, and damage
  that builds up to amorphization. Today it means IMSIL, Crystal-TRIM, MARLOWE
  or commercial TCAD.
- **Ion codes that are open** are either copyleft (RustBCA, iradina, IM3D) or
  use stopping tables taken from SRIM (OpenTRIM).
- **Low-energy electron transport** (eV–50 keV, for e-beam lithography
  proximity functions, SEM and EBIC) has one permissive open code, Nebula, and
  it has stalled. The reference physics lives inside Geant4 or AGPL codes.

lindhard aims to be **permissive, embeddable, fast, reproducible and
validated**:

- MIT license, with no copyleft dependency anywhere in the build.
- A library first; the CLI and later Python/WASM bindings sit on top of it.
- Every ion runs on its own seeded random stream, so results do not depend on
  the thread count.
- Every model is checked against analytic limits, independent codes and
  published measurements ([`docs/validation.md`](docs/validation.md)).

[`docs/prior-art.md`](docs/prior-art.md) is the survey this project starts
from. [`docs/history.md`](docs/history.md) traces the method back to the
1947 Los Alamos neutron histories.

The name honours Jens Lindhard. His LSS theory is the standard theory of ion
range, he did the foundational work on channeling, and his dielectric function
is the starting point for inelastic electron scattering. Between them those
cover both halves of this project.

## Layout

| Path | What |
|---|---|
| [`lindhard/`](lindhard/README.md) | The library: materials, potentials, stopping, transport, tallies |
| [`lindhard-cli/`](lindhard-cli/README.md) | `lindhard` binary: TOML in, JSON/CSV out |
| [`lindhard-py/`](lindhard-py/README.md) | `lindhard` Python package (pyo3 + maturin): NumPy arrays out, same TOML schema |
| [`examples/`](examples/README.md) | Example CLI inputs |
| [`validation/`](validation/README.md) | Validation harness: oracle runner, experimental datasets, results |
| [`docs/`](docs/README.md) | Design, physics and project documentation (index below) |

## Documentation

- [The lindhard book](https://2amlogic.github.io/lindhard/) (source in
  [`book/`](book/)): the physics manual, one page per model, and the user
  guide, including a [first run](book/src/guide/first-run.md) from start to
  finish.
- [`docs/architecture.md`](docs/architecture.md): design and module plan.
- [`docs/cli.md`](docs/cli.md): the `lindhard` command, input schema and outputs.
- [`docs/stopping-models.md`](docs/stopping-models.md): electronic stopping
  models, their validity ranges and sources.
- [`docs/validation.md`](docs/validation.md): validation plan and current results.
- [`docs/benchmarks.md`](docs/benchmarks.md): how performance is measured, and
  first numbers.
- [`docs/data-provenance.md`](docs/data-provenance.md): where every dataset
  comes from.
- [`docs/prior-art.md`](docs/prior-art.md) and
  [`docs/history.md`](docs/history.md): the survey and the lineage.

## Contributing

Read [`CONTRIBUTING.md`](CONTRIBUTING.md) first. It sets out the **clean-room
rule**: what may be learned from other codes, what may be ported, and what may
never enter this tree.

## License

MIT. See [`LICENSE`](LICENSE).
