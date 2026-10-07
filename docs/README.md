# Documentation

The physics manual (one page per model: equations, assumptions, validity,
references) and the user guide (install, a first run, the input and output
reference) are an mdBook under [`../book/`](../book/), published at
<https://2amlogic.github.io/lindhard/>. Build it locally with
`mdbook build book`. The book includes `cli.md` and `stopping-models.md`
from this directory rather than copying them, so edit them here.

| Document | What |
|---|---|
| [`architecture.md`](architecture.md) | Design intent: crate shape, module plan, and which milestone delivers each part |
| [`cli.md`](cli.md) | The `lindhard` command: subcommands, TOML input schema, output files, reproducibility |
| [`stopping-models.md`](stopping-models.md) | Electronic stopping models, their validity ranges, sources, and the terms declined or deferred |
| [`stopping-data-format.md`](stopping-data-format.md) | The CSV format for measured stopping points (per-point citation and uncertainty), its validation rules, and the model comparison |
| [`validation.md`](validation.md) | The three validation levels (analytic, oracle, experiment), current results and known deviations |
| [`benchmarks.md`](benchmarks.md) | How to run the benches, what they measure, and first numbers |
| [`data-provenance.md`](data-provenance.md) | One row per dataset in the tree: origin, citation and terms |
| [`prior-art.md`](prior-art.md) | Survey of open and closed ion and electron transport codes |
| [`history.md`](history.md) | A historical note, from the 1947 neutron histories to lindhard |

The project rules are in [`../CONTRIBUTING.md`](../CONTRIBUTING.md) and the
roadmap in [`../WORK_PLAN.md`](../WORK_PLAN.md).
