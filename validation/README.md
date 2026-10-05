# Validation harness

The three levels of [`docs/validation.md`](../docs/validation.md). One command
regenerates the results tables there:

```sh
validation/run.sh [--experiments] [--oracles]
```

| Path | Level | What |
|---|---|---|
| `../lindhard/tests/validation/` | 1 | Analytic and internal checks; a `cargo test` target, run in CI |
| `oracles/run.py`, `oracles/problems.json` | 2 | Runs `lindhard` and user-installed oracles (RustBCA, OpenTRIM) on matched problems |
| `oracles/summaries/` | 2 | Committed comparison summaries (our metrics only) |
| `oracle-runs/` | 2, 3 | Raw runs, regenerated locally; **gitignored** |
| `data/` | 3 | Published measurements with citations ([schema](data/README.md)) |
| `experiments/run.py` | 3 | Runs `lindhard` on each dataset; writes `experiments/results.json` |
| `lib/lindhard_cli.py` | 2, 3 | Builds and drives the `lindhard` command |
| `update_docs.py` | all | Splices the tables into `docs/validation.md` between marker comments |

The Python scripts use the standard library only (Python 3.9 or newer).

## Clean-room rules that apply here

From [`CONTRIBUTING.md`](../CONTRIBUTING.md):

- Oracles are run **unmodified**, installed by the user, and never vendored.
  Adapters write an oracle's input from its *published documentation* (paper,
  manual), never from its source code or its example inputs and test
  fixtures. RustBCA is GPL (Tier B); OpenTRIM is MIT (Tier A), but its
  SRIM-2013 stopping tables are Tier C and are neither used for matched runs
  nor copied.
- Raw oracle output stays in `oracle-runs/`. Commit only our summary metrics,
  with the oracle's name and version.
- No SRIM/TRIM output anywhere, including as reference data. SRIM may be run
  locally as a familiarity check; its output never becomes a fixture.
- Every experimental dataset is cited to its paper and has a row in
  [`docs/data-provenance.md`](../docs/data-provenance.md).
