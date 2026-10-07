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
| `oracles/summaries/` | 2 | Committed comparison summaries (scalar summary metrics only; [format](oracles/summaries/README.md)) |
| `oracle-runs/` | 2, 3 | Raw runs, regenerated locally; **gitignored** |
| `data/` | 3 | Published measurements with citations ([schema](data/README.md)); `data/digitize/` holds the scripts that digitized figures |
| `experiments/run.py` | 3 | Checks each dataset's provenance, runs `lindhard` on it with the defaults and with the stopping inputs varied; writes `experiments/results.json` |
| `lib/lindhard_cli.py` | 2, 3 | Builds and drives the `lindhard` command |
| `update_docs.py` | all | Splices the tables into `docs/validation.md` between marker comments |
| `check_manual_coverage.py` | docs | Fails if a variant of a model-selecting enum is not named in the physics manual (`book/src/models/`); run in CI |
| `book_walkthrough.py` | docs | Runs `examples/b_5keV_si.toml` once and writes (or with `--check`, compares) the output excerpts the user guide quotes; run in CI |

The harness scripts use the standard library only (Python 3.9 or newer). The
one-off digitizing scripts in `data/digitize/` also need numpy and Pillow;
the harness never runs them.

## Clean-room rules that apply here

From [`CONTRIBUTING.md`](../CONTRIBUTING.md):

- Oracles are run **unmodified**, installed by the user, and never vendored.
  Adapters write an oracle's input from its *published documentation* (paper,
  manual), never from its source code or its example inputs and test
  fixtures. RustBCA is GPL (Tier B); OpenTRIM is MIT (Tier A), but its
  SRIM-2013 stopping tables are Tier C and are neither used for matched runs
  nor copied.
- Raw oracle output, and the input files written for an oracle, stay in
  `oracle-runs/`. Commit only scalar summary metrics, with the oracle's name,
  version and upstream commit. Install oracles outside this tree.
- No SRIM/TRIM output anywhere, including as reference data. SRIM may be run
  locally as a familiarity check; its output never becomes a fixture.
- Every experimental dataset is cited to its paper and has a row in
  [`docs/data-provenance.md`](../docs/data-provenance.md).
