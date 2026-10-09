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
| `oracles/run_electron.py`, `oracles/electron_problems.json`, `oracles/geant4_microelec/` | 2 | Electron comparison: runs `lindhard` and user-installed Nebula and Geant4 MicroElec (through our own Geant4 application) on matched electron problems ([docs](../docs/validation.md#electron-oracles-nebula-and-geant4-microelec-150)) |
| `oracles/summaries/` | 2 | Committed comparison summaries (scalar summary metrics only; [format](oracles/summaries/README.md)) |
| `oracle-runs/` | 2, 3 | Raw runs, regenerated locally; **gitignored** |
| `data/` | 3 | Published measurements with citations ([schema](data/README.md)); `data/digitize/` holds the scripts that digitized figures |
| `experiments/run.py` | 3 | Checks each dataset's provenance, runs `lindhard` on it with the defaults and with the stopping inputs varied; writes `experiments/results.json` |
| `experiments/backscatter.py`, `experiments/backscatter/` | 3 | Electron backscatter coefficient vs measurements (#148): checks `data/backscatter/`, runs the committed `[electron]` inputs at 1 to 30 keV and the elastic-correction sensitivity; writes `experiments/backscatter_results.json` (`run.sh --backscatter`) |
| `experiments/se_yield.py` | 3 | Secondary-electron yield δ(E) of Al and Cu against the measured sets of `data/se_yield/` (#149): checks them, runs `lindhard` electron runs, writes `experiments/se_yield_results.json`, evaluates the initial bounds (`--check`) and prints the tables spliced into `docs/validation.md` |
| `lib/lindhard_cli.py` | 2, 3 | Builds and drives the `lindhard` command |
| `update_docs.py` | all | Splices the tables into `docs/validation.md` between marker comments |
| `check_manual_coverage.py` | docs | Fails if a variant of a model-selecting enum is not named in the physics manual (`book/src/models/`); run in CI |
| `book_walkthrough.py` | docs | Runs `examples/b_5keV_si.toml` once and writes (or with `--check`, compares) the output excerpts the user guide quotes; run in CI |

The harness scripts use the standard library only (Python 3.9 or newer). The
one-off digitizing scripts in `data/digitize/` also need numpy and Pillow;
the harness never runs them.

## Unit tests

Three stdlib `unittest` files guard the harness itself (level-3 provenance
enforcement, secondary-electron yield evaluation, and electron table-coverage
summing). They need no build, simulation or third-party package and run in
under a second. CI runs them in the `validation-docs` job, and
`validation/run.sh` runs them first:

```sh
python3 validation/experiments/test_run.py
python3 validation/experiments/test_se_yield.py
python3 validation/oracles/test_run_electron.py
```

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
