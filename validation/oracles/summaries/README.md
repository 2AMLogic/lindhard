# Oracle comparison summaries

One JSON file per (oracle, problem), written by `../run.py` and committed by
hand after review. They hold **scalar summary metrics only**, never an
oracle's tables, curves, particle lists or input files
(`CONTRIBUTING.md`, Oracles); those stay in the gitignored
`validation/oracle-runs/`.

Format `lindhard-oracle-summary/2`:

| Field | Contents |
|---|---|
| `problem` | Id of the matched problem in `../problems.json` |
| `oracle`, `oracle_version`, `oracle_commit`, `oracle_source`, `oracle_build`, `oracle_license` | Which program, at which upstream commit, how it was built, and its licence tier |
| `lindhard_version`, `date`, `host` | Our side, the run date, the machine and its load average at start |
| `ions`, `threads` | Histories and worker threads on each side |
| `lindhard`, `oracle_values` | Each side's metrics of the problem (`rp_nm`, `drp_nm`, `backscatter`, `sputter_yield`, `ions_per_s`, `ions_per_s_marginal`) with standard errors; `null` where a code has no comparable quantity |
| `comparison` | lindhard relative to the oracle: `*_rel_diff` for moments and yields, `backscatter_abs_diff`, each with `*_se` and `*_z` (difference / combined standard error), and `ions_per_s_ratio` / `ions_per_s_marginal_ratio` (> 1: lindhard faster) |
| `timing` | How the two speeds were measured |
| `lindhard_settings`, `oracle_settings` | The physics each side ran with, as our own summary of the options (not the oracle's input file) |
| `matched`, `mismatches` | Choices made identical, and every choice that could not be matched |

Format `/1` (no commit, no oracle values) was never committed.

## Electron oracles: `lindhard-electron-<problem>.json`, `nebula-<problem>.json`, `geant4_microelec-<problem>.json`

Written by `../run_electron.py` on the problems of `../electron_problems.json`
(#150) and committed by hand after review; `validation/update_docs.py` builds
the electron table of `docs/validation.md` from them.

Format `lindhard-electron-run/1` (`lindhard-electron-<problem>.json`, our side
alone): `problem`, `lindhard_version`, `date`, `host`, `histories`, `batches`,
`values` and `std_err` of the four metrics (`eta`, `delta`,
`primary_depth_nm`, `r50_nm`; standard errors by batch means),
`primaries_stopped`, `wall_s`, `materials` (the cstool commit and the SHA-256
of the parameter and ELF files read at run time; no values from them) and
`lindhard_settings` (the run's own model and transport metadata).

Format `lindhard-oracle-electron-summary/1` (`<oracle>-<problem>.json`): the
fields of `lindhard-oracle-summary/2` above where they apply (`oracle*`,
`lindhard_version`, `date`, `host`, `histories`, `lindhard`,
`oracle_values`, `oracle_settings`, `matched`, `mismatches`), plus
`comparison` (`eta_abs_diff`, `delta_abs_diff`, `primary_depth_rel_diff`,
`r50_rel_diff`, each with `*_se` and `*_z`), `tolerance` (the #150 rule and,
where it applies, pass or fail per metric), `wall_s` and, for Nebula,
`materials`. A metric a code has no value for is `null`.

## Electron benchmarks: `bench-electron-<problem>.json`

Written by `../bench_electron.py` on the same problems (#152) and committed by
hand after review; `bench_electron.py --update-docs` builds the electron
tables of `docs/benchmarks.md` from them. Timings only, no physics metric.

Format `lindhard-electron-bench/1`: `problem`, `date`, `lindhard_version`,
`build` (profile, rustc), `host` (operating system, architecture, CPU model,
logical CPUs and physical cores; no host name or other identity),
`load_average` at the start and end of the problem, `materials` (as in the
electron summaries above: the cstool commit and SHA-256s, no values),
`settings` (energy, element, seed, the lindhard physics block of
`../electron_problems.json`, the tally grid), `lindhard` (the output of
`lindhard/examples/electron_scaling.rs`: histories, seed, chunk size,
repeats, the table build times, and per thread count the transport wall
times, their median, electrons/s, speedup, parallel efficiency and the
SHA-256 of the tally report, with `deterministic` true only if every digest
is equal), `nebula` (`null`, or Nebula's version, commit, build, threads,
process wall time, its own `Simulation` and material-loading times and
electrons/s), `comparison` (`null`, or the electrons/s ratio at Nebula's
thread count, > 1: lindhard faster) and `caveats`.

## Level-3 context: `rustbca-sputter_ar_<target>.json`

Format `lindhard-oracle-sputter-summary/1`, written by
`validation/experiments/run.py --rustbca` (RUSTBCA_BIN set; `--rustbca-target T`
limits it to one target) and committed by hand: RustBCA's Ar sputter yield of
Cu (#69), Si, Ag and Au (#70), with weak collisions 0 and 3, at the energies
of the measured level-3 data (`validation/data/sputtering/`), run through
this directory's RustBCA adapter with the settings of `ar_1keV_cu_ed_es`
(`E_d = E_s` = lindhard's tabulated cohesive energy of the target). `rows`
holds per energy only scalar yields (`rustbca_k0`, `rustbca_k3`, Poisson
standard errors) next to lindhard's; `oracle_settings` and `mismatches` are
as in the level-2 summaries. Context for docs/validation.md section 3, not a
level-2 comparison: `update_docs.py` leaves them out of the level-2 table.
