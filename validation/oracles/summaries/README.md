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
