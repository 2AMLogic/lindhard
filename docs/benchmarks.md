# Benchmarks

Performance is a tracked number. The benches are [criterion](https://crates.io/crates/criterion)
micro-benchmarks and an ions/s throughput benchmark. They are **local only**:
they never run under `cargo test`. CI only checks that they compile
(`cargo clippy --workspace --all-targets` builds every bench target).

## Running

```sh
cargo bench -p lindhard                      # everything (several minutes)
cargo bench -p lindhard --bench scattering   # one file
cargo bench -p lindhard --bench history -- 'history_'   # filter by name
cargo bench -p lindhard --bench history -- --quick      # fast, rough
cargo bench --no-run                         # compile only
```

Close other work first: on a loaded machine the numbers are inflated and
noisy, and the thread-scaling benchmark is the most sensitive to this.
All inputs are deterministic (fixed seeds, fixed ion counts, one shared
scattering table), so run-to-run differences are timing noise only.

## What is measured

| Bench file | Group | Measures |
|---|---|---|
| `scattering.rs` | `theta_1024_points` | Table lookup vs magic formula vs Gauss-Mehler quadrature for 1024 (eps, beta) points |
| `stopping.rs` | `stopping_256_energies` | `ElectronicStopping::stopping` over a 256-point log energy sweep for Lindhard-Scharff, Oen-Robinson, the equipartition mix and Bethe-Bloch, plus the Oen-Robinson local loss |
| `history.rs` | `history_100_ions_single_thread` | One history at a time (`Bca::history`) for B 5 keV, As 50 keV into Si and Ar 1 keV into Cu, 7 degrees off normal |
| `history.rs` | `cascade_and_free_path_B_5keV_Si_1000_ions` | Full cascade vs ions only (`follow_recoils`), constant vs energy-dependent free path |
| `history.rs` | `tally_overhead_B_5keV_Si_1000_ions` | The light `SummaryTally` vs the full `IonTally` (histograms, moments, damage, escapes) on the same histories |
| `history.rs` | `throughput_<problem>_10k_ions` | 10^4 ions through `Bca::run` on explicit rayon pools of 1, 2, 4, ... up to the machine's parallelism |

Free-path selection is measured through the engine (constant vs
energy-dependent convention), since the selection code is private to
`Bca`. Targets use illustrative `E_d` = 15 eV, cutoffs 5 eV (primary) and 2 eV
(recoil), and the ZBL universal potential with Lindhard-Scharff stopping; these
are benchmark settings, not recommendations.

## Reading the results

* Criterion prints `time:` as `[low estimate high]` (a confidence interval)
  and, where a throughput is declared, `thrpt:` in elements per second. For
  `history.rs` and `stopping.rs` an element is one ion or one evaluation, so
  **`thrpt` is ions/s** (or evaluations/s).
* Compare **ratios within one run** (cascade vs ions only, 1 thread vs N
  threads); absolute numbers belong to one machine.
* Thread scaling: ideal is `ions/s(N) = N * ions/s(1)`. Because the engine is
  deterministic across thread counts (the stream is keyed on the primary
  index, chunk size is fixed), the work is identical at every N; a shortfall
  is scheduling, memory or a loaded machine, not different physics.
* Cost per ion varies by orders of magnitude between problems (cascade size
  grows with the primary energy and the mass ratio), so ions/s is only
  comparable between codes on **matched problems and settings** (cascades on
  or off, cutoffs, mean-free-path convention).

## Sample numbers

Taken at commit `9aed2ab` (clean-room tree plus this bench set), 2026-10-05,
Apple M3 Ultra (28 logical CPUs), macOS, rustc 1.98.1, `cargo bench -- --quick`.
**The machine was heavily loaded (load average above 24) during these runs**, so
treat them as order-of-magnitude and rerun on a quiet machine before quoting.
They are not assertions and nothing checks them.

| Benchmark | Result |
|---|---|
| theta, table lookup | about 35 ns per point |
| theta, magic formula | about 1.05 us per point |
| theta, quadrature | about 2.75 us per point |
| Lindhard-Scharff stopping | about 73 ns per evaluation |
| Oen-Robinson / equipartition mix | about 66 ns per evaluation |
| Bethe-Bloch (proton in Si) | about 16 us per evaluation |
| One history, B 5 keV Si (cascade) | about 4.1 k ions/s, 1 thread |
| One history, As 50 keV Si (cascade) | about 0.41 k ions/s, 1 thread |
| One history, Ar 1 keV Cu (cascade) | about 23 k ions/s, 1 thread |
| B 5 keV Si, ions only | about 20 k ions/s (about 5x the cascade) |
| B 5 keV Si, `IonTally` vs `SummaryTally` | about 5% slower |
| B 5 keV Si, 10^4 ions, 1 thread | about 4.5 k ions/s |
| B 5 keV Si, 10^4 ions, 28 threads | about 50 k ions/s (loaded machine) |

Observations: the lookup table is about 30 times cheaper than the magic
formula, so it is not a hot spot; with cascades on, the cost is dominated by
the number of collisions, which the ions-only row shows is mostly recoil
cascade work. Bethe-Bloch is about 200 times the cost of Lindhard-Scharff per
evaluation, which is irrelevant for keV implants but would matter if it were
called per collision at MeV energies.

## Not yet covered

* **Comparison against oracles (RustBCA, OpenTRIM)** needs the oracle runner
  of the validation harness and matched inputs; tracked as the remainder of
  issue #9. Rules for running oracles and what may be committed (our numbers
  and ratios only, never oracle output or source) are in
  [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §Oracles.
* **A profile (flamegraph) of the 10^4-ion run**, with the top hot spots named,
  to decide whether SoA/SIMD work is worth doing. Also the remainder of #9.
