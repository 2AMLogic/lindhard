# lindhard-cli

The `lindhard` binary: one TOML input in, `summary.json` and CSV profiles out.
It is a thin front end over the [`lindhard`](https://github.com/2AMLogic/lindhard/blob/main/lindhard/README.md) library.

| Subcommand | What |
|---|---|
| `lindhard check <input.toml>` | Parse and validate the input, no transport |
| `lindhard run <input.toml> --out <dir>` | Run it and write `summary.json` and the CSV profiles to `<dir>`; `--ions`, `--seed` and `--threads` override the input (threads never change the results) |

From the repository root:

```sh
cargo run --release -p lindhard-cli -- run examples/b_5keV_si.toml --out out/b_5keV_si
```

The input schema, the output files and the reproducibility guarantees are in
[`docs/cli.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/cli.md). Example inputs are in
[`examples/`](https://github.com/2AMLogic/lindhard/blob/main/examples/README.md).

| Source | What |
|---|---|
| `src/main.rs` | Argument parsing, input loading and the run |
| `src/output.rs` | `summary.json` and the CSV profiles |
| `src/tally.rs` | The tally the CLI runs, built on the library's tallies |
| `build.rs` | Stamps `git describe` into `--version` |
| `tests/` | End-to-end tests of the binary |
