# Work Log

Chronological record of notable decisions and merges.

## 2026-10-07

- `lindhard-py` (#29): Python bindings with pyo3 and maturin. The classes
  wrap the library's own input structs (`lindhard::input`), so the TOML schema
  has one definition and `Run.from_toml()` / `Run.to_toml()` round-trip with
  the CLI. To make results identical by construction, `lindhard-cli` gained a
  library target (`sim::simulate`, the tally and the summary/CSV writers) that
  both the binary and the bindings call. `Run.run()` releases the GIL. The
  `numpy` crate is BSD-2-Clause and `target-lexicon` (a pyo3 build
  dependency) is Apache-2.0 WITH LLVM-exception; both are added to the
  `deny.toml` allow list.

## 2026-10-04

- Project started: a clean-room, MIT-licensed Rust engine for ion
  implantation and low-energy electron transport. The repo stays private until
  M0. The starting survey is in `docs/prior-art.md`.
- Milestone order: M0 amorphous ion core → M1 electron engine → M2
  crystalline implant → M3 reach (operator may reorder).
- Made public the same day (operator), ahead of the original "after M0"
  plan: private-repo Actions minutes were blocked by the org's billing
  limit, and public repos run standard runners free. A leak-pattern scan of
  the tree, the full history and the issues was clean before the flip. CI now
  runs on Linux x86-64, Linux arm64 and macOS, as fasterhenry's does.
