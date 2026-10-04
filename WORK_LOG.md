# Work Log

Chronological record of notable decisions and merges.

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
