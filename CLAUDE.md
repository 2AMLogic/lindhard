# lindhard

A clean-room, MIT-licensed Rust engine for Monte Carlo transport of ions and
electrons in matter. Before writing code, read `CONTRIBUTING.md`, in
particular the **clean-room license tiers** (what may be ported, what may only
be read about in papers, and what may never enter the tree) and the
**disclosure** section. Then read `docs/architecture.md`.

- Physics comes from published papers. Every model and constant cites its
  source in a doc comment.
- Every dataset gets a row in `docs/data-provenance.md`. No SRIM-derived
  tables, ever.
- Determinism across thread counts is a tested invariant. Do not break it.
- This repo will go public after M0. Keep it free of anything that is not
  already published physics, and of any description of downstream
  applications.
