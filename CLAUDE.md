# lindhard

A clean-room, MIT-licensed Rust engine for Monte Carlo transport of ions and
electrons in matter. Before writing code, read `CONTRIBUTING.md`, in
particular the **clean-room license tiers** (what may be ported, what may only
be read about in papers, and what may never enter the tree) and the
**disclosure** section. Then read `docs/architecture.md`.

- Physics comes from published papers. Every model and constant cites its
  source in a doc comment.
- **Data is not licensed; code is** ([operator ruling, 2026-11-10](https://github.com/2AMLogic/lindhard/issues/629)). The goal is to add to the
  digital commons: lindhard aims to be a complete open replacement for
  SRIM/TRIM. Published numbers from any source (papers, the SRIM book's
  tables, srim.org plots and curves, SRIM/TRIM output, ICRU, NIST, other
  codes' data files) may enter the tree as data. Every dataset gets a row in
  `docs/data-provenance.md` citing its source (edition, table and page; or
  program, version and file) and labelling its origin.
- The clean room is about **code**: never read, port or transcribe code from
  Tier B or Tier C programs (`CONTRIBUTING.md`).
- Keep validation honest: a comparison of lindhard against data that came
  from the same code is consistency, not independent validation, and a model
  fitted to a code's output is a named calibration, not physics. Page images
  of books stay out of the tree; commit the numbers.
- Determinism across thread counts is a tested invariant. Do not break it.
- This repo is public. Keep it free of anything that is not already
  published physics, and of any description of downstream applications.

<!-- BEGIN LOOM ORCHESTRATION -->
This repository uses [Loom](https://github.com/rjwalters/loom) for AI-powered development orchestration — see the Loom repository for the full guide (roles, labels, worktrees, configuration). When installed, Loom also writes a locally-substituted copy of that guide to `.loom/CLAUDE.md`.

Work is coordinated through `loom:` labels on issues and pull requests, and the same roles run either under `loom-daemon` or by hand in an attended session — daemon mode is optional. Create the labels once with `.loom/scripts/sync-labels.sh` (an install ships `.github/labels.yml` but does not create the labels on the forge). A pull request ready for review carries `loom:review-requested`; Judge reviews it and applies `loom:pr` (approved) or `loom:changes-requested`; Doctor fixes a `loom:changes-requested` pull request and returns it to `loom:review-requested`. Only a `loom:pr` pull request gets merged, and always via this repo's merge script (`.loom/scripts/merge-pr.sh`) — never a raw forge merge command such as `gh pr merge`. Full state machine: `.loom/docs/label-state-machine.md`.
<!-- END LOOM ORCHESTRATION -->

<!-- BEGIN REPO-SKILLS -->
This repository has [Repo Skills](https://github.com/rjwalters/repo) v0.25.10 installed —
general repository hygiene and environment commands invoked as `/repo:<command>`. Run
`/repo:help` for the command list, or see `.claude/skills/repo/SKILL.md` for the full
guide. Hygiene commands apply safe, reversible fixes by default and report each
change; run with `--ask` to review first, and `--prune` to allow irreversible
removals. Managed by `install.sh` — edit outside the markers only.
<!-- END REPO-SKILLS -->
