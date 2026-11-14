# Contributing

## The one hard rule: clean room for code

lindhard is released under MIT. **Software code is licensed; data is not**
([operator ruling, 2026-11-10](https://github.com/2AMLogic/lindhard/issues/629)). So the clean room governs code: the code, comments and test
harnesses in this tree must be either our own work from published methods, or
ported from a permissively licensed source with its notice kept. We learn from
other codes, but **how** we may learn from their code depends on each code's
license. Data (measured or computed numbers, from any source) may enter the
tree with provenance; see "Data" below. The goal is to add to the digital
commons, and lindhard aims to be a complete open replacement for SRIM/TRIM.

### Tier A: permissive. Read, learn, and port with attribution

| Code | License | Useful for |
|---|---|---|
| OpenTRIM (ir2-lab) | MIT | Architecture, 3D grid tallies, damage bookkeeping, GUI/HDF5 output ideas |
| DISPLATH (Julia) | MIT | BCA on explicit atom positions (a route to crystalline targets) |
| Nebula + cstool (TU Delft) | BSD-3-Clause | Low-energy electron MC, dielectric-function cross sections, GPU layout |
| ESPNN | MIT | A stopping model trained on experimental data |
| PyEPICS | BSD-3-Clause | EEDL/EADL/EPDL data access |

Porting from a Tier A source is allowed. Keep the upstream copyright notice,
add the source to `THIRD_PARTY_LICENSES.md`, and put a comment at the port
site naming the upstream file and commit. OpenTRIM's electronic stopping
tables come from SRIM-2013; as data they are admissible like any SRIM table
(see "Data"), cited to SRIM-2013.

### Tier B: copyleft (GPL, LGPL, AGPL). Papers, manuals and outputs only

RustBCA, iradina, IM3D, libdEdx, CATIMA, JIBAL, Corteo, PenRed, EGSnrc.

Code derived from these could never ship under MIT. You may use:

- **The papers and manuals.** Most of these codes are well documented in
  journals (RustBCA in JOSS, iradina in NIM B, CATIMA and libdEdx in their
  docs). Implement from the paper.
- **Their input/output formats**, as a convenience. A file format is not code.
- **The programs themselves as oracles.** Run them unmodified on your own
  machine and compare our results against theirs (see "Oracles" below).

You may not have their source open while writing code here, or transcribe it
from memory. Numbers they publish or ship (tables, constants) are data: they
may enter under "Data", cited to the program, version and file.

### Tier C: closed or non-commercial. Literature and data; never code

SRIM/TRIM (and every table SRIM produces), SDTrimSP, TRIDYN/TRI3DYN, IMSIL,
Crystal-TRIM, MARLOWE, IIS, CASINO, PHITS, FLUKA, MCNP, ELSEPA, DPASS/PASS,
CasP.

Use their published papers and manuals; never their code. Their output and
their printed tables are data, and may be committed with provenance (see
"Data"). That includes SRIM's stopping, range and straggling tables, the
tables printed in the SRIM book, and the curves on the srim.org plots.

**Label what came from where.** A dataset taken from a code's output names
the code and its version. When lindhard uses such a dataset as model input
(for example a SRIM stopping table selected as the stopping model), a
comparison against that same code is a consistency check, not independent
validation, and the report says so. Our physics models stay implementations
of the published physics: if one is fitted or tuned to a code's output, it
becomes a named, versioned calibration with its fit recorded (as for the
sputtering tuning sets in `docs/data-provenance.md`), never a silent change.

**The SRIM book.** J. F. Ziegler, J. P. Biersack and M. D. Ziegler, *SRIM:
The Stopping and Range of Ions in Matter* (SRIM Co., 2015; ISBN
978-0-9654207-1-6; earlier printings 1984 to 2008). Page and equation
numbers cited in this tree are from the 2015 printing. Treat it as Tier C
literature:

- **Citable:** the text, the derivations, and the published formulas and
  constants (for example the universal screening function and the
  magic-formula constants). Cite them by edition, page and equation, and
  check every digit against the page itself, not against memory.
- **Data:** the printed stopping, range and straggling tables, the SRIM
  curves on its plots, and the fitted electronic-stopping coefficient sets.
  These may be entered into the tree as data and used as model input. Enter
  each value twice, independently, and diff the two entries; cite edition,
  table and page.
- **Never:** the bundled software and any program listing (see "Mixed
  sources").
- **Page images stay out of the tree.** Commit the numbers and their
  citations; keep scans or photos locally (for example under the gitignored
  `validation/oracle-runs/`).

**The srim.org experimental-data compilation** (operator ruling, 2026-10-28).
srim.org publishes plots of measured stopping powers (about 22,000 points)
and a bibliography of the roughly 2,300 experimental papers behind them
(`PlotCit.htm` and its per-system citation files). The measured points and the
citations are facts reported by those papers, not SRIM output. They may enter
the tree under these conditions:

- Each point is cited to its **original paper** (the citation the plot's
  legend names) and flagged as *digitized from the srim.org compilation plot*.
  It carries a digitization uncertainty, and it is replaced by the paper's own
  value once the paper has been read.
- SRIM's own computed values on the plots (the fitted curve, the
  Experiment/Theory ratio panel, and the normalizations of the heavy-ion
  aggregate plots) are data too ([operator ruling, 2026-11-10](https://github.com/2AMLogic/lindhard/issues/629), superseding the 2026-10-28
  intermediate-only restriction recorded on #486). A digitized curve is
  labelled as SRIM output, never as a measurement. A measured point recovered
  through a SRIM ratio or normalization records the ratio used and is flagged
  as reconstructed.
- The site's history and commentary are welcome too, as facts restated in our
  own words with citations. The site's files (images, PDFs, page text) stay in
  the operator's archive rather than being copied into the tree.
- This data is released under CC0 (`validation/data/stopping/compilation/LICENSE`),
  separate from the MIT code. Every point keeps its citation regardless.

**Mixed sources.** Some published reports bundle a program listing with the
paper. The listing is Tier C even though the report is citable. Exclude the
listing pages before any OCR or text extraction, not after.

- Moller and Eckstein, *TRIDYN*, report IPP 9/64 (1988): read only the report
  body, PDF pp. 1-47. Never open Appendix 1 (the TRIDYN program listing, PDF
  pp. 48-86).

### Oracles

Comparing against third-party programs is encouraged, and the harness under
`validation/` drives them. Their raw outputs are regenerated locally into
`validation/oracle-runs/`, which is gitignored, to keep the tree small. Commit
our own **summary metrics** of a comparison (for example "Rp within 3% of
oracle X at version Y"). When an oracle's output is wanted as data (a fixture,
a selectable table), commit it as a dataset with a provenance row naming the
program and version. This includes SRIM: run it locally, or use the tables in
its book, naming the SRIM version or the book's edition and page.

### Data

Every number that enters the tree as data needs a provenance entry in
[`docs/data-provenance.md`](docs/data-provenance.md). Its origin must be
stated and cited: computed from a published formula, an experimental
measurement cited to its paper, or a table or output from a named source (book
edition, table and page; database and version; program, version and file).
"It came from a well-known table" is not provenance; "SRIM book, 2015
printing, table and page" is.

### CI enforces part of this

The `clean-room` job greps every file for copyleft license notices and fails
the build if it finds one. It catches copied code files, not paraphrase.
Paraphrase is on you.

## Disclosure

This repository is **public**. Everything in its history, its issues and its
PRs is published the moment it is pushed or posted. This project implements
*published* physics. If you think you have devised a method that is genuinely
new (one not found in the literature), stop and raise it with the operator
**before** committing it, filing an issue about it, or describing it in a PR.
Do not describe the operator's downstream applications in this repo either.

## Releasing

Releases are cut by the operator; the workflow in
[`.github/workflows/release.yml`](.github/workflows/release.yml) does the rest.
Both crates share one version.

1. Set `version` under `[workspace.package]` in `Cargo.toml` and the `version`
   of the `lindhard` entry under `[workspace.dependencies]` to the new
   version, and refresh `Cargo.lock`.
2. In `CHANGELOG.md`, rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD`
   and start a fresh `## [Unreleased]` above it. The section is the release
   notes.
3. Merge that change, then tag the merge commit with an annotated tag and
   push it: `git tag -a vX.Y.Z -m "lindhard X.Y.Z" && git push origin vX.Y.Z`.
   The tag must be annotated: `lindhard-cli/build.rs` stamps `git describe`
   into `--version` and the run report, and `git describe` ignores
   lightweight tags.
4. The workflow's preflight fails, before any build, if the tag differs from
   the workspace version or the changelog section is missing. Otherwise it
   builds `lindhard-<version>-<target>` archives (`.tar.gz`, `.zip` on
   Windows) for x86_64 and aarch64 Linux (static musl), x86_64 and aarch64
   macOS, and x86_64 Windows, then creates the GitHub Release with the
   archives and `SHA256SUMS`.
5. The `publish` job waits for approval of the protected `release`
   environment, then publishes `lindhard` and `lindhard-cli` to crates.io.
   Approve it only after checking the release assets.

Any pull request that changes `.github/workflows/release.yml` runs the
workflow as a dry run in PR CI: all five builds, `SHA256SUMS`, and the
archive-contents check, with the archives uploaded as workflow artifacts.
The workflow can also be run by hand from the Actions tab
(`workflow_dispatch`), which is likewise a dry run. Only a pushed `v*` tag
creates a release or publishes anything; on every other event the `release`
and `publish` jobs are skipped, and a pull request (including one from a
fork) gets a read-only token and no environment secrets. Release builds skip
the Rust build cache, so tag-built binaries come from a clean build.

One-time setup: create the `release` environment with required reviewers and
give it the secret `CARGO_REGISTRY_TOKEN` (a crates.io API token). The publish
job fails with a clear message if the secret is missing.

## Clippy pin and canary

The required clippy gate in [`ci.yml`](.github/workflows/ci.yml) runs on an exact
Rust release, `CLIPPY_TOOLCHAIN` (a full `1.x.y`), so a new Rust stable cannot
turn unrelated pull requests red. `cargo fmt` and `cargo test` stay on floating
stable, and `rust-toolchain.toml` is unchanged. The pin is independent of the
MSRV and of the `package` job's Cargo requirement.

The root `package.json` scripts run that same pinned gate: `npm run lint`
and `npm run check:ci` invoke `cargo +1.x.y clippy --workspace --all-targets
-- -D warnings`, so they need the pinned toolchain installed
(`rustup toolchain install 1.x.y --profile minimal --component clippy`).
Without it they fail loudly with rustup's "toolchain ... is not installed"
rather than quietly falling back to floating stable. `cargo fmt`,
`cargo test`, `cargo doc` and `cargo build` in those scripts stay on
floating stable, matching what CI gates for each.

Drift is reported by the weekly [`clippy-canary`](.github/workflows/clippy-canary.yml)
workflow (also runnable by `workflow_dispatch`). It runs the same
`cargo clippy --workspace --all-targets -- -D warnings` on `stable` and `beta`,
is not a pull-request workflow, and is not a required check.

**Reading the canary.** A failing leg opens, reopens or comments on the single
tracking issue (identified by the hidden `<!-- clippy-canary-tracking -->` marker
in its body). Each comment names the failing channel(s), the commit SHA and the
run URL. A fully green run comments "recovered" and closes the issue. If the run
itself is red with a reporting error (API failure, or more than one issue
carrying the marker), fix that first: remove the marker from all but one issue,
then re-run. Never create a second tracking issue by hand.

**Responding to drift.**

1. Reproduce locally with the channel the canary named, for example
   `cargo +stable clippy --workspace --all-targets -- -D warnings`.
2. Fix the new lints (or, with a stated reason, add a narrowly scoped `#[allow]`).
   Fixes that pass on both the current pin and the newer toolchain can land first.
3. Bump the pin in a pull request: set `CLIPPY_TOOLCHAIN` in `ci.yml` to the
   current stable release (check it exists with
   `rustup toolchain install 1.x.y --profile minimal --component clippy`, or at
   `https://static.rust-lang.org/dist/channel-rust-1.x.y.toml`), together with
   the lint fixes from step 2. Update the `cargo +1.x.y clippy` pin in the root
   [`package.json`](package.json) in the same commit: the `workflow-checks` job
   (`test_clippy_canary.py`) fails the bump pull request until the two agree.
4. Before merging the bump, confirm with that exact toolchain:
   `cargo +1.x.y clippy --workspace --all-targets -- -D warnings`, and that the
   PR's `Rust` jobs are green on all three operating systems.
5. After the bump merges and the next canary is green, the tracking issue closes
   itself.

To check the reporting path live, dispatch the canary with `simulate-failure`
set to true, confirm the tracking issue is created or updated, then dispatch a
normal run and confirm the recovery comment and closure. The reporting logic and
the workflow wiring have offline tests:
`python3 .github/scripts/test_clippy_canary.py` (needs PyYAML).

## Dependency advisories

The weekly [`advisories`](.github/workflows/advisories.yml) workflow (also
runnable by `workflow_dispatch`) runs `cargo deny check advisories` against the
RustSec database, with the policy in the `[advisories]` table of `deny.toml`. It
is not a pull-request workflow and not a required check, on purpose: the
database changes without any commit here, so a PR-time check would fail
unrelated PRs on a newly published advisory. `ci.yml` runs only
`cargo deny check licenses sources`.

The workflow pins its own cargo-deny (`CARGO_DENY_VERSION` in `advisories.yml`),
newer than the one in `ci.yml`: 0.18.2 cannot load the current database
(`unsupported CVSS version: 4.0`). Keep `deny.toml` parseable by both.

**Reading the report.** A failing run opens, reopens or comments on the single
tracking issue (identified by the hidden `<!-- advisories-tracking -->` marker in
its body); each comment carries the commit SHA and run URL, and the run log lists
the advisory IDs. A green run comments "recovered" and closes the issue. A red
run with a reporting error (API failure, or more than one issue carrying the
marker) needs that fixed first: remove the marker from all but one issue, then
re-run. Never create a second tracking issue by hand.

**Responding.** Reproduce with `cargo deny check advisories`, then update the
affected crate (`cargo update -p <crate>`), or replace it. If an advisory truly
does not apply, add it to `ignore` in `deny.toml` in a pull request, with the
advisory ID and a comment giving the reason and what would let the entry go.
Never ignore without both.

To check the reporting path live, dispatch the workflow with `simulate-failure`
set to true, confirm the tracking issue is created, dispatch it again and confirm
no second issue appears, then dispatch a normal run and confirm the recovery
comment and closure. The reporting logic and the workflow wiring have offline
tests: `python3 .github/scripts/test_advisories_report.py` (needs PyYAML). They
are not run by `ci.yml`.

## Everything else

- Rust stable, `cargo fmt`, `cargo test` green, `cargo clippy -- -D warnings` clean on
  the pinned clippy toolchain (see "Clippy pin and canary" above), and
  `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` clean.
- The `#[ignore]`d statistical and slow tests are not in the PR gate. The weekly
  [`statistical`](.github/workflows/statistical.yml) workflow runs them with
  `cargo test -p lindhard --release -- --ignored` (also by hand via
  `workflow_dispatch`).
- `#![forbid(unsafe_code)]` in the library stays. SIMD goes through safe
  crates (`wide`), not intrinsics.
- No C/Fortran dependencies in the default build.
- Every dependency must be permissively licensed and come from crates.io.
  `cargo deny check licenses sources` enforces both in CI, and `deny.toml`
  holds the allowlist. Known-vulnerable, yanked and unmaintained crates are
  reported weekly, not gated per PR (see "Dependency advisories" above).
- **Reproducibility is a feature.** Every random draw comes from a counter-based
  stream keyed on (seed, particle index), so a run gives the same bits on 1
  thread or 64. A PR that breaks this needs a test showing why it must. The invariant is
  tested by [`lindhard/tests/determinism.rs`](lindhard/tests/determinism.rs)
  (bit-identical tallies on 1, 2 and 8 threads); the helpers are in
  `lindhard::rng`.
- Physics constants and model choices cite their source in a doc comment.
- User-visible changes (behaviour, input/output formats, CLI flags, public API,
  packaging) add a line under `## [Unreleased]` in [`CHANGELOG.md`](CHANGELOG.md)
  ([Keep a Changelog](https://keepachangelog.com/) style). Internal refactors,
  tests and docs-only changes need no entry.
- The minimum supported Rust version is `rust-version` in the root
  `Cargo.toml`, and the CI `msrv` job checks the workspace on exactly that
  toolchain; raise it deliberately, in its own `CHANGELOG.md` entry, never as
  a side effect of using a newer language or std feature.
- One PR per issue. PRs are reviewed by Loom's Judge and merged by Champion,
  using merge commits.
