# Work Plan

Prioritized roadmap of upcoming work, maintained by the Guide role.

<!-- Maintained automatically by the Guide triage agent. Manual edits are fine but may be overwritten. -->

## Milestones

- **M0: amorphous ion core, validated.** Materials, potentials, scattering
  integral, electronic stopping, 1D layered BCA with full cascades, damage and
  range tallies, TOML/JSON CLI, a validation harness against analytic limits
  and oracles, and benchmarks.
- **M1: low-energy electron engine.** Partial-wave Mott elastic scattering,
  dielectric-function inelastic scattering, secondary electrons, layered
  resist/oxide/metal stacks, energy-deposition maps and PSF extraction.
  Validated against Nebula, Geant4 MicroElec and published η/δ/PSF data.
- **M2: crystalline implant.** Lattice targets, thermal vibration,
  tilt/twist, screen oxide, and dose-dependent damage leading to dechanneling
  and amorphization. Validated against published SIMS profiles.
- **M3: reach.** Dynamic composition, 2D/3D geometry, Python wheels, WASM.

Alongside the milestones run release engineering and docs, the permissive
stopping-data library, and performance work (release engineering #30 and
performance #36 closed in October 2026). Original order: M0 first. M1 Phase 1
can start alongside it (it needs only materials and the RNG). M2 follows M0's
engine: M2 Phase 1 may run alongside M0's remaining data verification and
validation datasets (#15, #51), which do not touch the code it builds on
(operator ruling, 2026-10-08, on #12).
The rest run as capacity allows.

## Where it stands

Status as of 2026-11-13, `main` @ `eff4bf6`. The engine work of every
milestone has merged. What is open is validation against published
measurements, data entry from printed sources, and reach (inputs, bindings).
Most of that waits on a source or a ruling, not on agent work.

| Area | Done and evidenced | Open |
|---|---|---|
| M0 amorphous ion core (#10) | Closed 2026-11-10; validated at levels 1 to 3 (`docs/validation.md`) | Recorded gaps: #15 (ASTM E521 `E_d`; C and Mn densities), #45 (Lenz-Jensen digits; Kr-C coefficients) |
| M1 electron engine (#11) | Engine, secondaries, PSF, exposure, charge tally (#546 delivered as a tally) | Gates #148 (η), #151 (PSF); #557 (depth-dose), #111, #162 wait on sources |
| M2 crystalline implant (#12) | Engine through step 23 (#513, #514) | SIMS validation #515 then #516; #251 (Oen-Robinson length) |
| M3 dynamic composition (#22) | Engine; validations #436, #437, #448, #465, #512; Ni datasets and Cu attribution (#543, #544) | #511 (parent of #543/#544) is still open; focused ion beams #554: #565, #561 |
| 2D/3D geometry (#25) | Voxel, CSG and mesh geometry in the library; voxel tallies; masked-edge and buried-edge validations (#435, #572) | No TOML input for a laterally varying target (`docs/cli.md`, "Direct Monte Carlo exposure"); #434 layout import; #121 |
| Bindings (#28) | Python package and WASM session and simulator | #469 real-browser check; the WASM session rejects `[[crystal]]` and `[dynamic]` |
| Stopping-data library (#33) | Epic closed 2026-11-13; loader, residual report, first primary-paper datasets (#483, PR #724) | #497 (heavy-ion plots), #486 (parent, parked on #497), #503 and #723 (papers), #487 (explorer) |
| SRIM/TRIM parity (#630) | Data policy #629; `stopping = "srim"` loader #632; gas targets #636; depth profiles #639; quick damage mode #640; PKA-spectrum primaries #641 | #637 (PR #650), #633, #658, #634, #638, #642, #643 |
| Photon channel (#548) | Part A implemented in PR #631 (open) | #562 until #631 merges; then #566, #567 |
| Damage outputs, optical data | #553 and #563 (versioned export), #552, #555 closed | Optical gaps listed in `docs/data-provenance.md` |

**0.1.0** (operator decision, 2026-11-09): the release is held until the M1
and M3 validation gates close (#148, #151, #511). It publishes a GitHub
release, the crate and the wheels (`CONTRIBUTING.md`, "Releasing"). The
website stays live but unannounced, so #469 does not gate it. Work merged
alongside ships with the tag and does not change the scope (ruling,
2026-11-10).

Changed in this revision (2026-11-13): #543 and #544 are closed, so #511 is
no longer waiting on them; Epic #33 and #483 are closed, so the stopping-data
library is no longer "after 0.1.0" as a whole; the SRIM/TRIM parity epic
(#630) was missing and is added; #549, #552, #553, #555, #563 and #572 are
closed and left the "alongside" table; the photon data question moved from
"which attenuation source" to the angular data in PR #631.

## Now

Can start today. Each item is one pull request.

| What | Done when | Issue |
|---|---|---|
| Structural clean-up already claimed: shared loaders, an `electron::beam` module, shared test and bench helpers, the stub `package.json` | Each issue's pull request merges with CI green | #750, #751, #752, #753, #755 |
| Opt-in Bohr energy-loss straggling in ion transport. `ion/stopping/straggling.rs` states that BCA transport does not sample it; the Bohr variance is already in the tree and cited | A `[physics]` key, off by default, with a thin-foil test: the transmitted energy variance matches the Bohr value within the stated statistical error; default results stay bit-identical; thread-count determinism test passes | none yet |
| TOML input for a laterally varying ion target built from the library's CSG solids (cuboid, cylinder, half-space) | `lindhard check` and `run` accept an example under `examples/`; `docs/cli.md` documents the tables; a test compares a one-cuboid target with the equivalent layer stack | none yet (Epic #25; lifts the block on #121) |
| `[[crystal]]` targets in the WASM session | `Simulation::new` accepts a crystal input; a `wasm-bindgen` test reproduces the native depth profile for `examples/b_5keV_si_crystal.toml` at the same seed | none yet (Epic #28) |

## Next

Unblocked by a named item above or by one merge.

| What | Done when | Waits on |
|---|---|---|
| `lindhard stopping-table` (#637) | PR #650 merged after a rebase (it conflicts with `main`) | A human merge: the pull request adds a dependency to `lindhard-cli/Cargo.toml`, which holds it |
| Lateral fate in `docs/cli.md` and `CHANGELOG.md` (#121) | The three criteria in #121 are met | The CSG-input item above |
| Close or restate #511 | The issue is closed, or its body says what remains after #543 and #544 | A re-read of the issue; both cited blockers are closed |
| Photon parts B and C (#566, #567) | Each issue's criteria | PR #631 merged (#562) |
| Straggling beyond Bohr (#638) | A Chu-type variant with paper-sourced coefficients, or a recorded decline | The Bohr item above for the transport hook, and a readable source |
| SRIM-format output writers and a TRIM.IN importer (#643, #642) | Each issue's criteria | Operator hold on both |

## Later

Direction, each with the question to answer first.

- **Electron and dynamic runs in the browser.** Which run drivers need a
  timer-free path first? `lindhard-wasm/README.md` records that the CLI
  driver's use of `std::time::Instant` panics on `wasm32-unknown-unknown`.
- **Layout cross-section import (#434).** Does `klayout-tools` export
  per-layer merged polygons in a versioned format? The issue forbids a second
  GDS or OASIS reader in this tree.
- **Stopping-data explorer in the book (#487).** Is the digitized compilation
  complete enough to show (#497)?
- **More wheel platforms.** `lindhard-py/README.md` lists musllinux and
  Windows arm64 as not covered. Is there a native runner or a cross build for
  each that the release workflow can use?
- **Kr-C magic-formula set.** Needs Table 4-1 of the SRIM book; checked
  against the quadrature like the other two sets.
- **Muffin-tin potential and a GPU backend.** Both are recorded deferrals
  with the conditions that lift them (`docs/muffin-tin-deferral.md`,
  `docs/gpu-feasibility.md`).
- **Default stopping model.** #630 adds SRIM-compatible options and changes
  no default; changing one is a separate decision.

## Blocked and waiting

Outside the Now table, no open issue can be started by an agent today. Every
other open issue waits on one of these.

| Input | Unblocks |
|---|---|
| A human merge of PR #650 (critical-file hold) | #637 |
| A ruling on the source of the form-factor and scattering-function data in PR #631 | #562, then #566, #567 |
| Scans of the SRIM book tables and pages (they stay out of the tree) | #633, #658, #634, the Kr-C set |
| A run on the machine that holds the srim.org mirror and the OCR tools | #497, then #486, #487 |
| Oen and Robinson, NIM 132, 647 (1976) | #251, then M2 accuracy |
| SIMS profiles from papers (B, P, As in Si; screen oxide) | #515, #516 |
| An openable measured resist PSF; a ruling on the CXRO tables for the PMMA ELF | #151 (release gate), #162 |
| A ruling on a Si core-region ELF source (Henke tables) and on running the #457 comparison | #148 (release gate) |
| A readable electron depth-dose or interface measurement | #557 |
| Oliva-Florio et al., Phys. Rev. B 35, 2198 (1987), or another measured yield against angle | #565; #561 is operator-held |
| Primary stopping papers (H, He, B, P, As in Si; He in Cu and Au) | #503, #723 |
| Chu (1976), Yang, O'Connor and Wang (1991), or the SRIM book's straggling text | #638 |
| Lenz (1932), Jensen (1932); Wilson, Haggmark, Biersack (1977); ASTM E521; a later CRC Handbook | the rest of #45 and #15 |
| A sourced high-l phase-shift method | #111 |
| An idle machine for benchmarks; a penn-full build measurement | #76, #177, #360 |
| A real-browser check of the simulator | #469, #709, then #467 |
| The post-merge failure-simulation run for the advisories gate | #199 |
| Operator rulings | #728 (shared blocks and formats), #739 (`statistical.yml` skip), #195 and #220 (guard reviews) |

## Not planned

- Reading, porting or transcribing code from SRIM/TRIM or any Tier B or
  Tier C program (`CONTRIBUTING.md`). Data from those sources is admitted
  with a provenance row (#629).
- A copyleft dependency anywhere in the build.
- A second GDS or OASIS reader (#434).
- Descriptions of downstream applications (`CLAUDE.md`).

<!-- guide:plan-body:start -->
## Operator Attention: Merge-Risk-Hold Pileup

Judge-approved PRs under a hold (`loom:blocked` or `loom:operator`) — only a human merge decision is missing.

- **#650**: feat: add stopping-table subcommand for stopping and CSDA range tables

## Very Important

`loom:very-important` (one counts per workspace); every stage takes these first.

- **#511**: [Epic #22] Single-element Ar → Cu and Ar → Ni yields: store a Ni dataset and explain the 0.67 Cu ratio without tuning

## Important

`loom:important`; taken before unlabelled work.

- **#148**: [Epic #11] Validation: backscatter coefficient η(E, Z) vs published measurements
- **#151**: [Epic #11] Validation: PSF fit vs a published resist PSF measurement
- **#486**: [Epic #33] Digitize the measured points from the srim.org stopping plots (cited, flagged, count-checked)
- **#497**: [Epic #33] Digitize heavy-ion aggregate stopping plots (STOP<Z>xx / STOPxx<Z>) with cross-view dedup
- **#637**: [Epic #630] `lindhard stopping-table`: stopping and range tables over an energy grid

## Ready

Human-approved issues ready for implementation (`loom:issue`).

- **#562**: [Photon channel #548] A: photon transport with cited attenuation data and per-layer photon dose
- **#773**: Finish the tests/common consolidation: 18 clusters of byte-identical test helpers remain
- **#785**: Share the 100 duplicated helper lines of tests/csg.rs and tests/mesh.rs via tests/common

## In Progress

Issues currently being built (`loom:building`).

- **#752**: Share one data-file loader across the input resolvers, and record the SRIM coefficient file's SHA-256 that is computed but never written
- **#753**: Move the electron beam model (BeamProfile, BeamSpread, beam serde macro) out of electron/transport.rs into its own electron::beam module
- **#758**: Collapse the duplicate output path in `lindhard exposure` into the shared OutputSet lifecycle
- **#763**: Root package.json clippy scripts run floating stable, not the pinned CI gate toolchain
- **#764**: Collapse the duplicated secondary-spawn tail in electron transport, and pass the Kieft-Bosch options as one value
- **#765**: THIRD_PARTY_LICENSES.md still points cstool's inelastic_axis attribution at lindhard-cli/src/electron.rs
- **#769**: [Epic #28] Accept `[[crystal]]` targets in the WASM simulation session
- **#778**: Share the oracle-harness scaffolding the two level-2 runners duplicate (Oracle, NotConfigured, _git, _round, sha256)

## PRs Awaiting Review

PRs waiting on Judge (`loom:review-requested`).

- **#770**: refactor(cli): publish exposure outputs through the shared OutputSet
- **#776**: build: run the pinned clippy from package.json and guard the pin
- **#783**: docs: repoint cstool inelastic_axis attribution at assemble.rs

## Approved (Awaiting Merge)

PRs that passed review and are queued for Champion auto-merge (`loom:pr`).

- **#650**: feat: add stopping-table subcommand for stopping and CSDA range tables
- **#759**: refactor(input): share one data-file loader and record the SRIM file's hash
- **#760**: refactor(electron): move the beam model out of transport.rs into electron::beam
- **#779**: feat(wasm): accept [[crystal]] targets in the simulation session
- **#780**: refactor(electron): share one secondary-spawn tail and one Kieft-Bosch options value
- **#784**: refactor(validation): share the level-2 oracle scaffolding through validation/lib

## Proposed

Issues carrying `loom:curated`.

- **#15**: [Epic #10] Verify element and constant data against primary sources *(curated)*
- **#45**: Verify Lenz-Jensen, Moliere magic and screening-length coefficients against primary sources (follow-up to #40) *(curated)*
- **#76**: Benchmarks: quiet-machine re-measure of oracle ions/s and a 1..N thread-scaling curve *(curated)*
- **#111**: [Epic #11] Mott elastic: WKB / asymptotic phase shifts for high-energy partial waves (follow-up to #17) *(curated)*
- **#151**: [Epic #11] Validation: PSF fit vs a published resist PSF measurement *(curated)*
- **#162**: Source the PMMA optical ELF (Ritsko 1978) via a legitimate copy; ask operator about CXRO tables *(curated)*
- **#177**: [Epic #11] Benchmarks follow-up: quiet-machine thread scaling and Nebula CPU timings (remainder of #152) *(curated)*
- **#199**: CI: add a RustSec advisories gate to cargo-deny *(curated)*
- **#360**: Measure the build time and peak RSS of a penn-full table with band-bottom rate refinement (operator-run) *(curated)*
- **#497**: [Epic #33] Digitize heavy-ion aggregate stopping plots (STOP<Z>xx / STOPxx<Z>) with cross-view dedup *(curated)*
- **#511**: [Epic #22] Single-element Ar → Cu and Ar → Ni yields: store a Ni dataset and explain the 0.67 Cu ratio without tuning *(curated)*
- **#515**: [Epic #12] Step 24: low-dose SIMS validation of channeled B, P, As in ⟨100⟩ Si vs tilt/twist and a screen-oxide case *(curated)*
- **#546**: Electron engine: deposited-charge tally (where electrons stop, net charge per layer and voxel) *(curated)*
- **#557**: Electron depth-dose and interface-enhancement validation against published measurements *(curated)*
- **#561**: [Epic #22] Focused ion beams, part 3: lateral implant and damage profiles, film cut-through fluence, trench comparison *(curated)*
- **#562**: [Photon channel #548] A: photon transport with cited attenuation data and per-layer photon dose *(curated)*
- **#565**: [Epic #22] Focused ion beams, part 2: validate sputter yield against incidence angle for one ion-target pair *(curated)*
- **#630**: Epic: lindhard as a comprehensive SRIM/TRIM replacement (feature-parity census and plan) *(curated)*
- **#633**: [Epic #630] Enter the SRIM book's stopping coefficient and effective-charge tables (double entry) *(curated)*
- **#634**: [Epic #630] Compound stopping: cores-and-bonds and phase corrections and a compound dictionary from the book *(curated)*
- **#637**: [Epic #630] `lindhard stopping-table`: stopping and range tables over an energy grid *(curated)*
- **#638**: [Epic #630] Energy-loss straggling beyond Bohr (verify transport use, add Chu-type model) *(curated)*
- **#642**: [Epic #630] TRIM.IN importer *(curated)*
- **#643**: [Epic #630] SRIM-format output writers *(curated)*
- **#658**: Check the low-energy form of the SRIM-data-sourced stopping against the book before entering coefficients *(curated)*
- **#763**: Root package.json clippy scripts run floating stable, not the pinned CI gate toolchain *(curated)*
- **#765**: THIRD_PARTY_LICENSES.md still points cstool's inelastic_axis attribution at lindhard-cli/src/electron.rs *(curated)*
- **#768**: [Epic #25] TOML input for a laterally varying ion target built from CSG solids *(curated)*
- **#769**: [Epic #28] Accept `[[crystal]]` targets in the WASM simulation session *(curated)*

## Proposed (Architect / Hermit)

- **#771**: Consolidate the three copies of Gaussian elimination with partial pivoting (and their divergent singularity tests) *(hermit)*
- **#777**: Trim never-called public API before 0.1.0 locks the surface: 8 items with zero references *(hermit)*
- **#787**: Share one thread-pool helper across the CLI and Python crates: ten copies of the threads == Some(0) guard, two of them unreachable *(hermit)*
- **#788**: Share the one histogram-CSV writer and the one element-symbol lookup: hist_rows exists twice, symbol() three times *(hermit)*
- **#789**: Remove the orphaned one-shot WASM range calculator: calculate()/range_profile lost its only caller when the page folded into the simulator *(hermit)*

## Epics

- **#11**: Epic: M1, low-energy electron engine
- **#12**: Epic: M2, crystalline implant (channeling, damage accumulation, amorphization)
- **#22**: Epic: M3, dynamic composition (fluence-dependent targets)
- **#25**: Epic: 2D/3D target geometry and layout cross-sections
- **#28**: Epic: Python and WASM bindings
- **#548**: Photon channel for electron runs: bremsstrahlung and characteristic x-rays, and the dose they carry beyond the electron range

## Backlog Balance

| Bucket | Count |
|------|-------|
| Operator merge-risk holds | 1 |
| Very important | 1 |
| Important | 5 |
| Ready (`loom:issue`) | 3 |
| In Progress (`loom:building`) | 8 |
| PRs awaiting review | 3 |
| Approved PRs awaiting merge | 6 |
| Curated | 29 |
| Architect / Hermit proposals | 5 |
| Active epics | 6 |
<!-- guide:plan-body:end -->

## Done

M0: #1 materials (PR #14), #2 deterministic RNG streams and parallel driver
(PR #43), #3 screened potentials and scattering integral (PR #40), #4
electronic stopping (PR #42), #5 BCA engine with cascades (PR #46), #6
tallies (PR #48), #7 CLI (PRs #47, #54), #41 density effect, with the
declined Barkas, shell and straggling terms documented (PR #49), #44
Biersack-Varelas join documented, not implemented (PR #60), #50
level-2 oracle adapters and summaries (PR #59), #55 user-supplied stopping
tables (PR #56), #64 opt-in TRIDYN weak collisions and a cascade
electronic-share check against the LNST partition (PR #67).

Since then, among others: M1 elastic validation (#93), hexagonal lattices
(#179), M2 steps 22 and 23 (#513, #514), M3 dynamic-target validations
(#436, #437, #448, #465, #512), sparse 3D voxel tallies and the masked-edge
validation (#433, #435), the WASM build and simulator session (#432, #479),
the GPU feasibility decision (#425, #464), the stopping compilation and
residual report (#430, #431, #485, #495, #496, #498, #499), the clippy pin
and canary (#336), hostile-input tests (#192), and the first SRIM-book
verification (#537).

Subsequent merges and closed issues are recorded in [WORK_LOG.md](WORK_LOG.md).
