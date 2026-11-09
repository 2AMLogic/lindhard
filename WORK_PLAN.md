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

## Current phase: close out and release 0.1.0

Operator decision, 2026-11-09: finish what is started and cut a first real
release before opening new fronts. Since the last refresh (2026-10-15), 54
issues closed. Merges have nearly stopped since 2026-10-31, because almost
every open issue now waits on an operator input (a source, a measurement, a
quiet machine, a browser check) rather than on agent work. The SRIM book
(Tier C literature, `CONTRIBUTING.md`) unblocked several of the source gaps
on 2026-11-08 and 2026-11-09 (#537, #540).

Where the milestones stand:

| Milestone | State | What is left |
|---|---|---|
| M0 amorphous ion core (#10) | Feature-complete and validated (levels 1 to 3) | #15 (constant-data verification: only inaccessible sources remain), #45 (screening and magic constants: narrowed by the SRIM book), #538 (0.8853 universal length) |
| M1 electron engine (#11) | Engine done, 24 sub-issues closed | Validation gates #148 (η), #151 (PSF); #111 and #162 wait on sources |
| M2 crystalline implant (#12) | Engine done through step 23 (#513, #514) | SIMS validation #515 then #516; #251 (Oen-Robinson length) affects channeled-range accuracy |
| M3 dynamic composition (#22) | Engine and two validations done | #511 (Ar to Cu/Ni yield ratio) |
| 2D/3D geometry (#25) | Voxel tallies and masked-edge validation done | #434 waits on a `klayout-tools` polygon export |
| Bindings (#28) | Python and WASM done | #469 browser check (operator) |
| Stopping-data library (#33) | Loader, residual report, compilation bibliography done | Digitization #486, #497 and seed datasets #483, #503: **after 0.1.0** |

### Close-out steps, in order

1. **Land the SRIM-book verifications.** #540 (in review); #538 needs only
   `EXPECTED_CRYSTAL_AARCH64_MACOS` regenerated on aarch64 macOS (the
   attended operator session can do it).
2. **Narrow #45 and #15 to what is truly inaccessible,** then close M0
   (#10). After #537 and #540, the ZBL and Moliere sets, both magic-formula
   sets and the nuclear stopping fit are verified. Left: the Lenz-Jensen
   polynomial digits and the Kr-C coefficients (Wilson, Haggmark and
   Biersack 1977). These are recorded gaps, not release blockers.
3. **M2 SIMS validation** (#515, then #516). Before it, try to close #251
   (Oen-Robinson length), since the channeled tail depends on it.
4. **Decide the 0.1.0 scope for M1 and M3:** ship them with their open
   validation gates (#148, #151, #511) listed as known deviations in the
   README, or hold the release for them. Default: ship and list them.
5. **Cut 0.1.0** (`CONTRIBUTING.md`, "Releasing"): CHANGELOG section,
   workspace version, annotated tag. Before tagging: the browser check of
   the WASM page (#469) if the site is announced with the release.

### Operator inputs that unblock work

| Input | Unblocks |
|---|---|
| macOS golden regeneration | #538 |
| Oen and Robinson, NIM 132, 647 (1976) | #251, then M2 accuracy |
| Lenz (1932) and Jensen (1932); Wilson, Haggmark, Biersack (1977) | the rest of #45 |
| SIMS profiles from openable papers (B, P, As in Si; screen oxide) | #515, #516 |
| Primary stopping papers for B, P, As in Si | #503, #483 |
| PMMA ELF (Ritsko 1978); a sourced high-l phase-shift method | #162; #111 |
| An idle machine for benchmarks; a penn-full build measurement | #76, #177; #360 |
| A real-browser check of the WASM simulator | #469 |
| Loom guard reviews | #195, #220 |

### After 0.1.0

The stopping-data library (#33: #486, #497, #483, #487), the public website
(#467), layout import (#434, once `klayout-tools` exports polygons), the
Kr-C magic-formula set (Table 4-1 of the SRIM book, p. 7-6; checked against
the quadrature like the other two), a RustSec advisory report (#199) and the
docs follow-up #121.

<!-- guide:plan-body:start -->
## Operator Attention: Merge-Risk-Hold Pileup

Judge-approved PRs stuck under a `loom:operator` merge-risk hold — implementation work is done, only a human merge decision is missing.

- **PR #538**: fix: use the Thomas-Fermi constant for the universal screening length *(needs the aarch64 macOS crystal golden set)*

## Operator Priority

Issues the operator starred (`loom:operator-priority`); land these first.

_None._

## Ready

Human-approved issues ready for implementation (`loom:issue`).

_None._

## In Progress

Issues currently being built (`loom:building`).

- **#536**: Universal screening length: use the 0.8853 prefactor printed in the SRIM book's Eq. 2-73 *(PR #538)*

## PRs Awaiting Review

PRs waiting on Judge (`loom:review-requested`).

- **PR #540**: docs: verify Molière constants and TRIM-chapter formulas against the SRIM book
- **PR #421**: docs: Guide document maintenance update *(2026-10-15; superseded by this refresh)*

## Approved (Awaiting Merge)

PRs that passed review and are queued for Champion auto-merge (`loom:pr`).

- **PR #539**: chore(deps): bump rizzma from 1.13.7 to 1.14.1 *(also `loom:blocked`)*

## Proposed

Issues carrying `loom:curated`.

- **#15**: [Epic #10] Verify element and constant data against primary sources *(curated, blocked)*
- **#45**: Verify Lenz-Jensen, Moliere magic and screening-length coefficients against primary sources *(curated, blocked; narrowed by #537, #540)*
- **#76**: Benchmarks: quiet-machine re-measure of oracle ions/s and a 1..N thread-scaling curve *(curated, blocked)*
- **#111**: [Epic #11] Mott elastic: WKB / asymptotic phase shifts for high-energy partial waves *(curated, blocked)*
- **#148**: [Epic #11] Validation: backscatter coefficient η(E, Z) vs published measurements *(curated, blocked)*
- **#151**: [Epic #11] Validation: PSF fit vs a published resist PSF measurement *(curated, blocked)*
- **#162**: Source the PMMA optical ELF (Ritsko 1978) via a legitimate copy *(curated, blocked)*
- **#177**: [Epic #11] Benchmarks follow-up: quiet-machine thread scaling and Nebula CPU timings *(curated, operator)*
- **#199**: CI: add a RustSec advisories gate to cargo-deny *(curated, blocked)*
- **#360**: Measure the build time and peak RSS of a penn-full table *(curated, blocked)*
- **#483**: [Epic #33] Seed measured stopping datasets from primary papers *(curated, blocked)*
- **#497**: [Epic #33] Digitize heavy-ion aggregate stopping plots *(curated, operator, important)*
- **#511**: [Epic #22] Ar → Cu and Ar → Ni yields: store a Ni dataset and explain the 0.67 Cu ratio *(curated, blocked)*
- **#515**: [Epic #12] Step 24: low-dose SIMS validation of channeled B, P, As in ⟨100⟩ Si *(curated, operator)*

## Proposed (Architect / Hermit)

_None._

## Epics

- **#10**: Epic: M0, amorphous ion core, validated
- **#11**: Epic: M1, low-energy electron engine
- **#12**: Epic: M2, crystalline implant (channeling, damage accumulation, amorphization)
- **#22**: Epic: M3, dynamic composition (fluence-dependent targets)
- **#25**: Epic: 2D/3D target geometry and layout cross-sections
- **#28**: Epic: Python and WASM bindings
- **#33**: Epic: Permissive stopping-power data library

## Backlog Balance

| Tier | Count |
|------|-------|
| Operator merge-risk holds | 1 |
| Operator priority | 0 |
| Ready (`loom:issue`) | 0 |
| In Progress (`loom:building`) | 1 |
| PRs awaiting review | 2 |
| Approved PRs awaiting merge | 1 |
| Curated | 14 |
| Architect / Hermit proposals | 0 |
| Active epics | 7 |
| Open issues, all | 34 |
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
