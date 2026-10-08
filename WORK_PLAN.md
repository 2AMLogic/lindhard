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
stopping-data library, and performance work. Order: M0 first. M1 Phase 1
can start alongside it (it needs only materials and the RNG). M2 follows M0's
engine: M2 Phase 1 may run alongside M0's remaining data verification and
validation datasets (#15, #51), which do not touch the code it builds on
(operator ruling, 2026-10-08, on #12).
The rest run as capacity allows.

<!-- guide:plan-body:start -->
## Operator Attention: Merge-Risk-Hold Pileup

Judge-approved PRs stuck under a `loom:operator` merge-risk hold — implementation work is done, only a human merge decision is missing.

- **#198**: Crystal step 19b: hexagonal lattices (wurtzite GaN, 4H/6H-SiC) (#179)

## Operator Priority

Issues the operator starred (`loom:operator-priority`); land these first.

_None._

## Ready

Human-approved issues ready for implementation (`loom:issue`).

- **#58**: CLI: an electronic-loss-off stopping choice, for like-for-like OpenTRIM comparisons (follow-up to #50)
- **#151**: [Epic #11] Validation: PSF fit vs a published resist PSF measurement
- **#179**: [Epic #12] Crystal step 19b: hexagonal lattices (wurtzite GaN, 4H/6H-SiC)

## In Progress

Issues currently being built (`loom:building`).

_None._

## PRs Awaiting Review

PRs waiting on Judge (`loom:review-requested`).

_None._

## Approved (Awaiting Merge)

PRs that passed review and are queued for Champion auto-merge (`loom:pr`).

- **#198**: Crystal step 19b: hexagonal lattices (wurtzite GaN, 4H/6H-SiC) (#179)

## Proposed

Issues carrying `loom:curated`.

- **#15**: [Epic #10] Verify element and constant data against primary sources *(curated)*
- **#45**: Verify Lenz-Jensen, Moliere magic and screening-length coefficients against primary sources (follow-up to #40) *(curated)*
- **#51**: Validation level 3: first published range datasets (B, P, As in amorphous Si) and stopping-input attribution (follow-up to #8) *(curated)*
- **#58**: CLI: an electronic-loss-off stopping choice, for like-for-like OpenTRIM comparisons (follow-up to #50) *(curated)*
- **#76**: Benchmarks: quiet-machine re-measure of oracle ions/s and a 1..N thread-scaling curve *(curated)*
- **#93**: [Epic #11] Elastic validation: total and transport cross sections vs published partial-wave values (C, Si, Cu, Au) *(curated)*
- **#148**: [Epic #11] Validation: backscatter coefficient η(E, Z) vs published measurements *(curated)*
- **#149**: [Epic #11] Validation: secondary-electron yield δ(E) vs published measurements *(curated)*
- **#151**: [Epic #11] Validation: PSF fit vs a published resist PSF measurement *(curated)*
- **#179**: [Epic #12] Crystal step 19b: hexagonal lattices (wurtzite GaN, 4H/6H-SiC) *(curated)*

## Proposed (Architect / Hermit)

- **#192**: Test: malformed and hostile input values must error, never panic or hang *(architect)*
- **#199**: CI: add a RustSec advisories gate to cargo-deny *(architect)*
- **#226**: Crystal flight: impact-parameter-dependent (Oen-Robinson) electronic loss for channeled ions (Epic #12, Phase 2) *(architect)*

## Epics

- **#10**: Epic: M0, amorphous ion core, validated
- **#11**: Epic: M1, low-energy electron engine
- **#12**: Epic: M2, crystalline implant (channeling, damage accumulation, amorphization)
- **#22**: Epic: M3, dynamic composition (fluence-dependent targets)
- **#25**: Epic: 2D/3D target geometry and layout cross-sections
- **#28**: Epic: Python and WASM bindings
- **#30**: Epic: Release engineering and documentation
- **#33**: Epic: Permissive stopping-power data library
- **#36**: Epic: Performance (profile-guided, SIMD, GPU feasibility)

## Backlog Balance

| Tier | Count |
|------|-------|
| Operator merge-risk holds | 1 |
| Operator priority | 0 |
| Ready (`loom:issue`) | 3 |
| In Progress (`loom:building`) | 0 |
| PRs awaiting review | 0 |
| Approved PRs awaiting merge | 1 |
| Curated | 10 |
| Architect / Hermit proposals | 3 |
| Active epics | 9 |
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

Subsequent merges and closed issues are recorded in [WORK_LOG.md](WORK_LOG.md).
