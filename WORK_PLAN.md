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
can start alongside it (it needs only materials and the RNG). M2 follows M0.
The rest run as capacity allows.

## Urgent

*No urgent issues.*

## Ready

Human-approved issues ready for implementation (`loom:issue`): the M0 set,
seeded at bootstrap (epic #10).

- #2 Deterministic per-particle RNG streams + parallel driver
- #3 Screened interatomic potentials + scattering integral
- #4 Electronic stopping models + Bragg additivity

Blocked on those: #5 BCA engine → #6 tallies, #7 CLI, #8 validation, #9 benchmarks.

## In Progress

*No issues currently being built.*

## Proposed

Phase 1 issues awaiting Champion approval (`loom:architect` + `loom:epic-phase`):
#15–#21, #23–#24, #26–#27, #29, #31–#32, #34–#35, #37. Later phases are
filed by Champion as each phase completes.

## Epics

| Epic | Scope | Phase 1 |
|---|---|---|
| #10 | M0: amorphous ion core, validated | #2–#4, #15 (#1 merged) |
| #11 | M1: low-energy electron engine | #16–#18 |
| #12 | M2: crystalline implant | #19–#21 |
| #22 | M3: dynamic composition | #23–#24 |
| #25 | 2D/3D geometry and layout cross-sections | #26–#27 |
| #28 | Python and WASM bindings | #29 |
| #30 | Release engineering and documentation | #31–#32 |
| #33 | Permissive stopping-data library (terms review needs an operator ruling) | #34–#35 |
| #36 | Performance (profile-guided, SIMD, GPU feasibility) | #37 |
