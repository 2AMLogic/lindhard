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

## Urgent

*No urgent issues.*

## Ready

Human-approved issues ready for implementation (`loom:issue`): the M0 set,
seeded at bootstrap (epic #10).

- #1 Materials, units and physical constants
- #2 Deterministic per-particle RNG streams + parallel driver
- #3 Screened interatomic potentials + scattering integral
- #4 Electronic stopping models + Bragg additivity

Blocked on those: #5 BCA engine → #6 tallies, #7 CLI, #8 validation, #9 benchmarks.

## In Progress

*No issues currently being built.*

## Proposed

*No proposed issues.*

## Epics

- #10 M0: amorphous ion core, validated
- #11 M1: low-energy electron engine
- #12 M2: crystalline implant
