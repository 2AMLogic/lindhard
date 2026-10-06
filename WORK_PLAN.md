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

Human-approved issues ready for implementation (`loom:issue`):

*None at present.* The remaining M0 items below are curated but blocked or
partly done, or are awaiting triage.

## In Progress

*No issues currently being built.*

## Proposed

Phase 1 issues awaiting Champion approval (`loom:architect` + `loom:epic-phase`):
#15–#21, #23–#24, #26–#27, #29, #31–#32, #34–#35, #37. Later phases are
filed by Champion as each phase completes.

M0 follow-ups awaiting triage (`loom:triage`):

- #58 CLI: an electronic-loss-off stopping choice, for like-for-like
  OpenTRIM comparisons (follow-up to #50)
- #61 Investigate low sputter yields: Ar 1 keV → Cu is 2.2x below RustBCA
  even with `E_d = E_s` (follow-up to #50). Investigated: the gap is
  electronic loss on low-energy recoils; the fix is #64
- #64 BCA: nuclear loss truncated at p_max while electronic loss is not,
  over-damping low-energy recoils (follow-up to #61)

Open M0 work (`loom:curated`), partly done:

- #8 Validation harness: level 1 in CI (PR #52) and level-2 oracles (#50,
  PR #59) are done; level 3 continues under #51.
- #9 Benchmarks: hot-path benches and thread scaling are done (PR #53), and
  oracle ions/s are recorded by #50; a quiet-machine re-measure and a
  profile of the hot spots remain.
- #45 Verify Lenz-Jensen, Moliere magic and screening-length coefficients
  against primary sources. PR #57 verified what open-access sources allow;
  the rest needs closed-access primary papers, so it is not in the ready
  queue.

## Done

M0: #1 materials (PR #14), #2 deterministic RNG streams and parallel driver
(PR #43), #3 screened potentials and scattering integral (PR #40), #4
electronic stopping (PR #42), #5 BCA engine with cascades (PR #46), #6
tallies (PR #48), #7 CLI (PRs #47, #54), #41 density effect, with the
declined Barkas, shell and straggling terms documented (PR #49), #44
Biersack-Varelas join documented, not implemented (PR #60), #50
level-2 oracle adapters and summaries (PR #59), #55 user-supplied stopping
tables (PR #56).

## Epics

| Epic | Scope | Phase 1 |
|---|---|---|
| #10 | M0: amorphous ion core, validated | Done: #1–#7, #41, #44, #50, #55. Partial: #8, #9. Open: #15, #45, #51, #58, #61 |
| #11 | M1: low-energy electron engine | #16–#18 |
| #12 | M2: crystalline implant | #19–#21 |
| #22 | M3: dynamic composition | #23–#24 |
| #25 | 2D/3D geometry and layout cross-sections | #26–#27 |
| #28 | Python and WASM bindings | #29 |
| #30 | Release engineering and documentation | #31–#32 |
| #33 | Permissive stopping-data library (terms review needs an operator ruling) | #34–#35 |
| #36 | Performance (profile-guided, SIMD, GPU feasibility) | #37 |
