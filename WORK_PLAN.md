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

Operator decision, 2026-11-09: cut a first real release, 0.1.0, with the
scope below. The release waits almost entirely on operator inputs (a source,
a measurement, a ruling) rather than on agent work, so agent capacity went to
new work from 2026-11-09: 27 issues (#546 to #572) on electron dose and
exposure, electrons on voxel targets, damage outputs, focused ion beams,
optical data and a photon channel. That work is listed under "Alongside the
close-out" below. **It does not change the 0.1.0 scope and does not gate the
release** (operator ruling, 2026-11-10): whatever has merged by the tag ships
with it, and its open validations are listed as such in
`docs/validation.md`.

Operator rulings, 2026-11-09:

- **Scope: hold 0.1.0 until the M1 and M3 validation gates close** (#148,
  #151, #511). They are release blockers, not known deviations.
- **Sources: the operator will obtain** the Oen and Robinson (1976) paper
  (#251) and the SIMS papers for #515, and scan the relevant pages, as for
  the SRIM book. Scans stay out of the tree.
- **Publishing: GitHub release, crates.io and PyPI.** The website stays
  live but unannounced, so the browser check (#469) does not gate the
  release.

Where the milestones stand:

| Milestone | State | What is left |
|---|---|---|
| M0 amorphous ion core (#10) | **Closed 2026-11-10.** Feature-complete and validated (levels 1 to 3) | Recorded gaps, not release blockers: #15 (ASTM E521 `E_d`; C and Mn densities), #45 (Lenz-Jensen polynomial digits; Kr-C screening coefficients) |
| M1 electron engine (#11) | Engine done, 24 sub-issues closed | Validation gates #148 (η), #151 (PSF); #111 and #162 wait on sources |
| M2 crystalline implant (#12) | Engine done through step 23 (#513, #514) | SIMS validation #515 then #516; #251 (Oen-Robinson length) affects channeled-range accuracy |
| M3 dynamic composition (#22) | Engine and two validations done | #511 (Ar to Cu/Ni yield ratio): #543, then #544. Focused ion beams (#554) run alongside |
| 2D/3D geometry (#25) | Voxel tallies and masked-edge validation done | #434 waits on a `klayout-tools` polygon export. Electrons on voxel targets (#549) run alongside |
| Bindings (#28) | Python and WASM done | #469 browser check (operator) |
| Stopping-data library (#33) | Loader, residual report, compilation bibliography done | Digitization #486, #497 and seed datasets #483, #503: **after 0.1.0** |

### Close-out steps, in order

1. ~~**Land the SRIM-book verifications.**~~ Done 2026-11-09: #540 and
   #538 (0.8853 universal length) merged.
2. ~~**Narrow #45 and #15 to what is truly inaccessible,** then close M0
   (#10).~~ Done 2026-11-10: both issues now open with their remaining
   scope, and #10 is closed. One of its success criteria (every provenance
   row verified against its source) is not met, and the close says so. Left
   in #45: the Lenz-Jensen polynomial digits and the Kr-C screening
   coefficients (Wilson, Haggmark and Biersack 1977). Left in #15: the ASTM
   E521 `E_d` values and the C and Mn densities. These are recorded gaps,
   not release blockers.
3. **M2 SIMS validation** (#515, then #516). Before it, close #251
   (Oen-Robinson length), since the channeled tail depends on it. Both wait
   on the papers the operator is obtaining.
4. **Close the M1 and M3 validation gates** (operator ruling: they block
   0.1.0):
   - #148 (backscatter η): parked by the operator's 2026-10-19 ruling (the
     0.05 tolerance stands; Si misses at 20 and 30 keV stay visible). #457
     delivered a source audit and comparison plan (`docs/validation.md`,
     "Si shell inputs, the 199 eV join and Au"). Executing it needs operator
     decisions: whether a source for the Si core region above 199 eV (the
     Henke tables, unopened) is admissible, and authorizing the comparison
     runs. Operator inputs.
   - #511 (Ar to Cu/Ni yield ratio): split on 2026-11-09 as the Builder
     recommended. #543 stores the Ni datasets (IPPJ-AM-32, PDF p. 103,
     digitized by the existing method); #544 then attributes the 0.67 Cu
     ratio without tuning. Agent work. A first Builder pass on #543
     (2026-11-09) stopped on glyph stacks at 60 to 150 eV, which lie below
     the stored range (196 eV to 10.2 keV); the issue was returned to the
     queue on 2026-11-10 with that noted.
   - #151 (PSF): needs an openable measured PSF for a stated resist stack
     and energy, and the PMMA optical ELF (#162; the open question there is
     whether the CXRO tables are acceptable). Operator inputs.
5. **Cut 0.1.0** (`CONTRIBUTING.md`, "Releasing"): CHANGELOG section,
   workspace version, annotated tag; publish the GitHub release, the crate
   to crates.io and the wheels to PyPI. The website is not announced with
   it, so #469 is not a prerequisite.

### Alongside the close-out

Opened 2026-11-09 and 2026-11-10. None of it gates 0.1.0. Each front carries
its own validation issue, and a front is not done until that closes.

| Front | Issues | Open validation or gap |
|---|---|---|
| Electron dose, charge and exposure | #546, #547, #550, #551 | #557 (depth-dose against published measurements), #555 (direct Monte Carlo check of the PSF convolution); #546 is delivered as a tally (#578) and still wants a comparison |
| Electrons on voxel targets (#549, Epic #25) | #570, #571 | #572 (position-dependent PSF and buried-edge validation) |
| Damage outputs (#553) | #558, #560 | #563 (versioned export with a Python round trip) |
| Focused ion beams (#554, Epic #22) | #559 | #565 (yield against incidence angle; needs a readable paper), #561 |
| Optical energy-loss functions (#552) | #556, #564, #568, #569 | Only amorphous Al2O3 has a sourced ELF so far (#577). W, Ta, Ti, Co, Ru, TiN, TaN and HSQ are recorded in `docs/data-provenance.md` as gaps: no openable source |
| Photon channel (#548) | #562, then #566, #567 | #562 waits on an operator ruling on the attenuation data source |

### Operator inputs that unblock work

| Input | Unblocks |
|---|---|
| Oen and Robinson, NIM 132, 647 (1976) (operator obtaining) | #251, then M2 accuracy |
| Lenz (1932) and Jensen (1932); Wilson, Haggmark, Biersack (1977) | the rest of #45 |
| ASTM E521 (any edition); a CRC Handbook edition later than the 80th | the rest of #15 |
| A ruling on the photon attenuation data source (NIST XCOM's terms reserve all rights; EPDL97/EPICS2017 is the proposed alternative) | #562, then #566, #567 |
| Oliva-Florio et al., Phys. Rev. B 35, 2198 (1987), or another measured yield against incidence angle | #565 |
| Openable optical data for W, Ta, Ti, Co, Ru, TiN, TaN and HSQ | the rest of #556, #564, #569 |
| SIMS profiles from papers (B, P, As in Si; screen oxide) (operator obtaining) | #515, #516 |
| An openable measured resist PSF; a ruling on the CXRO tables for the PMMA ELF | #151 (release gate), #162 |
| A ruling on a Si core-region ELF source (Henke tables) and authorization of the #457 comparison runs | #148 (release gate) |
| Primary stopping papers for B, P, As in Si | #503, #483 |
| A sourced high-l phase-shift method | #111 |
| An idle machine for benchmarks; a penn-full build measurement | #76, #177; #360 |
| A real-browser check of the WASM simulator (needed only to announce the website) | #469 |
| Loom guard reviews | #195, #220 |

### After 0.1.0

The stopping-data library (#33: #486, #497, #483, #487), the public website
(#467), layout import (#434, once `klayout-tools` exports polygons), the
Kr-C magic-formula set (Table 4-1 of the SRIM book, p. 7-6; checked against
the quadrature like the other two), a RustSec advisory report (#199) and the
docs follow-up #121.

<!-- guide:plan-body:start -->
## Operator Attention: Merge-Risk-Hold Pileup

Judge-approved PRs under a hold (`loom:blocked` or `loom:operator`) — only a human merge decision is missing.

- **#650**: feat: add stopping-table subcommand for stopping and CSDA range tables

## Very Important

`loom:very-important` (one counts per workspace); every stage takes these first.

_None._

## Important

`loom:important`; taken before unlabelled work.

- **#486**: [Epic #33] Digitize the measured points from the srim.org stopping plots (cited, flagged, count-checked)
- **#497**: [Epic #33] Digitize heavy-ion aggregate stopping plots (STOP<Z>xx / STOPxx<Z>) with cross-view dedup
- **#637**: [Epic #630] `lindhard stopping-table`: stopping and range tables over an energy grid

## Ready

Human-approved issues ready for implementation (`loom:issue`).

- **#562**: [Photon channel #548] A: photon transport with cited attenuation data and per-layer photon dose

## In Progress

Issues currently being built (`loom:building`).

- **#730**: Use common::with_threads for the 48 hand-built rayon pools in lindhard/tests

## PRs Awaiting Review

PRs waiting on Judge (`loom:review-requested`).

- **#731**: refactor(tests): call common::with_threads for the hand-built rayon pools

## Approved (Awaiting Merge)

PRs that passed review and are queued for Champion auto-merge (`loom:pr`).

- **#650**: feat: add stopping-table subcommand for stopping and CSDA range tables

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

## Proposed (Architect / Hermit)

_None._

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
| Very important | 0 |
| Important | 3 |
| Ready (`loom:issue`) | 1 |
| In Progress (`loom:building`) | 1 |
| PRs awaiting review | 1 |
| Approved PRs awaiting merge | 1 |
| Curated | 25 |
| Architect / Hermit proposals | 0 |
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
