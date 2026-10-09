# Work Log

Chronological record of notable decisions and merges.

### 2026-10-09

- **PR #286**: Table cache: verify cached table bytes and key the elastic model choice
- **PR #281**: docs: correct the PSF summary and CSV description in cli.md
- **PR #283**: Derive Debug on PsfOutcome (fix red CI on main)
- **PR #278**: Enable unreachable_pub and missing_debug_implementations in [workspace.lints]
- **PR #275**: Crystal: validate the off-axis thermal tail and rescope the 7/22 criterion (#225)
- **PR #274**: feat: CSG primitives and boolean combinations on the Geometry trait
- **Issue #282** (closed): CI red on main: PsfOutcome lacks Debug under missing_debug_implementations
- **Issue #271** (closed): CLI: expose cubic crystal targets and orientation through shared TOML input
- **Issue #268** (closed): Electron diagnostics: count capped secondary tracks and affected histories
- **Issue #265** (closed): CI: run the validation Python unit tests (three test files never run on PRs)
- **Issue #264** (closed): Refactor: split lindhard/src/geometry.rs into stack and voxel modules (no behaviour change)
- **Issue #260** (closed): Refactor: split lindhard/src/input.rs into schema and resolve modules (no behaviour change)
- **Issue #259** (closed): CI: scheduled release-mode run of the #[ignore]d statistical and slow tests
- **Issue #254** (closed): Correct the PSF summary and CSV documentation
- **Issue #253** (closed): Electron validation: report actual cross-section table endpoint use
- **Issue #249** (closed): Table cache: verify cached table bytes and key the elastic model choice (follow-up to #248)
- **Issue #225** (closed): Crystal: off-axis channeling tail grows with thermal vibration (7°/22°, 30°/17°): validate and rescope the 7/22 dRp criterion
- **Issue #222** (closed): Enable unreachable_pub and missing_debug_implementations in [workspace.lints]
- **Issue #194** (closed): [Epic #25] CSG primitives and boolean combinations on the Geometry trait
- **Issue #169** (closed): Backscatter validation: evaluate C, Si, Au and rerun with full Penn + DHFS (follow-up to #148)
- **Issue #131** (closed): Elastic solver: start-radius rule r_t exp(-60/|kappa|) is too close to the turning point for |kappa| of several hundred
- **Issue #123** (closed): tally::electron: state (or change) the energy reference of deposits when secondaries are on
- **Issue #81** (closed): docs: clean-room warning — cited IPP 9/64 PDF contains the TRIDYN program listing (Appendix 1)
- **PR #276**: CLI: expose cubic crystal targets and orientation through shared TOML input
- **PR #272**: Backscatter validation: DHFS for all five targets, add Si (#169)
- **PR #270**: Electron diagnostics: count capped secondary tracks and affected histories
- **PR #269**: Split geometry.rs into stack and voxel modules (no behaviour change)
- **PR #267**: Run the validation Python unit tests in CI
- **PR #266**: docs: warn that IPP 9/64 bundles the TRIDYN listing (Tier C)
- **PR #263**: Split input.rs into schema, resolve and resolved modules (no behaviour change)
- **PR #262**: CI: scheduled release-mode run of the #[ignore]d statistical and slow tests
- **PR #261**: Use the WKB start radius for every elastic potential (#131)
- **PR #258**: Document the band-bottom reference of electron deposits (#123)
- **PR #257**: Report electron cross-section table endpoint use (#253)
- **PR #246**: Build hydrogenic K-shell optical ELFs for the inner-shell channels (#135)
- **PR #252**: test: validate the crystal Oen-Robinson local loss and flag its unverified constants
- **PR #248**: Electron CLI: reuse built cross-section tables across runs (--table-cache)
- **PR #244**: docs(electron): trace the Al/Au single-pole Penn delta overestimate to the low-energy IMFP
- **PR #243**: Cited electron band defaults (Si, SiO2) and a narrowed band-parameter gap (#115)
- **PR #240**: Add Si optical ELF digitized from Yang et al. (2019) (#125)
- **PR #239**: Enter the Salvat et al. (1987) DHFS screening table (Z = 1..92)
- **PR #238**: feat(dynamic): sputter erosion and surface recession (#234)
- **Issue #234** (closed): Dynamic composition: sputter erosion and surface recession in the fluence loop (Epic #22, Phase 2)
- **Issue #226** (closed): Crystal flight: impact-parameter-dependent (Oen-Robinson) electronic loss for channeled ions (Epic #12, Phase 2)
- **Issue #205** (closed): validation: Au secondary-electron yield delta_max overestimate (3.42 vs 1.468)
- **Issue #204** (closed): Au secondary-electron yield with the default (single-pole Penn) model is about 2.3x the measured δ_max
- **Issue #173** (closed): Al secondary-electron yield with the default (single-pole Penn) model is about 5x the measured δ_max
- **Issue #168** (closed): Electron CLI: reuse built cross-section tables across runs, so the full Penn model can enter the #150 oracle comparison
- **Issue #130** (closed): Enter the Salvat et al. (1987) DHFS screening table (Z=1..92): the paper is now reachable in an open repository
- **Issue #127** (closed): Fix clippy nonminimal_bool in ion/scattering.rs:50 (blocks clippy -D warnings)
- **Issue #125** (closed): Source Si valence-region optical ELF data (6 to 30 eV) for #98
- **Issue #115** (closed): Electron band parameters: cited per-material defaults for work function, Fermi energy, affinity and band gap
- **PR #235**: ci: tag-triggered release workflow
- **Issue #184** (closed): Add #![forbid(unsafe_code)] to lindhard-cli and lindhard-py
- **Issue #231** (closed): Guard telemetry disposition: retain scope protection for /dev/shm cleanup
- **Issue #215** (closed): Guard decision review: stash-scope:create-redirect
- **Issue #232** (closed): Release workflow: tag-triggered CLI binaries, checksums and crates.io publish (Epic #30, Phase 2)
- **Issue #214** (closed): Guard decision review: worktree-write-confinement-unresolved-var

### 2026-10-08

- **PR #228**: ci: stubtest the lindhard-py .pyi stubs against the compiled extension
- **Issue #227** (closed): CI: stubtest the lindhard-py .pyi stubs against the compiled extension
- **PR #224**: feat(crystal): thermal vibration in the lattice collision search (#181)
- **PR #223**: Declare lint policy once in [workspace.lints] (#218)
- **Issue #218** (closed): Declare lint policy once in [workspace.lints] and opt all crates in
- **Issue #181** (closed): [Epic #12] Crystal step 21b: thermal vibrations in the lattice search
- **PR #219**: feat: ship es-sputter-ar-v1 opt-in E_s tuning set with fit record (#80)
- **PR #216**: test: pin fixed-seed golden results for cross-platform reproducibility
- **PR #212**: ci: add rustdoc job and fix rustdoc warnings (#191, #119)
- **PR #211**: fix(cli): encode material names in electron table CSV (#207)
- **PR #210**: docs: replace stale crate-doc status paragraph with a stable module map
- **PR #209**: ci: check the declared MSRV in a dedicated msrv job
- **PR #208**: feat(bca): crystal flight model with per-region partner selection
- **PR #206**: fix(cli): remove stale optional CSVs when reusing an output directory
- **PR #203**: validation: Au secondary-electron yield delta(E), bounds and barrier inputs (#149)
- **PR #202**: feat(input): opt-in E_s tuning-set plumbing (#80, plumbing only)
- **PR #197**: feat(geometry): triangle-mesh solids (STL/OBJ) with a BVH
- **PR #196**: docs(provenance): terms record for IAEA and NIST stopping-power data (#34)
- **PR #189**: feat(crystal): unit-cell lattice neighbour search (#20)
- **PR #187**: feat(crystal): cubic lattice model and wafer/beam orientation (#19)
- **Issue #213** (closed): Test: pin fixed-seed golden results so the CI OS matrix checks cross-platform reproducibility
- **Issue #207** (closed): CLI: escape material identifiers in electron table CSV exports
- **Issue #200** (closed): Docs: replace stale 'implemented so far' crate docs in lib.rs with a stable module map
- **Issue #191** (closed): CI: build rustdoc with warnings denied so doc-comment breakage cannot land
- **Issue #190** (closed): CI: check the declared MSRV (rust-version 1.85) in a dedicated job
- **Issue #185** (closed): CLI: remove stale optional CSVs when reusing an output directory
- **Issue #180** (closed): [Epic #12] Crystal step 20b: crystal flight model in the BCA engine
- **Issue #159** (closed): Guard refinement: resolve ephemeral TMPDIR writes in worktree sessions
- **Issue #158** (closed): Guard telemetry disposition: retain stash scope protection for lint command
- **Issue #119** (closed): rustdoc: private intra-doc link in Lenz-Jensen docs fails cargo doc -D warnings
- **Issue #80** (closed): Opt-in experiment-tuned mode: phenomenological correction factors behind an on/off flag
- **Issue #34** (closed): [Epic #33] Data-terms review: IAEA stopping database and NIST SRD (PSTAR/ASTAR/ESTAR)
- **Issue #27** (closed): [Epic #25] Triangle-mesh solids (STL/OBJ) with a BVH
- **Issue #20** (closed): [Epic #12] Crystal step 20a: lattice neighbor search (candidate collision partners along a path)
- **Issue #19** (closed): [Epic #12] Crystal step 19a: cubic lattice model and wafer/beam orientation
- **PR #186**: feat(crystal): Debye thermal-displacement model and sampling (#21)
- **PR #183**: docs: record resist PSF validation gap and sources tried
- **PR #182**: docs(WORK_PLAN): M2 Phase 1 may run alongside M0's remaining verification
- **PR #178**: validation: backscatter eta for glassy C and Au (C, Al, Cu, Au compared)
- **PR #176**: test: harden elastic reference harness (part of #93)
- **PR #175**: validation: secondary-electron yield δ(E) vs published measurements (#149)
- **PR #174**: Electron engine benchmarks: electrons/s, thread scaling, determinism (#152)
- **PR #160**: feat(tally): radial PSF extraction and double/triple-Gaussian fits (#146)
- **PR #155**: Verify E_s against Kittel 8th ed. Table 1 (p. 50); correct Na and K (part of #15)
- **Issue #152** (closed): [Epic #11] Benchmarks: electrons/s and thread scaling vs Nebula (CPU)
- **Issue #146** (closed): [Epic #11] PSF extraction: radial point-spread function and double/triple-Gaussian fit
- **Issue #120** (closed): WORK_PLAN.md: epic #11 children list is stale after the M1 restructure
- **Issue #21** (closed): [Epic #12] Crystal step 21a: Debye thermal-displacement model
- **Issue #8** (closed): Validation harness: analytic checks, oracle runner (RustBCA/OpenTRIM), first experimental range datasets
- **PR #171**: chore: adopt the 2AMLogic Renovate preset

### 2026-10-07

- **PR #170**: Electron oracle comparison vs Nebula and Geant4 MicroElec on matched problems (#150)
- **PR #167**: validation: backscatter coefficient eta(E, Z) vs published measurements (#148)
- **PR #166**: chore: resync installed Loom surfaces to 0.19.874
- **PR #164**: chore: resync installed Repo Skills to 0.21.2
- **PR #163**: chore: repo hygiene fixes (2026-10-07)
- **PR #161**: docs: PMMA optical ELF source could not be opened; gap recorded (#147)
- **PR #154**: feat(cli): [electron] run mode with electron_summary.json and CSV output
- **PR #153**: docs: record that the JSP 2004 reprint has no sigma_el/sigma_tr1 table for C, Si, Cu, Au (Part of #93)
- **PR #144**: Inelastic validation: full Penn IMFPs vs TPP 2011 for Al and Cu (#99)
- **PR #143**: feat(electron): Mermin dielectric function and MELF oscillator fit (#96)
- **PR #142**: Elastic validation harness vs partial-wave values (#93): placeholder, no reference values
- **PR #141**: chore(py): remove unused serde_json dependency from lindhard-py
- **PR #140**: validation: sputter sensitivity to doubtful sets; Fig. 310 tick wording (#78)
- **PR #139**: docs: #61 sputter-yield investigation summary and reproduction recipe
- **PR #138**: Fluence stepping loop for dynamic composition (#24)
- **PR #136**: electron::inelastic::table: energy-loss CrossSectionTable and q sampler (#94)
- **PR #134**: electron: shell-resolved inner-shell channels and Born-Ochkur exchange for the Penn model
- **PR #133**: electron::inelastic: full Penn algorithm over Lindhard dielectric functions
- **PR #132**: feat(electron::elastic): exchange and correlation-polarization corrections (#91)
- **PR #129**: docs: defer the muffin-tin potential, record the gap (#92)
- **PR #128**: data: cross-check Mn density (NBS Monograph 25); record Kittel and E521 access attempts (#15)
- **PR #126**: Electron transport: Frohlich phonon scattering and polaron trapping channels
- **PR #124**: docs: move Salvat DHFS and elastic-solver rows into the provenance inventory table
- **PR #122**: electron::elastic::table: energy-grid elastic CrossSectionTables (#90)
- **PR #118**: feat: single-pole Penn DIIMFP, IMFP and stopping from an OpticalElf
- **PR #117**: loom: enable per-worktree cargo target dirs (#114)
- **PR #116**: electron: secondary electrons, step barrier and interface refraction (#101)
- **PR #113**: tally::electron: 3D deposition grids, eta/delta yields, escape spectra (#103)
- **PR #110**: data: check Cr, Au, Ra, Ac densities against NBS and AEC reports (#15)
- **PR #109**: electron::elastic: radial Dirac partial-wave solver (#17)
- **PR #108**: electron::transport: event loop in layered stacks on CrossSectionTable inputs
- **PR #107**: Commit Al and Cu optical ELF datasets with sum-rule test (#98)
- **PR #106**: Commit EADL2017 subshell binding energies (Z=1..92)
- **PR #105**: data: check element densities against the X-Ray Data Booklet (CRC 80th ed.); record E_d access gaps
- **PR #104**: Voxel-grid geometry: Geometry trait, VoxelGrid and exact 3D DDA crossing
- **PR #88**: perf: collision hot path - sin/cos angle chain, cached stopping coefficient, reused buffers (#37)
- **PR #86**: Stopping-data format, loader and model comparison (#35)
- **PR #85**: docs: add the lindhard book (physics manual and user guide) with CI build and Pages deploy
- **PR #84**: lindhard-py: pyo3 bindings with NumPy outputs
- **PR #83**: Composition-depth grid with volume relaxation (ion::dynamic)
- **PR #82**: Make both crates publish-ready: allowlists, CHANGELOG, CI package job
- **Issue #150** (closed): [Epic #11] Validation: oracle comparison vs Nebula and Geant4 MicroElec on matched problems
- **Issue #147** (closed): [Epic #11] Commit PMMA optical ELF (Ritsko et al. 1978) as a cited measurement
- **Issue #145** (closed): [Epic #11] CLI: [electron] input section and electron report output
- **Issue #137** (closed): Remove unused serde_json dependency from lindhard-py
- **Issue #114** (closed): Loom config: enable per-worktree cargo target dirs (cargo.perWorktreeTargetDir)
- **Issue #112** (closed): docs/data-provenance.md: Salvat DHFS and elastic-solver rows sit under Open questions instead of the inventory table
- **Issue #103** (closed): [Epic #11] Electron transport, step 4: electron tallies — 3D deposition grid, BSE/SE yields and spectra
- **Issue #102** (closed): [Epic #11] Electron transport, step 3: insulator losses — Fröhlich phonon scattering and polaron trapping
- **Issue #101** (closed): [Epic #11] Electron transport, step 2: secondary-electron generation, surface barrier and interface refraction
- **Issue #100** (closed): [Epic #11] Electron transport, step 1: event-by-event loop in layered stacks on CrossSectionTable inputs
- **Issue #99** (closed): [Epic #11] Inelastic validation: IMFPs vs Tanuma–Powell–Penn for Al and Cu
- **Issue #98** (closed): [Epic #11] Commit first optical ELF datasets (Si, Al, Cu) as cited measurements
- **Issue #97** (closed): [Epic #11] Penn inelastic, step 5: inner-shell ionization channel and slow-electron exchange correction
- **Issue #96** (closed): [Epic #11] Penn inelastic, step 4: Mermin-ELF oscillator fit to optical data
- **Issue #95** (closed): [Epic #11] Penn inelastic, step 3: full Penn algorithm (Lindhard-function extension)
- **Issue #94** (closed): [Epic #11] Penn inelastic, step 2: energy-loss and angle sampling tables into CrossSectionTable
- **Issue #92** (closed): [Epic #11] Mott elastic, step 4: muffin-tin potential for condensed targets (or a documented deferral)
- **Issue #91** (closed): [Epic #11] Mott elastic, step 3: exchange and correlation-polarization corrections
- **Issue #90** (closed): [Epic #11] Mott elastic, step 2: energy-grid elastic tables into CrossSectionTable
- **Issue #78** (closed): Validation level 3: report the sensitivity of the sputter ratios to doubtful sets; Fig. 310 tick wording; checker gaps (follow-up to #70)
- **Issue #74** (closed): Commit EADL2017 subshell binding energies once their redistribution terms are ruled on
- **Issue #61** (closed): Investigate low sputter yields: Ar 1 keV -> Cu is 2.2x below RustBCA even with E_d = E_s (follow-up to #50)
- **Issue #37** (closed): [Epic #36] Profile-guided hot-path optimization
- **Issue #35** (closed): [Epic #33] Experimental stopping-data format and loader
- **Issue #32** (closed): [Epic #30] Physics manual + user guide (mdBook)
- **Issue #31** (closed): [Epic #30] Publish-readiness: packageable crates, cargo publish --dry-run in CI, CHANGELOG
- **Issue #29** (closed): [Epic #28] lindhard-py: pyo3 bindings with NumPy outputs
- **Issue #26** (closed): [Epic #25] Voxel-grid geometry
- **Issue #24** (closed): [Epic #22] Fluence stepping loop
- **Issue #23** (closed): [Epic #22] Composition–depth grid with volume relaxation
- **Issue #18** (closed): [Epic #11] Penn inelastic, step 1: single-pole Penn DIIMFP, IMFP and stopping from an OpticalElf
- **Issue #17** (closed): [Epic #11] Mott elastic, step 1: radial Dirac partial-wave solver for a screened potential

### 2026-10-06

- **PR #77**: feat: measured Ar -> Si, Ag, Au sputter yields with per-target comparison and generated interpretation
- **PR #75**: feat(electron): validated electron::data types for ELF, subshell binding energies and cross-section caches
- **PR #73**: docs: oracle ions/s comparison and profile (#9)
- **PR #72**: data: verify constants and atomic weights against CODATA 2022 and CIAAW; record E_s, density, E_d gaps
- **PR #71**: Level 3: measured Ar -> Cu sputter yields with independent double-check (#69)
- **PR #68**: docs: reflect #64 (opt-in weak collisions) in README and WORK_PLAN
- **PR #67**: bca: optional TRIDYN weak collisions beyond p_max, with an LNST partition check (#64)
- **PR #66**: docs: refresh README status and WORK_PLAN (post #62, #61)
- **Issue #70** (closed): Validation level 3: measured Ar sputtering yields for Si, Ag and Au (follow-up to #69)
- **Issue #69** (closed): Validation level 3: measured Ar sputtering yields (Cu, Si, Ag, Au) from open compilations
- **Issue #64** (closed): BCA: nuclear loss truncated at p_max while electronic loss is not, which over-damps low-energy recoils
- **Issue #16** (closed): [Epic #11] Electron interaction data model + optical ELF provenance
- **Issue #9** (closed): Benchmarks: hot-path criterion benches + ions/s vs oracles

### 2026-10-05

- **PR #63**: docs: repo hygiene pass (README status, stale refs, WORK_PLAN, directory READMEs)
- **PR #62**: Validation level 3: B-in-a-Si range data, provenance enforcement, stopping-input attribution
- **PR #60**: Docs: Biersack-Varelas join documented, not implemented (#44)
- **PR #59**: Validation level 2: RustBCA and OpenTRIM oracle adapters and first committed summaries
- **PR #57**: docs: verify screening lengths and impulse-approximation citation against primary sources
- **PR #56**: CLI: user-supplied stopping tables via [stopping]
- **PR #54**: CLI: tally outputs (range, damage, sputtering, escapes)
- **PR #53**: Add hot-path benches, thread-scaling throughput and benchmarks doc
- **PR #52**: Validation harness: level-1 checks in CI, oracle and experiment scaffolding (part of #8)
- **PR #49**: Stopping: single-oscillator density effect; document declined Barkas/shell/straggling terms
- **PR #48**: Tallies: range moments, Pearson IV / dual-Pearson, damage (NRT vs cascade), sputter/backscatter
- **PR #47**: CLI: TOML input, check/run, JSON/CSV output with full physics metadata
- **PR #46**: Amorphous 1D-layered BCA engine with full recoil cascades
- **PR #43**: Deterministic per-particle RNG streams + parallel driver
- **PR #42**: Electronic stopping models from published formulas (LS, OR, Bethe-Bloch, joined) + Bragg additivity
- **PR #40**: Screened interatomic potentials and scattering integral (#3)
- **PR #39**: docs: historical note, from the 1947 neutron histories to lindhard
- **PR #38**: WORK_PLAN: full roadmap (9 epics, Phase 1 issues filed)
- **PR #14**: Materials, units and physical constants
- **Issue #55** (closed): CLI: user-supplied stopping tables
- **Issue #50** (closed): Validation level 2: RustBCA and OpenTRIM oracle adapters and first committed summaries (follow-up to #8)
- **Issue #44** (closed): Biersack-Varelas interpolation joining low- and high-energy stopping (deferred from #4)
- **Issue #41** (closed): Stopping: add Barkas, shell and density-effect terms and a Chu/Yang-type straggling correction from cleanly sourced formulas
- **Issue #7** (closed): CLI: TOML input, JSON/CSV output with full physics metadata
- **Issue #6** (closed): Tallies: range moments, Pearson IV / dual-Pearson, damage (NRT vs full cascade), sputter/backscatter
- **Issue #5** (closed): Amorphous 1D-layered BCA engine with full recoil cascades
- **Issue #4** (closed): Electronic stopping models from published formulas (LS, OR, Bethe-Bloch, Biersack-Varelas) + Bragg additivity
- **Issue #3** (closed): Screened interatomic potentials + scattering integral (quadrature, magic formula, lookup tables)
- **Issue #2** (closed): Deterministic per-particle RNG streams + parallel driver
- **Issue #1** (closed): Materials, units and physical constants

### 2026-10-04

- **PR #13**: Public from 2026-10-04: drop private-until-M0 wording, full CI matrix

## 2026-10-07

- `lindhard-py` (#29): Python bindings with pyo3 and maturin. The classes
  wrap the library's own input structs (`lindhard::input`), so the TOML schema
  has one definition and `Run.from_toml()` / `Run.to_toml()` round-trip with
  the CLI. To make results identical by construction, `lindhard-cli` gained a
  library target (`sim::simulate`, the tally and the summary/CSV writers) that
  both the binary and the bindings call. `Run.run()` releases the GIL. The
  `numpy` crate is BSD-2-Clause and `target-lexicon` (a pyo3 build
  dependency) is Apache-2.0 WITH LLVM-exception; both are added to the
  `deny.toml` allow list.

## 2026-10-04

- Project started: a clean-room, MIT-licensed Rust engine for ion
  implantation and low-energy electron transport. The repo stays private until
  M0. The starting survey is in `docs/prior-art.md`.
- Milestone order: M0 amorphous ion core → M1 electron engine → M2
  crystalline implant → M3 reach (operator may reorder).
- Made public the same day (operator), ahead of the original "after M0"
  plan: private-repo Actions minutes were blocked by the org's billing
  limit, and public repos run standard runners free. A leak-pattern scan of
  the tree, the full history and the issues was clean before the flip. CI now
  runs on Linux x86-64, Linux arm64 and macOS, as fasterhenry's does.
