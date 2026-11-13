# Work Log

Chronological record of notable decisions and merges.

### 2026-11-13

- **PR #738**: refactor(geometry): share the ownership walk between mesh and CSG
- **Issue #735** (closed): Share the ownership walk between MeshGeometry and CsgGeometry (about 110 duplicated lines)
- **PR #737**: refactor(serde): share one beam_serde_module macro between the two exposure beams
- **Issue #733** (closed): Deduplicate the two exposure beam serde modules (beam_serde / exposure_beam)
- **PR #736**: test: require the four elastic reference fields and drop empty-fixture paths
- **Issue #734** (closed): Remove placeholder-fixture scaffolding and stale #93 markers from the elastic reference test
- **PR #731**: refactor(tests): call common::with_threads for the hand-built rayon pools
- **Issue #730** (closed): Use common::with_threads for the 48 hand-built rayon pools in lindhard/tests
- **Issue #33** (closed): Epic: Permissive stopping-power data library

### 2026-11-12

- **PR #724**: data: add Madsen 1953 proton stopping in Cu and Au from the paper's tables
- **Issue #483** (closed): [Epic #33] Seed measured stopping datasets from primary papers (first batch: H/He/B/P/As in Si, H/He in Cu and Au)
- **PR #721**: fix(sw): never store a redirected response in the offline cache
- **PR #720**: Bump Repo Skills 0.25.10 -> 0.26.1
- **PR #719**: feat(simulator): 3D trajectory view; fix the Setup column layout
- **PR #718**: feat(site): fold the range calculator into the simulator and run an example on load
- **PR #717**: refactor(exposure): share the electron transport beam profile and add rectangle
- **Issue #671** (closed): Decide: fold the range calculator into the simulator, and run the default case on page load?
- **Issue #608** (closed): Decide whether to collapse exposure::BeamProfile into transport::BeamProfile

### 2026-11-11

- **PR #716**: feat(site): restructure the landing page around a hero, trust rows and examples
- **PR #714**: feat(gallery): generate the quick-start figure from one wasm session
- **Issue #713** (closed): Landing hero: generate the quick-start figure (trajectories over depth profile) in gallery.py
- **Issue #665** (closed): Landing page: hero with a real figure, tools and physics as cards, a scannable trust section
- **PR #711**: test: cover the layer stack drawing at tilt 0 and with a long layer name
- **PR #710**: fix: reset hidden species on run start and make empty trajectory bounds 1:1
- **PR #706**: test: consolidate duplicated test helpers into tests/common
- **PR #705**: feat(simulator): trim plot axes, integer count ticks, species colours, radial plot, Rp marker
- **PR #704**: feat(simulator): readable trajectory cascades with axes, 1:1 scale and species toggles
- **PR #703**: test: assert stopping-residuals tables fit at 768 px
- **PR #700**: test: assert Setup's visible controls fit 1440x900, skipping closed details
- **PR #687**: feat(site): shared header, footer and design tokens across all site surfaces
- **PR #685**: feat(simulator): two-column desktop layout, single run bar, result tiles and plot grid
- **Issue #712** (closed): Simulator empty state: layer drawing may not span its row (flex-grow in metres)
- **Issue #708** (closed): Simulator trajectory view: reset hidden species per run; fix empty-segment lateralScale
- **Issue #699** (closed): Consolidate duplicated integration-test helpers into tests/common (with_threads, rel, simpson_log, close)
- **Issue #698** (closed): Simulator Setup: one-layer target ends 1333 px down at 1440x900; fit the 900 px viewport
- **Issue #696** (closed): docs: check stopping-residuals tables at 1440/768/390 px after #688
- **Issue #692** (closed): Simulator: "Add element" button throws ReferenceError; 1440 px single-screen layout not met
- **Issue #677** (closed): Simulator trajectory view: readable cascades, true aspect ratio, real axes
- **Issue #676** (closed): Simulator plots: axes that follow the data, integer count ticks, consistent species colours
- **Issue #667** (closed): Simulator: sharp, responsive canvases; axes that follow the data; a readable trajectory view
- **Issue #666** (closed): Simulator: one workspace on desktop, a single Run button, result tiles, and plots in a grid
- **Issue #664** (closed): One shared header, footer and visual system across the landing page, the book and the browser tools
- **PR #693**: fix: export el for the simulator's Add element button; add e2e tests
- **PR #701**: docs: regenerate screen_oxide gallery table and curve after ZBL constant change
- **PR #697**: docs: turn the exposure output-files table into a list and drop the last wide-table exception
- **Issue #695** (closed): Follow-on from PR #686: screen_oxide gallery numbers drift from the committed results
- **Issue #691** (closed): Simulator Setup: fix Add element ReferenceError; finish one-screen layout at 1440 px
- **Issue #663** (closed): Book: tables that need content changes to fit (output and input guides, stopping residuals, validity)
- **PR #690**: docs: reflow example band tables and shorten Pearson IV excerpt
- **PR #688**: docs: readable titles and summary table for the stopping-residuals page
- **PR #686**: docs(book): collapsed sidebar, readability stylesheet, dark-aware gallery figures
- **PR #683**: feat(simulator): size canvases from CSS box and pixel ratio, redraw on resize and theme change
- **Issue #675** (closed): Book: gallery and first-run code lines overflow the page
- **Issue #674** (closed): Simulator canvases: sharp at any pixel density, redraw on resize and theme change
- **Issue #673** (closed): Book: stopping-residuals page — shorter headers, summary table, readable headings
- **Issue #669** (closed): Book readability: collapsible sidebar with the user guide first, line length and heading rhythm, code blocks, gallery figures in dark mode
- **PR #689**: feat: compact simulator setup form with layer stack drawing
- **PR #684**: docs: turn long-prose book tables into lists and add header to coefficient tables
- **PR #682**: feat(wasm-calculator): fix dark mode, ticks, freeze, form and escaped-ion count
- **PR #681**: fix(wasm): label unfollowed-recoil energy and correct damage-plot note in quick mode
- **PR #680**: feat(book): style tables and display equations to fit their column
- **PR #679**: fix: link the shipped simulator from the landing page and fail the build if it is missing
- **Issue #672** (closed): Book: hand-written tables to definition lists (output, input, validity, potentials)
- **Issue #670** (closed): Simulator quick mode: the damage plot contradicts its note, and the energy partition reads as wrong
- **Issue #668** (closed): Simulator Setup: fields sized to their content, clearer layer controls, and a drawn layer stack
- **Issue #662** (closed): Book: tables do not fit their content, and equations pan the page sideways on a phone
- **Issue #661** (closed): Range calculator: unreadable in dark mode, and other small defects
- **Issue #660** (closed): Landing page says the trajectory simulator is "planned" although it is published
- **PR #617**: validation: attribute the pure-Cu sputter shortfall and compare engine Ni with stored sets
- **Issue #544** (closed): [Epic #22] Attribute the Ar → Cu 0.67 yield ratio and add the engine-Ni row, without tuning (split from #511, part 2)

### 2026-11-10

- **PR #657**: feat: add service worker for offline calculator use
- **PR #656**: fix(web): give Si and SiO2 presets an illustrative Si displacement energy
- **PR #655**: test(site): Playwright browser tests of the site, run in CI before and after publishing
- **PR #652**: feat: add stopping = "srim" model reading a cited coefficient file
- **PR #649**: feat: report electronic, lattice and recoil energy by depth
- **PR #648**: feat: tally NRT and Kinchin-Pease displacements by PKA depth
- **PR #647**: feat: weighted primary-energy spectrum for [beam]
- **PR #646**: feat: add ideal-gas material spec with a gas phase flag
- **PR #645**: docs: record that BCA transport does not sample straggling
- **Issue #654** (closed): Simulator: adding a Si, SiO2 or GaAs preset layer leaves the setup unrunnable until a displacement energy is typed
- **Issue #653** (closed): Simulator page does not work offline after its first load
- **Issue #641** (closed): [Epic #630] PKA-spectrum primaries (neutron-recoil style cascades)
- **Issue #640** (closed): [Epic #630] Quick-calculation mode: K-P/NRT damage without cascades (verify, expose)
- **Issue #639** (closed): [Epic #630] Ionization, phonon and energy-to-recoil depth profiles (verify, add missing)
- **Issue #636** (closed): [Epic #630] Gas targets: density from pressure and temperature, gas-phase flag
- **Issue #632** (closed): [Epic #630] `stopping = "srim"`: SRIM-book electronic stopping as a selectable model
- **PR #635**: docs: data is not licensed, code is; admit cited data from any source (#629)
- **Issue #629** (closed): Data policy: admit SRIM data (book and srim.org), including computed tables; rewrite the CLAUDE.md rule and README positioning
- **PR #627**: refactor: consolidate duplicated [f64; 3] helpers into a vec3 module
- **Issue #625** (closed): Simplify: consolidate duplicated 3-vector helpers (dot/cross/norm/normalize) into one module
- **PR #584**: feat: add per-layer dose, depth-dose and interface transmission tallies for electrons
- **Issue #547** (closed): Electron engine: dose to named buried layers, depth-dose near interfaces, and transmission through a stack
- **PR #624**: refactor: split ElectronInput::resolve_in into per-stage helpers
- **PR #623**: feat: direct Monte Carlo exposure mode compared with the PSF convolution
- **PR #621**: test: add density-free Al2O3 sum-rule integrals and gap rows
- **PR #585**: feat: add multi-shot exposure dose, iso-dose contours and lindhard exposure command
- **Issue #622** (closed): Simplify: split 630-line ElectronInput::resolve_in into per-stage helpers
- **Issue #618** (closed): Al2O3 optical ELF: sum rules and IMFP comparison
- **Issue #555** (closed): Direct Monte Carlo exposure of a small shot list, to check the PSF convolution on non-uniform targets
- **Issue #552** (closed): Sourced optical energy-loss functions and band parameters for more metals, nitrides, dielectrics and resists
- **Issue #551** (closed): Pattern exposure from the point-spread function: dose from a shot list, iso-dose contours, direct Monte Carlo check
- **PR #615**: feat: position-dependent PSF scan over voxel targets
- **PR #614**: feat: electron voxel damage map table and damage_export for electron runs
- **PR #612**: feat: report pure-Ni yield against stored Ar -> Ni sets
- **PR #611**: docs: record SiCOH and SiCN optical ELF as gaps
- **PR #610**: refactor: share Gaussian and disc beam offset helpers across transport and exposure
- **PR #609**: refactor: unify electron Face with geometry::Face via re-export
- **PR #607**: docs: record HfO2 optical ELF as a gap with routes tried
- **PR #604**: docs: record Si3N4 optical ELF as a gap with routes tried
- **PR #603**: test: per-species (He, Ne, Ga, Xe) 2-D dynamic profile runs
- **PR #601**: feat: versioned damage-export file with provenance header and NumPy round trip
- **PR #600**: docs: record SiO2 optical ELF as a gap with routes tried
- **PR #599**: docs: record full Al2O3 Table 9 second read against the scan
- **PR #596**: test: share synthetic table builder and bare-Si helper across electron tests
- **PR #595**: feat: lateral implant/damage profiles and film cut-through for the 2-D dynamic run
- **PR #593**: data: add measured Ar -> Ni sputter yields (IPPJ-AM-32 p. 103) with NIFS cross-check
- **PR #592**: docs: new fronts run alongside the 0.1.0 close-out; close-out step 2 done
- **PR #591**: fix: set empty species/defect fields in electron voxel_report
- **PR #586**: feat: split electron deposition into ionisation and lattice, add pairs and knock-on report
- **PR #579**: feat: electron beam landing point, lateral profile and divergence
- **PR #577**: data: add amorphous Al2O3 optical ELF (Hagemann 1975, Table 9)
- **Issue #602** (closed): Electron voxel damage map: input table and Python/CLI export
- **Issue #597** (closed): Regenerate the #465 Cu-Ni diagnostic now that Ar -> Ni sets are stored
- **Issue #572** (closed): [Epic #25] Position-dependent electron PSF report and buried-edge validation (part of #549)
- **Issue #568** (closed): Sourced optical ELFs: SiO2, Si3N4, SiCOH, SiCN, HfO2, Al2O3
- **Issue #549** (closed): Electron transport on 2-D/3-D voxel targets with a finite beam, and a position-dependent point-spread report
- **Issue #606** (closed): Collapse duplicate BeamProfile in tally::exposure into electron::transport::BeamProfile
- **Issue #605** (closed): Unify duplicate Face enum: electron::transport re-exports geometry::Face
- **Issue #598** (closed): data: second read of Al2O3 optical ELF rows above 20.5 eV against the DESY SR-74/7 scan
- **Issue #590** (closed): main does not compile: tally/electron.rs VoxelReport/VoxelRecord missing fields after #560/#570 merge
- **Issue #589** (closed): Consolidate duplicated synthetic-table and bare-Si helpers in electron tests
- **Issue #571** (closed): [Epic #25] Electron beam with a finite lateral profile, landing position and divergence (part of #549)
- **Issue #563** (closed): [Part of #553] Versioned damage-map export (ion and electron) with provenance header and Python round trip
- **Issue #558** (closed): [Part of #553] Electron runs: ionisation vs lattice energy, e-h pairs, and knock-on threshold report
- **Issue #553** (closed): Damage outputs: ionisation dose and knock-on threshold for electrons, per-species defect maps for ions, and a documented export
- **Issue #543** (closed): [Epic #22] Store Ar → Ni sputter-yield datasets from IPPJ-AM-32 (split from #511, part 1)
- **PR #583**: feat: add a 2-D dynamic run with finite beam, height-field surface and redeposition
- **PR #582**: feat: run electron transport on any Geometry with a per-voxel deposition tally
- **PR #581**: feat: per-element defect maps in the voxel tally, written by the CLI
- **PR #575**: docs: record TiN and TaN optical-data gaps in the provenance inventory
- **PR #573**: docs: record W, Ta, Ti, Co and Ru optical-data sources as gaps
- **Issue #588** (closed): CI flake: Wheel and CLI (linux-x86_64) fails pulling manylinux image from quay.io (504)
- **Issue #570** (closed): [Epic #25] Electron transport on any Geometry (voxel targets) with a 3-D energy-deposition tally (part of #549)
- **Issue #564** (closed): Sourced optical ELFs: TaN and TiN
- **Issue #560** (closed): [Part of #553] Ion runs: per-element vacancy, interstitial and replacement maps per voxel, written by the CLI
- **Issue #559** (closed): [Epic #22] Focused ion beams, part 1: 2-D dynamic voxel target with finite beam, local-angle yield and redeposition
- **Issue #556** (closed): Sourced optical ELFs: metals Co, Ru, W, Ta, Ti
- **Issue #10** (closed): Epic: M0, amorphous ion core, validated
- **PR #580**: feat: exposure statistics for small volumes via compound-Poisson tally
- **PR #578**: feat: tally stopped electrons and holes per layer and cell with an exact charge balance
- **PR #574**: docs: record HSQ as a sourced-ELF gap in the optical inventory
- **Issue #569** (closed): Sourced optical ELFs: PMMA and hydrogen silsesquioxane (HSQ)
- **Issue #550** (closed): Exposure statistics: distribution of energy deposited in a small volume by a finite number of electrons

### 2026-11-09

- **PR #542**: docs: record operator rulings on 0.1.0 scope and publishing
- **PR #539**: chore(deps): bump rizzma from 1.13.7 to 1.14.1 in the cargo-minor-patch group across 1 directory
- **Issue #536** (closed): Universal screening length: use the 0.8853 prefactor printed in the SRIM book's Eq. 2-73
- **PR #538**: fix: use the Thomas-Fermi constant for the universal screening length
- **PR #541**: docs: roadmap refresh — close out and release 0.1.0
- **PR #540**: docs: verify Molière constants and TRIM-chapter formulas against the SRIM book
- **PR #537**: docs: verify ZBL constants against the SRIM book, Ch. 2
- **PR #535**: docs: tidy SRIM book wording (follow-up to #534)

### 2026-11-08

- **PR #534**: docs: allow SRIM as a comparison reference, never as model input

### 2026-11-07

- **PR #527**: ci: report RustSec advisories weekly via a tracking issue
- **PR #533**: Bump Repo Skills 0.25.1 -> 0.25.10

### 2026-11-06

- **Issue #531** (closed): Consolidate copy-pasted Si/SiO2 Material test fixtures into tests/common
- **PR #532**: test: share Si and SiO2 material fixtures via tests/common
- **Issue #529** (closed): Dedupe the Si-Si ZBL ScatteringTable test helper copied across 15 integration tests
- **PR #530**: test: share Si-Si ZBL table helper across integration tests
- **Issue #192** (closed): Test: malformed and hostile input values must error, never panic or hang
- **PR #528**: test(input): add subprocess-contained malformed-input corpus; cap elastic grid size
- **Issue #157** (closed): Auditor Capability Request: Rust and Python toolchains unavailable on auditor host
- **Issue #193** (closed): Auditor process proposal: check current main for conflicts before review handoff
- **Issue #525** (closed): Add the H and He water/polystyrene compound points to the stopping residual report
- **PR #526**: feat(validation): add water and polystyrene compound points to the stopping residual report

### 2026-11-02

- **Issue #522** (closed): lindhard.org: public demo site (landing, in-browser WASM calculator, gallery, validation) on Cloudflare Pages
- **PR #523**: feat(site): Cloudflare Pages build script, CI job and validation/status/about pages
- **PR #524**: chore(deps): bump sha2 from 0.11.0 to 0.11.1 in the cargo-minor-patch group across 1 directory

### 2026-11-01

- **PR #519**: feat: add live distribution plots and CSV/JSON/PNG exports to the simulator Run view
- **PR #521**: Bump Repo Skills 0.24.7 -> 0.25.1

### 2026-10-31

- **Issue #512** (closed): [Epic #22] Run the Sb-in-Si high-fluence profile comparison against the stored Titze et al. SIMS data (criterion b)
- **PR #520**: validation: compare the fluence-dependent engine with the Sb-in-Si implant profile
- **Issue #513** (closed): [Epic #12] Step 22: impact-parameter-dependent electronic stopping (Oen–Robinson local term), calibrated per ion/crystal from primary sources
- **PR #517**: test: record channeled vs random electronic loss for B in Si
- **Issue #514** (closed): [Epic #12] Step 23: amorphous overlayers (native / screen oxide) over a crystalline substrate with exact interface handling
- **PR #518**: test: check crystal interface identity and per-layer bookkeeping in mixed targets
- **PR #510**: feat: add live counters, trajectory view and ion cap to the simulator Run view
- **PR #506**: feat(wasm): simulator page with layered-target setup, TOML round-trip and worker run

### 2026-10-30

- **Issue #498** (closed): [Epic #33] Digitize compound-target stopping plots (water, polystyrene, hydrocarbons, polymers)
- **PR #508**: data(stopping): count-check report and hand check for the srim.org compound figures
- **Issue #496** (closed): [Epic #33] Digitize per-pair stopping plots for Li and heavier ions
- **PR #509**: data(stopping): per-plot count check for the Li per-pair plots
- **Issue #495** (closed): [Epic #33] Digitize measured projected-range plots (B, P, As in Si; He in air)
- **PR #507**: data(ranges): count-check report and hand check for the srim.org range figures

### 2026-10-29

- **Issue #479** (closed): lindhard-wasm: Simulation session for the full engine (layered compound targets, chunked runs, IonReport and trajectories as JSON)
- **PR #482**: feat: add stepwise Simulation session to lindhard-wasm
- **Issue #470** (closed): Project landing page on the Pages site with validation status
- **PR #480**: feat: landing page over the book root, with a relative-link check
- **Issue #488** (closed): Supply primary papers to complete the first batch of measured stopping datasets
- **PR #505**: feat(data): add per-paper digitized stopping datasets for H and He in Si and Cu
- **Issue #502** (closed): Fix stale srim_org_h_he.csv name in stopping compilation README
- **PR #504**: docs: point compilation README at renamed srim_org_stopping.csv
- **PR #500**: Recover hidden srim.org points; add Li, ranges and compound figures (#486)
- **Issue #485** (closed): [Epic #33] Preserve the srim.org experimental stopping bibliography (~2,300 papers) as CC0 data
- **PR #494**: feat: add srim.org experimental-paper bibliography for pairs and compounds
- **Issue #499** (closed): [Epic #33] Admit the digitized-compilation origin in the stopping loader and residual report
- **PR #501**: feat(stopping): admit the digitized-compilation origin in loader and report
- **PR #493**: Digitize the srim.org H and He stopping compilation (7,241 measured points)
- **Issue #489** (closed): [Epic #33] Book chapter: history and literature of ion stopping, with a page on James F. Ziegler (1935–2024)
- **PR #492**: docs: add history-and-literature chapter and Ziegler page to the book
- **PR #484**: docs(clean-room): allow measured points from the srim.org compilation (operator ruling)
- **Issue #179** (closed): [Epic #12] Crystal step 19b: hexagonal lattices (wurtzite GaN, 4H/6H-SiC)
- **PR #491**: test: derive 4H-SiC density tolerance from measured spread
- **PR #490**: docs: record first-batch measured-stopping search and unreadable papers

### 2026-10-28

- **Issue #156** (closed): Electron CLI: subshell binding-energy input key once inner-shell channels reach the transport loop
- **PR #481**: feat: add inner_shells key for electron materials in lindhard run
- **Issue #336** (closed): CI: pin the clippy gate and move latest-stable lint drift to a scheduled canary
- **Issue #444** (closed): Verify the WASM range-calculator page in a real browser (follow-up to #432)
- **Issue #93** (closed): [Epic #11] Elastic validation: total and transport cross sections vs published partial-wave values (C, Si, Cu, Au)
- **Issue #301** (closed): Au and Cu secondary-electron yield excess of the Mermin model: the candidates #242 could not test
- **Issue #432** (closed): [Epic #28] WASM build and static range-calculator page
- **PR #439**: feat: add lindhard-wasm crate and static WASM range-calculator page
- **Issue #431** (closed): [Epic #33] Stopping residual report: one command, published in the docs
- **PR #443**: feat: generate the stopping-residual report page and check it in CI
- **PR #475**: chore(deps): bump the github-actions group across 1 directory with 2 updates

### 2026-10-27

- **Issue #474** (closed): Remove one-off #131 start-rule measurement harness from elastic tests
- **PR #478**: test(elastic): remove the one-off start-rule measurement harness
- **Issue #464** (closed): Benchmark lindhard-gpu-proto on real GPU hardware and add the rows to gpu-feasibility.md
- **PR #477**: docs(gpu-feasibility): add Apple M3 Ultra (Metal) rows, restate defer (#464)

### 2026-10-26

- **PR #476**: chore(deps): bump the cargo-minor-patch group across 1 directory with 2 updates

### 2026-10-24

- **Issue #465** (closed): [Epic #22] Diagnose why the dynamic engine is low against the Cu–Ni measurements (no tuning)
- **PR #473**: docs: diagnose the Cu-Ni shortfall with static pure-element and fixed-composition runs
- **Issue #468** (closed): Opt-in ion trajectory recorder in the engine (per-history collision points)
- **PR #472**: feat: add opt-in trajectory recorder tally for ion and recoil paths
- **Issue #466** (closed): [Epic #22] Source a measured high-fluence implant depth profile from an openable publication
- **PR #471**: docs: add measured high-fluence Sb-in-Si SIMS profile with two-read digitization
- **Issue #36** (closed): Epic: Performance (profile-guided, SIMD, GPU feasibility)

### 2026-10-23

- **Issue #462** (closed): CI: ubuntu-latest moved to 26.04; x86_64 Linux golden digests in crystal_amorphous_identity fail (and macOS job hits 30 min timeout)
- **PR #463**: ci: raise rust job timeout to 45 min for macOS headroom
- **Issue #460** (closed): CI: Rust job on main fails linking lindhard-py (unable to find -lpython3.12)
- **PR #461**: ci: pin the x86_64 Rust test leg to ubuntu-24.04
- **Issue #425** (closed): [Epic #36] GPU feasibility study: written decision with a prototype benchmark
- **PR #427**: docs: GPU feasibility decision (defer) with an out-of-workspace wgpu prototype

### 2026-10-21

- **PR #459**: Bump Repo Skills 0.24.3 -> 0.24.7

### 2026-10-20

- **Issue #448** (closed): [Epic #22] Validate dynamic Cu–Ni targets against Ho and Betz measurements
- **PR #456**: test: compare the fluence-dependent engine with measured Cu-Ni yields and surface fractions
- **Issue #457** (closed): Audit Si shell inputs and remaining Au backscatter discrepancy; prepare comparison plan
- **PR #458**: docs(validation): Si shell input audit, 199 eV join review, Au audit and comparison plan
- **Issue #454** (closed): Consolidate 14 copies of the 53-bit uniform helper into rng.rs
- **PR #455**: refactor: share one 53-bit uniform helper from rng.rs
- **Issue #452** (closed): [Epic #11] Elastic validation follow-up: test whether the potential difference explains the σ_el gap vs NIST SRD 64
- **PR #453**: test: report in-tree potential sensitivity of elastic gaps vs NIST SRD 64
- **PR #447**: Preserve Ho and Betz Cu–Ni experimental data with provenance
- **PR #451**: test: validate elastic cross sections against NIST SRD 64 and report all 16 gaps
- **Issue #449** (closed): Diagnose the silicon backscatter excess at 20–30 keV (#148)
- **PR #450**: feat(validation): Si backscatter with inner-shell channels (library-level example)

### 2026-10-16

- **Issue #437** (closed): [Epic #22] Validation: fluence-dependent sputter yield and high-fluence implant profile
- **Issue #435** (closed): [Epic #25] Validation: masked-implant lateral edge profile vs the analytic result
- **PR #446**: test: add masked-implant edge-profile validation on a voxel target
- **PR #445**: docs: record the fluence-dependent validation data gap
- **Issue #433** (closed): [Epic #25] Sparse 3D voxel tallies with CSV and VTK export
- **PR #440**: feat: add a sparse 3D voxel tally with CSV and VTK export
- **Issue #436** (closed): [Epic #22] Implanted-species retention limit (saturation and re-emission)
- **PR #442**: feat: add a per-species retention limit with re-emission to the fluence loop
- **PR #441**: docs: record klayout-tools pin and the missing polygon export in 0.7.0
- **Issue #30** (closed): Epic: Release engineering and documentation
- **Issue #430** (closed): [Epic #33] Per-pair stopping correction fits with bootstrap uncertainty and held-out residuals
- **PR #438**: feat: per-pair stopping correction fits with bootstrap band and held-out residuals
- **PR #423**: Bump Repo Skills 0.24.1 -> 0.24.3
- **Issue #424** (closed): [Epic #36] Re-profile after hot-path work and decide SoA/SIMD with numbers
- **PR #428**: perf: compute the collision azimuth with a short series instead of libm sin and cos
- **Issue #426** (closed): [Epic #30] Examples gallery: screen-oxide implant, sputtering yield vs angle, electron deposition in resist on Si
- **PR #429**: docs: add examples gallery with screen-oxide, sputter-angle and resist pages
- **Issue #217** (closed): Decide on pure-Rust libm: crystal BCA golden run diverges on aarch64 macOS
- **PR #422**: test: pin a labelled aarch64-macOS crystal golden set


### 2026-10-15

- **PR #419**: docs: record Au phonon energy-loss inputs as a provenance gap
- **Issue #313** (closed): SE yield candidates (#301 item 2): acoustic-phonon energy loss per event, or the gap
- **PR #417**: Bump Repo Skills 0.21.2 -> 0.24.1
- **Issue #51** (closed): Validation level 3: first published range datasets (B, P, As in amorphous Si) and stopping-input attribution (follow-up to #8)

### 2026-10-13

- **PR #415**: Dependabot: render dependabot.yml from the org policy
- **PR #414**: chore(deps): bump the github-actions group across 1 directory with 6 updates
- **PR #413**: Remove Renovate config: Dependabot takes over
- **PR #412**: chore(deps): bump criterion from 0.7.0 to 0.8.2
- **PR #411**: Dependabot: render dependabot.yml from the org policy
- **PR #410**: fix(clippy): satisfy beta needless_range_loop in BVH traversal and tests
- **Issue #409** (closed): CI: clippy canary failing on latest stable or beta Rust

### 2026-10-11

- **PR #407**: validation: fingerprint executable, input and tables for ion-run reuse (#406)
- **PR #404**: fix: bound the adaptive sputter rerun loop and report precision status
- **Issue #406** (closed): Ion validation: fingerprint executable and referenced data before resuming runs
- **Issue #403** (closed): Validation: bound adaptive sputter retries and handle zero-yield observations
- **Issue #287** (closed): Re-run the secondary-yield tables and the electron oracle summaries with the band-bottom inelastic tables of #241

### 2026-10-11

- **PR #401**: docs(crystal): record opened 4H-SiC density measurements; 4H gap stands
- **PR #400**: Python: expose electron transport through a TOML-first ElectronRun API
- **PR #399**: Python: expose fluence-dependent ion runs through run_dynamic (#397)
- **PR #396**: perf(rng): bound live partial tallies in the particle driver
- **PR #395**: docs(cli): distinguish Git dirty scope from provenance watch scope
- **PR #394**: Backscatter (#148 part 2): full-Penn rerun of the eta comparison
- **PR #392**: validation: run the cutoff-band-bottom δ(E) row to 1 keV for Al and Si (#380)
- **PR #390**: fix(cli): refresh git_describe after incremental source and ref changes
- **PR #389**: docs(oracles): update Nebula inelastic-model mismatch text for #241/#339 (#386)
- **PR #387**: refactor(electron): split inelastic/table.rs into density, refine, build and sampler modules
- **PR #384**: docs(benchmarks): mark electron bench digests stale after band-bottom changes (#366)
- **PR #383**: validation: re-run the lindhard side of the electron oracles (#365)
- **PR #382**: docs(validation): attribute the Si δ(E) rise of #364 by bisect (#379)
- **PR #325**: fix(table-cache): lock each entry across lookup, build and publication
- **PR #87**: feat: [physics] stopping = "none" for nuclear-only runs and like-for-like oracle checks
- **Issue #398** (closed): Python: expose electron transport through a TOML-first ElectronRun API
- **Issue #397** (closed): Python: expose fluence-dependent ion runs through run_dynamic
- **Issue #393** (closed): Bound temporary tally storage in the deterministic particle driver
- **Issue #391** (closed): docs(cli): distinguish Git dirty scope from provenance invalidation scope
- **Issue #388** (closed): Build provenance: refresh git_describe after incremental source and ref changes
- **Issue #386** (closed): run_electron.py: Nebula mismatches text still describes pre-#241 single-pole Penn ("model Fermi energy 0")
- **Issue #380** (closed): Complete the cutoff-band-bottom δ(E) row at the regenerated build (Al above 400 eV, Cu/Au above 1000 eV, Si above 400 eV)
- **Issue #379** (closed): Si secondary-electron yield rose 77–187 % after #241/#339/#342 (found regenerating #364)
- **Issue #375** (closed): Refactor: split lindhard/src/electron/inelastic/table.rs into density, refine, build and sampler modules (no behaviour change)
- **Issue #366** (closed): Refresh the six electron benchmark digests after the band-bottom changes
- **Issue #365** (closed): Regenerate the six lindhard electron-oracle summaries after the band-bottom changes
- **Issue #326** (closed): CI: clippy chunks_exact_to_as_chunks in lindhard/src/geometry/mesh.rs fails on current stable
- **Issue #305** (closed): Electron table cache: coordinate concurrent builders per entry
- **Issue #89** (closed): Validation: regenerate OpenTRIM and RustBCA summaries for the *_nuclear level-2 problems (follow-up to #58)
- **Issue #58** (closed): CLI: an electronic-loss-off stopping choice, for like-for-like OpenTRIM comparisons (follow-up to #50)
- **PR #381**: docs(provenance): final 4H-SiC density source search; 6H precision density (Part of #179)
- **PR #378**: validation: regenerate the δ(E) tables on band-bottom inelastic tables (#364)
- **PR #377**: refactor(input): split input/electron.rs into schema, resolve, resolved and tests modules
- **PR #376**: fix(inelastic): open the full-Penn row window below elf_min (#361)
- **Issue #374** (closed): Refactor: split lindhard/src/input/electron.rs into schema, resolve and resolved modules (no behaviour change)
- **Issue #364** (closed): Regenerate Al, Cu, and Au secondary-electron yield tables after the band-bottom changes
- **Issue #361** (closed): Full Penn tables: the loss window opens at the lowest ELF energy, but the full-Penn DIIMFP is not zero below it

### 2026-10-10

- **PR #373**: refactor(electron): split data.rs into optical, binding, endf6 and table modules (#349)
- **PR #370**: validation: preserve and surface event-cap truncation diagnostics
- **Issue #369** (closed): Electron validation: preserve and surface event-cap truncation diagnostics
- **Issue #349** (closed): Refactor: split lindhard/src/electron/data.rs into optical, binding, endf6 and table modules (no behaviour change)
- **PR #371**: feat(full_penn): report ω_p exhaustion on the direct DIIMFP and IMFP paths (#368)
- **PR #367**: docs(validation): knot-derived width cap does not bound the Mermin f-sum (#337, negative result)
- **PR #363**: fix(inelastic): start the rate refinement at the model's rate onset (#357)
- **PR #362**: Full Penn: bound the ω_p and profile bisections; Si table builds in 224 s (#298)
- **PR #359**: docs: stable README status, trimmed architecture module table
- **PR #358**: fix(inelastic): stop loss-density bisection at floating-point resolution (#345)
- **PR #355**: feat(inelastic): refine band-bottom tables in lindhard run (#339)
- **PR #354**: renovate: move onto the 2AMLogic org preset
- **PR #352**: feat(crystal): source GaN cell and measured density; record the 4H-SiC density gap
- **PR #299**: Full Penn: make one inelastic table build finish in minutes (#256)
- **PR #236**: ci: build, test and publish Python wheels (abi3, trusted publishing)
- **PR #198**: Crystal step 19b: hexagonal lattices (wurtzite GaN, 4H/6H-SiC) (#179)
- **Issue #368** (closed): Report ω_p exhaustion on the direct full-Penn DIIMFP and IMFP paths (follow-up to #298)
- **Issue #357** (closed): Rate refinement fails on an optical ELF whose onset lies inside the 2-20 eV window
- **Issue #356** (closed): Dependency Dashboard
- **Issue #350** (closed): Docs: replace the README status narrative and trim architecture.md module-table rows to stable summaries
- **Issue #345** (closed): Mermin inelastic table: Si row at T = 1000 eV exhausts MAX_NODES in the loss-density refinement
- **Issue #339** (closed): Resolve the inelastic rate from E_F upward on the band-bottom axis: extra table rows above the Fermi level
- **Issue #337** (closed): Mermin fit: measure a knot-derived upper bound on the oscillator width against the #311 wide oscillators
- **Issue #298** (closed): Full Penn: isolated (q, ω) points where one loss-function evaluation takes seconds (Si table build does not finish)
- **Issue #291** (closed): Electron table grid: the default 10 points per decade under-resolves the inelastic rate just above the Fermi level on the band-bottom axis
- **Issue #256** (closed): Full Penn: make one inelastic table build finish in usable time (split from #169)
- **Issue #233** (closed): Python wheels: maturin multi-platform builds and PyPI trusted publishing (Epic #28, Phase 2)
- **PR #351**: feat(inelastic): rate refinement of band-bottom tables above E_F (part of #339)
- **PR #348**: PSF parameter CSV: export convergence and iteration diagnostics
- **PR #346**: ci: pin the clippy gate and add a stable/beta clippy canary (#336)
- **PR #344**: fix(electron): open a Mermin table row's loss window at the rate's lower limit (#342)
- **PR #343**: fix(electron): store inelastic table rows with no loss window above the ELF table as empty (#340)
- **PR #341**: docs(electron): measure the band-bottom table grid above E_F (#291)
- **PR #338**: docs(electron): measure Mermin fit weighting options against the wide oscillator (#311)
- **PR #335**: backscatter: fingerprint optical_elf contents before reusing runs
- **PR #334**: crystal: density evidence per material; 4H-SiC/GaN density source gap (Refs #179)
- **PR #333**: Expose the Mermin fit options per material in the electron input (#306)
- **PR #332**: Plot backscatter eta(E) validation figure with rizzma (#279)
- **PR #330**: validation: Cu seed 2 across hosts — not isolated; record hosts and table hashes (#319)
- **PR #328**: fix(py): share CLI output publication in RunResult.write (#327)
- **Issue #347** (closed): PSF parameter CSV: export convergence and iteration diagnostics
- **Issue #342** (closed): Si Mermin rate table has an artificial zero row at 15.5 eV because the loss interval starts below the ELF table
- **Issue #340** (closed): Si Mermin inelastic table refused for a band-bottom row within 0.5 eV above E_F (negative loss-density integral)
- **Issue #331** (closed): Backscatter run reuse ignores optical-ELF content changes
- **Issue #327** (closed): Python RunResult.write leaves prior output files in place while replacing a run
- **Issue #319** (closed): Cu secondary-yield baseline differs by 5 sigma across machines with identical committed inputs
- **Issue #311** (closed): Mermin wide-loss oscillator fit is biased toward the dense low-energy ELF rows
- **Issue #306** (closed): Expose Mermin fit controls per electron material
- **Issue #279** (closed): Plot simulation and validation output with rizzma

### 2026-10-10

- **Issue #149** (closed): [Epic #11] Validation: secondary-electron yield δ(E) vs published measurements
- **Issue #241** (closed): Electron inelastic tables: build rows on the band-bottom axis with the band's Fermi energy (cstool convention) instead of relying on the E - E_F clamp
- **Issue #321** (closed): CLI: invalidate prior summaries when a rerun starts replacing outputs
- **PR #289**: Electron inelastic tables: band-bottom axis with the band's Fermi energy (#241)
- **PR #322**: CLI: invalidate prior summaries when a rerun starts replacing outputs
- **PR #324**: Backscatter (#148): Si in the fast-secondary sensitivity of eta
- **Issue #300** (closed): Cu Mermin inelastic mean free path is 15 to 41 % longer than TPP 2011 and than the single-pole model on the same ELF
- **Issue #307** (closed): Mermin fit with 10 or 16 oscillators puts a sub-meV oscillator between ELF knots (f-sum x7000 for Cu)
- **Issue #312** (closed): SE yield candidates (#301 item 1): acoustic-phonon scattering for Cu, or the gap
- **Issue #314** (closed): SE yield candidates (#301 item 3): a second Au ELF, licence verdict first, or the gap
- **Issue #315** (closed): SE yield candidates (#301 item 4): low-energy elastic corrections (correlation-polarization, muffin-tin), or the gap
- **PR #308**: validation: Cu Mermin IMFP excess is the default 3-oscillator fit (#300)
- **PR #310**: fix: bound Mermin fit widths below by the local knot spacing (#307)
- **PR #316**: feat(electron): secondary binding switch; run #301 item 5 (Azzolini B), split items 1-4
- **PR #317**: docs(provenance): Cu acoustic-phonon inputs from Arat 2019 and the gap (#312)
- **PR #318**: Au ELF from Hagemann Table 6 as a second ELF; licence verdicts; au-elf-table6 run (#314)
- **PR #320**: docs(provenance): record the #301 item 4 elastic rows as gaps (#315)


### 2026-10-10

- **PR #304**: feat(electron): sample inner-shell ionisation channels in the transport (#273)
- **PR #284**: Screen-oxide substrate channeling tail characterization (#280)
- **Issue #280** (closed): Crystal validation: characterize screen-oxide effects on substrate channeling tails
- **Issue #273** (closed): Electron transport: sample inner-shell ionisation channels in the transport loop

### 2026-10-09

- **PR #302**: validation: no tested candidate explains the Au and Cu Mermin secondary-yield excess (#242)
- **Issue #242** (closed): Au secondary-electron yield stays 77 % above the measured δ_max with the Mermin model (residual of #173)
- **PR #297**: cli: write dynamic summary last (#292)
- **PR #296**: fix: crystal local loss, explain the fixed-direction excess and assert the direction average (#250)
- **PR #295**: validation: delta(E) on DHFS, Si curves, band text (Part of #149)
- **PR #294**: Backscatter (#148): fast-secondary sensitivity of eta for Al, Cu and Au
- **PR #293**: feat: hydrogenic L-subshell ELFs for the inner-shell builder (#245)
- **PR #288**: Ion beam: connect divergence sampling to transport and shared TOML
- **PR #247**: CLI: expose the radial PSF tally in the [electron] section
- **Issue #292** (closed): Dynamic CLI: publish the summary only after required CSV writes succeed
- **Issue #250** (closed): Crystal flight: local/nonlocal electronic loss ratio is 5-8 % above the amorphous value at the same p_max in random directions
- **Issue #245** (closed): Penn inelastic: hydrogenic L-subshell ELFs for the inner-shell builder (follow-up to #135)
- **Issue #285** (closed): Ion beam: connect existing divergence sampling to transport and shared TOML
- **Issue #165** (closed): CLI: expose the radial PSF tally in the [electron] section
- **Issue #201** (closed): Pilot E_s tuning set: fit and held-out evaluation (follow-up to #80)
- **Issue #79** (closed): bca: opt-in time-integral path length for nonlocal electronic loss (TRIM tau = p tan(theta/2))
- **Issue #135** (closed): Penn inelastic: build per-shell optical ELFs for ShellResolvedChannels from a published source
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
