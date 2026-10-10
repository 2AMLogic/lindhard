# Changelog

All notable changes to `lindhard` and `lindhard-cli` are recorded here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the crates use [Semantic Versioning](https://semver.org/) (both crates share
one version).

## [Unreleased]

### Fixed

- Python `RunResult.write()` now follows the command's output lifecycle: it
  removes the previous `summary.json` before replacing any CSV, removes a
  stale `ions.csv` when per-ion output is off, and writes the summary last.
  Both callers share `lindhard_cli::publish::publish_static_ion` (#327).
- A run into a reused output directory removes its mode's previous summary
  before replacing any CSV, and publishes the new summary through a
  temporary file and rename, so a failed rerun no longer leaves an old summary
  beside replaced CSVs (#321).
- Dynamic runs write `dynamic_summary.json` last, after
  `dynamic_steps.csv` and the composition CSV, like ion and electron runs
  (#292). A failed CSV write no longer leaves a new summary.

### Added

- Inner-shell ionisation channels in the electron transport (#273), library
  only (no CLI input yet; that is #156). `build_shell_channel_tables` builds
  the valence table and one `ShellChannelTable` per shell of a
  `ShellResolvedChannels` (single-pole Penn per channel, no exchange);
  `Transport::with_inner_shells` adds them to a layer. A channel is chosen
  by its inverse IMFP and the loss drawn from that channel's table; a shell
  event under the Kieft-Bosch model liberates an electron of energy
  `E_F + ω - B` and leaves `B - E_F` in the solid. New `ElectronTally::inner_shell`
  hook, `SummaryTally::inner_shell_events` / `inner_shell_loss_ev`, and
  `LayerMetadata::inner_shells` (serialized only when a layer has shells).
  `ShellChannelTable` files carry the shell, its binding energy and that
  energy's provenance, with their own format version
  (`SHELL_CHANNEL_FORMAT_VERSION`); the `CrossSectionTable` cache format is
  unchanged. Layers without shells give bit-identical results.
- `BcaTally::partner` (#250): a tally hook, a no-op by default, that reports
  every collision partner's impact parameter and the number of partners of
  its collision step, before the collision changes the particle. The
  amorphous model reports its hard and weak collisions, the crystal flight
  model every simultaneous partner. Results are unchanged.
- `FullPenn::diimfp_grid` and `electron::inelastic::DiimfpGrid` (#256): the
  full Penn DIIMFP tabulated once for the energies of one table. The `q`
  integrand of each loss node is stored as the panels of its adaptive
  quadrature and integrated over each energy's momentum window; between
  nodes the window integrals are interpolated in `ln p` against `ln ω`,
  with nodes added until the interpolation is within the model tolerance
  (relative, or relative to the row's mean density in `ln ω` where the
  DIIMFP is below it). Built in parallel, bit-identical on any thread count.
  `DiimfpGrid::unresolved_cells` and `unresolved_loss_range_ev` report the
  cells left failing the check and their loss range; the grid does not
  answer inside those cells, so the table builder uses the direct model
  there rather than an unchecked interpolation. Near such a cell a row's
  loss-density panels are not split below the grid's narrowest cell width
  (1e-4 in `ln ω`): the model switch at the cell edges is a step that the
  width-invariant split test would otherwise bisect to floating-point
  resolution and pad up to `MAX_NODES` direct evaluations. Rows away from
  unresolved cells, and models without a grid, are unchanged. New
  `table::loss_density_node_count` reports a row's node count. New example
  `penn_full_build_time` (build time and accuracy checks).
- Hydrogenic L-subshell ELFs (#245): `hydrogenic_shell_elf` and
  `hydrogenic_shell_elfs` now also build L1 (2s) and L2 / L3 (2p), from the
  bound-free cross sections of Karzas and Latter, Astrophys. J. Suppl. 6, 167
  (1961), eqs. (36) and (37), reduced for `n = 2` and checked against the
  Gaunt factors of that paper's Table 1. L2 and L3 use the same 2p formula,
  each with its own binding energy and occupancy. New public
  `hydrogenic_2s_oscillator_strength_density_per_ev` and
  `hydrogenic_2p_oscillator_strength_density_per_ev`. K-shell results are
  unchanged; the M shell and above are still an error (with a new message).
- Ion beam divergence in transport (#285): the optional `[beam.divergence]`
  input table (`model = "gaussian"` with `sigma_deg` per plane, or
  `"uniform-cone"` with `half_angle_deg`; widths in `[0, 10]` degrees; static
  ion runs only), the opt-in `Bca::with_divergence` engine setter,
  `Bca::primary_direction`, `Bca::divergence_metadata`, and `Beam`-level
  Python arguments `divergence_model` / `divergence_deg`. Each primary's
  direction is sampled once about the nominal direction on a dedicated
  segment of its own stream (word `2^65`) and conditioned on pointing into the
  target (bounded rejection, error `BcaError::BeamDivergence` on exhaustion).
  Without it results are unchanged. The summary gains
  `physics.beam_divergence` (no `format.version` bump). `Bca::history` and
  `history_in` now return the new `HistoryError` (wrapping `StoppingError`)
  and `BcaError` gains the `BeamDivergence` variant. An input capability, not
  a validated channeling prediction. Example:
  `examples/b_5keV_si_crystal_divergence.toml`.
- `geometry::CsgGeometry`, `geometry::Csg` and `geometry::Primitive`:
  constructive-solid-geometry targets of box, capped-cylinder and half-space
  primitives combined by nestable union, intersection and difference, one
  region per top-level solid, each tagged with a material. Rays are
  classified by span combination (Roth 1982), with the overlap,
  surface-ownership and tolerance rules of the mesh target, so a CSG box gives
  the same events as the equivalent `VoxelGrid` and `MeshGeometry` (tested).
  Every top-level solid must be bounded (a half-space only inside an
  intersection or difference that bounds it), else the new
  `GeometryError::CsgUnbounded`; a bad primitive or an empty operator is the
  new `GeometryError::CsgInvalid`. `GeometryError` is not `non_exhaustive`, so
  an exhaustive `match` on it must add the two variants. Not yet reachable
  from the CLI or TOML input (#194).
- Electron event-cap diagnostics (#268): how many secondary tracks the
  collision cap (`max_events`) cut off, and how many primary histories had
  the primary or any secondary cut off (each history counted once). New
  public `tally::EventCapCounts` (`secondary_tracks`, `affected_histories`)
  and field `ElectronReport::event_caps` (so `results.event_caps` in
  `electron_summary.json`); reports written before it read back with zeros.
  New public field `SummaryTally::secondaries_event_capped` in
  `electron::transport`. Both are new public struct fields, so code that
  builds `ElectronReport` or `SummaryTally` with a struct literal must add
  them (`SummaryTally` derives `Default`). `FateCounts` still counts one
  fate per primary, and the `event_cap_ev` energy accounting, the random
  number sequence and every other result are unchanged.
- Electron table-coverage diagnostics (#253): per layer and per channel
  (elastic, inelastic), how many rate evaluations of the transport read the
  cross-section table inside its energy grid (endpoints included) and how
  many used the constant continuation below or above it, with the grid
  bounds, primaries and secondaries together. New `ElectronTally::table_lookup`
  hook (no-op default, so existing tallies compile unchanged),
  `electron::transport::{TableChannel, GridCoverage}`,
  `Transport::layer_tables`, and `tally::table_coverage`
  (`TableCoverageTally`, `LayerTableCoverage`, `TableCoverageCounts`).
  `ElectronReport::table_coverage` (and so `results.table_coverage` in
  `electron_summary.json`) carries the counts; reports written before it
  read back with an empty list. `validation/oracles/run_electron.py` sums
  them over its batches into the lindhard summary. These are evaluation
  counts, not collision counts. Counting draws no random number and changes
  no result.
- Validation and run metadata for the local Oen-Robinson electronic loss
  (`ElectronicLoss::EquipartitionLsOr`) in the crystal flight model (#226).
  The model was already in place: each lattice partner within `p_max` takes
  the local loss at its closest approach. New public
  `CrystalMetadata::electronic_constants_unverified`, `true` for crystal runs
  that use the local loss while the Oen-Robinson constants are not verified
  against the 1976 paper. It follows the new constant
  `ion::stopping::oen_robinson::OR_CONSTANTS_UNVERIFIED`, and a `false`
  value is left out of the serialised metadata. The new tests are in
  `tests/crystal_electronic.rs` and measure `R = E_local / E_nonlocal`
  (Ar 20 keV into Si). Along <110> and <100>, R is 0.26 and 0.60 times the
  off-axis value. Along the fixed off-axis direction 30/17 the crystal R is
  7.6 % above the amorphous one at the same `p_max`, outside the 5 % that
  #226 asked for. #250 traced this to the directions, not to the counting:
  averaged over all beam directions the crystal R is 0.990 of the amorphous
  one (the rule of angular averages of Lindhard, Mat. Fys. Medd. Dan. Vid.
  Selsk. 34, no. 14 (1965), section 5), and the channeling directions lie
  below that average. The 5 % criterion is now asserted on the direction
  average, and the 30/17 quotient is a recorded value. The module docs no
  longer call impact-parameter-dependent stopping in crystals a later step.
  Results are unchanged.
- Cited electron band defaults (#115): `electron::boundary::BAND_DEFAULTS`,
  `band_defaults`, `BandDefaults`, `BandKind`, `CitedValue`, `BandFill`,
  `BandDefaults::complete` and `BandStructure::from_defaults`. The table
  cites the Si band gap (1.1 eV) and electron affinity (4.05 eV) and the
  SiO2 band gap (9 eV) to Robertson and Wallace (2015), by table or figure
  and page. Parameters with no source that could be opened (the work
  function and Fermi energy of Al, Cu, Au and W, the valence band width of
  Si and SiO2, and the affinity of SiO2) have no number. `complete` takes
  them from the caller with their provenance.
- Electron table cache (#168): `lindhard run --table-cache DIR` reads the
  elastic and inelastic cross-section tables from `DIR` when it holds them
  for exactly this run's physics, grid and executable, and stores the tables
  it builds otherwise, so a series of runs builds its tables once. Entries are
  named by the SHA-256 of a key document covering every input of the build;
  a file that does not match its key is refused, never silently reused.
  Outputs are bit-identical with and without the cache; the summary records
  each table's `source` (`built` or `cache`) and its cache file (`path`,
  `sha256`, `key_sha256`) under `physics.materials`.
  `validation/oracles/run_electron.py` uses it across batches.
- Table cache integrity (#249): each entry now also stores the SHA-256 of
  the table file (`<kind>-<hash>.sha256`), and a lookup refuses a table whose
  bytes do not match it (or whose hash file is missing), naming the file,
  before parsing it. The elastic key includes `electron.elastic.model`
  (`ElasticChoice::model`, `ElasticModelChoice::label`), and the key version
  is 2, so entries written before this change are rebuilt.
- Inner-shell ELFs built in the engine (#135): `electron::inelastic::shell_elf`
  builds the optical ELF of a K shell from Stobbe's hydrogenic
  photoionization formula (dV2022 eq. (2), edge and occupancy from a
  `SubshellBindingTable`), ready for `ShellResolvedChannels::new`. New public
  `hydrogenic_shell_elf`, `hydrogenic_shell_elfs`,
  `hydrogenic_k_oscillator_strength_density_per_ev` and `ShellElfGrid`.
  Other subshells are still caller-supplied.
- Radial PSF tally in the electron CLI (#165): `[electron.tally.psf]` (depth
  slab, log radial bins, `fits`, `normalization`) writes
  `electron_psf_profile.csv` and `electron_psf_parameters.csv` and adds
  `results.psf` and `files.psf_profile` / `files.psf_parameters` to
  `electron_summary.json`. A failing fit is reported in `results.psf.fit_errors`
  without stopping the run. `ResolvedElectron` gains `psf_fits` and
  `psf_normalization`; `ElectronTallySpec` gains `psf`. Without the table,
  output is unchanged.
- Sputter erosion in dynamic runs (#234): `[dynamic] erosion = true` removes
  sputtered atoms from the front of the target (slab 0 first) instead of the
  slab where they were displaced, and the surface recedes. New public fields
  `DynamicConfig::erosion`, `StepRecord::{recession_m, recession_total_m}`,
  `Particle::origin_layer` and `DynamicRun::recession_m()`; `recession_nm` in
  `dynamic_summary.json` totals (erosion on only). With erosion off, results
  are bit-identical to before.
- The Salvat et al. (1987) DHFS screening table (Table I, Z = 1..92, #130):
  `SalvatDhfs::for_element` now returns the published coefficients (two-term
  potentials for the asterisked rows, `SalvatDhfs::from_two_terms`), so
  `potential = "salvat-dhfs"` runs instead of failing.
- Release workflow (`.github/workflows/release.yml`): tag-triggered `lindhard`
  CLI archives for five targets, `SHA256SUMS`, a GitHub Release, and a
  crates.io publish gated on a protected environment; see "Releasing" in
  `CONTRIBUTING.md`.
- Thermal vibration in the crystal flight model (#181): `CrystalTarget::thermal`
  / `with_thermal(Thermal { temperature_k, debye_temperature_k,
  include_zero_point })` displaces every lattice site the particle meets by an
  uncorrelated Gaussian of the Debye amplitude (`ion::crystal::debye`), drawn
  once per encounter from a dedicated segment of the history's random stream.
  Static crystals (`thermal: None`) and amorphous runs are bit-identical to
  before; `T = 0` without the zero-point term reproduces the static lattice bit
  for bit. The target temperature and the per-species amplitudes are in the
  new run metadata, `Bca::crystal_metadata` (`CrystalMetadata`, serialisable).
  B 5 keV along <110> in Si: the tail beyond twice the amorphous Rp falls from
  0.876 (0 K) to 0.833 (300 K) and 0.758 (600 K). At 7 degrees tilt, 22 degrees
  twist the dRp ratio to amorphous is 1.42 (B 5 keV) and 1.81 (As 30 keV) at
  300 K: vibration does not remove that tail.
- `Bca::with_crystal` and `CrystalTarget` (`ion::bca::crystal`): the crystal
  flight model of the BCA engine (#180). The partner model is chosen per
  region (amorphous, or sites of a `Lattice` in a given `Orientation`);
  amorphous regions and runs are bit-identical to before (a test compares
  complete reports). Partners come from the lattice neighbour search, nearest
  by path distance, with simultaneous partners merged by momentum balance
  (criterion after DISPLATH, MIT, see `THIRD_PARTY_LICENSES.md`); the lattice
  is static and perfect, with a random translation per history. B 5 keV along
  <110> reaches 6.1 times the amorphous Rp; at 7 degrees tilt and 22 degrees
  twist a static lattice keeps a channeling tail (dRp 1.25 times the
  amorphous value for B 5 keV, 1.88 for As 30 keV; Rp within 10 %). Criterion
  bench `crystal_flight`: about 17 times slower than amorphous at 7/22.
- `[physics] tuning = "none" | "<set>"`: opt-in phenomenological
  surface-binding-energy multiplier sets for static single-element ion runs
  (`input::TuningSet`). `"none"` and omission are bit-for-bit unchanged;
  `summary.json` gains `physics.tuning` (set, version, original and effective
  `E_s`) only when a set is used (#80).
- First shipped tuning set, `es-sputter-ar-v1` (#80): per-element `E_s`
  multipliers for Ar sputtering of Si, Cu, Ag and Au, a phenomenological
  calibration (not a published model value) fitted to measured yields only, by
  a deterministic grid search on training sets with a held-out evaluation
  (`validation/experiments/run.py --fit-tuning`, record in
  `validation/experiments/tuning_results.json` and `docs/data-provenance.md`).
  A tuning set now also names the beam species and energy range it was fitted
  for: other beams and unlisted target elements are rejected, an energy
  outside the range or a tilted beam is warned about. Untuned results stay the
  primary level-3 comparison; default physics is unchanged. Held-out error
  falls for all four targets, but the Ag factor sits at the fit grid's lower
  bound and the factors absorb the whole yield deficit, whatever its cause.
- `geometry::MeshGeometry` and `geometry::TriMesh`: targets of closed triangle
  solids, each tagged with a material, loaded from STL (ASCII or binary) or
  OBJ with a watertightness check (every edge shared by two consistently
  oriented triangles, no degenerate triangles, positive volume; a bad mesh is
  rejected with an error naming an edge). Ray queries use a binned-SAH BVH and
  a watertight ray-triangle test; the target implements `Geometry` with
  documented overlap, surface-ownership and tolerance rules (#27).
- `ion::crystal::search`: `LatticeSearch`, the lattice sites within `p_max` of
  a path segment found by walking unit cells (no global atom list), each with
  its impact parameter, path distance, species and cell/basis index, ordered
  along the path with a documented closed boundary rule. Not yet used by the
  transport engine (#20).
- `ion::crystal::debye`: Debye-model RMS thermal displacement (1D and 3D, with
  the zero-point term) and per-atom Gaussian displacement sampling on the
  particle's RNG stream; Debye temperatures of Si, Ge, GaAs and 3C-SiC as
  cited defaults (#21).
- `lindhard::ion::crystal` (#19): the cubic crystal data model, not yet
  used by the transport engine. `Lattice` with diamond and zincblende
  constructors and cited presets (Si, Ge, GaAs, 3C-SiC; lattice constants
  with their temperature), Miller-index plane normals and number densities;
  `Orientation`, the lab-to-crystal rotation from tilt, twist and wafer
  rotation about an explicit in-plane reference direction; and `Divergence`,
  Gaussian or uniform-cone beam divergence drawn from the per-particle random
  stream. Conventions and a figure: `docs/crystal-orientation.md`.
- `ElectronReport.stopping_points.primaries` (`PrimaryStoppingPoints`): the
  depth and radial moments of the primaries alone that fell below the stopping
  threshold, i.e. the penetration depth of stopped primaries; written in the
  CLI's `electron_summary.json` (#150).
- Electron code-to-code comparison harness (#150):
  `validation/oracles/run_electron.py` runs lindhard, Nebula and Geant4
  MicroElec (through our own application in
  `validation/oracles/geant4_microelec/`) on 1, 5 and 20 keV electrons into Si
  and Cu, and `docs/validation.md` tabulates η, δ, the primary penetration
  depth and the 50 %-energy radius from the committed scalar summaries.
- Electron engine benchmarks (#152): Criterion benches
  (`cargo bench -p lindhard --bench electron`) for the transport rate at 1, 5
  and 20 keV into Si and Cu and, separately, the cross-section table build;
  the `electron_scaling` example (wall time, electrons/s and parallel
  efficiency on 1..N threads, with a tally-report digest check at every
  count); and `validation/oracles/bench_electron.py`, which runs them (and
  Nebula's CPU build, when configured) on the #150 matched problems and
  writes the scalar `bench-electron-*.json` summaries and the tables of
  `docs/benchmarks.md`.
- Documented gap for the PMMA optical ELF (#147): Ritsko et al. (1978) and the
  Henke tail could not be opened on 2026-10-07, so no dataset is committed;
  `docs/data-provenance.md` records what was tried and what is needed.
- Documented gap for the resist PSF validation (#151): no measured PSF of a
  stated stack and energy, and neither the PMMA nor the Si optical ELF, could
  be opened on 2026-10-08, so no comparison is run; `docs/validation.md`
  fixes the tolerance rule and limitations for when the inputs exist, and
  `docs/data-provenance.md` records the sources tried.
- Electron run mode in `lindhard-cli`: an input with an `[electron]` table
  (beam, transport cutoff and escape rule, Kieft-Bosch secondaries, step
  barrier, Mott elastic with optional exchange and correlation-polarization
  corrections, single-pole Penn, full Penn or Mermin-ELF inelastic, per-material
  optical ELF files and band, phonon and polaron parameters, tally grids) runs
  the electron transport with the full electron tally and writes
  `electron_summary.json` (the `ElectronReport`, every model choice and every
  data provenance, with file SHA-256s) and `electron_escape_spectra.csv`,
  `electron_deposition_{cylindrical,cartesian}.csv` and `electron_tables.csv`.
  Data without a provenance is refused. Schema types:
  `lindhard::input::electron`; `--histories` is an alias of `--ions`. Example:
  `examples/electron/e_10keV_si.toml` (synthetic ELF and band data). Library
  plumbing for it: `electron::inelastic::table::build_inelastic_table_for_model`
  (energy-loss tables for any `PennInelastic` model) and
  `electron::elastic::table::AtomicElastic::compute_corrected` /
  `ThomasFermiYukawa::yukawa` (elastic tables with the corrections).
- `tally::psf`: the radial point-spread function of a pencil beam. Setting
  `ElectronTallyConfig::psf` (a depth slab and log-spaced radial bins from
  nanometres to the backscatter range) makes `FullElectronTally` report a
  `RadialProfile` in `DepositionReport::psf`, with per-bin standard errors
  from per-history accumulation. `fit_psf` fits the normalised double
  Gaussian (Mao et al. (2025) eq. (1), after Chang (1975)) or triple
  Gaussian (Rosa Figueiro (2015) eq. (72)) by weighted least squares on the
  bin energies, with `E` fixed to the slab energy by default or free, and
  returns the parameters, covariance, reduced χ² and residuals. `PsfReport`
  exports the profile and fits as CSV and derives `Serialize`/`Deserialize`
  for callers to write JSON. The CLI does not expose it
  yet.
- Fluence-dependent targets: `ion::dynamic::DynamicRun`, a fluence stepping
  loop on `CompositionGrid` with fixed or adaptive steps (bounded relative
  composition change per step, with reject and retry), and `InventoryTally`,
  the adapter from transport events to inventory deltas. Each step covers a
  range of the global primary indices (`Bca::run_range`,
  `rng::run_particles_range`; `run` is unchanged and bit-identical), so a run
  is reproducible end to end and identical at any thread count. The CLI reads
  an optional `[dynamic]` table and writes `dynamic_summary.json`,
  `dynamic_steps.csv` (cumulative yields per step) and
  `dynamic_composition.csv` (slab profile per step); inputs without it run as
  before. `CompositionGrid::seed_energies` supplies the energies of elements
  that only enter the target during the run. Example:
  `examples/dynamic/as_1keV_si_film.toml`.
- `electron::inelastic::ShellResolvedChannels`: a valence channel and
  inner-shell ionization channels, each from its own caller-supplied optical
  ELF (one loss function per shell, as in de Vera et al., Int. J. Mol. Sci.
  23, 6121 (2022), eqs. (1)-(2) and (32)), with secondary energy `ω - B`,
  channel-resolved inverse mean free paths and DIIMFPs, and channel sampling.
  No partition of a single total ELF is applied.
- `electron::inelastic::ExchangeCorrection`: optional Born-Ochkur exchange for
  `SinglePolePenn` below a stated energy, off by default. The denominator is
  `T' - W` with `W = ω - B` the emitted energy (de Vera et al. (2022), eq.
  (32)) and the loss limit `ω <= (T' + B)/2`
  (`SinglePolePenn::diimfp_with_binding_per_m_ev`,
  `imfp_and_stopping_with_binding`, `born_ochkur_factor`); the plain model
  and the valence channel use `B = 0`.

### Changed

- **Electron results change for every `lindhard run` whose materials have a
  `band`** (#241). The inelastic table of a material with a band is now
  built on the band-bottom energy axis the transport reads it on, with the
  model's Fermi energy taken from the band: the row at `E` is the model's at
  `E - E_F`, kinematics on `E`, losses below `E - E_F`. `E_F` is the band's
  minimum excitation energy, the new `BandStructure::min_excitation_ev` (the
  Fermi energy of a metal, the conduction-band bottom `W_v + E_g` of an
  insulator). This is the convention of the table compiler of Nebula
  (`compile_full_imfp_icdf` and its caller, cstool commit `0c739eb`,
  BSD-3-Clause). Before, every table was on the model's own axis with
  `[electron.inelastic] fermi_energy_ev` (default 0), its losses reached `E`,
  and the Kieft-Bosch secondary model clamped those beyond `E - E_F`, leaving
  the primary at the Fermi level (#173 measured 78 to 82 % of the Al events
  5 to 20 eV above the Fermi level in the clamp). Tables of materials
  without a band, and so runs without bands, are unchanged bit for bit, as
  are their cache keys.
  - `[electron.inelastic] fermi_energy_ev` now applies only to materials
    without a band, and a nonzero value is refused if any material has one
    (it would count the Fermi energy twice). Migration: remove the key from
    inputs that give a `band`.
  - Library: new `electron::inelastic::table::EnergyAxis`
    (`ModelFermiLevel`, the default, and `BandBottom`), new public field
    `InelasticTableOptions::axis` with `InelasticTableOptions::with_axis`
    (code that builds `InelasticTableOptions` with a struct literal must add
    the field; `InelasticTableOptions::new` sets the default), and new
    `PennInelastic::fermi_energy_ev`. A table built with the default axis is
    the table built before. The transport is unchanged: it still clamps, and
    the clamp no longer acts on a band-bottom table above its first row with
    losses.
  - The diagnostic example `lindhard-cli/examples/inelastic_low_energy.rs`
    loses its `fermi-reference` mode (now what `lindhard run` builds) and
    gains `legacy`, the table as built before.
  - The δ(E) tables of `docs/validation.md` (secondary-electron yield) and
    the committed electron oracle summaries were produced before this
    change and have not been re-run (#287).
- Faster full Penn (`penn-full`) inelastic tables (#256). A table now reads
  its rows' DIIMFP from `FullPenn::diimfp_grid` instead of the nested
  integrals, and the plasmon term and the `ω_p` integral of the model are
  cheaper (a sign test instead of a plasmon root solve per bisection step;
  the knot intervals' error checked in groups, every knot kept). The Al
  table on the #169 grid (51 eV to 30 keV, 20 points per decade) builds in
  203 s on two threads; before, one row did not finish in 15 minutes. The
  model's IMFP and stopping power change by at most 3e-8 and 2e-7 (Al, 57
  energies); table rows change within the model tolerance (the grid
  interpolation, `docs/validation.md`), so full-Penn tables built before
  and after differ in their last digits. Tables of the single-pole and
  Mermin models are unchanged. The full-Penn table's `model` string now
  says that its DIIMFP comes from the loss grid. The Si table on the #168
  grid still does not finish (isolated slow points of the `ω_p` integral,
  older than this change; #298). On the band-bottom axis (#241) the loss
  grid is built for the rows above the Fermi level only, the rows that get
  losses.
- Crystal off-axis channeling tail (#225). The 7°/22° orientation is
  2.6-2.7° from a {100} and a {110} plane, so it is no longer held to the
  #180 random-direction bound "dRp within 15 % of amorphous". Its ignored
  tests are renamed `near_planar_7_22_static_lattice` and
  `near_planar_7_22_at_300_k` and are now recorded-value regression checks;
  tolerances are three seed-to-seed standard deviations.
  - A new ignored test, `off_axis_tail_grows_with_vibration_under_both_losses`,
    records that at 7°/22° and 30°/17° the tail beyond twice the amorphous
    Rp grows at 300 K. It does so under both `NonLocal` and
    `EquipartitionLsOr`, by 5-9 σ at 16000 ions.
  - The 30°/17° check at 300 K now also asserts the 90th percentile within
    10 % of amorphous, and bounds dRp at its recorded value plus 0.12.
  - The `ion::bca::crystal` module docs record the literature search. No
    measured profile at matched conditions was usable. One published
    simulation (Bratchenko et al. 2009) found the same sign of thermal
    feeding-in, so the effect is consistent with it but not validated
    against measurement.
  - No transport code changed.
- `lindhard run` with `electron.elastic.potential = "salvat-dhfs"` now
  applies the `exchange` and `correlation_polarization` corrections the input
  asks for, solved on the DHFS Poisson density (Salvat et al. 1987,
  Eq. (12)). Before, the DHFS elastic table was built without them and the
  settings were silently ignored; DHFS runs with either correction on change
  results (#169). The stand-in path is unchanged.
- Elastic cross sections change for potentials without a polarization tail
  (static `Yukawa`, `SalvatDhfs`, `SquareWell`, and `CorrectedPotential` with
  exchange only): the radial Dirac solver now starts the outward integration
  by the WKB criterion of #91 for every potential, instead of
  `r_t exp(-60/|kappa|)`. For Au (Salvat DHFS) the change is below `1e-11`
  relative in `sigma_el` and `sigma_tr1` up to 50 keV, and `2.4e-5` in
  `sigma_tr1` at 100 keV, where the old rule was off by up to `2.5e-4` rad at
  `|kappa| = 779`. Breaking: `ScreenedPotential::long_range()` is removed; it
  no longer selected anything (#131).
- Dynamic runs: `surface_nm` is now the cumulative surface recession (always 0
  with `erosion = false`). Breaking for struct-literal construction: new public
  fields on `Particle`, `DynamicConfig` and `StepRecord` (#234).
- Lint policy is declared once in `[workspace.lints]` (`unsafe_code = "deny"`) and
  every crate opts in with `[lints] workspace = true`, so `lindhard-cli` and
  `lindhard-py` no longer lack the guard the library has (#218, supersedes #184).
- `lindhard run` into a reused output directory now removes reserved optional
  files the current run does not produce (`ions.csv` without `tally.per_ion`;
  the electron deposition CSVs without their grid), so they cannot be mistaken
  for current output; unrelated files are left alone, a failed removal is an
  error naming the path, and the summary is written last (#185).
- Collision hot path (about 2.2 to 2.7 times the ions/s, see
  `docs/benchmarks.md`): the scattering angle is carried as `tan(theta/2)` and
  sines and cosines instead of angles, the Lindhard-Scharff coefficient is
  computed once per ion and target instead of at every flight, and
  `Bca::run` reuses its particle stack and scratch buffers across histories.
  Results are statistically unchanged but not bit-identical to 0.0.1 (the
  floating-point operations differ); they remain bit-identical across thread
  counts.
- Default elemental densities now follow the X-Ray Data Booklet, LBNL/PUB-490
  Rev. 3 (2009), Table 5-2 (from the CRC Handbook, 80th ed.). 37 values
  changed, by up to 5 % (for example Sr 2.64 to 2.54, Cs 1.93 to 1.873, Ag
  10.49 to 10.50 g/cm³), which changes the number density of a pure-element
  target that relies on the default. The full list is in
  `docs/data-provenance.md`. Si, Cu and Au are unchanged.

### Added

- Inelastic electron tables (`electron::inelastic::table`, #94):
  `build_inelastic_table` builds a `CrossSectionTable` with
  `SamplingAxis::InelasticEnergyLoss` (inverse mean free path and the inverse
  CDF of the energy loss) from the single-pole Penn DIIMFP on a caller-chosen
  energy grid, with `stopping_power_ev_per_m` recovering `S(E)` from it to
  1e-3; `MomentumTransferSampler` draws the momentum transfer given `(E, W)`
  and gives the deflection cosine.
- Optical ELF dataset for Si (`validation/data/optical/si_elf_yang2019.toml`,
  #125): the bulk ELF of Yang et al., Phys. Rev. B 100, 245209 (2019), from
  REELS, 0.5 to 199 eV, digitized from their Fig. 7 (script and second read in
  `validation/data/digitize/`). The sum-rule test now checks a nonconductor's
  perfect-screening sum against `1 - 1/eps1(0)` (Si: 0.9135 against 0.9143);
  Si N_eff is 7.66 of 14 because the table stops before the K shell, pinned.
  The TPP 2011 IMFP check gains Si (99.5 eV to 9.9 keV, within 3.6 %).
- Optical ELF datasets for Al and Cu (`validation/data/optical/`, Hagemann,
  Gudat and Kunz 1975, read through the CC0 refractiveindex.info
  transcription and checked against the scanned tables of DESY report
  SR-74/7), with a sum-rule test gated on the linear interpolation
  `OpticalElf::elf()` serves, integrated exactly. Both materials miss a 5 %
  check there (Al N_eff +13.8 % and P_eff +9.4 %, Cu N_eff +7.3 %); the
  failures are pinned and explained in the test as the chord overshoot of
  log-spaced knots, with power-law segments kept as a labelled comparison
  (Cu passes, Al N_eff +7.5 % against the authors' own +4.5 %) (#98 in
  part). Si is a documented gap.
- `electron::phonon`: opt-in insulator channels for the electron loop.
  `FrohlichPhonon` (Fröhlich LO-phonon emission and absorption inverse mean
  free paths and angular distribution, after Llacer and Garwin (1969), as
  restated by Ding et al., Sci. Technol. Adv. Mater. 22, 932 (2021)), with a
  cited SiO₂ preset, and `PolaronTrapping` (`C exp(-γE)`, after Ganachaud and
  Mokrani (1995); `C`, `γ` are caller inputs). `Transport::with_insulator_channels`
  switches them on per layer (all layers default to none); `LayerMetadata`
  gains `phonon` and `polaron`. New `Fate::PolaronTrapped`, `PhononEvent`,
  `ElectronTally::phonon` and `ElectronTally::polaron_trapped` hooks, and
  `SummaryTally` counters for phonon events and trapped energy.
- `electron::elastic` exchange and correlation-polarization corrections
  (#91), both off by default: `solve_corrected` and `CorrectedPotential` add
  Furness-McCarthy local exchange and the Salvat (2003) correlation-polarization
  potential (Perdew-Zunger LDA correlation joined to a Buckingham
  `-alpha_d/(2(r^2+d^2)^2)` tail, Seltzer's cutoff or a caller-set `b_pol^2`)
  to a static potential at one energy. The atomic electron density is an
  explicit input (`ElectronDensity`, implemented for `SalvatDhfs` and `Yukawa`
  by their Poisson densities) and the dipole polarizability is a caller input
  with a required citation. `ElasticResult` gains a `corrections` field
  (`CorrectionMetadata`) recording which corrections, model identifiers and
  inputs were used; `solve` reports both off. `ScreenedPotential` gains
  `long_range()` (default `false`); with it the solver places the start of the
  outward integration by a WKB criterion. With both corrections off results
  are bit-identical to before.
- `geometry::Geometry`, the engine-facing target trait, implemented by `Stack`
  and the new `geometry::VoxelGrid` (regular 3D grid of material indices,
  periodic or vacuum boundaries per axis, exact 3D DDA traversal that
  truncates the free path at a material change; `Geometry::flight` also
  reports the region a flight without an event ends in, so the particle's
  region follows it across same-material voxel faces). `Bca::new` now takes any
  `&dyn Geometry` (a `&Stack` still works) and `Bca::with_entry_point` picks
  the incident point. `Face::Side` and `EnergyBudget::lateral` account for
  escapes through the lateral faces of a voxel grid; `Face` moved to
  `geometry` (still re-exported from `ion::bca`) and gained that variant;
  `kinematics::refract_out_normal` applies the surface barrier along an
  arbitrary face normal. `Particle::layer` is now the region index (layer or
  flat voxel index). Stack runs are unchanged.
- Secondary electrons and surface barriers in `electron::transport` (Kieft
  and Bosch; equations from Verduin's thesis, TU Delft 2017; kinematics and
  step ported from Nebula and the band model from cstool, BSD-3). New
  `TransportConfig` fields, all off by default so existing runs are unchanged
  bit for bit: `secondaries` (`SecondaryModel::KieftBosch`: one secondary of
  energy `E_F + W - B` per inelastic event, Ivanchenko direction model,
  pushed on the history's stack and transported in full), `boundary`
  (`BoundaryModel::StepBarrier`: quantum-mechanical transmission, refraction
  and reflection at every face by the inner-potential step, with the primary
  incident from vacuum) and `cutoff_reference` (cutoff from the band bottom
  or the vacuum level). Per-layer band parameters (`electron::boundary::
  BandStructure`: Fermi energy and work function, or valence band width, band
  gap and electron affinity, with a required provenance) go through
  `Transport::with_band_structures`; no material defaults are built in (see
  `docs/data-provenance.md`). New `ElectronTally` hooks with no-op defaults
  (`secondary`, `barrier`, `reflected`, `begin_secondary`, `end_secondary`),
  new `SummaryTally` counters, and the models and band parameters in
  `RunMetadata`. `FullElectronTally` implements the new hooks and balances
  with either model on: with a secondary model only the part of an inelastic
  loss that no secondary carries away is deposited at the event, and the
  budget gains a `fermi_sea_ev` source (the energy conduction electrons
  already had) and a `barrier_ev` term (the sum of `-ΔU` over face
  transmissions), so `incident + fermi_sea = deposited + escaped + trapped +
  barrier`. Stops are told from trapped electrons by the per-layer stopping
  threshold (new `Transport::stopping_thresholds_ev`, recorded as
  `ElectronTallyMetadata::stopping_threshold_ev`), and an electron cut off by
  the event cap is booked from its own last state, primary or secondary.
  Reports of runs with both models off are unchanged apart from the three new
  fields.
- `electron::transport`: the event-by-event electron loop (Kieft and Bosch,
  J. Phys. D 41, 215310 (2008)) over layered stacks on elastic and inelastic
  `CrossSectionTable`s. `Transport::run` goes through `rng::run_particles`
  (bit-identical at any thread count) and returns the tally with `RunMetadata`
  recording the energy cutoff and escape rule; per-event results go through
  the `ElectronTally` hook trait.
- `tally::electron`: `FullElectronTally`, an `ElectronTally` that records
  energy deposition per layer and on optional Cartesian and cylindrical
  (r-z) grids, backscatter (eta) and secondary (delta) yields split at a
  configurable energy (50 eV by default, after Chen et al., Sci. Rep. 12,
  18201 (2022); the split and its source are recorded in the report
  metadata), escape energy and polar-angle spectra per face,
  energy-weighted generation-volume moments, and the energy balance
  (deposited + escaped + trapped = incident). Its plain-data
  `ElectronReport` is serializable and bit-identical at any thread count.
  Example: `cargo run -p lindhard --example electron_tally`.
- `electron::elastic`: radial Dirac partial-wave solver for a screened central
  potential at one energy (phase shifts, differential cross section, Sherman
  function, `sigma_el`, `sigma_tr1`), with `Yukawa`, `SquareWell` and
  `SalvatDhfs` potentials. The Salvat et al. (1987) coefficient table is not
  yet in the tree (`SalvatDhfs::for_element` returns an error).
- `electron::elastic::table`: elastic `CrossSectionTable`s for any `Material`
  on a caller-chosen energy grid (default 10 eV to 50 keV, 20 points per
  decade). Elements combine by independent-atom additivity; the DCS is held
  in exact Legendre form and the shared probability grid is refined until the
  `sigma_tr1` recovered from the stored inverse CDF matches the solver's (to
  1e-3 at every grid energy, tested). The potential comes from a
  `PotentialSource`; `ThomasFermiYukawa` is a stand-in until the DHFS table is
  available, and the table's model and provenance strings say so.
- `SubshellBindingTable::eadl2017()`: the committed Z=1..92 subshell binding
  energies and occupancies (EADL as distributed in EPICS2017; D. E. Cullen,
  IAEA-NDS-224 Rev. 1), embedded in the library with its credit notice, under
  the operator ruling linked from `docs/data-provenance.md`. Its coverage test
  runs in the default test run.
- `electron::inelastic`: the single-pole Penn model (`SinglePolePenn`), built
  from a user-supplied `OpticalElf`, with the momentum-dependent loss
  function, the DIIMFP, the inelastic mean free path and the stopping power
  (nonrelativistic kinematics, stated integration tolerance); a sum-rule
  report for any optical ELF (`SumRuleReport`: f-sum `N_eff(W)`, `P_eff`,
  mean excitation energy); nonrelativistic Bethe stopping
  (`inelastic::bethe`); and an analytic Drude-Lorentz ELF with closed-form
  sum rules (`DrudeLorentz`) as a synthetic fixture. Example
  `penn_sum_rules` and bench `penn`.
- `electron::inelastic::FullPenn`: the full Penn algorithm, the optical ELF
  expanded over Lindhard dielectric functions (`LindhardGas`, with its
  plasmon and limiting-form tests), with the loss function, DIIMFP, IMFP and
  stopping power; it agrees with the single pole to 1e-3 at 10 to 50 keV on
  the synthetic Drude fixture and differs below 1 keV as documented in the
  module docs. `PennAlgorithm` / `PennInelastic` select either per material
  and give the identity string to record as a table's `model`. Example
  `penn_full_vs_single_pole`, bench group
  `penn_single_pole_vs_full_8_energies`. Penn (1987) and Lindhard (1954)
  were not opened; see `docs/data-provenance.md`.
- `electron::inelastic::MerminPenn`, `MerminGas`, `fit_mermin_oscillators`:
  the Mermin relaxation-time dielectric function (formula as in de Vera et
  al., Int. J. Mol. Sci. 23, 6121 (2022), eq. (4); γ -> 0 reproduces the
  Lindhard function to 1e-8), a deterministic MELF-GOS fit of oscillators
  (non-negative amplitudes, energies, widths) to an `OpticalElf` that reports
  the parameters, residuals, f-sum and `P_eff` (`MerminFit`), and an inelastic
  model (loss function, DIIMFP, IMFP, stopping power) on the fit.
  `PennAlgorithm::Mermin` (label `mermin-melf`) selects it beside the single
  pole and the full Penn algorithm, `PennInelastic::try_new` / `mermin` build
  it, and `model_identity` records the choice and the fitted parameters.
  Mermin (1970) was not opened and the threshold step of the MELF is not
  implemented; see `docs/data-provenance.md`.
- `Bca::history_in` and `HistoryBuffers`: `Bca::history` with caller-owned
  working memory, which allocates nothing in steady state.
- `ScatteringTable::half_angle_tan`, `kinematics::rotate_sc` and
  `kinematics::lab_projectile_sc`.
- `ElectronicStopping::sqrt_energy_coefficient` (defaulted, so existing
  implementations are unaffected), implemented by `LindhardScharff`.
- Packaging metadata for both crates: `include` allowlists, crate-local
  README, LICENSE and CHANGELOG in each archive, and a CI `package` job that
  runs `cargo package` and `cargo publish --dry-run`.
