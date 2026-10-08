# Changelog

All notable changes to `lindhard` and `lindhard-cli` are recorded here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the crates use [Semantic Versioning](https://semver.org/) (both crates share
one version).

## [Unreleased]

### Added

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
