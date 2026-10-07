# Changelog

All notable changes to `lindhard` and `lindhard-cli` are recorded here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the crates use [Semantic Versioning](https://semver.org/) (both crates share
one version).

## [Unreleased]

### Added

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
- `Bca::history_in` and `HistoryBuffers`: `Bca::history` with caller-owned
  working memory, which allocates nothing in steady state.
- `ScatteringTable::half_angle_tan`, `kinematics::rotate_sc` and
  `kinematics::lab_projectile_sc`.
- `ElectronicStopping::sqrt_energy_coefficient` (defaulted, so existing
  implementations are unaffected), implemented by `LindhardScharff`.
- Packaging metadata for both crates: `include` allowlists, crate-local
  README, LICENSE and CHANGELOG in each archive, and a CI `package` job that
  runs `cargo package` and `cargo publish --dry-run`.
