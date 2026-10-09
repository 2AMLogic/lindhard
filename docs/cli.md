# The `lindhard` command

```text
lindhard check input.toml                  # parse and validate, no transport
lindhard run input.toml --out dir/         # run, write dir/summary.json and CSVs
                                           # (dir/electron_summary.json for [electron])
lindhard run input.toml --out dir/ --ions 200 --seed 7 --threads 4
lindhard run electron.toml --out dir/ --table-cache cache/   # reuse built electron tables
lindhard --version                         # crate version and git describe
```

An input with a top-level `[electron]` table is an electron run (see
"Electron runs" below); every other input is an ion run.

`--ions` and `--seed` override `run.ions` and `run.seed` (`--ions`, or its
alias `--histories`, overrides `run.histories` of an electron run); the override is what
the output echoes. `--threads` overrides `run.threads` and never changes the
results. `--table-cache DIR` (electron runs only) reads the cross-section
tables from `DIR` when it holds them for exactly this run's physics, grid and
build, and stores the tables it builds otherwise; it never changes the results
(see "Cross-section table cache" below). Invalid input exits non-zero with a message naming the offending key
(`target.layers[0].thickness_nm: -5 nm must be finite and positive`), and
writes no output. Warnings (for example a beam energy outside the stopping
model's advisory range) go to stderr and do not stop the run.

Examples: [`../examples/`](../examples/).

## Input (TOML)

<!-- Rendered in the book (book/src/guide/input.md and output.md) between the anchors in this file. -->
<!-- ANCHOR: input -->
Units are in the key names: energies in eV (`_ev`), lengths in nm (`_nm`),
angles in degrees (`_deg`), densities in g/cm³. Every table rejects unknown
keys. The schema types are `lindhard::input` (shared with future front ends).

### `[beam]`

| Key | Default | Meaning |
|---|---|---|
| `ion` | required | Element symbol of the projectile (case-sensitive, `"As"`) |
| `mass_amu` | standard atomic weight | Projectile mass, u |
| `energy_ev` | required | Incident energy, eV |
| `tilt_deg` | 0 | Polar angle from the surface normal, `[0, 90)` |
| `azimuth_deg` | 0 | Azimuth of the incidence plane |

### `[materials.<name>]`

Named materials, in the form of `lindhard::material::MaterialSpec`:

```toml
[materials.SiO2]
density_g_cm3 = 2.2                    # required for compounds
elements = [
  { symbol = "Si", atom_fraction = 1.0 },
  { symbol = "O", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },
]
```

Each element takes `symbol` or `z`, `atom_fraction` or `mass_fraction`
(relative weights, normalised), and optional `e_d_ev`, `e_b_ev`, `e_s_ev`.
Defaults for the energies come from the element table where one is
tabulated; see `lindhard/src/material.rs`.

### `[target]`

```toml
[target]
substrate = "Si"            # optional semi-infinite substrate

[[target.layers]]           # finite layers, front to back
material = "SiO2"
thickness_nm = 10.0
```

A `material` (in a layer or as the substrate) is either a key of
`[materials]`, an element symbol (the pure element at its tabulated density),
or an inline material table with the same keys as `[materials.<name>]`. A
`[materials]` key wins over an element symbol of the same name. Without a
substrate the target has a back face and particles can be transmitted.

### `[physics]`

| Key | Default | Choices |
|---|---|---|
| `potential` | `"zbl"` | `zbl`, `kr-c`, `moliere`, `lenz-jensen` |
| `screening_length` | paired with the potential | `universal`, `firsov`, `lindhard` |
| `stopping` | `"lindhard-scharff"` | `lindhard-scharff`, `bethe-bloch`, `equipartition-ls-or` |
| `free_path` | `"constant"` | `constant`, `energy-dependent` |
| `min_cm_angle_deg` | none | Required with, and only with, `energy-dependent` |
| `weak_collisions` | 0 | `0` to `3`: weak collisions beyond `p_max` per collision step (Moller and Eckstein, IPP 9/64 (1988)); `constant` free path only. See the `ion::bca` docs, "Weak collisions" |
| `primary_cutoff_ev` | required | The primary stops below this energy |
| `recoil_cutoff_ev` | required | Recoils stop below this; keep it below the smallest `E_s` |
| `follow_recoils` | `true` | Full cascades |
| `primary_surface_binding_ev` | 0 | Surface barrier for the beam species |
| `tuning` | `"none"` | Opt-in phenomenological calibration: the name of a versioned factor set (see below) |

`[physics.energies.<symbol>]` sets `e_d_ev`, `e_b_ev` and/or `e_s_ev` for that
element in every layer that contains it, after (so overriding) the
material's own values. An element with no tabulated default and no value set
is an error naming the layer and the key to set.

**`tuning` (phenomenological calibration, not a published model choice).**
`"none"` or omission leaves the physics and the echoed input exactly as
without the key. A named set scales the surface binding energy `E_s` by a
per-element factor fitted to measured data. Rules: the factor multiplies the
*resolved* `E_s` (an explicit `[physics.energies]` or material value if given,
else the elemental default), once per layer, after overrides; the global
element table and the collision algorithm are untouched. The pilot supports
static ion runs on single-element layers only, with a beam species the set
was fitted for and target elements the set lists; compounds, `[dynamic]`
targets, other beams, unlisted elements and unknown set names are rejected
with a `physics.tuning` error. A beam energy outside the set's fitted range,
or a tilted beam, runs with a warning (an extrapolation). `summary.json` then
has `physics.tuning` with the set, its version and provenance, and per layer
the original `E_s`, the factor and the effective `E_s`. Tuned results must be
reported next to, never in place of, untuned ones.

Shipped sets (fit record and held-out scores: `docs/data-provenance.md`,
"Tuning factor sets"; a new version ships under a new name):

| Set | Beam | Elements | Fitted energies | Fitted under |
|---|---|---|---|---|
| `es-sputter-ar-v1` | Ar | Si, Cu, Ag, Au | 196 to 10020 eV, normal incidence | the matched level-3 sputter settings (`docs/validation.md`, section 3) |

The factors are a calibration of yields under those settings; other settings
(potential, `E_d`, cutoffs, weak collisions) were not part of the fit, and
the set does not claim better accuracy for them.

### `[stopping]`

Optional. Supplies user stopping tables for the electronic stopping of
particular (ion, target element) pairs.

```toml
[stopping]
tables = ["tables/b_in_si.toml", "tables/p_in_si.toml"]
```

| Key | Default | Meaning |
|---|---|---|
| `tables` | none | Paths of table files, one per pair |

Absent, the input means what it always meant (`format.version` is unchanged
and nothing is echoed).

**Table file.** The `lindhard::ion::stopping::table::StoppingTable` format:

```toml
provenance = "Author, Journal vol, page (year), Table N"   # required
ion_z = 5
ion_mass_amu = 11.0093            # optional; default: standard atomic weight
target_z = 14
energy_ev = [1.0e2, 1.0e3, 1.0e4]                 # strictly increasing, eV
stopping_ev_1e15_cm2 = [10.0, 30.0, 60.0]         # eV 1e-15 cm^2 per atom
```

Interpolation is piecewise linear in ln S versus ln E. `provenance` is
mandatory (data without an origin is not admitted). Do not use SRIM- or
ICRU-derived tables in anything committed to a repository; a table is the
user's own data and its terms are the user's concern.

**Paths.** A relative path resolves against the directory of the input file
(not the current directory). The echoed input keeps the path as written.

**Composition with `[physics] stopping`.** A table replaces the
`[physics] stopping` model for exactly the pair it declares (`ion_z`,
`target_z`), including recoils of that species when `follow_recoils` is on.
Every other pair uses the `[physics] stopping` model. A pair with a table is
never silently served by the model: a query outside the table's energy range,
or for a different ion mass, is an error that stops the run. Nothing is
extrapolated. Declare each pair once. A table for a pair that cannot occur
in the run (including a table for a target element's recoils when
`follow_recoils = false`) is accepted with a warning that it is unused.

**Recoil species.** Tables are keyed by (`ion_z`, `target_z`) only, so with
`follow_recoils = true` a table whose `ion_z` is a target element also serves
every recoil of that element. Recoils carry the standard atomic weight and are
followed down to `physics.recoil_cutoff_ev`, so such a table is checked up
front against both: its `ion_mass_amu` must be the standard weight, and it
must start at or below `recoil_cutoff_ev`; either failure is an error. A
consequence is that an isotopic self-ion beam (e.g. `beam.mass_amu = 27.9769`
for Si into Si) cannot take a table for its own pair while recoils are
followed: drop the table, use the standard weight, or set
`follow_recoils = false`. A recoil-species table that ends below the largest
energy the beam can transfer to that element warns. Tables cannot be combined with
`stopping = "equipartition-ls-or"` (that mode carries its own
Lindhard-Scharff/Oen-Robinson loss and would ignore them).

**Errors** name the field (`stopping.tables[0]`): an unknown key in
`[stopping]`, a missing or unreadable file, invalid table contents (including
a missing `provenance`), a duplicate pair, a table whose ion mass differs from
the beam ion's, a table whose range does not contain the beam energy, and the
recoil-species checks above. A
table that starts above `physics.primary_cutoff_ev` warns, because the run
fails if a projectile slows below it.

**Provenance in the output.** Tables are user data, so the run records them.
`summary.json` has `physics.stopping_tables`, one entry per table: `path` (as
written), `resolved_path` (absolute where possible), `sha256` of the file's
bytes, the table's `provenance` string, `ion_z`, `ion_mass_amu`, `target_z`
and the energy range. `physics.models` lists each as `user-table` with the path
and provenance as its source. The key is absent without `[stopping]`.

### `[dynamic]` (optional)

Makes the run fluence-dependent: the target composition is updated as the
fluence builds up (sputter erosion, build-up of implanted atoms). Without the
table nothing changes: the run, its output files and its bytes are those of a
static run. The model, its conventions and its limits are in
`lindhard::ion::dynamic` and the book chapter on dynamic composition.

`run.ions` is the number of histories of the **whole** run and
`fluence_cm2` the fluence they represent, so each ion stands for
`fluence_cm2 / run.ions` ions/cm². The ions are delivered in steps; after each
step the grid is updated from that step's events (a recoil is subtracted where
it is created and added where it stops, a stopped beam ion is added, an
escaped atom is a loss; the substrate is an immutable reservoir) and relaxed.
Primary `i` of the run always uses the random stream `(seed, i)`, whatever the
step sizes and thread count.

| Key | Default | Meaning |
|---|---|---|
| `fluence_cm2` | required | Total fluence of the run, ions/cm² |
| `ions_per_step` | required | Ions per step; with `max_change`, the largest and first step |
| `max_change` | absent (fixed steps) | Adaptive steps: largest relative composition change of a slab per step, the largest absolute change in atoms/m² of one element in one slab, divided by that slab's atoms/m². A larger step is discarded and retried from the same first ion with fewer ions, which consumes no ions of the run; the step doubles again after a step below half the bound |
| `min_ions_per_step` | 1 | Adaptive only: the smallest step. At this size a step is accepted whatever its change, and a removal beyond what a slab holds is capped at what it holds (`clamped` column). A fixed run that removes more than a slab holds fails, naming the slab: use smaller steps or `max_change` |
| `slab_nm` | one slab per layer | Split each finite layer into equal slabs at most this thick; the composition is tracked per slab |
| `relaxation` | `"ideal-mixing"` | How thickness follows inventory: `"ideal-mixing"` (additive atomic volumes) or `"fixed-number-density"` |
| `number_density_cm3` | none | Total atom density, atoms/cm³; required with `"fixed-number-density"` |
| `atomic_volume_nm3.<Sym>` | elemental solid volume from the element table | Atomic volume, nm³/atom, per element (ideal mixing). Required for an element with no tabulated solid density (a gas) |
| `energies.<Sym>` | `[physics.energies.<Sym>]`, then element defaults | `e_d_ev`, `e_b_ev`, `e_s_ev` of an element that enters the target during the run (the beam species, for example). Elements already in a layer keep that layer's energies |
| `erosion` | `false` | Sputter erosion: sputtered atoms are removed from the front of the target (slab 0 first, then deeper slabs) instead of from the slab where they were displaced, and the surface recedes. Must be a boolean |

With `erosion = false` the front surface stays at `x = 0`: swelling moves the
interior interfaces and the back face of the slabs, not the front surface (the
`surface_nm` column is that fixed frame, always 0). With `erosion = true` the
lost thickness is removed from the front and the grid is re-anchored so the
current surface is again `x = 0`; `surface_nm` is then the cumulative recession
`R` in nm and a depth `x` in the output is `x + R` in the original frame. The
recession of a step is the volume of the removed atoms per area under the
chosen `relaxation` (`sum Z removed_Z v_Z`, or `sum removed / n`), and removal
equals the sputtered counts per element, so no atom is created or lost. If a
step sputters more of an element than the slabs hold, the excess is not
removed and the element is counted in the `clamped` column. With a substrate,
atoms sputtered from the substrate remove nothing from the slabs, so the
recession falls short of `Y F / n` once the film is thin.
`dynamic_summary.json` totals gain `recession_nm` only with erosion on. Depths
in the output are measured from the front surface of that step. A dynamic run
needs `E_d` for every element that can occur, including the beam species. The Python bindings run static inputs only.

### `[run]`

| Key | Default | Meaning |
|---|---|---|
| `ions` | required | Number of primary histories |
| `seed` | required | Run seed |
| `threads` | all cores | Worker threads; does not affect results and is not echoed |

### `[tally]`

| Key | Default | Meaning |
|---|---|---|
| `depth_bin_nm` | 1 | Bin width of the stopped-primary depth profile |
| `depth_bins` | 1000 | Number of bins; the last also collects everything deeper |
| `per_ion` | `true` | Write `ions.csv` |
| `lateral_bin_nm` | 1 | Bin width of the lateral and radial profiles of stopped primaries |
| `lateral_bins` | 100 | Bins per side of the beam axis: `y` and `z` span `[-lateral_bins * lateral_bin_nm, +lateral_bins * lateral_bin_nm)` nm, the radial distance `[0, lateral_bins * lateral_bin_nm)` nm |
| `escape_energy_max_ev` | beam energy | Upper edge of the escape-energy spectra, eV (from 0) |
| `escape_energy_bins` | 100 | Escape-energy bins |
| `escape_polar_bins` | 30 | Polar-angle bins over `[0, 90)` degrees from the outward surface normal |
| `dual_pearson` | `false` | Also fit a dual-Pearson profile to the depth histogram |

The depth grid (`depth_bin_nm`, `depth_bins`, from the front face) is shared by
the range histogram, the dual-Pearson fit and the defect profiles. Particles
outside any grid are counted in explicit underflow and overflow entries, never
dropped. Every count and per-ion value is for the same incident ions.

### Electron runs (`[electron]`)

An input with an `[electron]` table runs the low-energy electron engine
(`lindhard::electron::transport`) instead of the ion BCA, with the full
electron tally (`lindhard::tally::FullElectronTally`). It exposes what the
library does and adds no physics; every choice and every data provenance is
written to the output, so a result can be reproduced from its own header. The
schema types are `lindhard::input::electron`. `[materials]` and `[target]`
are the ion run's tables (above); `[beam]`, `[physics]`, `[stopping]`,
`[tally]` and `[dynamic]` are not accepted. Example:
[`../examples/electron/e_10keV_si.toml`](../examples/electron/e_10keV_si.toml).

```toml
[electron.beam]
energy_ev = 10000.0

[electron.transport]
cutoff_ev = 1.0
cutoff_reference = "vacuum-level"
secondaries = "kieft-bosch"
boundary = "step-barrier"

[electron.elastic]
potential = "thomas-fermi-yukawa"

[electron.inelastic]
model = "penn-single-pole"

[electron.materials.Si]
optical_elf = "si_elf.toml"
band = { kind = "insulator", valence_band_width_ev = 10.0, band_gap_ev = 2.0, affinity_ev = 3.0, provenance = "..." }

[target]
substrate = "Si"

[run]
histories = 1000
seed = 1
```

**`[electron.beam]`**

| Key | Default | Meaning |
|---|---|---|
| `energy_ev` | required | Kinetic energy of the primaries, eV: the vacuum energy with `boundary = "step-barrier"`, the energy inside the first layer otherwise |
| `tilt_deg` | 0 | Polar angle from the surface normal, `[0, 90)` |
| `azimuth_deg` | 0 | Azimuth of the incidence plane |

Primaries start on the front face at `y = z = 0` (just outside it with the
step barrier).

**`[electron.transport]`** (`lindhard::electron::transport::TransportConfig`)

| Key | Default | Choices |
|---|---|---|
| `cutoff_ev` | required | An electron stops below this energy |
| `cutoff_reference` | `"band-bottom"` | `band-bottom`; `vacuum-level` (the threshold is `U + cutoff`: electrons that can no longer leave are not followed) |
| `escape_rule` | `"both-faces"` | `both-faces`; `front-only` (the back face absorbs) |
| `max_events` | 10000000 | Collision and reflection cap per electron |
| `secondaries` | `"off"` | `off`; `kieft-bosch` (Kieft and Bosch 2008) |
| `instantaneous_momentum`, `momentum_conservation` | `true` | Options of `kieft-bosch`; an error with `off` |
| `boundary` | `"transparent"` | `transparent`; `step-barrier` (inner-potential step with quantum transmission and refraction) |
| `quantum_transmission`, `refraction` | `true` | Options of `step-barrier`; an error with `transparent` |

**`[electron.elastic]`**

| Key | Default | Choices |
|---|---|---|
| `model` | `"mott"` | `mott`: Mott cross sections from radial-Dirac partial waves, independent-atom additivity (`electron::elastic::table`) |
| `potential` | required | `thomas-fermi-yukawa`: the Thomas-Fermi Yukawa **stand-in**; `salvat-dhfs`: the Salvat et al. (1987) DHFS potentials (Table I coefficients, Z = 1..92; `data-provenance.md`) |
| `exchange` | `false` | Furness-McCarthy exchange correction |
| `correlation_polarization` | absent (off) | A table: `polarizability.<Sym> = { bohr3 = ..., source = "..." }` for every target element (the source is required), optional `b_pol_squared` (absent: Seltzer's rule, which needs every table energy above 50 eV) and `outer_radius_bohr` (50) |

The corrections are solved per grid energy with the chosen potential's own
Poisson density: the stand-in's, or the DHFS density of Salvat et al. (1987)
Eq. (12) (`AtomicElastic::compute_corrected`); the elastic table's `model` and
`provenance` strings name them and every polarizability with its source.

**`[electron.inelastic]`**

| Key | Default | Choices |
|---|---|---|
| `model` | `"penn-single-pole"` | `penn-single-pole`, `penn-full`, `mermin-melf` (`electron::inelastic::PennAlgorithm`). The full Penn and Mermin models integrate numerically and build tables far more slowly. The single-pole model's mean free path is much longer than the other two below about 30 eV (Al: up to 23 times), which inflates the secondary yield; see `electron::inelastic::penn`, "Low energies" (#173) |
| `fermi_energy_ev` | 0 | Fermi energy of the model, eV. It is not the band's: the transport reads table rows at the electron's energy above the band bottom, so setting it to the band's Fermi energy counts that energy twice; see `electron::transport`, "Energy reference of the inelastic table" (#173) |

**`[electron.tables]`**: one log-spaced energy grid shared by the elastic and
inelastic tables of every material.

| Key | Default | Meaning |
|---|---|---|
| `min_energy_ev` | 10 | Lowest grid energy, eV |
| `max_energy_ev` | the beam energy (plus the largest inner potential with the step barrier) | Highest grid energy, eV |
| `points_per_decade` | 20 | Minimum points per decade |

The transport holds the rates of the first and last rows beyond the grid;
a grid that starts above the lowest stopping threshold or ends below the
largest possible energy warns.

**`[electron.materials.<name>]`**: the electron data of each material the
target uses, keyed by the name the target gives it (a `[materials]` key or an
element symbol; an inline target material is an error in an electron run).
Every name the target uses needs an entry; an unused entry warns.

| Key | Default | Meaning |
|---|---|---|
| `optical_elf` | required | Path of an optical ELF file, relative to the input file's directory, in the `lindhard::electron::data::OpticalElf` TOML form (`material`, `provenance`, `energy_ev`, `elf`). It is read with that type's loader, so **a file without a provenance is refused**, as is any invalid table |
| `band` | none | Band parameters, required with `kieft-bosch`, `step-barrier` or `vacuum-level`: `{ kind = "metal", fermi_ev, work_function_ev, provenance }`, `{ kind = "insulator", valence_band_width_ev, band_gap_ev, affinity_ev, provenance }` or `{ kind = "free-electron-metal", valence_electrons_per_atom, work_function_ev, provenance }` (`lindhard::electron::boundary::BandStructure`; a blank provenance is refused) |
| `phonon` | none (off) | Fröhlich LO-phonon channel, polar insulators only: `{ hbar_omega_ev, eps_static, eps_high_frequency, temperature_k, provenance }`, or `{ preset = "sio2-63mev" \| "sio2-153mev", temperature_k }` (the library's cited SiO₂ values) |
| `polaron` | none (off) | Polaron trapping `C exp(-γE)`: `{ c_per_nm, gamma_per_ev, provenance }` |

No optical or band data of any real material is committed
([`data-provenance.md`](data-provenance.md)); the data files are the user's,
and their terms are the user's concern. Subshell binding-energy tables
(`SubshellBindingTable`) have no key yet: the transport loop does not use
inner-shell channels, so there is nothing to feed them to (see "extending"
below).

**`[electron.tally]`** (`lindhard::tally::ElectronTallyConfig`)

| Key | Default | Meaning |
|---|---|---|
| `se_bse_split_ev` | 50 | Escaping electrons below it are slow (secondary), at or above it fast (backscattered) |
| `escape_energy_max_ev` | beam energy | Upper edge of the escape-energy spectra (from 0), eV |
| `escape_energy_bins` | 100 | Escape-energy bins |
| `escape_polar_max_deg` | 90 | Upper edge of the polar-angle spectra (from 0, at most 180), degrees from the outward normal |
| `escape_polar_bins` | 18 | Polar-angle bins |
| `cartesian` | none | Deposition grid `{ x, y, z }`, each `{ lo_nm, hi_nm, bins }` (`x` is depth) |
| `cylindrical` | none | Deposition grid `{ r, depth }` about the beam axis, each `{ lo_nm, hi_nm, bins }` (`r.lo_nm >= 0`) |

**`[run]`** of an electron run: `histories` (required), `seed` (required) and
`threads` (all cores; not echoed, never changes results).

<!-- ANCHOR_END: input -->

## Output

<!-- ANCHOR: output -->
### `summary.json`

| Key | Content |
|---|---|
| `format` | `{"name": "lindhard-summary", "version": 1}` |
| `software` | Crate `version` and `git_describe` of the binary |
| `input` | The input as run: defaults filled in, CLI overrides applied, `run.threads` removed. Deserializes back to the same `lindhard::input::Input`, so a run can be reproduced from its own header |
| `physics.models` | Every model in use: `role`, `name`, `citation` |
| `physics.stopping_tables` | Only with `[stopping]`: path, SHA-256, provenance and range of each user table (see `[stopping]`) |
| `physics.engine` | Cutoffs, free path, weak collisions, electronic-loss mode, seed, chunk size as passed to the engine |
| `physics.scattering_table` | Angle-table grid and its measured interpolation error |
| `physics.target` | Each layer: extent (nm; `back_nm` is `null` for a substrate), atom density, and the fully resolved material (every `E_d`, `E_b`, `E_s`) |
| `results.histories` | Primaries run |
| `results.primaries` | `stopped`, `backscattered`, `transmitted`; mean and standard deviation of the stopped-primary depth, nm (`stopped_depth_mean_nm`, `stopped_depth_std_nm`; the same numbers as `results.range.depth.mean_nm` and `std_dev_nm`, the projected range and straggle, kept under their original keys) |
| `results.recoils` | Atoms `displaced`, `sputtered` (left through the front face), `transmitted` |
| `results.yields` | The above per incident ion |
| `results.energy_budget_ev_per_ion` | Where the incident energy went, per ion, and the largest per-history relative bookkeeping residual |
| `results.range` | Where the beam particles came to rest, lengths in nm. `depth`: `n`, `mean_nm` (projected range `Rp`), `std_dev_nm` (straggle), `skewness`, `kurtosis` (`beta`, Gaussian 3) and the standard error of each (`null` below two stopped primaries). `pearson_iv`: the Pearson IV density with those moments (`m`, `nu`, `a_nm`, `lambda_nm`), or `null` with the reason in `pearson_iv_error`. `dual_pearson` (or `dual_pearson_error`): only with `tally.dual_pearson = true`; head fraction, head and tail components, chi-square of the fit and of the single Pearson IV. `lateral_y`, `lateral_z`, `radial`: moments of the lateral positions. `layers`: `stopped` and `depth` moments by the layer where the particle stopped |
| `results.damage` | `nrt`: Norgett-Robinson-Torrens and Kinchin-Pease displacement estimates from the primary knock-on atom damage energies (`pka_count`, `pka_energy_ev`, `damage_energy_ev`, `nrt_displacements`, `kinchin_pease_displacements`). `cascade`: defects counted event by event in the simulated cascades (`displacements`, `replacements`, `vacancies`, `interstitials`). The two are different quantities and are reported separately; see `lindhard::ion::damage` for the conventions, including the approximation used for compounds. Totals over all ions, with `per_ion` and a `layers` breakdown |
| `results.sputtering` | `yield_per_ion` and `by_element`: target atoms leaving the front face, with `count`, `per_ion` and `mean_energy_ev` for each element |
| `results.escapes` | `backscatter_coefficient`, `transmission_coefficient`, `energy_reflection_coefficient` (energy carried out of the front face by the beam particles, as a fraction of the incident energy), and `species`: every species (beam first) through the `front` and `back` face, with `count`, `per_ion` and `mean_energy_ev` |
| `files` | Names of the other files written (`null` if not written) |
| `run` | `threads`, `table_build_s`, `transport_s`, `ions_per_s` |

### `depth_profile.csv`

`depth_lo_nm,depth_hi_nm,stopped_primaries,fraction_per_nm`: stopped primaries
per depth bin, and that count per incident ion per nm. The last row's upper
edge is `inf` (overflow) and its density is empty.

### `lateral_profile.csv`

`quantity,lo_nm,hi_nm,count,per_ion_per_nm`: stopped primaries by `y`, `z` and
radial distance (`quantity` is `y`, `z` or `radial`), for the grid set by
`lateral_bin_nm` and `lateral_bins`. Each quantity ends with an `underflow`
row (`lo_nm = -inf`) and an `overflow` row (`hi_nm = inf`) with empty density.

### `damage_profile.csv`

`depth_lo_nm,depth_hi_nm,vacancies,interstitials,replacements,` then the same
three per incident ion per nm: the cascade defect counts by depth, on the
depth grid. `vacancies` is displacements minus replacements. The last row,
with `depth_hi_nm = inf`, holds everything deeper than the grid (empty
densities). The totals equal `results.damage.cascade`.

### `escape_spectra.csv`

`species_z,symbol,beam,face,spectrum,lo,hi,count,per_ion_per_unit`: the energy
(`spectrum = energy_ev`, `lo`/`hi` in eV) and polar-angle (`polar_deg`, degrees
from the outward normal) spectra of every species leaving each face (`front`
or `back`). Each spectrum ends with `-inf` and `inf` rows for entries outside
the grid, with empty densities. The density is per incident ion per eV or per
degree.

### `ions.csv`

`index,fate,x_nm,y_nm,z_nm,energy_ev,dir_x,dir_y,dir_z,layer`: the final
state of every primary, by history index. `fate` is `stopped`,
`backscattered` or `transmitted`; for escaped primaries the position is on the
face and the energy and direction are outside the target. `x` is depth.

### Dynamic runs (`[dynamic]`)

A run with a `[dynamic]` table writes `dynamic_summary.json`,
`dynamic_steps.csv` and `dynamic_composition.csv` instead of the static files.

`dynamic_summary.json` (format `lindhard-dynamic-summary`): the echoed
`input` (with `dynamic`), `physics`-model list `models`, the `species`,
`totals` (ions, steps, rejected attempts, fluence, slab count and thickness
before and after, yields per ion), `files` and the timing `run` object.

`dynamic_steps.csv`: one row per accepted step; step 0 is the initial target.
`step`, `first_index` (global index of the step's first ion), `ions`,
`ions_done`, `fluence_cm2` (delivered so far), `attempts` (more than 1 if the
adaptive bound rejected the step), `max_change`, `clamped`, `removed_slabs`
(slabs that emptied), `n_slabs`, `surface_nm` (cumulative surface recession
in nm; always 0 with `erosion = false`), `thickness_nm` (total of the finite
slabs), then cumulative counts since the start: `cum_backscattered`,
`cum_transmitted`, `cum_stopped_in_target`, `cum_stopped_in_substrate`, `cum_sputtered`, `sputter_yield` (atoms per ion so
far), `cum_recoils_transmitted` and `cum_sputtered_<Sym>` per element.

`dynamic_composition.csv`: the slab profile after every step (step 0 is the
initial target), one row per step and slab: `step`, `slab` (0 is the front),
`front_nm`, `back_nm`, `thickness_nm`, then per element
`atoms_per_cm2_<Sym>` and `fraction_<Sym>` (atom fraction). Slabs that
emptied are gone from later steps.

### Electron runs: `electron_summary.json` and `electron_*.csv`

An electron run writes these instead of the ion files.

`electron_summary.json` (format `{"name": "lindhard-electron-summary",
"version": 1}`):

| Key | Content |
|---|---|
| `software` | As for an ion run |
| `input` | The input as run: defaults filled in (including `tables.max_energy_ev`, `tally.escape_energy_max_ev` and the secondary and barrier options), CLI overrides applied, `run.threads` removed. Deserializes to `lindhard::input::electron::ElectronInput` |
| `physics.models` | Every model in use: `role`, `name`, `citation` (transport loop, elastic model and potential, corrections, inelastic model, secondaries, barrier, phonon and polaron channels, SE/BSE split) |
| `physics.transport` | The engine's `RunMetadata`: cutoff and its reference, escape rule, event cap, secondary and boundary models, seed, histories, chunk size, the primary, and per layer its extent (m), the `model` and `provenance` strings of both tables, the band parameters, phonon and polaron channels with their provenance |
| `physics.target` | Each layer: extent (nm), atom density and the resolved material |
| `physics.materials` | Each material: the ELF file (`path`, `resolved_path`, `sha256`, its `material` and `provenance`, energy range and point count), `band`, `phonon`, `polaron`, and for `elastic_table` and `inelastic_table` their `model`, `material`, `provenance`, cache `format_version`, energy range and grid sizes, `source` (`"built"` or `"cache"`) and `cache` (`null` without `--table-cache`, else the table file's `path`, `sha256` and `key_sha256`) |
| `results` | The `ElectronReport` (`lindhard::tally::ElectronReport`), lengths in m and energies in eV, summed over all histories unless named per primary: `histories`, `metadata` (split and its source, cutoff, stopping thresholds, tally settings), `fates` of the primaries, `budget` (the energy balance and its `relative_imbalance`; deposits are measured from the band bottom, so with secondaries in a layer with a Fermi energy `deposited_ev` includes the Fermi-sea energy of liberated conduction electrons and can exceed the energy imparted, which is `incident_ev - escaped_ev = deposited_ev + trapped_ev + barrier_ev - fermi_sea_ev - phonon_absorbed_ev`), `yields` (`backscatter_eta`, `secondary_delta`, `total_sigma`, transmitted), `front` and `back` (counts, energies, slow and fast classes), `deposition` (`per_layer_ev`; for each grid its binning, `inside_ev` and `outside_ev`), `generation_volume`, `stopping_points` (all electrons that fell below the stopping threshold, and under `primaries` the primaries alone: the penetration depth of stopped primaries), `table_coverage` (see below). The histograms and grid cells are in the CSV files, not here |
| `files` | Names of the CSV files (`null` if not written) |
| `run` | `threads`, `table_build_s`, `transport_s`, `histories_per_s` |

`results.table_coverage` is a numerical diagnostic: one entry per layer
(`layer`), and for its `elastic` and `inelastic` table the grid bounds
`energy_min_ev` and `energy_max_ev` and the counts `below` (`E <
energy_min_ev`: the first row's rate and distribution were used), `within`
(both bounds included: interpolated, or a row read exactly) and `above` (`E >
energy_max_ev`: the last row's were used), over primaries and secondaries.
They count **rate evaluations, not collisions**: the transport evaluates both
tables of the electron's layer once before every free flight, including
flights cut short at a layer face, the flight after a face reflection and
flights with zero total rate, so the totals exceed the number of elastic and
inelastic events, and a layer's elastic and inelastic totals are equal. They
are not fractions of path length or of deposited energy either. Nonzero
`below` or `above` counts say that part of the transport used the constant
continuation of a table beyond its grid (`[electron.tables]`); they do not say
how much that changed the result. Summaries written before the key existed
lack it.

`electron_escape_spectra.csv`: `face,spectrum,class,lo,hi,count,per_primary_per_unit`.
For each face (`front`, `back`): the energy spectrum of all escaping electrons
(`spectrum = energy_ev`, `class = all`, eV) and the polar-angle spectrum of
each class (`polar_deg`, `slow` or `fast`, degrees from the outward normal).
Each spectrum ends with `-inf` and `inf` rows for entries outside the grid,
with empty densities. The density is per primary per eV or per degree.

`electron_deposition_cylindrical.csv` (with `tally.cylindrical`):
`ir,ix,r_lo_nm,r_hi_nm,depth_lo_nm,depth_hi_nm,energy_ev,ev_per_primary_per_nm3`,
one row per cell. `electron_deposition_cartesian.csv` (with
`tally.cartesian`): `ix,iy,iz,x_lo_nm,x_hi_nm,y_lo_nm,y_hi_nm,z_lo_nm,z_hi_nm,energy_ev,ev_per_primary_per_nm3`.
Energy deposited outside a grid is `outside_ev` in the summary. Like `budget.deposited_ev`, the cell energies are
measured from the band bottom and, with secondaries, include Fermi-sea energy
the beam did not supply.

`electron_tables.csv`: `material,energy_ev,elastic_inverse_mfp_per_nm,inelastic_inverse_mfp_per_nm,inelastic_mean_loss_ev,inelastic_stopping_ev_per_nm`,
the tables the run used, per material and grid energy (the stopping power is
`λ⁻¹ ⟨W⟩` of the stored loss distribution).

### Cross-section table cache (`--table-cache`)

Building the elastic and inelastic tables is the slow part of a short
electron run (tens of seconds with `penn-single-pole`, far longer with
`penn-full`; see `docs/validation.md`, "Electron oracles"). With
`--table-cache DIR`, `lindhard run` looks each table up in `DIR` (created if
missing) and builds and stores only the ones it does not find, so a series of
runs that differ only in seed, history count, tallies or transport settings
builds its tables once.

Each entry is two files, named by the SHA-256 of a **key document**:
`<kind>-<sha256>.toml`, the table in the versioned cache form of
`lindhard::electron::data::CrossSectionTable`, read back with that type's
loader, and `<kind>-<sha256>.key.json`, the key itself. The key spells out
every input the table depends on:

- the key schema version and the table cache `format_version`;
- the build: crate version, `git describe --always --dirty`, and the SHA-256
  of the running `lindhard` executable, so any rebuild that changes the code
  (an uncommitted edit included) misses, and a rebuilt binary never reuses an
  older binary's table;
- the table kind and the exact energy grid (`electron.tables`, after the
  defaults are filled in);
- the material's name, composition and density;
- elastic: the potential, the exchange and correlation-polarization
  corrections with all their inputs, the starting probability grid and the
  refinement tolerance;
- inelastic: the model (`electron.inelastic.model`), its Fermi energy, the
  SHA-256 and provenance of the optical ELF file, and the material's band
  parameters.

Every `f64` is written in shortest round-trip form, so a change in the last
bit of any number is a different key. A lookup must find the stored key
equal, byte for byte, to the run's own (a mismatch under the same hash means
the file was edited, and is an error naming the differing field); the table
file is then loaded and validated by the library, and its axis and energy grid
are checked against the run. A file found under the run's key that fails any
of these checks is an error, never a silent rebuild. A table of another cache
`format_version` can never be found, since the version is part of the key.
Writes go to a temporary file renamed into place, the table before its key.

A cached table is the built one bit for bit (the TOML cache form round-trips
every `f64`), so outputs do not depend on whether the tables were built or
read, nor on the thread count. Only `physics.materials.*_table.source` and
`.cache`, and the timings in `run`, differ. `lindhard-cli/tests/examples.rs`
checks a building run, a storing run and reading runs on 1 and 8 threads
against each other. Remove the directory to reclaim space; nothing else
prunes it.

## Reusing an output directory

`--out` may name an existing directory; the run overwrites the files it
writes and creates the directory if needed. The CLI also owns the reserved
optional file names of the run's mode. After a successful run, an optional
file the run did not produce is removed if present: `ions.csv` (without
`tally.per_ion`), and `electron_deposition_cartesian.csv` or
`electron_deposition_cylindrical.csv` (without the matching deposition grid).
A missing file is not an error; a failed removal is, and names the path. The
summary is written last and lists only files that exist. Other files in the
directory are never touched, and no cleanup happens between ion, electron and
dynamic runs. Do not keep your own data under a reserved name.

## Reproducibility

Everything in `summary.json` except the trailing `run` object, and every CSV
file, is a function of the input and the binary only: byte-identical for the
same input and seed at any thread count (for a dynamic run: the same three
files, with `run` the only thread-dependent part of `dynamic_summary.json`).
For an electron run the same holds for `electron_summary.json` (apart from
`run`) and every `electron_*.csv`: the tables are built bit-identically on any
thread count, and histories run in chunks of a fixed size (16, recorded as
`physics.transport.chunk_size`) merged in chunk order.
`lindhard-cli/tests/examples.rs` checks this on 1 and 4 threads (1, 2 and 8
for the dynamic example; 1 and 4 for the electron example, and with
`--table-cache` 1 and 8). Floats are
written in shortest round-trip form.

## Compatibility and extension

Consumers must ignore keys they do not know. New results are added as new
keys, never by changing existing ones:

- New tallies become new objects under `results`, with their settings as new
  keys under `[tally]` and their profiles as new CSV files listed under
  `files`. `results.range`, `results.damage`, `results.sputtering` and
  `results.escapes` were added this way, without a version bump.
- New model choices become new values of the existing `[physics]` keys, or
  new keys with defaults, so existing inputs keep their meaning.
- `format.version` is bumped only when an existing key is removed or changes
  meaning.

The electron schema and its output follow the same rules, with
`lindhard-electron-summary` counting its versions separately:

- A new library model becomes a new value of an existing key
  (`electron.inelastic.model`, `electron.elastic.potential`,
  `electron.transport.secondaries`, `electron.transport.boundary`, a `band`
  `kind`, a `phonon` `preset`) and a new `physics.models` entry; existing
  values keep their meaning and defaults never change.
- New per-material data (for example a subshell binding-energy table, once
  inner-shell channels reach the transport loop, or a precomputed
  cross-section cache) becomes a new optional key of
  `[electron.materials.<name>]` that names a file. It must be read with the
  `lindhard::electron::data` loader of its type, so data without a
  provenance is refused, and recorded under `physics.materials` with its
  path, SHA-256 and provenance.
- The cross-section table cache is the one exception to the rule above, by
  decision (#168): it is a command-line flag (`--table-cache DIR`), not an
  input key, because it never changes a result (as with `--threads`, the
  input and its echo stay the same whether or not tables are reused), and a
  content-addressed directory keyed on every input of the build cannot name
  a stale file the way a hand-written path can. Its tables are still read
  with the `lindhard::electron::data` loader and echoed under
  `physics.materials` with their path, SHA-256 and provenance.
- A new tally becomes a key under `[electron.tally]`, an object under
  `results` and, for profiles, a new CSV file listed under `files`.
- Every table keeps `deny_unknown_fields`, and every default is echoed, so
  an input written today still means the same thing.
<!-- ANCHOR_END: output -->
