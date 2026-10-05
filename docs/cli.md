# The `lindhard` command

```text
lindhard check input.toml                  # parse and validate, no transport
lindhard run input.toml --out dir/         # run, write dir/summary.json and CSVs
lindhard run input.toml --out dir/ --ions 200 --seed 7 --threads 4
lindhard --version                         # crate version and git describe
```

`--ions` and `--seed` override `run.ions` and `run.seed`; the override is what
the output echoes. `--threads` overrides `run.threads` and never changes the
results. Invalid input exits non-zero with a message naming the offending key
(`target.layers[0].thickness_nm: -5 nm must be finite and positive`), and
writes no output. Warnings (for example a beam energy outside the stopping
model's advisory range) go to stderr and do not stop the run.

Examples: [`../examples/`](../examples/).

## Input (TOML)

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
| `primary_cutoff_ev` | required | The primary stops below this energy |
| `recoil_cutoff_ev` | required | Recoils stop below this; keep it below the smallest `E_s` |
| `follow_recoils` | `true` | Full cascades |
| `primary_surface_binding_ev` | 0 | Surface barrier for the beam species |

`[physics.energies.<symbol>]` sets `e_d_ev`, `e_b_ev` and/or `e_s_ev` for that
element in every layer that contains it, after (so overriding) the
material's own values. An element with no tabulated default and no value set
is an error naming the layer and the key to set.

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

## Output

### `summary.json`

| Key | Content |
|---|---|
| `format` | `{"name": "lindhard-summary", "version": 1}` |
| `software` | Crate `version` and `git_describe` of the binary |
| `input` | The input as run: defaults filled in, CLI overrides applied, `run.threads` removed. Deserializes back to the same `lindhard::input::Input`, so a run can be reproduced from its own header |
| `physics.models` | Every model in use: `role`, `name`, `citation` |
| `physics.stopping_tables` | Only with `[stopping]`: path, SHA-256, provenance and range of each user table (see `[stopping]`) |
| `physics.engine` | Cutoffs, free path, electronic-loss mode, seed, chunk size as passed to the engine |
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

## Reproducibility

Everything in `summary.json` except the trailing `run` object, and every CSV
file, is a function of the input and the binary only: byte-identical for the
same input and seed at any thread count. `lindhard-cli/tests/examples.rs`
checks this on 1 and 4 threads. Floats are written in shortest round-trip
form.

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
