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

## Output

### `summary.json`

| Key | Content |
|---|---|
| `format` | `{"name": "lindhard-summary", "version": 1}` |
| `software` | Crate `version` and `git_describe` of the binary |
| `input` | The input as run: defaults filled in, CLI overrides applied, `run.threads` removed. Deserializes back to the same `lindhard::input::Input`, so a run can be reproduced from its own header |
| `physics.models` | Every model in use: `role`, `name`, `citation` |
| `physics.engine` | Cutoffs, free path, electronic-loss mode, seed, chunk size as passed to the engine |
| `physics.scattering_table` | Angle-table grid and its measured interpolation error |
| `physics.target` | Each layer: extent (nm; `back_nm` is `null` for a substrate), atom density, and the fully resolved material (every `E_d`, `E_b`, `E_s`) |
| `results.histories` | Primaries run |
| `results.primaries` | `stopped`, `backscattered`, `transmitted`; mean and standard deviation of the stopped-primary depth, nm |
| `results.recoils` | Atoms `displaced`, `sputtered` (left through the front face), `transmitted` |
| `results.yields` | The above per incident ion |
| `results.energy_budget_ev_per_ion` | Where the incident energy went, per ion, and the largest per-history relative bookkeeping residual |
| `files` | Names of the other files written (`null` if not written) |
| `run` | `threads`, `table_build_s`, `transport_s`, `ions_per_s` |

### `depth_profile.csv`

`depth_lo_nm,depth_hi_nm,stopped_primaries,fraction_per_nm`: stopped primaries
per depth bin, and that count per incident ion per nm. The last row's upper
edge is `inf` (overflow) and its density is empty.

### `ions.csv`

`index,fate,x_nm,y_nm,z_nm,energy_ev,dir_x,dir_y,dir_z,layer`: the final
state of every primary, by history index. `fate` is `stopped`,
`backscattered` or `transmitted`; for escaped primaries the position is on the
face and the energy and direction are outside the target. `x` is depth.

## Reproducibility

Everything in `summary.json` except the trailing `run` object, and both CSV
files, is a function of the input and the binary only: byte-identical for the
same input and seed at any thread count. `lindhard-cli/tests/examples.rs`
checks this on 1 and 4 threads. Floats are written in shortest round-trip
form.

## Compatibility and extension

Consumers must ignore keys they do not know. New results are added as new
keys, never by changing existing ones:

- New tallies (range moments and fits, damage, sputtering by species and
  energy, lateral profiles) become new objects under `results`, for example
  `results.range`, `results.damage`, `results.sputtering`, with their settings
  as new keys under `[tally]` and their profiles as new CSV files listed
  under `files`.
- New model choices become new values of the existing `[physics]` keys, or
  new keys with defaults, so existing inputs keep their meaning.
- `format.version` is bumped only when an existing key is removed or changes
  meaning.
