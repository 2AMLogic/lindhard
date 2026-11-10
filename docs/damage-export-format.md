# Damage export format

One text file that carries a damage map from either engine: positions, units,
per-channel energies, optional per-species defect counts, and the provenance
of the run. It is the voxel-tally CSV of [`cli.md`](cli.md#voxel_tallycsv-and-voxel_tallyvtk)
(same columns, same row order) under a banner and a header that identifies the
run.

**Status.** The format is for exchange between this engine and the tools that
read its maps. It is not a stability promise beyond its `format_version`: a
change to the header keys, the columns or their meaning bumps the version, and
a reader must refuse a version it does not know.

## Layout

```
# lindhard damage export v1
# format_version = 1
# engine = ion
# lindhard_version = 0.0.1
# seed = 9
# input_sha256 = <64 hex digits>
# units = length=m; energy=eV; counts=1
# origin_m = 0e0,-1e-8,-1e-8
# spacing_m = 4e-9,1e-8,1e-8
# dims = 10,2,2
# histories = 100
# outside = <implanted>,<displacements>,<electronic_ev>,<lattice_ev>,<rest_ev>
# species = 0:B,1:Si,2:O
# outside_defects = <3 integers per species>
i,j,k,x_m,y_m,z_m,implanted,displacements,electronic_ev,lattice_ev,rest_ev,deposited_ev,vacancies_0_B,...
```

Header lines are `# key = value`; the first line is the banner. Keys:

| Key | Meaning |
|---|---|
| `format_version` | Integer, currently `1`. |
| `engine` | `ion` or `electron`. |
| `lindhard_version` | Version of the writer. |
| `seed` | Root seed of the run. |
| `input_sha256` | SHA-256 (hex) of the resolved input echo, serialized as compact JSON. The ion CLI and Python bindings compute it from the echoed input of `summary.json`. |
| `units` | The fixed string `length=m; energy=eV; counts=1`. |
| `origin_m`, `spacing_m`, `dims` | Grid: corner of voxel `[0,0,0]`, voxel edges (m) and voxel counts. Axis `x` is depth. |
| `histories` | Primary histories tallied. |
| `outside` | Sums of the events outside the grid (same five fields as the columns). |
| `species`, `outside_defects` | Present only with per-species defects (ion runs): species index and element symbol, and the defects outside the grid (vacancies, interstitials, replacements per species). |

## Columns

One row per occupied voxel, in increasing flat index `i + nx (j + ny k)`
(x fastest). Unlisted voxels are empty.

| Column | Unit | Meaning |
|---|---|---|
| `i`, `j`, `k` | - | Voxel indices. |
| `x_m`, `y_m`, `z_m` | m | Voxel centre. Derived from the grid; a reader may ignore them. |
| `implanted` | count | Beam atoms that came to rest in the voxel (zero for electrons). |
| `displacements` | count | Target atoms displaced (zero for electrons). |
| `electronic_ev` | eV | Energy deposited into electronic excitation: for ions the electronic loss, for electrons the **ionisation** channel (inelastic events). |
| `lattice_ev` | eV | Energy left in the lattice: for ions the subthreshold and binding energy, for electrons LO-phonon emission and polaron trapping. |
| `rest_ev` | eV | Kinetic energy of a particle that stopped below its cutoff. |
| `deposited_ev` | eV | The sum of the three energy channels. |
| `vacancies_<n>_<sym>`, `interstitials_<n>_<sym>`, `replacements_<n>_<sym>` | count | Per species `n` (0 is the beam, then target elements in order of first appearance). Ion runs only. |

Floats are written as the shortest decimal text that parses back to the same
`f64` bits, so reading a written file returns the report unchanged. The
electron engine has no species axis: its files have no `species` header and no
defect columns. An empty tally is a header and a column row only.

## Reading and writing

- Rust: `VoxelReport::write_export_csv` / `to_export_csv` take an
  `ExportProvenance`; `VoxelReport::read_export_csv` returns the report and
  the provenance, and refuses a missing banner, version or provenance key.
  `FullElectronTally::voxel_report` gives the electron map to write.
- CLI: `lindhard run` writes `damage_export.csv` next to the other outputs
  when the input has a `[tally.voxel]` table (ion runs) or an
  `[electron.tally.voxel]` table (electron runs, `engine = electron`).
- Python: `RunResult.damage_export()`, `ElectronRunResult.damage_export()` and `lindhard.read_damage_export(path)`
  return the same dict of NumPy arrays ([`python.md`](python.md#damage-export)).

Not in this format: any cylindrical electron grid, and the NRT model damage
estimate (a model value, not a tally; see `cli.md`).
