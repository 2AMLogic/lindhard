# Voxel tallies and VTK export

The 1D tallies of [Output files](output.md) report results against depth. For
a target that is not laterally uniform (a voxel grid, a triangle mesh or a CSG
solid), the library's `lindhard::tally::VoxelTally` bins the same events in
3D. It is a library tally: attach it to a run through
`Bca::run(|| VoxelTally::new(grid))` and read it with `VoxelTally::report()`.
The command-line driver writes it when the input has a `[tally.voxel]` table
([CLI reference](https://github.com/2AMLogic/lindhard/blob/main/docs/cli.md)):
`voxel_tally.csv` and `voxel_tally.vtk` in the output directory.

## What is tallied

The tally uses its own regular grid (`VoxelTallyGrid`: origin, voxel edge
lengths and voxel counts), independent of the target geometry. For a voxel
grid target, `VoxelTallyGrid::from_voxel_grid` gives the target's own voxels.
Each voxel owns the half-open box `[o + i h, o + (i + 1) h)` on every axis, and
`x` is the depth.

| Column | Meaning |
|---|---|
| `implanted` | Beam particles that came to rest in the voxel |
| `displacements` | Target atoms displaced in the voxel (recoils created, transfer above `E_d`) |
| `electronic_ev` | Electronic energy loss, nonlocal and local |
| `lattice_ev` | Nuclear energy left in the lattice: subthreshold and weak transfers, and the binding energy of displaced atoms |
| `rest_ev` | Kinetic energy of particles when they stop below their cutoff |
| `deposited_ev` | `electronic_ev + lattice_ev + rest_ev` |

* **Per-element defects.** A tally built with
  `VoxelTally::with_species(&bca.species_z())` (the driver does this) also
  counts, for each species index of the run, the `vacancies_<i>_<symbol>` at
  the displacement site, the `interstitials_<i>_<symbol>` where the displaced
  atom came to rest without filling a site, and the
  `replacements_<i>_<symbol>`, with the definitions of
  [Displacement damage](../models/damage.md#cascade-counts) and the same
  replacement rule as the 1D tally. A replacement is credited to the species
  of the atom it replaced, so `vacancies + replacements = displacements` per
  species and voxel. Beam atoms are `implanted`, not interstitials, and a
  displaced atom that leaves the target leaves a vacancy and no interstitial.
  Summed over voxels, `outside` and species the counts equal the run's
  cascade totals exactly (integer equality, `lindhard/tests/voxel_tally.rs`).
  With `follow_recoils = false` there are no replacements. The aggregate
  `displacements` column is unchanged.
* **Damage** is the simulated cascade displacement count, the voxel version of
  the cascade counts in [Displacement damage](../models/damage.md). The
  NRT / Kinchin-Pease model estimate is not tallied per voxel, because it needs
  the Lindhard partition with the mean `Z` and `A` of each material; the two
  are never mixed.
* **Deposited energy** is split by channel so that each can be checked against
  the engine's energy budget. The work done against the surface barrier and
  the energy carried out of the target are not deposited in any voxel and are
  not tallied.
* The nonlocal electronic loss of a free flight is spread uniformly along the
  straight segment, as the engine defines it, so a segment crossing voxel
  faces is split in proportion to its length in each voxel. All other events
  go to the voxel that contains them. Events outside the grid are summed
  separately (`outside` in the report and the CSV header).

## Sparse storage and determinism

Only voxels that received an event are stored (a `BTreeMap` on the flat voxel
index `i + nx (j + ny k)`), so a grid of 10¹² voxels costs only its occupied
voxels, and iteration, merging and export always follow the same order.
Counts are integers and energies are `f64` sums merged in chunk order, so for a
fixed chunk size the report is bit-identical at any thread count. The
conservation tests check that the totals over all voxels equal the engine's
own 1D tallies of the same run for a laterally uniform target
(`lindhard/tests/voxel_tally.rs`).

## CSV export

`VoxelReport::write_csv` writes a `# key = value` header that carries the grid,
then one row per occupied voxel in flat-index order, with the voxel centre in
metres and the energies in eV. Floats are written in the shortest form that
reads back to the same bits, and `VoxelReport::read_csv` returns the report
unchanged.

The format is version 1 (the banner `# lindhard voxel tally v1`) for a tally
without species and version 2 with them: v2 adds the header lines
`# species = 0:B,1:Si,2:O` and `# outside_defects = ...` (vacancies,
interstitials, replacements per species) and the per-species columns after
`deposited_ev`. The VTK title line carries `species=` and the same columns
become cell scalars. Both readers accept both versions; v1 files and the
checked-in examples below are unchanged. The file
[`examples/voxel-tally/example.csv`](https://github.com/2AMLogic/lindhard/blob/main/examples/voxel-tally/example.csv)
is a small synthetic report (made-up numbers, not a run):

```text
# lindhard voxel tally v1
# origin_m = 0e0,0e0,0e0
# spacing_m = 5e-9,5e-9,5e-9
# dims = 4,2,2
# histories = 100
# outside = 0,1,3.5e0,0e0,0e0
i,j,k,x_m,y_m,z_m,implanted,displacements,electronic_ev,lattice_ev,rest_ev,deposited_ev
0,0,0,2.5e-9,2.5e-9,2.5e-9,0,12,1.8005e3,9.625e1,0e0,1.89675e3
1,0,0,7.500000000000001e-9,2.5e-9,2.5e-9,3,41,5.210125e3,3.105e2,4.75e0,5.525375e3
...
```

## VTK export

`VoxelReport::write_vtk` writes a legacy ASCII VTK 3.0 file (Kitware, *The VTK
User's Guide*, "VTK File Formats"), which needs no library: an
`UNSTRUCTURED_GRID` with one `VTK_VOXEL` cell (type 11) per **occupied** voxel,
so the file is as sparse as the tally, and the columns above as cell scalars.
ParaView and VisIt open it directly. The grid and the history count are in the
title line. The file does not carry the `outside` sums; use the CSV for those.
`VoxelReport::read_vtk` reads back what `write_vtk` writes, and the checked-in
[`examples/voxel-tally/example.vtk`](https://github.com/2AMLogic/lindhard/blob/main/examples/voxel-tally/example.vtk)
is the same synthetic report as the CSV example. Both example files are
compared with the exporters' output in the test suite, and read back.
