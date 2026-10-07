# Stopping-data format

A plain-text format for measured electronic stopping points, one point per
line, each with its own uncertainty, method and citation. The loader is
`lindhard::ion::stopping::dataset::StoppingDataset`, and `dataset::compare`
compares any `ElectronicStopping` model against a dataset.

**No real data is shipped in this format yet.** Which published measurements
may enter the tree is decided in issue #34. Until that is settled, the only
datasets in the repository are synthetic ones generated at test time from
our own models (see [`data-provenance.md`](data-provenance.md)).

## Example

```text
# format: lindhard-stopping-data/1
# name: example
# description: free text, optional
ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation
5,,Si,1.0e4,12.5,0.3,transmission,"Author, Journal vol, page (year), Table 2"
5,11.009,SiO2,2.0e4,10.1,0,backscattering,"Database name, record 1234"
```

(The numbers above are placeholders, not measurements.)

## Header

Lines before the column row that start with `#` are the header. A header line
of the form `# key: value`, where `key` is lower-case letters, digits and
underscores, is a metadata entry. Any other `#` line is a comment.

| Key | Required | Meaning |
|---|---|---|
| `format` | yes | Must be `lindhard-stopping-data/1` |
| `name` | yes | Non-empty dataset name |
| `description` | no | Free text |
| anything else | no | Kept, in file order, by `StoppingDataset::metadata` |

A key given twice is an error.

## Columns

The first non-blank line that does not start with `#` is the column row. It
must list exactly these columns, in this order. The units are part of the
names, so a file in another unit (MeV cm²/mg, keV/nm, energy in keV) cannot be
read by accident.

| Column | Unit | Rule |
|---|---|---|
| `ion_z` | | Atomic number of the projectile, 1 to 92 |
| `ion_mass_amu` | u | Projectile mass. Empty means the standard atomic weight. Otherwise finite and positive. Points with different masses are different systems |
| `target` | | Element symbol (`Si`) or formula (`SiO2`, `Al2O3`, `Si0.5Ge0.5`): case-sensitive symbols, each with an optional positive count. No parentheses, and no element repeated. A formula with one element (`O2`) is an element target |
| `energy_ev` | eV | Projectile kinetic energy in the laboratory frame. Finite and positive |
| `s_e_ev_1e15_cm2` | eV·10⁻¹⁵ cm² per atom | Electronic stopping cross section. For a compound, per **average atom**: the stopping per molecule divided by the number of atoms in the formula. Finite and positive |
| `sigma_ev_1e15_cm2` | eV·10⁻¹⁵ cm² per atom | One-standard-deviation uncertainty of `s_e`. Finite and non-negative; `0` means none was stated. Not optional |
| `method` | | Measurement method, free text. Must not be empty |
| `citation` | | Source of the point, free text: a paper, or a database name and record identifier. Must not be empty |

## Lines and quoting

- Fields are separated by commas. Whitespace around a field is ignored.
- A field containing a comma is wrapped in double quotes; `""` inside quotes
  is a literal quote. A quoted field cannot span lines.
- After the column row, blank lines and lines starting with `#` are ignored.
- A leading byte-order mark and CRLF line endings are accepted.
- The same system may have several points at the same energy (repeated
  measurements).

## Validation

Every rule above is checked at load time. A violation is an error naming the
1-based line (`StoppingError::InvalidDataset { line, reason }`), and nothing
is loaded. In particular a row with an empty `citation` or `method` is
rejected: data without a stated origin is not admitted. A file with no
column row, or with no points, is `StoppingError::InvalidTable`.

## Comparison

`compare(model, correction, dataset)` evaluates the model at every point:

- element targets with `ElectronicStopping::stopping`;
- compound targets with `bragg::bragg_cross_section_per_atom` and the given
  `CompoundCorrection` (Bragg and Kleeman, Phil. Mag. 10, 318 (1905)). The
  correction is not applied to element targets.

The material built for a compound carries the formula as its name and a
placeholder density (`PLACEHOLDER_DENSITY_KG_M3`), since a dataset does not
record one and the stopping per atom under Bragg's rule does not depend on it.
A correction used with `compare` must not depend on the density.

Per point it returns, with the sign convention **positive when the model is
above the measurement**:

- `relative = model / measured − 1`;
- `absolute = model − measured`, eV·10⁻¹⁵ cm²;
- `pull = (model − measured) / sigma`, absent when sigma is 0.

Summary statistics are given over all points and per ion–target system (same
ion Z, same ion mass, same target as written; in order of first appearance):
the count `n`, the number skipped, the mean and the root mean square of
`relative`, and the largest `|relative|`. With no evaluated points the three
statistics are NaN.

A point where the model returns an error (for example `NotApplicable` outside
its formula's range) is listed as skipped with that error. It is never dropped
silently.

All sums run sequentially in dataset order, so the results are bitwise the
same whatever the caller's thread count.
