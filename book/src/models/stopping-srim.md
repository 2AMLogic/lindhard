# SRIM-data-sourced stopping

Code: `lindhard/src/ion/stopping/srim.rs` (`SrimStopping`, `SrimRow`).

## Model

The electronic stopping cross section per target atom is evaluated from a
five-coefficient fit whose coefficients are **data read from a file**, not
constants in the source. With \\( E \\) the projectile energy per atomic mass
unit in keV/u and \\( S \\) in eV·10⁻¹⁵ cm² per atom:

\\[ S = A_1 E^{0.45}, \qquad E < 10\ \text{keV/u} \\]

\\[ S_\text{low} = A_2 E^{0.45}, \qquad
   S_\text{high} = \frac{A_3}{E}\ln\\!\left(1 + \frac{A_4}{E} + A_5 E\right),
   \qquad S = \frac{S_\text{low}\,S_\text{high}}{S_\text{low} + S_\text{high}}
   \qquad E \ge 10\ \text{keV/u} \\]

Each row of the file carries the coefficients for one (ion \\( Z_1 \\), target
element \\( Z_2 \\)) pair, the energy span the row is valid over, and a
mandatory citation (edition, table, page). All loss is nonlocal.

## Selecting it

| TOML | Rust |
|---|---|
| `stopping = "srim"` (with `[stopping] srim_data = "file.toml"`) | `StoppingChoice::Srim` |

`srim_data` is required with this choice and rejected with any other. The
file format is `lindhard-srim-coefficients/1` (see
[`docs/cli.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/cli.md)).
Library-only: `SrimStopping::from_toml_str` and `from_toml_file`. It is not
the default.

## Assumptions

- The five-term form above describes the stopping of the pair over the row's
  whole energy span; the coefficients, not the form, carry the pair
  dependence.
- The projectile is characterised by its mass only through \\( E/M \\); no
  effective-charge or isotope correction is applied.
- Compounds are built from element rows by
  [Bragg additivity](stopping-bragg.md).

## Validity

Exactly each row's `e_min_kev_per_u` to `e_max_kev_per_u`. A query outside the
row's span, or for a pair with no row, is an error (`OutOfTableRange`,
`NoSrimCoefficients`) that stops the run; nothing is extrapolated. The
advisory validity range is the union of the ion's row spans.

## Verification status

**The functional form is unverified.** It was written from the published
description of the fit (`A1 E^0.45` at low energy, the harmonic combination
above 10 keV/u) and has **not been checked against the book page**; the
check is part of entering the full tables (#633). The loader's format version
is to be bumped if the form changes. Until then, do not rely on absolute
values from this model.

**Consistency, not independent validation.** The coefficients are data taken
from the SRIM literature, so the output is SRIM-data-sourced: agreement
between this model and SRIM, or any dataset fitted with the same form, shows
that the data were read and evaluated correctly, not that the stopping is
right. This change ships no book values; the unit-test fixture is synthetic
numbers made up for tests. See
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)
("SRIM coefficient file format").

## References

- H. H. Andersen and J. F. Ziegler, *Hydrogen Stopping Powers and Ranges in
  All Elements* (Pergamon, 1977).
- J. F. Ziegler, J. P. Biersack and U. Littmark, *The Stopping and Range of
  Ions in Solids* (Pergamon, 1985), ch. 3.
- The citation on each coefficient row is the reference for that row's
  values.
