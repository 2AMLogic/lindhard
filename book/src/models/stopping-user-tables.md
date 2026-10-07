# User stopping tables

Code: `lindhard/src/ion/stopping/table.rs` (`StoppingTable`,
`TableOverride`).

## Model

A user table gives \\( S_e(E) \\) for one projectile (atomic number and
mass) in one target element, as data the user supplies. Between the
tabulated points the cross section is interpolated piecewise linearly in
\\( \ln S \\) against \\( \ln E \\), which keeps the monotonicity of each
segment of the data. Nothing is extrapolated or clamped: a query outside
the table's energy range is an error (`OutOfTableRange`) that stops the run.

```toml
provenance = "Author, Journal vol, page (year), Table N"   # required
ion_z = 5
ion_mass_amu = 11.0093            # optional; default: standard atomic weight
target_z = 14
energy_ev = [1.0e2, 1.0e3, 1.0e4]                 # strictly increasing, eV
stopping_ev_1e15_cm2 = [10.0, 30.0, 60.0]         # eV 1e-15 cm^2 per atom
```

(The numbers are a format illustration, not data.)

## Provenance is mandatory

A table without a non-empty `provenance` string is rejected at load time:
data without an origin is not admitted. The run records each table in
`summary.json` (`physics.stopping_tables`: the path as written, the resolved
path, the SHA-256 of the file's bytes, the provenance string, the pair and
the energy range), and lists it under `physics.models` as `user-table` with
its provenance as the source.

**Do not load SRIM- or ICRU-derived tables into anything committed to a
repository.** A table is the user's own data and its terms are the user's
concern; this project's own tree never contains such tables.

## Selecting it

Declare the files under `[stopping] tables` (see
[TOML input reference](../guide/input.md)). A table replaces the
`[physics] stopping` model for exactly the (`ion_z`, `target_z`) pair it
declares, including recoils of that species when recoils are followed; every
other pair uses the model. A table cannot be combined with
`stopping = "equipartition-ls-or"`, which carries its own Lindhard-Scharff
and Oen-Robinson loss.

## Assumptions

- The table is tied to the projectile mass it was declared for. A query for
  a mass that differs by more than a small relative tolerance (tight enough
  to separate neighbouring isotopes) is an error (`TableMassMismatch`): the
  same energy at a different mass is a different speed, and no energy-axis
  conversion is attempted.
- Compounds are built from element tables by
  [Bragg additivity](stopping-bragg.md); a compound table is not a supported
  input.

## Validity

Exactly the table's energy range. The CLI checks up front that the range
contains the beam energy, and for a recoil species that it starts at or
below `recoil_cutoff_ev`; it warns if a table starts above
`primary_cutoff_ev`, because the run fails if a projectile slows below it.

## Verification status

The loader, the interpolation and the range and mass checks are unit-tested.
The tests use tables generated at test time from our own
[Lindhard-Scharff](stopping-lindhard-scharff.md) model; no table is
committed ([`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)).

## References

The table's own `provenance` string is its reference. The format and the
rules are this project's.
