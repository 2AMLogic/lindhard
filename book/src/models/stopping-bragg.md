# Bragg additivity

Code: `lindhard/src/ion/stopping/bragg.rs` (`bragg_cross_section_per_atom`,
`CompoundCorrection`).

## Model

The stopping cross section of a compound or mixture, per average atom, is
the atom-fraction-weighted sum of the element cross sections (Bragg and
Kleeman 1905):

\\[ S_\mathrm{compound} = f \sum_j x_j S_j, \qquad -\frac{dE}{dx} = N f \sum_j x_j S_j, \\]

with \\( x_j \\) the atom fractions, \\( S_j \\) the element cross sections
from the chosen model (or a [user table](stopping-user-tables.md) for that
pair), \\( N \\) the total atom density and \\( f \\) a per-compound
correction factor.

## Selecting it

Always on: every layer's electronic loss goes through the Bragg sum, which
for a pure element is just that element's cross section. `summary.json`
lists it under `physics.models` as `bragg-additivity`. The CLI applies no
correction (\\( f = 1 \\)). The library hook `CompoundCorrection` (a
constant or energy-dependent factor, finite and non-negative) is there for
chemical or phase corrections, such as cores-and-bonds schemes, supplied by
the caller from a cited source.

## Assumptions

- Each atom stops the ion independently of its chemical environment.

## Validity

Wherever the element models apply. Bragg additivity ignores the chemical
and physical state of the target, and no correction for it is applied
unless one is supplied.

## Verification status

The sum and the rejection of invalid correction factors are unit-tested.

## References

- W. H. Bragg and R. Kleeman, Phil. Mag. 10, 318 (1905).
