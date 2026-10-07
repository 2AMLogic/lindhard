# Electronic stopping

Code: `lindhard/src/ion/stopping/`.

A moving atom loses energy to the target electrons as well as to the nuclei.
In the BCA the electronic loss is handled separately from the
[collisions](scattering.md): continuously along each free flight (nonlocal),
at the collisions as a function of the distance of closest approach (local),
or as a mix of the two ([BCA transport](bca.md), "Electronic loss").

## Conventions

- Every model returns the **stopping cross section per target atom**,
  \\( S_e = -\frac{1}{N}\frac{dE}{dx} \\), in J m² internally and
  eV·10⁻¹⁵ cm² at the boundary.
- Energies are the laboratory kinetic energy of the moving atom.
- Every model is a pure function of (ion, target \\( Z_2 \\), energy): no
  random numbers, no state.
- Each model exposes an **advisory** validity range
  (`ElectronicStopping::validity`). The CLI warns on stderr when the beam
  energy lies outside it; it does not refuse the run.
- Compounds and mixtures combine the element cross sections by
  [Bragg additivity](stopping-bragg.md).

## Choosing a model

| `[physics] stopping` | Rust | Loss mode | Page |
|---|---|---|---|
| `stopping = "lindhard-scharff"` (default) | `StoppingChoice::LindhardScharff` | all nonlocal | [Lindhard-Scharff](stopping-lindhard-scharff.md) |
| `stopping = "bethe-bloch"` | `StoppingChoice::BetheBloch` | all nonlocal | [Bethe-Bloch](stopping-bethe-bloch.md) |
| `stopping = "equipartition-ls-or"` | `StoppingChoice::EquipartitionLsOr` | half nonlocal (LS), half local (Oen-Robinson) | [Oen-Robinson](stopping-oen-robinson.md) |
| `stopping = "none"` | `StoppingChoice::None` | none: nuclear loss only | below |

A `[stopping]` table replaces the chosen model for the one (ion, target
element) pair it declares ([User stopping tables](stopping-user-tables.md)).
Energy-loss straggling is described on its [own page](straggling.md). It is
a library function that the BCA engine does not call at present: the
electronic loss along each flight is deterministic, \\( N S_e(E)\\, s \\).

## No electronic stopping

`stopping = "none"` (`StoppingChoice::None`, model
`lindhard::ion::stopping::none::NoStopping`) sets \\( S_e = 0 \\) for every
ion, target and energy, so the moving atoms lose energy only in nuclear
collisions. This is the \\( k = 0 \\) limit of the range theory of Lindhard,
Scharff and Schiott (Mat. Fys. Medd. Dan. Vid. Selsk. 33 (14), 1963), where
the electronic stopping coefficient vanishes.

It is **not a physical model** of any target. It is there for two uses:
like-for-like transport comparisons against codes run with their electronic
stopping switched off (the `*_nuclear` level-2 problems of the validation
report), and nuclear-only studies such as the level-1 `range.si5k_si.*` and
damage checks, which use the same type. It has no validity range, so no
warning is raised for it, and the electronic terms of the energy budget are
exactly zero. `[stopping]` tables may still be given; they serve the pairs
they declare.

The validity ranges side by side, and the terms that are deliberately not
implemented, are on [Validity ranges and declined terms](stopping-validity.md).
