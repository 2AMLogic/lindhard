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

A `[stopping]` table replaces the chosen model for the one (ion, target
element) pair it declares ([User stopping tables](stopping-user-tables.md)).
Energy-loss straggling is described on its [own page](straggling.md). It is
a library function that the BCA engine does not call at present: the
electronic loss along each flight is deterministic, \\( N S_e(E)\\, s \\).

The validity ranges side by side, and the terms that are deliberately not
implemented, are on [Validity ranges and declined terms](stopping-validity.md).
