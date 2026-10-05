# Electronic stopping models: validity ranges

Code: `lindhard/src/ion/stopping/`. All models return the stopping cross section
per atom in J m² (`to_ev_1e15_cm2` converts to eV·10⁻¹⁵ cm²) and expose the
ranges below at run time through `ElectronicStopping::validity`. Ranges are
advisory; the models do not refuse energies outside them (Bethe-Bloch returns
`NotApplicable` where its bracket is not positive, and user tables return
`OutOfTableRange`). `v0 = α c` is the Bohr velocity, `Z1` the ion charge number.

| Model | Module | Intended range | Source |
|---|---|---|---|
| Lindhard-Scharff | `lindhard_scharff` | `v < v0 Z1^(2/3)` (E/A below about 25 keV · Z1^(4/3)); stopping ∝ `v` | Lindhard & Scharff, Phys. Rev. 124, 128 (1961) |
| Oen-Robinson (local) | `oen_robinson` | same as LS; impact-averaged value equals LS | Oen & Robinson, NIM 132, 647 (1976) |
| Equipartition LS/OR | `mix` | same as LS | as above |
| Bethe-Bloch | `bethe` | `v >= 3 v0 Z1^(2/3)` up to 1 GeV/u (no density effect) | Bethe 1930/32; Bloch 1933; Fano 1963 |
| User table | `table` | exactly the table's energy range | the table's own `provenance` |
| Bragg additivity | `bragg` | where the element models apply; ignores chemical state unless a correction is supplied | Bragg & Kleeman 1905 |
| Bohr straggling | `straggling` | high energy, `v >> v0 Z1^(2/3)`; overestimates below about 1 MeV/u | Bohr 1948 |

## What is deliberately not implemented

The clean-room rules (`CONTRIBUTING.md`) rule out tabulated coefficient sets
from SRIM/ZBL, ICRU reports and similar. As a result:

* **Barkas (`L1`) term** (Ashley, Ritchie, Brandt 1972): needs a tabulated
  function. Omitted. Expect errors of a few percent at 1 to 10 MeV for protons
  in the Bethe-Bloch model.
* **Shell correction `C/Z2`** and **density effect `δ`**: the published
  parameter sets are tables. Both default to 0; callers may set constants from
  a cited source (`BetheBloch::shell_over_z`, `density_delta`).
* **Mean excitation energy `I`**: defaults to the Bloch rule `10 eV · Z2`
  (rough); supply a cited measured value with `with_mean_excitation_ev`.
* **Biersack-Varelas interpolation joining the low- and high-energy
  regimes**: deferred to a follow-up issue. No verifiable published equation
  for the full interpolation was available, and an unsourced substitute is
  not allowed. Fitted ZBL/Biersack-Varelas coefficients are likewise not used.
* **Straggling**: Chu and Yang-O'Connor-Wang corrections (fits/tables) and the
  Lindhard-Scharff low-velocity correction (not verified) are omitted; Bohr
  plus an optional relativistic factor only.
* **Heavy-ion effective charge**: the Barkas empirical form is available in
  `bethe` and is not applied by default.

## Deviations from the issue text

* The Lindhard-Scharff `k_L` quoted in the issue omits the factor
  `ξ_e = Z1^(1/6)`. It is included here, because with it the reduced form
  `k_L ε^(1/2)` agrees with the dimensional LS expression to better than 1 %
  for all tested ion/target pairs (and without it Z1 > 1 disagrees by `Z1^(1/6)`).
* The Oen-Robinson local-loss constants are flagged unverified in
  `data-provenance.md`.
