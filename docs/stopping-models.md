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
| Single-oscillator density effect (opt-in) | `bethe::density_effect_single_oscillator` | `βγ > I/ħω_p` (about 5.6 for Si, beyond the Bethe-Bloch 1 GeV/u range); zero below, so no effect in range; underestimates `δ` through the transition | Fermi 1940; Sternheimer 1952; Fano 1963 (specialisation ours) |
| Bohr straggling | `straggling` | high energy, `v >> v0 Z1^(2/3)`; overestimates below about 1 MeV/u | Bohr 1948 |

## What is deliberately not implemented

The clean-room rules (`CONTRIBUTING.md`) rule out tabulated coefficient sets
from SRIM/ZBL, ICRU reports and similar. As a result:

* **Barkas (`L1`) term**: declined (see "Literature findings" below).
* **Shell correction `C/Z2`**: declined. Callers may set a constant from a
  cited source (`BetheBloch::shell_over_z`); default 0.
* **Density effect `δ`**: only the single-oscillator form is implemented
  (opt-in, `BetheBloch::with_density_effect`); the tabulated multi-oscillator
  parameter sets are declined. `density_delta` can still be set as a constant.
* **Mean excitation energy `I`**: defaults to the Bloch rule `10 eV · Z2`
  (rough); supply a cited measured value with `with_mean_excitation_ev`.
* **Biersack-Varelas interpolation joining the low- and high-energy
  regimes**: deferred to a follow-up issue. No verifiable published equation
  for the full interpolation was available, and an unsourced substitute is
  not allowed. Fitted ZBL/Biersack-Varelas coefficients are likewise not used.
* **Straggling**: Chu and Yang-O'Connor-Wang corrections declined (see
  below) and the Lindhard-Scharff low-velocity correction (not verified) is
  omitted; Bohr plus an optional relativistic factor only.
* **Heavy-ion effective charge**: the Barkas empirical form is available in
  `bethe` and is not applied by default.

## Literature findings (issue #41)

Each of the four omitted items was reviewed for a closed form whose
coefficients come from a paper itself. The review was from the contributor's
knowledge of the literature; no paper text was available to check digits
against, so every statement below about a paper's content is "as recalled" and
the declines are on that basis. Nothing was taken from ICRU/SRIM/NIST tables.

* **Barkas `L1`: declined.** Ashley, Ritchie & Brandt (Phys. Rev. B 5, 2393
  (1972)) give `L1 = F(b / x^(1/2)) / (Z2^(1/2) x^(3/2))` with `F` a function
  that is tabulated (and later fitted to it) rather than closed-form. Lindhard's
  1976 treatment of the Barkas effect (NIM 132, 1) gives asymptotic forms from a
  harmonic-oscillator model, whose numerical coefficients we could not
  reproduce with confidence from memory. Jackson & McCarthy (Phys. Rev. B 6,
  4131 (1972)) and later fits carry fitted coefficients that cannot be
  verified here. No form with verifiable coefficients exists to us; an
  invented one is not acceptable.
* **Shell correction: declined.** Walske and Bichsel shell corrections come from
  hydrogenic-shell calculations published as tables in `η = βγ` and `I`;
  the compact empirical forms (e.g. the ICRU 49 expression) are fits with
  coefficients tuned to data and tables, which Tier C of `CONTRIBUTING.md`
  excludes. A first-principles computation (hydrogenic shell sums) is a research
  task of its own, not a closed form.
* **Density effect: partly implemented.** The Sternheimer parameter sets
  (`x0, x1, a, m, C̄`) are tables, which we do not ingest. What is admissible is
  the dielectric formulation itself, computed from a model; we implement its
  one-oscillator specialisation (see `density_effect_single_oscillator`).
  Its high-energy limit is exact, but it switches on too late to matter in
  range. A multi-oscillator version needs per-material oscillator strengths,
  which are tabulated material data. Callers wanting a realistic `δ` in the
  range 1 to 100 GeV/u should supply one from a cited source via
  `density_delta`.
* **Chu / Yang-O'Connor-Wang straggling: declined.** Chu's correction
  (Phys. Rev. A 13, 2057 (1976)) is built from Hartree-Fock-Slater charge
  densities and is published as tables; Yang, O'Connor & Wang (NIM B 61, 149
  (1991)) is a fit whose coefficients are fitted to data tables. Neither has
  a closed form with paper-sourced coefficients that we could verify.

## Deviations from the issue text

* The Lindhard-Scharff `k_L` quoted in the issue omits the factor
  `ξ_e = Z1^(1/6)`. It is included here, because with it the reduced form
  `k_L ε^(1/2)` agrees with the dimensional LS expression to better than 1 %
  for all tested ion/target pairs (and without it Z1 > 1 disagrees by `Z1^(1/6)`).
* The Oen-Robinson local-loss constants are flagged unverified in
  `data-provenance.md`.
