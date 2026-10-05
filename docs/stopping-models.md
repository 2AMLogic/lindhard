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
  regimes**: documented but not implemented (issue #44, closed; see
  "Biersack-Varelas interpolation" below). Fitted ZBL/Biersack-Varelas
  coefficients are not used.
* **Straggling**: Chu and Yang-O'Connor-Wang corrections declined (see
  below) and the Lindhard-Scharff low-velocity correction (not verified) is
  omitted; Bohr plus an optional relativistic factor only.
* **Heavy-ion effective charge**: the Barkas empirical form is available in
  `bethe` and is not applied by default.

## Biersack-Varelas interpolation (issue #44, closed)

**Status: documented, not implemented.**

### The joining form

The interpolation joins a low-energy and a high-energy stopping branch
harmonically:

    1/S = 1/S_low + 1/S_high     (equivalently S = S_low S_high / (S_low + S_high))

Two open-access sources were read (rendered pages) and show this form:

* M. V. Moro, PhD thesis, Universidade de Sao Paulo (2017),
  doi:10.11606/T.43.2017.tde-18092017-095345, eq. (3.7), printed p. 40. It
  gives `s_low = A1 E^0.45` and `s_high = (A2/E) ln(1 + A3/E + A4 E)` (E in
  reduced units as defined there), attributes the form to Varelas and
  Biersack (1970), refined by Andersen and Ziegler (1977), and says A1 to A4
  are fitted to experimental data.
* P. de Vera et al., arXiv:2608.15368 (2026), eq. (9), p. 3: the same
  `s = s_low s_high / (s_low + s_high)`, citing Varelas and Biersack, and
  Andersen and Ziegler.

NISTIR 4999 (Berger 1992), section 3.5.2, p. 8, also mentions the "fitting
formula of Varelas and Biersack (1970)" used with ICRU 49 coefficients (no
equation given).

The original papers were **not consulted** (closed access; Unpaywall and
Semantic Scholar report no open copy, ScienceDirect returned 403):

* C. Varelas and J. P. Biersack, Nucl. Instrum. Methods 79, 213 (1970),
  doi:10.1016/0029-554X(70)90141-2.
* J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174, 257 (1980),
  doi:10.1016/0029-554X(80)90440-1.

Citation correction: "Biersack and Varelas, NIM 194, 93 (1982)" is not a
Biersack-Varelas paper. Crossref resolves NIM 194, 93-100 (1982) to Biersack
and Ziegler, "Refined universal potentials in atomic collisions"
(doi:10.1016/0029-554X(82)90496-7). It is not a source for this join.

### Why it is not implemented

The published `s_high = (A2/E) ln(1 + A3/E + A4 E)` is a fully fitted
functional form (de Vera et al. eqs. (10)-(11), following ICRU 49; A1 to A4
are all fitting parameters) built so that it never crosses zero. For positive
coefficients the `1 +` keeps the logarithm's argument above 1, so `s_high > 0`
at every E whatever A3 is. The harmonic join tends to `s_low` at low energy
even with A3 = 0: `s_high` then tends to the finite value `A2 A4` while
`s_low = A1 E^0.45` goes to 0. The fitted `A3/E` term only sets how fast
`s_high` grows as E goes to 0. The coefficients come from the
Andersen-Ziegler / ICRU 49 fits, which the clean-room policy
(`CONTRIBUTING.md`) does not allow.

lindhard's only sourced high-energy branch is `BetheBloch`, which has no such
structure: its bracket crosses zero inside the crossover region. As `S_high` goes to 0+, the join
`S_low S_high / (S_low + S_high)` goes to 0, not to `S_low`. Evaluated with
`LindhardScharff` and `BetheBloch` defaults (Bloch term on, `I = 10 eV Z2`),
in units of 1e-15 eV cm^2/atom (computed by a Python mirror of the formulas,
not by the Rust models):

| ion / target | E/u | LS | Bethe-Bloch | harmonic join |
|---|---|---|---|---|
| H in Si | 80 keV | 27.0 | < 0 (NotApplicable) | undefined |
| H in Si | 100 keV | 30.2 | 6.7 | 5.5 |
| H in Si | 225 keV | 45.3 | 16.8 | 12.3 |
| P in Si | 225 keV | 463 | < 0 | undefined |
| P in Si | 500 keV | 690 | 399 | 253 |
| He in Au | 225 keV | 114 | < 0 | undefined |
| He in Au | 500 keV | 170 | 18.5 | 16.7 |

Bethe-Bloch is negative for H in Si below roughly 90 keV, for P in Si up to
about 400 keV/u, and for He in Au up to about 450 keV/u. Neither workaround
is acceptable. Falling back to `S_low` below the Bethe zero makes S jump from
about 0 to the full LS value (about 27 for H in Si near 90 keV), which breaks
continuity. Propagating `NotApplicable` leaves the model undefined over the
whole low-energy regime, so ions slowing down through it cannot be
transported. Restricting the join to where Bethe-Bloch is positive violates
the low-energy limit and gives badly wrong stopping near the zero. A fix
would need a high branch that stays positive at low E. The options are the
fitted A1 to A4 (Tier C data), or a regularisation of our own, such as
grafting the `ln(1 + ...)` structure onto Bethe by identifying A2 and A4 with
Bethe quantities and setting A3 = 0. That identification is our own
construction, which neither source makes. Neither option is admissible.

### When it could be revisited

* A published, coefficient-free high-energy branch that stays positive below
  the Bethe zero (for example a Bethe form with a closed-form, non-tabulated
  shell correction; the review of the omitted correction terms under #41,
  now closed, is recorded below) is found and can be cited to an equation.
* The A1 to A4 coefficients become available from a source whose licence
  permits use, or are derived from our own fits to data that may be used.
* The primary papers are read and show a high branch that stays positive
  without fitted coefficients, or that ties the `ln(1 + ...)` form to Bethe
  quantities in a way that can be cited.

## Literature findings (issue #41, closed)

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
