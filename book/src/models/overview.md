# Conventions and the clean-room position

## Units

The library works in SI internally. At the input and output boundary the
units are in the key names: energies in eV (`_ev`), lengths in nm (`_nm`),
angles in degrees (`_deg`), densities in g/cm³. Electronic stopping cross
sections are reported per target atom, in the customary
eV·10⁻¹⁵ cm² (`to_ev_1e15_cm2` converts from J m²).

Throughout the manual \\( Z_1, M_1 \\) are the atomic number and mass of the
moving particle, \\( Z_2, M_2 \\) those of the target atom, \\( E \\) the
laboratory kinetic energy, \\( N \\) the atom density, \\( a_0 \\) the Bohr
radius, \\( v_0 = \alpha c \\) the Bohr velocity, and
\\( e^2 \\) stands for \\( e^2 / 4\pi\varepsilon_0 \\). Physical constants are
CODATA 2022 (Mohr et al. 2025).

## Reduced variables

The collision models use the reduced variables of Lindhard, Scharff and
Schiøtt (1963). With a screening length \\( a \\) (see
[Screening lengths](screening-lengths.md)):

\\[ x = \frac{r}{a}, \qquad \beta = \frac{b}{a}, \qquad
\varepsilon = \frac{a\\, E_\mathrm{cm}}{Z_1 Z_2 e^2}, \qquad
E_\mathrm{cm} = E \frac{M_2}{M_1 + M_2}, \\]

with \\( r \\) the separation and \\( b \\) the impact parameter. In these
variables the scattering angle depends only on \\( (\varepsilon, \beta) \\)
and on the screening function.

## Determinism

Every model is a pure function of its arguments. Random numbers come from a
counter-based stream (ChaCha8) keyed on the run seed and the index of the
primary history, so a run gives the same bits on one thread or many. A model
page that introduces randomness says which draws it makes and in which order.

## The clean-room position

lindhard implements published physics from the papers. It never uses the
code of closed or copyleft programs (the tiers are set out in
`CONTRIBUTING.md`). Data is different: data is not licensed, code is
(operator ruling, 2026-11-10, lindhard#629), and lindhard aims to be a
complete open replacement for SRIM/TRIM. In particular:

- **Tabulated stopping is admissible data**, SRIM, ICRU and ZBL tables
  included, as long as it carries its provenance. The electronic stopping
  models are computed from closed-form physics by default; a table, from a
  user or committed to the tree, is a separately named stopping source, and a
  comparison against the code it came from is a consistency check, not
  independent validation.
- The ZBL *universal screening function* (eight published coefficients) and
  the ZBL reduced nuclear stopping fit are used. The electronic stopping
  coefficient sets of the same book are not used yet; they are admissible
  data.
- Terms whose only published coefficients are tabulated (the Barkas term,
  shell corrections, multi-oscillator density-effect parameters, Chu and
  Yang-O'Connor-Wang straggling) were declined under the former no-tables
  rule. They are now admissible and stay listed in
  [Validity ranges and declined terms](stopping-validity.md) until each is
  implemented.

Every coefficient that enters the code as a number has a row in
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)
saying where it comes from and how far it has been verified. Where a value
has only been checked against a secondary source, its page says so.

## How a model page is laid out

Each page has the same sections: what the model is, the equations, the
assumptions, the validity range, how to select it (the TOML key and the Rust
type), its verification status, and its references. New models start from
the template `book/src/models/_template.md` in the repository.

## References

- J. Lindhard, M. Scharff, H. E. Schiøtt, Mat. Fys. Medd. Dan. Vid. Selsk.
  33 (14) (1963).
- P. J. Mohr, D. B. Newell, B. N. Taylor, E. Tiesinga, Rev. Mod. Phys. 97,
  025002 (2025), doi:10.1103/RevModPhys.97.025002 (CODATA 2022).
- D. J. Bernstein, *ChaCha, a variant of Salsa20* (2008); J. K. Salmon,
  M. A. Moraes, R. O. Dror, D. E. Shaw, Proc. SC'11 (2011) (counter-based
  random streams).
