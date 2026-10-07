# Interatomic potentials

Code: `lindhard/src/ion/potential.rs` (`Screening`, `Potential`).

## Model

Two atoms at separation \\( r \\) interact through a screened Coulomb
potential

\\[ V(r) = \frac{Z_1 Z_2 e^2}{r}\\, \phi\\!\left(\frac{r}{a}\right), \\]

where \\( \phi \\) is a universal screening function (\\( \phi(0) = 1 \\),
decreasing to 0) and \\( a \\) a screening length that depends on
\\( Z_1, Z_2 \\) ([Screening lengths](screening-lengths.md)). The scattering
code works in \\( x = r/a \\), where \\( \phi \\) is the only input; the
screening length only converts to and from SI.

Three of the four screening functions are sums of exponentials,

\\[ \phi(x) = \sum_i c_i\\, e^{-b_i x}, \\]

and the fourth is a polynomial times an exponential.

## The screening functions

| Function | TOML (`[physics]`) | Rust | Default length |
|---|---|---|---|
| ZBL universal | `potential = "zbl"` (default) | `PotentialChoice::Zbl`, `Screening::ZblUniversal` | universal |
| Kr-C | `potential = "kr-c"` | `PotentialChoice::KrC`, `Screening::KrC` | Firsov |
| Molière | `potential = "moliere"` | `PotentialChoice::Moliere`, `Screening::Moliere` | Firsov |
| Lenz-Jensen | `potential = "lenz-jensen"` | `PotentialChoice::LenzJensen`, `Screening::LenzJensen` | Lindhard |

The default length is the one each function was introduced or fitted with
(`Screening::default_length`); `screening_length` overrides it.

### ZBL universal

Four exponentials fitted by Ziegler, Biersack and Littmark (1985, ch. 2) to
Hartree-Fock-Slater solid-state pair potentials:

| \\( c_i \\) | 0.1818 | 0.5099 | 0.2802 | 0.02817 |
|---|---|---|---|---|
| \\( b_i \\) | 3.2 | 0.9423 | 0.4029 | 0.2016 |

These are the commonly printed four-digit rounding of the published values.

### Kr-C

Three exponentials fitted to the Hartree-Fock Kr-Kr pair interaction (Wilson,
Haggmark and Biersack 1977):

| \\( c_i \\) | 0.190945 | 0.473674 | 0.335381 |
|---|---|---|---|
| \\( b_i \\) | 0.278544 | 0.637174 | 1.919249 |

### Molière

Molière's three-exponential approximation to the Thomas-Fermi screening
function (Molière 1947):

| \\( c_i \\) | 0.35 | 0.55 | 0.10 |
|---|---|---|---|
| \\( b_i \\) | 0.3 | 1.2 | 6.0 |

### Lenz-Jensen

A polynomial-times-exponential approximation to the Thomas-Fermi-Jensen
statistical model (Lenz 1932; Jensen 1932), in the form printed by Möller
(2017, p. 11, eq. (28)):

\\[ \phi(x) = \left(1 + q + 0.3344\\, q^2 + 0.0485\\, q^3 + 0.002647\\, q^4\right) e^{-q},
\qquad q = \sqrt{9.67\\, x}. \\]

## Assumptions

- Pair potentials: the interaction of two atoms does not depend on any third
  atom.
- The screening function is universal: \\( Z_1 \\) and \\( Z_2 \\) enter
  only through the screening length.
- The potentials are purely repulsive. There is no attractive well, so they
  are not meant for energies comparable with chemical binding (a few eV),
  where the cutoffs and binding energies of the [BCA](bca.md) take over.

## Validity

Screened Coulomb potentials describe the repulsive part of the interaction,
from the close collisions of keV ions down to the tens of eV of cascade
atoms. The ZBL universal function was fitted to pair potentials over a broad
set of ion-target pairs and is the usual default. It is not exact for any
particular pair: the computed ranges of B in amorphous Si run long against
the measurement, and the validation attributes the offset below about 5 keV
mainly to the nuclear stopping (the ZBL function being too soft for B on Si)
([`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md),
level 3).

## Verification status

From [`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md):
the ZBL, Kr-C and Molière coefficients have been cross-checked against an
independent MIT-licensed implementation (`ir2-lab/screened_coulomb`, commit
f84c3c8), a secondary source; the Lenz-Jensen set agrees with a textbook
tabulation (Möller 2017) but has not been checked against the primary
papers, which were not accessible.

## References

- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions
  in Solids* (Pergamon, New York, 1985), ch. 2.
- W. D. Wilson, L. G. Haggmark, J. P. Biersack, Phys. Rev. B 15, 2458 (1977).
- G. Molière, Z. Naturforsch. A 2, 133 (1947).
- W. Lenz, Z. Phys. 77, 713 (1932), doi:10.1007/BF01342150.
- H. Jensen, Z. Phys. 77, 722 (1932), doi:10.1007/BF01342151.
- W. Möller, *Fundamentals of Ion-Solid Interaction*, HZDR-073
  (Helmholtz-Zentrum Dresden-Rossendorf, 2017), p. 11, eqs. (26)-(30),
  <https://www.hzdr.de/publications/PublDoc-10091.pdf>.
