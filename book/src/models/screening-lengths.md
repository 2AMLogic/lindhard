# Screening lengths

Code: `lindhard/src/ion/potential.rs` (`ScreeningLength`).

## Model

The screening length \\( a(Z_1, Z_2) \\) scales the separation in the
[screening function](potentials.md), \\( x = r/a \\). Three forms are
implemented. Two of them use the Thomas-Fermi constant

\\[ C_\mathrm{TF} = \left(\frac{9\pi^2}{128}\right)^{1/3} = 0.8853, \\]

in the closed form of Firsov (1958, p. 535), printed as 0.8853 by Lindhard,
Scharff and Schiøtt (1963, p. 8).

| Length | Formula | TOML (`[physics]`) | Rust |
|---|---|---|---|
| Universal | \\( a_U = 0.8853\\, a_0 / (Z_1^{0.23} + Z_2^{0.23}) \\) | `screening_length = "universal"` | `LengthChoice::Universal`, `ScreeningLength::Universal` |
| Firsov | \\( a_F = 0.8853\\, a_0 / (Z_1^{1/2} + Z_2^{1/2})^{2/3} \\) | `screening_length = "firsov"` | `LengthChoice::Firsov`, `ScreeningLength::Firsov` |
| Lindhard (Thomas-Fermi) | \\( a_L = 0.8853\\, a_0 / (Z_1^{2/3} + Z_2^{2/3})^{1/2} \\) | `screening_length = "lindhard"` | `LengthChoice::Lindhard`, `ScreeningLength::Lindhard` |

Without `screening_length` the length paired with the potential is used:
universal for ZBL, Firsov for Kr-C and Molière, Lindhard for Lenz-Jensen.
The run echoes the length it used in `summary.json`
(`input.physics.screening_length`).

## Assumptions

- The universal length is the one the ZBL screening function was fitted
  with; the other two come from Thomas-Fermi scaling arguments for the
  two-atom system.
- Mixing a screening function with a length other than its own is allowed
  (it is a model choice), but the fitted functions were fitted with their own
  length.

## Validity

As for the [potentials](potentials.md): screened Coulomb collisions from
cascade energies upward.

## Verification status

From [`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md):
the Firsov and Lindhard lengths have been verified against scans of the
primary papers (page references above). The universal length uses the same
constant: the SRIM book (Ziegler, Biersack and Ziegler, 2015, Ch. 2) confirms
the exponent 0.23 and prints the prefactor as 0.8853 in Eq. 2-73 (p. 2-30),
while the labels of its Figs. 2-16 and 2-17 (p. 2-32) print 0.8854, a 10⁻⁴
relative difference. The code follows the equation and uses the Thomas-Fermi
constant, 0.88534.

## References

- O. B. Firsov, Sov. Phys. JETP 6, 534 (1958), pp. 535-536.
- J. Lindhard, M. Scharff, H. E. Schiøtt, Mat. Fys. Medd. Dan. Vid. Selsk.
  33 (14) (1963), p. 8.
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions
  in Solids* (Pergamon, New York, 1985), ch. 2.
- J. F. Ziegler, J. P. Biersack, M. D. Ziegler, *SRIM: The Stopping and Range
  of Ions in Matter* (SRIM Co., 2015), ch. 2, pp. 2-29 to 2-32.
