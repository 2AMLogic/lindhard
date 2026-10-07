# Lindhard-Scharff

Code: `lindhard/src/ion/stopping/lindhard_scharff.rs` (`LindhardScharff`).

## Model

At ion velocities below \\( v_0 Z_1^{2/3} \\) the electronic stopping is
proportional to the velocity (Lindhard and Scharff 1961). In the reduced
units of Lindhard, Scharff and Schiøtt (1963):

\\[ a = 0.8853\\, a_0 \left(Z_1^{2/3} + Z_2^{2/3}\right)^{-1/2}, \qquad
\varepsilon = \frac{E\\, a\\, M_2}{Z_1 Z_2 e^2 (M_1 + M_2)}, \qquad
\rho = N x\\, 4\pi a^2 \frac{M_1 M_2}{(M_1 + M_2)^2}, \\]

the reduced electronic stopping is

\\[ \left(\frac{d\varepsilon}{d\rho}\right)_e = k_L\\, \varepsilon^{1/2}, \qquad
k_L = \xi_e\\, \frac{0.0793\\, Z_1^{1/2} Z_2^{1/2} (A_1 + A_2)^{3/2}}
{\left(Z_1^{2/3} + Z_2^{2/3}\right)^{3/4} A_1^{3/2} A_2^{1/2}}, \qquad
\xi_e = Z_1^{1/6}, \\]

with \\( A_1, A_2 \\) the masses in u. The equivalent dimensional form is

\\[ S_e = \xi_e\\, \frac{8\pi e^2 a_0\\, Z_1 Z_2}{\left(Z_1^{2/3} + Z_2^{2/3}\right)^{3/2}}\\, \frac{v}{v_0}. \\]

The factor \\( \xi_e \approx Z_1^{1/6} \\) (order 1 to 2) is part of the
Lindhard-Scharff result. A test checks the reduced form, with the rounded
0.0793, against the dimensional form to 1 % for all tested ion-target pairs.

A per-element multiplicative correction \\( f(Z_2) \\) can be attached in
the library (`LindhardScharff::with_correction`, default 1; it must be finite
and non-negative, and 0 switches the element's stopping off). The CLI does
not expose it; use a [user table](stopping-user-tables.md) for measured
stopping instead.

## Selecting it

`stopping = "lindhard-scharff"` in `[physics]` (the default),
`StoppingChoice::LindhardScharff`. All of the loss is applied nonlocally,
along the free flights. It is also the nonlocal half of
`stopping = "equipartition-ls-or"` ([Oen-Robinson](stopping-oen-robinson.md)).

## Assumptions

- A free-electron-gas picture with Thomas-Fermi scaling; the stopping is
  smooth in \\( Z_1 \\) and \\( Z_2 \\) and has no shell structure (the
  measured \\( Z_1 \\) oscillations are not reproduced).
- No dependence on the chemical or physical state of the target.

## Validity

\\( v < v_0 Z_1^{2/3} \\), that is \\( E/A_1 \\) below about
\\( 25\ \mathrm{keV} \cdot Z_1^{4/3} \\). This is the regime of keV-range
implantation and of cascade atoms. Above it the stopping peaks and turns
over, which this model does not describe; see
[Bethe-Bloch](stopping-bethe-bloch.md) for high velocities and
[Validity ranges and declined terms](stopping-validity.md) for why no
interpolation joins the two.

The validation pages report a known offset: for B in Si the computed ranges
run long against the measurement, and at 10 to 20 keV part of the offset is
attributed to LS stopping being too small for that pair
([`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md),
level 3).

## Verification status

The constants 0.8853, 0.0793 and \\( \xi_e = Z_1^{1/6} \\) are checked
against the dimensional form in a unit test (1 %); they have not yet been
checked digit by digit against the papers
([`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)).

## References

- J. Lindhard and M. Scharff, Phys. Rev. 124, 128 (1961).
- J. Lindhard, M. Scharff, H. E. Schiøtt, Mat. Fys. Medd. Dan. Vid. Selsk.
  33 (14) (1963).
