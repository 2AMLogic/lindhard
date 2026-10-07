# The scattering integral

Code: `lindhard/src/ion/scattering.rs` (`theta_quadrature`, `theta_magic`,
`ScatteringTable`).

## Model

A binary collision in a central [potential](potentials.md) is classical
elastic scattering. In the reduced variables
\\( x = r/a \\), \\( \beta = b/a \\),
\\( \varepsilon = a E_\mathrm{cm} / (Z_1 Z_2 e^2) \\), the radial motion in
the centre-of-mass frame obeys

\\[ G(x) = 1 - \frac{\phi(x)}{\varepsilon x} - \frac{\beta^2}{x^2}, \\]

with the distance of closest approach \\( x_0 \\) the root of
\\( G(x_0) = 0 \\). The centre-of-mass scattering angle is (Goldstein 1980;
Ziegler, Biersack and Littmark 1985, ch. 2)

\\[ \theta = \pi - 2\beta \int_{x_0}^{\infty} \frac{dx}{x^2 \sqrt{G(x)}}. \\]

### Distance of closest approach

\\( G \\) is strictly increasing, so the root is unique and lies in
\\( [\beta,\\, (1/\varepsilon + \sqrt{1/\varepsilon^2 + 4\beta^2})/2] \\).
Newton's method runs inside that bracket, and a step that leaves it is
replaced by bisection.

### Gauss-Mehler quadrature

With \\( u = x_0/x \\) and then \\( u = \cos t \\) the integrable
endpoint singularity is removed, leaving

\\[ \int_0^{\pi/2} \frac{dt}{\sqrt{H(\cos t)}}, \qquad
H(u) = \frac{\beta^2}{x_0^2} + \frac{\phi(x_0) - u\\, \phi(x_0/u)}{\varepsilon x_0 (1 - u^2)}, \\]

which the midpoint rule in \\( t \\) (the Gauss-Mehler, or Gauss-Chebyshev
of the first kind, nodes; Mendenhall and Weller 1991) integrates with
exponential convergence. The default is 64 nodes. \\( H \\) is evaluated in
this cancellation-free form.

### Precomputed angle table (the transport hot path)

The BCA does not integrate at every collision. Once per run and per
\\( (Z_1, Z_2, \text{potential}) \\) pair it builds a table of
\\( y = \ln\tan(\theta/2) \\) on a grid uniform in \\( \ln\varepsilon \\)
and \\( \ln\beta \\), by quadrature, and interpolates bilinearly in those
coordinates. \\( y \\) is close to linear in \\( \ln\beta \\) at both ends,
which is what makes bilinear interpolation accurate, and an error
\\( \delta y \\) in the table is at most \\( \delta y \\) radians in
\\( \theta \\). At build time the interpolation error is measured against
direct quadrature at cell centres on a deterministic subset of cells; the
grid and the largest absolute and relative errors found are written to
`summary.json` under `physics.scattering_table`. These are empirical bounds
from sampling, not proofs. The run grid is \\( 10^{-6} \le \varepsilon \le 10^{4} \\),
\\( 10^{-5} \le \beta \le 10^{2} \\), 32 points per decade
(`lindhard::input::TABLE_SPEC`). Beyond the largest tabulated
\\( \beta \\) the angle is taken as 0 (the constant-free-path
\\( p_\mathrm{max} \\) is about 20 screening lengths or less at solid
densities); below the smallest it is clamped; and at a reduced energy outside
the table the engine falls back to direct quadrature (deterministic, only
slower).

### Magic formula (cross-check only)

The Biersack-Haggmark "magic formula" is kept as a fast approximate
cross-check, not used in transport:

\\[ \cos\frac{\theta}{2} = \frac{\beta + \rho + \Delta}{x_0 + \rho}, \qquad
\Delta = \frac{A (x_0 - \beta)}{1 + G}, \\]

\\[ A = 2\alpha\varepsilon\beta^{b}, \quad \alpha = 1 + C_1 \varepsilon^{-1/2}, \quad
b = \frac{C_2 + \varepsilon^{1/2}}{C_3 + \varepsilon^{1/2}}, \quad
G = \frac{\gamma}{\sqrt{1 + A^2} - A}, \quad \gamma = \frac{C_4 + \varepsilon}{C_5 + \varepsilon}, \\]

with \\( \rho \\) the radius of curvature of the trajectory at closest
approach. Constant sets exist for the ZBL universal function (Ziegler,
Biersack and Littmark 1985) and for Molière (Biersack and Haggmark 1980).

### Nuclear stopping

The reduced nuclear stopping cross section follows from the angle,

\\[ s_n(\varepsilon) = 2\varepsilon \int_0^\infty \sin^2\frac{\theta}{2}\\, \beta\\, d\beta, \\]

integrated over \\( \ln\beta \\) by Simpson's rule. It is used in tests and
validation, not in transport (the BCA samples the collisions themselves).
Two references are checked against it: the ZBL universal fit
\\( s_n = \ln(1 + 1.1383\varepsilon) / [2(\varepsilon + 0.01321\varepsilon^{0.21226} + 0.19593\varepsilon^{0.5})] \\)
for \\( \varepsilon \le 30 \\) (Ziegler, Biersack and Littmark 1985), and the
high-energy limit for a sum of exponentials from the first-order impulse
approximation of Lindhard, Nielsen and Scharff (1968, eq. (3.4)).

## Assumptions

- Classical, elastic, binary collisions in a central potential; no
  inelastic energy loss inside the collision itself (electronic loss is
  handled separately, see [Electronic stopping](stopping.md)).
- The angle depends only on \\( (\varepsilon, \beta) \\) and the screening
  function; \\( Z_1, Z_2 \\) and the screening length enter through the
  reduced variables.

## Validity

The classical treatment holds when the collision is well localised compared
with the screening length, which is the case for heavy particles from
cascade energies (a few eV) to well above the MeV range.

## Verification status

The quadrature is checked against the small-angle perturbation formula of
Lindhard, Nielsen and Scharff (1968) and against the ZBL nuclear stopping
fit; the magic formula is checked against the quadrature for both constant
sets. The ZBL magic-formula constants were cross-checked against
`ir2-lab/screened_coulomb` (MIT); the Molière set is unverified against the
paper and is checked only against the quadrature. See
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md).

## References

- H. Goldstein, *Classical Mechanics*, 2nd ed. (Addison-Wesley, 1980).
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions
  in Solids* (Pergamon, New York, 1985), ch. 2.
- M. H. Mendenhall and R. A. Weller, Nucl. Instrum. Methods B 58, 11 (1991).
- J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174, 257 (1980),
  doi:10.1016/0029-554X(80)90440-1.
- J. Lindhard, V. Nielsen, M. Scharff, Mat. Fys. Medd. Dan. Vid. Selsk.
  36 (10) (1968), eqs. (3.3)-(3.4).
