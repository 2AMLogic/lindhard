# Oen-Robinson and the equipartition mix

Code: `lindhard/src/ion/stopping/oen_robinson.rs` (`OenRobinson`),
`lindhard/src/ion/stopping/mix.rs` (`EquipartitionMix`).

## Model

Oen and Robinson (1976) make the electronic loss *local*: it is taken at
each collision and depends on how close the two atoms came. A collision with
distance of closest approach \\( r_\mathrm{min} \\) loses

\\[ \Delta E_e(E, r_\mathrm{min}) = S_\mathrm{LS}(E)\\, \frac{c^2}{2\pi a^2}\\,
\exp\\!\left(-\frac{c\\, r_\mathrm{min}}{a}\right), \qquad c = 0.3, \\]

with \\( S_\mathrm{LS} \\) the [Lindhard-Scharff](stopping-lindhard-scharff.md)
cross section and \\( a \\) the Firsov screening length
\\( 0.8853\\, a_0 (Z_1^{1/2} + Z_2^{1/2})^{-2/3} \\). The prefactor
normalises \\( \int \Delta E_e\\, 2\pi p\\, dp = S_\mathrm{LS} \\) when
\\( r_\mathrm{min} \\) is identified with the impact parameter \\( p \\);
that identification is the approximation of the original paper, so the
impact-averaged stopping equals the LS value.

### Equipartition mix

The equipartition mix splits the Lindhard-Scharff loss in two equal halves:
half nonlocal, continuous along the free flight, and half local, at each
collision, by the Oen-Robinson formula evaluated at the distance of closest
approach. This combination is the one described in the BCA literature that
builds on Oen and Robinson (1976). Averaged over impact parameter the mix
equals the LS stopping (a unit test checks this to 10⁻¹⁴).

## Selecting it

`stopping = "equipartition-ls-or"` in `[physics]`, which is
`StoppingChoice::EquipartitionLsOr` in the input and
`ElectronicLoss::EquipartitionLsOr` in the engine (see
[BCA transport](bca.md)). The stopping model passed to the engine is then not
used: this mode carries its own LS and Oen-Robinson losses, so it cannot be
combined with `[stopping]` tables. The local losses are reported in the
energy budget as `electronic_local`.

## Assumptions

- \\( r_\mathrm{min} \approx p \\) in the normalisation (Oen and Robinson
  1976).
- The local part is only sampled at the collisions the BCA makes, out to
  \\( p_\mathrm{max} \\) (to \\( p_\mathrm{max}\sqrt{K + 1} \\) with \\( K \\)
  weak collisions). Where that radius is not large compared with
  \\( a / 0.3 \\), the decay length of the local loss, the total electronic
  stopping comes out somewhat below the LS value.

## Validity

The same as [Lindhard-Scharff](stopping-lindhard-scharff.md): velocities
below \\( v_0 Z_1^{2/3} \\).

## Verification status

The constants \\( c = 0.3 \\) and the choice of the Firsov length were
entered from the contributor's recollection of the paper and are **not
verified** against it; spot-check them before relying on the local loss
([`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)).
The averaging identity with LS is tested.

## References

- O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132, 647 (1976).
- O. B. Firsov, Sov. Phys. JETP 6, 534 (1958) (screening length).
- J. Lindhard and M. Scharff, Phys. Rev. 124, 128 (1961).
