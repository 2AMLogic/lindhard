# Energy-loss straggling

Code: `lindhard/src/ion/stopping/straggling.rs` (`StragglingModel`,
`variance_per_atom`).

## Model

The energy loss of an ion along a path fluctuates. For a thick target in the
free-electron picture, the variance of the energy loss per unit areal
density is (Bohr 1948)

\\[ \Omega_B^2 = 4\pi z^2 Z_2 (e^2)^2 \quad \text{(J}^2\\,\text{m}^2 \text{ per target atom)}, \\]

independent of energy in the non-relativistic regime, with \\( z \\) the
projectile charge from an [effective-charge model](stopping-bethe-bloch.md).
The engine computes it for the bare charge \\( z = Z_1 \\) regardless of the
mean loss's effective charge, so a `bethe-bloch` mean loss configured with a
non-bare charge samples against a bare-charge variance — a known
inconsistency between mean loss and variance, not a fitted or corrected
model.
The relativistic variant multiplies it by
\\( (1 - \beta^2/2)/(1 - \beta^2) \\) (Bethe and Livingston 1937; Fano 1963).
For a compound the variances add by atom fraction,
\\( \sum_j x_j \Omega_j^2 \\).

| Model | Rust |
|---|---|
| Bohr, non-relativistic (default) | `StragglingModel::Bohr` |
| Bohr with the relativistic factor | `StragglingModel::BohrRelativistic` |

## Selecting it

Ion BCA transport samples electronic straggling when the input sets
`[physics] straggling = "bohr"` (or `"bohr-relativistic"`); the default
`"none"` keeps the deterministic mean loss and draws nothing. Every flight
still loses its mean; the flight also adds its mean loss and its Bohr
variance `n Omega^2 s` to the particle's pending straggling segment. When
the segment's accumulated mean loss reaches 4 of its accumulated standard
deviations — or the particle stops below its cutoff, escapes or reflects —
one deviation is drawn from a Gaussian of the segment's accumulated
standard deviation and applied, clamped to `[-mu_acc, E]`. Because Bohr's
variance is linear in path length, the accumulated variance of a segment of
any number of flights is exactly the Bohr variance of that segment, so the
one draw reproduces the published variance per unit path length whatever
the flight convention; and at the closing criterion the Gaussian sits 4
standard deviations above zero, so the probability of a draw below it is
`Phi(-4) ~ 3.2e-5` and the clamp — which keeps every history's energy exact
and non-negative — is essentially inactive. The mean loss stays unbiased.
Because the bare-charge Bohr variance does not depend on energy, the
engine precomputes it per layer alongside its mean-loss coefficient (the
relativistic variant keeps its energy dependence and is evaluated per
flight).

The segment, not the flight, is the sampling unit because Bohr's variance
is a many-collisions result: a Gaussian describes the loss of a path that
sits many standard deviations above zero, which a single free flight is not
(the constant free path flies `lambda = N^(-1/3)`, about 0.21 nm in carbon;
a crystal region flies the lattice, `search_length_nm` defaulting to the
lattice constant, 0.543 nm in Si, with no long-flight option). An earlier
implementation drew per flight and truncated each draw to `[0, E]`, which
left-censored those short-flight Gaussians: the mean electronic loss
measured about 2.5x too high in amorphous carbon, 2.8x in amorphous silicon
and 3.7x on crystal Si lattice flights, with the transmitted variance
depressed to about 0.42 of the Bohr value. The per-segment cadence removes
that bias; the figures are history, not current behaviour. The deviation is
tallied on the flight that closes the segment — localized at the segment's
end rather than spread along it — and layer boundaries do not close a
segment, so thin layers do not reintroduce the short-segment regime.

The draw comes from the history's own random stream, so results stay
deterministic per `(seed, index)` at any thread count, and a run with the
option on is a different trajectory sequence than a run without it:
comparisons between the two must be statistical over many histories, never
per history. The library entry points (`variance_per_atom`,
`variance_per_atom_material`) remain available independently.

## Assumptions

- Free, stationary target electrons; every electron of the target
  contributes, whatever its binding.

## Validity

High energy, \\( v \gg v_0 Z_1^{2/3} \\). Below about 1 MeV/u Bohr's value
is an overestimate, because the bound electrons contribute less than free
ones; `lindhard check` warns when the beam energy sits below this floor. The corrections that reduce it at intermediate energy (Chu;
Yang-O'Connor-Wang) are declined, and the Lindhard-Scharff low-velocity
correction is omitted (not verified); see
[Validity ranges and declined terms](stopping-validity.md).

## Verification status

Not verified against the original papers
([`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)).

## References

- N. Bohr, Mat. Fys. Medd. Dan. Vid. Selsk. 18 (8) (1948).
- H. A. Bethe and M. S. Livingston, Rev. Mod. Phys. 9, 245 (1937).
- U. Fano, Ann. Rev. Nucl. Sci. 13, 1 (1963).
- W. K. Chu, Phys. Rev. A 13, 2057 (1976) (declined correction).
- Q. Yang, D. J. O'Connor, Z. Wang, Nucl. Instrum. Methods B 61, 149 (1991)
  (declined correction).
