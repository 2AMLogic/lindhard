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
`"none"` keeps the deterministic mean loss and draws nothing. Each free
flight's nonlocal loss is then drawn from a Gaussian centred on the mean loss
it would take anyway, with standard deviation
`sqrt(n Omega^2 s)`: the Bragg-added variance above
(`variance_per_atom_material`) scaled by the flight's areal atom density
`n s`. The draw comes from the history's own random stream, so results stay
deterministic per `(seed, index)` at any thread count, and a run with the
option on is a different trajectory sequence than a run without it:
comparisons between the two must be statistical over many histories, never
per history.

The Gaussian is unbounded below, so a draw can exceed the ion's energy or
come out negative. The sampled loss is truncated by clamping to `[0, E]`,
which biases the mean loss upward wherever the per-flight standard deviation
is comparable to the per-flight mean loss — with the constant free path
(`lambda = N^(-1/3)`) at MeV energies it always is, and the sampled variance
of a transmitted spectrum then falls well below the Bohr value (measured
about 0.42 for a 1 MeV proton in a 100 nm carbon foil; see the
`ion::stopping::straggling` module docs). With long flights — the
energy-dependent free path at a small `min_cm_angle_deg` — the truncation is
inactive and the transmitted variance matches Bohr to about two percent in
the same foil. The library entry points
(`variance_per_atom`, `variance_per_atom_material`) remain available
independently.

## Assumptions

- Free, stationary target electrons; every electron of the target
  contributes, whatever its binding.

## Validity

High energy, \\( v \gg v_0 Z_1^{2/3} \\). Below about 1 MeV/u Bohr's value
is an overestimate, because the bound electrons contribute less than free
ones. The corrections that reduce it at intermediate energy (Chu;
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
