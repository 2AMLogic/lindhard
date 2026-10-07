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

Library only (`variance_per_atom`, `variance_per_atom_material`). The BCA
engine does not use it at present and there is no TOML key: the electronic
loss along each free flight is deterministic.

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
