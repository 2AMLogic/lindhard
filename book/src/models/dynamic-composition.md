# Dynamic composition (bookkeeping only)

Code: `lindhard/src/ion/dynamic.rs` (`CompositionGrid`, `Relaxation`).

**Status.** This is the start of milestone M3 (targets that change with
fluence). Only the bookkeeping exists: a grid of slabs whose composition
callers change with explicit inventory deltas. The fluence loop, and the
adapter that turns transport tallies into deltas, are not written yet, so no
CLI run uses it and there is no TOML key.

## Model

The finite layers of a target are held as slabs. Each slab stores an
**areal inventory** \\( A_i \\) (atoms/m²) per element, kept sorted by
\\( Z \\) so every derived quantity is deterministic; a slab built from a
layer of atom density \\( n_i \\) and thickness \\( t \\) starts with
\\( A_i = n_i t \\). An optional semi-infinite substrate backs the grid and
never changes.

Inventories change only through `CompositionGrid::apply`. A retained primary
adds one atom of its species where it came to rest; a recoil subtracts an
atom where it was created and adds it where it stops (an escaping recoil is
a loss). The update either succeeds completely or leaves the grid untouched.

### Volume relaxation

After every update each slab's thickness is recomputed from its inventory by
one of two stated conventions:

| Convention | Rust | Thickness |
|---|---|---|
| Ideal mixing of atomic volumes | `Relaxation::IdealMixing` | \\( t = \sum_i A_i v_i \\), with reference atomic volumes \\( v_i \\) |
| Fixed total number density | `Relaxation::FixedNumberDensity` | \\( t = \sum_i A_i / n_\mathrm{mix} \\) |

For an element in its own solid, `atomic_volume_from_density` gives
\\( v = M / (N_A \rho) \\). After relaxation the slab boundaries are rebuilt
from the front surface, which stays at \\( x = 0 \\): erosion or swelling
moves the interior interfaces and the back face.

## Assumptions

- Ideal mixing is additivity of atomic volumes (a Vegard-type rule applied to
  atoms); it ignores chemistry, voids and amorphisation swelling, and a
  compound's real density generally differs from its prediction.
- A fixed number density reproduces a known compound density but cannot
  respond to composition, and one value serves every slab.
- Neither is an equation of state with a pressure or phase model from the
  literature; both are **this crate's stated conventions**, and the module
  claims no published source for them.

## Validity

Bookkeeping for compositions that stay close to the phase whose volumes or
density the caller supplies. The caller must give the volume of every
species that occurs (none is inferred, in particular not for gases).

## Verification status

Conservation, transactional updates and the relaxation conventions are
covered by `lindhard/tests/dynamic.rs`. There is no comparison with a
measured dynamic profile yet, because no fluence driver exists.

## References

No published source is claimed for the two relaxation conventions (see
above). The atomic volume uses the Avogadro constant of CODATA 2022 (P. J.
Mohr, D. B. Newell, B. N. Taylor, E. Tiesinga, Rev. Mod. Phys. 97, 025002
(2025)) and the element densities of `lindhard/src/elements.rs`, whose
sources are in
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md).
