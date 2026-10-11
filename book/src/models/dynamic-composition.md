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
from the front surface at \\( x = 0 \\). Swelling moves the interior
interfaces and the back face. With sputter erosion on, the sputtered atoms
of each element are removed from the front slabs (the loss they would
otherwise cause in the slab where they were displaced is cancelled, so
nothing is removed twice), the thickness they occupied under the chosen
convention is the recession of that step, and the grid is re-anchored so the
surface is again \\( x = 0 \\); a depth plus the cumulative recession is
the depth in the original frame. Erosion is off by default and then changes
nothing.

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
covered by `lindhard/tests/dynamic.rs`. The fluence stepping loop
(`DynamicRun`, with the tally-to-delta adapter) is covered by
`lindhard/tests/dynamic_run.rs`: the low-fluence limit equals the static
engine, refining the step size converges, steps continue one global random
stream, results are bit-identical on 1, 2 and 8 threads, adaptive rejection
consumes no indices, and inventory follows the event conventions. There is no
comparison with a measured dynamic profile yet.

## The fluence loop

`DynamicRun` delivers the beam's primaries in steps. A step runs `n` primaries
on the current target, scales the integer atom counts of the events by the
fluence one primary stands for, applies them to the grid and relaxes it. The
target is held fixed within a step, so the step must be small enough that the
composition does not change much inside it: the adaptive policy bounds the
largest relative change of a slab per step and retries a too-large step with
fewer ions from the same first primary. The front surface stays at `x = 0`;
the time series reports the interface depths and total thickness measured
from it, not a receding surface. The step-size policy is this crate's own
design, not a published scheme. In the CLI it is the `[dynamic]` table
(`docs/cli.md`).

## References

No published source is claimed for the two relaxation conventions (see
above). The atomic volume uses the Avogadro constant of CODATA 2022 (P. J.
Mohr, D. B. Newell, B. N. Taylor, E. Tiesinga, Rev. Mod. Phys. 97, 025002
(2025)) and the element densities of `lindhard/src/elements.rs`, whose
sources are in
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md).
