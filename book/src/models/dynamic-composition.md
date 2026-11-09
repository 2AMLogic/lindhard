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

### Retention limit

Optionally, each implanted species `Z` gets a **saturation limit**: the largest
atomic fraction \\( f_Z \\) of a slab it may make up. Dynamic-composition codes
(TRIDYN, SDTrimSP; see References) describe a saturation concentration of the
implanted species above which further atoms are not retained and are released
again. The update here is this crate's own statement of that idea; no numeric
limit is built in, it is always given by the caller. Per accepted step, each
slab is looked at with the composition it would have after the step. With
\\( A_Z \\) of the limited species and \\( O \\) atoms/m² of everything else,
the slab may hold at most

\\[ A_\mathrm{max} = \frac{f_Z}{1 - f_Z}\, O \\]

(that is \\( A_Z / (A_Z + O) \le f_Z \\)), and the excess
\\( A_Z - A_\mathrm{max} \\) is **re-emitted**: removed from the slab through
the same inventory deltas as any other change, and it leaves the target. A slab
holding only the limited species loses all of it. With several limited species
in one slab they are met in ascending \\( Z \\), each seeing the amounts left by
the earlier ones, so the limit of a later species is met exactly and an earlier
one can end marginally below its own after a later removal.

The atom balance of each element over a step is
\\( \Delta A_Z^\mathrm{total} = \text{transport net} - \text{re-emitted} \\),
where the transport net is the retained primaries and the moved, sputtered and
transmitted recoils under the conventions above. The re-emission is a pure
function of the merged integer tally and the grid and draws no random numbers,
so results stay bit-identical at any thread count. The model is off by
default and then changes nothing (the CLI adds no columns either).

## Assumptions

- The retention limit is a bookkeeping cap on a slab's composition. It has no
  diffusion, no temperature or time dependence, no trapping chemistry and no
  record of where re-emitted atoms go; the saturation fraction is a free
  parameter, not a prediction. A limit reached within one step is applied
  to the step's result, so steps small enough that no slab jumps far past the
  limit give a better approximation (the adaptive policy sees the removal).

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
consumes no indices, and inventory follows the event conventions. The
retention limit is checked for the analytic plateau (a saturated slab holds
exactly \\( f_Z \\) and none exceeds it), for the atom balance
(implanted − sputtered − re-emitted = retained), for bit-identity at 1, 2 and 8
threads, and for leaving results unchanged when off or out of reach. There is no
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

- W. Möller and W. Eckstein, Nucl. Instrum. Methods B 2, 814 (1984), and
  report IPP 9/64 (1988): TRIDYN, dynamic composition changes. Read for the
  idea of a saturation concentration only; nothing was ported.
- A. Mutzke et al., *SDTrimSP Version 6.00*, report IPP 2019-02 (2019): the
  same, for the saturation of implanted species.

No published source is claimed for the two relaxation conventions or for the
exact form of the retention update (see above). The atomic volume uses the Avogadro constant of CODATA 2022 (P. J.
Mohr, D. B. Newell, B. N. Taylor, E. Tiesinga, Rev. Mod. Phys. 97, 025002
(2025)) and the element densities of `lindhard/src/elements.rs`, whose
sources are in
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md).
