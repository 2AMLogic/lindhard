# BCA transport and free-flight conventions

Code: `lindhard/src/ion/bca/` (`Bca`, `BcaConfig`, `MeanFreePath`,
`ElectronicLoss`; `kinematics.rs`).

## Model

The target is a stack of homogeneous, structureless (amorphous) layers, the
first starting at depth \\( x = 0 \\); a semi-infinite substrate may close
it. A moving atom alternates straight free flights and binary elastic
collisions with target atoms, losing energy to electrons along each flight.
This is the amorphous-target BCA of Biersack and Haggmark (1980), within the
general BCA framework of Robinson and Torrens (1974) and the treatment of
Eckstein (1991). One history:

1. The primary enters at the origin of the front face with the beam's
   direction.
2. **Free flight** of length \\( \tau\lambda \\) (\\( \lambda \\) and the
   distribution of \\( \tau \\) depend on the free-path convention below),
   with the nonlocal electronic loss \\( N S_e(E)\\, s \\) taken at the energy
   at the start of the segment (with [Bragg additivity](stopping-bragg.md)
   over the layer composition).
3. **Collision** with a partner drawn by stoichiometry, at impact parameter
   \\( p = p_\mathrm{max}\sqrt{R} \\) (uniform over a disc) and a uniform
   azimuth. The centre-of-mass angle comes from the
   [scattering table](scattering.md).
4. **Recoil.** If the transfer \\( T \\) exceeds the partner's displacement
   energy \\( E_d \\), the atom is displaced with energy \\( T - E_b \\) (the
   lattice binding \\( E_b \\) stays in the lattice); otherwise \\( T \\)
   stays in the lattice at the site (Biersack and Haggmark 1980; Eckstein
   1991). Displaced atoms are followed in turn as full cascades.
5. A particle **stops** when its energy falls below its cutoff, or **leaves**
   through the front or back face if it can overcome the surface barrier;
   otherwise it is reflected specularly back into the target.

### Kinematics

Classical elastic two-body kinematics (Goldstein 1980; Robinson and Torrens
1974), with \\( \mu = M_1/M_2 \\):

\\[ \tan\psi = \frac{\sin\theta}{\cos\theta + \mu}, \qquad
\varphi = \frac{\pi - \theta}{2}, \qquad
T = \gamma E \sin^2\frac{\theta}{2}, \qquad \gamma = \frac{4 M_1 M_2}{(M_1 + M_2)^2}, \\]

where \\( \psi \\) is the projectile's laboratory deflection and
\\( \varphi \\) the recoil's, on the opposite azimuth.

### Surface barrier

A planar barrier (Sigmund 1969; Eckstein 1991): a particle of energy
\\( E \\) reaching a face with direction cosine \\( c \\) to the normal
escapes if \\( E c^2 > E_s \\). Outside, its energy is \\( E - E_s \\), the
parallel momentum is unchanged, and the normal component satisfies
\\( E' c'^2 = E c^2 - E_s \\) (refraction away from the normal). Target
atoms use the surface binding energy \\( E_s \\) of the material at the
face; the beam species uses `primary_surface_binding_ev` (default 0).

## Free-path conventions

| TOML (`[physics]`) | Rust | Flight length | Impact parameter |
|---|---|---|---|
| `free_path = "constant"` (default) | `FreePathChoice::Constant`, `MeanFreePath::Constant` | fixed, \\( l = N^{-1/3} \\) | \\( p_\mathrm{max} = (\pi N^{2/3})^{-1/2} \\) |
| `free_path = "energy-dependent"` | `FreePathChoice::EnergyDependent`, `MeanFreePath::EnergyDependent` | exponential, mean \\( \lambda(E) \\) | \\( p_i(E) \\), capped at the constant \\( p_\mathrm{max} \\) |

### Constant

The flight length is the mean interatomic distance \\( l = N^{-1/3} \\) and
impact parameters are uniform over a disc of radius
\\( p_\mathrm{max} = (\pi N^{2/3})^{-1/2} \\), so
\\( N\pi p_\mathrm{max}^2 l = 1 \\): each flight sweeps exactly one atom's
worth of target (Biersack and Haggmark 1980). So that every primary does
not make its first collision at the same depth, the primary's first flight
is \\( R\\, l \\) with \\( R \\) uniform in \\( [0, 1) \\); recoils start at
an atom site and fly a full \\( l \\). Collisions with
\\( p > p_\mathrm{max} \\) are dropped unless weak collisions add them.

### Energy-dependent

Collisions that deflect by less than a minimum centre-of-mass angle
\\( \theta_\mathrm{min} \\) (`min_cm_angle_deg`, required with this
convention and only with it) are neglected. For each element \\( i \\),
\\( p_i(E) \\) is the impact parameter at which the angle equals
\\( \theta_\mathrm{min} \\), capped at the constant-convention
\\( p_\mathrm{max} \\). The mean free path is

\\[ \lambda = \frac{1}{N\pi \sum_i x_i p_i^2}, \\]

flight lengths are exponential with that mean, the partner is drawn with
probability \\( x_i p_i^2 / \sum_j x_j p_j^2 \\), and
\\( p = p_i\sqrt{R} \\). This is the standard cross-section cut-off
treatment of a Poisson collision process (Eckstein 1991). It reduces to the
constant convention, with exponential flight lengths, when every
\\( p_i \\) hits the cap. The nuclear loss of the neglected small-angle
collisions is dropped, so choose \\( \theta_\mathrm{min} \\) small.

### Layer boundaries

A flight is drawn as a dimensionless number of mean free paths \\( \tau \\).
When it reaches an interface it is truncated there, the electronic loss for
the truncated length is applied, the material is switched, and the flight
continues in the new layer with the unused part \\( \tau - s/\lambda \\) (it
is not redrawn). For exponential paths this is exact by memorylessness; for
the constant path it means that splitting one layer into two of the same
material changes nothing but rounding (a test checks this).

### Weak collisions (optional)

With the constant free path every flight ends in one collision with
\\( p \le p_\mathrm{max} \\), and the nuclear loss of collisions beyond
\\( p_\mathrm{max} \\) is dropped while the electronic loss of the flight is
charged in full. At cascade-tail energies that dropped part is large, so
the electronic share of a cascade comes out too high.
`weak_collisions = K` (0 to 3) adds the weak collisions of Möller and
Eckstein (1988): before the hard collision, \\( K \\) collisions with
partners at

\\[ p_k = p_\mathrm{max}\sqrt{k + R}, \qquad k = 1, \dots, K, \\]

each uniform over an annulus of area \\( \pi p_\mathrm{max}^2 \\), with its
own partner and azimuth. A weak collision deflects the particle and takes its
transfer, but never makes a recoil (its \\( T \\) stays in the lattice). A
weak collision whose partner would lie in front of the front surface is
skipped. Each collision is evaluated at the energy left after the previous
one, which differs from the report, where all collisions of a step use the
starting energy (see the `ion::bca` module docs for the measured effect).
Weak collisions are available with the constant free path only. **Default
0**, because there is no single published convention and the option changes
sputter yields substantially (see
[`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md)).

## Electronic loss

| Mode | Rust | Selected by |
|---|---|---|
| All loss continuous along the flight, from the chosen stopping model | `ElectronicLoss::NonLocal` | `stopping = "lindhard-scharff"`, `stopping = "bethe-bloch"`, or user tables |
| Half Lindhard-Scharff along the flight, half Oen-Robinson at each collision | `ElectronicLoss::EquipartitionLsOr` | `stopping = "equipartition-ls-or"` |

The local half is evaluated at the distance of closest approach of every
collision, weak ones included; see
[Oen-Robinson](stopping-oen-robinson.md).

## Cutoffs and energy bookkeeping

`primary_cutoff_ev` and `recoil_cutoff_ev` are required: they have no
defensible default. A recoil cutoff above the surface binding energies
suppresses sputtering, so keep it below the smallest \\( E_s \\) when
sputtering matters (Eckstein 1991). With `follow_recoils = false` a
displaced atom stops where it was created.

Every energy change is subtracted from the particle and added to exactly one
field of the energy budget (`results.energy_budget_ev_per_ion`), so each
history conserves energy up to rounding; the largest per-history relative
residual is reported. Electronic losses larger than the particle's energy
are clamped to it.

## Randomness and determinism

Each primary has its own random stream, keyed on the run seed and the
primary's index; its recoils draw from the same stream, and cascades are
followed from an explicit stack in a fixed order. Histories are grouped in
chunks of a fixed size (64, independent of the thread count), which fixes
the summation order of the tallies. The output is therefore bit-identical on
any number of threads (`lindhard/tests/determinism.rs`).

## Assumptions

- Amorphous, structureless targets: no crystal structure (no channeling).
- The target does not change with fluence (no dynamic composition).
- No refraction of the incident beam at the entrance surface (negligible at
  keV energies).
- Binary collisions only; many-body effects at low energy are represented
  only through the cutoffs, \\( E_d \\), \\( E_b \\) and \\( E_s \\).

## Validity

The BCA is meant for energies well above the binding energies of the
target. It follows cascade atoms down to a few eV, where its results depend on the cutoffs and binding
energies chosen. Comparisons with other codes and with measured ranges and
sputter yields, with the known deviations, are in
[`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md).

## References

- J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174, 257 (1980),
  doi:10.1016/0029-554X(80)90440-1.
- M. T. Robinson and I. M. Torrens, Phys. Rev. B 9, 5008 (1974).
- W. Eckstein, *Computer Simulation of Ion-Solid Interactions* (Springer,
  Berlin, 1991).
- H. Goldstein, *Classical Mechanics*, 2nd ed. (Addison-Wesley, 1980),
  sec. 3.11.
- P. Sigmund, Phys. Rev. 184, 383 (1969).
- W. Möller and W. Eckstein, *TRIDYN - Binary collision simulation of atomic
  collisions and dynamic composition changes in solids*, report IPP 9/64,
  Max-Planck-Institut für Plasmaphysik, Garching (1988), p. 14, eq. (26).
- O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132, 647 (1976).
