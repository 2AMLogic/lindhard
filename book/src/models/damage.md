# Displacement damage

Code: `lindhard/src/ion/damage.rs` (models), `lindhard/src/tally/` (counts).

lindhard reports two different quantities side by side and never adds one
to the other or substitutes one for the other:

- **Model estimates** of the number of displacements made by a primary
  knock-on atom (PKA) of a given energy: NRT and Kinchin-Pease, from the
  Lindhard partition of the PKA energy (`results.damage.nrt`).
- **Event counts** from the full-cascade simulation: vacancies,
  interstitials and replacements counted collision by collision as the BCA
  follows every recoil (`results.damage.cascade`).

## Damage energy (Lindhard partition)

Of a PKA's energy \\( T \\), the part eventually given to atomic motion is
(Lindhard, Nielsen, Scharff and Thomsen 1963)

\\[ T_\mathrm{dam} = \frac{T}{1 + k\\, g(\varepsilon)}, \qquad
g(\varepsilon) = 3.4008\\, \varepsilon^{1/6} + 0.40244\\, \varepsilon^{3/4} + \varepsilon, \\]

with the analytic fit \\( g \\) of Robinson (1970) as adopted by the NRT
standard, and

\\[ \varepsilon = \frac{T a A_2}{Z_1 Z_2 e^2 (A_1 + A_2)}, \qquad
a = 0.8853\\, a_0 \left(Z_1^{2/3} + Z_2^{2/3}\right)^{-1/2}, \\]

\\[ k = \frac{0.0793\\, Z_1^{2/3} Z_2^{1/2} (A_1 + A_2)^{3/2}}
{\left(Z_1^{2/3} + Z_2^{2/3}\right)^{3/4} A_1^{3/2} A_2^{1/2}}, \\]

the Lindhard reduced energy and the
[Lindhard-Scharff](stopping-lindhard-scharff.md) electronic stopping
coefficient of the PKA (\\( Z_1, A_1 \\)) in the target (\\( Z_2, A_2 \\)).
For a self-ion these reduce to the forms written in the NRT standard,
\\( \varepsilon = T / (86.931\\, Z^{7/3}) \\) (\\( T \\) in eV) and
\\( k = 0.1337\\, Z^{1/6} (Z/A)^{1/2} \\); a test checks this.

## NRT and Kinchin-Pease

Norgett, Robinson and Torrens (1975):

\\[ N_\mathrm{NRT} = \begin{cases}
0 & T_\mathrm{dam} < E_d, \\\\
1 & E_d \le T_\mathrm{dam} < 2E_d/0.8, \\\\
0.8\\, T_\mathrm{dam} / (2 E_d) & T_\mathrm{dam} \ge 2E_d/0.8,
\end{cases} \\]

with \\( E_d \\) the displacement threshold and 0.8 the displacement
efficiency the authors took from BCA simulations. The older Kinchin and
Pease (1955) count is 0, 1, or \\( T/(2E_d) \\) with thresholds \\( E_d \\)
and \\( 2E_d \\); the tallies evaluate it on the damage energy (the usual
"modified" form).

**NRT is not a cascade count.** Stoller et al. (2013) discuss how BCA
full-cascade vacancy counts and NRT estimates differ, and recommend that a
displacement dose comparable with the NRT standard be computed from the
damage energy with the NRT formula, not taken from the full-cascade vacancy
count. NRT also overestimates the number of defects that survive in-cascade
recombination; it is a standard exposure unit, not a prediction of surviving
defects.

## Cascade counts

During transport a target atom is **displaced** when it receives
\\( T > E_d \\) ([BCA transport](bca.md), step 4). The counts are defined in
`lindhard::tally::CascadeDefects`:

- **Replacement**: a moving atom comes to rest, immediately after the
  collision in which it displaced an atom of its own element, at that atom's
  site, so it fills the site. A particle stops when its energy falls below
  its cutoff, so this count depends on the cutoffs (a recoil cutoff near
  \\( E_d \\) gives the most replacements); with `follow_recoils = false` it
  is always 0.
- **Vacancies** = displacements minus replacements.
- **Interstitials**: recoils that came to rest in the target and did not
  fill a site. Implanted beam particles are not counted here (they are in
  the range tally).

These are kept in total, per layer and as a depth profile
(`damage_profile.csv`).

## Selecting it

Always on: both are computed in every run. The inputs are the
per-element energies: \\( E_d \\) (`e_d_ev`), \\( E_b \\) (`e_b_ev`) and
\\( E_s \\) (`e_s_ev`), set per material or for an element in every layer
with `[physics.energies.<symbol>]`. There are no model variants to choose.

## Assumptions

- **Compounds (this crate's own convention).** NRT and the Lindhard partition
  are defined for a monatomic target. For a layer with several elements the
  tally uses the PKA's own \\( Z_1, A_1 \\), the atom-fraction-weighted mean
  \\( Z_2, A_2 \\) of the layer the PKA starts in (non-integer in general),
  and the PKA element's own \\( E_d \\) in that layer. No published source is
  claimed for this averaging; it reduces exactly to standard NRT for a
  monatomic layer. Treat NRT numbers for compounds as an exposure index
  comparable within lindhard, not as a value following any compound-target
  standard.
- The cascade counts depend directly on \\( E_d \\), \\( E_b \\) and the
  recoil cutoff, and on the amorphous-target BCA itself (no recombination,
  no crystal structure).

## Validity

NRT is a convention for damage accounting, defined for monatomic targets.
The cascade counts are BCA counts of displacement events, not of surviving
defects.

## Verification status

The self-ion reductions of \\( \varepsilon \\) and \\( k \\) are checked in a
test. The default \\( E_d \\) values of the element table follow the ASTM
E521 convention and are not yet verified against the standard
([`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md));
Si has no default, so an input must choose one.

## References

- J. Lindhard, V. Nielsen, M. Scharff, P. V. Thomsen, "Integral equations
  governing radiation effects", Mat. Fys. Medd. Dan. Vid. Selsk. 33 (10)
  (1963).
- M. T. Robinson, in *Nuclear Fusion Reactors* (British Nuclear Energy
  Society, London, 1970), p. 364.
- M. J. Norgett, M. T. Robinson, I. M. Torrens, Nucl. Eng. Des. 33, 50
  (1975).
- G. H. Kinchin and R. S. Pease, Rep. Prog. Phys. 18, 1 (1955).
- R. E. Stoller, M. B. Toloczko, G. S. Was, A. G. Certain, S. Dwaraknath,
  F. A. Garner, Nucl. Instrum. Methods B 310, 75 (2013).
