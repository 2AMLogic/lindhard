# Validation plan

A model counts as done only when it passes the checks below, at three levels.
The tolerances are starting points that the M0 issues will refine.

## 1. Analytic and internal checks (CI, every PR)

- Integrating the scattering integral by quadrature reproduces the ZBL
  universal nuclear stopping formula to < 0.5%.
- The magic formula and quadrature agree within the magic formula's published
  accuracy.
- Reduced range ρ(ε) agrees with LSS theory in the regime where nuclear
  stopping dominates.
- Energy is conserved per history: deposited + escaped + bound equals the
  incident energy.
- Determinism: same seed gives bit-identical tallies on 1, 2 and N threads.
- Electrons: elastic total cross sections from our partial-wave Mott solution
  match published Mott values at spot energies. Inelastic: the dielectric model
  satisfies the f-sum rule and the perfect-screening sum rule.

## 2. Code-to-code oracles (local harness, summaries committed)

Programs are run unmodified on matched problems, under the rules in
[`../CONTRIBUTING.md`](../CONTRIBUTING.md) §Oracles:

| Oracle | Compared on |
|---|---|
| RustBCA | Amorphous ranges, sputter yields, reflection coefficients, ions/s |
| OpenTRIM | Ranges, damage profiles, ions/s |
| iradina / IM3D | 3D geometry cases (M3) |
| SRIM (via Wine, local only) | Range moments only, as a familiarity check; never a fixture |
| Nebula | Electron PSFs, SE/BSE yields (M1) |
| Geant4 MicroElec | Electron energy deposition in Si/SiO₂ (M1) |
| PENELOPE / penEasy | Electron depth–dose above ~1 keV (M1) |

## 3. Experiment (the real bar)

- **Ion ranges:** published SIMS and RBS depth profiles of B, BF₂, P, As and
  Sb in amorphous (pre-amorphized) Si, and in Ge and SiC. Cite each dataset in
  [`data-provenance.md`](data-provenance.md).
- **Channeling (M2):** published SIMS profiles in crystalline Si as a function
  of tilt and twist, dose and screen oxide.
- **Electronic stopping:** IAEA stopping database experimental points, with
  per-system residual statistics reported.
- **Sputtering:** published yields (e.g. Ar → Si, Ar → Cu) vs. energy and angle.
- **Electrons (M1):** published backscatter coefficients η(E, Z), SE yields
  δ(E), and resist-exposure PSF measurements.

## Reporting

`validation/` will hold the harness, and its committed output is a table of
summary metrics per release. A release note states which validation set it
passes.
