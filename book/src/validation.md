# Validation

What lindhard has been compared with, against which published reference, and
how far off it is. The full record, with every row, tolerance and
investigation, is
[`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md);
this page summarises it. Where lindhard disagrees with a reference the
disagreement is reported, not hidden. Provenance of every number in the tree:
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md).

The numbers below are the level 1 results (`full` statistics, `lindhard`
0.0.1) as of 2026-11-01. Level 1 runs in CI on every pull request.

## Level 1: analytic limits and internal consistency

Each row compares the engine with an equation or limit that can be evaluated
independently. "Deviation" is the measured difference; "tolerance" is the
bound the check asserts (three standard errors plus solver uncertainty plus a
stated systematic allowance for Monte Carlo rows).

| Compared | Against | Deviation | Tolerance |
|---|---|---|---|
| ZBL nuclear stopping by quadrature, 1e-4 ≤ ε ≤ 0.3 | ZBL universal fit (Ziegler, Biersack, Littmark 1985, ch. 2) | 0.388 % | 0.5 % |
| Nuclear stopping by quadrature, 1e3 ≤ ε ≤ 1e6 | Impulse-approximation asymptote | 0.110 % | 0.2 % |
| Quadrature scattering angle, small deflection | Rutherford formula | 0.009 % | 0.2 % |
| Projected range, 50 keV P into Si, 100 000 ions | First-moment equation of Lindhard, Scharff and Schiøtt (LSS) | 0.009 % | 0.889 % |
| Projected range, 50 keV As into Si, 100 000 ions | LSS first-moment equation | 0.020 % | 0.849 % |
| Reduced path, 5 keV Si into Si, nuclear stopping only | LSS first-moment equation | 0.438 % | 0.914 % |
| Lindhard partition solver | Lindhard, Nielsen, Scharff and Thomsen (1963) analytic solution for power-law scattering | 0.206 % | 1.0 % |
| Per-history energy budget, 900 Ar histories | Energy conservation | 2.7e-15 | 1e-9 |
| Masked-implant lateral straggle | One-dimensional lateral straggle | 0.148 % | 2.179 % |

Results are byte-identical on 1, 2 and 8 threads (`engine.determinism`).

## Known deviations

These are measured and not asserted away:

- **ZBL universal fit, 0.3 < ε < 1000:** the fit runs up to 5.5 % low against
  the quadrature, because its high-energy branch drops the screening constant.
  The quadrature agrees with the exact asymptote to 0.1 %.
- **Range of light ions in a semi-infinite target:** 10 keV B into Si comes out
  3.97 % above the infinite-medium equation, because backscattered ions leave
  the target. A difference of definition, not of transport.
- **Range against measurement:** for B in amorphous Si the default stopping
  inputs (Lindhard–Scharff electronic, ZBL nuclear) give a projected range 24
  to 36 % too long. Because the transport reproduces the LSS equation for the
  same inputs, the offset comes from those inputs.
- **Sputter yield against displacement energy:** Ar at 1 keV on Si gives 0.30
  with E_d = 15 eV and 0.49 with E_d = E_s. Set E_d ≤ E_s when low-energy
  yields matter.
- **Dynamic composition:** Cu–Ni fluence-dependent targets are compared with
  measurement, and the engine does **not** reproduce the measured profile.
  They are compared, not validated.

## Levels 2 and 3

Level 2 compares with other codes through oracle adapters (summaries are
committed; the third-party outputs are not). Level 3 compares with
measurement: stopping data, sputter yields, backscatter coefficients,
secondary-electron yields and inelastic mean free paths. Each comparison
has its data source, residuals and verdict in
[`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md),
and the page also says where a comparison is missing. The
[residuals against stopping measurements](models/stopping-residuals.md) are
in the manual.
