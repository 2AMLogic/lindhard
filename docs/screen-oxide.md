# Screen-oxide characterization

A numerical characterization of the implemented model (#280): how an amorphous
SiO2 overlayer changes what reaches a crystalline Si substrate, and the
channeling tail there. It is **not** an experimental validation. Agreement
with measured SIMS profiles is a separate deliverable
([`validation.md`](validation.md), section 3).

## Reproduce

```sh
LINDHARD_SCREEN_OXIDE_N=3000 \
LINDHARD_SCREEN_OXIDE_OUT=validation/results/screen_oxide_baseline.json \
  cargo test --release -p lindhard --test crystal_screen_oxide -- --ignored --nocapture
```

The sweep (`lindhard/tests/crystal_screen_oxide.rs`, `screen_oxide_sweep`) is
`#[ignore]`d like the other statistical crystal tests. It takes about 10
minutes on one core at 3000 primaries per configuration (40 configurations),
prints one summary line per configuration and writes the JSON. The committed
`validation/results/screen_oxide_baseline.json` is its output at the settings
below; the header of the file records the model and inputs, and every
crystalline record carries `Bca::crystal_metadata`. The ordinary (non-ignored)
test `small_screen_oxide_case_conserves_energy_and_is_thread_independent`
runs 150 histories per case and checks energy conservation, that every primary
falls in exactly one class, and bit-identical tallies on 1, 2 and 8 threads
(same platform), for the crystal and for its amorphous control.

## Setup

| Item | Choice |
|---|---|
| Ions | B 5 keV; As 30 keV |
| Directions | aligned: tilt 45 deg, twist 0 (<110> of a (100) wafer, reference [010]); off-axis control: tilt 30 deg, twist 17 deg |
| Stack | SiO2 of T = 0 (bare), 1, 2, 5, 10 nm on Si; SiO2 amorphous, 2200 kg/m3, E_d = 15 eV |
| Substrate | crystalline Si (crystal flight model, 300 K, Debye 640 K with zero point) and, as the matched control, amorphous Si at the crystal density in the same stack |
| Electronic loss | `ElectronicLoss::NonLocal` (engine default), Lindhard-Scharff. No local loss; the crystal local/nonlocal ratio deviation of #250 is neither addressed nor changed here |
| Statistics | 3000 primaries per configuration, seed 1 for all, recoils followed, cutoffs 5 eV / 2 eV |

Coordinates: `x = 0` is the entrance surface (front of the oxide, or of the
bare crystal); the interface is at `x = T`. Results are given as "from
surface" and "from interface" (`x - T`). The bare crystal and the bare
amorphous control use the same ion, direction and seed. "Entry" into the
substrate is the first inward arrival at `x = T` (the incident beam for
T = 0); energy and angles are taken there. Histories that stop in the oxide,
or leave through the front before reaching the substrate, are reported
separately (`stopped_in_oxide`, `backscattered_before_entry`) and are in no
substrate statistic, so a lost transmitted history is not read as
dechanneling. Histories that enter and then return to the oxide are counted
(`entered_then_returned_to_oxide`, `entered_then_backscattered`) and are also
excluded from the substrate depth statistics.

The tail metric is the fraction of substrate-stopped primaries deeper than 2
(and 3) times the bare-amorphous Rp of the same ion and direction, measured
below the interface; `fraction_of_incident` gives the same count over all
incident primaries. Uncertainties are one standard error, `sigma/sqrt(n)` for
means and `sqrt(f(1-f)/n)` for fractions. Seeds are shared between
configurations, but the histories diverge in the oxide, so differences are
compared in quadrature as if independent.

## Measured baselines (3000 primaries)

Aligned <110>, tail fraction deeper than 2 Rp(amorphous, bare) below the
interface, crystal vs amorphous control (one standard error):

| Ion | T (nm) | enter | E at entry (eV) | angle to beam (deg) | crystal tail | amorphous tail |
|---|---|---|---|---|---|---|
| B 5 keV | 0 | 1.000 | 5000 | 0 | 0.833 +- 0.007 | 0.056 +- 0.005 |
| B 5 keV | 1 | 0.989 | 4798 | 7.6 | 0.549 +- 0.009 | 0.048 +- 0.004 |
| B 5 keV | 2 | 0.974 | 4581 | 12.4 | 0.384 +- 0.009 | 0.039 +- 0.004 |
| B 5 keV | 5 | 0.917 | 3970 | 21.7 | 0.173 +- 0.008 | 0.020 +- 0.003 |
| B 5 keV | 10 | 0.764 | 3073 | 31.1 | 0.083 +- 0.006 | 0.009 +- 0.002 |
| As 30 keV | 0 | 1.000 | 30000 | 0 | 0.862 +- 0.006 | 0.021 +- 0.003 |
| As 30 keV | 1 | 1.000 | 28014 | 4.6 | 0.622 +- 0.009 | 0.016 +- 0.002 |
| As 30 keV | 2 | 1.000 | 25999 | 7.3 | 0.456 +- 0.009 | 0.013 +- 0.002 |
| As 30 keV | 5 | 0.986 | 20414 | 13.0 | 0.230 +- 0.008 | 0.005 +- 0.001 |
| As 30 keV | 10 | 0.897 | 12857 | 19.8 | 0.093 +- 0.006 | 0.000 +- 0.000 |

Off-axis control (30 deg / 17 deg), crystal tail vs amorphous control:

| Ion | T (nm) | crystal tail | amorphous tail |
|---|---|---|---|
| B 5 keV | 0 | 0.069 +- 0.005 | 0.041 +- 0.004 |
| B 5 keV | 2 | 0.085 +- 0.005 | 0.033 +- 0.004 |
| B 5 keV | 10 | 0.068 +- 0.005 | 0.010 +- 0.002 |
| As 30 keV | 0 | 0.028 +- 0.003 | 0.014 +- 0.002 |
| As 30 keV | 2 | 0.051 +- 0.004 | 0.008 +- 0.002 |
| As 30 keV | 10 | 0.050 +- 0.004 | 0.000 +- 0.000 |

All the 40 records, including projected range from surface and interface,
entry-angle statistics to the surface normal and to the beam, and the 3 Rp
tail, are in the JSON. Energy was conserved in every history (relative
residual < 1e-9).

## What the results support

- In this model a thin amorphous screen reduces the aligned-direction
  substrate tail: in these runs the crystal tail fell from 0.83 (B) and 0.86
  (As) bare to 0.08 and 0.09 at 10 nm, beyond the statistical errors. The
  decrease was observed, not asserted: the tests assert only conservation,
  accounting and determinism.
- The suppression comes with a large angular broadening at entry (mean angle
  to the incident axis 7.6 deg at 1 nm to 31 deg at 10 nm for B) and with
  entry-energy loss (B 5000 to 3073 eV at 10 nm). Substrate entry is reduced
  only by oxide stopping and reflection (B: 0.76 enter at 10 nm; As: 0.90).
  The tail fractions above are conditional on entry, so they are not
  explained by lost histories.
- The substrate tail does not reach the amorphous control at 10 nm: aligned
  crystal tails remain about ten times (B) the control or well above it (As),
  as the beam keeps some collimation after the oxide.

## Deviations and caveats to report

- Off-axis, a screen does **not** reduce the crystal tail: the crystal tail
  is above the matched amorphous control at every thickness, and for B at
  1 to 5 nm it is larger than that of the bare crystal (for example B 2 nm:
  0.085 +- 0.005 against 0.069 +- 0.005 bare, about 2 standard errors of
  the difference; As 2 nm 0.051 +- 0.004 against 0.028 +- 0.003). The oxide
  moves a nominally random beam into a distribution of angles that includes
  channeling directions. No fault is claimed; it follows from the model and
  is a candidate for separate investigation (not made here).
- The tail threshold is fixed by the bare amorphous control, but the beam
  arrives with less energy after the oxide, so the amorphous control's own
  tail also falls with T. Read the crystal tail against the matched control
  of the same row, not against the bare row.
- Electronic loss is the nonlocal Lindhard-Scharff default with no
  equipartition local term. The known crystal local/nonlocal loss-ratio
  deviation (#250) is unchanged, so absolute tail fractions and ranges, not
  only their trends, depend on that choice. The sweep does not vary it.
- Oxide stopping uses the same amorphous Lindhard-Scharff model and a
  nominal 2200 kg/m3, with no dependence of the electronic stopping on the
  oxide's O content beyond the Bragg sum. No oxide-specific stopping data
  were used or fitted.
- The Debye temperature is the one cited default of Si. No beam divergence,
  no dose, no accumulating damage (a static oxide, fresh crystal), no
  thermal expansion and no interface roughness.
- Statistics: 3000 primaries and one seed. Tail uncertainties of a few
  tenths of a percent are statistical only; they do not include model or
  seed-to-seed variation.

## Needs experimental validation

Everything above is a property of the implemented model. None of it has been
compared with measurement. In particular, the oxide thickness at which a
real tail is suppressed, the magnitude of suppression, the off-axis
tail enhancement and the sensitivity to electronic loss need low-dose SIMS
profiles with the tilt, twist, reference direction, divergence, oxide
thickness and dose all stated (see [`validation.md`](validation.md),
"Channeling (M2)"). Dose-dependent models should not be calibrated before
that comparison.
