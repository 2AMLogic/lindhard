# Cu–Ni alloy measurements

These source tables make the published numbers usable without access to the
papers. They are reference inputs. The single-element sputter-yield harness does
not load this directory; the Cu–Ni comparison with the fluence-dependent engine
(#448) is `validation/experiments/cuni_dynamic.py`, with its results and verdicts
in [the validation plan](../../../docs/validation.md#fluence-dependent-targets-vs-measurement-437-448-cuni-compared-not-validated).
That comparison does not validate the engine: the verdicts there are failures
and some conversions remain unresolved.

## Ho et al. (1976)

P. S. Ho, J. E. Lewis, H. S. Wildman and J. K. Howard, “Auger study of preferred
sputtering on binary alloy surfaces”, *Surface Science* **57**, 393–405 (1976),
[doi:10.1016/0039-6028(76)90191-6](https://doi.org/10.1016/0039-6028(76)90191-6).

| File | Printed location | Coverage |
|---|---|---|
| [ho1976_table1.csv](ho1976_table1.csv) | Table 1, p. 401 (PDF p. 9) | Five Cu–Ni compositions at 0.5 keV Ar+ |
| [ho1976_table2.csv](ho1976_table2.csv) | Table 2, p. 402 (PDF p. 10) | Four energies at 38.4 at% Cu |
| [ho1976.json](ho1976.json) | Source metadata | Conditions, column meanings, extraction checks and limitations |

All nine rows are transcribed as printed. Energy in Table 1 comes from its
footnote; composition in Table 2 comes from its caption. `reported_pm` columns
hold the printed ± amounts, **not an assumed standard deviation**. Intensities
are ions cm⁻² s⁻¹. Erosion rates are Å/min, not atoms/ion. `r0_over_rinf` and
`tau_s` are the authors' fitted/extrapolated AES results, rather than raw
detector readings. The numerical facts are stored with attribution; the paper,
page images and article text are not redistributed.

### Conditions needed to use the tables

Section 3 (pp. 397–398) describes homogeneous, coevaporated Cu–Ni films about
5000 Å thick on thermally grown SiO₂. Their grains have random orientation.
The ion gun is approximately **70° from the sample surface, or 20° from its
normal**. The deposition substrate was at room temperature; an irradiation
temperature is not specified in that procedure. Argon pressure is 5 × 10⁻⁵ Torr.
The authors measured ion intensity with a Faraday cup and crater depth using
Talystep and laser interferometry, favoring the latter.

The initial oxide is about 20 Å thick and changes the transient behavior.
The authors fit the later, oxide-cleared portion of the AES signal and
extrapolate to obtain `R0`; `Rinf` is the equilibrium Cu/Ni signal ratio.
Cu 920 eV and Ni 716 eV are the preferred AES lines; the Ni signal requires
correction for Cu overlap. Composition calibration and extrapolation are
described on pp. 398–402. Treat these as dependencies of the measurement.
AES samples a finite depth, so its inferred composition is not automatically
the composition of the first simulation bin.

### Printed inconsistency and selection

At 38.4 at% Cu and 0.5 keV, Table 1 (`t1-r3`) prints `R0/Rinf = 1.7 ± 0.2`
and `tau = 141 ± 7 s`; Table 2 (`t2-r1`) prints `1.6 ± 0.2` and `92 ± 8 s`.
The latter pair also appears in Table 1's 71.9 at% row. Both tables have been
checked visually against the printed pages: this is not an OCR correction.
Preserve their separate row identities; do not average them or silently pick
one as a correction. Exclude the ambiguous 0.5 keV / 38.4 at% pair from an
initial benchmark until resolved.

Table 2's 2.0 keV row (`t2-r4`) is a useful initial candidate. The discussion
(pp. 403–404) also cautions that the AES escape depth can exceed the altered
layer at 0.5 keV; at energies of 1 keV and above the authors consider the
altered layer adequately thick for this measurement.

### Derived quantities: keep the assumptions explicit

With the authors' proportional AES calibration, bulk Cu fraction `c`, and
`q = R0/Rinf`, the inferred equilibrium Cu fraction is

```text
m = c / (q * (1 - c) + c)
```

This follows by taking `R0 ∝ c/(1-c)` and `Rinf ∝ m/(1-m)` with the same
calibration factor. It depends on the oxide-cleared extrapolation, calibration
and sampling depth; it is **not a separately measured table column**. For
`t2-r4`, the formula gives approximately 0.2683. This number is our calculation,
not a value printed in the paper. The CSVs retain only the source quantities;
propagate the calibration and composition uncertainties before comparing.

Total yield would be `Y = N * V / J` for alloy number density `N` in atoms/cm³,
erosion speed `V` in cm/s (`Å/min × 1e-8 / 60`), and ion flux `J` per unit
sample area. Establish the density and the area convention of the reported
ion intensity before making that conversion; do not silently apply an angle
correction. No total-yield column has been invented. The individual alloy
yields and altered-layer thickness in Figs. 4–5 are inferred through the
authors' kinetic model and are not independent measured yield fixtures.

## Betz (1980): total yields

G. Betz, “Alloy sputtering”, *Surface Science* **92**, 283–309 (1980),
[doi:10.1016/0039-6028(80)90258-7](https://doi.org/10.1016/0039-6028(80)90258-7),
has measured Cu–Ni total yields in Fig. 13, p. 304 (PDF p. 22). All **22
visible measurement markers** are in [betz1980_fig13.csv](betz1980_fig13.csv),
with source metadata in [betz1980.json](betz1980.json). The fitted line and
extrapolated pure-element endpoints are excluded.

The first reading was an operator-supplied approximate coordinate table. A
second manual reading used a separate 300 dpi rendering of the original PDF.
Both readings, the pixel centres and axis calibration are retained in
[the digitization record](../digitize/betz1980_fig13.json). The second reading
sets the CSV values, rounded to 0.1 at% Cu and 0.01 atoms/ion. It is a manual
marker read, not an automated detection result. Reproduce the CSV without
needing the paper or any image package:

```sh
python3 validation/data/digitize/betz1980_fig13.py --check
# Omit --check to regenerate the CSV from the recorded coordinates.
```

Estimated digitization bounds are ±1 at% Cu and ±0.03 atoms/ion for rows 1–15,
and ±1.5 at% Cu and ±0.05 atoms/ion for crowded rows 16–22. All points agree
between readings within the sum of their respective digitization bounds.
These are reading bounds, **not experimental uncertainties or 1-sigma errors**.
The paper gives no pointwise error bars or complete uncertainty budget for
Fig. 13. Its 3% beam-current stability is only one contribution; full
experimental uncertainty remains unknown and is represented as null in metadata.

The 2 keV Ar+ experiment uses **normal incidence**, unlike Ho's 20° from
normal. Cu–Ni samples are about 2000 Å thick, sputter-deposited on glass.
Section 3 describes inventory divided by fluence to the film/substrate
interface, located at a 50% AES signal decrease. Section 3.2 states a different
preparation route for Cu–Ni than for the evaporated alloys and does not detail
its separate inventory calibration. A steady-state simulation comparison
needs to account for this film-integrated measurement. Do not combine Betz
and Ho into one matched experimental row.

## Related papers already examined

G. Betz and G. K. Wehner, “Sputtering of Multicomponent Materials”, in
R. Behrisch (ed.), *Sputtering by Particle Bombardment II*, Topics in Applied
Physics **52**, Springer (1983), pp. 11–90, reproduces these earlier results
(Figs. 2.10 and 2.13(b), pp. 35 and 43). It is a source map, not an independent
replicate. The supplied 1982 Betz book review has no experimental table.
V. I. Shapovalov, *Materials* **16**, 3258 (2023),
[doi:10.3390/ma16083258](https://doi.org/10.3390/ma16083258), concerns reactive
sputtering; its chemistry differs from this inert-ion Cu–Ni comparison.

These alloy tables do not supply a high-fluence implant depth profile. The Ho
`t2-r4` and Betz comparisons for both volume-relaxation conventions and the
fixed tolerances are in
[the validation plan](../../../docs/validation.md#fluence-dependent-targets-vs-measurement-437-448-cuni-compared-not-validated);
the other Ho rows have not been compared.
