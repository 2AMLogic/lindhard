# High-fluence implant depth profile

One measured implant depth profile at a fluence high enough that sputtering
changes its shape. Reference input only; no engine comparison has been run.

## Titze et al. (2021), Fig. 3, 50 keV Sb into Si at 5×10¹⁶ ions/cm²

M. Titze, J. L. Pacheco, T. Byers, S. B. Van Deusen, D. L. Perry, D. Weathers and
E. S. Bielejec, Sandia report SAND2021-13560J (author manuscript, 25 October 2021),
public full text at <https://www.osti.gov/servlets/purl/1830531> (accessed 2026-10-24).

| File | Printed location | Coverage |
|---|---|---|
| [titze2021_fig3_5e16.csv](titze2021_fig3_5e16.csv) | Fig. 3, p. 3 | 53 SIMS marker centres, depth in nm and normalized Sb concentration |
| [titze2021.json](titze2021.json) | Source metadata | Conditions, printed text values, limits |

The digitization record (two reads, calibration) is
[`../digitize/titze2021_fig3.json`](../digitize/titze2021_fig3.json); reproduce the
CSV and the cross-check with `python3 validation/data/digitize/titze2021_fig3.py --check`.

The concentration is normalized to the curve maximum by the authors. The profile
gives shape and depth, not retained dose. The authors report that the 5×10¹⁶
profile peaks near 15 nm against about 27.5 nm below their 6×10¹⁵ ions/cm²
saturation threshold, with a skewed shape. Digitization bounds are 0.3 nm and
0.02 (normalized), estimated, not statistical. Markers hidden behind other
series (about 10 to 19 nm) are not in the file. The paper's SRIM and TRIDYN
curves are code output and are not stored. The paper and its images are not
redistributed.
