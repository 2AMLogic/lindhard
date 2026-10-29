# srim.org stopping-plot digitizer

Reads the measured points off srim.org's per-system stopping plots
(`StopPlots/STOPzzZZ.gif`, ion Z then target Z) and writes
[`../../stopping/compilation/`](../../stopping/compilation/README.md). Pure
Python 3 with no packages; needs ImageMagick (`magick`) to decode the GIFs,
Tesseract for the axis labels and the plot's totals line, and `pdftotext`
(poppler) for the citation lists. Issue #486.

The input is a local mirror of the site (the operator keeps one; it is not in
this repository and is never committed):

```sh
python3 -I run_all.py <mirror>/site <work-dir> ../../stopping/compilation/srim_org_h_he.csv
python3 -I bibliography.py <mirror>/site ../../stopping/compilation/srim_org_h_he.csv \
    ../../stopping/compilation/srim_org_bibliography.csv
python3 -I digitize.py <mirror>/site/StopPlots/STOP0114.gif "<mirror>/site/Citations/Cit 01-14.pdf"  # one plot, JSON
```

`<work-dir>` gets one JSON per plot (calibration, legend rows, every point
with its pixel, and the per-colour count check) for review; a plot already
there is not redone.

## What it does, per plot

1. **Frame and grid** (`digitize.py`). The top panel's box is the same on
   every plot; grid lines are the cyan rows and columns spanning it.
2. **Axes.** x is log, 1 to 10^5 keV/amu, with decades found from the grid
   pattern and checked against the OCR'd labels. y is linear from 0 on the
   bottom edge, so each OCR'd grid label votes for the step per grid line and
   the commonest step wins (a misread label is outvoted).
3. **Legend** (`legocr.py`). Rows sit on a fixed grid at the top right. Each
   row's colour is the commonest ink in its number columns. The four numbers
   (total, points below and at/above 1 MeV/amu, old citation number) are
   read with templates of the plot font's digits (`digits.json`, learned from
   the H→Si legend and then from every row that passed both checks below),
   by a best-path search that tolerates grid lines across the glyphs. Checks:
   total = below + above, and the sum over rows against the plot's own
   "Data Sets = N (M pts)" line. A legend that matches that line is
   *complete*.
4. **Citation keys** (`citpdf.py`). The old citation number is the number in
   the system's citation list (`Citations/Cit zz-ZZ.pdf`), which gives the key
   (e.g. `1981-Sant`) and the reference. Numbers repeat across the
   compilation, so a number found only in another system's list is used just
   when that entry names this plot's ion (`run_all.key_for`).
5. **Markers.** Hydrogen plots draw `*H` or `*D`: an asterisk (several raster
   variants, found by its 5-pixel bars and diagonals, tolerating cover by
   other ink) followed by a one-pixel label stroke. Helium plots draw `He`
   only (`symdetect.find_H`: two 8-pixel strokes 5 apart with a crossbar,
   and ink after it). Rejected: the black curve (`on_curve`), thick title
   text, the blue annotation block, the legend box.
6. **Anchor.** The plotted point is the centre of the label text: 4.72 px
   right of the asterisk centre (fitted on 1981-Sant, D→Si, 19 points at
   0.1 MeV steps, scatter 0.27 px) and 7.0 px right of the H's left stroke for
   `He` (1974-Bagl, 2 MeV, three plots).
7. **Experiment/Theory panel** (`residual.py`). Under each plot the same
   points are drawn divided by SRIM's curve, on a fixed scale (1.0 at row
   786, 0.05 per 39.5 px). Stacks on top spread apart here by their
   few-percent differences. The curve is traced in the top panel (a
   least-jump path through the black runs of each column), each lower marker
   becomes ratio × curve at its energy, and it is paired with a top marker of
   its colour (within 1.5 px in x and 1 % or 3 px in stopping, then 4 % for
   what is left on both sides, where markers on the curve blur its trace).
   The panel's blue title block is skipped, as blue is also a marker colour.
   A pair checks
   the top reading (`lower_check_rel` in the JSON); a lower marker with no
   partner is a point hidden on top and is added, flagged `panel=lower`,
   with the pairs' scatter on that plot as its uncertainty. The curve itself
   is never written out.
8. **Output.** Energy, stopping, uncertainties, the candidate keys for the
   point's colour and an attribution grade (see the data README).

## Limits

- Stacked markers cover each other in both panels, so dense regions still
  lose points; recall on the complete-legend plots is 84 % (61 % from the top panel
  alone). Where a complete legend's colour count is exceeded, the excess
  lower-only points nearest a top point are left out (`run_all.cap_lower`). The count check
  per colour is the guard against false points.
- Ions heavier than helium, compound targets (`SRIM/Compounds`) and the
  range plots are not handled yet.
