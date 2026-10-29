# srim.org stopping-plot digitizer

Reads the measured points off srim.org's per-system stopping plots
(`StopPlots/STOPzzZZ.gif`, ion Z then target Z, and the two with data that
exist only under `SRIM/SRIMPICS/STOP01/`) and writes
[`../../stopping/compilation/`](../../stopping/compilation/README.md).
`figures.py` does the same for the site's one-off figures, the measured
ranges and the stopping in compounds ([below](#the-one-off-figures-figurespy)). Pure
Python 3 with no packages; needs ImageMagick (`magick`) to decode the GIFs,
Tesseract for the axis labels and the plot's totals line, and `pdftotext`
(poppler) for the citation lists. Issue #486.

The input is a local mirror of the site (the operator keeps one; it is not in
this repository and is never committed):

```sh
python3 -I run_all.py <mirror>/site <work-dir> ../../stopping/compilation/srim_org_stopping.csv
python3 -I bibliography.py <mirror>/site ../../stopping/compilation/srim_org_stopping.csv \
    ../../stopping/compilation/srim_org_bibliography.csv
python3 -I digitize.py <mirror>/site/StopPlots/STOP0114.gif "<mirror>/site/Citations/Cit 01-14.pdf"  # one plot, JSON
python3 -I convert_to_loader.py ../../stopping/compilation/srim_org_stopping.csv \
    ../../stopping/compilation/srim_org_bibliography.csv ../../stopping/compilation/srim_org_h_he_loader.csv
python3 -I figures.py <mirror>/site <work-dir> --ranges=../../ranges/compilation/srim_org_ranges.csv \
    --compounds=../../stopping/compilation/srim_org_compounds.csv
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
   and ink after it). Lithium plots draw `Li` (`symdetect.find_L`: the L's
   corner, an 8-pixel stroke meeting a 6-pixel bar with no crossbar, and the
   i's stem 10 to 12 px to the right, as the scaled font varies by a pixel).
   Rejected: the black curve (`on_curve`), thick title text, the blue title
   block (its "L" in "Lithium" looks like a marker), the blue annotation
   block, the legend box.
6. **Anchor.** The plotted point is the centre of the label text: 4.72 px
   right of the asterisk centre (fitted on 1981-Sant, D→Si, 19 points at
   0.1 MeV steps, scatter 0.27 px) and 7.0 px right of the H's left stroke for
   `He` (1974-Bagl, 2 MeV, three plots) and `Li`. For `Li` that offset was
   checked on 1975-Neuw2 (11 energies, 80 to 840 keV, three plots): it puts
   them on round lab energies to 0.14 ± 0.15 px.
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
  lose points; recall on the complete-legend plots is 84 % (H and He 84 %,
  Li 85 %; 64 % from the top panel alone). Where a complete legend's colour count is exceeded, the excess
  lower-only points nearest a top point are left out (`run_all.cap_lower`). The count check
  per colour is the guard against false points.
- Ions heavier than lithium are out of scope: the site plots them only in
  the aggregates `STOPZZxx` and `STOPxxZZ`, scaled to a reference target or
  ion by SRIM's own ratio, so the printed value is not the measurement.

## The one-off figures (`figures.py`)

The four range figures (`SRIM/SRIMPICS/RANGES/`) and the compound figures
(`SRIM/Compounds/`) each have their own layout, so each is described once in
`figures.json`: two axis labels per axis, read off the figure; the legend
rows as printed (transcribed by eye: the legend font is too small for
Tesseract); and one or more clean examples of each marker. The rest is
measured and checked:

- **Axes.** The labels are snapped to grid lines, and every other grid line
  must then land on a round value (worst 1.2 px over all nine figures).
- **Legend.** On the range figures the rows sit on a fixed pitch and take
  their colours from the plotting program's 20-colour cycle; the colours
  found must agree with the cycle on 90 % of rows. The compound figures give
  their rows explicitly: a letter per data set (H in water), a colour per
  data set, or a colour per target with one reference printed under the
  plot.
- **Keys.** A printed key is accepted only if a citation list on the site
  has it and its entry names the ion and the target; one confusable digit or
  author letter may be corrected to the only key that does.
- **Markers.** Templates are the examples' own ink plus every isolated mark
  of about their size. Letters share a size, so an isolated letter becomes a
  template of the example letter it overlaps best. Where every marker
  touches a curve (He in hydrocarbon gases), the H of `He` is drawn by hand
  in its twelve raster variants. A match is ranked by the share of the
  template it covers, less the ink of its colour around it that it does not
  explain, so an E is not also read as an F; the best match wins where
  several fit one place. The point is the centre of the marker's label.
- **Seeds and rejects.** In the H-in-water cluster and the H-in-polystyrene
  tail, markers were read by eye at 10× and listed in `figures.json`: a seed
  is placed by the best template match within 4 px of the reading, or by eye
  where nothing fits; a reject removes a match read as part of a neighbour.
  Every such point says how it was placed (`how`, `note`).

The checks that the figures pass: H in water 60 of 62 legend points, He in
polystyrene 26 of 26 with energies on round values (99.8, 199.6, 249.8
keV/amu), He in hydrocarbon gases 18 points per gas on 300 to 2000 keV in
100 keV steps (299 to 1997 keV read). The counts and key statuses per
figure are in the two data READMEs.

Not handled: the low-resolution "He into Water" and "He into Ethelene Gas"
(illegible legends, 5-pixel markers in overlapping chains), and the figures
that hold no measured points as printed (listed in the stopping data
README).
