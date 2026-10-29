# Measured stopping powers digitized from the srim.org compilation

srim.org publishes, for each ion–target system, a plot of the stopping powers
measured by the papers in its bibliography. This directory holds those
measured points for **hydrogen, helium and lithium ions** (H, D, He and Li;
174 ion–target systems covering 73 target elements; 11,706 points, 1.1 keV/amu
to 65 MeV/amu), read back from the plots
by [`../../digitize/srim_org/`](../../digitize/srim_org/README.md), under the
rules in [`CONTRIBUTING.md`](../../../../CONTRIBUTING.md) ("The srim.org
experimental-data compilation") and origin 4 of
[`docs/data-provenance.md`](../../../../docs/data-provenance.md).

These are measurements, not SRIM output. SRIM's fitted curve is never stored.
It is read only as an intermediate: the Experiment/Theory panel under each
plot shows the same points divided by that curve, and a point hidden under
other markers on top is recovered there as ratio × curve and flagged
(`panel` = `lower`, operator ruling of 2026-10-28, #486). No plot image, PDF
or page text from the site is stored here.

**License:** CC0 1.0 ([`LICENSE`](LICENSE)), separate from the MIT code. Cite
the original papers.

## Files

- `srim_org_stopping.csv`: one row per digitized point.
- `srim_org_h_he_loader.csv`: the H, D and He points of `origin = digitized-compilation`
  in the stopping loader's format (`docs/stopping-data-format.md`), made by
  [`convert_to_loader.py`](../../digitize/srim_org/convert_to_loader.py). Rows
  with several candidate papers keep every key, joined with `|`. Li points and
  the panel-recovered points (`origin = residual-reconstructed`) are not in it.
- `srim_org_bibliography.csv`: the papers behind the citation keys (key, old
  citation number, authors, year, title, journal reference), parsed from the
  compilation's per-system citation lists. Parsed automatically; check a
  reference against the paper before relying on its wording.

These files sit in a subdirectory, so the stopping report
(`lindhard-plots/src/stopping_report.rs`), which reads only the CSVs at the
top of `stopping/`, does not load them. Its source guard refuses any dataset
that mentions srim unless every point has `origin = digitized-compilation`,
which `srim_org_h_he_loader.csv` has. Moving that file to the top of
`stopping/` would add it to the report and change the committed page; that is
a separate change.

### `srim_org_stopping.csv` columns

| Column | Meaning |
|---|---|
| `ion_z`, `ion` | Ion. `ion` is `D` where the marker's label reads D (deuteron); `H` otherwise for Z = 1. |
| `target_z`, `target` | Target element. |
| `energy_kev_per_amu` | Energy per mass unit, as plotted (log axis). |
| `stopping_ev_cm2_per_1e15_atoms` | Stopping cross section as plotted, eV/(10^15 atoms/cm²). The compilation plots total measured stopping; below about 10 keV/amu that includes a nuclear part. |
| `u_energy_rel` | Relative digitization uncertainty of the energy (anchor scatter 0.27 px plus half a pixel, 0.8 %). |
| `u_stopping` | Digitization uncertainty of the stopping, in the same unit: half a pixel for a point read on top; for a `lower` point, the scatter between the two panels' readings of the plot's paired points (0.7 to 2.1 px, median 0.4 % of the value). |
| `citation_keys` | Candidate papers for the point: every legend row printed in the point's colour. `ref<N>` where the old citation number `N` could not be matched to a paper (see below). |
| `old_refs` | The legend's old citation numbers for those rows, in the same order. |
| `attribution` | How sure the paper is (below). |
| `label` | The marker label read: `H`, `D`, `He`, `Li`, or `?` if covered. |
| `panel` | Where the point was seen: `top` (top panel only), `both` (top, and confirmed by its Experiment/Theory marker), or `lower` (Experiment/Theory panel only; value is ratio × curve). |
| `plot`, `colour`, `pixel_x`, `pixel_y` | Where the point was read: the plot file name on the site, the marker colour and the pixel, for review against the operator's archive. |
| `origin`, `method` | `digitized-compilation`, or `residual-reconstructed` for a `lower` point; and the digitizer version. |

### Attribution

Each plot colours its data sets from a palette of 18 colours, so a colour can
belong to more than one paper, and plots with many data sets list only some
of them in the legend ("Not All Data Shown").

- `unique` (4,493 points): the legend lists every data set on the plot (its
  row counts sum to the plot's own "Data Sets = N (M pts)" line), and one row
  has this colour. The paper is known.
- `shared-colour` (758): the legend is complete, but two or more rows share
  the colour. `citation_keys` lists them; the point belongs to one of them.
  Between H and D rows the label usually decides it.
- `partial-legend` (6,455): the plot lists only some data sets, so the point
  may also belong to one that is not listed. Use these for broad checks only.

## How good it is

- **Values.** Against Bader et al., Phys. Rev. 103, 32 (1956), Table I
  (`../bader1956_p_cu.csv`), the ten 1956-Bade points on the H→Cu plot agree
  within 0.2 % in stopping at the same energies (50 to 600 keV), and the
  three points recovered from the Experiment/Theory panel (74, 248 and
  346 keV) within 0.3 %.
- **Two panels.** 6,405 points are seen in both panels; there the
  reconstructed value agrees with the top reading to +0.16 % on average,
  scatter 0.6 %.
- **Energy anchor.** A marker's point is the centre of its label text, not
  the asterisk. That offset was fitted on 1981-Sant (D→Si, 0.2 to 2.0 MeV in
  0.1 MeV steps, 19 points; 4.72 ± 0.27 px) and, for He labels, on 1974-Bagl
  (2 MeV He, the same position on the Rh, Hf and W plots). `Li` labels are
  also two characters from the left of the first letter and take the He
  offset; on 1975-Neuw2 (Li, 11 energies 80 to 840 keV, the same pixels on
  the Be, Cu and Ta plots) that offset puts the points on round lab
  energies to within 0.14 ± 0.15 px.
- **Completeness.** Markers stack where data are dense, and many are covered
  by others; the Experiment/Theory panel spreads a stack by its few-percent
  differences. On the 171 plots with complete legends, 5,248 of the 6,229
  listed points were found (84 %, alike for H, He and Li; 61 % for H and He
  from the top panel alone). Where a colour has more points than its legend
  lists, the excess lower-only points nearest a top point (likely the same
  marker) are left out (26). Three top-panel points on those plots still
  exceed their colour's count, which bounds false detections at about
  0.1 %.
- **Unresolved papers.** 1,192 points carry a `ref<N>` candidate: the legend's
  number was not in that system's citation list, and the entry with that
  number elsewhere does not mention the ion. Old citation numbers are not
  unique across the compilation, so these are not guessed.

Every value here is a stand-in until the paper itself is read: a point read
from its paper replaces the digitized one (issue #483). The compound plots
and the range plots are not yet digitized (#486).

**Not digitized, by design:** ions heavier than lithium. srim.org has no
per-target plot for them, only the aggregate plots `STOPZZxx` (one ion into
all targets) and `STOPxxZZ` (all ions into one target). Those plot every
point scaled to a reference target or ion ("Normalized Tgt = Cu",
"Normalized Ion = Al") by SRIM's own ratio between the two, so the value
printed is not the measurement, and undoing the scaling would need SRIM's
tables.
