# Measured stopping powers digitized from the srim.org compilation

srim.org publishes, for each ion–target system, a plot of the stopping powers
measured by the papers in its bibliography. This directory holds those
measured points for **hydrogen, helium and lithium ions** (H, D, He and Li;
176 ion–target systems covering 74 target elements; 11,771 points, 1.1 keV/amu
to 65 MeV/amu), read back from the plots
by [`../../digitize/srim_org/`](../../digitize/srim_org/README.md), under the
rules in [`CONTRIBUTING.md`](../../../../CONTRIBUTING.md) ("The srim.org
experimental-data compilation") and origin 4 of
[`docs/data-provenance.md`](../../../../docs/data-provenance.md).

These are measurements, not SRIM output. SRIM's fitted curve is never stored.
It is read only as an intermediate: the Experiment/Theory panel under each
plot shows the same points divided by that curve, and a point hidden under
other markers on top is recovered there as ratio × curve and flagged
(`panel` = `lower`; operator ruling of 2026-10-28,
[recorded on #486](https://github.com/2AMLogic/lindhard/issues/486#issuecomment-6488202598)). No plot image, PDF
or page text from the site is stored here.

**License:** CC0 1.0 ([`LICENSE`](LICENSE)), separate from the MIT code. Cite
the original papers.

## Files

- `srim_org_stopping.csv`: one row per digitized point of the per-system
  plots (`StopPlots/`, plus the two data-bearing plots that exist only under
  `SRIM/SRIMPICS/STOP01/`: H in B and H in P).
- `srim_org_compounds.csv`: stopping in compounds, from the one-off figures
  of the site's Compounds page (H and He in water and polystyrene, He in five
  hydrocarbon gases, Li in three polymers); see
  [below](#srim_org_compoundscsv). The measured ranges of the same site are in
  [`../../ranges/compilation/`](../../ranges/compilation/README.md).
- `srim_org_h_he_loader.csv`: the H, D and He points of `origin = digitized-compilation`
  in the stopping loader's format (`docs/stopping-data-format.md`), made by
  [`convert_to_loader.py`](../../digitize/srim_org/convert_to_loader.py). Rows
  with several candidate papers keep every key, joined with `|`. Li points and
  the panel-recovered points (`origin = residual-reconstructed`) are not in it.
- `srim_org_bibliography.csv`: the papers behind the citation keys (key, old
  citation number, authors, year, title, journal reference), parsed from the
  compilation's per-system citation lists. Parsed automatically; check a
  reference against the paper before relying on its wording.

- `bibliography.csv`: the experimental stopping and range papers that srim.org's
  citation files record, one row per paper (1230 rows), with a compiler's
  comment, the ion–target pairs and compounds each is listed under, and a
  Crossref DOI where exactly one candidate matches. Described in
  [Experimental-paper bibliography](#experimental-paper-bibliography) below.
- `fetch_bibliography.py`, `test_fetch_bibliography.py`: the script that builds
  it, and its parser test.
- `bibliography-unparsed.txt`: entries the parser could not split, kept
  verbatim (two, the same paper listed in two compound files).

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
| `plot`, `colour`, `pixel_x`, `pixel_y` | Where the point was read: the plot's path on the site, the marker colour and the pixel, for review against the operator's archive. |
| `origin`, `method` | `digitized-compilation`, or `residual-reconstructed` for a `lower` point; and the digitizer version. |

### Attribution

Each plot colours its data sets from a palette of 18 colours, so a colour can
belong to more than one paper, and plots with many data sets list only some
of them in the legend ("Not All Data Shown").

- `unique` (4,558 points): the legend lists every data set on the plot (its
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
- **Two panels.** 6,450 points are seen in both panels; there the
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
  differences. On the 122 plots whose legends are complete and list points,
  5,313 of the 6,295 listed points were found (84 %: H and He 84 %, Li 85 %;
  64 % from the top panel alone). Where a colour has more points than its legend
  lists, the excess lower-only points nearest a top point (likely the same
  marker) are left out (26). Three top-panel points on those plots still
  exceed their colour's count, which bounds false detections at about
  0.1 %.
- **Unresolved papers.** 1,192 points carry a `ref<N>` candidate: the legend's
  number was not in that system's citation list, and the entry with that
  number elsewhere does not mention the ion. Old citation numbers are not
  unique across the compilation, so these are not guessed.

Every value here is a stand-in until the paper itself is read: a point read
from its paper replaces the digitized one (issue #483).

**Not digitized, by design:** ions heavier than lithium. srim.org has no
per-target plot for them, only the aggregate plots `STOPZZxx` (one ion into
all targets) and `STOPxxZZ` (all ions into one target). Those plot every
point scaled to a reference target or ion ("Normalized Tgt = Cu",
"Normalized Ion = Al") by SRIM's own ratio between the two, so the value
printed is not the measurement, and undoing the scaling would need SRIM's
tables.

Also not digitized, from the rest of the site:

- `SRIM/SRIMPICS/STOPHE*` (ratios of theory to experiment), `IPLOTS` (mean
  ionization potentials fitted through Bethe's formula), `STOPNUMB` (SRIM's
  own stopping components) and the History figures: curves and SRIM-derived
  values, no measured points as printed.
- `SRIM/SRIMPICS/STOP01` to `STOP03` repeat the `StopPlots/` files byte for
  byte and add 50 plots; 48 of those have no data points.
- Compounds page: "Compound Fig1" (lines fitted to hydrocarbon data, no
  points), "Compound Fig2" (a scanned figure from a paper whose few measured
  points have no stated energy), and "Core and Bond Stopping Values" (SRIM's
  own core-and-bond model values, not measurements).
- "He into Water" and "He into Ethelene Gas": 573 × 392 and 471 × 322 pixel
  images whose legends cannot be read (the keys are a few pixels tall) and
  whose 5-pixel letter markers sit in overlapping chains; a reading would be
  mostly guesswork.

## `srim_org_compounds.csv`

Digitized by `../../digitize/srim_org/figures.py` from the one-off figures of
the Compounds page, each described once in `figures.json` (axes, legend rows
as printed, marker examples). 263 points:

| Figure | Points | Legend | Keys |
|---|---|---|---|
| H in water (letters A to H) | 60 of 62 | A–F condensed water (the legend says "Water Solid"), G–H vapour | 8 printed; 6 exact, 1 corrected, 1 found in the water list but its entry does not name the system |
| H in polystyrene (`H*`) | 36 of 40 | four colours | 4 printed, all resolved |
| He in polystyrene (`He`) | 26 of 26 | three colours | 3 printed, all resolved |
| Li in polyimide, formvar, polysulphone (`Li`) | 51 | colour = target; no counts | none: the figure cites F. Munnik et al., J. Appl. Phys. 85, 3934 (1999) |
| He in five hydrocarbon gases (`He`) | 90 (18 per gas, 300 to 2000 keV in 100 keV steps) | colour = gas; no counts | none: the site says "All of the data is from the papers of D. Powers (Baylor Univ., USA)"; `candidate_key` 1978-Chau for the three gases its entry names |

Units differ by figure and are kept as printed: energy in keV/amu, MeV or
keV; stopping per molecule or per atom as the axis says (`energy_unit`,
`value_unit`).

Columns: `figure`, `quantity`, `ion_z`, `ion`, `target`, `phase`, `energy`,
`energy_unit`, `value`, `value_unit`; `u_energy_rel` (relative) and
`u_value` (absolute), read-off only: half a pixel and the anchor's
uncertainty (0.7 px; 1.0 px for the hand-drawn hydrocarbon template; 2.4 px
for Li, whose energies read 1.05 % below round MeV values that are not
assumed), plus 1.5 px for a point placed by eye; `symbol`, `colour`;
`key_printed`, `citation_key`, `key_status`, `reference` (one entry per
candidate legend row, `|`-joined in the same order); `candidate_key`; `how`
(`auto`, or `seeded`: read by eye and placed by template match within 4 px,
or `eye`: placed by eye where no template fits), `note`; `pixel_x`,
`pixel_y`; `origin`, `method`.

`key_status`: `exact`; `corrected` (one confusable digit or author letter,
to the only key in the list that names the system); `unconfirmed` (in the
figure's own per-compound citation list, but the entry does not name the
system); `key-mismatch` (the key exists, for another system); `not-found`;
`not-printed`.

Where markers overlap (the H-in-water cluster near 40 keV/amu, the
high-energy tail of H in polystyrene), points were checked by eye at 10×.
Four H-in-water points keep an uncertain letter (`symbol` `C|D` or
`C|D|F`), so their paper is one of those rows; one more is read as D from
its lower half (noted).

## Experimental-paper bibliography

`bibliography.csv` covers every ion–target pair (202 citation files, from
`PlotCit.htm`) and every compound (148 files, from `SRIM/Compounds.htm`), not
only the H and He plots, fetched 2026-10-29. Titles, authors, journal
references and years are facts about published papers; the compiler's comment
says what each paper measured. The citation PDFs and page text are not copied
into the tree. The license is the CC0 1.0 one above.

`fetch_bibliography.py` fetches the files (at most one request per second, into
a scratch directory outside the tree), runs `pdftotext -layout` and parses the
entries. The `Cit <Z>-All` and `Cit All-<Z>` aggregates repeat per-pair
entries and serve only as a cross-check.

| column | meaning |
|---|---|
| `citation_key` | the compiler's key as printed (for example `1962-Ewin`); not unique, since 73 keys are used for more than one paper |
| `year`, `authors`, `title`, `journal` | as printed; title text order comes from the PDF layout and can be scrambled in a few entries |
| `comment` | the compiler's note: quantity (S stopping, R range, d prefix for straggling), energy range, ions and targets |
| `old_ref` | the old reference number; empty in the compound files, which do not print one |
| `pairs` | ion–target pairs the paper is listed under, as `<ion Z>-<target Z>` (`01-14` is H in Si), joined with `;` |
| `target` | compound names the paper is listed under, joined with `;`; empty for papers listed only under elemental pairs |
| `doi` | Crossref DOI, only when exactly one candidate matches the normalized title, the year, the journal and, where printed, the first page (275 of 1230 filled); never guessed |

Checks: per file, the `Comment :` lines, the entry-start lines and the parsed
plus unparsed entries must agree (all 350 files pass). No measured value is
in `bibliography.csv`. Rebuild with `python3 fetch_bibliography.py` (about 40
minutes, needs poppler's `pdftotext`; `--no-fetch` reuses the scratch
directory). Test the parser with
`python3 -B -m unittest test_fetch_bibliography`.

### Relation to `srim_org_bibliography.csv`

The two files are **complementary, not a superset pair**, and neither
supersedes the other. `srim_org_bibliography.csv` (267 rows) is keyed to the
citation lists behind the digitized H/He plots and carries `old_ref` as the
join to `old_refs` in `srim_org_stopping.csv`. `bibliography.csv` is built from all
350 citation files and has a different parser. Compared on `citation_key`:

- 262 of the 267 keys also occur in `bibliography.csv`; year agrees on all 262.
  Normalized title agrees on 253 (the rest differ in apostrophes, truncation
  or layout order); normalized journal reference on 203 (the older file keeps
  trailing page-footer text on some rows). 25 of the 262 keys name more than
  one paper in `bibliography.csv`, so match on key plus title.
- Five keys of the older file are **absent** from `bibliography.csv`:
  `1961-Niel`, `1966-Palm`, `1968-Powe`, `1980-Sant`, `2000-Berg`. Their
  titles do not appear there under other keys either. Keep
  `srim_org_bibliography.csv` for resolving the digitized points' citation
  keys.

The comparison was made once by script on the committed CSVs; neither file's
rows were checked by hand against the papers.
