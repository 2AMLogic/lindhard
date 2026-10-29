# Measured stopping powers digitized from the srim.org compilation

srim.org publishes, for each ion–target system, a plot of the stopping powers
measured by the papers in its bibliography. This directory holds those
measured points for **hydrogen and helium ions** (H, D and He into 132
targets, 7,241 points, 1.1 keV/amu to 65 MeV/amu), read back from the plots
by [`../../digitize/srim_org/`](../../digitize/srim_org/README.md), under the
rules in [`CONTRIBUTING.md`](../../../../CONTRIBUTING.md) ("The srim.org
experimental-data compilation") and origin 4 of
[`docs/data-provenance.md`](../../../../docs/data-provenance.md).

These are measurements, not SRIM output. The fitted SRIM curve on each plot and
the Experiment/Theory panel under it were never read. No plot image, PDF or
page text from the site is stored here.

**License:** CC0 1.0 ([`LICENSE`](LICENSE)), separate from the MIT code. Cite
the original papers.

## Files

- `srim_org_h_he.csv`: one row per digitized point.
- `srim_org_bibliography.csv`: the papers behind the citation keys (key, old
  citation number, authors, year, title, journal reference), parsed from the
  compilation's per-system citation lists. Parsed automatically; check a
  reference against the paper before relying on its wording.

These files sit in a subdirectory, so the stopping report
(`lindhard-plots/src/stopping_report.rs`), which reads only the CSVs at the
top of `stopping/`, does not load them. Its source guard refuses any
dataset that mentions srim. Wiring this compilation into the report needs
that guard to accept origin 4 explicitly; that is a separate change.

### `srim_org_h_he.csv` columns

| Column | Meaning |
|---|---|
| `ion_z`, `ion` | Ion. `ion` is `D` where the marker's label reads D (deuteron); `H` otherwise for Z = 1. |
| `target_z`, `target` | Target element. |
| `energy_kev_per_amu` | Energy per mass unit, as plotted (log axis). |
| `stopping_ev_cm2_per_1e15_atoms` | Stopping cross section as plotted, eV/(10^15 atoms/cm²). The compilation plots total measured stopping; below about 10 keV/amu that includes a nuclear part. |
| `u_energy_rel` | Relative digitization uncertainty of the energy (anchor scatter 0.27 px plus half a pixel, 0.8 %). |
| `u_stopping` | Digitization uncertainty of the stopping, half a pixel, in the same unit. |
| `citation_keys` | Candidate papers for the point: every legend row printed in the point's colour. `ref<N>` where the old citation number `N` could not be matched to a paper (see below). |
| `old_refs` | The legend's old citation numbers for those rows, in the same order. |
| `attribution` | How sure the paper is (below). |
| `label` | The marker label read: `H`, `D`, `He`, or `?` if covered. |
| `plot`, `colour`, `pixel_x`, `pixel_y` | Where the point was read: the plot file name on the site, the marker colour and the pixel, for review against the operator's archive. |
| `origin`, `method` | Always `digitized-compilation` and the digitizer version. |

### Attribution

Each plot colours its data sets from a palette of 18 colours, so a colour can
belong to more than one paper, and plots with many data sets list only some
of them in the legend ("Not All Data Shown").

- `unique` (3,433 points): the legend lists every data set on the plot (its
  row counts sum to the plot's own "Data Sets = N (M pts)" line), and one row
  has this colour. The paper is known.
- `shared-colour` (636): the legend is complete, but two or more rows share
  the colour. `citation_keys` lists them; the point belongs to one of them.
  Between H and D rows the label usually decides it.
- `partial-legend` (3,172): the plot lists only some data sets, so the point
  may also belong to one that is not listed. Use these for broad checks only.

## How good it is

- **Values.** Against Bader et al., Phys. Rev. 103, 32 (1956), Table I
  (`../bader1956_p_cu.csv`), the ten 1956-Bade points on the H→Cu plot agree
  within 0.2 % in stopping at the same energies (50 to 550 keV).
- **Energy anchor.** A marker's point is the centre of its label text, not
  the asterisk. That offset was fitted on 1981-Sant (D→Si, 0.2 to 2.0 MeV in
  0.1 MeV steps, 19 points; 4.72 ± 0.27 px) and, for He labels, on 1974-Bagl
  (2 MeV He, the same position on the Rh, Hf and W plots).
- **Completeness.** Markers stack where data are dense, and many are covered
  by others. On the 99 plots with complete legends, 4,069 of the 6,680 listed
  points were found (61 %); isolated points are found, the bottom of a dense
  stack is not. Six points on those plots exceed their colour's legend count,
  which bounds false detections at about 0.1 %.
- **Unresolved papers.** 726 points carry a `ref<N>` candidate: the legend's
  number was not in that system's citation list, and the entry with that
  number elsewhere does not mention the ion. Old citation numbers are not
  unique across the compilation, so these are not guessed.

Every value here is a stand-in until the paper itself is read: a point read
from its paper replaces the digitized one (issue #483). Ions heavier than
helium, the compound plots and the range plots are not yet digitized (#486).
