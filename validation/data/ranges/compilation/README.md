# Measured ranges digitized from the srim.org figures

srim.org shows four figures of measured ranges against SRIM's projected-range
curve (`SRIM/SRIMPICS/RANGES/`): B, P and As in Si, and He in air (STP). This
directory holds the measured points of those figures, read back by
[`../../digitize/srim_org/figures.py`](../../digitize/srim_org/README.md)
under the rules in [`CONTRIBUTING.md`](../../../../CONTRIBUTING.md) ("The
srim.org experimental-data compilation") and origin 4 of
[`docs/data-provenance.md`](../../../../docs/data-provenance.md).

These are measurements, not SRIM output. SRIM's curve is not read at all, and
no image, PDF or page text from the site is stored here.

**License:** CC0 1.0 ([`LICENSE`](LICENSE)), separate from the MIT code. Cite
the original papers.

## `srim_org_ranges.csv`

221 points, one row each. The columns are those of
[`srim_org_compounds.csv`](../../stopping/compilation/README.md#srim_org_compoundscsv);
here `quantity` is `range`, and `value` is the figure's "Target depth" in nm
(Si) or cm (air). The papers measured ranges of different kinds (projected,
mean, extrapolated), and the figures do not say which; read the paper before
comparing a point with a projected range.

| Figure | Points found / listed | Energy | Range |
|---|---|---|---|
| B in Si | 84 / 200 | 18 keV to 2.5 MeV | 60 nm to 3.1 µm |
| P in Si | 49 / 58 | 30 keV to 2.7 MeV | 40 nm to 2.5 µm |
| As in Si | 18 / 37 | 19 keV to 1.2 MeV | 21 to 735 nm |
| He in air (STP) | 70 / 121 | 29 keV to 8.7 MeV | 0.04 to 8.7 cm |

Markers stack where many papers measured the same energies, and the hidden
ones cannot be recovered (these figures have no Experiment/Theory panel).
On an overlay of B in Si the markers left out are nearly all under others.

### Attribution, and the scrambled legends

The legend rows are coloured from the plotting program's 20-colour cycle, so
a colour belongs to every row 20 places apart (and two entries of the cycle
are the same red, and two the same magenta). A point lists every row of its
colour, `|`-joined, with each row's printed key, resolved key, status and
reference in the same order. Only on P in Si, As in Si and He in air does
a colour often belong to one row (78 points); on B in Si (57 rows) none does.

The printed keys of the Si figures are often not papers about that system:
1913-Tayl is He in air and 1958-Kepl is about muons. The legend font is
small, and its 9 and 8 (and 1 and 7, 3 and 5) look alike, so `key_printed`
is the transcription as read (`1852-Simm` may print 1952). Each key is checked
against the ion's citation list (`Citations/Cit zz-All.pdf`) and accepted
only if its entry names the ion and Si; `corrected` takes the one key that
does, one confusable digit or author letter away (8 for 6, 7 for 1). The
rest keep their printed key with status `key-mismatch` (the key exists, for
another system) or `not-found`. On He in air all ten keys resolve against
the site's air list (`SRIM/Compounds/Citations/Air-He.pdf`).

| Figure | Legend rows | exact | corrected | key-mismatch | not-found |
|---|---|---|---|---|---|
| B in Si | 57 | 31 | 9 | 11 | 6 |
| P in Si | 24 | 4 | 3 | 12 | 5 |
| As in Si | 16 | 2 | 1 | 12 | 1 |
| He in air | 10 | 8 | 2 | 0 | 0 |

## How good it is

- **Axes.** Every grid line lands within 1.2 px of a round value.
- **Anchor.** The point is the centre of the marker's label. On B and P in
  Si, markers at round energies sit 0.5 and 0.7 px from them, so 0.7 px is
  added to the energy uncertainty (`u_energy_rel`, about 0.6 %).
- **Uncertainty.** `u_energy_rel` and `u_value` are read-off only (half a
  pixel and the anchor). They say nothing about the paper's own error, or
  about which range it measured.

## Count check and hand check

`count_check.csv` is the per-colour count check of each figure: the legend
rows printed in that colour, their printed total, the points found, and a
status. `match` and `mismatch` are per colour; `ambiguous` means several
legend rows share the colour, so only the group total can be checked (every
colour on B in Si). The last row of each figure is its total against the
legend total. Nothing is forced to match: a `mismatch` is a count the
digitizer could not reach, almost always markers hidden under others.

`HANDCHECK.md` records ten points of B in Si read by eye at 8x on
2026-10-30 against the recorded values: all within 1.7 % in depth and 0.7 %
in energy (about 1.5 px at worst), with the axis references checked against
the tick labels. On a rerun of `figures.py` the same day, from the same
mirror, `srim_org_ranges.csv` came back byte-identical.

The issue that asked for this data (#495) wanted one file per (ion, target,
citation key). That is not possible honestly on these figures: a colour
belongs to every legend row 20 places apart, so most points have several
candidate papers (the `|`-joined columns above). One CSV with the candidate
list per point is what the figures support; splitting by paper is the work of
#483, after the papers are read.
