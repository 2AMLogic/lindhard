#!/usr/bin/env python3
"""Digitize the measured Ar -> Cu, Si, Ag and Au sputtering yields of the
compilation N. Matsunami, Y. Yamamura, Y. Itikawa, N. Itoh, Y. Kazumata,
S. Miyagawa, K. Morita, R. Shimizu and H. Tawara, "Energy dependence of the
yields of ion-induced sputtering of monatomic solids", report IPPJ-AM-32,
Institute of Plasma Physics, Nagoya University (1983); published as At. Data
Nucl. Data Tables 31, 1 (1984), doi:10.1016/0092-640X(84)90016-0.

This is the extraction behind `validation/data/sputtering/ar_<target>_sputter_*.json`
(provenance: docs/data-provenance.md, operator decision #69; Cu in #69, Si, Ag
and Au in #70). It is not part of the harness and needs numpy, scipy and
Pillow, which the harness does not. Neither the report nor any page image is
in this tree: fetch the report yourself and extract the page raster of the
figure with poppler (PDF p. 118 Ar -> Cu, p. 47 Ar -> Si, p. 191 Ar -> Ag,
p. 262 Ar -> Au):

    curl -O http://dpc.nifs.ac.jp/IPPJ-AM/IPPJ-AM-32.pdf
    pdfimages -f 118 -l 118 -png IPPJ-AM-32.pdf am32
    validation/data/digitize/matsunami1983_ipp_am32.py am32-000.png
    validation/data/digitize/matsunami1983_ipp_am32.py am32-000.png --write validation/data/sputtering
    pdfimages -f 191 -l 191 -png IPPJ-AM-32.pdf am32ag
    validation/data/digitize/matsunami1983_ipp_am32.py am32ag-000.png --target Ag --write validation/data/sputtering \
        --record validation/data/digitize/secondread_ar_ag_am32.json

Each page is a lossless 1-bit CCITT scan, 3307 x 4677 px at 400 ppi, so every
reader gets the same pixels and the pixel coordinates below are
reproducible. Every figure has the same log-log frame (50 eV to 10^6 eV,
10^-3 to 50 atoms/ion) and plots each reference's measured points with its
own letter (legend in the figure) and the compilation's empirical fit as a
solid line. Only the letters are read; the fitted line is never digitized or
stored (it is masked out, step 3). Everything that depends on the figure
(frame threshold, legend box, template scale, curve anchors, symbols,
dropped symbols, references) is in that figure's table below (CU_*, SI_*,
AG_*, AU_*, collected in FIGURES); the method is shared.

Method (also summarized in each dataset's `extraction`):

1. Axes. The frame is slightly rotated in the scan (about 0.1-0.5 deg), so
   the pixel -> (log10 E, log10 Y) map is a 2-D affine fit, not two 1-D fits.
   The frame sides are located from rows and columns that are mostly ink
   (the column threshold is per figure: 40 % of the page height on p. 118,
   20 % on the others, whose vertical sides are fainter or more tilted).
   Every tick on all four sides (minor ticks at 2..9 x 10^n, long ticks at
   10^n and 5 x 10^n) is found as a run of ink perpendicular to the fitted
   frame line. A run reaching the 80 px scan limit is a frame line, not a
   tick. Each long tick takes the nearest 10^n or 5 x 10^n on the straight
   line through the frame corners' values; a line fitted to the long ticks
   then gives every tick its value (nearest 1..9 x 10^n), long ticks that
   miss their value and ticks off the comb by more than 0.02 decade being
   dropped and printed. The six affine coefficients are least-squares fits
   to all kept ticks. The script stops if a side keeps fewer ticks than the
   figure's minimum (35 on p. 118, 30 on p. 47, whose bottom axis is faint,
   35 on the others), and prints the residuals (decades).
2. Templates. Each legend letter (the column of letters beside the reference
   names) is cut out of the page. The plotted letters are drawn smaller than
   the legend letters, so the templates are resampled by the figure's
   `scale`. That factor is a constant in the figure's table, measured once
   as the ratio of the ink bounding boxes of isolated plotted glyphs to their
   legend glyphs (Cu: on T, V, A, I; Si, Ag, Au: the median over 7, 8 and 11
   isolated glyphs, among them the figure's `isolated_check` glyphs). The
   script does not recompute it; it prints the template coverage of the
   `isolated_check` glyphs at both template sizes.
3. Fitted curve. The compilation's empirical curve crosses many symbols. It
   is located from anchor points on the curve (hand-read for Cu; for Si, Ag
   and Au taken about every 60 px where the curve's ink is a single thin run,
   with one hand-read point in the crowded 200 eV region of Ag), snapped to
   the thin ink run nearest a spline through them, and treated as "don't
   care" (neither ink nor background) when templates are matched. It is not
   stored.
4. Symbols. Many letters overlap. Each symbol was identified by eye from 4x
   to 8x zooms and character dumps of the raster, helped by the template
   coverage of every legend letter at the spot (for Cu also by a
   template-matching cluster solver); the identifications are the figure's seeds
   (letter, approximate centre in pixels, status, note, glyph box). Identity
   rules: a glyph is kept only where its distinguishing strokes are visible
   (e.g. the H crossbar sits at 31-42 % of the glyph height and the A
   crossbar at 51-64 %; R and B differ only by B's bottom bar; E and F by
   E's bottom bar; the slashed O; the G inner bar). "infer" marks a glyph
   whose identity rests on an argument stated in the seed note; the figure's
   dropped list holds every symbol in range that was not stored, and why.
   Nothing is guessed.
5. Centres. The plotted symbols are centred on the data point: in x the glyph
   centres of the evidently round-energy points (150, 200, 250, 300, 350 eV
   and so on) fall within about 1 % of the round values, whereas lower-left
   anchoring would put them 6-10 % low. The same centring is assumed in y.
   Each centre is refined twice, by the maximum of template coverage with
   the scaled template (1 px tolerance) and with the unscaled legend glyph
   (2 px tolerance); the value is their mean and half their difference is
   the centring repeatability.
6. Values. E and Y come from the affine map. Only points with
   196 eV <= E <= 10.2 keV are stored (about 0.2 to 10 keV, with 2 % slack
   for nominal 200 eV and 10 keV points). Uncertainties are digitizing
   uncertainties only (the compilation gives no measurement errors): in each
   log coordinate, the quadrature sum of the calibration residual rms, the
   centring repeatability and the symbol-placement scatter, the latter
   measured as the rms deviation of the round-energy points from their
   nominal energies; floored at 2 %.
7. Second read (Si, Ag and Au; `--record`). Every stored point is located a
   second time without templates: the sixth seed field is a box around the
   glyph read by eye from the zooms and dumps, and the centre is the centre
   of the bounding box of the ink inside it (`second_read()`; the curve's and
   neighbours' ink is not removed, so where they cross the box the centre is
   that of the box as read). The record lists the
   manual / template ratios in E and Y per point against the stored
   uncertainty u and against sqrt(2) u; each dataset's `crosscheck` gets a
   `second_read` summary. The Cu seeds have no boxes (Cu was cross-checked
   against NIFS-DATA-23 in #69).

Reference letters, years and citations are those of the figure legend and of
the report's "References for graphs" (PDF pp. 281-286). The compilation
states (pp. 1-2) that it holds normal-incidence, room-temperature data, that
single-crystal targets are excluded, that most specimens were thinned
polycrystalline foils or evaporated films, and that energies are the total
ion energy in eV and yields are atoms per ion.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.signal import fftconvolve

E_MIN_EV, E_MAX_EV = 196.0, 10200.0
FRAME_E, FRAME_Y = (50.0, 1e6), (1e-3, 50.0)  # values at the frame sides, the same in every figure

COMPILATION = "Matsunami et al., IPPJ-AM-32 (1983); At. Data Nucl. Data Tables 31, 1 (1984)"
URL = "http://dpc.nifs.ac.jp/IPPJ-AM/IPPJ-AM-32.pdf"
TARGET_STATE = (
    "not stated in the compilation for this set; the compilation excludes single-crystal targets and states "
    "that most specimens were thinned polycrystalline foils or evaporated films (IPPJ-AM-32, p. 2)"
)
RELIABILITY = (
    "Included in the compilation's Ar -> {target} figure, so it passed the compilation's selection (IPPJ-AM-32, pp. 1-2: "
    "absolute or normalized-relative yields, room temperature, no dose effects, sufficient beam current; doubtful "
    "sets excluded). The compilation does not rank the sets it shows."
)

# ======================================================================
# Ar -> Cu, PDF p. 118 (printed p. 113); issue #69
# ======================================================================
# Step 3: hand-read points on the fitted curve (x, y in px), snapped to ink.
CU_CURVE_ANCHORS = [
    (700, 2026), (780, 1924), (840, 1858), (900, 1798), (960, 1745), (1020, 1697), (1080, 1653),
    (1140, 1612), (1200, 1581), (1320, 1515), (1380, 1486), (1440, 1462), (1560, 1403), (1600, 1396),
    (1660, 1385), (1720, 1370), (1780, 1355), (1830, 1347), (1860, 1344), (1920, 1339), (1980, 1335),
    (2100, 1333),
]

# Step 4: every symbol stored. (letter, x px, y px, status, note)
CU_SEEDS = [
    ("T", 880, 1776, "keep", "overlaps the H below"),
    ("H", 881, 1780, "infer", "H and R are indistinguishable under the overlapping T (R = H plus a top bar); "
     "R (Oechsner 1973) appears only near 1 keV, and the H series continues at 250-500 eV"),
    ("H", 934, 1728, "keep", ""),
    ("T", 935, 1755, "keep", "overlaps the H above"),
    ("T", 975, 1804, "keep", "isolated"),
    ("F", 996, 1671, "keep", "overlaps H(1012,1663)"),
    ("H", 1012, 1663, "keep", ""),
    ("N", 1012, 1710, "keep", ""),
    ("T", 1012, 1755, "keep", ""),
    ("H", 1044, 1634, "keep", ""),
    ("T", 1045, 1670, "keep", "bar joins the next T"),
    ("T", 1071, 1671, "keep", ""),
    ("A", 1071, 1755, "keep", "isolated"),
    ("H", 1095, 1609, "keep", ""),
    ("N", 1095, 1643, "keep", "best template fit (0.98, next 0.96)"),
    ("T", 1117, 1611, "keep", ""),
    ("T", 1137, 1559, "keep", ""),
    ("T", 1157, 1569, "keep", ""),
    ("N", 1175, 1576, "keep", ""),
    ("F", 1213, 1567, "keep", "no bottom bar (not E), no bowl (not P)"),
    ("T", 1218, 1516, "keep", ""),
    ("N", 1245, 1550, "keep", ""),
    ("O", 1256, 1608, "keep", "slashed O; overlaps the T below"),
    ("T", 1256, 1620, "keep", "bar inside the O"),
    ("R", 1268, 1484, "keep", "no bottom bar at rows 39-44 (not B)"),
    ("J", 1268, 1571, "keep", ""),
    ("T", 1288, 1516, "keep", ""),
    ("T", 1350, 1470, "keep", "overlaps the F"),
    ("J", 1352, 1505, "keep", ""),
    ("F", 1365, 1457, "keep", ""),
    ("T", 1415, 1452, "keep", "top bar 26 px wide (T), not 16 px (I)"),
    ("V", 1416, 1475, "keep", ""),
    ("J", 1419, 1463, "infer", "top of the J stem visible right of the T stem; the hook lies under the curve and the V"),
    ("O", 1416, 1531, "keep", "isolated"),
    ("A", 1459, 1433, "keep", "crossbar at 51-64 % of the height"),
    ("J", 1470, 1430, "keep", "stem and curled hook end visible"),
    ("F", 1478, 1394, "keep", ""),
    ("F", 1582, 1356, "keep", ""),
    ("B", 1630, 1457, "keep", "isolated"),
    ("O", 1674, 1348, "keep", "overlaps the F above"),
    ("F", 1677, 1326, "keep", ""),
    ("D", 1725, 1316, "keep", "overlaps the O"),
    ("O", 1741, 1313, "keep", ""),
    ("I", 1754, 1433, "keep", "isolated"),
    ("G", 1792, 1348, "keep", "10 keV stack; inner bar visible"),
    ("M", 1792, 1420, "keep", "10 keV stack; central V visible"),
    ("W", 1792, 1444, "keep", "10 keV stack; inverted V from the bottom visible"),
]

# Symbols in the energy range that were not stored (x px, y px, reason).
CU_DROPPED = [
    (934, 1777, "250 eV, Y ~ 1.1: N, M or H overlapped by a T stem and the curve; not resolved"),
    (977, 1680, "300 eV, Y ~ 1.6-1.9: two overlapping glyphs (V, M or A); no template fits"),
    (1138, 1599, "600 eV, Y ~ 2.2: R or B; B's bottom bar would lie under the curve"),
    (1257, 1512, "1 keV, Y ~ 3: B overlapped by other glyphs; not resolved"),
    (1525, 1430, "3-4 keV, Y ~ 4.5-5.5: four or more overlapping glyphs (G, U, S, H, M?); not resolved"),
    (1585, 1420, "4-5 keV, Y ~ 5: several overlapping glyphs; not resolved"),
    (1631, 1348, "5 keV, Y ~ 6.4: B, or F + D (their union is B); not resolved"),
    (1631, 1377, "5 keV, Y ~ 5.7: C or L partly under the curve; not resolved"),
    (1793, 1279, "10 keV stack, Y ~ 8.5: B or O over a second glyph; not resolved"),
    (1793, 1299, "10 keV stack, Y ~ 7.8: O, G or D; not resolved"),
    (1792, 1404, "10 keV stack, Y ~ 5.1: S, F or P; not resolved"),
]

# AM-32 letter -> (dataset id, original reference, DOI or None, NIFS-DATA-23 letter or None; the last is
# informational only: yamamura1995_nifs23.py keeps its own map)
CU_REFERENCES = {
    "A": ("ar_cu_sputter_keywell1955", "F. Keywell, Phys. Rev. 97, 1611 (1955)", "10.1103/PhysRev.97.1611", "A"),
    "B": ("ar_cu_sputter_guseva1960", "M. I. Guseva, Sov. Phys. Solid State 1, 1410 (1960)", None, "B"),
    "D": ("ar_cu_sputter_yonts1960", "O. C. Yonts, C. E. Normand and D. E. Harrison Jr., J. Appl. Phys. 31, 447 (1960)",
          "10.1063/1.1735605", "D"),
    "F": ("ar_cu_sputter_bader1961", "M. Bader, F. C. Witteborn and T. W. Snouse, NASA Technical Report R-105 (1961)",
          None, "Z"),
    "G": ("ar_cu_sputter_fert1961", "C. Fert, N. Colombie, B. Fagot and P. V. Chuong, Ionic Bombardment, Bellevue "
          "(1961), p. 67 (as cited by the compilation)", None, "F"),
    "H": ("ar_cu_sputter_laegreid1961", "N. Laegreid and G. K. Wehner, Trans. 6th Natl. Vacuum Symp. (Pergamon, 1959), "
          "p. 164; J. Appl. Phys. 32, 365 (1961)", "10.1063/1.1736012", "G"),
    "I": ("ar_cu_sputter_perovic1961", "B. Perovic and B. Cobic, Proc. 5th Int. Conf. Ionization Phenomena in Gases "
          "(Munich, 1961), p. 1165", None, "H"),
    "J": ("ar_cu_sputter_southern1963", "A. L. Southern, W. R. Willis and M. T. Robinson, J. Appl. Phys. 34, 153 (1963)",
          "10.1063/1.1729057", "J"),
    "M": ("ar_cu_sputter_dupp1966", "G. Dupp and A. Scharmann, Z. Phys. 192, 284 (1966)", "10.1007/BF01325803", "L"),
    "N": ("ar_cu_sputter_weijsenfeld1967", "C. H. Weijsenfeld, Philips Res. Rep. Suppl. No. 2 (1967)", None, "P"),
    "O": ("ar_cu_sputter_koshkin1969", "V. K. Koshkin, J. A. Rysov, I. I. Shkarban and B. M. Gourmin, Proc. 9th Int. "
          "Conf. Phenomena in Ionized Gases (Bucharest, 1969), p. 92", None, "N"),
    "R": ("ar_cu_sputter_oechsner1973", "H. Oechsner, Z. Phys. 261, 37 (1973)", "10.1007/BF01402280", "U"),
    # Fourth author is Sukenobu per Crossref (10.3131/jvsj.20.161); the compilation misprints it.
    "T": ("ar_cu_sputter_akaishi1977", "K. Akaishi, A. Miyahara, Z. Kabeya, S. Sukenobu, M. Komizo and T. Gotoh, "
          "J. Vac. Soc. Japan 20, 161 (1977)", "10.3131/jvsj.20.161", "O"),
    "V": ("ar_cu_sputter_bohdansky1980", "J. Bohdansky, J. Nucl. Mater. 93-94, 44 (1980)",
          "10.1016/0022-3115(80)90302-5", "Q"),
    "W": ("ar_cu_sputter_okajima1981", "Y. Okajima, Jpn. J. Appl. Phys. 20, 2313 (1981)", "10.1143/JJAP.20.2313", "V"),
}
CU_LEGEND_NAMES = {
    "A": "KEYWELL (1955)", "B": "GUSEVA (1960)", "D": "YONTS, NORMAND, HARRISON (1960)",
    "F": "BADER, WITTEBORN, SNOUSE (1961)", "G": "FERT, COLOMBIE, FAGOT (1961)", "H": "LAEGREID, WEHNER (1959, 1961)",
    "I": "PEROVIC, COBIC (1961)", "J": "SOUTHERN, WILLIS, ROBINSON (1963)", "M": "DUPP, SCHARMANN (1966)",
    "N": "WEIJSENFELD (1967)", "O": "KOSHKIN, RYSOV, SHKARBAN (1969)", "R": "OECHSNER (1973)",
    "T": "AKAISHI, MIYAHARA, KABEYA (1977)", "V": "BOHDANSKY (1980)", "W": "OKAJIMA (1981)",
}
CU = {
    "page": 118,
    "frame_col_frac": 0.4,
    "min_ticks": 35,
    "legend_letters": "ABCDEFGHIJKLMNOPQRSTUVW",
    "legend_box": (1480, 1600, 1600, 3400),  # x0, x1, y0, y1 of the legend letter column
    "legend_mask": (1620, 1495),  # legend block: y >= 1620 and x >= 1495 (no data there in range)
    "scale": (0.89, 0.957),  # plotted / legend glyph size (step 2)
    "isolated_check": [("T", 975, 1804), ("V", 812, 1800), ("A", 1071, 1755), ("I", 1754, 1433)],
    "curve_anchors": CU_CURVE_ANCHORS,
    "seeds": CU_SEEDS,
    "dropped": CU_DROPPED,
    "references": CU_REFERENCES,
    "legend_names": CU_LEGEND_NAMES,
    "target_state": TARGET_STATE,
}

# Crossref was queried (2026-10-06) for every original reference below with a
# journal citation; a DOI is given only where Crossref returns the cited
# journal, volume and first page, and author names follow Crossref where the
# compilation misspells them. "No DOI found" means Crossref has no record.
LAEGREID = ("N. Laegreid and G. K. Wehner, Trans. 6th Natl. Vacuum Symp. (Pergamon, 1959), p. 164; J. Appl. Phys. 32, "
            "365 (1961)", "10.1063/1.1736012")
NOT_STATED = TARGET_STATE
SI_STATE = (
    "not stated in the compilation for this set, nor whether the surface was crystalline or amorphized by the Ar "
    "dose; the compilation excludes single-crystal targets and states that most specimens were thinned "
    "polycrystalline foils or evaporated films (IPPJ-AM-32, p. 2)"
)

# ======================================================================
# Ar -> Si, PDF p. 47 (printed p. 42); issue #70. Legend A-G; in 0.2-10 keV
# only A, E, F and G have points (B, C, D lie above 20 keV).
# ======================================================================
SI_CURVE_ANCHORS = [
    (864, 2358), (906, 2268), (960, 2172), (1029, 2069), (1080, 2006), (1164, 1916), (1200, 1883), (1260, 1830),
    (1320, 1787), (1380, 1746), (1440, 1711), (1500, 1678), (1560, 1648), (1620, 1623), (1680, 1601), (1740, 1581),
    (1800, 1564), (1860, 1548), (1920, 1536), (1980, 1526), (2040, 1518), (2100, 1514), (2160, 1514),
]
# (letter, x px, y px, status, note, glyph box x0, y0, x1, y1 for the second read)
SI_SEEDS = [
    ("A", 1015, 2031, "keep", "crossbar at 50-60 % of the height", (1001, 2010, 1029, 2052)),
    ("A", 1015, 2090, "keep", "crossed by the curve; crossbar at 51-60 %", (1001, 2068, 1030, 2113)),
    ("A", 1066, 2030, "keep", "crossed by the curve; crossbar at 50-60 %", (1052, 2008, 1080, 2052)),
    ("A", 1110, 1957, "keep", "crossbar at 53 %; overlaps the F below", (1096, 1935, 1124, 1978)),
    ("F", 1110, 1987, "keep", "top bar (rows 1966-1970) across the A's legs, short middle bar at 37 %, stem to y 2008; "
     "no right stem below the A", (1097, 1966, 1123, 2008)),
    ("A", 1148, 1908, "keep", "curve crosses its right leg", (1134, 1885, 1162, 1930)),
    ("A", 1180, 1950, "keep", "isolated", (1165, 1928, 1194, 1972)),
    ("F", 1231, 1856, "keep", "top bar 10 px above the A's top, middle bar at 41 %, no bottom bar at y 1876-1878 (not E)",
     (1219, 1834, 1245, 1878)),
    ("A", 1233, 1866, "keep", "rounded top, crossbar at 51-60 %; under the F", (1219, 1844, 1246, 1889)),
    ("A", 1275, 1826, "keep", "crossed by the curve; crossbar at 50-60 %", (1260, 1803, 1290, 1848)),
    ("F", 1329, 1840, "keep", "isolated", (1317, 1817, 1341, 1863)),
    ("E", 1372, 1830, "keep", "isolated; bottom bar", (1359, 1808, 1386, 1853)),
    ("F", 1396, 1764, "keep", "top bar touches the curve; no bottom bar", (1384, 1742, 1411, 1786)),
    ("F", 1496, 1754, "keep", "isolated", (1483, 1732, 1508, 1777)),
    ("F", 1560, 1700, "keep", "isolated", (1546, 1678, 1573, 1723)),
    ("F", 1654, 1682, "keep", "isolated", (1642, 1659, 1667, 1704)),
    ("G", 1656, 1626, "keep", "top under the curve; inner bar visible (not C)", (1642, 1604, 1671, 1647)),
    ("G", 1780, 1606, "keep", "isolated; inner bar", (1767, 1584, 1794, 1628)),
    ("G", 1943, 1606, "keep", "isolated; inner bar", (1929, 1584, 1957, 1628)),
]
SI_DROPPED = []  # every glyph in 196 eV - 10.2 keV was identified
SI_REFERENCES = {
    "A": ("ar_si_sputter_laegreid1961", *LAEGREID, "A"),
    "E": ("ar_si_sputter_poate1976", "J. M. Poate, W. L. Brown, R. Homer, W. M. Augustyniak, J. W. Mayer, K. N. Tu and "
          "W. F. van der Weg, Nucl. Instrum. Methods 132, 345 (1976)", "10.1016/0029-554x(76)90756-4", "K"),
    "F": ("ar_si_sputter_coburn1977", "J. W. Coburn, H. F. Winters and T. J. Chuang, J. Appl. Phys. 48, 3532 (1977)",
          "10.1063/1.324150", "F"),
    "G": ("ar_si_sputter_kang1979", "S. T. Kang, R. Shimizu and T. Okutani, Jpn. J. Appl. Phys. 18, 1717 (1979)",
          "10.1143/jjap.18.1717", "G"),
}
SI_LEGEND_NAMES = {
    "A": "LAEGREID, WEHNER (1959, 1961)", "E": "POATE, BROWN, HOMER (1976)", "F": "COBURN, WINTERS, CHUANG (1977)",
    "G": "KANG, SHIMIZU, OKUTANI (1979)",
}
SI = {
    "page": 47,
    "frame_col_frac": 0.2,
    "min_ticks": 30,  # the bottom axis is faint: 32 ticks found
    "legend_letters": "ABCDEFG",
    "legend_box": (850, 905, 2520, 3140),
    "legend_mask": (2510, 850),
    "scale": (1.0, 0.957),
    "isolated_check": [("A", 1180, 1950), ("G", 1780, 1606), ("F", 1560, 1700), ("E", 1372, 1830)],
    "curve_anchors": SI_CURVE_ANCHORS,
    "seeds": SI_SEEDS,
    "dropped": SI_DROPPED,
    "references": SI_REFERENCES,
    "legend_names": SI_LEGEND_NAMES,
    "target_state": SI_STATE,
}

# ======================================================================
# Ar -> Ag, PDF p. 191 (printed p. 186); issue #70. Legend A-O.
# ======================================================================
AG_CURVE_ANCHORS = [
    (723, 2035), (780, 1950), (861, 1853), (1000, 1713), (1080, 1650), (1140, 1608), (1200, 1570), (1260, 1532),
    (1320, 1498), (1380, 1466), (1440, 1437), (1500, 1411), (1560, 1384), (1620, 1359), (1680, 1337), (1740, 1316),
    (1800, 1298), (1860, 1281), (1920, 1266), (1980, 1254), (2040, 1242), (2100, 1233), (2160, 1226),
]
# On this page the H crossbar sits at 35-47 % of the height, the F middle bar
# at 31-48 % (shorter than the glyph), the A crossbar at 50-62 % under a
# rounded top; F and H drawn on one spot therefore look like an H with a top
# bar, told apart by the top bar not joining the right stem or by the second,
# shorter bar.
AG_SEEDS = [
    ("B", 1017, 1732, "keep", "crossed by the curve; middle bar visible", (1002, 1709, 1032, 1754)),
    ("F", 1017, 1655, "keep", "drawn on the H: flat top bar and the short middle bar at 45-50 % (rows 1653-1655)",
     (1002, 1633, 1028, 1677)),
    ("H", 1017, 1655, "keep", "drawn with the F: right stem from top to bottom and crossbar at 34-43 %",
     (1002, 1633, 1032, 1677)),
    ("B", 1061, 1652, "keep", "rounded top right corner and middle bar at 39-50 % (not D); bottom bar under the curve",
     (1047, 1630, 1075, 1674)),
    ("B", 1061, 1668, "infer", "B or D (its middle lies under the curve): only the top, right stem and bottom bar are "
     "visible; D (Almen and Bruce 1961A) appears only at 10 keV and above in this figure, while the B series "
     "(Koedam, 40-240 eV per the original's title) runs through this energy", (1047, 1646, 1075, 1690)),
    ("B", 1040, 1688, "infer", "B or D (its middle lies under the curve), same argument as the B above",
     (1027, 1666, 1053, 1710)),
    ("F", 1044, 1603, "keep", "isolated", (1031, 1580, 1058, 1626)),
    ("F", 1069, 1554, "keep", "isolated", (1056, 1532, 1082, 1576)),
    ("F", 1093, 1587, "keep", "top bar touches the H's stem", (1080, 1565, 1106, 1608)),
    ("H", 1115, 1566, "keep", "crossbar at 27-41 %, both stems", (1101, 1544, 1129, 1588)),
    ("F", 1112, 1577, "keep", "top bar on the H's crossbar, short middle bar (rows 1571-1577), stem to y 1597 below "
     "the H", (1101, 1556, 1126, 1597)),
    ("F", 1152, 1560, "keep", "no right stem below the next F's stem (not H or A)", (1139, 1537, 1163, 1582)),
    ("F", 1169, 1547, "keep", "no right stem below y 1548 (not H or A)", (1156, 1524, 1182, 1570)),
    ("F", 1187, 1526, "keep", "drawn on the H: its flat top bar ends at x 1194, apart from the H's right stem",
     (1173, 1504, 1196, 1548)),
    ("H", 1187, 1526, "keep", "drawn with the F: right stem (x 1196-1201) separate from the F's top bar", (1173, 1504, 1202, 1548)),
    ("H", 1240, 1474, "keep", "isolated; crossbar at 33-44 %", (1226, 1451, 1253, 1496)),
    ("M", 1240, 1527, "keep", "crossed by the curve; diagonals visible", (1226, 1506, 1253, 1549)),
    ("H", 1286, 1465, "keep", "crossbar at 39-45 %; stems to y 1487", (1270, 1443, 1301, 1487)),
    ("F", 1283, 1473, "keep", "top bar 8 px below the H's top, short middle bar, stem to y 1495 below the H",
     (1270, 1451, 1296, 1495)),
    ("A", 1307, 1403, "keep", "isolated", (1291, 1381, 1323, 1425)),
    ("J", 1410, 1322, "keep", "isolated", (1396, 1300, 1424, 1344)),
    ("L", 1421, 1394, "keep", "isolated", (1407, 1372, 1436, 1415)),
    ("M", 1410, 1460, "keep", "under the curve; the left diagonal and the V's bottom at mid-height are visible, no "
     "crossbar (not H) and no diagonal to the lower right (not N)", (1395, 1439, 1424, 1481)),
    ("A", 1493, 1333, "keep", "crossbar at 52-64 %", (1479, 1311, 1507, 1355)),
    ("M", 1510, 1373, "keep", "under the A", (1497, 1352, 1524, 1395)),
    ("A", 1622, 1318, "keep", "isolated", (1608, 1295, 1637, 1340)),
    ("C", 1674, 1291, "keep", "isolated", (1660, 1268, 1688, 1313)),
    ("A", 1711, 1266, "keep", "isolated", (1696, 1243, 1726, 1288)),
    ("A", 1770, 1268, "keep", "touches the C", (1758, 1247, 1782, 1290)),
    ("C", 1797, 1246, "keep", "", (1785, 1225, 1810, 1268)),
    ("A", 1824, 1224, "keep", "overlaps the next A", (1810, 1203, 1838, 1245)),
    ("A", 1846, 1226, "keep", "overlaps the previous A", (1832, 1205, 1860, 1248)),
    ("C", 1878, 1244, "keep", "isolated", (1865, 1221, 1891, 1267)),
    ("O", 1963, 1156, "keep", "top of the 10 keV stack; slashed O", (1949, 1134, 1977, 1177)),
]
AG_DROPPED = [
    (1963, 1195, "10 keV stack below the O, Y ~ 12: an O, D or B (bars at y 1197-1201, a slash at 1202-1211); not resolved"),
    (1963, 1219, "10 keV stack, Y ~ 11: B or D (bars at y 1211-1219 and 1235-1239); not resolved"),
    (1963, 1234, "10 keV stack bottom, Y ~ 10: B, D or E partly under the curve; not resolved"),
]
AG_REFERENCES = {
    "A": ("ar_ag_sputter_keywell1955", "F. Keywell, Phys. Rev. 97, 1611 (1955)", "10.1103/PhysRev.97.1611", "A"),
    "B": ("ar_ag_sputter_koedam1958", "M. Koedam, Physica 24, 692 (1958)", "10.1016/s0031-8914(58)80083-x", None),
    "C": ("ar_ag_sputter_guseva1960", "M. I. Guseva, Sov. Phys. Solid State 1, 1410 (1960)", None, "B"),
    "F": ("ar_ag_sputter_laegreid1961", *LAEGREID, "E"),
    "H": ("ar_ag_sputter_wehner1961", "G. K. Wehner and D. Rosenberg (1961), cited by the compilation as J. Appl. Phys. "
          "32, 887 (1961); Crossref gives that paper's title as \"Mercury Ion Beam Sputtering of Metals at Energies "
          "4-15 kev\", which does not match these Ar points at 0.2-0.6 keV, so the publication holding them is "
          "uncertain (the compilation also lists G. K. Wehner, R. V. Stuart and D. Rosenberg, General Mills Report "
          "No. 2243 (1961))", None, "G"),
    "J": ("ar_ag_sputter_benninghoven1969", "A. Benninghoven, Z. Angew. Phys. 27, 51 (1969)", None, "O"),
    "L": ("ar_ag_sputter_oechsner1973", "H. Oechsner, Z. Phys. 261, 37 (1973)", "10.1007/BF01402280", "L"),
    "M": ("ar_ag_sputter_smith1975", "J. N. Smith, C. H. Meyer and J. K. Layton, J. Appl. Phys. 46, 4291 (1975) (the "
          "legend's \"(1975)\"; of the compilation's two 1975 entries for these authors this is the one on Ar+ "
          "sputtering of Ag, per its Crossref title)", "10.1063/1.321449", "I"),
    "O": ("ar_ag_sputter_okajima1981", "Y. Okajima, Jpn. J. Appl. Phys. 20, 2313 (1981)", "10.1143/JJAP.20.2313", "M"),
}
AG_LEGEND_NAMES = {
    "A": "KEYWELL (1955)", "B": "KOEDAM (1958)", "C": "GUSEVA (1960)", "F": "LAEGREID, WEHNER (1959, 1961)",
    "H": "WEHNER, ROSENBERG (1961)", "J": "BENNINGHOVEN (1969)", "L": "OECHSNER (1973)", "M": "SMITH, MEYER, LAYTON (1975)",
    "O": "OKAJIMA (1981)",
}
AG = {
    "page": 191,
    "frame_col_frac": 0.2,
    "min_ticks": 35,
    "legend_letters": "ABCDEFGHIJKLMNO",
    "legend_box": (1440, 1500, 2225, 3390),
    "legend_mask": (2225, 1440),
    "scale": (0.905, 0.989),
    "isolated_check": [("A", 1622, 1318), ("C", 1674, 1291), ("J", 1410, 1322), ("L", 1421, 1394), ("H", 1240, 1474)],
    "curve_anchors": AG_CURVE_ANCHORS,
    "seeds": AG_SEEDS,
    "dropped": AG_DROPPED,
    "references": AG_REFERENCES,
    "legend_names": AG_LEGEND_NAMES,
    "target_state": {**{L: NOT_STATED for L in "ACFHJLMO"},
                     "B": "polycrystalline, per the original's title (Crossref): \"Sputtering of a polycristalline silver "
                          "surface bombarded with monoenergetic argon ions of low energy (40-240 eV)\"",
                     "L": "polycrystalline, per the original's title (Crossref; oblique-incidence study of polycrystalline "
                          "targets)"},
}

# ======================================================================
# Ar -> Au, PDF p. 262 (printed p. 257); issue #70. Legend A-T.
# ======================================================================
AU_CURVE_ANCHORS = [
    (600, 2042), (660, 1951), (720, 1872), (801, 1781), (840, 1742), (900, 1690), (960, 1636), (1020, 1591),
    (1080, 1550), (1140, 1510), (1200, 1476), (1260, 1442), (1320, 1412), (1380, 1382), (1440, 1357), (1515, 1329),
    (1560, 1312), (1620, 1292), (1680, 1272), (1740, 1257), (1800, 1242), (1860, 1230), (1920, 1218), (1980, 1210),
]
# On this page T's top bar is about 23-27 px wide and I's top serif about
# 15 px; T (Yamashita 1982) and S (Yamashita 1980) are often drawn on one spot.
AU_SEEDS = [
    ("B", 857, 1654, "keep", "isolated", (842, 1631, 872, 1676)),
    ("B", 886, 1606, "keep", "isolated", (871, 1584, 900, 1628)),
    ("S", 910, 1736, "keep", "isolated", (896, 1714, 923, 1758)),
    ("F", 954, 1535, "keep", "top bar 13 px above the B's top with a square end, no bottom bar at y 1552-1557 (not E)",
     (941, 1513, 965, 1557)),
    ("B", 955, 1547, "keep", "under the F; middle bar visible", (941, 1526, 968, 1569)),
    ("S", 988, 1764, "keep", "isolated", (974, 1742, 1003, 1786)),
    ("B", 1022, 1500, "keep", "isolated", (1007, 1478, 1036, 1523)),
    ("B", 1074, 1455, "keep", "isolated", (1061, 1431, 1087, 1478)),
    ("S", 1074, 1606, "keep", "isolated", (1061, 1583, 1086, 1628)),
    ("F", 1116, 1430, "keep", "isolated", (1102, 1407, 1129, 1452)),
    ("T", 1170, 1540, "keep", "top bar 22 px wide (T, not I); over the S", (1157, 1518, 1182, 1562)),
    ("S", 1170, 1556, "keep", "under the T", (1157, 1534, 1184, 1578)),
    ("F", 1186, 1414, "keep", "isolated", (1173, 1392, 1198, 1437)),
    ("C", 1237, 1637, "keep", "isolated", (1224, 1615, 1250, 1659)),
    ("T", 1237, 1480, "infer", "centred stem over the S; the left half of the top bar lies under the curve, so T and "
     "I (Sletten, whose points lie near 10 keV and above in this figure) cannot be told apart by the bar width; T by the "
     "T-over-S pairs at 0.75, 1.5 and 2 keV (the Yamashita sets)", (1225, 1458, 1250, 1502)),
    ("S", 1238, 1513, "keep", "under the T", (1225, 1491, 1251, 1534)),
    ("G", 1240, 1274, "keep", "isolated; inner bar", (1226, 1252, 1253, 1297)),
    ("J", 1250, 1332, "keep", "hook on the F stack", (1238, 1310, 1264, 1352)),
    ("F", 1237, 1394, "keep", "1 keV stack, middle glyph: middle bar at y 1389-1393, no bottom bar at y 1411-1416",
     (1225, 1372, 1249, 1416)),
    ("F", 1237, 1402, "keep", "1 keV stack, bottom glyph: middle bar at y 1397-1401, no bottom bar at y 1419-1424",
     (1225, 1380, 1249, 1424)),
    ("S", 1289, 1479, "keep", "isolated", (1275, 1458, 1303, 1500)),
    ("E", 1335, 1294, "keep", "isolated", (1321, 1272, 1349, 1316)),
    ("T", 1334, 1458, "keep", "top bar 23 px wide (T, not I) and centred stem, drawn on the S", (1321, 1436, 1346, 1479)),
    ("S", 1335, 1458, "keep", "drawn with the T", (1321, 1436, 1349, 1479)),
    ("S", 1371, 1447, "keep", "isolated", (1358, 1425, 1384, 1469)),
    ("E", 1402, 1252, "keep", "isolated", (1389, 1230, 1416, 1274)),
    ("T", 1402, 1405, "keep", "top bar at least 19 px wide outside the curve (T, not I); over the S",
     (1390, 1383, 1415, 1427)),
    ("S", 1402, 1430, "keep", "under the T", (1388, 1409, 1416, 1452)),
    ("C", 1402, 1534, "keep", "isolated", (1388, 1511, 1415, 1556)),
    ("P", 1403, 1168, "keep", "isolated", (1389, 1145, 1417, 1191)),
    ("E", 1458, 1224, "keep", "isolated", (1444, 1202, 1471, 1246)),
    ("E", 1500, 1206, "keep", "isolated", (1487, 1183, 1513, 1228)),
    ("C", 1500, 1349, "keep", "top under the curve; open right side, no inner bar (not G or O)", (1486, 1327, 1513, 1370)),
    ("E", 1536, 1194, "keep", "isolated", (1523, 1172, 1549, 1217)),
    ("E", 1566, 1190, "keep", "isolated", (1555, 1168, 1578, 1212)),
    ("C", 1567, 1255, "keep", "isolated", (1554, 1233, 1580, 1277)),
    ("E", 1595, 1180, "keep", "three bars; its top bar meets the next E's stem", (1582, 1158, 1608, 1201)),
    ("E", 1620, 1170, "keep", "three bars", (1608, 1148, 1633, 1191)),
    ("C", 1620, 1220, "keep", "isolated", (1607, 1198, 1634, 1241)),
    ("C", 1662, 1160, "keep", "isolated", (1650, 1138, 1675, 1183)),
    ("Q", 1664, 1032, "keep", "isolated; tail at the lower right", (1651, 1010, 1680, 1054)),
    ("C", 1768, 1160, "keep", "10 keV group, left: open right side at mid-height, no inner bar (not G or O)",
     (1754, 1137, 1781, 1183)),
    ("A", 1782, 1152, "keep", "10 keV group: rounded top, crossbar at 51-62 %", (1767, 1129, 1796, 1173)),
    ("I", 1780, 1306, "keep", "isolated; top serif 17 px wide", (1771, 1283, 1788, 1328)),
]
AU_DROPPED = [
    (1238, 1360, "1 keV stack, top glyph (Y ~ 3.5): E, F or P (its bottom bar or bowl would lie in the merged bars and "
     "under the J's hook); not resolved"),
    (1712, 1185, "7.4 keV, Y ~ 7.4: a rounded A-like top with a crossbar at 51-62 % over a bottom bar (A drawn on a D or "
     "B?); not resolved"),
]
AU_REFERENCES = {
    "A": ("ar_au_sputter_almen1961", "O. Almen and G. Bruce, Nucl. Instrum. Methods 11, 257 (1961)",
          "10.1016/0029-554x(61)90026-x", "A"),
    "B": ("ar_au_sputter_laegreid1961", *LAEGREID, "B"),
    "C": ("ar_au_sputter_patterson1962", "H. Patterson and D. H. Tomlin, Proc. R. Soc. London A 265, 474 (1962)",
          "10.1098/rspa.1962.0037", "D"),
    "E": ("ar_au_sputter_robinson1967", "M. T. Robinson and A. L. Southern, J. Appl. Phys. 38, 2969 (1967)",
          "10.1063/1.1710034", "M"),
    "F": ("ar_au_sputter_weijsenfeld1967", "C. H. Weijsenfeld, Philips Res. Rep. Suppl. No. 2 (1967)", None, "N"),
    "G": ("ar_au_sputter_benninghoven1969", "A. Benninghoven, Z. Angew. Phys. 27, 51 (1969)", None, "Y"),
    "I": ("ar_au_sputter_sletten1972", "G. Sletten and P. Knudsen, Nucl. Instrum. Methods 102, 459 (1972)",
          "10.1016/0029-554x(72)90633-7", "L"),
    "J": ("ar_au_sputter_oechsner1973", "H. Oechsner, Z. Phys. 261, 37 (1973)", "10.1007/BF01402280", "U"),
    "P": ("ar_au_sputter_holloway1977", "P. H. Holloway, Surf. Sci. 66, 479 (1977)", "10.1016/0039-6028(77)90033-4", "V"),
    "Q": ("ar_au_sputter_szymonski1978", "M. Szymonski, R. S. Bhattacharya, H. Overeijnder and A. E. de Vries, J. Phys. "
          "D: Appl. Phys. 11, 751 (1978) (Crossref title: \"Sputtering of an AgAu alloy by bombardment with 6 keV Xe+ "
          "ions\"; the compilation plots this point as Ar -> Au, and whether the paper reports that yield was not "
          "checked)", "10.1088/0022-3727/11/5/018", "Z"),
    "S": ("ar_au_sputter_yamashita1980", "M. Yamashita, S. Baba and A. Kinbara, Proc. 4th Symp. Ion Sources and Ion "
          "Application Technology (Tokyo, 1980), p. 311", None, "T"),
    "T": ("ar_au_sputter_yamashita1982", "M. Yamashita, S. Baba and A. Kinbara, J. Vac. Soc. Japan 25, 249 (1982)",
          "10.3131/jvsj.25.249", "^"),
}
AU_LEGEND_NAMES = {
    "A": "ALMEN, BRUCE (1961A)", "B": "LAEGREID, WEHNER (1959, 1961)", "C": "PATTERSON, TOMLIN (1962)",
    "E": "ROBINSON, SOUTHERN (1967)", "F": "WEIJSENFELD (1967)", "G": "BENNINGHOVEN (1969)",
    "I": "SLETTEN, KNUDSEN (1972)", "J": "OECHSNER (1973)", "P": "HOLLOWAY (1977)",
    "Q": "SZYMONSKI, BHATTACHARYA, OVEREIJNDER (1978)", "S": "YAMASHITA, BABA, KINBARA (1980)",
    "T": "YAMASHITA, BABA, KINBARA (1982)",
}
AU = {
    "page": 262,
    "frame_col_frac": 0.2,
    "min_ticks": 35,
    "legend_letters": "ABCDEFGHIJKLMNOPQRST",
    "legend_box": (1415, 1475, 1650, 3280),
    "legend_mask": (1650, 1415),
    "scale": (0.933, 1.0),
    "isolated_check": [("E", 1458, 1224), ("C", 1567, 1255), ("S", 910, 1736), ("B", 1022, 1500), ("G", 1240, 1274)],
    "curve_anchors": AU_CURVE_ANCHORS,
    "seeds": AU_SEEDS,
    "dropped": AU_DROPPED,
    "references": AU_REFERENCES,
    "legend_names": AU_LEGEND_NAMES,
    "target_state": {**{L: NOT_STATED for L in "ABCFGIJPQST"},
                     "E": "monocrystalline, per the original's title (Crossref): \"Sputtering Experiments with 1- to "
                          "5-keV Ar+ Ions. II. Monocrystalline Targets of Al, Cu, and Au\"; the compilation shows the set "
                          "although it states that it excludes single-crystal targets (p. 2); orientation not recorded "
                          "here",
                     "J": "polycrystalline, per the original's title (Crossref; oblique-incidence study of polycrystalline "
                          "targets)"},
}

FIGURES = {"Cu": CU, "Si": SI, "Ag": AG, "Au": AU}


GRID = np.log10(np.array([k * 10.0**n for n in range(-4, 8) for k in range(1, 10)]))
LONG = np.log10(np.array([k * 10.0**n for n in range(-4, 8) for k in (1, 5)]))  # long ticks: 10^n and 5 x 10^n


# --- 1. axes ---

def frame(ink, col_frac):
    """Approximate frame sides: rows with ink over 40 % of the page width, and
    columns with ink over `col_frac` of the page height (per figure)."""
    h, w = ink.shape
    rows, cols = ink.sum(1), ink.sum(0)
    ys = [i for i in range(h) if rows[i] > 0.4 * w]
    xs = [j for j in range(w) if cols[j] > col_frac * h]
    return (np.mean([y for y in ys if y < h / 2]), np.mean([y for y in ys if y > h / 2]),
            np.mean([x for x in xs if x < w / 2]), np.mean([x for x in xs if x > w / 2]))


def side_line(ink, horizontal, approx, a, b, half=14):
    """Centre line of one frame side as perp = s * along + c."""
    pts = []
    for t in range(int(a), int(b), 20):
        seg = ink[int(approx) - half:int(approx) + half, t] if horizontal else ink[t, int(approx) - half:int(approx) + half]
        idx = np.nonzero(seg)[0]
        if 3 <= len(idx) <= 14 and idx.max() - idx.min() < 14:
            pts.append((t, int(approx) - half + idx.mean()))
    pts = np.array(pts)
    s, c = np.polyfit(pts[:, 0], pts[:, 1], 1)
    keep = np.abs(pts[:, 1] - (s * pts[:, 0] + c)) < 2.5
    return np.polyfit(pts[keep, 0], pts[keep, 1], 1)


def side_ticks(ink, horizontal, line, a, b, inward, start=6, minlen=15):
    """Ticks protruding inward from a side: [(along, perp, length)], split ticks merged."""
    s, c = line
    lens = []
    for t in range(int(a), int(b)):
        p0, n = s * t + c, 0
        while n < 80:
            p = int(round(p0 + inward * (start + n)))
            if not (ink[p, t] if horizontal else ink[t, p]):
                break
            n += 1
        lens.append(n)
    lens = np.array(lens)
    out, i = [], 0
    while i < len(lens):
        if lens[i] >= minlen:
            j = i
            while j + 1 < len(lens) and lens[j + 1] >= minlen:
                j += 1
            seg = np.arange(i, j + 1)
            full = seg[lens[seg] >= 0.7 * lens[seg].max()]
            cen = a + full.mean()
            out.append([cen, s * cen + c, int(lens[seg].max()), j - i + 1])
            i = j + 1
        else:
            i += 1
    merged = []
    for t in out:
        if merged and t[0] - merged[-1][0] < 10:
            p = merged[-1]
            w0, w1 = p[3], t[3]
            merged[-1] = [(p[0] * w0 + t[0] * w1) / (w0 + w1), (p[1] * w0 + t[1] * w1) / (w0 + w1), max(p[2], t[2]), w0 + w1]
        else:
            merged.append(t)
    return merged


def assign(name, ticks, end_a, val_a, end_b, val_b):
    """Value (log10) of each tick. A run that reaches the 80 px scan limit is a
    frame line, not a tick, and is skipped. Long ticks (at least 0.8 of the
    longest) mark 10^n and 5 x 10^n: each takes the nearest such value on the
    straight line through the frame corners' values (val_a at end_a, val_b at
    end_b), then a line is fitted to the long ticks, the worst dropped while
    it misses its value by more than 0.02 decade. Every tick then takes the
    nearest 1..9 x 10^n on that line; a tick more than 0.02 decade from it is
    not a tick of the comb (e.g. a symbol stroke) and is dropped. Both
    rejections are printed."""
    ticks = [t for t in ticks if t[2] < 80]
    lmax = max(t[2] for t in ticks)
    longs = [t for t in ticks if t[2] >= 0.8 * lmax]
    xs = np.array([t[0] for t in longs])
    guess = val_a + (val_b - val_a) * (xs - end_a) / (end_b - end_a)
    vals = LONG[np.abs(guess[:, None] - LONG[None, :]).argmin(1)]
    keep = np.ones(len(xs), bool)
    while True:
        p = np.polyfit(xs[keep], vals[keep], 1)
        miss = np.where(keep, np.abs(np.polyval(p, xs) - vals), 0.0)
        if miss.max() <= 0.02:
            break
        print(f"  {name} axis: long tick at {xs[miss.argmax()]:.0f} px dropped (misses its value by {miss.max():.3f} decade)")
        keep[miss.argmax()] = False
    out = []
    for t in ticks:
        v = np.polyval(p, t[0])
        g = GRID[np.argmin(np.abs(GRID - v))]
        if abs(v - g) > 0.02:
            print(f"  {name} axis: tick at {t[0]:.0f} px dropped ({v - g:+.3f} decade from the comb)")
            continue
        out.append((t[0], t[1], g))
    return out


def calibrate(ink, fig):
    t, b, l, r = frame(ink, fig["frame_col_frac"])
    lines = {
        "bottom": (True, side_line(ink, True, b, l + 30, r - 30), -1),
        "top": (True, side_line(ink, True, t, l + 30, r - 30), +1),
        "left": (False, side_line(ink, False, l, t + 30, b - 30), +1),
        "right": (False, side_line(ink, False, r, t + 30, b - 30), -1),
    }
    (e_lo, e_hi), (y_lo, y_hi) = FRAME_E, FRAME_Y
    X, Y = [], []
    for name, (hz, line, inward) in lines.items():
        a, bb = (l + 10, r - 10) if hz else (t + 10, b - 10)
        tk = side_ticks(ink, hz, line, a, bb, inward)
        if hz:  # x axis: the frame runs from 50 eV (left) to 10^6 eV (right)
            got = assign(name, tk, l, math.log10(e_lo), r, math.log10(e_hi))
            X += [(al, pp, v) for al, pp, v in got]
        else:  # y axis: from 10^-3 (bottom) to 50 (top)
            got = assign(name, tk, b, math.log10(y_lo), t, math.log10(y_hi))
            Y += [(pp, al, v) for al, pp, v in got]
        if len(got) < fig["min_ticks"]:
            sys.exit(f"{name} axis: {len(got)} ticks, expected at least {fig['min_ticks']}; wrong page or raster?")
    X, Y = np.array(X), np.array(Y)
    A = np.c_[np.ones(len(X)), X[:, 0], X[:, 1]]
    B = np.c_[np.ones(len(Y)), Y[:, 0], Y[:, 1]]
    cx = np.linalg.lstsq(A, X[:, 2], rcond=None)[0]
    cy = np.linalg.lstsq(B, Y[:, 2], rcond=None)[0]
    rx, ry = A @ cx - X[:, 2], B @ cy - Y[:, 2]
    print(f"axes: {len(X)} x ticks, {len(Y)} y ticks; residual rms {rx.std():.4f} / {ry.std():.4f} decade, "
          f"max {np.abs(rx).max():.4f} / {np.abs(ry).max():.4f} decade")
    print(f"  log10 E = {cx[0]:.6f} + {cx[1]:.7e} x + {cx[2]:.3e} y;  log10 Y = {cy[0]:.6f} + {cy[1]:.3e} x + {cy[2]:.7e} y")
    return cx, cy, float(rx.std()), float(ry.std()), len(X) + len(Y)


# --- 2. templates, 3. curve ---

def templates(ink, fig, sx, sy):
    x0, x1, y0, y1 = fig["legend_box"]
    letters = fig["legend_letters"]
    sub = ink[y0:y1, x0:x1]
    rows, runs, cur = sub.sum(1), [], None
    for i, v in enumerate(rows):
        if v > 0 and cur is None:
            cur = i
        if v == 0 and cur is not None:
            runs.append((cur, i - 1))
            cur = None
    if len(runs) != len(letters):
        sys.exit(f"legend: found {len(runs)} letters, expected {len(letters)}")
    G = {}
    for (a, b), L in zip(runs, letters):
        g = sub[a:b + 1]
        cs = np.nonzero(g.sum(0))[0]
        g = g[:, cs.min():cs.max() + 1]
        if sx != 1.0 or sy != 1.0:
            h, w = g.shape
            img = Image.fromarray((g * 255).astype(np.uint8)).resize((round(w * sx), round(h * sy)), Image.BOX)
            g = np.asarray(img) >= 128
        G[L] = g
    return G


def curve_mask(ink, anchors):
    a = np.array(anchors, float)
    cs = CubicSpline(a[:, 0], a[:, 1])
    xs, ys = [], []
    for x in range(int(a[0, 0]), int(a[-1, 0])):
        yp, s = float(cs(x)), abs(float(cs(x, 1)))
        maxrun, best, i, col = 9 * math.sqrt(1 + s * s) + 2, None, int(yp - 15), ink[:, x]
        while i < int(yp + 15):
            if col[i]:
                j = i
                while col[j + 1]:
                    j += 1
                c = (i + j) / 2
                if j - i + 1 <= maxrun and abs(c - yp) <= 4 and (best is None or abs(c - yp) < abs(best - yp)):
                    best = c
                i = j + 1
            else:
                i += 1
        if best is not None:
            xs.append(x)
            ys.append(best)
    xs, ys = np.array(xs, float), np.array(ys)
    keep = np.ones(len(xs), bool)
    for _ in range(4):
        sp = UnivariateSpline(xs[keep], ys[keep], s=keep.sum() * 1.5)
        keep = np.abs(ys - sp(xs)) < 2.5
    mask = np.zeros_like(ink)
    d = sp.derivative()
    for x in range(int(a[0, 0]), int(a[-1, 0])):
        y, h = float(sp(x)), 6.5 * math.sqrt(1 + float(d(x)) ** 2)
        mask[int(math.floor(y - h)):int(math.ceil(y + h)) + 1, x] = True
    print(f"fitted curve masked: {keep.sum()} of {len(xs)} snapped columns, rms {np.std(ys[keep] - sp(xs[keep])):.2f} px")
    return mask


# --- 5. centres ---

def refine(ink, dc, g, x, y, tol, r=5):
    """Centre maximizing template coverage (ink dilated by `tol`, curve = don't care);
    sub-pixel: mean of the positions within 0.005 of the maximum."""
    inkd = ndi.binary_dilation(ink, iterations=tol) | dc
    gh, gw = g.shape
    win = inkd[y - r - gh:y + r + gh + 1, x - r - gw:x + r + gw + 1].astype(float)
    c = fftconvolve(win, g[::-1, ::-1].astype(float), mode="same") / g.sum()
    cy0, cx0 = r + gh, r + gw
    sub = c[cy0 - r:cy0 + r + 1, cx0 - r:cx0 + r + 1]
    best = sub.max()
    yy, xx = np.nonzero(sub >= best - 0.005)
    # 'same' convolution centres the kernel at (gh // 2, gw // 2); move to the geometric glyph centre
    off_y, off_x = (gh - 1) / 2 - gh // 2, (gw - 1) / 2 - gw // 2
    return x - r + xx.mean() + off_x, y - r + yy.mean() + off_y, float(best)


def load(raster, target):
    """Ink of the page, the figure's table, axis maps, templates and curve mask."""
    fig = FIGURES[target]
    ink = np.asarray(Image.open(raster).convert("L")) < 128
    if ink.shape != (4677, 3307):
        sys.exit(f"raster is {ink.shape[1]} x {ink.shape[0]} px; expected 3307 x 4677 (pdfimages, PDF p. {fig['page']})")
    cal = calibrate(ink, fig)
    G = templates(ink, fig, *fig["scale"])
    G1 = templates(ink, fig, 1.0, 1.0)
    data = ink.copy()
    my, mx = fig["legend_mask"]
    data[my:, mx:] = False
    dc = curve_mask(data, fig["curve_anchors"])
    return fig, cal, G, G1, data, dc


def second_read(data, seed):
    """Manual centre (step 7): the centre of the bounding box of the ink inside
    the glyph box read by eye, or None if the seed has no box."""
    if len(seed) < 6 or seed[5] is None:
        return None
    x0, y0, x1, y1 = seed[5]
    ys, xs = np.nonzero(data[y0:y1 + 1, x0:x1 + 1])
    return x0 + (xs.min() + xs.max()) / 2, y0 + (ys.min() + ys.max()) / 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raster", help="PNG of the figure's PDF page from `pdfimages -f P -l P -png`")
    ap.add_argument("--target", choices=sorted(FIGURES), default="Cu", help="which figure (default Cu, PDF p. 118)")
    ap.add_argument("--write", metavar="DIR", help="write one dataset JSON per reference into DIR")
    ap.add_argument("--record", metavar="JSON", help="write the manual-versus-template record (Si, Ag, Au)")
    ap.add_argument("--added", default="2026-10-06", help="date for the datasets' `added` field")
    args = ap.parse_args()

    fig, (cx, cy, rms_x, rms_y, n_ticks), G, G1, data, dc = load(args.raster, args.target)

    def to_log(x, y):
        return cx[0] + cx[1] * x + cx[2] * y, cy[0] + cy[1] * x + cy[2] * y

    print("template check on isolated glyphs (coverage, scaled / legend size): " + ", ".join(
        f"{L} {refine(data, dc, G[L], x, y, 1)[2]:.3f} / {refine(data, dc, G1[L], x, y, 1)[2]:.3f}"
        for L, x, y in fig["isolated_check"]))

    pts = []
    for seed in fig["seeds"]:
        L, x, y, status, note = seed[:5]
        xa, ya, ca = refine(data, dc, G[L], x, y, 1)
        xb, yb, cb = refine(data, dc, G1[L], x, y, 2)
        xm, ym = (xa + xb) / 2, (ya + yb) / 2
        le, ly = to_log(xm, ym)
        # centring repeatability, in decades
        rep_x = abs(to_log(xa, ya)[0] - to_log(xb, yb)[0]) / 2
        rep_y = abs(to_log(xa, ya)[1] - to_log(xb, yb)[1]) / 2
        man = second_read(data, seed)
        pts.append(dict(letter=L, x=xm, y=ym, log_e=le, log_y=ly, cov=(ca, cb), rep=(rep_x, rep_y), status=status,
                        note=note, manual=None if man is None else to_log(*man), box=seed[5] if len(seed) > 5 else None))

    # Symbol-placement scatter: deviation of evidently round-energy points from the nearest round energy.
    rounds = np.log10([200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 1000, 1500, 2000, 2500, 10000])
    dev = []
    for p in pts:
        d = rounds - p["log_e"]
        k = np.argmin(np.abs(d))
        if abs(d[k]) < 0.012:  # within ~3 %: treated as a nominal round energy
            dev.append(-d[k])
    dev = np.array(dev)
    place = float(np.sqrt(np.mean(dev**2)))
    print(f"round-energy check: {len(dev)} points within 3 % of a round energy; mean offset {dev.mean():+.4f} decade "
          f"({100 * (10**dev.mean() - 1):+.1f} %), rms {place:.4f} decade ({100 * (10**place - 1):.1f} %)")

    print("\nletter  x px    y px     E eV     Y      cov(A/B)     unc E  unc Y  status  reference")
    rows = []
    for p in sorted(pts, key=lambda p: (p["letter"], p["log_e"])):
        e, yv = 10 ** p["log_e"], 10 ** p["log_y"]
        ux = math.hypot(math.hypot(rms_x, p["rep"][0]), place)
        uy = math.hypot(math.hypot(rms_y, p["rep"][1]), place)
        ue = max(0.02, math.log(10) * ux)
        uyr = max(0.02, math.log(10) * uy)
        inrange = E_MIN_EV <= e <= E_MAX_EV
        print(f"  {p['letter']}   {p['x']:7.1f} {p['y']:7.1f} {e:8.0f} {yv:6.3f}  {p['cov'][0]:.3f}/{p['cov'][1]:.3f}  "
              f"{100 * ue:4.1f}%  {100 * uyr:4.1f}%  {p['status']:6s}{'' if inrange else ' (out of range)'}  "
              f"{fig['legend_names'][p['letter']]}")
        if inrange:
            rows.append((p, e, yv, ue, uyr))
    print("\nnot stored (in range):")
    for x, y, why in fig["dropped"]:
        le, ly = to_log(x, y)
        print(f"  ({x}, {y}) E ~ {10**le:.0f} eV: {why}")

    # Step 7: manual (box) versus automatic (template) centres of every stored point.
    second = []
    for p, e, yv, ue, uyr in rows:
        if p["manual"] is None:
            continue
        me, my = 10 ** p["manual"][0], 10 ** p["manual"][1]
        comb_e, comb_y = math.sqrt(2) * ue, math.sqrt(2) * uyr
        lr_e, lr_y = abs(math.log(me / e)), abs(math.log(my / yv))
        verdict = ("within the digitizing uncertainty" if lr_e <= ue and lr_y <= uyr else
                   "within the combined uncertainty" if lr_e <= comb_e and lr_y <= comb_y else
                   "beyond the combined uncertainty")
        second.append({"letter": p["letter"], "id": fig["references"][p["letter"]][0], "energy_ev": float(f"{e:.4g}"),
                       "yield": float(f"{yv:.4g}"), "box_px": list(p["box"]), "energy_manual_ev": float(f"{me:.4g}"),
                       "yield_manual": float(f"{my:.4g}"), "energy_ratio": round(me / e, 4),
                       "yield_ratio": round(my / yv, 4), "unc_rel": [round(ue, 3), round(uyr, 3)], "verdict": verdict})
    if second:
        print("\nsecond read, manual box centre / template centre:")
        for s in second:
            print(f"  {s['letter']} {s['energy_ev']:8.0f} eV  Y {s['yield']:6.3f}: E ratio {s['energy_ratio']:.3f}, "
                  f"Y ratio {s['yield_ratio']:.3f}  {s['verdict']}")
        lr = np.array([math.log(s["yield_ratio"]) for s in second])
        le = np.array([math.log(s["energy_ratio"]) for s in second])
        print(f"  {len(second)} points: Y ratio rms {lr.std():.4f} (mean {math.exp(lr.mean()):.4f}), max |ln| "
              f"{np.abs(lr).max():.4f}; E ratio rms {le.std():.4f} (mean {math.exp(le.mean()):.4f}), max |ln| "
              f"{np.abs(le).max():.4f}")
    if args.record:
        rec = {
            "format": "lindhard-digitize-secondread/1",
            "figure": f"IPPJ-AM-32 (1983), Ar -> {args.target}, PDF p. {fig['page']} (matsunami1983_ipp_am32.py)",
            "automatic": "template-coverage centre (mean of the scaled and the legend-size template), as stored",
            "manual": "centre of the bounding box of the ink inside a glyph box read by eye from zoomed rasters "
                      "(no template, same axis calibration)",
            "criterion": "|ln(manual / automatic)| in E and in Y against the stored digitizing uncertainty u "
                         "(within u) and against sqrt(2) u (within the combined uncertainty)",
            "points": second,
        }
        Path(args.record).write_text(json.dumps(rec, indent=2) + "\n")
        print(f"wrote {args.record}")

    if args.write:
        out = Path(args.write)
        out.mkdir(parents=True, exist_ok=True)
        for L in sorted({p["letter"] for p, *_ in rows}):
            ident, ref, doi, nifs = fig["references"][L]
            mine = sorted((r for r in rows if r[0]["letter"] == L), key=lambda r: r[1])
            points = []
            for p, e, yv, ue, uyr in mine:
                pt = {"energy_ev": float(f"{e:.4g}"), "energy_unc_rel": round(ue, 3),
                      "yield": float(f"{yv:.4g}"), "yield_unc_rel": round(uyr, 3)}
                if p["status"] == "infer":
                    pt["note"] = "symbol identity inferred: " + p["note"]
                points.append(pt)
            state = fig["target_state"]
            d = {
                "kind": "sputter_yield",
                "id": ident,
                "ion": "Ar",
                "mass_amu": None,
                "target": args.target,
                "target_state": state[L] if isinstance(state, dict) else state,
                "incidence_deg": 0.0,
                "points": points,
                "original_reference": ref,
                "original_doi": doi,
                "compilation": COMPILATION,
                "compilation_figure": f"Ar -> {args.target}",
                "compilation_pdf_page": fig["page"],
                "compilation_symbol": f"{L} ({fig['legend_names'][L]})",
                "compilation_table_ref": "References for graphs, PDF pp. 281-286 (the figure legend gives the letter)",
                "url": URL,
                "crosscheck": None,  # filled in by yamamura1995_nifs23.py --compare
                "extraction": (
                    f"Digitized from the 400 ppi scan of PDF p. {fig['page']} with validation/data/digitize/"
                    "matsunami1983_ipp_am32.py: 2-D affine axis calibration from all "
                    f"~{round(n_ticks, -1)} tick marks of the log-log "
                    f"frame (residual rms {rms_x:.4f} decade in E, {rms_y:.4f} in Y); symbols identified by eye from "
                    "zoomed rasters (listed in the script) and centred by template matching against the figure's legend "
                    "glyphs, two template sizes; the fitted curve masked out and not stored. Points with 196 eV <= E <= "
                    "10.2 keV only; symbols in that range that could not be identified are listed in the script and not "
                    "stored. Uncertainties are digitizing only (calibration, centring repeatability and the "
                    f"{100 * (10**place - 1):.1f} % rms scatter of round-energy points), floored at 2 %; the compilation "
                    "gives no measurement uncertainty. Energy is the total ion energy, yield in atoms per ion "
                    "(IPPJ-AM-32, p. 6)."
                ),
                "reliability_note": RELIABILITY.format(target=args.target),
                "terms": "Facts, cited; stored under operator decision (#69) despite the compilation's cover note. "
                         "No figure image and no fitted curve stored.",
                "added": args.added,
            }
            path = out / f"{ident}.json"
            if path.exists():  # keep a crosscheck, and any point flag, written earlier by the second read
                prev = json.loads(path.read_text())
                d["crosscheck"] = prev.get("crosscheck")
                flags = {(q["energy_ev"], q["yield"]): q["flag"] for q in prev.get("points", []) if "flag" in q}
                for pt in points:
                    if (pt["energy_ev"], pt["yield"]) in flags:
                        pt["flag"] = flags[pt["energy_ev"], pt["yield"]]
            mine2 = [s for s in second if s["letter"] == L]
            if mine2:
                cc = dict(d["crosscheck"] or {})
                cc["second_read"] = (
                    f"manual glyph-box centre against the template centre (validation/data/digitize/"
                    f"secondread_ar_{args.target.lower()}_am32.json): "
                    + ", ".join(f"Y ratio {s['yield_ratio']:.3f}, E ratio {s['energy_ratio']:.3f} at "
                                f"{s['energy_ev']:.0f} eV ({s['verdict']})" for s in mine2))
                d["crosscheck"] = cc
            path.write_text(json.dumps(d, indent=2) + "\n")
            print(f"wrote {path} ({len(points)} points)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
