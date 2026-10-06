#!/usr/bin/env python3
"""Digitize the measured Ar -> Cu sputtering yields of the compilation
N. Matsunami, Y. Yamamura, Y. Itikawa, N. Itoh, Y. Kazumata, S. Miyagawa,
K. Morita, R. Shimizu and H. Tawara, "Energy dependence of the yields of
ion-induced sputtering of monatomic solids", report IPPJ-AM-32, Institute of
Plasma Physics, Nagoya University (1983); published as At. Data Nucl. Data
Tables 31, 1 (1984), doi:10.1016/0092-640X(84)90016-0.

This is the extraction behind `validation/data/sputtering/ar_cu_sputter_*.json`
(provenance: docs/data-provenance.md, operator decision #69). It is not part
of the harness and needs numpy, scipy and Pillow, which the harness does not.
Neither the report nor any page image is in this tree: fetch the report
yourself and extract the page raster of the Ar -> Cu figure with poppler:

    curl -O http://dpc.nifs.ac.jp/IPPJ-AM/IPPJ-AM-32.pdf
    pdfimages -f 118 -l 118 -png IPPJ-AM-32.pdf am32
    validation/data/digitize/matsunami1983_ipp_am32.py am32-000.png
    validation/data/digitize/matsunami1983_ipp_am32.py am32-000.png --write validation/data/sputtering

The page (PDF p. 118, printed p. 113) is a lossless 1-bit CCITT scan,
3307 x 4677 px at 400 ppi, so every reader gets the same pixels and the
pixel coordinates below are reproducible. The figure plots each reference's
measured points with its own letter (legend in the figure) and the
compilation's empirical fit as a solid line. Only the letters are read; the
fitted line is never digitized or stored (it is masked out, step 3).

Method (also summarized in each dataset's `extraction`):

1. Axes. The frame is slightly rotated in the scan (about 0.17 deg), so the
   pixel -> (log10 E, log10 Y) map is a 2-D affine fit, not two 1-D fits.
   Every tick on all four sides (minor ticks at 2..9 x 10^n, long ticks at
   10^n and 5 x 10^n; about 40 per side) is found as a run of ink
   perpendicular to the fitted frame line, assigned its value from the long
   ticks, and the six affine coefficients are least-squares fits to all of
   them. The script stops if a side has fewer than 35 ticks, and prints the
   residuals (decades).
2. Templates. Each legend letter (the column A..W beside the reference names)
   is cut out of the page. The plotted letters are drawn smaller than the
   legend letters: isolated plotted glyphs are 0.89 times as wide and 0.957
   times as tall (measured on T, V, A, I; the check prints their template
   coverage), so the templates are resampled by that factor.
3. Fitted curve. The compilation's empirical curve crosses many symbols. It
   is located from hand-read anchor points snapped to the thin ink run
   nearest a spline through them, and treated as "don't care" (neither ink
   nor background) when templates are matched. It is not stored.
4. Symbols. Many letters overlap. Each symbol was identified by eye from 4x
   to 6x zooms and ASCII dumps of the raster, helped by a template-matching
   cluster solver; the identifications are the SEEDS list below (letter,
   approximate centre in pixels, status). Identity rules: a glyph is kept
   only where its distinguishing strokes are visible (e.g. the H crossbar
   sits at 31-42 % of the glyph height and the A crossbar at 51-64 %; R and
   B differ only by B's bottom bar; the slashed O; the G inner bar).
   "infer" marks a glyph whose identity rests on an argument stated in the
   seed note; DROPPED lists every symbol in range that was not stored, and
   why. Nothing is guessed.
5. Centres. The plotted symbols are centred on the data point: in x the glyph
   centres of the evidently round-energy points (150, 200, 250, 300, 350 eV
   and so on) fall within about 1 % of the round values, whereas lower-left
   anchoring would put them 6-10 % low. The same centring is assumed in y.
   Each centre is refined twice, by the maximum of template coverage with
   the scaled template (1 px tolerance) and with the unscaled legend glyph
   (2 px tolerance); the value is their mean and half their difference is
   the centring repeatability.
6. Values. E and Y come from the affine map. Only points with
   196 eV <= E <= 10.2 keV are stored (the issue's "about 0.2 to 10 keV",
   with 2 % slack for nominal 200 eV and 10 keV points). Uncertainties are
   digitizing uncertainties only (the compilation gives no measurement
   errors): in each log coordinate, the quadrature sum of the calibration
   residual rms, the centring repeatability and the symbol-placement scatter,
   the latter measured as the rms deviation of the round-energy points from
   their nominal energies; floored at 2 %.

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
LEGEND_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVW"
LEGEND_BOX = (1480, 1600, 1600, 3400)  # x0, x1, y0, y1 of the legend letter column
LEGEND_MASK = (1620, 1495)  # legend block: y >= 1620 and x >= 1495 (no data there in range)
SCALE_X, SCALE_Y = 0.89, 0.957  # plotted / legend glyph size (step 2)
ISOLATED_CHECK = [("T", 975, 1804), ("V", 812, 1800), ("A", 1071, 1755), ("I", 1754, 1433)]

# Step 3: hand-read points on the fitted curve (x, y in px), snapped to ink.
CURVE_ANCHORS = [
    (700, 2026), (780, 1924), (840, 1858), (900, 1798), (960, 1745), (1020, 1697), (1080, 1653),
    (1140, 1612), (1200, 1581), (1320, 1515), (1380, 1486), (1440, 1462), (1560, 1403), (1600, 1396),
    (1660, 1385), (1720, 1370), (1780, 1355), (1830, 1347), (1860, 1344), (1920, 1339), (1980, 1335),
    (2100, 1333),
]

# Step 4: every symbol stored. (letter, x px, y px, status, note)
SEEDS = [
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
DROPPED = [
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

COMPILATION = "Matsunami et al., IPPJ-AM-32 (1983); At. Data Nucl. Data Tables 31, 1 (1984)"
URL = "http://dpc.nifs.ac.jp/IPPJ-AM/IPPJ-AM-32.pdf"
TARGET_STATE = (
    "not stated in the compilation for this set; the compilation excludes single-crystal targets and states "
    "that most specimens were thinned polycrystalline foils or evaporated films (IPPJ-AM-32, p. 2)"
)
RELIABILITY = (
    "Included in the compilation's Ar -> Cu figure, so it passed the compilation's selection (IPPJ-AM-32, pp. 1-2: "
    "absolute or normalized-relative yields, room temperature, no dose effects, sufficient beam current; doubtful "
    "sets excluded). The compilation does not rank the sets it shows."
)
# AM-32 letter -> (dataset id, original reference, DOI or None, NIFS-DATA-23 Fig. 120 letter or None)
REFERENCES = {
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
    "T": ("ar_cu_sputter_akaishi1977", "K. Akaishi, A. Miyahara, Z. Kabeya, S. Skenobu, M. Komizo and T. Gotoh, "
          "J. Vac. Soc. Japan 20, 161 (1977)", "10.3131/jvsj.20.161", "O"),
    "V": ("ar_cu_sputter_bohdansky1980", "J. Bohdansky, J. Nucl. Mater. 93-94, 44 (1980)",
          "10.1016/0022-3115(80)90302-5", "Q"),
    "W": ("ar_cu_sputter_okajima1981", "Y. Okajima, Jpn. J. Appl. Phys. 20, 2313 (1981)", "10.1143/JJAP.20.2313", "V"),
}
LEGEND_NAMES = {
    "A": "KEYWELL (1955)", "B": "GUSEVA (1960)", "D": "YONTS, NORMAND, HARRISON (1960)",
    "F": "BADER, WITTEBORN, SNOUSE (1961)", "G": "FERT, COLOMBIE, FAGOT (1961)", "H": "LAEGREID, WEHNER (1959, 1961)",
    "I": "PEROVIC, COBIC (1961)", "J": "SOUTHERN, WILLIS, ROBINSON (1963)", "M": "DUPP, SCHARMANN (1966)",
    "N": "WEIJSENFELD (1967)", "O": "KOSHKIN, RYSOV, SHKARBAN (1969)", "R": "OECHSNER (1973)",
    "T": "AKAISHI, MIYAHARA, KABEYA (1977)", "V": "BOHDANSKY (1980)", "W": "OKAJIMA (1981)",
}
GRID = np.log10(np.array([k * 10.0**n for n in range(-4, 8) for k in range(1, 10)]))


# --- 1. axes ---

def frame(ink):
    h, w = ink.shape
    rows, cols = ink.sum(1), ink.sum(0)
    ys = [i for i in range(h) if rows[i] > 0.4 * w]
    xs = [j for j in range(w) if cols[j] > 0.4 * h]
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


def assign(ticks, first_long, sign):
    """Value (log10) of each tick: long ticks alternate 10^n, 5 x 10^n from `first_long`."""
    ticks = sorted(ticks, key=lambda t: sign * t[0])
    lmax = max(t[2] for t in ticks)
    longs = [t for t in ticks if t[2] >= 0.8 * lmax]
    vals, v = [], first_long
    for _ in longs:
        vals.append(v)
        v = round(v + (math.log10(5) if abs(v - round(v)) < 1e-9 else 1 - math.log10(5)), 6)
    p = np.polyfit([t[0] for t in longs], vals, 1)
    return [(t[0], t[1], GRID[np.argmin(np.abs(GRID - np.polyval(p, t[0])))]) for t in ticks]


def calibrate(ink):
    t, b, l, r = frame(ink)
    lines = {
        "bottom": (True, side_line(ink, True, b, l + 30, r - 30), -1),
        "top": (True, side_line(ink, True, t, l + 30, r - 30), +1),
        "left": (False, side_line(ink, False, l, t + 30, b - 30), +1),
        "right": (False, side_line(ink, False, r, t + 30, b - 30), -1),
    }
    X, Y = [], []
    for name, (hz, line, inward) in lines.items():
        a, bb = (l + 10, r - 10) if hz else (t + 10, b - 10)
        tk = side_ticks(ink, hz, line, a, bb, inward)
        if len(tk) < 35:
            sys.exit(f"{name} axis: found {len(tk)} ticks, expected about 40; wrong page or raster?")
        if hz:  # x axis: the first long tick from the left is 100 eV
            X += [(al, pp, v) for al, pp, v in assign(tk, 2.0, +1)]
        else:  # y axis: the first long tick from the bottom is 5e-3
            Y += [(pp, al, v) for al, pp, v in assign(tk, math.log10(5e-3), -1)]
    X, Y = np.array(X), np.array(Y)
    A = np.c_[np.ones(len(X)), X[:, 0], X[:, 1]]
    B = np.c_[np.ones(len(Y)), Y[:, 0], Y[:, 1]]
    cx = np.linalg.lstsq(A, X[:, 2], rcond=None)[0]
    cy = np.linalg.lstsq(B, Y[:, 2], rcond=None)[0]
    rx, ry = A @ cx - X[:, 2], B @ cy - Y[:, 2]
    print(f"axes: {len(X)} x ticks, {len(Y)} y ticks; residual rms {rx.std():.4f} / {ry.std():.4f} decade, "
          f"max {np.abs(rx).max():.4f} / {np.abs(ry).max():.4f} decade")
    print(f"  log10 E = {cx[0]:.6f} + {cx[1]:.7e} x + {cx[2]:.3e} y;  log10 Y = {cy[0]:.6f} + {cy[1]:.3e} x + {cy[2]:.7e} y")
    return cx, cy, float(rx.std()), float(ry.std())


# --- 2. templates, 3. curve ---

def templates(ink, sx, sy):
    x0, x1, y0, y1 = LEGEND_BOX
    sub = ink[y0:y1, x0:x1]
    rows, runs, cur = sub.sum(1), [], None
    for i, v in enumerate(rows):
        if v > 0 and cur is None:
            cur = i
        if v == 0 and cur is not None:
            runs.append((cur, i - 1))
            cur = None
    if len(runs) != len(LEGEND_LETTERS):
        sys.exit(f"legend: found {len(runs)} letters, expected {len(LEGEND_LETTERS)}")
    G = {}
    for (a, b), L in zip(runs, LEGEND_LETTERS):
        g = sub[a:b + 1]
        cs = np.nonzero(g.sum(0))[0]
        g = g[:, cs.min():cs.max() + 1]
        if sx != 1.0 or sy != 1.0:
            h, w = g.shape
            img = Image.fromarray((g * 255).astype(np.uint8)).resize((round(w * sx), round(h * sy)), Image.BOX)
            g = np.asarray(img) >= 128
        G[L] = g
    return G


def curve_mask(ink):
    a = np.array(CURVE_ANCHORS, float)
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raster", help="PNG of PDF p. 118 from `pdfimages -f 118 -l 118 -png`")
    ap.add_argument("--write", metavar="DIR", help="write one dataset JSON per reference into DIR")
    ap.add_argument("--added", default="2026-10-06", help="date for the datasets' `added` field")
    args = ap.parse_args()

    ink = np.asarray(Image.open(args.raster).convert("L")) < 128
    if ink.shape != (4677, 3307):
        sys.exit(f"raster is {ink.shape[1]} x {ink.shape[0]} px; expected 3307 x 4677 (pdfimages, PDF p. 118)")
    cx, cy, rms_x, rms_y = calibrate(ink)

    def to_log(x, y):
        return cx[0] + cx[1] * x + cx[2] * y, cy[0] + cy[1] * x + cy[2] * y

    G = templates(ink, SCALE_X, SCALE_Y)
    G1 = templates(ink, 1.0, 1.0)
    data = ink.copy()
    data[LEGEND_MASK[0]:, LEGEND_MASK[1]:] = False
    dc = curve_mask(data)
    print("template check on isolated glyphs (coverage, scaled / legend size): " + ", ".join(
        f"{L} {refine(data, dc, G[L], x, y, 1)[2]:.3f} / {refine(data, dc, G1[L], x, y, 1)[2]:.3f}"
        for L, x, y in ISOLATED_CHECK))

    pts = []
    for L, x, y, status, note in SEEDS:
        xa, ya, ca = refine(data, dc, G[L], x, y, 1)
        xb, yb, cb = refine(data, dc, G1[L], x, y, 2)
        xm, ym = (xa + xb) / 2, (ya + yb) / 2
        le, ly = to_log(xm, ym)
        # centring repeatability, in decades
        rep_x = abs(to_log(xa, ya)[0] - to_log(xb, yb)[0]) / 2
        rep_y = abs(to_log(xa, ya)[1] - to_log(xb, yb)[1]) / 2
        pts.append(dict(letter=L, x=xm, y=ym, log_e=le, log_y=ly, cov=(ca, cb), rep=(rep_x, rep_y), status=status, note=note))

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
              f"{LEGEND_NAMES[p['letter']]}")
        if inrange:
            rows.append((p, e, yv, ue, uyr))
    print("\nnot stored (in range):")
    for x, y, why in DROPPED:
        le, ly = to_log(x, y)
        print(f"  ({x}, {y}) E ~ {10**le:.0f} eV: {why}")

    if args.write:
        out = Path(args.write)
        out.mkdir(parents=True, exist_ok=True)
        for L in sorted({p["letter"] for p, *_ in rows}):
            ident, ref, doi, nifs = REFERENCES[L]
            mine = sorted((r for r in rows if r[0]["letter"] == L), key=lambda r: r[1])
            points = []
            for p, e, yv, ue, uyr in mine:
                pt = {"energy_ev": float(f"{e:.4g}"), "energy_unc_rel": round(ue, 3),
                      "yield": float(f"{yv:.4g}"), "yield_unc_rel": round(uyr, 3)}
                if p["status"] == "infer":
                    pt["note"] = "symbol identity inferred: " + p["note"]
                points.append(pt)
            d = {
                "kind": "sputter_yield",
                "id": ident,
                "ion": "Ar",
                "mass_amu": None,
                "target": "Cu",
                "target_state": TARGET_STATE,
                "incidence_deg": 0.0,
                "points": points,
                "original_reference": ref,
                "original_doi": doi,
                "compilation": COMPILATION,
                "compilation_figure": "Ar -> Cu",
                "compilation_pdf_page": 118,
                "compilation_symbol": f"{L} ({LEGEND_NAMES[L]})",
                "compilation_table_ref": "References for graphs, PDF pp. 281-286 (the figure legend gives the letter)",
                "url": URL,
                "crosscheck": None,  # filled in by yamamura1995_nifs23.py --compare
                "extraction": (
                    "Digitized from the 400 ppi scan of PDF p. 118 with validation/data/digitize/"
                    "matsunami1983_ipp_am32.py: 2-D affine axis calibration from all ~160 tick marks of the log-log "
                    f"frame (residual rms {rms_x:.4f} decade in E, {rms_y:.4f} in Y); symbols identified by eye from "
                    "zoomed rasters (listed in the script) and centred by template matching against the figure's legend "
                    "glyphs, two template sizes; the fitted curve masked out and not stored. Points with 196 eV <= E <= "
                    "10.2 keV only; symbols in that range that could not be identified are listed in the script and not "
                    "stored. Uncertainties are digitizing only (calibration, centring repeatability and the "
                    f"{100 * (10**place - 1):.1f} % rms scatter of round-energy points), floored at 2 %; the compilation "
                    "gives no measurement uncertainty. Energy is the total ion energy, yield in atoms per ion "
                    "(IPPJ-AM-32, p. 6)."
                ),
                "reliability_note": RELIABILITY,
                "terms": "Facts, cited; stored under operator decision (#69) despite the compilation's cover note. "
                         "No figure image and no fitted curve stored.",
                "added": args.added,
            }
            path = out / f"{ident}.json"
            if path.exists():  # keep a crosscheck written earlier by the second read
                d["crosscheck"] = json.loads(path.read_text()).get("crosscheck")
            path.write_text(json.dumps(d, indent=2) + "\n")
            print(f"wrote {path} ({len(points)} points)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
