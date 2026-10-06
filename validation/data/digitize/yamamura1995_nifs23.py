#!/usr/bin/env python3
"""Second, independent read of the measured Ar -> Cu, Si, Ag and Au
sputtering yields, from Y. Yamamura and H. Tawara, "Energy dependence of
ion-induced sputtering yields from monatomic solids at normal incidence",
report NIFS-DATA-23, National Institute for Fusion Science (1995); published
as At. Data Nucl. Data Tables 62, 149 (1996), doi:10.1006/adnd.1996.0005.
Figures: Ar -> Cu Fig. 120, PDF p. 49 (upper right; #69); Ar -> Si Fig. 41,
PDF p. 29 (lower left), Ar -> Ag Fig. 209, PDF p. 71 (lower left), Ar -> Au
Fig. 310, PDF p. 96 (lower right; #70). The pages were located by reading
the figure titles.

This is the double-check of `matsunami1983_ipp_am32.py` (issues #69, #70). It
shares no code and no constants with that script: its own axis calibration
(a log-comb fit, not long-tick assignment; per-side 1-D fits blended across
the frame, not a 2-D affine fit), its own template-free centring (ink bounding
box and ink centroid of a hand-drawn box; no legend templates, no
cross-correlation). Everything figure-specific (panel, frame threshold,
letters, symbols, stacks) is in that figure's table (CU, SI, AG, AU, collected
in FIGURES). Neither the report nor any page image is in this tree; fetch it
yourself:

    curl -O https://www.nifs.ac.jp/report/NIFS-DATA-023.pdf
    pdfimages -f 49 -l 49 -png NIFS-DATA-023.pdf nifs
    validation/data/digitize/yamamura1995_nifs23.py nifs-000.png
    validation/data/digitize/yamamura1995_nifs23.py nifs-000.png --compare validation/data/sputtering \
        --write validation/data/digitize/crosscheck_ar_cu_nifs23.json
    pdfimages -f 96 -l 96 -png NIFS-DATA-023.pdf nifsau
    validation/data/digitize/yamamura1995_nifs23.py nifsau-000.png --target Au \
        --compare validation/data/sputtering --write validation/data/digitize/crosscheck_ar_au_nifs23.json

Each page is a lossless 1-bit scan at 200 ppi (1603, 1596 or 1594 x 2314 px,
checked per figure). The symbols are about 9 x 12 px, a quarter of the linear
size of IPPJ-AM-32's, so many overlap beyond reading.

Method:

1. Axes. Each frame side's inner edge is fitted with a straight line; every
   tick (a 3-20 px run of ink beyond that edge) is located; an ideal
   logarithmic comb (log10 of k x 10^n, k = 1..9) is fitted to the tick
   positions by a grid search over offset and scale, then by least squares
   on the matched ticks. The scans are sheared or rotated (on p. 49 the left
   edge leans about 0.5 deg, the bottom edge 0.02 deg), so log10 E is the
   bottom-axis fit and the top-axis fit blended linearly in y, and log10 Y the
   left and right fits blended linearly in x. The right side of Fig. 310 has
   no ticks in the scan; there the right-side fit is the left-side fit moved
   along the mean slope of the fitted top and bottom edges (the panel is
   rotated by about 0.6 deg). Using the left fit alone instead read all nine
   compared Au points 4-9 % high, which is how the rotation was noticed; the
   correction uses only the frame geometry. Residuals are printed.
2. Symbols. Each figure uses its own letters (legend in each figure table's
   comment). Each symbol in a figure's seeds was read by eye from 6x to 10x
   zooms of the raster. The search was guided by where IPPJ-AM-32 puts the
   same reference's point (stated, not hidden), but each position below is
   read and then measured in this raster only. Legibility classes:
   "isolated" (the glyph stands alone or its outline is unambiguous),
   "coincident" (two references' symbols overlap at one spot in both
   compilations; compared to the common centre), "cluster" (inside a stack
   where no single glyph can be separated at 200 ppi; reported, not used for
   the agreement statistics). Points labelled ACAT (calculated, the star
   symbol) are never read; sets that IPPJ-AM-32 does not show are not read.
3. Centres. Two template-free estimates per symbol, from the ink inside its
   hand-drawn box in the seeds (`centres()`): the centre of the ink's
   bounding box, and the centroid of the ink pixels. The value is their mean,
   half their difference the centring repeatability.
4. Uncertainty (digitizing only): calibration residual rms, centring
   repeatability, and one pixel of symbol placement (about 0.0085 decade in E
   and 0.0076 in Y), in quadrature.

`--compare DIR` pairs every symbol with the stored IPPJ-AM-32 points of the
same reference (nearest within 0.05 decade in E and 0.1 in Y) and reports the
ratio of the two reads against their combined uncertainty; `--write` keeps
that record as JSON and writes each dataset's `crosscheck` summary (keeping a
`second_read` entry written by matsunami1983_ipp_am32.py). A read beyond the
combined uncertainty is called explained by coincident symbols only by the
rule in `coinc()`: the NIFS read must lie within the range of the IPPJ-AM-32
reads of the references drawn on that spot. A read beyond twice the combined
uncertainty that the rule does not explain gets
`"flag": "disagrees_between_compilations"` on its point in the dataset (none
does at present).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image

# Fig. 120 (Ar -> Cu), PDF p. 49, upper right. Legend: A Keywell, B Guseva, D Yonts, F Fert, G Laegreid,
# H Perovic, J Southern, L Dupp, N Koshkin, O Akaishi, P Weijsenfeld 1967, Q Bohdansky, U Oechsner,
# V Okajima, Z Bader, ...; no point of these references is marked ACAT; the "#" set (Eckstein 1983) is not
# in IPPJ-AM-32 and is not read.
# Dataset id (validation/data/sputtering) -> this figure's letter for the same reference.
CU_LETTER_OF = {
    "ar_cu_sputter_keywell1955": "A", "ar_cu_sputter_guseva1960": "B", "ar_cu_sputter_yonts1960": "D",
    "ar_cu_sputter_fert1961": "F", "ar_cu_sputter_laegreid1961": "G", "ar_cu_sputter_perovic1961": "H",
    "ar_cu_sputter_southern1963": "J", "ar_cu_sputter_dupp1966": "L", "ar_cu_sputter_koshkin1969": "N",
    "ar_cu_sputter_akaishi1977": "O", "ar_cu_sputter_weijsenfeld1967": "P", "ar_cu_sputter_bohdansky1980": "Q",
    "ar_cu_sputter_oechsner1973": "U", "ar_cu_sputter_okajima1981": "V", "ar_cu_sputter_bader1961": "Z",
}

# Glyphs read from this raster: (letters, legibility, box x0, y0, x1, y1 in px, note). The box encloses
# the glyph's ink and excludes the fitted curve and neighbouring glyphs.
CU_SEEDS = [
    ("GO", "coincident", 1092, 411, 1101, 425, "about 200 eV: G and O drawn on one spot (also coincident in IPPJ-AM-32)"),
    ("O", "isolated", 1103, 407, 1113, 420, "about 250 eV"),
    ("O", "isolated", 1111, 418, 1123, 433, "about 300 eV, below the curve; its top bar touches the curve"),
    ("A", "isolated", 1134, 407, 1143, 421, "about 450 eV, below the curve"),
    ("ON", "coincident", 1172, 376, 1183, 389, "about 1 keV: O and N on one spot (also coincident in IPPJ-AM-32)"),
    ("B", "isolated", 1255, 338, 1266, 351, "about 5 keV, below the curve"),
    ("DN", "coincident", 1277, 303, 1287, 317, "7.5-8 keV: D and N overlap (also in IPPJ-AM-32)"),
]
# Stacks where no single glyph can be separated at 200 ppi (x0, y0, x1, y1 px).
CU_CLUSTERS = [
    (1100, 380, 1139, 418, "350-450 eV stack"),
    (1138, 355, 1170, 395, "500-800 eV stack"),
    (1165, 340, 1200, 378, "0.8-1.5 keV stack"),
    (1193, 325, 1235, 372, "1.5-2.7 keV stack"),
    (1240, 300, 1300, 345, "3-10 keV stack"),
]
CU = {"figure": "Fig. 120", "page": 49, "am32_page": 118, "shape": (2314, 1603), "panel": (900, 1560, 140, 840),
      "frame_frac": 0.5, "letter_of": CU_LETTER_OF, "seeds": CU_SEEDS, "clusters": CU_CLUSTERS}

# Fig. 41 (Ar -> Si), PDF p. 29, lower left. Legend: A Laegreid, B Wehner, C Eernisse, D Sommerfeldt,
# E Andersen, F Coburn, G Kang, H Tung-Ti Tu, I Morgan, J Oostra, K Poate; the star is ACAT (calculated)
# and is not read.
SI_LETTER_OF = {"ar_si_sputter_laegreid1961": "A", "ar_si_sputter_coburn1977": "F", "ar_si_sputter_kang1979": "G",
                "ar_si_sputter_poate1976": "K"}
SI_SEEDS = [
    ("F", "isolated", 438, 1451, 448, 1466, "about 1.5 keV; the curve crosses its top bar"),
    ("G", "isolated", 531, 1419, 545, 1430, "about 10 keV, right of an I; its top lies under the curve, so the box "
     "starts below the curve"),
]
SI_CLUSTERS = [
    (320, 1440, 432, 1550, "0.2-1 keV stack of Wehner B, Laegreid A, Coburn F, Poate K and ACAT stars"),
    (432, 1395, 530, 1456, "1.5-5 keV stack of Coburn F, Kang G, Wehner B, Tung-Ti Tu H and Morgan I"),
]
SI = {"figure": "Fig. 41", "page": 29, "am32_page": 47, "shape": (2314, 1603), "panel": (150, 800, 1150, 1865),
      "frame_frac": 0.5, "letter_of": SI_LETTER_OF, "seeds": SI_SEEDS, "clusters": SI_CLUSTERS}

# Fig. 209 (Ar -> Ag), PDF p. 71, lower left. Legend: A Keywell, B Guseva, C Almen, D Fert, E Laegreid,
# F Perovic, G Wehner, H Ramer, I Smith, J Brawn (1958), K Andersen, L Oechsner, M Okajima, N Braun,
# O Benninghoven; star ACAT. Koedam (IPPJ-AM-32's B) has no symbol here.
AG_LETTER_OF = {"ar_ag_sputter_keywell1955": "A", "ar_ag_sputter_guseva1960": "B", "ar_ag_sputter_laegreid1961": "E",
                "ar_ag_sputter_wehner1961": "G", "ar_ag_sputter_smith1975": "I", "ar_ag_sputter_oechsner1973": "L",
                "ar_ag_sputter_okajima1981": "M", "ar_ag_sputter_benninghoven1969": "O"}
AG_SEEDS = [
    ("A", "isolated", 393, 1342, 406, 1357, "about 650 eV; the curve runs under its feet"),
    ("O", "isolated", 416, 1322, 428, 1339, "about 1 keV, above the curve"),
    ("A", "isolated", 432, 1326, 447, 1341, "about 1.4 keV; the box stops at the curve"),
    ("B", "isolated", 471, 1316, 481, 1331, "about 3 keV; its lower right touches the next glyph"),
]
AG_CLUSTERS = [
    (300, 1335, 386, 1395, "0.2-0.6 keV stack of Laegreid E, Wehner G, Smith I and ACAT stars"),
    (455, 1280, 545, 1340, "2.5-10 keV stack of Keywell A, Guseva B, Okajima M and others"),
]
AG = {"figure": "Fig. 209", "page": 71, "am32_page": 191, "shape": (2314, 1596), "panel": (145, 795, 1150, 1865),
      "frame_frac": 0.3, "letter_of": AG_LETTER_OF, "seeds": AG_SEEDS, "clusters": AG_CLUSTERS}

# Fig. 310 (Ar -> Au), PDF p. 96, lower right. Legend: A Almen, B Laegreid, C Wehner, D Patterson,
# E Colombie, F Weijsenfeld (1966), G Nenadovic, H Andersen, I Wittmaack, J Eernisse, K Colligon,
# L Sletten, M Robinson, N Weijsenfeld (1967), O Fitch, P Fitch, Q Chenecrojian, R Brauer, S Oliva-Florio,
# T Yamashita (1980), U Oechsner, V Holloway, W Ato, X Braun, Y Benninghoven, Z Szymonski, and a filled
# triangle (written "^" here) for Yamashita (1982); no ACAT points. The right axis of this scan has no
# ticks (right_ticks False).
AU_LETTER_OF = {"ar_au_sputter_almen1961": "A", "ar_au_sputter_laegreid1961": "B", "ar_au_sputter_patterson1962": "D",
                "ar_au_sputter_robinson1967": "M", "ar_au_sputter_weijsenfeld1967": "N",
                "ar_au_sputter_benninghoven1969": "Y", "ar_au_sputter_sletten1972": "L", "ar_au_sputter_oechsner1973": "U",
                "ar_au_sputter_holloway1977": "V", "ar_au_sputter_szymonski1978": "Z", "ar_au_sputter_yamashita1980": "T",
                "ar_au_sputter_yamashita1982": "^"}
AU_SEEDS = [
    ("T", "isolated", 1091, 1460, 1101, 1474, "about 250 eV"),
    ("T", "isolated", 1109, 1467, 1118, 1481, "about 350 eV"),
    ("T", "isolated", 1128, 1428, 1137, 1443, "about 500 eV"),
    ("D", "isolated", 1161, 1436, 1170, 1450, "about 1 keV"),
    ("Y", "isolated", 1161, 1347, 1170, 1360, "about 1 keV, above another glyph"),
    ("V", "isolated", 1195, 1321, 1206, 1335, "about 2 keV"),
    ("D", "isolated", 1195, 1410, 1206, 1424, "about 2 keV"),
    ("Z", "isolated", 1251, 1287, 1262, 1301, "about 6 keV"),
    ("L", "isolated", 1276, 1353, 1285, 1368, "about 10 keV"),
]
AU_CLUSTERS = [
    (1050, 1378, 1145, 1445, "0.2-0.8 keV stack of Laegreid B and Weijsenfeld N"),
    (1140, 1370, 1215, 1432, "0.75-2 keV stack of Yamashita T and triangle, Weijsenfeld N and Oechsner U"),
    (1170, 1270, 1300, 1370, "2-10 keV stack of Robinson M, Patterson D, Almen A and others"),
]
AU = {"figure": "Fig. 310", "page": 96, "am32_page": 262, "shape": (2314, 1594), "panel": (890, 1550, 1160, 1875),
      "frame_frac": 0.3, "right_ticks": False, "letter_of": AU_LETTER_OF, "seeds": AU_SEEDS, "clusters": AU_CLUSTERS}

FIGURES = {"Cu": CU, "Si": SI, "Ag": AG, "Au": AU}
COMB = np.log10(np.array([k * 10.0**n for n in range(-4, 8) for k in range(1, 10)]))


# --- 1. axes: log-comb fit ---

class Axes:
    def __init__(self, ink, fig):
        x0, x1, y0, y1 = fig["panel"]
        f = fig["frame_frac"]
        sub = ink[y0:y1, x0:x1]
        rows, cols = sub.sum(1), sub.sum(0)
        ys = [i + y0 for i in range(len(rows)) if rows[i] > f * (x1 - x0)]
        xs = [j + x0 for j in range(len(cols)) if cols[j] > f * (y1 - y0)]
        ym, xm = (y0 + y1) / 2, (x0 + x1) / 2
        self.top, self.bottom = float(np.mean([y for y in ys if y < ym])), float(np.mean([y for y in ys if y > ym]))
        self.left, self.right = float(np.mean([x for x in xs if x < xm])), float(np.mean([x for x in xs if x > xm]))
        self.ink = ink
        span_x, span_y = self.right - self.left, self.bottom - self.top
        res = []
        self.edge_slope = {}
        # x: frame spans 10^1 .. 10^6 eV; y: 10^-3 .. 10^2 atoms/ion
        self.fx_bottom, r = self._side(True, self.bottom, -1, 5 / span_x, 1.0, self.left)
        res.append(("bottom", r))
        self.fx_top, r = self._side(True, self.top, +1, 5 / span_x, 1.0, self.left)
        res.append(("top", r))
        self.fy_left, r = self._side(False, self.left, +1, -5 / span_y, -3.0, self.bottom)
        res.append(("left", r))
        if fig.get("right_ticks", True):
            self.fy_right, r = self._side(False, self.right, -1, -5 / span_y, -3.0, self.bottom)
            res.append(("right", r))
        else:
            # No ticks on the right side of this scan: the right-side log10 Y fit is the
            # left-side fit moved along the mean slope of the fitted top and bottom edges.
            tilt = (self.edge_slope[True, -1] + self.edge_slope[True, +1]) / 2
            self.fy_right = (self.fy_left[0] - self.fy_left[1] * tilt * span_x, self.fy_left[1])
            print(f"axis right: no ticks; left-side fit moved along the frame tilt ({tilt:+.4f} px/px)")
        allr = np.concatenate([r for _, r in res])
        self.rms = float(allr.std())
        for name, r in res:
            print(f"axis {name}: {len(r)} ticks on the comb, residual rms {r.std():.4f}, max {np.abs(r).max():.4f} decade")
        if min(len(r) for _, r in res) < 12:
            sys.exit("too few ticks on a side; wrong page or raster?")

    def _side(self, horizontal, approx, inward, b0, v_end, t_end):
        ink = self.ink
        lo, hi = (int(self.left) + 3, int(self.right) - 2) if horizontal else (int(self.top) + 3, int(self.bottom) - 2)
        tt, ee = [], []
        for t in range(lo, hi):
            line = [ink[p, t] if horizontal else ink[t, p] for p in range(int(approx) - 4, int(approx) + 5)]
            idx = [i for i, v in enumerate(line) if v]
            if idx:
                tt.append(t)
                ee.append(int(approx) - 4 + (max(idx) if inward > 0 else min(idx)))
        tt, ee = np.array(tt), np.array(ee)
        p = np.polyfit(tt, ee, 1)
        for _ in range(3):
            k = np.abs(ee - np.polyval(p, tt)) <= 1.5
            p = np.polyfit(tt[k], ee[k], 1)
        self.edge_slope[horizontal, inward] = float(p[0])
        lens = []
        for t in range(lo, hi):
            e, n = np.polyval(p, t), 0
            while n < 25:
                q = int(round(e + inward * (1 + n)))
                if not (ink[q, t] if horizontal else ink[t, q]):
                    break
                n += 1
            lens.append(n)
        lens = np.array(lens)
        ticks, i = [], 0
        while i < len(lens):
            if 3 <= lens[i] <= 20:
                j = i
                while j + 1 < len(lens) and 3 <= lens[j + 1] <= 20:
                    j += 1
                if j - i + 1 <= 4:
                    ticks.append(lo + (i + j) / 2)
                i = j + 1
            else:
                i += 1
        ts = np.array(ticks)
        best = None
        for b in np.linspace(b0 * 0.95, b0 * 1.05, 401):
            for da in np.linspace(-0.06, 0.06, 121):
                v = v_end + da + b * (ts - t_end)
                s = np.exp(-((np.abs(v[:, None] - COMB[None, :]).min(1) / 0.006) ** 2)).sum()
                if best is None or s > best[0]:
                    best = (s, v_end + da - b * t_end, b)
        _, a, b = best
        v = a + b * ts
        tgt = COMB[np.abs(v[:, None] - COMB[None, :]).argmin(1)]
        m = np.abs(v - tgt) < 0.015
        b2, a2 = np.polyfit(ts[m], tgt[m], 1)
        return (a2, b2), a2 + b2 * ts[m] - tgt[m]

    def to_pixel(self, le, ly):
        x, y = (self.left + self.right) / 2, (self.top + self.bottom) / 2
        for _ in range(30):
            a, b = self.to_log(x, y)
            x += (le - a) / self.fx_bottom[1]
            y += (ly - b) / self.fy_left[1]
        return x, y

    def to_log(self, x, y):
        w = (y - self.top) / (self.bottom - self.top)
        le = w * (self.fx_bottom[0] + self.fx_bottom[1] * x) + (1 - w) * (self.fx_top[0] + self.fx_top[1] * x)
        u = (x - self.left) / (self.right - self.left)
        ly = (1 - u) * (self.fy_left[0] + self.fy_left[1] * y) + u * (self.fy_right[0] + self.fy_right[1] * y)
        return le, ly


# --- 2./3. symbols ---

def centres(ink, x0, y0, x1, y1):
    """Two template-free centres of the ink in a box: bounding-box centre and ink centroid."""
    ys, xs = np.nonzero(ink[y0:y1 + 1, x0:x1 + 1])
    bb = (x0 + (xs.min() + xs.max()) / 2, y0 + (ys.min() + ys.max()) / 2)
    cg = (x0 + xs.mean(), y0 + ys.mean())
    return bb, cg


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raster", help="PNG of the figure's PDF page from `pdfimages -f P -l P -png`")
    ap.add_argument("--target", choices=sorted(FIGURES), default="Cu", help="which figure (default Cu, Fig. 120)")
    ap.add_argument("--compare", metavar="DIR", help="pair with the stored IPPJ-AM-32 datasets in DIR")
    ap.add_argument("--write", metavar="JSON", help="write the comparison record (needs --compare)")
    args = ap.parse_args()
    fig = FIGURES[args.target]
    t = args.target.lower()

    ink = np.asarray(Image.open(args.raster).convert("L")) < 128
    if ink.shape != fig["shape"]:
        sys.exit(f"raster is {ink.shape[1]} x {ink.shape[0]} px; expected {fig['shape'][1]} x {fig['shape'][0]} "
                 f"(pdfimages, PDF p. {fig['page']})")
    ax = Axes(ink, fig)
    px_e = abs(ax.fx_bottom[1])
    px_y = abs(ax.fy_left[1])

    reads = []
    print("\nletters  centre px        E eV     Y     unc E  unc Y  legibility")
    for letters, leg, x0, y0, x1, y1, note in fig["seeds"]:
        (xa, ya), (xb, yb) = centres(ink, x0, y0, x1, y1)
        le, ly = ax.to_log((xa + xb) / 2, (ya + yb) / 2)
        la, lya = ax.to_log(xa, ya)
        lb, lyb = ax.to_log(xb, yb)
        ue = math.log(10) * math.sqrt(ax.rms**2 + ((la - lb) / 2) ** 2 + px_e**2)
        uy = math.log(10) * math.sqrt(ax.rms**2 + ((lya - lyb) / 2) ** 2 + px_y**2)
        r = dict(letters=letters, x=(xa + xb) / 2, y=(ya + yb) / 2, energy_ev=10**le, yield_=10**ly,
                 energy_unc_rel=ue, yield_unc_rel=uy, legibility=leg, note=note)
        reads.append(r)
        print(f"  {letters:3s}  ({r['x']:6.1f}, {r['y']:5.1f})  {10**le:8.0f} {10**ly:6.3f}  {100 * ue:4.1f}%  "
              f"{100 * uy:4.1f}%  {leg}: {note}")

    if not args.compare:
        return 0
    letter_of = fig["letter_of"]
    rows = []
    for f in sorted(Path(args.compare).glob(f"ar_{t}_sputter_*.json")):
        d = json.loads(f.read_text())
        L = letter_of.get(d["id"])
        for p in d["points"]:
            cands = [r for r in reads if L is not None and L in r["letters"]
                     and abs(math.log10(r["energy_ev"] / p["energy_ev"])) < 0.05
                     and abs(math.log10(r["yield_"] / p["yield"])) < 0.1]
            row = {"id": d["id"], "nifs_letter": L, "energy_ev": p["energy_ev"], "yield_am32": p["yield"]}
            if L is None:
                row.update(verdict=f"reference not in NIFS-DATA-23 {fig['figure']}; not compared")
            elif not cands:
                x, y = ax.to_pixel(math.log10(p["energy_ev"]), math.log10(p["yield"]))
                inside = [c[4] for c in fig["clusters"] if c[0] <= x <= c[2] and c[1] <= y <= c[3]]
                row.update(verdict=("inside an unreadable NIFS stack (" + inside[0] + "); not compared") if inside
                           else f"no separable symbol of this reference found in NIFS-DATA-23 {fig['figure']}")
            else:
                r = min(cands, key=lambda r: abs(math.log(r["energy_ev"] / p["energy_ev"])))
                ratio = r["yield_"] / p["yield"]
                comb = math.hypot(p["yield_unc_rel"], r["yield_unc_rel"])
                row.update(energy_nifs_ev=round(r["energy_ev"], 1), energy_ratio=round(r["energy_ev"] / p["energy_ev"], 4),
                           yield_nifs=round(r["yield_"], 4), yield_ratio=round(ratio, 4),
                           combined_unc_rel=round(comb, 4), legibility=r["legibility"])
                if abs(math.log(ratio)) <= comb:
                    row["verdict"] = "agree within the combined uncertainty"
                elif abs(math.log(ratio)) <= 2 * comb:
                    row["verdict"] = "agree within twice the combined uncertainty"
                else:
                    row["verdict"] = "disagree"
            rows.append(row)
    print("\nid                               E eV    Y AM-32  Y NIFS   ratio  comb.unc  verdict")
    for row in rows:
        print(f"  {row['id']:30s} {row['energy_ev']:7.0f}  {row['yield_am32']:6.3f}  "
              + (f"{row['yield_nifs']:6.3f}  {row['yield_ratio']:6.3f}  {100 * row['combined_unc_rel']:5.1f}%  "
                 if "yield_nifs" in row else " " * 32)
              + row["verdict"])
    used = [r for r in rows if r.get("legibility") in ("isolated", "coincident")]
    if used:
        lr = np.array([math.log(r["yield_ratio"]) for r in used])
        le = np.array([math.log(r["energy_ratio"]) for r in used])
        print(f"\ncompared {len(used)} points: Y ratio NIFS/AM-32 mean {math.exp(lr.mean()):.3f}, "
              f"rms of ln ratio {lr.std():.3f}; E ratio mean {math.exp(le.mean()):.3f}, rms {le.std():.3f}")

    def coinc(r):
        if r["legibility"] != "coincident":
            return ""
        # A coincident blob explains a split only if its read lies within the
        # range of the AM-32 reads of the references drawn on that spot.
        ys = [q["yield_am32"] for q in rows if q.get("yield_nifs") == r["yield_nifs"]]
        if r["verdict"].startswith("agree within the") or min(ys) <= r["yield_nifs"] <= max(ys):
            return ", symbols coincident in both figures"
        return (", symbols coincident in both figures, but the NIFS read lies outside the AM-32 reads "
                "of the coincident references, so coincidence does not explain it")

    if args.write:
        rec = {
            "format": "lindhard-digitize-crosscheck/1",
            "figure_a": f"IPPJ-AM-32 (1983), Ar -> {args.target}, PDF p. {fig['am32_page']} (matsunami1983_ipp_am32.py)",
            "figure_b": f"NIFS-DATA-23 (1995), {fig['figure']}, PDF p. {fig['page']} (yamamura1995_nifs23.py)",
            "calibration_b_rms_decade": round(ax.rms, 4),
            "pairing": "same reference, nearest NIFS symbol within 0.05 decade in E and 0.1 in Y",
            "criterion": "|ln(Y_b / Y_a)| against the combined digitizing uncertainty (quadrature)",
            "nifs_reads": [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in reads],
            "pairs": rows,
        }
        Path(args.write).write_text(json.dumps(rec, indent=2) + "\n")
        print(f"wrote {args.write}")
        # Summarize per dataset into its `crosscheck` field.
        for f in sorted(Path(args.compare).glob(f"ar_{t}_sputter_*.json")):
            d = json.loads(f.read_text())
            mine = [r for r in rows if r["id"] == d["id"]]
            comp = [r for r in mine if "yield_ratio" in r]
            L = letter_of.get(d["id"])
            if comp:
                agreement = (f"{len(comp)} of {len(mine)} points compared; Y ratio NIFS/AM-32 "
                             + ", ".join(f"{r['yield_ratio']:.3f} at {r['energy_ev']:.0f} eV ({r['verdict']}"
                                         + coinc(r)
                                         + ")" for r in comp)
                             + "; the other points lie in stacks not separable at 200 ppi")
            elif L is None:
                agreement = f"this reference has no symbol in {fig['figure']}; not compared"
            else:
                agreement = (f"present in {fig['figure']} (letter {L}), but none of its {len(mine)} points "
                             "is separable from neighbouring symbols at 200 ppi; not compared")
            # A read beyond twice the combined uncertainty that coincidence does not explain is flagged;
            # validation/experiments/run.py leaves flagged points out of its statistics.
            bad = {r["energy_ev"] for r in comp
                   if r["verdict"] == "disagree" and coinc(r) != ", symbols coincident in both figures"}
            for pt in d["points"]:
                if pt["energy_ev"] in bad:
                    pt["flag"] = "disagrees_between_compilations"
                elif pt.get("flag") == "disagrees_between_compilations":
                    del pt["flag"]
            old = d.get("crosscheck") or {}
            d["crosscheck"] = {
                "compilation": "Yamamura and Tawara, NIFS-DATA-23 (1995); At. Data Nucl. Data Tables 62, 149 (1996)",
                "figure": fig["figure"],
                "pdf_page": fig["page"],
                "symbol": L,
                "agreement": agreement + f" (validation/data/digitize/crosscheck_ar_{t}_nifs23.json)",
            }
            if "second_read" in old:  # written by matsunami1983_ipp_am32.py (Si, Ag, Au)
                d["crosscheck"]["second_read"] = old["second_read"]
            f.write_text(json.dumps(d, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
