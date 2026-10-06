#!/usr/bin/env python3
"""Second, independent read of the measured Ar -> Cu sputtering yields, from
Y. Yamamura and H. Tawara, "Energy dependence of ion-induced sputtering
yields from monatomic solids at normal incidence", report NIFS-DATA-23,
National Institute for Fusion Science (1995); published as At. Data Nucl.
Data Tables 62, 149 (1996), doi:10.1006/adnd.1996.0005. Fig. 120, PDF p. 49.

This is the double-check of `matsunami1983_ipp_am32.py` (issue #69). It shares
no code and no constants with that script: its own axis calibration (a
log-comb fit, not long-tick assignment; per-side 1-D fits blended across the
frame, not a 2-D affine fit), its own template-free centring (ink bounding box
and ink centroid of a hand-drawn box; no legend templates, no cross-correlation). Neither the
report nor any page image is in this tree; fetch it yourself:

    curl -O https://www.nifs.ac.jp/report/NIFS-DATA-023.pdf
    pdfimages -f 49 -l 49 -png NIFS-DATA-023.pdf nifs
    validation/data/digitize/yamamura1995_nifs23.py nifs-000.png
    validation/data/digitize/yamamura1995_nifs23.py nifs-000.png --compare validation/data/sputtering \
        --write validation/data/digitize/crosscheck_ar_cu_nifs23.json

The page is a lossless 1-bit scan, 1603 x 2314 px at 200 ppi; Fig. 120 is the
upper-right panel. Its symbols are about 9 x 12 px, a quarter of the linear
size of IPPJ-AM-32's, so many overlap beyond reading.

Method:

1. Axes. Each frame side's inner edge is fitted with a straight line; every
   tick (a 3-20 px run of ink beyond that edge) is located; an ideal
   logarithmic comb (log10 of k x 10^n, k = 1..9) is fitted to the tick
   positions by a grid search over offset and scale, then by least squares
   on the matched ticks. The scan is sheared (the left edge leans about
   0.5 deg, the bottom edge 0.02 deg), so log10 E is the bottom-axis fit and
   the top-axis fit blended linearly in y, and log10 Y the left and right fits
   blended linearly in x. Residuals are printed.
2. Symbols. This figure uses its own letters (legend: A Keywell, B Guseva,
   D Yonts, F Fert, G Laegreid, H Perovic, J Southern, L Dupp, N Koshkin,
   O Akaishi, P Weijsenfeld 1967, Q Bohdansky, U Oechsner, V Okajima,
   Z Bader, ...). Each symbol in SEEDS was read by eye from 8x zooms of this
   raster. The search was guided by where IPPJ-AM-32 puts the same
   reference's point (stated, not hidden), but each position below is read
   and then measured in this raster only. Legibility classes: "isolated" (the
   glyph stands alone or its outline is unambiguous), "coincident" (two
   references' symbols overlap at one spot in both compilations; compared to
   the common centre), "cluster" (inside a stack where no single glyph can be
   separated at 200 ppi; reported, not used for the agreement statistics).
   The figure marks no point of these references as calculated (ACAT); the
   "#" set (Eckstein 1983) is not in IPPJ-AM-32 and is not read.
3. Centres. Two template-free estimates per symbol, from the ink inside its
   hand-drawn box in SEEDS (`centres()`): the centre of the ink's bounding box,
   and the centroid of the ink pixels. The value is their mean, half their
   difference the centring repeatability.
4. Uncertainty (digitizing only): calibration residual rms, centring
   repeatability, and one pixel of symbol placement (0.0085 decade in E,
   0.0076 in Y), in quadrature.

`--compare DIR` pairs every symbol with the stored IPPJ-AM-32 points of the
same reference (nearest within 0.05 decade in E and 0.1 in Y) and reports the
ratio of the two reads against their combined uncertainty; `--write` keeps
that record as JSON.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image

PANEL = (900, 1560, 140, 840)  # x0, x1, y0, y1 of Fig. 120 on the page
# Dataset id (validation/data/sputtering) -> this figure's letter for the same reference.
LETTER_OF = {
    "ar_cu_sputter_keywell1955": "A", "ar_cu_sputter_guseva1960": "B", "ar_cu_sputter_yonts1960": "D",
    "ar_cu_sputter_fert1961": "F", "ar_cu_sputter_laegreid1961": "G", "ar_cu_sputter_perovic1961": "H",
    "ar_cu_sputter_southern1963": "J", "ar_cu_sputter_dupp1966": "L", "ar_cu_sputter_koshkin1969": "N",
    "ar_cu_sputter_akaishi1977": "O", "ar_cu_sputter_weijsenfeld1967": "P", "ar_cu_sputter_bohdansky1980": "Q",
    "ar_cu_sputter_oechsner1973": "U", "ar_cu_sputter_okajima1981": "V", "ar_cu_sputter_bader1961": "Z",
}

# Glyphs read from this raster: (letters, legibility, box x0, y0, x1, y1 in px, note). The box encloses
# the glyph's ink and excludes the fitted curve and neighbouring glyphs.
SEEDS = [
    ("GO", "coincident", 1092, 411, 1101, 425, "about 200 eV: G and O drawn on one spot (also coincident in IPPJ-AM-32)"),
    ("O", "isolated", 1103, 407, 1113, 420, "about 250 eV"),
    ("O", "isolated", 1111, 418, 1123, 433, "about 300 eV, below the curve; its top bar touches the curve"),
    ("A", "isolated", 1134, 407, 1143, 421, "about 450 eV, below the curve"),
    ("ON", "coincident", 1172, 376, 1183, 389, "about 1 keV: O and N on one spot (also coincident in IPPJ-AM-32)"),
    ("B", "isolated", 1255, 338, 1266, 351, "about 5 keV, below the curve"),
    ("DN", "coincident", 1277, 303, 1287, 317, "7.5-8 keV: D and N overlap (also in IPPJ-AM-32)"),
]
# Stacks where no single glyph can be separated at 200 ppi (x0, y0, x1, y1 px).
CLUSTERS = [
    (1100, 380, 1139, 418, "350-450 eV stack"),
    (1138, 355, 1170, 395, "500-800 eV stack"),
    (1165, 340, 1200, 378, "0.8-1.5 keV stack"),
    (1193, 325, 1235, 372, "1.5-2.7 keV stack"),
    (1240, 300, 1300, 345, "3-10 keV stack"),
]
COMB = np.log10(np.array([k * 10.0**n for n in range(-4, 8) for k in range(1, 10)]))


# --- 1. axes: log-comb fit ---

class Axes:
    def __init__(self, ink):
        x0, x1, y0, y1 = PANEL
        sub = ink[y0:y1, x0:x1]
        rows, cols = sub.sum(1), sub.sum(0)
        ys = [i + y0 for i in range(len(rows)) if rows[i] > 0.5 * (x1 - x0)]
        xs = [j + x0 for j in range(len(cols)) if cols[j] > 0.5 * (y1 - y0)]
        self.top, self.bottom = float(np.mean([y for y in ys if y < 500])), float(np.mean([y for y in ys if y > 500]))
        self.left, self.right = float(np.mean([x for x in xs if x < 1200])), float(np.mean([x for x in xs if x > 1200]))
        self.ink = ink
        span_x, span_y = self.right - self.left, self.bottom - self.top
        res = []
        # x: frame spans 10^1 .. 10^6 eV; y: 10^-3 .. 10^2 atoms/ion
        self.fx_bottom, r = self._side(True, self.bottom, -1, 5 / span_x, 1.0, self.left)
        res.append(("bottom", r))
        self.fx_top, r = self._side(True, self.top, +1, 5 / span_x, 1.0, self.left)
        res.append(("top", r))
        self.fy_left, r = self._side(False, self.left, +1, -5 / span_y, -3.0, self.bottom)
        res.append(("left", r))
        self.fy_right, r = self._side(False, self.right, -1, -5 / span_y, -3.0, self.bottom)
        res.append(("right", r))
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
    ap.add_argument("raster", help="PNG of PDF p. 49 from `pdfimages -f 49 -l 49 -png`")
    ap.add_argument("--compare", metavar="DIR", help="pair with the stored IPPJ-AM-32 datasets in DIR")
    ap.add_argument("--write", metavar="JSON", help="write the comparison record (needs --compare)")
    args = ap.parse_args()

    ink = np.asarray(Image.open(args.raster).convert("L")) < 128
    if ink.shape != (2314, 1603):
        sys.exit(f"raster is {ink.shape[1]} x {ink.shape[0]} px; expected 1603 x 2314 (pdfimages, PDF p. 49)")
    ax = Axes(ink)
    px_e = abs(ax.fx_bottom[1])
    px_y = abs(ax.fy_left[1])

    reads = []
    print("\nletters  centre px        E eV     Y     unc E  unc Y  legibility")
    for letters, leg, x0, y0, x1, y1, note in SEEDS:
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
    rows = []
    for f in sorted(Path(args.compare).glob("ar_cu_sputter_*.json")):
        d = json.loads(f.read_text())
        L = LETTER_OF.get(d["id"])
        for p in d["points"]:
            cands = [r for r in reads if L in r["letters"]
                     and abs(math.log10(r["energy_ev"] / p["energy_ev"])) < 0.05
                     and abs(math.log10(r["yield_"] / p["yield"])) < 0.1]
            row = {"id": d["id"], "nifs_letter": L, "energy_ev": p["energy_ev"], "yield_am32": p["yield"]}
            if not cands:
                x, y = ax.to_pixel(math.log10(p["energy_ev"]), math.log10(p["yield"]))
                inside = [c[4] for c in CLUSTERS if c[0] <= x <= c[2] and c[1] <= y <= c[3]]
                row.update(verdict=("inside an unreadable NIFS stack (" + inside[0] + "); not compared") if inside
                           else "no separable symbol of this reference found in NIFS-DATA-23 Fig. 120")
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
    if args.write:
        rec = {
            "format": "lindhard-digitize-crosscheck/1",
            "figure_a": "IPPJ-AM-32 (1983), Ar -> Cu, PDF p. 118 (matsunami1983_ipp_am32.py)",
            "figure_b": "NIFS-DATA-23 (1995), Fig. 120, PDF p. 49 (yamamura1995_nifs23.py)",
            "calibration_b_rms_decade": round(ax.rms, 4),
            "pairing": "same reference, nearest NIFS symbol within 0.05 decade in E and 0.1 in Y",
            "criterion": "|ln(Y_b / Y_a)| against the combined digitizing uncertainty (quadrature)",
            "nifs_reads": [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in reads],
            "pairs": rows,
        }
        Path(args.write).write_text(json.dumps(rec, indent=2) + "\n")
        print(f"wrote {args.write}")
        # Summarize per dataset into its `crosscheck` field.
        for f in sorted(Path(args.compare).glob("ar_cu_sputter_*.json")):
            d = json.loads(f.read_text())
            mine = [r for r in rows if r["id"] == d["id"]]
            comp = [r for r in mine if "yield_ratio" in r]
            if comp:
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
                agreement = (f"{len(comp)} of {len(mine)} points compared; Y ratio NIFS/AM-32 "
                             + ", ".join(f"{r['yield_ratio']:.3f} at {r['energy_ev']:.0f} eV ({r['verdict']}"
                                         + coinc(r)
                                         + ")" for r in comp)
                             + "; the other points lie in stacks not separable at 200 ppi")
            else:
                agreement = (f"present in Fig. 120 (letter {LETTER_OF[d['id']]}), but none of its {len(mine)} points "
                             "is separable from neighbouring symbols at 200 ppi; not compared")
            d["crosscheck"] = {
                "compilation": "Yamamura and Tawara, NIFS-DATA-23 (1995); At. Data Nucl. Data Tables 62, 149 (1996)",
                "figure": "Fig. 120",
                "pdf_page": 49,
                "symbol": LETTER_OF[d["id"]],
                "agreement": agreement + " (validation/data/digitize/crosscheck_ar_cu_nifs23.json)",
            }
            f.write_text(json.dumps(d, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
