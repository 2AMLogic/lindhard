#!/usr/bin/env python3
"""Digitize the bulk optical energy-loss function (ELF) of Si of
L. H. Yang, K. Tokesi, J. Toth, B. Da, H. M. Li and Z. J. Ding, "Optical
properties of silicon and germanium determined by high-precision analysis of
reflection electron energy loss spectroscopy spectra", Phys. Rev. B 100,
245209 (2019), doi:10.1103/PhysRevB.100.245209, Fig. 7 (the "Present" curve,
the ELF averaged over the 3, 4 and 5 keV REELS analyses), with a second read
from Fig. 6(b) of the same paper.

This is the extraction behind `validation/data/optical/si_elf_yang2019.toml`
(provenance: docs/data-provenance.md, issue #125). It is not part of the
harness and needs numpy and Pillow, which the harness does not. Neither the
article nor any figure is in this tree: fetch the article yourself (the
published version, also deposited at the Library of the Hungarian Academy of
Sciences, http://real.mtak.hu/105943/1/PhysRevB.100.245209.pdf, 17 pages,
3,694,542 bytes, SHA-256
99a5f2cc92615f2f483f42e122b810d92c96244b6a3cec6663ddddbe71028648) and extract
the two figure rasters of PDF page 10 with poppler:

    pdfimages -f 10 -l 10 -png PhysRevB.100.245209.pdf p10
    # p10-000.png: Fig. 6 (1800 x 671 px), p10-001.png: Fig. 7 (1020 x 773 px)
    validation/data/digitize/yang2019_si_elf.py p10-001.png p10-000.png
    validation/data/digitize/yang2019_si_elf.py p10-001.png p10-000.png \\
        --write validation/data/optical/si_elf_yang2019.toml \\
        --record validation/data/digitize/secondread_si_yang2019.json

Both rasters are embedded losslessly at 300 ppi with flat colours (pure red
(255, 0, 0) curves with anti-aliased edges, black frames), so every reader
gets the same pixels and the pixel coordinates below are reproducible.

Method (also stated in the dataset's provenance):

1. Axes. Each frame side is a band of rows or columns that is mostly black.
   Every tick is a run of black pixels perpendicular to its frame side,
   measured from the frame outward or inward as the figure draws it; the
   longest ticks are the labelled (major) ones, and a tick's position is the
   length-weighted centre of its rows or columns. A least-squares line maps
   pixels to energy (linear) and to ELF (linear in the main panel, log10 in
   the inset). Ticks that coincide with a frame side are not used. Where the
   curves run over some majors (main panel x beyond 25 eV; Fig. 6(b) y at
   10^-2 and 10^-3), the clean majors place the comb of all ticks (every
   2 eV; k x 10^n) and every clean tick on the comb enters the fit. The script
   stops if a side gives an unexpected number of majors, and prints the
   residuals of each fit and the value it assigns to the frame sides.
2. Curve. A pixel's red weight is (R - max(G, B)) / 255, which is 1 inside
   the red stroke, falls off over its anti-aliased edge and is 0 for the
   black, blue and grey symbols of the other data sets (the red curve is drawn
   on top of them). In each pixel column the red run with the largest total
   weight is taken (grown by its fringe), and its weighted centre row is the
   curve at that column.
3. Panels and joins. Below 3 eV the main panel's stroke lies on the bottom
   frame (the plot area clips it there, so its centre is biased upward by up
   to 0.008 in ELF), so the knots from 0.5 to 2.75 eV come from the log inset,
   whose first column right of its frame is at 0.46 eV. The linear main panel
   (18 px/eV, 159 px per unit of ELF) gives the knots from 3 to 24.75 eV; the
   script checks that its stroke is clear of the frame there. The log inset
   (2.25 px/eV, 72 px/decade) gives the knots from 25 eV to 199 eV, the last
   full pixel column before the right frame (the axis ends at 199.5 eV; the
   paper states its range as 0-200 eV). The join at 25 eV is where the inset's
   left frame meets the main panel's curve; both panels show the curve from
   3 to 50 eV, and their agreement is printed (and recorded, 20 to 49.5 eV).
4. Knots. Every 0.25 eV from 0.5 to 24.75 eV, every 1 eV from 25 to 96 eV
   and from 104 to 199 eV, and every 0.25 eV across the Si L2,3 edge from 96
   to 104 eV. Each knot value is interpolated linearly between the two
   nearest column reads (in log10 ELF in the inset); nothing is smoothed or
   extrapolated. The script prints how far linear interpolation between the
   knots departs from the column reads (inset: median 0.15 %, at most 2.2 %;
   main panel: median 0.0007, at most 0.014 in ELF).
5. Second read (`--record`). Fig. 6(b) of the same paper plots the ELFs from
   the 3, 4 and 5 keV spectra separately on a log axis (3.73 px/eV,
   107 px/decade); the 5 keV curve is red and drawn on top. The paper states
   that the three are "almost the same in the whole energy loss region"
   (p. 11), and Fig. 7 is their average. The same method (steps 1 and 2) reads
   that red curve, independently calibrated, and it is compared with the
   committed knots at fixed energies.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

PROVENANCE = (
    "Bulk energy-loss function Im[-1/eps] of Si, 0.5 to 199 eV, determined from reflection "
    "electron energy loss spectra (REELS) measured at 3, 4 and 5 keV by a reverse Monte Carlo "
    "analysis (full Penn algorithm for the bulk term, Ritchie-Howie surface term), averaged over "
    "the three primary energies: L. H. Yang, K. Tokesi, J. Toth, B. Da, H. M. Li and Z. J. Ding, "
    "Phys. Rev. B 100, 245209 (2019), doi:10.1103/PhysRevB.100.245209, Fig. 7 (red 'Present' "
    "curve; linear main panel 0-50 eV, log inset 0-200 eV). Read from the published version "
    "deposited at http://real.mtak.hu/105943/1/PhysRevB.100.245209.pdf (SHA-256 "
    "99a5f2cc92615f2f483f42e122b810d92c96244b6a3cec6663ddddbe71028648), 2026-10-09. The paper "
    "prints no ELF table, so the curve was DIGITIZED by "
    "validation/data/digitize/yang2019_si_elf.py: axes from the tick marks, the red stroke's "
    "weighted centre per pixel column; knots every 0.25 eV to 24.75 eV (log inset below 3 eV, "
    "where the main panel's stroke lies on the axis; main panel from 3 eV), then from the log "
    "inset every 1 eV to 199 eV (every 0.25 eV across the L2,3 edge, 96 to 104 eV), joined at "
    "25 eV; values interpolated between column reads, not smoothed. Digitizing uncertainty "
    "(half a pixel): about 0.003 in ELF from 3 to 25 eV, about 1.6 % from 25 eV and below 3 eV "
    "(more on the steep flanks below 3 eV and at the edge). Second read: the 5 keV curve of "
    "Fig. 6(b) of the same paper, digitized the same way at 32 energies "
    "(validation/data/digitize/secondread_si_yang2019.json). Caveats: an ELF from the inversion "
    "of electron spectra through a dielectric model, not a direct optical measurement; ends at "
    "199 eV, so the Si K shell (about 1.84 keV) and the L-shell tail above 199 eV are absent. "
    "Operator ruling #34: committed as cited facts. Issue #125."
)

FIG7_SIZE = (1020, 773)
FIG6_SIZE = (1800, 671)

LOW_EV = 0.5  # first knot (first inset column right of its frame: 0.46 eV)
MAIN_EV = 3.0  # main panel from here (below, its stroke sits on the axis)
JOIN_EV = 25.0  # log inset from here
TOP_EV = 199.0  # last knot (the axis ends at 199.5 eV)


def masks(path: str, size: tuple[int, int]):
    im = np.asarray(Image.open(path).convert("RGB")).astype(int)
    h, w, _ = im.shape
    if (w, h) != size:
        sys.exit(f"{path}: expected a {size[0]} x {size[1]} raster, got {w} x {h}")
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    black = (r < 128) & (g < 128) & (b < 128)
    red = np.clip(r - np.maximum(g, b), 0, 255) / 255.0
    return black, red


def runs(flags):
    """Index groups of consecutive True entries."""
    out, cur = [], []
    for i, f in enumerate(flags):
        if f:
            if cur and i != cur[-1] + 1:
                out.append(cur)
                cur = []
            cur.append(i)
    if cur:
        out.append(cur)
    return out


def frame_lines(profile, frac, n):
    """Centres of the bands of rows (or columns) whose black count exceeds frac * n."""
    return [float(np.mean(g)) for g in runs(profile > frac * n)]


def ticks(black, axis, start, step, span, lo, hi):
    """Ticks perpendicular to a frame side.

    `axis` 'y' scans rows lo..hi (ticks on a vertical side), 'x' scans columns
    lo..hi (ticks on a horizontal side). For each row/column the tick length
    is the run of black pixels from `start` in direction `step` (+1 or -1), up
    to `span`. Returns (centre, length) per tick, centre weighted by length.
    """
    lengths = np.zeros(hi - lo)
    for k, i in enumerate(range(lo, hi)):
        n = 0
        while n < span:
            p = start + step * n
            if (black[i, p] if axis == "y" else black[p, i]):
                n += 1
            else:
                break
        lengths[k] = n
    out = []
    for gp in runs(lengths >= 2):
        wts = lengths[gp]
        out.append((lo + float(np.dot(gp, wts) / wts.sum()), float(wts.max()), len(gp)))
    return out


def majors(found, expected, what):
    """The longest ticks (within 2 px of the longest), narrow ones only."""
    found = [t for t in found if t[2] <= 4]
    longest = max(t[1] for t in found)
    maj = sorted(c for c, length, _ in found if length >= longest - 2)
    if len(maj) != expected:
        sys.exit(f"{what}: expected {expected} major ticks, found {maj}")
    return maj


def comb(found, maj, maj_values, candidates, min_len, tol, what):
    """Every tick of length >= min_len (narrow ones) whose value on the line
    through the majors lies within `tol` of one of `candidates`; returns the
    ticks and those candidate values, for the final fit."""
    a = np.polyfit(maj, maj_values, 1)
    cand = np.asarray(candidates, dtype=float)
    px, vals, dropped = [], [], []
    for c, length, width in found:
        if length < min_len or width > 4:
            continue
        u = np.polyval(a, c)
        k = int(np.argmin(np.abs(cand - u)))
        if abs(cand[k] - u) <= tol:
            px.append(c)
            vals.append(float(cand[k]))
        else:
            dropped.append(round(c, 1))
    if len(set(vals)) != len(vals):
        sys.exit(f"{what}: two ticks on one comb value {vals}")
    print(f"{what}: {len(px)} ticks on the comb ({min(vals):g} to {max(vals):g}), "
          f"off-comb dropped {dropped}")
    return px, vals


def log_comb(lo_decade, hi_decade):
    """log10 of k * 10^n, k = 1..9, n = lo..hi-1, and 10^hi."""
    return [n + np.log10(k) for n in range(lo_decade, hi_decade) for k in range(1, 10)] + [hi_decade]


def fit(pixels, values, what, unit):
    a = np.polyfit(pixels, values, 1)
    res = np.abs(np.polyval(a, pixels) - np.asarray(values)).max()
    print(f"{what}: ticks {[round(p, 2) for p in pixels]}, {1 / abs(a[0]):.3f} px/{unit}, "
          f"max residual {res:.4f} {unit}")
    return a


def column_reads(red, cols, rows):
    """Per column: (column, centre row, last row, total weight) of the heaviest
    red run within `rows` (a (lo, hi) pair, or a function of the column)."""
    out = []
    for x in cols:
        lo, hi = rows(x) if callable(rows) else rows
        w = red[lo:hi, x]
        best = None
        for gp in runs(w > 0.25):
            # Grow the run by its anti-aliased fringe (weight > 0).
            a, b = gp[0], gp[-1]
            while a > 0 and w[a - 1] > 0:
                a -= 1
            while b < len(w) - 1 and w[b + 1] > 0:
                b += 1
            seg = np.arange(a, b + 1)
            tot = w[seg].sum()
            if best is None or tot > best[0]:
                best = (tot, seg)
        if best is None:
            continue
        tot, seg = best
        out.append((x, lo + float(np.dot(seg, w[seg]) / tot), lo + int(seg[-1]), float(tot)))
    return out


def read_fig7(path):
    black, red = masks(path, FIG7_SIZE)
    h, w = black.shape
    rows_full = frame_lines(black.sum(1), 0.6, w)  # top and main bottom
    cols_full = frame_lines(black.sum(0), 0.6, h)  # main left and right
    if len(rows_full) != 2 or len(cols_full) != 2:
        sys.exit(f"Fig. 7 frame: rows {rows_full}, columns {cols_full}")
    top, bottom = rows_full
    left, right = cols_full
    # Inset frame: a column with black over 30 % of the height and a row with
    # black over 40 % of the width, other than the main frame.
    inset_left = [c for c in frame_lines(black.sum(0), 0.3, h) if left + 5 < c < right - 5]
    inset_bottom = [r for r in frame_lines(black.sum(1), 0.4, w) if top + 5 < r < bottom - 5]
    if len(inset_left) != 1 or len(inset_bottom) != 1:
        sys.exit(f"Fig. 7 inset frame: columns {inset_left}, rows {inset_bottom}")
    il, ib = inset_left[0], inset_bottom[0]
    left_in = int(np.ceil(left + 1.5))  # first column right of the 3 px frame
    bottom_in = int(np.floor(bottom - 1.0))  # first row above the 2 px frame
    il_out = int(np.floor(il - 1.0))
    ib_in = int(np.floor(ib - 1.0))
    print(f"Fig. 7 frame: main left {left}, right {right}, top {top}, bottom {bottom}; "
          f"inset left {il}, bottom {ib}")

    # Main panel: y ticks inward from the left side, x ticks inward from the
    # bottom. Beyond about 25 eV the red curve and the symbols run over the
    # bottom ticks, so only the 10 and 20 eV majors are clean: they place the
    # comb of all ticks (every 2 eV), and every clean tick on it is fitted.
    my = majors(ticks(black, "y", left_in, +1, 20, int(top) + 4, bottom_in - 3), 4, "main y")
    xt = ticks(black, "x", bottom_in, -1, 20, left_in + 3, int(right) - 3)
    mx, mxv = comb(xt, majors(xt, 2, "main x"), [10, 20], np.arange(0, 51, 2), 5, 0.2, "Fig. 7 main x")
    ay = fit(my, [4, 3, 2, 1], "Fig. 7 main y (ELF)", "ELF")
    ax = fit(mx, mxv, "Fig. 7 main x (eV)", "eV")
    # Inset: log y ticks outward to the left of its left side (down to its
    # bottom side, 10^-4), x ticks inward from its bottom side.
    iy = majors(ticks(black, "y", il_out, -1, 20, int(top) + 4, int(ib) + 2), 5, "inset y")
    ix = majors(ticks(black, "x", ib_in, -1, 20, int(il) + 3, int(right) - 3), 3, "inset x")
    by = fit(iy, [0, -1, -2, -3, -4], "Fig. 7 inset y (log10 ELF)", "decade")
    bx = fit(ix, [50, 100, 150], "Fig. 7 inset x (eV)", "eV")
    print(f"  checks: main x at left frame {np.polyval(ax, left):.3f} eV (0), at right frame "
          f"{np.polyval(ax, right):.3f} eV (50); main ELF at bottom frame "
          f"{np.polyval(ay, bottom):.4f} (0); inset x at its left frame {np.polyval(bx, il):.3f} eV (0)")

    # Main panel reads: left of the inset the whole height; under the inset
    # only below ELF = 1 (the legend's red sample line is at about 1.5).
    below_1 = int(np.polyval(np.polyfit([4, 3, 2, 1], my, 1), 1.0))

    def main_rows(x):
        return (int(top) + 2, bottom_in + 2) if x < il - 1 else (below_1, bottom_in + 2)

    main = column_reads(red, range(left_in, int(right) - 1), main_rows)
    main_e = np.array([np.polyval(ax, t[0]) for t in main])
    main_v = np.array([np.polyval(ay, t[1]) for t in main])
    # A stroke that reaches the frame (row bottom_in + 1 or below) is clipped
    # by the plot area there, so its centre is biased upward.
    touch_e = [e for e, t in zip(main_e, main) if t[2] >= bottom_in + 1 and e < 25]
    inset = column_reads(red, range(int(il) + 2, int(right) - 1), (int(top) + 2, ib_in))
    ins_e = np.array([np.polyval(bx, t[0]) for t in inset])
    ins_v = np.array([np.polyval(by, t[1]) for t in inset])  # log10 ELF
    print(f"  main panel: {len(main)} columns read; below 25 eV the stroke reaches the bottom frame up to "
          f"{max(touch_e):.2f} eV; inset: {len(inset)} columns, {ins_e.min():.2f} to {ins_e.max():.2f} eV")
    return dict(main_e=main_e, main_v=main_v, ins_e=ins_e, ins_v=ins_v, touch_ev=max(touch_e),
                main_px_per_ev=1 / ax[0], main_px_per_elf=-1 / ay[0],
                ins_px_per_ev=1 / bx[0], ins_px_per_decade=-1 / by[0])


def knot_energies():
    e = list(np.arange(LOW_EV, MAIN_EV - 0.01, 0.25))
    e += list(np.arange(MAIN_EV, JOIN_EV - 0.01, 0.25))
    e += list(np.arange(JOIN_EV, 96.0 - 0.01, 1.0))
    e += list(np.arange(96.0, 104.0 - 0.01, 0.25))
    e += list(np.arange(104.0, TOP_EV + 0.01, 1.0))
    return np.round(np.array(e), 4)


def knots(d):
    e = knot_energies()
    main = (e >= MAIN_EV) & (e < JOIN_EV)
    if d["touch_ev"] >= MAIN_EV - 0.1:  # the two columns around MAIN_EV must be clear
        sys.exit(f"main-panel stroke touches the frame at {d['touch_ev']:.2f} eV, inside the main-panel knots")
    if e[~main].min() < d["ins_e"].min() or e[~main].max() > d["ins_e"].max():
        sys.exit("inset knots outside the columns read")
    v = np.empty_like(e)
    v[main] = np.interp(e[main], d["main_e"], d["main_v"])
    v[~main] = 10 ** np.interp(e[~main], d["ins_e"], d["ins_v"])
    if (v <= 0).any():
        sys.exit(f"non-positive ELF read at {e[v <= 0]}")
    print(f"\n{len(e)} knots: {int((e < MAIN_EV).sum())} from the inset below {MAIN_EV} eV, "
          f"{int(main.sum())} from the main panel, {int((e >= JOIN_EV).sum())} from the inset from {JOIN_EV} eV")
    return e, v


def report(d, e, v):
    i = int(np.argmax(v))
    print(f"plasmon peak: ELF {v[i]:.3f} at {e[i]:.2f} eV (knots); "
          f"column max {d['main_v'].max():.3f} at {d['main_e'][np.argmax(d['main_v'])]:.2f} eV")
    # L2,3 edge: minimum before and maximum after 99 eV.
    pre = (e > 90) & (e < 100)
    post = (e > 100) & (e < 120)
    print(f"L2,3 edge: minimum {v[pre].min():.2e} at {e[pre][np.argmin(v[pre])]:.2f} eV, "
          f"maximum {v[post].max():.3e} at {e[post][np.argmax(v[post])]:.2f} eV")
    for lo_e, hi_e in ((0.5, 1.0), (1.0, 3.0), (3.0, 6.0)):
        sel = (e >= lo_e) & (e <= hi_e)
        print(f"ELF {lo_e}-{hi_e} eV: {v[sel].min():.4f} to {v[sel].max():.4f}")
    # Below 3 eV: main panel against the inset (the main stroke sits on the axis).
    for e0 in (1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0):
        print(f"  {e0:4.1f} eV: main {np.interp(e0, d['main_e'], d['main_v']):+.4f}, "
              f"inset {10 ** np.interp(e0, d['ins_e'], d['ins_v']):.4f}")
    # Panel overlap: both panels from 20 to 50 eV.
    ov = (d["main_e"] >= 20) & (d["main_e"] <= 49.5)
    me, mv = d["main_e"][ov], d["main_v"][ov]
    iv = 10 ** np.interp(me, d["ins_e"], d["ins_v"])
    ratio = mv / iv
    print(f"panel overlap 20-49.5 eV ({ov.sum()} main columns): main/inset median {np.median(ratio):.3f}, "
          f"range {ratio.min():.3f} to {ratio.max():.3f}")
    rows = []
    for e0 in (20.0, 22.0, 24.0, 25.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 49.0):
        m = float(np.interp(e0, d["main_e"], d["main_v"]))
        n = float(10 ** np.interp(e0, d["ins_e"], d["ins_v"]))
        rows.append({"energy_ev": e0, "elf_main_panel": round(m, 4), "elf_inset": round(n, 4),
                     "ratio_main_over_inset": round(m / n, 3)})
        print(f"  {e0:5.1f} eV: main {m:.4f}, inset {n:.4f}, ratio {m / n:.3f}")
    # Representation: linear interpolation between knots against every column read.
    lin = np.interp(d["ins_e"][d["ins_e"] >= JOIN_EV], e, v)
    col = 10 ** d["ins_v"][d["ins_e"] >= JOIN_EV]
    rel = np.abs(lin / col - 1)
    print(f"knots vs column reads, inset 25-199 eV: max |rel| {rel.max():.3f}, median {np.median(rel):.4f}")
    sel = (d["main_e"] >= MAIN_EV) & (d["main_e"] < JOIN_EV)
    lin = np.interp(d["main_e"][sel], e, v)
    ab = np.abs(lin - d["main_v"][sel])
    print(f"knots vs column reads, main {MAIN_EV}-{JOIN_EV} eV: max |abs| {ab.max():.4f}, median {np.median(ab):.4f}")
    return rows, float(np.median(ratio)), float(ratio.min()), float(ratio.max())


def read_fig6(path):
    black, red = masks(path, FIG6_SIZE)
    h, w = black.shape
    # Panel (b) is the right half; its frame is the black rows and columns there.
    half = w // 2
    rows_b = frame_lines(black[:, half:].sum(1), 0.6, w - half)
    cols_b = [c for c in frame_lines(black.sum(0), 0.6, h) if c > half]
    if len(rows_b) != 2 or len(cols_b) != 2:
        sys.exit(f"Fig. 6(b) frame: rows {rows_b}, columns {cols_b}")
    top, bottom = rows_b
    left, right = cols_b
    left_in = int(np.ceil(left + 1.0))
    bottom_in = int(np.floor(bottom - 1.0))
    print(f"\nFig. 6(b) frame: left {left}, right {right}, top {top}, bottom {bottom}")
    # y: inward from the left side. The steep low-energy branch of the curves
    # runs over the 10^-2 and 10^-3 majors, so the clean 10^1, 10^0 and 10^-1
    # majors place the log comb (k x 10^n) and every clean tick on it is
    # fitted. x: inward from the bottom side, 20 .. 180 eV (0 and 200 eV are
    # the frames).
    yt = ticks(black, "y", left_in, +1, 20, int(top) + 4, bottom_in - 3)
    ly, lyv = comb(yt, majors(yt, 3, "Fig. 6(b) y"), [1, 0, -1], log_comb(-4, 1), 4, 0.02, "Fig. 6(b) y")
    lx = majors(ticks(black, "x", bottom_in, -1, 20, left_in + 3, int(right) - 3), 9, "Fig. 6(b) x")
    ay = fit(ly, lyv, "Fig. 6(b) y (log10 ELF)", "decade")
    ax = fit(lx, [20, 40, 60, 80, 100, 120, 140, 160, 180], "Fig. 6(b) x (eV)", "eV")
    print(f"  checks: x at left frame {np.polyval(ax, left):.3f} eV (0), at right frame "
          f"{np.polyval(ax, right):.3f} eV (200); log10 ELF at bottom frame {np.polyval(ay, bottom):.3f} (-4)")
    # The legend's red sample line (lower right, at about 2e-3 beyond 130 eV)
    # is excluded: there only rows above 5e-3 are read (the curve is near 2e-2).
    inv_ay = np.polyfit(lyv, ly, 1)
    legend_top = int(np.polyval(inv_ay, np.log10(5e-3)))

    def rows(x):
        return (int(top) + 2, legend_top if np.polyval(ax, x) > 130 else bottom_in)

    reads = column_reads(red, range(left_in + 1, int(right) - 1), rows)
    e = np.array([np.polyval(ax, t[0]) for t in reads])
    v = np.array([np.polyval(ay, t[1]) for t in reads])
    print(f"  {len(reads)} columns, {e.min():.2f} to {e.max():.2f} eV, {1 / ax[0]:.3f} px/eV, "
          f"{-1 / ay[0]:.2f} px/decade")
    return e, v, 1 / ax[0], -1 / ay[0]


SECOND_READ_EV = (3.0, 5.0, 8.0, 10.0, 12.0, 14.0, 15.0, 16.0, 16.5, 17.0, 18.0, 19.0, 20.0, 22.0,
                  25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 95.0, 97.0, 102.0, 105.0, 110.0,
                  120.0, 140.0, 160.0, 180.0, 198.0)


def second_read(d, e, v, f6):
    e6, v6, px_ev6, px_dec6 = f6
    # Digitizing uncertainty of a read (half a pixel in each axis, through the
    # local slope), per figure, as a relative uncertainty of the ELF.
    def unc7(x):
        h = 0.5
        if x < JOIN_EV:
            slope = (np.interp(x + h, e, v) - np.interp(x - h, e, v)) / (2 * h)  # ELF/eV
            val = np.interp(x, e, v)
            return float(np.hypot(0.5 / d["main_px_per_elf"], slope * 0.5 / d["main_px_per_ev"]) / val)
        lv = np.log10(np.maximum(v, 1e-12))
        slope = (np.interp(x + h, e, lv) - np.interp(x - h, e, lv)) / (2 * h)  # decade/eV
        return float(np.log(10) * np.hypot(0.5 / d["ins_px_per_decade"], slope * 0.5 / d["ins_px_per_ev"]))

    def unc6(x):
        h = 0.5
        slope = (np.interp(x + h, e6, v6) - np.interp(x - h, e6, v6)) / (2 * h)
        return float(np.log(10) * np.hypot(0.5 / px_dec6, slope * 0.5 / px_ev6))

    pts = []
    print("\nSecond read, Fig. 6(b) 5 keV curve against the committed Fig. 7 knots:")
    print("  E/eV   Fig.7      Fig.6(b)   ratio   u7     u6     verdict")
    for x in SECOND_READ_EV:
        a = float(np.interp(x, e, v))
        b = float(10 ** np.interp(x, e6, v6))
        u7, u6 = unc7(x), unc6(x)
        u = float(np.hypot(u7, u6))
        r = b / a
        verdict = ("within the combined digitizing uncertainty" if abs(np.log(r)) <= u
                   else "within twice the combined digitizing uncertainty" if abs(np.log(r)) <= 2 * u
                   else "outside twice the combined digitizing uncertainty")
        print(f"  {x:6.1f} {a:.4e} {b:.4e} {r:6.3f}  {u7:.3f}  {u6:.3f}  {verdict}")
        pts.append({"energy_ev": x, "elf_fig7_committed": round(a, 5), "elf_fig6b_5kev": round(b, 5),
                    "ratio_fig6b_over_fig7": round(r, 3), "unc_rel_fig7": round(u7, 3),
                    "unc_rel_fig6b": round(u6, 3), "verdict": verdict})
    n_in = sum(p["verdict"].startswith("within the") for p in pts)
    n_2 = sum(p["verdict"].startswith("within twice") for p in pts)
    print(f"  {n_in} of {len(pts)} within the combined uncertainty, {n_2} more within twice it")
    return pts


def write_toml(path, e, v):
    lines = ['material = "Si"', f'provenance = "{PROVENANCE}"', "energy_ev = ["]
    lines += [f"  {x:.4f}," for x in e]
    lines += ["]", "elf = ["]
    lines += [f"  {x:.4e}," for x in v]
    lines += ["]"]
    Path(path).write_text("\n".join(lines) + "\n")
    print(f"wrote {path}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fig7", help="Fig. 7 raster (p10-001.png)")
    p.add_argument("fig6", nargs="?", help="Fig. 6 raster (p10-000.png), for the second read")
    p.add_argument("--write", help="write the OpticalElf TOML here")
    p.add_argument("--record", help="write the second-read record (JSON) here; needs fig6")
    args = p.parse_args()

    d = read_fig7(args.fig7)
    e, v = knots(d)
    overlap, med, rmin, rmax = report(d, e, v)
    rec = None
    if args.fig6:
        pts = second_read(d, e, v, read_fig6(args.fig6))
        rec = {
            "format": "lindhard-digitize-secondread/1",
            "figure": "Yang et al., Phys. Rev. B 100, 245209 (2019), Fig. 7 (PDF p. 10), red 'Present' curve "
                      "(yang2019_si_elf.py)",
            "automatic": "Fig. 7: weighted centre of the red stroke per pixel column, main panel below 25 eV, "
                         "log inset from 25 eV; the committed knots (si_elf_yang2019.toml), linearly interpolated",
            "second_read": "Fig. 6(b) of the same paper (PDF p. 10), the red 5 keV curve (drawn on top of the "
                           "3 and 4 keV curves), read by the same per-column method on its own, independently "
                           "calibrated log axes (3.73 px/eV, 107 px/decade). Fig. 7 is the average over 3, 4 and "
                           "5 keV, which the paper calls 'almost the same in the whole energy loss region' (p. 11), "
                           "so the two reads differ by the digitizing and by that spread",
            "uncertainty": "unc_rel_*: half a pixel in each axis of that figure, through the local slope; "
                           "the criterion is |ln(ratio)| against sqrt(unc7^2 + unc6^2)",
            "panel_overlap": {
                "what": "Fig. 7 main (linear) panel against its log inset, every main-panel column from 20 to "
                        "49.5 eV (the join is at 25 eV)",
                "ratio_median": round(med, 3), "ratio_min": round(rmin, 3), "ratio_max": round(rmax, 3),
                "points": overlap,
            },
            "points": pts,
        }
    if args.write:
        write_toml(args.write, e, v)
    if args.record:
        if rec is None:
            sys.exit("--record needs the Fig. 6 raster")
        Path(args.record).write_text(json.dumps(rec, indent=2) + "\n")
        print(f"wrote {args.record}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
