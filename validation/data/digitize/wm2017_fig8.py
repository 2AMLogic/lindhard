#!/usr/bin/env python3
"""Digitize the low-energy B-in-a-Si ranges of Wittmaack and Mutzke,
J. Appl. Phys. 121, 105104 (2017), Fig. 8 (doi:10.1063/1.4978016).

This is the extraction behind `validation/data/ranges/b_*keV_asi_wach1982.json`
(provenance: docs/data-provenance.md). It is not part of the harness and needs
numpy and Pillow, which the harness does not. Neither the article nor the
figure is in this tree: fetch the published version yourself from the
Helmholtz Munich repository
(https://push-zb.helmholtz-munich.de/frontdoor.php?source_opus=50831) and
extract the figure raster with poppler:

    pdfimages -f 7 -l 7 -png <article.pdf> fig8
    validation/data/digitize/wm2017_fig8.py fig8-000.png

Method (also stated in each dataset's `extraction`):

1. Axes. The embedded raster is 952 x 733 px at 300 ppi. log10(E) and
   log10(Rp) are linear in pixels; each is a least-squares line through the
   centres of the major (decade) tick marks, found as the longest ticks on the
   bottom and right axes.
2. Circles ("raw experimental data", filled black, 15 px across in the
   legend). The blue fit line crosses every marker, so a plain centroid of
   the black pixels is biased. Each centre is found by matching a filled disk
   of radius R against black-or-blue pixels, minus twice the black pixels in
   the ring R..R+3, on a 0.5 px then 0.1 px grid; R = 7.0 and 7.5 px, and the
   value is their mean.
3. Second read. The blue "selected data (mean) x 0.5" triangles at the same
   energies, doubled, by bounding-box centre and by area centroid.

Only the points below 25 keV are used: the paper's Table I lists no set other
than Wach and Wittmaack (a-Si, 0 deg) there.
"""

from __future__ import annotations

import sys
from collections import deque

import numpy as np
from PIL import Image

NOMINAL_KEV = (1, 2, 3, 5, 10, 20)


def main(path: str) -> int:
    im = np.asarray(Image.open(path).convert("RGB")).astype(int)
    h, w, _ = im.shape
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    black = (r < 110) & (g < 110) & (b < 110)
    blue = (b > 140) & (r < 120) & (g < 140) & ~black

    # 1. Axis calibration from the major ticks (longest tick marks).
    def tick_centres(profile, lo, hi):
        groups, cur = [], []
        for i in range(lo, hi):
            if profile[i] >= 1:
                if cur and i != cur[-1] + 1:
                    groups.append(cur)
                    cur = []
                cur.append(i)
        if cur:
            groups.append(cur)
        # Centre weighted by the black-pixel count of each row/column, so an
        # anti-aliased edge row does not shift it by half a pixel.
        ticks = [
            (sum(i * profile[i] for i in gp) / sum(profile[i] for i in gp), max(profile[i] for i in gp))
            for gp in groups
        ]
        longest = max(length for _, length in ticks)
        return [c for c, length in ticks if length >= longest - 2]

    # Inner edges of the (2 to 3 px wide) bottom and right frame lines.
    frame_bottom = min(i for i in range(h // 2, h) if black[i, :].sum() > 0.5 * w)
    frame_right = min(j for j in range(w // 2, w) if black[:, j].sum() > 0.5 * h)
    xmaj = tick_centres(black[frame_bottom - 26 : frame_bottom, :].sum(0), 115, frame_right - 3)
    ymaj = tick_centres(black[:, frame_right - 24 : frame_right].sum(1), 6, frame_bottom)
    if len(xmaj) != 5 or len(ymaj) != 4:
        sys.exit(f"unexpected major ticks: x {xmaj}, y {ymaj}")
    ax = np.polyfit(xmaj, [0, 1, 2, 3, 4], 1)  # log10(E / keV)
    ay = np.polyfit(ymaj, [4, 3, 2, 1], 1)  # log10(Rp / nm)
    print(f"x major ticks {[round(float(c), 2) for c in xmaj]}, {1 / ax[0]:.2f} px/decade, max residual "
          f"{np.abs(np.polyval(ax, xmaj) - [0, 1, 2, 3, 4]).max():.4f} decade")
    print(f"y major ticks {[round(float(c), 2) for c in ymaj]}, {-1 / ay[0]:.2f} px/decade, max residual "
          f"{np.abs(np.polyval(ay, ymaj) - [4, 3, 2, 1]).max():.4f} decade")

    # 2. Circles by disk-template matching.
    cover = black | blue
    yy, xx = np.mgrid[0:h, 0:w]

    def score(cx, cy, rad):
        d = np.hypot(xx - cx, yy - cy)
        return cover[d <= rad].sum() - 2 * black[(d > rad) & (d <= rad + 3)].sum()

    def x_of(e_kev):
        return (np.log10(e_kev) - ax[1]) / ax[0]

    circles = {}
    for e in NOMINAL_KEV:
        x0 = x_of(e)
        col = black[300:600, int(round(x0)) - 2 : int(round(x0)) + 3].sum(1)
        y0 = float(np.median([i + 300 for i in range(len(col)) if col[i] >= 4]))
        reads = []
        for rad in (7.0, 7.5):
            best = max(
                (score(cx, cy, rad), cx, cy)
                for cx in np.arange(x0 - 8, x0 + 8.01, 0.5)
                for cy in np.arange(y0 - 8, y0 + 8.01, 0.5)
            )
            _, cx, cy = best
            best = max(
                (score(cx2, cy2, rad), cx2, cy2)
                for cx2 in np.arange(cx - 0.5, cx + 0.51, 0.1)
                for cy2 in np.arange(cy - 0.5, cy + 0.51, 0.1)
            )
            _, cx, cy = best
            reads.append((10 ** np.polyval(ax, cx), 10 ** np.polyval(ay, cy)))
        circles[e] = reads

    # 3. Triangles (selected means x 0.5), connected components of blue.
    seen = np.zeros_like(blue)
    triangles = []
    for y0 in range(300, frame_bottom - 10):
        for x0 in range(140, int(x_of(25))):
            if blue[y0, x0] and not seen[y0, x0]:
                q, pts = deque([(y0, x0)]), []
                seen[y0, x0] = True
                while q:
                    y, x = q.popleft()
                    pts.append((y, x))
                    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        y2, x2 = y + dy, x + dx
                        if 0 <= y2 < h and 0 <= x2 < w and blue[y2, x2] and not seen[y2, x2]:
                            seen[y2, x2] = True
                            q.append((y2, x2))
                p = np.array(pts)
                ph = p[:, 0].max() - p[:, 0].min() + 1
                pw = p[:, 1].max() - p[:, 1].min() + 1
                if 60 < len(p) < 200 and 10 <= pw <= 18 and 10 <= ph <= 18:
                    bx = (p[:, 1].max() + p[:, 1].min()) / 2
                    by = (p[:, 0].max() + p[:, 0].min()) / 2
                    triangles.append(
                        (10 ** np.polyval(ax, bx), 2 * 10 ** np.polyval(ay, by), 2 * 10 ** np.polyval(ay, p[:, 0].mean()))
                    )

    print("\nE nominal | circle E, Rp (R=7.0) | circle E, Rp (R=7.5) | Rp mean | triangle x2 (box, centroid)")
    for e in NOMINAL_KEV:
        (e1, r1), (e2, r2) = circles[e]
        tri = min(triangles, key=lambda t: abs(np.log(t[0] / e)))
        print(f"{e:>3} keV | {e1:6.3f} {r1:7.3f} | {e2:6.3f} {r2:7.3f} | {(r1 + r2) / 2:7.2f} nm "
              f"| {tri[1]:7.3f} {tri[2]:7.3f}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
