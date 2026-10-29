"""Digitize the measured points on one srim.org stopping plot (top panel only).

Usage: python3 -I digitize.py <STOPzzZZ.gif> [<Cit zz-ZZ.pdf>]

Reads only the top panel's coloured markers ("*H"/"*D" for hydrogen, "He"
for helium). SRIM's fitted curve (black) and the lower Experiment/Theory
panel are never read as data. Prints JSON: axis calibration, the legend rows
(counts, old citation numbers, keys joined from the citation PDF), points
(energy keV/amu, stopping, uncertainties, colour, candidate keys) and
per-colour count checks. See README.md.
"""
import json, math, os, re, subprocess, sys, tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ppm

WHITE, BLACK, GRID = (255, 255, 255), (0, 0, 0), (85, 255, 255)


def frame(w, h, r):
    """Top panel box: (x0, x1, y0, y1) from the long black frame lines."""
    # Count black pixels rather than the longest run: data labels printed
    # over a frame line break it into pieces.
    rows = [y for y in range(h) if sum(1 for x in range(w) if r[y][x] == BLACK) > w * 0.5]
    groups = []
    for y in rows:
        if groups and y - groups[-1][-1] <= 1:
            groups[-1].append(y)
        else:
            groups.append([y])
    top, bottom = groups[0][-1], groups[1][0]
    cols = [x for x in range(w) if sum(1 for y in range(top, bottom) if r[y][x] == BLACK) > (bottom - top) * 0.8]
    return cols[0] + (1 if cols[1] == cols[0] + 1 else 0), cols[-1] - (1 if cols[-2] == cols[-1] - 1 else 0), top, bottom


def gridlines(r, x0, x1, y0, y1):
    vx = [x for x in range(x0 + 1, x1) if sum(1 for y in range(y0 + 1, y1) if r[y][x] == GRID) > (y1 - y0) * 0.6]
    hy = [y for y in range(y0 + 1, y1) if sum(1 for x in range(x0 + 1, x1) if r[y][x] == GRID) > (x1 - x0) * 0.6]
    return vx, hy


def ocr(gif, box, scale=4, whitelist="0123456789.", colour=BLACK):
    x, y, w, h = box
    hexc = "#%02X%02X%02X" % colour
    with tempfile.TemporaryDirectory() as d:
        png = os.path.join(d, "c.png")
        subprocess.run(["magick", gif, "-crop", f"{w}x{h}+{x}+{y}", "+repage", "-fill", "white", "+opaque", hexc,
                        "-fill", "black", "-opaque", hexc, "-scale", f"{scale * 100}%",
                        "-bordercolor", "white", "-border", "20", png], check=True)
        out = subprocess.run(["tesseract", png, "-", "--psm", "7", "-c", f"tessedit_char_whitelist={whitelist}"],
                             capture_output=True, text=True).stdout
    return out.strip()


def calibrate_x(gif, r, x0, x1, y1, vx):
    """Log axis: decade lines are those followed by the 2,4,6,8 minor pattern."""
    xs = [x0] + vx + [x1]
    best = None
    # Try each candidate decade width D from consecutive line spacing.
    for i in range(len(xs) - 1):
        for j in range(i + 1, len(xs)):
            D = xs[j] - xs[i]
            if D < 60:
                continue
            pred = [xs[i] + D * k + D * math.log10(m) for k in range(-12, 12) for m in (1, 2, 4, 6, 8)]
            pred = [p for p in pred if x0 - 2 <= p <= x1 + 2]
            score = sum(1 for p in pred if min(abs(p - q) for q in xs) <= 1.5)
            if best is None or (score, -abs(len(pred) - len(xs))) > best[0]:
                best = ((score, -abs(len(pred) - len(xs))), xs[i], D)
    _, anchor, D = best
    decades = [anchor + D * k for k in range(-12, 12) if x0 - 2 <= anchor + D * k <= x1 + 2]
    # Value of the leftmost decade line from its tick label (e.g. "1", "10", or "10" + superscript n).
    labels = []
    for xd in decades:
        txt = ocr(gif, (int(xd) - 22, y1 + 4, 44, 26))
        labels.append(txt)
    exps = []
    for k, txt in enumerate(labels):
        t = txt.replace(".", "")
        if t == "1":
            exps.append(0 - k)
        elif t.startswith("10") and len(t) == 3:
            exps.append(int(t[2]) - k)
        elif t == "10":
            exps.append(1 - k)
    if not exps:
        raise ValueError(f"x tick labels unreadable: {labels}")
    e0 = max(set(exps), key=exps.count)
    agree = exps.count(e0)
    return {"x_decade_px": D, "x_first_decade_px": decades[0], "x_first_decade_log10": e0,
            "x_labels_ocr": labels, "x_label_agreement": f"{agree}/{len(exps)}"}


def calibrate_y(gif, r, x0, y0, y1, hy):
    """Linear stopping axis. Every plot's axis starts at 0 on the frame's
    bottom edge and has evenly spaced grid lines, so gridline n (counted up
    from the bottom) is worth n * step. Each readable label votes for a step;
    the commonest wins, and the fit then uses every grid line."""
    from collections import Counter
    lines = sorted(set(hy) | {y1}, reverse=True)  # bottom first
    gaps = sorted(b - a for a, b in zip(sorted(lines), sorted(lines)[1:]))
    g = gaps[len(gaps) // 2]
    pairs, votes = [], Counter()
    for y in lines:
        n = round((y1 - y) / g)
        txt = ocr(gif, (x0 - 62, y - 11, 58, 22))
        try:
            v = float(txt)
        except ValueError:
            continue
        pairs.append((y, v))
        if n > 0 and v > 0:
            votes[float("%.3g" % (v / n))] += 1
    if not votes:
        raise ValueError(f"y tick labels unreadable: {pairs}")
    step, nvote = votes.most_common(1)[0]
    pts = [(y, step * round((y1 - y) / g)) for y in lines]
    n = len(pts)
    my = sum(p for p, _ in pts) / n
    mv = sum(v for _, v in pts) / n
    a = sum((p - my) * (v - mv) for p, v in pts) / sum((p - my) ** 2 for p, _ in pts)
    b = mv - a * my
    agree = sum(1 for y, v in pairs if abs(v - step * round((y1 - y) / g)) < 1e-6 * max(1, v))
    return {"y_per_px": a, "y_at_px0": b, "y_ticks": pairs, "y_step": step,
            "y_label_agreement": f"{agree}/{len(pairs)}"}


PLUS = [(dx, 0) for dx in range(-2, 3)] + [(0, dy) for dy in (-2, -1, 1, 2)]
DIAG = [(-1, -1), (1, -1), (-1, 1), (1, 1), (-2, -1), (2, -1), (-2, 1), (2, 1)]


def find_stars(r, x0, x1, y0, y1):
    """Data asterisks. The glyph is rasterized in several mirror variants, but
    every variant has a full 5-pixel horizontal and vertical bar through the
    point plus some diagonal pixels. Pixels covered by another colour count as
    hidden, not missing; only white where the bar should be counts against.
    Every marker colour present on the plus is tested, since the centre pixel
    is often covered by another marker's label or by the black curve."""
    hits = []
    for y in range(y0 + 3, y1 - 4):
        for x in range(x0 + 3, x1 - 3):
            for c in {r[y + dy][x + dx] for dx, dy in PLUS}:
                if c == WHITE:
                    continue
                if (is_star(r, x, y, c, x1) or is_small_star(r, x, y, c, x1)) and not in_thick_text(r, x, y, c, x1) and not on_curve(r, x, y, c):
                    on = sum(1 for dx, dy in PLUS if r[y + dy][x + dx] == c)
                    diag = sum(1 for dx, dy in DIAG if r[y + dy][x + dx] == c)
                    hits.append((x, y, c, on + diag))
    hits.sort(key=lambda t: -t[3])
    kept = []
    for h in hits:
        if all(h[2] != k[2] or abs(h[0] - k[0]) > 2 or abs(h[1] - k[1]) > 2 for k in kept):
            kept.append(h)
    return kept


MIN_ON, MAX_MISS, MIN_DIAG, LABEL_MIN_ON = 7, 0, 2, 3
# The plotted point is the centre of the two-character symbol ("*H", "*D"),
# not the asterisk: fitted on 1981-Sant (D in Si, 0.2-2.0 MeV in 0.1 MeV
# steps, 19 points) as 4.72 +/- 0.27 px right of the asterisk centre.
ANCHOR_DX = 4.72
# "He" markers have no asterisk; the point is again the centre of the
# string. Fitted on 1974-Bagl (2 MeV He; Rh, Hf and W plots agree): 7.17 px
# right of the H's left stroke.
HE_ANCHOR_DX = 7.0
# The title block ("Ion = ...", "Target = ...") is drawn in pure blue over
# the top of the panel.
TITLE_PX = 60
# Pairing the Experiment/Theory panel with the top panel: a lower point is the
# same as a top point of its colour within 1.5 px in x and within this many
# top-panel pixels in stopping (about three times the scatter of ratio x
# curve against a direct reading, 0.7-1.6 px on the plots checked).
LOWER_PAIR_PX, LOWER_DEFAULT_PX = 3.0, 1.6
# Below the title, the annotations ("Plotted are Ion Names", V_F, <I>) are
# pure blue text in a fixed box at the top left; "He" in them reads as a marker.
ANNOT_BLUE, ANNOT_BOX = (0, 0, 255), (300, 255)


def is_star(r, x, y, c, x1):
    on = miss = 0
    for dx, dy in PLUS:
        p = r[y + dy][x + dx]
        if p == c:
            on += 1
        elif p == WHITE:
            miss += 1
    if on < MIN_ON or miss > MAX_MISS:
        return False
    if sum(1 for dx, dy in DIAG if r[y + dy][x + dx] == c) < MIN_DIAG:
        return False
    # Long straight runs are grid lines or text strokes, not an asterisk.
    if r[y][x - 3] == c and r[y][x + 3] == c:
        return False
    if r[y - 3][x] == c and r[y + 3][x] == c:
        return False
    return has_label(r, x, y, c, x1)


def on_curve(r, x, y, c):
    """Black is also the SRIM curve and the annotation text. A black marker
    stands alone: no black left of it, or above or below it."""
    if c != BLACK:
        return False
    ring = [(dx, dy) for dx in range(-7, -3) for dy in range(-6, 7)]
    ring += [(dx, dy) for dx in range(-3, 4) for dy in (-7, -6, -5, 5, 6, 7)]
    return sum(1 for dx, dy in ring if r[y + dy][x + dx] == BLACK) > 3


def in_thick_text(r, x, y, c, x1):
    """Title and annotation text is a large font with strokes several pixels
    thick; marker and label strokes are one pixel. A 4x4 solid block of the
    colour near the candidate means it sits on text, not a marker."""
    for by in range(y - 5, y + 3):
        for bx in range(x - 5, min(x + 8, x1 - 4)):
            if all(r[by + j][bx + i] == c for i in range(4) for j in range(4)):
                return True
    return False


def has_label(r, x, y, c, x1):
    """The ion label follows every marker: its first stroke is a vertical bar
    4-8 px to the right spanning rows -4..+4 (H, D, L all start with one)."""
    for k in range(4, 9):
        xs = x + k
        if xs >= x1:
            continue
        col = [r[y + dy][xs] for dy in range(-4, 5)]
        # Another marker's ink may cover part of the stroke; white may not.
        on = sum(1 for p in col if p == c)
        if on < LABEL_MIN_ON or sum(1 for p in col if p == WHITE) > 2:
            continue
        # One pixel wide: only crossbars and the D's bars touch either side.
        if all(sum(1 for dy in range(-4, 5) if r[y + dy][xs + s] == c) <= 5 for s in (-1, 1)):
            return True
    return False


def label_letter(r, x, y, c, x1):
    """'H' or 'D' for the label after a marker at (x, y), '?' if unreadable.
    H has a crossbar on the marker's row between its strokes; D has none, and
    bars on its top and bottom rows instead."""
    for k in range(4, 9):
        xs = x + k
        if xs + 5 >= x1 or r[y][xs] != c:
            continue
        col = [r[y + dy][xs] for dy in range(-4, 5)]
        if sum(1 for p in col if p == c) < LABEL_MIN_ON:
            continue
        bar = [r[y][xs + i] for i in range(1, 5)]
        top = [r[y - 4][xs + i] for i in range(0, 4)] + [r[y - 3][xs + i] for i in range(0, 4)]
        if bar.count(c) >= 3:
            return "H"
        if bar.count(WHITE) >= 3 and top.count(c) >= 3:
            return "D"
        return "?"
    return "?"


def is_small_star(r, x, y, c, x1):
    """The squat variant: a 5-px horizontal bar with diagonal pixels above
    and below the centre but no full vertical bar."""
    bar = [r[y][x + dx] for dx in range(-2, 3)]
    if sum(1 for p in bar if p == c) < 4 or WHITE in bar:
        return False
    if r[y][x - 3] == c or r[y][x + 3] == c:
        return False
    above = [r[y - 1][x + dx] for dx in (-1, 0, 1)]
    below = [r[y + 1][x + dx] for dx in (-1, 0, 1)]
    if above.count(c) + below.count(c) < 2:
        return False
    if 3 - above.count(WHITE) < 1 or 3 - below.count(WHITE) < 1 or 6 - above.count(WHITE) - below.count(WHITE) < 3:
        return False
    return has_label(r, x, y, c, x1)


def plot_totals(gif):
    """The plot's own "Data Sets =N (Mpts)" line, black text at a fixed spot."""
    png = subprocess.run(["magick", gif, "-crop", "280x20+470+432", "-fill", "white", "+opaque", "black",
                          "-fill", "black", "-opaque", "black", "-scale", "300%", "-bordercolor", "white",
                          "-border", "20", "png:-"], capture_output=True, check=True).stdout
    txt = subprocess.run(["tesseract", "stdin", "-", "--psm", "7"], input=png, capture_output=True).stdout.decode()
    # Tesseract mangles "pts" ("pis", "pit}", "s"); the digits are reliable,
    # except that a 0 may come back as O; anything else unreadable is None.
    m = re.search(r"=\s*([0-9OoIl]+)\s*\(\s*([0-9OoIl]+?)\s*(?:S?p|s\b|\))", txt)
    if not m:
        return None, None, txt.strip()
    fix = lambda t: int(t.replace("O", "0").replace("o", "0").replace("I", "1").replace("l", "1"))
    return fix(m.group(1)), fix(m.group(2)), txt.strip()


def lower_marks(r, x0, x1, z1, rows):
    """Markers in the Experiment/Theory panel, found as in the top panel.
    Excluded: its title block (blue, like one marker colour), the blue
    annotation box, and the magenta
    "5%"/"10%" labels at the right edge."""
    import residual
    y0, y1 = residual.PANEL_TOP, residual.PANEL_BOTTOM
    if z1 == 1:
        marks = find_stars(r, x0, x1, y0, y1)
    else:
        import symdetect
        colours = {tuple(int(row["colour"][i:i + 2], 16) for i in (1, 3, 5)) for row in rows}
        marks = [(hx, hy_ + 4, c, 0) for hx, hy_, c, _ in symdetect.find_H(r, x0, x1, y0, y1, colours)
                 if not on_curve(r, hx, hy_ + 4, c)]
    return [m for m in marks if m[1] > y0 + residual.TITLE_PX and m[0] < x1 - 30
            and not (m[2] == ANNOT_BLUE and m[0] < ANNOT_BOX[0] and m[1] < y0 + ANNOT_BOX[1] - 24)
            and not (m[2] == ANNOT_BLUE and m[0] < x0 + residual.TITLE_BOX[0] and m[1] < y0 + residual.TITLE_BOX[1])]


def lower_points(r, x0, x1, y0, y1, cx, cy, anchor, lower, pts, rows, z1):
    """Add the points seen only in the Experiment/Theory panel to pts, as
    ratio x SRIM curve, and mark top points also seen there.

    The curve is read off the top panel as an intermediate only. Points seen
    in both panels measure how well ratio x curve reproduces a direct reading;
    that scatter, per plot, is the uncertainty given to points seen only below."""
    import residual
    D = cx["x_decade_px"]
    curve = residual.trace_curve(r, x0, x1, y0, y1)
    to_e = lambda xf: 10 ** (cx["x_first_decade_log10"] + (xf - cx["x_first_decade_px"]) / D)
    ypx = abs(cy["y_per_px"])

    def curve_s(xf):
        yy = residual.curve_at(curve, xf)
        return None if yy is None else cy["y_per_px"] * yy + cy["y_at_px0"]

    top = [{"x": p["px"][0] + anchor, "s": p["stopping_ev_cm2_per_1e15"], "colour": p["colour"]} for p in pts]
    bottom = [{"x": m[0] + anchor, "ratio": residual.ratio(m[1]), "colour": "#%02X%02X%02X" % m[2], "m": m}
              for m in lower]
    pairs, extra = residual.merge(top, bottom, curve_s, tol_px=1.5,
                                  tol_s=lambda s: max(0.01 * s, LOWER_PAIR_PX * ypx))
    for p in pts:
        p["panel"] = "top"
    dpx = []
    for i, j in pairs:
        pts[i]["panel"] = "both"
        pts[i]["lower_check_rel"] = round(bottom[j]["s"] / top[i]["s"] - 1, 4)
        dpx.append((bottom[j]["s"] - top[i]["s"]) / ypx)
    if len(dpx) >= 10:
        m = sum(dpx) / len(dpx)
        sd = (sum((v - m) ** 2 for v in dpx) / (len(dpx) - 1)) ** 0.5
    else:
        sd = LOWER_DEFAULT_PX
    for j in extra:
        b = bottom[j]
        mx, my, c, q = b["m"]
        e = to_e(b["x"])
        col = b["colour"]
        pts.append({"energy_kev_per_amu": round(e, 4 - int(math.floor(math.log10(e)))),
                    "stopping_ev_cm2_per_1e15": round(b["s"], 4),
                    "u_energy_rel": round(math.log(10) / D * 0.6, 4),
                    "u_stopping": round(sd * ypx, 4),
                    "colour": col, "px": [mx, my], "match": q,
                    "label": label_letter(r, mx, my, c, x1) if z1 == 1 else "He",
                    "panel": "lower", "ratio": round(b["ratio"], 4),
                    "candidates": [row["key"] or "ref%d" % row["ref"] for row in rows if row["colour"] == col]})
    return {"pairs": len(pairs), "lower_only": len(extra), "pair_sd_px": round(sd, 3)}


def digitize(gif, pdf=None):
    """One per-pair stopping plot -> dict with calibration, legend and points."""
    import legocr, citpdf
    w, h, r = ppm.load(gif)
    x0, x1, y0, y1 = frame(w, h, r)
    vx, hy = gridlines(r, x0, x1, y0, y1)
    cx = calibrate_x(gif, r, x0, x1, y1, vx)
    cy = calibrate_y(gif, r, x0, y0, y1, hy)
    entries = citpdf.parse(pdf) if pdf and os.path.exists(pdf) else []
    byref = {e["ref"]: e for e in entries if e["ref"]}
    rows = [row for row in legocr.read_legend(r, x0, x1, y0, y1, vx, hy, set(byref)) if row["total"] is not None]
    for row in rows:
        e = byref.get(row["ref"])
        row["key"] = e["key"] if e else None
    # The legend box: right of x1-200, down to the last legend row.
    leg_bottom = max((row["top"] + 8 for row in rows), default=y0)
    z1 = int(os.path.basename(gif)[4:6])
    if z1 == 1:
        marks = find_stars(r, x0, x1, y0, y1)
        anchor = ANCHOR_DX
    else:
        import symdetect
        colours = {tuple(int(row["colour"][i:i + 2], 16) for i in (1, 3, 5)) for row in rows}
        # Below the title block; black only where it is not the curve.
        marks = [(hx, hy_ + 4, c, 0) for hx, hy_, c, _ in symdetect.find_H(r, x0, x1, y0, y1, colours)
                 if hy_ > y0 + TITLE_PX and not on_curve(r, hx, hy_ + 4, c)]
        anchor = HE_ANCHOR_DX
    stars = [s for s in marks if not (s[0] >= x1 - 205 and s[1] <= leg_bottom + 2)
             and not (s[2] == ANNOT_BLUE and s[0] < ANNOT_BOX[0] and s[1] < ANNOT_BOX[1])]
    lower = lower_marks(r, x0, x1, z1, rows)
    D = cx["x_decade_px"]
    pts = []
    for (px, py, c, q) in stars:
        col = "#%02X%02X%02X" % c
        e = 10 ** (cx["x_first_decade_log10"] + (px + anchor - cx["x_first_decade_px"]) / D)
        s = cy["y_per_px"] * py + cy["y_at_px0"]
        pts.append({"energy_kev_per_amu": round(e, 4 - int(math.floor(math.log10(e)))),
                    "stopping_ev_cm2_per_1e15": round(s, 4),
                    "u_energy_rel": round(math.log(10) / D * 0.6, 4),  # anchor scatter 0.27 px plus half a pixel
                    "u_stopping": round(abs(cy["y_per_px"]) * 0.5, 4),
                    "colour": col, "px": [px, py], "match": q,
                    "label": label_letter(r, px, py, c, x1) if z1 == 1 else "He",
                    "candidates": [row["key"] or "ref%d" % row["ref"] for row in rows if row["colour"] == col]})
    lower_stats = lower_points(r, x0, x1, y0, y1, cx, cy, anchor, lower, pts, rows, z1)
    counts = defaultdict(int)
    for p in pts:
        counts[p["colour"]] += 1
    check = []
    for col in sorted({row["colour"] for row in rows} | set(counts)):
        check.append({"colour": col, "keys": [row["key"] or "ref%d" % row["ref"] for row in rows if row["colour"] == col],
                      "legend_total": sum(row["total"] for row in rows if row["colour"] == col),
                      "found": counts.get(col, 0)})
    n_sets, n_pts, _ = plot_totals(gif)
    legend_sum = sum(row["total"] for row in rows)
    return {"plot": os.path.basename(gif), "frame": [x0, x1, y0, y1], **cx, **cy,
            "lower_panel": lower_stats,
            "plot_data_sets": n_sets, "plot_points": n_pts, "legend_rows": len(rows), "legend_points": legend_sum,
            "legend_complete": n_sets == len(rows) and n_pts == legend_sum,
            "legend": rows, "count_check": check, "n_points": len(pts), "points": pts}


def main():
    gif = sys.argv[1]
    pdf = sys.argv[2] if len(sys.argv) > 2 else None
    print(json.dumps(digitize(gif, pdf), indent=1))


if __name__ == "__main__":
    main()
