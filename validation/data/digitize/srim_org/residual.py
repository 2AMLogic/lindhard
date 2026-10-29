"""The lower panel of each plot: the same measured points divided by SRIM's
curve (Experiment/Theory). Markers that sit in a stack in the top panel are
spread apart here by their few-percent differences, so a point hidden on top
can still be read below, and recovered as ratio x curve.

SRIM's curve is read here only as that intermediate. It is never written out,
and no point is made from it alone (operator ruling 2026-10-28, #486).
"""
WHITE, BLACK = (255, 255, 255), (0, 0, 0)

# The lower panel is the same box on every plot: rows 549..1022, the x axis of
# the top panel. Ratio 1.0 is the black line at row 786; the magenta guide
# lines at 0.90, 0.95, 1.05 and 1.10 are 39.5 px apart.
PANEL_TOP, PANEL_BOTTOM = 549, 1022
ONE_ROW, ROWS_PER_005 = 786, 39.5
# Below the panel's top edge, the title block and "STOP-2011" (as above).
TITLE_PX = 60
# The blue "Ion = ... / Target = ..." title reaches lower than that at the
# left, in the same blue as one marker colour: (width from the frame, height).
TITLE_BOX = (460, 90)


def ratio(y):
    return 1 + (ONE_ROW - y) * 0.05 / ROWS_PER_005


def ratio_per_px():
    return 0.05 / ROWS_PER_005


def black_runs(r, x, y0, y1):
    """(start, end) of each vertical run of black in column x, inside the
    panel (frame rows excluded)."""
    runs, start = [], None
    for y in range(y0 + 2, y1 - 1):
        if r[y][x] == BLACK:
            if start is None:
                start = y
        elif start is not None:
            runs.append((start, y - 1))
            start = None
    if start is not None:
        runs.append((start, y1 - 2))
    return runs


def trace_curve(r, x0, x1, y0, y1):
    """y of the curve's centre line in each column of the top panel, {x: y}.

    The curve is the one black line that crosses the whole panel; black text
    ("STOP-2011", "Data Sets") and black data markers also make runs. A path
    through the columns is chosen that minimises the jumps between columns,
    with a penalty per column skipped, so short text runs cannot capture it.
    Columns where the curve is covered by other ink are interpolated."""
    cols = list(range(x0 + 2, x1 - 1))
    cand = {x: [((a + b) / 2, b - a + 1) for a, b in black_runs(r, x, y0, y1) if b - a + 1 >= 2] for x in cols}
    SKIP = 6.0
    # best[x][i] = (cost, prev); dynamic programming over columns.
    best = {}
    prev_states = []  # (x, i, y, cost)
    for x in cols:
        cur = {}
        for i, (y, n) in enumerate(cand[x]):
            c0 = None
            for px, pi, py, pc in prev_states:
                gap = x - px - 1
                c = pc + abs(y - py) / (gap + 1) + SKIP * gap
                if c0 is None or c < c0[0]:
                    c0 = (c, (px, pi))
            if c0 is None:
                c0 = (SKIP * (x - cols[0]), None)
            cur[i] = c0
        best[x] = cur
        # keep states of the last few columns so a covered column can be skipped
        prev_states = [s for s in prev_states if x - s[0] < 25] + [(x, i, cand[x][i][0], cur[i][0]) for i in cur]
    # End state: cheapest among the last columns, counting the skip to the edge.
    ends = [(c + SKIP * (cols[-1] - x), x, i) for x in cols[-25:] for i, (c, _) in best[x].items()]
    if not ends:
        return {}
    _, x, i = min(ends)
    path = {}
    while True:
        path[x] = cand[x][i][0]
        prev = best[x][i][1]
        if prev is None:
            break
        x, i = prev
    xs = sorted(path)
    out = dict(path)
    for a, b in zip(xs, xs[1:]):
        for x in range(a + 1, b):
            out[x] = path[a] + (path[b] - path[a]) * (x - a) / (b - a)
    return out


def curve_at(curve, xf):
    """Curve row at a fractional column (linear between columns)."""
    a = int(xf)
    if a not in curve or a + 1 not in curve:
        return curve.get(round(xf))
    return curve[a] + (curve[a + 1] - curve[a]) * (xf - a)


def merge(top, bottom, curve_s, tol_px, tol_s, loose_rel=0.04):
    """Pair lower-panel points with top-panel points of the same colour.

    top: [{x, s, colour, ...}] with x the anchor column (fractional) and s the
    stopping read on top. bottom: [{x, ratio, colour, ...}]. curve_s(x) gives
    the curve's stopping at a column. Each lower point's stopping is
    ratio x curve; it pairs with the nearest unpaired top point of its colour
    within tol_px in x and tol_s(top stopping) in stopping. A second pass pairs
    what is left on both sides within loose_rel in stopping: where markers
    cover the curve its trace can be off by a few pixels, and an unpaired top
    and lower marker of one colour at one energy are the same point. Returns
    (pairs, extra): pairs of (top index, bottom index) and the indices of
    lower points with no top partner (points hidden on top)."""
    for b in bottom:
        c = curve_s(b["x"])
        if c is not None:
            b["s"] = b["ratio"] * c
    used_t, used_b, pairs = set(), set(), []
    for tol_of in (tol_s, lambda s: loose_rel * s):
        cands = []
        for j, b in enumerate(bottom):
            if "s" not in b or j in used_b:
                continue
            for i, t in enumerate(top):
                if i in used_t or t["colour"] != b["colour"] or abs(t["x"] - b["x"]) > tol_px:
                    continue
                ds, tol = abs(b["s"] - t["s"]), tol_of(t["s"])
                if ds <= tol:
                    cands.append((abs(t["x"] - b["x"]) / tol_px + ds / tol, i, j))
        cands.sort()
        for _, i, j in cands:
            if i in used_t or j in used_b:
                continue
            used_t.add(i)
            used_b.add(j)
            pairs.append((i, j))
    extra = [j for j, b in enumerate(bottom) if j not in used_b and "s" in b]
    return pairs, extra
