#!/usr/bin/env python3
"""Fold `perf script -F ip,sym` stacks into a small flame-graph SVG.

Standard library only; written for lindhard (MIT), no third-party code.

    perf script -i perf.data -F ip,sym | python3 perf_to_svg.py out.svg [title]

Also prints, to stdout, the top self-time and top inclusive-time frames, which
is what docs/benchmarks.md quotes. Frames narrower than MIN_FRAC of the total
are dropped from the drawing (not from the tables) to keep the file small.
"""
import re
import sys
from collections import defaultdict
from html import escape

MIN_FRAC = 0.004
WIDTH, ROW, PAD = 1200, 16, 24


def clean(sym):
    sym = re.sub(r"::h[0-9a-f]{16}$", "", sym)
    sym = re.sub(r"\+0x[0-9a-f]+$", "", sym)
    return sym


def read_stacks(f):
    stack = []
    for line in f:
        line = line.rstrip("\n")
        if not line.strip():
            if stack:
                yield list(reversed(stack))  # root first
            stack = []
            continue
        parts = line.split(None, 1)
        stack.append(clean(parts[1]) if len(parts) == 2 else "[unknown]")
    if stack:
        yield list(reversed(stack))


def main():
    out = sys.argv[1]
    title = sys.argv[2] if len(sys.argv) > 2 else "flame graph"
    tree = {"n": 0, "kids": {}}
    selfc, incl, total = defaultdict(int), defaultdict(int), 0
    for st in read_stacks(sys.stdin):
        total += 1
        selfc[st[-1]] += 1
        for s in set(st):
            incl[s] += 1
        node = tree
        node["n"] += 1
        for s in st:
            node = node["kids"].setdefault(s, {"n": 0, "kids": {}})
            node["n"] += 1
    if not total:
        sys.exit("no samples on stdin")

    rects, depth_max = [], 0

    def walk(node, name, x, depth):
        nonlocal depth_max
        w = node["n"] / total
        if w < MIN_FRAC:
            return
        depth_max = max(depth_max, depth)
        rects.append((x, depth, w, name, node["n"]))
        cx = x
        for k, v in sorted(node["kids"].items()):
            walk(v, k, cx, depth + 1)
            cx += v["n"] / total

    walk(tree, "all", 0.0, 0)
    height = (depth_max + 1) * ROW + 2 * PAD
    w_px = WIDTH - 20
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
        f'viewBox="0 0 {WIDTH} {height}" font-family="monospace" font-size="11">',
        f'<rect width="100%" height="100%" fill="#f8f8f8"/>',
        f'<text x="10" y="16" font-size="13">{escape(title)} ({total} samples)</text>',
    ]
    for x, d, w, name, n in rects:
        y = height - PAD - (d + 1) * ROW
        px, pw = 10 + x * w_px, w * w_px
        hue = (sum(map(ord, name)) * 7) % 50 + 5
        label = ""
        chars = int(pw / 6.6)
        if chars >= 4:
            t = name if len(name) <= chars else name[: chars - 2] + ".."
            label = f'<text x="{px + 2:.1f}" y="{y + 11}">{escape(t)}</text>'
        svg.append(
            f'<g><title>{escape(name)} ({n} samples, {100 * n / total:.1f}%)</title>'
            f'<rect x="{px:.2f}" y="{y}" width="{max(pw - 0.3, 0.3):.2f}" height="{ROW - 1}" '
            f'fill="hsl({hue},80%,62%)"/>{label}</g>'
        )
    svg.append("</svg>")
    with open(out, "w") as f:
        f.write("\n".join(svg) + "\n")

    def top(counter, label):
        print(f"{label} (of {total} samples)")
        for s, c in sorted(counter.items(), key=lambda kv: -kv[1])[:12]:
            print(f"  {100 * c / total:5.1f}%  {s}")

    top(selfc, "self time")
    top(incl, "inclusive time")


main()
