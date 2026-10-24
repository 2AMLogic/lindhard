#!/usr/bin/env python3
"""Recompute Titze et al. Fig. 3 (5e16 ions/cm2) numbers from recorded marker centres.

Standard library only. No PDF or image required: source hash, render command,
axis calibration, first-read and second-read marker centres are in the JSON.
This performs calibration and the two-read cross-check, not image detection.
"""

import argparse
import csv
import io
import json
from pathlib import Path


def calibrate(record):
    (x0, d0), (x1, d1) = record["axis_calibration"]["x_px_at_depth_nm"]
    (y0, c0), (y1, c1) = record["axis_calibration"]["y_px_at_normalized_concentration"]

    def convert(point):
        x, y = point
        return d0 + (x - x0) / (x1 - x0) * (d1 - d0), c0 + (y - y0) / (y1 - y0) * (c1 - c0)

    return convert


def render_csv(record):
    convert = calibrate(record)
    bounds = record["bounds"]
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow([
        "row_id", "depth_nm", "normalized_sb_concentration",
        "depth_bound_nm", "concentration_bound",
    ])
    points = sorted(record["first_read"]["marker_px"])
    for index, point in enumerate(points, start=1):
        depth, conc = convert(point)
        writer.writerow([
            f"fig3-5e16-p{index:02d}", f"{depth:.2f}", f"{conc:.3f}",
            bounds["depth_nm"], bounds["normalized_concentration"],
        ])
    return output.getvalue()


def cross_check(record):
    """Each second-read marker against the nearest first-read marker in depth."""
    convert = calibrate(record)
    bounds = record["bounds"]
    first = [convert(p) for p in record["first_read"]["marker_px"]]
    rows = []
    for point in record["second_read"]["marker_px"]:
        depth, conc = convert(point)
        near = min(first, key=lambda f: abs(f[0] - depth))
        ok = (abs(near[0] - depth) <= 2 * bounds["depth_nm"]
              and abs(near[1] - conc) <= 2 * bounds["normalized_concentration"])
        rows.append((depth, conc, near[0], near[1], ok))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check instead of writing CSV")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    record = json.loads((here / "titze2021_fig3.json").read_text())
    target = here.parent / "implant" / "titze2021_fig3_5e16.csv"
    rendered = render_csv(record)
    if args.check:
        if target.read_text() != rendered:
            raise SystemExit("Titze Fig. 3 CSV differs from recorded digitization")
        rows = cross_check(record)
        bad = [r for r in rows if not r[4]]
        if bad:
            raise SystemExit(f"second read disagrees with first read at {len(bad)} markers")
        print(f"Titze Fig. 3: {len(record['first_read']['marker_px'])}-point CSV matches; "
              f"{len(rows)} second-read markers agree with the first read within the summed bounds")
    else:
        target.write_text(rendered)


if __name__ == "__main__":
    main()
