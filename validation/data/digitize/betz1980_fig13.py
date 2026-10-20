#!/usr/bin/env python3
"""Recompute Betz Fig. 13 numbers from recorded manual marker centres.

Standard library only. No PDF or image required: source hash, render command,
axis calibration, first reading and second-reading centres are in the JSON.
This performs calibration, not automated image/marker detection.
"""

import argparse
import csv
import io
import json
from pathlib import Path


def render_csv(record):
    axes = record["second_read"]["axis_calibration"]
    x0, x1 = axes["x_endpoints_px"]
    c0, c1 = axes["x_values_at_percent"]
    ys = axes["left_y_ticks_px"]
    values = axes["left_y_values_atoms_per_ion"]
    mean_y, mean_v = sum(ys) / len(ys), sum(values) / len(values)
    slope = sum((v - mean_v) * (y - mean_y) for v, y in zip(values, ys)) / sum(
        (v - mean_v) ** 2 for v in values
    )
    intercept = mean_y - slope * mean_v
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow([
        "row_id", "bulk_cu_at_percent", "total_yield_atoms_per_ion",
        "cu_digitization_bound_at_percent", "yield_digitization_bound_atoms_per_ion",
    ])
    for point in record["points"]:
        x, y = point["second_read_marker_px"]
        fraction = (x - x0) / (x1 - x0)
        composition = c0 + fraction * (c1 - c0)
        skew = fraction * axes["right_minus_left_y_shift_px"]
        measured_yield = (y - skew - intercept) / slope
        writer.writerow([
            point["row_id"], f"{composition:.1f}", f"{measured_yield:.2f}",
            point["digitization_bound_cu_at_percent"],
            point["digitization_bound_yield_atoms_per_ion"],
        ])
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check instead of writing CSV")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    record = json.loads((here / "betz1980_fig13.json").read_text())
    target = here.parent / "alloys" / "betz1980_fig13.csv"
    rendered = render_csv(record)
    if args.check:
        if target.read_text() != rendered:
            raise SystemExit("Betz Fig. 13 CSV differs from recorded digitization")
        print("Betz Fig. 13: 22-point CSV matches recorded digitization")
    else:
        target.write_text(rendered)


if __name__ == "__main__":
    main()
