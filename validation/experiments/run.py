#!/usr/bin/env python3
"""Level-3 comparison with published range measurements
(docs/validation.md section 3; data and schema: validation/data/README.md).

For every `validation/data/ranges/*.json`, runs `lindhard` on the same ion,
energy, tilt and target with the default physics below, and writes
`validation/experiments/results.json` with the measured and computed moments
and their difference in units of the measurement uncertainty.

Usage:
    validation/experiments/run.py [--ions N]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

DATA = HERE.parent / "data" / "ranges"
RESULTS = HERE / "results.json"

# The physics every dataset is run with: the M0 defaults. E_d does not
# affect primary-only range runs. Recorded in results.json.
PHYSICS = {
    "potential": "zbl",
    "stopping": "lindhard-scharff",
    "free_path": "constant",
    "primary_cutoff_ev": 5.0,
    "recoil_cutoff_ev": 2.0,
    "follow_recoils": False,
}


def problem(d: dict, ions: int) -> dict:
    beam = {"ion": d["ion"], "energy_ev": d["energy_ev"], "tilt_deg": d.get("tilt_deg", 0.0)}
    if d.get("mass_amu"):
        beam["mass_amu"] = d["mass_amu"]
    rp = d["measured"]["rp_nm"]
    return {
        "id": d["id"],
        "beam": beam,
        "target": {"substrate": d["target"]},
        "physics": {**PHYSICS, "energies": {d["target"]: {"e_d_ev": 15.0}}},
        "run": {"ions": ions, "seed": 1},
        "tally": {"depth_bin_nm": max(rp / 100.0, 0.1), "depth_bins": 1000},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ions", type=int, default=20000)
    args = ap.parse_args()

    files = sorted(DATA.glob("*.json")) if DATA.is_dir() else []
    if not files:
        print(f"no datasets in {DATA.relative_to(lindhard_cli.REPO)}; see validation/data/README.md")
        return 0
    binary = lindhard_cli.lindhard_binary()
    rows = []
    for f in files:
        d = json.loads(f.read_text())
        m = d["measured"]
        ours = lindhard_cli.run(problem(d, args.ions), lindhard_cli.RUNS / "experiments" / d["id"], binary)
        row = {
            "id": d["id"],
            "case": f"{d['ion']} {d['energy_ev'] / 1e3:g} keV -> {d['target']} ({m.get('method', '?')})",
            "citation": d["citation"],
            "rp_measured_nm": m["rp_nm"],
            "rp_unc_nm": m["rp_unc_nm"],
            "rp_lindhard_nm": ours["rp_nm"],
            "rp_rel_diff": ours["rp_nm"] / m["rp_nm"] - 1.0,
            "rp_diff_sigma": (ours["rp_nm"] - m["rp_nm"]) / m["rp_unc_nm"],
            "drp_measured_nm": m["drp_nm"],
            "drp_unc_nm": m["drp_unc_nm"],
            "drp_lindhard_nm": ours["drp_nm"],
            "drp_rel_diff": ours["drp_nm"] / m["drp_nm"] - 1.0,
            "drp_diff_sigma": (ours["drp_nm"] - m["drp_nm"]) / m["drp_unc_nm"],
        }
        rows.append(row)
        print(
            f"{row['case']}: Rp {ours['rp_nm']:.1f} vs {m['rp_nm']:.1f} nm ({100 * row['rp_rel_diff']:+.1f} %), "
            f"dRp {ours['drp_nm']:.1f} vs {m['drp_nm']:.1f} nm ({100 * row['drp_rel_diff']:+.1f} %)"
        )
    RESULTS.write_text(
        json.dumps(
            {
                "format": "lindhard-experiment-results/1",
                "lindhard_version": lindhard_cli.version(binary),
                "ions": args.ions,
                "physics": PHYSICS,
                "datasets": rows,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {RESULTS.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
