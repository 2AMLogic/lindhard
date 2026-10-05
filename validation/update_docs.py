#!/usr/bin/env python3
"""Rewrite the results tables in docs/validation.md between their marker
comments, so a rerun is idempotent:

    <!-- validation:level1:begin --> ... <!-- validation:level1:end -->
    <!-- validation:level2:begin --> ... <!-- validation:level2:end -->
    <!-- validation:level3:begin --> ... <!-- validation:level3:end -->

Level 1 comes from the Markdown the harness writes (LINDHARD_VALIDATION_OUT),
level 2 from the committed oracle summaries, level 3 from
validation/experiments/results.json. Called by validation/run.sh.

Usage: validation/update_docs.py [--level1 FILE]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "validation.md"
SUMMARIES = ROOT / "validation" / "oracles" / "summaries"
RESULTS = ROOT / "validation" / "experiments" / "results.json"


def splice(text: str, name: str, body: str) -> str:
    begin, end = f"<!-- validation:{name}:begin -->", f"<!-- validation:{name}:end -->"
    pat = re.compile(re.escape(begin) + r".*?" + re.escape(end), re.S)
    if not pat.search(text):
        sys.exit(f"error: markers for {name} not found in {DOC}")
    return pat.sub(lambda _: f"{begin}\n{body.strip()}\n{end}", text)


def pct(x) -> str:
    return "-" if x is None else f"{100 * x:+.1f} %"


def level2() -> str:
    files = sorted(SUMMARIES.glob("*.json")) if SUMMARIES.is_dir() else []
    if not files:
        return "_No oracle summaries committed yet; see \"Running the oracles\" above._"
    lines = [
        "| Problem | Oracle (version) | Rp | ΔRp | Backscatter (abs.) | Sputter yield | Speed ratio | Date |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for f in files:
        s = json.loads(f.read_text())
        c = s["comparison"]
        bs = c.get("backscatter_abs_diff")
        sp = c.get("ions_per_s_ratio")
        lines.append(
            f"| `{s['problem']}` | {s['oracle']} ({s['oracle_version']}) | {pct(c.get('rp_rel_diff'))} "
            f"| {pct(c.get('drp_rel_diff'))} | {'-' if bs is None else f'{bs:+.4f}'} "
            f"| {pct(c.get('sputter_yield_rel_diff'))} | {'-' if sp is None else f'{sp:.2f}x'} | {s['date']} |"
        )
    lines.append("")
    lines.append("Differences are lindhard relative to the oracle; speed ratio > 1 means lindhard is faster.")
    return "\n".join(lines)


def level3() -> str:
    if not RESULTS.exists():
        return "_No experimental datasets yet; see validation/data/README.md._"
    r = json.loads(RESULTS.read_text())
    lines = [
        f"{r['lindhard_version']}, {r['ions']} ions per case, physics {r['physics']['potential']} + {r['physics']['stopping']}.",
        "",
        "| Case | Rp measured (nm) | Rp lindhard (nm) | Diff. | Diff. / unc. | ΔRp measured (nm) | ΔRp lindhard (nm) | Diff. | Source |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for d in r["datasets"]:
        lines.append(
            f"| {d['case']} | {d['rp_measured_nm']:.1f} ± {d['rp_unc_nm']:.1f} | {d['rp_lindhard_nm']:.1f} "
            f"| {pct(d['rp_rel_diff'])} | {d['rp_diff_sigma']:+.1f} | {d['drp_measured_nm']:.1f} ± {d['drp_unc_nm']:.1f} "
            f"| {d['drp_lindhard_nm']:.1f} | {pct(d['drp_rel_diff'])} | {d['citation']} |"
        )
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level1", type=Path, help="Markdown table written by the level-1 harness")
    args = ap.parse_args()
    text = DOC.read_text()
    if args.level1:
        text = splice(text, "level1", args.level1.read_text())
    text = splice(text, "level2", level2())
    text = splice(text, "level3", level3())
    DOC.write_text(text)
    print(f"updated {DOC.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
