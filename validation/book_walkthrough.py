#!/usr/bin/env python3
"""Run the user-guide walkthrough (5 keV B into Si) and write, or check, the
command output that book/src/guide/first-run.md quotes.

The guide includes the files under book/src/guide/walkthrough/ with mdBook's
{{#include}}, so every number it quotes comes from an actual run of the
`lindhard` binary on examples/b_5keV_si.toml, not from memory. This script
regenerates those files; with --check it regenerates them in memory and
exits 1 if any differs from the committed copy (CI).

One run of the example (10 000 ions; a few seconds in a release build). The
results do not depend on the thread count, so --threads only sets the speed.

Usage: validation/book_walkthrough.py --lindhard PATH [--threads N] [--check]

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "book" / "src" / "guide" / "walkthrough"
EXAMPLE = "examples/b_5keV_si.toml"
RUN_OUT = "out/b_5keV_si"
CSV_HEAD_LINES = 8


def run(cmd: list[str], cwd: Path) -> str:
    """The command's terminal output: stdout and stderr interleaved, as a
    user sees it (`lindhard run` reports its one-line summary on stderr)."""
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        raise SystemExit(f"book_walkthrough: {' '.join(cmd)} exited {proc.returncode}")
    return proc.stdout


def excerpt(obj: object) -> str:
    return json.dumps(obj, indent=2) + "\n"


def generate(lindhard: Path, threads: int | None) -> dict[str, str]:
    files: dict[str, str] = {}
    with tempfile.TemporaryDirectory(prefix="lindhard-walkthrough-") as tmp:
        work = Path(tmp)
        # Run from a scratch directory that sees the repository's examples/
        # under the same relative path, so the commands and the paths they
        # print are exactly the ones the guide shows.
        os.symlink(ROOT / "examples", work / "examples")
        files["check.txt"] = run([str(lindhard), "check", EXAMPLE], work)
        cmd = [str(lindhard), "run", EXAMPLE, "--out", RUN_OUT]
        if threads is not None:
            cmd += ["--threads", str(threads)]
        files["run.txt"] = run(cmd, work)
        out = work / RUN_OUT
        files["ls.txt"] = "".join(f"{p.name}\n" for p in sorted(out.iterdir()))
        s = json.loads((out / "summary.json").read_text(encoding="utf-8"))
        r = s["results"]
        files["primaries.json"] = excerpt({"primaries": r["primaries"]})
        files["range_depth.json"] = excerpt(
            {
                "depth": r["range"]["depth"],
                "pearson_iv": r["range"]["pearson_iv"],
                "pearson_iv_error": r["range"].get("pearson_iv_error"),
            }
        )
        files["dual_pearson.json"] = excerpt(
            {
                "head_fraction": r["range"]["dual_pearson"]["head_fraction"],
                "head": {
                    k: r["range"]["dual_pearson"]["head"][k]
                    for k in ("mean_nm", "std_dev_nm")
                },
                "tail": {
                    k: r["range"]["dual_pearson"]["tail"][k]
                    for k in ("mean_nm", "std_dev_nm")
                },
                "chi_square": r["range"]["dual_pearson"]["chi_square"],
            }
        )
        files["damage_per_ion.json"] = excerpt({"per_ion": r["damage"]["per_ion"]})
        files["escapes.json"] = excerpt(
            {
                "sputtering": {"yield_per_ion": r["sputtering"]["yield_per_ion"]},
                "escapes": {
                    k: r["escapes"][k]
                    for k in (
                        "backscatter_coefficient",
                        "transmission_coefficient",
                        "energy_reflection_coefficient",
                    )
                },
            }
        )
        files["energy_budget.json"] = excerpt(
            {"energy_budget_ev_per_ion": r["energy_budget_ev_per_ion"]}
        )
        files["models.json"] = excerpt(
            {
                "models": [
                    {"role": m["role"], "name": m["name"]}
                    for m in s["physics"]["models"]
                ]
            }
        )
        lines = (out / "depth_profile.csv").read_text(encoding="utf-8").splitlines(True)
        files["depth_profile_head.csv"] = "".join(lines[:CSV_HEAD_LINES])
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--lindhard", required=True, type=Path, help="path to the lindhard binary"
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=None,
        help="worker threads (does not change results)",
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="compare with the committed files, write nothing",
    )
    args = ap.parse_args(argv)
    files = generate(args.lindhard.resolve(), args.threads)
    if args.check:
        stale = [
            name
            for name, text in files.items()
            if not (OUT_DIR / name).is_file()
            or (OUT_DIR / name).read_text(encoding="utf-8") != text
        ]
        extra = (
            sorted(
                p.name for p in OUT_DIR.glob("*") if p.is_file() and p.name not in files
            )
            if OUT_DIR.is_dir()
            else []
        )
        if stale or extra:
            for name in stale:
                print(
                    f"book_walkthrough: {OUT_DIR.relative_to(ROOT)}/{name} differs from the run",
                    file=sys.stderr,
                )
            for name in extra:
                print(
                    f"book_walkthrough: {OUT_DIR.relative_to(ROOT)}/{name} is not generated",
                    file=sys.stderr,
                )
            print(
                "Regenerate with: cargo build --release -p lindhard-cli && "
                "python3 validation/book_walkthrough.py --lindhard target/release/lindhard",
                file=sys.stderr,
            )
            return 1
        print(f"walkthrough output: {len(files)} files match the run")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (OUT_DIR / name).write_text(text, encoding="utf-8")
    print(f"wrote {len(files)} files to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
