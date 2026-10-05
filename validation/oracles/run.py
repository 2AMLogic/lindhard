#!/usr/bin/env python3
"""Level-2 code-to-code comparison (docs/validation.md, CONTRIBUTING.md
"Oracles").

Runs `lindhard` and every *configured* oracle on the matched problems in
`problems.json`, then writes one summary per (oracle, problem) to
`validation/oracles/summaries/`. Only our own comparison metrics go into a
summary (relative differences, the oracle's name and version, our own
values); the oracle's raw output stays in the gitignored
`validation/oracle-runs/`.

Oracles are third-party programs that the user installs; they are never
vendored, and nothing in CI needs them. Each one is configured by an
environment variable naming its executable:

    RUSTBCA_BIN   RustBCA (GPL: Tier B, run unmodified as an oracle only)
    OPENTRIM_BIN  OpenTRIM (MIT; its SRIM-2013 stopping tables are Tier C and
                  must not be used for these runs or copied anywhere)

An oracle whose variable is unset is skipped with a message, and the script
still exits 0.

Usage:
    validation/oracles/run.py                 # all problems, all configured oracles
    validation/oracles/run.py --problem b_5keV_si
    validation/oracles/run.py --lindhard-only # exercise the lindhard side only
    validation/oracles/run.py --ions 2000     # override the history count
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

SUMMARIES = HERE / "summaries"

# Follow-up that implements the oracle adapters (input writers and output
# parsers from each oracle's published documentation).
ADAPTER_ISSUE = "#50"


class NotConfigured(Exception):
    pass


class Oracle:
    """One external program. Subclasses implement `run`."""

    name = "?"
    env = "?"
    license_note = ""

    def binary(self) -> Path:
        b = os.environ.get(self.env)
        if not b:
            raise NotConfigured(f"{self.name}: not configured (set {self.env}); skipped")
        p = Path(b)
        if not p.exists():
            raise NotConfigured(f"{self.name}: {self.env}={b} does not exist; skipped")
        return p

    def version(self, binary: Path) -> str:
        """`$<NAME>_VERSION` if set (some oracles print no version), else
        the first line of `--version`."""
        override = os.environ.get(self.env.replace("_BIN", "_VERSION"))
        if override:
            return override
        try:
            out = subprocess.run(
                [str(binary), "--version"], capture_output=True, text=True, timeout=30
            )
            line = (out.stdout or out.stderr).strip().splitlines()
            return line[0] if line else "unknown"
        except (OSError, subprocess.TimeoutExpired):
            return "unknown"

    def run(self, problem: dict, workdir: Path, binary: Path) -> dict:
        """Run the oracle unmodified on `problem` in `workdir`; return a dict
        with the keys of `problem['metrics']` and `mismatches` (a list of
        physics choices the oracle could not match)."""
        raise NotImplementedError


class RustBca(Oracle):
    name = "RustBCA"
    env = "RUSTBCA_BIN"
    license_note = "GPL (Tier B): run unmodified; source, tables and fixtures never read or copied"

    def run(self, problem, workdir, binary):
        # The input writer must be written from RustBCA's published
        # documentation (its JOSS paper and user manual), never from its
        # source or example inputs. Not done yet.
        raise NotImplementedError(
            f"{self.name}: input writer / output parser not implemented yet ({ADAPTER_ISSUE})"
        )


class OpenTrim(Oracle):
    name = "OpenTRIM"
    env = "OPENTRIM_BIN"
    license_note = "MIT (Tier A); its SRIM-2013 stopping tables are Tier C and are not used"

    def run(self, problem, workdir, binary):
        raise NotImplementedError(
            f"{self.name}: input writer / output parser not implemented yet ({ADAPTER_ISSUE})"
        )


ORACLES = [RustBca(), OpenTrim()]


def rel(a: float | None, b: float | None) -> float | None:
    if a is None or b is None or b == 0:
        return None
    return a / b - 1.0


def compare(ours: dict, theirs: dict, metrics: list[str]) -> dict:
    """Our summary metrics of a comparison. Relative differences for
    moments and yields, an absolute one for the backscatter fraction, a ratio
    for speed. No raw oracle values."""
    out = {}
    for m in metrics:
        if m == "backscatter":
            if theirs.get(m) is not None:
                out["backscatter_abs_diff"] = ours[m] - theirs[m]
        elif m == "ions_per_s":
            r = rel(ours[m], theirs.get(m))
            if r is not None:
                out["ions_per_s_ratio"] = r + 1.0
        else:
            r = rel(ours[m], theirs.get(m))
            if r is not None:
                out[f"{m.removesuffix('_nm')}_rel_diff"] = r
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", action="append", help="problem id (repeatable)")
    ap.add_argument("--ions", type=int, help="override run.ions")
    ap.add_argument("--lindhard-only", action="store_true", help="run lindhard only and print its metrics")
    args = ap.parse_args()

    spec = json.loads((HERE / "problems.json").read_text())
    problems = [p for p in spec["problems"] if not args.problem or p["id"] in args.problem]
    if not problems:
        lindhard_cli.die(f"no problem matches {args.problem}")
    if args.ions:
        for p in problems:
            p["run"]["ions"] = args.ions

    # Which oracles can run at all? Report before doing any work.
    ready = []
    if not args.lindhard_only:
        for o in ORACLES:
            try:
                b = o.binary()
                ready.append((o, b, o.version(b)))
                print(f"{o.name}: {b} ({ready[-1][2]})")
            except NotConfigured as e:
                print(e)
        if not ready:
            print("no oracle configured; nothing to compare (use --lindhard-only to exercise the lindhard side)")
            return 0

    binary = lindhard_cli.lindhard_binary()
    ours_version = lindhard_cli.version(binary)
    today = _dt.date.today().isoformat()
    SUMMARIES.mkdir(exist_ok=True)
    for p in problems:
        ours = lindhard_cli.run(p, lindhard_cli.RUNS / "lindhard" / p["id"], binary)
        print(f"lindhard {p['id']}: " + ", ".join(f"{m} {ours[m]:.4g}" for m in p["metrics"] if ours[m] is not None))
        for o, b, v in ready:
            work = lindhard_cli.RUNS / o.name.lower() / p["id"]
            work.mkdir(parents=True, exist_ok=True)
            try:
                theirs = o.run(p, work, b)
            except NotImplementedError as e:
                print(f"  {e}; skipped")
                continue
            summary = {
                "format": "lindhard-oracle-summary/1",
                "problem": p["id"],
                "oracle": o.name,
                "oracle_version": v,
                "oracle_license": o.license_note,
                "lindhard_version": ours_version,
                "date": today,
                "ions": {"lindhard": ours["ions"], "oracle": theirs.get("ions")},
                "lindhard": {m: ours[m] for m in p["metrics"]},
                "comparison": compare(ours, theirs, p["metrics"]),
                "mismatches": theirs.get("mismatches", []),
            }
            path = SUMMARIES / f"{o.name.lower()}-{p['id']}.json"
            path.write_text(json.dumps(summary, indent=2) + "\n")
            print(f"  {o.name}: wrote {path.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
