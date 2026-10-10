#!/usr/bin/env python3
"""Level-3 comparison of the electron backscatter coefficient eta(E, Z) with
published measurements (#148; docs/validation.md, "Backscatter coefficient";
data and schema: validation/data/README.md).

1. Checks every `validation/data/backscatter/*.json` (one measured set per
   file) and that `docs/data-provenance.md` has a row naming it; any failure
   stops the run (exit 2). `--check` does only this.
2. Groups the measured points per target and comparison energy: a set
   contributes its point nearest to the energy if that point is within 2 %
   of it (the grouping tolerance of the sputter comparison), at most one
   point per set; the group keeps the median, min and max over the sets.
   Sets are never averaged at ingest; the median is taken here.
3. For every target with a committed input
   (`validation/experiments/backscatter/eta_<target>.toml`) runs `lindhard`
   through the CLI's `[electron]` input at every comparison energy, changing
   only `electron.beam.energy_ev` and `run.histories`; eta and its binomial
   standard error sqrt(eta (1 - eta) / N) come from `electron_summary.json`
   (with secondaries off each primary escapes at most once).
4. Sensitivity: at SENSITIVITY_ENERGIES_KEV the same input is rerun with the
   two elastic corrections switched off one at a time and together.

5. Fast-secondary sensitivity (#148 part 1): Al, Cu, Au and Si rerun with Kieft-Bosch
   secondaries and the step barrier, band inputs from se_yield.py; `--secondaries-only`
   runs just these and merges them into the committed results, after rerunning
   the baseline at the same points and checking it reproduces bit for bit
   (`--targets Si` limits it to some targets, `--verify-baseline-only` does only
   the reproduction check, writing nothing).

A target without a committed input is not run; its measured groups are
still reported. All five targets have one (Si since #169).

Writes `validation/experiments/backscatter_results.json`, which
`validation/update_docs.py` turns into the tables in docs/validation.md.

Usage:
    validation/experiments/backscatter.py [--histories N] [--threads N] [--reuse]
    validation/experiments/backscatter.py --secondaries-only [--targets Si,Au] [--threads N] [--reuse]
    validation/experiments/backscatter.py --verify-baseline-only [--targets Al,Cu] [--threads N] [--reuse]
    validation/experiments/backscatter.py --check
"""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402
import se_yield  # noqa: E402  (BAND, METAL_BAND_PROVENANCE: the delta(E) band inputs, read only)

DATA = HERE.parent / "data" / "backscatter"
INPUTS = HERE / "backscatter"
RESULTS = HERE / "backscatter_results.json"
PROVENANCE = lindhard_cli.REPO / "docs" / "data-provenance.md"

TARGETS = ["C", "Al", "Si", "Cu", "Au"]
# Comparison energies, keV: the energies most measured sets share (1-30 keV,
# the issue's range).
ENERGIES_KEV = [1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 30.0]
GROUP_TOLERANCE = 0.02
# The issue's initial tolerance (#148): at E >= 5 keV, |eta_sim - median| <= 0.05
# absolute for every element; below 5 keV reported without pass/fail.
# Loosening it needs operator sign-off on the issue.
PASS_MIN_KEV = 5.0
TOLERANCE = 0.05
SENSITIVITY_ENERGIES_KEV = [1.0, 10.0]
# (id, label, exchange, correlation-polarization). The first is the baseline
# of the main table (both corrections on, as in the committed inputs).
VARIANTS = [
    ("both", "exchange + polarization (baseline)", True, True),
    ("exchange", "exchange only", True, False),
    ("polarization", "polarization only", False, True),
    ("none", "no corrections", False, False),
]

# Fast-secondary sensitivity (#148, part 1): the baseline input rerun with
# `secondaries = "kieft-bosch"` and a step barrier at the mid work function,
# with exactly the band inputs of the delta(E) runs (se_yield.BAND and
# se_yield.BAND_PROVENANCE): free-electron metal for Al, Cu and Au; for Si the
# insulator band (gap and affinity: BAND_DEFAULTS, #115; valence-band width:
# Chelikowsky and Cohen 1974, added by #149), whose barrier is the electron
# affinity. C has no band data, so it is a stated gap, not run. Everything else
# (elastic model, optical ELF, tables, 50 eV band-bottom cutoff, seed, primaries)
# is the baseline's.
SECONDARY_TARGETS = ["Al", "Cu", "Si", "Au"]
SECONDARY_ENERGIES_KEV = [1.0, 5.0, 10.0, 30.0]
SECONDARY_GAPS = {
    "C": "no band data committed",
}

# --- dataset checks ---------------------------------------------------------

REQUIRED_TEXT = [
    "id", "target", "target_state", "original_reference", "compilation",
    "compilation_set", "extraction", "terms", "added",
]


def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def check_dataset(path: Path, d: dict, provenance: str) -> list[str]:
    """Every problem with one backscatter dataset, as messages; empty if it passes."""
    errs = []
    if d.get("kind") != "backscatter_coefficient":
        errs.append("`kind` must be \"backscatter_coefficient\"")
    for k in REQUIRED_TEXT:
        if not isinstance(d.get(k), str) or not d[k].strip():
            errs.append(f"`{k}` missing or empty")
    if d.get("id") != path.stem:
        errs.append(f"`id` {d.get('id')!r} differs from the file name {path.stem!r}")
    if not any(isinstance(d.get(k), str) and d[k].strip() for k in ("original_doi", "url")):
        errs.append("needs an `original_doi` or a `url` for the source that was read")
    if not (isinstance(d.get("z"), int) and d["z"] > 0):
        errs.append("`z` must be a positive integer")
    if not _num(d.get("incidence_deg")):
        errs.append("`incidence_deg` must be a number (normal incidence = 0)")
    pts = d.get("points")
    if not isinstance(pts, list) or not pts:
        errs.append("`points` missing or empty")
    else:
        for i, p in enumerate(pts):
            if not (_num(p.get("energy_ev")) and p["energy_ev"] > 0):
                errs.append(f"points[{i}].energy_ev must be a positive number")
            if not (_num(p.get("eta")) and 0 < p["eta"] < 1):
                errs.append(f"points[{i}].eta must be a number in (0, 1)")
            if not (_num(p.get("eta_unc_abs")) and p["eta_unc_abs"] > 0):
                errs.append(f"points[{i}].eta_unc_abs must be a positive number")
    if isinstance(d.get("id"), str) and d["id"] not in provenance:
        errs.append(f"no row in docs/data-provenance.md names `{d['id']}`")
    return errs


def load_datasets(directory: Path = DATA) -> list[dict]:
    """All backscatter datasets, checked. Exits 2 on any failure."""
    provenance = PROVENANCE.read_text()
    out, failed = [], False
    for f in sorted(directory.glob("eta_*.json")):
        rel = f.relative_to(lindhard_cli.REPO) if f.is_relative_to(lindhard_cli.REPO) else f
        try:
            d = json.loads(f.read_text())
        except json.JSONDecodeError as e:
            print(f"error: {rel}: not valid JSON: {e}", file=sys.stderr)
            failed = True
            continue
        errs = check_dataset(f, d, provenance)
        for e in errs:
            print(f"error: {rel}: {e}", file=sys.stderr)
        failed |= bool(errs)
        out.append(d)
    if failed:
        print("error: backscatter dataset check failed (rules: validation/data/README.md)", file=sys.stderr)
        sys.exit(2)
    return out


# --- grouping ---------------------------------------------------------------


def measured_group(datasets: list[dict], target: str, e_kev: float) -> dict:
    """The measured points of `target` within GROUP_TOLERANCE of `e_kev`: at
    most one per set (its nearest point), normal incidence only."""
    e_ev = 1000.0 * e_kev
    members = []
    for d in datasets:
        if d["target"] != target or d["incidence_deg"] != 0.0:
            continue
        near = [p for p in d["points"] if abs(p["energy_ev"] / e_ev - 1.0) <= GROUP_TOLERANCE]
        if near:
            p = min(near, key=lambda p: abs(p["energy_ev"] - e_ev))
            members.append({"id": d["id"], "energy_ev": p["energy_ev"], "eta": p["eta"]})
    etas = [m["eta"] for m in members]
    return {
        "energy_kev": e_kev,
        "sets": len(members),
        "median": statistics.median(etas) if etas else None,
        "min": min(etas) if etas else None,
        "max": max(etas) if etas else None,
        "members": sorted(members, key=lambda m: m["id"]),
    }


# --- runs -------------------------------------------------------------------


def _replace_once(text: str, pattern: str, repl) -> str:
    new, n = re.subn(pattern, repl, text, flags=re.M)
    if n != 1:
        raise SystemExit(f"error: expected one match of {pattern!r} in the committed input, found {n}")
    return new


def variant_input(base: str, base_dir: Path, e_kev: float, histories: int,
                  exchange: bool, polarization: bool) -> str:
    """The committed input with the energy, history count and corrections set,
    and its data paths made absolute (the copy runs in validation/oracle-runs/)."""
    t = _replace_once(base, r"^energy_ev = .*$", f"energy_ev = {1000.0 * e_kev!r}")
    t = _replace_once(t, r"^histories = .*$", f"histories = {histories}")
    t = _replace_once(t, r"^exchange = .*$", f"exchange = {'true' if exchange else 'false'}")
    t, n = re.subn(
        r'^optical_elf = "(.*)"$',
        lambda m: f'optical_elf = {json.dumps(str((base_dir / m.group(1)).resolve()))}',
        t, flags=re.M,
    )
    if n < 1:
        raise SystemExit("error: the committed input names no optical_elf")
    if not polarization:
        # Drop the [electron.elastic.correlation_polarization] table: its header
        # and every line up to the next table.
        t, n = re.subn(r"^\[electron\.elastic\.correlation_polarization\]\n(?:(?!\[).*\n)*", "", t, flags=re.M)
        if n != 1:
            raise SystemExit("error: the committed input has no correlation_polarization table")
    return t


def secondary_input(base: str, base_dir: Path, target: str, e_kev: float, histories: int) -> str:
    """The baseline input of `target` at `e_kev` with Kieft-Bosch secondaries on,
    the step barrier (metals: the mid work function; Si: its electron affinity) and
    the delta(E) band inputs (se_yield.band_line)."""
    t = variant_input(base, base_dir, e_kev, histories, True, True)
    t = _replace_once(t, r'^secondaries = ".*"$', 'secondaries = "kieft-bosch"')
    t = _replace_once(t, r'^boundary = ".*"$', 'boundary = "step-barrier"')
    line = se_yield.band_line(target, "mid")
    return _replace_once(t, r"^(optical_elf = .*)$", lambda m: m.group(1) + "\n" + line)


def reusable(workdir: Path, text: str, version: str) -> bool:
    """True if `workdir` holds a finished run of exactly `text` by the binary
    whose `--version` is `version` ("lindhard X.Y.Z (<git describe>)"): the
    stored input is byte-identical, its summary exists and names the same
    git describe. Results do not depend on the thread count (a tested
    invariant), so such a run is the run."""
    inp, summary = workdir / "input.toml", workdir / "out" / "electron_summary.json"
    if not (inp.is_file() and summary.is_file() and inp.read_text() == text):
        return False
    try:
        described = json.loads(summary.read_text())["software"]["git_describe"]
    except (json.JSONDecodeError, KeyError, TypeError):
        return False
    return isinstance(described, str) and bool(described) and version.endswith(f"({described})")


def run_one(binary: Path, text: str, workdir: Path, threads: int | None,
            reuse_version: str | None = None) -> dict:
    out = workdir / "out"
    if reuse_version is not None and reusable(workdir, text, reuse_version):
        print(f"  reusing {workdir.name}", file=sys.stderr)
    else:
        workdir.mkdir(parents=True, exist_ok=True)
        inp = workdir / "input.toml"
        inp.write_text(text)
        cmd = [str(binary), "run", str(inp), "--out", str(out)]
        if threads:
            cmd += ["--threads", str(threads)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            lindhard_cli.die(f"lindhard failed on {inp}:\n{proc.stderr}")
    s = json.loads((out / "electron_summary.json").read_text())
    r = s["results"]
    n = r["histories"]
    eta = r["yields"]["backscatter_eta"]
    fast = r["front"]["fast"]
    elastic = s["physics"]["transport"]["layers"][0]["elastic_model"]
    return {
        "version": s["software"]["git_describe"],
        "histories": n,
        "eta": eta,
        "eta_se": math.sqrt(eta * (1.0 - eta) / n),
        "secondary_delta": r["yields"]["secondary_delta"],
        "fast_count": fast["count"],
        "relative_imbalance": r["budget"]["relative_imbalance"],
        "elastic_model": elastic,
        "table_build_s": s["run"]["table_build_s"],
    }


def run_secondaries(binary: Path, entry: dict, histories: int, threads: int | None, reuse,
                    verify_baseline: bool, baseline_only: bool = False) -> None:
    """Fill entry["secondary_sensitivity"] for a target with band inputs.

    With `verify_baseline` (the `--secondaries-only` mode, which keeps the
    committed baseline) the baseline run at each energy is redone first with
    this binary and must reproduce the committed eta and delta bit for bit;
    otherwise the comparison would be against a baseline this binary does not
    give, and the script stops."""
    t = entry["target"]
    inp = INPUTS / f"eta_{t.lower()}.toml"
    base = inp.read_text()
    out = {}
    for e in SECONDARY_ENERGIES_KEV:
        if verify_baseline:
            print(f"{t} {e:g} keV baseline (reproduction check)", file=sys.stderr)
            text = variant_input(base, inp.parent, e, histories, True, True)
            again = run_one(binary, text, lindhard_cli.RUNS / "backscatter" / f"{t}_{e:g}keV_both", threads, reuse)
            committed = entry["runs"][f"{e:g}"]
            for key in ("histories", "eta", "secondary_delta"):
                if again[key] != committed[key]:
                    lindhard_cli.die(f"{t} {e:g} keV: baseline {key} {again[key]!r} does not reproduce the "
                                     f"committed {committed[key]!r} ({committed['version']}); regenerate the "
                                     "whole table instead")
        if baseline_only:
            continue
        print(f"{t} {e:g} keV secondaries", file=sys.stderr)
        text = secondary_input(base, inp.parent, t, e, histories)
        run = run_one(binary, text, lindhard_cli.RUNS / "backscatter" / f"{t}_{e:g}keV_secondaries", threads, reuse)
        # With secondaries on a primary can yield several fast electrons, so the
        # binomial error of the baseline does not apply. The summary has no
        # per-primary multiplicity, so this is the Poisson estimate sqrt(n)/N:
        # it exceeds the binomial error if no primary yields more than one, and
        # can understate the error only to the extent that some yield several.
        run["eta_se"] = math.sqrt(run["fast_count"]) / run["histories"]
        run["eta_se_kind"] = "poisson"
        out[f"{e:g}"] = run
    if not baseline_only:
        entry["secondary_sensitivity"] = out


def secondary_meta(binary: Path) -> dict:
    return {
        "lindhard": lindhard_cli.version(binary),
        "energies_kev": SECONDARY_ENERGIES_KEV,
        "secondaries": "kieft-bosch",
        "boundary": "step-barrier",
        "work_function": "mid",
        "band_inputs": {t: se_yield.BAND[t] for t in SECONDARY_TARGETS},
        "band_provenance": {t: se_yield.BAND_PROVENANCE[t] for t in SECONDARY_TARGETS},
        "gaps": SECONDARY_GAPS,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="check the datasets only")
    ap.add_argument("--histories", type=int, default=100000)
    ap.add_argument("--threads", type=int, default=None,
                    help="threads per run (results do not depend on it)")
    ap.add_argument("--secondaries-only", action="store_true",
                    help="run only the fast-secondary sensitivity and merge it into the committed "
                         "backscatter_results.json (baseline entries stay as committed)")
    ap.add_argument("--targets", default=None,
                    help="with --secondaries-only / --verify-baseline-only: comma-separated targets to "
                         "run (default: all of SECONDARY_TARGETS); the others keep their committed entries")
    ap.add_argument("--verify-baseline-only", action="store_true",
                    help="only rerun the baseline at the secondary-sensitivity points and check it "
                         "reproduces the committed eta and delta bit for bit; writes nothing")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse a finished run under validation/oracle-runs/backscatter/ whose input is "
                         "byte-identical and whose summary names this binary's version (resumes an "
                         "interrupted sweep)")
    args = ap.parse_args()
    datasets = load_datasets()
    if args.check:
        print(f"{len(datasets)} backscatter datasets OK")
        return 0
    binary = lindhard_cli.lindhard_binary()
    reuse = lindhard_cli.version(binary) if args.reuse else None
    if args.secondaries_only or args.verify_baseline_only:
        res = json.loads(RESULTS.read_text())
        if res["histories"] != args.histories:
            lindhard_cli.die(f"--histories {args.histories} differs from the committed {res['histories']}")
        chosen = args.targets.split(",") if args.targets else SECONDARY_TARGETS
        bad = [t for t in chosen if t not in SECONDARY_TARGETS]
        if bad:
            lindhard_cli.die(f"--targets {bad} not in {SECONDARY_TARGETS}")
        for entry in res["targets"]:
            if entry["target"] in chosen:
                run_secondaries(binary, entry, args.histories, args.threads, reuse, True,
                                baseline_only=args.verify_baseline_only)
        if args.verify_baseline_only:
            print("baseline reproduces the committed eta and delta bit for bit")
            return 0
        # Each run keeps the `version` of the binary that produced it; the
        # top-level entry names the binary of this invocation.
        res["secondary_sensitivity"] = secondary_meta(binary)
        res["secondary_sensitivity"]["baseline_reproduced"] = True
        RESULTS.write_text(json.dumps(res, indent=1) + "\n")
        print(f"updated {RESULTS.relative_to(lindhard_cli.REPO)}")
        return 0
    targets = []
    for t in TARGETS:
        groups = [measured_group(datasets, t, e) for e in ENERGIES_KEV]
        inp = INPUTS / f"eta_{t.lower()}.toml"
        entry = {"target": t, "input": None, "groups": groups, "runs": {}, "sensitivity": {}}
        if inp.is_file():
            entry["input"] = str(inp.relative_to(lindhard_cli.REPO))
            base = inp.read_text()
            for e in ENERGIES_KEV:
                print(f"{t} {e:g} keV", file=sys.stderr)
                text = variant_input(base, inp.parent, e, args.histories, True, True)
                entry["runs"][f"{e:g}"] = run_one(
                    binary, text, lindhard_cli.RUNS / "backscatter" / f"{t}_{e:g}keV_both", args.threads, reuse)
            for e in SENSITIVITY_ENERGIES_KEV:
                row = {"both": entry["runs"][f"{e:g}"]}
                for vid, _label, ex, pol in VARIANTS[1:]:
                    print(f"{t} {e:g} keV {vid}", file=sys.stderr)
                    text = variant_input(base, inp.parent, e, args.histories, ex, pol)
                    row[vid] = run_one(
                        binary, text, lindhard_cli.RUNS / "backscatter" / f"{t}_{e:g}keV_{vid}", args.threads, reuse)
                entry["sensitivity"][f"{e:g}"] = row
        if t in SECONDARY_TARGETS and inp.is_file():
            run_secondaries(binary, entry, args.histories, args.threads, reuse, False)
        targets.append(entry)
    RESULTS.write_text(json.dumps({
        "format": "lindhard-backscatter-results/1",
        "lindhard": lindhard_cli.version(binary),
        "histories": args.histories,
        "seed": 1,
        "energies_kev": ENERGIES_KEV,
        "group_tolerance": GROUP_TOLERANCE,
        "pass_min_kev": PASS_MIN_KEV,
        "tolerance": TOLERANCE,
        "variants": [{"id": v[0], "label": v[1], "exchange": v[2], "polarization": v[3]} for v in VARIANTS],
        "secondary_sensitivity": secondary_meta(binary),
        "targets": targets,
    }, indent=1) + "\n")
    print(f"wrote {RESULTS.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
