#!/usr/bin/env python3
"""Rewrite the results tables in docs/validation.md between their marker
comments, so a rerun is idempotent:

    <!-- validation:level1:begin --> ... <!-- validation:level1:end -->
    <!-- validation:level2:begin --> ... <!-- validation:level2:end -->
    <!-- validation:level3:begin --> ... <!-- validation:level3:end -->
    <!-- validation:level3-sputter:begin --> ... <!-- validation:level3-sputter:end -->

Level 1 comes from the Markdown the harness writes (LINDHARD_VALIDATION_OUT),
level 2 from the committed oracle summaries, level 3 from
validation/experiments/results.json (the sputter-yield table also reads the
optional RustBCA context summary, validation/oracles/summaries/rustbca-sputter_ar_cu.json). Called by validation/run.sh.

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
        "| Problem | Oracle (version, commit) | Rp | ΔRp | Backscatter (abs.) | Sputter yield "
        "| Speed ratio (end-to-end) | Speed ratio (marginal) | Date |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    def sig(c: dict, key: str) -> str:
        z = c.get(f"{key}_z")
        return "" if z is None else f" ({abs(z):.1f}σ)"

    def ratio(x) -> str:
        return "-" if x is None else f"{x:.2f}x"

    for f in files:
        s = json.loads(f.read_text())
        if s.get("format") != "lindhard-oracle-summary/2":
            continue  # e.g. the level-3 sputtering context (rustbca-sputter_ar_cu.json)
        c = s["comparison"]
        bs = c.get("backscatter_abs_diff")
        commit = s.get("oracle_commit", "")[:7]
        oracle = f"{s['oracle']} ({s['oracle_version']}{', ' + commit if commit else ''})"
        lines.append(
            f"| `{s['problem']}` | {oracle} | {pct(c.get('rp_rel_diff'))}{sig(c, 'rp')} "
            f"| {pct(c.get('drp_rel_diff'))}{sig(c, 'drp')} "
            f"| {'-' if bs is None else f'{bs:+.4f}'}{sig(c, 'backscatter')} "
            f"| {pct(c.get('sputter_yield_rel_diff'))}{sig(c, 'sputter_yield')} "
            f"| {ratio(c.get('ions_per_s_ratio'))} | {ratio(c.get('ions_per_s_marginal_ratio'))} | {s['date']} |"
        )
    lines.append("")
    lines.append(
        "Differences are lindhard relative to the oracle, with the difference in units of its combined "
        "standard error in brackets (both runs' statistics; the sputter-yield error assumes Poisson counts). "
        "A speed ratio > 1 means lindhard is faster: end-to-end is process wall clock at the run's ion "
        "count, marginal removes fixed setup costs (see each summary's `timing`). Every row's settings, "
        "both sides' values and the full list of mismatches are in its file under "
        "`validation/oracles/summaries/`."
    )
    return "\n".join(lines)


def level3() -> str:
    if not RESULTS.exists():
        return "_No experimental datasets yet; see validation/data/README.md._"
    r = json.loads(RESULTS.read_text())

    def mm(x, unc=None) -> str:
        if x is None:
            return "-"
        return f"{x:.2f}" + ("" if unc is None else f" ± {unc:.2f}")

    lines = [
        f"{r['lindhard_version']}, {r['ions']} ions per case, physics {r['physics']['potential']} + "
        f"{r['physics']['stopping']} (the defaults).",
        "",
        "| Case | Rp measured (nm) | Rp lindhard (nm) | Diff. | Diff. / σ | ΔRp measured (nm) "
        "| ΔRp lindhard (nm) | Diff. | Source |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for d in r["datasets"]:
        lines.append(
            f"| {d['case']} | {mm(d['rp_measured_nm'], d['rp_unc_nm'])} | {d['rp_lindhard_nm']:.2f} "
            f"| {pct(d['rp_rel_diff'])} | {d['rp_diff_sigma']:+.1f} | {mm(d['drp_measured_nm'], d['drp_unc_nm'])} "
            f"| {d['drp_lindhard_nm']:.2f} | {pct(d['drp_rel_diff'])} | {d.get('short_citation', d['citation'])} |"
        )
    lines.append("")
    lines.append(
        "σ combines the measurement uncertainty with the Monte Carlo standard error of lindhard's Rp. "
        "Datasets, their extraction and uncertainty: `validation/data/ranges/`; provenance: "
        "[`data-provenance.md`](data-provenance.md)."
    )
    a = r.get("attribution")
    if not a or not a.get("rows"):
        return "\n".join(lines)
    variants = a["variants"]
    fams = list(a["k_families"])
    fam_label = {"zbl": "k needed, ZBL", "krc_aL": "k needed, Kr-C (Lindhard a)"}
    lines += [
        "",
        "**Stopping-input attribution.** lindhard's Rp relative to the measurement, with the nuclear "
        "or the electronic stopping varied separately (existing models only; LS = Lindhard-Scharff, "
        "k LS = LS times k, through a generated user table; a = screening length). "
        "\"k needed\" is the LS factor at which the potential's Rp would meet the measurement "
        "(range from the measurement uncertainty). A value without a star lies between the probed k "
        "(1, 1.46 and 2 for ZBL; 1 and 1.46 for Kr-C) and is a log-log interpolation. "
        "**A starred value (*) lies outside the probed k: it extrapolates the nearest pair of "
        "runs, assuming Rp follows a power law in k there, and no run verifies it.** "
        "It is a sensitivity estimate, not a result, and the measurement-uncertainty range does not "
        "include the extrapolation uncertainty.",
        "",
        "| Case | Rp measured (nm) | "
        + " | ".join(v["label"] for v in variants)
        + " | "
        + " | ".join(fam_label.get(f, f"k needed, {f}") for f in fams)
        + " |",
        "|---|---|" + "---|" * (len(variants) + len(fams)),
    ]
    for row in a["rows"]:

        def kn(f: str) -> str:
            k = row["k_needed"][f]
            if k["k"] is None:
                return "-"
            def v(key: str) -> str:
                return f"{k[key]:.2f}" + ("*" if k.get(key + "_extrapolated") else "")

            return f"{v('k')} ({v('k_low')}-{v('k_high')})"

        lines.append(
            f"| {row['case']} | {mm(row['rp_measured_nm'], row['rp_unc_nm'])} | "
            + " | ".join(pct(row["ratio"][v["id"]] - 1.0) for v in variants)
            + " | "
            + " | ".join(kn(f) for f in fams)
            + " |"
        )
    lines.append("")
    lines.append(
        f"Control: a k = 1 table reproduces the built-in Lindhard-Scharff Rp to "
        f"{a['table_control_max_rel_diff']:.1e} (relative), so the k tables change only the magnitude. "
        "Per-run values: `validation/experiments/results.json`."
    )
    return "\n".join(lines)


def level3_sputter() -> str:
    if not RESULTS.exists():
        return "_No experimental datasets yet; see validation/data/README.md._"
    sp = json.loads(RESULTS.read_text()).get("sputtering")
    if not sp or not sp.get("rows"):
        return "_No sputter-yield datasets yet; see validation/data/README.md._"
    rb_path = SUMMARIES / "rustbca-sputter_ar_cu.json"
    rb = json.loads(rb_path.read_text()) if rb_path.exists() else None
    rb_at = {r["energy_ev"]: r for r in rb["rows"]} if rb else {}

    def y(row, k) -> str:
        v = row[f"lindhard_k{k}"]
        return f"{v:.2f} ({row[f'lindhard_k{k}_over_median']:.2f}){'' if not row[f'lindhard_k{k}_in_band'] else ' in'}"

    lines = [
        f"lindhard, {sp['ions_per_run']} ions per run, seed {sp['seed']}; matched settings of `ar_1keV_cu_ed_es` "
        f"(ZBL, Lindhard-Scharff nonlocal, constant free path, cutoffs {sp['physics']['primary_cutoff_ev']:g} / "
        f"{sp['physics']['recoil_cutoff_ev']:g} eV, recoils followed, E_b = 0, E_d = E_s = "
        + ", ".join(f"{v:g} eV ({k})" for k, v in sp["e_s_ev"].items())
        + ", lindhard's tabulated cohesive energy). Measured energies within "
        f"{100 * sp['energy_merge_rel']:g} % share one run.",
        "",
        "| Ar → Cu, E (eV) | Points (sets) | Measured min..max (median) | lindhard K = 0 (/median) "
        "| lindhard K = 3 (/median) | RustBCA K = 0 / K = 3 |",
        "|---|---|---|---|---|---|",
    ]
    for row in sp["rows"]:
        if (row["ion"], row["target"]) != ("Ar", "Cu"):
            continue
        r = rb_at.get(row["energy_ev"])
        rbs = f"{r['rustbca_k0']:.2f} / {r['rustbca_k3']:.2f}" if r else "-"
        lines.append(
            f"| {row['energy_ev']:g} | {row['n_points']} ({len(row['sets'])}) | {row['yield_min']:.2f}..{row['yield_max']:.2f} "
            f"({row['yield_median']:.2f}) | {y(row, 0)} | {y(row, 3)} | {rbs} |"
        )
    k0 = [row["lindhard_k0_over_median"] for row in sp["rows"]]
    k3 = [row["lindhard_k3_over_median"] for row in sp["rows"]]
    se_max = max(row[f"lindhard_k{k}_se"] / row[f"lindhard_k{k}"] for row in sp["rows"] for k in (0, 3))
    inb0 = sum(row["lindhard_k0_in_band"] for row in sp["rows"])
    inb3 = sum(row["lindhard_k3_in_band"] for row in sp["rows"])
    lines += [
        "",
        f"(/median): lindhard over the measured median; \"in\": inside the measured min..max. K = 0 lies in the "
        f"band at {inb0} of {len(k0)} energies (ratio {min(k0):.2f} to {max(k0):.2f}), K = 3 at {inb3} "
        f"(ratio {min(k3):.2f} to {max(k3):.2f}). Statistical errors (Poisson) are at most {100 * se_max:.2f} % of the yield. "
        + (f"RustBCA {rb['oracle_version']} ({rb['ions']} ions, same settings through the level-2 adapter; summary "
           "`validation/oracles/summaries/rustbca-sputter_ar_cu.json`, context only)." if rb else
           "RustBCA: not run (set RUSTBCA_BIN and run `validation/experiments/run.py --rustbca`).")
        + " Measured points: `validation/data/sputtering/`; provenance: [`data-provenance.md`](data-provenance.md).",
    ]
    c = sp.get("control")
    if c:
        lines += ["", f"Control: the 1 keV run gives {c['ar_cu_1keV_k0']:.4f}, the level-2 `ar_1keV_cu_ed_es` value "
                      f"{c['level2_ar_1keV_cu_ed_es']:.4f} ({c['z']:+.1f} σ)."]
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
    text = splice(text, "level3-sputter", level3_sputter())
    DOC.write_text(text)
    print(f"updated {DOC.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
