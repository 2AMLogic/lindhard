#!/usr/bin/env python3
"""Rewrite the results tables in docs/validation.md between their marker
comments, so a rerun is idempotent:

    <!-- validation:level1:begin --> ... <!-- validation:level1:end -->
    <!-- validation:level2:begin --> ... <!-- validation:level2:end -->
    <!-- validation:electron-oracles:begin --> ... <!-- validation:electron-oracles:end -->
    <!-- validation:level3:begin --> ... <!-- validation:level3:end -->
    <!-- validation:level3-sputter:begin --> ... <!-- validation:level3-sputter:end -->
    <!-- validation:level3-sputter-crosscheck:begin --> ... <!-- validation:level3-sputter-crosscheck:end -->
    <!-- validation:level3-sputter-summary:begin --> ... <!-- validation:level3-sputter-summary:end -->

Level 1 comes from the Markdown the harness writes (LINDHARD_VALIDATION_OUT),
level 2 from the committed oracle summaries, level 3 from
validation/experiments/results.json. The sputter-yield blocks also read the
optional RustBCA context summaries
(validation/oracles/summaries/rustbca-sputter_ar_<target>.json) and the
committed digitizing cross-check records (validation/data/digitize/
crosscheck_ar_*_nifs23.json, secondread_ar_*_am32.json). Every number in the
sputter-yield interpretation, and every comparative word in it ("consistent
with one common value", "monotonic", "flat", "between"), is computed here by a
rule stated in the code next to it (#70). The summary block ends with a
sensitivity table: the same comparison with the sets named in
SENSITIVITY_SCENARIOS left out, regrouped from the committed datasets (#78).
Called by validation/run.sh.

Usage: validation/update_docs.py [--level1 FILE] [--check]

--check regenerates every block in memory and exits 1 if docs/validation.md
differs, without writing it (CI).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "validation.md"
SUMMARIES = ROOT / "validation" / "oracles" / "summaries"
RESULTS = ROOT / "validation" / "experiments" / "results.json"
SPUTTER_DATA = ROOT / "validation" / "data" / "sputtering"

sys.path.insert(0, str(ROOT / "validation" / "experiments"))
import run as experiments  # noqa: E402  (energy_groups: the one grouping rule of the level-3 sputter runs)


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
            continue  # e.g. the level-3 sputtering context (rustbca-sputter_ar_<target>.json)
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


ELECTRON_PROBLEMS = ROOT / "validation" / "oracles" / "electron_problems.json"
ELECTRON_ORACLES = (("Nebula", "nebula"), ("Geant4 MicroElec", "geant4_microelec"))


def electron_oracles() -> str:
    """The matched-problem table of the electron oracles (#150), from the
    committed `lindhard-electron-<problem>.json` and `<oracle>-<problem>.json`
    summaries that validation/oracles/run_electron.py writes. One row per code
    and problem: each metric with its batch-means standard error, and on an
    oracle's row the difference lindhard minus oracle (absolute for the
    yields, relative for the lengths) with its combined standard error and z.
    The tolerance column applies the rule stored in the summary itself."""
    if not ELECTRON_PROBLEMS.exists():
        return "_No electron problems defined._"
    spec = json.loads(ELECTRON_PROBLEMS.read_text())
    lines = [
        "| Problem | Code (version, commit) | Histories | η | δ | Primary depth (nm) | r50 (nm) | Tolerance (vs Nebula) |",
        "|---|---|---|---|---|---|---|---|",
    ]

    def val(x, se, digits: int) -> str:
        if x is None:
            return "n/a"
        return f"{x:.{digits}f} ± {se:.{digits}f}" if se is not None else f"{x:.{digits}f}"

    def diff(c: dict, metric: str) -> str:
        if metric in ("eta", "delta"):
            d, se, z = c.get(f"{metric}_abs_diff"), c.get(f"{metric}_abs_diff_se"), c.get(f"{metric}_z")
            if d is None:
                return ""
            s = f"{d:+.3f}" + (f" ± {se:.3f}" if se is not None else "")
        else:
            key = metric.removesuffix("_nm")
            d, se, z = c.get(f"{key}_rel_diff"), c.get(f"{key}_rel_diff_se"), c.get(f"{key}_z")
            if d is None:
                return ""
            s = f"{100 * d:+.1f} %" + (f" ± {100 * se:.1f} %" if se is not None else "")
        return f"<br>Δ {s}" + (f" ({abs(z):.1f}σ)" if z is not None else "")

    any_row = False
    for p in spec["problems"]:
        pid = p["id"]
        ours_path = SUMMARIES / f"lindhard-electron-{pid}.json"
        if not ours_path.exists():
            lines.append(f"| `{pid}` | lindhard | - | not run | not run | not run | not run | - |")
            continue
        any_row = True
        o = json.loads(ours_path.read_text())
        v, se = o["values"], o["std_err"]
        lines.append(
            f"| `{pid}` | lindhard ({o['lindhard_version']}) | {o['histories']} "
            f"| {val(v['eta'], se['eta'], 3)} | {val(v['delta'], se['delta'], 3)} "
            f"| {val(v['primary_depth_nm'], se['primary_depth_nm'], 1)} | {val(v['r50_nm'], se['r50_nm'], 1)} | - |"
        )
        for name, slug in ELECTRON_ORACLES:
            path = SUMMARIES / f"{slug}-{pid}.json"
            if not path.exists():
                why = ("not applicable (MicroElec: silicon only)"
                       if slug == "geant4_microelec" and not p.get("geant4_material") else "**not run**")
                lines.append(f"| `{pid}` | {name} | - | {why} | | | | - |")
                continue
            s = json.loads(path.read_text())
            if s.get("format") != "lindhard-oracle-electron-summary/1" or s.get("lindhard_version") != o["lindhard_version"]:
                sys.exit(f"error: {path.name} is not a lindhard-oracle-electron-summary/1 of the same lindhard run "
                         f"as {ours_path.name}; rerun validation/oracles/run_electron.py for {pid}")
            ov, ose, c = s["oracle_values"], s["oracle_values"]["std_err"], s["comparison"]
            t = s["tolerance"]
            if t.get("applies"):
                verdict = []
                for m, label in (("eta", "η"), ("r50_nm", "r50")):
                    ok = t["checks"][m]["pass"]
                    verdict.append(f"{label} " + ("n/a" if ok is None else "pass" if ok else "**FAIL**"))
                tol = ", ".join(verdict)
            else:
                tol = "reported only"
            commit = s.get("oracle_commit", "")[:7]
            lines.append(
                f"| `{pid}` | {name} ({s['oracle_version']}{', ' + commit if commit else ''}) | {s['histories']['oracle']} "
                f"| {val(ov['eta'], ose['eta'], 3)}{diff(c, 'eta')} "
                f"| {val(ov['delta'], ose['delta'], 3)}{diff(c, 'delta')} "
                f"| {val(ov['primary_depth_nm'], ose['primary_depth_nm'], 1)}{diff(c, 'primary_depth_nm')} "
                f"| {val(ov['r50_nm'], ose['r50_nm'], 1)}{diff(c, 'r50_nm')} | {tol} |"
            )
    if not any_row:
        return "_No electron oracle summaries committed yet; see \"Electron oracles\" above._"
    lines.append("")
    lines.append(
        "Values are pooled over all histories, ± the batch-means standard error "
        f"({spec['batches']} batches). Δ is lindhard minus the oracle: absolute for η and δ, relative "
        "for the lengths, ± the combined standard error, with the difference in units of it in brackets. "
        "η and δ split the front-face escapes at 50 eV (vacuum energy). r50 is the radius of the cylinder "
        "about the beam axis that holds half the energy deposited. Tolerances (vs Nebula at 5 and 20 keV "
        f"only): |Δη| ≤ {spec['tolerances']['eta_abs']} and |Δr50| ≤ {100 * spec['tolerances']['r50_rel']:.0f} %. "
        "Every summary's settings and its list of differing inputs are in its file under "
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


BANDS = ((0.0, 600.0, "0.2-0.6 keV"), (600.0, 2000.0, "0.6-2 keV"), (2000.0, math.inf, "2-10 keV"))  # lo <= E < hi
FLAT_BAND_RATIO = 1.15  # "flat over energy": the band geometric means agree within this factor
CROSSCHECKS = ROOT / "validation" / "data" / "digitize"


def gmean(xs) -> float:
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def sputter_data():
    """(sputtering section of results.json, targets sorted by mass, RustBCA summary per target or None)."""
    if not RESULTS.exists():
        return None, [], {}
    sp = json.loads(RESULTS.read_text()).get("sputtering")
    if not sp or not sp.get("rows"):
        return None, [], {}
    mass = sp.get("mass_amu", {})
    targets = sorted({r["target"] for r in sp["rows"]}, key=lambda t: (mass.get(t, 0.0), t))
    rb = {}
    for t in targets:
        f = SUMMARIES / f"rustbca-sputter_ar_{t.lower()}.json"
        rb[t] = json.loads(f.read_text()) if f.exists() else None
    return sp, targets, rb


def rb_ratios(sp, t, rb, k):
    """RustBCA K=k over the measured median at every energy of target t (None if not run)."""
    if not rb.get(t):
        return None
    at = {r["energy_ev"]: r for r in rb[t]["rows"]}
    rows = [r for r in sp["rows"] if r["target"] == t]
    if any(r["energy_ev"] not in at for r in rows):
        return None
    return [at[r["energy_ev"]][f"rustbca_k{k}"] / r["yield_median"] for r in rows]


def level3_sputter() -> str:
    sp, targets, rb = sputter_data()
    if not sp:
        return "_No sputter-yield datasets yet; see validation/data/README.md._"

    def y(row, k) -> str:
        v = row[f"lindhard_k{k}"]
        return f"{v:.2f} ({row[f'lindhard_k{k}_over_median']:.2f}){'' if not row[f'lindhard_k{k}_in_band'] else ' in'}"

    ions = sorted({row[f"lindhard_k{k}_ions"] for row in sp["rows"] for k in (0, 3) if f"lindhard_k{k}_ions" in row})
    lines = [
        f"lindhard, at least {sp['ions_per_run']} ions per run"
        + (f" (up to {ions[-1]} where the Poisson error would exceed 3 %)" if ions and ions[-1] > sp["ions_per_run"] else "")
        + f", seed {sp['seed']}; matched settings of `ar_1keV_cu_ed_es` "
        f"(ZBL, Lindhard-Scharff nonlocal, constant free path, cutoffs {sp['physics']['primary_cutoff_ev']:g} / "
        f"{sp['physics']['recoil_cutoff_ev']:g} eV, recoils followed, E_b = 0, E_d = E_s, with E_s lindhard's tabulated "
        "cohesive energy read back from a probe run: "
        + ", ".join(f"{v:g} eV ({k.split('->')[1]})" for k, v in sp["e_s_ev"].items())
        + f"). Measured energies within {100 * sp['energy_merge_rel']:g} % share one run.",
    ]
    for t in targets:
        rows = [r for r in sp["rows"] if r["target"] == t]
        at = {r["energy_ev"]: r for r in rb[t]["rows"]} if rb.get(t) else {}
        lines += [
            "",
            f"| Ar → {t}, E (eV) | Points (sets) | Measured min..max (median) | lindhard K = 0 (/median) "
            "| lindhard K = 3 (/median) | RustBCA K = 0 / K = 3 |",
            "|---|---|---|---|---|---|",
        ]
        for row in rows:
            r = at.get(row["energy_ev"])
            rbs = f"{r['rustbca_k0']:.2f} / {r['rustbca_k3']:.2f}" if r else "-"
            lines.append(
                f"| {row['energy_ev']:g} | {row['n_points']} ({len(row['sets'])}) | {row['yield_min']:.2f}..{row['yield_max']:.2f} "
                f"({row['yield_median']:.2f}) | {y(row, 0)} | {y(row, 3)} | {rbs} |"
            )
        k0 = [row["lindhard_k0_over_median"] for row in rows]
        k3 = [row["lindhard_k3_over_median"] for row in rows]
        inb0 = sum(row["lindhard_k0_in_band"] for row in rows)
        inb3 = sum(row["lindhard_k3_in_band"] for row in rows)
        f = f"validation/oracles/summaries/rustbca-sputter_ar_{t.lower()}.json"
        lines += [
            "",
            f"Ar → {t}: K = 0 lies in the measured band at {inb0} of {len(k0)} energies (ratio {min(k0):.2f} to "
            f"{max(k0):.2f}), K = 3 at {inb3} (ratio {min(k3):.2f} to {max(k3):.2f}). "
            + (f"RustBCA {rb[t]['oracle_version']} ({rb[t]['ions']} ions, same settings through the level-2 adapter; "
               f"summary `{f}`, context only)." if rb.get(t) else
               "RustBCA: not run (set RUSTBCA_BIN and run `validation/experiments/run.py --rustbca`)."),
        ]
    se_max = max(row[f"lindhard_k{k}_se"] / row[f"lindhard_k{k}"] for row in sp["rows"] for k in (0, 3))
    lines += [
        "",
        f"(/median): lindhard over the measured median; \"in\": inside the measured min..max. Statistical errors "
        f"(Poisson) are at most {100 * se_max:.2f} % of the yield. Measured points: `validation/data/sputtering/`; "
        "provenance: [`data-provenance.md`](data-provenance.md).",
    ]
    c = sp.get("control")
    if c:
        lines += ["", f"Control: the 1 keV Ar → Cu run gives {c['ar_cu_1keV_k0']:.4f}, the level-2 `ar_1keV_cu_ed_es` value "
                      f"{c['level2_ar_1keV_cu_ed_es']:.4f} ({c['z']:+.1f} σ)."]
    return "\n".join(lines)


def level3_sputter_crosscheck() -> str:
    """Per-point tables of the NIFS-DATA-23 reads and per-target summaries of the
    manual second reads, from the committed records; the coincidence rule of
    yamamura1995_nifs23.py (`coinc()`) is replayed here on every pair."""
    sp, targets, _ = sputter_data()
    if not sp:
        return "_No sputter-yield datasets yet._"
    lines = []
    for t in targets:
        f = CROSSCHECKS / f"crosscheck_ar_{t.lower()}_nifs23.json"
        if not f.exists():
            continue
        rec = json.loads(f.read_text())
        pairs = [p for p in rec["pairs"] if "yield_ratio" in p]
        npts = len(rec["pairs"])
        lines += [
            f"**Ar → {t}**, {rec['figure_b'].split(' (')[0]}: {len(pairs)} of {npts} stored points compared "
            f"(NIFS calibration rms {rec['calibration_b_rms_decade']:.4f} decade).",
            "",
            "| Set (NIFS letter) | E (eV) | Y, IPPJ-AM-32 | Y, NIFS-DATA-23 | Y ratio | E ratio | Combined unc. | Symbols "
            "| Verdict |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        n1 = n2 = nexp = nunexp = ndis = 0
        for p in pairs:
            lr = abs(math.log(p["yield_ratio"]))
            le = abs(math.log(p["energy_ratio"]))
            u = p["combined_unc_rel"]
            worst = max(lr, le)
            if worst <= u:
                v = "within 1σ"
                n1 += 1
            else:
                what = "Y" if lr >= le else "E"
                tier = "within 2σ" if worst <= 2 * u else "beyond 2σ"
                ys = [q["yield_am32"] for q in rec["pairs"] if q.get("yield_nifs") == p["yield_nifs"]]
                explained = p["legibility"] == "coincident" and min(ys) <= p["yield_nifs"] <= max(ys)
                v = f"{'' if what == 'Y' else 'E '}{tier}, " + ("explained (coincident)" if explained else "unexplained")
                n2 += tier == "within 2σ"
                ndis += tier != "within 2σ"
                nexp += explained
                nunexp += not explained
            m = re.match(r"ar_\w+?_sputter_([a-z]+)(\d{4})", p["id"])
            lines.append(
                f"| {m.group(1).capitalize()} {m.group(2)} ({p['nifs_letter']}) | {p['energy_ev']:g} | {p['yield_am32']:.3f} "
                f"| {p['yield_nifs']:.3f} | {p['yield_ratio']:.3f} | {p['energy_ratio']:.3f} | {100 * u:.1f} % "
                f"| {p['legibility']} | {v} |"
            )
        if pairs:
            lr = [math.log(p["yield_ratio"]) for p in pairs]
            le = [math.log(p["energy_ratio"]) for p in pairs]
            mu, me = sum(lr) / len(lr), sum(le) / len(le)
            sr = math.sqrt(sum((x - mu) ** 2 for x in lr) / len(lr))
            se = math.sqrt(sum((x - me) ** 2 for x in le) / len(le))
            lines += [
                "",
                f"Mean Y ratio {math.exp(mu):.3f} (rms {100 * sr:.1f} %), mean E ratio {math.exp(me):.3f} (rms {100 * se:.1f} %); "
                f"{n1} within 1σ, {n2} within 2σ ({nexp} explained by coincident symbols under the rule that the NIFS read "
                f"lies within the IPPJ-AM-32 reads drawn on that spot, {nunexp} not), {ndis} beyond 2σ; σ is the "
                f"combined digitizing uncertainty of the two yield reads, used for both ratios. "
                f"Record: `{f.relative_to(ROOT)}`.",
            ]
        s = CROSSCHECKS / f"secondread_ar_{t.lower()}_am32.json"
        if s.exists():
            sr2 = json.loads(s.read_text())["points"]
            ly = [math.log(p["yield_ratio"]) for p in sr2]
            le2 = [math.log(p["energy_ratio"]) for p in sr2]
            within = sum(p["verdict"] == "within the digitizing uncertainty" for p in sr2)
            lines += [
                "",
                f"Manual second read (glyph boxes against template centres, `{s.relative_to(ROOT)}`): {within} of "
                f"{len(sr2)} points within the stored digitizing uncertainty; rms of ln(manual/template) "
                f"{math.sqrt(sum(x * x for x in ly) / len(ly)):.4f} in Y and {math.sqrt(sum(x * x for x in le2) / len(le2)):.4f} "
                f"in E, largest {max(abs(x) for x in ly + le2):.4f}.",
            ]
        lines.append("")
    return "\n".join(lines).strip()


# --- sensitivity of the summary to doubtful sets (#78) ---
#
# Each scenario leaves out existing datasets (by `id`, the file stem in validation/data/sputtering/)
# of one target, for a doubt that is already documented: `caveat` names where (the dataset's own
# `target_state` or `original_reference`, or the "Titles that name another system" paragraph of
# docs/validation.md, section 3). These are sensitivity scenarios only: every stored set stays in
# the baseline comparison, nothing is deleted, and whether an original reports the plotted Ar yield
# is not settled here (the originals were not read). Inclusion is decided by this table only,
# never by matching words in the datasets' prose.
CAVEAT_DOC = "docs/validation.md"
CAVEAT_SOURCES = ("target_state", "original_reference", CAVEAT_DOC)
SENSITIVITY_SCENARIOS = (
    {"target": "Si", "exclude": ("ar_si_sputter_poate1976",), "caveat": CAVEAT_DOC,
     "reason": "cited paper's Crossref title names PtSi and NiSi"},
    {"target": "Ag", "exclude": ("ar_ag_sputter_wehner1961",), "caveat": "original_reference",
     "reason": "cited paper's Crossref title is Hg+ at 4-15 keV"},
    {"target": "Ag", "exclude": ("ar_ag_sputter_okajima1981",), "caveat": CAVEAT_DOC,
     "reason": "cited paper's Crossref title names O2+"},
    {"target": "Au", "exclude": ("ar_au_sputter_robinson1967",), "caveat": "target_state",
     "reason": "monocrystalline target per the original's title; a target-state doubt, not an attribution one"},
    {"target": "Au", "exclude": ("ar_au_sputter_szymonski1978",), "caveat": "original_reference",
     "reason": "cited paper's Crossref title is 6 keV Xe+ on an AgAu alloy"},
    {"target": "Au", "exclude": ("ar_au_sputter_holloway1977",), "caveat": CAVEAT_DOC,
     "reason": "cited paper's Crossref title names Cr in Au"},
    {"target": "Au", "exclude": ("ar_au_sputter_robinson1967", "ar_au_sputter_szymonski1978",
                                 "ar_au_sputter_holloway1977"), "caveat": "target_state, original_reference, " + CAVEAT_DOC,
     "reason": "the three Au rows above together"},
)


def load_sputter_sets(directory: Path = SPUTTER_DATA) -> list[tuple[Path, dict]]:
    """The committed sputter-yield datasets, read only (run.py --check validates them)."""
    files = sorted(directory.glob("*.json")) if directory.is_dir() else []
    return [(f, json.loads(f.read_text())) for f in files]


def check_scenarios(scenarios, datasets) -> list[str]:
    """Every problem with the scenario table, as actionable messages; empty if it is usable."""
    by_id = {d.get("id"): d for _, d in datasets}
    errs, seen = [], {}
    for i, s in enumerate(scenarios):
        where = f"SENSITIVITY_SCENARIOS[{i}] (target {s.get('target')!r})"
        ids = list(s.get("exclude") or ())
        if not ids:
            errs.append(f"{where}: `exclude` is empty; the baseline is added automatically, list at least one dataset id")
        for x in sorted({x for x in ids if ids.count(x) > 1}):
            errs.append(f"{where}: dataset id {x!r} is listed more than once in `exclude`")
        for x in ids:
            d = by_id.get(x)
            if d is None:
                errs.append(f"{where}: unknown dataset id {x!r}; ids are the file stems of "
                            "validation/data/sputtering/*.json")
            elif d.get("target") != s.get("target"):
                errs.append(f"{where}: dataset {x!r} is an {d.get('ion')} -> {d.get('target')} set, "
                            f"not {s.get('target')}; move it to a scenario of its own target")
        for c in (c.strip() for c in str(s.get("caveat", "")).split(",")):
            if c not in CAVEAT_SOURCES:
                errs.append(f"{where}: `caveat` {c!r} is not one of {', '.join(CAVEAT_SOURCES)}")
            elif c != CAVEAT_DOC and len(ids) == 1 and ids[0] in by_id and not str(by_id[ids[0]].get(c, "")).strip():
                errs.append(f"{where}: `caveat` names the dataset field {c!r}, which is empty in {ids[0]!r}")
        if not str(s.get("reason", "")).strip():
            errs.append(f"{where}: `reason` is empty")
        key = (s.get("target"), frozenset(ids))
        if ids and key in seen:
            errs.append(f"{where}: same target and exclusions as SENSITIVITY_SCENARIOS[{seen[key]}]")
        seen.setdefault(key, i)
    return errs


def code_yield(curve: dict, e: float):
    """A code's yield at energy `e` from its committed runs `curve` (energy -> yield):
    (yield, "exact") at a run energy; (yield, "interpolated") between the two bracketing run
    energies, linear in log E - log Y (positive yields only); (None, why) otherwise. Never
    extrapolates. An interpolated yield approximates the committed code curve; it is not a run."""
    if e in curve:
        y = curve[e]
        return (y, "exact") if y > 0 else (None, "nonpositive yield")
    lo = [x for x in curve if x < e]
    hi = [x for x in curve if x > e]
    if not lo or not hi:
        return None, "outside the committed run energies"
    e0, e1 = max(lo), min(hi)
    y0, y1 = curve[e0], curve[e1]
    if y0 <= 0 or y1 <= 0:
        return None, "nonpositive yield"
    f = math.log(e / e0) / math.log(e1 / e0)
    return math.exp(math.log(y0) + f * math.log(y1 / y0)), "interpolated"


def scenario_comparison(datasets, ion: str, target: str, exclude, curves: dict) -> dict:
    """One scenario: the datasets of `target` minus `exclude`, grouped by the unchanged rule of
    run.energy_groups (2 % grouping, incidence and flag rules, median), then each code's yield
    (`curves`: name -> {energy: yield}, or None if that code was not run) over the new median.
    The inputs are not modified."""
    kept = [(p, d) for p, d in datasets if d["ion"] == ion and d["target"] == target and d["id"] not in set(exclude)]
    groups = [g for g in experiments.energy_groups(kept) if (g["ion"], g["target"]) == (ion, target)]
    out = {
        "groups": groups,
        "n_groups": len(groups),
        "n_sets": len({s for g in groups for s in g["sets"]}),
        "n_points": sum(g["n_points"] for g in groups),
        "codes": {},
    }
    for name, curve in curves.items():
        if curve is None:
            out["codes"][name] = None
            continue
        ratios, missing, interpolated = [], [], 0
        for g in groups:
            y, how = code_yield(curve, g["energy_ev"])
            if y is None:
                missing.append((g["energy_ev"], how))
            else:
                ratios.append(y / g["yield_median"])
                interpolated += how == "interpolated"
        out["codes"][name] = {"ratios": ratios, "missing": missing, "interpolated": interpolated}
    return out


def ratio_summary(c, n_groups: int) -> str:
    """Geometric mean (range) of one code's ratios, marked partial if any group is unsupported,
    'unavailable' if there is nothing to compare; never a fabricated number."""
    if c is None:
        return "unavailable (not run)"
    if not c["ratios"]:
        return "unavailable (no energies)" if n_groups == 0 else "unavailable (no supported energy)"
    s = f"{gmean(c['ratios']):.2f} ({min(c['ratios']):.2f}-{max(c['ratios']):.2f})"
    if c["missing"]:
        s += f", partial: {len(c['ratios'])} of {n_groups} energies"
    return s


def sputter_sensitivity(sp, targets, rb, datasets, scenarios=SENSITIVITY_SCENARIOS) -> list[str]:
    """The generated sensitivity table (#78); [] if no scenario applies to a committed target."""
    errs = check_scenarios(scenarios, datasets)
    if errs:
        sys.exit("error: " + "\n       ".join(errs))
    ion = "Ar"
    lines = [
        "**Sensitivity to doubtful sets (scenarios, not a new baseline).** The table above keeps every stored set. "
        "Below, the sets named in each row are left out *before* the measured points are regrouped by the same rule "
        "(energies within 2 %, flagged points left out, median), so the median and the representative energy of a "
        "group can change and a group whose only set is left out disappears. lindhard's and RustBCA's yields are the "
        "committed runs: exact at a run energy, otherwise interpolated linearly in log E - log Y between the two "
        "bracketing run energies (an approximation to the committed code curve, not a new run and not an "
        "uncertainty); nothing is extrapolated, and an energy outside the runs is counted as unsupported. The "
        "reasons are the caveats stored in each dataset's `target_state` or `original_reference`, or listed under "
        "\"Titles that name another system\" above; none of the originals was read, so none is settled here.",
        "",
        "| Target | Left out | Caveat (where stored) | Energies (sets, points) | Interpolated | Unsupported "
        "| lindhard K = 0 / median | RustBCA K = 0 / median |",
        "|---|---|---|---|---|---|---|---|",
    ]
    rows = 0
    for t in targets:
        sc = [s for s in scenarios if s["target"] == t]
        if not sc:
            continue
        curves = {"lindhard": {r["energy_ev"]: r["lindhard_k0"] for r in sp["rows"] if r["target"] == t},
                  "rustbca": {r["energy_ev"]: r["rustbca_k0"] for r in rb[t]["rows"]} if rb.get(t) else None}
        base = scenario_comparison(datasets, ion, t, (), curves)
        committed = [r for r in sp["rows"] if r["target"] == t]
        # The no-exclusion scenario must reproduce the committed baseline: same groups, same
        # medians, and the same printed geometric means; otherwise results.json is stale.
        if [(g["energy_ev"], g["yield_median"], g["n_points"]) for g in base["groups"]] != \
                [(r["energy_ev"], r["yield_median"], r["n_points"]) for r in committed]:
            sys.exit(f"error: regrouping the committed Ar -> {t} datasets does not reproduce the groups in "
                     f"{RESULTS.relative_to(ROOT)}; rerun validation/experiments/run.py")
        printed = f"{gmean([r['lindhard_k0_over_median'] for r in committed]):.2f}"
        if ratio_summary(base["codes"]["lindhard"], base["n_groups"]).split(" ")[0] != printed:
            sys.exit(f"error: the Ar -> {t} no-exclusion scenario does not reproduce the printed K = 0 ratio {printed}")
        for s in ({"exclude": (), "caveat": "-", "reason": "baseline"}, *sc):
            r = base if not s["exclude"] else scenario_comparison(datasets, ion, t, s["exclude"], curves)
            interp = r["codes"]["lindhard"]["interpolated"]
            rb_c = r["codes"]["rustbca"]
            left = "none (baseline)" if not s["exclude"] else ", ".join(f"`{x}`" for x in s["exclude"])
            why = "-" if not s["exclude"] else f"{s['reason']} ({s['caveat']})"
            interp_txt = f"{interp}" + ("" if rb_c is None or rb_c["interpolated"] == interp
                                        else f" (RustBCA {rb_c['interpolated']})")
            missing = [f"{name} at {e:g} eV ({how})" for name, c in r["codes"].items() if c for e, how in c["missing"]]
            lines.append(
                f"| {t} | {left} | {why} | {r['n_groups']} ({r['n_sets']}, {r['n_points']}) | {interp_txt} "
                f"| {'; '.join(missing) if missing else 'none'} "
                f"| {ratio_summary(r['codes']['lindhard'], r['n_groups'])} | {ratio_summary(rb_c, r['n_groups'])} |"
            )
            rows += 1
    if not rows:
        return []
    lines += [
        "",
        "Geometric mean (range) over the regrouped energies, as in the table above; \"Interpolated\" counts the "
        "energies whose code yield is interpolated rather than a committed run. The baseline rows are recomputed "
        "by this regrouping and reproduce the table above to its printed precision (update_docs.py stops "
        "otherwise). Scenarios: `SENSITIVITY_SCENARIOS` in `validation/update_docs.py`; semantics: "
        "`validation/data/README.md`.",
    ]
    return lines


def level3_sputter_summary() -> str:
    """The per-target summary table and the interpretation, every number and
    every comparative word computed here from results.json and the RustBCA
    summaries (#70)."""
    sp, targets, rb = sputter_data()
    if not sp or "mass_amu" not in sp:
        return "_Rerun validation/experiments/run.py to record the masses._"
    m1 = sp["mass_amu"]["Ar"]
    es = {k.split("->")[1]: v for k, v in sp["e_s_ev"].items()}
    by = {t: [r for r in sp["rows"] if r["target"] == t] for t in targets}

    def band_means(rows, key):
        out = []
        for lo, hi, _ in BANDS:
            sel = [r[key] for r in rows if lo <= r["energy_ev"] < hi]
            out.append((gmean(sel), len(sel)) if sel else (None, 0))
        return out

    def rng(xs) -> str:
        return f"{gmean(xs):.2f} ({min(xs):.2f}-{max(xs):.2f})"

    lines = [
        "| Target | M2/M1 | E_s = E_d (eV) | Energies (sets, points) | lindhard K = 0 / median | lindhard K = 3 / median "
        "| RustBCA K = 0 / median | RustBCA K = 3 / median | lindhard K = 0 by band: "
        + " / ".join(b[2] for b in BANDS) + " |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    stats = {}
    for t in targets:
        rows = by[t]
        k0 = [r["lindhard_k0_over_median"] for r in rows]
        k3 = [r["lindhard_k3_over_median"] for r in rows]
        r0, r3 = rb_ratios(sp, t, rb, 0), rb_ratios(sp, t, rb, 3)
        bm = band_means(rows, "lindhard_k0_over_median")
        nsets = len({s for r in rows for s in r["sets"]})
        npts = sum(r["n_points"] for r in rows)
        stats[t] = dict(k0=k0, k3=k3, r0=r0, r3=r3, bm=bm, mass=sp["mass_amu"][t] / m1, es=es[t])
        lines.append(
            f"| {t} | {sp['mass_amu'][t] / m1:.2f} | {es[t]:g} | {len(rows)} ({nsets}, {npts}) | {rng(k0)} | {rng(k3)} "
            f"| {rng(r0) if r0 else '-'} | {rng(r3) if r3 else '-'} | "
            + " / ".join(f"{g:.2f} ({n})" if g else "-" for g, n in bm) + " |"
        )
    lines += [
        "",
        "Geometric mean over the merged energies of each target's ratio to the measured median, with its range in "
        "brackets; band columns give the geometric mean and, in brackets, the number of energies in the band.",
        "",
        "**Interpretation (generated from the table; no tuning, no default changed).**",
        "",
    ]
    g0 = {t: gmean(stats[t]["k0"]) for t in targets}
    spread = {}
    for t in targets:
        multi = [r for r in by[t] if len(r["sets"]) > 1]
        spread[t] = (max(r["yield_max"] / r["yield_min"] for r in multi), len(multi)) if multi else (None, 0)
    lines.append(
        "- **Measured scatter.** At energies where more than one set has a point, the sets differ by up to a factor "
        + ", ".join((f"{spread[t][0]:.2f} ({t}, {spread[t][1]} of {len(by[t])} energies)" if spread[t][0]
                     else f"- ({t}, no such energy)") for t in targets)
        + "; elsewhere the median is a single set."
    )
    # Is the K = 0 ratio one common, target-independent value? Each target's geometric mean has a
    # 95 % interval from its own energy-to-energy scatter (Student t, n - 1 degrees of freedom);
    # the ratio is called consistent with one common value only if all intervals share a point.
    tq = [12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.36, 2.31, 2.26, 2.23, 2.20, 2.18, 2.16, 2.14, 2.13, 2.12, 2.11, 2.10,
          2.09, 2.09, 2.08, 2.07, 2.07, 2.06, 2.06, 2.06, 2.05, 2.05, 2.05, 2.04]  # t(0.975, df = 1..30)
    iv = {}
    for t in targets:
        lr = [math.log(x) for x in stats[t]["k0"]]
        n = len(lr)
        mu = sum(lr) / n
        sd = math.sqrt(sum((x - mu) ** 2 for x in lr) / (n - 1)) if n > 1 else float("inf")
        half = (tq[min(n - 1, 30) - 1] if n > 1 else float("inf")) * sd / math.sqrt(n)
        iv[t] = (math.exp(mu - half), math.exp(mu + half), sd)
    common = max(v[0] for v in iv.values()) <= min(v[1] for v in iv.values())
    lo_t, hi_t = min(targets, key=g0.get), max(targets, key=g0.get)
    lines.append(
        f"- **Across targets the K = 0 ratio is "
        + ("consistent with one common value" if common else "not one common value")
        + f".** Its geometric means are "
        + ", ".join(f"{g0[t]:.2f} ({t})" for t in targets)
        + f", a factor {g0[hi_t] / g0[lo_t]:.2f} between {lo_t} and {hi_t}. Within each target the energy-to-energy "
        "scatter of ln(ratio) is "
        + ", ".join(f"{iv[t][2]:.2f} ({t})" for t in targets)
        + "; the 95 % intervals of the means that this scatter implies are "
        + ", ".join(f"{iv[t][0]:.2f}-{iv[t][1]:.2f} ({t})" for t in targets)
        + (", and they all overlap." if common else ", and they do not all overlap.")
        + " These intervals treat energies as independent, which the shared laboratories and the one-set-per-energy "
        "medians are not, so they understate the uncertainty."
    )

    def order(key):
        return [t for t in sorted(targets, key=lambda t: stats[t][key])]

    def monotone(seq):
        vals = [g0[t] for t in seq]
        if all(a < b for a, b in zip(vals, vals[1:])):
            return "rising"
        if all(a > b for a, b in zip(vals, vals[1:])):
            return "falling"
        return None

    mm, me = monotone(order("mass")), monotone(order("es"))
    p_chance = 2 / math.factorial(len(targets))

    def trend(name, key, mono):
        seq = order(key)
        s = f"in order of {name} ({', '.join(f'{t} {stats[t][key]:.2f}' for t in seq)}) the means run " + \
            ", ".join(f"{g0[t]:.2f}" for t in seq)
        return s + (f", monotonically {mono}" if mono else ", not monotonically")

    lines.append(
        f"- **No trend can be claimed from {len(targets)} targets.** "
        + (lambda s: s[0].upper() + s[1:])(trend("M2/M1", "mass", mm)) + "; " + trend("E_s", "es", me) + ". "
        + f"With {len(targets)} targets a monotonic order in a chosen variable arises by chance with probability "
        f"{p_chance:.2f} even when there is no trend, so "
        + ("even the monotonic order above would be weak evidence" if (mm or me) else "the data support no trend in either")
        + "; no physical cause is identified, and the target dependence is **unexplained**."
    )
    flat = {t: None for t in targets}
    parts = []
    for t in targets:
        bm = [g for g, _ in stats[t]["bm"] if g]
        if len(bm) >= 2:
            ratio = max(bm) / min(bm)
            flat[t] = ratio <= FLAT_BAND_RATIO
            parts.append(f"{t} {ratio:.2f}")
    flats = [t for t in targets if flat[t]]
    notflat = [t for t in targets if flat[t] is False]
    lines.append(
        f"- **Energy dependence.** The ratio of the largest to the smallest band mean is " + ", ".join(parts)
        + f"; it is flat over energy (band means within a factor {FLAT_BAND_RATIO}) for "
        + (", ".join(flats) if flats else "no target")
        + (f", and not for {', '.join(notflat)}" if notflat else "")
        + ". Where it is not flat, a single factor does not describe the deficit."
    )
    below = {t: sum(r["lindhard_k0"] < r["yield_min"] for r in by[t]) for t in targets}
    above = {t: sum(r["lindhard_k0"] > r["yield_max"] for r in by[t]) for t in targets}
    inb = {t: sum(r["lindhard_k0_in_band"] for r in by[t]) for t in targets}
    lines.append(
        "- **Against the measured band**, lindhard K = 0 is below it at "
        + ", ".join(f"{below[t]} of {len(by[t])} energies ({t})" for t in targets)
        + ", inside it at " + ", ".join(f"{inb[t]} ({t})" for t in targets)
        + " and above it at " + ", ".join(f"{above[t]} ({t})" for t in targets) + "."
    )
    k3lower = sum(r["lindhard_k3"] < r["lindhard_k0"] for r in sp["rows"])
    further = [t for t in targets if abs(math.log(gmean(stats[t]["k3"]))) > abs(math.log(g0[t]))]
    lines.append(
        f"- **Weak collisions (K = 3)** give a lower lindhard yield than K = 0 at {k3lower} of {len(sp['rows'])} "
        "energies; their geometric mean ratio to the median is "
        + ", ".join(f"{gmean(stats[t]['k3']):.2f} ({t})" for t in targets)
        + ", further from the measured median than K = 0 for "
        + (", ".join(further) if further else "no target")
        + ". No default is changed here (#61)."
    )
    have_rb = [t for t in targets if stats[t]["r0"]]
    if have_rb:
        frac, desc = {}, []
        for t in have_rb:
            gl, gr = g0[t], gmean(stats[t]["r0"])
            if gl < 1 < gr:
                frac[t] = math.log(1 / gl) / math.log(gr / gl)
                desc.append(f"{t}: lindhard {gl:.2f}, RustBCA {gr:.2f}, the median {100 * frac[t]:.0f} % of the log "
                            "distance from lindhard toward RustBCA")
            else:
                desc.append(f"{t}: lindhard {gl:.2f}, RustBCA {gr:.2f}, the median not between them")
        rabove = {t: sum(r > 1 for r in stats[t]["r0"]) for t in have_rb}
        nearer = [t for t in have_rb if abs(math.log(gmean(stats[t]["r0"]))) < abs(math.log(g0[t]))]
        lr = {t: gmean([l / r * 1.0 for l, r in zip(stats[t]["k0"], stats[t]["r0"])]) for t in have_rb}
        lines.append(
            "- **RustBCA (context, K = 0)**: " + "; ".join(desc) + ". RustBCA lies above the measured median at "
            + ", ".join(f"{rabove[t]} of {len(stats[t]['r0'])} energies ({t})" for t in have_rb)
            + ". The measured median lies between the two codes' geometric means for "
            + (", ".join(frac) if frac else "no target")
            + ("" if len(frac) == len(have_rb) else f", and not for {', '.join(t for t in have_rb if t not in frac)}")
            + "; RustBCA's geometric mean is the nearer to the median for "
            + (", ".join(nearer) if nearer else "no target")
            + ". lindhard is "
            + ", ".join(f"{lr[t]:.2f} ({t})" for t in have_rb)
            + " of RustBCA (geometric mean over energies): the gap of #61 is present for "
            + ("every target" if all(v < 1 for v in lr.values()) else ", ".join(t for t in have_rb if lr[t] < 1))
            + "; it stays **unexplained**."
        )
    sens = sputter_sensitivity(sp, targets, rb, load_sputter_sets())
    if sens:
        lines += [""] + sens
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level1", type=Path, help="Markdown table written by the level-1 harness")
    ap.add_argument("--check", action="store_true",
                    help="regenerate in memory and exit 1 if docs/validation.md differs (CI)")
    args = ap.parse_args()
    old = DOC.read_text()
    text = old
    if args.level1:
        text = splice(text, "level1", args.level1.read_text())
    text = splice(text, "level2", level2())
    text = splice(text, "electron-oracles", electron_oracles())
    text = splice(text, "level3", level3())
    text = splice(text, "level3-sputter", level3_sputter())
    text = splice(text, "level3-sputter-crosscheck", level3_sputter_crosscheck())
    text = splice(text, "level3-sputter-summary", level3_sputter_summary())
    if args.check:
        if text != old:
            print(f"error: {DOC.relative_to(ROOT)} is not what validation/update_docs.py generates from the committed "
                  "results; run it and commit the result", file=sys.stderr)
            return 1
        print(f"{DOC.relative_to(ROOT)} matches the committed results")
        return 0
    DOC.write_text(text)
    print(f"updated {DOC.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
