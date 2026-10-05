#!/usr/bin/env python3
"""Level-3 comparison with published range measurements
(docs/validation.md section 3; data and schema: validation/data/README.md).

For every `validation/data/ranges/*.json`:

1. checks the dataset's provenance fields and that `docs/data-provenance.md`
   has a row naming it (any failure stops the run, exit 2);
2. runs `lindhard` on the same ion, energy, tilt and target with the default
   physics below and records measured and computed moments;
3. reruns it with the nuclear and the electronic stopping varied separately
   (the stopping-input attribution, docs/validation.md section 3), using only
   lindhard's existing models: the CLI's interatomic potentials and screening
   lengths for the nuclear part, and the Lindhard-Scharff stopping scaled by
   a constant factor k (a user stopping table, generated here) for the
   electronic part. Nothing is fitted or tuned; the defaults do not change.

Writes `validation/experiments/results.json`.

Usage:
    validation/experiments/run.py [--ions N]   check, run, write results
    validation/experiments/run.py --check      check the datasets only
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

DATA = HERE.parent / "data" / "ranges"
RESULTS = HERE / "results.json"
PROVENANCE = lindhard_cli.REPO / "docs" / "data-provenance.md"

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

# Stopping-input attribution. Each variant changes either the nuclear part
# (potential and screening length, all existing CLI choices; docs/cli.md) or
# the electronic part (k times the Lindhard-Scharff stopping; or the same LS
# magnitude split into local Oen-Robinson and nonlocal halves), not both at
# once, except the last, which reproduces the combination that Wittmaack
# and Mutzke, J. Appl. Phys. 121, 105104 (2017), report as best for B in Si
# (Kr-C with the Lindhard screening length, k = 1.46). k = 1.46 and 2 are the
# values quoted in that paper and in the PR #46 review; they are probes, not
# fits, and no default changes.
VARIANTS = [
    {"id": "zbl_k1", "part": "default", "potential": "zbl", "k": 1.0, "label": "ZBL, LS"},
    {"id": "zbl_k1.46", "part": "electronic", "potential": "zbl", "k": 1.46, "label": "ZBL, 1.46 LS"},
    {"id": "zbl_k2", "part": "electronic", "potential": "zbl", "k": 2.0, "label": "ZBL, 2 LS"},
    {
        "id": "zbl_eqLSOR",
        "part": "electronic form",
        "potential": "zbl",
        "stopping": "equipartition-ls-or",
        "k": 1.0,
        "label": "ZBL, LS/Oen-Robinson",
    },
    {"id": "krc_k1", "part": "nuclear", "potential": "kr-c", "k": 1.0, "label": "Kr-C (Firsov a), LS"},
    {"id": "moliere_k1", "part": "nuclear", "potential": "moliere", "k": 1.0, "label": "Moliere (Firsov a), LS"},
    {
        "id": "krc_aL_k1",
        "part": "nuclear",
        "potential": "kr-c",
        "screening_length": "lindhard",
        "k": 1.0,
        "label": "Kr-C (Lindhard a), LS",
    },
    {
        "id": "krc_aL_k1.46",
        "part": "both",
        "potential": "kr-c",
        "screening_length": "lindhard",
        "k": 1.46,
        "label": "Kr-C (Lindhard a), 1.46 LS",
    },
]
# Families over which the k needed to match each measurement is interpolated.
K_FAMILIES = {
    "zbl": ["zbl_k1", "zbl_k1.46", "zbl_k2"],
    "krc_aL": ["krc_aL_k1", "krc_aL_k1.46"],
}

# --- dataset checks (the provenance hooks of validation/data/README.md) ---

REQUIRED_TEXT = [
    "id",
    "ion",
    "target",
    "target_state",
    "citation",
    "source_location",
    "extraction",
    "terms",
    "added",
]


def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def check_dataset(path: Path, d: dict, provenance: str) -> list[str]:
    """Every problem with one dataset file, as messages; empty if it passes."""
    errs = []
    for k in REQUIRED_TEXT:
        if not isinstance(d.get(k), str) or not d[k].strip():
            errs.append(f"`{k}` missing or empty")
    if "short_citation" in d and (not isinstance(d["short_citation"], str) or not d["short_citation"].strip()):
        errs.append("`short_citation`, if given, must be a non-empty string")
    if not any(isinstance(d.get(k), str) and d[k].strip() for k in ("doi", "url")):
        errs.append("needs a `doi` or a `url` for the source that was read")
    if d.get("id") != path.stem:
        errs.append(f"`id` {d.get('id')!r} differs from the file name {path.stem!r}")
    if not (_num(d.get("energy_ev")) and d["energy_ev"] > 0):
        errs.append("`energy_ev` must be a positive number")
    if "tilt_deg" in d and not _num(d["tilt_deg"]):
        errs.append("`tilt_deg` must be a number")
    m = d.get("measured")
    if not isinstance(m, dict):
        errs.append("`measured` missing")
    else:
        if not isinstance(m.get("method"), str) or not m["method"].strip():
            errs.append("`measured.method` missing or empty")
        for k in ("rp_nm", "rp_unc_nm"):
            if not (_num(m.get(k)) and m[k] > 0):
                errs.append(f"`measured.{k}` must be a positive number")
        drp, unc = m.get("drp_nm"), m.get("drp_unc_nm")
        if drp is None:
            if unc is not None:
                errs.append("`measured.drp_unc_nm` given without `drp_nm`")
        elif not (_num(drp) and drp > 0 and _num(unc) and unc > 0):
            errs.append("`measured.drp_nm` and `drp_unc_nm` must both be positive numbers, or both null")
    if isinstance(d.get("id"), str) and d["id"] not in provenance:
        errs.append(f"no row in docs/data-provenance.md names `{d['id']}`")
    return errs


def load_datasets() -> list[tuple[Path, dict]]:
    """All datasets, checked, ordered by ion and energy. Exits 2 on any failure."""
    files = sorted(DATA.glob("*.json")) if DATA.is_dir() else []
    provenance = PROVENANCE.read_text()
    out, failed = [], False
    for f in files:
        rel = f.relative_to(lindhard_cli.REPO)
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
        out.append((f, d))
    if failed:
        print("error: dataset provenance check failed (rules: validation/data/README.md)", file=sys.stderr)
        sys.exit(2)
    return sorted(out, key=lambda fd: (fd[1]["ion"], fd[1]["energy_ev"], fd[1]["id"]))


# --- Lindhard-Scharff tables for the electronic-stopping probes ---

# The constants and formula of lindhard/src/ion/stopping/lindhard_scharff.rs
# (Lindhard and Scharff, Phys. Rev. 124, 128 (1961); Lindhard, Scharff and
# Schiott, Mat. Fys. Medd. 33 (14) (1963)) and lindhard/src/constants.rs
# (CODATA 2022). A k = 1 table must reproduce the built-in model; run() checks
# that it does.
ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878188e-12
BOHR_RADIUS = 5.29177210544e-11
COULOMB_E2 = ELEMENTARY_CHARGE**2 / (4.0 * math.pi * VACUUM_PERMITTIVITY)
# Z and lindhard's standard atomic weight (lindhard/src/elements.rs) of the
# species the datasets use.
ELEMENTS = {"B": (5, 10.81), "Si": (14, 28.085), "P": (15, 30.973761998), "As": (33, 74.921595)}


def ls_stopping(z1: int, m1: float, z2: int, m2: float, e_ev: float) -> float:
    """Lindhard-Scharff electronic stopping cross section, eV 1e-15 cm^2."""
    s23 = z1 ** (2.0 / 3.0) + z2 ** (2.0 / 3.0)
    a = 0.8853 * BOHR_RADIUS / math.sqrt(s23)
    eps = e_ev * ELEMENTARY_CHARGE * a * m2 / (z1 * z2 * COULOMB_E2 * (m1 + m2))
    k_l = (
        0.0793 * z1 ** (1.0 / 6.0) * math.sqrt(z1 * z2) * (m1 + m2) ** 1.5
        / (s23**0.75 * m1**1.5 * math.sqrt(m2))
    )
    to_si = 4.0 * math.pi * a * z1 * z2 * COULOMB_E2 * m1 / (m1 + m2)  # J m^2
    return k_l * math.sqrt(eps) * to_si / ELEMENTARY_CHARGE / 1e-19


def ls_table(d: dict, k: float) -> Path:
    """Write a k x LS table for the dataset's (ion, target) pair; return its path."""
    z1, w1 = ELEMENTS[d["ion"]]
    z2, m2 = ELEMENTS[d["target"]]
    m1 = d.get("mass_amu") or w1
    e_max = 1.5 * d["energy_ev"]
    n = 80
    energies = [math.exp(math.log(1.0) + (math.log(e_max) - math.log(1.0)) * i / n) for i in range(n + 1)]
    lines = [
        f'provenance = "k = {k:g} times the Lindhard-Scharff stopping (Phys. Rev. 124, 128 (1961)), '
        f'computed by validation/experiments/run.py for the level-3 attribution; not data"',
        f"ion_z = {z1}",
        f"ion_mass_amu = {m1!r}",
        f"target_z = {z2}",
        "energy_ev = [" + ", ".join(repr(e) for e in energies) + "]",
        "stopping_ev_1e15_cm2 = [" + ", ".join(repr(k * ls_stopping(z1, m1, z2, m2, e)) for e in energies) + "]",
    ]
    path = lindhard_cli.RUNS / "experiments" / "tables" / f"{d['ion']}_{m1:g}_{d['target']}_{e_max:g}eV_k{k:g}.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def problem(d: dict, ions: int, variant: dict | None = None, force_table: bool = False) -> dict:
    beam = {"ion": d["ion"], "energy_ev": d["energy_ev"], "tilt_deg": d.get("tilt_deg", 0.0)}
    if d.get("mass_amu"):
        beam["mass_amu"] = d["mass_amu"]
    physics = dict(PHYSICS)
    vid = ""
    if variant:
        physics["potential"] = variant["potential"]
        for key in ("screening_length", "stopping"):
            if key in variant:
                physics[key] = variant[key]
        vid = "." + variant["id"]
    p = {
        "id": d["id"] + vid + (".table" if force_table else ""),
        "beam": beam,
        "target": {"substrate": d["target"]},
        "physics": {**physics, "energies": {d["target"]: {"e_d_ev": 15.0}}},
        "run": {"ions": ions, "seed": 1},
        "tally": {"depth_bin_nm": max(d["measured"]["rp_nm"] / 100.0, 0.1), "depth_bins": 1000},
    }
    k = variant["k"] if variant else 1.0
    if k != 1.0 or force_table:
        p["stopping_tables"] = [ls_table(d, k)]
    return p


def k_needed(points: list[tuple[float, float]], target: float) -> float | None:
    """k at which Rp(k) meets `target`, interpolating ln Rp linearly in ln k
    between the probed k. Outside the probed Rp range the nearest pair is
    extrapolated, which assumes a power law Rp ~ k^s that no run verifies;
    `is_extrapolated` flags those results."""
    pts = sorted(points)
    if len(pts) < 2:
        return None
    lt = math.log(target)
    pairs = list(zip(pts, pts[1:]))
    for (k0, r0), (k1, r1) in pairs:
        if min(r0, r1) <= target <= max(r0, r1):
            break
    else:
        (k0, r0), (k1, r1) = pairs[0] if target > pts[0][1] else pairs[-1]
    slope = (math.log(r1) - math.log(r0)) / (math.log(k1) - math.log(k0))
    return math.exp(math.log(k0) + (lt - math.log(r0)) / slope)


def is_extrapolated(points: list[tuple[float, float]], k: float | None) -> bool:
    """True if `k` lies outside the probed k, so `k_needed` extrapolated to it."""
    ks = [p[0] for p in points]
    return k is not None and not min(ks) <= k <= max(ks)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ions", type=int, default=20000)
    ap.add_argument("--check", action="store_true", help="check the datasets' provenance fields only")
    args = ap.parse_args()

    datasets = load_datasets()
    if args.check:
        print(f"{len(datasets)} dataset(s) pass the provenance check")
        return 0
    if not datasets:
        print(f"no datasets in {DATA.relative_to(lindhard_cli.REPO)}; see validation/data/README.md")
        return 0
    binary = lindhard_cli.lindhard_binary()
    runs = lindhard_cli.RUNS / "experiments"
    rows, attribution, control = [], [], []
    for _, d in datasets:
        m = d["measured"]
        ours = lindhard_cli.run(problem(d, args.ions), runs / d["id"], binary)
        se = ours["drp_nm"] / math.sqrt(ours["ions"] * (1.0 - ours["backscatter"]))
        drp = m.get("drp_nm")
        row = {
            "id": d["id"],
            "case": f"{d['ion']} {d['energy_ev'] / 1e3:g} keV -> {d['target']} ({m.get('method', '?').split(' ')[0]})",
            "citation": d["citation"],
            "short_citation": d.get("short_citation") or d["citation"],
            "rp_measured_nm": m["rp_nm"],
            "rp_unc_nm": m["rp_unc_nm"],
            "rp_lindhard_nm": ours["rp_nm"],
            "rp_lindhard_se_nm": se,
            "rp_rel_diff": ours["rp_nm"] / m["rp_nm"] - 1.0,
            "rp_diff_sigma": (ours["rp_nm"] - m["rp_nm"]) / math.hypot(m["rp_unc_nm"], se),
            "drp_measured_nm": drp,
            "drp_unc_nm": m.get("drp_unc_nm"),
            "drp_lindhard_nm": ours["drp_nm"],
            "drp_rel_diff": None if drp is None else ours["drp_nm"] / drp - 1.0,
            "drp_diff_sigma": None if drp is None else (ours["drp_nm"] - drp) / m["drp_unc_nm"],
            "backscatter": ours["backscatter"],
        }
        rows.append(row)
        print(
            f"{row['case']}: Rp {ours['rp_nm']:.2f} vs {m['rp_nm']:.2f} nm ({100 * row['rp_rel_diff']:+.1f} %)"
            + ("" if drp is None else f", dRp {ours['drp_nm']:.1f} vs {drp:.1f} nm ({100 * row['drp_rel_diff']:+.1f} %)")
        )

        # Control: a k = 1 table must reproduce the built-in Lindhard-Scharff model.
        tab = lindhard_cli.run(problem(d, args.ions, force_table=True), runs / (d["id"] + ".table"), binary)
        control.append(abs(tab["rp_nm"] / ours["rp_nm"] - 1.0))

        rp = {}
        for v in VARIANTS:
            if v["id"] == "zbl_k1":
                rp[v["id"]] = ours["rp_nm"]
                continue
            r = lindhard_cli.run(problem(d, args.ions, v), runs / f"{d['id']}.{v['id']}", binary)
            rp[v["id"]] = r["rp_nm"]
        kn = {}
        for fam, ids in K_FAMILIES.items():
            pts = [(next(v["k"] for v in VARIANTS if v["id"] == i), rp[i]) for i in ids]
            kn[fam] = {
                "k": k_needed(pts, m["rp_nm"]),
                "k_low": k_needed(pts, m["rp_nm"] + m["rp_unc_nm"]),
                "k_high": k_needed(pts, m["rp_nm"] - m["rp_unc_nm"]),
            }
            for key in ("k", "k_low", "k_high"):
                kn[fam][key + "_extrapolated"] = is_extrapolated(pts, kn[fam][key])
        attribution.append(
            {
                "id": d["id"],
                "case": row["case"],
                "rp_measured_nm": m["rp_nm"],
                "rp_unc_nm": m["rp_unc_nm"],
                "rp_nm": rp,
                "ratio": {k: v / m["rp_nm"] for k, v in rp.items()},
                "k_needed": kn,
            }
        )
        print(
            "  attribution: "
            + ", ".join(f"{v['id']} {rp[v['id']] / m['rp_nm']:.3f}" for v in VARIANTS)
            + "; k needed: "
            + ", ".join(f"{f} {kn[f]['k']:.2f}" for f in K_FAMILIES)
        )

    worst = max(control)
    if worst > 2e-3:
        lindhard_cli.die(
            f"a k = 1 Lindhard-Scharff table changed Rp by {100 * worst:.3f} %: run.py's LS formula "
            "no longer matches lindhard's; fix ls_stopping() before trusting the attribution"
        )
    RESULTS.write_text(
        json.dumps(
            {
                "format": "lindhard-experiment-results/2",
                "lindhard_version": lindhard_cli.version(binary),
                "ions": args.ions,
                "physics": PHYSICS,
                "datasets": rows,
                "attribution": {
                    "variants": VARIANTS,
                    "k_families": K_FAMILIES,
                    "table_control_max_rel_diff": worst,
                    "rows": attribution,
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(f"k = 1 table control: Rp reproduced to {worst:.1e} (relative)")
    print(f"wrote {RESULTS.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
