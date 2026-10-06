#!/usr/bin/env python3
"""Level-3 comparison with published measurements: projected ranges and
sputter yields (docs/validation.md section 3; data and schema:
validation/data/README.md).

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

For every `validation/data/sputtering/*.json` (measured sputter yields):

1. checks it the same way (`check_sputter_dataset`);
2. merges the measured energies of each ion-target pair into groups within
   2 % and runs `lindhard` once per group with the matched settings of the
   level-2 problem `ar_1keV_cu_ed_es` (SPUTTER_PHYSICS below), with
   `weak_collisions` = 0 and 3; records the measured band (min, median, max
   over every stored point in the group) next to both yields. Each run uses
   at least `--ions` ions, and more where the yield is low: a run whose
   Poisson error exceeds 3 % is repeated with the ion count scaled to bring
   it to 2.5 % (rounded up to 1000; deterministic, seed 1). Nothing is
   fitted to the measurements and no default changes.

Writes `validation/experiments/results.json`.

Usage:
    validation/experiments/run.py [--ions N]   check, run, write results
    validation/experiments/run.py --check      check the datasets only
    validation/experiments/run.py --rustbca [--rustbca-target T ...]
                                               also run RustBCA (RUSTBCA_BIN, level-2
                                               adapter) at the sputtering energies of each
                                               target (default: every target), for context;
                                               writes the summaries
                                               validation/oracles/summaries/rustbca-sputter_ar_<t>.json
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
SPUTTER_DATA = HERE.parent / "data" / "sputtering"
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


# --- sputter yields (validation/data/README.md, "Sputter-yield datasets") ---

SPUTTER_REQUIRED_TEXT = [
    "id",
    "ion",
    "target",
    "target_state",
    "original_reference",
    "compilation",
    "compilation_figure",
    "extraction",
    "terms",
    "added",
]

# Matched settings: those of the level-2 problem `ar_1keV_cu_ed_es`
# (validation/oracles/problems.json): ZBL, Lindhard-Scharff all nonlocal,
# constant free path, cutoffs 2 eV (primary) and 1 eV (recoils), recoils
# followed, E_b = 0, and E_d = E_s with E_s lindhard's tabulated cohesive
# energy (read back from a probe run, never hardcoded). E_d <= E_s is
# Eckstein's convention when yields matter (Computer Simulation of Ion-Solid
# Interactions, 1991); #61 showed E_d is inert once E_d <= E_s.
SPUTTER_PHYSICS = {
    "potential": "zbl",
    "stopping": "lindhard-scharff",
    "free_path": "constant",
    "primary_cutoff_ev": 2.0,
    "recoil_cutoff_ev": 1.0,
    "follow_recoils": True,
}
SPUTTER_WEAK = (0, 3)  # weak collisions per step: lindhard's default, and TRIDYN's maximum (#64)
MERGE_REL = 0.02  # measured energies within 2 % share one lindhard run
MAX_POISSON_REL = 0.03  # rerun with more ions above this relative Poisson error (#70)
TARGET_POISSON_REL = 0.025  # ... aiming at this one
AVOGADRO = 6.02214076e23  # exact (SI 2019)
PROBE_E_D_EV = 1.0  # E_d of the E_s read-back probe only (sputter_problem); never in a recorded run


def check_sputter_dataset(path: Path, d: dict, provenance: str) -> list[str]:
    """Every problem with one sputter-yield dataset, as messages; empty if it passes."""
    errs = []
    if d.get("kind") != "sputter_yield":
        errs.append("`kind` must be \"sputter_yield\"")
    for k in SPUTTER_REQUIRED_TEXT:
        if not isinstance(d.get(k), str) or not d[k].strip():
            errs.append(f"`{k}` missing or empty")
    if not any(isinstance(d.get(k), str) and d[k].strip() for k in ("original_doi", "url")):
        errs.append("needs an `original_doi` or a `url` for the source that was read")
    if d.get("id") != path.stem:
        errs.append(f"`id` {d.get('id')!r} differs from the file name {path.stem!r}")
    if not (_num(d.get("compilation_pdf_page")) and d["compilation_pdf_page"] > 0):
        errs.append("`compilation_pdf_page` must be a positive number")
    if not _num(d.get("incidence_deg")):
        errs.append("`incidence_deg` missing or not a number")
    pts = d.get("points")
    if not isinstance(pts, list) or not pts:
        errs.append("`points` missing or empty")
    else:
        for i, p in enumerate(pts):
            if not isinstance(p, dict):
                errs.append(f"`points[{i}]` is not an object")
                continue
            for k in ("energy_ev", "yield", "energy_unc_rel", "yield_unc_rel"):
                if not (_num(p.get(k)) and p[k] > 0):
                    errs.append(f"`points[{i}].{k}` must be a positive number")
    if isinstance(d.get("id"), str) and d["id"] not in provenance:
        errs.append(f"no row in docs/data-provenance.md names `{d['id']}`")
    return errs


def load_sputter_datasets() -> list[tuple[Path, dict]]:
    """All sputter-yield datasets, checked. Exits 2 on any failure."""
    files = sorted(SPUTTER_DATA.glob("*.json")) if SPUTTER_DATA.is_dir() else []
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
        errs = check_sputter_dataset(f, d, provenance)
        for e in errs:
            print(f"error: {rel}: {e}", file=sys.stderr)
        failed |= bool(errs)
        out.append((f, d))
    if failed:
        print("error: dataset provenance check failed (rules: validation/data/README.md)", file=sys.stderr)
        sys.exit(2)
    return out


def energy_groups(datasets: list[tuple[Path, dict]]) -> list[dict]:
    """Measured points of each (ion, target), normal incidence, grouped by
    energy: a group starts at its lowest energy and takes every point up to
    2 % above it; its run energy is the geometric mean of its points. A point
    flagged `disagrees_between_compilations` is left out."""
    pts = {}
    for _, d in datasets:
        if d["incidence_deg"] != 0.0:
            continue
        for p in d["points"]:
            if p.get("flag") == "disagrees_between_compilations":
                continue  # a second read disagreed beyond twice the combined uncertainty (data README)
            pts.setdefault((d["ion"], d["target"]), []).append((p["energy_ev"], p["yield"], d["id"]))
    groups = []
    for (ion, target), xs in sorted(pts.items()):
        xs.sort()
        cur = []
        for x in xs:
            if cur and x[0] > cur[0][0] * (1 + MERGE_REL):
                groups.append((ion, target, cur))
                cur = []
            cur.append(x)
        if cur:
            groups.append((ion, target, cur))
    out = []
    for ion, target, g in groups:
        ys = sorted(y for _, y, _ in g)
        n = len(ys)
        med = ys[n // 2] if n % 2 else 0.5 * (ys[n // 2 - 1] + ys[n // 2])
        e = math.exp(sum(math.log(x[0]) for x in g) / n)
        out.append({
            "ion": ion,
            "target": target,
            "energy_ev": float(f"{e:.4g}"),
            "measured_energies_ev": [x[0] for x in g],
            "n_points": n,
            "sets": sorted({x[2] for x in g}),
            "yield_min": ys[0],
            "yield_median": med,
            "yield_max": ys[-1],
        })
    return out


def sputter_problem(ion: str, target: str, energy_ev: float, weak: int, ions: int, e_s_ev: float | None) -> dict:
    """The matched lindhard problem at one energy. `e_s_ev` None: a probe
    with the target's default E_s, to read it back; the probe sets E_d
    (PROBE_E_D_EV) because some elements (Si) have no default E_d, and its
    yield is not used."""
    energies = ({target: {"e_d_ev": PROBE_E_D_EV}} if e_s_ev is None
                else {target: {"e_d_ev": e_s_ev, "e_b_ev": 0.0, "e_s_ev": e_s_ev}})
    physics = dict(SPUTTER_PHYSICS)
    if weak:
        physics["weak_collisions"] = weak
    return {
        "id": f"sputter_{ion.lower()}_{target.lower()}_{energy_ev:g}eV_k{weak}" + ("" if e_s_ev is not None else "_probe"),
        "beam": {"ion": ion, "energy_ev": energy_ev, "tilt_deg": 0.0},
        "target": {"substrate": target},
        "physics": {**physics, "energies": energies},
        "run": {"ions": ions, "seed": 1},
        "tally": {"depth_bin_nm": 0.1, "depth_bins": 100},
    }


def tabulated_e_s(ion: str, target: str, binary: Path) -> float:
    """lindhard's default E_s of `target` (its tabulated cohesive energy), from a probe run's echo."""
    return probe_target(ion, target, binary)[0]


def probe_target(ion: str, target: str, binary: Path) -> tuple[float, float]:
    """lindhard's default E_s of `target` and its atomic mass (amu), both from a
    probe run's echo; the mass is density x N_A / atom density of the echoed
    pure-element target, i.e. the weight lindhard uses."""
    p = sputter_problem(ion, target, 1000.0, 0, 10, None)
    r = lindhard_cli.run(p, lindhard_cli.RUNS / "experiments" / p["id"], binary)
    (el,) = [e for e in r["target_elements"] if e["symbol"] == target]
    if not el.get("e_s_ev"):
        lindhard_cli.die(f"lindhard has no tabulated E_s for {target}; the matched settings need one")
    layer = json.loads((lindhard_cli.RUNS / "experiments" / p["id"] / "out" / "summary.json").read_text())[
        "physics"]["target"][0]
    mass = layer["material"]["density_g_cm3"] * AVOGADRO / layer["atom_density_per_cm3"]
    return el["e_s_ev"], float(f"{mass:.6g}")


def ion_mass_amu(symbol: str) -> float:
    """The ion's mass from the oracle adapter's element table (validation/oracles/run.py,
    ELEMENTS, the standard atomic weights of lindhard/src/elements.rs); the CLI echo
    does not carry the ion mass."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("oracles_run", HERE.parent / "oracles" / "run.py")
    orc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orc)
    return orc.ELEMENTS[symbol][1]


def run_sputter(ions: int, binary: Path, groups: list[dict]) -> dict:
    runs = lindhard_cli.RUNS / "experiments"
    e_s, mass = {}, {}
    rows = []
    for g in groups:
        key = (g["ion"], g["target"])
        if key not in e_s:
            e_s[key], mass[g["target"]] = probe_target(*key, binary)
        row = dict(g)
        for k in SPUTTER_WEAK:
            n = ions
            while True:
                p = sputter_problem(g["ion"], g["target"], g["energy_ev"], k, n, e_s[key])
                r = lindhard_cli.run(p, runs / p["id"], binary)
                rel = r["sputter_yield_se"] / r["sputter_yield"] if r["sputter_yield"] > 0 else math.inf
                if rel <= MAX_POISSON_REL:
                    break
                n = 1000 * math.ceil(n * (rel / TARGET_POISSON_REL) ** 2 / 1000)
            row[f"lindhard_k{k}_ions"] = n
            row[f"lindhard_k{k}"] = r["sputter_yield"]
            row[f"lindhard_k{k}_se"] = float(f"{r['sputter_yield_se']:.4g}")
            row[f"lindhard_k{k}_in_band"] = g["yield_min"] <= r["sputter_yield"] <= g["yield_max"]
            row[f"lindhard_k{k}_over_median"] = float(f"{r['sputter_yield'] / g['yield_median']:.4g}")
        rows.append(row)
        print(f"{g['ion']} {g['energy_ev']:g} eV -> {g['target']}: measured {g['yield_min']:.3g}..{g['yield_max']:.3g} "
              f"(median {g['yield_median']:.3g}, {g['n_points']} pt), lindhard K=0 {row['lindhard_k0']:.3f}, "
              f"K=3 {row['lindhard_k3']:.3f}")
    # Sanity: the 1 keV problem must reproduce the committed level-2 value of ar_1keV_cu_ed_es.
    summ = lindhard_cli.REPO / "validation" / "oracles" / "summaries" / "rustbca-ar_1keV_cu_ed_es.json"
    control = None
    if ("Ar", "Cu") in e_s and summ.exists():
        ref = json.loads(summ.read_text())["lindhard"]
        p = sputter_problem("Ar", "Cu", 1000.0, 0, 20000, e_s[("Ar", "Cu")])
        r = lindhard_cli.run(p, runs / p["id"], binary)
        z = (r["sputter_yield"] - ref["sputter_yield"]) / math.hypot(r["sputter_yield_se"], ref["std_err"]["sputter_yield"])
        control = {"ar_cu_1keV_k0": r["sputter_yield"], "level2_ar_1keV_cu_ed_es": ref["sputter_yield"], "z": round(z, 3)}
        print(f"control: Ar 1 keV -> Cu, K=0: {r['sputter_yield']:.4f} (level-2 ar_1keV_cu_ed_es: {ref['sputter_yield']:.4f})")
        if abs(z) > 3:
            lindhard_cli.die("the matched 1 keV problem no longer reproduces ar_1keV_cu_ed_es; the settings differ")
    return {
        "physics": SPUTTER_PHYSICS,
        "weak_collisions": list(SPUTTER_WEAK),
        "e_s_ev": {f"{i}->{t}": v for (i, t), v in sorted(e_s.items())},
        # Target masses: lindhard's, from the probe echo (probe_target). The ion
        # mass is the oracle adapter's table entry, which mirrors
        # lindhard/src/elements.rs (validation/oracles/run.py, ELEMENTS).
        "mass_amu": {**{i: ion_mass_amu(i) for i, _ in sorted(e_s)}, **dict(sorted(mass.items()))},
        "e_d_ev": "equal to e_s_ev",
        "e_b_ev": 0.0,
        "energy_merge_rel": MERGE_REL,
        "ions_per_run": ions,
        "ions_rule": f"at least {ions}; rerun with more where the Poisson error exceeds {100 * MAX_POISSON_REL:g} % "
                     "(per-row `lindhard_k*_ions`)",
        "seed": 1,
        "std_err": "Poisson, sqrt(sputtered) / ions; underestimates the true error (correlated bursts)",
        "control": control,
        "rows": rows,
    }


def run_rustbca_context(ions: int, binary: Path, groups: list[dict], targets: list[str] | None = None) -> None:
    """RustBCA at the same energies and settings through the level-2 adapter
    (validation/oracles/run.py), for context only; writes one summary of
    scalar yields per target, never RustBCA output (CONTRIBUTING.md, Oracles).
    `targets` None: every target with sputter-yield data."""
    import importlib.util
    import datetime as _dt

    spec = importlib.util.spec_from_file_location("oracles_run", HERE.parent / "oracles" / "run.py")
    orc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orc)
    rb = orc.RustBca()
    try:
        rb_bin = rb.binary()
    except orc.NotConfigured as e:
        print(e)
        return
    for target in targets or sorted({g["target"] for g in groups if g["ion"] == "Ar"}):
        rows, settings, mismatches = [], None, None
        e_s = tabulated_e_s("Ar", target, binary)
        for g in groups:
            if (g["ion"], g["target"]) != ("Ar", target):
                continue
            row = {"energy_ev": g["energy_ev"]}
            for k in SPUTTER_WEAK:
                p = sputter_problem("Ar", target, g["energy_ev"], k, ions, e_s)
                # lindhard_settings() reads the lindhard run from RUNS/lindhard/<id>.
                ours = lindhard_cli.run(p, lindhard_cli.RUNS / "lindhard" / p["id"], binary)
                matched = orc.lindhard_settings(p)
                work = lindhard_cli.RUNS / "rustbca" / p["id"]
                work.mkdir(parents=True, exist_ok=True)
                t = rb.run(p, work, rb_bin, matched)
                row[f"lindhard_k{k}"] = ours["sputter_yield"]
                row[f"rustbca_k{k}"] = t["sputter_yield"]
                row[f"rustbca_k{k}_se"] = float(f"{t['std_err']['sputter_yield']:.4g}")
                settings, mismatches = t["settings"], t["mismatches"]
            print(f"RustBCA Ar {g['energy_ev']:g} eV -> {target}: K=0 {row['rustbca_k0']:.3f}, K=3 {row['rustbca_k3']:.3f}")
            rows.append(row)
        if not rows:
            print(f"no Ar -> {target} sputter-yield data; nothing to run")
            continue
        out = {
            "format": "lindhard-oracle-sputter-summary/1",
            "problem": f"Ar -> {target} sputter yield at the level-3 measured energies, settings of ar_1keV_cu_ed_es",
            "oracle": rb.name,
            "oracle_version": rb.version(rb_bin),
            "oracle_commit": rb.commit(),
            "oracle_source": rb.source_url,
            "oracle_license": rb.license_note,
            "lindhard_version": lindhard_cli.version(binary),
            "date": _dt.date.today().isoformat(),
            "ions": ions,
            "oracle_settings": orc._round(settings),
            "mismatches": mismatches,
            "rows": orc._round(rows),
        }
        path = HERE.parent / "oracles" / "summaries" / f"rustbca-sputter_ar_{target.lower()}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {path.relative_to(lindhard_cli.REPO)}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ions", type=int, default=20000)
    ap.add_argument("--check", action="store_true", help="check the datasets' provenance fields only")
    ap.add_argument("--rustbca", action="store_true", help="also run RustBCA at the sputtering energies (context)")
    ap.add_argument("--rustbca-target", action="append", metavar="T",
                    help="with --rustbca: only this target (repeatable; default every target)")
    args = ap.parse_args()

    datasets = load_datasets()
    sputter = load_sputter_datasets()
    if args.check:
        npts = sum(len(d["points"]) for _, d in sputter)
        print(f"{len(datasets)} range dataset(s) and {len(sputter)} sputter-yield dataset(s) ({npts} points) "
              "pass the provenance check")
        return 0
    if not datasets and not sputter:
        print("no datasets; see validation/data/README.md")
        return 0
    binary = lindhard_cli.lindhard_binary()
    groups = energy_groups(sputter)
    if args.rustbca:
        run_rustbca_context(args.ions, binary, groups, args.rustbca_target)
        return 0
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

    sputtering = run_sputter(args.ions, binary, groups) if groups else None
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
                "sputtering": sputtering,
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
