#!/usr/bin/env python3
"""Level-3 comparison of the simulated secondary-electron yield delta(E)
(escaping electrons below 50 eV, per incident electron, normal incidence,
100 eV to 5 keV) with the measurements compiled by D. C. Joy (issue #149;
data and schema: validation/data/README.md, "Secondary-electron yields").

Standard library only.

    se_yield.py --datasets         datasets and provenance rows only (validation/run.sh)
    se_yield.py --check            the above, and the initial bounds against the
                                   committed results (exit 1 if a bound fails)
    se_yield.py --run [--material Al Cu Au] [--configs ID ...] [--histories N]
                                   [--energies E ...]
                                   run `lindhard` (LINDHARD_BIN or a release
                                   build) and write results into
                                   validation/experiments/se_yield_results.json
                                   (existing runs of other keys are kept)
    se_yield.py --markdown         print the tables spliced into
                                   docs/validation.md

Reference statistics (the "measured median", written down before looking at
any simulated value; the compared-range clause was added afterwards, see
WINDOW_EV and docs/validation.md):

* every stored data set is one original measurement; they are never merged;
* a data set *resolves* the maximum if it has at least 5 points, its
  largest yield is not at its lowest or highest energy (the first such point
  in energy order if the maximum is tied), and that point lies within the
  compared range, 100 eV to 5 keV (WINDOW_EV). Its delta_max, E_max are that
  tabulated point: nothing is fitted or interpolated;
* the measured median is the median over the resolving data sets of delta_max
  and, separately, of E_max (median of two = their mean). The spread (min and
  max over the same sets) is reported with it. Data sets that do not resolve
  the maximum are listed, not used.

Initial bounds (issue #149, gated for the default configuration of every
material that was run): E_max within a factor 2 of the measured median's, and
delta_max within 50 % of the measured median's. E_max is the grid energy of the
largest simulated delta, so its resolution is the energy grid. Every other
point is reported, not gated. Loosening a bound needs operator sign-off.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

REPO = lindhard_cli.REPO
DATA = HERE.parent / "data" / "se_yield"
OPTICAL = HERE.parent / "data" / "optical"
RESULTS = HERE / "se_yield_results.json"
PROVENANCE = REPO / "docs" / "data-provenance.md"

MATERIALS = ["Al", "Cu", "Si", "Au"]
# Materials that can be run: an optical ELF with provenance is committed. The Si ELF (Yang et al. 2019, #125)
# ends at 199 eV: the Si K shell and the L-shell tail above 199 eV are absent (docs/data-provenance.md).
OPTICAL_ELF = {"Al": "al_elf_hagemann1975.toml", "Cu": "cu_elf_hagemann1975.toml",
               "Si": "si_elf_yang2019.toml", "Au": "au_elf_hagemann1975.toml"}
MIN_POINTS = 5
# The energy range the issue compares (normal incidence, 100 eV to 5 keV). A data set whose interior maximum
# lies outside it is not measuring the low-energy maximum of delta(E) and does not resolve it (amended in #149
# before review: as first written the rule let se_al_czaja1966, measured only from 10 to 40 keV, count with a
# "maximum" at 15 keV; docs/validation.md gives the effect of the amendment).
WINDOW_EV = (100.0, 5000.0)

# Initial bounds (issue #149).
EMAX_FACTOR = 2.0
DMAX_REL = 0.5

ENERGIES_EV = [100, 150, 200, 300, 400, 600, 800, 1000, 1500, 2000, 3000, 5000]
SLOW_ENERGIES_EV = [100, 200, 400, 800, 1500, 3000]  # Mermin: table builds take minutes
# Full Penn: not run for the committed results. One Al run at 200 eV (the smallest table) had not finished
# its table build after 45 minutes on two threads of a shared host (2026-10-08) and was stopped.
PENN_FULL_ENERGIES_EV = [200, 800]

# Band inputs. Metals (Al, Cu, Au): caller-supplied. #115 added cited band defaults
# (lindhard::electron::boundary::BAND_DEFAULTS) but found no openable source for the work function or Fermi
# energy of Al, Cu or Au, so these stay as they are. Si: an insulator band; the band gap and the electron
# affinity are the cited BAND_DEFAULTS values of #115 (Robertson and Wallace 2015), and the valence band width,
# which BAND_DEFAULTS lacks, is the computed Gamma25' - Gamma1 separation of Chelikowsky and Cohen (1974), whose
# Table IV lists two photoemission measurements that agree with it within their stated errors. See
# docs/data-provenance.md, "Band inputs of the SE-yield validation runs", "Cited band defaults" and
# "Per-material barrier parameters".
BAND = {
    "Al": {
        "kind": "free-electron-metal",
        "valence_electrons_per_atom": 3.0,
        "work_function_ev": {"low": 4.06, "mid": 4.16, "high": 4.26},
    },
    "Cu": {
        "kind": "free-electron-metal",
        "valence_electrons_per_atom": 1.0,
        "work_function_ev": {"low": 4.53, "mid": 4.815, "high": 5.10},
    },
    "Si": {
        "kind": "insulator",
        "valence_band_width_ev": 12.36,
        "band_gap_ev": 1.1,
        "affinity_ev": 4.05,
    },
    "Au": {
        "kind": "free-electron-metal",
        "valence_electrons_per_atom": 1.0,
        "work_function_ev": {"low": 5.10, "mid": 5.285, "high": 5.47},
    },
}
METAL_BAND_PROVENANCE = (
    "Caller-supplied for the delta(E) validation (#149); #115 found no openable source for these, so they stay: "
    "work function from the range printed "
    "for the element in the Wikipedia 'Work function' table, revision 1368612509 of 2026-08-10, which cites "
    "CRC Handbook of Chemistry and Physics (2008), p. 12-124 (the Handbook itself was not opened); "
    "'low'/'high' are the table's endpoints and 'mid' their mean; valence electrons per atom from the "
    "ground-state configuration (Wikipedia 'Electron configuration', Al [Ne]3s2 3p1 -> 3, Cu [Ar]3d10 4s1 -> 1; Au [Xe]4f14 5d10 6s1 -> 1, from the table of anomalous configurations of the revision 1378968605 of 2026-10-07 of that article); "
    "Fermi energy by the library's free-electron formula (Verduin 2017, Eq. 3.133) from the density of the element table."
)
SI_BAND_PROVENANCE = (
    "Si insulator band for the delta(E) validation (#149). Band gap 1.1 eV and electron affinity 4.05 eV: the cited "
    "lindhard::electron::boundary::BAND_DEFAULTS values (#115), J. Robertson and R. M. Wallace, Mater. Sci. Eng. R 88, "
    "1-41 (2015), doi:10.1016/j.mser.2014.11.001, authors' manuscript hdl:1810/246441, Table 1, p. 43 (gap) and Fig. 42, "
    "p. 63 (affinity). Valence band width 12.36 eV: Gamma1 at -12.36 eV below Gamma25' (the valence band top, 0.00 eV) "
    "in the energy-dependent nonlocal pseudopotential calculation of J. R. Chelikowsky and M. L. Cohen, 'Electronic "
    "structure of silicon', Phys. Rev. B 10, 5095 (1974), doi:10.1103/PhysRevB.10.5095, Table II (and the 'Nonlocal' "
    "column of Table IV), read in the preprint LBL-3127 (escholarship.org/uc/item/4jf90532, preprint pp. 38 and 40), "
    "opened 2026-10-09; a published computed value. Its Table IV lists, for the same level, the photoemission values "
    "-12.4 +- 0.6 eV (W. D. Grobman and D. E. Eastman, Phys. Rev. Lett. 29, 1508 (1972)) and -12.5 +- 0.6 eV (L. Ley et "
    "al., Phys. Rev. Lett. 29, 1088 (1972)); those two papers were not opened as published (closed access; the Ley et "
    "al. preprint LBL-688 gives Gamma1 relative to the Fermi level only)."
)
BAND_PROVENANCE = {"Al": METAL_BAND_PROVENANCE, "Cu": METAL_BAND_PROVENANCE, "Au": METAL_BAND_PROVENANCE,
                   "Si": SI_BAND_PROVENANCE}

# Elastic scattering of every delta(E) run (#149 remaining scope 1a). Mott partial waves on the DHFS potential of
# Salvat et al. (1987) (#130) with the Furness-McCarthy exchange correction, as in
# validation/experiments/backscatter/eta_*.toml. The correlation-polarization correction of those inputs is OFF
# here: with Seltzer's cutoff rule, the only one with a source in the library (Salvat 2003, Eq. (5)), it needs
# every table energy above 50 eV, and delta(E) needs elastic tables down to 5 eV; no cited b_pol^2 for lower
# energies is available, and choosing one would be a tuned constant.
ELASTIC_POTENTIAL = "salvat-dhfs"
ELASTIC_EXCHANGE = True
ELASTIC_ID = "mott/salvat-dhfs/exchange"

# id -> (label, inelastic model, work function key or None for barrier off, cutoff reference, energies).
# "Barrier on" is the inner-potential step of the material's band: work function above the Fermi level for the
# metals (the key picks the value in BAND), the electron affinity above the conduction band for Si (one cited
# value, so the key does not apply to it).
CONFIGS = {
    "default": ("single-pole Penn, barrier on (metals: mid work function; Si: electron affinity), vacuum-level cutoff", "penn-single-pole", "mid", "vacuum-level", ENERGIES_EV),
    "phi-low": ("as default, work function at the low end of the cited range (metals only)", "penn-single-pole", "low", "vacuum-level", ENERGIES_EV),
    "phi-high": ("as default, work function at the high end of the cited range (metals only)", "penn-single-pole", "high", "vacuum-level", ENERGIES_EV),
    "barrier-off": ("as default, transparent boundary (no barrier), vacuum-level cutoff", "penn-single-pole", None, "vacuum-level", ENERGIES_EV),
    "cutoff-band-bottom": ("as default, cutoff measured from the band bottom, 1 eV above the Fermi level (below the vacuum level)", "penn-single-pole", "mid", "band-bottom", ENERGIES_EV),
    "penn-full": ("full Penn, barrier on, vacuum-level cutoff", "penn-full", "mid", "vacuum-level", PENN_FULL_ENERGIES_EV),
    "mermin": ("Mermin (MELF), barrier on, vacuum-level cutoff", "mermin-melf", "mid", "vacuum-level", SLOW_ENERGIES_EV),
}
DEFAULT_CONFIG = "default"
# The work-function range configurations apply to the metals only (Si's barrier is its one cited affinity).
METAL_ONLY_CONFIGS = ("phi-low", "phi-high")


def applies(material: str, cfg_id: str) -> bool:
    """Whether configuration `cfg_id` is defined for `material`."""
    return not (cfg_id in METAL_ONLY_CONFIGS and BAND[material]["kind"] != "free-electron-metal")


# Threshold of the `cutoff-band-bottom` configuration, eV above the band bottom. The library refuses a
# band-bottom threshold at or below the Fermi energy with secondaries on, so this is the Fermi energy the
# library computes from the band inputs (printed as `band.model.fermi_ev` in every run's metadata: Al
# 11.6555 eV, Cu 7.0445 eV, Au 5.5269 eV; Si, an insulator, W_v + E_g/2 = 12.91 eV) plus 1 eV, rounded up to
# 0.01 eV. It lies below the vacuum level, so unlike the default (threshold 1 eV above the vacuum level) every
# electron that can leave is followed. `run_one` checks it against the run's own Fermi energy.
BAND_BOTTOM_CUTOFF_EV = {"Al": 12.66, "Cu": 8.05, "Si": 13.91, "Au": 6.53}


# --- reference data -----------------------------------------------------------------------------

def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


REQUIRED_TEXT = ["id", "target", "target_state", "original_reference", "compilation", "compilation_location",
                 "extraction", "terms", "added", "url"]


def check_dataset(path: Path, d: dict, provenance: str) -> list[str]:
    errs = []
    for k in REQUIRED_TEXT:
        if not isinstance(d.get(k), str) or not d[k].strip():
            errs.append(f"`{k}` missing or empty")
    if d.get("kind") != "se_yield":
        errs.append("`kind` must be 'se_yield'")
    if d.get("id") != path.stem:
        errs.append(f"`id` {d.get('id')!r} differs from the file name {path.stem!r}")
    if d.get("target") not in MATERIALS:
        errs.append(f"`target` {d.get('target')!r} is not one of {MATERIALS}")
    if "incidence_deg" not in d:
        errs.append("`incidence_deg` missing (null = not stated by the source)")
    elif d["incidence_deg"] is not None and not _num(d["incidence_deg"]):
        errs.append("`incidence_deg` must be a number or null")
    elif d["incidence_deg"] is None and not (isinstance(d.get("incidence_note"), str) and d["incidence_note"].strip()):
        errs.append("`incidence_deg` is null without an `incidence_note`")
    if not isinstance(d.get("compilation_reference_number"), int):
        errs.append("`compilation_reference_number` must be an integer")
    pts = d.get("points")
    if not isinstance(pts, list) or not pts:
        errs.append("no points")
    else:
        for i, p in enumerate(pts):
            for k in ("energy_ev", "yield"):
                if not (_num(p.get(k)) and p[k] > 0):
                    errs.append(f"`points[{i}].{k}` must be a positive number")
    if isinstance(d.get("id"), str) and d["id"] not in provenance:
        errs.append(f"no row in docs/data-provenance.md names `{d['id']}`")
    return errs


def load_datasets() -> list[dict]:
    provenance = PROVENANCE.read_text()
    out, failed = [], False
    for f in sorted(DATA.glob("*.json")) if DATA.is_dir() else []:
        rel = f.relative_to(REPO)
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
        print("error: SE-yield dataset check failed (rules: validation/data/README.md)", file=sys.stderr)
        sys.exit(2)
    return out


def dataset_peak(d: dict):
    """(E_max_ev, delta_max) if the data set resolves its maximum, else None."""
    pts = sorted((p["energy_ev"], p["yield"]) for p in d["points"])
    if len(pts) < MIN_POINTS:
        return None
    m = max(y for _, y in pts)
    i = next(k for k, (_, y) in enumerate(pts) if y == m)
    if i == 0 or i == len(pts) - 1:
        return None
    if not WINDOW_EV[0] <= pts[i][0] <= WINDOW_EV[1]:
        return None
    return pts[i]


def reference_stats(datasets: list[dict], material: str) -> dict:
    sets = [d for d in datasets if d["target"] == material]
    res, unres = [], []
    for d in sorted(sets, key=lambda d: d["id"]):
        pk = dataset_peak(d)
        (res if pk else unres).append((d, pk))
    out = {"material": material, "n_sets": len(sets), "resolving": [], "not_resolving": [d["id"] for d, _ in unres]}
    for d, pk in res:
        out["resolving"].append({"id": d["id"], "emax_ev": pk[0], "dmax": pk[1]})
    if res:
        e = [pk[0] for _, pk in res]
        y = [pk[1] for _, pk in res]
        out.update(emax_median_ev=statistics.median(e), emax_min_ev=min(e), emax_max_ev=max(e),
                   dmax_median=statistics.median(y), dmax_min=min(y), dmax_max=max(y))
    return out


def curve_at(datasets: list[dict], material: str, e_lo=100.0, e_hi=5000.0) -> list[tuple[float, float, str]]:
    return sorted((p["energy_ev"], p["yield"], d["id"]) for d in datasets if d["target"] == material
                  for p in d["points"] if e_lo <= p["energy_ev"] <= e_hi)


# --- simulation -------------------------------------------------------------------------------

def cutoff_ev(material: str, cfg_id: str) -> float:
    """The `cutoff_ev` of a run: 1 eV above the vacuum level, or BAND_BOTTOM_CUTOFF_EV for band-bottom."""
    return BAND_BOTTOM_CUTOFF_EV[material] if CONFIGS[cfg_id][3] == "band-bottom" else 1.0


def band_line(material: str, phi: str | None) -> str:
    """The `band = { ... }` entry of `material`'s run input. `phi` picks a metal's work function; with the
    barrier off (`None`) the band still supplies the secondary-electron binding (Kieft-Bosch needs it), so the
    mid work function is used there and is not applied as a step."""
    band = BAND[material]
    prov = BAND_PROVENANCE[material].replace('"', "'")
    if band["kind"] == "free-electron-metal":
        w = band["work_function_ev"][phi or "mid"]
        return (f'band = {{ kind = "free-electron-metal", valence_electrons_per_atom = {band["valence_electrons_per_atom"]!r}, '
                f'work_function_ev = {w!r}, provenance = "{prov}" }}')
    if band["kind"] == "insulator":
        return (f'band = {{ kind = "insulator", valence_band_width_ev = {band["valence_band_width_ev"]!r}, '
                f'band_gap_ev = {band["band_gap_ev"]!r}, affinity_ev = {band["affinity_ev"]!r}, provenance = "{prov}" }}')
    raise ValueError(f"unknown band kind {band['kind']!r} for {material}")


def make_input(material: str, cfg_id: str, energy_ev: float, histories: int, seed: int, elf_name: str) -> str:
    _label, model, phi, cutoff_ref, _energies = CONFIGS[cfg_id]
    if not applies(material, cfg_id):
        raise ValueError(f"configuration {cfg_id!r} does not apply to {material}")
    lines = [
        "# Generated by validation/experiments/se_yield.py; do not edit.",
        "[electron.beam]",
        f"energy_ev = {float(energy_ev)!r}",
        "[electron.transport]",
        f"cutoff_ev = {cutoff_ev(material, cfg_id)!r}",
        f'cutoff_reference = "{cutoff_ref}"',
        'secondaries = "kieft-bosch"',
        f'boundary = "{"step-barrier" if phi else "transparent"}"',
        "[electron.elastic]",
        'model = "mott"',
        f'potential = "{ELASTIC_POTENTIAL}"',
        f"exchange = {'true' if ELASTIC_EXCHANGE else 'false'}",
        "[electron.inelastic]",
        f'model = "{model}"',
        "[electron.tables]",
        "min_energy_ev = 5.0",
        "points_per_decade = 10.0",
        f"[electron.materials.{material}]",
        f'optical_elf = "{elf_name}"',
        band_line(material, phi),
        "[target]",
        f'substrate = "{material}"',
        "[run]",
        f"histories = {int(histories)}",
        f"seed = {int(seed)}",
        "threads = 2",
        "",
    ]
    return "\n".join(lines)


def band_fermi_ev(model: dict) -> float:
    """The Fermi energy above the band bottom of a run's band metadata (`physics.materials[].band.model`): printed
    for a metal; for an insulator the library's `BandStructure::fermi_ev`, W_v + E_g/2 (mid-gap), which the
    metadata does not print."""
    if "fermi_ev" in model:
        return model["fermi_ev"]
    if model.get("kind") == "insulator":
        return model["valence_band_width_ev"] + 0.5 * model["band_gap_ev"]
    raise KeyError(f"no Fermi energy in band model {model!r}")


def run_one(binary: Path, material: str, cfg_id: str, energy_ev: float, histories: int, seed: int) -> dict:
    elf = OPTICAL / OPTICAL_ELF[material]
    with tempfile.TemporaryDirectory(prefix="se_yield_", dir=lindhard_cli.RUNS if lindhard_cli.RUNS.is_dir() else None) as td:
        td = Path(td)
        (td / elf.name).write_bytes(elf.read_bytes())
        (td / "input.toml").write_text(make_input(material, cfg_id, energy_ev, histories, seed, elf.name))
        proc = subprocess.run([str(binary), "run", str(td / "input.toml"), "--out", str(td / "out")],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            sys.exit(f"error: lindhard failed ({material} {cfg_id} {energy_ev} eV):\n{proc.stderr}")
        s = json.loads((td / "out" / "electron_summary.json").read_text())
    r = s["results"]
    n = r["histories"]
    fermi = band_fermi_ev(s["physics"]["materials"][0]["band"]["model"])
    cut = cutoff_ev(material, cfg_id)
    # The 1e-9 eV slack absorbs the rounding of W_v + E_g/2 (Si: 12.36 + 0.55 prints as 12.910000000000002).
    if CONFIGS[cfg_id][3] == "band-bottom" and not -1e-9 <= cut - (fermi + 1.0) < 0.01:
        sys.exit(f"error: BAND_BOTTOM_CUTOFF_EV[{material}] = {cut} is not the run's Fermi energy {fermi} + 1 eV")
    slow, fast = r["front"]["slow"], r["front"]["fast"]
    return {
        "material": material, "config": cfg_id, "energy_ev": float(energy_ev), "histories": n, "seed": seed,
        "elastic": ELASTIC_ID,
        "delta": slow["per_primary"], "eta": fast["per_primary"],
        # sqrt(N_slow)/N: the Poisson floor of the statistical error; it ignores the correlation of
        # electrons of one cascade and so understates it.
        "delta_poisson_floor": math.sqrt(slow["count"]) / n if n else None,
        "energy_balance_relative_imbalance": r["budget"]["relative_imbalance"],
        "fermi_ev": fermi, "cutoff_ev": cut, "cutoff_reference": CONFIGS[cfg_id][3],
    }


def do_run(args) -> int:
    binary = lindhard_cli.lindhard_binary()
    prev = json.loads(RESULTS.read_text()) if RESULTS.exists() else {"runs": []}
    # Runs with another elastic model are dropped, not mixed into the tables (#149: the DHFS rerun replaces
    # the Thomas-Fermi Yukawa stand-in runs).
    keep = {(r["material"], r["config"], r["energy_ev"]): r for r in prev["runs"] if r.get("elastic") == ELASTIC_ID}
    for m in args.material:
        for cid in args.configs:
            if not applies(m, cid):
                continue
            for e in CONFIGS[cid][4]:
                if args.energies and e not in args.energies:
                    continue
                r = run_one(binary, m, cid, e, args.histories, 1)
                keep[(m, cid, float(e))] = r
                print(f"{m} {cid} {e} eV: delta {r['delta']:.4f} eta {r['eta']:.4f}", flush=True)
                RESULTS.write_text(json.dumps(assemble(binary, keep.values()), indent=1, sort_keys=True) + "\n")
    return 0


def assemble(binary: Path, runs) -> dict:
    return {
        "format": {"name": "lindhard-se-yield-results", "version": 1},
        "software": lindhard_cli.version(binary),
        "settings": {
            "se_split_ev": 50, "incidence": "normal (the simulation's; the compilation does not state its sets')",
            "configs": {k: {"label": v[0], "inelastic": v[1], "work_function": v[2], "cutoff_reference": v[3]} for k, v in CONFIGS.items()},
            "band_bottom_cutoff_ev": BAND_BOTTOM_CUTOFF_EV,
            "elastic": ("Mott partial waves on the Salvat et al. (1987) DHFS potential with the Furness-McCarthy exchange "
                        "correction; correlation-polarization off (Seltzer's b_pol^2 rule needs table energies above "
                        "50 eV, the delta(E) tables start at 5 eV)"),
            "elastic_id": ELASTIC_ID,
            "band_inputs": BAND, "band_provenance": BAND_PROVENANCE,
        },
        "runs": sorted(runs, key=lambda r: (r["material"], r["config"], r["energy_ev"])),
    }


# --- reporting and bounds ---------------------------------------------------------------------------

def load_results():
    return json.loads(RESULTS.read_text()) if RESULTS.exists() else None


def sim_curve(results, material: str, cfg: str):
    return sorted((r["energy_ev"], r["delta"]) for r in results["runs"] if r["material"] == material and r["config"] == cfg)


def sim_peak(curve):
    if not curve:
        return None
    m = max(y for _, y in curve)
    return next((e, y) for e, y in curve if y == m)


def evaluate_bounds(ref: dict, curve) -> dict | None:
    """The two initial bounds for one material's default curve against its reference statistics."""
    pk = sim_peak(curve)
    if pk is None or "emax_median_ev" not in ref:
        return None
    e_ratio = pk[0] / ref["emax_median_ev"]
    d_rel = pk[1] / ref["dmax_median"] - 1.0
    return {
        "emax_ev": pk[0], "dmax": pk[1], "emax_ratio": e_ratio, "dmax_rel": d_rel,
        "emax_ok": 1.0 / EMAX_FACTOR <= e_ratio <= EMAX_FACTOR,
        "dmax_ok": abs(d_rel) <= DMAX_REL,
    }


def fmt(x, nd=3):
    return "-" if x is None else f"{x:.{nd}g}"


def markdown(datasets, results) -> str:
    out = []
    out.append("**Reference data** (Joy's compilation, one row per original measurement; the maximum is the tabulated point, "
               "not a fit; a set resolves the maximum if it has at least 5 points and the maximum is interior to it and "
               "lies in the compared range, 100 eV to 5 keV):\n")
    out.append("| Material | Sets | Resolving | Measured median E_max (eV) [min, max] | Measured median δ_max [min, max] | Not resolving (listed, not used) |")
    out.append("|---|---|---|---|---|---|")
    for m in MATERIALS:
        s = reference_stats(datasets, m)
        if "emax_median_ev" in s:
            em = f"{s['emax_median_ev']:.0f} [{s['emax_min_ev']:.0f}, {s['emax_max_ev']:.0f}]"
            dm = f"{s['dmax_median']:.3f} [{s['dmax_min']:.3f}, {s['dmax_max']:.3f}]"
        else:
            em = dm = "-"
        res = ", ".join(f"`{r['id']}` ({r['emax_ev']:.0f} eV, {r['dmax']:.3f})" for r in s["resolving"]) or "none"
        out.append(f"| {m} | {s['n_sets']} | {len(s['resolving'])}: {res} | {em} | {dm} | {len(s['not_resolving'])} sets |")
    out.append("")
    if not results:
        out.append("_No simulated results are committed yet._")
        return "\n".join(out)
    out.append(f"**Simulated δ(E)** (`{results['software']}`; seed 1; histories per run in the results file; "
               "δ = electrons escaping the front face below 50 eV per primary). Statistical error: the Poisson floor "
               "√N_slow/N is in the results file and understates the true error by the cascade correlation.\n")
    for m in MATERIALS:
        cfgs = [c for c in CONFIGS if sim_curve(results, m, c)]
        if not cfgs:
            out.append(f"**{m}**: not run (the reason is given in the text below the tables).\n")
            continue
        energies = sorted({e for c in cfgs for e, _ in sim_curve(results, m, c)})
        out.append(f"**{m}**\n")
        out.append("| Configuration | " + " | ".join(f"{e:g}" for e in energies) + " | E_max (eV) | δ_max |")
        out.append("|---|" + "---|" * (len(energies) + 2))
        for c in cfgs:
            cur = dict(sim_curve(results, m, c))
            pk = sim_peak(sorted(cur.items()))
            cells = " | ".join(f"{cur[e]:.3f}" if e in cur else "" for e in energies)
            missing = [e for e in CONFIGS[c][4] if float(e) not in cur]
            ran = sorted(cur)
            if missing and pk[0] in (ran[0], ran[-1]):
                # The largest value of an incomplete curve sits at the end of what ran: no maximum to report.
                note = "incomplete: not run at " + ", ".join(f"{e:g}" for e in missing) + " eV"
                out.append(f"| `{c}`: {CONFIGS[c][0]} | {cells} | {note} | - |")
            elif missing:
                note = "; not run at " + ", ".join(f"{e:g}" for e in missing) + " eV"
                out.append(f"| `{c}`: {CONFIGS[c][0]} | {cells} | {pk[0]:g}{note} | {pk[1]:.3f} |")
            else:
                out.append(f"| `{c}`: {CONFIGS[c][0]} | {cells} | {pk[0]:g} | {pk[1]:.3f} |")
        out.append("")
        unrun = [c for c in CONFIGS if c not in cfgs and applies(m, c)]
        if unrun:
            out.append("Not run for " + m + ": " + "; ".join(f"`{c}` ({CONFIGS[c][0]})" for c in unrun)
                       + ". The reason is given in the text below the tables.\n")
        na = [c for c in CONFIGS if not applies(m, c)]
        if na:
            out.append("Not defined for " + m + ": " + ", ".join(f"`{c}`" for c in na)
                       + " (work-function ranges of the metals; the barrier of " + m + " is its cited electron affinity).\n")
    out.append("**Initial bounds** (default configuration; gated by `validation/experiments/se_yield.py --check`):\n")
    out.append("| Material | Simulated E_max (eV) | Measured median E_max (eV) | Ratio | Within factor 2 | Simulated δ_max | Measured median δ_max | Deviation | Within 50 % |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for m in MATERIALS:
        ref = reference_stats(datasets, m)
        b = evaluate_bounds(ref, sim_curve(results, m, DEFAULT_CONFIG))
        if b is None:
            out.append(f"| {m} | not run | {fmt(ref.get('emax_median_ev'), 4)} | - | - | not run | {fmt(ref.get('dmax_median'))} | - | - |")
        else:
            out.append(f"| {m} | {b['emax_ev']:g} | {ref['emax_median_ev']:g} | {b['emax_ratio']:.2f} | {'yes' if b['emax_ok'] else '**NO**'} "
                       f"| {b['dmax']:.3f} | {ref['dmax_median']:.3f} | {100 * b['dmax_rel']:+.0f} % | {'yes' if b['dmax_ok'] else '**NO**'} |")
    return "\n".join(out)


def do_check(strict_bounds: bool = True, bounds: bool = True) -> int:
    datasets = load_datasets()
    if not datasets:
        print("error: no SE-yield datasets", file=sys.stderr)
        return 2
    print(f"{len(datasets)} SE-yield data sets pass the provenance check")
    results = load_results()
    if not results or not bounds:
        return 0
    failed = False
    for m in MATERIALS:
        b = evaluate_bounds(reference_stats(datasets, m), sim_curve(results, m, DEFAULT_CONFIG))
        if b is None:
            continue
        for name, ok in (("E_max within a factor of 2 of the measured median", b["emax_ok"]),
                         ("delta_max within 50 % of the measured median", b["dmax_ok"])):
            print(f"{m}: {name}: {'ok' if ok else 'FAILS'} (E_max {b['emax_ev']:g} eV, ratio {b['emax_ratio']:.2f}; "
                  f"delta_max {b['dmax']:.3f}, {100 * b['dmax_rel']:+.0f} %)")
            failed |= not ok
    return 1 if (failed and strict_bounds) else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--datasets", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--material", nargs="+", default=list(OPTICAL_ELF), choices=list(OPTICAL_ELF))
    ap.add_argument("--configs", nargs="+", default=list(CONFIGS), choices=list(CONFIGS))
    ap.add_argument("--histories", type=int, default=2000)
    ap.add_argument("--results", type=Path, default=None,
                    help="with --run: read and write this results file instead of the committed one (to run two "
                         "workers side by side; merge the files by hand afterwards)")
    ap.add_argument("--energies", nargs="+", type=float, help="only these grid energies, eV (to resume a partial run)")
    args = ap.parse_args()
    if args.run:
        if args.results:
            global RESULTS
            RESULTS = args.results
        return do_run(args)
    if args.markdown:
        print(markdown(load_datasets(), load_results()))
        return 0
    return do_check(strict_bounds=not args.datasets, bounds=not args.datasets)


if __name__ == "__main__":
    sys.exit(main())
