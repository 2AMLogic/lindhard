#!/usr/bin/env python3
"""Level-3 comparison of the fluence-dependent engine with a measured
high-fluence implant profile: 50 keV Sb into Si at 5e16 ions/cm2 (criterion (b)
of #437, #512; Titze et al., SAND2021-13560J, Fig. 3; data and conditions:
validation/data/implant/README.md, titze2021.json).

The measured profile is the SIMS Sb concentration normalized by the authors to
its own maximum, so it gives shape and depth only. Retained dose, absolute
concentration and surface recession are not in the data and are not compared.
The metric is fixed here before the first full run:

* Peak depth. The engine profile is the seed-mean Sb number concentration
  (atoms/cm3) of the final step on a fixed grid of cells `BASE["slab_nm"]` wide
  (the engine's own slabs change thickness with their content and differ
  between seeds; the thickness of the engine slab at the peak is reported); the
  peak is the centre of the maximum cell, so its resolution is the cell width,
  and the cells within `PLATEAU` of the maximum are reported with it. It is
  compared with the printed peak (about 15 nm, titze2021.json) and, separately,
  with the depth of the largest stored marker (the stored 10-19 nm region has
  no markers, so the peak depth cannot be read closer from the data). Pass:
  relative difference within `TOL` of the reference.
* Shape residual. The engine profile divided by its own maximum is
  interpolated linearly in depth between cell centres (flat before the first
  centre) to each stored marker; the residual is the root-mean-square
  difference from the stored normalized concentration over all markers (the
  maximum absolute difference is reported too). Pass: RMS within `TOL`, in
  units of the maximum. The unmeasured 10-19 nm region has no markers, so it
  contributes no residual.

`TOL` is the 25 % preset tolerance of docs/validation.md. Both conventions
(`ideal-mixing`, `fixed-number-density`) are run, with fixed seeds, no `tuning`
and no fitted parameter; pass or fail is reported either way. Sensitivities
(beam angle, native oxide, surface binding, displacement energy, slab width)
are reported side by side and never used to pick a run.

Writes `validation/experiments/sb_si_dynamic_results.json`, which
`validation/update_docs.py` turns into a section of docs/validation.md.

Usage:
    validation/experiments/sb_si_dynamic.py [--threads N]   # full statistics
    validation/experiments/sb_si_dynamic.py --reduced       # CI variant, writes nothing
    validation/experiments/sb_si_dynamic.py --check         # data and committed results only
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
sys.path.insert(0, str(HERE))
import lindhard_cli  # noqa: E402
from cuni_dynamic import mean_se, rel_diff, run_lindhard, verdict  # noqa: E402

DATA = HERE.parent / "data" / "implant"
RESULTS = HERE / "sb_si_dynamic_results.json"
RUNDIR = lindhard_cli.RUNS / "sb-si-dynamic"  # gitignored

N_A = 6.02214076e23  # 1/mol (2019 SI)

TOL = 0.25  # the preset tolerance (docs/validation.md)
PLATEAU = 0.95  # slabs at or above this fraction of the maximum are reported with the peak
PRINTED_PEAK_NM = 15.0  # titze2021.json: "peak near 15 nm"

# Beam and target of the measurement (titze2021.json): 50 keV Sb, 5e16 ions/cm2,
# wafer at 7 degrees to the beam, no azimuth stated, temperature and flux not
# stated (no dynamic annealing exists in the engine).
ENERGY_EV = 50000.0
TILT_DEG = 7.0
FLUENCE_CM2 = 5.0e16
OXIDE_NM = 2.9  # measured native oxide, +- 1.2 nm
OXIDE_DENSITY_G_CM3 = 2.2  # vitreous silica handbook density (examples/as_50keV_si_sio2.toml)
OXIDE_O_ENERGIES = {"e_d_ev": 20.0, "e_s_ev": 2.0}  # illustrative, as in that example

# Matched sputter convention of docs/validation.md section 3 and cuni_dynamic.py:
# E_b = 0 and E_d = E_s, E_s the tabulated cohesive energy (lindhard/src/elements.rs).
# The cutoffs are those of the examples for 50 keV implants.
PHYSICS = {
    "potential": "zbl",
    "stopping": "lindhard-scharff",
    "free_path": "constant",
    "primary_cutoff_ev": 5.0,
    "recoil_cutoff_ev": 2.0,
}
ELEMENTS = {  # symbol: (E_s eV, atomic mass u, solid density g/cm3), as in elements.rs
    "Si": (4.63, 28.085, 2.329),
    "Sb": (2.75, 121.760, 6.691),
}
CONVENTIONS = ["ideal-mixing", "fixed-number-density"]

# Numerical settings of the baseline; the sensitivity runs change one at a time.
# The finite Si layer is far deeper than the profile plus the recession.
BASE = {"slab_nm": 1.0, "max_change": 0.1, "ions_per_step": 1000, "film_nm": 150.0, "tilt_deg": TILT_DEG,
        "oxide": False, "es_scale": 1.0, "ed_ev": None}
SENSITIVITY = [  # (label, change, conventions it applies to)
    ("baseline, seed 1", {}, CONVENTIONS),  # the reference of the single-seed rows below
    ("tilt=0deg", {"tilt_deg": 0.0}, CONVENTIONS),
    ("oxide=2.9nm SiO2", {"oxide": True}, ["ideal-mixing"]),  # one number density cannot describe two materials
    ("es=0.5x", {"es_scale": 0.5}, CONVENTIONS),
    ("es=2x", {"es_scale": 2.0}, CONVENTIONS),
    ("ed=15eV", {"ed_ev": 15.0}, CONVENTIONS),
    ("slab=0.5nm", {"slab_nm": 0.5}, CONVENTIONS),
]
SEEDS = [1, 2, 3, 4]
IONS = 20000  # per seed, standing for FLUENCE_CM2
IONS_REDUCED = 4000


# --- data -------------------------------------------------------------------


def load_markers() -> list[dict]:
    with open(DATA / "titze2021_fig3_5e16.csv", newline="") as f:
        return [{k: (v if k == "row_id" else float(v)) for k, v in r.items()} for r in csv.DictReader(f)]


def stored_peak(markers: list[dict]) -> dict:
    """The stored peak region: the largest marker, and the printed peak depth."""
    top = max(markers, key=lambda m: m["normalized_sb_concentration"])
    depths = sorted(m["depth_nm"] for m in markers)
    gap = max(zip(depths, depths[1:]), key=lambda p: p[1] - p[0])
    return {"marker_max_depth_nm": top["depth_nm"], "marker_max_value": top["normalized_sb_concentration"],
            "printed_peak_nm": PRINTED_PEAK_NM, "largest_gap_nm": list(gap)}


# --- inputs -----------------------------------------------------------------


def si_number_density_cm3() -> float:
    _, mass, rho = ELEMENTS["Si"]
    return rho * N_A / mass


def oxide_o_volume_nm3() -> float:
    """Atomic volume of O in the oxide: (V_SiO2 - V_Si) / 2 per formula unit, V_SiO2 from
    `OXIDE_DENSITY_G_CM3` (60.0843 u) and V_Si the elemental solid volume."""
    v_sio2 = (ELEMENTS["Si"][1] + 2 * 15.999) / (N_A * OXIDE_DENSITY_G_CM3) * 1e21
    v_si = ELEMENTS["Si"][1] / (N_A * ELEMENTS["Si"][2]) * 1e21
    return (v_sio2 - v_si) / 2


def build_input(convention: str, ions: int, seed: int, tilt_deg: float = TILT_DEG, oxide: bool = False,
                es_scale: float = 1.0, ed_ev: float | None = None, slab_nm: float = BASE["slab_nm"],
                max_change: float = BASE["max_change"], ions_per_step: int = BASE["ions_per_step"],
                film_nm: float = BASE["film_nm"]) -> str:
    def energies(sym: str) -> tuple[float, float]:
        es = ELEMENTS[sym][0] * es_scale
        return (es if ed_ev is None else ed_ev), es

    lines = ["# Generated by validation/experiments/sb_si_dynamic.py; do not edit.",
             "[beam]", 'ion = "Sb"', f"energy_ev = {ENERGY_EV!r}", f"tilt_deg = {tilt_deg!r}", ""]
    if oxide:
        o_e = OXIDE_O_ENERGIES
        lines += ["[materials.SiO2]", f"density_g_cm3 = {OXIDE_DENSITY_G_CM3!r}", "elements = [",
                  '  { symbol = "Si", atom_fraction = 1.0 },',
                  f'  {{ symbol = "O", atom_fraction = 2.0, e_d_ev = {o_e["e_d_ev"]!r}, e_s_ev = {o_e["e_s_ev"]!r} }},',
                  "]", ""]
    lines += ["[target]", 'substrate = "Si"', ""]
    if oxide:
        lines += ["[[target.layers]]", 'material = "SiO2"', f"thickness_nm = {OXIDE_NM!r}", ""]
    lines += ["[[target.layers]]", 'material = "Si"', f"thickness_nm = {film_nm!r}", "", "[physics]"]
    for k, v in PHYSICS.items():
        lines.append(f"{k} = {json.dumps(v)}" if isinstance(v, str) else f"{k} = {v!r}")
    ed, es = energies("Si")
    lines += ["", "[physics.energies.Si]", f"e_d_ev = {ed!r}", "e_b_ev = 0.0", f"e_s_ev = {es!r}"]
    lines += ["", "[run]", f"ions = {ions}", f"seed = {seed}", "",
              "[dynamic]", f"fluence_cm2 = {FLUENCE_CM2!r}", f"ions_per_step = {ions_per_step}",
              f"max_change = {max_change!r}", f"slab_nm = {slab_nm!r}", "erosion = true",
              f'relaxation = "{convention}"']
    if convention == "fixed-number-density":
        lines.append(f"number_density_cm3 = {si_number_density_cm3()!r}")
    if oxide:  # O has no elemental solid density: the volume that keeps SiO2 at its density
        lines += ["", "[dynamic.atomic_volume_nm3]", f"O = {oxide_o_volume_nm3()!r}", "",
                  "[dynamic.energies.O]", f'e_d_ev = {OXIDE_O_ENERGIES["e_d_ev"]!r}', "e_b_ev = 0.0",
                  f'e_s_ev = {OXIDE_O_ENERGIES["e_s_ev"]!r}']  # O recoils enter the Si layer
    ed, es = energies("Sb")
    lines += ["", "[dynamic.energies.Sb]", f"e_d_ev = {ed!r}", "e_b_ev = 0.0", f"e_s_ev = {es!r}"]
    return "\n".join(lines) + "\n"


# --- analysis ---------------------------------------------------------------


def read_final_profile(out: Path) -> dict:
    """Sb concentration (atoms/cm3) per slab at the final step, depth from the
    current surface, and the cumulative recession."""
    with open(out / "dynamic_steps.csv", newline="") as f:
        steps = list(csv.DictReader(f))
    last = steps[-1]["step"]
    slabs = []
    with open(out / "dynamic_composition.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["step"] == last:
                t = float(r["thickness_nm"])
                slabs.append({"front_nm": float(r["front_nm"]), "back_nm": float(r["back_nm"]), "thickness_nm": t,
                              "sb_cm3": float(r["atoms_per_cm2_Sb"]) / (t * 1e-7)})
    slabs.sort(key=lambda s: s["front_nm"])
    return {"slabs": slabs, "recession_nm": float(steps[-1]["surface_nm"]),
            "sputter_yield": float(steps[-1]["sputter_yield"]), "steps": len(steps) - 1,
            "ions_done": float(steps[-1]["ions_done"])}


def regrid(slabs: list[dict], extent_nm: float, cell_nm: float = None) -> list[dict]:
    """The slab profile as the mean concentration in fixed cells of `cell_nm`
    from the surface: the engine's slabs change thickness with their content (and
    differ between seeds), so seeds are only comparable on a common grid. A cell
    holds the overlap-weighted concentration of the slabs it covers; depth beyond
    the last slab counts as zero."""
    cell = cell_nm or BASE["slab_nm"]
    n = int(math.ceil(extent_nm / cell))
    out = []
    for i in range(n):
        lo, hi = i * cell, (i + 1) * cell
        atoms = sum(s["sb_cm3"] * max(0.0, min(hi, s["back_nm"]) - max(lo, s["front_nm"])) for s in slabs)
        out.append({"centre_nm": 0.5 * (lo + hi), "thickness_nm": cell, "sb_cm3": atoms / cell})
    return out


def mean_profile(runs: list[list[dict]]) -> list[dict]:
    """Seed mean on the common grid, out to the deepest slab of any seed."""
    extent = max(s["back_nm"] for slabs in runs for s in slabs)
    grids = [regrid(slabs, extent) for slabs in runs]
    return [{"centre_nm": g[0]["centre_nm"], "thickness_nm": g[0]["thickness_nm"],
             "sb_cm3": sum(x["sb_cm3"] for x in col) / len(col)}
            for col in zip(*grids) for g in [col]]


def slab_at(slabs: list[dict], depth_nm: float) -> float:
    """Thickness of the engine slab that contains `depth_nm`."""
    for s in slabs:
        if s["front_nm"] <= depth_nm < s["back_nm"]:
            return s["thickness_nm"]
    return math.nan


def peak(profile: list[dict]) -> dict:
    top = max(s["sb_cm3"] for s in profile)
    if top <= 0.0:
        raise ValueError("no Sb in the engine profile")
    best = max(profile, key=lambda s: s["sb_cm3"])
    plateau = [s["centre_nm"] for s in profile if s["sb_cm3"] >= PLATEAU * top]
    return {"depth_nm": best["centre_nm"], "resolution_nm": best["thickness_nm"], "max_cm3": top,
            "plateau_nm": [min(plateau), max(plateau)]}


def interpolate(profile: list[dict], depth_nm: float) -> float:
    """Linear interpolation of the max-normalized engine profile between slab
    centres; flat before the first centre, zero beyond the last (the profile
    ends in the unimplanted deep film)."""
    top = max(s["sb_cm3"] for s in profile)
    xs, ys = [s["centre_nm"] for s in profile], [s["sb_cm3"] / top for s in profile]
    if depth_nm <= xs[0]:
        return ys[0]
    if depth_nm > xs[-1]:
        return 0.0
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if x0 <= depth_nm <= x1:
            return y0 + (y1 - y0) * (depth_nm - x0) / (x1 - x0)
    return 0.0


def shape_residual(profile: list[dict], markers: list[dict]) -> dict:
    diffs = [interpolate(profile, m["depth_nm"]) - m["normalized_sb_concentration"] for m in markers]
    worst = max(range(len(diffs)), key=lambda i: abs(diffs[i]))
    return {"rms": math.sqrt(sum(d * d for d in diffs) / len(diffs)), "max_abs": abs(diffs[worst]),
            "max_abs_at_depth_nm": markers[worst]["depth_nm"], "markers": len(diffs)}


def compare(profile: list[dict], markers: list[dict]) -> dict:
    pk = peak(profile)
    ref = stored_peak(markers)
    res = shape_residual(profile, markers)
    out = {"peak": pk, "tolerance": TOL, "peak_depth": {}, "shape_residual": {**res, "verdict": verdict(res["rms"], TOL)}}
    for name, depth in (("printed_peak", ref["printed_peak_nm"]), ("marker_maximum", ref["marker_max_depth_nm"])):
        rel = rel_diff(pk["depth_nm"], depth)
        out["peak_depth"][name] = {"reference_nm": depth, "relative_difference": rel, "verdict": verdict(rel, TOL)}
    # Peak inside the unmeasured gap: the data cannot place it closer.
    lo, hi = ref["largest_gap_nm"]
    out["peak_in_unmeasured_gap"] = lo < pk["depth_nm"] < hi
    return out


# --- the runs ---------------------------------------------------------------


class Runner:
    """Runs and analyses one input at a time; the analysed seed result is kept
    in `analysis.json` next to its input and reused when the input text and the
    executable are the same (`reuse`)."""

    def __init__(self, binary: Path, threads: int, root: Path | None = None, reuse: bool = False):
        self.binary, self.threads, self.reuse = binary, threads, reuse
        self.root = root or RUNDIR
        self.binary_id = lindhard_cli.binary_sha256(binary)

    def one(self, tag: str, convention: str, ions: int, seed: int, **over) -> dict:
        text = build_input(convention, ions, seed, **over)
        key = hashlib.sha256((self.binary_id + text).encode()).hexdigest()
        cache = self.root / tag / "analysis.json"
        if self.reuse and cache.exists():
            cached = json.loads(cache.read_text())
            if cached["key"] == key:
                return cached["result"]
        out = run_lindhard(text, self.root / tag, self.binary, self.threads)
        res = read_final_profile(out)
        res["input_sha256"] = hashlib.sha256(text.encode()).hexdigest()
        cache.write_text(json.dumps({"key": key, "result": res}))
        for name in ("dynamic_composition.csv", "dynamic_steps.csv"):
            (out / name).unlink()
        return res

    def case(self, tag: str, convention: str, ions: int, seeds: list[int], markers: list[dict], **over) -> dict:
        runs = [self.one(f"{tag}_s{s}", convention, ions, s, **over) for s in seeds]
        profile = mean_profile([r["slabs"] for r in runs])
        res = compare(profile, markers)
        res["seeds"] = seeds
        res["ions_per_seed"] = ions
        res["recession_nm"] = mean_se([r["recession_nm"] for r in runs])[0]
        res["sputter_yield"] = mean_se([r["sputter_yield"] for r in runs])[0]
        res["peak_per_seed_nm"] = [peak(mean_profile([r["slabs"]]))["depth_nm"] for r in runs]
        res["engine_slab_thickness_at_peak_nm"] = [slab_at(r["slabs"], res["peak"]["depth_nm"]) for r in runs]
        res["profile"] = [{"centre_nm": s["centre_nm"], "sb_cm3": s["sb_cm3"]} for s in profile]
        res["input_sha256"] = [r["input_sha256"] for r in runs]
        return res


def command_line(argv_extra: str = "") -> str:
    return "python3 validation/experiments/sb_si_dynamic.py" + argv_extra


def settings() -> dict:
    return {
        "physics": PHYSICS, "base_numerics": BASE, "energy_ev": ENERGY_EV, "tilt_deg": TILT_DEG,
        "fluence_cm2": FLUENCE_CM2, "seeds": SEEDS, "ions_per_seed": IONS,
        "surface_binding_ev": {s: e[0] for s, e in ELEMENTS.items()},
        "displacement_energy_ev": "E_d = E_s (matched sputter convention)",
        "oxide": {"thickness_nm": OXIDE_NM, "density_g_cm3": OXIDE_DENSITY_G_CM3, "o_energies": OXIDE_O_ENERGIES},
        "tolerance": TOL, "plateau_fraction": PLATEAU, "printed_peak_nm": PRINTED_PEAK_NM,
        "si_number_density_cm3": si_number_density_cm3(),
        "sensitivity": [{"label": lab, "change": ch, "conventions": cv} for lab, ch, cv in SENSITIVITY],
        "not_modelled": ["temperature and beam flux (not stated in the paper): no dynamic annealing",
                         "azimuth of the 7 degree tilt (not stated)", "crystal channelling (amorphous target)"],
    }


def full_run(threads: int, reuse: bool) -> dict:
    binary = lindhard_cli.lindhard_binary()
    runner = Runner(binary, threads, reuse=reuse)
    markers = load_markers()
    res = {
        "format": {"name": "lindhard-sb-si-dynamic-results", "version": 1},
        "issue": "#512",
        "reproduce": command_line(),
        "binary_sha256": lindhard_cli.binary_sha256(binary),
        "software": lindhard_cli.version(binary),
        "settings": settings(),
        "stored": stored_peak(markers),
        "runs": {}, "sensitivity": {},
    }
    for conv in CONVENTIONS:
        res["runs"][conv] = runner.case(f"base_{conv}", conv, IONS, SEEDS, markers)
        for label, change, convs in SENSITIVITY:
            if conv in convs:
                tag = "sens_" + "".join(ch if ch.isalnum() else "-" for ch in label) + f"_{conv}"
                res["sensitivity"].setdefault(conv, {})[label] = runner.case(tag, conv, IONS, SEEDS[:1], markers, **change)
    return res


def reduced_run(threads: int, root: Path | None = None) -> dict:
    """CI variant: both conventions, one seed, a fifth of the ions. Writes nothing
    under the committed results."""
    runner = Runner(lindhard_cli.lindhard_binary(), threads, root=root)
    markers = load_markers()
    return {conv: runner.case(f"reduced_{conv}", conv, IONS_REDUCED, [1], markers) for conv in CONVENTIONS}


# --- checks -----------------------------------------------------------------


def check_data(results: bool = True) -> list[str]:
    errs = []
    m = load_markers()
    if len(m) != 53:
        errs.append(f"expected 53 Titze markers, found {len(m)}")
    elif [x["depth_nm"] for x in m] != sorted(x["depth_nm"] for x in m):
        errs.append("Titze markers are not in depth order")
    if results and RESULTS.exists():
        r = json.loads(RESULTS.read_text())
        if r["settings"]["tolerance"] != TOL:
            errs.append("committed results were produced with a different tolerance")
        if r["stored"] != json.loads(json.dumps(stored_peak(m))):
            errs.append("committed stored-peak values differ from the data files; rerun the harness")
        if set(r["runs"]) != set(CONVENTIONS):
            errs.append("committed results do not hold both relaxation conventions")
    return errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=2, help="worker threads of each run (never changes results)")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse the analysed result of a run whose input and executable are unchanged")
    ap.add_argument("--reduced", action="store_true", help="CI variant; prints a summary and writes nothing")
    ap.add_argument("--check", action="store_true", help="check the data files and committed results only")
    args = ap.parse_args()
    errs = check_data(results=args.check)
    if errs:
        for e in errs:
            print(f"error: {e}", file=sys.stderr)
        return 1 if args.check else 2
    if args.check:
        return 0
    if args.reduced:
        with tempfile.TemporaryDirectory() as d:
            res = reduced_run(args.threads, Path(d))
        for conv, r in res.items():
            r.pop("profile")
        print(json.dumps(res, indent=1))
        return 0
    res = full_run(args.threads, args.reuse)
    RESULTS.write_text(json.dumps(res, indent=1) + "\n")
    print(f"wrote {RESULTS.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
