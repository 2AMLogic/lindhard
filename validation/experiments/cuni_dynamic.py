#!/usr/bin/env python3
"""Level-3 comparison of the fluence-dependent engine with measured Cu-Ni
observables (#448; docs/validation.md, "Fluence-dependent targets vs
measurement"; data: validation/data/alloys/README.md).

Two separate experiments, never merged:

* Ho et al. (1976), Table 2 row `t2-r4`: 2 keV Ar+ on Cu 38.4 at% / Ni, about
  20 degrees from normal. Observables: the total alloy yield (from the erosion
  rate, Y = N V / J) and the equilibrium Cu surface fraction (from the fitted
  AES ratio R0/Rinf, mapped through a finite-depth AES weighting). The
  conflicting 0.5 keV / 38.4 at% pair (`t1-r3`, `t2-r1`) is excluded.
* Betz (1980), Fig. 13: 2 keV Ar+, normal incidence, total yield of 22
  compositions, a film-integrated measurement compared here with the steady-
  state yield, with the transient bias of the engine's own run reported.

Each case is run for both volume-relaxation conventions (`ideal-mixing`,
`fixed-number-density`) and for two treatments of the implanted argon
(retained as a solid-argon volume, or re-emitted at once). Tolerances are the
ones fixed in advance (docs/validation.md): 20 % relative on the surface
fraction, 30 % relative on the total yield. They are not tuned or relaxed, and
no model parameter is fitted. Pass or fail is reported per observable,
convention and argon treatment.

Writes `validation/experiments/cuni_dynamic_results.json`, which
`validation/update_docs.py` turns into the tables of docs/validation.md.

Usage:
    validation/experiments/cuni_dynamic.py [--threads N]           # full statistics
    validation/experiments/cuni_dynamic.py --reduced               # CI variant, writes nothing
    validation/experiments/cuni_dynamic.py --check                 # data and committed results only
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

DATA = HERE.parent / "data" / "alloys"
RESULTS = HERE / "cuni_dynamic_results.json"
RUNDIR = lindhard_cli.RUNS / "cuni-dynamic"  # gitignored

N_A = 6.02214076e23  # 1/mol (2019 SI)

# Element data as in lindhard/src/elements.rs (`test_cuni_dynamic.py` checks the
# resulting initial number density against the engine's own).
ELEMENTS = {
    "Cu": {"mass": 63.546, "rho": 8.96, "e_s_ev": 3.49},
    "Ni": {"mass": 58.6934, "rho": 8.902, "e_s_ev": 4.44},
}

# Fixed in advance (docs/validation.md): not tuned, not relaxed.
TOL_SURFACE_FRACTION = 0.20
TOL_TOTAL_YIELD = 0.30

# Matched sputter settings of the single-element comparison (run.py,
# SPUTTER_*): ZBL, Lindhard-Scharff, constant free path, cutoffs 2 / 1 eV,
# recoils followed, E_b = 0 and E_d = E_s with E_s the tabulated cohesive energy.
PHYSICS = {
    "potential": "zbl",
    "stopping": "lindhard-scharff",
    "free_path": "constant",
    "primary_cutoff_ev": 2.0,
    "recoil_cutoff_ev": 1.0,
}
# Argon in the target: a noble gas has no tabulated solid density or binding.
# Assumptions, stated rather than sourced: atomic volume of solid argon
# (density 1.77 g/cm3, 20 K: CRC Handbook), no surface binding, E_d placeholder
# equal to Cu's. The two treatments below bracket its effect.
AR_VOLUME_NM3 = 39.948 / (1.77 * N_A) * 1e21
AR_ENERGIES = {"e_d_ev": 3.49, "e_b_ev": 0.0, "e_s_ev": 0.0}
AR_MODES = {
    "ar-retained": None,  # solid-argon volume, no retention limit
    "ar-re-emitted": 1e-6,  # retention limit so small that argon leaves after each step
}
CONVENTIONS = ["ideal-mixing", "fixed-number-density"]

# Numerical settings of the baseline; the sensitivity runs change one at a time.
BASE = {"slab_nm": 0.5, "max_change": 0.1, "ions_per_step": 1000}
SENSITIVITY = [
    ("slab_nm=0.25", {"slab_nm": 0.25}),
    ("slab_nm=1.0", {"slab_nm": 1.0}),
    ("max_change=0.05", {"max_change": 0.05}),
    ("ions=2x", {"ions_scale": 2}),
]
FILM_NM = 40.0  # eroded depth stays well below this at the fluences used
FLUENCE_CM2 = 3.0e16  # about 6 tau of the Ho transient (tau = 50 s, J = 9.4e13)
SEEDS_HO = [1, 2, 3, 4]
IONS_HO = 60000
IONS_BETZ = 30000
BLOCKS = 5  # blocks of the steady window for the standard errors
WINDOW_START = 0.5  # fraction of the run after which the steady window starts

# AES weighting. Seah and Dench (Surf. Interface Anal. 1, 2 (1979)) for elements:
# lambda_ML = 538 / E^2 + 0.41 (a E)^(1/2), E in eV, a the monolayer thickness in
# nm; the attenuation length is taken as this inelastic mean free path at normal
# emission. The analyzer geometry of Ho et al. is not recorded here, so the
# result is also given for 0.5x and 2x the length and for the first slab alone.
AES_LINES_EV = {"Cu": 920.0, "Ni": 716.0}  # Ho et al. 1976, Section 3
AES_SCALES = [0.5, 1.0, 2.0]

HO_ROW = "t2-r4"
HO_EXCLUDED = ["t1-r3", "t2-r1"]

# Diagnostics (#465): static, composition-fixed runs that locate the shortfall.
SPUTTER_DATA = HERE.parent / "data" / "sputtering"
DIAG_ENERGY_EV = 2000.0
DIAG_SEEDS = [1, 2, 3, 4]
DIAG_IONS = 100000  # per seed
DIAG_WINDOWS = {"within_2_percent": 0.02, "within_25_percent": 0.25}  # of the energy, around 2 keV
DIAG_COMPOSITIONS = [0.20, 0.24, 0.268, 0.28, 0.32, 0.34, 0.36, 0.384, 0.40]  # Cu fraction, uniform slab


# --- conversions ------------------------------------------------------------


def atomic_volume_nm3(sym: str) -> float:
    e = ELEMENTS[sym]
    return e["mass"] / (N_A * e["rho"]) * 1e21


def alloy_number_density_cm3(cu_fraction: float) -> float:
    """Total atom density of the ideal-mixing (additive atomic volume) alloy."""
    v = cu_fraction * atomic_volume_nm3("Cu") + (1 - cu_fraction) * atomic_volume_nm3("Ni")
    return 1e21 / v


def alloy_density_g_cm3(cu_fraction: float) -> float:
    m = cu_fraction * ELEMENTS["Cu"]["mass"] + (1 - cu_fraction) * ELEMENTS["Ni"]["mass"]
    return m / N_A * alloy_number_density_cm3(cu_fraction)


def surface_fraction_from_ratio(cu_fraction: float, q: float) -> float:
    """m = c / (q (1 - c) + c) with q = R0/Rinf (alloys/README.md)."""
    return cu_fraction / (q * (1 - cu_fraction) + cu_fraction)


def yield_from_erosion(n_cm3: float, rate_angstrom_min: float, flux_cm2_s: float) -> float:
    """Y = N V / J, V in cm/s; `flux_cm2_s` is per unit sample area."""
    return n_cm3 * rate_angstrom_min * 1e-8 / 60.0 / flux_cm2_s


def seah_dench_nm(sym: str, energy_ev: float) -> float:
    a_nm = atomic_volume_nm3(sym) ** (1 / 3)
    ml = 538.0 / energy_ev**2 + 0.41 * math.sqrt(a_nm * energy_ev)
    return ml * a_nm


def rel_diff(sim: float, ref: float) -> float:
    return (sim - ref) / ref


# --- data -------------------------------------------------------------------


def load_ho() -> dict:
    rows = {}
    for name in ("ho1976_table1.csv", "ho1976_table2.csv"):
        with open(DATA / name, newline="") as f:
            for r in csv.DictReader(f):
                rows[r["row_id"]] = {k: (v if k == "row_id" else float(v)) for k, v in r.items()}
    return rows


def load_betz() -> list[dict]:
    with open(DATA / "betz1980_fig13.csv", newline="") as f:
        return [{k: (v if k == "row_id" else float(v)) for k, v in r.items()} for r in csv.DictReader(f)]


def ho_experiment(row: dict) -> dict:
    """The measured observables of one Ho row and their conversions. Intervals
    use the printed +- amounts; those are not standard deviations."""
    c = row["bulk_cu_at_percent"] / 100
    n = alloy_number_density_cm3(c)
    out = {
        "row_id": row["row_id"],
        "cu_fraction": c,
        "energy_keV": row["energy_keV"],
        "number_density_cm3": n,
        "mass_density_g_cm3": alloy_density_g_cm3(c),
        "erosion_rate_angstrom_min": row["erosion_rate_angstrom_min"],
        "erosion_rate_pm": row["erosion_rate_reported_pm_angstrom_min"],
        "ion_intensity_cm2_s": row["ion_intensity_cm2_s"],
        "q": row["r0_over_rinf"],
        "q_pm": row["r0_over_rinf_reported_pm"],
    }
    cos20 = math.cos(math.radians(20.0))
    conv = {
        # Intensity per unit sample area (no angle correction).
        "flux_per_sample_area": 1.0,
        # Intensity per area normal to the beam: sample-area flux = J cos(20 deg).
        "flux_per_beam_normal_area": 1.0 / cos20,
    }
    y = {}
    for name, factor in conv.items():
        c0 = yield_from_erosion(n, row["erosion_rate_angstrom_min"], row["ion_intensity_cm2_s"]) * factor
        d = yield_from_erosion(n, row["erosion_rate_reported_pm_angstrom_min"], row["ion_intensity_cm2_s"]) * factor
        y[name] = {"yield": c0, "lo": c0 - d, "hi": c0 + d}
    out["yield_conventions"] = y
    q, dq = row["r0_over_rinf"], row["r0_over_rinf_reported_pm"]
    out["surface_cu_fraction"] = {
        "value": surface_fraction_from_ratio(c, q),
        # larger q (stronger depletion) gives the smaller m
        "lo": surface_fraction_from_ratio(c, q + dq),
        "hi": surface_fraction_from_ratio(c, q - dq),
    }
    return out


# --- inputs -----------------------------------------------------------------


def build_input(cu_fraction: float, tilt_deg: float, energy_ev: float, relaxation: str, ar_mode: str,
                ions: int, seed: int, fluence_cm2: float = FLUENCE_CM2, slab_nm: float = BASE["slab_nm"],
                max_change: float = BASE["max_change"], ions_per_step: int = BASE["ions_per_step"]) -> str:
    c = cu_fraction
    lines = [
        "# Generated by validation/experiments/cuni_dynamic.py; do not edit.",
        "[beam]", 'ion = "Ar"', f"energy_ev = {energy_ev!r}", f"tilt_deg = {tilt_deg!r}", "",
        "[materials.CuNi]", f"density_g_cm3 = {alloy_density_g_cm3(c)!r}",
        "elements = [",
        f'  {{ symbol = "Cu", atom_fraction = {c!r} }},',
        f'  {{ symbol = "Ni", atom_fraction = {1 - c!r} }},',
        "]", "",
        "[target]", 'substrate = "CuNi"', "",
        "[[target.layers]]", 'material = "CuNi"', f"thickness_nm = {FILM_NM!r}", "",
        "[physics]",
    ]
    for k, v in PHYSICS.items():
        lines.append(f"{k} = {json.dumps(v)}" if isinstance(v, str) else f"{k} = {v!r}")
    for sym, e in ELEMENTS.items():
        lines += ["", f"[physics.energies.{sym}]", f"e_d_ev = {e['e_s_ev']!r}", "e_b_ev = 0.0", f"e_s_ev = {e['e_s_ev']!r}"]
    lines += ["", "[run]", f"ions = {ions}", f"seed = {seed}", "",
              "[dynamic]", f"fluence_cm2 = {fluence_cm2!r}", f"ions_per_step = {ions_per_step}",
              f"max_change = {max_change!r}", f"slab_nm = {slab_nm!r}", "erosion = true",
              f'relaxation = "{relaxation}"']
    if relaxation == "fixed-number-density":
        lines.append(f"number_density_cm3 = {alloy_number_density_cm3(c)!r}")
    lines += ["", "[dynamic.energies.Ar]"] + [f"{k} = {v!r}" for k, v in AR_ENERGIES.items()]
    if relaxation == "ideal-mixing":
        lines += ["", "[dynamic.atomic_volume_nm3]", f"Ar = {AR_VOLUME_NM3!r}"]
    limit = AR_MODES[ar_mode]
    if limit is not None:
        lines += ["", "[dynamic.retention_limit]", f"Ar = {limit!r}"]
    return "\n".join(lines) + "\n"


def run_lindhard(text: str, workdir: Path, binary: Path, threads: int) -> Path:
    """Run `lindhard` on an input; return the output directory."""
    workdir.mkdir(parents=True, exist_ok=True)
    inp, out = workdir / "input.toml", workdir / "out"
    inp.write_text(text)
    proc = subprocess.run([str(binary), "run", str(inp), "--out", str(out), "--threads", str(threads)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        lindhard_cli.die(f"lindhard failed in {workdir}:\n{proc.stderr}")
    return out


# --- analysis ---------------------------------------------------------------


def _num(v: str) -> float:
    """A CSV cell as a float; an empty cell (not applicable) is NaN."""
    return float(v) if v != "" else math.nan


def read_steps(out: Path) -> list[dict]:
    with open(out / "dynamic_steps.csv", newline="") as f:
        return [{k: _num(v) for k, v in r.items()} for r in csv.DictReader(f)]


def read_profiles(out: Path, steps: set[int]) -> dict[int, list[dict]]:
    """Slab rows of the given steps only (the file is large: filter on the step
    column before parsing a row)."""
    prof: dict[int, list[dict]] = {s: [] for s in steps}
    with open(out / "dynamic_composition.csv") as f:
        header = f.readline().strip().split(",")
        for line in f:
            s = int(line[: line.index(",")])
            if s in prof:
                prof[s].append({k: _num(v) for k, v in zip(header, line.strip().split(","))})
    return prof


def aes_surface_fraction(slabs: list[dict], scale: float | None) -> float:
    """Cu fraction of Cu + Ni seen by an AES signal attenuated as exp(-x / lambda_e)
    in depth x (lambda_e: `seah_dench_nm` of the Cu 920 eV and Ni 716 eV lines times
    `scale`), calibrated on the uniform bulk: the same construction as the authors'
    R0 -> m (alloys/README.md). `scale` None: the first slab alone."""
    if scale is None:
        s0 = min(slabs, key=lambda s: s["slab"])
        cu, ni = s0["atoms_per_cm2_Cu"], s0["atoms_per_cm2_Ni"]
        return cu / (cu + ni)
    eff = {}
    for sym in ("Cu", "Ni"):
        lam = scale * seah_dench_nm(sym, AES_LINES_EV[sym])
        tot = 0.0
        for s in slabs:
            a, x1, x2, t = s[f"atoms_per_cm2_{sym}"], s["front_nm"], s["back_nm"], s["thickness_nm"]
            g = math.exp(-x1 / lam) if t < 1e-9 * lam else (math.exp(-x1 / lam) - math.exp(-x2 / lam)) * lam / t
            tot += a * g
        eff[sym] = tot / lam
    return eff["Cu"] / (eff["Cu"] + eff["Ni"])


def mean_se(xs: list[float]) -> tuple[float, float | None]:
    n = len(xs)
    m = sum(xs) / n
    if n < 2:
        return m, None
    return m, math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1) / n)


def alloy_sputtered(row: dict) -> float:
    return row["cum_sputtered_Cu"] + row["cum_sputtered_Ni"]


def analyze(out: Path, ions: int, aes_scales=AES_SCALES) -> dict:
    """Steady-window yield and surface fractions with block standard errors.

    The yield is the number of Cu and Ni atoms sputtered per ion (sputtered
    argon is re-emitted gas, not erosion of the alloy). The window is the
    second half of the run, in `BLOCKS` blocks of equal ions."""
    steps = read_steps(out)
    done = [r["ions_done"] for r in steps]
    edges = [ions * (WINDOW_START + (1 - WINDOW_START) * i / BLOCKS) for i in range(BLOCKS + 1)]

    def row_at(n_ions: float) -> dict:
        best = steps[0]
        for r in steps:
            if r["ions_done"] <= n_ions + 1e-9:
                best = r
        return best

    bounds = [row_at(e) for e in edges]
    block_y = [(alloy_sputtered(b) - alloy_sputtered(a)) / (b["ions_done"] - a["ions_done"])
               for a, b in zip(bounds, bounds[1:])]
    ymean, yse = mean_se(block_y)
    # Drift: yields over three consecutive windows of the run's quarter points.
    q = [row_at(ions * f) for f in (0.25, 0.5, 0.75, 1.0)]
    windows = [(alloy_sputtered(b) - alloy_sputtered(a)) / (b["ions_done"] - a["ions_done"]) for a, b in zip(q, q[1:])]
    # Transient: intercept b of the line cum(F) = Y_ss F + b over the steady window.
    last = bounds[-1]
    intercept_sputtered_atoms = alloy_sputtered(last) - ymean * last["ions_done"]

    # Surface fractions at the steps of the steady window.
    in_window = [r for r in steps if bounds[0]["ions_done"] < r["ions_done"] <= bounds[-1]["ions_done"]]
    prof = read_profiles(out, {int(r["step"]) for r in in_window})
    modes = {"first_slab": None, **{f"aes_x{s:g}": s for s in aes_scales}}
    per_block = {k: [] for k in modes}
    ar_block = []
    for a, b in zip(bounds, bounds[1:]):
        sel = [r for r in in_window if a["ions_done"] < r["ions_done"] <= b["ions_done"]]
        wsum = sum(r["ions"] for r in sel)
        for k, sc in modes.items():
            per_block[k].append(sum(aes_surface_fraction(prof[int(r["step"])], sc) * r["ions"] for r in sel) / wsum)
        ar = 0.0
        for r in sel:
            s0 = min(prof[int(r["step"])], key=lambda s: s["slab"])
            ar += s0["fraction_Ar"] * r["ions"]
        ar_block.append(ar / wsum)
    surf = {}
    for k, xs in per_block.items():
        mean, se = mean_se(xs)
        surf[k] = {"value": mean, "block_se": se}
    return {
        "ions": ions,
        "yield": {"value": ymean, "block_se": yse, "blocks": block_y},
        "yield_windows_quarters_2_to_4": windows,
        "transient_intercept_sputtered_atoms": intercept_sputtered_atoms,
        "surface_cu_fraction": surf,
        "front_slab_ar_fraction": mean_se(ar_block)[0],
        "total_attempts_rejected": int(sum(r["attempts"] - 1 for r in steps)),
        "steps": len(steps) - 1,
        "recession_nm": steps[-1]["surface_nm"],
    }


def combine(runs: list[dict]) -> dict:
    """Seeds -> one estimate: mean and the standard error of the mean over seeds
    (the block standard error if there is one seed)."""
    def comb(vals, blocks):
        m, se = mean_se(vals)
        if se is None:
            se = blocks[0]
        return {"value": m, "se": se, "per_seed": vals}

    out = {
        "seeds": len(runs),
        "yield": comb([r["yield"]["value"] for r in runs], [runs[0]["yield"]["block_se"]]),
        "surface_cu_fraction": {
            k: comb([r["surface_cu_fraction"][k]["value"] for r in runs],
                    [runs[0]["surface_cu_fraction"][k]["block_se"]]) for k in runs[0]["surface_cu_fraction"]
        },
        "yield_windows_quarters_2_to_4": [mean_se([r["yield_windows_quarters_2_to_4"][i] for r in runs])[0]
                                          for i in range(3)],
        "transient_intercept_sputtered_atoms": mean_se([r["transient_intercept_sputtered_atoms"] for r in runs])[0],
        "front_slab_ar_fraction": mean_se([r["front_slab_ar_fraction"] for r in runs])[0],
        "recession_nm": mean_se([r["recession_nm"] for r in runs])[0],
    }
    return out


def verdict(rel: float, tol: float) -> str:
    return "pass" if abs(rel) <= tol else "fail"


# --- the runs ---------------------------------------------------------------


class Runner:
    """Runs and analyses one input at a time. The analysed result of a run is
    kept in `analysis.json` next to its input and reused when the input text and
    the executable are the same (`reuse`); the large CSV files are deleted after
    the analysis unless `keep_output`."""

    def __init__(self, binary: Path, threads: int, root: Path | None = None, reuse: bool = False,
                 keep_output: bool = False):
        self.binary, self.threads, self.reuse, self.keep_output = binary, threads, reuse, keep_output
        self.root = root or RUNDIR
        self.binary_id = lindhard_cli.binary_sha256(binary)

    def one(self, tag: str, c: float, tilt: float, relaxation: str, ar_mode: str, ions: int, seed: int,
            **over) -> dict:
        scale = over.pop("ions_scale", 1)
        ions *= scale
        text = build_input(c, tilt, 2000.0, relaxation, ar_mode, ions, seed, **over)
        key = hashlib.sha256((self.binary_id + text).encode()).hexdigest()
        cache = self.root / tag / "analysis.json"
        if self.reuse and cache.exists():
            cached = json.loads(cache.read_text())
            if cached["key"] == key:
                return cached["result"]
        out = run_lindhard(text, self.root / tag, self.binary, self.threads)
        res = analyze(out, ions)
        res["input_sha256"] = hashlib.sha256(text.encode()).hexdigest()
        cache.write_text(json.dumps({"key": key, "result": res}))
        if not self.keep_output:
            for name in ("dynamic_composition.csv", "dynamic_steps.csv"):
                (out / name).unlink()
        return res


def ho_case(runner: Runner, row: dict, seeds: list[int], ions: int, sensitivity: bool) -> dict:
    exp = ho_experiment(row)
    c = exp["cu_fraction"]
    res = {"experiment": exp, "runs": {}, "sensitivity": {}}
    for conv in CONVENTIONS:
        for ar in AR_MODES:
            key = f"{conv}/{ar}"
            runs = [runner.one(f"ho_{row['row_id']}_{conv}_{ar}_s{s}", c, 20.0, conv, ar, ions, s) for s in seeds]
            res["runs"][key] = combine(runs)
            if sensitivity:
                for label, over in SENSITIVITY:
                    r = runner.one(f"ho_{row['row_id']}_{conv}_{ar}_{label}", c, 20.0, conv, ar, ions, seeds[0], **over)
                    res["sensitivity"].setdefault(key, {})[label] = r
    res["comparison"] = compare_ho(exp, res["runs"])
    return res


def compare_ho(exp: dict, runs: dict) -> dict:
    out = {}
    m_ref = exp["surface_cu_fraction"]["value"]
    for key, r in runs.items():
        entry = {"yield": {}, "surface_cu_fraction": {}}
        for cname, y in exp["yield_conventions"].items():
            rel = rel_diff(r["yield"]["value"], y["yield"])
            entry["yield"][cname] = {"relative_difference": rel, "tolerance": TOL_TOTAL_YIELD,
                                     "verdict": verdict(rel, TOL_TOTAL_YIELD),
                                     "interval_overlaps_printed_pm": y["lo"] <= r["yield"]["value"] <= y["hi"]}
        for k, v in r["surface_cu_fraction"].items():
            rel = rel_diff(v["value"], m_ref)
            entry["surface_cu_fraction"][k] = {"relative_difference": rel, "tolerance": TOL_SURFACE_FRACTION,
                                               "verdict": verdict(rel, TOL_SURFACE_FRACTION)}
        out[key] = entry
    return out


def betz_case(runner: Runner, points: list[dict], ions: int) -> dict:
    film_nm = 200.0  # about 2000 angstrom (betz1980.json)
    out = {"points": [], "film_nm": film_nm}
    for p in points:
        c = p["bulk_cu_at_percent"] / 100
        entry = {"row_id": p["row_id"], "cu_fraction": c, "measured_yield": p["total_yield_atoms_per_ion"],
                 "digitization_bound": p["yield_digitization_bound_atoms_per_ion"], "runs": {}}
        for conv in CONVENTIONS:
            for ar in AR_MODES:
                r = runner.one(f"betz_{p['row_id']}_{conv}_{ar}", c, 0.0, conv, ar, ions, 1)
                inventory = film_nm * 1e-7 * alloy_number_density_cm3(c)  # atoms/cm2 of the film
                fluence_per_ion = FLUENCE_CM2 / r["ions"]
                bias = r["transient_intercept_sputtered_atoms"] * fluence_per_ion / inventory
                rel = rel_diff(r["yield"]["value"], p["total_yield_atoms_per_ion"])
                entry["runs"][f"{conv}/{ar}"] = {
                    "yield": r["yield"]["value"], "block_se": r["yield"]["block_se"],
                    "relative_difference": rel, "verdict": verdict(rel, TOL_TOTAL_YIELD),
                    "transient_bias_of_film_yield": bias,
                    "yield_windows_quarters_2_to_4": r["yield_windows_quarters_2_to_4"],
                    "input_sha256": r["input_sha256"],
                }
        out["points"].append(entry)
    return out


PARTS = ("ho1976", "betz1980", "diagnostics")


def command_line(argv_extra: str = "") -> str:
    return "python3 validation/experiments/cuni_dynamic.py" + argv_extra


def full_run(threads: int, parts: list[str], reuse: bool) -> dict:
    """The committed results. `parts` ("ho1976", "betz1980") selects what to
    (re)run; a part not selected is kept from the committed file, which then has
    to come from the same binary."""
    binary = lindhard_cli.lindhard_binary()
    runner = Runner(binary, threads, reuse=reuse)
    old = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    sha = lindhard_cli.binary_sha256(binary)
    if old and set(parts) != set(PARTS) and old.get("binary_sha256") != sha:
        lindhard_cli.die("the committed results come from another executable; run every part")
    res = {
        "format": {"name": "lindhard-cuni-dynamic-results", "version": 1},
        "issue": "#448, #465",
        "reproduce": command_line(),
        "binary_sha256": sha,
        "software": lindhard_cli.version(binary),
        "settings": {
            "physics": PHYSICS, "base_numerics": BASE, "film_nm": FILM_NM, "fluence_cm2": FLUENCE_CM2,
            "ions_ho_per_seed": IONS_HO, "seeds_ho": SEEDS_HO, "ions_betz": IONS_BETZ, "seed_betz": 1,
            "blocks": BLOCKS, "steady_window_start": WINDOW_START,
            "argon": {"volume_nm3": AR_VOLUME_NM3, "energies": AR_ENERGIES, "modes": AR_MODES},
            "surface_binding_ev": {s: e["e_s_ev"] for s, e in ELEMENTS.items()},
            "tolerances": {"surface_cu_fraction": TOL_SURFACE_FRACTION, "total_yield": TOL_TOTAL_YIELD},
            "aes": {"lines_ev": AES_LINES_EV, "scales": AES_SCALES,
                    "attenuation_length_nm": {s: seah_dench_nm(s, e) for s, e in AES_LINES_EV.items()}},
        },
        "ho_excluded_rows": HO_EXCLUDED,
    }
    for part in PARTS:
        if part in parts:
            if part == "ho1976":
                res[part] = ho_case(runner, load_ho()[HO_ROW], SEEDS_HO, IONS_HO, sensitivity=True)
            elif part == "betz1980":
                res[part] = betz_case(runner, load_betz(), IONS_BETZ)
            else:
                ho_exp = res["ho1976"]["experiment"] if "ho1976" in res else ho_experiment(load_ho()[HO_ROW])
                ho_run = res.get("ho1976") or old.get("ho1976")
                res[part] = diagnostics(binary, threads, ho_exp, load_betz(), {"ho1976": ho_run} if ho_run else None)
        elif part in old:
            res[part] = old[part]
    return res


def reduced_run(threads: int, root: Path | None = None) -> dict:
    """CI variant: Ho `t2-r4`, both conventions, argon re-emitted, one seed, a
    sixth of the ions. Reads and writes nothing under the committed results."""
    binary = lindhard_cli.lindhard_binary()
    runner = Runner(binary, threads, root=root, keep_output=True)
    exp = ho_experiment(load_ho()[HO_ROW])
    runs = {}
    for conv in CONVENTIONS:
        runs[f"{conv}/ar-re-emitted"] = combine([runner.one(f"reduced_{conv}", exp["cu_fraction"], 20.0, conv,
                                                            "ar-re-emitted", IONS_HO // 6, 1)])
    return {"experiment": exp, "runs": runs, "comparison": compare_ho(exp, runs)}


# --- diagnostics (#465) --------------------------------------------------------


def build_static_input(cu_fraction: float, tilt_deg: float, energy_ev: float, ions: int, seed: int,
                       e_s_ev: dict | None = None) -> str:
    """A static (fixed-composition, no erosion) Ar+ run on a uniform Cu-Ni slab
    or a pure element, with the matched sputter settings of the dynamic runs.
    `e_s_ev` overrides the surface binding energies (mechanism test only)."""
    c = cu_fraction
    es = {sym: (e_s_ev or {}).get(sym, e["e_s_ev"]) for sym, e in ELEMENTS.items()}
    lines = ["# Generated by validation/experiments/cuni_dynamic.py (diagnostics); do not edit.",
             "[beam]", 'ion = "Ar"', f"energy_ev = {energy_ev!r}", f"tilt_deg = {tilt_deg!r}", ""]
    if c in (0.0, 1.0):
        sym = "Cu" if c == 1.0 else "Ni"
        lines += ["[target]", f'substrate = "{sym}"', ""]
        syms = [sym]
    else:
        lines += ["[materials.CuNi]", f"density_g_cm3 = {alloy_density_g_cm3(c)!r}", "elements = [",
                  f'  {{ symbol = "Cu", atom_fraction = {c!r} }},',
                  f'  {{ symbol = "Ni", atom_fraction = {1 - c!r} }},', "]", "",
                  "[target]", 'substrate = "CuNi"', ""]
        syms = ["Cu", "Ni"]
    lines.append("[physics]")
    for k, v in PHYSICS.items():
        lines.append(f"{k} = {json.dumps(v)}" if isinstance(v, str) else f"{k} = {v!r}")
    lines.append("follow_recoils = true")
    for sym in syms:
        lines += ["", f"[physics.energies.{sym}]", f"e_d_ev = {es[sym]!r}", "e_b_ev = 0.0", f"e_s_ev = {es[sym]!r}"]
    lines += ["", "[run]", f"ions = {ions}", f"seed = {seed}"]
    return "\n".join(lines) + "\n"


def static_yields(binary: Path, threads: int, cu_fraction: float, tilt_deg: float, energy_ev: float = DIAG_ENERGY_EV,
                  e_s_ev: dict | None = None, seeds=DIAG_SEEDS, ions: int = DIAG_IONS) -> dict:
    """Per-element sputter yields (atoms/ion) of a static run, mean and standard
    error of the mean over seeds, and the Cu fraction of the sputtered flux."""
    per = {"Cu": [], "Ni": []}
    for seed in seeds:
        text = build_static_input(cu_fraction, tilt_deg, energy_ev, ions, seed, e_s_ev)
        with tempfile.TemporaryDirectory() as d:
            out = run_lindhard(text, Path(d), binary, threads)
            by = {e["symbol"]: e["per_ion"] for e in json.loads((out / "summary.json").read_text())["results"]["sputtering"]["by_element"]}
        for sym in per:
            per[sym].append(by.get(sym, 0.0))
    tot = [a + b for a, b in zip(per["Cu"], per["Ni"])]
    flux_cu = [a / t for a, t in zip(per["Cu"], tot)]
    res = {"cu_fraction": cu_fraction, "tilt_deg": tilt_deg, "energy_ev": energy_ev, "seeds": len(seeds),
           "ions_per_seed": ions}
    for name, xs in (("y_cu", per["Cu"]), ("y_ni", per["Ni"]), ("y_total", tot), ("cu_fraction_of_flux", flux_cu)):
        m, se = mean_se(xs)
        res[name] = {"value": m, "se": se}
    return res


def stored_cu_band(energy_ev: float, window: float) -> dict:
    """The stored measured Ar -> Cu points (normal incidence, flagged points left
    out as in run.py) within `window` (relative) of `energy_ev`."""
    pts = []
    for f in sorted(SPUTTER_DATA.glob("ar_cu_sputter_*.json")):
        d = json.loads(f.read_text())
        if d["incidence_deg"] != 0.0:
            continue
        for p in d["points"]:
            if p.get("flag") == "disagrees_between_compilations":
                continue
            if abs(p["energy_ev"] / energy_ev - 1) <= window:
                pts.append((p["yield"], d["id"]))
    ys = sorted(y for y, _ in pts)
    n = len(ys)
    med = ys[n // 2] if n % 2 else 0.5 * (ys[n // 2 - 1] + ys[n // 2])
    return {"points": n, "sets": sorted({i for _, i in pts}), "min": ys[0], "median": med, "max": ys[-1]}


def count_stored_sets(prefix: str) -> int:
    return len(list(SPUTTER_DATA.glob(f"{prefix}_sputter_*.json")))


def betz_interpolated_yield(cu_percent: float, points: list[dict]) -> dict:
    """Linear interpolation of the Betz Fig. 13 markers to a composition: the
    two bracketing markers and the interpolated value (the markers scatter, so
    this is indicative, with the digitization bound only)."""
    xs = sorted(points, key=lambda p: p["bulk_cu_at_percent"])
    lo = max((p for p in xs if p["bulk_cu_at_percent"] <= cu_percent), key=lambda p: p["bulk_cu_at_percent"])
    hi = min((p for p in xs if p["bulk_cu_at_percent"] >= cu_percent), key=lambda p: p["bulk_cu_at_percent"])
    if hi is lo:
        v = lo["total_yield_atoms_per_ion"]
    else:
        t = (cu_percent - lo["bulk_cu_at_percent"]) / (hi["bulk_cu_at_percent"] - lo["bulk_cu_at_percent"])
        v = lo["total_yield_atoms_per_ion"] * (1 - t) + hi["total_yield_atoms_per_ion"] * t
    return {"value": v, "bracket_at_percent": [lo["bulk_cu_at_percent"], hi["bulk_cu_at_percent"]],
            "bracket_yields": [lo["total_yield_atoms_per_ion"], hi["total_yield_atoms_per_ion"]]}


def required_surface_fraction(grid: list[dict], target_flux_fraction: float) -> float | None:
    """Uniform-slab Cu fraction at which the sputtered flux has the bulk Cu
    fraction (the steady-state condition), by linear interpolation of the grid."""
    g = sorted(grid, key=lambda r: r["cu_fraction"])
    for a, b in zip(g, g[1:]):
        fa, fb = a["cu_fraction_of_flux"]["value"], b["cu_fraction_of_flux"]["value"]
        if (fa - target_flux_fraction) * (fb - target_flux_fraction) <= 0 and fa != fb:
            t = (target_flux_fraction - fa) / (fb - fa)
            return a["cu_fraction"] + t * (b["cu_fraction"] - a["cu_fraction"])
    return None


def diagnostics(binary: Path, threads: int, ho_exp: dict, betz_points: list[dict], committed: dict | None) -> dict:
    """Where the Cu-Ni shortfall comes from (#465): single elements against the
    stored data, then the binary at fixed composition, then the data. No
    parameter is fitted; the only variants are mechanism tests, labelled as such."""
    c_bulk = ho_exp["cu_fraction"]
    out = {"energy_ev": DIAG_ENERGY_EV, "seeds": DIAG_SEEDS, "ions_per_seed": DIAG_IONS,
           "stored_data": {"ar_cu_sets": count_stored_sets("ar_cu"), "ar_ni_sets": count_stored_sets("ar_ni"),
                           "cu_band_normal_incidence": {k: stored_cu_band(DIAG_ENERGY_EV, w) for k, w in DIAG_WINDOWS.items()}}}
    # 1. single elements, static.
    pure = {}
    for tilt in (0.0, 20.0):
        for sym, c in (("Cu", 1.0), ("Ni", 0.0)):
            pure[f"{sym}_{tilt:g}deg"] = static_yields(binary, threads, c, tilt)
    out["pure_elements"] = pure
    band = out["stored_data"]["cu_band_normal_incidence"]["within_2_percent"]
    y_cu0 = pure["Cu_0deg"]["y_total"]["value"]
    out["pure_cu_over_stored"] = {"engine_0deg": y_cu0, "stored_median": band["median"],
                                  "ratio_to_median": y_cu0 / band["median"],
                                  "ratio_to_min": y_cu0 / band["min"], "ratio_to_max": y_cu0 / band["max"]}
    out["pure_ni_over_stored"] = None  # no stored Ar -> Ni set (#465)
    # 2. the binary at fixed composition (uniform slab; the real surface has a depth profile).
    for tilt, tag in ((20.0, "ho_20deg"), (0.0, "betz_0deg")):
        yc, yn = pure[f"Cu_{tilt:g}deg"]["y_total"]["value"], pure[f"Ni_{tilt:g}deg"]["y_total"]["value"]
        grid = [static_yields(binary, threads, c, tilt) for c in DIAG_COMPOSITIONS]
        for r in grid:
            c = r["cu_fraction"]
            r["additive_expectation"] = c * yc + (1 - c) * yn  # unweighted mixture of pure-element yields
            r["total_over_additive"] = r["y_total"]["value"] / r["additive_expectation"]
            r["preferential_factor_cu_over_ni"] = (r["y_cu"]["value"] / r["y_ni"]["value"]) / (c / (1 - c))
        entry = {"grid": grid, "required_surface_fraction_equal_flux": required_surface_fraction(grid, c_bulk)}
        mech = {}
        for c in (0.268, c_bulk):
            mech[f"equal_binding_cu_{c:g}"] = static_yields(binary, threads, c, tilt,
                                                            e_s_ev={"Cu": ELEMENTS["Cu"]["e_s_ev"], "Ni": ELEMENTS["Cu"]["e_s_ev"]})
        entry["mechanism_test_equal_binding_3.49_eV"] = mech
        out[tag] = entry
    # Committed dynamic engine results for the same case, for the comparison.
    if committed:
        runs = committed["ho1976"]["runs"]
        out["dynamic_reference"] = {k: {"yield": v["yield"]["value"],
                                        "first_slab_cu": v["surface_cu_fraction"]["first_slab"]["value"]}
                                    for k, v in runs.items()}
    # 3. the data.
    cos20 = math.cos(math.radians(20.0))
    betz_at_ho = betz_interpolated_yield(ho_exp["cu_fraction"] * 100, betz_points)
    yc = ho_exp["yield_conventions"]
    out["data_consistency"] = {
        "betz_normal_incidence_at_ho_bulk": betz_at_ho,
        "ho_yield_flux_per_sample_area": yc["flux_per_sample_area"]["yield"],
        "ho_yield_flux_per_beam_normal_area": yc["flux_per_beam_normal_area"]["yield"],
        "ho_over_betz_flux_per_sample_area": yc["flux_per_sample_area"]["yield"] / betz_at_ho["value"],
        "ho_over_betz_flux_per_beam_normal_area": yc["flux_per_beam_normal_area"]["yield"] / betz_at_ho["value"],
        "cos20": cos20,
        "engine_tilt_20_over_0_pure_cu": pure["Cu_20deg"]["y_total"]["value"] / y_cu0,
        "ho_alloy_number_density_cm3": ho_exp["number_density_cm3"],
        "ho_alloy_mass_density_g_cm3": ho_exp["mass_density_g_cm3"],
    }
    return out


# --- checks -----------------------------------------------------------------


def check_data(results: bool = True) -> list[str]:
    errs = []
    ho = load_ho()
    if len(ho) != 9:
        errs.append(f"expected 9 Ho rows, found {len(ho)}")
    if HO_ROW not in ho:
        errs.append(f"{HO_ROW} missing")
    betz = load_betz()
    if len(betz) != 22:
        errs.append(f"expected 22 Betz markers, found {len(betz)}")
    if results and RESULTS.exists():
        r = json.loads(RESULTS.read_text())
        if r.get("ho_excluded_rows") != HO_EXCLUDED:
            errs.append("committed results do not record the excluded Ho rows")
        if len(r.get("betz1980", {}).get("points", [])) != 22:
            errs.append("committed results do not hold all 22 Betz points")
        if r["settings"]["tolerances"] != {"surface_cu_fraction": TOL_SURFACE_FRACTION, "total_yield": TOL_TOTAL_YIELD}:
            errs.append("committed results were produced with different tolerances")
        committed_ho = r["ho1976"]["experiment"]
        fresh = ho_experiment(ho[HO_ROW])
        if committed_ho != json.loads(json.dumps(fresh)):
            errs.append("committed Ho conversions differ from the data files; rerun the harness")
    return errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=2, help="worker threads of each run (never changes results)")
    ap.add_argument("--part", choices=list(PARTS), action="append",
                    help="run only this part (repeatable); the other is kept from the committed results")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse the analysed result of a run whose input and executable are unchanged")
    ap.add_argument("--reduced", action="store_true", help="CI variant; prints a summary and writes nothing")
    ap.add_argument("--check", action="store_true", help="check the data files and committed results only")
    args = ap.parse_args()
    if args.check:
        errs = check_data()
        for e in errs:
            print(f"error: {e}", file=sys.stderr)
        return 1 if errs else 0
    errs = check_data(results=False)
    if errs:
        for e in errs:
            print(f"error: {e}", file=sys.stderr)
        return 2
    if args.reduced:
        with tempfile.TemporaryDirectory() as d:
            res = reduced_run(args.threads, Path(d))
        print(json.dumps(res["comparison"], indent=1))
        return 0
    res = full_run(args.threads, args.part or list(PARTS), args.reuse)
    RESULTS.write_text(json.dumps(res, indent=1) + "\n")
    print(f"wrote {RESULTS.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
