#!/usr/bin/env python3
"""One-change-at-a-time runs for the Au and Cu secondary-electron yield excess
of the Mermin model (issue #242, the residual of #173 under #149).

Every run is the `mermin` configuration of se_yield.py at 800 eV (Mermin
inelastic model, Mott elastic scattering on the DHFS potential with exchange,
step barrier at the mid work function, cutoff 1 eV above the vacuum level,
Kieft-Bosch secondaries) with exactly one thing changed. Nothing here is
gated and nothing changes a default: the table says how far each candidate
moves delta, and in which direction.

Standard library only.

    se_yield_candidates.py --run [--material Au Cu] [--candidates ID ...]
                                   [--histories N] [--seeds S ...]
                                   run and write se_yield_candidates_results.json
                                   (runs of other keys are kept)
    se_yield_candidates.py --markdown   print the table spliced into docs/validation.md
    se_yield_candidates.py --check      the committed results are complete and the
                                        block in docs/validation.md is the --markdown
                                        output (exit 1 if not)

Statistics: SEEDS independent runs of HISTORIES primaries per row. delta is
the mean over the seeds; its error is the standard error of that mean (batch
means, one batch per seed), which includes the correlation of the electrons
of one cascade that the Poisson floor of se_yield.py leaves out. The change
against the baseline is the mean of the per-seed differences (the same seeds
are used for every row), with the standard error of those differences.

`--run` needs the `lindhard` binary (LINDHARD_BIN or a release build) and,
for the acoustic-phonon row, the example
`lindhard-cli/examples/acoustic_phonon_elastic.rs`
(LINDHARD_ACOUSTIC_EXAMPLE or a release build of it).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import se_yield as sy  # noqa: E402

lindhard_cli = sy.lindhard_cli
REPO = sy.REPO
RESULTS = HERE / "se_yield_candidates_results.json"
DOC = REPO / "docs" / "validation.md"
BEGIN = "<!-- validation:level3-se-yield-candidates:begin -->"
END = "<!-- validation:level3-se-yield-candidates:end -->"

MATERIALS = ["Au", "Cu"]
ENERGY_EV = 800.0  # the reference energy of #242
BASE_CONFIG = "mermin"  # of se_yield.CONFIGS
HISTORIES = 2000
SEEDS = list(range(1, 11))

# Fermi energies printed in Table 1 of S. Tanuma, C. J. Powell and D. R. Penn, Surf. Interface Anal. 43, 689
# (2011), doi:10.1002/sia.3522, read in the authors' manuscript at NIMS MDR (doi:10.48505/nims.3238,
# https://mdr.nims.go.jp/pid/f9222d77-030b-44e8-b8d5-5661ab007778, SHA-256 859c35f5...6c63aaa, p. 31 of 60),
# opened 2026-10-09. The manuscript calls E_F "a parameter used in the IMFP calculations" (p. 5) and does not
# say where the values come from, so they are a second published value for a sensitivity run, not a measured
# Fermi energy, and they are not adopted as a band input (docs/data-provenance.md).
FERMI_TPP2011_EV = {"Cu": 8.7, "Au": 9.0}
FERMI_TPP2011_PROVENANCE = (
    "Sensitivity run of #242, not a band default: Fermi energy printed in Table 1 of Tanuma, Powell and Penn, "
    "Surf. Interface Anal. 43, 689 (2011), doi:10.1002/sia.3522 (authors' manuscript, NIMS MDR, "
    "doi:10.48505/nims.3238, p. 31), which does not give its source; work function: the mid value of the "
    "delta(E) runs (se_yield.py, METAL_BAND_PROVENANCE)."
)

# The elastic model of the stand-in row: the Thomas-Fermi Yukawa potential the delta(E) runs used before #149's
# DHFS rerun, with the same exchange correction.
STAND_IN_POTENTIAL = "thomas-fermi-yukawa"

# The materials of this script that have a row in Verduin (2017), Table 3.2 (Al, Si, Au, SiO2), which the
# example reads its acoustic-phonon parameters from. Cu has none.
ACOUSTIC_MATERIALS = ("Au",)

# id -> (candidate of #242, what is changed against the baseline).
CANDIDATES = {
    "baseline": ("-", "nothing: the `mermin` configuration of the δ(E) tables"),
    "elastic-stand-in": ("(a) elastic potential",
                         "Thomas-Fermi Yukawa stand-in potential instead of DHFS (exchange on in both)"),
    "acoustic-phonon": ("(b) quasi-elastic scattering",
                        "below 100 eV the Mott rows are replaced by the acoustic-phonon mean free path and angle "
                        "of Verduin (2017), Table 3.2, mixed linearly up to 200 eV; no energy loss"),
    "phi-low": ("(c) band", "work function at the low end of its cited range"),
    "phi-high": ("(c) band", "work function at the high end of its cited range"),
    "fermi-tpp2011": ("(c) band", "Fermi energy of TPP 2011, Table 1, instead of the free-electron value "
                                  "(one valence electron)"),
    "barrier-off": ("context", "transparent boundary: no barrier at all (not a candidate; the largest effect the "
                               "barrier can have)"),
}
BASELINE = "baseline"


def applies(material: str, cand: str) -> bool:
    return cand != "acoustic-phonon" or material in ACOUSTIC_MATERIALS


def _replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise ValueError(f"expected exactly one {old!r} in the generated input")
    return text.replace(old, new)


def _replace_band(text: str, new_line: str) -> str:
    lines = text.split("\n")
    idx = [i for i, l in enumerate(lines) if l.startswith("band = ")]
    if len(idx) != 1:
        raise ValueError("expected exactly one band line in the generated input")
    lines[idx[0]] = new_line
    return "\n".join(lines)


def make_input(material: str, cand: str, histories: int, seed: int, elf_name: str) -> str:
    """The input of one run: the baseline input with the one change of `cand`. The acoustic-phonon row uses the
    baseline input unchanged; its change is made by the example that runs it."""
    if cand not in CANDIDATES:
        raise ValueError(f"unknown candidate {cand!r}")
    if not applies(material, cand):
        raise ValueError(f"candidate {cand!r} does not apply to {material}")
    text = sy.make_input(material, BASE_CONFIG, ENERGY_EV, histories, seed, elf_name)
    if cand == "elastic-stand-in":
        text = _replace_once(text, f'potential = "{sy.ELASTIC_POTENTIAL}"', f'potential = "{STAND_IN_POTENTIAL}"')
    elif cand in ("phi-low", "phi-high"):
        text = _replace_band(text, sy.band_line(material, cand.split("-")[1]))
    elif cand == "fermi-tpp2011":
        w = sy.BAND[material]["work_function_ev"]["mid"]
        text = _replace_band(text, f'band = {{ kind = "metal", fermi_ev = {FERMI_TPP2011_EV[material]!r}, '
                                   f'work_function_ev = {w!r}, provenance = "{FERMI_TPP2011_PROVENANCE}" }}')
    elif cand == "barrier-off":
        text = _replace_once(text, 'boundary = "step-barrier"', 'boundary = "transparent"')
    return text


# --- running --------------------------------------------------------------------------------------

def acoustic_example() -> Path:
    """`$LINDHARD_ACOUSTIC_EXAMPLE`, or the release build of the example."""
    env = os.environ.get("LINDHARD_ACOUSTIC_EXAMPLE")
    if env:
        return Path(env)
    proc = subprocess.run(
        ["cargo", "build", "--release", "-p", "lindhard-cli", "--example", "acoustic_phonon_elastic",
         "--message-format=json"], cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    if proc.returncode != 0:
        lindhard_cli.die("building the example acoustic_phonon_elastic failed; run cargo by hand to see why")
    for line in proc.stdout.splitlines():
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            continue
        if (msg.get("reason") == "compiler-artifact" and msg.get("executable")
                and msg["target"]["name"] == "acoustic_phonon_elastic"):
            return Path(msg["executable"])
    lindhard_cli.die("could not find the acoustic_phonon_elastic executable in cargo's build output")


def run_one(binary: Path, example, cache: Path, material: str, cand: str, histories: int, seed: int) -> dict:
    elf = sy.OPTICAL / sy.OPTICAL_ELF[material]
    base = lindhard_cli.RUNS if lindhard_cli.RUNS.is_dir() else None
    with tempfile.TemporaryDirectory(prefix="se_yield_candidates_", dir=base) as td:
        td = Path(td)
        shutil.copy(elf, td / elf.name)
        (td / "input.toml").write_text(make_input(material, cand, histories, seed, elf.name))
        out = {"material": material, "candidate": cand, "energy_ev": ENERGY_EV, "seed": seed}
        if cand == "acoustic-phonon":
            proc = subprocess.run([str(example), str(td / "input.toml"), "acoustic"], capture_output=True, text=True)
            if proc.returncode != 0:
                sys.exit(f"error: the example failed ({material} {cand} seed {seed}):\n{proc.stderr}")
            r = json.loads(proc.stdout)
            out.update(histories=r["histories"], slow=r["slow"], delta=r["delta"], eta=r["eta"],
                       elastic_below_100_ev_per_primary=r["elastic_below_100_ev_per_primary"])
            if seed == SEEDS[0]:
                # The control: the example with the input's own tables.
                proc = subprocess.run([str(example), str(td / "input.toml"), "mott"], capture_output=True, text=True)
                if proc.returncode != 0:
                    sys.exit(f"error: the example failed ({material} mott seed {seed}):\n{proc.stderr}")
                c = json.loads(proc.stdout)
                out["control"] = {"delta": c["delta"], "eta": c["eta"],
                                  "elastic_below_100_ev_per_primary": c["elastic_below_100_ev_per_primary"]}
            return out
        proc = subprocess.run([str(binary), "run", str(td / "input.toml"), "--out", str(td / "out"),
                               "--table-cache", str(cache)], capture_output=True, text=True)
        if proc.returncode != 0:
            sys.exit(f"error: lindhard failed ({material} {cand} seed {seed}):\n{proc.stderr}")
        s = json.loads((td / "out" / "electron_summary.json").read_text())
    r = s["results"]
    slow, fast = r["front"]["slow"], r["front"]["fast"]
    out.update(histories=r["histories"], slow=slow["count"], delta=slow["per_primary"], eta=fast["per_primary"])
    return out


def do_run(args) -> int:
    binary = lindhard_cli.lindhard_binary()
    example = acoustic_example() if "acoustic-phonon" in args.candidates else None
    lindhard_cli.RUNS.mkdir(parents=True, exist_ok=True)
    cache = lindhard_cli.RUNS / "se_yield_candidates_table_cache"
    prev = json.loads(RESULTS.read_text()) if RESULTS.exists() else {"runs": []}
    keep = {(r["material"], r["candidate"], r["seed"]): r for r in prev["runs"]}
    for m in args.material:
        for cand in args.candidates:
            if not applies(m, cand):
                continue
            for seed in args.seeds:
                r = run_one(binary, example, cache, m, cand, args.histories, seed)
                keep[(m, cand, seed)] = r
                print(f"{m} {cand} seed {seed}: delta {r['delta']:.4f} eta {r['eta']:.4f}", flush=True)
            RESULTS.write_text(json.dumps(assemble(binary, keep.values()), indent=1, sort_keys=True) + "\n")
    return 0


def assemble(binary: Path, runs) -> dict:
    return {
        "format": {"name": "lindhard-se-yield-candidates-results", "version": 1},
        "software": lindhard_cli.version(binary),
        "settings": {
            "energy_ev": ENERGY_EV, "base_config": BASE_CONFIG, "se_split_ev": 50,
            "elastic_baseline": sy.ELASTIC_ID, "stand_in_potential": STAND_IN_POTENTIAL,
            "fermi_tpp2011_ev": FERMI_TPP2011_EV, "fermi_tpp2011_provenance": FERMI_TPP2011_PROVENANCE,
            "acoustic_materials": list(ACOUSTIC_MATERIALS),
            "candidates": {k: {"candidate": v[0], "change": v[1]} for k, v in CANDIDATES.items()},
        },
        "runs": sorted(runs, key=lambda r: (r["material"], r["candidate"], r["seed"])),
    }


# --- statistics and reporting ---------------------------------------------------------------------

def mean_se(values) -> tuple[float, float]:
    """The mean and its standard error (batch means); the error of a single value is not defined (nan)."""
    values = list(values)
    if len(values) < 2:
        return (values[0] if values else math.nan, math.nan)
    return statistics.fmean(values), statistics.stdev(values) / math.sqrt(len(values))


def by_seed(results, material: str, cand: str, key: str = "delta") -> dict:
    return {r["seed"]: r[key] for r in results["runs"] if r["material"] == material and r["candidate"] == cand}


def summarize(results, material: str, cand: str):
    """delta, eta and the paired change of delta against the baseline for one row, or None if it was not run."""
    d = by_seed(results, material, cand)
    if not d:
        return None
    base = by_seed(results, material, BASELINE)
    common = sorted(set(d) & set(base))
    delta, delta_se = mean_se(d.values())
    eta, eta_se = mean_se(by_seed(results, material, cand, "eta").values())
    out = {"seeds": len(d), "delta": delta, "delta_se": delta_se, "eta": eta, "eta_se": eta_se,
           "histories": sum(r["histories"] for r in results["runs"]
                            if r["material"] == material and r["candidate"] == cand)}
    if cand != BASELINE and common:
        diff, diff_se = mean_se(d[s] - base[s] for s in common)
        b = statistics.fmean(base[s] for s in common)
        out.update(change=diff, change_se=diff_se, change_rel=diff / b)
    return out


def _pm(x, se, nd=3) -> str:
    return f"{x:.{nd}f}" if math.isnan(se) else f"{x:.{nd}f} ± {se:.{nd}f}"


def markdown(datasets, results) -> str:
    if not results:
        return "_No candidate runs are committed yet._"
    out = []
    seeds = sorted({r["seed"] for r in results["runs"]})
    hist = sorted({r["histories"] for r in results["runs"]})
    out.append(f"**One change at a time, {ENERGY_EV:g} eV** (`{results['software']}`; baseline: the `mermin` "
               f"configuration; {len(seeds)} seeds ({seeds[0]} to {seeds[-1]}) of "
               f"{', '.join(str(h) for h in hist)} primaries per row; δ ± the standard error of the mean over "
               "the seeds; the change is the mean of the per-seed differences from the baseline ± its standard "
               "error, and in per cent of the baseline):\n")
    out.append("| Row | Candidate | What is changed | " + " | ".join(f"{m} δ | {m} change | {m} η" for m in MATERIALS) + " |")
    out.append("|---|---|---|" + "---|---|---|" * len(MATERIALS))
    for cand, (group, change) in CANDIDATES.items():
        cells = []
        for m in MATERIALS:
            if not applies(m, cand):
                cells += ["not testable", "-", "-"]
                continue
            s = summarize(results, m, cand)
            if s is None:
                cells += ["not run", "-", "-"]
                continue
            cells.append(_pm(s["delta"], s["delta_se"]))
            cells.append("-" if "change" not in s
                         else f"{s['change']:+.3f} ± {s['change_se']:.3f} ({100 * s['change_rel']:+.0f} %)")
            cells.append(_pm(s["eta"], s["eta_se"]))
        out.append(f"| `{cand}` | {group} | {change} | " + " | ".join(cells) + " |")
    out.append("")
    out.append("**Against the measurements** (measured median δ_max of the reference table above; the baseline "
               "is the simulated δ at 800 eV, not the maximum of its curve):\n")
    out.append("| Material | Measured median δ_max [min, max] | Baseline δ (800 eV) | Excess | Largest reduction by a "
               "candidate row | Excess with it |")
    out.append("|---|---|---|---|---|---|")
    for m in MATERIALS:
        ref = sy.reference_stats(datasets, m)
        base = summarize(results, m, BASELINE)
        if base is None or "dmax_median" not in ref:
            continue
        rows = [(c, summarize(results, m, c)) for c in CANDIDATES
                if c != BASELINE and CANDIDATES[c][0] != "context" and applies(m, c)]
        rows = [(c, s) for c, s in rows if s is not None and "change" in s]
        best = min(rows, key=lambda cs: cs[1]["change"], default=None)
        med = ref["dmax_median"]
        line = (f"| {m} | {med:.3f} [{ref['dmax_min']:.3f}, {ref['dmax_max']:.3f}] | {base['delta']:.3f} | "
                f"{100 * (base['delta'] / med - 1):+.0f} % | ")
        if best is None or best[1]["change"] >= 0:
            line += "none lowers δ | - |"
        else:
            line += (f"`{best[0]}`: {best[1]['change']:+.3f} | "
                     f"{100 * (best[1]['delta'] / med - 1):+.0f} % |")
        out.append(line)
    ctl = [r for r in results["runs"] if r["candidate"] == "acoustic-phonon" and "control" in r]
    if ctl:
        out.append("")
        for r in ctl:
            b = by_seed(results, r["material"], BASELINE).get(r["seed"])
            a = by_seed(results, r["material"], "acoustic-phonon", "elastic_below_100_ev_per_primary")
            out.append(f"Control of the `acoustic-phonon` row ({r['material']}, seed {r['seed']}): the example with "
                       f"the input's own tables gives δ {r['control']['delta']:.4f}; `lindhard run` gives "
                       f"{'not run' if b is None else format(b, '.4f')}. Elastic collisions below 100 eV per "
                       f"primary: {r['control']['elastic_below_100_ev_per_primary']:.0f} with the Mott rows, "
                       f"{statistics.fmean(a.values()):.0f} with the acoustic-phonon rows.")
    return "\n".join(out)


def doc_block() -> str | None:
    text = DOC.read_text()
    if BEGIN not in text or END not in text:
        return None
    return text.split(BEGIN, 1)[1].split(END, 1)[0].strip("\n")


def missing_runs(results) -> list[str]:
    have = {(r["material"], r["candidate"], r["seed"]): r for r in results["runs"]}
    out = []
    for m in MATERIALS:
        for cand in CANDIDATES:
            if not applies(m, cand):
                continue
            for seed in SEEDS:
                r = have.get((m, cand, seed))
                if r is None:
                    out.append(f"{m} {cand} seed {seed}: not run")
                elif r["histories"] != HISTORIES:
                    out.append(f"{m} {cand} seed {seed}: {r['histories']} histories, expected {HISTORIES}")
    return out


def control_errors(results) -> list[str]:
    """The example's control run must reproduce `lindhard run` exactly (same input, same seed)."""
    out = []
    for r in results["runs"]:
        if r["candidate"] != "acoustic-phonon" or "control" not in r:
            continue
        b = by_seed(results, r["material"], BASELINE).get(r["seed"])
        if b is not None and b != r["control"]["delta"]:
            out.append(f"{r['material']} seed {r['seed']}: the example's control delta {r['control']['delta']} "
                       f"differs from the baseline run's {b}")
    return out


def do_check() -> int:
    results = json.loads(RESULTS.read_text()) if RESULTS.exists() else None
    if not results:
        print("error: no committed candidate results", file=sys.stderr)
        return 1
    errs = missing_runs(results) + control_errors(results)
    block = doc_block()
    if block is None:
        errs.append(f"docs/validation.md has no {BEGIN} block")
    elif block != markdown(sy.load_datasets(), results):
        errs.append("the block in docs/validation.md differs from `se_yield_candidates.py --markdown`")
    for e in errs:
        print(f"error: {e}", file=sys.stderr)
    if not errs:
        n = len(results["runs"])
        print(f"{n} candidate runs are complete and docs/validation.md matches them")
    return 1 if errs else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--material", nargs="+", default=MATERIALS, choices=MATERIALS)
    ap.add_argument("--candidates", nargs="+", default=list(CANDIDATES), choices=list(CANDIDATES))
    ap.add_argument("--histories", type=int, default=HISTORIES)
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--results", type=Path, default=None,
                    help="with --run or --markdown: this results file instead of the committed one (screening)")
    args = ap.parse_args()
    if args.results:
        global RESULTS
        RESULTS = args.results
    if args.run:
        return do_run(args)
    if args.markdown:
        print(markdown(sy.load_datasets(), json.loads(RESULTS.read_text()) if RESULTS.exists() else None))
        return 0
    return do_check()


if __name__ == "__main__":
    sys.exit(main())
