#!/usr/bin/env python3
"""Level-2 code-to-code comparison (docs/validation.md, CONTRIBUTING.md
"Oracles").

Runs `lindhard` and every *configured* oracle on the matched problems in
`problems.json`, then writes one summary per (oracle, problem) to
`validation/oracles/summaries/`. A summary holds scalar summary metrics only
(range moments, yields, throughput, their differences and standard errors),
the oracle's name, version, commit and the settings each side was run with.
The oracle's raw output and the input file written for it stay in the
gitignored `validation/oracle-runs/`.

Oracles are third-party programs that the user installs outside this tree;
they are never vendored, and nothing in CI needs them. Each one is configured
by environment variables:

    RUSTBCA_BIN    RustBCA executable (GPL: Tier B, run unmodified as an
                   oracle only; adapter written from its wiki, never from its
                   source or example inputs)
    RUSTBCA_SRC    optional: the git clone it was built from, for the
                   recorded version and commit (`git describe`, `rev-parse`)
    OPENTRIM_BIN   OpenTRIM `opentrim` CLI (MIT; its SRIM-2013/SRIM-1996
                   stopping tables are Tier C and are never selected)
    OPENTRIM_SRC   optional: the git clone it was built from
    H5DUMP_BIN     optional: `h5dump` (HDF5 tools) to read OpenTRIM's output;
                   default: `h5dump` on PATH
    <NAME>_VERSION, <NAME>_COMMIT, <NAME>_BUILD
                   optional overrides of the recorded version, commit and
                   build description (compiler, profile)

An oracle whose executable is not configured is skipped with a message, and
the script still exits 0.

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
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

SUMMARIES = HERE / "summaries"
FORMAT = "lindhard-oracle-summary/2"

# Atomic numbers and the masses `lindhard` uses by default for the elements of
# the matched problems: the standard atomic weights of lindhard/src/elements.rs
# (IUPAC/CIAAW 2021; row "Standard atomic weights" in docs/data-provenance.md).
# `matched_settings` cross-checks each target mass against lindhard's own
# echoed density and atom density, so a drift fails loudly.
ELEMENTS = {
    "B": (5, 10.81),
    "Si": (14, 28.085),
    "Ar": (18, 39.95),
    "Cu": (29, 63.546),
    "As": (33, 74.921595),
}
AVOGADRO = 6.02214076e23  # exact (SI 2019)


class NotConfigured(Exception):
    pass


def _git(src: str | None, *args: str) -> str | None:
    if not src:
        return None
    try:
        out = subprocess.run(["git", "-C", src, *args], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


# ---------------------------------------------------------------------------
# lindhard side


def lindhard_settings(problem: dict) -> dict:
    """What `lindhard` actually ran with, read back from the input echo of its
    `summary.json` (docs/cli.md, "Output"): the resolved target elements with
    E_d/E_b/E_s, the atom density, cutoffs and thread count. Every oracle input
    is written from these values, not from a second copy of the defaults."""
    s = json.loads((lindhard_cli.RUNS / "lindhard" / problem["id"] / "out" / "summary.json").read_text())
    eng = s["physics"]["engine"]
    layer = s["physics"]["target"][0]
    mat = layer["material"]
    n_cm3 = layer["atom_density_per_cm3"]
    elements = []
    mean_mass = 0.0
    for e in mat["elements"]:
        z, m = ELEMENTS[e["symbol"]]
        mean_mass += e["atom_fraction"] * m
        elements.append({**e, "z": z, "mass_amu": m})
    # rho = n <M> / N_A: check the mass table against lindhard's own numbers.
    rho = n_cm3 * mean_mass / AVOGADRO
    if abs(rho / mat["density_g_cm3"] - 1.0) > 1e-6:
        lindhard_cli.die(
            f"{problem['id']}: ELEMENTS masses disagree with lindhard's echo "
            f"({rho:.6g} vs {mat['density_g_cm3']} g/cm3); update ELEMENTS from lindhard/src/elements.rs"
        )
    beam = problem["beam"]
    z1, m1 = ELEMENTS[beam["ion"]]
    return {
        "beam": {"symbol": beam["ion"], "z": z1, "mass_amu": beam.get("mass_amu", m1),
                 "energy_ev": beam["energy_ev"], "tilt_deg": beam.get("tilt_deg", 0.0)},
        "material": mat["name"],
        "density_g_cm3": mat["density_g_cm3"],
        "atom_density_per_nm3": n_cm3 * 1e-21,
        "elements": elements,
        "primary_cutoff_ev": eng["primary_cutoff_ev"],
        "recoil_cutoff_ev": eng["recoil_cutoff_ev"],
        "follow_recoils": eng["follow_recoils"],
        "primary_surface_binding_ev": eng["primary_surface_binding_ev"],
        "free_path": eng["free_path"],
        # Weak collisions per collision step (absent from echoes written
        # before #64, which had none).
        "weak_collisions": int(eng.get("weak_collisions", 0)),
        "electronic_loss": eng["electronic_loss"],
        "potential": problem["physics"]["potential"],
        "stopping": problem["physics"]["stopping"],
        "threads": s["run"]["threads"],
        "ions_per_s_transport": s["run"]["ions_per_s"],
        "range_std_err": {
            "rp_nm": s["results"]["range"]["depth"]["mean_std_err_nm"] if s["results"]["range"]["depth"] else None,
            "drp_nm": s["results"]["range"]["depth"]["std_dev_std_err_nm"] if s["results"]["range"]["depth"] else None,
        },
    }


# ---------------------------------------------------------------------------
# Oracles


class Oracle:
    """One external program. Subclasses implement `run`."""

    name = "?"
    env = "?"
    source_url = "?"
    license_note = ""

    def _var(self, suffix: str) -> str | None:
        return os.environ.get(self.env.replace("_BIN", suffix))

    def binary(self) -> Path:
        b = os.environ.get(self.env)
        if not b:
            raise NotConfigured(f"{self.name}: not configured (set {self.env}); skipped")
        p = Path(b)
        if not p.exists():
            raise NotConfigured(f"{self.name}: {self.env}={b} does not exist; skipped")
        self.check(p)
        return p

    def check(self, binary: Path) -> None:
        """Extra prerequisites; raise NotConfigured if one is missing."""

    def version(self, binary: Path) -> str:
        """`$<NAME>_VERSION`, else `git describe` of `$<NAME>_SRC`, else the
        first line of `--version`."""
        v = self._var("_VERSION") or _git(self._var("_SRC"), "describe", "--tags", "--always")
        if v:
            return v
        try:
            out = subprocess.run([str(binary), "--version"], capture_output=True, text=True, timeout=30)
            line = (out.stdout or out.stderr).strip().splitlines()
            return line[0] if line else "unknown"
        except (OSError, subprocess.TimeoutExpired):
            return "unknown"

    def commit(self) -> str:
        return self._var("_COMMIT") or _git(self._var("_SRC"), "rev-parse", "HEAD") or "unknown"

    def build(self) -> str:
        return self._var("_BUILD") or "unknown"

    def run(self, problem: dict, workdir: Path, binary: Path, matched: dict) -> dict:
        """Run the oracle unmodified on `problem` in `workdir`. Return a dict
        with the keys of `problem['metrics']` (None where the oracle has no
        comparable quantity), `ions`, `ions_per_s` (end-to-end wall clock),
        `std_err` (per metric), `settings` (the oracle options used, as our
        own summary), `matched` and `mismatches` (lists of strings)."""
        raise NotImplementedError


def _timed(cmd: list[str], cwd: Path, log: Path, env: dict | None = None) -> float:
    t0 = time.perf_counter()
    with open(log, "w") as f:
        proc = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, env=env)
    dt = time.perf_counter() - t0
    if proc.returncode != 0:
        lindhard_cli.die(f"{cmd[0]} failed (exit {proc.returncode}); see {log}")
    return dt


def _moments(xs: list[float]) -> tuple[float | None, float | None, int]:
    n = len(xs)
    if n < 2:
        return None, None, n
    mean = sum(xs) / n
    var = sum((x - mean) ** 2 for x in xs) / (n - 1)
    return mean, math.sqrt(var), n


def _toml(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        r = repr(v)
        # The RustBCA wiki ("Input File Fields", NOTE): floats need a
        # digit after the point, `10.0` not `10.`.
        return r if any(c in r for c in ".eEn") else r + ".0"
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, list):
        return "[" + ", ".join(_toml(x) for x in v) + "]"
    raise TypeError(v)


class RustBca(Oracle):
    name = "RustBCA"
    env = "RUSTBCA_BIN"
    source_url = "https://github.com/lcpp-org/RustBCA"
    license_note = (
        "GPL-3.0 (Tier B): run unmodified as an oracle; adapter written from the RustBCA wiki only; "
        "source, example inputs, tables and fixtures never read or copied"
    )

    def run(self, problem, workdir, binary, matched):
        # Input format: RustBCA wiki, page "Standalone Code: Input File"
        # (https://github.com/lcpp-org/RustBCA/wiki/Standalone-Code:-Input-File,
        # wiki commit ee59b62), sections "Options", "Material Parameters",
        # "Particle Parameters" and "Geometry Input / Mesh0D". The table names
        # `[options]` and `[material_parameters]` are named on that page; the
        # particle and geometry tables follow the same naming of its section
        # headings (`[particle_parameters]`, `[geometry_input]`). The geometry
        # mode is the documented command-line argument `0D` (same page,
        # second paragraph). Every key below is one documented there.
        b = matched["beam"]
        els = matched["elements"]
        threads = os.cpu_count() or 1
        name = "rb_"
        tilt = math.radians(b["tilt_deg"])
        # "Particle Parameters", dir: x is the depth axis and the x component
        # may not be exactly 1.0, so normal incidence is a 1e-5 rad tilt.
        tilt = max(tilt, 1e-5)
        options = {
            "name": name,
            "track_trajectories": False,
            "track_recoils": matched["follow_recoils"],
            "track_recoil_trajectories": False,
            "track_displacements": False,
            "track_energy_losses": False,
            # "weak_collision_order": weak collisions per step, partners in
            # the annuli sqrt(k + [0, 1)) p_max: the same rings as lindhard's
            # `weak_collisions` (Moller and Eckstein, IPP 9/64 (1988), eq. (26)).
            "weak_collision_order": matched["weak_collisions"],
            "suppress_deep_recoils": False,
            "high_energy_free_flight_paths": False,
            "num_threads": threads,
            # "num_chunks": the page advises raising it when memory use is
            # high; results do not depend on it (checked), run time does
            # (about 4x for 1e5-ion cascade runs), so one chunk per 2000 ions.
            "num_chunks": max(1, int(problem["run"]["ions"]) // 2000),
            "seed": int(problem["run"]["seed"]),
            # "electronic_stopping_mode": LOW_ENERGY_NONLOCAL is
            # Lindhard-Scharff, all loss nonlocal (lindhard's `nonlocal`).
            "electronic_stopping_mode": "LOW_ENERGY_NONLOCAL",
            # "mean_free_path_model": LIQUID is the constant free path
            # dx = n^(-1/3), lindhard's `constant` free path.
            "mean_free_path_model": "LIQUID",
        }
        lines = ["# Written by validation/oracles/run.py from documented keys; do not edit.", "", "[options]"]
        lines += [f"{k} = {_toml(v)}" for k, v in options.items()]
        # NxN arrays of enums ("interaction_potential", "scattering_integral",
        # "root_finder"); a single species interaction set is [[...]].
        lines += [
            'interaction_potential = [["ZBL"]]',
            'scattering_integral = [[{"GAUSS_MEHLER" = {n_points = 10}}]]',
            'root_finder = [[{"NEWTON" = {max_iterations = 100, tolerance = 1e-6}}]]',
            "",
            "[material_parameters]",
            'energy_unit = "EV"',
            'mass_unit = "AMU"',
            f"Eb = {_toml([float(e['e_b_ev']) for e in els])}",
            f"Es = {_toml([float(e['e_s_ev']) for e in els])}",
            f"Ec = {_toml([float(matched['recoil_cutoff_ev']) for _ in els])}",
            f"Ed = {_toml([float(e['e_d_ev']) for e in els])}",
            f"Z = {_toml([e['z'] for e in els])}",
            f"m = {_toml([float(e['mass_amu']) for e in els])}",
            'surface_binding_model = "TARGET"',
            'bulk_binding_model = "INDIVIDUAL"',
            "",
            "[particle_parameters]",
            'length_unit = "NM"',
            'energy_unit = "EV"',
            'mass_unit = "AMU"',
            f"N = [{int(problem['run']['ions'])}]",
            f"m = [{_toml(float(b['mass_amu']))}]",
            f"Z = [{b['z']}]",
            f"E = [{_toml(float(b['energy_ev']))}]",
            f"Ec = [{_toml(float(matched['primary_cutoff_ev']))}]",
            f"Es = [{_toml(float(matched['primary_surface_binding_ev']))}]",
            "pos = [[0.0, 0.0, 0.0]]",
            f"dir = [[{_toml(math.cos(tilt))}, {_toml(math.sin(tilt))}, 0.0]]",
            "",
            "[geometry_input]",
            'length_unit = "NM"',
            "electronic_stopping_correction_factor = 1.0",
            f"densities = {_toml([matched['atom_density_per_nm3'] * e['atom_fraction'] for e in els])}",
            "",
        ]
        inp = workdir / "input.toml"
        inp.write_text("\n".join(lines))
        for old in workdir.glob(f"{name}*.output"):
            old.unlink()
        wall = _timed([str(binary), "0D", str(inp)], workdir, workdir / "stdout.log")

        # Output format: RustBCA wiki, page "Standalone Code: Output Files",
        # sections "[name]summary.output", "[name]deposited.output"
        # (`M, Z, x, y, z, collisions`, x the depth) and
        # "[name]reflected.output" (one row per reflected ion).
        n = int(problem["run"]["ions"])
        depths = []
        for line in (workdir / f"{name}deposited.output").read_text().splitlines():
            f = line.split(",")
            if len(f) >= 3 and f[0].strip():
                depths.append(float(f[2]))
        reflected = sum(1 for ln in (workdir / f"{name}reflected.output").read_text().splitlines() if ln.strip())
        sputtered = None
        if matched["follow_recoils"]:
            sputtered = sum(1 for ln in (workdir / f"{name}sputtered.output").read_text().splitlines() if ln.strip())
        rp, drp, ns = _moments(depths)
        bs = reflected / n
        out = {
            "ions": n,
            "rp_nm": rp,
            "drp_nm": drp,
            "backscatter": bs,
            "sputter_yield": None if sputtered is None else sputtered / n,
            "ions_per_s": n / wall,
            "std_err": {
                "rp_nm": drp / math.sqrt(ns) if drp else None,
                "drp_nm": drp / math.sqrt(2 * (ns - 1)) if drp else None,
                "backscatter": math.sqrt(bs * (1 - bs) / n),
                "sputter_yield": None if sputtered is None else math.sqrt(sputtered) / n,
            },
            "settings": {
                "geometry": "0D (Mesh0D: semi-infinite homogeneous target from x = 0)",
                "options": options,
                "interaction_potential": "ZBL",
                "scattering_integral": "GAUSS_MEHLER, n_points = 10",
                "root_finder": "NEWTON, max_iterations = 100, tolerance = 1e-6",
                "material": {
                    "Z": [e["z"] for e in els],
                    "m_amu": [e["mass_amu"] for e in els],
                    "Eb_ev": [e["e_b_ev"] for e in els],
                    "Es_ev": [e["e_s_ev"] for e in els],
                    "Ed_ev": [e["e_d_ev"] for e in els],
                    "Ec_ev": [matched["recoil_cutoff_ev"] for _ in els],
                    "surface_binding_model": "TARGET (planar)",
                    "bulk_binding_model": "INDIVIDUAL",
                },
                "particle": {
                    "Z": b["z"], "m_amu": b["mass_amu"], "E_ev": b["energy_ev"],
                    "Ec_ev": matched["primary_cutoff_ev"], "Es_ev": matched["primary_surface_binding_ev"],
                    "dir": "normal incidence, 1e-5 rad off the x axis (documented requirement)",
                },
                "densities_per_nm3": [matched["atom_density_per_nm3"] * e["atom_fraction"] for e in els],
                "threads": threads,
                "timing": "end-to-end process wall clock, including setup and writing the particle lists",
            },
            "matched": [
                "potential ZBL (universal screening length)",
                "electronic stopping Lindhard-Scharff, all nonlocal (LOW_ENERGY_NONLOCAL)",
                "constant free path n^(-1/3) (LIQUID)",
                f"weak collisions per step: {matched['weak_collisions']} (weak_collision_order)",
                "primary cutoff Ec and E_s of the beam ion",
                "recoil cutoff Ec, E_b, E_s and E_d of each target element",
                "planar surface barrier (surface_binding_model TARGET, planar by default)",
                "recoils followed" if matched["follow_recoils"] else "recoils not followed (track_recoils = false)",
                "number density and masses as echoed by lindhard",
            ],
            "mismatches": [
                "random-number streams differ (different generators); agreement is statistical only",
                "impact-parameter limit: not stated in the RustBCA documentation; lindhard uses N pi p_max^2 l = 1 "
                "(Biersack-Haggmark 1980)",
                "first free flight: lindhard draws the primary's first flight uniformly in [0, l); RustBCA's "
                "choice is not documented",
                "scattering angle: RustBCA Gauss-Mehler quadrature (10 points) per collision; lindhard a "
                "Gauss-Mehler table interpolated in (eps, beta), max error 3.8 mrad",
                "Lindhard-Scharff: the RustBCA documentation cites Lindhard and Scharff (1961) without stating "
                "the formula; constants (e.g. Firsov vs Lindhard screening length in k_L) may differ",
                "normal incidence is 1e-5 rad off-normal in RustBCA (it rejects an exact x direction)",
            ],
        }
        if matched["weak_collisions"]:
            out["mismatches"].append(
                "weak collisions: the RustBCA input page places the partners in the same annuli but does not say "
                "whether they make recoils, which energy each is evaluated at, or whether they carry local "
                "electronic loss; lindhard: no recoils, each at the energy left after the previous one "
                "(ion::bca module docs)"
            )
        if matched["follow_recoils"]:
            out["mismatches"].append(
                "E_d: the RustBCA input page describes Ed as a filter of the displacement list output, its BCA "
                "page as the threshold for removing an atom from its site. Its sputter yield is the same with "
                "E_d = 30 eV and E_d = E_s (compare ar_1keV_cu and ar_1keV_cu_ed_es), so it does not act as a "
                "displacement threshold there; lindhard uses E_d as the displacement criterion "
                "(Biersack-Haggmark), which controls low-energy sputtering (docs/validation.md, "
                "'Sputter yield vs E_d')"
            )
        return out


class OpenTrim(Oracle):
    name = "OpenTRIM"
    env = "OPENTRIM_BIN"
    source_url = "https://github.com/ir2-lab/OpenTRIM"
    license_note = (
        "MIT (Tier A); its SRIM-2013 and SRIM-1996 electronic stopping (Tier C data) is never selected; "
        "DPASS (Tier C) is not used either"
    )

    def h5dump(self) -> str | None:
        return os.environ.get("H5DUMP_BIN") or shutil.which("h5dump")

    def check(self, binary):
        if not self.h5dump():
            raise NotConfigured(f"{self.name}: needs h5dump (HDF5 tools) to read its output; set H5DUMP_BIN; skipped")

    def version(self, binary):
        v = self._var("_VERSION")
        if v:
            return v
        try:
            out = subprocess.run([str(binary), "-v"], capture_output=True, text=True, timeout=30).stdout
        except (OSError, subprocess.TimeoutExpired):
            out = ""
        m = re.search(r"version\s+(\S+)", out)
        ver = m.group(1) if m else "unknown"
        desc = _git(self._var("_SRC"), "describe", "--tags", "--always")
        return f"{ver} ({desc})" if desc else ver

    def _read(self, h5: Path, path: str) -> list[float]:
        out = subprocess.run(
            [self.h5dump(), "-d", path, "-y", "-w", "0", str(h5)], capture_output=True, text=True
        )
        if out.returncode != 0:
            lindhard_cli.die(f"h5dump {path} failed: {out.stderr.strip()}")
        m = re.search(r"DATA \{\s*(.*?)\s*\}", out.stdout, re.S)
        if not m:
            lindhard_cli.die(f"h5dump {path}: no DATA block")
        body = m.group(1).strip()
        if body.startswith("("):  # scalar printed as (0): value
            body = body.split(":", 1)[1]
        return [float(x) for x in re.split(r"[,\s]+", body) if x]

    def run(self, problem, workdir, binary, matched):
        # Input format: OpenTRIM documentation, "The JSON configuration
        # string" (doc/json_config.md and the generated option reference
        # doc/include/options.dox.md, https://ir2-lab.gitlab.io/opentrim/json_config.html):
        # groups Simulation, Transport, IonBeam, Target, Output, Run and
        # UserTally. Output: "Tallies & events" (doc/tallies.md), sections
        # "User tallies" (IonStop / IonExit events, bins x and atom_id,
        # `data` = events per source ion) and the HDF5 layout
        # (doc/h5file.md, `/run_info/total_ion_count`).
        b = matched["beam"]
        els = matched["elements"]
        threads = os.cpu_count() or 1
        dx = 0.05  # nm, stop-depth bin of the user tally
        depth = 400.0 if not matched["follow_recoils"] else 50.0
        width = 100.0
        # "Fixed flight path": l0 = flight_path_const * R_at with
        # R_at = (4 pi N / 3)^(-1/3) (doc/flightpath.md). lindhard's constant
        # flight is N^(-1/3), so flight_path_const = (4 pi / 3)^(1/3); the
        # impact parameter disc p_max = (N pi l0)^(-1/2) is then lindhard's.
        fp_const = (4.0 * math.pi / 3.0) ** (1.0 / 3.0)
        # Transport/min_energy is one cutoff for the beam ion and recoils.
        min_energy = matched["primary_cutoff_ev"]
        composition = []
        for e in els:
            composition.append({
                "element": {"symbol": e["symbol"], "atomic_number": e["z"], "atomic_mass": e["mass_amu"]},
                "X": e["atom_fraction"],
                "Ed": e["e_d_ev"],
                # El has a documented minimum of 0.001 eV; lindhard's E_b = 0.
                "El": max(e["e_b_ev"], 0.001),
                "Es": e["e_s_ev"],
                # doc/damage.md: "The value of E_r is typically set equal to E_d".
                "Er": e["e_d_ev"],
                "Rc": 0.5,
            })
        cfg = {
            "Simulation": {
                "simulation_type": "FullCascade" if matched["follow_recoils"] else "IonsOnly",
                "screening_type": "ZBL",
                # Values: Off | SRIM96 | SRIM13 | DPASS. SRIM96/SRIM13 are SRIM
                # tables and DPASS is Tier C (CONTRIBUTING.md); none is a
                # published formula we can match, so electronic loss is off.
                "electronic_stopping": "Off",
                "electronic_straggling": "Off",
                "defect_recombination": False,
            },
            "Transport": {
                "min_energy": min_energy,
                "flight_path_type": "Constant",
                "flight_path_const": fp_const,
            },
            "IonBeam": {
                "ion": {"symbol": b["symbol"], "atomic_number": b["z"], "atomic_mass": b["mass_amu"]},
                "energy_distribution": {"type": "SingleValue", "center": b["energy_ev"], "fwhm": 1.0},
                "spatial_distribution": {"geometry": "Surface", "type": "SingleValue",
                                         "center": [0.0, width / 2, width / 2], "fwhm": 1.0},
                "angular_distribution": {"type": "SingleValue",
                                         "center": [math.cos(math.radians(b["tilt_deg"])),
                                                    math.sin(math.radians(b["tilt_deg"])), 0.0],
                                         "fwhm": 1.0},
            },
            "Target": {
                "size": [depth, width, width],
                "origin": [0.0, 0.0, 0.0],
                "cell_count": [1, 1, 1],
                "periodic_bc": [0, 1, 1],
                "materials": [{"id": "M", "density": matched["density_g_cm3"], "composition": composition}],
                "regions": [{"id": "R1", "material_id": "M", "origin": [0.0, 0.0, 0.0], "size": [depth, width, width]}],
            },
            "Output": {
                "title": problem["id"],
                "outfilename": "ot",
                "store_exit_events": False,
                "store_pka_events": False,
                "store_damage_events": False,
                "store_ion_track_events": False,
                "store_dedx": False,
            },
            "Run": {"max_no_ions": int(problem["run"]["ions"]), "threads": threads, "seed": int(problem["run"]["seed"])},
            "UserTally": [
                {"id": "stop", "description": "beam ions at rest, by depth", "event": "IonStop",
                 "bins": {"x": [round(i * dx, 6) for i in range(int(round(depth / dx)) + 1)], "atom_id": [0, 1]}},
                {"id": "exit", "description": "ions leaving through the front face, by species", "event": "IonExit",
                 "bins": {"x": [-1.0, 1e-3], "atom_id": [0, 1, 2]}},
            ],
        }
        inp = workdir / "config.json"
        inp.write_text(json.dumps(cfg, indent=1))
        h5 = workdir / "ot.h5"
        if h5.exists():
            h5.unlink()
        wall = _timed([str(binary), "-f", str(inp)], workdir, workdir / "stdout.log")

        n = int(round(self._read(h5, "/run_info/total_ion_count")[0]))
        stop = self._read(h5, "/user_tally/stop/data")  # per ion, per depth bin
        exit_ = self._read(h5, "/user_tally/exit/data")  # [beam, target species 1]
        w = sum(stop)
        rp = drp = None
        if w > 0:
            centres = [(i + 0.5) * dx for i in range(len(stop))]
            rp = sum(c * p for c, p in zip(centres, stop)) / w
            # Sheppard's correction removes the dx^2/12 of the binning.
            drp = math.sqrt(max(sum((c - rp) ** 2 * p for c, p in zip(centres, stop)) / w - dx * dx / 12, 0.0))
        ns = w * n
        bs = exit_[0]
        mism = [
            "electronic stopping: OpenTRIM offers Off, SRIM96, SRIM13 or DPASS. SRIM tables and DPASS (Tier C) "
            "are excluded, so it runs with electronic loss Off, while lindhard runs Lindhard-Scharff (its CLI has "
            "no electronic-loss-off choice). Ranges and backscatter are NOT like-for-like; the differences are "
            "dominated by the missing electronic loss",
            "random-number streams differ (different generators); agreement is statistical only",
            "one energy cutoff (Transport/min_energy) for the beam ion and recoils, set to lindhard's primary "
            "cutoff",
            "E_b = 0 in lindhard; OpenTRIM's El has a minimum of 0.001 eV",
            "replacement energy Er set to E_d (OpenTRIM doc/damage.md); lindhard's replacement rule is its own",
            "finite box (depth x 100 nm x 100 nm, laterally periodic) instead of a semi-infinite target",
            "ZBL scattering from OpenTRIM's interpolation tables; lindhard a Gauss-Mehler table in (eps, beta)",
            "normal incidence along +x; first-flight convention not documented",
        ]
        if matched["follow_recoils"]:
            mism.append(
                "no surface binding: OpenTRIM 1.2 stores Es but applies no surface barrier (its TODO.md lists "
                "'Handle surface effects (sputtering etc.)'), so recoils leaving the front face are not a "
                "sputter yield; the sputter yield is not compared"
            )
        if matched["weak_collisions"]:
            mism.append(
                f"weak collisions: lindhard runs {matched['weak_collisions']} per collision step; the adapter "
                "maps no OpenTRIM option to them"
            )
        return {
            "ions": n,
            "rp_nm": rp,
            "drp_nm": drp,
            "backscatter": bs,
            "sputter_yield": None,
            "ions_per_s": n / wall,
            "std_err": {
                "rp_nm": drp / math.sqrt(ns) if drp else None,
                "drp_nm": drp / math.sqrt(2 * (ns - 1)) if drp else None,
                "backscatter": math.sqrt(bs * (1 - bs) / n),
                "sputter_yield": None,
            },
            "settings": {
                "Simulation": cfg["Simulation"],
                "Transport": cfg["Transport"],
                "target": {"size_nm": [depth, width, width], "periodic_bc": [0, 1, 1],
                           "density_g_cm3": matched["density_g_cm3"], "composition": composition},
                "tallies": f"IonStop by depth (bins {dx} nm, atom_id 0) and IonExit through x = 0 by species",
                "threads": threads,
                "timing": "end-to-end process wall clock, including setup and writing the HDF5 file",
            },
            "matched": [
                "potential ZBL",
                "constant free path N^(-1/3) (flight_path_const = (4 pi/3)^(1/3) R_at) and its p_max disc",
                "E_d of each target element",
                "primary cutoff",
                "recoils followed" if matched["follow_recoils"] else "recoils not followed (IonsOnly)",
                "density and masses as echoed by lindhard",
            ],
            "mismatches": mism,
        }


ORACLES = [RustBca(), OpenTrim()]


# ---------------------------------------------------------------------------
# Comparison


def timing_problem(problem: dict, factor: int) -> dict | None:
    """The same problem at `factor` times the ion count, for timing only."""
    if factor <= 1:
        return None
    big = json.loads(json.dumps(problem))
    big["run"]["ions"] = int(problem["run"]["ions"]) * factor
    return big


def marginal(n1: int, t1: float, n2: int, t2: float) -> float | None:
    """Ions per second with the fixed (setup) cost removed."""
    return (n2 - n1) / (t2 - t1) if t2 > t1 and n2 > n1 else None


def compare(ours: dict, ours_se: dict, theirs: dict, metrics: list[str]) -> dict:
    """Relative differences for moments and yields, an absolute one for the
    backscatter fraction, a ratio for speed; each difference with its
    standard error (independent runs, errors added in quadrature) and the
    difference in units of it (`*_z`)."""
    out = {}
    tse = theirs.get("std_err", {})
    for m in metrics:
        a, b = ours.get(m), theirs.get(m)
        if a is None or b is None:
            continue
        if m == "ions_per_s":
            out["ions_per_s_ratio"] = a / b
            am, bm = ours.get("ions_per_s_marginal"), theirs.get("ions_per_s_marginal")
            if am and bm:
                out["ions_per_s_marginal_ratio"] = am / bm
            continue
        sa, sb = ours_se.get(m), tse.get(m)
        if m == "backscatter":
            out["backscatter_abs_diff"] = a - b
            if sa is not None and sb is not None:
                se = math.hypot(sa, sb)
                out["backscatter_abs_diff_se"] = se
                out["backscatter_z"] = (a - b) / se if se > 0 else None
            continue
        if b == 0:
            continue
        key = m.removesuffix("_nm")
        out[f"{key}_rel_diff"] = a / b - 1.0
        if sa is not None and sb is not None:
            se = (a / b) * math.hypot(sa / a if a else 0.0, sb / b)
            out[f"{key}_rel_diff_se"] = se
            out[f"{key}_z"] = (a / b - 1.0) / se if se > 0 else None
    return out


def _round(x):
    if isinstance(x, float):
        return float(f"{x:.6g}")
    if isinstance(x, dict):
        return {k: _round(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_round(v) for v in x]
    return x


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", action="append", help="problem id (repeatable)")
    ap.add_argument("--ions", type=int, help="override run.ions")
    ap.add_argument("--lindhard-only", action="store_true", help="run lindhard only and print its metrics")
    ap.add_argument(
        "--timing-factor", type=int, default=5,
        help="also time every code at this many times run.ions, for the marginal (startup-free) ions/s; "
             "0 or 1 disables (default 5)",
    )
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
    host = f"{platform.system()} {platform.machine()}, {os.cpu_count()} logical CPUs"
    if hasattr(os, "getloadavg"):
        # Timings are only as good as the machine was quiet: record the load.
        host += ", load average {:.1f} at start".format(os.getloadavg()[0])
    SUMMARIES.mkdir(exist_ok=True)
    for p in problems:
        t0 = time.perf_counter()
        ours = lindhard_cli.run(p, lindhard_cli.RUNS / "lindhard" / p["id"], binary)
        wall = time.perf_counter() - t0
        matched = lindhard_settings(p)
        ours_se = {
            **matched["range_std_err"],
            "backscatter": math.sqrt(ours["backscatter"] * (1 - ours["backscatter"]) / ours["ions"]),
            "sputter_yield": math.sqrt(ours["sputter_yield"] * ours["ions"]) / ours["ions"],
        }
        # Like-for-like speed: end-to-end process wall clock on both sides.
        ours_cmp = {**ours, "ions_per_s": ours["ions"] / wall}
        big = timing_problem(p, args.timing_factor)
        if big:
            t0 = time.perf_counter()
            lindhard_cli.run(big, lindhard_cli.RUNS / "lindhard-timing" / p["id"], binary)
            ours_cmp["ions_per_s_marginal"] = marginal(ours["ions"], wall, big["run"]["ions"], time.perf_counter() - t0)
        print(f"lindhard {p['id']}: " + ", ".join(f"{m} {ours_cmp[m]:.4g}" for m in p["metrics"] if ours_cmp[m] is not None))
        for o, b, v in ready:
            work = lindhard_cli.RUNS / o.name.lower() / p["id"]
            work.mkdir(parents=True, exist_ok=True)
            try:
                theirs = o.run(p, work, b, matched)
                if big:
                    w2 = lindhard_cli.RUNS / f"{o.name.lower()}-timing" / p["id"]
                    w2.mkdir(parents=True, exist_ok=True)
                    t_big = o.run(big, w2, b, matched)
                    theirs["ions_per_s_marginal"] = marginal(
                        theirs["ions"], theirs["ions"] / theirs["ions_per_s"],
                        t_big["ions"], t_big["ions"] / t_big["ions_per_s"],
                    )
            except NotImplementedError as e:
                print(f"  {e}; skipped")
                continue
            print(f"  {o.name}: " + ", ".join(
                f"{m} {theirs[m]:.4g}" for m in p["metrics"] if theirs.get(m) is not None))
            summary = {
                "format": FORMAT,
                "problem": p["id"],
                "oracle": o.name,
                "oracle_version": v,
                "oracle_commit": o.commit(),
                "oracle_source": o.source_url,
                "oracle_build": o.build(),
                "oracle_license": o.license_note,
                "lindhard_version": ours_version,
                "date": today,
                "host": host,
                "ions": {"lindhard": ours["ions"], "oracle": theirs.get("ions")},
                "threads": {"lindhard": matched["threads"], "oracle": theirs["settings"].get("threads")},
                "lindhard": {
                    **{m: ours_cmp[m] for m in p["metrics"]},
                    "std_err": {m: ours_se.get(m) for m in p["metrics"] if m != "ions_per_s"},
                    "ions_per_s_marginal": ours_cmp.get("ions_per_s_marginal"),
                    "ions_per_s_transport_only": matched["ions_per_s_transport"],
                },
                "oracle_values": {
                    **{m: theirs.get(m) for m in p["metrics"]},
                    "ions_per_s_marginal": theirs.get("ions_per_s_marginal"),
                    "std_err": {m: theirs["std_err"].get(m) for m in p["metrics"] if m != "ions_per_s"},
                },
                "comparison": compare(ours_cmp, ours_se, theirs, p["metrics"]),
                "timing": {
                    "ions_per_s": "end-to-end process wall clock at the run's ion count (setup, table "
                                  "builds and output writing included)",
                    "ions_per_s_marginal": None if not big else (
                        f"(N2 - N1) / (t2 - t1) between runs of N1 = {p['run']['ions']} and "
                        f"N2 = {big['run']['ions']} ions: throughput with fixed setup costs removed"),
                    "threads": "every code runs on all logical CPUs of the host",
                },
                "lindhard_settings": {
                    k: matched[k]
                    for k in ("potential", "stopping", "free_path", "weak_collisions", "electronic_loss",
                              "primary_cutoff_ev",
                              "recoil_cutoff_ev", "follow_recoils", "primary_surface_binding_ev",
                              "density_g_cm3", "atom_density_per_nm3")
                } | {"elements": [{k: e[k] for k in ("symbol", "atom_fraction", "e_d_ev", "e_b_ev", "e_s_ev")}
                                  for e in matched["elements"]]},
                "oracle_settings": theirs["settings"],
                "matched": theirs.get("matched", []),
                "mismatches": theirs.get("mismatches", []),
            }
            path = SUMMARIES / f"{o.name.lower()}-{p['id']}.json"
            path.write_text(json.dumps(_round(summary), indent=2, ensure_ascii=False) + "\n")
            print(f"  {o.name}: wrote {path.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
