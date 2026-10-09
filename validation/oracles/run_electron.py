#!/usr/bin/env python3
"""Electron code-to-code comparison (#150; docs/validation.md, "Electron
oracles: Nebula and Geant4 MicroElec"; CONTRIBUTING.md, "Oracles").

Runs `lindhard` (its `[electron]` mode) and every *configured* electron oracle
on the matched problems of `electron_problems.json`: electrons at normal
incidence into bulk Si and bulk Cu at 1, 5 and 20 keV. Compares four scalar
outputs: the backscatter yield eta, the secondary yield delta, the mean depth
at which stopped primaries came to rest, and the radius containing 50 % of the
deposited energy (r50). Writes one summary per (code, problem) to
`validation/oracles/summaries/` (`lindhard-electron-<problem>.json` for our
side alone, `<oracle>-<problem>.json` for each comparison): scalar metrics
with statistical errors, versions, commits, the settings each side ran with,
and the inputs that differ. Raw output and every input file written stay in
the gitignored `validation/oracle-runs/electron/`.

Statistical errors: every code's histories are split into `batches`
independent groups (lindhard: separate runs with seeds seed, seed + 1, ...;
the oracles: primary index modulo `batches`). Each metric is evaluated on the
pooled histories, and its standard error is the standard deviation of the
batch values divided by sqrt(batches) (the batch-means method).

The oracles are third-party programs the user installs **outside this tree**,
unmodified; nothing in CI needs them. Configuration, by environment variable:

    CSTOOL_SRC      git clone of Nebula's cstool (BSD-3-Clause, Tier A),
                    https://github.com/Nebula-simulator/cstool. Required for
                    the lindhard side too: the optical ELF and band parameters
                    of each material are read at run time from its
                    `data/materials/<material>.yaml` and the ELF file it
                    names, so that lindhard and Nebula run on the same
                    dielectric and band inputs. These files cite no source
                    (docs/data-provenance.md), so they are never copied into
                    this tree; the summaries record their commit and SHA-256.
    NEBULA_BIN      Nebula's `nebula_cpu_edep` executable (BSD-3-Clause,
                    Tier A), https://github.com/Nebula-simulator/nebula
    NEBULA_SRC      optional: the clone it was built from (version, commit)
    NEBULA_MATERIALS directory holding the `.mat` files that cstool compiled
                    from CSTOOL_SRC's `data/materials/*.yaml`
    GEANT4_MICROELEC_BIN
                    `lindhard_microelec`, our application in
                    `geant4_microelec/`, built against a Geant4 installation
                    outside this tree (Geant4 MicroElec: oracle and papers
                    only)
    GEANT4_SRC      optional: the Geant4 clone it was built from
    <NAME>_VERSION, <NAME>_COMMIT, <NAME>_BUILD
                    optional overrides (NAME = NEBULA, GEANT4)

An oracle that is not configured is skipped with a message; the script still
exits 0. Without CSTOOL_SRC nothing can run and the script says so.

Usage:
    validation/oracles/run_electron.py                     # all problems, all configured oracles
    validation/oracles/run_electron.py --problem e_5keV_si
    validation/oracles/run_electron.py --lindhard-only
    validation/oracles/run_electron.py --scale 0.1         # a quick run at 10 % of the histories
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import math
import os
import platform
import re
import struct
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import lindhard_cli  # noqa: E402

SUMMARIES = HERE / "summaries"
PROBLEMS = HERE / "electron_problems.json"
RUNS = lindhard_cli.RUNS / "electron"
# Shared by every problem and batch: the key of each entry covers everything
# the table depends on (lindhard `--table-cache`).
TABLE_CACHE = RUNS / "lindhard" / "table-cache"
FORMAT = "lindhard-oracle-electron-summary/1"
LINDHARD_FORMAT = "lindhard-electron-run/1"
METRICS = ("eta", "delta", "primary_depth_nm", "r50_nm")
# Absolute differences for the yields, relative ones for the lengths.
ABSOLUTE = ("eta", "delta")


class NotConfigured(Exception):
    pass


def _git(src: str | None, *args: str) -> str | None:
    if not src:
        return None
    try:
        out = subprocess.run(["git", "-C", src, *args], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return (out.stdout.strip() or None) if out.returncode == 0 else None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Material parameters, read from the user's cstool clone


def _quantity(text: str, unit: str) -> float:
    """`7.27 eV` -> 7.27; refuses any other unit (cstool documentation,
    "Parameter file format": quantities carry their unit)."""
    m = re.fullmatch(r"\s*([-+0-9.eE]+)\s*" + re.escape(unit) + r"\s*(#.*)?", text)
    if not m:
        lindhard_cli.die(f"cannot read {text!r} as a quantity in {unit}")
    return float(m.group(1))


def cstool_material(src: Path, name: str) -> dict:
    """Band structure and optical file of a cstool parameter file, read with a
    minimal parser of the documented layout (Nebula documentation, "cstool
    parameter file format", sections "Band structure" and "Optical")."""
    path = src / "data" / "materials" / name
    if not path.is_file():
        lindhard_cli.die(f"{path} not found (CSTOOL_SRC)")
    section = None
    band: dict[str, str] = {}
    df = None
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith((" ", "\t")):
            section = line.split(":", 1)[0].strip()
            continue
        key, _, value = line.strip().partition(":")
        if section == "band_structure":
            band[key.strip()] = value.strip()
        elif section == "optical" and key.strip() == "df_file":
            df = value.strip()
    if not band or not df:
        lindhard_cli.die(f"{path}: no band_structure or optical.df_file")
    model = band["model"]
    if model == "metal":
        lband = {
            "kind": "metal",
            "fermi_ev": _quantity(band["fermi"], "eV"),
            "work_function_ev": _quantity(band["work_function"], "eV"),
        }
        inner = lband["fermi_ev"] + lband["work_function_ev"]
    elif model in ("semiconductor", "insulator"):
        lband = {
            "kind": "insulator",
            "valence_band_width_ev": _quantity(band["valence"], "eV"),
            "band_gap_ev": _quantity(band["band_gap"], "eV"),
            "affinity_ev": _quantity(band["affinity"], "eV"),
        }
        inner = lband["valence_band_width_ev"] + lband["band_gap_ev"] + lband["affinity_ev"]
    else:
        lindhard_cli.die(f"{path}: unknown band model {model!r}")
    df_path = path.parent / df
    return {
        "yaml": path,
        "yaml_sha256": sha256(path),
        "df": df_path,
        "df_sha256": sha256(df_path),
        "band": lband,
        "inner_potential_ev": inner,
    }


def read_df(path: Path) -> tuple[list[float], list[float]]:
    """A cstool optical data file (Nebula documentation, "cstool optical data
    file format"): a first line of binding energies ending in -1, then
    `energy_eV ELF` pairs, ended by `-1 -1`."""
    lines = [ln.split() for ln in path.read_text().splitlines() if ln.strip()]
    energy, elf = [], []
    for row in lines[1:]:
        e, f = float(row[0]), float(row[1])
        if e < 0:
            break
        energy.append(e)
        elf.append(f)
    return energy, elf


# ---------------------------------------------------------------------------
# Statistics


def r50(hist: list[float], outside: float, rmax: float, rbins: int) -> float | None:
    """Radius containing half the deposited energy, from a linear radial
    histogram over [0, rmax) plus the energy beyond it; linear within the bin
    where the cumulative sum crosses one half."""
    total = sum(hist) + outside
    if total <= 0:
        return None
    half = 0.5 * total
    width = rmax / rbins
    acc = 0.0
    for i, e in enumerate(hist):
        if acc + e >= half and e > 0:
            return (i + (half - acc) / e) * width
        acc += e
    return None  # more than half beyond rmax


def batch_se(values: list[float | None]) -> float | None:
    v = [x for x in values if x is not None]
    if len(v) < 2 or len(v) != len(values):
        return None
    m = sum(v) / len(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1) / len(v))


class Accumulator:
    """Pooled and per-batch tallies of one code on one problem."""

    def __init__(self, batches: int, rmax: float, rbins: int):
        self.b, self.rmax, self.rbins = batches, rmax, rbins
        self.n = [0] * batches
        self.fast = [0] * batches
        self.slow = [0] * batches
        self.stops = [0] * batches
        self.depth_sum = [0.0] * batches
        self.hist = [[0.0] * (rbins + 1) for _ in range(batches)]  # last: beyond rmax
        self.has_depth = True

    def deposit(self, batch: int, r_nm: float, e: float) -> None:
        i = int(r_nm / self.rmax * self.rbins)
        self.hist[batch][i if 0 <= i < self.rbins else self.rbins] += e

    def metrics(self, which: range | list[int]) -> dict:
        n = sum(self.n[i] for i in which)
        stops = sum(self.stops[i] for i in which)
        h = [sum(self.hist[i][k] for i in which) for k in range(self.rbins + 1)]
        return {
            "eta": sum(self.fast[i] for i in which) / n if n else None,
            "delta": sum(self.slow[i] for i in which) / n if n else None,
            "primary_depth_nm": (sum(self.depth_sum[i] for i in which) / stops
                                 if self.has_depth and stops else None),
            "r50_nm": r50(h[:-1], h[-1], self.rmax, self.rbins),
        }

    def result(self) -> dict:
        pooled = self.metrics(range(self.b))
        per = [self.metrics([i]) for i in range(self.b)]
        return {
            **pooled,
            "histories": sum(self.n),
            "primaries_stopped": sum(self.stops) if self.has_depth else None,
            "std_err": {m: (batch_se([p[m] for p in per]) if pooled[m] is not None else None) for m in METRICS},
            "batch_values": {m: [p[m] for p in per] for m in METRICS},
        }


# ---------------------------------------------------------------------------
# lindhard side


def lindhard_input(problem: dict, spec: dict, mat: dict, elf_name: str, seed: int, histories: int) -> str:
    lc = spec["lindhard"]
    band = dict(mat["band"])
    band["provenance"] = (
        f"oracle matching only: cstool data/materials/{problem['cstool_material']} "
        f"(sha256 {mat['yaml_sha256'][:16]}), the parameters Nebula runs with; cited to no source, "
        "never committed (docs/data-provenance.md)"
    )

    def inline(d: dict) -> str:
        return "{ " + ", ".join(f"{k} = {lindhard_cli._toml_value(v)}" for k, v in d.items()) + " }"

    return "\n".join([
        "# Written by validation/oracles/run_electron.py; do not edit.",
        "[electron.beam]",
        f"energy_ev = {problem['energy_ev']!r}",
        "[electron.transport]",
        f"cutoff_ev = {lc['cutoff_ev']!r}",
        f"cutoff_reference = \"{lc['cutoff_reference']}\"",
        f"secondaries = \"{lc['secondaries']}\"",
        f"boundary = \"{lc['boundary']}\"",
        "[electron.elastic]",
        "model = \"mott\"",
        f"potential = \"{lc['elastic_potential']}\"",
        f"exchange = {lindhard_cli._toml_value(lc['exchange'])}",
        "[electron.inelastic]",
        f"model = \"{lc['inelastic_model']}\"",
        "[electron.tables]",
        f"min_energy_ev = {lc['min_energy_ev']!r}",
        f"points_per_decade = {lc['points_per_decade']!r}",
        f"[electron.materials.{problem['element']}]",
        f"optical_elf = \"{elf_name}\"",
        f"band = {inline(band)}",
        "[electron.tally]",
        f"se_bse_split_ev = {spec['split_ev']!r}",
        "[electron.tally.cylindrical]",
        f"r = {{ lo_nm = 0.0, hi_nm = {problem['rmax_nm']!r}, bins = {problem['rbins']} }}",
        # One depth bin deep enough for everything, so that the grid's
        # outside_ev is the energy beyond rmax only.
        "depth = { lo_nm = 0.0, hi_nm = 100000000.0, bins = 1 }",
        "[target]",
        f"substrate = \"{problem['element']}\"",
        "[run]",
        f"histories = {histories}",
        f"seed = {seed}",
        "",
    ])


def write_elf(problem: dict, mat: dict, work: Path) -> str:
    energy, elf = read_df(mat["df"])
    name = f"{problem['element'].lower()}_elf_cstool.toml"
    prov = (
        f"oracle matching only: the ELF file cstool data/materials/{problem['cstool_material']} names "
        f"({mat['df'].name}, sha256 {mat['df_sha256'][:16]}), the optical data Nebula runs with; its "
        "origin is not stated (docs/data-provenance.md), so it is never committed"
    )
    (work / name).write_text(
        f"material = \"{problem['element']}\"\nprovenance = {json.dumps(prov)}\n"
        f"energy_ev = {json.dumps(energy)}\nelf = {json.dumps(elf)}\n"
    )
    return name


REDACTED = "read at run time from the cstool parameter file named in `materials`; not committed"


def redact_settings(settings: dict) -> dict:
    """Remove every value that comes from cstool's material files from the
    run metadata kept in a summary: the band parameters echoed in the input,
    and the table's upper energy (the beam energy plus the inner potential).
    Those numbers cite no source (docs/data-provenance.md), so they never
    enter the tree; the summaries record the files' SHA-256 instead."""
    s = json.loads(json.dumps(settings))
    el = s.get("input", {}).get("electron", {})
    for m in el.get("materials", {}).values():
        band = m.get("band")
        if isinstance(band, dict):
            m["band"] = {k: (v if k in ("kind", "provenance") else REDACTED) for k, v in band.items()}
    if "max_energy_ev" in el.get("tables", {}):
        el["tables"]["max_energy_ev"] = REDACTED
    return s


def run_lindhard(problem: dict, spec: dict, mat: dict, binary: Path, histories: int) -> tuple[dict, dict]:
    b = spec["batches"]
    work = RUNS / "lindhard" / problem["id"]
    work.mkdir(parents=True, exist_ok=True)
    elf_name = write_elf(problem, mat, work)
    acc = Accumulator(b, problem["rmax_nm"], problem["rbins"])
    per_batch = histories // b
    settings = None
    t0 = time.perf_counter()
    for i in range(b):
        inp = work / f"input_{i}.toml"
        inp.write_text(lindhard_input(problem, spec, mat, elf_name, spec["seed"] + i, per_batch))
        out = work / f"out_{i}"
        # The batches differ only in their seed, so the first one builds the
        # cross-section tables and the others read them from the table cache
        # (docs/cli.md, "Cross-section table cache"). The cache is keyed on
        # the physics, the grid and the lindhard executable, so a rebuilt
        # lindhard or a changed problem never reuses a stale table.
        proc = subprocess.run(
            [str(binary), "run", str(inp), "--out", str(out), "--table-cache", str(TABLE_CACHE)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            lindhard_cli.die(f"lindhard failed on {problem['id']} batch {i}:\n{proc.stderr}")
        tables_line = next((l for l in proc.stderr.splitlines() if l.startswith("tables: ")), None)
        # Log the build of batch 0, and any later batch that did not reuse it.
        if tables_line and (i == 0 or ", 0 built" not in tables_line):
            print(f"  lindhard {problem['id']} batch {i}: {tables_line}")
        s = json.loads((out / "electron_summary.json").read_text())
        r = s["results"]
        acc.n[i] = r["histories"]
        acc.fast[i] = r["front"]["fast"]["count"]
        acc.slow[i] = r["front"]["slow"]["count"]
        prim = r["stopping_points"].get("primaries")
        if prim is None:
            lindhard_cli.die("this lindhard has no results.stopping_points.primaries; rebuild it")
        acc.stops[i] = prim["stopped"]
        if prim["depth"]:
            acc.depth_sum[i] = prim["depth"]["mean"] * 1e9 * prim["stopped"]
        rows = (out / "electron_deposition_cylindrical.csv").read_text().splitlines()
        head = rows[0].split(",")
        ir, ie = head.index("ir"), head.index("energy_ev")
        for row in rows[1:]:
            f = row.split(",")
            acc.hist[i][int(f[ir])] += float(f[ie])
        acc.hist[i][-1] += r["deposition"]["cylindrical"]["outside_ev"]
        if settings is None:
            settings = {
                "software": s["software"],
                "models": s["physics"]["models"],
                "transport": {k: v for k, v in s["physics"]["transport"].items() if k not in ("layers",)},
                "input": s["input"],
                "threads": s["run"]["threads"],
            }
    res = acc.result()
    res["wall_s"] = time.perf_counter() - t0
    return res, redact_settings(settings)


# ---------------------------------------------------------------------------
# Oracles


class Oracle:
    name = "?"
    env = "?"
    prefix = "?"
    source_url = "?"
    license_note = ""

    def _var(self, suffix: str) -> str | None:
        return os.environ.get(self.prefix + suffix)

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
        pass

    def version(self, binary: Path) -> str:
        return self._var("_VERSION") or _git(self._var("_SRC"), "describe", "--tags", "--always") or "unknown"

    def commit(self) -> str:
        return self._var("_COMMIT") or _git(self._var("_SRC"), "rev-parse", "HEAD") or "unknown"

    def build(self) -> str:
        return self._var("_BUILD") or "unknown"

    def applies(self, problem: dict) -> str | None:
        """None if the oracle can run the problem, else the reason it cannot."""
        return None


class Nebula(Oracle):
    name = "Nebula"
    env = "NEBULA_BIN"
    prefix = "NEBULA"
    source_url = "https://github.com/Nebula-simulator/nebula"
    license_note = (
        "BSD-3-Clause (Tier A); run unmodified (CPU build, nebula_cpu_edep) with material files compiled by its "
        "cstool (BSD-3-Clause) from cstool's own parameter files; cstool calls ELSEPA (Tier C: run unmodified, "
        "locally, only inside this toolchain; none of its output is committed) for the Mott cross sections. "
        "Only scalar summaries are committed"
    )
    # Geometry (Nebula documentation, "Nebula geometry format"; the layout of
    # its REELS tutorial): the sample-vacuum surface at z = 0 (material 0
    # below, vacuum -123 above, normal +z), a detector (-126) plane at
    # z = +10 nm, a terminator (-127) plane deep below; lateral half width L.
    # No side mirrors: an electron leaving the simulation domain is removed,
    # which with L much larger than any range here never happens.
    HALF_WIDTH_NM = 1.0e5
    DEPTH_NM = 1.0e5
    DETECTOR_Z_NM = 10.0
    START_Z_NM = 5.0

    def check(self, binary):
        if not os.environ.get("NEBULA_MATERIALS"):
            raise NotConfigured("Nebula: set NEBULA_MATERIALS to the directory of the cstool-compiled .mat files; skipped")

    def version(self, binary):
        v = self._var("_VERSION")
        if v:
            return v
        desc = _git(self._var("_SRC"), "describe", "--tags", "--always")
        return desc or "unknown"

    def geometry(self) -> str:
        L, D = self.HALF_WIDTH_NM, self.DEPTH_NM
        rows = []

        def square(m_in: int, m_out: int, z: float) -> None:
            a, b, c, d = (-L, -L, z), (L, -L, z), (L, L, z), (-L, L, z)
            # Counter-clockwise seen from +z: normal +z, towards m_out.
            for tri in ((a, b, c), (a, c, d)):
                rows.append(f"{m_in} {m_out} " + " ".join(f"{v:.9g}" for p in tri for v in p))

        square(0, -123, 0.0)
        square(-126, -126, self.DETECTOR_Z_NM)
        square(-127, -127, -D)
        return "\n".join(rows) + "\n"

    def run(self, problem, spec, binary, histories):
        work = RUNS / "nebula" / problem["id"]
        work.mkdir(parents=True, exist_ok=True)
        tri = work / "geometry.tri"
        tri.write_text(self.geometry())
        # Electron file format (Nebula documentation, "Nebula electron file
        # format"): x, y, z (nm), direction, energy (eV) as float32, then two
        # int32 tags, native endianness; secondaries inherit the tags. The
        # first tag is the primary's index, so every detected electron and
        # every deposit can be traced to its history (and batch).
        pri = work / "primaries.pri"
        rec = struct.Struct("=7f2i")
        with open(pri, "wb") as f:
            for i in range(histories):
                f.write(rec.pack(0.0, 0.0, self.START_Z_NM, 0.0, 0.0, -1.0, problem["energy_ev"], i, 0))
        mat = Path(os.environ["NEBULA_MATERIALS"]) / problem["nebula_material"]
        if not mat.is_file():
            lindhard_cli.die(f"Nebula material file {mat} not found (NEBULA_MATERIALS)")
        det = work / "detected.bin"
        log = work / "stderr.log"
        seed = spec["seed"]
        # Options (Nebula documentation, "Running Nebula"): seed,
        # energy-threshold (default 0: every electron that can still reach
        # the vacuum is followed), and for nebula_cpu_edep detect-filename
        # and deposit-filename (default: standard output).
        cmd = [str(binary), f"--seed={seed}", f"--detect-filename={det}", str(tri), str(pri), str(mat)]
        b = spec["batches"]
        acc = Accumulator(b, problem["rmax_nm"], problem["rbins"])
        acc.has_depth = False
        for i in range(histories):
            acc.n[i % b] += 1
        t0 = time.perf_counter()
        # The deposit stream: one record per scattering event that changed an
        # electron's energy, 5 float32 (position after the event, nm; energy
        # before it, eV; energy lost, eV) and the 2 int32 tags. Nebula's
        # documentation names this output but not its layout; the layout is
        # read from cpu_energydep.cpp (Nebula commit a50a8e8, BSD-3, Tier A),
        # whose deposit callback writes {x, y, z, E_before, E_before -
        # E_after} and the two tags.
        drec = struct.Struct("=5f2i")
        with open(log, "w") as err:
            proc = subprocess.Popen(cmd, cwd=work, stdout=subprocess.PIPE, stderr=err)
            buf = b""
            hypot = math.hypot
            while True:
                chunk = proc.stdout.read(drec.size * 65536)
                if not chunk:
                    break
                buf += chunk
                usable = len(buf) - len(buf) % drec.size
                for x, y, _z, _e, de, px, _py in drec.iter_unpack(buf[:usable]):
                    acc.deposit(px % b, hypot(x, y), de)
                buf = buf[usable:]
            proc.wait()
        wall = time.perf_counter() - t0
        if proc.returncode != 0:
            lindhard_cli.die(f"Nebula failed on {problem['id']} (exit {proc.returncode}); see {log}")
        version_line = next((ln for ln in log.read_text().splitlines() if ln.startswith("This is Nebula version")), "")
        for x, y, z, dx, dy, dz, e, px, _py in rec.iter_unpack(det.read_bytes()):
            if e >= spec["split_ev"]:
                acc.fast[px % b] += 1
            else:
                acc.slow[px % b] += 1
        res = acc.result()
        res["wall_s"] = wall
        res["settings"] = {
            "executable": "nebula_cpu_edep (CPU, all hardware threads; Nebula's documentation: results are "
                          "reproducible only on one thread, so agreement is statistical)",
            "nebula_reports": version_line.replace("This is ", "") or None,
            "physics": "as compiled in (Nebula physics_config.h at the commit above): full Penn inelastic "
                       "with secondaries, Kieft elastic (ELSEPA Mott above 200 eV, acoustic phonon below "
                       "100 eV, interpolated between, with acoustic-phonon and atomic-recoil losses), quantum "
                       "transmission and refraction at the surface, no empirical interface absorption",
            "material_file": problem["nebula_material"],
            "material_sha256": sha256(mat),
            "energy_threshold_ev": 0.0,
            "seed": seed,
            "geometry": f"surface z = 0, detector z = +{self.DETECTOR_Z_NM} nm, terminator z = -{self.DEPTH_NM} nm, "
                        f"half width {self.HALF_WIDTH_NM} nm; primaries from z = +{self.START_Z_NM} nm along -z",
            "deposit": "energy lost per scattering event at the event position (all electrons), from the "
                       "nebula_cpu_edep deposit stream",
        }
        return res


class Geant4MicroElec(Oracle):
    name = "Geant4 MicroElec"
    env = "GEANT4_MICROELEC_BIN"
    prefix = "GEANT4"
    source_url = "https://github.com/Geant4/geant4"
    license_note = (
        "Geant4 Software License; MicroElec run unmodified as an oracle (papers and outputs only), through our "
        "own application validation/oracles/geant4_microelec/ written from the Geant4 user guides"
    )

    def applies(self, problem):
        if not problem.get("geant4_material"):
            return ("not applicable: the Geant4 Physics Reference Manual 11.4 ('The MicroElec extension for "
                    "microelectronics applications') states that the MicroElec models are valid for silicon only")
        return None

    def run(self, problem, spec, binary, histories):
        work = RUNS / "geant4_microelec" / problem["id"]
        work.mkdir(parents=True, exist_ok=True)
        for old in work.glob("*_*.txt"):
            old.unlink()
        b = spec["batches"]
        threads = os.cpu_count() or 1
        cmd = [str(binary), "--energy-ev", repr(problem["energy_ev"]), "--material", problem["geant4_material"],
               "--events", str(histories), "--threads", str(threads), "--seed", str(spec["seed"]),
               "--out", str(work), "--rmax-nm", repr(problem["rmax_nm"]), "--rbins", str(problem["rbins"]),
               "--batches", str(b)]
        t0 = time.perf_counter()
        with open(work / "stdout.log", "w") as f:
            proc = subprocess.run(cmd, cwd=work, stdout=f, stderr=subprocess.STDOUT)
        wall = time.perf_counter() - t0
        if proc.returncode != 0:
            lindhard_cli.die(f"Geant4 MicroElec failed on {problem['id']} (exit {proc.returncode}); "
                             f"see {work / 'stdout.log'}")
        acc = Accumulator(b, problem["rmax_nm"], problem["rbins"])
        for ev in work.glob("events_*.txt"):
            for line in ev.read_text().splitlines():
                if line.startswith("#"):
                    continue
                i, fast, slow, stopped, depth = line.split()
                k = int(i) % b
                acc.n[k] += 1
                acc.fast[k] += int(fast)
                acc.slow[k] += int(slow)
                if stopped == "1":
                    acc.stops[k] += 1
                    acc.depth_sum[k] += float(depth)
        for dep in work.glob("deposit_*.txt"):
            for line in dep.read_text().splitlines():
                if line.startswith("#"):
                    continue
                k, i, e = line.split()
                acc.hist[int(k)][int(i)] += float(e)
        if sum(acc.n) != histories:
            lindhard_cli.die(f"Geant4 MicroElec: {sum(acc.n)} events tallied, {histories} expected")
        res = acc.result()
        res["wall_s"] = wall
        res["settings"] = {
            "application": "validation/oracles/geant4_microelec/lindhard_microelec.cc",
            "physics": "G4EmStandardPhysics_option4 with G4EmParameters::AddMicroElec(\"Target\") on the slab's "
                       "region; default cut 1 nm",
            "material": problem["geant4_material"],
            "geometry": "slab 0 <= z <= 100 um, half width 100 um, in G4_Galactic; primaries from z = -1 nm along +z",
            "threads": threads,
            "seed": spec["seed"],
            "deposit": "G4Step::GetTotalEnergyDeposit at the post-step point, in the slab",
            "escape": "electrons crossing z = 0 out of the slab, at their kinetic energy there (no surface barrier "
                      "in this setup)",
        }
        return res


ORACLES = [Nebula(), Geant4MicroElec()]


def mismatches(oracle: str, problem: dict) -> tuple[list[str], list[str]]:
    """Inputs made identical, and every input or definition that differs."""
    if oracle == "Nebula":
        matched = [
            "optical ELF: the same cstool ELF file on both sides (lindhard reads it at run time)",
            "band parameters (Fermi energy and work function, or valence width, gap and affinity) from the same "
            "cstool parameter file, so the same inner potential",
            "Kieft-Bosch secondary model; step barrier with quantum transmission and refraction",
            "stopping at the vacuum level (Nebula energy-threshold 0; lindhard vacuum-level cutoff 0.01 eV)",
            "SE/BSE split 50 eV on the vacuum energy, E = 50 eV counted as backscattered",
            "normal incidence, primaries entering at the beam axis",
        ]
        mism = [
            "inelastic model: Nebula runs the full Penn algorithm (Nebula physics_config.h, full_penn; tables "
            "by cstool) with inner-shell ionization energies from the ENDF/B photo-atomic data cstool downloads "
            "(cstool/endf/) and cstool's Fermi correction (energy loss below K - E_F, "
            "cstool/dielectric_function/compile.py); lindhard runs the single-pole Penn approximation with its "
            "model Fermi energy 0 and no inner-shell channels, the loss clamp W <= E - E_F being applied by its "
            "transport. lindhard's own full Penn model was not used: one Si table build at 5 keV did not finish "
            "within 60 min on this host, against about 20 s for the single-pole tables (docs/validation.md)",
            "elastic model: Nebula uses Mott cross sections from ELSEPA as cstool calls it (cstool/mott/mott.py "
            "and elsepa_input.py: Dirac-Fock electron density, Fermi nuclear charge distribution, Furness-McCarthy "
            "exchange, LDA correlation-polarization, muffin-tin model for Z other than 1, 7 and 8, full "
            "partial-wave calculation) above 200 eV and acoustic-phonon scattering below 100 eV, interpolated "
            "between (Nebula physics/kieft/elastic.h), with acoustic-phonon and atomic-recoil energy losses "
            "(physics_config.h); lindhard uses Mott cross sections of its Thomas-Fermi Yukawa stand-in potential, "
            "free atom, with Furness-McCarthy exchange, no correlation-polarization and no phonon or recoil losses "
            "(the Salvat DHFS coefficients and the muffin-tin option are documented gaps)",
            "energy deposition: Nebula's deposit stream gives the energy each electron loses at each event, so a "
            "secondary's energy is counted where its parent lost it and again along its own track, and the "
            "remaining energy of an electron dropped below the vacuum level is not in it; lindhard deposits only "
            "what no secondary carries away, plus the remaining energy where an electron stops. r50 is "
            "therefore a loss-weighted radius for Nebula and a deposit-weighted one for lindhard",
            "primary penetration depth: Nebula's outputs do not distinguish primaries from secondaries, so it has "
            "no value for this metric",
            "random-number streams differ; agreement is statistical only",
        ]
    else:
        matched = [
            "SE/BSE split 50 eV, E = 50 eV counted as backscattered",
            "normal incidence, primaries entering at the beam axis",
            "same target element and density basis (G4_Si; lindhard's element table)",
        ]
        mism = [
            "every cross section differs: MicroElec uses its own dielectric-formalism inelastic and partial-wave "
            "elastic models for Si (Valentin et al. 2012); lindhard the single-pole Penn model on cstool's ELF "
            "and its Thomas-Fermi Yukawa Mott stand-in",
            "low-energy end: MicroElec kills electrons below 16.7 eV and deposits their energy locally (Physics "
            "Reference Manual 11.4); lindhard follows electrons to the vacuum level. delta is therefore not "
            "comparable in kind",
            "surface: no surface barrier in the Geant4 setup (electrons escape at their kinetic energy inside); "
            "lindhard applies the step barrier from cstool's band parameters",
            "random-number streams differ; agreement is statistical only",
        ]
    return matched, mism


def compare(ours: dict, theirs: dict) -> dict:
    out = {}
    for m in METRICS:
        a, b = ours.get(m), theirs.get(m)
        if a is None or b is None:
            continue
        sa, sb = ours["std_err"].get(m), theirs["std_err"].get(m)
        if m in ABSOLUTE:
            d = a - b
            se = math.hypot(sa, sb) if sa is not None and sb is not None else None
            out[f"{m}_abs_diff"] = d
        else:
            if b == 0:
                continue
            d = a / b - 1.0
            se = (a / b) * math.hypot(sa / a, sb / b) if sa is not None and sb is not None and a else None
            out[f"{m.removesuffix('_nm')}_rel_diff"] = d
        key = m.removesuffix("_nm")
        if se is not None:
            out[f"{key}_{'abs' if m in ABSOLUTE else 'rel'}_diff_se"] = se
            out[f"{key}_z"] = d / se if se > 0 else None
    return out


def tolerance(spec: dict, problem: dict, oracle: str, comparison: dict) -> dict:
    t = spec["tolerances"]
    if oracle != t["reference_oracle"] or problem["energy_ev"] not in t["energies_ev"]:
        return {"applies": False, "note": t["report_only"]}
    checks = {}
    d = comparison.get("eta_abs_diff")
    checks["eta"] = {"limit_abs": t["eta_abs"], "diff": d, "pass": None if d is None else abs(d) <= t["eta_abs"]}
    d = comparison.get("r50_rel_diff")
    checks["r50_nm"] = {"limit_rel": t["r50_rel"], "diff": d, "pass": None if d is None else abs(d) <= t["r50_rel"]}
    return {"applies": True, "source": t["source"], "checks": checks}


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
    ap.add_argument("--lindhard-only", action="store_true", help="run the lindhard side only")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply every history count (quick runs)")
    args = ap.parse_args()

    spec = json.loads(PROBLEMS.read_text())
    problems = [p for p in spec["problems"] if not args.problem or p["id"] in args.problem]
    if not problems:
        lindhard_cli.die(f"no problem matches {args.problem}")
    src = os.environ.get("CSTOOL_SRC")
    if not src:
        print("CSTOOL_SRC not set: the material parameters of every problem come from cstool's files "
              "(see the header of this script); nothing run")
        return 0
    cstool = Path(src)

    ready = []
    if not args.lindhard_only:
        for o in ORACLES:
            try:
                b = o.binary()
                ready.append((o, b, o.version(b)))
                print(f"{o.name}: {b} ({ready[-1][2]})")
            except NotConfigured as e:
                print(e)

    binary = lindhard_cli.lindhard_binary()
    ours_version = lindhard_cli.version(binary)
    today = _dt.date.today().isoformat()
    host = f"{platform.system()} {platform.machine()}, {os.cpu_count()} logical CPUs"
    if hasattr(os, "getloadavg"):
        host += ", load average {:.1f} at start".format(os.getloadavg()[0])
    materials_note = {
        "cstool_source": "https://github.com/Nebula-simulator/cstool",
        "cstool_commit": _git(src, "rev-parse", "HEAD") or "unknown",
    }
    SUMMARIES.mkdir(exist_ok=True)
    for p in problems:
        n = max(spec["batches"], int(round(p["histories"] * args.scale)))
        n -= n % spec["batches"]
        mat = cstool_material(cstool, p["cstool_material"])
        materials = {**materials_note, "parameter_file": f"data/materials/{p['cstool_material']}",
                     "parameter_file_sha256": mat["yaml_sha256"],
                     "elf_file": f"data/materials/{mat['df'].relative_to(mat['yaml'].parent)}",
                     "elf_file_sha256": mat["df_sha256"]}
        ours, lsettings = run_lindhard(p, spec, mat, binary, n)
        print(f"lindhard {p['id']}: " + ", ".join(
            f"{m} {ours[m]:.4g} ± {ours['std_err'][m]:.2g}" for m in METRICS
            if ours[m] is not None and ours["std_err"][m] is not None))
        lind = {
            "format": LINDHARD_FORMAT,
            "problem": p["id"],
            "lindhard_version": ours_version,
            "date": today,
            "host": host,
            "histories": ours["histories"],
            "batches": spec["batches"],
            "values": {m: ours[m] for m in METRICS},
            "std_err": ours["std_err"],
            "primaries_stopped": ours["primaries_stopped"],
            "wall_s": ours["wall_s"],
            "materials": materials,
            "lindhard_settings": lsettings,
        }
        path = SUMMARIES / f"lindhard-electron-{p['id']}.json"
        path.write_text(json.dumps(_round(lind), indent=2, ensure_ascii=False) + "\n")
        print(f"  wrote {path.relative_to(lindhard_cli.REPO)}")
        for o, b, v in ready:
            why = o.applies(p)
            if why:
                print(f"  {o.name}: {why}")
                continue
            theirs = o.run(p, spec, b, n)
            print(f"  {o.name}: " + ", ".join(
                f"{m} {theirs[m]:.4g} ± {theirs['std_err'][m]:.2g}" for m in METRICS
                if theirs[m] is not None and theirs["std_err"][m] is not None))
            comparison = compare(ours, theirs)
            matched, mism = mismatches(o.name, p)
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
                "histories": {"lindhard": ours["histories"], "oracle": theirs["histories"]},
                "batches": spec["batches"],
                "lindhard": {**{m: ours[m] for m in METRICS}, "std_err": ours["std_err"]},
                "oracle_values": {**{m: theirs[m] for m in METRICS}, "std_err": theirs["std_err"]},
                "comparison": comparison,
                "tolerance": tolerance(spec, p, o.name, comparison),
                "wall_s": {"lindhard": ours["wall_s"], "oracle": theirs["wall_s"]},
                "materials": materials if o.name == "Nebula" else None,
                "oracle_settings": theirs["settings"],
                "matched": matched,
                "mismatches": mism,
            }
            slug = o.name.lower().replace(" ", "_")
            path = SUMMARIES / f"{slug}-{p['id']}.json"
            path.write_text(json.dumps(_round(summary), indent=2, ensure_ascii=False) + "\n")
            print(f"  wrote {path.relative_to(lindhard_cli.REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
