#!/usr/bin/env python3
"""Tests of the input generation of backscatter.py (#148): the baseline
variants and the fast-secondary sensitivity inputs. Standard library only;
no simulation is run:

    python3 validation/experiments/test_backscatter.py
"""

from __future__ import annotations

import json
import stat
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import backscatter as bs  # noqa: E402
import se_yield  # noqa: E402


def committed(target: str) -> tuple[str, Path]:
    inp = bs.INPUTS / f"eta_{target.lower()}.toml"
    return inp.read_text(), inp.parent


class SecondaryInputs(unittest.TestCase):
    def test_band_provenance_per_material(self):
        # Al, Cu and Au use the single free-electron-metal string; Si is an
        # insulator with its own provenance (#149). The results JSON records
        # the per-target dict under "band_provenance".
        self.assertIsInstance(se_yield.METAL_BAND_PROVENANCE, str)
        for t in bs.SECONDARY_TARGETS:
            if t == "Si":
                self.assertEqual(se_yield.BAND[t]["kind"], "insulator")
            else:
                self.assertEqual(se_yield.BAND[t]["kind"], "free-electron-metal")
                self.assertIs(se_yield.BAND_PROVENANCE[t], se_yield.METAL_BAND_PROVENANCE)
        committed_meta = json.loads(bs.RESULTS.read_text())["secondary_sensitivity"]
        self.assertEqual(committed_meta["band_provenance"],
                         {t: se_yield.BAND_PROVENANCE[t] for t in bs.SECONDARY_TARGETS})

    def test_targets_and_gaps(self):
        self.assertEqual(bs.SECONDARY_TARGETS, ["Al", "Cu", "Si", "Au"])
        self.assertEqual(set(bs.SECONDARY_GAPS), {"C"})
        for t in bs.SECONDARY_TARGETS:
            self.assertIn(t, se_yield.BAND)
        for e in (1.0, 5.0, 10.0, 30.0):
            self.assertIn(e, bs.SECONDARY_ENERGIES_KEV)

    def test_only_secondaries_boundary_and_band_differ_from_the_baseline(self):
        for t in bs.SECONDARY_TARGETS:
            base, d = committed(t)
            plain = tomllib.loads(bs.variant_input(base, d, 5.0, 1000, True, True))
            sec = tomllib.loads(bs.secondary_input(base, d, t, 5.0, 1000))
            band = sec["electron"]["materials"][t].pop("band")
            expected = tomllib.loads(se_yield.band_line(t, "mid"))["band"]
            self.assertEqual(band, expected)
            self.assertEqual(band["provenance"], se_yield.BAND_PROVENANCE[t].replace('"', "'"))
            if t == "Si":
                self.assertEqual(band["kind"], "insulator")
            else:
                self.assertEqual(band["kind"], "free-electron-metal")
                self.assertEqual(band["valence_electrons_per_atom"],
                                 se_yield.BAND[t]["valence_electrons_per_atom"])
                self.assertEqual(band["work_function_ev"], se_yield.BAND[t]["work_function_ev"]["mid"])
            tr = sec["electron"]["transport"]
            self.assertEqual(tr.pop("secondaries"), "kieft-bosch")
            self.assertEqual(tr.pop("boundary"), "step-barrier")
            ptr = plain["electron"]["transport"]
            self.assertEqual(ptr.pop("secondaries"), "off")
            self.assertEqual(ptr.pop("boundary"), "transparent")
            self.assertEqual(sec, plain)
            self.assertEqual(sec["electron"]["beam"]["energy_ev"], 5000.0)
            self.assertEqual(sec["run"]["histories"], 1000)


class ReuseManifest(unittest.TestCase):
    VERSION = "lindhard 0.0.0 (v1)"

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.elf = self.root / "a.elf"
        self.elf.write_bytes(b"data-1")
        self.text = f'[x]\noptical_elf = {json.dumps(str(self.elf))}\n'
        self.work = self.root / "run"
        self.binary = self.root / "fake"

    def fake(self, code: int = 0, write_summary: bool = True, pre: str = ""):
        summary = '{"software":{"git_describe":"v1"},"results":{"histories":100,' \
                  '"yields":{"backscatter_eta":0.5,"secondary_delta":0.0},"front":{"fast":{"count":1}},' \
                  '"budget":{"relative_imbalance":0.0}},"physics":{"transport":{"layers":[{"elastic_model":"m"}]}},' \
                  '"run":{"table_build_s":0.0}}'
        body = f'mkdir -p "$4"\n' + (f"echo '{summary}' > \"$4/electron_summary.json\"\n" if write_summary else "")
        self.binary.write_text(f"#!/bin/sh\n{pre}{body}exit {code}\n")
        self.binary.chmod(self.binary.stat().st_mode | stat.S_IXUSR)

    def run_it(self, reuse=VERSION):
        return bs.run_one(self.binary, self.text, self.work, None, reuse)

    def test_unchanged_reused(self):
        self.fake()
        self.run_it()
        self.assertTrue(bs.reusable(self.work, self.text, self.VERSION))
        self.binary.unlink()  # a rerun would fail now
        self.assertEqual(self.run_it()["histories"], 100)

    def test_edited_elf_not_reusable(self):
        self.fake()
        self.run_it()
        self.elf.write_bytes(b"data-2")
        self.assertFalse(bs.reusable(self.work, self.text, self.VERSION))

    def test_missing_elf_not_reusable(self):
        self.fake()
        self.run_it()
        self.elf.unlink()
        self.assertFalse(bs.reusable(self.work, self.text, self.VERSION))

    def test_missing_malformed_stale_manifest(self):
        self.fake()
        self.run_it()
        m = self.work / bs.MANIFEST
        good = m.read_text()
        for bad in ("{not json", "[]", '{"version": 99, "optical_elf": {}}', '{"version": 1}',
                    json.dumps({"version": 1, "optical_elf": {str(self.elf): "0" * 64}})):
            m.write_text(bad)
            self.assertFalse(bs.reusable(self.work, self.text, self.VERSION), bad)
        m.unlink()
        self.assertFalse(bs.reusable(self.work, self.text, self.VERSION))
        m.write_text(good)
        self.assertTrue(bs.reusable(self.work, self.text, self.VERSION))

    def test_failed_rerun_leaves_no_manifest(self):
        self.fake()
        self.run_it()
        self.elf.write_bytes(b"data-2")
        self.fake(code=1)
        with self.assertRaises(SystemExit):
            self.run_it()
        self.assertFalse((self.work / bs.MANIFEST).exists())
        self.assertFalse(bs.reusable(self.work, self.text, self.VERSION))

    def test_rerun_after_edit_rewrites_manifest(self):
        self.fake()
        self.run_it()
        self.elf.write_bytes(b"data-2")
        self.run_it()
        self.assertTrue(bs.reusable(self.work, self.text, self.VERSION))

    def test_elf_edited_during_run_writes_no_manifest(self):
        # The CLI "reads" data-1, the file becomes data-2 mid-run, the CLI
        # exits 0: the result must not be recorded as a run of data-2.
        self.fake(pre=f"printf data-2 > {json.dumps(str(self.elf))}\n")
        with self.assertRaises(SystemExit):
            self.run_it()
        self.assertEqual(self.elf.read_bytes(), b"data-2")
        self.assertFalse((self.work / bs.MANIFEST).exists())
        self.assertFalse(bs.reusable(self.work, self.text, self.VERSION))

    def test_manifest_holds_pre_run_prints(self):
        self.fake()
        self.run_it()
        m = json.loads((self.work / bs.MANIFEST).read_text())
        self.assertEqual(m["optical_elf"], bs.elf_fingerprints(self.text))


if __name__ == "__main__":
    unittest.main()
