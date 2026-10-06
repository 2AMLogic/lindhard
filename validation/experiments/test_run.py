#!/usr/bin/env python3
"""Negative tests of the level-3 provenance enforcement in run.py
(validation/data/README.md, "Enforced"). Standard library only; run by
validation/run.sh, or directly:

    python3 validation/experiments/test_run.py
"""

from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run  # noqa: E402


class SputterDatasetChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        files = sorted(run.SPUTTER_DATA.glob("*.json"))
        if not files:
            raise unittest.SkipTest("no sputter-yield datasets")
        cls.path = files[0]
        cls.good = json.loads(cls.path.read_text())
        cls.provenance = run.PROVENANCE.read_text()

    def errs(self, d, path=None, provenance=None):
        return run.check_sputter_dataset(path or self.path, d, self.provenance if provenance is None else provenance)

    def mutated(self):
        return copy.deepcopy(self.good)

    def test_valid_dataset_passes(self):
        self.assertEqual(self.errs(self.good), [])

    def test_missing_original_reference(self):
        d = self.mutated()
        del d["original_reference"]
        self.assertIn("`original_reference` missing or empty", self.errs(d))

    def test_zero_yield(self):
        d = self.mutated()
        d["points"][0]["yield"] = 0.0
        self.assertIn("`points[0].yield` must be a positive number", self.errs(d))

    def test_zero_uncertainty(self):
        d = self.mutated()
        d["points"][0]["yield_unc_rel"] = 0.0
        self.assertIn("`points[0].yield_unc_rel` must be a positive number", self.errs(d))

    def test_id_differs_from_file_name(self):
        d = self.mutated()
        d["id"] = d["id"] + "_x"
        self.assertTrue(any("differs from the file name" in e for e in self.errs(d)))

    def test_missing_incidence(self):
        d = self.mutated()
        del d["incidence_deg"]
        self.assertIn("`incidence_deg` missing or not a number", self.errs(d))

    def test_empty_points(self):
        d = self.mutated()
        d["points"] = []
        self.assertIn("`points` missing or empty", self.errs(d))

    def test_no_provenance_row(self):
        provenance = self.provenance.replace(self.good["id"], "")
        self.assertIn(f"no row in docs/data-provenance.md names `{self.good['id']}`", self.errs(self.good, provenance=provenance))

    def test_wrong_kind(self):
        d = self.mutated()
        d["kind"] = "range"
        self.assertTrue(any("`kind`" in e for e in self.errs(d)))


class RangeDatasetRegression(unittest.TestCase):
    def test_ranges_still_pass(self):
        provenance = run.PROVENANCE.read_text()
        for f in sorted(run.DATA.glob("*.json")):
            self.assertEqual(run.check_dataset(f, json.loads(f.read_text()), provenance), [], f.name)


class EnergyGroups(unittest.TestCase):
    def test_groups_merge_within_two_percent(self):
        ds = [(Path("a.json"), {"ion": "Ar", "target": "Cu", "incidence_deg": 0.0, "id": "a",
                                "points": [{"energy_ev": 1000.0, "yield": 2.0}, {"energy_ev": 1015.0, "yield": 3.0},
                                           {"energy_ev": 1100.0, "yield": 2.5}]})]
        g = run.energy_groups(ds)
        self.assertEqual([x["n_points"] for x in g], [2, 1])
        self.assertEqual((g[0]["yield_min"], g[0]["yield_median"], g[0]["yield_max"]), (2.0, 2.5, 3.0))


if __name__ == "__main__":
    unittest.main()
