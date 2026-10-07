#!/usr/bin/env python3
"""Negative tests of the level-3 provenance enforcement in run.py
(validation/data/README.md, "Enforced"), and tests of the sputter
sensitivity scenarios in validation/update_docs.py (README, "Sensitivity
scenarios"). Standard library only; run by
validation/run.sh, or directly:

    python3 validation/experiments/test_run.py
"""

from __future__ import annotations

import copy
import json
import math
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

    def test_targets_are_grouped_separately(self):
        pts = [{"energy_ev": 1000.0, "yield": 2.0}]
        ds = [(Path("a.json"), {"ion": "Ar", "target": "Cu", "incidence_deg": 0.0, "id": "a", "points": pts}),
              (Path("b.json"), {"ion": "Ar", "target": "Au", "incidence_deg": 0.0, "id": "b", "points": pts})]
        g = run.energy_groups(ds)
        self.assertEqual([(x["target"], x["n_points"]) for x in g], [("Au", 1), ("Cu", 1)])

    def test_flagged_point_is_left_out(self):
        ds = [(Path("a.json"), {"ion": "Ar", "target": "Cu", "incidence_deg": 0.0, "id": "a",
                                "points": [{"energy_ev": 1000.0, "yield": 2.0},
                                           {"energy_ev": 1005.0, "yield": 9.0,
                                            "flag": "disagrees_between_compilations"}]})]
        g = run.energy_groups(ds)
        self.assertEqual((len(g), g[0]["n_points"], g[0]["yield_max"]), (1, 1, 2.0))


class SputterProblem(unittest.TestCase):
    def test_matched_problem_sets_e_d_equal_to_e_s(self):
        p = run.sputter_problem("Ar", "Si", 500.0, 0, 20000, 4.63)
        self.assertEqual(p["physics"]["energies"], {"Si": {"e_d_ev": 4.63, "e_b_ev": 0.0, "e_s_ev": 4.63}})
        self.assertNotIn("weak_collisions", p["physics"])
        self.assertEqual(run.sputter_problem("Ar", "Si", 500.0, 3, 20000, 4.63)["physics"]["weak_collisions"], 3)

    def test_probe_leaves_e_s_to_lindhard(self):
        # Si has no default E_d, so the probe sets one; it must not set E_s, which is what it reads back.
        p = run.sputter_problem("Ar", "Si", 1000.0, 0, 10, None)
        self.assertEqual(p["physics"]["energies"], {"Si": {"e_d_ev": run.PROBE_E_D_EV}})
        self.assertTrue(p["id"].endswith("_probe"))


class AllSputterDatasets(unittest.TestCase):
    def test_every_dataset_passes_and_names_its_target(self):
        provenance = run.PROVENANCE.read_text()
        for f in sorted(run.SPUTTER_DATA.glob("*.json")):
            d = json.loads(f.read_text())
            self.assertEqual(run.check_sputter_dataset(f, d, provenance), [], f.name)
            self.assertTrue(f.stem.startswith(f"{d['ion'].lower()}_{d['target'].lower()}_sputter_"), f.name)
            self.assertEqual(d["incidence_deg"], 0.0, f.name)


sys.path.insert(0, str(HERE.parent))
import update_docs  # noqa: E402


def ds(id_, target, points, ion="Ar"):
    """A minimal sputter dataset for the scenario tests."""
    return (Path(f"{id_}.json"), {"id": id_, "ion": ion, "target": target, "incidence_deg": 0.0,
                                  "points": [{"energy_ev": e, "yield": y} for e, y in points]})


class CodeYield(unittest.TestCase):
    """update_docs.code_yield: committed code yields at exact or bracketed energies (#78)."""

    curve = {100.0: 1.0, 1000.0: 100.0}

    def test_exact_energy(self):
        self.assertEqual(update_docs.code_yield(self.curve, 1000.0), (100.0, "exact"))

    def test_interior_log_log(self):
        # log Y = 2 log E - 4 between the two points: Y(200) = 200**2 / 1e4 = 4.
        y, how = update_docs.code_yield(self.curve, 200.0)
        self.assertEqual(how, "interpolated")
        self.assertAlmostEqual(y, 4.0, places=12)

    def test_no_extrapolation_below_or_above(self):
        self.assertEqual(update_docs.code_yield(self.curve, 99.0), (None, "outside the committed run energies"))
        self.assertEqual(update_docs.code_yield(self.curve, 1001.0), (None, "outside the committed run energies"))

    def test_nonpositive_yields(self):
        self.assertEqual(update_docs.code_yield({100.0: 0.0, 1000.0: 1.0}, 100.0), (None, "nonpositive yield"))
        self.assertEqual(update_docs.code_yield({100.0: 0.0, 1000.0: 1.0}, 500.0), (None, "nonpositive yield"))
        self.assertEqual(update_docs.code_yield({100.0: 1.0, 1000.0: -1.0}, 500.0), (None, "nonpositive yield"))


class Scenarios(unittest.TestCase):
    """update_docs.scenario_comparison: filter, then the unchanged energy_groups rule (#78)."""

    def setUp(self):
        # a: 1000 eV (2.0), 1015 eV (4.0); b: 1010 eV (3.0), 2000 eV (5.0).
        # Together: one group {1000, 1010, 1015} (median 3.0) and {2000} (5.0).
        self.data = [ds("a", "Cu", [(1000.0, 2.0), (1015.0, 4.0)]), ds("b", "Cu", [(1010.0, 3.0), (2000.0, 5.0)]),
                     ds("c", "Au", [(1000.0, 9.0)])]
        e3 = float(f"{(1000.0 * 1010.0 * 1015.0) ** (1 / 3):.4g}")
        # Code yields at the baseline group energies (and one further run energy above).
        self.curves = {"lindhard": {e3: 1.5, 2000.0: 2.5, 4000.0: 4.0}, "rustbca": None}
        self.e3 = e3

    def test_baseline_groups(self):
        r = update_docs.scenario_comparison(self.data, "Ar", "Cu", (), self.curves)
        self.assertEqual([(g["energy_ev"], g["n_points"], g["yield_median"]) for g in r["groups"]],
                         [(self.e3, 3, 3.0), (2000.0, 1, 5.0)])
        c = r["codes"]["lindhard"]
        self.assertEqual((c["ratios"], c["missing"], c["interpolated"]), ([0.5, 0.5], [], 0))

    def test_excluding_sole_contributor_removes_group(self):
        r = update_docs.scenario_comparison(self.data, "Ar", "Cu", ("b",), self.curves)
        # Left: {1000, 1015} -> median 3.0, energy sqrt(1000 * 1015); the 2000 eV group is gone.
        e2 = float(f"{(1000.0 * 1015.0) ** 0.5:.4g}")
        self.assertEqual([(g["energy_ev"], g["n_points"], g["yield_median"]) for g in r["groups"]], [(e2, 2, 3.0)])
        self.assertEqual((r["n_groups"], r["n_sets"], r["n_points"]), (1, 1, 2))

    def test_shifted_energy_is_interpolated(self):
        r = update_docs.scenario_comparison(self.data, "Ar", "Cu", ("a",), self.curves)
        # Left: {1010} (3.0) and {2000} (5.0); 1010 eV lies between the runs at e3 and 2000 eV.
        self.assertEqual([(g["energy_ev"], g["yield_median"]) for g in r["groups"]], [(1010.0, 3.0), (2000.0, 5.0)])
        f = math.log(1010.0 / self.e3) / math.log(2000.0 / self.e3)
        y = math.exp(math.log(1.5) + f * math.log(2.5 / 1.5))
        c = r["codes"]["lindhard"]
        self.assertEqual((c["interpolated"], c["missing"]), (1, []))
        self.assertAlmostEqual(c["ratios"][0], y / 3.0, places=12)
        self.assertAlmostEqual(c["ratios"][1], 0.5, places=12)

    def test_shift_below_runs_is_unsupported_not_dropped(self):
        curves = {"lindhard": {1010.0: 1.0, 2000.0: 2.5}, "rustbca": None}
        r = update_docs.scenario_comparison(self.data, "Ar", "Cu", ("b",), curves)
        c = r["codes"]["lindhard"]
        self.assertEqual((len(c["ratios"]), len(c["missing"])), (0, 1))
        self.assertEqual(update_docs.ratio_summary(c, r["n_groups"]), "unavailable (no supported energy)")

    def test_partial_is_marked(self):
        curves = {"lindhard": {2000.0: 2.5, 4000.0: 4.0}, "rustbca": None}
        r = update_docs.scenario_comparison(self.data, "Ar", "Cu", (), curves)
        self.assertEqual(update_docs.ratio_summary(r["codes"]["lindhard"], r["n_groups"]),
                         "0.50 (0.50-0.50), partial: 1 of 2 energies")

    def test_absent_rustbca_and_empty_scenario_are_unavailable(self):
        r = update_docs.scenario_comparison(self.data, "Ar", "Au", ("c",), {"lindhard": {1000.0: 1.0}, "rustbca": None})
        self.assertEqual(r["n_groups"], 0)
        self.assertEqual(update_docs.ratio_summary(r["codes"]["lindhard"], 0), "unavailable (no energies)")
        self.assertEqual(update_docs.ratio_summary(r["codes"]["rustbca"], 0), "unavailable (not run)")

    def test_inputs_not_mutated(self):
        before = copy.deepcopy(self.data)
        update_docs.scenario_comparison(self.data, "Ar", "Cu", ("a",), self.curves)
        self.assertEqual(self.data, before)

    def test_flag_rule_kept(self):
        data = [ds("a", "Cu", [(1000.0, 2.0)]), ds("b", "Cu", [(1005.0, 9.0)])]
        data[1][1]["points"][0]["flag"] = "disagrees_between_compilations"
        r = update_docs.scenario_comparison(data, "Ar", "Cu", (), {"lindhard": {1000.0: 1.0}})
        self.assertEqual([(g["n_points"], g["yield_median"]) for g in r["groups"]], [(1, 2.0)])


class ScenarioTable(unittest.TestCase):
    """update_docs.check_scenarios: unknown, duplicate and wrong-target ids are named (#78)."""

    data = [ds("ar_cu_sputter_x1", "Cu", [(1000.0, 2.0)]), ds("ar_au_sputter_y1", "Au", [(1000.0, 2.0)])]

    def errs(self, **kw):
        s = {"target": "Cu", "exclude": ("ar_cu_sputter_x1",), "caveat": update_docs.CAVEAT_DOC, "reason": "r"}
        s.update(kw)
        return update_docs.check_scenarios([s], self.data)

    def test_valid(self):
        self.assertEqual(self.errs(), [])

    def test_unknown_id(self):
        self.assertTrue(any("unknown dataset id 'ar_cu_sputter_nope'" in e for e in self.errs(exclude=("ar_cu_sputter_nope",))))

    def test_duplicate_id(self):
        self.assertTrue(any("more than once" in e for e in self.errs(exclude=("ar_cu_sputter_x1", "ar_cu_sputter_x1"))))

    def test_wrong_target(self):
        e = self.errs(exclude=("ar_au_sputter_y1",))
        self.assertTrue(any("is an Ar -> Au set, not Cu" in x for x in e), e)

    def test_empty_exclude_and_bad_caveat(self):
        self.assertTrue(any("`exclude` is empty" in e for e in self.errs(exclude=())))
        self.assertTrue(any("`caveat` 'somewhere'" in e for e in self.errs(caveat="somewhere")))
        self.assertTrue(any("which is empty in" in e for e in self.errs(caveat="target_state")))

    def test_duplicate_scenario(self):
        s = {"target": "Cu", "exclude": ("ar_cu_sputter_x1",), "caveat": update_docs.CAVEAT_DOC, "reason": "r"}
        self.assertTrue(any("same target and exclusions" in e for e in update_docs.check_scenarios([s, s], self.data)))

    def test_committed_scenarios_are_valid(self):
        self.assertEqual(update_docs.check_scenarios(update_docs.SENSITIVITY_SCENARIOS, update_docs.load_sputter_sets()), [])


class CommittedBaseline(unittest.TestCase):
    """The no-exclusion scenario reproduces the committed results.json rows and printed means."""

    def test_no_exclusion_reproduces_results(self):
        sp, targets, rb = update_docs.sputter_data()
        if not sp:
            raise unittest.SkipTest("no committed sputter results")
        data = update_docs.load_sputter_sets()
        before = copy.deepcopy(data)
        for t in targets:
            rows = [r for r in sp["rows"] if r["target"] == t]
            curves = {"lindhard": {r["energy_ev"]: r["lindhard_k0"] for r in rows}}
            r = update_docs.scenario_comparison(data, "Ar", t, (), curves)
            self.assertEqual([(g["energy_ev"], g["yield_median"], g["n_points"], g["sets"]) for g in r["groups"]],
                             [(x["energy_ev"], x["yield_median"], x["n_points"], x["sets"]) for x in rows], t)
            got = update_docs.ratio_summary(r["codes"]["lindhard"], r["n_groups"])
            k0 = [x["lindhard_k0_over_median"] for x in rows]
            want = f"{update_docs.gmean(k0):.2f} ({min(k0):.2f}-{max(k0):.2f})"
            self.assertEqual(got, want, t)
        self.assertEqual(data, before)


if __name__ == "__main__":
    unittest.main()
