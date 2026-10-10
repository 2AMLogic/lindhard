#!/usr/bin/env python3
"""Tests of the one-change-at-a-time runs of se_yield_candidates.py (#242):
each candidate changes one thing, the statistics are batch means over the
seeds, and the committed results and the block in docs/validation.md agree.
Standard library only; run by validation/run.sh, or directly:

    python3 validation/experiments/test_se_yield_candidates.py
"""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import se_yield as sy  # noqa: E402
import se_yield_candidates as sc  # noqa: E402


def changed_lines(material, cand):
    base = sc.make_input(material, sc.BASELINE, 2000, 1, "elf.toml").split("\n")
    other = sc.make_input(material, cand, 2000, 1, "elf.toml").split("\n")
    assert len(base) == len(other)
    return [(a, b) for a, b in zip(base, other) if a != b]


class OneChangeAtATime(unittest.TestCase):
    def test_baseline_is_the_mermin_configuration_at_800_ev(self):
        for m in sc.MATERIALS:
            self.assertEqual(sc.make_input(m, sc.BASELINE, 2000, 1, "elf.toml"),
                             sy.make_input(m, "mermin", 800.0, 2000, 1, "elf.toml"))

    def test_every_candidate_changes_exactly_one_line(self):
        expect = {"elastic-stand-in": "potential = ", "elastic-no-exchange": "exchange = ", "phi-low": "band = ",
                  "phi-high": "band = ", "fermi-tpp2011": "band = ", "barrier-off": "boundary = ",
                  "binding-azzolini": "band = "}
        for m in sc.MATERIALS:
            for cand, start in expect.items():
                diff = changed_lines(m, cand)
                self.assertEqual(len(diff), 1, (m, cand, diff))
                self.assertTrue(diff[0][0].startswith(start) and diff[0][1].startswith(start), (m, cand, diff))

    def test_the_acoustic_row_runs_the_baseline_input(self):
        # Its change is the elastic table the example swaps in, not the input.
        self.assertEqual(changed_lines("Au", "acoustic-phonon"), [])

    def test_the_acoustic_row_is_only_for_materials_of_the_published_table(self):
        self.assertTrue(sc.applies("Au", "acoustic-phonon"))
        self.assertFalse(sc.applies("Cu", "acoustic-phonon"))
        with self.assertRaises(ValueError):
            sc.make_input("Cu", "acoustic-phonon", 2000, 1, "elf.toml")

    def test_stand_in_keeps_the_exchange_correction(self):
        text = sc.make_input("Au", "elastic-stand-in", 2000, 1, "elf.toml")
        self.assertIn('potential = "thomas-fermi-yukawa"', text)
        self.assertIn("exchange = true", text)

    def test_work_function_rows_use_the_ends_of_the_cited_range(self):
        for m in sc.MATERIALS:
            w = sy.BAND[m]["work_function_ev"]
            self.assertIn(f"work_function_ev = {w['low']!r}", sc.make_input(m, "phi-low", 2000, 1, "elf.toml"))
            self.assertIn(f"work_function_ev = {w['high']!r}", sc.make_input(m, "phi-high", 2000, 1, "elf.toml"))

    def test_fermi_row_is_the_printed_value_with_the_mid_work_function(self):
        self.assertEqual(sc.FERMI_TPP2011_EV, {"Cu": 8.7, "Au": 9.0})
        for m in sc.MATERIALS:
            text = sc.make_input(m, "fermi-tpp2011", 2000, 1, "elf.toml")
            self.assertIn('kind = "metal"', text)
            self.assertIn(f"fermi_ev = {sc.FERMI_TPP2011_EV[m]!r}", text)
            self.assertIn(f"work_function_ev = {sy.BAND[m]['work_function_ev']['mid']!r}", text)
            self.assertIn("doi:10.1002/sia.3522", text)

    def test_binding_row_adds_only_the_printed_binding_to_the_baseline_band(self):
        self.assertEqual(sc.AZZOLINI_BINDING_EV, {"Cu": 7.726, "Au": 9.226})
        for m in sc.MATERIALS:
            (base, new), = changed_lines(m, "binding-azzolini")
            b = sc.AZZOLINI_BINDING_EV[m]
            self.assertIn(f"valence_binding_ev = {b!r}", new)
            self.assertIn("arXiv:1809.00859v1", new)
            # Everything else of the baseline band is kept: kind, valence count, mid work function.
            self.assertEqual(new.replace(f"valence_binding_ev = {b!r}, ", "").replace(
                f"{sc.AZZOLINI_PROVENANCE}; ", ""), base)

    def test_no_exchange_keeps_the_dhfs_potential(self):
        text = sc.make_input("Au", "elastic-no-exchange", 2000, 1, "elf.toml")
        self.assertIn('potential = "salvat-dhfs"', text)
        self.assertIn("exchange = false", text)

    def test_the_pre_149_row_is_both_elastic_changes_and_is_not_a_candidate(self):
        diff = changed_lines("Cu", "elastic-pre-149")
        self.assertEqual(sorted(b for _, b in diff), ["exchange = false", 'potential = "thomas-fermi-yukawa"'])
        self.assertEqual(sc.CANDIDATES["elastic-pre-149"][0], "context")

    def test_barrier_off_keeps_the_band(self):
        text = sc.make_input("Cu", "barrier-off", 2000, 1, "elf.toml")
        self.assertIn('boundary = "transparent"', text)
        self.assertIn("band = ", text)

    def test_seed_and_histories_reach_the_input(self):
        text = sc.make_input("Au", "phi-low", 123, 7, "elf.toml")
        self.assertIn("histories = 123", text)
        self.assertIn("seed = 7", text)

    def test_unknown_candidate(self):
        with self.assertRaises(ValueError):
            sc.make_input("Au", "tuned", 2000, 1, "elf.toml")


def fake(runs):
    return {"software": "test", "runs": [
        {"material": m, "candidate": c, "seed": s, "histories": 2000, "delta": d, "eta": 0.5, "energy_ev": 800.0}
        for (m, c, s, d) in runs]}


class Statistics(unittest.TestCase):
    def test_mean_and_standard_error(self):
        m, se = sc.mean_se([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(m, 2.5)
        self.assertAlmostEqual(se, math.sqrt(5.0 / 3.0) / 2.0)

    def test_one_value_has_no_error(self):
        m, se = sc.mean_se([2.0])
        self.assertEqual(m, 2.0)
        self.assertTrue(math.isnan(se))

    def test_change_is_paired_by_seed(self):
        # The candidate is the baseline plus 0.5 on every seed: the change has no scatter although both rows do.
        r = fake([("Au", "baseline", s, 2.0 + 0.1 * s) for s in (1, 2, 3)]
                 + [("Au", "phi-low", s, 2.5 + 0.1 * s) for s in (1, 2, 3)])
        s = sc.summarize(r, "Au", "phi-low")
        self.assertAlmostEqual(s["change"], 0.5)
        self.assertAlmostEqual(s["change_se"], 0.0)
        self.assertAlmostEqual(s["change_rel"], 0.5 / 2.2)
        self.assertGreater(s["delta_se"], 0.0)
        self.assertNotIn("change", sc.summarize(r, "Au", "baseline"))
        self.assertIsNone(sc.summarize(r, "Cu", "baseline"))

    def test_missing_runs_are_named(self):
        r = fake([("Au", "baseline", 1, 2.0)])
        miss = sc.missing_runs(r)
        self.assertIn("Au baseline seed 2: not run", miss)
        self.assertNotIn("Au baseline seed 1: not run", miss)
        self.assertFalse(any(x.startswith("Cu acoustic-phonon") for x in miss))

    def test_control_must_reproduce_the_baseline(self):
        r = fake([("Au", "baseline", 1, 2.0), ("Au", "acoustic-phonon", 1, 3.0)])
        r["runs"][1]["control"] = {"delta": 2.0, "eta": 0.5, "elastic_below_100_ev_per_primary": 1.0}
        self.assertEqual(sc.control_errors(r), [])
        r["runs"][1]["control"]["delta"] = 2.001
        self.assertEqual(len(sc.control_errors(r)), 1)


class Committed(unittest.TestCase):
    def setUp(self):
        if not sc.RESULTS.exists():
            self.skipTest("no committed candidate results")
        self.results = json.loads(sc.RESULTS.read_text())

    def test_results_are_complete(self):
        self.assertEqual(sc.missing_runs(self.results), [])
        self.assertEqual(sc.control_errors(self.results), [])

    def test_every_run_is_at_the_reference_energy(self):
        self.assertTrue(all(r["energy_ev"] == sc.ENERGY_EV for r in self.results["runs"]))

    def test_docs_block_is_the_markdown_output(self):
        self.assertEqual(sc.doc_block(), sc.markdown(sy.load_datasets(), self.results))


if __name__ == "__main__":
    unittest.main()
