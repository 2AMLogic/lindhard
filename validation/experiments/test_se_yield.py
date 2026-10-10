#!/usr/bin/env python3
"""Tests of the SE-yield reference statistics, the dataset checks and the
initial bounds in se_yield.py (#149). Standard library only; run by
validation/run.sh, or directly:

    python3 validation/experiments/test_se_yield.py
"""

from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import se_yield as sy  # noqa: E402


def ds(ident, pts, target="Al"):
    return {"id": ident, "target": target, "points": [{"energy_ev": e, "yield": y} for e, y in pts]}


class ReferenceRule(unittest.TestCase):
    def test_interior_maximum_resolves(self):
        d = ds("a", [(100, 1), (200, 2), (300, 3), (400, 2.5), (500, 2)])
        self.assertEqual(sy.dataset_peak(d), (300, 3))

    def test_maximum_at_an_end_does_not_resolve(self):
        self.assertIsNone(sy.dataset_peak(ds("a", [(100, 5), (200, 4), (300, 3), (400, 2), (500, 1)])))
        self.assertIsNone(sy.dataset_peak(ds("a", [(100, 1), (200, 2), (300, 3), (400, 4), (500, 5)])))

    def test_interior_maximum_outside_the_compared_range_does_not_resolve(self):
        d = ds("a", [(10000, 0.15), (15000, 0.27), (20000, 0.1), (25000, 0.06), (30000, 0.05)])
        self.assertIsNone(sy.dataset_peak(d))
        d = ds("a", [(20, 0.5), (50, 1.0), (80, 1.2), (150, 1.1), (300, 0.9)])
        self.assertIsNone(sy.dataset_peak(d))
        d = ds("a", [(50, 0.5), (100, 1.2), (150, 1.1), (300, 0.9), (400, 0.8)])
        self.assertEqual(sy.dataset_peak(d), (100, 1.2))

    def test_too_few_points_do_not_resolve(self):
        self.assertIsNone(sy.dataset_peak(ds("a", [(100, 1), (200, 3), (300, 1)])))

    def test_tie_takes_the_lowest_energy(self):
        d = ds("a", [(100, 1), (200, 3), (300, 3), (400, 2), (500, 1)])
        self.assertEqual(sy.dataset_peak(d), (200, 3))

    def test_unsorted_input_is_sorted(self):
        d = ds("a", [(500, 2), (100, 1), (300, 3), (200, 2), (400, 2.5)])
        self.assertEqual(sy.dataset_peak(d), (300, 3))

    def test_median_of_two_is_their_mean_and_sets_stay_separate(self):
        a = ds("a", [(100, 1), (200, 2), (300, 3), (400, 2), (500, 1)])
        b = ds("b", [(100, 1), (400, 4), (500, 1), (600, 0.5), (700, 0.2)])
        r = sy.reference_stats([a, b, ds("c", [(100, 1)])], "Al")
        self.assertEqual(r["n_sets"], 3)
        self.assertEqual([x["id"] for x in r["resolving"]], ["a", "b"])
        self.assertEqual(r["dmax_median"], 3.5)
        self.assertEqual(r["emax_median_ev"], 350)
        self.assertEqual(r["not_resolving"], ["c"])


class Bounds(unittest.TestCase):
    ref = {"emax_median_ev": 400.0, "dmax_median": 1.0}

    def test_inside(self):
        b = sy.evaluate_bounds(self.ref, [(100, 0.5), (200, 1.4), (400, 1.0)])
        self.assertTrue(b["emax_ok"] and b["dmax_ok"])
        self.assertEqual(b["emax_ev"], 200)

    def test_energy_factor_edges(self):
        self.assertTrue(sy.evaluate_bounds(self.ref, [(200, 1.0), (300, 0.5)])["emax_ok"])
        self.assertTrue(sy.evaluate_bounds(self.ref, [(800, 1.0), (300, 0.5)])["emax_ok"])
        self.assertFalse(sy.evaluate_bounds(self.ref, [(100, 1.0), (300, 0.5)])["emax_ok"])
        self.assertFalse(sy.evaluate_bounds(self.ref, [(1000, 1.0), (300, 0.5)])["emax_ok"])

    def test_yield_deviation(self):
        self.assertTrue(sy.evaluate_bounds(self.ref, [(400, 1.5)])["dmax_ok"])
        self.assertFalse(sy.evaluate_bounds(self.ref, [(400, 1.51)])["dmax_ok"])
        self.assertFalse(sy.evaluate_bounds(self.ref, [(400, 0.49)])["dmax_ok"])

    def test_nothing_to_compare(self):
        self.assertIsNone(sy.evaluate_bounds({}, [(400, 1.0)]))
        self.assertIsNone(sy.evaluate_bounds(self.ref, []))


class RunInputs(unittest.TestCase):
    """The generated run inputs (#149 remaining scope 1: DHFS elastic, the Si band)."""

    def inp(self, material, cfg="default"):
        return sy.make_input(material, cfg, 400.0, 10, 1, "elf.toml")

    def test_every_material_runs_on_dhfs_with_exchange(self):
        for m in sy.OPTICAL_ELF:
            t = self.inp(m)
            self.assertIn('potential = "salvat-dhfs"', t, m)
            self.assertIn("exchange = true", t, m)
            self.assertNotIn("thomas-fermi-yukawa", t, m)
            # Seltzer's polarization cutoff needs every table energy above 50 eV; the tables start at 5 eV.
            self.assertNotIn("correlation_polarization", t, m)
            self.assertIn("min_energy_ev = 5.0", t, m)

    def test_si_is_run_with_its_elf_and_an_insulator_band(self):
        self.assertEqual(sy.OPTICAL_ELF["Si"], "si_elf_yang2019.toml")
        self.assertTrue((sy.OPTICAL / sy.OPTICAL_ELF["Si"]).is_file())
        t = self.inp("Si")
        self.assertIn('band = { kind = "insulator", valence_band_width_ev = 12.36, band_gap_ev = 1.1, '
                      'affinity_ev = 4.05, provenance = "', t)
        self.assertIn('boundary = "step-barrier"', t)
        self.assertNotIn("work_function_ev", t)

    def test_si_gap_and_affinity_are_the_cited_band_defaults(self):
        # lindhard/src/electron/boundary.rs, BAND_DEFAULTS (#115): Si 1.1 eV gap, 4.05 eV affinity.
        src = (sy.REPO / "lindhard" / "src" / "electron" / "boundary.rs").read_text()
        si = src[src.index('material: "Si"'):src.index('material: "SiO2"')]
        self.assertIn("value_ev: 1.1,", si)
        self.assertIn("value_ev: 4.05,", si)
        self.assertEqual((sy.BAND["Si"]["band_gap_ev"], sy.BAND["Si"]["affinity_ev"]), (1.1, 4.05))

    def test_every_band_value_names_its_source(self):
        for m, b in sy.BAND.items():
            p = sy.BAND_PROVENANCE[m]
            self.assertTrue(p.strip(), m)
        p = sy.BAND_PROVENANCE["Si"]
        for s in ("Robertson and R. M. Wallace", "Chelikowsky and M. L. Cohen", "10.1103/PhysRevB.10.5095",
                  "Table II", "12.36"):
            self.assertIn(s, p)
        prov = sy.PROVENANCE.read_text()
        self.assertTrue("Chelikowsky" in prov and "12.36" in prov,
                        "docs/data-provenance.md has no row for the Si valence band width")

    def test_work_function_configs_are_metal_only(self):
        self.assertFalse(sy.applies("Si", "phi-low"))
        self.assertFalse(sy.applies("Si", "phi-high"))
        self.assertTrue(sy.applies("Si", "default"))
        self.assertTrue(sy.applies("Al", "phi-low"))
        with self.assertRaises(ValueError):
            self.inp("Si", "phi-low")

    def test_barrier_off_keeps_the_band_for_the_secondaries(self):
        t = self.inp("Si", "barrier-off")
        self.assertIn('boundary = "transparent"', t)
        self.assertIn('kind = "insulator"', t)
        t = self.inp("Cu", "barrier-off")
        self.assertIn("work_function_ev = 4.815", t)

    def test_band_bottom_cutoff_is_one_ev_above_the_fermi_level(self):
        si = sy.BAND["Si"]
        fermi = sy.band_fermi_ev({"kind": "insulator", "valence_band_width_ev": si["valence_band_width_ev"],
                                  "band_gap_ev": si["band_gap_ev"], "affinity_ev": si["affinity_ev"]})
        self.assertAlmostEqual(fermi, 12.91, places=9)
        self.assertTrue(-1e-9 <= sy.BAND_BOTTOM_CUTOFF_EV["Si"] - (fermi + 1.0) < 0.01)
        self.assertEqual(sy.band_fermi_ev({"kind": "metal", "fermi_ev": 7.0}), 7.0)
        t = self.inp("Si", "cutoff-band-bottom")
        self.assertIn("cutoff_ev = 13.91", t)
        self.assertIn('cutoff_reference = "band-bottom"', t)


class DatasetChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = sorted(sy.DATA.glob("*.json"))[0]
        cls.good = json.loads(cls.path.read_text())
        cls.prov = sy.PROVENANCE.read_text()

    def errs(self, d):
        return sy.check_dataset(self.path, d, self.prov)

    def test_committed_file_passes(self):
        self.assertEqual(self.errs(self.good), [])

    def test_missing_incidence_note(self):
        d = copy.deepcopy(self.good)
        d["incidence_note"] = ""
        self.assertIn("`incidence_deg` is null without an `incidence_note`", self.errs(d))

    def test_zero_yield(self):
        d = copy.deepcopy(self.good)
        d["points"][0]["yield"] = 0
        self.assertIn("`points[0].yield` must be a positive number", self.errs(d))

    def test_unknown_in_provenance(self):
        d = copy.deepcopy(self.good)
        self.assertTrue(any("no row in docs/data-provenance.md" in e for e in sy.check_dataset(self.path, d, "")))

    def test_every_committed_set_is_one_original_measurement(self):
        seen = set()
        for f in sy.DATA.glob("*.json"):
            d = json.loads(f.read_text())
            key = (d["target"], d["compilation_reference_number"], d["id"].endswith("amorphous_joy1996"))
            self.assertNotIn(key, seen, f.name)
            seen.add(key)


if __name__ == "__main__":
    unittest.main()
