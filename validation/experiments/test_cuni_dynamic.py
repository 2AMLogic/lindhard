#!/usr/bin/env python3
"""Tests of cuni_dynamic.py (#448): the conversions, the committed data and
results, the generated inputs, and a reduced-statistics run of the real CLI
(Ho `t2-r4`, both relaxation conventions, one seed, a sixth of the ions).
Standard library only. The reduced run needs the `lindhard` binary
(`$LINDHARD_BIN`, else cargo builds it); it fails, it is not skipped, without:

    python3 validation/experiments/test_cuni_dynamic.py
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cuni_dynamic as cd  # noqa: E402


class Conversions(unittest.TestCase):
    def test_surface_fraction_formula_matches_the_readme_value(self):
        # alloys/README.md: about 0.2683 for t2-r4 (c = 0.384, q = 1.7)
        self.assertAlmostEqual(cd.surface_fraction_from_ratio(0.384, 1.7), 0.2683, places=4)
        # q = 1 means no change of the AES ratio: m = c
        self.assertAlmostEqual(cd.surface_fraction_from_ratio(0.384, 1.0), 0.384)

    def test_yield_from_erosion(self):
        # 1 angstrom/min at 1e13 ions/cm2/s in a 1e23 /cm3 solid:
        # V = 1e-8/60 cm/s, Y = 1e23 * 1.6667e-10 / 1e13 = 1.6667
        self.assertAlmostEqual(cd.yield_from_erosion(1e23, 1.0, 1e13), 1e23 * 1e-8 / 60 / 1e13)
        self.assertAlmostEqual(cd.yield_from_erosion(1e23, 60.0, 1e13), 100.0)

    def test_number_density_is_between_the_elements(self):
        n_cu = 1e21 / cd.atomic_volume_nm3("Cu")
        n_ni = 1e21 / cd.atomic_volume_nm3("Ni")
        self.assertAlmostEqual(cd.alloy_number_density_cm3(1.0), n_cu)
        self.assertAlmostEqual(cd.alloy_number_density_cm3(0.0), n_ni)
        self.assertTrue(min(n_cu, n_ni) < cd.alloy_number_density_cm3(0.384) < max(n_cu, n_ni))

    def test_ho_t2r4_conversion(self):
        e = cd.ho_experiment(cd.load_ho()[cd.HO_ROW])
        a = e["yield_conventions"]["flux_per_sample_area"]["yield"]
        b = e["yield_conventions"]["flux_per_beam_normal_area"]["yield"]
        self.assertAlmostEqual(a, e["number_density_cm3"] * 23.0e-8 / 60 / 9.4e13)
        self.assertAlmostEqual(b / a, 1 / math.cos(math.radians(20.0)))
        s = e["surface_cu_fraction"]
        self.assertLess(s["lo"], s["value"])
        self.assertLess(s["value"], s["hi"])

    def test_seah_dench_lengths_are_of_the_order_of_a_nanometre(self):
        for sym, ev in cd.AES_LINES_EV.items():
            self.assertTrue(0.8 < cd.seah_dench_nm(sym, ev) < 2.0, sym)

    def test_verdict_uses_the_fixed_tolerances(self):
        self.assertEqual((cd.TOL_SURFACE_FRACTION, cd.TOL_TOTAL_YIELD), (0.20, 0.30))
        self.assertEqual(cd.verdict(-0.30, 0.30), "pass")
        self.assertEqual(cd.verdict(0.31, 0.30), "fail")


class Aes(unittest.TestCase):
    @staticmethod
    def slab(i, x1, x2, cu, ni):
        return {"slab": i, "front_nm": x1, "back_nm": x2, "thickness_nm": x2 - x1,
                "atoms_per_cm2_Cu": cu, "atoms_per_cm2_Ni": ni}

    def test_uniform_bulk_returns_the_bulk_fraction(self):
        slabs = [self.slab(i, 0.5 * i, 0.5 * (i + 1), 0.384, 0.616) for i in range(80)]
        for scale in (0.5, 1.0, 2.0):
            self.assertAlmostEqual(cd.aes_surface_fraction(slabs, scale), 0.384, places=3)
        self.assertAlmostEqual(cd.aes_surface_fraction(slabs, None), 0.384)

    def test_depleted_front_slab_is_seen_less_with_a_longer_length(self):
        slabs = [self.slab(0, 0.0, 0.5, 0.2, 0.8)] + [self.slab(i, 0.5 * i, 0.5 * (i + 1), 0.384, 0.616)
                                                      for i in range(1, 80)]
        first = cd.aes_surface_fraction(slabs, None)
        near = cd.aes_surface_fraction(slabs, 0.5)
        far = cd.aes_surface_fraction(slabs, 2.0)
        self.assertAlmostEqual(first, 0.2)
        self.assertTrue(first < near < far < 0.384)

    def test_empty_thickness_slab_does_not_divide_by_zero(self):
        slabs = [self.slab(0, 0.0, 0.0, 1e-3, 1e-3), self.slab(1, 0.0, 5.0, 0.384, 0.616)]
        self.assertTrue(0 < cd.aes_surface_fraction(slabs, 1.0) < 1)


class Inputs(unittest.TestCase):
    def parse(self, relaxation, ar_mode, tilt=20.0):
        return tomllib.loads(cd.build_input(0.384, tilt, 2000.0, relaxation, ar_mode, 1000, 7))

    def test_inputs_parse_and_carry_the_case(self):
        t = self.parse("ideal-mixing", "ar-retained")
        self.assertEqual(t["beam"]["ion"], "Ar")
        self.assertEqual(t["beam"]["tilt_deg"], 20.0)
        self.assertEqual(t["run"]["seed"], 7)
        self.assertTrue(t["dynamic"]["erosion"])
        self.assertEqual(t["dynamic"]["relaxation"], "ideal-mixing")
        self.assertIn("Ar", t["dynamic"]["atomic_volume_nm3"])
        self.assertNotIn("retention_limit", t["dynamic"])

    def test_fixed_number_density_and_argon_modes(self):
        t = self.parse("fixed-number-density", "ar-re-emitted", tilt=0.0)
        self.assertEqual(t["dynamic"]["number_density_cm3"], cd.alloy_number_density_cm3(0.384))
        self.assertNotIn("atomic_volume_nm3", t["dynamic"])
        self.assertEqual(t["dynamic"]["retention_limit"]["Ar"], cd.AR_MODES["ar-re-emitted"])

    def test_only_the_expected_keys_differ_between_conventions(self):
        a, b = self.parse("ideal-mixing", "ar-retained"), self.parse("fixed-number-density", "ar-retained")
        for k in ("beam", "target", "materials", "physics", "run"):
            self.assertEqual(a[k], b[k], k)

    def test_no_tuning_and_matched_energies(self):
        t = self.parse("ideal-mixing", "ar-retained")
        self.assertNotIn("tuning", t["physics"])
        for sym in ("Cu", "Ni"):
            e = t["physics"]["energies"][sym]
            self.assertEqual((e["e_d_ev"], e["e_b_ev"], e["e_s_ev"]), (e["e_s_ev"], 0.0, cd.ELEMENTS[sym]["e_s_ev"]))


class Committed(unittest.TestCase):
    def test_data_and_results_are_consistent(self):
        self.assertEqual(cd.check_data(), [])

    def test_excluded_pair_is_not_used(self):
        r = json.loads(cd.RESULTS.read_text())
        self.assertEqual(r["ho1976"]["experiment"]["row_id"], "t2-r4")
        self.assertEqual(r["ho_excluded_rows"], ["t1-r3", "t2-r1"])

    def test_unknown_betz_uncertainty_is_not_zero(self):
        meta = json.loads((cd.DATA / "betz1980.json").read_text())
        self.assertIsNone(meta["uncertainty"]["experimental_yield_uncertainty"])
        r = json.loads(cd.RESULTS.read_text())
        for p in r["betz1980"]["points"]:
            self.assertGreater(p["digitization_bound"], 0)
            self.assertNotIn("experimental_uncertainty", p)

    def test_committed_verdicts_follow_the_fixed_tolerances(self):
        r = json.loads(cd.RESULTS.read_text())
        for entry in r["ho1976"]["comparison"].values():
            for group, tol in (("yield", cd.TOL_TOTAL_YIELD), ("surface_cu_fraction", cd.TOL_SURFACE_FRACTION)):
                for c in entry[group].values():
                    self.assertEqual(c["tolerance"], tol)
                    self.assertEqual(c["verdict"], cd.verdict(c["relative_difference"], tol))
        for p in r["betz1980"]["points"]:
            for c in p["runs"].values():
                self.assertEqual(c["verdict"], cd.verdict(c["relative_difference"], cd.TOL_TOTAL_YIELD))


class Diagnostics(unittest.TestCase):
    """The #465 diagnosis: inputs, interpolation and the committed findings."""

    def test_static_inputs_parse_and_have_no_dynamic_block(self):
        pure = tomllib.loads(cd.build_static_input(1.0, 20.0, 2000.0, 1000, 3))
        self.assertEqual(pure["target"]["substrate"], "Cu")
        self.assertEqual(set(pure["physics"]["energies"]), {"Cu"})
        alloy = tomllib.loads(cd.build_static_input(0.384, 0.0, 2000.0, 1000, 3))
        self.assertEqual(set(alloy["physics"]["energies"]), {"Cu", "Ni"})
        self.assertEqual(alloy["materials"]["CuNi"]["density_g_cm3"], cd.alloy_density_g_cm3(0.384))
        for t in (pure, alloy):
            self.assertNotIn("dynamic", t)
            self.assertNotIn("tuning", t["physics"])
            self.assertTrue(t["physics"]["follow_recoils"])

    def test_surface_binding_override_is_for_the_mechanism_test_only(self):
        t = tomllib.loads(cd.build_static_input(0.3, 0.0, 2000.0, 10, 1, e_s_ev={"Ni": 3.49}))
        self.assertEqual(t["physics"]["energies"]["Ni"]["e_s_ev"], 3.49)
        self.assertEqual(t["physics"]["energies"]["Cu"]["e_s_ev"], cd.ELEMENTS["Cu"]["e_s_ev"])

    def test_stored_cu_band_and_the_missing_ni_data(self):
        band = cd.stored_cu_band(2000.0, 0.02)
        self.assertGreater(band["points"], 0)
        self.assertLessEqual(band["min"], band["median"])
        self.assertLessEqual(band["median"], band["max"])
        self.assertEqual(cd.count_stored_sets("ar_ni"), 0)
        self.assertGreater(cd.count_stored_sets("ar_cu"), 0)

    def test_betz_interpolation_stays_between_the_bracketing_markers(self):
        r = cd.betz_interpolated_yield(38.4, cd.load_betz())
        lo, hi = sorted(r["bracket_yields"])
        self.assertLessEqual(lo, r["value"])
        self.assertLessEqual(r["value"], hi)
        self.assertLessEqual(r["bracket_at_percent"][0], 38.4)
        self.assertLessEqual(38.4, r["bracket_at_percent"][1])

    def test_required_surface_fraction_interpolates_the_flux_fraction(self):
        grid = [{"cu_fraction": c, "cu_fraction_of_flux": {"value": f}} for c, f in ((0.2, 0.3), (0.3, 0.5))]
        self.assertAlmostEqual(cd.required_surface_fraction(grid, 0.4), 0.25)
        self.assertIsNone(cd.required_surface_fraction(grid, 0.9))

    def test_committed_findings_are_reported_not_tuned(self):
        r = json.loads(cd.RESULTS.read_text())
        d = r["diagnostics"]
        self.assertEqual(d["stored_data"]["ar_ni_sets"], 0)
        self.assertIsNone(d["pure_ni_over_stored"])
        self.assertEqual(r["settings"]["surface_binding_ev"], {s: e["e_s_ev"] for s, e in cd.ELEMENTS.items()})
        for tag in ("ho_20deg", "betz_0deg"):
            self.assertEqual(len(d[tag]["grid"]), len(cd.DIAG_COMPOSITIONS))
            for g in d[tag]["grid"]:
                # alloy yield is additive in the pure-element yields at fixed composition
                self.assertAlmostEqual(g["total_over_additive"], 1.0, delta=0.02)


class StaticRun(unittest.TestCase):
    """The real CLI, one seed, a few thousand ions (about 1 s)."""

    def test_static_yields_have_the_expected_structure(self):
        binary = cd.lindhard_cli.lindhard_binary()
        r = cd.static_yields(binary, 2, 0.384, 20.0, seeds=[1], ions=5000)
        self.assertAlmostEqual(r["y_total"]["value"], r["y_cu"]["value"] + r["y_ni"]["value"])
        self.assertTrue(1.0 < r["y_total"]["value"] < 6.0)
        # Cu, the more weakly bound element, is over-represented in the flux
        self.assertGreater(r["cu_fraction_of_flux"]["value"], 0.384)


class ReducedRun(unittest.TestCase):
    """The real CLI at a sixth of the ions (about 10 s)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.res = cd.reduced_run(2, Path(cls.tmp.name))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_both_conventions_run_and_are_physical(self):
        self.assertEqual(set(self.res["runs"]), {f"{c}/ar-re-emitted" for c in cd.CONVENTIONS})
        for r in self.res["runs"].values():
            self.assertTrue(1.0 < r["yield"]["value"] < 6.0)
            self.assertTrue(0.0 < r["yield"]["se"] < 0.5)
            for v in r["surface_cu_fraction"].values():
                self.assertTrue(0.0 < v["value"] < 0.6)

    def test_surface_is_depleted_in_the_preferentially_sputtered_cu(self):
        # Cu has the lower binding energy and is the preferentially sputtered
        # component: the surface (first slab) holds less Cu than the bulk.
        for r in self.res["runs"].values():
            self.assertLess(r["surface_cu_fraction"]["first_slab"]["value"], 0.384)

    def test_comparison_reports_both_flux_conventions_and_all_depth_mappings(self):
        for entry in self.res["comparison"].values():
            self.assertEqual(set(entry["yield"]), {"flux_per_sample_area", "flux_per_beam_normal_area"})
            self.assertEqual(set(entry["surface_cu_fraction"]), {"first_slab", "aes_x0.5", "aes_x1", "aes_x2"})

    def test_initial_density_matches_the_engines(self):
        # The harness's number density is the one the engine builds from the
        # material (first row of the composition profile, step 0).
        import csv
        out = Path(self.tmp.name) / "reduced_ideal-mixing" / "out"
        with open(out / "dynamic_composition.csv", newline="") as f:
            row = next(r for r in csv.DictReader(f) if r["step"] == "0" and r["slab"] == "0")
        n = (float(row["atoms_per_cm2_Cu"]) + float(row["atoms_per_cm2_Ni"])) / (float(row["thickness_nm"]) * 1e-7)
        self.assertAlmostEqual(n / cd.alloy_number_density_cm3(0.384), 1.0, places=6)


if __name__ == "__main__":
    unittest.main()
