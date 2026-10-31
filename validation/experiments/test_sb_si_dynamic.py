#!/usr/bin/env python3
"""Tests of sb_si_dynamic.py (#512): the metric, the data loading, the generated
inputs, the committed results, and a reduced-statistics run of the real CLI
(both relaxation conventions, one seed). Standard library only. The reduced run
needs the `lindhard` binary (`$LINDHARD_BIN`, else cargo builds it); it fails,
it is not skipped, without:

    python3 validation/experiments/test_sb_si_dynamic.py
"""

from __future__ import annotations

import json
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sb_si_dynamic as sd  # noqa: E402


def flat_profile(values, width=1.0):
    return [{"centre_nm": width * (i + 0.5), "thickness_nm": width, "sb_cm3": v} for i, v in enumerate(values)]


class Metric(unittest.TestCase):
    def test_peak_is_the_maximum_slab_at_its_centre_with_the_slab_width(self):
        pk = sd.peak(flat_profile([1.0, 4.0, 10.0, 6.0, 1.0]))
        self.assertEqual(pk["depth_nm"], 2.5)
        self.assertEqual(pk["resolution_nm"], 1.0)
        self.assertEqual(pk["max_cm3"], 10.0)

    def test_plateau_lists_the_slabs_near_the_maximum(self):
        pk = sd.peak(flat_profile([1.0, 9.8, 10.0, 9.4, 1.0]))
        self.assertEqual(pk["plateau_nm"], [1.5, 2.5])

    def test_peak_in_the_surface_slab(self):
        pk = sd.peak(flat_profile([10.0, 5.0, 1.0]))
        self.assertEqual(pk["depth_nm"], 0.5)

    def test_empty_engine_profile_is_an_error(self):
        with self.assertRaises(ValueError):
            sd.peak(flat_profile([0.0, 0.0]))

    def test_interpolation_is_normalized_by_the_engines_own_maximum(self):
        # the absolute scale of the engine profile is irrelevant (retained dose is untested)
        a = flat_profile([1.0, 4.0, 10.0, 6.0])
        b = flat_profile([1e20, 4e20, 1e21, 6e20])
        for z in (0.2, 1.0, 2.5, 3.1, 3.9):
            self.assertAlmostEqual(sd.interpolate(a, z), sd.interpolate(b, z))
        self.assertAlmostEqual(sd.interpolate(a, 2.5), 1.0)
        self.assertAlmostEqual(sd.interpolate(a, 2.0), 0.7)  # halfway between 0.4 and 1.0

    def test_interpolation_edges(self):
        p = flat_profile([2.0, 4.0])
        self.assertAlmostEqual(sd.interpolate(p, 0.1), 0.5)  # flat before the first centre
        self.assertEqual(sd.interpolate(p, 5.0), 0.0)  # beyond the extent of the profile

    def test_residual_is_zero_for_the_identical_shape(self):
        markers = [{"depth_nm": z, "normalized_sb_concentration": v} for z, v in ((0.5, 0.2), (1.5, 0.6), (2.5, 1.0))]
        r = sd.shape_residual(flat_profile([2.0, 6.0, 10.0]), markers)
        self.assertAlmostEqual(r["rms"], 0.0)
        self.assertAlmostEqual(r["max_abs"], 0.0)

    def test_residual_rms_and_worst_marker(self):
        markers = [{"depth_nm": 0.5, "normalized_sb_concentration": 0.5},
                   {"depth_nm": 2.5, "normalized_sb_concentration": 0.0}]
        r = sd.shape_residual(flat_profile([0.5, 0.0, 1.0]), markers)  # engine 0.5 and 1.0
        self.assertAlmostEqual(r["rms"], (0.5 * (0.0 + 1.0)) ** 0.5)
        self.assertAlmostEqual(r["max_abs"], 1.0)
        self.assertEqual(r["max_abs_at_depth_nm"], 2.5)

    def test_verdicts_use_the_preset_tolerance(self):
        self.assertEqual(sd.TOL, 0.25)
        pk15 = flat_profile([0.0] * 14 + [1.0])  # peak at 14.5 nm
        m = sd.load_markers()
        c = sd.compare(pk15, m)
        self.assertEqual(c["peak_depth"]["printed_peak"]["verdict"], "pass")
        # 14.5 nm against the largest marker at 19.09 nm: -24 %, inside the tolerance
        self.assertEqual(c["peak_depth"]["marker_maximum"]["verdict"], "pass")
        shallow = sd.compare(flat_profile([1.0]), m)  # peak at 0.5 nm
        self.assertEqual(shallow["peak_depth"]["printed_peak"]["verdict"], "fail")
        self.assertFalse(shallow["peak_in_unmeasured_gap"])
        self.assertTrue(c["peak_in_unmeasured_gap"])


class Regrid(unittest.TestCase):
    def test_regrid_conserves_atoms_per_area_and_handles_uneven_slabs(self):
        slabs = [{"front_nm": 0.0, "back_nm": 0.4, "thickness_nm": 0.4, "sb_cm3": 10.0},
                 {"front_nm": 0.4, "back_nm": 2.5, "thickness_nm": 2.1, "sb_cm3": 20.0}]
        g = sd.regrid(slabs, 2.5, 1.0)
        self.assertEqual(len(g), 3)
        self.assertAlmostEqual(g[0]["sb_cm3"], 10.0 * 0.4 + 20.0 * 0.6)
        self.assertAlmostEqual(g[1]["sb_cm3"], 20.0)
        self.assertAlmostEqual(g[2]["sb_cm3"], 20.0 * 0.5)  # half a cell beyond the last slab counts as zero
        self.assertAlmostEqual(sum(c["sb_cm3"] for c in g), 10.0 * 0.4 + 20.0 * 2.1)

    def test_seed_mean_is_taken_on_the_common_grid(self):
        a = [{"front_nm": 0.0, "back_nm": 2.0, "thickness_nm": 2.0, "sb_cm3": 4.0}]
        b = [{"front_nm": 0.0, "back_nm": 1.0, "thickness_nm": 1.0, "sb_cm3": 2.0},
             {"front_nm": 1.0, "back_nm": 2.0, "thickness_nm": 1.0, "sb_cm3": 6.0}]
        m = sd.mean_profile([a, b])
        self.assertEqual([c["sb_cm3"] for c in m], [3.0, 5.0])
        self.assertEqual([c["centre_nm"] for c in m], [0.5, 1.5])

    def test_slab_at(self):
        slabs = [{"front_nm": 0.0, "back_nm": 0.4, "thickness_nm": 0.4}, {"front_nm": 0.4, "back_nm": 2.5, "thickness_nm": 2.1}]
        self.assertEqual(sd.slab_at(slabs, 1.0), 2.1)


class Data(unittest.TestCase):
    def test_markers_are_the_stored_53_in_depth_order(self):
        m = sd.load_markers()
        self.assertEqual(len(m), 53)
        self.assertEqual([x["depth_nm"] for x in m], sorted(x["depth_nm"] for x in m))
        self.assertEqual(max(x["normalized_sb_concentration"] for x in m), 0.978)

    def test_stored_peak_and_the_unmeasured_gap(self):
        s = sd.stored_peak(sd.load_markers())
        self.assertEqual(s["marker_max_depth_nm"], 19.09)
        self.assertEqual(s["printed_peak_nm"], 15.0)
        self.assertEqual(s["largest_gap_nm"], [13.47, 19.09])

    def test_conditions_match_the_source_metadata(self):
        meta = json.loads((sd.DATA / "titze2021.json").read_text())
        self.assertEqual(meta["conditions"]["energy_keV"] * 1000.0, sd.ENERGY_EV)
        self.assertEqual(meta["conditions"]["fluence_stored"], sd.FLUENCE_CM2)
        self.assertEqual(meta["conditions"]["projectile"], "Sb")

    def test_check_data_passes(self):
        self.assertEqual(sd.check_data(), [])


class Inputs(unittest.TestCase):
    def test_base_input_parses_and_has_no_tuning_or_oxide(self):
        t = tomllib.loads(sd.build_input("ideal-mixing", 1000, 3))
        self.assertEqual(t["beam"]["ion"], "Sb")
        self.assertEqual(t["beam"]["tilt_deg"], 7.0)
        self.assertEqual([layer["material"] for layer in t["target"]["layers"]], ["Si"])
        self.assertNotIn("tuning", t["physics"])
        self.assertEqual(t["dynamic"]["fluence_cm2"], 5.0e16)
        self.assertTrue(t["dynamic"]["erosion"])
        self.assertNotIn("number_density_cm3", t["dynamic"])
        self.assertEqual(t["physics"]["energies"]["Si"]["e_d_ev"], t["physics"]["energies"]["Si"]["e_s_ev"])

    def test_fixed_number_density_input_carries_the_silicon_density(self):
        t = tomllib.loads(sd.build_input("fixed-number-density", 1000, 1))
        self.assertAlmostEqual(t["dynamic"]["number_density_cm3"] / 4.994e22, 1.0, places=3)

    def test_oxide_variant_has_the_measured_thickness_on_top(self):
        t = tomllib.loads(sd.build_input("ideal-mixing", 1000, 1, oxide=True))
        layers = t["target"]["layers"]
        self.assertEqual([(layer["material"], layer["thickness_nm"]) for layer in layers],
                         [("SiO2", 2.9), ("Si", sd.BASE["film_nm"])])

    def test_oxide_oxygen_volume_keeps_the_oxide_density(self):
        t = tomllib.loads(sd.build_input("ideal-mixing", 1000, 1, oxide=True))
        v_o = t["dynamic"]["atomic_volume_nm3"]["O"]
        self.assertTrue(0.008 < v_o < 0.02)
        v_si = 28.085 / (6.02214076e23 * 2.329) * 1e21
        self.assertAlmostEqual((v_si + 2 * v_o) * 6.02214076e23 * 2.2 / 1e21, 60.083, places=2)

    def test_sensitivity_changes_only_what_it_names(self):
        base = tomllib.loads(sd.build_input("ideal-mixing", 1000, 1))
        t = tomllib.loads(sd.build_input("ideal-mixing", 1000, 1, es_scale=2.0))
        self.assertEqual(t["physics"]["energies"]["Si"]["e_s_ev"], 2 * base["physics"]["energies"]["Si"]["e_s_ev"])
        self.assertEqual(t["beam"], base["beam"])
        self.assertEqual(t["dynamic"]["slab_nm"], base["dynamic"]["slab_nm"])

    def test_oxide_is_only_run_with_ideal_mixing(self):
        convs = {lab: cv for lab, _, cv in sd.SENSITIVITY}
        self.assertEqual(convs["oxide=2.9nm SiO2"], ["ideal-mixing"])


class CommittedResults(unittest.TestCase):
    def setUp(self):
        if not sd.RESULTS.exists():
            self.skipTest("no committed results")
        self.r = json.loads(sd.RESULTS.read_text())

    def test_verdicts_follow_the_fixed_tolerance(self):
        for conv, run in self.r["runs"].items():
            self.assertEqual(run["tolerance"], sd.TOL)
            for c in run["peak_depth"].values():
                self.assertEqual(c["verdict"], sd.verdict(c["relative_difference"], sd.TOL))
            self.assertEqual(run["shape_residual"]["verdict"], sd.verdict(run["shape_residual"]["rms"], sd.TOL))

    def test_no_tuning_and_sensitivities_are_reported_for_both_conventions(self):
        self.assertEqual(set(self.r["sensitivity"]), set(sd.CONVENTIONS))
        self.assertNotIn("tuning", json.dumps(self.r["settings"]["physics"]))


class ReducedRun(unittest.TestCase):
    """The real CLI at a fifth of the ions."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.res = sd.reduced_run(2, Path(cls.tmp.name))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_both_conventions_run_and_the_peak_is_inside_the_film(self):
        self.assertEqual(set(self.res), set(sd.CONVENTIONS))
        for r in self.res.values():
            self.assertTrue(2.0 < r["peak"]["depth_nm"] < 60.0)
            self.assertGreater(r["peak"]["max_cm3"], 0.0)
            self.assertEqual(r["shape_residual"]["markers"], 53)
            self.assertGreater(r["recession_nm"], 0.0)
            self.assertGreater(r["sputter_yield"], 0.0)

    def test_profile_is_a_buried_peak_not_a_surface_spike(self):
        for r in self.res.values():
            self.assertGreater(r["peak"]["depth_nm"], 1.5)


if __name__ == "__main__":
    unittest.main()
