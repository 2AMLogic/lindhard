#!/usr/bin/env python3
"""Tests of the input generation of backscatter.py (#148): the baseline
variants and the fast-secondary sensitivity inputs. Standard library only;
no simulation is run:

    python3 validation/experiments/test_backscatter.py
"""

from __future__ import annotations

import json
import sys
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


if __name__ == "__main__":
    unittest.main()
