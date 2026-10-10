#!/usr/bin/env python3
"""Unit tests of the report conversion in run_electron.py that need no
simulation: summing `results.table_coverage` over batches (#253).

    python3 validation/oracles/test_run_electron.py
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_electron  # noqa: E402


def layer(i, el, inel, bounds=((10.0, 1000.0), (20.0, 2000.0))):
    def ch(c, b):
        return {"energy_min_ev": b[0], "energy_max_ev": b[1], "below": c[0], "within": c[1], "above": c[2]}

    return {"layer": i, "elastic": ch(el, bounds[0]), "inelastic": ch(inel, bounds[1])}


class TableCoverage(unittest.TestCase):
    def test_batches_sum_per_layer_and_channel(self):
        acc = run_electron.Accumulator(2, 100.0, 10)
        acc.add_table_coverage([layer(0, (1, 2, 3), (0, 6, 0)), layer(1, (0, 4, 0), (4, 0, 0))])
        acc.add_table_coverage([layer(0, (10, 20, 30), (0, 60, 0)), layer(1, (0, 1, 0), (0, 0, 1))])
        cov = acc.result()["table_coverage"]
        self.assertEqual(cov[0]["elastic"]["below"], 11)
        self.assertEqual(cov[0]["elastic"]["within"], 22)
        self.assertEqual(cov[0]["elastic"]["above"], 33)
        self.assertEqual(cov[0]["inelastic"]["within"], 66)
        self.assertEqual(cov[1]["inelastic"], {"energy_min_ev": 20.0, "energy_max_ev": 2000.0,
                                               "below": 4, "within": 0, "above": 1})
        self.assertEqual(cov[1]["layer"], 1)

    def test_the_first_batch_is_not_modified(self):
        first = [layer(0, (1, 2, 3), (0, 6, 0))]
        acc = run_electron.Accumulator(2, 100.0, 10)
        acc.add_table_coverage(first)
        acc.add_table_coverage([layer(0, (1, 1, 1), (1, 1, 1))])
        self.assertEqual(first[0]["elastic"]["below"], 1)

    def test_a_summary_without_it_gives_none(self):
        for parts in ([None, [layer(0, (1, 1, 1), (1, 1, 1))]], [[layer(0, (1, 1, 1), (1, 1, 1))], None]):
            acc = run_electron.Accumulator(2, 100.0, 10)
            for p in parts:
                acc.add_table_coverage(p)
            self.assertIsNone(acc.result()["table_coverage"])

    def test_bounds_must_agree(self):
        acc = run_electron.Accumulator(2, 100.0, 10)
        acc.add_table_coverage([layer(0, (1, 1, 1), (1, 1, 1))])
        with self.assertRaises(SystemExit):
            acc.add_table_coverage([layer(0, (1, 1, 1), (1, 1, 1), bounds=((10.0, 999.0), (20.0, 2000.0)))])


if __name__ == "__main__":
    unittest.main()
