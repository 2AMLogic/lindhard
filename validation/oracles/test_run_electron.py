#!/usr/bin/env python3
"""Unit tests of run_electron.py that need no simulation: summing
`results.table_coverage` over batches (#253), and the identity each electron
oracle records in a committed summary (#778).

    python3 validation/oracles/test_run_electron.py
"""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import oracle  # noqa: E402
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


def caps(prim=0, tracks=0, aff=0, ev=0.0):
    return {"primary_event_capped": prim, "secondary_tracks": tracks, "affected_histories": aff, "event_cap_ev": ev}


def results(prim=0, tracks=0, aff=0, ev=0.0):
    return {"fates": {"event_capped": prim}, "event_caps": {"secondary_tracks": tracks, "affected_histories": aff},
            "budget": {"event_cap_ev": ev}}


class EventCaps(unittest.TestCase):
    def pooled(self, rs, ns):
        acc = run_electron.Accumulator(len(rs), 100.0, 10)
        for i, (r, n) in enumerate(zip(rs, ns)):
            acc.n[i] = n
            if r is not None:
                acc.add_event_caps(i, r)
        return acc.result()["event_cap_diagnostics"]

    def test_one_capped_secondary_of_a_completed_primary(self):
        d = self.pooled([results(0, 1, 1, 3.5)], [100])
        self.assertEqual((d["primary_event_capped"], d["secondary_tracks"], d["affected_histories"]), (0, 1, 1))
        self.assertEqual(d["status"], "truncated")
        self.assertIn("TRUNCATED", run_electron.truncation_note(d))

    def test_affected_histories_is_not_inflated_to_the_track_count(self):
        d = self.pooled([results(0, 7, 2, 1.0)], [50])
        self.assertEqual((d["secondary_tracks"], d["affected_histories"]), (7, 2))

    def test_unequal_batches_pool_counts_energy_and_fraction(self):
        d = self.pooled([results(1, 3, 2, 1.5), results(0, 0, 0, 0.0), results(0, 2, 1, 0.25)], [100, 300, 600])
        self.assertEqual((d["primary_event_capped"], d["secondary_tracks"], d["affected_histories"]), (1, 5, 3))
        self.assertEqual(d["event_cap_ev"], 1.75)
        self.assertAlmostEqual(d["affected_fraction"], 3 / 1000)
        # not the mean of the batch fractions, (0.02 + 0 + 1/600) / 3
        self.assertNotAlmostEqual(d["affected_fraction"], (0.02 + 1 / 600) / 3)

    def test_explicit_zero_is_uncapped(self):
        d = self.pooled([results(), results()], [10, 10])
        self.assertEqual(d["status"], "uncapped")
        self.assertEqual(d["affected_fraction"], 0.0)
        self.assertIn("uncapped", run_electron.truncation_note(d))

    def test_missing_diagnostics_are_unknown_not_zero(self):
        d = self.pooled([{"fates": {}, "budget": {}}], [10])
        self.assertEqual(d["status"], "unknown")
        self.assertTrue(all(d[f] is None for f in run_electron.CAP_FIELDS))
        self.assertIsNone(d["affected_fraction"])
        # an oracle records nothing
        self.assertIsNone(run_electron.Accumulator(1, 1.0, 1).result()["event_cap_diagnostics"])
        self.assertIn("unknown", run_electron.truncation_note(None))

    def test_mixed_old_and_new_batches_stay_unknown(self):
        d = self.pooled([results(), {"fates": {}, "budget": {}}], [10, 10])
        self.assertEqual(d["status"], "unknown")
        self.assertIsNone(d["affected_histories"])
        self.assertIsNone(d["affected_fraction"])
        d = self.pooled([results(0, 1, 1, 2.0), {"fates": {}, "budget": {}}], [10, 10])
        self.assertEqual(d["status"], "truncated")
        self.assertIsNone(d["affected_histories"])


class OracleIdentity(unittest.TestCase):
    """Both electron oracles stop at `git describe`: neither `nebula_cpu_edep`
    nor our `lindhard_microelec` answers `--version`, so the shared base's
    last fallback (`validation/lib/oracle.py`) must stay overridden here or a
    usage error would be recorded as the oracle's version."""

    def chain(self, o, env):
        ran = []

        def fake(cmd, **kwargs):
            ran.append(cmd)
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="usage: ...\n")

        with mock.patch.dict(oracle.os.environ, env, clear=True):
            with mock.patch.object(oracle.subprocess, "run", fake):
                return o.version(Path("/bin/oracle")), ran

    def test_the_binary_is_never_asked_for_its_version(self):
        for o in (run_electron.Nebula(), run_electron.Geant4MicroElec()):
            with self.subTest(oracle=o.name):
                v, ran = self.chain(o, {f"{o.prefix}_SRC": "/src"})
                self.assertEqual(v, "unknown")
                self.assertEqual([c[0] for c in ran], ["git"])

    def test_the_documented_variables_are_the_ones_read(self):
        for o, prefix in ((run_electron.Nebula(), "NEBULA"),
                          (run_electron.Geant4MicroElec(), "GEANT4")):
            with self.subTest(oracle=o.name):
                env = {f"{prefix}_VERSION": "v1", f"{prefix}_COMMIT": "c1", f"{prefix}_BUILD": "b1"}
                with mock.patch.dict(oracle.os.environ, env, clear=True):
                    self.assertEqual(o.version(Path("/bin/oracle")), "v1")
                    self.assertEqual(o.commit(), "c1")
                    self.assertEqual(o.build(), "b1")

    def test_geant4_is_silicon_only_and_nebula_is_not_restricted(self):
        si = {"geant4_material": "G4_Si"}
        self.assertIsNone(run_electron.Geant4MicroElec().applies(si))
        self.assertIn("silicon only", run_electron.Geant4MicroElec().applies({}))
        self.assertIsNone(run_electron.Nebula().applies({}))


if __name__ == "__main__":
    unittest.main()
