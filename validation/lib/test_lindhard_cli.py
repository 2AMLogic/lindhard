#!/usr/bin/env python3
"""Tests of the completion-manifest resume contract in lindhard_cli.run
(validation/README.md, "Resuming a sweep"). Temporary files and a mocked
subprocess only; standard library only:

    python3 validation/lib/test_lindhard_cli.py
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lindhard_cli as lc  # noqa: E402

SUMMARY = {
    "software": {"git_describe": "v0-dirty"},
    "results": {
        "histories": 10,
        "recoils": {"sputtered": 4},
        "primaries": {"stopped_depth_mean_nm": 1.0, "stopped_depth_std_nm": 0.1},
        "yields": {"backscattered_per_ion": 0.1, "sputtered_per_ion": 0.4},
    },
    "physics": {"target": [{"material": {"elements": []}}]},
    "run": {"ions_per_s": 1.0},
}


def problem(tables=()):
    return {
        "id": "p",
        "beam": {"ion": "Ar"},
        "target": {"substrate": "Si"},
        "physics": {},
        "run": {"ions": 10},
        "stopping_tables": list(tables),
    }


class Fixture(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = Path(self._tmp.name)
        self.exe = self.tmp / "lindhard"
        self.exe.write_bytes(b"build-1")
        self.work = self.tmp / "w"
        self.calls = 0
        self.fail = False
        self.summary = SUMMARY

    def fake(self, argv, **kw):
        self.calls += 1
        if self.fail:
            return mock.Mock(returncode=1, stderr="boom")
        out = Path(argv[argv.index("--out") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "summary.json").write_text(json.dumps(self.summary))
        return mock.Mock(returncode=0, stderr="")

    def go(self, p=None, exe_bytes=None):
        if exe_bytes is not None:
            self.exe.write_bytes(exe_bytes)
        bid = lc.binary_sha256(self.exe)
        with mock.patch.object(lc.subprocess, "run", self.fake):
            return lc.run(p or problem(), self.work, self.exe, binary_id=bid)


class ResumeContract(Fixture):
    def test_identical_reuses_without_launching(self):
        self.go()
        self.go()
        self.assertEqual(self.calls, 1)

    def test_same_describe_different_executable_reruns(self):
        self.go()
        self.go(exe_bytes=b"build-2")  # summary git_describe is unchanged
        self.assertEqual(self.calls, 2)

    def test_changed_input_reruns(self):
        self.go()
        p = problem()
        p["run"]["ions"] = 11
        self.go(p)
        self.assertEqual(self.calls, 2)

    def test_changed_stopping_table_reruns(self):
        t = self.tmp / "t.toml"
        t.write_text("a")
        p = problem([t])
        self.go(p)
        self.go(p)
        self.assertEqual(self.calls, 1)
        t.write_text("b")  # same path, new bytes
        self.go(p)
        self.assertEqual(self.calls, 2)

    def test_missing_stopping_table_is_an_error(self):
        with self.assertRaises(SystemExit), mock.patch.object(lc.subprocess, "run", self.fake):
            lc.run(problem([self.tmp / "nope.toml"]), self.work, self.exe, binary_id="x")
        self.assertEqual(self.calls, 0)

    def test_unreadable_binary_is_an_error(self):
        with self.assertRaises(SystemExit):
            lc.binary_sha256(self.tmp / "absent")

    def test_missing_malformed_obsolete_manifest_is_a_miss(self):
        self.go()
        m = self.work / lc.MANIFEST_NAME
        good = json.loads(m.read_text())
        n = 1
        for bad in (None, "{not json", "[]", json.dumps({**good, "manifest_version": 0}),
                    json.dumps({**good, "complete": False}),
                    json.dumps({k: v for k, v in good.items() if k != "binary_sha256"})):
            if bad is None:
                m.unlink()
            else:
                m.write_text(bad)
            self.go()
            n += 1
            self.assertEqual(self.calls, n)

    def test_incomplete_or_altered_summary_is_a_miss(self):
        self.go()
        s = self.work / "out" / "summary.json"
        s.write_text(json.dumps({"software": SUMMARY["software"]}))
        self.go()
        self.assertEqual(self.calls, 2)
        s.write_text(s.read_text() + " ")  # bytes differ from the manifest
        self.go()
        self.assertEqual(self.calls, 3)

    def test_failed_rerun_leaves_no_manifest(self):
        self.go()
        self.fail = True
        with self.assertRaises(SystemExit):
            self.go(exe_bytes=b"build-2")
        self.assertFalse((self.work / lc.MANIFEST_NAME).exists())
        self.fail = False
        self.go(exe_bytes=b"build-1")  # old outputs must not be certified
        self.assertEqual(self.calls, 3)

    def test_unreadable_summary_after_run_publishes_nothing(self):
        self.summary = {"software": {}}
        with self.assertRaises(SystemExit):
            self.go()
        self.assertFalse((self.work / lc.MANIFEST_NAME).exists())

    def test_without_binary_id_always_runs(self):
        with mock.patch.object(lc.subprocess, "run", self.fake):
            lc.run(problem(), self.work, self.exe)
            lc.run(problem(), self.work, self.exe)
        self.assertEqual(self.calls, 2)
        self.assertFalse((self.work / lc.MANIFEST_NAME).exists())


if __name__ == "__main__":
    unittest.main()
