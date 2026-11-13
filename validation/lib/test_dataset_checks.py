#!/usr/bin/env python3
"""Tests of the shared level-3 dataset checks (validation/data/README.md).

Temporary files only; standard library only:

    python3 validation/lib/test_dataset_checks.py
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dataset_checks as dc  # noqa: E402

REQUIRED = ["id", "citation", "added"]
PROVENANCE = "| `set_a` | a paper, table 1 | measured |\n"


def good(stem="set_a"):
    return {"id": stem, "citation": "A. Author, J. Phys. 1 (1970)", "added": "2026-01-01"}


class IsNumber(unittest.TestCase):
    def test_finite_numbers(self):
        for x in (0, -1, 1.5, 1e300):
            self.assertTrue(dc.is_number(x), x)

    def test_not_numbers(self):
        for x in (True, False, None, "1", float("nan"), float("inf"), float("-inf"), [1]):
            self.assertFalse(dc.is_number(x), x)


class CheckCommon(unittest.TestCase):
    def errs(self, d, stem="set_a", provenance=PROVENANCE):
        return dc.check_common(Path(f"/data/{stem}.json"), d, REQUIRED, provenance)

    def test_good_dataset_passes(self):
        self.assertEqual(self.errs(good()), [])

    def test_missing_required_key(self):
        d = good()
        del d["citation"]
        self.assertIn("`citation` missing or empty", self.errs(d))

    def test_blank_required_key(self):
        d = good()
        d["citation"] = "   "
        self.assertIn("`citation` missing or empty", self.errs(d))

    def test_non_string_required_key(self):
        d = good()
        d["citation"] = 7
        self.assertIn("`citation` missing or empty", self.errs(d))

    def test_id_differs_from_file_name(self):
        d = good()
        d["id"] = "other"
        self.assertIn("`id` 'other' differs from the file name 'set_a'", self.errs(d))

    def test_id_not_in_provenance(self):
        self.assertIn("no row in docs/data-provenance.md names `set_a`", self.errs(good(), provenance=""))

    def test_missing_id_is_reported_once(self):
        """A missing `id` is a missing required key and a name mismatch, but
        there is no id to look for in docs/data-provenance.md."""
        d = good()
        del d["id"]
        errs = self.errs(d)
        self.assertIn("`id` missing or empty", errs)
        self.assertTrue(any("differs from the file name" in e for e in errs), errs)
        self.assertFalse(any("no row in docs/data-provenance.md" in e for e in errs), errs)


class LoadChecked(unittest.TestCase):
    """`load_checked` against a temporary data directory and provenance file."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)  # stands in for the repository root
        self.data = self.repo / "validation" / "data" / "test"
        self.data.mkdir(parents=True)
        self.provenance = self.repo / "docs" / "data-provenance.md"
        self.provenance.parent.mkdir(parents=True)
        self.provenance.write_text(PROVENANCE)
        for name, value in (("REPO", self.repo), ("PROVENANCE", self.provenance)):
            patch = mock.patch.object(dc, name, value)
            patch.start()
            self.addCleanup(patch.stop)

    def write(self, stem, obj):
        f = self.data / f"{stem}.json"
        f.write_text(obj if isinstance(obj, str) else json.dumps(obj))
        return f

    @staticmethod
    def check(path, d, provenance):
        return dc.check_common(path, d, REQUIRED, provenance)

    def load(self, files):
        """`load_checked` with stderr captured: (result-or-exit-code, stderr)."""
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            try:
                return dc.load_checked(files, self.check, "test dataset"), err.getvalue()
            except SystemExit as e:
                return e.code, err.getvalue()

    def test_returns_path_and_dict_in_order(self):
        self.provenance.write_text(PROVENANCE + "| `set_b` | a paper, table 2 | measured |\n")
        files = [self.write("set_a", good()), self.write("set_b", good("set_b"))]
        out, err = self.load(files)
        self.assertEqual(err, "")
        self.assertEqual([p for p, _ in out], files)
        self.assertEqual([d["id"] for _, d in out], ["set_a", "set_b"])

    def test_no_files_is_not_a_failure(self):
        self.assertEqual(self.load([]), ([], ""))

    def test_failing_dataset_exits_2_naming_the_label(self):
        d = good()
        del d["citation"]
        f = self.write("set_a", d)
        code, err = self.load([f])
        self.assertEqual(code, 2)
        self.assertIn("error: validation/data/test/set_a.json: `citation` missing or empty", err)
        self.assertIn("error: test dataset check failed (rules: validation/data/README.md)", err)

    def test_invalid_json_is_reported_and_fails(self):
        f = self.write("set_a", "{not json")
        code, err = self.load([f])
        self.assertEqual(code, 2)
        self.assertIn("error: validation/data/test/set_a.json: not valid JSON:", err)

    def test_one_bad_file_fails_the_whole_load(self):
        """Every file is still read and reported before the exit."""
        self.provenance.write_text(PROVENANCE + "| `set_b` | a paper, table 2 | measured |\n")
        bad = good()
        del bad["added"]
        files = [self.write("set_a", bad), self.write("set_b", good("set_b"))]
        code, err = self.load(files)
        self.assertEqual(code, 2)
        self.assertIn("error: validation/data/test/set_a.json: `added` missing or empty", err)
        self.assertNotIn("set_b.json", err)

    def test_path_outside_the_repo_is_printed_whole(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        f = Path(tmp.name) / "set_a.json"
        d = good()
        del d["citation"]
        f.write_text(json.dumps(d))
        code, err = self.load([f])
        self.assertEqual(code, 2)
        self.assertIn(f"error: {f}: `citation` missing or empty", err)


if __name__ == "__main__":
    unittest.main()
