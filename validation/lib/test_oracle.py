#!/usr/bin/env python3
"""Tests of the oracle-harness scaffolding the two level-2 runners share
(#778): the one environment-variable convention, the `version()` fallback
chain, and the float rounding that writes the committed summaries. Mocked
subprocess and temporary files only; standard library only:

    python3 validation/lib/test_oracle.py
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oracle  # noqa: E402


class Derived(oracle.Oracle):
    name = "Derived"
    env = "DERIVED_BIN"


class Explicit(oracle.Oracle):
    name = "Explicit"
    env = "EXPLICIT_TOOL_BIN"
    prefix = "EXPLICIT"


def completed(code=0, out="", err=""):
    return subprocess.CompletedProcess([], code, stdout=out, stderr=err)


class Variables(unittest.TestCase):
    """One convention: `<PREFIX>_VERSION` etc., `<PREFIX>` from `env` unless
    the oracle says otherwise."""

    def resolved(self, o, suffix):
        """The variable name `_var` reads, recovered from the environment it
        is given."""
        with mock.patch.dict(oracle.os.environ, {}, clear=True):
            oracle.os.environ["SENTINEL"] = "x"
            for name in ("DERIVED_VERSION", "DERIVED_SRC", "EXPLICIT_VERSION", "EXPLICIT_SRC",
                         "EXPLICIT_TOOL_VERSION", "EXPLICIT_TOOL_SRC"):
                oracle.os.environ[name] = name
            seen = o._var(suffix)
        return seen

    def test_the_prefix_defaults_to_env_without_its_bin_suffix(self):
        self.assertEqual(self.resolved(Derived(), "_VERSION"), "DERIVED_VERSION")
        self.assertEqual(self.resolved(Derived(), "_SRC"), "DERIVED_SRC")

    def test_an_explicit_prefix_wins_over_the_derived_one(self):
        # The Geant4 case: GEANT4_MICROELEC_BIN but GEANT4_VERSION.
        self.assertEqual(self.resolved(Explicit(), "_VERSION"), "EXPLICIT_VERSION")
        self.assertEqual(self.resolved(Explicit(), "_SRC"), "EXPLICIT_SRC")

    def test_an_unset_variable_is_none(self):
        with mock.patch.dict(oracle.os.environ, {}, clear=True):
            self.assertIsNone(Derived()._var("_BUILD"))


class Binary(unittest.TestCase):
    def test_an_unconfigured_oracle_is_not_configured_not_an_error(self):
        with mock.patch.dict(oracle.os.environ, {}, clear=True):
            with self.assertRaises(oracle.NotConfigured):
                Derived().binary()

    def test_a_configured_path_that_does_not_exist_is_not_configured(self):
        with mock.patch.dict(oracle.os.environ, {"DERIVED_BIN": "/no/such/binary"}, clear=True):
            with self.assertRaises(oracle.NotConfigured):
                Derived().binary()

    def test_check_runs_on_a_resolved_binary(self):
        seen = []

        class Checked(Derived):
            def check(self, binary):
                seen.append(binary)

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "tool"
            p.write_text("")
            with mock.patch.dict(oracle.os.environ, {"DERIVED_BIN": str(p)}, clear=True):
                self.assertEqual(Checked().binary(), p)
        self.assertEqual(seen, [p])


class Version(unittest.TestCase):
    """`$<PREFIX>_VERSION` -> `git describe` of `$<PREFIX>_SRC` ->
    `<binary> --version` -> "unknown"."""

    def version(self, env, runner):
        with mock.patch.dict(oracle.os.environ, env, clear=True):
            with mock.patch.object(oracle.subprocess, "run", runner):
                return Derived().version(Path("/bin/tool"))

    def test_the_explicit_variable_wins_and_runs_nothing(self):
        def never(*a, **k):
            raise AssertionError(f"ran {a[0]}")

        self.assertEqual(self.version({"DERIVED_VERSION": "1.2.3", "DERIVED_SRC": "/src"}, never), "1.2.3")

    def test_git_describe_of_src_comes_before_the_binary(self):
        cmds = []

        def fake(cmd, **k):
            cmds.append(cmd)
            return completed(out="v0.9-3-gabc\n")

        self.assertEqual(self.version({"DERIVED_SRC": "/src"}, fake), "v0.9-3-gabc")
        self.assertEqual(cmds, [["git", "-C", "/src", "describe", "--tags", "--always"]])

    def test_the_binary_is_asked_last_and_only_its_first_line_is_kept(self):
        def fake(cmd, **k):
            if cmd[0] == "git":
                return completed(1)
            return completed(out="tool 2.0\nbuilt with gcc\n")

        self.assertEqual(self.version({"DERIVED_SRC": "/src"}, fake), "tool 2.0")

    def test_the_binary_may_answer_on_stderr(self):
        self.assertEqual(self.version({}, lambda cmd, **k: completed(out="", err="tool 2.0\n")), "tool 2.0")

    def test_everything_silent_is_unknown(self):
        self.assertEqual(self.version({}, lambda cmd, **k: completed(out="", err="")), "unknown")

    def test_a_binary_that_cannot_be_run_is_unknown_not_a_crash(self):
        def boom(cmd, **k):
            raise OSError("no such file")

        self.assertEqual(self.version({}, boom), "unknown")

    def test_a_binary_that_hangs_is_unknown_not_a_crash(self):
        def slow(cmd, **k):
            raise subprocess.TimeoutExpired(cmd, 30)

        self.assertEqual(self.version({}, slow), "unknown")


class Commit(unittest.TestCase):
    def test_the_explicit_variable_wins(self):
        with mock.patch.dict(oracle.os.environ, {"DERIVED_COMMIT": "deadbeef"}, clear=True):
            self.assertEqual(Derived().commit(), "deadbeef")

    def test_without_src_it_is_unknown(self):
        with mock.patch.dict(oracle.os.environ, {}, clear=True):
            self.assertEqual(Derived().commit(), "unknown")
            self.assertEqual(Derived().build(), "unknown")


class GitOutput(unittest.TestCase):
    def test_no_source_directory_is_none(self):
        self.assertIsNone(oracle.git_output(None, "rev-parse", "HEAD"))
        self.assertIsNone(oracle.git_output("", "rev-parse", "HEAD"))

    def test_a_directory_that_is_not_a_repository_is_none(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(oracle.git_output(d, "rev-parse", "HEAD"))

    def test_empty_output_is_none_not_an_empty_string(self):
        with mock.patch.object(oracle.subprocess, "run", lambda cmd, **k: completed(out="\n")):
            self.assertIsNone(oracle.git_output("/src", "describe"))

    def test_a_missing_or_slow_git_is_none(self):
        for exc in (OSError("no git"), subprocess.TimeoutExpired(["git"], 30)):
            def boom(cmd, _e=exc, **k):
                raise _e

            with mock.patch.object(oracle.subprocess, "run", boom):
                self.assertIsNone(oracle.git_output("/src", "describe"))


class RoundFloats(unittest.TestCase):
    """The committed summaries are rounded to six significant digits so a
    rerun that only changes the last bits leaves them byte-identical."""

    def test_floats_keep_six_significant_digits(self):
        self.assertEqual(oracle.round_floats(1.2345678901), 1.23457)
        self.assertEqual(oracle.round_floats(1.2345678901e-9), 1.23457e-9)

    def test_it_recurses_through_dicts_and_lists(self):
        got = oracle.round_floats({"a": [1.0000004999, {"b": (0.1 + 0.2)}]})
        self.assertEqual(got, {"a": [1.0, {"b": 0.3}]})

    def test_non_floats_are_returned_unchanged(self):
        for x in (7, True, False, None, "1.23456789"):
            self.assertIs(oracle.round_floats(x), x)
        self.assertEqual(oracle.round_floats({"xs": [], "d": {}, "s": "1.2345678"}),
                         {"xs": [], "d": {}, "s": "1.2345678"})

    def test_an_integer_does_not_become_a_float(self):
        self.assertIsInstance(oracle.round_floats({"n": 1000000})["n"], int)

    def test_the_committed_summaries_are_already_rounded(self):
        """Round-tripping each committed level-2 summary through the shared
        rounding reproduces it byte for byte."""
        import json

        summaries = Path(__file__).resolve().parents[1] / "oracles" / "summaries"
        files = sorted(summaries.glob("*.json"))
        self.assertTrue(files, f"no committed summaries under {summaries}")
        for f in files:
            raw = f.read_bytes()
            again = json.dumps(oracle.round_floats(json.loads(raw)), indent=2, ensure_ascii=False) + "\n"
            self.assertEqual(again.encode(), raw, f"{f.name} is not what round_floats writes")


if __name__ == "__main__":
    unittest.main()
