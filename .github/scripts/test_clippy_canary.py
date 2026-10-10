#!/usr/bin/env python3
"""Static workflow checks and reporting-path tests for the clippy pin/canary.

Needs PyYAML (CI installs it); otherwise standard library only.
Run: python3 .github/scripts/test_clippy_canary.py
"""
import os
import re
import sys
import unittest

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import clippy_canary_report as R  # noqa: E402

WF = os.path.join(HERE, "..", "workflows")


def load(name):
    with open(os.path.join(WF, name)) as f:
        text = f.read()
    return text, yaml.safe_load(text)


def triggers(doc):
    return doc.get("on", doc.get(True))  # PyYAML parses bare `on` as True


class StaticWorkflowTests(unittest.TestCase):
    def test_ci_clippy_uses_exact_pin_and_fmt_test_stay_floating(self):
        text, doc = load("ci.yml")
        job = doc["jobs"]["rust"]
        pin = job["env"]["CLIPPY_TOOLCHAIN"]
        self.assertRegex(pin, r"^\d+\.\d+\.\d+$", "pin must be a full release")
        steps = {s.get("name"): s for s in job["steps"]}
        clippy = [s for s in job["steps"] if "clippy" in s.get("run", "")
                  and "cargo" in s.get("run", "")]
        self.assertEqual(len(clippy), 1)
        self.assertIn('cargo +"$CLIPPY_TOOLCHAIN" clippy --workspace --all-targets -- -D warnings',
                      clippy[0]["run"])
        installs = [s for s in job["steps"] if "rustup toolchain install" in s.get("run", "")]
        self.assertEqual(len(installs), 1)
        self.assertIn("$CLIPPY_TOOLCHAIN", installs[0]["run"])
        self.assertIn("--component clippy", installs[0]["run"])
        stable = [s for s in job["steps"] if s.get("uses", "").endswith("rust-toolchain@stable")]
        self.assertEqual(len(stable), 1)
        self.assertEqual(stable[0]["with"]["components"], "rustfmt")
        runs = [s.get("run", "") for s in job["steps"]]
        self.assertIn("cargo fmt --all --check", runs)
        self.assertIn("cargo test --workspace", runs)
        for r in ("cargo fmt --all --check", "cargo test --workspace"):
            self.assertNotIn("+", r)
        # nothing else in ci.yml runs clippy
        self.assertEqual(len(re.findall(r"cargo[^\n]*clippy", text)), 1)

    def test_ci_does_not_reference_canary(self):
        text, _ = load("ci.yml")
        self.assertNotIn("clippy-canary.yml", text.replace("# Latest-stable and beta drift is reported by clippy-canary.yml.", ""))
        self.assertNotIn("workflow_call", text)

    def test_canary_triggers(self):
        _, doc = load("clippy-canary.yml")
        t = triggers(doc)
        self.assertIn("schedule", t)
        self.assertIn("workflow_dispatch", t)
        self.assertNotIn("pull_request", t)
        self.assertNotIn("pull_request_target", t)
        self.assertNotIn("push", t)
        inp = t["workflow_dispatch"]["inputs"]["simulate-failure"]
        self.assertEqual(inp["type"], "boolean")
        self.assertIs(inp["default"], False)

    def test_canary_matrix_and_command(self):
        text, doc = load("clippy-canary.yml")
        job = doc["jobs"]["clippy"]
        self.assertEqual(job["strategy"]["matrix"]["channel"], ["stable", "beta"])
        runs = [s.get("run", "") for s in job["steps"]]
        self.assertIn("cargo clippy --workspace --all-targets -- -D warnings", runs)
        sim = [s for s in job["steps"] if "simulate-failure" in str(s.get("if", ""))]
        self.assertEqual(len(sim), 1)
        self.assertIn("workflow_dispatch", sim[0]["if"])

    def test_canary_permissions(self):
        _, doc = load("clippy-canary.yml")
        self.assertEqual(doc["permissions"], {"contents": "read"})
        self.assertNotIn("permissions", doc["jobs"]["clippy"])
        self.assertEqual(doc["jobs"]["report"]["permissions"],
                         {"contents": "read", "issues": "write"})


class FakeApi:
    """In-memory issues API; records calls. `fail_on` is a method/path prefix."""

    def __init__(self, issues=(), fail_on=None):
        self.issues = [dict(i) for i in issues]
        self.calls, self.comments, self.fail_on = [], [], fail_on

    def request(self, method, path, body=None):
        self.calls.append((method, path))
        if self.fail_on and f"{method} {path}".startswith(self.fail_on):
            raise R.ReportError(f"{method} {path} failed: HTTP 500")
        if method == "GET":
            return self.issues if "page=1" in path else []
        if method == "POST" and path == "/issues":
            self.issues.append({"number": 99, "state": "open", **body})
            return {}
        n = int(path.split("/")[2])
        issue = next(i for i in self.issues if i["number"] == n)
        if path.endswith("/comments"):
            self.comments.append((n, body["body"]))
        else:
            issue["state"] = body["state"]
        return {}

    def writes(self):
        return [c for c in self.calls if c[0] != "GET"]


def marked(n, state):
    return {"number": n, "state": state, "body": f"{R.MARKER}\nx"}


SHA, URL = "abc123", "https://example.invalid/run/1"
OK, BAD = "success", "failure"


class ReportingTests(unittest.TestCase):
    def run_report(self, api, stable=OK, beta=OK):
        return R.report(api, {"stable": stable, "beta": beta}, SHA, URL)

    def test_failure_no_marker_creates_one_issue(self):
        api = FakeApi([{"number": 5, "state": "open", "body": "unrelated"}])
        code, _ = self.run_report(api, stable=BAD)
        self.assertEqual(code, 1)
        self.assertEqual(api.writes(), [("POST", "/issues")])
        body = api.issues[-1]["body"]
        self.assertIn(R.MARKER, body)
        self.assertIn("stable", body)
        self.assertIn(SHA, body)
        self.assertIn(URL, body)

    def test_marker_quoted_mid_body_is_not_canonical(self):
        quoting = {"number": 336, "state": "open",
                   "body": f"Track failures with `{R.MARKER}` in the body.\n"}
        api = FakeApi([quoting])
        self.assertEqual(R.find_tracking_issues(api), [])
        self.run_report(api, stable=BAD)
        self.assertEqual(api.writes(), [("POST", "/issues")])
        self.assertEqual(api.comments, [])
        self.assertEqual(api.issues[0]["state"], "open")

    def test_marker_quoted_mid_body_not_closed_on_green(self):
        quoting = {"number": 336, "state": "open",
                   "body": f"Example:\n\n{R.MARKER}\n"}
        api = FakeApi([quoting])
        self.run_report(api)
        self.assertEqual(api.writes(), [])

    def test_marker_at_start_is_canonical(self):
        anchored = {"number": 4, "state": "open",
                    "body": f"﻿ \n{R.MARKER}\n\nbody text"}
        api = FakeApi([anchored, marked(8, "closed")])
        nums = [i["number"] for i in R.find_tracking_issues(api)]
        self.assertEqual(nums, [4, 8])

    def test_created_issue_body_is_recognised(self):
        api = FakeApi()
        self.run_report(api, beta=BAD)
        self.assertTrue(R.has_marker(api.issues[-1]["body"]))

    def test_pull_request_with_marker_is_ignored(self):
        pr = {**marked(7, "open"), "pull_request": {}}
        api = FakeApi([pr])
        self.run_report(api, beta=BAD)
        self.assertEqual(api.writes(), [("POST", "/issues")])

    def test_failure_one_open_marker_comments_only(self):
        api = FakeApi([marked(3, "open")])
        self.run_report(api, beta=BAD)
        self.assertEqual(api.writes(), [("POST", "/issues/3/comments")])
        text = api.comments[0][1]
        self.assertIn("beta", text)
        self.assertNotIn("stable", text)
        self.assertIn(SHA, text)
        self.assertIn(URL, text)

    def test_failure_closed_marker_reopens_then_comments(self):
        api = FakeApi([marked(3, "closed")])
        self.run_report(api, stable=BAD)
        self.assertEqual(api.writes(),
                         [("PATCH", "/issues/3"), ("POST", "/issues/3/comments")])
        self.assertEqual(api.issues[0]["state"], "open")

    def test_both_failures_listed_in_one_comment(self):
        api = FakeApi([marked(3, "open")])
        self.run_report(api, stable=BAD, beta=BAD)
        self.assertEqual(len(api.comments), 1)
        self.assertIn("stable, beta", api.comments[0][1])

    def test_missing_result_counts_as_failure(self):
        api = FakeApi([marked(3, "open")])
        code, _ = R.report(api, {"stable": OK}, SHA, URL)
        self.assertEqual(code, 1)
        self.assertIn("beta", api.comments[0][1])

    def test_multiple_markers_fail_without_writes(self):
        for kwargs in ({"stable": BAD}, {}):
            api = FakeApi([marked(3, "open"), marked(4, "closed")])
            with self.assertRaises(R.ReportError) as cm:
                self.run_report(api, **kwargs)
            self.assertIn("#3", str(cm.exception))
            self.assertEqual(api.writes(), [])

    def test_green_open_issue_recovers_and_closes(self):
        api = FakeApi([marked(3, "open")])
        code, _ = self.run_report(api)
        self.assertEqual(code, 0)
        self.assertEqual(api.writes(),
                         [("POST", "/issues/3/comments"), ("PATCH", "/issues/3")])
        self.assertIn("recovered", api.comments[0][1])
        self.assertEqual(api.issues[0]["state"], "closed")

    def test_green_no_issue_or_closed_issue_is_noop(self):
        for issues in ([], [marked(3, "closed")]):
            api = FakeApi(issues)
            code, _ = self.run_report(api)
            self.assertEqual(code, 0)
            self.assertEqual(api.writes(), [])

    def test_api_failure_raises_and_creates_nothing_else(self):
        api = FakeApi([], fail_on="POST /issues")
        with self.assertRaises(R.ReportError):
            self.run_report(api, stable=BAD)
        self.assertEqual(api.writes(), [("POST", "/issues")])
        api = FakeApi(fail_on="GET")
        with self.assertRaises(R.ReportError):
            self.run_report(api, stable=BAD)
        self.assertEqual(api.writes(), [])

    def test_main_exit_codes(self):
        import tempfile
        from unittest import mock
        with tempfile.TemporaryDirectory() as d:
            for c, v in (("stable", OK), ("beta", BAD)):
                with open(os.path.join(d, f"{c}.txt"), "w") as f:
                    f.write(v + "\n")
            env = {"GITHUB_REPOSITORY": "o/r", "GITHUB_TOKEN": "t",
                   "RESULTS_DIR": d, "SHA": SHA, "RUN_URL": URL}
            with mock.patch.dict(os.environ, env), \
                 mock.patch.object(R, "GitHubApi", lambda *a: FakeApi([marked(3, "open")])):
                self.assertEqual(R.main(), 1)  # lint failed, reporting worked
            with mock.patch.dict(os.environ, env), \
                 mock.patch.object(R, "GitHubApi", lambda *a: FakeApi(fail_on="GET")):
                self.assertEqual(R.main(), 2)  # reporting failed


if __name__ == "__main__":
    unittest.main()
