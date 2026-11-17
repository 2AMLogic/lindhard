#!/usr/bin/env python3
"""Static workflow checks and reporting-path tests for the advisories workflow.

Needs PyYAML (like test_clippy_canary.py); otherwise standard library only.
Run: python3 .github/scripts/test_advisories_report.py
"""
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import advisories_report as R  # noqa: E402

WF = os.path.join(HERE, "..", "workflows")


def load(name):
    with open(os.path.join(WF, name)) as f:
        text = f.read()
    return text, yaml.safe_load(text)


def triggers(doc):
    return doc.get("on", doc.get(True))  # PyYAML parses bare `on` as True


OWN_UNIT_TEST_RUN = "python3 .github/scripts/test_advisories_report.py"

# Every spelling of a PR-time RustSec check that must never appear in
# `ci.yml` (CONTRIBUTING.md, "Dependency advisories"): the `cargo deny
# check advisories` family keys on "advisories", while `cargo audit`,
# the hyphenated `cargo-audit` binary and the `rustsec/audit-check`
# action contain no such word (issues #880, #882). Common spellings
# only; this list is final (CONTRIBUTING.md, "Dependency advisories").
BANNED_TOKENS = ("advisories", "cargo audit", "cargo-audit", "rustsec/",
                 "audit-check")


def no_advisories_check(doc):
    """True when a parsed `ci.yml` runs no PR-time RustSec (advisories) check.

    The exact `run` of this suite's own unit-test step in workflow-checks
    (issue #871) is exempted first; every remaining step, and all job and
    workflow level `env`, must not mention any of BANNED_TOKENS at all,
    so env-var or flag-splitting indirection cannot hide the check
    (issue #876) and the `cargo audit` / `cargo-audit` (issue #882) /
    `rustsec/audit-check` spellings cannot pass without the word
    "advisories" (issue #880). Operates on a
    deep copy; the caller's doc is left untouched.
    """
    doc = json.loads(json.dumps(doc))
    for job in doc.get("jobs", {}).values():
        job["steps"] = [s for s in job.get("steps", [])
                        if s.get("run", "").strip() != OWN_UNIT_TEST_RUN]
    text = json.dumps(doc)
    return not any(tok in text for tok in BANNED_TOKENS)


class StaticWorkflowTests(unittest.TestCase):
    def test_triggers_are_schedule_and_dispatch_only(self):
        _, doc = load("advisories.yml")
        t = triggers(doc)
        self.assertEqual(set(t), {"schedule", "workflow_dispatch"})
        inp = t["workflow_dispatch"]["inputs"]["simulate-failure"]
        self.assertEqual(inp["type"], "boolean")
        self.assertIs(inp["default"], False)

    def test_ci_runs_no_advisories_check(self):
        # ci.yml runs only `cargo deny check licenses sources`: a PR-time
        # RustSec check would fail unrelated PRs on a newly published
        # advisory (CONTRIBUTING.md, "Dependency advisories"). The exact
        # `run` of this suite's own unit-test step in workflow-checks
        # (issue #871) is exempted; every other step, and any job or
        # workflow level `env`, must not mention any BANNED_TOKENS at all,
        # so env-var or flag-splitting indirection cannot hide the check
        # (issue #876) and `cargo audit` / `cargo-audit` (issue #882) /
        # `rustsec/audit-check`, which lack the word "advisories", cannot
        # pass either (issue #880).
        _, doc = load("ci.yml")
        self.assertTrue(no_advisories_check(doc))

    def test_no_advisories_check_rejects_indirection(self):
        # The predicate must catch an advisories check hidden outside the
        # `run` string, with `cargo deny check` split up, or chained onto
        # the exempted own command (issue #876).
        _, doc = load("ci.yml")

        def injected(mutate):
            d = json.loads(json.dumps(doc))
            mutate(d)
            return no_advisories_check(d)

        def env_indirected(d):
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "deny", "env": {"DENY_ARGS": "advisories"},
                 "run": "cargo deny check $DENY_ARGS --show-stats"})

        def split_flags(d):
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "deny",
                 "run": "cargo deny --manifest-path Cargo.toml check advisories"})

        def own_chained_with_check(d):
            d["jobs"]["workflow-checks"]["steps"].append(
                {"name": "chained",
                 "run": OWN_UNIT_TEST_RUN + " && cargo deny check advisories"})

        def job_env_indirected(d):
            d["jobs"]["dependency-policy"]["env"]["DENY_ARGS"] = "advisories"
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "deny", "run": "cargo deny check $DENY_ARGS"})

        def uses_audit_action(d):
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "audit", "uses": "rustsec/audit-check@v2"})

        def cargo_audit_run(d):
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "audit", "run": "cargo audit --deny warnings"})

        def cargo_audit_dashed_run(d):
            # The hyphenated binary spelling (issue #882). The bare token
            # match also flags e.g. `cargo install cargo-audit`; that is
            # acceptable — the guard is defence-in-depth, not a sandbox
            # (CONTRIBUTING.md, "Dependency advisories").
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "audit", "run": "cargo-audit audit --deny warnings"})

        def audit_check_ref(d):
            d["jobs"]["dependency-policy"]["steps"].append(
                {"name": "audit", "uses": "acme/audit-check@v3"})

        for mutate in (env_indirected, split_flags, own_chained_with_check,
                       job_env_indirected, uses_audit_action, cargo_audit_run,
                       cargo_audit_dashed_run, audit_check_ref):
            self.assertFalse(injected(mutate), mutate.__name__)

    def test_no_advisories_check_accepts_own_step_anywhere(self):
        # The exemption is the exact `run`, not the step's name or job: a
        # copy of the unit-test step in another job stays allowed.
        _, doc = load("ci.yml")
        d = json.loads(json.dumps(doc))
        own = [s for s in d["jobs"]["workflow-checks"]["steps"]
               if s.get("run", "").strip() == OWN_UNIT_TEST_RUN]
        d["jobs"]["rust"]["steps"].extend(own)
        self.assertTrue(no_advisories_check(d))

    def test_permissions(self):
        _, doc = load("advisories.yml")
        self.assertEqual(doc["permissions"], {"contents": "read"})
        self.assertNotIn("permissions", doc["jobs"]["advisories"])
        self.assertEqual(doc["jobs"]["report"]["permissions"],
                         {"contents": "read", "issues": "write"})

    def test_check_command_and_pin(self):
        _, doc = load("advisories.yml")
        job = doc["jobs"]["advisories"]
        runs = [s.get("run", "") for s in job["steps"]]
        self.assertIn("cargo deny check advisories --show-stats", runs)
        self.assertRegex(job["env"]["CARGO_DENY_VERSION"], r"^\d+\.\d+\.\d+$")
        self.assertRegex(job["env"]["CARGO_DENY_SHA256"], r"^[0-9a-f]{64}$")

    def test_result_recorded_even_on_failure(self):
        _, doc = load("advisories.yml")
        rec = [s for s in doc["jobs"]["advisories"]["steps"]
               if "result/advisories.txt" in s.get("run", "")]
        self.assertEqual(len(rec), 1)
        self.assertEqual(rec[0]["if"], "always()")
        self.assertEqual(doc["jobs"]["report"]["if"], "always()")


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
OK, BAD = {"advisories": "success"}, {"advisories": "failure"}


class ReportingTests(unittest.TestCase):
    def test_failure_no_marker_creates_exactly_one_issue(self):
        api = FakeApi([{"number": 5, "state": "open", "body": "unrelated"}])
        code, _ = R.report(api, BAD, SHA, URL)
        self.assertEqual(code, 1)
        self.assertEqual(api.writes(), [("POST", "/issues")])
        body = api.issues[-1]["body"]
        self.assertTrue(R.has_marker(body))
        self.assertIn(SHA, body)
        self.assertIn(URL, body)

    def test_failure_open_issue_comments_without_creating(self):
        api = FakeApi([marked(3, "open")])
        R.report(api, BAD, SHA, URL)
        self.assertEqual(api.writes(), [("POST", "/issues/3/comments")])
        self.assertIn(SHA, api.comments[0][1])
        self.assertIn(URL, api.comments[0][1])

    def test_failure_closed_issue_reopens_then_comments(self):
        api = FakeApi([marked(3, "closed")])
        R.report(api, BAD, SHA, URL)
        self.assertEqual(api.writes(),
                         [("PATCH", "/issues/3"), ("POST", "/issues/3/comments")])
        self.assertEqual(api.issues[0]["state"], "open")

    def test_green_open_issue_recovers_and_closes(self):
        api = FakeApi([marked(3, "open")])
        code, _ = R.report(api, OK, SHA, URL)
        self.assertEqual(code, 0)
        self.assertEqual(api.writes(),
                         [("POST", "/issues/3/comments"), ("PATCH", "/issues/3")])
        self.assertIn("recovered", api.comments[0][1])
        self.assertEqual(api.issues[0]["state"], "closed")

    def test_green_no_issue_or_closed_issue_is_noop(self):
        for issues in ([], [marked(3, "closed")]):
            api = FakeApi(issues)
            code, _ = R.report(api, OK, SHA, URL)
            self.assertEqual(code, 0)
            self.assertEqual(api.writes(), [])

    def test_two_marker_matches_error_without_writes(self):
        for results in (BAD, OK):
            api = FakeApi([marked(3, "open"), marked(4, "closed")])
            with self.assertRaises(R.ReportError) as cm:
                R.report(api, results, SHA, URL)
            self.assertIn("#3", str(cm.exception))
            self.assertEqual(api.writes(), [])

    def test_marker_quoted_mid_body_is_not_canonical(self):
        quoting = {"number": 336, "state": "open",
                   "body": f"Track failures with `{R.MARKER}` in the body.\n"}
        api = FakeApi([quoting])
        R.report(api, BAD, SHA, URL)
        self.assertEqual(api.writes(), [("POST", "/issues")])

    def test_pull_request_with_marker_is_ignored(self):
        api = FakeApi([{**marked(7, "open"), "pull_request": {}}])
        R.report(api, BAD, SHA, URL)
        self.assertEqual(api.writes(), [("POST", "/issues")])

    def test_missing_result_counts_as_failure(self):
        api = FakeApi([marked(3, "open")])
        code, _ = R.report(api, {}, SHA, URL)
        self.assertEqual(code, 1)
        self.assertEqual(len(api.comments), 1)
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(R.read_results(d), {})

    def test_api_error_raises(self):
        api = FakeApi([], fail_on="POST /issues")
        with self.assertRaises(R.ReportError):
            R.report(api, BAD, SHA, URL)
        api = FakeApi(fail_on="GET")
        with self.assertRaises(R.ReportError):
            R.report(api, BAD, SHA, URL)
        self.assertEqual(api.writes(), [])

    def test_main_exit_codes(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "advisories.txt"), "w") as f:
                f.write("failure\n")
            env = {"GITHUB_REPOSITORY": "o/r", "GITHUB_TOKEN": "t",
                   "RESULTS_DIR": d, "SHA": SHA, "RUN_URL": URL}
            with mock.patch.dict(os.environ, env), \
                 mock.patch.object(R, "GitHubApi", lambda *a: FakeApi([marked(3, "open")])):
                self.assertEqual(R.main(), 1)  # check failed, reporting worked
            with mock.patch.dict(os.environ, env), \
                 mock.patch.object(R, "GitHubApi", lambda *a: FakeApi(fail_on="GET")):
                self.assertEqual(R.main(), 2)  # reporting failed (e.g. no issues: write)


if __name__ == "__main__":
    unittest.main()
