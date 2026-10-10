#!/usr/bin/env python3
"""Report the clippy canary result through one marker-identified GitHub issue.

Called by .github/workflows/clippy-canary.yml after the stable/beta matrix.
Reads one `<channel>.txt` file per channel from RESULTS_DIR (content: the job
status, `success` or anything else). A missing file counts as a failure.

Lifecycle (standard library only, no third-party actions):
  * any channel failed: find the canonical issue by MARKER; create it only if
    none exists, else reopen it if closed and add one comment (channels, SHA,
    run URL). More than one marker match is an error (exit 1), never a guess.
  * all green: comment "recovered" and close the canonical issue if it is
    open; no canonical issue, or one already closed, is a no-op.
  * a reporting/API failure exits non-zero, so the scheduled run stays red and
    a lint failure is never converted into success.
A failed channel makes the script exit 1 even when reporting succeeded, so the
run itself is red whenever the lint is.
"""
import json
import os
import sys
import urllib.error
import urllib.request

MARKER = "<!-- clippy-canary-tracking -->"
TITLE = "CI: clippy canary failing on latest stable or beta Rust"
CHANNELS = ("stable", "beta")


class ReportError(Exception):
    pass


class GitHubApi:
    """Minimal REST client. Any HTTP/network error raises ReportError."""

    def __init__(self, repo, token, base="https://api.github.com"):
        self.repo, self.token, self.base = repo, token, base

    def request(self, method, path, body=None):
        url = f"{self.base}/repos/{self.repo}{path}"
        data = None if body is None else json.dumps(body).encode()
        req = urllib.request.Request(url, data=data, method=method)
        req.add_header("Authorization", f"Bearer {self.token}")
        req.add_header("Accept", "application/vnd.github+json")
        req.add_header("X-GitHub-Api-Version", "2022-11-28")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read() or b"null")
        except urllib.error.HTTPError as e:
            detail = e.read().decode(errors="replace")
            raise ReportError(f"{method} {path} failed: HTTP {e.code}: {detail}")
        except (urllib.error.URLError, OSError) as e:
            raise ReportError(f"{method} {path} failed: {e}")


def has_marker(body):
    """True when MARKER opens the body, where report() itself places it.

    Only leading whitespace and a UTF-8 BOM are tolerated before it. A marker
    quoted further down (e.g. an issue that documents this script) must not
    make that issue the canonical tracking issue.
    """
    return (body or "").lstrip("﻿ \t\r\n").startswith(MARKER)


def find_tracking_issues(api):
    """All issues (open and closed, not PRs) whose body starts with MARKER."""
    found, page = [], 1
    while True:
        batch = api.request("GET", f"/issues?state=all&per_page=100&page={page}")
        for item in batch:
            if "pull_request" not in item and has_marker(item.get("body")):
                found.append(item)
        if len(batch) < 100:
            return found
        page += 1


def report(api, results, sha, run_url):
    """Apply the lifecycle. `results` maps channel -> status string.

    Returns (exit_code, messages). Raises ReportError on API failure.
    """
    failed = [c for c in CHANNELS if results.get(c) != "success"]
    issues = find_tracking_issues(api)
    if len(issues) > 1:
        nums = ", ".join(f"#{i['number']}" for i in issues)
        raise ReportError(
            f"multiple issues carry {MARKER} ({nums}); refusing to guess. "
            "Close or edit all but one, then re-run."
        )
    issue = issues[0] if issues else None
    msgs = []
    if failed:
        names = ", ".join(failed)
        comment = (
            f"Clippy canary failed on: {names}.\n\n"
            f"- Commit: {sha}\n- Run: {run_url}\n"
        )
        if issue is None:
            body = (
                f"{MARKER}\n\nThe scheduled clippy canary "
                "(`.github/workflows/clippy-canary.yml`) runs "
                "`cargo clippy --workspace --all-targets -- -D warnings` on "
                "latest stable and beta Rust. It is not a required check. See "
                "CONTRIBUTING.md, \"Clippy pin and canary\", for how to respond.\n\n"
                + comment
            )
            api.request("POST", "/issues", {"title": TITLE, "body": body})
            msgs.append("created tracking issue")
        else:
            if issue["state"] != "open":
                api.request("PATCH", f"/issues/{issue['number']}", {"state": "open"})
                msgs.append(f"reopened #{issue['number']}")
            api.request("POST", f"/issues/{issue['number']}/comments", {"body": comment})
            msgs.append(f"commented on #{issue['number']}")
        return 1, msgs
    if issue is not None and issue["state"] == "open":
        n = issue["number"]
        api.request(
            "POST",
            f"/issues/{n}/comments",
            {"body": f"Clippy canary recovered: stable and beta are green.\n\n"
                     f"- Commit: {sha}\n- Run: {run_url}\n"},
        )
        api.request("PATCH", f"/issues/{n}", {"state": "closed"})
        msgs.append(f"commented and closed #{n}")
    else:
        msgs.append("green; no open tracking issue, nothing to do")
    return 0, msgs


def read_results(results_dir):
    out = {}
    for c in CHANNELS:
        path = os.path.join(results_dir, f"{c}.txt")
        if os.path.isfile(path):
            with open(path) as f:
                out[c] = f.read().strip()
    return out


def main():
    try:
        api = GitHubApi(os.environ["GITHUB_REPOSITORY"], os.environ["GITHUB_TOKEN"])
        results = read_results(os.environ["RESULTS_DIR"])
        code, msgs = report(api, results, os.environ["SHA"], os.environ["RUN_URL"])
    except (ReportError, KeyError) as e:
        print(f"::error::clippy canary reporting failed: {e}", file=sys.stderr)
        return 2
    for m in msgs:
        print(m)
    if code:
        print("::error::clippy canary failed; see the tracking issue", file=sys.stderr)
    return code


if __name__ == "__main__":
    sys.exit(main())
