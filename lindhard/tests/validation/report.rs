//! One row of the validation table, and the table printers.

use std::fmt::Write as _;

/// Outcome of a check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Status {
    /// Measured value within its tolerance.
    Pass,
    /// Measured value outside its tolerance: the harness exits non-zero.
    Fail,
    /// Reported, not asserted (a known deviation or a context number). The
    /// note says why it is not asserted.
    Info,
}

impl Status {
    fn label(self) -> &'static str {
        match self {
            Status::Pass => "pass",
            Status::Fail => "FAIL",
            Status::Info => "info",
        }
    }
}

/// One row of the table.
#[derive(Debug, Clone)]
pub(crate) struct Check {
    /// Stable identifier, used by the command-line filter.
    pub id: &'static str,
    /// What is compared with what.
    pub what: String,
    /// The measured value, formatted.
    pub value: String,
    /// The tolerance it is held to, formatted ("-" for `Info`).
    pub tolerance: String,
    /// Outcome.
    pub status: Status,
    /// Reference, caveat or pointer to the issue that tracks a deviation.
    pub note: String,
}

impl Check {
    /// A check that passes when `measured <= tol`.
    pub(crate) fn at_most(
        id: &'static str,
        what: impl Into<String>,
        measured: f64,
        tol: f64,
        fmt: fn(f64) -> String,
        note: impl Into<String>,
    ) -> Self {
        Self {
            id,
            what: what.into(),
            value: fmt(measured),
            tolerance: format!("<= {}", fmt(tol)),
            // NaN fails.
            status: if measured <= tol {
                Status::Pass
            } else {
                Status::Fail
            },
            note: note.into(),
        }
    }

    /// A boolean check.
    pub(crate) fn holds(
        id: &'static str,
        what: impl Into<String>,
        ok: bool,
        value: impl Into<String>,
        note: impl Into<String>,
    ) -> Self {
        Self {
            id,
            what: what.into(),
            value: value.into(),
            tolerance: "exact".into(),
            status: if ok { Status::Pass } else { Status::Fail },
            note: note.into(),
        }
    }

    /// A reported, unasserted number.
    pub(crate) fn info(
        id: &'static str,
        what: impl Into<String>,
        value: impl Into<String>,
        note: impl Into<String>,
    ) -> Self {
        Self {
            id,
            what: what.into(),
            value: value.into(),
            tolerance: "-".into(),
            status: Status::Info,
            note: note.into(),
        }
    }
}

/// Relative deviation as a percentage with three significant figures.
pub(crate) fn pct(x: f64) -> String {
    format!("{:.3}%", 100.0 * x)
}

/// Signed percentage.
pub(crate) fn spct(x: f64) -> String {
    format!("{:+.2}%", 100.0 * x)
}

/// Scientific notation, two decimals.
pub(crate) fn sci(x: f64) -> String {
    format!("{x:.2e}")
}

/// Plain number, four significant decimals.
pub(crate) fn num(x: f64) -> String {
    format!("{x:.4}")
}

/// Plain-text table for the terminal.
pub(crate) fn text_table(rows: &[Check]) -> String {
    let mut s = String::new();
    let w_id = rows.iter().map(|r| r.id.len()).max().unwrap_or(2).max(2);
    let w_val = rows.iter().map(|r| r.value.len()).max().unwrap_or(5).max(5);
    let w_tol = rows
        .iter()
        .map(|r| r.tolerance.len())
        .max()
        .unwrap_or(9)
        .max(9);
    let _ = writeln!(
        s,
        "{:<w_id$}  {:<4}  {:>w_val$}  {:<w_tol$}  what",
        "id", "", "value", "tolerance"
    );
    for r in rows {
        let _ = writeln!(
            s,
            "{:<w_id$}  {:<4}  {:>w_val$}  {:<w_tol$}  {}",
            r.id,
            r.status.label(),
            r.value,
            r.tolerance,
            r.what
        );
    }
    s
}

/// Markdown table for `docs/validation.md`. Contains no timings or other
/// machine-dependent text, so a rerun on the same platform is idempotent.
pub(crate) fn markdown_table(rows: &[Check]) -> String {
    let esc = |t: &str| t.replace('|', "\\|");
    let mut s = String::new();
    s.push_str("| Check | What | Value | Tolerance | Result | Notes |\n");
    s.push_str("|---|---|---|---|---|---|\n");
    for r in rows {
        let _ = writeln!(
            s,
            "| `{}` | {} | {} | {} | {} | {} |",
            r.id,
            esc(&r.what),
            esc(&r.value),
            esc(&r.tolerance),
            r.status.label(),
            esc(&r.note)
        );
    }
    s
}
