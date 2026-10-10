//! The output lifecycle shared by the ion, electron and dynamic run modes.
//!
//! A run replaces its CSVs in place and publishes its summary last, so the
//! summary certifies the completed output set. [`OutputSet::begin`] removes
//! the mode's previous summary right before the first CSV is touched, so a
//! failure partway through never leaves an old summary beside replaced CSVs;
//! [`OutputSet::publish`] writes the new summary to a temporary sibling and
//! renames it into place, so a final summary is never truncated.
//!
//! Only the files the caller names are touched. Concurrent writers to one
//! directory are not supported.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

/// An output directory in the middle of a run, for one mode's summary file.
#[derive(Debug)]
pub struct OutputSet {
    out: PathBuf,
    summary_name: String,
    summary: String,
}

impl OutputSet {
    /// Invalidates the previous summary `summary_name` in `out` (a missing
    /// file is fine; any other failure aborts, naming the path, before any
    /// output changes). `summary` is the already-serialized new summary, held
    /// until [`publish`](Self::publish). The directory must exist.
    pub fn begin(out: &Path, summary_name: &str, summary: String) -> Result<Self> {
        let p = out.join(summary_name);
        match std::fs::remove_file(&p) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                return Err(e).with_context(|| format!("removing previous {}", p.display()));
            }
            _ => {}
        }
        Ok(Self {
            out: out.to_path_buf(),
            summary_name: summary_name.to_string(),
            summary,
        })
    }

    /// Writes (overwrites) the required file `name`.
    pub fn write(&self, name: &str, text: String) -> Result<()> {
        let p = self.out.join(name);
        std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
    }

    /// Writes the optional file `name` if `text` is `Some`; otherwise removes
    /// the file left by an earlier run. A missing file is fine; any other
    /// failure is an error naming the path.
    pub fn reconcile_optional(&self, name: &str, text: Option<String>) -> Result<()> {
        let p = self.out.join(name);
        match text {
            Some(t) => self.write(name, t),
            None => match std::fs::remove_file(&p) {
                Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                    Err(e).with_context(|| format!("removing stale {}", p.display()))
                }
                _ => Ok(()),
            },
        }
    }

    /// Publishes the summary: written to a temporary sibling, then renamed
    /// into place. The temporary file is removed on failure.
    pub fn publish(self) -> Result<()> {
        let tmp = self.out.join(format!(".{}.tmp", self.summary_name));
        let dst = self.out.join(&self.summary_name);
        let result = std::fs::write(&tmp, &self.summary)
            .with_context(|| format!("writing {}", tmp.display()))
            .and_then(|()| {
                std::fs::rename(&tmp, &dst).with_context(|| format!("writing {}", dst.display()))
            });
        if result.is_err() {
            let _ = std::fs::remove_file(&tmp);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dir(name: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("lindhard-publish-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn begin_removes_old_summary_and_publish_writes_new() {
        let d = dir("ok");
        std::fs::write(d.join("s.json"), "old").unwrap();
        let o = OutputSet::begin(&d, "s.json", "new".into()).unwrap();
        assert!(!d.join("s.json").exists());
        o.write("a.csv", "x".into()).unwrap();
        o.publish().unwrap();
        assert_eq!(std::fs::read_to_string(d.join("s.json")).unwrap(), "new");
        assert!(!d.join(".s.json.tmp").exists());
    }

    #[test]
    fn begin_failure_names_path_and_changes_nothing() {
        let d = dir("begin-fail");
        std::fs::create_dir_all(d.join("s.json/inner")).unwrap();
        std::fs::write(d.join("a.csv"), "keep").unwrap();
        let e = OutputSet::begin(&d, "s.json", "new".into()).err().unwrap();
        assert!(format!("{e:#}").contains("s.json"));
        assert_eq!(std::fs::read_to_string(d.join("a.csv")).unwrap(), "keep");
    }

    #[test]
    fn failed_publish_leaves_no_final_summary_or_temp() {
        let d = dir("publish-fail");
        let o = OutputSet::begin(&d, "s.json", "new".into()).unwrap();
        // A directory at the destination makes the rename fail.
        std::fs::create_dir_all(d.join("s.json/inner")).unwrap();
        assert!(o.publish().is_err());
        assert!(!d.join(".s.json.tmp").exists());
        assert!(d.join("s.json").is_dir());
    }
}
