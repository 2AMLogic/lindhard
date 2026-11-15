//! The one data-file loader both input resolvers use.
//!
//! An input may name files on disk: a user stopping table
//! (`stopping.tables[i]`), an SRIM coefficient file (`stopping.srim_data`),
//! an optical ELF or a subshell-binding table of an electron material. They
//! are read the same way — joined with the input file's directory, read,
//! checked for UTF-8, hashed — and recorded the same way, so that a result
//! identifies the exact bytes it was produced from (`CONTRIBUTING.md`:
//! provenance of every dataset).

use std::path::{Path, PathBuf};

use serde::Serialize;
use sha2::{Digest, Sha256};

use super::{invalid, InputError};

/// A data file a run read, identified for provenance.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DataFile {
    /// The path as written in the input.
    pub path: String,
    /// Where it was read from (absolute where possible).
    pub resolved_path: PathBuf,
    /// SHA-256 of the file's bytes, lowercase hex.
    pub sha256: String,
}

/// Lowercase hex SHA-256 of `bytes`.
pub(crate) fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// Read the data file `path` relative to `base_dir`, returning its text and
/// its provenance record. `field` is the dotted input key, for error
/// messages.
pub(crate) fn read_data_file(
    base_dir: &Path,
    field: &str,
    path: &str,
) -> Result<(String, DataFile), InputError> {
    let joined = base_dir.join(path);
    let bytes = std::fs::read(&joined)
        .map_err(|e| invalid(field, format!("cannot read {}: {e}", joined.display())))?;
    let text = String::from_utf8(bytes.clone()).map_err(|_| {
        invalid(
            field,
            format!("{} is not valid UTF-8 text", joined.display()),
        )
    })?;
    Ok((
        text,
        DataFile {
            path: path.to_string(),
            resolved_path: std::fs::canonicalize(&joined).unwrap_or(joined),
            sha256: sha256_hex(&bytes),
        },
    ))
}
