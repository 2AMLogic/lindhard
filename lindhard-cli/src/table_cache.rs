//! Reuse of built electron cross-section tables across runs
//! (`lindhard run --table-cache <DIR>`; `docs/cli.md`, "Cross-section table
//! cache").
//!
//! A table is stored in the versioned cache form of
//! [`lindhard::electron::data::CrossSectionTable`] and read back with that
//! type's loader. Its file name is the SHA-256 of a **key document**: a JSON
//! text that spells out every input the table build depends on. The key is
//! written next to the table (`<kind>-<hash>.key.json`), and on a lookup the
//! stored key must equal the run's key byte for byte, so a hash collision or
//! an edited key file is an error that names the differing field, never a
//! silent reuse. The SHA-256 of the table file's bytes is stored beside both
//! (`<kind>-<hash>.sha256`), and a table whose bytes do not hash to it is an
//! error naming the file, so an edited or corrupted number anywhere in the
//! table (not only in its key, axis or grid) is refused before it is parsed.
//!
//! What the key holds, and why it is conservative:
//!
//! - the key schema version ([`KEY_VERSION`]) and the cache format version
//!   ([`CACHE_FORMAT_VERSION`]);
//! - the build: crate version, `git describe --always --dirty`, and the
//!   SHA-256 of the running executable. Any rebuild that changes the code
//!   (including an uncommitted edit, which `git describe` alone would not
//!   show) changes the executable and so misses. Two builds of the same
//!   source may also miss; a miss only costs a rebuild;
//! - the table kind and the full incident-energy grid (every `f64` written in
//!   its shortest round-trip form, so distinct grids give distinct keys);
//! - the material: its name and its full `Debug` form (composition, atom
//!   fractions, density, exactly);
//! - elastic: the model (`electron.elastic.model`), the atomic potential, the
//!   corrections with all their inputs
//!   (exchange, per-element correlation-polarization), the probability grid
//!   and the refinement tolerance;
//! - inelastic: the Penn algorithm, the model's Fermi energy, the SHA-256 and
//!   provenance of the optical ELF file, and the band parameters of the
//!   material (the table of a material with a band is built on the
//!   band-bottom axis with the band's Fermi energy, `electron::inelastic_axis`),
//!   and, for `mermin-melf`, the options of the material's oscillator fit
//!   (`[electron.materials.<name>.mermin_fit]`, #306).
//!
//! A cache file of another format version can never be found under a
//! current key (the version is in the key). A file that **is** found under
//! the run's key is loaded with the library loader, whose errors (an
//! unsupported version, a failed validation) are returned, not treated as a
//! miss; after loading, its sampling axis and energy grid are checked
//! against the run.
//!
//! Writes go to a temporary file in the same directory and are renamed into
//! place, the table first, then its hash, then its key, so the key is the
//! commit marker: a key file on disk always has its table and hash beside it.
//! A key whose table is missing is rebuilt; a key and table whose hash file is
//! missing is an error, never a silent reuse.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use lindhard::electron::data::{CrossSectionTable, SamplingAxis, CACHE_FORMAT_VERSION};
use lindhard::electron::elastic::table::{default_probability_grid, DEFAULT_REFINE_TOLERANCE};
use lindhard::electron::inelastic::PennAlgorithm;
use lindhard::input::electron::{ResolvedElectron, ResolvedElectronMaterial};
use serde::Serialize;
use sha2::{Digest, Sha256};

/// Version of the key document's schema. Bump it when a field is added,
/// removed or changes meaning.
///
/// 3: `inelastic.mermin_fit`, the Mermin fit options (#306).
pub const KEY_VERSION: u32 = 3;

/// Lowercase hex SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// The two tables of a material.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TableKind {
    /// The elastic (polar-angle) table.
    Elastic,
    /// The inelastic (energy-loss) table.
    Inelastic,
}

impl TableKind {
    /// The stable label (also the file-name prefix).
    pub fn label(self) -> &'static str {
        match self {
            Self::Elastic => "elastic",
            Self::Inelastic => "inelastic",
        }
    }

    fn axis(self) -> SamplingAxis {
        match self {
            Self::Elastic => SamplingAxis::ElasticPolarAngle,
            Self::Inelastic => SamplingAxis::InelasticEnergyLoss,
        }
    }
}

/// Identity of the program that builds tables.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct BuildId {
    /// Crate version.
    pub lindhard_version: String,
    /// `git describe --always --dirty` of the build (`"unknown"` outside a
    /// checkout).
    pub git_describe: String,
    /// SHA-256 of the running executable.
    pub executable_sha256: String,
}

impl BuildId {
    /// The identity of the running executable.
    pub fn current() -> Result<Self> {
        let exe = std::env::current_exe().context("locating the running executable")?;
        let bytes = std::fs::read(&exe).with_context(|| {
            format!(
                "reading the running executable {} to key the table cache",
                exe.display()
            )
        })?;
        Ok(Self {
            lindhard_version: lindhard::VERSION.to_string(),
            git_describe: env!("LINDHARD_GIT_DESCRIBE").to_string(),
            executable_sha256: sha256_hex(&bytes),
        })
    }
}

#[derive(Serialize)]
struct ElasticKey {
    /// `electron.elastic.model`.
    model: &'static str,
    potential: &'static str,
    /// `Debug` of the resolved elastic choice: potential, exchange flag and
    /// every per-element correlation-polarization input.
    choice: String,
    corrections: String,
    /// The starting probability grid and the refinement tolerance of
    /// `combine` (the stored grid is the refined one).
    probability: Vec<f64>,
    refine_tolerance: f64,
}

#[derive(Serialize)]
struct InelasticKey {
    algorithm: &'static str,
    fermi_energy_ev: f64,
    optical_elf_sha256: String,
    optical_elf_provenance: String,
    /// `Debug` of the material's band parameters (`None` without a band).
    band: String,
    /// `Debug` of the material's resolved
    /// [`lindhard::electron::inelastic::MerminFitOptions`] with
    /// `mermin-melf`; `null` with the other models, which do not fit.
    mermin_fit: Option<String>,
}

#[derive(Serialize)]
struct Key<'a> {
    key_version: u32,
    table_format_version: u32,
    build: &'a BuildId,
    kind: TableKind,
    material_name: &'a str,
    /// `Debug` of the [`lindhard::material::Material`]: every field, with
    /// `f64` in shortest round-trip form.
    material: String,
    energy_ev: &'a [f64],
    elastic: Option<ElasticKey>,
    inelastic: Option<InelasticKey>,
}

/// The key document of one table, as compact JSON. Field order is fixed by
/// the struct, and `serde_json` writes each `f64` in its shortest round-trip
/// form, so equal inputs give equal text and different inputs different
/// text.
pub fn key_text(
    build: &BuildId,
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
    kind: TableKind,
) -> String {
    let (elastic, inelastic) = match kind {
        TableKind::Elastic => (
            Some(ElasticKey {
                model: r.elastic.model.label(),
                potential: r.elastic.potential.label(),
                choice: format!("{:?}", r.elastic),
                corrections: r.elastic.corrections_description(),
                probability: default_probability_grid(),
                refine_tolerance: DEFAULT_REFINE_TOLERANCE,
            }),
            None,
        ),
        TableKind::Inelastic => (
            None,
            Some(InelasticKey {
                algorithm: r.inelastic.label(),
                fermi_energy_ev: r.inelastic_fermi_ev,
                optical_elf_sha256: m.optical_elf_file.sha256.clone(),
                optical_elf_provenance: m.optical_elf.provenance().to_string(),
                band: format!("{:?}", m.band),
                mermin_fit: (r.inelastic == PennAlgorithm::Mermin)
                    .then(|| format!("{:?}", m.mermin_fit)),
            }),
        ),
    };
    let key = Key {
        key_version: KEY_VERSION,
        table_format_version: CACHE_FORMAT_VERSION,
        build,
        kind,
        material_name: &m.name,
        material: format!("{:?}", m.material),
        energy_ev: &r.table_energy_ev,
        elastic,
        inelastic,
    };
    serde_json::to_string(&key).expect("the key serializes")
}

/// The first field (as a dotted path) where two JSON documents differ, or
/// `None` if they are equal.
fn first_difference(a: &serde_json::Value, b: &serde_json::Value, path: &str) -> Option<String> {
    use serde_json::Value;
    match (a, b) {
        (Value::Object(x), Value::Object(y)) => {
            for (k, va) in x {
                let p = if path.is_empty() {
                    k.clone()
                } else {
                    format!("{path}.{k}")
                };
                match y.get(k) {
                    Some(vb) => {
                        if let Some(d) = first_difference(va, vb, &p) {
                            return Some(d);
                        }
                    }
                    None => return Some(p),
                }
            }
            y.keys()
                .find(|k| !x.contains_key(*k))
                .map(|k| format!("{path}.{k}").trim_start_matches('.').to_string())
        }
        _ if a == b => None,
        _ => Some(if path.is_empty() {
            "(document)".into()
        } else {
            path.into()
        }),
    }
}

/// Where a table the run used came from.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TableSource {
    /// Built by this run.
    Built,
    /// Read from the table cache.
    Cache,
}

/// The cache file of a table, as echoed in `physics.materials`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct CacheFile {
    /// The table file.
    pub path: PathBuf,
    /// SHA-256 of the table file's bytes.
    pub sha256: String,
    /// SHA-256 of the key document (the file-name stem).
    pub key_sha256: String,
}

/// A table and where it came from.
#[derive(Debug)]
pub struct CachedTable {
    /// The table.
    pub table: CrossSectionTable,
    /// Built or read.
    pub source: TableSource,
    /// The cache file, when a cache is in use.
    pub cache: Option<CacheFile>,
}

/// The files of one cache entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EntryPaths {
    /// `<kind>-<key hash>.toml`: the table.
    pub table: PathBuf,
    /// `<kind>-<key hash>.sha256`: the lowercase hex SHA-256 of the table
    /// file's bytes.
    pub hash: PathBuf,
    /// `<kind>-<key hash>.key.json`: the key document.
    pub key: PathBuf,
}

/// A directory of cached tables.
#[derive(Debug, Clone)]
pub struct TableCache {
    dir: PathBuf,
    build: BuildId,
}

fn write_atomically(path: &Path, bytes: &[u8]) -> Result<()> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .expect("cache file names are ASCII");
    let tmp = path.with_file_name(format!(".{name}.tmp.{}", std::process::id()));
    std::fs::write(&tmp, bytes).with_context(|| format!("writing {}", tmp.display()))?;
    std::fs::rename(&tmp, path)
        .with_context(|| format!("moving {} to {}", tmp.display(), path.display()))
}

impl TableCache {
    /// A cache in `dir` (created if missing), keyed on the running
    /// executable.
    pub fn open(dir: impl Into<PathBuf>) -> Result<Self> {
        Self::with_build(dir, BuildId::current()?)
    }

    /// A cache in `dir` keyed on an explicit build identity.
    pub fn with_build(dir: impl Into<PathBuf>, build: BuildId) -> Result<Self> {
        let dir = dir.into();
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("creating the table cache {}", dir.display()))?;
        let dir = std::fs::canonicalize(&dir).unwrap_or(dir);
        Ok(Self { dir, build })
    }

    /// The cache directory.
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    /// The build identity in every key.
    pub fn build(&self) -> &BuildId {
        &self.build
    }

    /// The files of the entry with key hash `key_sha256`.
    pub fn paths(&self, kind: TableKind, key_sha256: &str) -> EntryPaths {
        let stem = format!("{}-{key_sha256}", kind.label());
        EntryPaths {
            table: self.dir.join(format!("{stem}.toml")),
            hash: self.dir.join(format!("{stem}.sha256")),
            key: self.dir.join(format!("{stem}.key.json")),
        }
    }

    /// The table of `kind` for material `m`: read from the cache when its
    /// key is there, otherwise built with `build` and stored.
    pub fn get_or_build(
        &self,
        r: &ResolvedElectron,
        m: &ResolvedElectronMaterial,
        kind: TableKind,
        build: impl FnOnce() -> Result<CrossSectionTable>,
    ) -> Result<CachedTable> {
        let key = key_text(&self.build, r, m, kind);
        let key_sha256 = sha256_hex(key.as_bytes());
        let paths = self.paths(kind, &key_sha256);
        if let Some(hit) = self.lookup(r, kind, &key, &key_sha256, &paths)? {
            return Ok(hit);
        }
        let table = build()?;
        let text = table
            .to_toml_string()
            .with_context(|| format!("serializing the {} table", kind.label()))?;
        let sha256 = sha256_hex(text.as_bytes());
        // Table, then hash, then key: the key is the commit marker.
        write_atomically(&paths.table, text.as_bytes())?;
        write_atomically(&paths.hash, format!("{sha256}\n").as_bytes())?;
        write_atomically(&paths.key, key.as_bytes())?;
        Ok(CachedTable {
            table,
            source: TableSource::Built,
            cache: Some(CacheFile {
                path: paths.table,
                sha256,
                key_sha256,
            }),
        })
    }

    fn lookup(
        &self,
        r: &ResolvedElectron,
        kind: TableKind,
        key: &str,
        key_sha256: &str,
        paths: &EntryPaths,
    ) -> Result<Option<CachedTable>> {
        let (table_path, key_path) = (paths.table.as_path(), paths.key.as_path());
        let stored = match std::fs::read_to_string(key_path) {
            Ok(s) => s,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(e).with_context(|| format!("reading {}", key_path.display())),
        };
        if stored != key {
            let field = match (
                serde_json::from_str::<serde_json::Value>(&stored),
                serde_json::from_str::<serde_json::Value>(key),
            ) {
                (Ok(a), Ok(b)) => {
                    first_difference(&b, &a, "").unwrap_or_else(|| "(formatting)".into())
                }
                _ => "(not a JSON key document)".into(),
            };
            bail!(
                "table cache {}: the stored key does not match this run (field `{field}` \
                 differs) although its hash does; the cache file was edited or corrupted. \
                 Remove it to rebuild",
                key_path.display()
            );
        }
        let bytes = match std::fs::read(table_path) {
            Ok(b) => b,
            // A key without its table (removed by hand): rebuild.
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(e).with_context(|| format!("reading {}", table_path.display())),
        };
        let stored_sha256 = match std::fs::read_to_string(&paths.hash) {
            Ok(s) => s,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => bail!(
                "table cache {}: the table's SHA-256 file {} is missing, so its bytes \
                 cannot be verified. Remove the entry to rebuild",
                table_path.display(),
                paths.hash.display()
            ),
            Err(e) => return Err(e).with_context(|| format!("reading {}", paths.hash.display())),
        };
        let sha256 = sha256_hex(&bytes);
        if stored_sha256.trim_end() != sha256 {
            bail!(
                "table cache {}: the file's SHA-256 is {sha256}, but {} records {}; \
                 the table was edited or corrupted. Remove the entry to rebuild",
                table_path.display(),
                paths.hash.display(),
                stored_sha256.trim_end()
            );
        }
        let text = String::from_utf8(bytes.clone())
            .with_context(|| format!("{} is not UTF-8 text", table_path.display()))?;
        let table = CrossSectionTable::from_toml_str(&text)
            .with_context(|| format!("loading the cached table {}", table_path.display()))?;
        if table.axis() != kind.axis() {
            bail!(
                "table cache {}: field `axis` is {:?}, this run needs {:?}",
                table_path.display(),
                table.axis(),
                kind.axis()
            );
        }
        if table.energy_ev() != r.table_energy_ev.as_slice() {
            bail!(
                "table cache {}: field `energy_ev` differs from this run's grid \
                 ({} energies cached, {} needed)",
                table_path.display(),
                table.energy_ev().len(),
                r.table_energy_ev.len()
            );
        }
        Ok(Some(CachedTable {
            table,
            source: TableSource::Cache,
            cache: Some(CacheFile {
                path: table_path.to_path_buf(),
                sha256,
                key_sha256: key_sha256.to_string(),
            }),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lindhard::input::electron::{ElectronInput, PotentialChoice};

    fn example() -> ResolvedElectron {
        example_with(|t| t.to_string())
    }

    /// The example with its text changed by `edit`.
    fn example_with(edit: impl Fn(&str) -> String) -> ResolvedElectron {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples/electron");
        let text = std::fs::read_to_string(dir.join("e_10keV_si.toml")).unwrap();
        ElectronInput::from_toml_str(&edit(&text))
            .unwrap()
            .resolve_in(&dir)
            .unwrap()
    }

    /// The example with `model = "mermin-melf"` and, if given, a
    /// `mermin_fit` table for its material.
    fn mermin_example(fit: Option<&str>) -> ResolvedElectron {
        example_with(|t| {
            let t = t.replace("model = \"penn-single-pole\"", "model = \"mermin-melf\"");
            assert!(t.contains("mermin-melf"), "the example names its model");
            match fit {
                None => t,
                Some(f) => format!("{t}\n[electron.materials.Si.mermin_fit]\n{f}\n"),
            }
        })
    }

    fn build() -> BuildId {
        BuildId {
            lindhard_version: "0.0.0".into(),
            git_describe: "test".into(),
            executable_sha256: "00".into(),
        }
    }

    fn keys(r: &ResolvedElectron) -> (String, String) {
        let m = &r.materials[0];
        (
            key_text(&build(), r, m, TableKind::Elastic),
            key_text(&build(), r, m, TableKind::Inelastic),
        )
    }

    /// Changes `r` with `f` and asserts which of the two keys change.
    fn assert_changes(
        name: &str,
        elastic: bool,
        inelastic: bool,
        f: impl Fn(&mut ResolvedElectron),
    ) {
        let base = example();
        let (e0, i0) = keys(&base);
        let mut r = base.clone();
        f(&mut r);
        let (e1, i1) = keys(&r);
        assert_eq!(
            e0 != e1,
            elastic,
            "{name}: elastic key changed = {}",
            e0 != e1
        );
        assert_eq!(
            i0 != i1,
            inelastic,
            "{name}: inelastic key changed = {}",
            i0 != i1
        );
    }

    #[test]
    fn key_is_deterministic_and_kinds_differ() {
        let (e, i) = keys(&example());
        assert_eq!((e.clone(), i.clone()), keys(&example()));
        assert_ne!(e, i);
        assert!(e.contains("\"table_format_version\""));
    }

    #[test]
    fn key_changes_with_the_inelastic_model() {
        assert_changes("model", false, true, |r| {
            r.inelastic = lindhard::electron::inelastic::PennAlgorithm::Full;
        });
    }

    /// `fermi_energy_ev` keys the tables of materials without a band (with a
    /// band it must be 0, and the band keys the table).
    #[test]
    fn key_changes_with_the_fermi_energy() {
        assert_changes("fermi energy", false, true, |r| r.inelastic_fermi_ev = 1.0);
        // Even in the last bit.
        assert_changes("fermi energy ulp", false, true, |r| {
            r.inelastic_fermi_ev = f64::from_bits(1);
        });
    }

    /// The Mermin fit options key the `mermin-melf` table (#306): any
    /// change of the options changes it; an explicit table equal to the
    /// defaults keys the same table as no table (the key holds the resolved
    /// options); and the other models, which do not fit, carry no options.
    #[test]
    fn key_changes_with_the_mermin_fit_options() {
        use lindhard::electron::inelastic::{FitWeighting, MerminFitOptions};
        let base = mermin_example(None);
        assert_eq!(base.materials[0].mermin_fit, MerminFitOptions::default());
        let (e0, i0) = keys(&base);
        assert!(
            i0.contains(&format!(
                "\"mermin_fit\":{}",
                serde_json::to_string(&format!("{:?}", MerminFitOptions::default())).unwrap()
            )),
            "{i0}"
        );
        let explicit = mermin_example(Some(
            "oscillators = 3\nweighting = \"relative\"\nrelative_floor = 0.01",
        ));
        assert_eq!(keys(&explicit), (e0.clone(), i0.clone()));
        for (name, fit) in [
            ("oscillators", "oscillators = 6"),
            ("weighting", "weighting = \"uniform\""),
            ("relative_floor", "relative_floor = 0.02"),
        ] {
            let (e1, i1) = keys(&mermin_example(Some(fit)));
            assert_eq!(e1, e0, "{name}: the elastic key changed");
            assert_ne!(i1, i0, "{name}: the inelastic key did not change");
        }
        // Even in the last bit of the floor.
        let mut r = base.clone();
        r.materials[0].mermin_fit.weighting = FitWeighting::Relative {
            floor: f64::from_bits(0.01f64.to_bits() + 1),
        };
        assert_ne!(keys(&r).1, i0);
        // Without the Mermin model there is no fit to key.
        let single = example();
        let (_, i) = keys(&single);
        assert!(i.contains("\"mermin_fit\":null"), "{i}");
        let mut r = single.clone();
        r.materials[0].mermin_fit.n_oscillators = 6;
        assert_eq!(keys(&r).1, i);
    }

    #[test]
    fn key_changes_with_the_grid() {
        assert_changes("grid point", true, true, |r| {
            let x = &mut r.table_energy_ev[3];
            *x = f64::from_bits(x.to_bits() + 1);
        });
        assert_changes("grid length", true, true, |r| {
            r.table_energy_ev.pop();
        });
    }

    #[test]
    fn key_changes_with_the_elf_file() {
        assert_changes("elf sha256", false, true, |r| {
            r.materials[0].optical_elf_file.sha256 = "ff".repeat(32);
        });
    }

    #[test]
    fn key_changes_with_the_potential_and_corrections() {
        assert_changes("potential", true, false, |r| {
            r.elastic.potential = PotentialChoice::SalvatDhfs;
        });
        assert_changes("exchange", true, false, |r| {
            r.elastic.exchange = !r.elastic.exchange;
        });
    }

    #[test]
    fn key_changes_with_the_material_and_band() {
        assert_changes("density", true, true, |r| {
            let m = &r.materials[0].material;
            let spec = lindhard::material::MaterialSpec {
                density_g_cm3: Some(2.33),
                ..lindhard::material::MaterialSpec::from(m.clone())
            };
            r.materials[0].material = spec.try_into().unwrap();
        });
        assert_changes("name", true, true, |r| r.materials[0].name = "Si2".into());
        assert_changes("band", false, true, |r| r.materials[0].band = None);
        // The inelastic table of a material with a band is built with the
        // band's Fermi energy (its minimum excitation energy, #241), so a
        // band that moves it keys a different table, even in the last bit.
        use lindhard::electron::boundary::{BandModel, BandStructure};
        let moved = |r: &mut ResolvedElectron, f: &dyn Fn(f64) -> f64| {
            let b = r.materials[0].band.clone().unwrap();
            let BandModel::Insulator {
                valence_band_width_ev,
                band_gap_ev,
                affinity_ev,
            } = b.model()
            else {
                panic!("the example's band is an insulator");
            };
            let model = BandModel::Insulator {
                valence_band_width_ev: f(valence_band_width_ev),
                band_gap_ev,
                affinity_ev,
            };
            r.materials[0].band = Some(BandStructure::new(model, b.provenance()).unwrap());
            assert_ne!(
                r.materials[0].band.as_ref().unwrap().min_excitation_ev(),
                b.min_excitation_ev()
            );
        };
        assert_changes("band fermi energy", false, true, |r| moved(r, &|v| v + 0.5));
        assert_changes("band fermi energy ulp", false, true, |r| {
            moved(r, &|v| f64::from_bits(v.to_bits() + 1))
        });
    }

    #[test]
    fn key_changes_with_the_build() {
        let r = example();
        let m = &r.materials[0];
        let a = key_text(&build(), &r, m, TableKind::Inelastic);
        for b in [
            BuildId {
                executable_sha256: "01".into(),
                ..build()
            },
            BuildId {
                git_describe: "test-dirty".into(),
                ..build()
            },
            BuildId {
                lindhard_version: "0.0.1".into(),
                ..build()
            },
        ] {
            assert_ne!(a, key_text(&b, &r, m, TableKind::Inelastic), "{b:?}");
        }
    }

    #[test]
    fn elastic_key_carries_the_model() {
        let r = example();
        let (e, i) = keys(&r);
        let label = r.elastic.model.label();
        assert_eq!(label, "mott");
        assert!(e.contains(&format!("\"model\":\"{label}\"")), "{e}");
        assert!(e.contains(&format!("\"key_version\":{KEY_VERSION}")), "{e}");
        assert!(!i.contains("\"mott\""), "{i}");
    }

    /// Every field of the inputs that shape a table is accounted for by a
    /// field of the key. The destructurings have no `..`, so a field added
    /// to `ElasticSpec`, `CorrelationPolarizationSpec`, `PolarizabilitySpec`,
    /// `InelasticSpec`, `TableGridSpec`, `ElectronMaterialSpec` or
    /// `MerminFitSpec` fails to compile here until it is classified (and, if
    /// it changes a table, added to the key).
    #[test]
    fn key_covers_every_table_input_field() {
        use lindhard::input::electron::{
            CorrelationPolarizationSpec, ElasticSpec, ElectronMaterialSpec, InelasticSpec,
            MerminFitSpec, PolarizabilitySpec, TableGridSpec,
        };
        let r = example();
        let (e, i) = keys(&r);
        let ElasticSpec {
            // `elastic.model`.
            model,
            // `elastic.potential` (and `elastic.choice`).
            potential,
            // `elastic.choice` (`exchange`) and `elastic.corrections`.
            exchange,
            // `elastic.choice` (every per-element input) and
            // `elastic.corrections`.
            correlation_polarization,
        } = &r.input.electron.elastic;
        assert!(e.contains(&format!("\"model\":\"{}\"", model.label())));
        assert!(e.contains(&format!("\"potential\":\"{}\"", potential.label())));
        assert!(e.contains(&format!("exchange: {exchange}")), "{e}");
        if let Some(CorrelationPolarizationSpec {
            // `elastic.choice` (`PolarizationCutoff`).
            b_pol_squared: _,
            // `elastic.choice` (`outer_radius`).
            outer_radius_bohr: _,
            // `elastic.choice` (per element).
            polarizability,
        }) = correlation_polarization
        {
            for PolarizabilitySpec {
                // `elastic.choice` (`polarizability`).
                bohr3: _,
                // `elastic.choice` (`polarizability_source`).
                source: _,
            } in polarizability.values()
            {}
        }
        let InelasticSpec {
            // `inelastic.algorithm`.
            model,
            // `inelastic.fermi_energy_ev`.
            fermi_energy_ev,
        } = &r.input.electron.inelastic;
        assert!(i.contains(&format!("\"algorithm\":\"{model}\"")), "{i}");
        assert!(i.contains(&format!(
            "\"fermi_energy_ev\":{}",
            serde_json::to_string(fermi_energy_ev).unwrap()
        )));
        // All three set the grid, which is `energy_ev` of both keys.
        let TableGridSpec {
            min_energy_ev: _,
            max_energy_ev: _,
            points_per_decade: _,
        } = &r.input.electron.tables;
        let grid = serde_json::to_string(&r.table_energy_ev).unwrap();
        assert!(e.contains(&format!("\"energy_ev\":{grid}")));
        assert!(i.contains(&format!("\"energy_ev\":{grid}")));
        // The per-material data, with a Mermin fit table so that every
        // field is present.
        let r = mermin_example(Some(
            "oscillators = 5\nweighting = \"relative\"\nrelative_floor = 0.03",
        ));
        let (_, i) = keys(&r);
        let ElectronMaterialSpec {
            // `inelastic.optical_elf_sha256` and `.optical_elf_provenance`.
            optical_elf: _,
            // `inelastic.band`.
            band,
            // Not table inputs: the phonon and polaron channels are separate
            // channels of the transport, not part of either table.
            phonon: _,
            polaron: _,
            // `inelastic.mermin_fit`.
            mermin_fit,
        } = &r.input.electron.materials["Si"];
        assert!(band.is_some());
        assert!(i.contains(&format!(
            "\"band\":{}",
            serde_json::to_string(&format!("{:?}", r.materials[0].band)).unwrap()
        )));
        let MerminFitSpec {
            // `inelastic.mermin_fit` (`n_oscillators`).
            oscillators,
            // `inelastic.mermin_fit` (`weighting`).
            weighting: _,
            // `inelastic.mermin_fit` (`floor` of `Relative`).
            relative_floor,
        } = mermin_fit.as_ref().unwrap();
        assert!(i.contains(&format!("n_oscillators: {oscillators}")), "{i}");
        assert!(
            i.contains(&format!(
                "Relative {{ floor: {:?} }}",
                relative_floor.unwrap()
            )),
            "{i}"
        );
    }

    #[test]
    fn first_difference_names_the_field() {
        let a: serde_json::Value = serde_json::from_str(r#"{"x":1,"y":{"z":2}}"#).unwrap();
        let b: serde_json::Value = serde_json::from_str(r#"{"x":1,"y":{"z":3}}"#).unwrap();
        assert_eq!(first_difference(&a, &b, "").as_deref(), Some("y.z"));
        assert_eq!(first_difference(&a, &a, ""), None);
    }
}
