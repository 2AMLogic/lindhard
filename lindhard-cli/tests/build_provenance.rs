//! Freshness of the `LINDHARD_GIT_DESCRIBE` stamp across incremental builds.
//!
//! `build.rs` is copied into a throwaway workspace with the same layout as this
//! one (`lindhard/`, `lindhard-cli/`), under a fresh git repository and its own
//! `CARGO_TARGET_DIR`. The fixture binary prints the stamp. Each scenario then
//! changes one thing and runs `cargo build` again WITHOUT `cargo clean`, so a
//! stale stamp (a build script that was not rerun) is observable.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};
use std::thread::sleep;
use std::time::Duration;

const BUILD_RS: &str = include_str!("../build.rs");

/// Throwaway directory under the system temp dir (not under this repository,
/// whose own `.git` would otherwise be found by the non-git scenario).
struct Tmp(PathBuf);

impl Tmp {
    fn new(tag: &str) -> Self {
        static N: AtomicU32 = AtomicU32::new(0);
        let d = std::env::temp_dir().join(format!(
            "lindhard-prov-{}-{}-{tag}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Tmp(d)
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn write(path: &Path, text: &str) {
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    std::fs::write(path, text).unwrap();
}

/// Lets filesystem mtimes advance past the previous build before the next edit.
fn tick() {
    sleep(Duration::from_millis(100));
}

fn git_cmd(dir: &Path, args: &[&str]) -> Command {
    let mut c = Command::new("git");
    c.current_dir(dir)
        .args(args)
        .env("GIT_CONFIG_GLOBAL", "/dev/null")
        .env("GIT_CONFIG_NOSYSTEM", "1")
        .env("GIT_OPTIONAL_LOCKS", "0")
        .env("GIT_AUTHOR_NAME", "t")
        .env("GIT_AUTHOR_EMAIL", "t@example.invalid")
        .env("GIT_COMMITTER_NAME", "t")
        .env("GIT_COMMITTER_EMAIL", "t@example.invalid");
    c
}

fn git(dir: &Path, args: &[&str]) -> String {
    let o = git_cmd(dir, args).output().unwrap();
    assert!(
        o.status.success(),
        "git {args:?}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    String::from_utf8(o.stdout).unwrap().trim().to_string()
}

fn git_available() -> bool {
    Command::new("git").arg("--version").output().is_ok()
}

/// Writes the fixture workspace files (no git).
fn write_workspace(root: &Path) {
    write(
        &root.join("Cargo.toml"),
        "[workspace]\nresolver = \"2\"\nmembers = [\"lindhard\", \"lindhard-cli\"]\n",
    );
    write(
        &root.join("lindhard/Cargo.toml"),
        "[package]\nname = \"fixture-lib\"\nversion = \"0.0.0\"\nedition = \"2021\"\n",
    );
    write(
        &root.join("lindhard/src/lib.rs"),
        "pub fn id() -> u32 {\n    1\n}\n",
    );
    write(
        &root.join("lindhard-cli/Cargo.toml"),
        "[package]\nname = \"fixture-cli\"\nversion = \"0.0.0\"\nedition = \"2021\"\nbuild = \"build.rs\"\n\n\
         [[bin]]\nname = \"stamp\"\npath = \"src/main.rs\"\n\n\
         [dependencies]\nfixture-lib = { path = \"../lindhard\" }\n",
    );
    write(&root.join("lindhard-cli/build.rs"), BUILD_RS);
    write(
        &root.join("lindhard-cli/src/main.rs"),
        "fn main() {\n    println!(\"{} {}\", env!(\"LINDHARD_GIT_DESCRIBE\"), fixture_lib::id());\n}\n",
    );
}

/// A git repository with one commit holding the fixture workspace.
fn git_fixture(root: &Path) {
    write_workspace(root);
    git(root, &["init", "-q", "-b", "main"]);
    // Fix the lockfile before committing so building never dirties the tree.
    let o = Command::new(cargo())
        .current_dir(root)
        .args(["generate-lockfile", "--offline"])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    git(root, &["add", "-A"]);
    git(root, &["commit", "-q", "-m", "init"]);
}

fn cargo() -> PathBuf {
    std::env::var_os("CARGO").map_or_else(|| "cargo".into(), PathBuf::from)
}

/// A checkout plus its private target directory; builds are incremental.
struct Checkout<'a> {
    dir: &'a Path,
    target: PathBuf,
}

impl<'a> Checkout<'a> {
    fn new(dir: &'a Path, target: PathBuf) -> Self {
        Checkout { dir, target }
    }

    fn build(&self, describe_override: Option<&str>) -> Output {
        let mut c = Command::new(cargo());
        c.current_dir(self.dir)
            .args(["build", "--offline", "--quiet", "-j", "2"])
            .env("CARGO_TARGET_DIR", &self.target)
            .env_remove("LINDHARD_GIT_DESCRIBE")
            .env_remove("RUSTC_WRAPPER");
        if let Some(v) = describe_override {
            c.env("LINDHARD_GIT_DESCRIBE", v);
        }
        let o = c.output().unwrap();
        assert!(
            o.status.success(),
            "cargo build failed:\n{}",
            String::from_utf8_lossy(&o.stderr)
        );
        o
    }

    /// Rebuilds and returns the stamp the binary reports.
    fn stamp_with(&self, describe_override: Option<&str>) -> String {
        self.build(describe_override);
        let o = Command::new(self.target.join("debug/stamp"))
            .output()
            .unwrap();
        assert!(o.status.success());
        let out = String::from_utf8(o.stdout).unwrap();
        out.split_whitespace().next().unwrap().to_string()
    }

    fn stamp(&self) -> String {
        self.stamp_with(None)
    }

    /// What `git describe --always --dirty` says right now.
    fn expected(&self) -> String {
        git(self.dir, &["describe", "--always", "--dirty"])
    }

    fn assert_fresh(&self, what: &str) -> String {
        let got = self.stamp();
        assert_eq!(got, self.expected(), "stale stamp after {what}");
        got
    }
}

fn tag(dir: &Path, name: &str) {
    git(dir, &["tag", "-a", "-m", name, name]);
}

fn commit_file(dir: &Path, file: &str, text: &str, msg: &str) {
    write(&dir.join(file), text);
    git(dir, &["add", "-A"]);
    git(dir, &["commit", "-q", "-m", msg]);
}

#[test]
fn unstaged_source_edit_and_restore_refresh_the_stamp() {
    if !git_available() {
        return;
    }
    let tmp = Tmp::new("edit");
    let repo = tmp.0.join("repo");
    git_fixture(&repo);
    let co = Checkout::new(&repo, tmp.0.join("target"));

    let clean = co.assert_fresh("initial build");
    assert!(!clean.ends_with("-dirty"), "{clean}");

    let lib = repo.join("lindhard/src/lib.rs");
    let original = std::fs::read_to_string(&lib).unwrap();
    tick();
    std::fs::write(&lib, original.replace("1\n", "2\n")).unwrap();
    let dirty = co.assert_fresh("unstaged library edit");
    assert_eq!(dirty, format!("{clean}-dirty"));

    tick();
    std::fs::write(&lib, &original).unwrap();
    assert_eq!(co.assert_fresh("restoring the library source"), clean);

    // The CLI's own source, and a manifest, count too.
    let main = repo.join("lindhard-cli/src/main.rs");
    let main_text = std::fs::read_to_string(&main).unwrap();
    tick();
    std::fs::write(&main, format!("{main_text}// edit\n")).unwrap();
    assert_eq!(co.assert_fresh("CLI source edit"), dirty);
    tick();
    std::fs::write(&main, &main_text).unwrap();
    assert_eq!(co.assert_fresh("CLI source restore"), clean);

    let manifest = repo.join("lindhard/Cargo.toml");
    let manifest_text = std::fs::read_to_string(&manifest).unwrap();
    tick();
    std::fs::write(&manifest, format!("{manifest_text}\n# edit\n")).unwrap();
    assert_eq!(co.assert_fresh("manifest edit"), dirty);
    tick();
    std::fs::write(&manifest, &manifest_text).unwrap();
    assert_eq!(co.assert_fresh("manifest restore"), clean);
}

#[test]
fn branch_advance_and_tag_changes_refresh_the_stamp() {
    if !git_available() {
        return;
    }
    let tmp = Tmp::new("refs");
    let repo = tmp.0.join("repo");
    git_fixture(&repo);
    let co = Checkout::new(&repo, tmp.0.join("target"));
    let first = co.assert_fresh("initial build");

    // The branch ref moves; .git/HEAD (a symbolic ref) is untouched.
    tick();
    commit_file(&repo, "NOTES", "a\n", "second");
    let second = co.assert_fresh("branch advance");
    assert_ne!(first, second);

    tick();
    tag(&repo, "v1");
    assert_eq!(co.assert_fresh("new tag"), "v1");

    tick();
    commit_file(&repo, "NOTES", "b\n", "third");
    assert!(co.assert_fresh("commit after tag").starts_with("v1-1-g"));

    // Packed-ref layout: loose refs and tags are folded into packed-refs.
    tick();
    git(&repo, &["pack-refs", "--all", "--prune"]);
    assert!(co.assert_fresh("pack-refs").starts_with("v1-1-g"));

    tick();
    tag(&repo, "v2");
    assert_eq!(co.assert_fresh("tag on top of packed refs"), "v2");

    tick();
    commit_file(&repo, "NOTES", "c\n", "fourth");
    assert!(co
        .assert_fresh("branch advance over packed-refs")
        .starts_with("v2-1-g"));

    tick();
    git(&repo, &["tag", "-d", "v2"]);
    assert!(co.assert_fresh("packed tag deletion").starts_with("v1-2-g"));

    // Moving a tag (delete + recreate) is a ref change too.
    tick();
    git(&repo, &["tag", "-d", "v1"]);
    tag(&repo, "v1");
    assert_eq!(co.assert_fresh("tag recreated at the tip"), "v1");

    // Switching branches changes HEAD itself.
    tick();
    git(&repo, &["checkout", "-q", "-b", "side"]);
    commit_file(&repo, "NOTES", "d\n", "side");
    assert!(co.assert_fresh("side branch commit").starts_with("v1-1-g"));
}

#[test]
fn linked_worktree_follows_shared_refs_and_its_own_branch() {
    if !git_available() {
        return;
    }
    let tmp = Tmp::new("wt");
    let repo = tmp.0.join("repo");
    git_fixture(&repo);
    let linked = tmp.0.join("linked");
    git(
        &repo,
        &[
            "worktree",
            "add",
            "-q",
            "-b",
            "other",
            linked.to_str().unwrap(),
        ],
    );
    let co = Checkout::new(&linked, tmp.0.join("target-linked"));
    let first = co.assert_fresh("initial linked build");

    // A tag made from the main checkout lives in the common git dir.
    tick();
    tag(&repo, "shared1");
    assert_eq!(co.assert_fresh("tag created elsewhere"), "shared1");

    // The linked worktree's branch ref is in the common dir as well.
    tick();
    commit_file(&linked, "NOTES", "a\n", "linked commit");
    assert!(co
        .assert_fresh("linked branch advance")
        .starts_with("shared1-1-g"));

    // Packed layout, then change the shared tag from the main checkout.
    tick();
    git(&repo, &["pack-refs", "--all", "--prune"]);
    assert!(co.assert_fresh("pack-refs").starts_with("shared1-1-g"));
    tick();
    git(&repo, &["tag", "-d", "shared1"]);
    let after = co.assert_fresh("shared tag deleted");
    assert_ne!(after, first);
    assert!(!after.starts_with("shared1"), "{after}");

    // An unstaged edit in the linked worktree.
    let lib = linked.join("lindhard/src/lib.rs");
    let text = std::fs::read_to_string(&lib).unwrap();
    tick();
    std::fs::write(&lib, text.replace("1\n", "3\n")).unwrap();
    assert_eq!(co.assert_fresh("linked edit"), format!("{after}-dirty"));
    tick();
    std::fs::write(&lib, &text).unwrap();
    assert_eq!(co.assert_fresh("linked restore"), after);
}

#[test]
fn environment_override_is_exact_and_authoritative() {
    if !git_available() {
        return;
    }
    let tmp = Tmp::new("override");
    let repo = tmp.0.join("repo");
    git_fixture(&repo);
    let co = Checkout::new(&repo, tmp.0.join("target"));
    let real = co.assert_fresh("initial build");

    assert_eq!(co.stamp_with(Some("rel-1.2.3")), "rel-1.2.3");
    // Whitespace is trimmed; a blank override falls back to git.
    assert_eq!(co.stamp_with(Some("  rel-9\n")), "rel-9");
    tick();
    commit_file(&repo, "NOTES", "x\n", "later");
    assert_eq!(
        co.stamp_with(Some("rel-9")),
        "rel-9",
        "git must not leak in"
    );
    assert_eq!(co.stamp_with(Some("  ")), co.expected());
    assert_ne!(co.expected(), real);
    // Dropping the override returns to the git description.
    assert_eq!(co.stamp_with(Some("pinned")), "pinned");
    assert_eq!(co.assert_fresh("override removed"), co.expected());
}

#[test]
fn exported_tree_without_git_builds_with_unknown_fallback() {
    let tmp = Tmp::new("nogit");
    let tree = tmp.0.join("export");
    write_workspace(&tree);
    let co = Checkout::new(&tree, tmp.0.join("target"));
    // `git` itself may exist; it just finds no repository above the export.
    assert_eq!(co.stamp(), "unknown");
    // Rebuilding incrementally stays valid and does not loop on missing paths.
    assert_eq!(co.stamp(), "unknown");
    assert_eq!(co.stamp_with(Some("from-env")), "from-env");
    assert_eq!(co.stamp(), "unknown");
}
