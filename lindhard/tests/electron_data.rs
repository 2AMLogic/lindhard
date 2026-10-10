//! Contract tests for `lindhard::electron::data`: validated constructors,
//! loaders that refuse data without provenance, serde reads that go through
//! the same validation, and the cross-section cache format.
//!
//! Every numeric fixture in this file is **synthetic**: made up for the tests,
//! not physical data, and labelled so in its provenance string. The one test
//! that reads real data is the coverage test of the committed EADL2017 table
//! (`eadl2017_binding_energies_cover_z_1_to_92`).

use lindhard::electron::data::{
    AtomBindings, CrossSectionTable, CrossSectionTableParts, ElectronDataError, OpticalElf,
    SamplingAxis, ShellBinding, ShellChannelTable, Subshell, SubshellBindingTable,
    CACHE_FORMAT_VERSION, SHELL_CHANNEL_FORMAT_VERSION,
};

const SYNTHETIC: &str = "synthetic test fixture, not physical data";

// ---------------------------------------------------------------------------
// OpticalElf
// ---------------------------------------------------------------------------

const ELF_TOML: &str = r#"
material = "synthetic-material"
provenance = "synthetic test fixture, not physical data"
energy_ev = [1.0, 2.0, 4.0, 10.0]
elf = [0.0, 0.0, 2.0, 0.5]
"#;

fn elf() -> OpticalElf {
    OpticalElf::from_toml_str(ELF_TOML).unwrap()
}

#[test]
fn elf_reproduces_knots_exactly() {
    let t = elf();
    for (e, y) in t.energy_ev().iter().zip(t.elf_values()) {
        assert_eq!(t.elf(*e).unwrap(), *y, "knot {e} eV");
    }
}

#[test]
fn elf_interpolates_linearly_between_knots() {
    let t = elf();
    assert_eq!(t.elf(3.0).unwrap(), 1.0); // midpoint of (2, 0) .. (4, 2)
    assert!((t.elf(7.0).unwrap() - 1.25).abs() < 1e-15); // midpoint of (4, 2) .. (10, 0.5)
    assert!((t.elf(2.5).unwrap() - 0.5).abs() < 1e-15);
}

#[test]
fn elf_handles_zero_values() {
    let t = elf();
    // A zero-valued segment stays exactly zero (log-log interpolation could not).
    for e in [1.0, 1.25, 1.5, 1.999, 2.0] {
        assert_eq!(t.elf(e).unwrap(), 0.0, "{e} eV");
    }
    // Interpolation never goes negative.
    for i in 0..=900 {
        let e = 1.0 + f64::from(i) * 0.01;
        assert!(t.elf(e).unwrap() >= 0.0);
    }
}

#[test]
fn elf_endpoints_are_in_range_and_beyond_is_an_error() {
    let t = elf();
    assert_eq!(t.energy_range_ev(), (1.0, 10.0));
    assert_eq!(t.elf(1.0).unwrap(), 0.0);
    assert_eq!(t.elf(10.0).unwrap(), 0.5);
    for e in [0.999_999, 10.000_001, 0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(t.elf(e), Err(ElectronDataError::OutOfRange { .. })),
            "{e} eV"
        );
    }
}

#[test]
fn elf_loader_refuses_missing_provenance() {
    let no = ELF_TOML.replace(
        "provenance = \"synthetic test fixture, not physical data\"\n",
        "",
    );
    assert_eq!(
        OpticalElf::from_toml_str(&no).unwrap_err(),
        ElectronDataError::MissingProvenance
    );
    assert_eq!(
        OpticalElf::new("Si", "", vec![1.0, 2.0], vec![0.0, 1.0]).unwrap_err(),
        ElectronDataError::MissingProvenance
    );
}

#[test]
fn elf_loader_refuses_whitespace_only_provenance() {
    for blank in ["", " ", "   \t  "] {
        let bad = ELF_TOML.replace("synthetic test fixture, not physical data", blank);
        assert_eq!(
            OpticalElf::from_toml_str(&bad).unwrap_err(),
            ElectronDataError::MissingProvenance,
            "{blank:?}"
        );
    }
}

#[test]
fn elf_loader_refuses_malformed_provenance() {
    // Control characters: a provenance is one line of text.
    let ctl = ELF_TOML.replace(
        "synthetic test fixture, not physical data",
        "first line\\nsecond line",
    );
    assert!(matches!(
        OpticalElf::from_toml_str(&ctl),
        Err(ElectronDataError::MalformedText {
            field: "provenance",
            ..
        })
    ));
    // Wrong type.
    let num = ELF_TOML.replace("\"synthetic test fixture, not physical data\"", "42");
    assert!(matches!(
        OpticalElf::from_toml_str(&num),
        Err(ElectronDataError::Parse(_))
    ));
    // A misspelt key is not silently ignored (and so not a missing provenance
    // that slips through under another name).
    let typo = ELF_TOML.replace("provenance =", "provenence =");
    assert!(matches!(
        OpticalElf::from_toml_str(&typo),
        Err(ElectronDataError::Parse(_))
    ));
}

#[test]
fn elf_requires_material_identity() {
    let no = ELF_TOML.replace("material = \"synthetic-material\"\n", "");
    assert_eq!(
        OpticalElf::from_toml_str(&no).unwrap_err(),
        ElectronDataError::MissingIdentity("material")
    );
}

#[test]
fn elf_rejects_bad_grids_and_values() {
    let cases = [
        (
            "[1.0, 2.0, 4.0, 10.0]",
            "[0.0, 0.0, 2.0]",
            "mismatched lengths",
        ),
        (
            "[1.0, 2.0, 2.0, 10.0]",
            "[0.0, 0.0, 2.0, 0.5]",
            "duplicate point",
        ),
        (
            "[1.0, 4.0, 2.0, 10.0]",
            "[0.0, 0.0, 2.0, 0.5]",
            "descending",
        ),
        (
            "[1.0, nan, 4.0, 10.0]",
            "[0.0, 0.0, 2.0, 0.5]",
            "nan energy",
        ),
        ("[1.0, 2.0, 4.0, inf]", "[0.0, 0.0, 2.0, 0.5]", "inf energy"),
        (
            "[0.0, 2.0, 4.0, 10.0]",
            "[0.0, 0.0, 2.0, 0.5]",
            "zero energy",
        ),
        (
            "[-1.0, 2.0, 4.0, 10.0]",
            "[0.0, 0.0, 2.0, 0.5]",
            "negative energy",
        ),
        (
            "[1.0, 2.0, 4.0, 10.0]",
            "[0.0, -1e-9, 2.0, 0.5]",
            "negative ELF",
        ),
        ("[1.0, 2.0, 4.0, 10.0]", "[0.0, nan, 2.0, 0.5]", "nan ELF"),
        ("[1.0, 2.0, 4.0, 10.0]", "[0.0, 0.0, inf, 0.5]", "inf ELF"),
        ("[1.0]", "[0.0]", "single point"),
        ("[]", "[]", "empty"),
    ];
    for (e, y, why) in cases {
        let bad = ELF_TOML
            .replace("[1.0, 2.0, 4.0, 10.0]", e)
            .replace("[0.0, 0.0, 2.0, 0.5]", y);
        assert!(
            matches!(
                OpticalElf::from_toml_str(&bad),
                Err(ElectronDataError::Invalid { .. })
            ),
            "{why}: {:?}",
            OpticalElf::from_toml_str(&bad)
        );
    }
}

#[test]
fn elf_serde_reads_go_through_validation() {
    // A derived deserializer must not bypass the checks, whatever the format.
    let bad_json =
        r#"{"material": "m", "provenance": " ", "energy_ev": [1.0, 2.0], "elf": [0.0, 1.0]}"#;
    assert!(serde_json::from_str::<OpticalElf>(bad_json).is_err());
    let descending = r#"{"material": "m", "provenance": "synthetic", "energy_ev": [2.0, 1.0], "elf": [0.0, 1.0]}"#;
    assert!(serde_json::from_str::<OpticalElf>(descending).is_err());
    let good = elf();
    let json = serde_json::to_string(&good).unwrap();
    assert_eq!(serde_json::from_str::<OpticalElf>(&json).unwrap(), good);
    let toml = toml::to_string(&good).unwrap();
    assert_eq!(OpticalElf::from_toml_str(&toml).unwrap(), good);
    assert_eq!(good.provenance(), SYNTHETIC);
    assert_eq!(good.material(), "synthetic-material");
}

// ---------------------------------------------------------------------------
// Subshell binding energies
// ---------------------------------------------------------------------------

/// One 80-column ENDF-6 card: six 11-column fields, MAT, MF, MT, sequence.
fn card(fields: [&str; 6], mat: u32, mf: u32, mt: u32) -> String {
    let mut s: String = fields.iter().map(|f| format!("{f:>11}")).collect();
    s.push_str(&format!("{mat:>4}{mf:>2}{mt:>3}{:>5}", 1));
    s
}

/// A File 28 section for one element. `shells` are `(SUBI, EBI, ELN)`; each
/// gets one synthetic radiative transition so the transition lines are
/// exercised too. All numbers are synthetic.
fn mf28_section(z: u32, shells: &[(&str, &str, &str)]) -> Vec<String> {
    let mat = z * 100;
    let za = format!("{}.0", z * 1000);
    let nss = shells.len().to_string();
    let mut out = vec![card([&za, "1.0", "0", "0", &nss, "0"], mat, 28, 533)];
    for (subi, ebi, eln) in shells {
        out.push(card([subi, "0.0", "0", "0", "12", "1"], mat, 28, 533));
        out.push(card([ebi, eln, "0.0", "0.0", "0.0", "0.0"], mat, 28, 533));
        out.push(card(
            ["2.0", "0.0", "1.0+1", "1.0", "0.0", "0.0"],
            mat,
            28,
            533,
        ));
    }
    // SEND, FEND: other MF/MT, ignored by the reader.
    out.push(card(["", "", "", "", "", ""], mat, 28, 0));
    out.push(card(["", "", "", "", "", ""], mat, 0, 0));
    out
}

/// A synthetic two-element ENDF-6 File 28 text, with a File 1 header card to
/// show that other files are skipped, and one element beyond Z = 92.
fn synthetic_mf28() -> String {
    let mut lines = vec![card(["1.0+3", "1.0", "-1", "0", "0", "0"], 100, 1, 451)];
    lines.extend(mf28_section(1, &[("1.0", "1.000000+1", "1.0")]));
    lines.extend(mf28_section(
        3,
        &[
            ("1.0", "5.000000+1", "2.0"),
            ("2.0", "5.0D+0", "0.5"),
            ("4.0", "4.0", "0.5"),
        ],
    ));
    lines.extend(mf28_section(93, &[("1.0", "1.0+5", "93.0")]));
    lines.join("\n")
}

#[test]
fn endf6_mf28_reader_reads_binding_energies_and_occupancies() {
    let t = SubshellBindingTable::from_endf6_mf28_str(&synthetic_mf28(), SYNTHETIC).unwrap();
    assert_eq!(t.provenance(), SYNTHETIC);
    // Z = 93 is beyond the element table and skipped.
    assert_eq!(t.atoms().iter().map(|a| a.z()).collect::<Vec<_>>(), [1, 3]);
    let li = t.atom(3).unwrap();
    assert_eq!(li.binding_energy_ev(Subshell::K), Some(50.0));
    assert_eq!(
        li.binding_energy_ev(Subshell::from_label("L1").unwrap()),
        Some(5.0)
    );
    assert_eq!(
        li.binding_energy_ev(Subshell::from_label("L3").unwrap()),
        Some(4.0)
    );
    assert_eq!(
        li.shell(Subshell::from_label("L3").unwrap())
            .unwrap()
            .occupancy(),
        0.5
    );
    // L2 is absent (not listed), so it has no binding energy.
    assert_eq!(
        li.binding_energy_ev(Subshell::from_label("L2").unwrap()),
        None
    );
    assert_eq!(t.atom(2), None);
}

#[test]
fn endf6_mf28_reader_requires_provenance() {
    for p in ["", "  ", "\t"] {
        assert_eq!(
            SubshellBindingTable::from_endf6_mf28_str(&synthetic_mf28(), p).unwrap_err(),
            ElectronDataError::MissingProvenance
        );
    }
}

#[test]
fn endf6_mf28_reader_rejects_malformed_records() {
    let good = synthetic_mf28();
    let rejects = |text: &str, expect: &str| {
        let err = SubshellBindingTable::from_endf6_mf28_str(text, SYNTHETIC)
            .unwrap_err()
            .to_string();
        assert!(err.contains(expect), "expected {expect:?} in {err:?}");
    };
    let edit = |from: &str, to: &str| {
        let out = good.replacen(&format!("{from:>11}"), &format!("{to:>11}"), 1);
        assert_ne!(out, good, "fixture edit {from:?} -> {to:?} did nothing");
        out
    };
    // NW inconsistent with NTR (first subshell LIST record: NW = 12, NTR = 1).
    let bad = good.replacen(
        &format!("{:>11}{:>11}", "12", "1"),
        &format!("{:>11}{:>11}", "18", "1"),
        1,
    );
    assert_ne!(bad, good);
    rejects(&bad, "expected NW = 6 (1 + NTR)");
    // A non-numeric field.
    rejects(&edit("1.000000+1", "1.0x0000+1"), "is not a number");
    // A non-positive binding energy.
    rejects(&edit("1.000000+1", "0.0"), "finite and positive");
    // A duplicate subshell (Li's L1 LIST record relabelled as K).
    let l1_head = format!("{:>11}{:>11}{:>11}", "2.0", "0.0", "0");
    let k_head = format!("{:>11}{:>11}{:>11}", "1.0", "0.0", "0");
    let dup = good.replacen(&l1_head, &k_head, 1);
    assert_ne!(dup, good);
    rejects(&dup, "twice");
    // Occupancies that do not make a neutral atom.
    rejects(&edit("0.5", "0.25"), "neutral atom");
    // Not an elemental ZA.
    rejects(&edit("3000.0", "3007.0"), "elemental ZA");
    // Truncated section: drop the last card of the Li section's last subshell.
    let lines: Vec<&str> = good.lines().collect();
    let li_last = lines
        .iter()
        .rposition(|l| l[66..70].trim() == "300" && &l[70..75] == "28533")
        .unwrap();
    let mut cut = lines.clone();
    cut.remove(li_last);
    rejects(&cut.join("\n"), "truncated");
    // No File 28 at all.
    assert!(matches!(
        SubshellBindingTable::from_endf6_mf28_str("not an ENDF file", SYNTHETIC),
        Err(ElectronDataError::Parse(_))
    ));
}

fn shell(d: u8, ebi: f64, occ: f64) -> ShellBinding {
    ShellBinding::new(Subshell::new(d).unwrap(), ebi, occ).unwrap()
}

#[test]
fn atom_validation() {
    // Valid.
    assert!(AtomBindings::new(2, vec![shell(1, 20.0, 2.0)]).is_ok());
    // Unknown Z.
    assert!(AtomBindings::new(0, vec![shell(1, 20.0, 2.0)]).is_err());
    assert!(AtomBindings::new(93, vec![shell(1, 20.0, 2.0)]).is_err());
    // No shells.
    assert!(AtomBindings::new(1, vec![]).is_err());
    // Duplicate subshell.
    let dup = AtomBindings::new(4, vec![shell(1, 100.0, 2.0), shell(1, 90.0, 2.0)]);
    assert!(format!("{}", dup.unwrap_err()).contains("twice"));
    // Out of designator order.
    assert!(AtomBindings::new(4, vec![shell(2, 10.0, 2.0), shell(1, 100.0, 2.0)]).is_err());
    // Not a neutral atom.
    assert!(AtomBindings::new(3, vec![shell(1, 50.0, 2.0)]).is_err());
    // Shell values.
    for (e, o) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (f64::INFINITY, 1.0),
        (1.0, 0.0),
        (1.0, 3.0),
        (1.0, f64::NAN),
    ] {
        assert!(ShellBinding::new(Subshell::K, e, o).is_err(), "{e} eV, {o}");
    }
    // Occupancy up to the capacity 2j + 1: 4 for L3 (2p3/2), not 5.
    let l3 = Subshell::from_label("L3").unwrap();
    assert!(ShellBinding::new(l3, 10.0, 4.0).is_ok());
    assert!(ShellBinding::new(l3, 10.0, 4.5).is_err());
}

#[test]
fn binding_table_validation_coverage_and_round_trip() {
    let h = AtomBindings::new(1, vec![shell(1, 10.0, 1.0)]).unwrap();
    let he = AtomBindings::new(2, vec![shell(1, 20.0, 2.0)]).unwrap();
    assert!(SubshellBindingTable::new(SYNTHETIC, vec![he.clone(), h.clone()]).is_err());
    assert!(SubshellBindingTable::new(SYNTHETIC, vec![h.clone(), h.clone()]).is_err());
    assert!(SubshellBindingTable::new(SYNTHETIC, vec![]).is_err());
    assert_eq!(
        SubshellBindingTable::new(" ", vec![h.clone()]).unwrap_err(),
        ElectronDataError::MissingProvenance
    );
    let t = SubshellBindingTable::new(SYNTHETIC, vec![h, he]).unwrap();
    assert!(t.require_coverage(1..=2).is_ok());
    assert_eq!(t.missing(1..=4), vec![3, 4]);
    assert!(t.require_coverage(1..=92).is_err());

    let text = t.to_toml_string().unwrap();
    assert_eq!(SubshellBindingTable::from_toml_str(&text).unwrap(), t);
    let json = serde_json::to_string(&t).unwrap();
    assert_eq!(
        serde_json::from_str::<SubshellBindingTable>(&json).unwrap(),
        t
    );

    // Reading from disk goes through the same validation.
    let no_prov = text.replace(&format!("provenance = \"{SYNTHETIC}\"\n"), "");
    assert_eq!(
        SubshellBindingTable::from_toml_str(&no_prov).unwrap_err(),
        ElectronDataError::MissingProvenance
    );
    let neg = text.replacen("binding_energy_ev = 10.0", "binding_energy_ev = -10.0", 1);
    assert_ne!(neg, text);
    assert!(SubshellBindingTable::from_toml_str(&neg).is_err());
    let bad_json = json.replace("\"subshell\":1", "\"subshell\":40");
    assert_ne!(bad_json, json);
    assert!(serde_json::from_str::<SubshellBindingTable>(&bad_json).is_err());
}

/// Coverage of the committed EADL2017 table: always runs. Z = 1..92 present
/// (occupancies sum to Z, enforced by the constructor), the K shell present
/// and rising with Z, and hydrogen's K binding energy within 0.1 % of half
/// the CODATA Hartree energy.
#[test]
fn eadl2017_binding_energies_cover_z_1_to_92() {
    let t = SubshellBindingTable::eadl2017();
    assert!(t.provenance().contains("IAEA-NDS-224"));
    check_eadl2017_coverage(&t);
}

/// The same checks on a local copy of the source file named by
/// `LINDHARD_EADL2017`, and the committed table must agree with it.
#[test]
fn eadl2017_committed_table_matches_local_source_file() {
    let Ok(path) = std::env::var("LINDHARD_EADL2017") else {
        eprintln!("skipped: set LINDHARD_EADL2017 to a local EADL2017.ALL to run");
        return;
    };
    let t = SubshellBindingTable::from_endf6_mf28_file(&path, "EADL2017 local copy").unwrap();
    check_eadl2017_coverage(&t);
    assert_eq!(t.atoms(), SubshellBindingTable::eadl2017().atoms());
}

fn check_eadl2017_coverage(t: &SubshellBindingTable) {
    t.require_coverage(1..=92).unwrap();
    assert_eq!(t.atoms().len(), 92);
    let mut previous_k = 0.0;
    for a in t.atoms() {
        let occupancy: f64 = a.shells().iter().map(|s| s.occupancy()).sum();
        assert!(
            (occupancy - f64::from(a.z())).abs() < 1e-9,
            "occupancies of Z = {} sum to {occupancy}",
            a.z()
        );
        let k = a
            .binding_energy_ev(Subshell::K)
            .unwrap_or_else(|| panic!("Z = {} has no K shell", a.z()));
        assert!(k > previous_k, "K binding not increasing at Z = {}", a.z());
        previous_k = k;
    }
    let rydberg_ev = lindhard::constants::HARTREE_ENERGY / 2.0 / lindhard::units::J_PER_EV;
    let h_k = t.atom(1).unwrap().binding_energy_ev(Subshell::K).unwrap();
    assert!((h_k / rydberg_ev - 1.0).abs() < 1e-3, "H K = {h_k} eV");
}

// ---------------------------------------------------------------------------
// Cross-section cache
// ---------------------------------------------------------------------------

fn inelastic_parts() -> CrossSectionTableParts {
    CrossSectionTableParts {
        model: "synthetic-model v0".into(),
        material: "synthetic-material".into(),
        provenance: SYNTHETIC.into(),
        axis: SamplingAxis::InelasticEnergyLoss,
        energy_ev: vec![10.0, 100.0, 1000.0],
        // Zero at 10 eV: below a (synthetic) threshold, no interaction.
        inverse_mfp_per_m: vec![0.0, 1.0e9, 2.0e9],
        probability: vec![0.0, 0.5, 1.0],
        quantiles: vec![vec![], vec![0.0, 20.0, 100.0], vec![1.0, 30.0, 500.0]],
    }
}

fn elastic_parts() -> CrossSectionTableParts {
    CrossSectionTableParts {
        axis: SamplingAxis::ElasticPolarAngle,
        inverse_mfp_per_m: vec![3.0e9, 2.0e9, 1.0e9],
        quantiles: vec![
            vec![0.0, 1.0, std::f64::consts::PI],
            vec![0.0, 0.5, 3.0],
            vec![0.0, 0.1, 2.0],
        ],
        ..inelastic_parts()
    }
}

fn inelastic() -> CrossSectionTable {
    CrossSectionTable::new(inelastic_parts()).unwrap()
}

#[test]
fn cache_round_trip_preserves_everything() {
    for t in [
        inelastic(),
        CrossSectionTable::new(elastic_parts()).unwrap(),
    ] {
        let text = t.to_toml_string().unwrap();
        assert!(text.contains(&format!("format_version = {CACHE_FORMAT_VERSION}")));
        assert!(text.contains(&format!("ordinate_unit = \"{}\"", t.axis().unit())));
        let back = CrossSectionTable::from_toml_str(&text).unwrap();
        assert_eq!(back, t);
        assert_eq!(back.model(), "synthetic-model v0");
        assert_eq!(back.material(), "synthetic-material");
        assert_eq!(back.provenance(), SYNTHETIC);
        assert_eq!(back.format_version(), CACHE_FORMAT_VERSION);
        assert_eq!(back.axis(), t.axis());
        assert_eq!(back.parts(), t.parts());
        // Bit-exact values.
        for (a, b) in back.inverse_mfp_per_m().iter().zip(t.inverse_mfp_per_m()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        let json = serde_json::to_string(&t).unwrap();
        assert_eq!(serde_json::from_str::<CrossSectionTable>(&json).unwrap(), t);
    }
}

#[test]
fn cache_file_round_trip() {
    let t = inelastic();
    let dir = std::env::temp_dir().join(format!("lindhard-electron-data-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("inelastic.toml");
    t.write_toml_file(&path).unwrap();
    assert_eq!(CrossSectionTable::from_toml_file(&path).unwrap(), t);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn elastic_and_inelastic_tables_are_distinguishable() {
    let el = CrossSectionTable::new(elastic_parts()).unwrap();
    let inel = inelastic();
    assert_eq!(el.axis(), SamplingAxis::ElasticPolarAngle);
    assert_eq!(inel.axis(), SamplingAxis::InelasticEnergyLoss);
    assert_eq!(el.axis().unit(), "rad");
    assert_eq!(inel.axis().unit(), "eV");
    // A file whose unit does not match its axis is rejected.
    let text = el
        .to_toml_string()
        .unwrap()
        .replace("ordinate_unit = \"rad\"", "ordinate_unit = \"eV\"");
    assert!(matches!(
        CrossSectionTable::from_toml_str(&text),
        Err(ElectronDataError::MalformedText {
            field: "ordinate_unit",
            ..
        })
    ));
    // An unknown axis is rejected.
    let text = el
        .to_toml_string()
        .unwrap()
        .replace("elastic_polar_angle", "polar_angle");
    assert!(CrossSectionTable::from_toml_str(&text).is_err());
}

#[test]
fn cache_rejects_unsupported_versions_explicitly() {
    let text = inelastic().to_toml_string().unwrap();
    let v2 = text.replace("format_version = 1", "format_version = 2");
    assert_eq!(
        CrossSectionTable::from_toml_str(&v2).unwrap_err(),
        ElectronDataError::UnsupportedVersion {
            found: Some(2),
            supported: CACHE_FORMAT_VERSION
        }
    );
    // A future version with fields this build does not know still reports
    // the version, not a confusing field error.
    let v2_new_field = format!("future_field = [1, 2]\n{v2}");
    assert!(matches!(
        CrossSectionTable::from_toml_str(&v2_new_field),
        Err(ElectronDataError::UnsupportedVersion { found: Some(2), .. })
    ));
    let none = text.replace("format_version = 1\n", "");
    assert_eq!(
        CrossSectionTable::from_toml_str(&none).unwrap_err(),
        ElectronDataError::UnsupportedVersion {
            found: None,
            supported: CACHE_FORMAT_VERSION
        }
    );
    // An unknown field in a version-1 file is an error, not ignored.
    let extra = format!("surprise = 1\n{text}");
    assert!(matches!(
        CrossSectionTable::from_toml_str(&extra),
        Err(ElectronDataError::Parse(_))
    ));
    // The serde path (any format) also rejects it.
    let json = serde_json::to_string(&inelastic())
        .unwrap()
        .replace("\"format_version\":1", "\"format_version\":7");
    assert!(serde_json::from_str::<CrossSectionTable>(&json).is_err());
}

#[test]
fn cache_requires_provenance_and_identity() {
    let text = inelastic().to_toml_string().unwrap();
    let no_prov = text.replace(&format!("provenance = \"{SYNTHETIC}\"\n"), "");
    assert_ne!(no_prov, text);
    assert_eq!(
        CrossSectionTable::from_toml_str(&no_prov).unwrap_err(),
        ElectronDataError::MissingProvenance
    );
    let blank = text.replace(SYNTHETIC, "   ");
    assert_eq!(
        CrossSectionTable::from_toml_str(&blank).unwrap_err(),
        ElectronDataError::MissingProvenance
    );
    let no_model = text.replace("synthetic-model v0", "");
    assert_eq!(
        CrossSectionTable::from_toml_str(&no_model).unwrap_err(),
        ElectronDataError::MissingIdentity("model")
    );
    let no_material = text.replace("synthetic-material", "");
    assert_eq!(
        CrossSectionTable::from_toml_str(&no_material).unwrap_err(),
        ElectronDataError::MissingIdentity("material")
    );
}

#[test]
fn cache_rejects_invalid_contents() {
    type Edit = fn(&mut CrossSectionTableParts);
    let cases: [(&str, Edit); 19] = [
        ("imfp length", |p| {
            p.inverse_mfp_per_m.pop();
        }),
        ("negative imfp", |p| p.inverse_mfp_per_m[1] = -1.0),
        ("nan imfp", |p| p.inverse_mfp_per_m[1] = f64::NAN),
        ("inf imfp", |p| p.inverse_mfp_per_m[2] = f64::INFINITY),
        ("descending energies", |p| p.energy_ev.swap(0, 1)),
        ("duplicate energies", |p| p.energy_ev[1] = p.energy_ev[0]),
        ("nan energy", |p| p.energy_ev[2] = f64::NAN),
        ("probability not from 0", |p| p.probability[0] = 0.1),
        ("probability not to 1", |p| p.probability[2] = 0.99),
        ("probability not increasing", |p| p.probability[1] = 1.0),
        ("probability nan", |p| p.probability[1] = f64::NAN),
        ("probability one point", |p| {
            p.probability = vec![1.0];
        }),
        ("ragged row", |p| {
            p.quantiles[1].pop();
        }),
        ("row count", |p| {
            p.quantiles.pop();
        }),
        ("zero-rate row not empty", |p| {
            p.quantiles[0] = vec![0.0, 1.0, 2.0]
        }),
        ("non-zero-rate row empty", |p| p.quantiles[1].clear()),
        ("decreasing row", |p| p.quantiles[2] = vec![1.0, 0.5, 500.0]),
        ("loss above incident energy", |p| p.quantiles[1][2] = 100.5),
        ("negative loss", |p| p.quantiles[2][0] = -1.0),
    ];
    for (why, edit) in cases {
        let mut p = inelastic_parts();
        edit(&mut p);
        assert!(
            matches!(
                CrossSectionTable::new(p.clone()),
                Err(ElectronDataError::Invalid { .. })
            ),
            "{why}"
        );
        // The same data read from a file is rejected the same way.
        let raw = format!(
            "format_version = 1\nmodel = \"m\"\nmaterial = \"x\"\nprovenance = \"synthetic\"\n\
             axis = \"inelastic_energy_loss\"\nordinate_unit = \"eV\"\n\
             energy_ev = {:?}\ninverse_mfp_per_m = {:?}\nprobability = {:?}\nquantiles = {:?}\n",
            p.energy_ev, p.inverse_mfp_per_m, p.probability, p.quantiles
        )
        .replace("NaN", "nan");
        assert!(
            CrossSectionTable::from_toml_str(&raw).is_err(),
            "{why} (file)"
        );
    }
    // Elastic angles beyond pi.
    let mut p = elastic_parts();
    p.quantiles[0][2] = 3.2;
    assert!(CrossSectionTable::new(p).is_err());
}

#[test]
fn inverse_cdf_sampling_contract() {
    let t = inelastic();
    // Exact at the stored probability points.
    assert_eq!(t.inverse_cdf(1, 0.0).unwrap(), 0.0);
    assert_eq!(t.inverse_cdf(1, 0.5).unwrap(), 20.0);
    assert_eq!(t.inverse_cdf(1, 1.0).unwrap(), 100.0);
    // Linear between them.
    assert_eq!(t.inverse_cdf(1, 0.25).unwrap(), 10.0);
    assert_eq!(t.inverse_cdf(2, 0.75).unwrap(), 265.0);
    // Outside [0, 1] or past the grid: errors.
    for u in [-0.1, 1.1, f64::NAN] {
        assert!(matches!(
            t.inverse_cdf(1, u),
            Err(ElectronDataError::OutOfRange { .. })
        ));
    }
    assert!(matches!(
        t.inverse_cdf(3, 0.5),
        Err(ElectronDataError::IndexOutOfRange { index: 3, len: 3 })
    ));
}

#[test]
fn zero_rate_rows_are_never_sampled() {
    let t = inelastic();
    assert_eq!(t.inverse_mfp_per_m()[0], 0.0);
    assert_eq!(t.quantiles(0), None);
    for u in [0.0, 0.5, 1.0] {
        assert_eq!(
            t.inverse_cdf(0, u).unwrap_err(),
            ElectronDataError::ZeroRate { energy_ev: 10.0 }
        );
    }
    // And they survive a round trip as empty rows.
    let back = CrossSectionTable::from_toml_str(&t.to_toml_string().unwrap()).unwrap();
    assert_eq!(back.quantiles(0), None);
}

// ---------------------------------------------------------------------------
// ShellChannelTable
// ---------------------------------------------------------------------------

/// A synthetic inner-shell loss table with binding energy 100 eV: zero rate
/// at 10 eV, losses in `[100, 100]` at 100 eV and `[100, 150]` above.
fn shell_loss_table(axis: SamplingAxis, lowest_loss: f64) -> CrossSectionTable {
    let energy_ev = vec![10.0, 100.0, 1000.0, 10_000.0];
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic shell channel".into(),
        material: "synthetic-material".into(),
        provenance: SYNTHETIC.into(),
        axis,
        energy_ev,
        inverse_mfp_per_m: vec![0.0, 1e7, 2e7, 2e7],
        probability: vec![0.0, 0.5, 1.0],
        quantiles: vec![
            vec![],
            vec![100.0, 100.0, 100.0],
            vec![lowest_loss, 125.0, 150.0],
            vec![lowest_loss, 125.0, 150.0],
        ],
    })
    .unwrap()
}

fn shell_channel() -> ShellChannelTable {
    ShellChannelTable::new(
        14,
        Subshell::from_label("L3").unwrap(),
        100.0,
        "synthetic binding energy, not physical data",
        shell_loss_table(SamplingAxis::InelasticEnergyLoss, 100.0),
    )
    .unwrap()
}

#[test]
fn shell_channel_table_round_trips_identity_binding_and_provenance() {
    let t = shell_channel();
    assert_eq!(t.format_version(), SHELL_CHANNEL_FORMAT_VERSION);
    let text = t.to_toml_string().unwrap();
    assert!(text.contains("subshell = \"L3\""), "{text}");
    let back = ShellChannelTable::from_toml_str(&text).unwrap();
    assert_eq!(back, t);
    assert_eq!((back.z(), back.subshell().label()), (14, "L3"));
    assert_eq!(back.binding_energy_ev(), 100.0);
    assert_eq!(
        back.binding_provenance(),
        "synthetic binding energy, not physical data"
    );
    assert_eq!(back.table().provenance(), SYNTHETIC);
    // Through a file too.
    let path = std::env::temp_dir().join(format!(
        "lindhard-shell-channel-{}.toml",
        std::process::id()
    ));
    t.write_toml_file(&path).unwrap();
    let from_file = ShellChannelTable::from_toml_file(&path);
    let _ = std::fs::remove_file(&path);
    assert_eq!(from_file.unwrap(), t);
}

#[test]
fn shell_channel_table_rejects_stale_or_incompatible_files() {
    let text = shell_channel().to_toml_string().unwrap();
    // The wrapper's own version is checked first.
    let stale = text.replacen(
        &format!("format_version = {SHELL_CHANNEL_FORMAT_VERSION}"),
        "format_version = 999",
        1,
    );
    assert!(matches!(
        ShellChannelTable::from_toml_str(&stale),
        Err(ElectronDataError::UnsupportedVersion {
            found: Some(999),
            ..
        })
    ));
    let missing = text.replacen(
        &format!("format_version = {SHELL_CHANNEL_FORMAT_VERSION}\n"),
        "",
        1,
    );
    assert!(matches!(
        ShellChannelTable::from_toml_str(&missing),
        Err(ElectronDataError::UnsupportedVersion { found: None, .. })
    ));
    // The nested cross-section table keeps its own version check, with the
    // typed error.
    let (head, tail) = text.split_at(text.find("[table]").unwrap());
    let inner = tail.replacen(
        &format!("format_version = {CACHE_FORMAT_VERSION}"),
        "format_version = 998",
        1,
    );
    assert!(matches!(
        ShellChannelTable::from_toml_str(&format!("{head}{inner}")),
        Err(ElectronDataError::UnsupportedVersion {
            found: Some(998),
            ..
        })
    ));
    // Unknown keys, a bad label, blank provenance, a missing table.
    let unknown = format!("fermi_shift_ev = 1.0\n{text}");
    assert!(matches!(
        ShellChannelTable::from_toml_str(&unknown),
        Err(ElectronDataError::Parse(_))
    ));
    let label = text.replacen("subshell = \"L3\"", "subshell = \"L9\"", 1);
    assert!(matches!(
        ShellChannelTable::from_toml_str(&label),
        Err(ElectronDataError::MalformedText {
            field: "subshell",
            ..
        })
    ));
    let blank = text.replacen(
        "binding_provenance = \"synthetic binding energy, not physical data\"",
        "binding_provenance = \"  \"",
        1,
    );
    assert_eq!(
        ShellChannelTable::from_toml_str(&blank),
        Err(ElectronDataError::MissingProvenance)
    );
    assert!(matches!(
        ShellChannelTable::from_toml_str(head),
        Err(ElectronDataError::Parse(_))
    ));
}

#[test]
fn shell_channel_table_rejects_invalid_metadata() {
    let l3 = Subshell::from_label("L3").unwrap();
    let ok = || shell_loss_table(SamplingAxis::InelasticEnergyLoss, 100.0);
    assert!(ShellChannelTable::new(0, l3, 100.0, SYNTHETIC, ok()).is_err());
    for b in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ShellChannelTable::new(14, l3, b, SYNTHETIC, ok()).is_err(),
            "{b}"
        );
    }
    assert_eq!(
        ShellChannelTable::new(14, l3, 100.0, "", ok()),
        Err(ElectronDataError::MissingProvenance)
    );
    // A loss below the binding energy, or a binding energy above some loss.
    let low = shell_loss_table(SamplingAxis::InelasticEnergyLoss, 99.0);
    assert!(ShellChannelTable::new(14, l3, 100.0, SYNTHETIC, low).is_err());
    assert!(ShellChannelTable::new(14, l3, 100.5, SYNTHETIC, ok()).is_err());
    // An elastic table is not a loss channel.
    let elastic = CrossSectionTable::new(CrossSectionTableParts {
        axis: SamplingAxis::ElasticPolarAngle,
        quantiles: vec![
            vec![],
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 2.0],
        ],
        ..ok().parts().clone()
    })
    .unwrap();
    assert!(ShellChannelTable::new(14, l3, 0.5, SYNTHETIC, elastic).is_err());
}
