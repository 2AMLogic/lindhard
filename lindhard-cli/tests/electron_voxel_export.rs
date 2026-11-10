//! `[electron.tally.voxel]` through the electron driver: the per-voxel
//! deposit plus `outside` is the energy budget's `deposited_ev`, the export
//! carries `engine = electron`, it is the same bytes for 1 and 4 threads, and
//! it is absent without the table.

use std::path::PathBuf;

use lindhard::input::electron::ElectronInput;
use lindhard::tally::VoxelReport;
use lindhard_cli::electron::{damage_export_csv, simulate_electron};

const VOXEL: &str = "\n[electron.tally.voxel]\norigin_nm = [0.0, -500.0, -500.0]\n\
spacing_nm = [300.0, 500.0, 500.0]\ndims = [5, 2, 2]\n";

fn setup(name: &str) -> (PathBuf, String) {
    let ex = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../examples/electron");
    let dir = std::env::temp_dir().join(format!("lindhard-evoxel-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(
        ex.join("synthetic_plasmon_elf.toml"),
        dir.join("synthetic_plasmon_elf.toml"),
    )
    .unwrap();
    let text = std::fs::read_to_string(ex.join("e_10keV_si.toml"))
        .unwrap()
        .replace("points_per_decade = 8.0", "points_per_decade = 3.0")
        .replace("histories = 1000", "histories = 30");
    (dir, text)
}

#[test]
fn electron_voxel_export_balances_is_deterministic_and_optional() {
    let (dir, text) = setup("a");
    let with = format!("{text}{VOXEL}");
    let r = ElectronInput::from_toml_str(&with)
        .unwrap()
        .resolve_in(&dir)
        .unwrap();
    let s1 = simulate_electron(&r, Some(1), None).unwrap();
    let s4 = simulate_electron(&r, Some(4), None).unwrap();
    let csv1 = damage_export_csv(&r, &s1).expect("export text");
    assert_eq!(damage_export_csv(&r, &s4).unwrap(), csv1);

    let (back, prov) = VoxelReport::read_export_csv(&csv1).unwrap();
    assert_eq!(prov.engine, "electron");
    assert_eq!(back.histories, 30);
    let total = back.grand_totals().deposited_ev();
    let budget = s1.report.budget.deposited_ev;
    assert!(budget > 0.0);
    assert!(
        (total - budget).abs() <= 1e-9 * budget,
        "voxels + outside {total} vs budget {budget}"
    );

    // A grid that misses the interaction volume puts everything outside.
    let far = with.replace(
        "origin_nm = [0.0, -500.0, -500.0]",
        "origin_nm = [1e6, 1e6, 1e6]",
    );
    let rf = ElectronInput::from_toml_str(&far)
        .unwrap()
        .resolve_in(&dir)
        .unwrap();
    let sf = simulate_electron(&rf, Some(1), None).unwrap();
    let (miss, _) = VoxelReport::read_export_csv(&damage_export_csv(&rf, &sf).unwrap()).unwrap();
    assert_eq!(miss.totals().deposited_ev(), 0.0);
    let budget = sf.report.budget.deposited_ev;
    assert!((miss.outside.deposited_ev() - budget).abs() <= 1e-9 * budget);

    // Without the table: no report, no export.
    let rn = ElectronInput::from_toml_str(&text)
        .unwrap()
        .resolve_in(&dir)
        .unwrap();
    let sn = simulate_electron(&rn, Some(1), None).unwrap();
    assert!(sn.voxel.is_none());
    assert!(damage_export_csv(&rn, &sn).is_none());
}
