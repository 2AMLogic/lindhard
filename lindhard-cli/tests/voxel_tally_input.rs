//! `[tally.voxel]` through the shared input and driver: the per-element
//! voxel tally of a TOML run sums to the run's cascade defect totals, is
//! published as CSV and VTK next to the other outputs, is independent of the
//! thread count, and bad grids fail with the field name.

use lindhard::input::Input;
use lindhard::tally::{SpeciesDefects, VoxelReport};
use lindhard_cli::output;
use lindhard_cli::publish::publish_static_ion;
use lindhard_cli::sim::simulate;

const BASE: &str = r#"
[beam]
ion = "B"
energy_ev = 5000.0

[materials.SiO2]
density_g_cm3 = 2.2
elements = [
  { symbol = "Si", atom_fraction = 1.0 },
  { symbol = "O", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },
]

[target]
substrate = "Si"
[[target.layers]]
material = "SiO2"
thickness_nm = 20.0

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 200
seed = 9

[tally]
depth_bin_nm = 2.0
depth_bins = 50
per_ion = false
"#;

const VOXEL: &str = r#"
[tally.voxel]
origin_nm = [0.0, -10.0, -10.0]
spacing_nm = [4.0, 10.0, 10.0]
dims = [10, 2, 2]
"#;

fn resolve(text: &str) -> lindhard::input::Resolved {
    Input::from_toml_str(text).unwrap().resolve().unwrap()
}

#[test]
fn voxel_tally_matches_the_cascade_totals_and_is_published() {
    let r = resolve(&format!("{BASE}{VOXEL}"));
    let s1 = simulate(&r, Some(1)).unwrap();
    let report = s1.tally.voxel.as_ref().expect("voxel tally").report();
    assert_eq!(report.species_z, vec![5, 14, 8]);
    let c = s1.report.damage.cascade;
    let mut sum = SpeciesDefects::default();
    for d in report.grand_defect_totals() {
        sum.add(&d);
    }
    assert_eq!(sum.vacancies, c.vacancies);
    assert_eq!(sum.interstitials, c.interstitials);
    assert_eq!(sum.replacements, c.replacements);
    assert!(c.displacements > 0);
    // Electronic loss reaches the tally through the CLI tally too.
    assert!(report.grand_totals().electronic_ev > 0.0);

    let s4 = simulate(&r, Some(4)).unwrap();
    assert_eq!(s4.tally.voxel.as_ref().unwrap().report(), report);

    // Publication: both files, and the summary names them.
    let dir = std::env::temp_dir().join(format!("lindhard-voxel-input-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let summary = "{}".to_string();
    let export = output::damage_export_csv(&r, &s1.tally).expect("export text");
    publish_static_ion(&dir, summary, &s1.tally, &s1.report, Some(export.clone())).unwrap();
    let on_disk = std::fs::read_to_string(dir.join(output::DAMAGE_EXPORT_FILE)).unwrap();
    assert_eq!(on_disk, export);
    let (back, prov) = VoxelReport::read_export_csv(&on_disk).unwrap();
    assert_eq!(back, report);
    assert_eq!((prov.engine.as_str(), prov.seed), ("ion", r.config.seed));
    assert_eq!(prov.input_sha256.len(), 64);
    // The export is a function of the input: other thread counts give the same bytes.
    assert_eq!(output::damage_export_csv(&r, &s4.tally).unwrap(), export);
    let csv = std::fs::read_to_string(dir.join(output::VOXEL_CSV_FILE)).unwrap();
    assert_eq!(VoxelReport::read_csv(&csv).unwrap(), report);
    VoxelReport::read_vtk(&std::fs::read_to_string(dir.join(output::VOXEL_VTK_FILE)).unwrap())
        .unwrap();

    // A rerun without the table removes the stale files.
    let plain = simulate(&resolve(BASE), Some(1)).unwrap();
    assert!(plain.tally.voxel.is_none());
    publish_static_ion(&dir, "{}".into(), &plain.tally, &plain.report, None).unwrap();
    assert!(!dir.join(output::VOXEL_CSV_FILE).exists());
    assert!(!dir.join(output::VOXEL_VTK_FILE).exists());
    assert!(!dir.join(output::DAMAGE_EXPORT_FILE).exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn bad_voxel_grids_name_the_field() {
    for (table, field) in [
        (
            "dims = [0, 1, 1]\nspacing_nm = [1.0, 1.0, 1.0]",
            "tally.voxel.dims",
        ),
        (
            "dims = [1, 1, 1]\nspacing_nm = [1.0, 0.0, 1.0]",
            "tally.voxel.spacing_nm",
        ),
    ] {
        let text = format!("{BASE}\n[tally.voxel]\norigin_nm = [0.0, 0.0, 0.0]\n{table}\n");
        let e = Input::from_toml_str(&text).unwrap().resolve().unwrap_err();
        assert!(e.to_string().contains(field), "{e}");
    }
}
