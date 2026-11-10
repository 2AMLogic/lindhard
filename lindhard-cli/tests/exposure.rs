//! `lindhard exposure`: files written, PSF read from a parameters CSV, and
//! the dose independent of the thread count.

use lindhard_cli::exposure::{
    psf_from_parameters_csv, run, ExposureArgs, CONTOURS_FILE, DOSE_FILE, FEATURES_FILE,
    SUMMARY_FILE,
};

const PARAMS: &str = "model,parameter,value,std_error\n\
double_gaussian,total_ev,1e4,\ndouble_gaussian,alpha_m,4e-8,\n\
double_gaussian,beta_m,1.2e-7,\ndouble_gaussian,eta,0.7,\n\
double_gaussian,chi2,3,\n";

fn args(dir: &std::path::Path, threads: Option<usize>) -> ExposureArgs {
    ExposureArgs {
        shots: dir.join("shots.csv"),
        psf_parameters: Some(dir.join("params.csv")),
        model: None,
        total_ev: None,
        alpha_m: None,
        beta_m: None,
        eta: None,
        gamma_m: None,
        nu: None,
        blur_m: 1e-8,
        subsamples: 8,
        x0_m: -2e-7,
        y0_m: -2e-7,
        dx_m: 2e-8,
        dy_m: 2e-8,
        nx: 30,
        ny: 30,
        threshold: None,
        threshold_fraction: Some(0.5),
        direct: None,
        direct_histories: 80_000,
        direct_replicas: 8,
        direct_seed: 1,
        slab_lo_m: 0.0,
        slab_hi_m: None,
        tolerance_sigmas: 5.0,
        conv_rel_err: 0.0,
        floor_fraction: 0.01,
        out: dir.join(format!("out{}", threads.unwrap_or(0))),
        threads,
    }
}

#[test]
fn writes_dose_contours_and_summary_identically_at_any_thread_count() {
    let dir = std::env::temp_dir().join(format!("lindhard-exposure-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("params.csv"), PARAMS).unwrap();
    std::fs::write(
        dir.join("shots.csv"),
        "x_m,y_m,weight,shape,size_m\n0,0,1,point,\n5e-8,0,2,disc,2e-8\n",
    )
    .unwrap();
    let mut dose = Vec::new();
    for t in [1, 2, 4] {
        let a = args(&dir, Some(t));
        run(&a).unwrap();
        for f in [DOSE_FILE, CONTOURS_FILE, FEATURES_FILE, SUMMARY_FILE] {
            assert!(a.out.join(f).exists(), "{f}");
        }
        dose.push(std::fs::read(a.out.join(DOSE_FILE)).unwrap());
    }
    assert!(dose[0] == dose[1] && dose[1] == dose[2]);
    assert_eq!(
        String::from_utf8_lossy(&dose[0]).lines().count(),
        1 + 30 * 30
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn psf_parameters_csv_picks_the_requested_fit() {
    let psf = psf_from_parameters_csv(PARAMS, None).unwrap();
    assert_eq!((psf.total_ev, psf.alpha_m, psf.eta), (1e4, 4e-8, 0.7));
    assert!(psf.gamma_m.is_none());
    assert!(psf_from_parameters_csv(PARAMS, Some("triple_gaussian")).is_err());
    assert!(psf_from_parameters_csv(PARAMS, Some("nope")).is_err());
}

/// `--direct`: the transport dose and the comparison are written, listed in
/// the summary, and identical for 1 and 4 threads.
#[test]
fn direct_mode_writes_dose_and_comparison_identically_at_any_thread_count() {
    use lindhard_cli::exposure::{COMPARISON_FILE, DIRECT_FILE};

    let ex = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../examples/electron");
    let dir = std::env::temp_dir().join(format!("lindhard-exposure-direct-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(
        ex.join("synthetic_plasmon_elf.toml"),
        dir.join("synthetic_plasmon_elf.toml"),
    )
    .unwrap();
    let toml = std::fs::read_to_string(ex.join("e_10keV_si.toml"))
        .unwrap()
        .replace("points_per_decade = 8.0", "points_per_decade = 3.0");
    std::fs::write(dir.join("target.toml"), toml).unwrap();
    std::fs::write(dir.join("params.csv"), PARAMS).unwrap();
    std::fs::write(dir.join("shots.csv"), "x_m,y_m,weight\n0,0,1\n5e-8,0,2\n").unwrap();

    let mut files = Vec::new();
    for threads in [1, 4] {
        let mut a = args(&dir, Some(threads));
        a.direct = Some(dir.join("target.toml"));
        a.direct_histories = 80;
        a.direct_replicas = 2;
        a.slab_hi_m = Some(5e-8);
        a.threshold_fraction = None;
        run(&a).unwrap();
        let read = |n: &str| std::fs::read_to_string(a.out.join(n)).unwrap();
        let direct = read(DIRECT_FILE);
        assert!(direct.starts_with("x_m,y_m,dose_ev_per_m2,std_err_ev_per_m2\n"));
        assert_eq!(direct.lines().count(), 1 + 30 * 30);
        let cmp = read(COMPARISON_FILE);
        assert!(cmp.starts_with("x_m,y_m,convolution_ev_per_m2,direct_ev_per_m2,"));
        let summary: serde_json::Value = serde_json::from_str(&read(SUMMARY_FILE)).unwrap();
        assert_eq!(summary["direct"]["replicas"], 2);
        assert_eq!(summary["direct"]["points"], 900);
        files.push((direct, cmp));
    }
    assert_eq!(files[0], files[1]);

    // Without --slab-hi-m the mode is refused.
    let mut a = args(&dir, Some(1));
    a.direct = Some(dir.join("target.toml"));
    assert!(run(&a).is_err());
}
