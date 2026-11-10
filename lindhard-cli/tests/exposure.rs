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
