//! Masked implant against the analytic edge profile (issue #435): a beam
//! that illuminates the half-space `y > 0` uniformly, behind an ideal mask
//! that stops everything at `y < 0`, implanted into a voxel target and
//! tallied with the sparse `VoxelTally` (issue #433). The lateral dose
//! profile at the mask edge is the superposition of the pencil beams that
//! enter through the open aperture. When one pencil beam's lateral
//! distribution is a Gaussian of standard deviation `s`, the sum over a
//! half-plane is `n(y) = n0 * Phi(y / s)` with `Phi` the standard normal
//! distribution function (an error function), and for an aperture of finite
//! width `a` it is `n0 * (Phi(y / s) - Phi((y - a) / s))`. This is the
//! convolution of a step with a Gaussian, derived here, not taken from a
//! table. `s` is the lateral straggle of the corresponding 1D profile.
//!
//! The mask is realised by running one pencil beam per aperture voxel, entered
//! at the voxel centre (`Bca::with_entry_point`), with its own seed, and
//! merging the voxel tallies. The reference straggle is the standard deviation
//! of the lateral position `y` of the stopped primaries of a separate
//! `IonTally` run on a semi-infinite layer stack (the layer engine, a
//! different geometry code from the voxel DDA), at a different seed. Only
//! the shape is compared: the amplitude is fixed by the number of implanted
//! atoms, which the tally conserves.
//!
//! A pencil beam's lateral distribution is not exactly Gaussian (its kurtosis
//! is reported), so the error function is an approximation and the rows say
//! how good it is.

use lindhard::geometry::{Boundary, Stack, VoxelGrid};
use lindhard::ion::bca::{Bca, BcaConfig, BcaTally, Beam};
use lindhard::ion::scattering::ScatteringTable;
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonTally, IonTallyConfig, VoxelTally, VoxelTallyGrid};

use crate::report::{num, pct, Check};

const NM: f64 = 1e-9;

/// Beam: 2 keV boron at normal incidence into silicon. Test parameters, not
/// data: `E_d` = 15 eV, `E_s` = 2 eV (as the other level-1 engine checks).
const ION_Z: u8 = 5;
const ENERGY_EV: f64 = 2.0e3;

fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    m
}

/// Standard normal distribution function, from the complementary error
/// function of Abramowitz and Stegun, *Handbook of Mathematical Functions*
/// (1964), eq. 7.1.26 (absolute error below 1.5e-7).
fn big_phi(t: f64) -> f64 {
    let x = t.abs() / std::f64::consts::SQRT_2;
    let u = 1.0 / (1.0 + 0.327_591_1 * x);
    let poly = u
        * (0.254_829_592
            + u * (-0.284_496_736
                + u * (1.421_413_741 + u * (-1.453_152_027 + u * 1.061_405_429))));
    let erfc = poly * (-x * x).exp();
    // Phi(t) = 1 - erfc(t / sqrt 2) / 2 for t >= 0, erfc(|t| / sqrt 2) / 2 below.
    if t >= 0.0 {
        1.0 - 0.5 * erfc
    } else {
        0.5 * erfc
    }
}

/// `G(t) = integral of Phi(u / s) du` up to a constant: `t Phi(t / s) +
/// s phi(t / s)`, with `phi` the standard normal density.
fn big_g(t: f64, s: f64) -> f64 {
    let z = t / s;
    t * big_phi(z) + s * (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// Expected fraction of the implanted atoms in the lateral bin `[y0, y1]`
/// for an aperture `[0, a]` and Gaussian pencil-beam straggle `s`: the bin
/// integral of `(Phi(y / s) - Phi((y - a) / s)) / a`.
fn expected_fraction(y0: f64, y1: f64, a: f64, s: f64) -> f64 {
    ((big_g(y1, s) - big_g(y0, s)) - (big_g(y1 - a, s) - big_g(y0 - a, s))) / a
}

pub(crate) fn checks(table: &ScatteringTable, quick: bool) -> Vec<Check> {
    let ls = LindhardScharff::new();
    let mut out = Vec::new();
    let ion = Ion::new(ION_Z).unwrap();

    // Reference: lateral straggle of the 1D (layer-engine) profile.
    let ref_count: u64 = if quick { 20_000 } else { 200_000 };
    let stack = Stack::semi_infinite(si());
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 0x1D_435;
    let bca = Bca::new(
        Beam::normal(ion, ENERGY_EV, ref_count),
        &stack,
        cfg,
        &ls,
        table,
    )
    .unwrap();
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 40.0 * NM, 40).unwrap(),
        lateral: Binning::new(-40.0 * NM, 40.0 * NM, 80).unwrap(),
        radial: Binning::new(0.0, 40.0 * NM, 40).unwrap(),
        escape_energy: Binning::new(0.0, ENERGY_EV, 10).unwrap(),
        escape_polar: Binning::new(0.0, 0.5 * std::f64::consts::PI, 9).unwrap(),
    };
    let proto = IonTally::new(&stack, &bca.species_z(), tc).unwrap();
    let rep = bca.run(|| proto.clone()).unwrap().report(false);
    let ly = rep.range.lateral_y.expect("lateral moments of the 1D run");
    let (s_ref, s_ref_err) = (ly.std_dev, ly.std_dev_std_err);

    // Masked implant on a voxel target. Depth: 40 x 1 nm (the range is a few
    // nm); lateral y: voxels of `h`, the mask edge at a voxel face, the
    // aperture `a` wide and margins of at least 5 s either side so the grid
    // sides do not clip the profile; z: one periodic voxel (only the y
    // marginal is tallied, and wrapping only relabels z).
    let h = 1.0 * NM;
    let margin = (6.0 * s_ref / h).ceil() * h;
    let a = (40.0 * NM / h).round() * h;
    let n_lo = (margin / h).round() as usize;
    let n_ap = (a / h).round() as usize;
    let ny = n_lo + n_ap + n_lo;
    let nx = 40;
    let grid = VoxelGrid::new(
        vec![si()],
        [nx, ny, 1],
        [NM, h, 20.0 * NM],
        vec![0; nx * ny],
        [Boundary::Vacuum, Boundary::Vacuum, Boundary::Periodic],
    )
    .unwrap();
    let tgrid = VoxelTallyGrid::from_voxel_grid(&grid);
    let y_origin = -margin; // y of the grid's low face; the mask edge is at y = 0
    let per_entry: u64 = if quick { 600 } else { 6_000 };
    let mut tally = VoxelTally::new(tgrid);
    for j in 0..n_ap {
        let mut cfg = BcaConfig::new(5.0, 2.0);
        cfg.seed = 0x435_0000 + j as u64;
        let yc = (n_lo + j) as f64 * h + 0.5 * h;
        let bca = Bca::new(
            Beam::normal(ion, ENERGY_EV, per_entry),
            &grid,
            cfg,
            &ls,
            table,
        )
        .unwrap()
        .with_entry_point([0.0, yc, 0.5 * 20.0 * NM])
        .unwrap();
        let t = bca.run(|| VoxelTally::new(tgrid)).unwrap();
        tally.merge(t);
    }
    let report = tally.report();

    // Implanted atoms per lateral voxel row.
    let mut rows = vec![0u64; ny];
    for v in &report.voxels {
        rows[v.index[1]] += v.cell.implanted;
    }
    let n_in: u64 = rows.iter().sum();
    let n_out = report.outside.implanted;
    let n_total = n_in + n_out;
    let histories = per_entry * n_ap as u64;
    out.push(Check::holds(
        "mask.conservation",
        format!(
            "masked implant: primaries at rest in the grid and outside it ({n_out} outside) do not exceed the {histories} histories"
        ),
        n_total <= histories,
        format!("{n_total} <= {histories}"),
        "a history comes to rest at most once; the rest escaped by backscatter",
    ));

    // Fit s to the observed rows by minimising Pearson chi-square over the
    // voxels with at least 25 expected atoms.
    let chi2 = |s: f64| -> (f64, usize) {
        let (mut c, mut n) = (0.0, 0usize);
        for (j, &obs) in rows.iter().enumerate() {
            let y0 = y_origin + j as f64 * h;
            let e = n_in as f64 * expected_fraction(y0, y0 + h, a, s);
            if e >= 25.0 {
                c += (obs as f64 - e).powi(2) / e;
                n += 1;
            }
        }
        (c, n)
    };
    let (mut lo, mut hi) = (0.2 * s_ref, 3.0 * s_ref);
    for _ in 0..80 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if chi2(m1).0 < chi2(m2).0 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let s_fit = 0.5 * (lo + hi);
    let (c_min, ndf_pts) = chi2(s_fit);
    // One-parameter error from the curvature of chi-square: delta s =
    // sqrt(2 / chi2'').
    let d = 0.01 * s_fit;
    let curv = (chi2(s_fit + d).0 - 2.0 * c_min + chi2(s_fit - d).0) / (d * d);
    let s_fit_err = (2.0 / curv).sqrt();
    let (c_ref, _) = chi2(s_ref);

    out.push(Check::info(
        "mask.diagnostics",
        format!(
            "B 2 keV into Si, {n_ap} pencil beams x {per_entry} ions, mask edge at y = 0; 1D lateral straggle s_ref = {} nm, kurtosis {}",
            num(s_ref / NM),
            num(ly.kurtosis)
        ),
        format!(
            "chi2/pts: fit {} / {ndf_pts}, s_ref {} / {ndf_pts}",
            num(c_min),
            num(c_ref)
        ),
        "chi-square of the row counts against the erf profile",
    ));
    // Measured straggle against the 1D straggle: a systematic allowance of
    // 1 % plus three combined standard errors.
    let comb = (s_fit_err.powi(2) + s_ref_err.powi(2)).sqrt();
    out.push(Check::at_most(
        "mask.straggle",
        "masked implant: lateral straggle fitted to the edge profile vs the 1D lateral straggle, |s_fit / s_ref - 1|",
        (s_fit / s_ref - 1.0).abs(),
        0.01 + 3.0 * comb / s_ref,
        pct,
        format!(
            "s_fit = {} +- {} nm, s_ref = {} +- {} nm (1 % + 3 standard errors)",
            num(s_fit / NM),
            num(s_fit_err / NM),
            num(s_ref / NM),
            num(s_ref_err / NM)
        ),
    ));

    // Mask position: shift the whole profile by d and minimise chi-square.
    let chi2_shift = |d: f64| -> f64 {
        let mut c = 0.0;
        for (j, &obs) in rows.iter().enumerate() {
            let y0 = y_origin + j as f64 * h - d;
            let e = n_in as f64 * expected_fraction(y0, y0 + h, a, s_fit);
            if e >= 25.0 {
                c += (obs as f64 - e).powi(2) / e;
            }
        }
        c
    };
    let (mut lo, mut hi) = (-2.0 * s_ref, 2.0 * s_ref);
    for _ in 0..80 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if chi2_shift(m1) < chi2_shift(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let d_fit = 0.5 * (lo + hi);
    let dd = 0.01 * s_ref;
    let curv =
        (chi2_shift(d_fit + dd) - 2.0 * chi2_shift(d_fit) + chi2_shift(d_fit - dd)) / (dd * dd);
    let d_err = (2.0 / curv).sqrt();
    out.push(Check::at_most(
        "mask.edge_position",
        "masked implant: fitted shift of the edge profile from the mask position y = 0, in standard errors",
        d_fit.abs() / d_err,
        3.0,
        num,
        format!(
            "shift {} +- {} nm; an ideal edge sits at the half-dose point of an even pencil-beam distribution",
            num(d_fit / NM),
            num(d_err / NM)
        ),
    ));

    // Distribution-level shape: largest gap between the measured and the
    // erf-predicted cumulative lateral distribution (no fitted parameter, s =
    // s_ref), as a fraction of the implanted atoms.
    let (mut cum_obs, mut cum_exp, mut dmax) = (0.0, 0.0, 0.0f64);
    for (j, &obs) in rows.iter().enumerate() {
        let y0 = y_origin + j as f64 * h;
        cum_obs += obs as f64 / n_in as f64;
        cum_exp += expected_fraction(y0, y0 + h, a, s_ref);
        dmax = dmax.max((cum_obs - cum_exp).abs());
    }
    out.push(Check::at_most(
        "mask.shape",
        "masked implant: largest gap between the measured and the erf-predicted cumulative lateral dose (s = s_ref, no fit), fraction of the implanted atoms",
        dmax,
        0.005 + 2.0 / (n_in as f64).sqrt(),
        num,
        "0.005 + 2 / sqrt(N); the pencil-beam distribution is only approximately Gaussian (kurtosis in mask.diagnostics)",
    ));
    out
}
