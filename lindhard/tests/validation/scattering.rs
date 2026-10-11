//! Scattering-integral checks (issue #3): nuclear stopping from the
//! quadrature against the ZBL universal fit and the exact high-energy
//! asymptote, the magic formula against the quadrature, and the Rutherford
//! limit.
//!
//! The reference functions are copies of the test-only helpers in
//! `lindhard/src/ion/scattering.rs` (same citations, same digits; see the
//! "ZBL universal nuclear stopping fit" row of `docs/data-provenance.md`).
//! They are copied, not shared, so the library keeps no test-only API.

use lindhard::ion::potential::Screening;
use lindhard::ion::scattering::{
    nuclear_stopping_reduced, theta_magic, theta_quadrature, MagicConstants, ScatteringTable,
};

use crate::report::{pct, sci, spct, Check};

/// ZBL universal reduced nuclear stopping fit, Ziegler, Biersack, Littmark,
/// *The Stopping and Range of Ions in Solids* (Pergamon, 1985), Ch. 2:
/// `ln(1 + 1.1383 eps) / (2 (eps + 0.01321 eps^0.21226 + 0.19593 eps^0.5))`
/// for `eps <= 30`, `ln(eps) / (2 eps)` above. The `eps <= 30` coefficients
/// were cross-checked against the MIT-licensed `ir2-lab/screened_coulomb`
/// (commit f84c3c8) in PR #40.
fn zbl_fit(eps: f64) -> f64 {
    if eps <= 30.0 {
        (1.0 + 1.1383 * eps).ln()
            / (2.0 * (eps + 0.01321 * eps.powf(0.21226) + 0.19593 * eps.sqrt()))
    } else {
        eps.ln() / (2.0 * eps)
    }
}

/// High-energy asymptote of the reduced nuclear stopping for a
/// sum-of-exponentials screening function, `s_n -> (S^2 ln eps + C)/(2 eps)`,
/// from the first-order momentum (impulse) approximation (Lindhard, Nielsen,
/// Scharff, Mat. Fys. Medd. Dan. Vid. Selsk. 36(10) (1968), p. 12,
/// eqs. (3.1)-(3.4); verified in #45) matched to Rutherford scattering; `C`
/// from Lommel integrals (Watson, *Theory of Bessel Functions*, 2nd ed.,
/// 1944, Sec. 5.11). Derivation in the doc comment of the same helper in
/// `lindhard/src/ion/scattering.rs`.
fn high_energy_asymptote(s: Screening, eps: f64) -> f64 {
    let terms = s
        .exponential_terms()
        .expect("sum-of-exponentials screening");
    let sum: f64 = terms.iter().map(|&(c, _)| c).sum();
    let l = |a: f64, b: f64| {
        if a == b {
            (0.5 * a).ln() + 0.5
        } else {
            (a * a * (0.5 * a).ln() - b * b * (0.5 * b).ln()) / (a * a - b * b)
        }
    };
    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
    let mut c = sum * sum * (std::f64::consts::LN_2 - sum.ln() - EULER_GAMMA);
    for &(ci, bi) in terms {
        for &(cj, bj) in terms {
            c -= ci * cj * l(bi, bj);
        }
    }
    (sum * sum * eps.ln() + c) / (2.0 * eps)
}

fn beta_min(eps: f64) -> f64 {
    1e-4f64.min(1e-3 / eps)
}

fn sn(eps: f64) -> f64 {
    nuclear_stopping_reduced(Screening::ZblUniversal, eps, beta_min(eps), 1e3, 40)
}

pub(crate) fn checks(table: &ScatteringTable) -> Vec<Check> {
    let mut out = Vec::new();

    // 1. Low energy: quadrature s_n vs the ZBL fit, eps in [1e-4, 0.3].
    let mut worst = (0.0f64, 0.0);
    for k in 0..=7 {
        let eps = 1e-4 * 10f64.powf(k as f64 / 2.0);
        let d = (sn(eps) / zbl_fit(eps) - 1.0).abs();
        if d > worst.0 {
            worst = (d, eps);
        }
    }
    out.push(Check::at_most(
        "scatter.sn_vs_zbl_fit.low",
        format!(
            "ZBL nuclear stopping, quadrature vs ZBL universal fit, 1e-4 <= eps <= 0.3 (worst at eps = {:.1e})",
            worst.1
        ),
        worst.0,
        5e-3,
        pct,
        "ZBL 1985 Ch. 2; the fit's own accuracy is the limit",
    ));

    // 2. Middle: reported, not asserted. Known deviation of the fit.
    let mut worst = (0.0f64, 0.0);
    for k in 0..=7 {
        let eps = 10f64.powf(-0.5 + k as f64 / 2.0);
        let d = sn(eps) / zbl_fit(eps) - 1.0;
        if d.abs() > worst.0.abs() {
            worst = (d, eps);
        }
    }
    out.push(Check::info(
        "scatter.sn_vs_zbl_fit.mid",
        format!(
            "ZBL nuclear stopping, quadrature vs ZBL universal fit, 0.3 <= eps <= 1e3 (worst at eps = {:.0})",
            worst.1
        ),
        spct(worst.0),
        "known: the fit's eps > 30 branch omits the screening constant C (it runs low by C / ln eps) and its eps <= 30 branch joins it; not a quadrature error (#3 amended AC, PR #40 review). Anchored instead by the asymptote row and the convergence test in scattering.rs",
    ));

    // 3. High energy: quadrature vs the exact asymptote, three screenings.
    let mut worst = 0.0f64;
    for s in [Screening::ZblUniversal, Screening::KrC, Screening::Moliere] {
        for &eps in &[1e3, 1e4, 1e5, 1e6] {
            let q = nuclear_stopping_reduced(s, eps, beta_min(eps), 1e3, 40);
            worst = worst.max((q / high_energy_asymptote(s, eps) - 1.0).abs());
        }
    }
    out.push(Check::at_most(
        "scatter.sn_vs_asymptote.high",
        "Nuclear stopping, quadrature vs impulse-approximation asymptote (S^2 ln eps + C)/(2 eps), 1e3 <= eps <= 1e6, ZBL/Kr-C/Moliere",
        worst,
        2e-3,
        pct,
        "the asymptote has its own O(1/eps) error, largest at eps = 1e3",
    ));

    // 4. Magic formula vs quadrature, worst |d cos(theta/2)|.
    for (id, s, k, tol, pubacc) in [
        (
            "scatter.magic_vs_quadrature.zbl",
            Screening::ZblUniversal,
            MagicConstants::ZBL,
            2e-2,
            "ZBL",
        ),
        (
            "scatter.magic_vs_quadrature.moliere",
            Screening::Moliere,
            MagicConstants::MOLIERE,
            1e-2,
            "Moliere",
        ),
    ] {
        let mut worst = 0.0f64;
        for &eps in &[1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3] {
            for j in 0..50 {
                let beta = 1e-3 * 10f64.powf(j as f64 / 10.0);
                let q = theta_quadrature(s, eps, beta);
                let m = theta_magic(s, &k, eps, beta);
                worst = worst.max(((0.5 * m).cos() - (0.5 * q).cos()).abs());
            }
        }
        out.push(Check::at_most(
            id,
            format!(
                "Magic formula ({pubacc} constants) vs quadrature, worst |d cos(theta/2)|, 1e-3 <= eps <= 1e3, 1e-3 <= beta <= 1e2"
            ),
            worst,
            tol,
            sci,
            "Biersack and Haggmark, NIM 174 (1980) 257; tolerance is the fit accuracy established in PR #40",
        ));
    }

    // 5. Rutherford limit: at eps beta >> 1 inside the screening radius,
    // tan(theta/2) = 1 / (2 eps beta).
    let mut worst = 0.0f64;
    for s in Screening::ALL {
        for &(eps, beta) in &[(1e6, 1e-3), (1e5, 1e-3), (1e6, 1e-4), (1e7, 1e-4)] {
            let th = theta_quadrature(s, eps, beta);
            let expect = 2.0 * (1.0 / (2.0 * eps * beta)).atan();
            worst = worst.max((th / expect - 1.0).abs());
        }
    }
    out.push(Check::at_most(
        "scatter.rutherford_limit",
        "Quadrature angle vs Rutherford tan(theta/2) = 1/(2 eps beta), beta <= 1e-3, eps >= 1e5, all four screenings",
        worst,
        2e-3,
        pct,
        "residual is screening: phi(x) differs from 1 by O(x) at x ~ beta",
    ));

    // 6. The engine's scattering table: measured interpolation error.
    out.push(Check::at_most(
        "scatter.table_interpolation",
        "Engine scattering table (CLI grid, 32/decade) vs quadrature, worst |d theta| at sampled cell centres, rad",
        table.max_abs_error(),
        1e-2,
        sci,
        "measured at build time; bound from the table tests of PR #40",
    ));
    out
}
