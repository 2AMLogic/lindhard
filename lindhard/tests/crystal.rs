//! Acceptance tests for the cubic lattice model and the wafer/beam
//! orientation (`lindhard::ion::crystal`, issue #19).

use std::f64::consts::{FRAC_1_SQRT_2, PI};

use lindhard::ion::bca::Beam;
use lindhard::ion::crystal::lattice::GAAS_DENSITY_G_CM3;
use lindhard::ion::crystal::{Divergence, Lattice, Orientation};
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::rng::{run_particles, stream};
use lindhard::units::g_cm3_to_kg_m3;

fn deg(x: f64) -> f64 {
    x.to_radians()
}

fn assert_close(got: [f64; 3], want: [f64; 3], tol: f64) {
    for i in 0..3 {
        assert!(
            (got[i] - want[i]).abs() <= tol,
            "component {i}: got {got:?}, want {want:?}"
        );
    }
}

/// At zero tilt the beam travels along the inward surface normal, so a wafer
/// cut on (100), (110) or (111) sends it straight down the <100>, <110> or
/// <111> axial channel: `[100]`, `[110]/sqrt 2`, `[111]/sqrt 3` in the
/// crystal frame. Twist and wafer rotation must not matter at zero tilt.
#[test]
fn zero_tilt_gives_the_low_index_channels() {
    let si = Lattice::silicon();
    let s3 = 1.0 / 3.0f64.sqrt();
    let cases = [
        ([1, 0, 0], [0, 1, 0], [1.0, 0.0, 0.0]),
        ([1, 1, 0], [0, 0, 1], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0]),
        ([1, 1, 1], [1, -1, 0], [s3, s3, s3]),
    ];
    for (hkl, uvw, want) in cases {
        for (twist, rot) in [(0.0, 0.0), (deg(22.0), deg(-40.0)), (1.3, 2.9)] {
            let o = Orientation::new(&si, hkl, uvw, 0.0, twist, rot).unwrap();
            assert_close(o.beam_crystal(), want, 1e-12);
        }
    }
}

/// Tilt 7°, twist 22° on (100) Si.
///
/// Derivation. With surface `(hkl) = (100)` the inward normal is
/// `n = [100]`. Take the in-plane reference `r = [010]` (zone law:
/// `1·0 + 0·1 + 0·0 = 0`); then `t = n x r = [001]`. The lab beam is
/// `(cos θ, sin θ cos φ, sin θ sin φ)` and the crystal-frame beam is
/// `d = cos θ n + sin θ cos(φ - ω) r + sin θ sin(φ - ω) t`, so with
/// `θ = 7°`, `φ = 22°`, `ω = 0`:
///
/// ```text
/// d = (cos 7°, sin 7° cos 22°, sin 7° sin 22°)
///   = (0.992546151641322, 0.11299528757190813, 0.04565305957483624).
/// ```
///
/// With the reference along a <110> instead, `r = [011]/√2` (zone law
/// `0 + 0 + 0 = 0`), `t = n x r = [0,-1,1]/√2`, and
///
/// ```text
/// d = (cos 7°, sin 7° (cos 22° - sin 22°)/√2, sin 7° (cos 22° + sin 22°)/√2)
///   = (0.992546151641322, 0.047618146076940114, 0.11218132209150039).
/// ```
///
/// The decimal values were evaluated separately (Python `math`) from the
/// closed forms, so the test checks the matrix composition, not a copy of
/// it. A wafer rotation `ω` enters only as `φ - ω`, checked last.
#[test]
fn tilt_7_twist_22_on_100_silicon() {
    let si = Lattice::silicon();
    let (th, ph) = (deg(7.0), deg(22.0));

    let o = Orientation::new(&si, [1, 0, 0], [0, 1, 0], th, ph, 0.0).unwrap();
    let closed = [th.cos(), th.sin() * ph.cos(), th.sin() * ph.sin()];
    assert_close(o.beam_crystal(), closed, 1e-12);
    assert_close(
        o.beam_crystal(),
        [0.992546151641322, 0.11299528757190813, 0.04565305957483624],
        1e-12,
    );

    let o = Orientation::new(&si, [1, 0, 0], [0, 1, 1], th, ph, 0.0).unwrap();
    let closed = [
        th.cos(),
        th.sin() * (ph.cos() - ph.sin()) * FRAC_1_SQRT_2,
        th.sin() * (ph.cos() + ph.sin()) * FRAC_1_SQRT_2,
    ];
    assert_close(o.beam_crystal(), closed, 1e-12);
    assert_close(
        o.beam_crystal(),
        [0.992546151641322, 0.047618146076940114, 0.11218132209150039],
        1e-12,
    );

    // Twist 22° with the wafer turned by 22° is twist 0 on an unturned wafer.
    let a = Orientation::new(&si, [1, 0, 0], [0, 1, 0], th, ph, ph).unwrap();
    let b = Orientation::new(&si, [1, 0, 0], [0, 1, 0], th, 0.0, 0.0).unwrap();
    assert_close(a.beam_crystal(), b.beam_crystal(), 1e-12);
}

/// The lab beam direction is the BCA engine's `Beam::direction` for the same
/// polar angle (tilt) and azimuth (twist), so one convention serves both.
#[test]
fn lab_beam_matches_the_bca_beam_convention() {
    let si = Lattice::silicon();
    let (th, ph) = (deg(7.0), deg(22.0));
    let o = Orientation::new(&si, [1, 0, 0], [0, 1, 1], th, ph, deg(45.0)).unwrap();
    let beam = Beam {
        ion: Ion::new(5).unwrap(),
        energy_ev: 1.0e3,
        polar_rad: th,
        azimuth_rad: ph,
        count: 1,
    };
    assert_eq!(o.beam_lab(), beam.direction());
}

/// The channel example of Nordlund, Djurabekova and Hobler, Phys. Rev. B 94,
/// 214109 (2016), Sec. II B: on a [001] surface, `θ = 45°, ϕ = 0°` is
/// `[011]` and `θ = 54.73°, ϕ = 45°` is a `<111>` axis. Our convention with
/// `r = [010]` gives `[011]/√2` and `[-111]/√3` (exactly, with
/// `θ = atan √2`).
#[test]
fn nordlund_channel_example() {
    let si = Lattice::silicon();
    let o = Orientation::new(&si, [0, 0, 1], [0, 1, 0], deg(45.0), 0.0, 0.0).unwrap();
    assert_close(o.beam_crystal(), [0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2], 1e-12);
    let th = 2.0f64.sqrt().atan();
    let o = Orientation::new(&si, [0, 0, 1], [0, 1, 0], th, deg(45.0), 0.0).unwrap();
    let s3 = 1.0 / 3.0f64.sqrt();
    assert_close(o.beam_crystal(), [-s3, s3, s3], 1e-12);
}

/// Atom number density from the lattice, `8 / a³`, against the `material`
/// module's `ρ N_A / M̄`.
///
/// * Si: the default element density 2.329 g/cm³ (printed to 4 significant
///   figures, half a unit in the last place is 2.1e-4 relative) and the
///   CODATA 2022 lattice parameter (relative uncertainty 1.6e-8, so `8/a³` to
///   about 5e-8). The atomic weight interval [28.084, 28.086] adds 3.6e-5.
///   Tolerance: 2.1e-4 + 3.6e-5 + 5e-8, rounded up to 2.5e-4.
/// * GaAs: the density 5.321 g/cm³ and `a = 5.652 Å`, both printed on the
///   same NBS page (Monograph 25, Section 3, p. 33, at 25 °C). Half a unit in
///   the last place: 9.4e-5 for the density and 8.8e-5 for `a`, so 2.7e-4 for
///   `8/a³`. Tolerance: 9.4e-5 + 2.7e-4, rounded up to 3.7e-4.
#[test]
fn number_density_matches_the_material_module() {
    let si = Lattice::silicon();
    let m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    let rel = si.atom_number_density() / m.atom_number_density() - 1.0;
    assert!(rel.abs() < 2.5e-4, "Si relative difference {rel:e}");

    let gaas = Lattice::gallium_arsenide();
    let m = Material::from_atom_fractions(
        &[(31, 1.0), (33, 1.0)],
        Some(g_cm3_to_kg_m3(GAAS_DENSITY_G_CM3)),
    )
    .unwrap();
    let rel = gaas.atom_number_density() / m.atom_number_density() - 1.0;
    assert!(rel.abs() < 3.7e-4, "GaAs relative difference {rel:e}");
    for z in [31, 33] {
        let rel = gaas.number_density_of(z) / m.number_density_of(z).unwrap() - 1.0;
        assert!(rel.abs() < 3.7e-4, "GaAs Z={z} relative difference {rel:e}");
    }
}

/// One-sample Kolmogorov-Smirnov test.
///
/// Statistic as defined in the NIST/SEMATECH e-Handbook of Statistical
/// Methods, section 1.3.5.16 "Kolmogorov-Smirnov Goodness-of-Fit Test"
/// (<https://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm>):
/// `D = max_i max(F(Y_i) - (i-1)/N, i/N - F(Y_i))` over the ordered sample.
/// The p-value uses Kolmogorov's limiting distribution,
/// `P(√N D <= x) -> L(x) = 1 - 2 Σ_{i>=1} (-1)^(i-1) exp(-2 i² x²)`,
/// G. Marsaglia, W. W. Tsang and J. Wang, "Evaluating Kolmogorov's
/// Distribution", J. Stat. Softw. 8(18) (2003), doi:10.18637/jss.v008.i18,
/// section 3; they put the error of the limit near `0.278/√N`, about 0.002
/// at the sample sizes used here.
fn ks_p_value(mut sample: Vec<f64>, cdf: impl Fn(f64) -> f64) -> (f64, f64) {
    sample.sort_by(f64::total_cmp);
    let n = sample.len() as f64;
    let mut d = 0.0f64;
    for (i, &y) in sample.iter().enumerate() {
        let f = cdf(y);
        d = d.max(f - i as f64 / n).max((i + 1) as f64 / n - f);
    }
    let x = n.sqrt() * d;
    let mut tail = 0.0;
    for i in 1..=100 {
        let k = i as f64;
        let term = (-2.0 * k * k * x * x).exp();
        tail += if i % 2 == 1 { term } else { -term };
    }
    (d, (2.0 * tail).clamp(0.0, 1.0))
}

const KS_N: u64 = 20_000;

/// Polar deflections for particles `0..KS_N`, one per particle stream,
/// measured back from the sampled directions: the angle from the central
/// direction via `atan2(|c x d|, c . d)`, accurate at small angles. This
/// checks the full path (sampling and rotation), not only the angle draw.
fn sampled_polar_angles(div: Divergence, seed: u64) -> Vec<f64> {
    let central = [0.8, 0.36, 0.48]; // unit
    (0..KS_N)
        .map(|i| {
            let d = div.sample_direction(central, &mut stream(seed, i));
            let c = [
                central[1] * d[2] - central[2] * d[1],
                central[2] * d[0] - central[0] * d[2],
                central[0] * d[1] - central[1] * d[0],
            ];
            let s = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
            let co = central[0] * d[0] + central[1] * d[1] + central[2] * d[2];
            s.atan2(co)
        })
        .collect()
}

/// The sampled divergence matches its target distribution by a KS test at
/// the 1 % level: the polar deflection against the Rayleigh (Gaussian) or
/// cone CDF, and the azimuth against the uniform CDF on `[0, 2π)`. A
/// mis-scaled target (σ off by 5 %) must be rejected, so the test has power.
#[test]
fn divergence_passes_ks_at_one_percent() {
    let gauss = Divergence::Gaussian {
        sigma_rad: deg(0.5),
    };
    let cone = Divergence::UniformCone {
        half_angle_rad: deg(2.0),
    };
    for (div, seed) in [(gauss, 101), (cone, 202)] {
        let polar = sampled_polar_angles(div, seed);
        let (d, p) = ks_p_value(polar, |t| div.polar_cdf(t));
        assert!(p > 0.01, "{div:?}: polar KS D = {d}, p = {p}");

        let azimuth: Vec<f64> = (0..KS_N)
            .map(|i| div.sample_deflection(&mut stream(seed, i)).1)
            .collect();
        let (d, p) = ks_p_value(azimuth, |a| (a / (2.0 * PI)).clamp(0.0, 1.0));
        assert!(p > 0.01, "{div:?}: azimuth KS D = {d}, p = {p}");
    }

    let wrong = Divergence::Gaussian {
        sigma_rad: deg(0.5) * 1.05,
    };
    let (d, p) = ks_p_value(sampled_polar_angles(gauss, 101), |t| wrong.polar_cdf(t));
    assert!(p < 0.01, "wrong sigma not rejected: D = {d}, p = {p}");
}

/// Divergence draws come from the particle's own stream, so the sampled
/// directions are bit-identical on 1 and 4 threads.
#[test]
fn divergence_is_thread_count_independent() {
    let div = Divergence::Gaussian { sigma_rad: 0.01 };
    let si = Lattice::silicon();
    let o = Orientation::new(&si, [1, 0, 0], [0, 1, 1], deg(7.0), deg(22.0), 0.0).unwrap();
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            run_particles(
                77,
                2_000,
                64,
                Vec::new,
                |v: &mut Vec<[u64; 3]>, rng, _i| {
                    let d = o.to_crystal(div.sample_direction(o.beam_lab(), rng));
                    v.push(d.map(f64::to_bits));
                },
                |a, b| a.extend(b),
            )
        })
    };
    assert_eq!(run(1), run(4));
}
