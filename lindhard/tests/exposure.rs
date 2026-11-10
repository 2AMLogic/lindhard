//! Multi-shot exposure from a Gaussian-sum PSF (`lindhard::exposure`).

use std::f64::consts::PI;

use lindhard::exposure::{parse_shots_csv, uniform_area_dose, Exposure, Grid, Shot, SpotShape};
use lindhard::tally::GaussianPsf;

const NM: f64 = 1e-9;

fn double() -> GaussianPsf {
    GaussianPsf::double(1.0e4, 40.0 * NM, 120.0 * NM, 0.7).unwrap()
}

fn triple() -> GaussianPsf {
    GaussianPsf::triple(1.0e4, 40.0 * NM, 120.0 * NM, 0.7, 400.0 * NM, 0.2).unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

#[test]
fn isolated_shot_is_the_psf() {
    for psf in [double(), triple()] {
        let e = Exposure::new(psf, vec![Shot::point(5.0 * NM, -3.0 * NM, 1.0)], 0.0).unwrap();
        for k in 0..20 {
            let (x, y) = (k as f64 * 30.0 * NM, k as f64 * 7.0 * NM);
            let r = ((x - 5.0 * NM).powi(2) + (y + 3.0 * NM).powi(2)).sqrt();
            let want = psf.density_ev_per_m2(r);
            assert!(rel(e.dose_at(x, y), want) < 1e-12, "k={k}");
        }
    }
}

#[test]
fn weight_scales_and_overlapping_shots_add() {
    let psf = double();
    let one = Exposure::new(psf, vec![Shot::point(0.0, 0.0, 1.0)], 0.0).unwrap();
    let two = Exposure::new(
        psf,
        vec![Shot::point(0.0, 0.0, 1.0), Shot::point(0.0, 0.0, 2.0)],
        0.0,
    )
    .unwrap();
    let (a, b) = (one.dose_at(20.0 * NM, 0.0), two.dose_at(20.0 * NM, 0.0));
    assert!(rel(b, 3.0 * a) < 1e-12);
}

#[test]
fn large_uniform_area_gives_area_dose_with_backscatter() {
    for psf in [double(), triple()] {
        let pitch = 20.0 * NM;
        let n = 100; // 2 µm square, more than 4 of the longest width each side
        let mut shots = Vec::new();
        for j in 0..n {
            for i in 0..n {
                shots.push(Shot::point(i as f64 * pitch, j as f64 * pitch, 1.0));
            }
        }
        let e = Exposure::new(psf, shots, 0.0).unwrap();
        let c = (n as f64 - 1.0) / 2.0 * pitch;
        let got = e.dose_at(c, c);
        let area = uniform_area_dose(&psf, 1.0 / (pitch * pitch));
        assert!(rel(got, area.total_ev_per_m2) < 1e-4, "{got} vs {area:?}");
        // The forward-only part is a strict fraction of the total.
        let share = 1.0 / (1.0 + 0.7 + if psf.gamma_m.is_some() { 0.2 } else { 0.0 });
        assert!(rel(area.forward_ev_per_m2, area.total_ev_per_m2 * share) < 1e-12);
        assert!(
            rel(
                area.forward_ev_per_m2 + area.long_range_ev_per_m2,
                area.total_ev_per_m2
            ) < 1e-12
        );
        assert!(got > 1.2 * area.forward_ev_per_m2);
    }
}

#[test]
fn zero_blur_is_the_unblurred_result_and_quadrature_matches_the_integral() {
    let psf = double();
    let (x, y) = (70.0 * NM, 10.0 * NM);
    let plain = Exposure::new(psf, vec![Shot::point(0.0, 0.0, 1.0)], 0.0).unwrap();
    let zero = Exposure::new(psf, vec![Shot::point(0.0, 0.0, 1.0)], 0.0).unwrap();
    assert_eq!(plain.dose_at(x, y), zero.dose_at(x, y));

    // Numerical ∫ f(|r - s|) exp(-s²/b²)/(π b²) d²s for b = 30 nm.
    let b = 30.0 * NM;
    let n = 400;
    let half = 5.0 * b;
    let h = 2.0 * half / n as f64;
    let mut sum = 0.0;
    for j in 0..n {
        for i in 0..n {
            let sx = -half + (i as f64 + 0.5) * h;
            let sy = -half + (j as f64 + 0.5) * h;
            let g = (-(sx * sx + sy * sy) / (b * b)).exp() / (PI * b * b);
            let r = ((x - sx).powi(2) + (y - sy).powi(2)).sqrt();
            sum += psf.density_ev_per_m2(r) * g * h * h;
        }
    }
    for blurred in [
        Exposure::new(psf, vec![Shot::point(0.0, 0.0, 1.0)], b).unwrap(),
        Exposure::new(
            psf,
            vec![Shot {
                blur_m: b,
                ..Shot::point(0.0, 0.0, 1.0)
            }],
            0.0,
        )
        .unwrap(),
    ] {
        assert!(rel(blurred.dose_at(x, y), sum) < 1e-6);
    }
    // Global and per-shot blur add in quadrature.
    let both = Exposure::new(
        psf,
        vec![Shot {
            blur_m: 24.0 * NM,
            ..Shot::point(0.0, 0.0, 1.0)
        }],
        18.0 * NM,
    )
    .unwrap();
    let direct = Exposure::new(psf, vec![Shot::point(0.0, 0.0, 1.0)], b).unwrap();
    assert!(rel(both.dose_at(x, y), direct.dose_at(x, y)) < 1e-12);
}

#[test]
fn shaped_spots_conserve_weight_in_the_far_field() {
    let psf = double();
    let (x, y) = (50.0 * NM, 0.0);
    // A 100 nm beam blur keeps the dose smooth over the spot (the tail of a
    // narrow Gaussian is too steep for a point-source comparison).
    let blur = 100.0 * NM;
    let point = Exposure::new(psf, vec![Shot::point(0.0, 0.0, 2.0)], blur)
        .unwrap()
        .dose_at(x, y);
    for shape in [
        SpotShape::Disc {
            radius_m: 15.0 * NM,
        },
        SpotShape::Rect {
            width_m: 20.0 * NM,
            height_m: 30.0 * NM,
        },
    ] {
        let e = Exposure::new(
            psf,
            vec![Shot {
                shape,
                ..Shot::point(0.0, 0.0, 2.0)
            }],
            blur,
        )
        .unwrap();
        assert!(rel(e.dose_at(x, y), point) < 1e-2, "{shape:?}");
    }
}

#[test]
fn grid_is_bit_identical_at_1_2_4_threads() {
    let psf = triple();
    let shots: Vec<Shot> = (0..50)
        .map(|k| {
            let t = k as f64;
            Shot {
                blur_m: (k % 3) as f64 * 5.0 * NM,
                shape: if k % 7 == 0 {
                    SpotShape::Disc {
                        radius_m: 10.0 * NM,
                    }
                } else {
                    SpotShape::Point
                },
                ..Shot::point(
                    (t * 37.0 % 400.0) * NM,
                    (t * 91.0 % 400.0) * NM,
                    1.0 + t / 10.0,
                )
            }
        })
        .collect();
    let e = Exposure::new(psf, shots, 8.0 * NM).unwrap();
    let grid = Grid {
        x0_m: -50.0 * NM,
        y0_m: -50.0 * NM,
        dx_m: 13.0 * NM,
        dy_m: 11.0 * NM,
        nx: 40,
        ny: 45,
    };
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| e.dose_on_grid(&grid).unwrap())
    };
    let base = run(1);
    for t in [2, 4] {
        let other = run(t);
        let same = base
            .dose_ev_per_m2
            .iter()
            .zip(&other.dose_ev_per_m2)
            .all(|(a, b)| a.to_bits() == b.to_bits());
        assert!(same, "threads={t}");
    }
}

#[test]
fn shot_outside_grid_still_contributes_and_bad_input_is_rejected() {
    let psf = double();
    let e = Exposure::new(psf, vec![Shot::point(1e-6, 0.0, 1.0)], 0.0).unwrap();
    assert!(e.dose_at(0.0, 0.0) > 0.0);
    assert!(Exposure::new(psf, vec![Shot::point(f64::NAN, 0.0, 1.0)], 0.0).is_err());
    assert!(Exposure::new(psf, vec![], -1.0).is_err());
    let g = Grid {
        x0_m: 0.0,
        y0_m: 0.0,
        dx_m: 0.0,
        dy_m: 1e-9,
        nx: 2,
        ny: 2,
    };
    assert!(e.dose_on_grid(&g).is_err());
}

#[test]
fn shot_csv_round_trip_and_errors() {
    let shots = parse_shots_csv(
        "# comment\nx_m,y_m,weight,blur_m,shape,size_m,size2_m\n\
         0,0,1,,point,,\n1e-8,0,2.5,5e-9,disc,2e-8,\n0,1e-8,1,,rect,3e-8,1e-8\n",
    )
    .unwrap();
    assert_eq!(shots.len(), 3);
    assert_eq!(shots[1].shape, SpotShape::Disc { radius_m: 2e-8 });
    assert_eq!(shots[1].blur_m, 5e-9);
    assert_eq!(
        shots[2].shape,
        SpotShape::Rect {
            width_m: 3e-8,
            height_m: 1e-8
        }
    );
    assert_eq!(
        parse_shots_csv("x_m,y_m,weight\n1,2,3\n").unwrap()[0].weight,
        3.0
    );
    assert!(parse_shots_csv("x_m,y_m\n1,2\n").is_err());
    assert!(parse_shots_csv("x_m,y_m,weight\n1,2\n").is_err());
    assert!(parse_shots_csv("x_m,y_m,weight,shape\n1,2,3,disc\n").is_err());
    assert!(parse_shots_csv("x_m,y_m,weight\n1,x,3\n").is_err());
}

#[test]
fn half_dose_contour_of_a_square_pad_has_the_pad_area() {
    use lindhard::contour::contours;
    // A straight edge of a uniform exposure has exactly half the interior
    // dose under every Gaussian term, so the half-dose contour is the pad
    // edge up to corner rounding.
    let psf = double();
    let (pitch, n) = (20.0 * NM, 60);
    let mut shots = Vec::new();
    for j in 0..n {
        for i in 0..n {
            shots.push(Shot::point(i as f64 * pitch, j as f64 * pitch, 1.0));
        }
    }
    let e = Exposure::new(psf, shots, 0.0).unwrap();
    let grid = Grid {
        x0_m: -300.0 * NM,
        y0_m: -300.0 * NM,
        dx_m: 16.0 * NM,
        dy_m: 16.0 * NM,
        nx: 110,
        ny: 110,
    };
    let dose = e.dose_on_grid(&grid).unwrap();
    let half = 0.5 * uniform_area_dose(&psf, 1.0 / (pitch * pitch)).total_ev_per_m2;
    let c = contours(&dose, half);
    assert_eq!(c.len(), 1);
    // The shots sit at cell centres of an n·pitch square.
    let side = n as f64 * pitch;
    let area = c[0].signed_area_m2.unwrap();
    assert!(rel(area, side * side) < 0.03, "{area} vs {}", side * side);
    let m = c[0].centroid_m.unwrap();
    let mid = (n as f64 - 1.0) / 2.0 * pitch;
    assert!((m[0] - mid).abs() < 2.0 * NM && (m[1] - mid).abs() < 2.0 * NM);
}
