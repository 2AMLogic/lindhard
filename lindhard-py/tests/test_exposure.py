"""Exposure dose and iso-dose contours."""

import lindhard as lh
import numpy as np
import pytest

PSF = dict(total_ev=1.0e4, alpha_m=4e-8, beta_m=1.2e-7, eta=0.7)


def test_isolated_shot_is_the_psf():
    x = np.array([0.0, 3e-8, 1e-7])
    d = lh.exposure_dose(np.array([[0.0, 0.0, 1.0]]), x, np.zeros(3), **PSF)
    a, b, e = PSF["alpha_m"], PSF["beta_m"], PSF["eta"]
    f = (np.exp(-(x / a) ** 2) / a**2 + e * np.exp(-(x / b) ** 2) / b**2) / (
        np.pi * (1 + e)
    )
    assert np.allclose(d, 1.0e4 * f, rtol=1e-12)


def test_grid_threads_and_contour():
    shots = "x_m,y_m,weight,shape,size_m\n0,0,1,disc,2e-8\n"
    args = (shots, -2e-7, -2e-7, 1e-8, 1e-8, 41, 41)
    g1 = lh.exposure_dose_grid(*args, **PSF, threads=1)
    g4 = lh.exposure_dose_grid(*args, **PSF, threads=4)
    assert g1.shape == (41, 41) and np.array_equal(g1, g4)
    c = lh.iso_contours(g1, -2e-7, -2e-7, 1e-8, 1e-8, 0.5 * g1.max())
    assert len(c) == 1 and c[0]["closed"] and c[0]["area_m2"] > 0


def test_area_dose_and_errors():
    a = lh.area_dose(**PSF, weight_per_m2=2.0)
    assert a["total_ev_per_m2"] == pytest.approx(2.0e4)
    assert a["forward_ev_per_m2"] == pytest.approx(2.0e4 / 1.7)
    with pytest.raises(lh.InputError):
        lh.area_dose(**PSF, gamma_m=1e-6)
    with pytest.raises(lh.InputError):
        lh.exposure_dose("x_m,y_m\n1,2\n", np.zeros(1), np.zeros(1), **PSF)


def _electron_run():
    from pathlib import Path

    case = Path(__file__).resolve().parents[2] / "examples" / "electron" / "e_10keV_si.toml"
    text = case.read_text().replace("points_per_decade = 8.0", "points_per_decade = 3.0")
    return lh.ElectronRun.from_toml(text, base_dir=case.parent)


def test_exposure_direct_shapes_threads_and_errors():
    run = _electron_run()
    shots = np.array([[0.0, 0.0, 1.0], [5e-8, 0.0, 2.0]])
    args = (shots, run, -2e-7, -2e-7, 2e-8, 2e-8, 20, 20)
    kw = dict(
        **PSF, slab_hi_m=5e-8, histories_per_shot=80, replicas=2, blur_m=1e-8
    )
    r1 = lh.exposure_direct(*args, **kw, threads=1)
    r4 = lh.exposure_direct(*args, **kw, threads=4)
    assert set(r1) == {
        "dose_ev_per_m2",
        "std_err_ev_per_m2",
        "convolution_ev_per_m2",
        "deviation_sigmas",
        "flagged",
    }
    for k, v in r1.items():
        assert v.shape == (20, 20)
        assert np.array_equal(v, r4[k], equal_nan=True)
    assert r1["flagged"].dtype == np.bool_
    with pytest.raises(lh.InputError):
        lh.exposure_direct(*args, **{**kw, "replicas": 1})
