"""Local regression, checked against cases with a known closed-form answer."""

from __future__ import annotations

import numpy as np
import pytest

from tricyclepy import fit_periodic_loess, loess_fit

RNG = np.random.default_rng(11)


def test_degree_two_reproduces_a_quadratic_exactly():
    """A local quadratic fit of an exact quadratic returns it, at any span."""
    x = np.linspace(-3.0, 4.0, 200)
    y = 1.5 - 2.0 * x + 0.75 * x ** 2
    for span in (0.2, 0.5, 0.9):
        fit = loess_fit(x, y, span=span, degree=2)
        assert np.allclose(fit.fitted, y, atol=1e-10)


def test_degree_one_reproduces_a_line_exactly():
    x = np.linspace(0.0, 10.0, 150)
    y = -4.0 + 3.0 * x
    fit = loess_fit(x, y, span=0.3, degree=1)
    assert np.allclose(fit.fitted, y, atol=1e-10)


def test_degree_zero_is_the_tricube_weighted_window_mean():
    x = np.sort(RNG.uniform(0, 10, 120))
    y = RNG.normal(size=120)
    span = 0.4
    fit = loess_fit(x, y, span=span, degree=0)

    q = int(np.floor(120 * span + 1e-9))
    i = 60
    d = np.abs(x - x[i])
    order = np.argsort(d, kind="stable")[:q]
    h = d[order].max()
    u = d[order] / h
    w = np.where(u < 1, (1 - u ** 3) ** 3, 0.0)
    expected = float((w * y[order]).sum() / w.sum())
    assert fit.fitted[i] == pytest.approx(expected, rel=1e-12)


def test_fit_does_not_depend_on_input_order():
    x = RNG.uniform(0, 5, 90)
    y = np.sin(x) + RNG.normal(0, 0.1, 90)
    a = loess_fit(x, y, span=0.4).fitted
    shuffle = RNG.permutation(90)
    b = loess_fit(x[shuffle], y[shuffle], span=0.4).fitted
    assert np.allclose(a[shuffle], b, atol=1e-12)


def test_predict_at_the_data_points_reproduces_the_fitted_values():
    x = np.sort(RNG.uniform(0, 8, 130))
    y = np.cos(x) + RNG.normal(0, 0.2, 130)
    fit = loess_fit(x, y, span=0.35)
    assert np.allclose(fit.predict(x), fit.fitted, atol=1e-12)


def test_span_above_one_widens_the_bandwidth():
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    narrow = loess_fit(x, y, span=1.0, degree=1).fitted
    wide = loess_fit(x, y, span=4.0, degree=1).fitted
    # A wider bandwidth flattens the fit, so it departs further from the data.
    assert np.abs(wide - y).max() > np.abs(narrow - y).max()


def test_constant_response_is_returned_unchanged():
    x = np.sort(RNG.uniform(0, 3, 50))
    y = np.full(50, 7.0)
    assert np.allclose(loess_fit(x, y, span=0.3).fitted, 7.0, atol=1e-12)


def test_degree_is_validated():
    with pytest.raises(ValueError):
        loess_fit(np.arange(10.0), np.arange(10.0), degree=3)


# --- the periodic wrapper --------------------------------------------------


def test_periodic_prediction_curve_closes_on_itself():
    theta = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    result = fit_periodic_loess(theta, np.sin(theta), span=0.3)
    # 0 and 2*pi are the same angle, so the curve must meet itself there.
    assert result.pred_y[0] == pytest.approx(result.pred_y[-1], abs=1e-9)


def test_periodic_fit_is_invariant_to_rotating_the_predictor():
    """Rotating every angle by the same amount must rotate the fit with it."""
    theta = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    y = np.cos(theta) + 0.3 * np.sin(3 * theta)
    base = fit_periodic_loess(theta, y, span=0.3).fitted

    shift = 90                                  # a quarter turn, on the grid
    rotated = np.mod(theta + shift * (2 * np.pi / 360), 2 * np.pi)
    moved = fit_periodic_loess(rotated, y, span=0.3).fitted
    assert np.allclose(base, moved, atol=1e-10)


def test_periodic_fit_beats_a_non_periodic_one_at_the_seam():
    """The seam is where a non-periodic fit fails, so compare it there only.

    Both fits get the same effective bandwidth: the periodic fit triples the
    data, so its span has to be divided by three to match.
    """
    theta = np.linspace(0, 2 * np.pi, 300, endpoint=False)
    y = np.cos(theta)
    periodic = fit_periodic_loess(theta, y, span=0.1).fitted
    plain = loess_fit(theta, y, span=0.3, degree=2).fitted
    seam = np.r_[0:5, 295:300]
    assert np.abs(periodic[seam] - y[seam]).max() < \
        np.abs(plain[seam] - y[seam]).max()


def test_periodic_fit_tracks_a_smooth_cycle_at_a_narrow_span():
    theta = np.linspace(0, 2 * np.pi, 600, endpoint=False)
    y = np.sin(theta)
    result = fit_periodic_loess(theta, y, span=0.05)
    assert np.abs(result.fitted - y).max() < 5e-3


def test_a_wide_span_shrinks_the_fit_towards_the_mean():
    """span applies to the tripled data, so 0.3 is nearly a whole period."""
    theta = np.linspace(0, 2 * np.pi, 600, endpoint=False)
    y = np.sin(theta)
    narrow = fit_periodic_loess(theta, y, span=0.05).fitted
    wide = fit_periodic_loess(theta, y, span=0.3).fitted
    assert np.abs(wide).max() < np.abs(narrow).max()
    assert np.abs(wide).max() < np.abs(y).max()


def test_rsquared_is_one_minus_rss_over_the_unrepeated_tss():
    theta = np.linspace(0, 2 * np.pi, 250, endpoint=False)
    y = np.sin(theta) + RNG.normal(0, 0.2, 250)
    result = fit_periodic_loess(theta, y, span=0.3)
    tss = float(((y - y.mean()) ** 2).sum())
    rss = float((result.residual ** 2).sum())
    assert result.rsquared == pytest.approx(1.0 - rss / tss, rel=1e-12)
    assert np.allclose(result.residual, y - result.fitted)


def test_prediction_grid_spans_the_full_circle():
    theta = RNG.uniform(0, 2 * np.pi, 120)
    result = fit_periodic_loess(theta, RNG.normal(size=120), length_out=37)
    assert result.pred_x.shape == (37,)
    assert result.pred_x[0] == 0.0
    assert result.pred_x[-1] == pytest.approx(2 * np.pi)


def test_theta_outside_the_circle_is_rejected():
    with pytest.raises(ValueError, match="between 0 and 2"):
        fit_periodic_loess(np.array([-0.1, 1.0]), np.array([0.0, 1.0]))


def test_mismatched_lengths_are_rejected():
    with pytest.raises(ValueError, match="same length"):
        fit_periodic_loess(np.array([0.1, 1.0]), np.array([0.0]))
