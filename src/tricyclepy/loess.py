"""Local regression, and the periodic loess fit tricycle builds on it.

Port of ``R/fit_periodic_loess.R`` from R ``tricycle`` 1.12.0.

What is reproduced, and what is not
-----------------------------------
R's ``stats::loess`` has two surfaces. ``surface="direct"`` evaluates the local
regression at every point. ``surface="interpolate"``, the default, evaluates it
at the vertices of a k-d tree and interpolates between them, which is faster and
approximate.

:func:`loess_fit` here is the direct surface. It agrees with
``loess(..., control = loess.control(surface = "direct"))`` to floating-point
noise. Against R's default it agrees only to the k-d tree interpolation error.
Both gaps are measured in ``equivalence/`` and reported; neither is assumed.

Only the pieces tricycle uses are implemented: one predictor, Gaussian family,
no robustness iterations, no surface interpolation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

__all__ = ["loess_fit", "LoessFit", "fit_periodic_loess", "PeriodicLoessResult"]


def _neighbourhood_size(n: int, span: float) -> int:
    """Number of points in a local neighbourhood, as R's C code computes it."""
    if span <= 0:
        raise ValueError("span must be positive")
    q = int(np.floor(n * span + 1e-9))
    return int(min(max(q, 2), n))


def _window_starts(xs: np.ndarray, at: np.ndarray, q: int) -> np.ndarray:
    """Start index of the window of the ``q`` nearest neighbours of each point.

    In one dimension the ``q`` nearest neighbours of a point are contiguous in
    sorted order, so a window is fully described by its start. The start is
    found by bisection on the "stop moving right" test, which is monotone in
    the start index.
    """
    n = xs.shape[0]
    if q >= n:
        return np.zeros(at.shape[0], dtype=np.intp)
    last = n - q
    lo = np.zeros(at.shape[0], dtype=np.intp)
    hi = np.full(at.shape[0], last, dtype=np.intp)
    # The leftmost distance falls and the rightmost distance rises as the start
    # moves right, so "the left edge is farther than the right edge" is true for
    # small starts and false for large ones. Bisect for the first false.
    for _ in range(int(np.ceil(np.log2(last + 1))) + 1):
        active = lo < hi
        mid = (lo + hi) // 2
        move_right = (at - xs[mid]) > (xs[mid + q - 1] - at)
        lo = np.where(active & move_right, mid + 1, lo)
        hi = np.where(active & ~move_right, mid, hi)

    # The crossing sits between lo-1 and lo. Take whichever window is tighter.
    below = np.maximum(lo - 1, 0)
    width = lambda start: np.maximum(at - xs[start], xs[start + q - 1] - at)  # noqa: E731
    return np.where(width(below) <= width(lo), below, lo)


def _local_fit_block(xs: np.ndarray, ys: np.ndarray, at: np.ndarray, q: int,
                     span: float, degree: int) -> np.ndarray:
    """Fit a local polynomial at every point of ``at``. ``xs`` is sorted."""
    n = xs.shape[0]
    q = int(min(q, n))
    starts = _window_starts(xs, at, q)
    idx = starts[:, None] + np.arange(q, dtype=np.intp)[None, :]
    dx = xs[idx] - at[:, None]
    d = np.abs(dx)
    h = d.max(axis=1)
    if span > 1.0:
        h = h * span            # one predictor, so span**(1/p) is span
    yw = ys[idx]

    safe_h = np.where(h > 0.0, h, 1.0)
    u = d / safe_h[:, None]
    w = np.where(u < 1.0, (1.0 - u ** 3) ** 3, 0.0)

    # Scaling the local basis by the bandwidth keeps the normal equations
    # well conditioned, and leaves the value at the centre unchanged.
    t = dx / safe_h[:, None]
    p = degree + 1
    moments = np.empty((at.shape[0], 2 * degree + 1), dtype=np.float64)
    rhs = np.empty((at.shape[0], p), dtype=np.float64)
    term = w.copy()                 # w * t**k, built up one power at a time
    wy = w * yw
    term_y = wy.copy()
    for k in range(2 * degree + 1):
        if k:
            term *= t
        moments[:, k] = term.sum(axis=1)
    for k in range(p):
        if k:
            term_y *= t
        rhs[:, k] = term_y.sum(axis=1)

    gram = np.empty((at.shape[0], p, p), dtype=np.float64)
    for a in range(p):
        for b in range(p):
            gram[:, a, b] = moments[:, a + b]

    out = np.empty(at.shape[0], dtype=np.float64)
    degenerate = (h <= 0.0) | (moments[:, 0] <= 0.0)
    good = ~degenerate
    if good.any():
        out[good] = np.linalg.solve(gram[good], rhs[good][..., None])[:, 0, 0]
    if degenerate.any():
        # No spread, or every weight zero: fall back to the window mean, which
        # is what a local fit degenerates to.
        out[degenerate] = yw[degenerate].mean(axis=1)
    return out


@dataclass(frozen=True)
class LoessFit:
    """A fitted direct-surface loess, with the data needed to predict again."""

    x: np.ndarray
    y: np.ndarray
    span: float
    degree: int
    #: Fitted values in the order the data were given. Empty when the caller
    #: only needed the fit at other points, as :func:`fit_periodic_loess` does.
    fitted: np.ndarray

    def predict(self, newx: np.ndarray) -> np.ndarray:
        """Evaluate the fit at new predictor values."""
        return _evaluate(self.x, self.y, np.asarray(newx, dtype=np.float64),
                         self.span, self.degree)


def _evaluate(xs: np.ndarray, ys: np.ndarray, at: np.ndarray, span: float,
              degree: int, block: int = 2048) -> np.ndarray:
    """Evaluate the local fit at every point of ``at``, in memory-bounded blocks.

    Each block holds a dense ``len(block) by q`` neighbourhood array, so ``block``
    trades peak memory against the number of passes. Results do not depend on it.
    """
    q = _neighbourhood_size(xs.shape[0], span)
    out = np.empty(at.shape[0], dtype=np.float64)
    for start in range(0, at.shape[0], block):
        stop = min(start + block, at.shape[0])
        out[start:stop] = _local_fit_block(xs, ys, at[start:stop], q, span, degree)
    return out


def loess_fit(x, y, span: float = 0.3, degree: int = 2) -> LoessFit:
    """Fit a direct-surface loess of ``y`` on a single predictor ``x``.

    Parameters
    ----------
    x, y
        Predictor and response, same length.
    span
        Fraction of points in each local neighbourhood. R's ``span``.
    degree
        Degree of the local polynomial, 0, 1 or 2. R's default is 2.

    Returns
    -------
    LoessFit
        Carries the fitted values and can predict at new points.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same length")
    if degree not in (0, 1, 2):
        raise ValueError("degree must be 0, 1 or 2")
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ys = y[order]
    fitted_sorted = _evaluate(xs, ys, xs, span, degree)
    fitted = np.empty_like(fitted_sorted)
    fitted[order] = fitted_sorted
    return LoessFit(x=xs, y=ys, span=float(span), degree=int(degree), fitted=fitted)


@dataclass(frozen=True)
class PeriodicLoessResult:
    """Output of :func:`fit_periodic_loess`, mirroring the R list."""

    fitted: np.ndarray
    residual: np.ndarray
    pred_x: np.ndarray
    pred_y: np.ndarray
    rsquared: float
    loess: LoessFit


def fit_periodic_loess(theta, y, span: float = 0.3, length_out: int = 200,
                       degree: int = 2) -> PeriodicLoessResult:
    """Fit a loess line of ``y`` against the circular predictor ``theta``.

    Circularity is handled the way R does it: the predictor is repeated at
    ``theta - 2*pi``, ``theta`` and ``theta + 2*pi``, the response is repeated
    three times, one ordinary loess is fitted, and only the middle third is
    returned.

    Parameters
    ----------
    theta
        Cell cycle position in radians, within ``[0, 2*pi]``.
    y
        Response, same length as ``theta``.
    span
        Smoothing parameter.
    length_out
        Number of points on the returned prediction curve.
    degree
        Degree of the local polynomial.

    Returns
    -------
    PeriodicLoessResult
        ``rsquared`` is ``1 - RSS / TSS`` where ``TSS`` is the total sum of
        squares of the original, unrepeated ``y``.
    """
    theta = np.asarray(theta, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if theta.shape != y.shape:
        raise ValueError("theta and y must have the same length")
    if theta.size == 0:
        raise ValueError("theta is empty")
    if theta.min() < 0 or theta.max() > 2 * np.pi:
        raise ValueError("theta must be between 0 and 2*pi")

    n = theta.shape[0]
    two_pi = 2.0 * np.pi
    x3 = np.concatenate([theta - two_pi, theta, theta + two_pi])
    y3 = np.concatenate([y, y, y])
    ss_total = float(np.sum((y - y.mean()) ** 2))

    # R fits on all 3n points and then discards two thirds of the fitted values.
    # Only the middle third and the prediction grid are ever used, so evaluate
    # the fit there and nowhere else. The fit itself is unchanged.
    order = np.argsort(x3, kind="mergesort")
    fit = LoessFit(x=x3[order], y=y3[order], span=float(span), degree=int(degree),
                   fitted=np.empty(0))
    fitted = _evaluate(fit.x, fit.y, theta, span, degree)
    residual = y - fitted
    # A constant response has no total sum of squares. R returns NaN there
    # rather than raising, and so does this.
    rsquared = (np.nan if ss_total == 0.0
                else 1.0 - float(np.sum(residual ** 2)) / ss_total)

    pred_x = np.linspace(0.0, two_pi, length_out)
    pred_y = fit.predict(pred_x)
    return PeriodicLoessResult(fitted=fitted, residual=residual, pred_x=pred_x,
                               pred_y=pred_y, rsquared=rsquared, loess=fit)
