"""Diagnostic for UMI datasets whose cells are mostly not cycling.

Port of ``R/diagnose_totalumi.R`` from R ``tricycle`` 1.12.0.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .loess import PeriodicLoessResult, fit_periodic_loess

__all__ = ["diagnose_total_umi", "TotalUmiDiagnosis"]


@dataclass(frozen=True)
class TotalUmiDiagnosis:
    """Result of :func:`diagnose_total_umi`."""

    fit: PeriodicLoessResult
    peak: float
    valley: float
    difference: float
    passed: bool
    fraction_between: float


def diagnose_total_umi(theta, total_umis, span: float = 0.3,
                       length_out: int = 200) -> TotalUmiDiagnosis:
    """Check whether the estimated cell cycle position is trustworthy.

    A cycling population has its highest total UMI count near the end of S and
    roughly half of that near M. The function fits ``log2(total_umis + 1)``
    against the position and compares the fitted peak on ``(0.75*pi, 1.25*pi)``
    with the fitted valley on ``(1.35*pi, 1.85*pi)``. A gap below 0.4 means the
    embedding centre is probably shifted and the positions should not be used.

    A second warning fires when more than half the cells sit between ``0.5*pi``
    and ``1.8*pi``, which points the same way.

    Parameters
    ----------
    theta
        Cell cycle position in radians, within ``[0, 2*pi]``.
    total_umis
        Total UMI count per cell, not log transformed.
    span, length_out
        Passed to :func:`~tricyclepy.fit_periodic_loess`.

    Returns
    -------
    TotalUmiDiagnosis
        ``passed`` is ``False`` when the peak-to-valley gap is below 0.4.
    """
    theta = np.asarray(theta, dtype=np.float64)
    total_umis = np.asarray(total_umis, dtype=np.float64)
    if theta.shape != total_umis.shape:
        raise ValueError("theta and total_umis must have the same length")
    if theta.min() < 0 or theta.max() > 2 * np.pi:
        raise ValueError("theta must be between 0 and 2*pi")

    fit = fit_periodic_loess(theta, np.log2(total_umis + 1.0), span=span,
                             length_out=length_out)
    pi = np.pi
    peak = float(fit.pred_y[(fit.pred_x > 0.75 * pi) & (fit.pred_x < 1.25 * pi)].max())
    valley = float(fit.pred_y[(fit.pred_x > 1.35 * pi) & (fit.pred_x < 1.85 * pi)].min())
    difference = peak - valley
    passed = difference >= 0.4
    if not passed:
        warnings.warn(
            "The difference of predicted log2(total UMIs) between the end of S "
            f"phase and M phase is {difference:.2f}, below 0.4. The cell cycle "
            "positions may be wrong because the embedding centre is shifted. "
            "This happens when most cells in the dataset are not cycling.",
            RuntimeWarning, stacklevel=2,
        )

    fraction_between = float(np.mean((theta > 0.5 * pi) & (theta < 1.8 * pi)))
    if fraction_between > 0.5:
        warnings.warn(
            "More than 50 percent of cells sit between 0.5*pi and 1.8*pi. The "
            "cell cycle positions may be wrong because the embedding centre is "
            "shifted.", RuntimeWarning, stacklevel=2,
        )
    return TotalUmiDiagnosis(fit=fit, peak=peak, valley=valley,
                             difference=difference, passed=passed,
                             fraction_between=fraction_between)
