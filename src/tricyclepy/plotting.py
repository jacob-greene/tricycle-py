"""Plot helpers, and the computations underneath them.

Port of ``R/plot_emb_circle_scale.R`` and ``R/plot_ccposition_den.R`` from R
``tricycle`` 1.12.0.

An angle wraps, so it must be drawn with a cyclic colormap whose two ends are
the same colour. :data:`CYCLIC_COLORS` is the palette tricycle uses.

``matplotlib`` is an optional dependency. The numeric parts of this module --
:func:`circular_density` above all -- work without it.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "CYCLIC_COLORS",
    "cyclic_colormap",
    "circular_density",
    "plot_embedding_circle_scale",
    "circle_scale_legend",
    "plot_ccposition_density",
]

#: The eight hues tricycle cycles through. The last is close to the first, so
#: the scale closes on itself.
CYCLIC_COLORS: Tuple[str, ...] = (
    "#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B",
    "#C4E416", "#4BBA0F", "#447D87", "#2C24E9",
)


def cyclic_colormap(name: str = "tricycle", n: int = 500):
    """Return the tricycle cyclic colormap as a matplotlib colormap."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name, list(CYCLIC_COLORS), N=n)


def circular_density(theta, bw: float = 30.0, n: int = 512,
                     at: Optional[Sequence[float]] = None
                     ) -> Tuple[np.ndarray, np.ndarray]:
    """Von Mises kernel density on the circle.

    This is ``circular::density.circular(circular(theta), bw = bw)`` with its
    default kernel and grid. ``bw`` is the von Mises concentration, not a
    standard deviation, so a larger ``bw`` gives a narrower kernel.

    Parameters
    ----------
    theta
        Angles in radians.
    bw
        Von Mises concentration. tricycle passes 30.
    n
        Grid size when ``at`` is not given.
    at
        Explicit evaluation grid in radians.

    Returns
    -------
    x : numpy.ndarray
        Evaluation grid in radians.
    y : numpy.ndarray
        Density at each grid point. Integrates to one over the circle.
    """
    # exp(bw*cos(.)) overflows for a large concentration, so factor out exp(bw)
    # against the exponentially scaled Bessel function: i0e(bw) == i0(bw)*exp(-bw).
    from scipy.special import i0e

    theta = np.asarray(theta, dtype=np.float64)
    if theta.size == 0:
        raise ValueError("theta is empty")
    grid = (np.linspace(0.0, 2.0 * np.pi, n) if at is None
            else np.asarray(at, dtype=np.float64))
    z = np.cos(grid[:, None] - theta[None, :])
    dens = np.exp(bw * (z - 1.0)).sum(axis=1)
    y = dens / (theta.size * 2.0 * np.pi * i0e(bw))
    return grid, y


def _as_arrays(adata, color_key, embedding_key):
    if hasattr(adata, "obsm"):
        emb = np.asarray(adata.obsm[embedding_key])[:, :2]
        colour = np.asarray(adata.obs[color_key], dtype=np.float64)
        return emb, colour
    raise TypeError("expected an AnnData with .obsm and .obs")


def plot_embedding_circle_scale(adata, *, embedding_key: str = "X_umap",
                                color_key: str = "tricyclePosition",
                                ax=None, point_size: float = 2.1,
                                alpha: float = 0.6, title: Optional[str] = None):
    """Scatter an embedding coloured by cell cycle position on a cyclic scale.

    Parameters
    ----------
    adata
        ``AnnData`` carrying the embedding in ``.obsm`` and the angle in
        ``.obs``.
    embedding_key
        Key in ``.obsm``. The first two columns are drawn.
    color_key
        Column of ``.obs`` holding the angle in radians.
    ax
        Axes to draw on. A new figure is made when this is ``None``.
    point_size, alpha, title
        Cosmetic.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    emb, colour = _as_arrays(adata, color_key, embedding_key)
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(emb[:, 0], emb[:, 1], c=colour, cmap=cyclic_colormap(),
               vmin=0.0, vmax=2.0 * np.pi, s=point_size, alpha=alpha,
               linewidths=0, rasterized=True)
    ax.set_xlabel(f"{embedding_key} 1")
    ax.set_ylabel(f"{embedding_key} 2")
    ax.set_title(title if title is not None else f"(n={emb.shape[0]})")
    return ax


def circle_scale_legend(ax=None, *, add_stage_label: bool = False,
                        g1_pos: float = 0.0, s_pos: float = 2.2,
                        g2m_pos: float = 3.9, n: int = 500):
    """Draw the cyclic colour legend as a ring.

    Parameters
    ----------
    ax
        Polar axes to draw on. A new polar figure is made when ``None``.
    add_stage_label
        Write the approximate stage names on the ring.
    g1_pos, s_pos, g2m_pos
        Angles in radians for those labels.
    n
        Number of wedges in the ring.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(2, 2), subplot_kw={"projection": "polar"})
    edges = np.linspace(0.0, 2.0 * np.pi, n + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    cmap = cyclic_colormap()
    ax.bar(centres, height=1.0, width=np.diff(edges), bottom=1.5,
           color=cmap(centres / (2.0 * np.pi)), linewidth=0)
    ax.set_yticks([])
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π"])
    ax.set_ylim(0, 3.2)
    if add_stage_label:
        for angle, label in ((g1_pos, "G1/G0"), (s_pos, "S"), (g2m_pos, "G2M")):
            ax.text(angle, 2.9, label, ha="center", va="center", fontsize=7)
    return ax


def plot_ccposition_density(theta, groups=None, *, bw: float = 30.0,
                            ax=None, palette: Optional[Sequence[str]] = None,
                            title: Optional[str] = None):
    """Draw circular densities of cell cycle position, overall and per group.

    Parameters
    ----------
    theta
        Angles in radians, within ``[0, 2*pi]``.
    groups
        One label per cell. A dashed black curve always shows all cells.
    bw
        Von Mises concentration passed to :func:`circular_density`.
    ax, palette, title
        Cosmetic.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    theta = np.asarray(theta, dtype=np.float64)
    if theta.min() < 0 or theta.max() > 2 * np.pi:
        raise ValueError("theta must be between 0 and 2*pi")
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))
    x, y = circular_density(theta, bw=bw)
    ax.plot(x, y, color="black", linestyle="--", linewidth=0.8, label="all")
    if groups is not None:
        groups = np.asarray(groups)
        levels = [lv for lv in dict.fromkeys(groups.tolist())]
        colours = palette if palette is not None else [
            f"C{i}" for i in range(len(levels))]
        for lv, colour in zip(levels, colours):
            sel = theta[groups == lv]
            if sel.size == 0:
                continue
            gx, gy = circular_density(sel, bw=bw)
            ax.plot(gx, gy, color=colour, linewidth=0.8, label=f"{lv} (n={sel.size})")
        ax.legend(fontsize=6, frameon=False)
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel("θ")
    ax.set_ylabel("Density")
    ax.set_title(title if title is not None else f"(n={theta.size})")
    return ax
