"""Build the README figure: one row per dataset, three panels per row.

Panels, left to right:

1. UMAP coloured by cell type.
2. UMAP coloured by the cell cycle position this package computes.
3. Python position against R ``tricycle`` position for the same cells, with
   the circular correlation coefficient.

Each input is an ``.h5ad`` that holds:

* ``X``: log-normalised expression, cells by genes, gene symbols as
  ``var_names``;
* ``obsm["X_umap"]``: a two-dimensional embedding;
* ``obs[<celltype key>]``: a categorical cell type label;
* ``obs[<R key>]``: ``tricyclePosition`` from R ``tricycle`` on the same
  matrix (see ``equivalence/run_r.R``).

The Python position is computed here, from ``X``, by this package. Nothing is
read from a previous Python run.

Usage::

    python docs/readme_figure.py \\
        --row "Human bone marrow" bmmc.h5ad \\
        --row "Human CD34+ HSPC" cd34.h5ad \\
        --out docs/readme_figure.png

The circular correlation is the Jammalamadaka-SenGupta coefficient
(Jammalamadaka and SenGupta 2001, *Topics in Circular Statistics*, eq. 8.2.2):

    r = sum(sin(a - mean_a) * sin(b - mean_b))
        / sqrt(sum(sin(a - mean_a)**2) * sum(sin(b - mean_b)**2))

where ``mean_a`` and ``mean_b`` are circular means. It is the statistic
``astropy.stats.circcorrcoef`` computes.
"""

from __future__ import annotations

import argparse

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

import tricyclepy as tp  # noqa: E402

TWO_PI = 2.0 * np.pi


def circular_mean(a: np.ndarray) -> float:
    return float(np.arctan2(np.sin(a).sum(), np.cos(a).sum()))


def circular_correlation(a, b) -> float:
    """Jammalamadaka-SenGupta circular correlation of two angle vectors."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    sa = np.sin(a - circular_mean(a))
    sb = np.sin(b - circular_mean(b))
    return float((sa * sb).sum() / np.sqrt((sa ** 2).sum() * (sb ** 2).sum()))


def circular_difference(a, b) -> np.ndarray:
    """Absolute angular difference, wrapped onto [0, pi]."""
    d = np.mod(np.asarray(a) - np.asarray(b), TWO_PI)
    return np.minimum(d, TWO_PI - d)


def style() -> None:
    rc = plt.rcParams
    rc["figure.dpi"] = 150
    rc["savefig.dpi"] = 150
    for side in ("top", "right", "bottom", "left"):
        rc[f"axes.spines.{side}"] = False
    rc["axes.grid"] = False
    rc["font.family"] = "sans-serif"
    rc["font.sans-serif"] = ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
    rc["font.size"] = 10
    rc["pdf.fonttype"] = 42


def celltype_palette(categories) -> dict:
    """One fixed colour per cell type, shared by every row."""
    base = [matplotlib.colors.rgb2hex(c) for c in matplotlib.colormaps["Paired"].colors]
    extra = [matplotlib.colors.rgb2hex(c) for c in matplotlib.colormaps["Dark2"].colors]
    colours = base + extra
    cats = sorted(categories)
    if len(cats) > len(colours):
        raise ValueError(f"{len(cats)} cell types; at most {len(colours)} supported")
    return dict(zip(cats, colours))


def umap_axes(ax) -> None:
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--row", nargs=2, action="append", required=True,
                   metavar=("LABEL", "H5AD"), help="row label and input file")
    p.add_argument("--celltype-key", default="celltype")
    p.add_argument("--r-key", default="tricyclePosition_R")
    p.add_argument("--species", default="human")
    p.add_argument("--gname-type", default="SYMBOL")
    p.add_argument("--point-size", type=float, default=2.0)
    p.add_argument("--seed", type=int, default=0, help="draw-order shuffle")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    style()
    rows = []
    for label, path in args.row:
        a = ad.read_h5ad(path)
        tp.estimate_cycle_position(a, species=args.species, gname_type=args.gname_type)
        rows.append((label, a))

    palette = celltype_palette(
        set().union(*(set(a.obs[args.celltype_key].astype(str)) for _, a in rows)))
    cmap = tp.cyclic_colormap()
    norm = Normalize(0.0, TWO_PI)
    rng = np.random.default_rng(args.seed)

    n = len(rows)
    fig = plt.figure(figsize=(15.0, 4.2 * n))
    gs = fig.add_gridspec(n, 6, width_ratios=[1.0, 0.62, 1.0, 0.04, 0.22, 1.0],
                          wspace=0.08, hspace=0.25)

    for i, (label, a) in enumerate(rows):
        xy = np.asarray(a.obsm["X_umap"])
        ct = a.obs[args.celltype_key].astype(str).to_numpy()
        py = a.obs["tricyclePosition"].to_numpy(np.float64)
        rr = a.obs[args.r_key].to_numpy(np.float64)
        order = rng.permutation(a.n_obs)
        r = circular_correlation(py, rr)
        dmax = float(circular_difference(py, rr).max())
        n_genes = a.uns["tricycleEmbedding"]["n_genes"]
        print(f"{label}: cells={a.n_obs} genes={n_genes} r={r:.15f} "
              f"max_circ_diff={dmax:.3e}")

        ax = fig.add_subplot(gs[i, 0])
        ax.scatter(xy[order, 0], xy[order, 1], s=args.point_size, linewidths=0,
                   c=[palette[c] for c in ct[order]], rasterized=True)
        umap_axes(ax)
        ax.set_title(f"{label}\n{a.n_obs:,} cells: cell type", loc="left")
        lax = fig.add_subplot(gs[i, 1])
        lax.axis("off")
        present = sorted(set(ct))
        handles = [plt.Line2D([], [], marker="o", ls="", ms=5, color=palette[c])
                   for c in present]
        lax.legend(handles, present, loc="center left", frameon=False,
                   fontsize=7.5, ncol=1 if len(present) <= 10 else 2,
                   handletextpad=0.2, columnspacing=0.6, borderaxespad=0.0)

        ax = fig.add_subplot(gs[i, 2])
        sc = ax.scatter(xy[order, 0], xy[order, 1], s=args.point_size, linewidths=0,
                        c=py[order], cmap=cmap, norm=norm, rasterized=True)
        umap_axes(ax)
        ax.set_title("cell cycle position (Python)", loc="left")
        cax = fig.add_subplot(gs[i, 3])
        cax.set_box_aspect(12)
        cb = fig.colorbar(sc, cax=cax, ticks=[0, np.pi, TWO_PI])
        cb.ax.set_yticklabels(["0", "π", "2π"])
        cb.outline.set_visible(False)

        ax = fig.add_subplot(gs[i, 5])
        ax.scatter(rr[order], py[order], s=args.point_size, linewidths=0,
                   c=py[order], cmap=cmap, norm=norm, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(0, TWO_PI)
        ax.set_ylim(0, TWO_PI)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_ticks([0, TWO_PI])
            axis.set_ticklabels(["0", "2π"])
        ax.spines["left"].set_visible(True)
        ax.spines["bottom"].set_visible(True)
        ax.set_xlabel("R tricycle position (rad)")
        ax.set_ylabel("Python position (rad)")
        ax.set_title("Python against R", loc="left")
        ax.text(0.04, 0.96, f"circular r = {r:.6f}\nmax difference {dmax:.1e} rad\n{n_genes} genes",
                transform=ax.transAxes, va="top", ha="left", fontsize=9)

    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
