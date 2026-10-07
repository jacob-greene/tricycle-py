"""Build the README figure: one row per dataset, five panels per row.

Panels, left to right:

1. Cell type annotation: the UMAP coloured by cell type, with the palette
   stored in the input. Named cell types are labelled at their median UMAP
   position.
2. Cell cycle stage (binned): the UMAP coloured by stage bin of the position
   this package computes (G1/G0, S, G2/M).
3. The tricycle embedding, coloured by stage bin. Rays from the origin mark
   the bin edges; the position is the angle of each cell about the origin.
4. The UMAP coloured by the continuous position, with a colour-wheel legend.
5. Python position against R ``tricycle`` position for the same cells, with
   the circular correlation coefficient and lines at the bin edges. This
   column is drawn only when every input carries an R position.

Each column has one legend, shared by every row.

Each input is an ``.h5ad`` that holds:

* ``X``: log-normalised expression, cells by genes, gene symbols as
  ``var_names``;
* ``obsm["X_umap"]``: a two-dimensional embedding;
* ``obs[<celltype key>]``: a categorical cell type label, with its colours in
  ``uns["<celltype key>_colors"]`` in category order (the scanpy convention);
* optionally ``obs[<R key>]``: ``tricyclePosition`` from R ``tricycle`` on the
  same matrix (see ``docs/readme_figure_r.R``).

The Python position is computed here, from ``X``, by this package. Nothing is
read from a previous Python run. ``docs/make_readme_figure.py`` downloads the
public data and calls this script; use it to rebuild the README figure.

Usage::

    python docs/readme_figure.py \\
        --row "Human bone marrow" bm.h5ad --annotate HSC,Mono,Ery,NaiveB \\
        --row "Human CD34+ HSPC" cd34.h5ad --annotate HSC,Mono,Ery,CLP \\
        --out docs/readme_figure.png

One legend must be true for every row, so one cell type palette is used for
all rows: the palette of the row named by ``--palette-row`` (the first row by
default). That row must carry every cell type of every row.

Stage bins follow the tricycle vignette (Zheng et al. 2022, *Genome Biology*
23:41): position 0.5π is about the start of S, π the start of G2/M, 1.5π the
middle of M, and 1.75π to 0.25π is G1/G0. To cover the whole circle, the
bins here are

* G1/G0: 1.75π to 0.5π, through 0 (the vignette's G1/G0, plus 0.25π to 0.5π
  before the start of S);
* S: 0.5π to π;
* G2/M: π to 1.75π.

The bins are a guide. The position itself is continuous.

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
from matplotlib import patheffects  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

import tricyclepy as tp  # noqa: E402

TWO_PI = 2.0 * np.pi

#: Bin edges in radians: start of S, start of G2/M, start of G1/G0.
S_START, G2M_START, G1_START = 0.5 * np.pi, np.pi, 1.75 * np.pi
EDGES = ((S_START, "0.5π"), (G2M_START, "π"), (G1_START, "1.75π"))
STAGES = ("G1/G0", "S", "G2/M")
#: Middle of each bin; the bin colour is the cyclic colour at that angle.
STAGE_MID = {"G1/G0": 0.125 * np.pi, "S": 0.75 * np.pi, "G2/M": 1.375 * np.pi}


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


def stage_bin(theta) -> np.ndarray:
    """Name the stage bin of each position (radians on [0, 2*pi))."""
    theta = np.mod(np.asarray(theta, dtype=np.float64), TWO_PI)
    out = np.full(theta.shape, "G1/G0", dtype=object)
    out[(theta >= S_START) & (theta < G2M_START)] = "S"
    out[(theta >= G2M_START) & (theta < G1_START)] = "G2/M"
    return out


def style() -> None:
    rc = plt.rcParams
    rc["figure.dpi"] = 150
    rc["savefig.dpi"] = 150
    for side in ("top", "right", "bottom", "left"):
        rc[f"axes.spines.{side}"] = False
    rc["axes.grid"] = False
    rc["font.family"] = "sans-serif"
    rc["font.sans-serif"] = ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
    rc["font.size"] = 14
    rc["axes.titlesize"] = 15
    rc["legend.fontsize"] = 14
    rc["legend.title_fontsize"] = 15
    rc["pdf.fonttype"] = 42


def stored_palette(a, key: str) -> dict:
    """Cell type to colour, from ``uns[key + '_colors']`` in category order."""
    cats = list(a.obs[key].cat.categories)
    colours = list(a.uns[f"{key}_colors"])
    if len(cats) != len(colours):
        raise ValueError(f"{len(cats)} categories but {len(colours)} colours")
    return dict(zip(cats, colours))


def umap_axes(ax) -> None:
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def dot_legend(ax, labels, colours, title) -> None:
    handles = [plt.Line2D([], [], marker="o", ls="", ms=10, color=c) for c in colours]
    ax.legend(handles, labels, title=title, loc="center left", frameon=False,
              ncol=1, handletextpad=0.3, borderaxespad=0.0, labelspacing=0.35,
              alignment="left")


SUPTITLE = ("Validation: blood cell differentiation originates from a quiescent "
            "hematopoietic stem cell (HSC) population")


def umap_arrows(ax) -> None:
    """Small arrows in the lower left corner that name the UMAP axes."""
    o, L = 0.03, 0.16
    kw = dict(xycoords="axes fraction", textcoords="axes fraction",
              arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))
    ax.annotate("", xy=(o + L, o), xytext=(o, o), **kw)
    ax.annotate("", xy=(o, o + L), xytext=(o, o), **kw)
    ax.text(o + L + 0.01, o, "UMAP 1", transform=ax.transAxes, ha="left",
            va="center", fontsize=11)
    ax.text(o, o + L + 0.01, "UMAP 2", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=11, rotation=90)


def pad_lower_left(ax, xy, frac=0.12) -> None:
    """Widen the limits so the corner arrows sit clear of the points."""
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = hi - lo
    ax.set_xlim(lo[0] - frac * span[0], hi[0] + 0.02 * span[0])
    ax.set_ylim(lo[1] - frac * span[1], hi[1] + 0.02 * span[1])


def label_cell_types(ax, xy, ct, names) -> None:
    """Write each named cell type at the median UMAP position of its cells."""
    for name in names:
        hit = ct == name
        if not hit.any():
            raise ValueError(f"no cell has cell type {name!r}")
        x, y = np.median(xy[hit], axis=0)
        ax.text(x, y, name, ha="center", va="center", fontsize=14,
                fontweight="bold", zorder=5,
                path_effects=[patheffects.withStroke(linewidth=3.5, foreground="white")])


def colour_wheel(fig, ax, cmap) -> None:
    """A ring legend: the colour at each angle is the colour of that θ."""
    ax.axis("off")
    pax = ax.inset_axes([0.12, 0.34, 0.76, 0.32], projection="polar")
    t = np.linspace(0.0, TWO_PI, 361)
    r = np.array([0.62, 1.0])
    T, R = np.meshgrid(t, r)
    pax.pcolormesh(T, R, T[:-1, :-1], cmap=cmap, norm=Normalize(0.0, TWO_PI),
                   shading="flat", rasterized=True)
    pax.set_theta_zero_location("E")
    pax.set_theta_direction(1)
    pax.set_ylim(0.0, 1.0)
    pax.set_yticks([])
    pax.set_xticks([0.0, S_START, G2M_START, G1_START])
    pax.set_xticklabels(["0", "0.5π", "π", "1.75π"], fontsize=12)
    pax.tick_params(axis="x", pad=2)
    pax.grid(False)
    pax.spines["polar"].set_visible(False)
    pax.set_title("θ (rad)", fontsize=14, pad=16)


def main(argv=None) -> dict:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--row", nargs=2, action="append", required=True,
                   metavar=("LABEL", "H5AD"), help="row label and input file")
    p.add_argument("--annotate", action="append", default=[],
                   metavar="TYPE,TYPE,...",
                   help="cell types to label in column 1; give once per row, in row order")
    p.add_argument("--celltype-key", default="celltype")
    p.add_argument("--palette-row", type=int, default=0,
                   help="row whose stored cell type palette every row uses")
    p.add_argument("--r-key", default="tricyclePosition_R")
    p.add_argument("--species", default="human")
    p.add_argument("--gname-type", default="SYMBOL")
    p.add_argument("--point-size", type=float, default=2.0)
    p.add_argument("--seed", type=int, default=0, help="draw-order shuffle")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    if args.annotate and len(args.annotate) != len(args.row):
        p.error("give --annotate once per --row, or not at all")
    annotate = [[s for s in x.split(",") if s] for x in args.annotate] or [[] for _ in args.row]

    style()
    rows = []
    for label, path in args.row:
        a = ad.read_h5ad(path)
        tp.estimate_cycle_position(a, species=args.species, gname_type=args.gname_type)
        rows.append((label, a))
    with_r = all(args.r_key in a.obs for _, a in rows)
    if not with_r:
        print(f"obs[{args.r_key!r}] is missing from at least one input: "
              "column 5 (Python against R) is left out")

    key = args.celltype_key
    palette = stored_palette(rows[args.palette_row][1], key)
    for label, a in rows:
        own = stored_palette(a, key)
        missing = set(own) - set(palette)
        if missing:
            raise ValueError(f"{label}: no shared colour for {sorted(missing)}")
        differ = sorted(c for c in own if own[c] != palette[c])
        if differ:
            print(f"{label}: own stored colour replaced by the shared palette "
                  f"for {len(differ)} cell types")
    cmap = tp.cyclic_colormap()
    norm = Normalize(0.0, TWO_PI)
    stage_colour = {s: matplotlib.colors.rgb2hex(cmap(norm(m))) for s, m in STAGE_MID.items()}
    rng = np.random.default_rng(args.seed)

    n = len(rows)
    # Panel columns 0, 2, 4, 6 (and 8 with R); legend columns 1, 3 and 7, each
    # shared by every row; column 5 is a spacer.
    widths = [1.0, 0.42, 1.0, 0.42, 1.0, 0.10, 1.0, 0.62] + ([1.0] if with_r else [])
    fig = plt.figure(figsize=(sum(widths) * 3.85, 4.9 * n + 1.3))
    gs = fig.add_gridspec(n, len(widths), width_ratios=widths, wspace=0.12, hspace=0.12)
    present = set()
    stats = {}
    top_axes = {}

    for i, (label, a) in enumerate(rows):
        xy = np.asarray(a.obsm["X_umap"])
        ct = a.obs[key].astype(str).to_numpy()
        present |= set(ct)
        py = a.obs["tricyclePosition"].to_numpy(np.float64)
        emb = np.asarray(a.obsm["tricycleEmbedding"])[:, :2]
        # the angle drawn in panel 3 is the position itself
        assert np.allclose(circular_difference(np.arctan2(emb[:, 1], emb[:, 0]), py), 0.0)
        stage = stage_bin(py)
        order = rng.permutation(a.n_obs)
        n_genes = a.uns["tricycleEmbedding"]["n_genes"]
        counts = {s: int((stage == s).sum()) for s in STAGES}
        st = dict(cells=a.n_obs, genes=n_genes, bins=counts)
        if with_r:
            rr = a.obs[args.r_key].to_numpy(np.float64)
            st.update(r=circular_correlation(py, rr),
                      max_circ_diff=float(circular_difference(py, rr).max()),
                      bin_agreement=float(np.mean(stage == stage_bin(rr))))
        stats[label] = st
        print(label + ": " + " ".join(f"{k}={v}" for k, v in st.items()))
        stage_c = np.array([stage_colour[s] for s in stage], dtype=object)
        top = i == 0

        # 1. cell type, with the row label on its left
        ax = fig.add_subplot(gs[i, 0])
        ax.scatter(xy[order, 0], xy[order, 1], s=args.point_size, linewidths=0,
                   c=[palette[c] for c in ct[order]], rasterized=True)
        umap_axes(ax)
        pad_lower_left(ax, xy)
        umap_arrows(ax)
        label_cell_types(ax, xy, ct, annotate[i])
        ax.text(-0.06, 0.5, f"{label}\n{a.n_obs:,} cells", transform=ax.transAxes,
                rotation=90, ha="right", va="center", fontsize=15, fontweight="bold")
        if top:
            ax.set_title("Cell type annotation", loc="left")

        # 2. stage bin on the UMAP
        ax = fig.add_subplot(gs[i, 2])
        ax.scatter(xy[order, 0], xy[order, 1], s=args.point_size, linewidths=0,
                   c=list(stage_c[order]), rasterized=True)
        umap_axes(ax)
        pad_lower_left(ax, xy)
        if top:
            ax.set_title("Cell cycle stage (binned)", loc="left")

        # 3. tricycle embedding; θ is the angle about the origin
        ax = fig.add_subplot(gs[i, 4])
        if top:
            top_axes["embedding"] = ax
        ax.scatter(emb[order, 0], emb[order, 1], s=args.point_size, linewidths=0,
                   c=list(stage_c[order]), rasterized=True)
        lo = np.minimum(emb.min(axis=0), 0.0)
        hi = np.maximum(emb.max(axis=0), 0.0)
        pad = 0.10 * (hi - lo)
        lo, hi = lo - pad, hi + pad
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_aspect("equal")
        ax.axhline(0, color="0.8", lw=0.6, zorder=0)
        ax.axvline(0, color="0.8", lw=0.6, zorder=0)

        def to_edge(angle):
            # distance from the origin to the axes edge along this angle
            c, s = np.cos(angle), np.sin(angle)
            t = [(hi[0] if c > 0 else lo[0]) / c if abs(c) > 1e-12 else np.inf,
                 (hi[1] if s > 0 else lo[1]) / s if abs(s) > 1e-12 else np.inf]
            return 0.97 * min(t), c, s

        t, c, s = to_edge(0.0)
        ax.annotate("", xy=(t, 0), xytext=(0, 0), zorder=3,
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.4))
        ax.text(0.97 * hi[0], 0.04 * (hi[1] - lo[1]), "θ = 0", ha="right",
                va="bottom", fontsize=13, zorder=5,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))
        for angle, text in EDGES:
            t, c, s = to_edge(angle)
            ax.plot([0, t * c], [0, t * s], color="black", lw=1.2, ls="--", zorder=3)
            ax.text(0.80 * t * c, 0.80 * t * s, text, fontsize=13, zorder=4,
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))
        # counter-clockwise arrow: θ increases from embedding 1 towards embedding 2
        rad = 0.55 * min(hi[0], hi[1])
        ax.annotate("", xy=(rad * np.cos(1.2), rad * np.sin(1.2)), xytext=(rad, 0.0),
                    zorder=4, arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2,
                                              connectionstyle="arc3,rad=0.4"))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("tricycle embedding 1")
        ax.set_ylabel("tricycle embedding 2")

        # 4. continuous position on the UMAP
        ax = fig.add_subplot(gs[i, 6])
        if top:
            top_axes["theta"] = ax
        ax.scatter(xy[order, 0], xy[order, 1], s=args.point_size, linewidths=0,
                   c=py[order], cmap=cmap, norm=norm, rasterized=True)
        umap_axes(ax)
        pad_lower_left(ax, xy)

        # 5. Python against R
        if not with_r:
            continue
        ax = fig.add_subplot(gs[i, 8])
        ax.scatter(rr[order], py[order], s=args.point_size, linewidths=0,
                   c=py[order], cmap=cmap, norm=norm, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(0, TWO_PI)
        ax.set_ylim(0, TWO_PI)
        for angle, _ in EDGES:
            ax.axvline(angle, color="0.4", lw=0.9, ls="--", zorder=0)
            ax.axhline(angle, color="0.4", lw=0.9, ls="--", zorder=0)
        # 2π is the top-right corner; its label would collide with 1.75π
        ticks = [0.0] + [e for e, _ in EDGES]
        names = ["0"] + [t for _, t in EDGES]
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_ticks(ticks)
            axis.set_ticklabels(names, fontsize=12)
        # bin names inside the panel, along the bottom edge of each bin
        for name, b0, b1, rot in (("G1/G0", 0.0, S_START, 0), ("S", S_START, G2M_START, 0),
                                  ("G2/M", G2M_START, G1_START, 0),
                                  ("G1/G0", G1_START, TWO_PI, 90)):
            ax.text(0.5 * (b0 + b1), 0.12 if rot == 0 else 0.25, name, rotation=rot,
                    ha="center", va="bottom", fontsize=12, color="0.15")
        ax.spines["left"].set_visible(True)
        ax.spines["bottom"].set_visible(True)
        ax.set_xlabel("R tricycle θ")
        ax.set_ylabel("Python θ")
        if top:
            ax.set_title("Consistent results between Python\nand R tricycle implementations",
                         loc="left")
        ax.text(0.04, 0.96, f"circular r = {st['r']:.6f}\nmax diff {st['max_circ_diff']:.1e} rad",
                transform=ax.transAxes, va="top", ha="left", fontsize=12,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))

    # One title over columns 3 and 4, at the height of the other titles.
    l = top_axes["embedding"].get_position().x0
    r = top_axes["theta"].get_position().x1
    y = top_axes["theta"].get_position().y1
    fig.text(0.5 * (l + r), y + 0.012, "Cell cycle staging\nfrom a continuous variable (θ)",
             ha="center", va="bottom", fontsize=plt.rcParams["axes.titlesize"])
    fig.suptitle(SUPTITLE, fontsize=19, fontweight="bold",
                 y=y + 0.11 if n == 1 else y + 0.075)

    # One legend per column, spanning every row.
    lax = fig.add_subplot(gs[:, 1])
    lax.axis("off")
    cts = [c for c in palette if c in present]
    dot_legend(lax, cts, [palette[c] for c in cts], "cell type")

    sax = fig.add_subplot(gs[:, 3])
    sax.axis("off")
    dot_legend(sax, ["G1/G0\n1.75π to 0.5π", "S\n0.5π to π", "G2/M\nπ to 1.75π"],
               [stage_colour[s] for s in STAGES], "stage bin")

    colour_wheel(fig, fig.add_subplot(gs[:, 7]), cmap)

    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print("wrote", args.out)
    return stats


if __name__ == "__main__":
    main()
