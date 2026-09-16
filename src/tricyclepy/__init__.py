"""tricycle-py: cell cycle position from single-cell RNA-seq, on AnnData.

A port of the R/Bioconductor package ``tricycle`` (Zheng, Hansen et al.) to
Python. The port is function by function against ``tricycle`` 1.12.0. The
reference embedding, the marker lists and the identifier maps are carried over
from the R package rather than regenerated, so the numbers match.

Orientation. R holds expression as genes by cells; ``AnnData`` holds it as
cells by genes. Every function here takes the ``AnnData`` orientation.

Input requirement. ``project_cycle_space`` and everything built on it need
log-transformed, library-size-normalised expression. Counts give wrong answers
without warning.

Cite Zheng SC et al., *Universal prediction of cell cycle position using
transfer learning*, Genome Biology 23:41 (2022),
doi:10.1186/s13059-021-02581-y. Cite Schwabe D et al., Molecular Systems
Biology 16:e9946 (2020) as well if you use :func:`estimate_schwabe_stage`.
"""

from __future__ import annotations

from .diagnostics import TotalUmiDiagnosis, diagnose_total_umi
from .loess import (
    LoessFit,
    PeriodicLoessResult,
    fit_periodic_loess,
    loess_fit,
)
from .plotting import (
    CYCLIC_COLORS,
    circle_scale_legend,
    circular_density,
    cyclic_colormap,
    plot_ccposition_density,
    plot_embedding_circle_scale,
)
from .projection import (
    estimate_cycle_position,
    project_cycle_space,
    theta_from_embedding,
)
from .refdata import (
    get_rotation,
    load_ensembl_to_symbol,
    load_go_cell_cycle_genes,
    load_neuro_ref,
    load_revelio_gene_list,
    provenance,
)
from .reference import CellCyclePCA, run_pca_cc_genes
from .schwabe import estimate_schwabe_stage

__version__ = "0.1.0"

#: The R package version this port was written against.
TRICYCLE_R_VERSION = "1.12.0"

__all__ = [
    "__version__",
    "TRICYCLE_R_VERSION",
    # core
    "project_cycle_space",
    "estimate_cycle_position",
    "theta_from_embedding",
    "estimate_schwabe_stage",
    # loess and diagnostics
    "loess_fit",
    "LoessFit",
    "fit_periodic_loess",
    "PeriodicLoessResult",
    "diagnose_total_umi",
    "TotalUmiDiagnosis",
    # reference building
    "run_pca_cc_genes",
    "CellCyclePCA",
    # shipped data
    "get_rotation",
    "load_neuro_ref",
    "load_revelio_gene_list",
    "load_ensembl_to_symbol",
    "load_go_cell_cycle_genes",
    "provenance",
    # plotting
    "CYCLIC_COLORS",
    "cyclic_colormap",
    "circular_density",
    "plot_embedding_circle_scale",
    "circle_scale_legend",
    "plot_ccposition_density",
]
