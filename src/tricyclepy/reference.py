"""Building a new reference projection from Gene Ontology cell cycle genes.

Port of ``R/run_pca_cc_genes.R`` from R ``tricycle`` 1.12.0, which wraps
``scater::runPCA`` on the genes annotated to GO:0007049.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp

from .projection import _get_matrix, _is_anndata, _resolve_gene_names
from .refdata import load_go_cell_cycle_genes

__all__ = ["run_pca_cc_genes", "CellCyclePCA"]

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CellCyclePCA:
    """PCA of the GO cell cycle genes, and the reference it yields."""

    coords: np.ndarray          #: (n_cells, n_components)
    rotation: np.ndarray        #: (n_top_genes, n_components)
    genes: np.ndarray           #: gene names, one per rotation row
    percent_var: np.ndarray     #: percent of total variance per component
    n_cc_genes: int             #: GO cell cycle genes present in the input

    def as_reference(self, n_components: int = 2):
        """Return ``(gene_names, rotation)`` ready for ``ref=`` in a projection."""
        return self.genes, self.rotation[:, :n_components]


def run_pca_cc_genes(
    adata,
    *,
    layer: Optional[str] = None,
    gname: Optional[Sequence] = None,
    gname_type: str = "ENSEMBL",
    species: str = "mouse",
    go_genes: Optional[Sequence[str]] = None,
    ntop: int = 500,
    n_components: int = 20,
) -> CellCyclePCA:
    """Run PCA on the Gene Ontology cell cycle genes to learn a new reference.

    The input must be library-size normalised log-expression. The gene set is
    GO:0007049 as shipped with the package. Within it, the ``ntop`` genes with
    the highest variance across cells are kept, then centred and decomposed.
    Variables are centred but not scaled, which is what ``scater::runPCA``
    does by default.

    Parameters
    ----------
    adata
        ``AnnData``, or a cells-by-genes array or sparse matrix.
    layer
        Layer holding log-expression. ``None`` uses ``adata.X``.
    gname, gname_type, species
        Gene identity, as in :func:`~tricyclepy.project_cycle_space`.
    go_genes
        Replacement cell cycle gene set.
    ntop
        Number of highest-variance genes to decompose.
    n_components
        Number of components to return.

    Returns
    -------
    CellCyclePCA

    Notes
    -----
    Component signs are arbitrary. A singular value decomposition may return
    a component and its negation with equal right. Two runs that differ only
    in sign describe the same reference, but they give cell cycle positions
    that differ by a reflection, a rotation, or both. Fix the sign before you
    compare positions from two separately learned references.
    """
    names = _resolve_gene_names(adata, gname, layer)
    mat, _ = _get_matrix(adata, layer)

    cc = (load_go_cell_cycle_genes(species=species, gname_type=gname_type)
          if go_genes is None else np.asarray(go_genes, dtype=object))
    cc_set = {g for g in cc if g is not None}
    keep = np.flatnonzero(np.asarray([n in cc_set for n in names]))
    if keep.size < 100:
        raise ValueError(
            f"Only {keep.size} Gene Ontology cell cycle genes found in the data. "
            "Check the data, gname and gname_type."
        )
    logger.info("%d out of %d Gene Ontology cell cycle genes found in your data.",
                keep.size, len(cc))

    sub = mat[:, keep]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    sub = np.asarray(sub, dtype=np.float64)

    variance = sub.var(axis=0, ddof=1)
    n_top = int(min(ntop, keep.size))
    # Highest variance first, ties broken by gene order, as R's order() does.
    # scater keeps that order in the result, so the rotation rows come out
    # ranked by variance rather than in gene order. Re-sorting here would give
    # the same subspace under a different row order, which is a silent mismatch
    # against R rather than a numerical one.
    top = np.argsort(-variance, kind="stable")[:n_top]
    block = sub[:, top]
    genes = np.asarray(names, dtype=object)[keep][top]

    centred = block - block.mean(axis=0, keepdims=True)
    n_comp = int(min(n_components, min(centred.shape)))
    _, s, vt = np.linalg.svd(centred, full_matrices=False)
    rotation = vt[:n_comp].T
    coords = centred @ rotation
    total_var = float((s ** 2).sum())
    percent_var = (s[:n_comp] ** 2) / total_var * 100.0 if total_var > 0 else s[:n_comp]
    return CellCyclePCA(coords=coords, rotation=rotation, genes=genes,
                        percent_var=percent_var, n_cc_genes=int(keep.size))
