"""Projection into the cell cycle space, and the cell cycle position angle.

Port of ``R/project_cycle_space.R`` and ``R/estimate_cycle_position.R`` from the
R package ``tricycle`` 1.12.0.

Orientation warning. R holds expression as genes by cells. ``AnnData`` holds it
as cells by genes. Every public function here takes the ``AnnData`` orientation.
"""

from __future__ import annotations

import logging
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp

from .refdata import get_rotation, map_ensembl_to_symbol

__all__ = ["project_cycle_space", "estimate_cycle_position", "theta_from_embedding"]

logger = logging.getLogger(__name__)

#: Cells per block when centring and projecting. Each block is centred densely,
#: exactly as R does, and every cell's result depends only on its own row. The
#: block size therefore does not change what is computed, but it can change how
#: BLAS associates the sums inside a row: measured at 1.7e-14 radians between
#: blocks of 20,000 and 1,000 on 8,627 cells, which is floating-point noise.
DEFAULT_CHUNK_SIZE = 20_000

Matrixlike = Union["AnnData", np.ndarray, sp.spmatrix]  # noqa: F821


def _is_anndata(x) -> bool:
    return hasattr(x, "obs") and hasattr(x, "var") and hasattr(x, "X")


def _get_matrix(x, layer: Optional[str]):
    """Return the cells-by-genes matrix and the gene names for ``x``."""
    if _is_anndata(x):
        mat = x.X if layer is None else x.layers[layer]
        return mat, np.asarray(x.var_names, dtype=object)
    if layer is not None:
        raise ValueError("layer= is only meaningful for an AnnData input")
    return x, None


def _resolve_gene_names(x, gname, layer) -> np.ndarray:
    mat, var_names = _get_matrix(x, layer)
    n_genes = mat.shape[1]
    if gname is not None:
        gname = np.asarray(gname, dtype=object)
        if gname.shape[0] != n_genes:
            raise ValueError(
                f"gname has {gname.shape[0]} entries but the matrix has {n_genes} genes"
            )
        return gname
    if var_names is None:
        raise ValueError("gname= is required when the input is a bare matrix")
    return var_names


def _match_genes(
    data_names: np.ndarray, ref_names: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Reproduce R's ``data.m[intersect(rownames(data.m), rownames(ref.m)), ]``.

    R's ``intersect`` keeps the first occurrence of each name, in the order the
    names appear in the data. Character subsetting of a matrix then takes the
    first row carrying that name. Both rules are reproduced here, so duplicated
    gene names give the same answer as R rather than a silently different one.
    """
    ref_first: Dict[object, int] = {}
    for j, name in enumerate(ref_names):
        if name is not None and name not in ref_first:
            ref_first[name] = j
    data_idx = []
    ref_idx = []
    seen = set()
    for i, name in enumerate(data_names):
        if name is None or name in seen:
            continue
        j = ref_first.get(name)
        if j is None:
            continue
        seen.add(name)
        data_idx.append(i)
        ref_idx.append(j)
    return np.asarray(data_idx, dtype=np.intp), np.asarray(ref_idx, dtype=np.intp)


def _column_means(mat, cols: np.ndarray) -> np.ndarray:
    """Mean of each selected column, over all cells. Never densifies the input."""
    n_cells = mat.shape[0]
    if sp.issparse(mat):
        sub = mat[:, cols]
        total = np.asarray(sub.sum(axis=0), dtype=np.float64).ravel()
        return total / n_cells
    return np.asarray(mat[:, cols], dtype=np.float64).mean(axis=0)


def _project(mat, cols: np.ndarray, rotation: np.ndarray,
             chunk_size: int) -> np.ndarray:
    """Centre the selected genes and multiply by the rotation, block by block."""
    n_cells = mat.shape[0]
    means = _column_means(mat, cols)
    out = np.empty((n_cells, rotation.shape[1]), dtype=np.float64)
    is_sparse = sp.issparse(mat)
    csr = mat.tocsr() if is_sparse and not sp.isspmatrix_csr(mat) else mat
    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        block = csr[start:stop, :][:, cols]
        block = block.toarray() if sp.issparse(block) else np.asarray(block)
        block = block.astype(np.float64, copy=False) - means
        out[start:stop] = block @ rotation
    return out


def project_cycle_space(
    adata,
    *,
    layer: Optional[str] = None,
    ref: Optional[Tuple[Sequence, np.ndarray]] = None,
    gname: Optional[Sequence] = None,
    gname_type: str = "ENSEMBL",
    species: str = "mouse",
    ensembl_to_symbol: Optional[Mapping[str, str]] = None,
    key_added: str = "tricycleEmbedding",
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    copy: bool = False,
):
    """Project log-expression into the pre-learned cell cycle embedding.

    The projection is the sum of gene-mean-centred log-expression weighted by
    the reference rotation. The input must be log-transformed and library-size
    normalised. Centring uses the mean over all cells in ``adata``, so the
    embedding of one cell depends on which other cells are present.

    Parameters
    ----------
    adata
        ``AnnData`` with cells in rows, or a cells-by-genes array or sparse
        matrix. A bare matrix requires ``gname``.
    layer
        Layer holding log-expression. ``None`` uses ``adata.X``.
    ref
        Custom reference as ``(gene_names, rotation)``, where ``rotation`` has
        one row per gene name and at least two columns. Supplying this ignores
        ``species`` and ``gname_type``, matching R's ``ref.m``.
    gname
        Gene names to match on, one per column. Overrides ``var_names``.
    gname_type
        ``"ENSEMBL"`` or ``"SYMBOL"``. Describes ``gname``.
    species
        ``"mouse"`` or ``"human"``.
    ensembl_to_symbol
        Replacement Ensembl-to-symbol map for the human Ensembl path. The
        shipped ``org.Hs.eg.db`` map is used when this is ``None``.
    key_added
        Key in ``adata.obsm`` to write. Ignored for a bare matrix.
    chunk_size
        Cells per block. Lower it to cut peak memory. The result changes only
        by floating-point re-association, far below any meaningful tolerance.
    copy
        Return a modified copy instead of writing in place.

    Returns
    -------
    numpy.ndarray or AnnData
        A bare-matrix input returns the ``(n_cells, 2)`` projection. An
        ``AnnData`` input returns ``None``, or the copy when ``copy=True``, with
        the projection in ``.obsm[key_added]`` and the rotation actually used in
        ``.uns[key_added]``.

    Notes
    -----
    Ported from ``.project_cycle_space`` and ``.calProjection`` in R
    ``tricycle`` 1.12.0.
    """
    names = _resolve_gene_names(adata, gname, layer)
    mat, _ = _get_matrix(adata, layer)

    if ref is None:
        ref_names, rotation = get_rotation(gname_type=gname_type, species=species)
        if species == "human" and gname_type == "ENSEMBL":
            logger.info(
                "Reference genes are mouse symbols upper-cased; mapping the human "
                "Ensembl identifiers to symbols before matching."
            )
            names = map_ensembl_to_symbol(
                names, species="human", mapping=ensembl_to_symbol
            )
    else:
        ref_names, rotation = ref
        ref_names = np.asarray(ref_names, dtype=object)
        rotation = np.asarray(rotation, dtype=np.float64)
        if rotation.ndim != 2 or rotation.shape[0] != ref_names.shape[0]:
            raise ValueError("ref rotation must have one row per reference gene name")

    data_idx, ref_idx = _match_genes(names, ref_names)
    if data_idx.size == 0:
        raise ValueError(
            "No reference genes found in the data. Check gname_type and species: "
            "the internal reference is keyed on mouse Ensembl ids, mouse symbols, "
            "or upper-cased symbols for human."
        )
    logger.info("The number of projection genes found in the new data is %d.",
                data_idx.size)

    projection = _project(mat, data_idx, rotation[ref_idx], chunk_size)

    if not _is_anndata(adata):
        return projection

    out = adata.copy() if copy else adata
    out.obsm[key_added] = projection
    out.uns[key_added] = {
        "genes": np.asarray(names, dtype=object)[data_idx],
        "rotation": rotation[ref_idx],
        "n_genes": int(data_idx.size),
    }
    return out if copy else None


def theta_from_embedding(embedding: np.ndarray, center_pc1: float = 0.0,
                         center_pc2: float = 0.0) -> np.ndarray:
    """Turn a two-column embedding into an angle on ``[0, 2*pi)``.

    This is R's ``.getTheta``, which calls ``circular::coord2rad``. That
    function is ``atan2(pc2, pc1)`` reduced modulo ``2*pi``, so PC1 is the zero
    axis and the angle runs counter-clockwise.
    """
    emb = np.asarray(embedding, dtype=np.float64)
    if emb.ndim != 2 or emb.shape[1] < 2:
        raise ValueError("embedding must have at least two columns")
    pc1 = emb[:, 0] - center_pc1
    pc2 = emb[:, 1] - center_pc2
    return np.mod(np.arctan2(pc2, pc1), 2.0 * np.pi)


def estimate_cycle_position(
    adata,
    *,
    layer: Optional[str] = None,
    dimred: str = "tricycleEmbedding",
    center_pc1: float = 0.0,
    center_pc2: float = 0.0,
    key_added: str = "tricyclePosition",
    copy: bool = False,
    **kwargs,
):
    """Assign each cell a cell cycle position in radians on ``[0, 2*pi)``.

    Reuses ``adata.obsm[dimred]`` when it exists, and otherwise projects first.
    That is R's behaviour and it matters: a stale embedding is reused silently
    on both sides.

    Parameters
    ----------
    adata
        ``AnnData``, or a cells-by-genes array or sparse matrix.
    layer, kwargs
        Passed to :func:`project_cycle_space` when a projection is needed.
    dimred
        Key in ``adata.obsm`` holding the embedding.
    center_pc1, center_pc2
        Origin of the angle in embedding coordinates.
    key_added
        Column in ``adata.obs`` to write.
    copy
        Return a modified copy instead of writing in place.

    Returns
    -------
    numpy.ndarray or AnnData
        A bare-matrix input returns the angle per cell. An ``AnnData`` input
        returns ``None``, or the copy when ``copy=True``.

    Notes
    -----
    tricycle's own reading of the scale: ``0.5*pi`` near the start of S,
    ``pi`` near the start of G2M, ``1.5*pi`` near mid M, and ``1.75*pi`` to
    ``0.25*pi`` G1/G0. The value is continuous; the stage names are a guide.
    """
    if not _is_anndata(adata):
        projection = project_cycle_space(adata, layer=layer, **kwargs)
        return theta_from_embedding(projection, center_pc1, center_pc2)

    out = adata.copy() if copy else adata
    if dimred not in out.obsm:
        logger.info("%s not found in .obsm; running project_cycle_space.", dimred)
        project_cycle_space(out, layer=layer, key_added=dimred, **kwargs)
    out.obs[key_added] = theta_from_embedding(
        np.asarray(out.obsm[dimred]), center_pc1, center_pc2
    )
    return out if copy else None
