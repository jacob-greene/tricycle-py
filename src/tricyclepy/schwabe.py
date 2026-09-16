"""Discrete five-stage cell cycle assignment, the Schwabe method.

Port of ``R/estimate_Schwabe_stage.R`` from R ``tricycle`` 1.12.0.

The assignment method is from Schwabe et al. 2020, not from the tricycle
authors. tricycle ships it for convenience and says so. Cite Schwabe et al. if
you use it.
"""

from __future__ import annotations

import logging
import warnings
from typing import Dict, Mapping, Optional, Sequence

import numpy as np
import scipy.sparse as sp

from .projection import _get_matrix, _is_anndata, _resolve_gene_names
from .refdata import load_ensembl_to_symbol, load_revelio_gene_list

__all__ = ["estimate_schwabe_stage"]

logger = logging.getLogger(__name__)


def _to_symbol(gname: np.ndarray, species: str,
               mapping: Optional[Mapping[str, str]]) -> np.ndarray:
    """R's ``.getSYMBOL``: map Ensembl to symbol, then upper-case the result."""
    table = load_ensembl_to_symbol(species) if mapping is None else mapping
    out = []
    for g in gname:
        s = table.get(g)
        out.append(None if s is None else s.upper())
    return np.asarray(out, dtype=object)


def _first_occurrence(names: np.ndarray) -> Dict[object, int]:
    first: Dict[object, int] = {}
    for i, name in enumerate(names):
        if name is not None and name not in first:
            first[name] = i
    return first


def _stage_score(block: np.ndarray, rows: np.ndarray, cor_thres: float,
                 stage: str, batch: object) -> np.ndarray:
    """One stage's raw score for one batch. ``block`` is genes by cells."""
    sub = block[rows]
    gene_sum = sub.sum(axis=1)
    keep = gene_sum > 0
    if int(keep.sum()) <= 3:
        raise ValueError(
            f"Less than 3 {stage} stage genes expressed in batch {batch}. "
            "Consider changing the batch setting."
        )
    sub = sub[keep]
    mean_v = sub.mean(axis=0)

    # Pearson correlation of each gene with the stage mean. A gene with no
    # variance within the batch gives NaN here, exactly as R's cor() gives NA.
    centred = sub - sub.mean(axis=1, keepdims=True)
    mc = mean_v - mean_v.mean()
    denom = np.sqrt((centred ** 2).sum(axis=1)) * np.sqrt((mc ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        cor_v = np.where(denom > 0, (centred @ mc) / denom, np.nan)

    selected = cor_v > cor_thres          # NaN compares False here...
    n_na = int(np.isnan(cor_v).sum())     # ...so count the NaNs separately.
    n_sel = int(selected.sum())
    if n_sel < 3:
        warnings.warn(
            f"Batch {batch} phase {stage} too few genes (<3).", RuntimeWarning,
            stacklevel=2,
        )
    logger.info("Batch %s phase %s gene: %d", batch, stage, n_sel)
    if n_na > 0 or n_sel == 0:
        # R indexes with a logical vector holding NA, which injects all-NA rows
        # and drives colMeans to NA for every cell in the batch. An empty
        # selection lands in the same place by a different route.
        return np.full(block.shape[1], np.nan)
    return sub[selected].mean(axis=0)


def _scale_columns(mat: np.ndarray) -> np.ndarray:
    """R's ``scale()``: centre by the column mean, divide by the column sd."""
    centred = mat - mat.mean(axis=0, keepdims=True)
    sd = mat.std(axis=0, ddof=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return centred / sd


def _assign(scores: np.ndarray, stage_names: Sequence[str],
            tolerance: float) -> np.ndarray:
    """Pick a stage per cell, or leave it unassigned."""
    n_stages = len(stage_names)
    out = np.full(scores.shape[0], None, dtype=object)
    for c in range(scores.shape[0]):
        s = scores[c]
        # R's order() sorts NA last whatever the direction, and breaks ties by
        # position. A stable sort of the negated finite values does both.
        finite = np.flatnonzero(~np.isnan(s))
        nan_idx = np.flatnonzero(np.isnan(s))
        order = np.concatenate([finite[np.argsort(-s[finite], kind="stable")], nan_idx])
        eta, eta2 = int(order[0]), int(order[1])
        gap = abs(eta - eta2)
        non_adjacent = 1 < gap < (n_stages - 1)
        diff = s[eta] - s[eta2]
        if non_adjacent:
            continue
        if np.isnan(diff):
            raise ValueError(
                "Stage scores contain NaN, so the tolerance test is undecidable. "
                "R raises 'missing value where TRUE/FALSE needed' here. This "
                "happens when a marker gene has no variance within a batch."
            )
        if diff < tolerance:
            continue
        out[c] = stage_names[eta]
    return out


def estimate_schwabe_stage(
    adata,
    *,
    layer: Optional[str] = None,
    batch: Optional[Sequence] = None,
    batch_key: Optional[str] = None,
    cycle_genes: Optional[Mapping[str, Sequence[str]]] = None,
    gname: Optional[Sequence] = None,
    gname_type: str = "ENSEMBL",
    species: str = "mouse",
    ensembl_to_symbol: Optional[Mapping[str, str]] = None,
    cor_thres: float = 0.2,
    tolerance: float = 0.3,
    key_added: str = "CCStage",
    copy: bool = False,
):
    """Assign each cell one of five discrete cell cycle stages.

    For each batch and stage the function averages the stage's marker genes,
    keeps only markers correlating above ``cor_thres`` with that average,
    turns the five stage averages into z-scores twice (across cells, then
    across stages within a cell), and takes the highest. A cell is left
    unassigned when its top two stages are not neighbours on the cycle, or
    when they are closer together than ``tolerance``.

    Parameters
    ----------
    adata
        ``AnnData``, or a cells-by-genes array or sparse matrix.
    layer
        Layer holding log-expression. ``None`` uses ``adata.X``.
    batch
        One batch label per cell. Stages are assigned within each batch.
    batch_key
        Column of ``adata.obs`` to use as ``batch``.
    cycle_genes
        Stage name to marker list. The five Revelio stages are used when this
        is ``None``. Stage order defines which stages count as neighbours.
    gname, gname_type, species, ensembl_to_symbol
        Gene identity, as in :func:`~tricyclepy.project_cycle_space`.
    cor_thres
        Correlation floor for keeping a marker.
    tolerance
        Smallest z-score gap between the top two stages that still assigns.
    key_added
        Column of ``adata.obs`` to write.
    copy
        Return a modified copy instead of writing in place.

    Returns
    -------
    numpy.ndarray or AnnData
        A bare-matrix input returns an object array of stage names with
        ``None`` for unassigned cells. An ``AnnData`` input returns ``None``,
        or the copy when ``copy=True``.

    Notes
    -----
    Identifier handling follows R exactly, asymmetries included. Mouse symbols
    are upper-cased; mouse Ensembl ids are mapped then upper-cased; human
    Ensembl ids are mapped then upper-cased; human symbols are used as given.
    """
    names = _resolve_gene_names(adata, gname, layer)
    mat, _ = _get_matrix(adata, layer)

    if batch_key is not None:
        if batch is not None:
            raise ValueError("pass batch= or batch_key=, not both")
        if not _is_anndata(adata):
            raise ValueError("batch_key= requires an AnnData input")
        batch = np.asarray(adata.obs[batch_key])
    if batch is None:
        batch = np.ones(mat.shape[0], dtype=int)
    batch = np.asarray(batch)
    if batch.shape[0] != mat.shape[0]:
        raise ValueError("batch has the wrong length")
    if np.any([b is None or (isinstance(b, float) and np.isnan(b)) for b in batch]):
        raise ValueError("batch may not contain missing values")

    if cycle_genes is None:
        if species == "mouse":
            names = (np.asarray([n.upper() for n in names], dtype=object)
                     if gname_type == "SYMBOL"
                     else _to_symbol(names, "mouse", ensembl_to_symbol))
        elif gname_type == "ENSEMBL":
            names = _to_symbol(names, "human", ensembl_to_symbol)
        cycle_genes = load_revelio_gene_list()

    # R: lapply(cycleGene.l, intersect, rownames(data.m)) then Reduce(union, .)
    first = _first_occurrence(names)
    per_stage: Dict[str, list] = {}
    for stage, genes in cycle_genes.items():
        seen = set()
        kept = []
        for g in genes:
            if g in first and g not in seen:
                seen.add(g)
                kept.append(g)
        per_stage[stage] = kept
    all_genes: list = []
    seen_all = set()
    for genes in per_stage.values():
        for g in genes:
            if g not in seen_all:
                seen_all.add(g)
                all_genes.append(g)
    if len(all_genes) < 30:
        raise ValueError(
            "Less than 30 cell cycle genes found. Not enough information to "
            "assign the 5 stages."
        )

    rows = np.asarray([first[g] for g in all_genes], dtype=np.intp)
    position = {g: k for k, g in enumerate(all_genes)}
    stage_rows = {s: np.asarray([position[g] for g in genes], dtype=np.intp)
                  for s, genes in per_stage.items()}

    sub = mat[:, rows]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    sub = np.asarray(sub, dtype=np.float64).T          # genes by cells, as in R

    stage_names = list(cycle_genes.keys())
    assigned = np.full(mat.shape[0], None, dtype=object)
    _, first_seen = np.unique(batch, return_index=True)
    for b in batch[np.sort(first_seen)]:
        idx = np.flatnonzero(batch == b)
        block = sub[:, idx]
        raw = np.column_stack([
            _stage_score(block, stage_rows[s], cor_thres, s, b) for s in stage_names
        ])
        scores = _scale_columns(_scale_columns(raw).T).T
        assigned[idx] = _assign(scores, stage_names, tolerance)

    if not _is_anndata(adata):
        return assigned

    import pandas as pd

    out = adata.copy() if copy else adata
    out.obs[key_added] = pd.Categorical(assigned, categories=stage_names)
    return out if copy else None
