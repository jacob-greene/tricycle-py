"""Projection and angle. The expectations are computed independently here."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import tricyclepy as tp

RNG = np.random.default_rng(20260916)


def _toy(n_cells=40, genes=("GA", "GB", "GC", "GD")):
    x = RNG.random((n_cells, len(genes))) * 4.0
    return ad.AnnData(X=x, var=pd.DataFrame(index=list(genes)),
                      obs=pd.DataFrame(index=[f"c{i}" for i in range(n_cells)]))


def _reference(genes, rotation):
    return np.asarray(genes, dtype=object), np.asarray(rotation, dtype=float)


def _expected(x, cols, rotation):
    """The definition, written out: centre each gene, then weight and sum."""
    block = np.asarray(x, dtype=float)[:, cols]
    return (block - block.mean(axis=0)) @ rotation


def test_projection_matches_the_written_out_definition():
    adata = _toy()
    ref = _reference(["GC", "GA"], [[0.5, -1.5], [2.0, 0.25]])
    out = tp.project_cycle_space(adata, ref=ref, copy=True)
    # Gene order follows the data, not the reference: GA is column 0, GC is 2.
    want = _expected(adata.X, [0, 2], np.array([[2.0, 0.25], [0.5, -1.5]]))
    assert np.allclose(out.obsm["tricycleEmbedding"], want, atol=0, rtol=1e-14)
    assert list(out.uns["tricycleEmbedding"]["genes"]) == ["GA", "GC"]
    assert out.uns["tricycleEmbedding"]["n_genes"] == 2


def test_gene_order_follows_the_data_not_the_reference():
    adata = _toy(genes=("GD", "GC", "GB", "GA"))
    ref = _reference(["GA", "GB", "GC", "GD"], np.arange(8.0).reshape(4, 2))
    out = tp.project_cycle_space(adata, ref=ref, copy=True)
    assert list(out.uns["tricycleEmbedding"]["genes"]) == ["GD", "GC", "GB", "GA"]


def test_duplicate_gene_names_take_the_first_occurrence():
    """R's `data.m[gene, ]` takes the first row carrying that name.

    The duplicate column must vary across cells. A constant column centres to
    exactly zero, so counting it twice would add nothing and the test would
    pass whether or not the rule is implemented.
    """
    adata = _toy(genes=("GA", "GB", "GA"))
    adata.X[:, 2] = RNG.random(adata.n_obs) * 50.0 + 10.0
    assert adata.X[:, 2].std() > 1.0           # or the test proves nothing
    ref = _reference(["GA"], [[1.0, 0.0]])
    out = tp.project_cycle_space(adata, ref=ref, copy=True)
    want = _expected(adata.X, [0], np.array([[1.0, 0.0]]))
    got = out.obsm["tricycleEmbedding"]
    assert np.allclose(got, want)
    # And it must not be the answer you would get by using both copies.
    both = _expected(adata.X, [0, 2], np.array([[1.0, 0.0], [1.0, 0.0]]))
    assert not np.allclose(got, both)


def test_duplicate_reference_names_take_the_first_reference_row():
    """The same first-occurrence rule applies on the reference side."""
    adata = _toy(genes=("GA", "GB"))
    ref = _reference(["GA", "GA"], [[1.0, 0.0], [0.0, 7.0]])
    out = tp.project_cycle_space(adata, ref=ref, copy=True)
    want = _expected(adata.X, [0], np.array([[1.0, 0.0]]))
    got = out.obsm["tricycleEmbedding"]
    assert np.allclose(got, want)
    last = _expected(adata.X, [0], np.array([[0.0, 7.0]]))
    assert not np.allclose(got, last)


def test_unmatched_genes_are_dropped_and_missing_names_are_skipped():
    adata = _toy(genes=("GA", "GB"))
    names = np.asarray([None, "GB"], dtype=object)
    ref = _reference(["GA", "GB"], [[1.0, 1.0], [2.0, -2.0]])
    out = tp.project_cycle_space(adata, ref=ref, gname=names, copy=True)
    want = _expected(adata.X, [1], np.array([[2.0, -2.0]]))
    assert np.allclose(out.obsm["tricycleEmbedding"], want)


def test_no_shared_genes_raises_with_a_useful_message():
    adata = _toy()
    with pytest.raises(ValueError, match="No reference genes"):
        tp.project_cycle_space(adata, ref=_reference(["ZZZ"], [[1.0, 0.0]]))


def test_sparse_and_dense_agree_and_chunking_does_not_matter():
    adata = _toy(n_cells=200)
    adata.X[adata.X < 2.0] = 0.0
    ref = _reference(["GA", "GB", "GC"], RNG.normal(size=(3, 2)))
    dense = tp.project_cycle_space(adata, ref=ref, copy=True)
    sparse_adata = ad.AnnData(X=sp.csr_matrix(adata.X), var=adata.var, obs=adata.obs)
    sparse = tp.project_cycle_space(sparse_adata, ref=ref, copy=True)
    chunked = tp.project_cycle_space(sparse_adata, ref=ref, chunk_size=7, copy=True)
    a = dense.obsm["tricycleEmbedding"]
    b = sparse.obsm["tricycleEmbedding"]
    c = chunked.obsm["tricycleEmbedding"]
    assert np.allclose(a, b, atol=1e-13, rtol=0)
    assert np.allclose(b, c, atol=1e-13, rtol=0)


def test_layer_is_used_when_named():
    adata = _toy()
    adata.layers["logcounts"] = adata.X * 3.0
    ref = _reference(["GA"], [[1.0, 0.0]])
    out = tp.project_cycle_space(adata, ref=ref, layer="logcounts", copy=True)
    want = _expected(adata.layers["logcounts"], [0], np.array([[1.0, 0.0]]))
    assert np.allclose(out.obsm["tricycleEmbedding"], want)


def test_bare_matrix_requires_gene_names():
    with pytest.raises(ValueError, match="gname"):
        tp.project_cycle_space(RNG.random((5, 3)), ref=_reference(["A"], [[1, 0]]))


def test_gname_length_is_checked():
    adata = _toy()
    with pytest.raises(ValueError, match="entries"):
        tp.project_cycle_space(adata, gname=["A", "B"],
                               ref=_reference(["A"], [[1.0, 0.0]]))


# --- the angle -------------------------------------------------------------


@pytest.mark.parametrize(
    "pc1,pc2,expected",
    [
        (1.0, 0.0, 0.0),                    # PC1 is the zero axis
        (0.0, 1.0, np.pi / 2),              # counter-clockwise
        (-1.0, 0.0, np.pi),
        (0.0, -1.0, 3 * np.pi / 2),         # wrapped onto [0, 2*pi)
        (1.0, -1.0, 7 * np.pi / 4),
        (0.0, 0.0, 0.0),                    # atan2(0, 0) is 0, as in R
    ],
)
def test_angle_convention(pc1, pc2, expected):
    got = tp.theta_from_embedding(np.array([[pc1, pc2]]))
    assert got[0] == pytest.approx(expected, abs=1e-15)


def test_angle_is_always_in_range():
    emb = RNG.normal(size=(500, 2))
    theta = tp.theta_from_embedding(emb)
    assert theta.min() >= 0.0
    assert theta.max() < 2 * np.pi


def test_angle_centre_shifts_the_origin():
    emb = np.array([[2.0, 3.0]])
    shifted = tp.theta_from_embedding(emb, center_pc1=2.0, center_pc2=2.0)
    assert shifted[0] == pytest.approx(np.pi / 2)


def test_estimate_cycle_position_reuses_an_existing_embedding():
    adata = _toy()
    adata.obsm["tricycleEmbedding"] = np.tile([0.0, 1.0], (adata.n_obs, 1))
    tp.estimate_cycle_position(adata)
    assert np.allclose(adata.obs["tricyclePosition"], np.pi / 2)


def test_estimate_cycle_position_projects_when_no_embedding_exists():
    adata = _toy()
    ref = _reference(["GA", "GB"], [[1.0, 0.0], [0.0, 1.0]])
    tp.estimate_cycle_position(adata, ref=ref)
    assert "tricycleEmbedding" in adata.obsm
    want = tp.theta_from_embedding(adata.obsm["tricycleEmbedding"])
    assert np.allclose(adata.obs["tricyclePosition"], want)


def test_bare_matrix_round_trip_gives_the_angle_directly():
    x = RNG.random((30, 2)) * 5
    ref = _reference(["A", "B"], [[1.0, 0.5], [-0.5, 1.0]])
    theta = tp.estimate_cycle_position(x, gname=["A", "B"], ref=ref)
    emb = _expected(x, [0, 1], ref[1])
    assert np.allclose(theta, np.mod(np.arctan2(emb[:, 1], emb[:, 0]), 2 * np.pi))


def test_human_ensembl_path_maps_through_symbols():
    """A human Ensembl input must reach the reference through its symbol."""
    genes = ["ENSG00000131747", "ENSG00000012048"]     # TOP2A, BRCA1
    adata = ad.AnnData(X=RNG.random((25, 2)) * 3,
                       var=pd.DataFrame(index=genes))
    tp.project_cycle_space(adata, species="human", gname_type="ENSEMBL")
    # TOP2A is in neuroRef, BRCA1 is not.
    assert list(adata.uns["tricycleEmbedding"]["genes"]) == ["TOP2A"]


def test_copy_leaves_the_input_untouched():
    adata = _toy()
    ref = _reference(["GA"], [[1.0, 0.0]])
    out = tp.project_cycle_space(adata, ref=ref, copy=True)
    assert "tricycleEmbedding" in out.obsm
    assert "tricycleEmbedding" not in adata.obsm
