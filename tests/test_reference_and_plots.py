"""Reference building, the UMI diagnostic, and the circular density."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tricyclepy as tp

RNG = np.random.default_rng(7)


def _go_backed_adata(n_cells=300, n_genes=260):
    """Cells over real GO:0007049 human symbols, so the gene filter has work."""
    go = [g for g in dict.fromkeys(tp.load_go_cell_cycle_genes("human", "SYMBOL"))
          if g and g != "NA"][:n_genes]
    extra = [f"NOTACYCLEGENE{i}" for i in range(40)]
    x = RNG.normal(2.0, 1.0, (n_cells, len(go) + len(extra)))
    return ad.AnnData(X=x, var=pd.DataFrame(index=go + extra)), go


def test_run_pca_cc_genes_keeps_only_cell_cycle_genes():
    adata, go = _go_backed_adata()
    pca = tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL",
                              ntop=50, n_components=5)
    assert pca.n_cc_genes == len(go)
    assert set(pca.genes) <= set(go)
    assert pca.genes.size == 50
    assert pca.coords.shape == (adata.n_obs, 5)
    assert pca.rotation.shape == (50, 5)


def test_run_pca_cc_genes_selects_the_highest_variance_genes():
    adata, go = _go_backed_adata()
    adata.X[:, :5] *= 10.0                         # make five genes dominate
    pca = tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL",
                              ntop=5, n_components=2)
    assert set(pca.genes) == set(go[:5])


def test_run_pca_cc_genes_returns_genes_ranked_by_variance():
    """scater keeps the variance ranking, so the rotation rows must too.

    Re-sorting the rows into gene order gives the same subspace under a
    different row order. That is a silent mismatch against R, not a numerical
    one, so it is pinned here.
    """
    adata, _ = _go_backed_adata()
    pca = tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL",
                              ntop=25, n_components=3)
    cols = adata.var_names.get_indexer(pca.genes)
    variance = np.asarray(adata.X)[:, cols].var(axis=0, ddof=1)
    assert np.all(np.diff(variance) <= 0)


def test_run_pca_cc_genes_matches_a_plain_svd_on_the_same_genes():
    adata, _ = _go_backed_adata()
    pca = tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL",
                              ntop=30, n_components=4)
    cols = adata.var_names.get_indexer(pca.genes)
    block = np.asarray(adata.X)[:, cols]
    centred = block - block.mean(axis=0)
    _, s, vt = np.linalg.svd(centred, full_matrices=False)
    assert np.allclose(np.abs(pca.rotation), np.abs(vt[:4].T), atol=1e-10)
    assert np.allclose(pca.coords, centred @ pca.rotation, atol=1e-10)
    assert pca.percent_var[0] == pytest.approx(s[0] ** 2 / (s ** 2).sum() * 100)


def test_run_pca_cc_genes_refuses_a_dataset_with_too_few_cycle_genes():
    adata = ad.AnnData(X=RNG.normal(size=(20, 5)),
                       var=pd.DataFrame(index=[f"X{i}" for i in range(5)]))
    with pytest.raises(ValueError, match="Gene Ontology cell cycle genes"):
        tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL")


def test_as_reference_round_trips_into_a_projection():
    adata, _ = _go_backed_adata()
    pca = tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL",
                              ntop=40, n_components=3)
    theta = tp.estimate_cycle_position(adata, ref=pca.as_reference(2), copy=True)
    assert theta.obs["tricyclePosition"].between(0, 2 * np.pi).all()


# --- diagnostics -----------------------------------------------------------


def test_diagnose_total_umi_passes_a_cycling_profile():
    theta = np.linspace(0, 2 * np.pi, 800, endpoint=False)
    # Peak near pi, trough near 1.6*pi, which is the shape it looks for.
    umis = 2 ** (12.0 + 0.8 * np.cos(theta - np.pi))
    result = tp.diagnose_total_umi(theta, umis)
    assert result.passed
    assert result.difference > 0.4


def test_diagnose_total_umi_warns_on_a_flat_profile():
    theta = np.linspace(0, 2 * np.pi, 600, endpoint=False)
    umis = np.full(600, 4000.0)
    with pytest.warns(RuntimeWarning, match="below 0.4"):
        result = tp.diagnose_total_umi(theta, umis)
    assert not result.passed


def test_diagnose_total_umi_warns_when_cells_pile_up_off_g1():
    theta = RNG.uniform(0.6 * np.pi, 1.7 * np.pi, 400)
    umis = 2 ** (12.0 + 0.8 * np.cos(theta - np.pi))
    with pytest.warns(RuntimeWarning, match="50 percent"):
        tp.diagnose_total_umi(theta, umis)


# --- circular density ------------------------------------------------------


def test_circular_density_integrates_to_one():
    theta = RNG.uniform(0, 2 * np.pi, 400)
    x, y = tp.circular_density(theta, bw=30.0, n=4096)
    assert np.trapezoid(y, x) == pytest.approx(1.0, rel=1e-6)


def test_circular_density_of_one_point_is_the_von_mises_kernel():
    from scipy.special import i0

    x, y = tp.circular_density(np.array([1.0]), bw=5.0, at=np.array([1.3]))
    want = np.exp(5.0 * np.cos(1.3 - 1.0)) / (2 * np.pi * i0(5.0))
    assert y[0] == pytest.approx(want, rel=1e-12)


def test_circular_density_survives_a_large_concentration():
    """exp(bw*cos) would overflow; the scaled Bessel form must not."""
    x, y = tp.circular_density(RNG.uniform(0, 2 * np.pi, 50), bw=700.0)
    assert np.isfinite(y).all()
    assert y.max() > 0


def test_circular_density_peaks_where_the_mass_is():
    theta = np.concatenate([RNG.normal(1.0, 0.05, 500) % (2 * np.pi)])
    x, y = tp.circular_density(theta, bw=30.0)
    assert x[int(np.argmax(y))] == pytest.approx(1.0, abs=0.1)


def test_cyclic_palette_closes_on_itself():
    first, last = tp.CYCLIC_COLORS[0], tp.CYCLIC_COLORS[-1]
    assert first != last                       # not literally the same string
    a = np.array([int(first[i:i + 2], 16) for i in (1, 3, 5)])
    b = np.array([int(last[i:i + 2], 16) for i in (1, 3, 5)])
    assert np.abs(a - b).max() <= 8            # but visually the same colour


def test_plot_helpers_draw_without_a_display():
    import matplotlib

    matplotlib.use("Agg")
    adata = ad.AnnData(X=RNG.normal(size=(50, 3)),
                       var=pd.DataFrame(index=list("ABC")))
    adata.obsm["X_umap"] = RNG.normal(size=(50, 2))
    adata.obs["tricyclePosition"] = RNG.uniform(0, 2 * np.pi, 50)
    assert tp.plot_embedding_circle_scale(adata) is not None
    assert tp.circle_scale_legend(add_stage_label=True) is not None
    assert tp.plot_ccposition_density(
        adata.obs["tricyclePosition"].to_numpy(),
        groups=np.where(np.arange(50) % 2 == 0, "x", "y")) is not None
