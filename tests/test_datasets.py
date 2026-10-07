"""The shipped example loads, and reproduces the position R computes on it."""

import numpy as np
import scipy.sparse as sp

import tricyclepy as tp


def _circ_diff(a, b):
    d = np.mod(np.asarray(a) - np.asarray(b), 2.0 * np.pi)
    return np.minimum(d, 2.0 * np.pi - d)


def test_neurosphere_example_shape_and_content():
    adata = tp.datasets.neurosphere_example()
    assert adata.shape == (400, 1500)
    assert sp.issparse(adata.X) and adata.X.dtype == np.float64
    assert adata.var_names.str.startswith("ENSMUSG").all()
    assert "symbol" in adata.var
    theta_r = adata.obs["tricyclePosition_R"].to_numpy()
    assert np.all((theta_r >= 0) & (theta_r < 2 * np.pi))


def test_neurosphere_example_matches_r_with_defaults():
    adata = tp.datasets.neurosphere_example()
    tp.estimate_cycle_position(adata)                 # mouse + ENSEMBL defaults
    assert adata.uns["tricycleEmbedding"]["n_genes"] == 500
    diff = _circ_diff(adata.obs["tricyclePosition"], adata.obs["tricyclePosition_R"])
    assert diff.max() < 1e-9


def test_neurosphere_example_returns_independent_copies():
    a = tp.datasets.neurosphere_example()
    a.X.data[:] = 0.0
    b = tp.datasets.neurosphere_example()
    assert b.X.data.any()
