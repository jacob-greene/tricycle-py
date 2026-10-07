"""Example data shipped with the package.

``neurosphere_example`` is the dataset R ``tricycle`` ships under the same
name and uses in its vignette: 400 mouse neurosphere cells, 1,500 genes,
log-normalised. It is carried over from ``tricycle`` 1.12.0 together with the
``tricyclePosition`` R computes on it, so this port can be checked against R
without installing R. The export is ``equivalence/extract_example_data.R``.
"""

from __future__ import annotations

from importlib import resources

import numpy as np

__all__ = ["neurosphere_example"]


def neurosphere_example():
    """Load tricycle's 400-cell mouse neurosphere example as ``AnnData``.

    Returns
    -------
    anndata.AnnData
        ``X`` holds log-normalised expression, cells by genes, as a sparse
        float64 matrix. ``var_names`` are mouse Ensembl ids and
        ``var["symbol"]`` the gene symbols. ``obs["tricyclePosition_R"]`` is
        the position R ``tricycle`` 1.12.0 computes with
        ``estimate_cycle_position(species = "mouse", gname.type = "ENSEMBL")``.

    Examples
    --------
    >>> import tricyclepy as tp
    >>> adata = tp.datasets.neurosphere_example()
    >>> tp.estimate_cycle_position(adata)          # mouse, Ensembl ids: the defaults
    >>> adata.obs["tricyclePosition"].head()
    """
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    path = resources.files(__package__).joinpath("data", "neurosphere_example.npz")
    with resources.as_file(path) as p, np.load(p, allow_pickle=False) as z:
        x = sp.csr_matrix((z["data"], z["indices"], z["indptr"]),
                          shape=tuple(z["shape"]))
        var = pd.DataFrame({"symbol": z["symbol"].astype(object)},
                           index=pd.Index(z["ensembl"].astype(object)))
        obs = pd.DataFrame({"tricyclePosition_R": z["tricyclePosition_R"]},
                           index=pd.Index(z["cells"].astype(object)))
    return ad.AnnData(X=x, obs=obs, var=var)
