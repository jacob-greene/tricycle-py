#!/usr/bin/env python3
"""Turn the HDF5 written by ``export_seurat.R`` into an ``.h5ad``.

Every fingerprint the R side recorded is re-checked here. A conversion that
drops, reorders or rescales anything fails loudly rather than producing an
``.h5ad`` that merely looks plausible.

Usage:
    python build_h5ad.py <input.h5> <output.h5ad>
"""

from __future__ import annotations

import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

#: Fingerprints are sums of float64 values, so they carry accumulated rounding.
#: This bound is far tighter than any difference a real conversion bug makes.
FINGERPRINT_TOL = 1e-9


def _strings(dset) -> np.ndarray:
    raw = dset[:]
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in raw],
                      dtype=object)


def build(in_path: str, out_path: str) -> ad.AnnData:
    with h5py.File(in_path, "r") as handle:
        shape = tuple(int(v) for v in handle["data/shape"][:])
        matrix = sp.csc_matrix(
            (handle["data/x"][:].astype(np.float64),
             handle["data/i"][:].astype(np.int64),
             handle["data/p"][:].astype(np.int64)),
            shape=shape,
        )                                             # genes x cells, as in R
        features = _strings(handle["features"])
        barcodes = _strings(handle["barcodes"])
        obs = {k: _strings(handle[f"obs/{k}"]) for k in handle["obs"].keys()}
        total_counts = (handle["total_counts"][:].astype(np.float64)
                        if "total_counts" in handle else None)
        gene_sums = handle["fingerprint/gene_sums"][:].astype(np.float64)
        cell_sums = handle["fingerprint/cell_sums"][:].astype(np.float64)
        nnz = int(handle["fingerprint/nnz"][0])

    if matrix.nnz != nnz:
        raise AssertionError(f"nnz mismatch: R {nnz}, Python {matrix.nnz}")
    got_gene = np.asarray(matrix.sum(axis=1)).ravel()
    got_cell = np.asarray(matrix.sum(axis=0)).ravel()
    for name, want, got in (("gene sums", gene_sums, got_gene),
                            ("cell sums", cell_sums, got_cell)):
        worst = float(np.abs(want - got).max())
        if worst > FINGERPRINT_TOL:
            raise AssertionError(f"{name} differ from R by up to {worst:.3e}")
        print(f"fingerprint {name:<11} max |R - Python| = {worst:.3e}")

    frame = pd.DataFrame(obs, index=pd.Index(barcodes, name=None))
    for column in frame.columns:
        converted = pd.to_numeric(frame[column], errors="coerce")
        if converted.notna().all():
            frame[column] = converted
    if total_counts is not None:
        frame["total_counts"] = total_counts

    adata = ad.AnnData(
        X=matrix.T.tocsr(),                           # cells x genes
        obs=frame,
        var=pd.DataFrame(index=pd.Index(features, name=None)),
    )
    adata.write_h5ad(out_path, compression="gzip")
    print(f"wrote {out_path}  cells={adata.n_obs} genes={adata.n_vars} "
          f"nnz={adata.X.nnz}")
    return adata


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    build(sys.argv[1], sys.argv[2])
