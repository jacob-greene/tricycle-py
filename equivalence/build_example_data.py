"""Pack the export of ``extract_example_data.R`` into the shipped ``.npz``.

Usage::

    python build_example_data.py <export dir> src/tricyclepy/data/neurosphere_example.npz

The values are parsed from ``%.17g`` text, so they round-trip exactly. The
matrix is stored cells by genes, the ``AnnData`` orientation.
"""

import gzip
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

src, out = sys.argv[1], sys.argv[2]
with gzip.open(f"{src}/logcounts.tsv.gz", "rt") as fh:
    n_genes, n_cells = map(int, fh.readline().split())
    trip = np.loadtxt(fh, dtype=np.float64, ndmin=2)
gi = trip[:, 0].astype(np.int64)
ci = trip[:, 1].astype(np.int64)
x = sp.csr_matrix((trip[:, 2], (ci, gi)), shape=(n_cells, n_genes))
x.sort_indices()
genes = pd.read_csv(f"{src}/genes.tsv", sep="\t", keep_default_na=False)
cells = pd.read_csv(f"{src}/cells.tsv", sep="\t")
assert len(genes) == n_genes and len(cells) == n_cells
np.savez_compressed(
    out,
    data=x.data, indices=x.indices.astype(np.int32), indptr=x.indptr.astype(np.int32),
    shape=np.array(x.shape, dtype=np.int64),
    ensembl=genes["ensembl"].to_numpy(str), symbol=genes["symbol"].to_numpy(str),
    cells=cells["cell"].to_numpy(str),
    tricyclePosition_R=cells["tricyclePosition_R"].to_numpy(np.float64),
)
print(f"wrote {out}: {n_cells} cells x {n_genes} genes, nnz={x.nnz}")
