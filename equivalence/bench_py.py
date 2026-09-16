#!/usr/bin/env python3
"""Time one tricycle-py function. Peak memory is measured outside, by
``/usr/bin/time``, so the two sides are measured the same way.

Usage:
    python bench_py.py <input.h5ad> <what> [replicate_factor]
    what: load | position | schwabe | loess_direct
"""

from __future__ import annotations

import sys
import time

import anndata as ad
import numpy as np
import scipy.sparse as sp

import tricyclepy as tp


def main() -> None:
    path, what = sys.argv[1], sys.argv[2]
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    adata = ad.read_h5ad(path)
    if reps > 1:
        x = sp.vstack([adata.X] * reps, format="csr")
        names = np.concatenate([[f"{r}_{b}" for b in adata.obs_names]
                                for r in range(reps)])
        adata = ad.AnnData(X=x, var=adata.var.copy())
        adata.obs_names = names
    print(f"cells\t{adata.n_obs}")

    if what == "load":
        print("seconds\t0")
        return

    if what == "position":
        t0 = time.perf_counter()
        tp.estimate_cycle_position(adata, species="human", gname_type="SYMBOL")
    elif what == "schwabe":
        t0 = time.perf_counter()
        tp.estimate_schwabe_stage(adata, species="human", gname_type="SYMBOL")
    elif what == "loess_direct":
        tp.estimate_cycle_position(adata, species="human", gname_type="SYMBOL")
        y = adata[:, "TOP2A"].X
        y = np.asarray(y.todense() if hasattr(y, "todense") else y).ravel()
        t0 = time.perf_counter()          # exclude the projection from this one
        tp.fit_periodic_loess(
            np.asarray(adata.obs["tricyclePosition"], dtype=float), y)
    else:
        raise SystemExit("unknown benchmark target")
    print(f"seconds\t{time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()
