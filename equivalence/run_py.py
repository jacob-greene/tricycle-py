#!/usr/bin/env python3
"""Run ``tricycle-py`` on the ``.h5ad`` and write its outputs as text.

Output file names and columns mirror ``run_r.R`` so ``compare.py`` can line them
up without guessing.

``--perturb`` deliberately breaks the Python side. It exists so the comparison
can be shown to fail when the two implementations really do differ. A
comparison that has never been made to fail is not evidence.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd

import tricyclepy as tp

G17 = "%.17g"

PERTURBATIONS = {
    "none": "no change",
    "drop-gene": "remove one reference gene from the projection",
    "flip-pc2": "negate the PC2 column of the reference rotation",
    "nudge-rotation": "add --epsilon to one reference gene's PC1 weight",
}


def _fmt(values) -> np.ndarray:
    return np.asarray([G17 % float(v) for v in np.asarray(values, dtype=float)],
                      dtype=object)


def _reference(species: str, gname_type: str, perturb: str, gene: str,
               epsilon: float):
    names, rotation = tp.get_rotation(gname_type=gname_type, species=species)
    if perturb == "none":
        return None
    names = np.asarray(names, dtype=object).copy()
    rotation = np.asarray(rotation, dtype=float).copy()
    if perturb == "flip-pc2":
        rotation[:, 1] *= -1.0
        return names, rotation
    where = np.flatnonzero(names == gene)
    if where.size == 0:
        raise SystemExit(f"perturbation gene {gene!r} is not in the reference")
    if perturb == "drop-gene":
        keep = np.ones(names.shape[0], dtype=bool)
        keep[where[0]] = False
        return names[keep], rotation[keep]
    if perturb == "nudge-rotation":
        rotation[where[0], 0] += epsilon
        return names, rotation
    raise SystemExit(f"unknown perturbation {perturb!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5ad")
    parser.add_argument("outdir")
    parser.add_argument("--species", default="human")
    parser.add_argument("--gname-type", default="SYMBOL")
    parser.add_argument("--batch-column", default=None)
    parser.add_argument("--perturb", default="none", choices=sorted(PERTURBATIONS))
    parser.add_argument("--perturb-gene", default="TOP2A")
    parser.add_argument("--epsilon", type=float, default=1e-6)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    ref = _reference(args.species, args.gname_type, args.perturb,
                     args.perturb_gene, args.epsilon)
    if args.perturb != "none":
        print(f"PERTURBED: {args.perturb} ({PERTURBATIONS[args.perturb]})")

    t0 = time.perf_counter()
    tp.estimate_cycle_position(adata, species=args.species,
                               gname_type=args.gname_type, ref=ref)
    t_position = time.perf_counter() - t0

    embedding = np.asarray(adata.obsm["tricycleEmbedding"])
    theta = np.asarray(adata.obs["tricyclePosition"], dtype=float)
    info = adata.uns["tricycleEmbedding"]
    genes = np.asarray(info["genes"], dtype=object)
    rotation = np.asarray(info["rotation"], dtype=float)

    out = args.outdir
    pd.DataFrame({"barcode": adata.obs_names, "theta": _fmt(theta),
                  "pc1": _fmt(embedding[:, 0]), "pc2": _fmt(embedding[:, 1])}
                 ).to_csv(f"{out}/py_position.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": genes, "pc1_rot": _fmt(rotation[:, 0]),
                  "pc2_rot": _fmt(rotation[:, 1])}
                 ).to_csv(f"{out}/py_projection_genes.tsv", sep="\t", index=False)

    columns = adata.var_names.get_indexer(genes)
    block = adata.X[:, columns]
    block = np.asarray(block.todense() if hasattr(block, "todense") else block,
                       dtype=float)
    pd.DataFrame({"gene": genes, "mean": _fmt(block.mean(axis=0))}
                 ).to_csv(f"{out}/py_gene_means.tsv", sep="\t", index=False)

    t0 = time.perf_counter()
    try:
        batch = (np.asarray(adata.obs[args.batch_column])
                 if args.batch_column else None)
        tp.estimate_schwabe_stage(adata, species=args.species,
                                  gname_type=args.gname_type, batch=batch)
        stages = np.asarray(adata.obs["CCStage"].astype(object))
        stages = np.asarray(["NA" if s is None or (isinstance(s, float) and np.isnan(s))
                             else s for s in stages], dtype=object)
        pd.DataFrame({"barcode": adata.obs_names, "stage": stages}
                     ).to_csv(f"{out}/py_stage.tsv", sep="\t", index=False)
    except Exception as exc:                                   # noqa: BLE001
        print(f"Schwabe failed: {exc}")
    t_stage = time.perf_counter() - t0

    probe = "TOP2A" if "TOP2A" in adata.var_names else adata.var_names[0]
    y = adata[:, probe].X
    y = np.asarray(y.todense() if hasattr(y, "todense") else y, dtype=float).ravel()
    t0 = time.perf_counter()
    fit = tp.fit_periodic_loess(theta, y)
    t_loess = time.perf_counter() - t0
    pd.DataFrame({"barcode": adata.obs_names, "y": _fmt(y),
                  "fitted": _fmt(fit.fitted)}
                 ).to_csv(f"{out}/py_loess.tsv", sep="\t", index=False)
    pd.DataFrame({"x": _fmt(fit.pred_x), "y": _fmt(fit.pred_y)}
                 ).to_csv(f"{out}/py_loess_pred.tsv", sep="\t", index=False)

    try:
        pca = tp.run_pca_cc_genes(adata, species=args.species,
                                  gname_type=args.gname_type, ntop=500,
                                  n_components=10)
        pd.DataFrame({"gene": pca.genes, "pc1": _fmt(pca.rotation[:, 0]),
                      "pc2": _fmt(pca.rotation[:, 1])}
                     ).to_csv(f"{out}/py_pca_rotation.tsv", sep="\t", index=False)
        pd.DataFrame({"percent_var": _fmt(pca.percent_var)}
                     ).to_csv(f"{out}/py_pca_percentvar.tsv", sep="\t", index=False)
    except Exception as exc:                                   # noqa: BLE001
        print(f"run_pca_cc_genes failed: {exc}")

    dx, dy = tp.circular_density(theta, bw=30.0)
    pd.DataFrame({"x": _fmt(dx), "y": _fmt(dy)}
                 ).to_csv(f"{out}/py_density.tsv", sep="\t", index=False)

    summary = {
        "probe_gene": probe,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_projection_genes": len(genes),
        "loess_rsquared_direct": G17 % fit.rsquared,
        "seconds_estimate_cycle_position": G17 % t_position,
        "seconds_estimate_schwabe_stage": G17 % t_stage,
        "seconds_fit_periodic_loess": G17 % t_loess,
        "tricyclepy_version": tp.__version__,
        "python_version": sys.version.split()[0],
        "perturbation": args.perturb,
    }
    if "total_counts" in adata.obs:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            diagnosis = tp.diagnose_total_umi(
                theta, np.asarray(adata.obs["total_counts"], dtype=float))
        summary["diagnose_peak"] = G17 % diagnosis.peak
        summary["diagnose_valley"] = G17 % diagnosis.valley
        summary["diagnose_difference"] = G17 % diagnosis.difference

    with open(f"{out}/py_summary.tsv", "w", encoding="utf-8") as handle:
        for key, value in summary.items():
            handle.write(f"{key}\t{value}\n")
    print(f"Python side done: {out}")


if __name__ == "__main__":
    main()
