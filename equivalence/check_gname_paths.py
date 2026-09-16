#!/usr/bin/env python3
"""Score the identifier code paths the main equivalence run does not reach.

Reads what ``check_gname_paths.R`` wrote: for each arm, the exact matrix R used
and the positions R produced. Runs ``tricycle-py`` on that same matrix and
compares the angles circularly, against the same bounds ``compare.py`` uses.

The human Ensembl and mouse Ensembl arms are the ones that matter most. R calls
``AnnotationDbi`` at run time there; this package reads a frozen extract of the
same query. If the extract had drifted, these arms would show it.

Usage:
    python check_gname_paths.py <dir_written_by_check_gname_paths.R> [report.md]
"""

from __future__ import annotations

import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import tricyclepy as tp
from compare import TOL_MAX_ARC, TOL_MAX_RAD, TOL_MEDIAN_RAD, circular_difference

ARMS = ("mouse_ensembl", "mouse_symbol", "human_ensembl")


def _strings(dset) -> np.ndarray:
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v)
                       for v in dset[:]], dtype=object)


def _load(path: str):
    with h5py.File(path, "r") as handle:
        shape = tuple(int(v) for v in handle["data/shape"][:])
        matrix = sp.csc_matrix(
            (handle["data/x"][:].astype(np.float64),
             handle["data/i"][:].astype(np.int64),
             handle["data/p"][:].astype(np.int64)), shape=shape)
        return matrix.T.tocsr(), _strings(handle["features"]), \
            _strings(handle["barcodes"])


def main(indir: str, report: str | None = None) -> int:
    rows = []
    failures = []
    for arm in ARMS:
        try:
            meta = dict(line.split("\t", 1) for line in
                        open(f"{indir}/{arm}_meta.tsv", encoding="utf-8")
                        .read().splitlines())
            r_pos = pd.read_csv(f"{indir}/{arm}_position.tsv", sep="\t", dtype=str)
        except FileNotFoundError:
            rows.append((arm, "-", "-", "-", "-", "not run on the R side"))
            continue

        x, genes, barcodes = _load(f"{indir}/{arm}.h5")
        if not np.array_equal(barcodes, r_pos["barcode"].to_numpy()):
            failures.append(f"{arm}: cell order differs")
        embedding = tp.project_cycle_space(
            x, gname=genes, species=meta["species"],
            gname_type=meta["gname_type"])
        py_theta = tp.theta_from_embedding(embedding)
        r_theta = r_pos["theta"].astype(np.float64).to_numpy()

        diff = circular_difference(r_theta, py_theta)
        radius = np.hypot(*(r_pos[c].astype(np.float64).to_numpy()
                            for c in ("pc1", "pc2")))
        arc = diff * radius
        median_d, max_d, max_arc = (float(np.median(diff)), float(diff.max()),
                                    float(arc.max()))

        # The projection gene count is its own check. A mapping that silently
        # lost genes would still give a small angular difference on whatever
        # survived on both sides, so agreement alone is not enough.
        matched = _matched_gene_count(x, genes, meta)
        if matched != int(meta["n_projection_genes"]):
            failures.append(
                f"{arm}: matched {matched} reference genes, R matched "
                f"{meta['n_projection_genes']}")
        for name, value, bound in (("median", median_d, TOL_MEDIAN_RAD),
                                   ("maximum", max_d, TOL_MAX_RAD),
                                   ("arc", max_arc, TOL_MAX_ARC)):
            if value > bound:
                failures.append(f"{arm}: {name} {value:.3e} exceeds {bound:.3e}")
        rows.append((arm, meta["n_cells"],
                     f"{matched} / {meta['n_projection_genes']}",
                     f"{median_d:.3e}", f"{max_d:.3e}", f"{max_arc:.3e}"))

    lines = ["# Identifier code paths, scored against R", "",
             "The main equivalence run exercises `species=\"human\"` with "
             "`gname_type=\"SYMBOL\"`. These arms exercise the others, "
             "including both paths where R queries `AnnotationDbi` at run time "
             "and this package reads a frozen extract instead.", "",
             "| Arm | Cells | Projection genes, Python / R | Median circular "
             "difference | Maximum | Maximum arc error |",
             "|---|---:|---|---|---|---|"]
    for row in rows:
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    lines += ["", f"Bounds are the ones `compare.py` fixes: median "
                  f"{TOL_MEDIAN_RAD:g} rad, maximum {TOL_MAX_RAD:g} rad, arc "
                  f"{TOL_MAX_ARC:g}.", ""]
    verdict = "PASS" if not failures else "FAIL"
    lines.insert(2, f"**Verdict: {verdict}.**" + ("" if not failures else
                 " " + "; ".join(failures)))
    text = "\n".join(lines) + "\n"
    if report:
        with open(report, "w", encoding="utf-8") as handle:
            handle.write(text)
    print(text)
    for failure in failures:
        print(f"failure: {failure}")
    return 0 if not failures else 1


def _matched_gene_count(x, genes, meta) -> int:
    """How many reference genes this package matched, counted independently."""
    from tricyclepy.projection import _match_genes
    from tricyclepy.refdata import get_rotation, map_ensembl_to_symbol

    ref_names, _ = get_rotation(gname_type=meta["gname_type"],
                                species=meta["species"])
    names = genes
    if meta["species"] == "human" and meta["gname_type"] == "ENSEMBL":
        names = map_ensembl_to_symbol(genes, species="human")
    data_idx, _ = _match_genes(names, ref_names)
    return int(data_idx.size)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
