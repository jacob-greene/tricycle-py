"""Download the public data and build the README figure from it.

This is the whole path from public files to ``docs/readme_figure.png``:

1. Download two ``.h5ad`` files (see ``SOURCES``) into a cache directory and
   check each against its md5 checksum. A cached file whose checksum matches
   is not downloaded again.
2. Keep only what the figure needs from each file: the log-normalised
   expression ``X``, the UMAP, and the cell type labels with their stored
   colours.
3. If R and the Bioconductor package ``tricycle`` are available, run
   ``docs/readme_figure_r.R`` on the same matrix to get R's position for
   column 5. Without R, the figure has four columns.
4. Call ``docs/readme_figure.py``, which computes the Python position with
   this package and draws the figure.

Usage::

    uv pip install "tricycle-py[figure] @ git+https://github.com/jacob-greene/tricycle-py"
    python docs/make_readme_figure.py --out docs/readme_figure.png

Options: ``--cache DIR`` (default ``~/.cache/tricycle-py-figure``), ``--no-r``
to skip R, ``--rscript PATH`` to name the ``Rscript`` to use. The download is
about 3.2 GB.

Data sources
------------

* CD34+ hematopoietic stem and progenitor cells: ``cd34_multiome_rna.h5ad``
  from Zenodo record 6383269 (Persad et al. 2022, doi:10.5281/zenodo.6383269,
  CC BY 4.0). The md5 is the one Zenodo publishes.
* T cell depleted bone marrow: ``preprocessed_t-cell-depleted-bm-rna.h5ad``,
  the file the Mellon tutorial downloads (Otto et al. 2024, *Nature Methods*;
  github.com/settylab/Mellon, ``notebooks/basic_tutorial.ipynb``). It holds
  the same expression as ``bm_multiome_rna.h5ad`` in the Zenodo record, plus
  1,188 B lineage cells, a UMAP and finer cell type labels. No checksum is
  published for it; the md5 here was measured on 2026-10-07.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import readme_figure  # noqa: E402

ZENODO = "https://zenodo.org/records/6383269/files/{}?download=1"

#: One row per dataset, in figure order.
SOURCES = [
    dict(label="Human bone marrow",
         file="preprocessed_t-cell-depleted-bm-rna.h5ad",
         url="https://fh-pi-setty-m-eco-public.s3.amazonaws.com/mellon-tutorial/"
             "preprocessed_t-cell-depleted-bm-rna.h5ad",
         md5="64909a0be6cd51fe79b7d70b5f19bf1c",
         annotate="HSC,Mono,Ery,NaiveB"),
    dict(label="Human CD34+ HSPC",
         file="cd34_multiome_rna.h5ad",
         url=ZENODO.format("cd34_multiome_rna.h5ad"),
         md5="4cd8d82adfe267f54e13d8a383918fd0",
         annotate="HSC,Mono,Ery,CLP"),
]

R_KEY = "tricyclePosition_R"


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url: str, dest: Path, md5: str) -> Path:
    """Download ``url`` to ``dest`` unless a file with the right md5 is there."""
    if dest.exists() and md5sum(dest) == md5:
        print(f"cached: {dest.name}")
        return dest
    part = dest.with_suffix(dest.suffix + ".part")
    print(f"downloading {dest.name} ...", flush=True)
    with urllib.request.urlopen(url) as r, open(part, "wb") as f:
        shutil.copyfileobj(r, f, length=1 << 22)
    got = md5sum(part)
    if got != md5:
        raise RuntimeError(f"{dest.name}: md5 {got}, expected {md5}")
    part.replace(dest)
    return dest


def slim(src: Path, dest: Path, key: str = "celltype") -> ad.AnnData:
    """Keep X, the UMAP and the cell type labels with their colours."""
    a = ad.read_h5ad(src)
    out = ad.AnnData(X=sp.csr_matrix(a.X), obs=pd.DataFrame(index=a.obs_names),
                     var=pd.DataFrame(index=a.var_names))
    out.obsm["X_umap"] = np.asarray(a.obsm["X_umap"], dtype=np.float64)
    out.obs[key] = a.obs[key].astype("category")
    out.uns[f"{key}_colors"] = list(a.uns[f"{key}_colors"])
    out.write_h5ad(dest)
    return out


def have_r(rscript: str) -> bool:
    if shutil.which(rscript) is None:
        return False
    ok = subprocess.run([rscript, "-e", "suppressMessages(library(tricycle))"],
                        capture_output=True)
    return ok.returncode == 0


def r_position(a: ad.AnnData, workdir: Path, rscript: str) -> np.ndarray:
    """Run R tricycle on ``a.X`` and return its position, in cell order."""
    workdir.mkdir(parents=True, exist_ok=True)
    # genes by cells, as R expects; 17 significant digits round-trip a double
    scipy.io.mmwrite(workdir / "logcounts.mtx", sp.csr_matrix(a.X, dtype=np.float64).T,
                     precision=17)
    pd.Series(a.var_names).to_csv(workdir / "genes.txt", index=False, header=False)
    pd.Series(a.obs_names).to_csv(workdir / "cells.txt", index=False, header=False)
    subprocess.run([rscript, str(HERE / "readme_figure_r.R"), str(workdir)], check=True)
    r = pd.read_csv(workdir / "r_position.tsv", sep="\t", index_col=0)
    if not r.index.equals(pd.Index(a.obs_names)):
        raise RuntimeError("R output is not in input cell order")
    return r["theta"].to_numpy(np.float64)


def main(argv=None) -> dict:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", default=str(HERE / "readme_figure.png"))
    p.add_argument("--cache", default=str(Path.home() / ".cache" / "tricycle-py-figure"))
    p.add_argument("--no-r", action="store_true", help="leave out column 5")
    p.add_argument("--rscript", default="Rscript")
    args = p.parse_args(argv)

    cache = Path(args.cache).expanduser()
    cache.mkdir(parents=True, exist_ok=True)
    use_r = not args.no_r and have_r(args.rscript)
    if not args.no_r and not use_r:
        print("R with the tricycle package was not found: column 5 is left out")

    fig_argv = []
    for s in SOURCES:
        raw = fetch(s["url"], cache / s["file"], s["md5"])
        stem = raw.name.removesuffix(".h5ad")
        dest = cache / f"{stem}.fig_input.h5ad"
        a = slim(raw, dest)
        if use_r:
            a.obs[R_KEY] = r_position(a, cache / f"{stem}.r", args.rscript)
            a.write_h5ad(dest)
        fig_argv += ["--row", s["label"], str(dest), "--annotate", s["annotate"]]
    return readme_figure.main(fig_argv + ["--r-key", R_KEY, "--out", args.out])


if __name__ == "__main__":
    main()
