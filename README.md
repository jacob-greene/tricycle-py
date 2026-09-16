# tricycle-py

Cell cycle position from single-cell RNA-seq, on `AnnData`.

`tricycle-py` is a Python port of the R/Bioconductor package
[`tricycle`](https://github.com/hansenlab/tricycle) by Zheng, Hansen and
colleagues. It infers a continuous cell cycle position, an angle called
`tricyclePosition`, by projecting log-normalised expression onto a fixed
reference cell cycle embedding learned from mouse neurosphere data.

The port targets `tricycle` 1.12.0 and is written function by function from the
R sources. The reference embedding, the marker gene lists and the gene
identifier maps are carried over from the R package rather than regenerated.

**Licence: GPL-3.** `tricycle` is GPL-3, so this port is too. See
[`ATTRIBUTION.md`](ATTRIBUTION.md) for the copyright notice, what the licence
obliges, and credit to the original authors.

**Input dtype.** Expression is accepted in any dtype, and the per-gene means
are accumulated in float64 whichever you give. A float32 matrix, which is the
common case in the scverse ecosystem, reaches the same answer as float64: both
agree with R at 8.882e-16 radians on the dataset above.

## Install

```bash
pip install tricycle-py          # or: uv pip install tricycle-py
```

Plotting helpers need `matplotlib`:

```bash
pip install "tricycle-py[plot]"
```

## Use

```python
import anndata as ad
import tricyclepy as tp

adata = ad.read_h5ad("my_cells.h5ad")     # X must be log-normalised

tp.estimate_cycle_position(adata, species="human", gname_type="SYMBOL")

adata.obs["tricyclePosition"]             # angle in radians, [0, 2*pi)
adata.obsm["tricycleEmbedding"]           # the two projected components
adata.uns["tricycleEmbedding"]["n_genes"] # reference genes actually matched
```

Two things decide whether the answer is right, and both are silent when wrong.

**The matrix must be log-transformed and library-size normalised.** The
projection is a weighted sum of gene-mean-centred log-expression. Raw counts
give numbers, and the numbers are wrong.

**`species` and `gname_type` must match your `var` index.** The internal
reference is keyed on mouse Ensembl ids, mouse symbols, or upper-cased symbols
for human. An `.h5ad` often carries Ensembl ids in `var_names` with symbols in
`var["gene_name"]`. Check which one you have. Use `gname=` to point at a
different column:

```python
tp.estimate_cycle_position(adata, species="human", gname_type="SYMBOL",
                           gname=adata.var["gene_name"].to_numpy())
```

### Discrete stages

```python
tp.estimate_schwabe_stage(adata, species="human", gname_type="SYMBOL",
                          batch_key="sample")
adata.obs["CCStage"]      # G1.S, S, G2, G2.M, M.G1, or unassigned
```

### Diagnostics and plots

```python
d = tp.diagnose_total_umi(adata.obs["tricyclePosition"], adata.obs["total_counts"])
d.passed, d.difference

ax = tp.plot_embedding_circle_scale(adata, embedding_key="X_umap")
tp.circle_scale_legend(add_stage_label=True)
```

An angle wraps, so it needs a cyclic colormap. `tp.cyclic_colormap()` is the
one tricycle uses.

### Learning a new reference

```python
pca = tp.run_pca_cc_genes(adata, species="human", gname_type="SYMBOL")
ref = pca.as_reference(n_components=2)
tp.estimate_cycle_position(adata, ref=ref, dimred="tricycleEmbedding2")
```

Component signs from a decomposition are arbitrary. Two references learned
separately can differ by a reflection, a rotation, or both, and the positions
they give will differ accordingly.

## API

| Python | R |
|---|---|
| `project_cycle_space` | `project_cycle_space` |
| `estimate_cycle_position` | `estimate_cycle_position` |
| `theta_from_embedding` | `.getTheta` |
| `estimate_schwabe_stage` | `estimate_Schwabe_stage` |
| `fit_periodic_loess` | `fit_periodic_loess` |
| `loess_fit` | `stats::loess`, direct surface |
| `diagnose_total_umi` | `diagnose_totalumi` |
| `run_pca_cc_genes` | `run_pca_cc_genes` |
| `plot_embedding_circle_scale` | `plot_emb_circle_scale` |
| `circle_scale_legend` | `circle_scale_legend` |
| `plot_ccposition_density` | `plot_ccposition_den` |
| `circular_density` | `circular::density.circular` |
| `get_rotation`, `load_neuro_ref` | `.getRotation`, `data(neuroRef)` |
| `load_revelio_gene_list` | `data(RevelioGeneList)` |

**Orientation.** R holds expression as genes by cells. `AnnData` holds it as
cells by genes. Every function here takes the `AnnData` orientation.

**Argument defaults match R**, including `gname_type="ENSEMBL"` and
`species="mouse"`. A script ported from R keeps its behaviour.

## Differences from the R package, stated plainly

| Area | Difference |
|---|---|
| `loess` surface | R's `stats::loess` defaults to `surface="interpolate"`, a k-d tree approximation. `loess_fit` implements the exact `surface="direct"` computation. It matches R's direct surface to floating-point noise and differs from R's default by the interpolation error. Both gaps are measured in `equivalence/`. |
| Gene identifier maps | R queries `org.Hs.eg.db` / `org.Mm.eg.db` at run time. This package ships a frozen extract of the same query. A newer annotation release will change R's answer and not this package's. `PROVENANCE.tsv` records the frozen version; `ensembl_to_symbol=` overrides it. |
| Version suffixes | Neither R nor this package strips an Ensembl version suffix. `ENSG00000141510.16` fails to map on both sides. Strip it yourself. |
| `diagnose_total_umi` | It calls `fit_periodic_loess`, so it inherits the surface difference above. Against R's *default* it differs by 1.444e-02, about 1.7 percent, which is R's k-d tree approximation error rather than a port difference. Against R's exact surface it agrees at 1.776e-15. |
| Plot output | Plots are `matplotlib` axes, not `ggplot2` objects. The computations underneath, including the circular density, are ported; the rendering is not pixel-identical. |
| PCA backend | `run_pca_cc_genes` uses a NumPy SVD instead of `scater::runPCA`. Gene selection, its variance ranking, centring and scaling all match; the measured rotation difference is 8.0e-16 after matching component signs. Signs themselves are arbitrary on both sides. |

## Equivalence with the R package

The claim that this port reproduces R `tricycle` is tested, not asserted. The
harness in [`equivalence/`](equivalence/) runs the real Bioconductor package and
this package on the same cells, the same genes and the same normalised matrix,
then compares the angles circularly.

Measured on human bone marrow single-cell RNA-seq, 8,627 cells, 17,226 genes,
461 of the 500 reference genes matched. R `tricycle` 1.12.0 against
`tricycle-py` 0.1.0.

| Statistic | Measured | Bound, fixed before measuring |
|---|---|---|
| Median circular difference in `tricyclePosition` | 0.000e+00 rad | 1e-9 rad |
| Maximum circular difference, over every cell | 8.882e-16 rad | 1e-6 rad |
| Maximum arc error, difference times projected radius | 1.148e-14 | 1e-9 |
| Cells with an identical `CCStage` call | 8,627 of 8,627 | all |
| `fit_periodic_loess` against R's direct surface | 4.219e-14 | 1e-8 |
| Reference rotation weights against R's | 0.000e+00, bit-identical | - |
| `circular_density` against `circular::density.circular` | 2.220e-16 | - |
| `diagnose_total_umi` against R's direct surface | 1.776e-15 | 1e-8 |
| `run_pca_cc_genes` rotation, after matching sign | 8.049e-16 | 1e-6 |

8.9e-16 radians is four units in the last place of a double, and the median is
exactly zero: most cells agree bit for bit. No cell inverted. What residual
there is concentrates on cells near the embedding centre, where `atan2` is ill
conditioned, which is why the arc bound is scored as well.

### All four identifier paths, not just the one

The run above uses `species="human"` with `gname_type="SYMBOL"`. The other
three combinations resolve reference gene names differently, and two of them
are where R queries `org.*.eg.db` at run time while this package reads a frozen
extract. Each was scored against R separately, on real data.

| Arm | Cells | Reference genes matched, Python / R | Median circular difference | Maximum |
|---|---:|---|---|---|
| mouse, Ensembl | 400 | 500 / 500 | 0.000e+00 rad | 1.776e-15 rad |
| mouse, symbol | 400 | 500 / 500 | 0.000e+00 rad | 1.776e-15 rad |
| human, Ensembl | 8,627 | 437 / 437 | 0.000e+00 rad | 8.882e-16 rad |

The gene count is checked as well as the angle. A mapping that silently lost
genes would still agree closely on whatever survived on both sides, so
agreement alone would not have caught it. Full report:
[`equivalence/RESULTS_gname_paths.md`](equivalence/RESULTS_gname_paths.md).

### Showing the comparison can fail

The comparison has been shown to fail. Nudging a single reference weight out of
461 by 1e-9 is caught; 1e-10 passes. The clean run sits eleven orders of
magnitude below that floor. See
[`equivalence/RESULTS.md`](equivalence/RESULTS.md) for the full report and
[`equivalence/README.md`](equivalence/README.md) for how to reproduce it.

## Speed and memory

Same cells, same machine, one process per measurement, peak resident set size
from `/usr/bin/time -v`. Peak memory includes loading the object, which costs
4,042 MB in R and 348 MB in Python before any function runs. Raw numbers,
including that baseline row, are in
[`equivalence/results/`](equivalence/results/).

| Cells | Function | R seconds | Python seconds | R peak MB | Python peak MB |
|---:|---|---:|---:|---:|---:|
| 8,627 | `estimate_cycle_position` | 2.41 | 0.47 | 4,178 | 583 |
| 8,627 | `estimate_schwabe_stage` | 2.44 | 0.64 | 4,484 | 474 |
| 8,627 | `fit_periodic_loess`, exact surface | 8.29 | 5.97 | 4,186 | 1,593 |
| 103,524 | `estimate_cycle_position` | 12.54 | 7.74 | 11,605 | 3,869 |
| 103,524 | `estimate_schwabe_stage` | 18.67 | 5.45 | 11,812 | 4,370 |

The projection never densifies the whole matrix. It takes the per-gene means in
one sparse pass, then centres and multiplies in blocks over only the matched
genes. Lower `chunk_size` to cut peak memory further. Every cell's result
depends only on its own row, so the block size does not change what is
computed, but it can change how BLAS associates the sums within a row:
measured at 1.7e-14 radians between blocks of 20,000 and 1,000.

R's own `fit_periodic_loess` takes 0.76 s at 8,627 cells using its k-d tree
approximation. The 8.29 s in the table is R computing the exact surface, which
is the like-for-like comparison.

## Tests

```bash
pip install "tricycle-py[test]"
pytest
```

The test suite runs without any external dataset. It pins the shipped reference
data, the gene matching rules including the duplicate and unmapped cases, the
angle convention, the loess fit against closed-form cases, and the Schwabe
assignment rules.

## Citation

Cite the original work, not this port:

> Zheng SC, Stein-O'Brien G, Augustin JJ, Slosberg J, Carosso GA, Winer B,
> Shin G, Bjornsson HT, Goff LA, Hansen KD. Universal prediction of cell cycle
> position using transfer learning. *Genome Biology* 23:41 (2022).
> doi:10.1186/s13059-021-02581-y

Add Schwabe et al. (2020) if you use `estimate_schwabe_stage`. See
[`CITATION.cff`](CITATION.cff).
