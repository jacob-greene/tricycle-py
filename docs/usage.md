# Using tricycle-py

The [README](../README.md) has the quick start. This page covers the details.

## Basic use

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

## Discrete stages

```python
tp.estimate_schwabe_stage(adata, species="human", gname_type="SYMBOL",
                          batch_key="sample")
adata.obs["CCStage"]      # G1.S, S, G2, G2.M, M.G1, or unassigned
```

## Diagnostics and plots

```python
d = tp.diagnose_total_umi(adata.obs["tricyclePosition"], adata.obs["total_counts"])
d.passed, d.difference

ax = tp.plot_embedding_circle_scale(adata, embedding_key="X_umap")
tp.circle_scale_legend(add_stage_label=True)
```

An angle wraps, so it needs a cyclic colormap. `tp.cyclic_colormap()` is the
one tricycle uses.

## Learning a new reference

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

