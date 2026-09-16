# Equivalence harness

This directory holds the evidence that `tricycle-py` reproduces the R package
`tricycle`, and the scripts that produce it. The measured outcome is in
[`RESULTS.md`](RESULTS.md).

Nothing here reimplements tricycle on the R side. It installs and runs the real
Bioconductor package, so the comparison is against the real thing.

## What it does

1. `extract_reference_data.R` dumps the data objects `tricycle` ships, and the
   `AnnotationDbi` queries its `ENSEMBL` paths depend on, into plain text.
   Those files ship inside the Python package. They are never regenerated.
   Doubles are written with `%.17g`, not `write.table`'s default 15 significant
   digits. The default perturbed 933 of the 1000 rotation weights by up to 33
   units in the last place, and that was the dominant term in the residual
   disagreement with R -- larger than summation order.
2. `export_seurat.R` exports one assay of a Seurat `.rds` to a plain HDF5
   container, and writes fingerprints of it.
3. `build_h5ad.py` turns that container into an `.h5ad` and re-checks every
   fingerprint. A conversion that drops, reorders or rescales anything fails
   here rather than surviving as a plausible-looking `.h5ad`.
4. `run_r.R` runs R `tricycle` on the `.rds`.
5. `run_py.py` runs `tricycle-py` on the `.h5ad`.
6. `compare.py` scores the two against a tolerance fixed before measuring, and
   writes a report.
7. `benchmark.sh`, `bench_r.R` and `bench_py.py` measure runtime and peak
   resident memory on both sides, each in its own process.
8. `check_gname_paths.R` and `check_gname_paths.py` score the identifier code
   paths the main run does not reach.

## Why an angle needs a circular comparison

`tricyclePosition` wraps. Two values a hair apart across zero differ by almost
`2*pi` under `abs(a - b)`. `compare.py` wraps the difference onto `[0, pi]`, so
0.001 and 6.282 read as 0.002 apart, which is what they are.

## Running it

The inputs are a Seurat `.rds` and the `.h5ad` built from it. They must hold the
same cells, the same genes and the same log-normalised values, which is what
steps 2 and 3 guarantee.

```bash
RDS=/path/to/your_rna.rds
OUT=/path/to/scratch

Rscript equivalence/export_seurat.R "$RDS" "$OUT/data.h5"
python  equivalence/build_h5ad.py    "$OUT/data.h5" "$OUT/data.h5ad"

Rscript equivalence/run_r.R  "$RDS"        "$OUT/r"  human SYMBOL sample
python  equivalence/run_py.py "$OUT/data.h5ad" "$OUT/py" \
        --species human --gname-type SYMBOL --batch-column sample

python equivalence/compare.py "$OUT/r" "$OUT/py" "$OUT/RESULTS.md"
```

`compare.py` exits 0 on pass and 1 on fail.

Pass the same `species`, `gname_type` and batch column to both sides. Giving
one side a different batch column changes only the Schwabe call, which then
looks like a port defect and is not one.

### A caution about environment modules

If R comes from an environment module, it may export a `PYTHONPATH`. That
shadows a virtual environment's packages with older ones and breaks the Python
side at import time. Run the Python commands with `env -u PYTHONPATH`, which is
what `benchmark.sh` does.

## Showing that the comparison can fail

A comparison that has never fired is not evidence. `run_py.py --perturb`
deliberately breaks the Python side, through the real code path rather than by
editing the comparison arithmetic.

| `--perturb` | What it changes |
|---|---|
| `drop-gene` | removes one gene from the reference |
| `flip-pc2` | negates the PC2 column of the reference rotation |
| `nudge-rotation` | adds `--epsilon` to one gene's PC1 weight |

```bash
python equivalence/run_py.py "$OUT/data.h5ad" "$OUT/py_bad" --perturb flip-pc2
python equivalence/compare.py "$OUT/r" "$OUT/py_bad" "$OUT/BAD.md"   # must exit 1
```

Measured sensitivity, nudging one weight out of 461 on a real dataset:

| Nudge to one rotation weight | Median circular difference | Verdict |
|---|---|---|
| none | 0.000e+00 rad | pass |
| 1e-10 | 1.609e-11 rad | pass |
| 1e-9 | 1.609e-10 rad | **fail** |
| 1e-8 | 1.609e-09 rad | **fail** |
| 1e-7 | 1.609e-08 rad | **fail** |

The comparison catches a change of 1e-9 to a single weight. The clean run sits
eleven orders of magnitude below that floor.

## The identifier code paths

`check_gname_paths.R` and `check_gname_paths.py` score the three
`species` / `gname_type` combinations the main run does not reach. The mouse
arms use `neurosphere_example`, the SCE `tricycle` ships, so they need no data
of yours. The human Ensembl arm relabels a human matrix to the Ensembl ids of
its symbols, keeping only symbols with exactly one Ensembl id so the
relabelling is a bijection.

```bash
Rscript equivalence/check_gname_paths.R "$OUT/gname" "$RDS"   # $RDS optional
python  equivalence/check_gname_paths.py "$OUT/gname" "$OUT/gname.md"
```

Both the angle and the number of reference genes matched are scored. A mapping
that silently dropped genes would still agree on whatever survived, so the
angle alone is not enough. Measured outcome:
[`RESULTS_gname_paths.md`](RESULTS_gname_paths.md).

## Anonymising the report

`compare.py --anonymise-barcodes` names cells by row index instead of barcode.
A barcode often carries a dataset or sample identifier, and the report is meant
to travel. The committed `RESULTS.md` was produced with this flag.

## Publishing a report

The committed `RESULTS.md` records the human bone marrow single-cell RNA-seq
run described in it. The dataset itself is not redistributed here. Re-run the
harness on your own data to produce your own report.
