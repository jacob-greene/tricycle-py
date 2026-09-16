#!/usr/bin/env Rscript
# Export one assay of a Seurat object to a plain HDF5 container.
#
# The equivalence harness needs both ecosystems to read the same numbers. The R
# side reads the .rds directly. The Python side reads an .h5ad built from this
# container by build_h5ad.py. Fingerprints written here are re-checked there, so
# a conversion bug cannot pass as agreement.
#
# Usage:
#   Rscript export_seurat.R <input.rds> <output.h5> [assay] [layer]

suppressMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(rhdf5)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: export_seurat.R <input.rds> <output.h5> [assay] [layer]")
src <- args[1]
out <- args[2]
assay_name <- if (length(args) >= 3) args[3] else "RNA"
layer_name <- if (length(args) >= 4) args[4] else "data"

put <- function(file, name, v) {
  mode <- if (is.character(v)) "character" else if (is.integer(v)) "integer" else "double"
  n <- length(v)
  if (mode == "character") {
    h5createDataset(file, name, dims = n, storage.mode = "character",
                    size = max(nchar(v)) + 1L, chunk = min(n, 4096L), level = 4)
  } else {
    h5createDataset(file, name, dims = n, storage.mode = mode,
                    chunk = min(n, 262144L), level = 4)
  }
  h5write(v, file, name)
}

o <- readRDS(src)
a <- o[[assay_name]]
d <- as(LayerData(a, layer_name), "CsparseMatrix")   # genes x cells
counts <- tryCatch(as(LayerData(a, "counts"), "CsparseMatrix"), error = function(e) NULL)
if (!is.null(counts)) stopifnot(identical(dimnames(d), dimnames(counts)))

if (file.exists(out)) file.remove(out)
h5createFile(out)
h5createGroup(out, "data"); h5createGroup(out, "obs"); h5createGroup(out, "fingerprint")
put(out, "data/i", d@i); put(out, "data/p", d@p); put(out, "data/x", d@x)
put(out, "data/shape", dim(d))
put(out, "features", rownames(d))
put(out, "barcodes", colnames(d))
if (!is.null(counts)) put(out, "total_counts", Matrix::colSums(counts))

md <- o@meta.data
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (is.numeric(v) || is.character(v) || is.factor(v) || is.logical(v)) {
    put(out, paste0("obs/", cn), as.character(v))
  }
}

# Fingerprints the Python side must reproduce before any comparison is scored.
put(out, "fingerprint/gene_sums", Matrix::rowSums(d))
put(out, "fingerprint/cell_sums", Matrix::colSums(d))
put(out, "fingerprint/nnz", length(d@x))
H5close()
cat(sprintf("wrote %s  genes=%d cells=%d nnz=%d\n",
            out, nrow(d), ncol(d), length(d@x)))
