#!/usr/bin/env Rscript
# Compare the identifier code paths the main run does not reach.
#
# The main equivalence run uses species="human", gname.type="SYMBOL". Three
# other combinations exist, and each resolves reference gene names differently.
# Two of them call AnnotationDbi at run time, which is exactly where this port
# substitutes a frozen extract. That substitution has to be measured, not
# assumed.
#
# This script writes, for each arm, the exact log-expression matrix R used and
# the cell cycle positions R produced. check_gname_paths.py reads both and
# compares. The matrix travels with the answer, so neither side can be
# comparing a different input.
#
# Arms:
#   mouse_ensembl  neurosphere_example, the mouse SCE tricycle ships
#   mouse_symbol   the same cells, keyed on mouse gene symbols
#   human_ensembl  human cells relabelled to the Ensembl ids of their symbols
#
# Usage:
#   Rscript check_gname_paths.R <outdir> [human_rds]

suppressMessages({
  library(Matrix); library(rhdf5); library(SingleCellExperiment)
  library(tricycle); library(AnnotationDbi); library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
human_rds <- if (length(args) >= 2) args[2] else NA_character_
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
g17 <- function(x) sprintf("%.17g", x)

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

write_arm <- function(arm, mat, species, gname_type) {
  sce <- SingleCellExperiment(assays = list(logcounts = mat))
  sce <- estimate_cycle_position(sce, species = species, gname.type = gname_type)
  emb <- reducedDim(sce, "tricycleEmbedding")
  rot <- attr(emb, "rotation")

  h5 <- file.path(outdir, paste0(arm, ".h5"))
  if (file.exists(h5)) file.remove(h5)
  h5createFile(h5); h5createGroup(h5, "data")
  d <- as(mat, "CsparseMatrix")
  put(h5, "data/i", d@i); put(h5, "data/p", d@p); put(h5, "data/x", d@x)
  put(h5, "data/shape", dim(d))
  put(h5, "features", rownames(d)); put(h5, "barcodes", colnames(d))
  H5close()

  write.table(
    data.frame(barcode = colnames(mat), theta = g17(sce$tricyclePosition),
               pc1 = g17(emb[, 1]), pc2 = g17(emb[, 2])),
    file.path(outdir, paste0(arm, "_position.tsv")), sep = "\t",
    row.names = FALSE, quote = FALSE)
  writeLines(c(paste0("arm\t", arm), paste0("species\t", species),
               paste0("gname_type\t", gname_type),
               paste0("n_cells\t", ncol(mat)), paste0("n_genes\t", nrow(mat)),
               paste0("n_projection_genes\t", nrow(rot))),
             file.path(outdir, paste0(arm, "_meta.tsv")))
  cat(sprintf("%s: %d cells, %d projection genes\n", arm, ncol(mat), nrow(rot)))
}

# --- mouse, from the SCE tricycle ships --------------------------------------
env <- new.env(parent = emptyenv())
data("neurosphere_example", package = "tricycle", envir = env)
neuro <- env[["neurosphere_example"]]
neuro_mat <- assay(neuro, "logcounts")
rownames(neuro_mat) <- rowData(neuro)$Accession        # mouse Ensembl ids
write_arm("mouse_ensembl", neuro_mat, "mouse", "ENSEMBL")

neuro_sym <- neuro_mat
rownames(neuro_sym) <- rowData(neuro)$Gene             # mouse symbols
keep <- !is.na(rownames(neuro_sym)) & !duplicated(rownames(neuro_sym))
write_arm("mouse_symbol", neuro_sym[keep, ], "mouse", "SYMBOL")

# --- human, relabelled to Ensembl so the AnnotationDbi path runs --------------
# Only symbols with exactly one Ensembl id are kept, so the relabelling is a
# bijection and nothing is decided by tie-breaking.
if (!is.na(human_rds)) {
  suppressMessages({library(Seurat); library(SeuratObject)})
  o <- readRDS(human_rds)
  hmat <- LayerData(o[["RNA"]], "data")
  map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = rownames(hmat), keytype = "SYMBOL", columns = "ENSEMBL"))
  map <- map[!is.na(map$ENSEMBL), ]
  counts <- table(map$SYMBOL)
  unique_sym <- names(counts)[counts == 1]
  map <- map[map$SYMBOL %in% unique_sym & !duplicated(map$ENSEMBL), ]
  hsub <- hmat[map$SYMBOL, , drop = FALSE]
  rownames(hsub) <- map$ENSEMBL
  write_arm("human_ensembl", hsub, "human", "ENSEMBL")
}

cat("R side done:", outdir, "\n")
