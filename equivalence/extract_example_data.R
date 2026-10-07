#!/usr/bin/env Rscript
# Export tricycle's own `neurosphere_example` so the Python package can ship it.
#
# The example is 400 mouse neurosphere cells, log-normalised, keyed on mouse
# Ensembl ids. It is the dataset the R vignette uses. The package ships it with
# the tricyclePosition that R computes on it, so a user can check this port
# against R without installing R.
#
# Usage:
#   Rscript extract_example_data.R <outdir>

suppressMessages({library(SingleCellExperiment); library(tricycle); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "neurosphere_example"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
g17 <- function(x) sprintf("%.17g", x)   # exact round trip for a double

env <- new.env(parent = emptyenv())
data("neurosphere_example", package = "tricycle", envir = env)
sce <- env[["neurosphere_example"]]
mat <- as(assay(sce, "logcounts"), "CsparseMatrix")       # genes x cells
ensembl <- rowData(sce)$Accession
symbol <- rowData(sce)$Gene
rownames(mat) <- ensembl

trip <- summary(mat)                                      # i, j, x; 1-based
con <- gzfile(file.path(outdir, "logcounts.tsv.gz"), "w")
writeLines(sprintf("%d\t%d", nrow(mat), ncol(mat)), con)  # header: genes cells
writeLines(sprintf("%d\t%d\t%s", trip$i - 1L, trip$j - 1L, g17(trip$x)), con)
close(con)

write.table(data.frame(ensembl = ensembl, symbol = symbol),
            file.path(outdir, "genes.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")

r <- estimate_cycle_position(sce, species = "mouse", gname.type = "ENSEMBL",
                             gname = ensembl)
write.table(data.frame(cell = colnames(sce),
                       tricyclePosition_R = g17(r$tricyclePosition)),
            file.path(outdir, "cells.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
cat(sprintf("wrote %s: %d genes, %d cells, %d non-zero; tricycle %s\n",
            outdir, nrow(mat), ncol(mat), length(trip$x),
            as.character(packageVersion("tricycle"))))
