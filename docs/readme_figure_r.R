#!/usr/bin/env Rscript
# Run R tricycle for column 5 of the README figure.
#
# docs/make_readme_figure.py writes <dir>/logcounts.mtx (genes by cells,
# log-normalised), <dir>/genes.txt and <dir>/cells.txt, then calls this
# script. It writes <dir>/r_position.tsv: cell, theta.
#
# Needs the Bioconductor packages tricycle and SingleCellExperiment.
#
# Usage: Rscript readme_figure_r.R <dir>

suppressMessages({
  library(Matrix); library(SingleCellExperiment); library(tricycle)
})

dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(dir)) stop("usage: readme_figure_r.R <dir>")
mat <- as(readMM(file.path(dir, "logcounts.mtx")), "CsparseMatrix")
rownames(mat) <- readLines(file.path(dir, "genes.txt"))
colnames(mat) <- readLines(file.path(dir, "cells.txt"))

sce <- SingleCellExperiment(assays = list(logcounts = mat))
sce <- estimate_cycle_position(sce, species = "human", gname.type = "SYMBOL")

write.table(data.frame(cell = colnames(mat), theta = sprintf("%.17g", sce$tricyclePosition)),
            file.path(dir, "r_position.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("R tricycle", as.character(packageVersion("tricycle")), ":", ncol(mat), "cells\n")
