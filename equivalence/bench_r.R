#!/usr/bin/env Rscript
# Time one tricycle function in R. Peak memory is measured outside, by
# /usr/bin/time, so the two sides are measured the same way.
#
# Usage:
#   Rscript bench_r.R <input.rds> <what> [replicate_factor]
#   what: load | position | schwabe | loess_direct | loess_interpolate

suppressMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(SingleCellExperiment); library(tricycle)
})

args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
what <- args[2]
reps <- if (length(args) >= 3) as.integer(args[3]) else 1L

o <- readRDS(src)
mat <- LayerData(o[["RNA"]], "data")
if (reps > 1L) {
  mat <- do.call(cbind, replicate(reps, mat, simplify = FALSE))
  colnames(mat) <- paste0(rep(seq_len(reps), each = ncol(mat) / reps), "_",
                          colnames(mat))
}
cat(sprintf("cells\t%d\n", ncol(mat)))

if (what == "load") {
  cat("seconds\t0\n")
  quit(save = "no")
}

t0 <- proc.time()[["elapsed"]]
if (what == "position") {
  sce <- SingleCellExperiment(assays = list(logcounts = mat))
  sce <- estimate_cycle_position(sce, species = "human", gname.type = "SYMBOL")
  invisible(sce$tricyclePosition)
} else if (what == "schwabe") {
  invisible(estimate_Schwabe_stage(mat, species = "human", gname.type = "SYMBOL"))
} else if (what %in% c("loess_direct", "loess_interpolate")) {
  sce <- SingleCellExperiment(assays = list(logcounts = mat))
  sce <- estimate_cycle_position(sce, species = "human", gname.type = "SYMBOL")
  y <- as.numeric(mat["TOP2A", ])
  t0 <- proc.time()[["elapsed"]]        # exclude the projection from this one
  if (what == "loess_direct") {
    invisible(fit_periodic_loess(sce$tricyclePosition, y,
                                 control = loess.control(surface = "direct")))
  } else {
    invisible(fit_periodic_loess(sce$tricyclePosition, y))
  }
} else {
  stop("unknown benchmark target")
}
cat(sprintf("seconds\t%.3f\n", proc.time()[["elapsed"]] - t0))
