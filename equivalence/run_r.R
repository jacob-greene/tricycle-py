#!/usr/bin/env Rscript
# Run the real R package tricycle on a Seurat .rds and write its outputs as text.
#
# Nothing here reimplements tricycle. Every number comes from the installed
# Bioconductor package, so the comparison is against the real thing.
#
# Usage:
#   Rscript run_r.R <input.rds> <outdir> [species] [gname_type] [batch_column]

suppressMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(SingleCellExperiment); library(tricycle)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run_r.R <input.rds> <outdir> [species] [gname_type] [batch_column]")
src <- args[1]
outdir <- args[2]
species <- if (length(args) >= 3) args[3] else "human"
gname_type <- if (length(args) >= 4) args[4] else "SYMBOL"
batch_col <- if (length(args) >= 5) args[5] else NA_character_
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

g17 <- function(x) sprintf("%.17g", x)   # exact round trip for a double

o <- readRDS(src)
mat <- LayerData(o[["RNA"]], "data")     # genes x cells, log-normalised
sce <- SingleCellExperiment(assays = list(logcounts = mat))
colData(sce)$barcode <- colnames(mat)

t0 <- proc.time()[["elapsed"]]
sce <- estimate_cycle_position(sce, species = species, gname.type = gname_type)
t_position <- proc.time()[["elapsed"]] - t0

emb <- reducedDim(sce, "tricycleEmbedding")
rot <- attr(emb, "rotation")
write.table(
  data.frame(barcode = colnames(mat),
             theta = g17(sce$tricyclePosition),
             pc1 = g17(emb[, 1]), pc2 = g17(emb[, 2])),
  file.path(outdir, "r_position.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(gene = rownames(rot), pc1_rot = g17(rot[, 1]), pc2_rot = g17(rot[, 2])),
  file.path(outdir, "r_projection_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Per-gene means over the matched genes: the quantity the projection centres by.
# The Python side must reproduce these before any angle is compared.
sub <- mat[rownames(rot), , drop = FALSE]
write.table(
  data.frame(gene = rownames(rot), mean = g17(Matrix::rowMeans(sub))),
  file.path(outdir, "r_gene_means.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

batch_v <- NULL
if (!is.na(batch_col)) batch_v <- as.character(o@meta.data[[batch_col]])
t0 <- proc.time()[["elapsed"]]
stage <- tryCatch(
  estimate_Schwabe_stage(mat, batch.v = batch_v, species = species,
                         gname.type = gname_type),
  error = function(e) { message("Schwabe failed: ", conditionMessage(e)); NULL })
t_stage <- proc.time()[["elapsed"]] - t0
if (!is.null(stage)) {
  write.table(data.frame(barcode = colnames(mat), stage = as.character(stage)),
              file.path(outdir, "r_stage.tsv"), sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
}

# Periodic loess, on both loess surfaces, for a gene that cycles strongly.
probe <- if ("TOP2A" %in% rownames(mat)) "TOP2A" else rownames(mat)[1]
y <- as.numeric(mat[probe, ])
theta <- sce$tricyclePosition
t0 <- proc.time()[["elapsed"]]
fit_interp <- fit_periodic_loess(theta, y)
t_loess <- proc.time()[["elapsed"]] - t0
fit_direct <- fit_periodic_loess(theta, y,
                                 control = loess.control(surface = "direct"))
write.table(
  data.frame(barcode = colnames(mat), y = g17(y),
             fitted_interpolate = g17(fit_interp$fitted),
             fitted_direct = g17(fit_direct$fitted)),
  file.path(outdir, "r_loess.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(x = g17(fit_direct$pred.df$x),
             y_interpolate = g17(fit_interp$pred.df$y),
             y_direct = g17(fit_direct$pred.df$y)),
  file.path(outdir, "r_loess_pred.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Reference building. scater::runPCA decides which genes enter and how they are
# centred; the decomposition's component signs are arbitrary on both sides.
pca_ok <- TRUE
tryCatch({
  go_sce <- run_pca_cc_genes(sce, species = species, gname.type = gname_type,
                             ntop = 500, ncomponents = 10)
  pca_rot <- attr(reducedDim(go_sce, "PCA"), "rotation")
  write.table(
    data.frame(gene = rownames(pca_rot),
               pc1 = g17(pca_rot[, 1]), pc2 = g17(pca_rot[, 2])),
    file.path(outdir, "r_pca_rotation.tsv"), sep = "\t", row.names = FALSE,
    quote = FALSE)
  write.table(
    data.frame(percent_var = g17(attr(reducedDim(go_sce, "PCA"), "percentVar"))),
    file.path(outdir, "r_pca_percentvar.tsv"), sep = "\t", row.names = FALSE,
    quote = FALSE)
}, error = function(e) {
  message("run_pca_cc_genes failed: ", conditionMessage(e)); pca_ok <<- FALSE
})

# Circular density, the computation under plot_ccposition_den.
dens <- circular::density.circular(circular::circular(theta), bw = 30)
write.table(data.frame(x = g17(as.numeric(dens$x)), y = g17(as.numeric(dens$y))),
            file.path(outdir, "r_density.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)

summary_lines <- c(
  paste0("probe_gene\t", probe),
  paste0("n_cells\t", ncol(mat)),
  paste0("n_genes\t", nrow(mat)),
  paste0("n_projection_genes\t", nrow(rot)),
  paste0("loess_rsquared_interpolate\t", g17(fit_interp$rsquared)),
  paste0("loess_rsquared_direct\t", g17(fit_direct$rsquared)),
  paste0("seconds_estimate_cycle_position\t", g17(t_position)),
  paste0("seconds_estimate_schwabe_stage\t", g17(t_stage)),
  paste0("seconds_fit_periodic_loess\t", g17(t_loess)),
  paste0("tricycle_version\t", as.character(packageVersion("tricycle"))),
  paste0("r_version\t", R.version.string))

# diagnose_totalumi needs raw counts, which are not always present.
#
# It is reported on BOTH loess surfaces. The function calls fit_periodic_loess,
# so on R's default it inherits the k-d tree approximation. Comparing the port's
# exact surface against R's approximate one is not like for like, and the gap it
# produces is R's approximation error, not a defect of the port. The direct row
# is the one to score; the interpolate row is the disclosure.
counts <- tryCatch(LayerData(o[["RNA"]], "counts"), error = function(e) NULL)
if (!is.null(counts)) {
  totalumis <- Matrix::colSums(counts)
  peak_valley <- function(pred) {
    pk <- max(pred$y[(pred$x > 0.75 * pi) & (pred$x < 1.25 * pi)])
    vl <- min(pred$y[(pred$x > 1.35 * pi) & (pred$x < 1.85 * pi)])
    c(pk, vl)
  }
  fit_i <- fit_periodic_loess(theta, log2(totalumis + 1))
  fit_d <- fit_periodic_loess(theta, log2(totalumis + 1),
                              control = loess.control(surface = "direct"))
  pv_i <- peak_valley(fit_i$pred.df)
  pv_d <- peak_valley(fit_d$pred.df)
  summary_lines <- c(summary_lines,
                     paste0("diagnose_peak\t", g17(pv_i[1])),
                     paste0("diagnose_valley\t", g17(pv_i[2])),
                     paste0("diagnose_difference\t", g17(pv_i[1] - pv_i[2])),
                     paste0("diagnose_peak_direct\t", g17(pv_d[1])),
                     paste0("diagnose_valley_direct\t", g17(pv_d[2])),
                     paste0("diagnose_difference_direct\t", g17(pv_d[1] - pv_d[2])))
}
writeLines(summary_lines, file.path(outdir, "r_summary.tsv"))
cat("R side done:", outdir, "\n")
