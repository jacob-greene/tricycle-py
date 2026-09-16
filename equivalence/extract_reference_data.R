# Extract the data objects tricycle (GPL-3) ships, plus the AnnotationDbi-derived
# identifier maps its ENSEMBL code paths depend on, into plain text.
suppressMessages({
  library(tricycle); library(AnnotationDbi)
  library(org.Hs.eg.db); library(org.Mm.eg.db)
})
out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

e <- new.env(); data("neuroRef", package = "tricycle", envir = e)
nr <- e[["neuroRef"]]
nr$rowname <- rownames(nr)
# write.table formats a double at 15 significant digits, which is lossy: the
# shipped weights then differ from R's by up to tens of units in the last
# place, and that difference propagates into every projected angle. %.17g is
# the shortest format that round-trips an IEEE 754 double exactly, so the
# shipped reference is bit-identical to the one R holds in memory.
g17 <- function(x) sprintf("%.17g", x)
write.table(data.frame(rowname = nr$rowname,
                       pc1.rot = g17(nr$pc1.rot), pc2.rot = g17(nr$pc2.rot),
                       ensembl = nr$ensembl, symbol = nr$symbol,
                       SYMBOL = nr$SYMBOL),
            file.path(out, "neuroRef.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

e2 <- new.env(); data("RevelioGeneList", package = "tricycle", envir = e2)
rg <- e2[["RevelioGeneList"]]
rl <- do.call(rbind, lapply(names(rg), function(n) data.frame(stage = n, gene = rg[[n]])))
write.table(rl, file.path(out, "RevelioGeneList.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ENSEMBL -> SYMBOL, exactly as AnnotationDbi::mapIds(..., multiVals = "first") gives it.
for (sp in c("human", "mouse")) {
  db <- if (sp == "human") org.Hs.eg.db else org.Mm.eg.db
  keys <- keys(db, keytype = "ENSEMBL")
  sym <- suppressMessages(AnnotationDbi::mapIds(db, keys = keys, column = "SYMBOL",
                                                keytype = "ENSEMBL", multiVals = "first"))
  df <- data.frame(ensembl = keys, symbol = unname(sym), stringsAsFactors = FALSE)
  df <- df[!is.na(df$symbol), ]
  gz <- gzfile(file.path(out, paste0("ensembl2symbol_", sp, ".tsv.gz")), "w")
  write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

  # GO:0007049 (cell cycle), the gene set run_pca_cc_genes subsets to.
  for (gt in c("SYMBOL", "ENSEMBL")) {
    g <- suppressMessages(AnnotationDbi::select(db, keytype = "GOALL", keys = "GO:0007049",
                                                columns = gt)[, gt])
    writeLines(as.character(g), file.path(out, paste0("go_cc_", sp, "_", tolower(gt), ".txt")))
  }
}

writeLines(c(
  paste0("tricycle\t", as.character(packageVersion("tricycle"))),
  paste0("org.Hs.eg.db\t", as.character(packageVersion("org.Hs.eg.db"))),
  paste0("org.Mm.eg.db\t", as.character(packageVersion("org.Mm.eg.db"))),
  paste0("R\t", R.version.string)), file.path(out, "PROVENANCE.tsv"))
cat("done\n")
