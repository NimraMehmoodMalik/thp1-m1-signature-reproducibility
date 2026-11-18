# =====================================================================
# 10_enrichment.R
# Over-representation analysis (hypergeometric test) of the gene classes against MSigDB Hallmark and GO
# biological process gene sets.
#
# Input : results/meta/conserved_vs_specific.csv (06_thp1_meta_classes.R)
# Output: results/meta/enrichment_<group>.csv   set, size, overlap, p, FDR (BH), fold enrichment
#         groups: conserved, weaker_in_THP1, not_reproduced_or_opposite, THP1_specific
#                 (THP1_specific = THP1_specific + THP1_strong_primary_subthreshold)
#
# Universe: all genes tested in both the primary and the THP-1 meta-analysis. Gene sets with 10-500 genes
# in the universe are tested. Set contents depend on the msigdbr version, which is printed.
# =====================================================================
source("R/00_config.R")
source("R/meta_functions.R")

cm <- read.csv(file.path(PATH$results, "meta", "conserved_vs_specific.csv"), colClasses = c(gene = "character"))
sets <- msigdb_sets(); cat("msigdbr version:", as.character(utils::packageVersion("msigdbr")), "\n")
universe <- cm$gene
sets <- lapply(sets, function(s) intersect(s, universe)); sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]

ora <- function(genes_in) {
  genes_in <- intersect(genes_in, universe); N <- length(universe); n <- length(genes_in)
  res <- do.call(rbind, lapply(names(sets), function(s) {
    K <- length(sets[[s]]); k <- length(intersect(genes_in, sets[[s]]))
    data.frame(set = s, size = K, overlap = k, p = phyper(k - 1, K, N - K, n, lower.tail = FALSE))
  }))
  res$FDR <- p.adjust(res$p, "BH"); res$fold <- round((res$overlap / n) / (res$size / N), 2)
  res[order(res$p), ]
}
groups <- list(conserved = cm$gene[cm$cls == "conserved"],
               weaker_in_THP1 = cm$gene[cm$cls == "weaker_in_THP1"],
               not_reproduced_or_opposite = cm$gene[cm$cls %in% c("not_reproduced_in_THP1", "opposite")],
               THP1_specific = cm$gene[cm$cls %in% c("THP1_specific", "THP1_strong_primary_subthreshold")])
for (g in names(groups)) {
  r <- ora(groups[[g]])
  write.csv(r, file.path(PATH$results, "meta", sprintf("enrichment_%s.csv", g)), row.names = FALSE)
  cat(sprintf("\n=== %s (n = %d genes): top 8 sets ===\n", g, length(groups[[g]])))
  top <- head(r, 8); top$p <- signif(top$p, 2); top$FDR <- signif(top$FDR, 2)
  print(top[, c("set", "size", "overlap", "fold", "FDR")], row.names = FALSE)
}
