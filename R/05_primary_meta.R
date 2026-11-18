# =====================================================================
# 05_primary_meta.R
# Random-effects meta-analysis of the primary human macrophage datasets (M1 versus ctrl), core signature,
# and leave-one-dataset-out validation.
#
# Input : results/de/ (04_differential_expression.R)
# Output: results/meta/meta_primary_M1.csv   one row per gene: RE estimate (beta, se, tau2, I2), Stouffer
#                                            FDR, per-dataset log2FC (lfc_*), core flag
#         results/meta/loo_summary.csv       core signature derived without one dataset, tested in that dataset
#
# Primary datasets: GSE61298, GSE57614, GSE162698, GSE157182; genes present in at least 3 of the 4.
# Core gene: Stouffer FDR < 0.05, |RE log2FC| >= 1, same direction in all datasets, nominal P < 0.05 in at
# least max(2, k - 1) datasets.
# =====================================================================
source("R/00_config.R")
source("R/meta_functions.R")

de   <- load_de()
prim <- DATASETS$accession[DATASETS$role == "primary"]
dir.create(file.path(PATH$results, "meta"), showWarnings = FALSE, recursive = TRUE)

# ---- primary meta-analysis --------------------------------------------------
mp <- run_meta(de, prim, 3, "RE")
mp$core <- is_core(mp)
write.csv(mp, file.path(PATH$results, "meta", "meta_primary_M1.csv"), row.names = FALSE)
cat("=== PRIMARY META (M1 vs ctrl) ===\n")
cat("genes tested:", nrow(mp), "| FDR<0.05:", sum(mp$FDR < 0.05, na.rm = TRUE),
    "| core:", sum(mp$core), " (up", sum(mp$core & mp$beta > 0), "/ down", sum(mp$core & mp$beta < 0), ")\n")
cat("median I2 among core genes:", round(median(mp$I2[mp$core], na.rm = TRUE), 1), "%\n")

# ---- leave-one-dataset-out --------------------------------------------------
loo <- do.call(rbind, lapply(prim, function(h) {
  m <- run_meta(de, setdiff(prim, h), 3, "RE"); m$core <- is_core(m)
  j <- merge(m, de[[de_key(h)]][, c("gene", "logFC", "P", "FDR")], by = "gene", suffixes = c("", ".h"))
  cj <- j[j$core, ]
  data.frame(held_out = h, core_in_others = sum(m$core), core_testable_in_held_out = nrow(cj),
             spearman_all = round(cor(j$beta, j$logFC, method = "spearman"), 3),
             sign_concordance_core = round(mean(sign(cj$beta) == sign(cj$logFC)), 3),
             pct_core_sig_same_sign = round(100 * mean(sign(cj$beta) == sign(cj$logFC) & cj$P < 0.05), 1))
}))
write.csv(loo, file.path(PATH$results, "meta", "loo_summary.csv"), row.names = FALSE)
cat("\n=== LEAVE-ONE-DATASET-OUT ===\n"); print(loo, row.names = FALSE)
