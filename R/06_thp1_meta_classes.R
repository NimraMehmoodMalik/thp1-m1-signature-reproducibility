# =====================================================================
# 06_thp1_meta_classes.R
# THP-1 meta-analysis (GSE159112, GSE234829, E-MTAB-5917) and classification of every gene by how its
# primary macrophage effect is reproduced in THP-1; readout genes; LPS-only THP-1 dataset as a directional check.
#
# Input : results/de/, results/meta/meta_primary_M1.csv (05_primary_meta.R)
# Output: results/meta/meta_thp1_M1.csv           THP-1 meta-analysis (model THP1_MODEL in 00_config.R)
#         results/meta/conserved_vs_specific.csv  primary and THP-1 estimates and class per gene
#         results/meta/readout_classes.csv        the four readout genes
#         results/meta/lps_only_check.csv         GSE181851 (LPS only) against the primary meta-analysis
#
# THP-1 side of a gene: "any" = Stouffer FDR < 0.05, same direction in all datasets and nominal P < 0.05 in
# both of the two or more datasets; "strong" = any and |RE log2FC| >= 1.
# Classes of primary core genes (same direction as primary):
#   conserved          THP-1 strong
#   weaker_in_THP1     THP-1 any but |log2FC| < 1
#   not_reproduced_in_THP1  THP-1 not significant
#   opposite           THP-1 significant in the opposite direction
# Other genes: THP1_specific (THP-1 strong, primary P > 0.05 and |log2FC| < 0.5),
#   THP1_strong_primary_subthreshold (THP-1 strong, primary not core but not absent), neither.
# =====================================================================
source("R/00_config.R")
source("R/meta_functions.R")

de   <- load_de()
thp  <- DATASETS$accession[DATASETS$role == "THP1"]
mp   <- read.csv(file.path(PATH$results, "meta", "meta_primary_M1.csv"), colClasses = c(gene = "character"))

mt <- run_meta(de, thp, 2, THP1_MODEL)
write.csv(mt, file.path(PATH$results, "meta", "meta_thp1_M1.csv"), row.names = FALSE)
cat("=== THP-1 META (", THP1_MODEL, ", k = ", length(thp), ") ===  genes: ", nrow(mt), " | FDR<0.05: ", sum(mt$FDR < 0.05), "\n", sep = "")

cm <- merge(mp[, c("gene", "beta", "se", "p", "FDR", "I2", "core")],
            mt[, c("gene", "beta", "se", "p", "FDR", "n_sig", "concordant")], by = "gene", suffixes = c("_prim", "_thp"))
cm$cls <- "neither"
ts_any    <- cm$FDR_thp < 0.05 & cm$concordant & cm$n_sig >= 2
ts_strong <- ts_any & abs(cm$beta_thp) >= 1
ps   <- cm$core
same <- sign(cm$beta_prim) == sign(cm$beta_thp)
cm$cls[ps & ts_strong & same]           <- "conserved"
cm$cls[ps & ts_any & !ts_strong & same] <- "weaker_in_THP1"
cm$cls[ps & !ts_any]                    <- "not_reproduced_in_THP1"
cm$cls[ps & ts_any & !same]             <- "opposite"
absent_prim <- cm$p_prim > 0.05 & abs(cm$beta_prim) < 0.5
cm$cls[!ps & ts_strong & absent_prim]   <- "THP1_specific"
cm$cls[!ps & ts_strong & !absent_prim]  <- "THP1_strong_primary_subthreshold"
write.csv(cm, file.path(PATH$results, "meta", "conserved_vs_specific.csv"), row.names = FALSE)

cat("\n=== PRIMARY vs THP-1: universe =", nrow(cm), "genes ===\n")
cat("Spearman (meta log2FC):", round(cor(cm$beta_prim, cm$beta_thp, method = "spearman"), 3),
    "| slope THP1~primary:", round(coef(lm(beta_thp ~ beta_prim, cm))[2], 2), "\n")
print(table(cm$cls))
ft <- fisher.test(table(ps, ts_strong & same))
cat("Fisher overlap (primary core vs THP-1 strong same-direction): OR =", round(ft$estimate, 1), " p =", signif(ft$p.value, 2), "\n")
cat("Share of primary core genes that are conserved:", round(100 * mean(cm$cls[ps] == "conserved"), 1), "% | conserved or weaker:",
    round(100 * mean(cm$cls[ps] %in% c("conserved", "weaker_in_THP1")), 1), "%\n")

# ---- readout genes ----------------------------------------------------------
ro <- cm[cm$gene %in% READOUTS, c("gene", "beta_prim", "FDR_prim", "I2", "beta_thp", "FDR_thp", "cls")]
ro$beta_prim <- round(ro$beta_prim, 2); ro$beta_thp <- round(ro$beta_thp, 2)
ro$FDR_prim <- signif(ro$FDR_prim, 2); ro$FDR_thp <- signif(ro$FDR_thp, 2); ro$I2 <- round(ro$I2)
write.csv(ro, file.path(PATH$results, "meta", "readout_classes.csv"), row.names = FALSE)
cat("\n=== READOUT GENES ===\n"); print(ro, row.names = FALSE)
cat("\nPer-dataset log2FC of the readouts, primary meta table:\n")
print(round(mp[mp$gene %in% READOUTS, c("gene", grep("^lfc_", names(mp), value = TRUE))][, -1], 2), digits = 3)

# ---- GSE181851 (THP-1, LPS only) against the primary meta-analysis ----------
l  <- de[["GSE181851|LPS-ctrl"]]
j  <- merge(mp[, c("gene", "beta", "core")], l[, c("gene", "logFC", "P")], by = "gene")
cj <- j[j$core, ]
lps <- data.frame(spearman_all = round(cor(j$beta, j$logFC, method = "spearman"), 3),
                  sign_concordance_core = round(mean(sign(cj$beta) == sign(cj$logFC)), 3),
                  pct_core_sig_same_sign = round(100 * mean(sign(cj$beta) == sign(cj$logFC) & cj$P < 0.05), 1))
write.csv(lps, file.path(PATH$results, "meta", "lps_only_check.csv"), row.names = FALSE)
cat("\n=== GSE181851 (THP-1, LPS only) vs primary meta ===\n"); print(lps, row.names = FALSE)
