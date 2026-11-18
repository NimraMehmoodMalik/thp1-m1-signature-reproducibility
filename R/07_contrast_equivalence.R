# =====================================================================
# 07_contrast_equivalence.R
# Direct contrast between the primary and the THP-1 effect of each primary core gene, and equivalence test.
#
# Input : results/meta/conserved_vs_specific.csv, results/meta/meta_thp1_M1.csv (06_thp1_meta_classes.R)
# Output: results/meta/contrast_equivalence.csv  06 table plus delta, se_diff and eq_class
#         results/meta/equivalence_summary.csv   equivalence classes for margins 0.5, 1 and 1.5 log2 and
#                                                THP-1 standard errors inflated by 1 and 2
#
# delta = primary - THP-1 meta log2FC, se_diff = sqrt(se_prim^2 + se_thp^2). Main analysis: margin 1 log2, factor 1.
# =====================================================================
source("R/00_config.R")
source("R/meta_functions.R")

cm <- read.csv(file.path(PATH$results, "meta", "conserved_vs_specific.csv"), colClasses = c(gene = "character"))
mt <- read.csv(file.path(PATH$results, "meta", "meta_thp1_M1.csv"), colClasses = c(gene = "character"))

cm$delta <- cm$beta_prim - cm$beta_thp; cm$se_diff <- sqrt(cm$se_prim^2 + cm$se_thp^2)
core_cls <- c("conserved", "weaker_in_THP1", "not_reproduced_in_THP1", "opposite")
iscore <- cm$cls %in% core_cls; nc <- sum(iscore)
cm$eq_class <- equivalence_class(cm, 1, 1)

summ <- do.call(rbind, lapply(c(0.5, 1, 1.5), function(m) do.call(rbind, lapply(c(1, 2), function(f) {
  k <- factor(equivalence_class(cm, m, f)[iscore], levels = c("equivalent", "inconclusive", "different"))
  data.frame(margin_log2 = m, se_thp_factor = f, t(as.matrix(table(k))), pct_equivalent = round(100 * mean(k == "equivalent"), 1))
}))))
write.csv(summ, file.path(PATH$results, "meta", "equivalence_summary.csv"), row.names = FALSE)
cat("=== EQUIVALENCE SUMMARY: primary core genes tested in THP-1 (n =", nc, ") ===\n"); print(summ, row.names = FALSE)
cat("\n=== Margin 1, factor 1, by threshold-based class ===\n")
print(table(class = cm$cls[iscore], eq_class = cm$eq_class[iscore]))

dsig <- cm$eq_class == "different"; eq0 <- (abs(cm$beta_prim) + 1.645 * cm$se_prim) < 0.5; sp <- cm$cls == "THP1_specific"
cat("\nTHP1_specific genes:", sum(sp), "| significantly different from primary:", sum(sp & dsig),
    "| primary effect equivalent to 0:", sum(sp & eq0), "| both (strict THP-1-specific):", sum(sp & dsig & eq0), "\n")

r <- cm[cm$gene %in% READOUTS, c("gene", "beta_prim", "beta_thp", "delta", "se_diff")]
r$CI90_low <- r$delta - 1.645 * r$se_diff; r$CI90_high <- r$delta + 1.645 * r$se_diff
r$eq_class <- cm$eq_class[cm$gene %in% READOUTS]
r[, 2:7] <- signif(r[, 2:7], 3)
cat("\n=== READOUTS: primary vs THP-1, difference with 90% CI ===\n"); print(r, row.names = FALSE)
cat("\nPer-dataset THP-1 log2FC of the readouts:\n")
print(round(mt[mt$gene %in% READOUTS, c("gene", grep("^lfc_", names(mt), value = TRUE))][, -1], 2), digits = 3)
write.csv(cm, file.path(PATH$results, "meta", "contrast_equivalence.csv"), row.names = FALSE)
