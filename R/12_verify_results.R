# =====================================================================
# 12_verify_results.R
# Checks the main results against the values reported in the manuscript and stops with an error if any differ.
#
# Input : results/ (scripts 04 to 10)
# Output: results/verification.csv   expected and observed value of each check
#
# The checks assume the main analysis settings of 00_config.R (THP1_MODEL = "RE", PSEUDOCOUNT = 0.1).
# =====================================================================
source("R/00_config.R")
meta <- function(f, ...) read.csv(file.path(PATH$results, "meta", f), colClasses = c(gene = "character"), ...)
ds   <- read.csv(file.path(PATH$results, "de_summary.csv"))
mp <- meta("meta_primary_M1.csv"); mt <- meta("meta_thp1_M1.csv"); cm <- meta("conserved_vs_specific.csv")
eq <- meta("contrast_equivalence.csv"); rc <- meta("robust_core.csv"); ch <- read.csv(file.path(PATH$metadata, "chen2019_LPS_IFNg_vs_M0.csv"), colClasses = c(Gene = "character"))
core_cls <- c("conserved", "weaker_in_THP1", "not_reproduced_in_THP1", "opposite")
e1 <- table(factor(eq$eq_class[eq$cls %in% core_cls], levels = c("equivalent", "inconclusive", "different")))
m  <- merge(mp[, c("gene", "beta")], ch[, c("Gene", "mu")], by.x = "gene", by.y = "Gene")

acc8 <- c("GSE61298", "GSE57614", "GSE162698", "GSE157182", "GSE159112", "GSE234829", "E-MTAB-5917", "GSE181851")
chk <- data.frame(check = c(
    paste("genes tested in DE:", acc8),
    "primary meta: genes tested", "primary meta: core genes", "primary meta: core up", "primary meta: core down",
    "THP-1 meta: genes tested", "primary vs THP-1: genes in both",
    "class conserved", "class weaker_in_THP1", "class not_reproduced_in_THP1", "class opposite", "class THP1_specific",
    "equivalent (margin 1 log2)", "inconclusive (margin 1 log2)", "different (margin 1 log2)",
    "robust core genes", "Spearman primary vs THP-1 (3 decimals)", "Spearman primary vs Chen 2019 (3 decimals)"),
  expected = c(19460, 19553, 11599, 16256, 12292, 11257, 21355, 14928,
               13196, 1741, 894, 847, 12234, 10434, 503, 365, 543, 29, 209, 130, 1072, 238, 1238, 0.573, 0.811),
  observed = c(ds$genes[match(acc8, ds$accession)],
               nrow(mp), sum(mp$core), sum(mp$core & mp$beta > 0), sum(mp$core & mp$beta < 0), nrow(mt), nrow(cm),
               sum(cm$cls == "conserved"), sum(cm$cls == "weaker_in_THP1"), sum(cm$cls == "not_reproduced_in_THP1"),
               sum(cm$cls == "opposite"), sum(cm$cls == "THP1_specific"), e1[["equivalent"]], e1[["inconclusive"]], e1[["different"]],
               nrow(rc), round(cor(cm$beta_prim, cm$beta_thp, method = "spearman"), 3), round(cor(m$beta, m$mu, method = "spearman"), 3)))
chk$ok <- abs(chk$expected - chk$observed) < 1e-9
write.csv(chk, file.path(PATH$results, "verification.csv"), row.names = FALSE)
print(chk, row.names = FALSE)
if (!all(chk$ok)) stop("Verification failed for: ", paste(chk$check[!chk$ok], collapse = "; "))
cat("\nAll", nrow(chk), "checks passed.\n")
