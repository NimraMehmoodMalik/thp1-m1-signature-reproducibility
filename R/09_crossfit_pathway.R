# =====================================================================
# 09_crossfit_pathway.R
# Three robustness analyses of the primary-versus-THP-1 comparison.
#   A) Cross-fitting (winner's curse). The primary core is selected on |log2FC| >= 1 in the same data that are
#      later compared with THP-1. The four primary datasets are split 2 + 2 (3 splits, both directions = 6 runs):
#      the core is defined in one half, the primary effect is re-estimated in the other half, and that estimate
#      is compared with THP-1.
#   B) Pathway-level equivalence. The unit of replication is the dataset (4 primary versus 3 THP-1), so the
#      uncertainty comes from the spread between datasets (Welch 90% CI of the difference of dataset means).
#   C) Sensitivity to the THP-1 meta-analysis model (fixed effect versus random effects).
#
# Input : results/de/, results/meta/contrast_equivalence.csv (07_contrast_equivalence.R)
# Output: results/meta/crossfit_runs.csv, crossfit_summary.csv, pathway_equivalence.csv,
#         pathway_scores.csv, thp1_model_sensitivity.csv
# =====================================================================
source("R/00_config.R")
source("R/meta_functions.R")

de   <- load_de()
prim <- DATASETS$accession[DATASETS$role == "primary"]
thp  <- DATASETS$accession[DATASETS$role == "THP1"]
outdir <- file.path(PATH$results, "meta")
cm <- read.csv(file.path(outdir, "contrast_equivalence.csv"), colClasses = c(gene = "character"))
core_cls <- c("conserved", "weaker_in_THP1", "not_reproduced_in_THP1", "opposite")
pct <- function(x, lv) round(100 * mean(x == lv), 1)
# difference / equivalence class from effect estimates and standard errors (margin in log2 units)
classify <- function(ba, sa, bb, sb, margin = 1)
  equivalence_class(data.frame(beta_prim = ba, se_prim = sa, beta_thp = bb, se_thp = sb), margin, 1)
cat("THP-1 meta-analysis model:", THP1_MODEL, "|", nrow(cm), "genes tested in both systems;",
    sum(cm$cls %in% core_cls), "primary core genes testable in THP-1\n")

# ---- A) cross-fitting ---------------------------------------------------------------
cat("\n=============== A) CROSS-FITTED CONTRAST ===============\n")
splits <- list(c("GSE61298", "GSE57614"), c("GSE61298", "GSE162698"), c("GSE61298", "GSE157182"))
thp_cols <- cm[, c("gene", "beta_thp", "se_thp", "FDR_thp", "n_sig", "concordant", "beta_prim", "se_prim")]
runs <- list()
for (a in splits) {
  b <- setdiff(prim, a)
  for (dir in 1:2) {
    sel <- if (dir == 1) a else b; est <- if (dir == 1) b else a
    ms <- run_meta(de, sel, 2, "RE"); ms$core <- is_core(ms)     # core defined in the selection half only
    me <- run_meta(de, est, 2, "RE")                              # effect re-estimated in the other half only
    S <- data.frame(gene = ms$gene, b_sel = ms$beta, core_sel = ms$core)
    E <- data.frame(gene = me$gene, b_est = me$beta, se_est = me$se)
    J <- merge(merge(S, E, by = "gene"), thp_cols, by = "gene")
    J <- J[is.finite(J$b_est) & is.finite(J$se_est) & is.finite(J$beta_thp) & is.finite(J$se_thp) & is.finite(J$b_sel), ]
    J$cls_cf    <- classify(J$b_est,     J$se_est,  J$beta_thp, J$se_thp)
    J$cls_naive <- classify(J$beta_prim, J$se_prim, J$beta_thp, J$se_thp)
    ts_any <- J$FDR_thp < 0.05 & J$concordant & J$n_sig >= 2
    J$cons_cf    <- ts_any & abs(J$beta_thp) >= 1 & sign(J$b_est)     == sign(J$beta_thp)
    J$cons_naive <- ts_any & abs(J$beta_thp) >= 1 & sign(J$beta_prim) == sign(J$beta_thp)
    K <- J[J$core_sel, ]
    runs[[length(runs) + 1]] <- data.frame(
      select_in = paste(sel, collapse = "+"), estimate_in = paste(est, collapse = "+"),
      n_core_selected = sum(ms$core), n_testable_in_THP1 = nrow(K),
      mean_abs_logFC_selection = round(mean(abs(K$b_sel)), 2),
      mean_logFC_estimation_same_sign = round(mean(sign(K$b_sel) * K$b_est), 2),
      pct_same_direction_in_estimation_half = round(100 * mean(sign(K$b_sel) == sign(K$b_est)), 1),
      conserved_pct_naive = round(100 * mean(K$cons_naive), 1), conserved_pct_crossfit = round(100 * mean(K$cons_cf), 1),
      equiv_pct_naive = pct(K$cls_naive, "equivalent"), inconcl_pct_naive = pct(K$cls_naive, "inconclusive"), diff_pct_naive = pct(K$cls_naive, "different"),
      equiv_pct_crossfit = pct(K$cls_cf, "equivalent"), inconcl_pct_crossfit = pct(K$cls_cf, "inconclusive"), diff_pct_crossfit = pct(K$cls_cf, "different"))
  }
}
runs <- do.call(rbind, runs)
write.csv(runs, file.path(outdir, "crossfit_runs.csv"), row.names = FALSE)
num <- sapply(runs, is.numeric)
summ <- data.frame(measure = names(runs)[num], mean = round(colMeans(runs[, num]), 1),
                   min = round(apply(runs[, num], 2, min), 1), max = round(apply(runs[, num], 2, max), 1))
write.csv(summ, file.path(outdir, "crossfit_summary.csv"), row.names = FALSE)
cat("\nOne row per run (selection half -> estimation half):\n"); print(runs, row.names = FALSE)
cat("\nMean, min and max over the 6 runs:\n"); print(summ, row.names = FALSE)

# ---- B) pathway-level equivalence ------------------------------------------------------
cat("\n=============== B) PATHWAY-LEVEL EQUIVALENCE ===============\n")
sets <- msigdb_sets(); cat("msigdbr version:", as.character(utils::packageVersion("msigdbr")), "\n")
want <- c("Interferon-gamma response" = "HALLMARK_INTERFERON_GAMMA_RESPONSE", "Interferon-alpha response" = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "TNF-alpha signaling via NF-kB" = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Inflammatory response" = "HALLMARK_INFLAMMATORY_RESPONSE",
          "Antigen processing and presentation" = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
          "Peptide antigen assembly with MHC class II" = "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX")
miss <- want[!want %in% names(sets)]; if (length(miss)) cat("Gene sets not found (skipped):", paste(miss, collapse = ", "), "\n")
want <- want[want %in% names(sets)]
dsets <- c(prim, thp)
common <- Reduce(intersect, lapply(dsets, function(d) de[[de_key(d)]]$gene))      # genes measured in all 7 datasets
cat("Genes measured in all 7 datasets:", length(common), "\n")
score <- function(d, genes) { t <- de[[de_key(d)]]; mean(t$logFC[match(genes, t$gene)], na.rm = TRUE) }
welch <- function(x, y, margin) {
  m1 <- mean(x); m2 <- mean(y); v1 <- var(x); v2 <- var(y); n1 <- length(x); n2 <- length(y)
  se <- sqrt(v1 / n1 + v2 / n2); df <- se^4 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
  d <- m1 - m2; lo <- d - qt(0.95, df) * se; hi <- d + qt(0.95, df) * se
  c(delta = d, lo = lo, hi = hi, df = df, equivalent = as.numeric(lo > -margin & hi < margin))
}
pw <- do.call(rbind, lapply(names(want), function(nm) {
  g <- intersect(sets[[want[nm]]], common); x <- sapply(prim, score, genes = g); y <- sapply(thp, score, genes = g)
  w <- welch(x, y, 1); w5 <- welch(x, y, 0.5); gg <- cm[cm$gene %in% g, ]
  data.frame(pathway = nm, n_genes = length(g), primary_mean = round(mean(x), 2), thp1_mean = round(mean(y), 2),
             thp1_over_primary = round(mean(y) / mean(x), 2), delta = round(w["delta"], 2),
             CI90_low = round(w["lo"], 2), CI90_high = round(w["hi"], 2), welch_df = round(w["df"], 1),
             equivalent_margin_1 = w["equivalent"] == 1, equivalent_margin_0.5 = w5["equivalent"] == 1,
             range_primary_datasets = paste(round(range(x), 2), collapse = " to "), range_thp1_datasets = paste(round(range(y), 2), collapse = " to "),
             pct_genes_same_direction = round(100 * mean(sign(gg$beta_prim) == sign(gg$beta_thp)), 1),
             spearman_genes = round(cor(gg$beta_prim, gg$beta_thp, method = "spearman"), 2), row.names = NULL)
}))
write.csv(pw, file.path(outdir, "pathway_equivalence.csv"), row.names = FALSE)
print(pw, row.names = FALSE)
ps <- do.call(rbind, lapply(names(want), function(nm) {
  g <- intersect(sets[[want[nm]]], common); data.frame(pathway = nm, t(round(sapply(dsets, score, genes = g), 2))) }))
write.csv(ps, file.path(outdir, "pathway_scores.csv"), row.names = FALSE)
cat("\nPer-dataset pathway scores (mean log2FC of the set's genes):\n"); print(ps, row.names = FALSE)

# ---- C) THP-1 meta-analysis model ------------------------------------------------------
cat("\n=============== C) THP-1 META-ANALYSIS MODEL ===============\n")
alt_model <- if (THP1_MODEL == "RE") "FE" else "RE"
malt <- run_meta(de, thp, 2, alt_model)
C <- merge(cm[, c("gene", "cls", "beta_prim", "se_prim", "beta_thp", "se_thp", "eq_class")],
           malt[, c("gene", "beta", "se", "tau2", "I2")], by = "gene", suffixes = c("", "_alt"))
C$eq_class_alt <- classify(C$beta_prim, C$se_prim, C$beta, C$se)
ic <- C$cls %in% core_cls
count3 <- function(x) c(equivalent = sum(x[ic] == "equivalent"), inconclusive = sum(x[ic] == "inconclusive"), different = sum(x[ic] == "different"))
tab <- data.frame(THP1_meta = c(paste(THP1_MODEL, "(main)"), paste(alt_model, "(sensitivity)")), rbind(count3(C$eq_class), count3(C$eq_class_alt)))
tab$pct_equivalent <- round(100 * tab$equivalent / (tab$equivalent + tab$inconclusive + tab$different), 1)
tab$pct_different  <- round(100 * tab$different  / (tab$equivalent + tab$inconclusive + tab$different), 1)
write.csv(tab, file.path(outdir, "thp1_model_sensitivity.csv"), row.names = FALSE); print(tab, row.names = FALSE)
cat("median THP-1 SE:", THP1_MODEL, round(median(C$se_thp[ic]), 2), "|", alt_model, round(median(C$se[ic]), 2),
    "| median I2 among core genes in THP-1:", round(median(C$I2[ic], na.rm = TRUE), 1), "%\n")
r <- C[C$gene %in% READOUTS, ]; r$delta <- r$beta_prim - r$beta; r$se_d <- sqrt(r$se_prim^2 + r$se^2)
cat("\nReadouts with the", alt_model, "THP-1 estimate:\n")
print(data.frame(gene = r$gene, beta_prim = round(r$beta_prim, 2), beta_thp_alt = round(r$beta, 2), delta = round(r$delta, 2),
                 CI90_low = round(r$delta - 1.645 * r$se_d, 2), CI90_high = round(r$delta + 1.645 * r$se_d, 2),
                 class_main = r$eq_class, class_alt = r$eq_class_alt), row.names = FALSE)
