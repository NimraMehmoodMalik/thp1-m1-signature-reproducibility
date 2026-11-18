# =====================================================================
# 08_sensitivity_robust_core.R
# Sensitivity of the conserved fraction to the analysis choices, robust core, benchmark against an
# independent published meta-analysis, and dataset-level ranking.
#
# Input : results/meta/ (05-07), results/de/, metadata/chen2019_LPS_IFNg_vs_M0.csv
# Output: results/meta/sensitivity.csv        conserved share at THP-1 effect-size thresholds 1, 0.75, 0.5
#         results/meta/robust_core.csv        core genes that remain core when GSE57614 is left out
#         results/meta/benchmark_chen.csv     direction and significance of the full and robust core in Chen et al. 2019
#         results/meta/dataset_ranking.csv    agreement of each dataset with the primary core
#
# GSE57614 is a two-colour Agilent experiment (log-ratio data), so the core is also derived without it.
# Benchmark: Chen HJ et al., Front Immunol 2019;10:2887 (doi:10.3389/fimmu.2019.02887), random-effects
# meta-analysis of in vitro macrophages, M-LPS+IFNg versus M0 (mu = effect, Z_pvalBH = BH-adjusted P).
# Dataset ranking: a primary dataset is compared with the core and meta-analysis of the OTHER three primary
# datasets, so that the comparison is not circular.
# =====================================================================
source("R/00_config.R")
source("R/meta_functions.R")

de   <- load_de()
prim <- DATASETS$accession[DATASETS$role == "primary"]
thp  <- DATASETS$accession[DATASETS$role == "THP1"]
mp <- read.csv(file.path(PATH$results, "meta", "meta_primary_M1.csv"), colClasses = c(gene = "character"), check.names = FALSE)
cm <- read.csv(file.path(PATH$results, "meta", "conserved_vs_specific.csv"), colClasses = c(gene = "character"), check.names = FALSE)
ch <- read.csv(file.path(PATH$metadata, "chen2019_LPS_IFNg_vs_M0.csv"), colClasses = c(Gene = "character"), check.names = FALSE)
outdir <- file.path(PATH$results, "meta")

# ---- A) conserved share at different THP-1 thresholds -------------------------
ps     <- cm$core
ts_any <- cm$FDR_thp < 0.05 & cm$concordant & cm$n_sig >= 2
same   <- sign(cm$beta_prim) == sign(cm$beta_thp)
sens <- do.call(rbind, lapply(c(1, 0.75, 0.5), function(t)
  data.frame(analysis = paste0("THP-1 |log2FC| >= ", t), primary_core = sum(ps),
             pct_conserved = round(100 * mean((ts_any & same & abs(cm$beta_thp) >= t)[ps]), 1),
             pct_opposite = round(100 * mean((ts_any & !same)[ps]), 1))))
write.csv(sens, file.path(outdir, "sensitivity.csv"), row.names = FALSE)
cat("=== A) SENSITIVITY OF THE CONSERVED SHARE ===\n"); print(sens, row.names = FALSE)

# ---- B) robust core --------------------------------------------------------------
m3 <- re_meta_common(de, setdiff(prim, "GSE57614"))
core4 <- mp$gene[mp$core]; core3 <- m3$gene[m3$core]
robust <- intersect(core4, core3)
write.csv(data.frame(gene = robust), file.path(outdir, "robust_core.csv"), row.names = FALSE)
cat("\n=== B) ROBUST CORE ===\nfull core:", length(core4), "| core without GSE57614:", length(core3),
    "| overlap:", length(robust), "| % of the full core retained:", round(100 * length(robust) / length(core4), 1), "\n")
cm$robust <- cm$gene %in% robust
lv <- c("conserved", "weaker_in_THP1", "not_reproduced_in_THP1", "opposite")
cls_tab <- rbind(full_core = table(factor(cm$cls[cm$core], levels = lv)), robust_core = table(factor(cm$cls[cm$robust], levels = lv)))
print(cls_tab)
cat("% conserved: full", round(100 * cls_tab[1, 1] / sum(cls_tab[1, ]), 1), "| robust", round(100 * cls_tab[2, 1] / sum(cls_tab[2, ]), 1), "\n")
eq <- read.csv(file.path(outdir, "contrast_equivalence.csv"), colClasses = c(gene = "character"), check.names = FALSE)
e1 <- eq[eq$gene %in% robust & eq$cls %in% lv, ]
cat("Robust core genes tested in THP-1:", nrow(e1), "\nequivalence classes (margin 1 log2):\n")
print(table(e1$eq_class)); print(round(100 * prop.table(table(e1$eq_class)), 1))

# ---- C) benchmark against Chen et al. 2019 ----------------------------------------
m <- merge(mp[, c("gene", "beta", "core")], ch[, c("Gene", "mu", "Z_pvalBH")], by.x = "gene", by.y = "Gene")
m$robust <- m$gene %in% robust; m$chen_sig <- m$Z_pvalBH < 0.05
bench <- function(sel, label) {
  d <- m[sel, ]
  data.frame(set = label, n_in_chen = nrow(d), same_direction_pct = round(100 * mean(sign(d$beta) == sign(d$mu)), 1),
             chen_sig_same_dir_pct = round(100 * mean(sign(d$beta) == sign(d$mu) & d$chen_sig), 1))
}
b <- rbind(bench(m$core, "full core"), bench(m$robust, "robust core"))
write.csv(b, file.path(outdir, "benchmark_chen.csv"), row.names = FALSE)
cat("\n=== C) BENCHMARK vs Chen et al. 2019 ===\ngenes in both:", nrow(m), "| Spearman (all common genes):",
    round(cor(m$beta, m$mu, method = "spearman"), 3), "\n")
print(b, row.names = FALSE)
ft <- fisher.test(table(core = m$core, chen_sig = m$chen_sig))
cat("Fisher OR (core vs Chen FDR<0.05):", round(ft$estimate, 1), " P =", signif(ft$p.value, 2), "\n")
print(m[m$gene %in% READOUTS, c("gene", "beta", "mu", "Z_pvalBH")], row.names = FALSE)

# ---- D) dataset-level ranking -----------------------------------------------------
ref_full <- setNames(mp$beta, mp$gene)
rows <- list()
add_row <- function(name, type, ref, cg = mp$gene[mp$core]) {
  t <- de[[name]]; g <- intersect(intersect(cg, names(ref)), t$gene); x <- ref[g]; y <- t$logFC[match(g, t$gene)]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  rows[[length(rows) + 1]] <<- data.frame(dataset = sub("\\|.*", "", name), type = type, core_genes_tested = length(x),
    spearman = round(cor(x, y, method = "spearman"), 3), same_direction_pct = round(100 * mean(sign(x) == sign(y)), 1),
    slope = round(unname(coef(lm(y ~ x))[2]), 2))
}
for (p in prim) { mo <- re_meta_common(de, setdiff(prim, p)); add_row(de_key(p), "primary (vs other 3)", setNames(mo$beta, mo$gene), cg = mo$gene[mo$core]) }
for (t in thp) add_row(de_key(t), "THP-1 (vs primary meta)", ref_full)
add_row("GSE181851|LPS-ctrl", "THP-1, LPS only (vs primary meta)", ref_full)
rk <- do.call(rbind, rows); rk <- rk[order(-rk$spearman), ]
write.csv(rk, file.path(outdir, "dataset_ranking.csv"), row.names = FALSE)
cat("\n=== D) DATASET-LEVEL RANKING (core selected without the tested dataset) ===\n"); print(rk, row.names = FALSE)
