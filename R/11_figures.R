# =====================================================================
# 11_figures.R
# Publication figures 2 to 7 (ggplot2) from the result tables of scripts 05 to 10.
#
# Input : results/meta/ (contrast_equivalence, meta_primary_M1, meta_thp1_M1, loo_summary, dataset_ranking,
#         robust_core, enrichment_<group>)
# Output: figures/Figure<n>_*.png (600 dpi) and figures/Figure<n>_*.pdf (vector)
#
# Colour-blind-safe Okabe-Ito palette. Each figure is wrapped in tryCatch, so one failure does not stop the others.
# =====================================================================
source("R/00_config.R")
require_packages(cran = c("ggplot2", "patchwork", "ggrepel", "scales"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ggrepel); library(scales) })
OUT <- PATH$figures; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rdf <- function(f) read.csv(file.path(PATH$results, "meta", f), check.names = FALSE, colClasses = c(gene = "character"))
cm <- rdf("contrast_equivalence.csv"); mp <- rdf("meta_primary_M1.csv"); mt <- rdf("meta_thp1_M1.csv")
loo <- read.csv(file.path(PATH$results, "meta", "loo_summary.csv")); rk <- read.csv(file.path(PATH$results, "meta", "dataset_ranking.csv"))
robust <- rdf("robust_core.csv")$gene

theme_set(theme_classic(base_size = 9, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 9.5), plot.tag = element_text(face = "bold", size = 12),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5)))
lab <- c(conserved = "Conserved", weaker_in_THP1 = "Weaker in THP-1", not_reproduced_in_THP1 = "Not reproduced",
         opposite = "Opposite", THP1_specific = "THP-1-specific", THP1_strong_primary_subthreshold = "THP-1-strong, primary sub-threshold", neither = "Neither")
pal <- c(conserved = "#009E73", weaker_in_THP1 = "#56B4E9", not_reproduced_in_THP1 = "#E69F00", opposite = "#CC79A7",
         THP1_specific = "#5B4B9A", THP1_strong_primary_subthreshold = "#B8B0DC", neither = "grey85")
core_cls <- c("conserved", "weaker_in_THP1", "not_reproduced_in_THP1", "opposite")
rdg <- c("CXCL8", "IL1B", "IL6", "TNF")
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, device = "png", type = "cairo")
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  cat("saved", name, "\n") }
run <- function(name, expr) tryCatch(expr, error = function(e) cat("!! figure", name, "failed:", conditionMessage(e), "\n"))

# shared numbers -------------------------------------------------------
cm$cls <- factor(cm$cls, levels = names(lab))
delta <- cm$beta_prim - cm$beta_thp; se_d <- sqrt(cm$se_prim^2 + cm$se_thp^2)
fdr_d <- p.adjust(2 * pnorm(-abs(delta / se_d)), "BH"); diff_ <- fdr_d < 0.05 & abs(delta) >= 1
isc <- cm$cls %in% core_cls

# ---- Figure 2: a) primary vs THP-1 scatter, b) tile matrix: threshold class x equivalence class ----------------
run("Figure2", {
  rho <- cor(cm$beta_prim, cm$beta_thp, method = "spearman")
  pa <- ggplot(cm, aes(beta_prim, beta_thp)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3) +
    geom_point(data = subset(cm, cls == "neither"), colour = "grey85", size = 0.35, stroke = 0) +
    geom_point(data = subset(cm, cls != "neither"), aes(colour = cls), size = 0.7, stroke = 0, alpha = 0.9) +
    geom_point(data = subset(cm, gene %in% rdg), shape = 21, size = 2.6, stroke = 0.5, fill = NA, colour = "black") +
    geom_label(data = transform(cm[match(rdg, cm$gene), ], lx = beta_prim + c(-0.2, 0.0, 0.45, 0.0), ly = beta_thp + c(1.4, 1.2, 0.0, -1.2),
                                hj = c(1, 0.5, 0, 0.5)), aes(x = lx, y = ly, label = gene, hjust = hj), fontface = "italic", size = 3,
               fill = alpha("white", 0.85), label.size = 0, label.padding = unit(0.6, "mm")) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.9, lineheight = 1.05,
             label = sprintf("Spearman ρ = %.2f\n%s genes tested in both systems", rho, format(nrow(cm), big.mark = ","))) +
    scale_colour_manual(values = pal, labels = lab[setdiff(names(lab), "neither")], name = NULL, breaks = setdiff(names(lab), "neither")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1), ncol = 2, byrow = TRUE)) +
    labs(x = "Primary macrophages, meta-analysed log2FC", y = "THP-1, meta-analysed log2FC", title = "M1 response: primary vs THP-1") +
    theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 7))
  dc <- droplevels(cm[isc, ]); dc$eq <- factor(dc$eq_class, levels = c("equivalent", "inconclusive", "different"))
  tb <- as.data.frame(table(cls = dc$cls, eq = dc$eq)); tb$pct <- ave(tb$Freq, tb$cls, FUN = function(x) 100 * x / sum(x))
  tot <- table(dc$cls); tb$row <- factor(paste0(lab[as.character(tb$cls)], "\n(n = ", tot[as.character(tb$cls)], ")"),
                                         levels = rev(paste0(lab[core_cls], "\n(n = ", tot[core_cls], ")")))
  pb <- ggplot(tb, aes(eq, row, fill = pct)) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", Freq, pct), colour = pct > 55), size = 2.9, lineheight = 0.95) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F3F7FB", high = "#0072B2", limits = c(0, 100), breaks = c(0, 50, 100), name = "% of class", guide = guide_colourbar(barwidth = unit(28, "mm"), barheight = unit(2.5, "mm"))) +
    scale_x_discrete(position = "top", labels = c("Equivalent", "Inconclusive", "Different")) +
    labs(x = NULL, y = NULL, title = "Threshold class vs equivalence test") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
  save_fig(pa + pb + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "a"), "Figure2_scatter_and_classes", 185, 100) })

# ---- Figure 3: a) conserved share vs THP-1 effect-size cut-off, b) equivalence vs margin (stacked area) -----------
run("Figure3", {
  ts_any <- cm$FDR_thp < 0.05 & cm$concordant & cm$n_sig >= 2; same <- sign(cm$beta_prim) == sign(cm$beta_thp)
  thr <- seq(0.25, 1.5, by = 0.05); rob <- cm$core & cm$gene %in% robust
  cv <- do.call(rbind, lapply(thr, function(t) { hit <- ts_any & same & abs(cm$beta_thp) >= t
    rbind(data.frame(thr = t, pct = 100 * mean(hit[cm$core]), set = "Full core"), data.frame(thr = t, pct = 100 * mean(hit[rob]), set = "Robust core (without GSE57614)")) }))
  pk <- cv[round(cv$thr, 2) %in% c(0.5, 0.75, 1), ]
  pa <- ggplot(cv, aes(thr, pct, colour = set)) + geom_line(linewidth = 0.7) +
    geom_point(data = pk, size = 1.8) +
    geom_text(data = subset(pk, set == "Robust core (without GSE57614)"), aes(label = sprintf("%.1f%%", pct)), vjust = -1.1, hjust = 0.4, size = 2.6, show.legend = FALSE) +
    geom_text(data = subset(pk, set == "Full core"), aes(label = sprintf("%.1f%%", pct)), vjust = 2.1, hjust = 0.6, size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = c("Full core" = "#0072B2", "Robust core (without GSE57614)" = "#D55E00"), name = NULL) +
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5)) + scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
    labs(x = "THP-1 effect-size cut-off, |log2FC|", y = "Primary core genes classed as conserved", title = "Conserved share vs cut-off") +
    theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0))
  mar <- seq(0.25, 2, by = 0.05)
  ar <- do.call(rbind, lapply(mar, function(m) { eqv <- !diff_ & (abs(delta) + 1.645 * se_d) < m
    data.frame(margin = m, Equivalent = 100 * mean(eqv[isc]), Different = 100 * mean(diff_[isc]), Inconclusive = 100 * mean((!eqv & !diff_)[isc])) }))
  al <- rbind(data.frame(margin = ar$margin, pct = ar$Different, cl = "Different"), data.frame(margin = ar$margin, pct = ar$Inconclusive, cl = "Inconclusive"),
              data.frame(margin = ar$margin, pct = ar$Equivalent, cl = "Equivalent"))
  al$cl <- factor(al$cl, levels = c("Different", "Inconclusive", "Equivalent"))
  at <- ar[round(ar$margin, 2) %in% c(0.5, 1, 1.5), ]
  pb <- ggplot(al, aes(margin, pct, fill = cl)) + geom_area(colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = c(0.5, 1, 1.5), linetype = "22", linewidth = 0.3) +
    geom_text(data = at, aes(x = margin, y = Equivalent + 5, label = sprintf("%.0f%%", Equivalent)), inherit.aes = FALSE, size = 2.8, colour = "black", fontface = "bold", hjust = -0.15) +
    scale_fill_manual(values = c(Different = "#E69F00", Inconclusive = "grey80", Equivalent = "#009E73"), name = NULL) +
    scale_x_continuous(breaks = c(0.5, 1, 1.5, 2), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = "Equivalence margin (log2)", y = "Primary core genes", title = "Equivalence vs margin") + theme(legend.position = "right")
  save_fig(pa + pb + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "a"), "Figure3_sensitivity", 185, 78) })

# ---- Figure 4: readout genes: a) per-dataset log2FC, b) forest plot of the primary-minus-THP-1 difference ------------------
run("Figure4", {
  pc <- grep("^lfc_", names(mp), value = TRUE); tc <- grep("^lfc_", names(mt), value = TRUE)
  mk <- function(tab, cols, grp) do.call(rbind, lapply(cols, function(cn) data.frame(gene = tab$gene, ds = sub("^lfc_", "", cn), lfc = tab[[cn]], grp = grp)))
  hm <- rbind(mk(mp[mp$gene %in% rdg, ], pc, "Primary"), mk(mt[mt$gene %in% rdg, ], tc, "THP-1"))
  hm$gene <- factor(hm$gene, levels = rev(rdg)); hm$grp <- factor(hm$grp, levels = c("Primary", "THP-1"))
  pa <- ggplot(hm, aes(ds, gene, fill = lfc)) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f", lfc), colour = abs(lfc) > 5.5), size = 2.8) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-6, 10), oob = squish, name = "log2FC") +
    facet_grid(. ~ grp, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = "Per-dataset response") + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(face = "italic"))
  fr <- cm[cm$gene %in% rdg, ]; fr$lo <- (fr$beta_prim - fr$beta_thp) - 1.645 * sqrt(fr$se_prim^2 + fr$se_thp^2)
  fr$hi <- (fr$beta_prim - fr$beta_thp) + 1.645 * sqrt(fr$se_prim^2 + fr$se_thp^2); fr$d <- fr$beta_prim - fr$beta_thp; fr$gene <- factor(fr$gene, levels = rev(rdg))
  pb <- ggplot(fr, aes(d, gene)) + annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#009E73", alpha = 0.12) +
    geom_vline(xintercept = 0, linewidth = 0.3) + geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.8, colour = "grey30") +
    geom_point(size = 2.4, colour = "#0072B2") + scale_x_continuous(breaks = seq(-2, 8, 2)) +
    labs(x = "Primary − THP-1 log2FC (90% CI)", y = NULL, title = "Difference between systems") + theme(axis.text.y = element_text(face = "italic"))
  save_fig(pa + pb + plot_layout(widths = c(1.5, 1)) + plot_annotation(tag_levels = "a"), "Figure4_readouts", 180, 65) })

# ---- Figure 5: validation: a) volcano of the primary meta-analysis, b) leave-one-dataset-out dumbbell -------------------
run("Figure5", {
  vo <- mp; vo$y <- pmin(-log10(pmax(vo$FDR, 1e-50)), 50); vo$grp <- ifelse(vo$core, "Core", "Other")
  pa <- ggplot(vo, aes(beta, y)) + geom_point(data = subset(vo, grp == "Other"), colour = "grey82", size = 0.3, stroke = 0) +
    geom_point(data = subset(vo, grp == "Core"), aes(colour = beta > 0), size = 0.5, stroke = 0) +
    geom_point(data = vo[match(rdg, vo$gene), ], shape = 21, size = 2.4, stroke = 0.5, fill = NA, colour = "black") +
    geom_label(data = transform(vo[match(rdg, vo$gene), ], lx = beta + c(0.3, -0.3, 0.3, -0.3), ly = y + c(-2.5, 0.8, 2.8, -2.8), hj = c(0, 1, 0, 1)),
               aes(x = lx, y = ly, label = gene, hjust = hj), fontface = "italic", size = 2.8, fill = alpha("white", 0.85), label.size = 0, label.padding = unit(0.5, "mm")) +
    geom_vline(xintercept = c(-1, 1), linetype = "22", linewidth = 0.3) + geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.3) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), labels = c(`TRUE` = "Up", `FALSE` = "Down"), name = "Core gene") +
    labs(x = "Primary meta-analysed log2FC", y = expression(-log[10]~FDR), title = sprintf("Primary M1 core signature (n = %s)", format(sum(mp$core), big.mark = ","))) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
  lg <- rbind(data.frame(ds = loo$held_out, v = 100 * loo$sign_concordance_core, m = "Same direction"), data.frame(ds = loo$held_out, v = loo$pct_core_sig_same_sign, m = "Same direction and P < 0.05"))
  lg$ds <- factor(lg$ds, levels = loo$held_out[order(loo$pct_core_sig_same_sign)])
  pb <- ggplot(lg, aes(v, ds)) + geom_line(aes(group = ds), colour = "grey70", linewidth = 0.8) + geom_point(aes(colour = m), size = 2.6) +
    geom_text(aes(label = sprintf("%.0f", v)), vjust = -1, size = 2.6) +
    scale_colour_manual(values = c("Same direction" = "#009E73", "Same direction and P < 0.05" = "#5B4B9A"), name = NULL) +
    scale_x_continuous(limits = c(50, 105), breaks = seq(50, 100, 10), labels = function(x) paste0(x, "%")) +
    labs(x = "Core genes replicated in the held-out dataset", y = NULL, title = "Leave-one-dataset-out validation") + theme(legend.position = "bottom")
  save_fig(pa + pb + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "a"), "Figure5_validation", 185, 80) })

# ---- Figure 6: dataset-level similarity to the primary core (lollipop) ---------------------------------------------
run("Figure6", {
  rk$grp <- ifelse(grepl("^primary", rk$type), "Primary (vs the other 3)", ifelse(grepl("LPS only", rk$type), "THP-1, LPS only", "THP-1, LPS + IFN-γ"))
  rk$dataset <- factor(rk$dataset, levels = rev(rk$dataset[order(-rk$spearman)]))
  p <- ggplot(rk, aes(spearman, dataset, colour = grp)) + geom_segment(aes(x = 0, xend = spearman, yend = dataset), linewidth = 0.8) +
    geom_point(size = 3) + geom_text(aes(label = sprintf("%.2f", spearman)), hjust = -0.5, size = 2.8, show.legend = FALSE) +
    scale_colour_manual(values = c("Primary (vs the other 3)" = "#009E73", "THP-1, LPS + IFN-γ" = "#5B4B9A", "THP-1, LPS only" = "#B8B0DC"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.08))) + labs(x = "Spearman correlation with the primary M1 core", y = NULL, title = "Which datasets reproduce the primary core?") +
    theme(legend.position = "bottom")
  save_fig(p, "Figure6_dataset_ranking", 110, 75) })

# ---- Figure 7: enrichment dot plot  ----------------------------------------------------------------
run("Figure7", {
  fl <- c(conserved = "Conserved", weaker_in_THP1 = "Weaker in THP-1", not_reproduced_or_opposite = "Not reproduced / opposite", THP1_specific = "THP-1-strong, not in primary core")
  tidy <- function(s) { s <- sub("^HALLMARK_", "", sub("^GOBP_", "", s)); s <- tolower(gsub("_", " ", s)); s <- paste0(toupper(substr(s, 1, 1)), substring(s, 2))
    s <- gsub("Tnfa", "TNF-α", s); s <- gsub("nfkb", "NF-κB", s); s <- gsub("Il2 stat5", "IL2-STAT5", s); s <- gsub("mhc class ii", "MHC class II", s)
    s <- gsub("mhc", "MHC", s); s <- gsub("Interferon gamma", "Interferon-γ", s); s <- gsub("Interferon alpha", "Interferon-α", s)
    ifelse(nchar(s) > 58, paste0(substr(s, 1, 56), "..."), s) }
  en <- do.call(rbind, lapply(names(fl), function(g) { f <- file.path(PATH$results, "meta", paste0("enrichment_", g, ".csv")); if (!file.exists(f)) return(NULL)
    d <- read.csv(f); d <- d[d$FDR < 0.05, ]; d <- head(d[order(d$p), ], 5); if (nrow(d) == 0) return(NULL); d$cl <- fl[g]; d }))
  en$lab <- tidy(en$set); en$cl <- factor(en$cl, levels = fl); en$lab <- factor(paste(en$lab, en$cl), levels = rev(paste(en$lab, en$cl)))
  p <- ggplot(en, aes(fold, lab, size = overlap, colour = -log10(FDR))) + geom_point() +
    facet_grid(cl ~ ., scales = "free_y", space = "free_y") + scale_y_discrete(labels = function(x) sub(" (Conserved|Weaker in THP-1|Not reproduced / opposite|THP-1-strong, not in primary core)$", "", x)) +
    scale_colour_viridis_c(option = "D", end = 0.9, name = expression(-log[10]~FDR)) + scale_size_continuous(range = c(1.5, 5), name = "Genes in set") + scale_x_continuous(expand = expansion(mult = c(0.12, 0.1))) +
    labs(x = "Fold enrichment", y = NULL, title = "Over-representation of Hallmark and GO:BP sets") + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5))
  save_fig(p, "Figure7_enrichment", 170, 125) })
cat("\nDone. Figures in: ", normalizePath(OUT), "\n", sep = "")
