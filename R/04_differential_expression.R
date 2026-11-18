# =====================================================================
# 04_differential_expression.R
# Per-dataset differential expression, stimulated (case) versus ctrl.
#
# Input : data/processed/matrices.rds (03_build_matrices.R)
# Output: results/de/de_<accession>_<case>-ctrl.csv       main analysis (pseudocount 0.1 for TPM/FPKM data)
#         results/de_pc1/de_<accession>_<case>-ctrl.csv   sensitivity analysis (pseudocount 1; TPM/FPKM data only)
#         results/de_summary.csv                          genes tested and readout statistics per dataset
#         columns of the DE tables: gene, logFC (log2), SE, P, FDR
#
# Methods
#   microarray and TPM/FPKM data: limma on log2 values, empirical Bayes with the robust option and, except for
#     E-MTAB-5917, the mean-variance trend (DATASETS$trend in 00_config.R).
#     TPM/FPKM genes are kept when log2(value + pseudocount) exceeds log2(1 + pseudocount) in at least
#     as many samples as the smallest group, or when the gene is a readout gene.
#   counts: edgeR filterByExpr (or readout gene), TMM normalisation, voom, limma with the robust option.
#   design: ~ 0 + group + covariates of the dataset (see DATASETS in 00_config.R); a covariate is dropped
#     if it has missing values, a single level, or makes the design rank deficient.
#   SE is the standard error of the log2 fold change: sqrt(s2.post) * stdev.unscaled.
# =====================================================================
source("R/00_config.R")
require_packages(bioc = c("limma", "edgeR"))
suppressPackageStartupMessages({ library(limma); library(edgeR) })

PSEUDOCOUNT_SENS <- 1
mats <- readRDS(file.path(PATH$processed, "matrices.rds"))
dir.create(file.path(PATH$results, "de"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATH$results, "de_pc1"), showWarnings = FALSE, recursive = TRUE)

make_design <- function(group, covars = list()) {
  group <- factor(group)
  covars <- Filter(function(v) !anyNA(v) && nlevels(factor(v)) > 1, covars)
  repeat {
    df <- data.frame(group = group)
    for (n in names(covars)) df[[n]] <- factor(covars[[n]])
    des <- model.matrix(~ 0 + ., df)
    colnames(des)[seq_len(nlevels(group))] <- levels(group)
    if (qr(des)$rank == ncol(des) || !length(covars)) break
    covars <- covars[-length(covars)]
  }
  list(des = des, covars = names(covars))
}

run_de <- function(expr, group, covars, kind, case, trend = TRUE, pseudocount = PSEUDOCOUNT) {
  des <- make_design(group, covars)$des
  if (kind == "counts") {
    y <- DGEList(expr)
    keep <- filterByExpr(y, design = des) | rownames(y) %in% READOUTS
    y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    fit <- lmFit(voom(y, des), des); trend <- FALSE
  } else {
    if (kind == "tpm") {
      expr <- log2(pmax(2^expr - 1, 0) + pseudocount)     # stored as log2(value + 1)
      keep <- rowSums(expr > log2(1 + pseudocount)) >= min(table(group)) | rownames(expr) %in% READOUTS
      expr <- expr[keep, , drop = FALSE]
    }
    fit <- lmFit(expr, des)
  }
  cm <- matrix(0, nrow = ncol(des), ncol = 1, dimnames = list(colnames(des), paste0(case, "-ctrl")))
  cm[case, 1] <- 1; cm["ctrl", 1] <- -1
  f2 <- eBayes(contrasts.fit(fit, cm), trend = trend, robust = TRUE)
  tt <- topTable(f2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, SE = sqrt(f2$s2.post) * f2$stdev.unscaled[, 1],
             P = tt$P.Value, FDR = tt$adj.P.Val)
}

summ <- list()
for (acc in names(mats)) {
  m <- mats[[acc]]
  runs <- list(de = PSEUDOCOUNT)
  if (m$kind == "tpm") runs$de_pc1 <- PSEUDOCOUNT_SENS
  for (r in names(runs)) {
    res <- run_de(m$expr, m$group, m$covars, m$kind, m$case, DATASETS$trend[DATASETS$accession == acc], runs[[r]])
    write.csv(res, file.path(PATH$results, r, sprintf("de_%s_%s-ctrl.csv", acc, m$case)), row.names = FALSE)
    if (r == "de") {
      ro <- res[match(READOUTS, res$gene), ]
      summ[[acc]] <- data.frame(accession = acc, model = paste(c("group", make_design(m$group, m$covars)$covars), collapse = " + "),
                                genes = nrow(res), sig_FDR05 = sum(res$FDR < 0.05, na.rm = TRUE),
                                t(setNames(round(ro$logFC, 2), paste0("logFC_", READOUTS))))
    }
  }
  message(sprintf("%-12s %6d genes tested", acc, summ[[acc]]$genes))
}
summ <- do.call(rbind, summ)
write.csv(summ, file.path(PATH$results, "de_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
