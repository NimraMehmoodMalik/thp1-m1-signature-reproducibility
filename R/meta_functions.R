# =====================================================================
# meta_functions.R
# Functions shared by the meta-analysis scripts 05 to 07 (sourced by them).
#
#   load_de()          reads the per-dataset DE tables written by 04_differential_expression.R
#   meta_one()         meta-analysis of one gene across datasets (random effects or fixed effect)
#   run_meta()         meta-analysis of all genes for a set of datasets
#   is_core()          core-gene rule for the primary macrophage signature
#   equivalence_class() difference / equivalence classification of primary versus THP-1 effects
# =====================================================================

# Named list of DE tables, keyed "<accession>|<case>-ctrl". CXCL8 is listed as IL8 in older annotations.
load_de <- function(dir = file.path(PATH$results, "de")) {
  de <- list()
  for (i in seq_len(nrow(DATASETS))) {
    acc <- DATASETS$accession[i]; cs <- DATASETS$case[i]
    t <- read.csv(file.path(dir, sprintf("de_%s_%s-ctrl.csv", acc, cs)), colClasses = c(gene = "character"))
    if ("IL8" %in% t$gene && !"CXCL8" %in% t$gene) t$gene[t$gene == "IL8"] <- "CXCL8"
    de[[paste0(acc, "|", cs, "-ctrl")]] <- t
  }
  de
}
de_key <- function(acc) paste0(acc, "|M1-ctrl")

# Meta-analysis of one gene. y: log2 fold changes, s: their standard errors.
#   model "RE": DerSimonian-Laird tau2 with the modified Hartung-Knapp standard error (t distribution, k - 1 df)
#   model "FE": inverse-variance fixed effect
meta_one <- function(y, s, model) {
  ok <- is.finite(y) & is.finite(s) & s > 0; y <- y[ok]; s <- s[ok]; k <- length(y)
  out <- c(beta = NA, se = NA, p = NA, tau2 = NA, I2 = NA, k = k, npos = sum(y > 0), nneg = sum(y < 0))
  if (k < 2) return(out)
  w <- 1 / s^2; mu <- sum(w * y) / sum(w); Q <- sum(w * (y - mu)^2)
  if (model == "FE") {
    out[c("beta", "se")] <- c(mu, sqrt(1 / sum(w))); out["p"] <- 2 * pnorm(-abs(mu / sqrt(1 / sum(w))))
    out["tau2"] <- 0
  } else {
    C <- sum(w) - sum(w^2) / sum(w); tau2 <- max(0, (Q - (k - 1)) / C)
    ws <- 1 / (s^2 + tau2); b <- sum(ws * y) / sum(ws)
    qhk <- sum(ws * (y - b)^2) / (k - 1)
    se <- sqrt(1 / sum(ws)) * sqrt(max(1, qhk))
    out[c("beta", "se")] <- c(b, se); out["p"] <- 2 * pt(-abs(b / se), df = k - 1); out["tau2"] <- tau2
  }
  out["I2"] <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  out
}

# Meta-analysis of all genes present in at least min_k of the datasets 'dsets'.
# Significance is the signed Stouffer combination of the per-dataset limma P-values (BH-adjusted);
# the effect size beta is the meta-analytic estimate. With k = 3-4 datasets the Hartung-Knapp test
# has only 2-3 degrees of freedom and cannot reach FDR < 0.05 across ~13,000 tests, so it is
# reported (p, FDR_hk) but not used for selection.
run_meta <- function(de, dsets, min_k, model) {
  genes <- Reduce(union, lapply(dsets, function(d) de[[de_key(d)]]$gene))
  pick <- function(col) sapply(dsets, function(d) { t <- de[[de_key(d)]]; t[[col]][match(genes, t$gene)] })
  Y <- pick("logFC"); S <- pick("SE"); P <- pick("P")
  rownames(Y) <- rownames(S) <- rownames(P) <- genes
  keep <- rowSums(is.finite(Y) & is.finite(S)) >= min_k
  Y <- Y[keep, , drop = FALSE]; S <- S[keep, , drop = FALSE]; P <- P[keep, , drop = FALSE]
  res <- t(vapply(seq_len(nrow(Y)), function(i) meta_one(Y[i, ], S[i, ], model), numeric(8)))
  res <- data.frame(gene = rownames(Y), res)
  res$FDR_hk <- p.adjust(res$p, "BH")
  Z <- sign(Y) * qnorm(pmax(P, 1e-300) / 2, lower.tail = FALSE)
  nn <- rowSums(is.finite(Z)); zs <- rowSums(Z, na.rm = TRUE) / sqrt(nn)
  res$p_stouffer <- 2 * pnorm(-abs(zs)); res$FDR <- p.adjust(res$p_stouffer, "BH")
  res$n_sig <- rowSums(P < 0.05, na.rm = TRUE)
  res$concordant <- res$npos == res$k | res$nneg == res$k
  cbind(res, setNames(as.data.frame(Y), paste0("lfc_", dsets)))
}

# Core gene of the primary signature: Stouffer FDR < 0.05, |RE log2FC| >= thr, same direction in every dataset,
# nominal P < 0.05 in at least max(2, k - 1) datasets.
is_core <- function(m, thr = 1)
  !is.na(m$FDR) & m$FDR < 0.05 & abs(m$beta) >= thr & m$concordant & m$n_sig >= pmax(2, m$k - 1)

# Difference / equivalence class of a gene between the primary and the THP-1 meta-analysis (table d with
# beta_prim, se_prim, beta_thp, se_thp). 'different': BH-FDR < 0.05 and |difference| >= 1 log2;
# 'equivalent': the 90% CI of the difference lies inside +/- margin (two one-sided tests at 5%);
# otherwise 'inconclusive'. f inflates the THP-1 standard error (sensitivity analysis).
equivalence_class <- function(d, margin, f = 1) {
  se_d <- sqrt(d$se_prim^2 + (f * d$se_thp)^2); delta <- d$beta_prim - d$beta_thp
  fdr <- p.adjust(2 * pnorm(-abs(delta / se_d)), "BH")
  dsig <- fdr < 0.05 & abs(delta) >= 1
  ifelse(dsig, "different", ifelse((abs(delta) + 1.645 * se_d) < margin, "equivalent", "inconclusive"))
}

# Random-effects (DerSimonian-Laird) estimate and Stouffer FDR for genes present in ALL the given datasets,
# without the Hartung-Knapp floor. Used for the robustness analyses (core without one dataset, cross-fitting).
re_meta_common <- function(de, dsets) {
  genes <- Reduce(intersect, lapply(dsets, function(d) de[[de_key(d)]]$gene))
  pick <- function(col) sapply(dsets, function(d) { t <- de[[de_key(d)]]; t[[col]][match(genes, t$gene)] })
  Y <- pick("logFC"); S <- pick("SE"); P <- pick("P"); k <- length(dsets)
  w <- 1 / S^2; mu <- rowSums(w * Y) / rowSums(w); Q <- rowSums(w * (Y - mu)^2)
  C <- rowSums(w) - rowSums(w^2) / rowSums(w); tau2 <- pmax(0, (Q - (k - 1)) / C)
  ws <- 1 / (S^2 + tau2); beta <- rowSums(ws * Y) / rowSums(ws)
  Z <- sign(Y) * qnorm(pmax(P, 1e-300) / 2, lower.tail = FALSE); pS <- 2 * pnorm(-abs(rowSums(Z) / sqrt(k)))
  fdr <- p.adjust(pS, "BH")
  core <- fdr < 0.05 & abs(beta) >= 1 & (rowSums(Y > 0) == k | rowSums(Y < 0) == k) & rowSums(P < 0.05) >= max(2, k - 1)
  data.frame(gene = genes, beta = beta, FDR = fdr, core = core)
}

# MSigDB gene sets (Hallmark and GO biological process) as a named list of gene symbols.
# Gene-set contents depend on the msigdbr version; the version is printed by the calling scripts.
msigdb_sets <- function() {
  require_packages(cran = "msigdbr")
  get_sets <- function(coll, sub = NULL) {
    r <- try(if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", collection = coll)
             else msigdbr::msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub), silent = TRUE)
    if (inherits(r, "try-error"))
      r <- if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", category = coll)
           else msigdbr::msigdbr(species = "Homo sapiens", category = coll, subcategory = sub)
    split(r$gene_symbol, r$gs_name)
  }
  c(get_sets("H"), get_sets("C5", "GO:BP"))
}
