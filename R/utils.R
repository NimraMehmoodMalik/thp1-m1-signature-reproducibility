# =====================================================================
# utils.R
# Helper functions shared by the analysis scripts (sourced by 00_config.R).
# =====================================================================

# Download one file to 'dest' unless it is already there. Returns TRUE when the file exists afterwards.
download_if_missing <- function(url, dest) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(dest) && file.size(dest) > 0) return(TRUE)
  ok <- tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) { message("  download failed: ", conditionMessage(e)); FALSE },
                 warning = function(w) { message("  download warning: ", conditionMessage(w)); FALSE })
  if (!ok && file.exists(dest)) file.remove(dest)
  ok && file.exists(dest)
}

# Map identifiers to HGNC gene symbols with org.Hs.eg.db (keytype "ENTREZID" or "ENSEMBL").
map_ids <- function(ids, keytype) {
  ids <- as.character(ids)
  r <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids), column = "SYMBOL",
                                              keytype = keytype, multiVals = "first"))
  unname(r[ids])
}

# Collapse rows that share a gene symbol by keeping the row with the highest mean value.
# Used for log-scale data (microarray probes, log2 FPKM/TPM).
collapse_by_symbol <- function(m, sym) {
  keep <- !is.na(sym) & nzchar(sym); m <- m[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowMeans(m, na.rm = TRUE), decreasing = TRUE); m <- m[o, , drop = FALSE]; sym <- sym[o]
  dup <- duplicated(sym); m <- m[!dup, , drop = FALSE]; rownames(m) <- sym[!dup]; m
}

# Collapse raw counts by summing rows that share a gene symbol.
sum_by_symbol <- function(m, sym) {
  keep <- !is.na(sym) & nzchar(sym)
  rowsum(m[keep, , drop = FALSE], group = sym[keep])
}

# One-line summary of a built matrix, used for the build report.
describe_matrix <- function(accession, kind, m, group) {
  data.frame(accession = accession, kind = kind, genes = nrow(m), samples = ncol(m),
             groups = paste(names(table(group)), table(group), collapse = "; "))
}
