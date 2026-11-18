# =====================================================================
# 03_build_matrices.R
# Builds one gene-symbol x sample expression matrix per dataset from the raw files.
#
# Input : data/raw/ (01_download_data.R), metadata/design_final.csv
# Output: data/processed/matrices.rds   list with one element per dataset:
#           kind, expr, group, covars, case
#         results/matrix_summary.csv   genes, samples and groups per dataset
#
# Processing rules
#   microarray (GSE61298, GSE57614): intensities are log2-transformed when they are on the linear scale,
#     quantile-normalised across all arrays of the series, annotated with the platform gene symbol
#     (first symbol of multi-symbol probes), and reduced to the probe with the highest mean per symbol.
#   microarray (E-MTAB-5917): normalised log2 values from expression.txt; probes are mapped to
#     symbols with hgu133plus2.db and reduced to the probe with the highest mean per symbol.
#   RNA-seq TPM / FPKM (GSE162698, GSE159112, GSE234829): stored as log2(value + 1); rows sharing a symbol
#     are reduced to the row with the highest mean. GSE234829 Ensembl gene ids are mapped with org.Hs.eg.db.
#   RNA-seq counts (GSE157182, GSE181851): counts of rows sharing a symbol are summed.
# =====================================================================
source("R/00_config.R")
require_packages(cran = c("readxl"),
                 bioc = c("GEOquery", "Biobase", "limma", "AnnotationDbi", "org.Hs.eg.db", "hgu133plus2.db"))

design <- read.csv(file.path(PATH$metadata, "design_final.csv"), colClasses = "character")
design[design == ""] <- NA

raw_file <- function(acc, pattern) {
  f <- list.files(file.path(PATH$raw, acc), pattern = pattern, full.names = TRUE)
  if (!length(f)) stop(acc, ": no raw file matching '", pattern, "' in ", file.path(PATH$raw, acc))
  f[1]
}
take_columns <- function(m, cols, acc) {
  if (!all(cols %in% colnames(m))) stop(acc, ": columns not found in the raw file: ", paste(setdiff(cols, colnames(m)), collapse = ", "))
  m[, cols, drop = FALSE]
}
first_symbol <- function(x) { x <- trimws(sub("///.*$", "", as.character(x))); x[x %in% c("", "---", "NA")] <- NA; x }
to_matrix <- function(df, cols) { m <- as.matrix(as.data.frame(lapply(df[, cols, drop = FALSE], as.numeric), check.names = FALSE)); dimnames(m) <- list(NULL, cols); m }

# ---- microarray series from GEO ---------------------------------------------
build_geo_array <- function(acc, platform_file, symbol_columns) {
  ex <- Biobase::exprs(GEOquery::getGEO(filename = raw_file(acc, "series_matrix"), GSEMatrix = TRUE, getGPL = FALSE))
  qx <- quantile(ex, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
  if (qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0)) { ex[ex <= 0] <- NA; ex <- log2(ex) }
  ex <- limma::normalizeBetweenArrays(ex, method = "quantile")
  plat <- GEOquery::Table(GEOquery::getGEO(filename = raw_file(acc, platform_file)))
  sc <- intersect(symbol_columns, colnames(plat))
  if (!length(sc)) stop(acc, ": no gene-symbol column in the platform table: ", paste(colnames(plat), collapse = ", "))
  sym <- first_symbol(plat[[sc[1]]][match(rownames(ex), as.character(plat$ID))])
  collapse_by_symbol(ex, sym)
}

# ---- build all datasets -----------------------------------------------------
mats <- list()
for (i in seq_len(nrow(DATASETS))) {
  acc <- DATASETS$accession[i]; kind <- DATASETS$kind[i]
  d <- design[design$accession == acc, ]
  message("Building ", acc, " (", DATASETS$platform[i], ") ...")

  full <- switch(acc,
    "GSE61298"  = build_geo_array(acc, "GPL15207", c("Gene Symbol")),
    "GSE57614"  = build_geo_array(acc, "GPL6480", c("Gene symbol")),
    "E-MTAB-5917" = {
      x <- read.delim(raw_file(acc, "expression"), check.names = FALSE, quote = "")
      m <- to_matrix(x, setdiff(colnames(x), colnames(x)[1])); rownames(m) <- x[[1]]
      sym <- unname(suppressMessages(AnnotationDbi::mapIds(hgu133plus2.db::hgu133plus2.db, keys = rownames(m),
                                                           column = "SYMBOL", keytype = "PROBEID", multiVals = "first")))
      collapse_by_symbol(m, sym)
    },
    "GSE162698" = {
      f <- raw_file(acc, "tpm")
      nc <- max(count.fields(f, sep = "\t", quote = "", comment.char = ""))
      x <- read.delim(f, header = FALSE, quote = "", comment.char = "", colClasses = "character",
                      col.names = paste0("V", seq_len(nc)), fill = TRUE)
      nm <- as.character(unlist(x[3, -(1:2)])); x <- x[-(1:6), ]
      cols <- intersect(d$data_column, nm)
      m <- to_matrix(setNames(x[, -(1:2)], nm), cols)
      rownames(m) <- NULL
      collapse_by_symbol(log2(m + 1), first_symbol(x[[2]]))
    },
    "GSE159112" = {
      x <- read.delim(raw_file(acc, "fpkm"), check.names = FALSE, quote = "", comment.char = "")
      m <- to_matrix(x, grep("_FPKM$", colnames(x), value = TRUE))
      collapse_by_symbol(log2(m + 1), first_symbol(x$Symbol))
    },
    "GSE234829" = {
      x <- as.data.frame(readxl::read_excel(raw_file(acc, "xlsx"), sheet = "gene_fpkm"), check.names = FALSE)
      cols <- grep("^WT_M[012]_", colnames(x), value = TRUE)   # wild-type arms only; CES1 knock-down arms are not used
      m <- to_matrix(x, cols)
      collapse_by_symbol(log2(m + 1), first_symbol(map_ids(sub("\\.[0-9]+$", "", x$gene_id), "ENSEMBL")))
    },
    "GSE157182" = {
      x <- as.data.frame(readxl::read_excel(raw_file(acc, "xlsx")), check.names = FALSE)
      cols <- setdiff(colnames(x), "Gene")
      sum_by_symbol(to_matrix(x, cols), first_symbol(x$Gene))
    },
    "GSE181851" = {
      x <- read.delim(raw_file(acc, "counts"), check.names = FALSE, quote = "", comment.char = "")
      cols <- grep("^(Unstim|LPS|L_PW|PW)[0-9]$", colnames(x), value = TRUE)
      sum_by_symbol(to_matrix(x, cols), first_symbol(x$gene_name))
    })

  e <- take_columns(full, d$data_column, acc)
  colnames(e) <- d$sample_id
  covars <- lapply(setNames(nm = DATASETS$covariates[i][nzchar(DATASETS$covariates[i])]), function(v) d[[v]])
  mats[[acc]] <- list(kind = kind, expr = e, group = d$group, covars = covars, case = DATASETS$case[i])
}

saveRDS(mats, file.path(PATH$processed, "matrices.rds"))

summ <- do.call(rbind, lapply(names(mats), function(a) describe_matrix(a, mats[[a]]$kind, mats[[a]]$expr, mats[[a]]$group)))
write.csv(summ, file.path(PATH$results, "matrix_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
cat("\nReadout genes present per dataset (before expression filtering):\n")
print(sapply(mats, function(m) READOUTS %in% rownames(m$expr), USE.NAMES = TRUE) |> `rownames<-`(READOUTS))
