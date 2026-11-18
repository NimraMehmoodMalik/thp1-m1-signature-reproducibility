# =====================================================================
# 02_check_design.R
# Validates the curated sample design (metadata/design_final.csv) against the downloaded input files.
#
# Input : metadata/design_final.csv, data/raw/ (output of 01_download_data.R)
# Output: results/design_summary.csv (samples per dataset and group)
#
# design_final.csv has one row per sample used in the analysis:
#   accession   dataset identifier (GEO series or ArrayExpress experiment)
#   role        primary | THP1 | THP1_LPS (see 00_config.R)
#   sample_id   GEO sample accession (GSM...) or ArrayExpress assay name (NUID-...)
#   sample_name sample title as given by the submitters
#   group       ctrl (unstimulated) | M1 (LPS + IFN-gamma) | M2 (IL-4) | LPS (LPS only, GSE181851)
#               E-MTAB-5917 also lists its other arms (IFNg, IL4, IL4_IL13, LPS_only; monocytes are not used); they enter the
#               model of that dataset (pooled variance) but not the M1 versus ctrl contrast.
#   lineage     GSE61298 only: GM-CSF or M-CSF differentiation
#   timepoint   GSE57614 only: 6h, 12h or 24h after stimulation
#   donor       GSE162698 only: donor number (the same donors were used in all arms)
#   medium      GSE159112 only: LN (normal) or LH (lipid-enriched) culture medium
#   replicate   replicate label taken from the sample title (A/B/C in GSE157182, used as a blocking covariate there)
#   data_column name of the sample's column in the raw expression file (GSM id for the GEO series matrices)
# The M1 versus ctrl contrast is the analysis contrast; M2 samples are listed for completeness.
# =====================================================================
source("R/00_config.R")

design <- read.csv(file.path(PATH$metadata, "design_final.csv"), colClasses = "character")

# ---- 1) structure --------------------------------------------------------
req <- c("accession", "role", "sample_id", "sample_name", "group", "lineage", "timepoint", "donor", "medium", "replicate", "data_column")
stopifnot(all(req %in% names(design)))
stopifnot(setequal(unique(design$accession), DATASETS$accession))
stopifnot(!anyDuplicated(paste(design$accession, design$sample_id)))
stopifnot(all(design$group %in% c("ctrl", "M1", "M2", "LPS", "IFNg", "IL4", "IL4_IL13", "LPS_only")))
chk_role <- merge(unique(design[, c("accession", "role")]), DATASETS[, c("accession", "role")], by = "accession")
stopifnot(all(chk_role$role.x == chk_role$role.y))

# ---- 2) expected number of samples per dataset and group -------------------
expected <- data.frame(
  accession = c("GSE61298", "GSE57614", "GSE162698", "GSE157182", "GSE159112", "GSE234829", "E-MTAB-5917", "GSE181851"),
  n_ctrl    = c(6, 9, 3, 3, 6, 4, 3, 3),
  n_case    = c(6, 9, 3, 3, 6, 4, 3, 3))   # M1 samples; LPS samples for GSE181851
obs <- do.call(rbind, lapply(split(design, design$accession), function(x)
  data.frame(accession = x$accession[1], n_ctrl = sum(x$group == "ctrl"), n_case = sum(x$group %in% c("M1", "LPS")))))
cmp <- merge(expected, obs, by = "accession", suffixes = c("_expected", "_found"))
if (any(cmp$n_ctrl_expected != cmp$n_ctrl_found | cmp$n_case_expected != cmp$n_case_found)) {
  print(cmp); stop("Sample counts in design_final.csv differ from the expected design.")
}

# ---- 3) sample identifiers present in the downloaded files -----------------
# Microarray series: GSM accessions are listed in the series matrix header.
for (gse in c("GSE61298", "GSE57614")) {
  f <- list.files(file.path(PATH$raw, gse), pattern = "series_matrix", full.names = TRUE)[1]
  hdr <- grep("^!Sample_geo_accession", readLines(gzfile(f)), value = TRUE)
  gsm <- gsub("\"", "", strsplit(hdr, "\t")[[1]][-1])
  ids <- design$sample_id[design$accession == gse]
  if (!all(ids %in% gsm)) stop(gse, ": sample ids not found in the series matrix: ", paste(setdiff(ids, gsm), collapse = ", "))
}
# E-MTAB-5917: one expression column per assay ("<assay name> normalised").
expr_cols <- strsplit(readLines(file.path(PATH$raw, "E-MTAB-5917", "expression.txt"), n = 1), "\t")[[1]]
ids <- design$sample_id[design$accession == "E-MTAB-5917"]
if (!all(paste(ids, "normalised") %in% expr_cols)) stop("E-MTAB-5917: assay names not found in expression.txt")

# ---- 4) summary ------------------------------------------------------------
summ <- as.data.frame.matrix(table(design$accession, design$group))
summ <- cbind(accession = rownames(summ), role = DATASETS$role[match(rownames(summ), DATASETS$accession)], summ)
write.csv(summ, file.path(PATH$results, "design_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
cat("\nDesign check passed: ", nrow(design), " samples in ", length(unique(design$accession)), " datasets.\n", sep = "")
