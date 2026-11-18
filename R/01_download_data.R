# =====================================================================
# 01_download_data.R
# Downloads the public input files for the eight datasets into data/raw/ and records a download log
# (file name, size, MD5 checksum) in data/raw/download_log.csv.
#
# Input : none (needs internet access)
# Output: data/raw/<accession>/...  and  data/raw/download_log.csv
#
# Notes
#  * GSE61298 and GSE57614 (microarray): series matrix and platform annotation via GEOquery.
#  * RNA-seq datasets: processed supplementary files from the GEO FTP site.
#  * E-MTAB-5917: the sample table (SDRF) and the processed (normalised) expression table
#    "expression.txt" from ArrayExpress. Both files can alternatively be placed manually in
#    data/raw/E-MTAB-5917/.
#  * Files that are already present are not downloaded again.
# =====================================================================
source("R/00_config.R")
require_packages(bioc = "GEOquery")

geo_suppl <- function(gse, file) {
  stub <- sub("[0-9]{1,3}$", "nnn", gse)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", stub, gse, file)
}

# ---- 1) microarray series (series matrix + annotation) ---------------
for (gse in c("GSE61298", "GSE57614")) {
  message("Microarray series: ", gse)
  dir.create(file.path(PATH$raw, gse), showWarnings = FALSE)
  gd <- file.path(PATH$raw, gse)
  tryCatch(GEOquery::getGEO(gse, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = gd),
           error = function(e) GEOquery::getGEO(gse, GSEMatrix = TRUE, destdir = gd))   # fall back to the submitter annotation
}

# ---- 2) RNA-seq supplementary files ------------------------------------
suppl <- data.frame(
  accession = c("GSE162698", "GSE159112", "GSE234829", "GSE157182", "GSE181851"),
  file      = c("GSE162698_tpm_summary_invitro_Macs.txt.gz",
                "GSE159112_gene_fpkm.txt.gz",
                "GSE234829_Processed_data_FPKM_06_06_23_O.A.xlsx",
                "GSE157182_rawcounts_Bertrams.xlsx",
                "GSE181851_gene_counts.txt.gz"),
  stringsAsFactors = FALSE)

for (i in seq_len(nrow(suppl))) {
  dest <- file.path(PATH$raw, suppl$accession[i], suppl$file[i])
  message("Supplementary file: ", suppl$file[i])
  if (!download_if_missing(geo_suppl(suppl$accession[i], suppl$file[i]), dest))
    warning("Could not download ", suppl$file[i], call. = FALSE)
}

# ---- 3) E-MTAB-5917 (ArrayExpress) -------------------------------------
am_dir  <- file.path(PATH$raw, "E-MTAB-5917")
am_base <- "https://www.ebi.ac.uk/biostudies/files/E-MTAB-5917/"
am_files <- c("E-MTAB-5917.sdrf.txt", "expression.txt")   # sample table, normalised expression table
for (f in am_files) {
  message("ArrayExpress file: ", f)
  download_if_missing(paste0(am_base, f), file.path(am_dir, f))
}

# ---- 4) download log and completeness check ----------------------------
files <- list.files(PATH$raw, recursive = TRUE, full.names = TRUE)
files <- files[basename(files) != "download_log.csv"]
log_tab <- data.frame(file = sub(paste0("^", PATH$raw, "/"), "", files),
                      bytes = file.size(files),
                      md5 = unname(tools::md5sum(files)))
write.csv(log_tab, file.path(PATH$raw, "download_log.csv"), row.names = FALSE)

expected <- c(file.path("GSE162698", suppl$file[1]), file.path("GSE159112", suppl$file[2]),
              file.path("GSE234829", suppl$file[3]), file.path("GSE157182", suppl$file[4]),
              file.path("GSE181851", suppl$file[5]),
              file.path("E-MTAB-5917", c("E-MTAB-5917.sdrf.txt", "expression.txt")))
missing <- expected[!file.exists(file.path(PATH$raw, expected))]
has_matrix <- vapply(c("GSE61298", "GSE57614"),
                     function(g) length(list.files(file.path(PATH$raw, g), pattern = "series_matrix")) > 0, logical(1))
if (length(missing) || !all(has_matrix)) {
  stop("Missing input files:\n  ", paste(c(missing, names(has_matrix)[!has_matrix]), collapse = "\n  "),
       "\nE-MTAB-5917 files are available at https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5917 ",
       "(SDRF file and 'expression.txt') and belong in ", am_dir, ".")
}
cat("\nAll input files present (", nrow(log_tab), " files). Log: ", file.path(PATH$raw, "download_log.csv"), "\n", sep = "")
