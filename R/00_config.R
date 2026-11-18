# =====================================================================
# 00_config.R
# Central settings for the whole analysis. Every other script sources this file.
#
# Project layout (all paths are relative to the project root):
#   R/          analysis scripts
#   metadata/   hand-curated inputs that are part of the repository
#   data/raw/   files downloaded from GEO / ArrayExpress (not versioned)
#   data/processed/  gene x sample matrices built from the raw files
#   results/    tables written by the analysis
#   figures/    figures written by the analysis
#
# Usage: scripts are run from the project root (the folder that contains "R/") in numerical
# order, e.g. source("R/01_download_data.R"), or all at once with source("run_all.R").
# =====================================================================

options(stringsAsFactors = FALSE, width = 120, timeout = 1200)

# ---- project root check ------------------------------------------------
if (!file.exists("R/00_config.R"))
  stop("Working directory is not the project root (the folder that contains 'R/' and 'metadata/').")

# ---- Analysis settings -------------------------------------------------
PSEUDOCOUNT <- 0.1   # offset added before log2 for FPKM/TPM data: log2(value + PSEUDOCOUNT).
                     # 0.1 is the main analysis; 1 is a sensitivity analysis (1 compresses the fold
                     # changes of lowly expressed genes such as IL6, TNF and IL1B).
THP1_MODEL  <- Sys.getenv("THP1_MODEL", "RE")   # meta-analysis model for the three THP-1 datasets: "RE" (random effects, main)
                     # or "FE" (fixed effect, sensitivity; set with Sys.setenv(THP1_MODEL = "FE") before running 06 onwards).
SEED        <- 1L    # seed for any random procedure
READOUTS    <- c("CXCL8", "IL6", "TNF", "IL1B")   # assay readout genes

# ---- Datasets ----------------------------------------------------------
# role: "primary" = primary human monocyte-derived macrophages (reference);
#       "THP1"    = THP-1 macrophages stimulated with LPS + IFN-gamma (compared with the reference);
#       "THP1_LPS"= THP-1, LPS only (directional check, not an M1 replicate).
# kind:       array  = normalised log2 microarray intensities
#             tpm    = RNA-seq TPM or FPKM (analysed as log2(value + pseudocount))
#             counts = RNA-seq raw counts (voom)
# covariates: columns of metadata/design_final.csv added to the model besides the stimulation group
#             (lineage, timepoint, donor and medium are balanced across ctrl and M1 within each dataset).
# case:       group compared with ctrl (M1 = LPS + IFN-gamma; LPS = LPS only).
# trend:      eBayes(trend = ) for array/tpm datasets (FALSE for E-MTAB-5917, whose 6-group model is
#             fitted with robust = TRUE only); NA for count datasets, which are analysed with voom.
DATASETS <- data.frame(
  accession  = c("GSE61298", "GSE57614", "GSE162698", "GSE157182",
                 "GSE159112", "GSE234829", "E-MTAB-5917", "GSE181851"),
  role       = c(rep("primary", 4), rep("THP1", 3), "THP1_LPS"),
  platform   = c("microarray", "microarray", "RNA-seq TPM", "RNA-seq counts",
                 "RNA-seq FPKM", "RNA-seq FPKM", "microarray", "RNA-seq counts"),
  kind       = c("array", "array", "tpm", "counts", "tpm", "tpm", "array", "counts"),
  covariates = c("lineage", "timepoint", "donor", "replicate", "medium", "", "", ""),
  case       = c(rep("M1", 7), "LPS"),
  trend      = c(TRUE, TRUE, TRUE, NA, TRUE, TRUE, FALSE, NA),   # eBayes(trend=): NA = not used (counts use voom)
  stringsAsFactors = FALSE
)

# ---- Paths -------------------------------------------------------------
PATH <- list(
  raw       = "data/raw",
  processed = "data/processed",
  metadata  = "metadata",
  results   = "results",
  figures   = "figures"
)
for (p in PATH[c("raw", "processed", "results", "figures")]) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# ---- Packages ----------------------------------------------------------
# Installs a package only if it is missing. Bioconductor packages need BiocManager.
require_packages <- function(cran = character(), bioc = character()) {
  miss_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_cran)) install.packages(miss_cran)
  miss_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)
  }
  invisible(TRUE)
}

source("R/utils.R")
