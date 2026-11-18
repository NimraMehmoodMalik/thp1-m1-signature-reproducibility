# =====================================================================
# run_all.R
# Runs the complete analysis in order. Start R in the project root.
# =====================================================================
steps <- c("R/01_download_data.R", "R/02_check_design.R", "R/03_build_matrices.R", "R/04_differential_expression.R",
           "R/05_primary_meta.R", "R/06_thp1_meta_classes.R", "R/07_contrast_equivalence.R",
           "R/08_sensitivity_robust_core.R", "R/09_crossfit_pathway.R", "R/10_enrichment.R", "R/11_figures.R", "R/12_verify_results.R")
for (s in steps) { message("\n===== ", s, " ====="); source(s, echo = FALSE) }
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
