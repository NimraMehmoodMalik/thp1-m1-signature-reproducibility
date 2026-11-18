# Conservation of the primary human macrophage M1 response in THP-1 cells

Code for the secondary analysis of public transcriptome data (GEO, ArrayExpress) that asks how much of the core
M1 signature (LPS + IFN-gamma versus unstimulated) of primary human monocyte-derived macrophages is reproduced in
THP-1 macrophages. The readout genes of the accompanying wet-lab assay are CXCL8, IL6, TNF and IL1B.

All input data are public and are downloaded by the first script; no raw data are stored in this repository.

## Run

Open the project folder in R (>= 4.3) and run, from the project root:

```r
source("run_all.R")
```

This runs scripts 01 to 12 in order and writes tables to `results/`, figures to `figures/` and the software
versions to `results/sessionInfo.txt`. The first run needs internet access and installs missing packages
(CRAN and Bioconductor); it takes roughly an hour, mostly for the downloads (about 360 MB) and the platform annotation.
The last script compares the main results with the values reported in the manuscript and stops if any differ
(`results/verification.csv`).

The R and package versions of the run are written to `results/sessionInfo.txt`; the reported results were obtained
with R 4.6.1 and Bioconductor 3.23.

## Scripts

| Script | Purpose |
|---|---|
| `R/00_config.R` | Settings (pseudocount, THP-1 meta-analysis model, seed, readout genes), dataset table, paths |
| `R/utils.R`, `R/meta_functions.R` | Helper functions; meta-analysis, core-gene and equivalence functions |
| `R/01_download_data.R` | Download the input files; `data/raw/download_log.csv` records file sizes and MD5 sums |
| `R/02_check_design.R` | Validate `metadata/design_final.csv` against the downloaded files |
| `R/03_build_matrices.R` | Gene-symbol expression matrix per dataset |
| `R/04_differential_expression.R` | Per-dataset limma / voom analysis, stimulated versus control |
| `R/05_primary_meta.R` | Random-effects meta-analysis of the four primary datasets, core signature, leave-one-dataset-out |
| `R/06_thp1_meta_classes.R` | THP-1 meta-analysis, conserved / weaker / not reproduced / opposite / THP-1-specific classes |
| `R/07_contrast_equivalence.R` | Direct primary-versus-THP-1 contrast and equivalence test |
| `R/08_sensitivity_robust_core.R` | Effect-size thresholds, robust core, benchmark against Chen et al. 2019, dataset ranking |
| `R/09_crossfit_pathway.R` | Cross-fitting, pathway-level equivalence, THP-1 meta-analysis model sensitivity |
| `R/10_enrichment.R` | Over-representation analysis (MSigDB Hallmark and GO:BP) |
| `R/11_figures.R` | Figures 2 to 7 |
| `R/12_verify_results.R` | Check of the headline numbers |

## Data

| Accession | Role | Platform | Contrast (n per arm) |
|---|---|---|---|
| GSE61298 | primary | microarray (Affymetrix PrimeView) | M1 vs ctrl (6 vs 6) |
| GSE57614 | primary | microarray (Agilent) | M1 vs ctrl (9 vs 9) |
| GSE162698 | primary | RNA-seq (TPM) | M1 vs ctrl (3 vs 3) |
| GSE157182 | primary | RNA-seq (counts) | M1 vs ctrl (3 vs 3) |
| GSE159112 | THP-1 | RNA-seq (FPKM) | M1 vs ctrl (6 vs 6) |
| GSE234829 | THP-1 | RNA-seq (FPKM), wild-type arms | M1 vs ctrl (4 vs 4) |
| E-MTAB-5917 | THP-1 | microarray (Affymetrix HG-U133 Plus 2) | M1 vs ctrl (3 vs 3) |
| GSE181851 | THP-1, LPS only | RNA-seq (counts) | LPS vs ctrl (3 vs 3), directional check only |

`metadata/design_final.csv` lists every sample used, its group and the covariates (lineage, timepoint, donor,
replicate, medium) of the model. In E-MTAB-5917 the control (mac) and M1 arms were harvested at different culture
times (1 and 3 days according to the SDRF), so the M1 effect in that dataset is confounded with culture time.
`metadata/chen2019_LPS_IFNg_vs_M0.csv` is the M-LPS+IFN-gamma versus M0 table of the meta-analysis by
Chen HJ et al., Front Immunol 2019;10:2887 (doi:10.3389/fimmu.2019.02887), used as an external benchmark.

## Methods in brief

- Per dataset: limma on log2 intensities (arrays) or log2(value + 0.1) (TPM / FPKM), voom for counts, empirical Bayes
  with the robust option. Covariates are included where they are balanced between arms.
- Primary meta-analysis: DerSimonian-Laird random effects with the modified Hartung-Knapp standard error;
  significance from the signed Stouffer combination of the per-dataset P values (BH). Core gene: FDR < 0.05,
  |log2FC| >= 1, same direction in all datasets and nominal P < 0.05 in at least max(2, k - 1) datasets.
- THP-1 meta-analysis: random effects (main), fixed effect as sensitivity (`THP1_MODEL` in `00_config.R`).
- Equivalence: the 90% confidence interval of the primary-minus-THP-1 difference is compared with a margin of 1 log2.

## Notes

- Gene-set results depend on the msigdbr version (26.1.1 for the reported results).
- GSE234829 symbols are mapped with org.Hs.eg.db; results can differ for a few genes with other annotation releases.
- A seed (`SEED` in `00_config.R`) is defined for completeness; no step of the analysis uses random numbers.

## Licence

MIT (see `LICENSE`). The input data keep the licences of their original repositories.
