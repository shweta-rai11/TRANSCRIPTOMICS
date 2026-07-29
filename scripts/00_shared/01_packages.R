#!/usr/bin/env Rscript
# =============================================================================
# 00_packages.R  —  ALL packages used across the sex-stratified RA biomarker
# pipeline. Installs anything missing (CRAN + Bioconductor + GitHub), then loads
# and prints versions. Run once to set up a fresh machine:  Rscript scripts/00_packages.R

# R version 4.4.2

# =============================================================================

# ---- packages grouped by pipeline stage ------------------------------------
cran <- c(
  # data wrangling & plotting
  "data.table", "dplyr", "reshape2",
  "ggplot2", "ggraph", "ggrepel", "patchwork", "scales",
  "pheatmap", "RColorBrewer", "VennDiagram", "ggVennDiagram", "corrplot",
  # networks
  "igraph",
  # WGCNA (CRAN, pulls Bioconductor deps)
  "WGCNA",
  # machine learning
  "glmnet", "randomForest", "e1071", "Boruta", "caret", "kernlab", "nnet",
  "pROC",
  # Mendelian randomization (ieugwasr is pulled in as a dependency)
  "MendelianRandomization", "ieugwasr"
)

bioc <- c(
  # GEO retrieval, expression objects, annotation
  "GEOquery", "Biobase", "AnnotationDbi", "org.Hs.eg.db",
  # normalization / batch / differential expression
  "limma", "sva", "edgeR", "DESeq2",
  # functional enrichment
  "clusterProfiler", "enrichplot", "GO.db",
  # protein-protein interaction networks
  "STRINGdb"
)

github <- c(
  TwoSampleMR = "MRCIEU/TwoSampleMR",   # two-sample MR
  IOBR        = "IOBR/IOBR"             # CIBERSORT / LM22 deconvolution
)


# ---- installer --------------------------------------------------------------
inst <- function(p) suppressWarnings(requireNamespace(p, quietly = TRUE))

message("Checking CRAN packages...")
miss <- cran[!vapply(cran, inst, logical(1))]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")

message("Checking Bioconductor packages...")
if (!inst("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
missb <- bioc[!vapply(bioc, inst, logical(1))]
if (length(missb)) BiocManager::install(missb, update = FALSE, ask = FALSE)

message("Checking GitHub packages...")
if (!inst("remotes")) install.packages("remotes", repos = "https://cloud.r-project.org")
for (pkg in names(github)) if (!inst(pkg)) remotes::install_github(github[[pkg]], upgrade = "never")

# ---- load all + report versions --------------------------------------------
all_pkgs <- c(cran, bioc, names(github))
ok <- vapply(all_pkgs, function(p)
  suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE)), logical(1))

cat("\n================ PACKAGE STATUS ================\n")
for (p in sort(all_pkgs)) {
  v <- if (ok[p]) as.character(packageVersion(p)) else "NOT INSTALLED"
  cat(sprintf("  %-24s %s\n", p, v))
}
cat(sprintf("\nR version: %s\n", getRversion()))
cat(sprintf("Installed OK: %d / %d\n", sum(ok), length(all_pkgs)))
if (any(!ok)) cat("MISSING:", paste(all_pkgs[!ok], collapse = ", "), "\n")

