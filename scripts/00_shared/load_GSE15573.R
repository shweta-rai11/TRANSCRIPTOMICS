#!/usr/bin/env Rscript
# =============================================================================
# load_GSE15573.R  —  SHARED LOADER for the external blood validation cohort,
# with VERIFIED sex labels.
#
# WHY THIS FILE EXISTS
#   GSE15573 is the external blood validation cohort for a SEX-STRATIFIED
#   biomarker study, and two of its 33 GEO sex labels are wrong.
#
#   Verification: XIST (X-inactive specific transcript, high in females) and
#   Y-linked genes were checked against the submitted `gender:ch1` label. On
#   GPL6102, RPS4Y1 separates the sexes cleanly and unambiguously - 24 samples
#   at 6.87-7.25 and 9 samples at 11.54-12.32, with nothing in between:
#
#       sample      GEO label   RPS4Y1   XIST    expression-based call
#       GSM389714   F           11.94    8.26    MALE
#       GSM389731   M            7.23    8.78    FEMALE
#
#   Both sit unambiguously in the opposite cluster. The pattern - one F that is
#   male and one M that is female - is the signature of a LABEL SWAP (a sample
#   mix-up, or a one-row offset in the submitted metadata table).
#
#   NOTE ON MARKER CHOICE: do NOT average a panel of Y genes on this platform.
#   On GPL6102 the F-vs-M deltas are RPS4Y1 +4.18, KDM5D +1.54, but DDX3Y +0.03,
#   UTY +0.07, USP9Y +0.01 - three dead probes. Averaging them dilutes a clean
#   signal into noise and makes the cohort look uncheckable. Use RPS4Y1.
#
# WHY IT MATTERS
#   The male stratum of this cohort is n = 9. One mislabelled sample is 11% of
#   it, and it also puts a male into the female stratum. Any sex-stratified
#   validation AUC computed on the raw GEO labels is computed on contaminated
#   strata.
#
# WHAT THIS DOES NOT DO
#   It does not modify data/raw/GSE15573_raw.rds. The raw file stays exactly as
#   downloaded from GEO; the correction is applied on load and is auditable.
#
# USAGE
#   source("scripts/00_shared/load_GSE15573.R")
#   g <- load_gse15573()          # list(eset, meta, corrections)
#   g$meta$sex                    # verified sex
#   g$meta$sex_geo                # original GEO label
#   g$meta$sex_corrected          # TRUE where they disagree
#
# ---- References -------------------------------------------------------------
#   Toro-Dominguez D, et al. Sex-related differences in gene expression.  (sex QC rationale)
#   Staedtler F, et al. Robust sex determination from expression data.
#   Davis S, Meltzer PS. GEOquery. Bioinformatics 2007;23(14):1846-1847.
# =============================================================================
suppressMessages({ library(Biobase) })

# Samples whose GEO sex label is contradicted by sex-chromosome expression.
# Verified 2026-07-27 against XIST + RPS4Y1 on GPL6102.
GSE15573_SEX_CORRECTIONS <- data.frame(
  sample     = c("GSM389714", "GSM389731"),
  sex_geo    = c("F",         "M"),
  sex_verified = c("M",       "F"),
  RPS4Y1     = c(11.94,        7.23),
  XIST       = c(8.26,         8.78),
  stringsAsFactors = FALSE)

load_gse15573 <- function(rds = "data/raw/GSE15573_raw.rds",
                          apply_correction = TRUE,
                          verify = TRUE, verbose = TRUE) {
  e <- readRDS(rds)
  if (is.list(e) && !is(e, "ExpressionSet")) e <- e[[1]]
  pd <- pData(e)

  meta <- data.frame(
    sample  = rownames(pd),
    group   = ifelse(grepl("rheumatoid", pd[["status:ch1"]], ignore.case = TRUE),
                     "RA", "HC"),
    sex_geo = toupper(substr(as.character(pd[["gender:ch1"]]), 1, 1)),
    age     = suppressWarnings(as.numeric(as.character(pd[["age:ch1"]]))),
    stringsAsFactors = FALSE)
  meta$sex <- meta$sex_geo
  meta$sex_corrected <- FALSE

  # ---- optional re-verification from the expression data itself -------------
  if (verify) {
    fd  <- fData(e)
    sym <- as.character(fd[[grep("^gene[ ._]?symbol$", colnames(fd),
                                 ignore.case = TRUE, value = TRUE)[1]]])
    ex  <- exprs(e)
    if (max(ex, na.rm = TRUE) > 100) ex <- log2(pmax(ex, 1))   # GSE15573 is LINEAR scale
    i <- which(sym == "RPS4Y1")
    if (length(i)) {
      y   <- colMeans(ex[i, , drop = FALSE], na.rm = TRUE)
      cut <- mean(range(y, na.rm = TRUE))                      # clean bimodal gap
      call <- ifelse(y[meta$sample] > cut, "M", "F")
      mism <- which(call != meta$sex_geo)
      if (verbose) {
        message(sprintf("GSE15573 sex verification (RPS4Y1, cut = %.2f): %d/%d discordant",
                        cut, length(mism), nrow(meta)))
        if (length(mism))
          for (k in mism)
            message(sprintf("   %s: GEO=%s  RPS4Y1=%.2f  ->  %s",
                            meta$sample[k], meta$sex_geo[k],
                            y[meta$sample[k]], call[k]))
      }
      meta$sex_expr <- unname(call)
    }
  }

  # ---- apply the curated correction ----------------------------------------
  if (apply_correction) {
    k <- match(GSE15573_SEX_CORRECTIONS$sample, meta$sample)
    ok <- !is.na(k)
    meta$sex[k[ok]]           <- GSE15573_SEX_CORRECTIONS$sex_verified[ok]
    meta$sex_corrected[k[ok]] <- TRUE
    if (verbose)
      message(sprintf("GSE15573: sex label corrected for %d sample(s): %s",
                      sum(ok), paste(GSE15573_SEX_CORRECTIONS$sample[ok], collapse = ", ")))
  }

  if (verbose)
    message(sprintf("GSE15573 loaded: %d samples | RA=%d HC=%d | F=%d M=%d%s",
                    nrow(meta), sum(meta$group == "RA"), sum(meta$group == "HC"),
                    sum(meta$sex == "F"), sum(meta$sex == "M"),
                    if (apply_correction) "  (sex labels VERIFIED)" else "  (RAW GEO labels)"))

  rownames(meta) <- meta$sample
  list(eset = e, meta = meta, corrections = GSE15573_SEX_CORRECTIONS)
}

# Run directly for a standalone QC report:  Rscript scripts/00_shared/load_GSE15573.R
if (sys.nframe() == 0) {
  g <- load_gse15573()
  cat("\n=== GSE15573 sex composition, before vs after correction ===\n")
  print(table(GEO = g$meta$sex_geo, verified = g$meta$sex))
  cat("\n=== per-sex, per-group counts USED FOR VALIDATION ===\n")
  print(table(sex = g$meta$sex, group = g$meta$group))
  dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
  write.csv(g$meta, "results/tables/GSE15573_verified_metadata.csv", row.names = FALSE)
  cat("\nWrote results/tables/GSE15573_verified_metadata.csv\n")
}
