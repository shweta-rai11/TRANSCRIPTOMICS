#!/usr/bin/env Rscript
# Shared loader for the external blood validation cohort (GSE15573), with 2 of 33 GEO sex labels corrected via RPS4Y1/XIST expression verification.
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

  # optional re-verification from the expression data itself
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

  # apply the curated correction
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
