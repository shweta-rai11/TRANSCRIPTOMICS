#!/usr/bin/env Rscript
# R version 4.4.2


# Loads, cleans and merges the two training datasets (GSE93272, GSE110169), splits 70:30 (train/frozen internal-validation holdout) before any cross-sample learning, then quantile-normalizes and ComBat batch-corrects the training set only

# loading libraries

library(GEOquery)
library(Biobase)
library(limma)
library(sva)
library(WGCNA)
library(data.table)
library(ggplot2)
library(caret)

# setting seed for reproducibility

options(stringsAsFactors = FALSE); 
set.seed(22); 
allowWGCNAThreads()

# setting up directories for raw, processed, figures, and tables
raw <- "data/raw"
proc <- "data/processed"
fig <- "results/figures"
tab <- "results/tables"

# creating directories 
for (d in c(raw, proc, fig, tab)) dir.create(d, showWarnings = FALSE, recursive = TRUE)


# ---- helpers ---------------------------------------------------------------
collapse_to_genes <- function(eset) {
  fd <- fData(eset)
  symbol <- grep("^gene[ ._]?symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  sym <- as.character(fd[[symbol]])
  ex  <- exprs(eset)

  # Keep only probes with usable gene symbols and remove ambiguous mappings.
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
  ex <- ex[keep, , drop = FALSE]; sym <- sym[keep]

  # Remove probes that are completely missing across all samples before collapsing.
  ok <- rowSums(is.na(ex)) < ncol(ex)
  ex <- ex[ok, , drop = FALSE]; sym <- sym[ok]

  if (nrow(ex) == 0) stop("No probes remain for gene-level collapsing.")

# collapseRows MaxMean rule (WGCNA) for probe->gene collapsing
  g <- WGCNA::collapseRows(
    ex,
    rowGroup = sym,
    rowID = rownames(ex),
    method = "MaxMean",
    connectivityBasedCollapsing = FALSE,
    connectivityPower = 1,
    selectFewestMissing = TRUE,
    thresholdCombine = NA
  )$datETcollapsed

  message(sprintf("   collapse [%s]: %d probes -> %d genes", symbol, nrow(ex), nrow(g)))
  g
}
log2_norm <- function(m, plot_dir = fig, prefix = "expr_distribution") {
  vals <- as.numeric(m)
  vals <- vals[is.finite(vals) & vals > 0]

  if (length(vals) < 20) {
    message("Not enough positive values for a distribution-based log2 check so leave the data unchanged.")
    return(m)
  }

  sample_vals <- sample(vals, min(5000, length(vals)))
  sample_df <- data.frame(value = sample_vals)

  p <- ggplot(sample_df, aes(x = value)) +
    geom_histogram(bins = 50, color = "black", fill = "steelblue", alpha = 0.7) +
    geom_density(color = "firebrick", linewidth = 0.9) +
    theme_minimal() +
    labs(title = "Expression value distribution",
         x = "Expression value",
         y = "Count")

  plot_file <- file.path(plot_dir, paste0(prefix, ".png"))
  ggsave(plot_file, p, width = 7, height = 4, dpi = 150)
  message(sprintf("Saved expression distribution plot: %s", plot_file))

  # Decide on scale (not normality): transform only if data are on a linear scale (max > 100)
  qx <- as.numeric(quantile(vals, c(0.25, 0.99), na.rm = TRUE))
  needs_log <- qx[2] > 100
  message(sprintf("log2 scale check: 99th pct = %.2f -> %s",
                  qx[2], if (needs_log) "linear scale, applying log2" else "already log-scaled, no transform"))

  if (needs_log) {
    m[m <= 0] <- NA
    m <- log2(m)
    m <- m[rowSums(is.na(m)) == 0, ]
  }

  m
}

summarize_norm_diagnostics <- function(m) {
  sample_medians <- apply(m, 2, median, na.rm = TRUE)
  sample_iqr <- apply(m, 2, IQR, na.rm = TRUE)
  data.frame(n_samples = ncol(m),
             n_genes = nrow(m),
             max_value = max(m, na.rm = TRUE),
             min_value = min(m, na.rm = TRUE),
             median_sd = sd(sample_medians, na.rm = TRUE),
             iqr_sd = sd(sample_iqr, na.rm = TRUE),
             median_range = diff(range(sample_medians, na.rm = TRUE)),
             iqr_range = diff(range(sample_iqr, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}



# ---- per-dataset processing -------------------------------------------------
harmonize_group <- function(dis) {
  dis <- as.character(dis)
  ifelse(grepl("sle|lupus", dis, ignore.case = TRUE), "SLE",
         ifelse(grepl("normal|healthy|control", dis, ignore.case = TRUE), "HC",
                ifelse(grepl("rheumatoid|(^|[^a-z])ra([^a-z]|$)", dis, ignore.case = TRUE), "RA", NA_character_)))
}

process <- function(rds, nm) {
  message("Processing ", nm, " (", rds, ")")
  e <- readRDS(file.path(raw, rds)); pd <- pData(e)
  clin <- grep(":ch1$", colnames(pd), value = TRUE)
  discol <- grep("disease|status", clin, value = TRUE, ignore.case = TRUE)[1]
  sexcol <- grep("gender|sex", clin, value = TRUE, ignore.case = TRUE)[1]
  agecol <- grep("age", clin, value = TRUE, ignore.case = TRUE)
  idcol  <- grep("individual|patient|subject", clin, value = TRUE, ignore.case = TRUE)
  batchcol <- grep("batch", clin, value = TRUE, ignore.case = TRUE)
  n <- nrow(pd)
  getnum <- function(col) if (col %in% colnames(pd)) suppressWarnings(as.numeric(as.character(pd[[col]]))) else rep(NA_real_, n)
  tx <- pmin(getnum("ifx_days:ch1"), getnum("mtx_days:ch1"), getnum("tcz_days:ch1"), na.rm = TRUE)
  tx[is.infinite(tx)] <- NA
  meta <- data.frame(sample = rownames(pd), dataset = nm,
                     group = harmonize_group(as.character(pd[[discol]])),
                     sex = { s <- toupper(substr(gsub("[^A-Za-z]", "", as.character(pd[[sexcol]])), 1, 1)); s[!s %in% c("F", "M")] <- NA; s },
                     age = if (length(agecol)) getnum(agecol[1]) else rep(NA_real_, n),
                     id = if (length(idcol)) as.character(pd[[idcol]]) else rownames(pd),
                     tx = tx,                                                        # min treatment-exposure days (for baseline dedup)
                     rin = getnum("rin:ch1"),                                        # RNA integrity (dedup tie-break)
                     batch = if (length(batchcol)) as.character(pd[[batchcol[1]]]) else NA_character_,
                     stringsAsFactors = FALSE)
  ex <- log2_norm(collapse_to_genes(e))
  list(ex = ex, meta = meta)
}

d93  <- process("GSE93272_raw.rds",  "GSE93272")
d110 <- process("GSE110169_raw.rds", "GSE110169")

# ---- keep RA and HC, drop SLE / unknown sex, sex-gate, dedup GSE93272 -----------
clean <- function(d, dedup = FALSE) {
  m <- d$meta
  ex <- d$ex
  keep <- m$group %in% c("RA", "HC") & !is.na(m$group) & !is.na(m$sex)   # RA+HC only, known sex and known group
  m <- m[keep, ]; ex <- ex[, m$sample, drop = FALSE]
  if (dedup) {                                  # baseline: min tx, tie-break max RIN
    ord <- order(m$id, ifelse(is.na(m$tx),Inf,m$tx), -ifelse(is.na(m$rin),-Inf,m$rin))
    m <- m[ord, ]; m <- m[!duplicated(m$id), ]; ex <- ex[, m$sample, drop = FALSE]
  }
  message(sprintf("   %s cleaned: %d samples (RA=%d, HC=%d; M=%d, F=%d)", m$dataset[1], nrow(m),
                  sum(m$group=="RA"), sum(m$group=="HC"), sum(m$sex=="M"), sum(m$sex=="F")))
  list(ex = ex, meta = m)
}
c93  <- clean(d93,  dedup = TRUE)
c110 <- clean(d110, dedup = FALSE)

# ---- intersect common genes + merge ----------------------------------------
common <- intersect(rownames(c93$ex), rownames(c110$ex))
message(sprintf("\nCommon genes: %d (GSE93272=%d, GSE110169=%d)",
                length(common), nrow(c93$ex), nrow(c110$ex)))
merged <- cbind(c93$ex[common, ], c110$ex[common, ])
meta <- rbind(c93$meta, c110$meta)
meta$batch_full <- paste(meta$dataset, meta$batch, sep = "_")
stopifnot(identical(colnames(merged), meta$sample))

# ---- 70:30 TRAIN / INTERNAL-VALIDATION SPLIT (leakage-safe) -----------------
# Split on merged pre-normalization data, before quantile normalization and ComBat learn cross-sample parameters, to avoid leaking the holdout into training
set.seed(70)                                          # reproducible 70:30 partition
split_frac <- 0.70                                    # 70:30 train:test ratio

# Stratified split via caret::createDataPartition (strata = dataset x group x sex)
strata    <- factor(with(meta, paste(dataset, group, sex, sep = "|")))
train_idx <- caret::createDataPartition(strata, p = split_frac, list = FALSE)
train_sel <- logical(nrow(meta))
train_sel[train_idx] <- TRUE

# Freeze the internal-validation holdout UNTOUCHED (pre-norm, pre-ComBat).
meta_test   <- meta[!train_sel, , drop = FALSE]
merged_test <- merged[, meta_test$sample, drop = FALSE]
saveRDS(list(expr_prenorm = as.matrix(merged_test),
             meta  = meta_test,
             genes = common,
             role  = "internal_validation_holdout_30pct",
             note  = paste("Held out BEFORE quantile normalisation and ComBat.",
                           "At test time, place on the TRAINING scale using parameters FROZEN",
                           "from training (e.g. quantile-normalise against the training reference;",
                           "reference-batch ComBat) - never re-estimate jointly with training.")),
        file.path(proc, "internal_val_holdout.rds"))
fwrite(meta_test, file.path(tab, "internal_val_holdout_meta.csv"))

# Restrict ALL downstream analysis to the 70% training set.
merged <- merged[, train_sel, drop = FALSE]
meta   <- meta[train_sel, , drop = FALSE]
stopifnot(identical(colnames(merged), meta$sample))
cat("\n===== 70:30 SPLIT (stratified by dataset x group x sex) =====\n")
cat(sprintf("Training set (analysis): %d samples (%.1f%%)\n", ncol(merged),  100 * ncol(merged)  / (ncol(merged) + nrow(meta_test))))
cat(sprintf("Internal-validation holdout (frozen): %d samples (%.1f%%) -> data/processed/internal_val_holdout.rds\n",
            nrow(meta_test), 100 * nrow(meta_test) / (ncol(merged) + nrow(meta_test))))
cat("Training composition:\n");            print(table(meta$dataset, meta$group, meta$sex))
cat("Holdout composition (untouched):\n"); print(table(meta_test$dataset, meta_test$group, meta_test$sex))

# ---- NORMALIZATION CHECK (before) ------------------------------------------
med_sd <- function(m, g) tapply(apply(m,2,median), g, sd)
cat("\n===== NORMALIZATION CHECK (merged) =====\n")
cat("per-sample median by dataset (mean):\n")
print(round(tapply(apply(merged,2,median), meta$dataset, mean), 3))
cat("SD of per-sample medians by dataset:\n"); print(round(med_sd(merged, meta$dataset), 3))
cat("overall SD of per-sample medians:", round(sd(apply(merged,2,median)),3), "\n")

before_diag <- summarize_norm_diagnostics(merged)
needs_qnorm <- before_diag$max_value > 100 || before_diag$median_sd > 0.5 || before_diag$iqr_sd > 0.5

if (needs_qnorm) {
  message("Quantile normalization will be applied because the training microarray data still show strong between-sample distribution differences.")
  # Quantile normalisation fitted on the TRAINING set only.
  qn_data <- limma::normalizeBetweenArrays(as.matrix(merged), method = "quantile")
  after_diag <- summarize_norm_diagnostics(qn_data)
  normalization_applied <- TRUE
} else {
  message("The merged data already look reasonably comparable, so quantile normalization is being skipped.")
  qn_data <- as.matrix(merged)
  after_diag <- before_diag
  normalization_applied <- FALSE
}

norm_diag <- rbind(cbind(stage = "before_qnorm", before_diag),
                   cbind(stage = "after_qnorm", after_diag))
fwrite(norm_diag, file.path(tab, "normalization_diagnostics.csv"))
cat("\n===== NORMALIZATION DECISION =====\n")
cat("Applied quantile normalization:", normalization_applied, "\n")
cat("Before -> median SD:", round(before_diag$median_sd, 4), "| IQR SD:", round(before_diag$iqr_sd, 4), "\n")
cat("After  -> median SD:", round(after_diag$median_sd, 4), "| IQR SD:", round(after_diag$iqr_sd, 4), "\n")

# ---- BATCH CORRECTION: ComBat (TRAINING set only) --------------------------
# ComBat estimated on the training set only; mod = ~group + sex protects biological signal from removal as batch
meta$group <- factor(meta$group, levels = c("HC","RA"))
meta$sex   <- factor(meta$sex,   levels = c("F","M"))
mod <- model.matrix(~ group + sex, data = meta)
cat("\nbatch_full x dataset:\n"); print(table(meta$batch_full, meta$dataset))
run_combat <- function(batch) ComBat(dat = qn_data, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
combat <- tryCatch(run_combat(meta$batch_full), error = function(e) {
  message("   combined batch failed (", conditionMessage(e), ") -> falling back to study-level batch")
  run_combat(meta$dataset) })
cat("ComBat done. batches used:", length(unique(meta$batch_full)), "(study+internal)\n")

# ---- SAVE -------------------------------------------------------------------
saveRDS(list(expr = combat,                       # final: normalized + batch-corrected
             expr_qnorm = qn_data,                # after quantile norm, before ComBat
             expr_prenorm = as.matrix(merged),    # merged, before normalization
             meta = meta, genes = common, role = "combined_training_70pct",
             split = "stratified 70:30 (dataset x group x sex); 30% held out untouched in internal_val_holdout.rds",
             normalization = "quantile (limma normalizeBetweenArrays) on merged common genes",
             normalization_applied = normalization_applied,
             normalization_diagnostics = norm_diag,
             batch_method = "ComBat (sva), batch=study+internal, protect group+sex"),
        file.path(proc, "combined_train.rds"))

# --- also save the BATCH-CORRECTED expression as standalone files -----------
saveRDS(combat, file.path(proc, "combined_expr_batchcorrected.rds"))
fwrite(data.table(gene = rownames(combat), as.data.table(combat)),
       file.path(proc, "combined_expr_batchcorrected.csv.gz"))
fwrite(meta, file.path(proc, "combined_meta.csv"))
cat("Saved batch-corrected expression: combined_expr_batchcorrected.rds / .csv.gz + combined_meta.csv\n")

cohort <- as.data.frame(with(meta, table(dataset, group, sex)))
fwrite(cohort, file.path(tab, "combined_cohort_summary.csv"))
cat("\n===== COMBINED COHORT =====\n"); print(table(meta$dataset, meta$group, meta$sex))
cat(sprintf("\nCombined: %d samples x %d genes | males=%d, females=%d\n",
            ncol(combat), nrow(combat), sum(meta$sex=="M"), sum(meta$sex=="F")))

cat(sprintf("\nNOTE: analysis used the 70%% TRAINING split only (%d samples). The 30%% internal-validation\n", ncol(combat)))
cat("holdout is frozen UNTOUCHED in data/processed/internal_val_holdout.rds (pre-norm, pre-ComBat).\n")
cat("Apply the frozen training normalisation/ComBat parameters to it at test time - do NOT re-fit jointly.\n")
cat("\nDONE (analysis). Figures are generated separately by 02_normalize_batch_figure.R\n")
