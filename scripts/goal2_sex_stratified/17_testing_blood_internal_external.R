#!/usr/bin/env Rscript
# =============================================================================
# 18b_mr_roc_validation.R
# -----------------------------------------------------------------------------
# Diagnostic ROC / AUC for the sex-stratified MR consensus biomarker panels
# (produced by 17b_mr_feature_selection.R) across THREE datasets:
#     1. TRAIN         : 70% discovery cohort (GSE93272 + GSE110169), the data
#                        the panel and model were built on.
#     2. INTERNAL TEST : the sealed 30% internal-validation holdout, projected
#                        onto the training scale with FROZEN quantile + ComBat
#                        parameters (04_apply_holdout.R) so no holdout sample
#                        influenced any training parameter (leakage-free).
#     3. EXTERNAL BLOOD: an independent RA PBMC cohort (GSE15573, Illumina),
#                        a different platform -> a true external test.
#   (Cross-tissue SYNOVIUM validation, GSE89408, is deferred: that raw data is
#    not present in data/raw. The dataset loader below is written so synovium
#    can be slotted in later without changing the analysis logic.)
#
# PANELS (three-method LASSO n RF n SVM-RFE consensus; see mr_fs_consensus_*.csv)
#     Female : BNIP2, NMI
#     Male   : CLSTN1, GABBR1, HLA-DMA, SSRP1
#   Females and males are analysed SEPARATELY throughout; the female panel is
#   only ever evaluated in females and the male panel only in males.
#
# TWO COMPLEMENTARY READOUTS PER DATASET (both requested)
#   (a) COMBINED PANEL MODEL: a multivariable logistic-regression classifier
#       (RA vs HC) fitted on the TRAINING data using the consensus genes. To
#       make the panel transferable across platforms, each gene is standardised
#       (z-score) WITHIN each dataset (and, for the training cross-validation,
#       within each CV training fold) before entering the model; this removes
#       platform-specific location/scale offsets and is unsupervised (it never
#       uses the RA/HC labels). AUC is reported as:
#         - TRAIN apparent (resubstitution) AND TRAIN 10-fold CV (honest, the
#           number to quote for the training cohort),
#         - INTERNAL TEST and EXTERNAL BLOOD from the LOCKED training model
#           applied once (honest out-of-sample).
#   (b) PER-GENE (univariate) ROC: each panel gene alone as a classifier. Its
#       orientation (which direction of expression marks RA) is FIXED on the
#       training data and that same orientation is applied to the internal and
#       external sets, so the per-gene external AUC reflects whether the
#       training-learned rule holds out of sample (a gene that reverses
#       direction is flagged as non-concordant, not silently re-oriented).
#
# All AUCs are reported with 95% DeLong confidence intervals. NOTE on power:
#   male internal test (6 RA / 7 HC) and male blood (small n) give very wide
#   CIs; these AUCs are indicative only and are flagged in the output.
#
# REPRODUCIBILITY: single seed (1234) for CV fold assignment. Package versions
#   printed at the end. Outputs:
#     results/tables/mr_roc_panel_auc.csv     (combined-panel AUC per dataset)
#     results/tables/mr_roc_pergene_auc.csv   (per-gene AUC per dataset)
#     data/processed/mr_roc_objects.rds        (ROC coordinates for the figures)
# =============================================================================

suppressMessages({library(Biobase); library(pROC); library(data.table)})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(procN, showWarnings=FALSE, recursive=TRUE); dir.create(tab, showWarnings=FALSE, recursive=TRUE)

ml <- readRDS(file.path(procN, "ml_features.rds"))
panel <- list(F = ml$female$consensus, M = ml$male$consensus)
cat(sprintf("Female panel (%d): %s\nMale panel (%d): %s\n",
            length(panel$F), paste(panel$F, collapse = ", "),
            length(panel$M), paste(panel$M, collapse = ", ")))

# -----------------------------------------------------------------------------
# Dataset loaders -> list(expr = genes x samples, group = factor(HC,RA),
#                         sex = c("F","M"), label).
# -----------------------------------------------------------------------------
load_train <- function() {
  o <- readRDS(file.path(proc, "combined_train.rds"))
  list(expr = o$expr, group = factor(o$meta$group, levels = c("HC", "RA")),
       sex = o$meta$sex, label = "Train")
}
load_internal <- function() {
  h <- readRDS(file.path(proc, "internal_val_holdout_processed.rds"))
  list(expr = h$expr, group = factor(h$meta$group, levels = c("HC", "RA")),
       sex = h$meta$sex, label = "Internal test")
}
# GSE15573 external blood: log2 if linear, collapse probes to gene by MaxMean.
load_blood <- function() {
  e <- readRDS("data/raw/GSE15573_raw.rds"); if (is.list(e)) e <- e[[1]]
  x <- exprs(e); if (max(x, na.rm = TRUE) > 50) x <- log2(x + 1)
  sym <- fData(e)[["Gene symbol"]]
  keep <- !is.na(sym) & sym != ""; x <- x[keep, ]; sym <- sym[keep]
  rmean <- rowMeans(x)
  best <- tapply(seq_along(sym), sym, function(ix) ix[which.max(rmean[ix])])
  xg <- x[unlist(best), ]; rownames(xg) <- names(best)
  p <- pData(e)
  grp <- ifelse(grepl("Rheumatoid|RA", p[["status:ch1"]], ignore.case = TRUE), "RA", "HC")
  sex <- ifelse(grepl("Female", p[["gender:ch1"]], ignore.case = TRUE), "F", "M")
  list(expr = xg, group = factor(grp, levels = c("HC", "RA")),
       sex = sex, label = "External blood")
}

datasets <- list(train = load_train(), internal = load_internal(), blood = load_blood())
for (d in datasets)
  cat(sprintf("  %-14s %d genes x %d samples | RA=%d HC=%d | F=%d M=%d\n",
              d$label, nrow(d$expr), ncol(d$expr), sum(d$group == "RA"),
              sum(d$group == "HC"), sum(d$sex == "F"), sum(d$sex == "M")))

# -----------------------------------------------------------------------------
# Helpers.
# -----------------------------------------------------------------------------
zscore <- function(M) {                              # z-score each gene (row); 0 if sd=0/NA
  t(apply(M, 1, function(v) { s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s }))
}
# DeLong CI when n is adequate; stratified BOOTSTRAP CI for small n (<20),
# where DeLong is unreliable and can print fake-precise 1.000-1.000 intervals.
auc_ci <- function(r) {
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < 20) { set.seed(GLOBAL_SEED)
    suppressWarnings(as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000)))
  } else as.numeric(ci.auc(r))
  c(auc = as.numeric(auc(r)), ci[c(1, 3)])
}

# subset one dataset to a sex and the panel genes present there
sex_panel <- function(d, sx, genes) {
  cols <- which(d$sex == sx)
  present <- genes[genes %in% rownames(d$expr)]
  list(cols = cols, present = present,
       X = d$expr[present, cols, drop = FALSE],
       y = factor(d$group[cols], levels = c("HC", "RA")))
}

# 10-fold stratified CV out-of-fold probabilities for the panel logistic model,
# with z-scoring re-estimated inside each training fold (honest).
cv_oof_prob <- function(Xz_raw, y, k = 10) {
  set.seed(GLOBAL_SEED)
  idx <- unlist(tapply(seq_along(y), y, function(ii) sample(ii)))
  fold <- integer(length(y)); fold[idx] <- rep_len(1:k, length(idx))
  prob <- rep(NA_real_, length(y))
  for (f in 1:k) {
    tr <- which(fold != f); te <- which(fold == f)
    if (length(unique(y[tr])) < 2) next
    mu <- colMeans(Xz_raw[tr, , drop = FALSE])
    sg <- apply(Xz_raw[tr, , drop = FALSE], 2, sd); sg[sg == 0 | is.na(sg)] <- 1
    Ztr <- scale(Xz_raw[tr, , drop = FALSE], center = mu, scale = sg)
    Zte <- scale(Xz_raw[te, , drop = FALSE], center = mu, scale = sg)
    df  <- data.frame(y = y[tr], Ztr, check.names = FALSE)
    fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial))
    prob[te] <- as.numeric(predict(fit, newdata = data.frame(Zte, check.names = FALSE),
                                   type = "response"))
  }
  prob
}

# -----------------------------------------------------------------------------
# Per-sex evaluation.
# -----------------------------------------------------------------------------
roc_store <- list(); panel_rows <- list(); gene_rows <- list()

eval_sex <- function(sx) {
  sexlab <- if (sx == "F") "Female" else "Male"
  genes  <- panel[[sx]]
  cat(sprintf("\n================  %s panel  ================\n", sexlab))

  ## ---- fit the LOCKED combined-panel model on TRAIN (z within train) -------
  tr <- sex_panel(datasets$train, sx, genes)
  Ztr <- as.data.frame(t(zscore(tr$X)))               # samples x genes, z-scored
  colnames(Ztr) <- make.names(tr$present)
  train_mod <- suppressWarnings(glm(tr$y ~ ., data = cbind(y = tr$y, Ztr),
                                    family = binomial))

  ## per-gene orientation fixed on TRAIN
  gene_dir <- sapply(tr$present, function(g) {
    v <- as.numeric(datasets$train$expr[g, tr$cols])
    r <- roc(tr$y, v, direction = "<", levels = c("HC", "RA"), quiet = TRUE)
    if (as.numeric(auc(r)) < 0.5) ">" else "<"
  })

  ## ---- loop datasets -------------------------------------------------------
  for (dn in names(datasets)) {
    d <- datasets[[dn]]; sp <- sex_panel(d, sx, genes)
    if (length(sp$cols) < 3 || length(unique(sp$y)) < 2) {
      cat(sprintf("  [%s] skipped (insufficient samples/classes)\n", d$label)); next }
    n_RA <- sum(sp$y == "RA"); n_HC <- sum(sp$y == "HC")
    miss <- setdiff(genes, sp$present)

    ## (a) combined panel model
    Zd <- zscore(d$expr[sp$present, sp$cols, drop = FALSE])   # genes x samples
    Zdf <- as.data.frame(t(Zd)); colnames(Zdf) <- make.names(sp$present)
    for (g in setdiff(make.names(genes), colnames(Zdf))) Zdf[[g]] <- 0  # missing gene -> mean (0)
    Zdf <- Zdf[, make.names(genes), drop = FALSE]

    if (dn == "train") {
      p_app <- as.numeric(predict(train_mod, type = "response"))
      r_app <- roc(sp$y, p_app, direction = "<", levels = c("HC", "RA"), quiet = TRUE)
      p_cv  <- cv_oof_prob(t(sp$X), sp$y, k = 10)
      r_cv  <- roc(sp$y[!is.na(p_cv)], p_cv[!is.na(p_cv)],
                   direction = "<", levels = c("HC", "RA"), quiet = TRUE)
      for (tagr in list(list("Train (apparent)", r_app), list("Train (10-fold CV)", r_cv))) {
        a <- auc_ci(tagr[[2]])
        panel_rows[[length(panel_rows) + 1]] <<- data.table(
          sex = sexlab, dataset = tagr[[1]], n = length(sp$y), n_RA = n_RA, n_HC = n_HC,
          AUC = round(a[1], 3), AUC_lo = round(a[2], 3), AUC_hi = round(a[3], 3),
          missing_genes = "")
        roc_store[[paste("panel", sexlab, tagr[[1]])]] <<- data.table(
          sex = sexlab, dataset = tagr[[1]], sens = tagr[[2]]$sensitivities,
          spec = tagr[[2]]$specificities, auc = a[1])
      }
      cat(sprintf("  [Train] panel AUC apparent=%.3f | 10-fold CV=%.3f\n",
                  as.numeric(auc(r_app)), as.numeric(auc(r_cv))))
    } else {
      p <- as.numeric(predict(train_mod, newdata = Zdf, type = "response"))
      r <- roc(sp$y, p, direction = "<", levels = c("HC", "RA"), quiet = TRUE)
      a <- auc_ci(r)
      panel_rows[[length(panel_rows) + 1]] <<- data.table(
        sex = sexlab, dataset = d$label, n = length(sp$y), n_RA = n_RA, n_HC = n_HC,
        AUC = round(a[1], 3), AUC_lo = round(a[2], 3), AUC_hi = round(a[3], 3),
        missing_genes = paste(miss, collapse = ";"))
      roc_store[[paste("panel", sexlab, d$label)]] <<- data.table(
        sex = sexlab, dataset = d$label, sens = r$sensitivities,
        spec = r$specificities, auc = a[1])
      cat(sprintf("  [%s] panel AUC=%.3f (95%% CI %.3f-%.3f)%s\n", d$label,
                  a[1], a[2], a[3],
                  if (length(miss)) paste0(" | missing: ", paste(miss, collapse = ",")) else ""))
    }

    ## (b) per-gene ROC, orientation fixed on TRAIN
    for (g in genes) {
      if (!g %in% sp$present) {
        gene_rows[[length(gene_rows) + 1]] <<- data.table(
          sex = sexlab, dataset = d$label, gene = g, present = FALSE,
          AUC = NA, AUC_lo = NA, AUC_hi = NA, concordant = NA); next }
      v   <- as.numeric(d$expr[g, sp$cols]); dir <- gene_dir[[g]]
      r   <- roc(sp$y, v, direction = dir, levels = c("HC", "RA"), quiet = TRUE)
      a   <- auc_ci(r)
      gene_rows[[length(gene_rows) + 1]] <<- data.table(
        sex = sexlab, dataset = d$label, gene = g, present = TRUE,
        AUC = round(a[1], 3), AUC_lo = round(a[2], 3), AUC_hi = round(a[3], 3),
        concordant = a[1] >= 0.5)
      roc_store[[paste("gene", sexlab, d$label, g)]] <<- data.table(
        sex = sexlab, dataset = d$label, gene = g, sens = r$sensitivities,
        spec = r$specificities, auc = a[1])
    }
  }
}

eval_sex("F"); eval_sex("M")

panel_tab <- rbindlist(panel_rows); gene_tab <- rbindlist(gene_rows)
fwrite(panel_tab, file.path(tab, "mr_roc_panel_auc.csv"))
fwrite(gene_tab,  file.path(tab, "mr_roc_pergene_auc.csv"))
saveRDS(list(roc = roc_store, panel = panel_tab, gene = gene_tab, panels = panel),
        file.path(procN, "mr_roc_objects.rds"))

cat("\n=======================  PANEL AUC  =======================\n"); print(panel_tab)
cat("\n=======================  PER-GENE AUC  ====================\n"); print(gene_tab)
cat("\nSaved mr_roc_panel_auc.csv, mr_roc_pergene_auc.csv, mr_roc_objects.rds\n")
cat("\n---- key package versions ----\n")
for (p in c("pROC", "Biobase", "data.table"))
  cat(sprintf("  %-12s %s\n", p, as.character(packageVersion(p))))
