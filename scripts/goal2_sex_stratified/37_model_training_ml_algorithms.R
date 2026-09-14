#!/usr/bin/env Rscript
# =============================================================================
# 37_model_training_ml_algorithms.R
# -----------------------------------------------------------------------------
# TRAINING PHASE for a five-algorithm diagnostic-model comparison, requested in
# addition to the elastic-net panel model already reported in section 2.10.
# Where section 2.10 fits ONE classifier (logistic regression, optionally
# elastic-net-penalised) on the sex-stratified consensus panel, this script
# fits FIVE standard classifier families on that SAME locked panel so the
# algorithms can be compared head-to-head under identical features:
#
#     Logistic regression (glm)      Support vector machine (RBF kernel)
#     k-nearest neighbours            Artificial neural network (single
#     Random forest                   hidden layer, nnet)
#
# FEATURES: the sex-specific consensus panel already selected and locked by
#   the LASSO n RF n SVM-RFE pipeline (12_feature_selection.R), read from
#   data/processed/new/ml_features.rds:
#       Female (6): C6orf136, ESYT1, GNL1, IKZF3, MED1, SMARCC2
#       Male   (6): ESYT1, HLA-DMA, INPP5B, MED1, SMARCC2, VPS52
#   No new feature selection is performed here -- reusing the locked panel
#   keeps this comparison about ALGORITHM choice, not panel choice, and avoids
#   re-introducing the feature-selection leakage that 14_/16_ already went to
#   considerable lengths to eliminate.
#
# TRAINING DATA: data/processed/combined_train.rds (blood, GSE93272+GSE110169,
#   the same 70% discovery cohort used throughout this chapter). Each sex is
#   modelled separately throughout (never pooled).
#
# PREPROCESSING: each panel gene is z-scored WITHIN the training set (own
#   mean/sd), matching the convention already justified in 18b/20_ ("each gene
#   is standardised within each dataset ... unsupervised, removes platform-
#   specific location/scale offsets"). The same within-dataset standardisation
#   is re-applied independently to each test set in 38_ rather than reusing the
#   training mean/sd, because one test set (synovium, RNA-seq log2-CPM) is on a
#   different platform and tissue from training (blood microarray) and a
#   frozen microarray scale would not transfer.
#
# HYPERPARAMETER TUNING: caret::train(), 5-fold cross-validation repeated 5
#   times (identical scheme for both sexes; unlike the outer resampling of
#   section 2.10's nested CV, this is tuning only, not the reported
#   performance estimate -- see 38_ for that), metric = ROC, seed 1234 fixed
#   before every model. Grids:
#       Logistic regression : none (glm has no hyperparameter; included as
#                              the common baseline against which the other
#                              four are compared)
#       SVM (svmRadial)     : C in {0.1,0.25,0.5,1,2,4,8,16}, sigma in
#                              {0.01,0.05,0.1,0.5,1}
#       KNN                 : k in {3,5,7,9,11,13,15}
#       Random forest       : mtry in {1,...,6} (all possible values at p=6),
#                              ntree fixed at 500
#       ANN (nnet)          : size in {1,3,5,7,9}, decay in
#                              {0,0.001,0.01,0.1,1}
#   The male stratum (n=38) is tuned with the same 5x5 scheme as the female
#   stratum (n=145) for comparability; this is smaller-sample than ideal for
#   a 5-fold split with 17 RA cases, and results are read with that in mind
#   (see evidence-tier flags in 38_).
#
# HONEST "TRAIN" ESTIMATE: for every algorithm, the pooled out-of-fold
#   predicted probabilities at the best tuning parameters (across all 5x5=25
#   resamples, averaged per sample) are kept and turned into an honest,
#   non-resubstitution "Train (resampled CV)" ROC in 38_/39_, exactly as
#   section 2.10 already does for the elastic-net panel. Resubstitution
#   (apparent) accuracy is not reported anywhere in this comparison.
#
# Outputs:
#   data/processed/new/ml_algo_models.rds  (fitted caret models, OOF preds,
#                                            best hyperparameters, panel, seed)
#   results/tables/ML_hyperparameter_tuning.csv
# =============================================================================
suppressMessages({
  library(caret); library(data.table); library(pROC)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
set.seed(GLOBAL_SEED)
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(procN, showWarnings = FALSE, recursive = TRUE)
dir.create(tab,   showWarnings = FALSE, recursive = TRUE)

ml <- readRDS(file.path(procN, "ml_features.rds"))
PANEL <- list(F = ml$female$consensus, M = ml$male$consensus)

o <- readRDS(file.path(proc, "combined_train.rds"))
expr <- o$expr; meta <- as.data.table(o$meta)

cat(sprintf("Female panel (%d): %s\n", length(PANEL$F), paste(PANEL$F, collapse = ", ")))
cat(sprintf("Male panel   (%d): %s\n", length(PANEL$M), paste(PANEL$M, collapse = ", ")))

## ---- within-dataset z-score (rows = genes) --------------------------------
zrows <- function(M) t(apply(M, 1, function(v) {
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
}))

## ---- algorithm grids (fixed, pre-specified; no post hoc adjustment) -------
ALGOS <- list(
  logreg = list(method = "glm",       grid = NULL, extra = list(family = "binomial"), label = "Logistic regression"),
  svm    = list(method = "svmRadial", grid = expand.grid(C = c(0.1,0.25,0.5,1,2,4,8,16), sigma = c(0.01,0.05,0.1,0.5,1)), extra = list(), label = "SVM (RBF)"),
  knn    = list(method = "knn",       grid = data.frame(k = seq(3, 15, 2)), extra = list(), label = "k-NN"),
  rf     = list(method = "rf",        grid = data.frame(mtry = 1:6), extra = list(ntree = 500), label = "Random forest"),
  ann    = list(method = "nnet",      grid = expand.grid(size = c(1,3,5,7,9), decay = c(0,0.001,0.01,0.1,1)), extra = list(trace = FALSE, maxit = 200, MaxNWts = 2000), label = "ANN (nnet)")
)

fit_tuned <- function(X, y, spec, seed = GLOBAL_SEED) {
  ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                        classProbs = TRUE, summaryFunction = twoClassSummary,
                        savePredictions = "final")
  set.seed(seed)
  args <- c(list(x = X, y = y, method = spec$method, trControl = ctrl,
                 metric = "ROC", tuneGrid = spec$grid), spec$extra)
  suppressWarnings(do.call(caret::train, args))
}

## pooled out-of-fold ROC at the best tuning parameters (honest "train" curve)
oof_roc <- function(m, y_full) {
  p <- as.data.table(m$pred)
  agg <- p[, .(prob = mean(RA)), by = rowIndex]
  setorder(agg, rowIndex)
  obs <- y_full[agg$rowIndex]
  roc(obs, agg$prob, direction = "<", levels = c("HC", "RA"), quiet = TRUE)
}

best_tune_string <- function(m, spec) {
  if (is.null(spec$grid)) return("(none)")
  bt <- m$bestTune
  paste(sprintf("%s=%s", names(bt), sapply(bt, format)), collapse = ", ")
}

run_sex <- function(sx) {
  sexlab <- if (sx == "F") "Female" else "Male"
  genes <- PANEL[[sx]]
  cols  <- meta$sample[meta$sex == sx]
  y     <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))
  Xz    <- t(zrows(expr[genes, cols, drop = FALSE]))
  colnames(Xz) <- make.names(genes)
  Xz <- as.data.frame(Xz, check.names = FALSE)

  cat(sprintf("\n================  %s (n=%d, RA=%d, HC=%d)  ================\n",
              sexlab, length(y), sum(y == "RA"), sum(y == "HC")))

  models <- list()
  rows <- list()
  for (a in names(ALGOS)) {
    spec <- ALGOS[[a]]
    cat(sprintf("  tuning %-22s ...", spec$label))
    m <- fit_tuned(Xz, y, spec)
    cvroc <- max(m$results$ROC, na.rm = TRUE)
    cvsd  <- m$results$ROCSD[which.max(m$results$ROC)]
    cat(sprintf(" done | CV ROC=%.3f +/- %.3f | best: %s\n",
                cvroc, cvsd, best_tune_string(m, spec)))
    models[[a]] <- list(model = m, oof = oof_roc(m, y), genes = genes, colnames = colnames(Xz))
    rows[[length(rows) + 1]] <- data.table(
      sex = sexlab, algorithm = spec$label, n_train = length(y),
      tuning_grid_size = if (is.null(spec$grid)) 0L else nrow(spec$grid),
      best_hyperparameters = best_tune_string(m, spec),
      cv_ROC_mean = round(cvroc, 3), cv_ROC_sd = round(cvsd, 3))
  }
  list(models = models, tuning = rbindlist(rows), genes = genes, colnames = colnames(Xz), y = y)
}

resF <- run_sex("F")
resM <- run_sex("M")

tuning_tbl <- rbindlist(list(resF$tuning, resM$tuning))
fwrite(tuning_tbl, file.path(tab, "ML_hyperparameter_tuning.csv"))

saveRDS(list(female = resF, male = resM, panel = PANEL, seed = GLOBAL_SEED,
             built = "goal2_sex_stratified/37_model_training_ml_algorithms.R"),
        file.path(procN, "ml_algo_models.rds"))

cat("\n====================  HYPERPARAMETER TUNING SUMMARY  ====================\n")
print(tuning_tbl, width = 200)
cat("\nSaved ML_hyperparameter_tuning.csv, ml_algo_models.rds\n")
