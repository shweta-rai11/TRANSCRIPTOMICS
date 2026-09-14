#!/usr/bin/env Rscript
# Training phase: tune five classifiers (logreg, SVM, KNN, RF, ANN) on the locked sex-specific consensus panel via repeated CV.
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

# within-dataset z-score (rows = genes)
zrows <- function(M) t(apply(M, 1, function(v) {
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
}))

# algorithm grids (fixed, pre-specified; no post hoc adjustment)
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

# pooled out-of-fold ROC at the best tuning parameters (honest "train" curve)
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
             built = "06_models/37_model_training_ml_algorithms.R"),
        file.path(procN, "ml_algo_models.rds"))

cat("\n====================  HYPERPARAMETER TUNING SUMMARY  ====================\n")
print(tuning_tbl, width = 200)
cat("\nSaved ML_hyperparameter_tuning.csv, ml_algo_models.rds\n")
