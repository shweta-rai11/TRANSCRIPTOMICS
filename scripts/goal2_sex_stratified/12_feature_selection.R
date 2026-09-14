#!/usr/bin/env Rscript
# Sex-stratified feature selection: LASSO + Random Forest + SVM-RFE consensus on MR-prioritised candidate genes, females and males analysed separately.

suppressMessages({
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(caret)          # 10-fold CV grid search for RF mtry tuning
  library(data.table)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
set.seed(GLOBAL_SEED)

proc  <- "data/processed"              # combined_train.rds lives here
procN <- "data/processed/new"          # NEW ml_features.rds output
tab   <- "results/tables"          # NEW FS_input read + mr_fs write
dir.create(tab,   showWarnings = FALSE, recursive = TRUE)
dir.create(procN, showWarnings = FALSE, recursive = TRUE)

# ---- load training cohort and MR candidate lists ---------------------------
o    <- readRDS(file.path(proc, "combined_train.rds"))
expr <- o$expr
meta <- o$meta

mr_female <- fread(file.path(tab, "FS_input_female.csv"))
mr_male   <- fread(file.path(tab, "FS_input_male.csv"))

# NEW symmetric disease-module MR panels (counts differ from the old 14/40)
cat(sprintf("NEW MR-prioritised candidates: female %d, male %d\n", nrow(mr_female), nrow(mr_male)))
mr_female <- mr_female[mr_female$gene %in% rownames(expr)]
mr_male   <- mr_male[mr_male$gene   %in% rownames(expr)]
stopifnot(nrow(mr_female) >= 2L, nrow(mr_male) >= 2L)

# -----------------------------------------------------------------------------
# SVM-RFE ranking: recursive elimination by squared linear-SVM weight.
# Removes the single least-important feature each round -> full ranking vector
# ordered from MOST to LEAST important.
# -----------------------------------------------------------------------------
svm_rfe_rank <- function(X, y, cost = 1) {
  feats   <- colnames(X)
  ranking <- character(0)                 # filled from least to most important
  while (length(feats) > 1) {
    m    <- svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = TRUE, cost = cost)
    w2   <- ((t(m$coefs) %*% m$SV)[1, ])^2
    drop <- names(sort(w2))[1]            # smallest weight = least important
    ranking <- c(drop, ranking)
    feats   <- setdiff(feats, drop)
  }
  c(feats, ranking)                        # most-important first
}

# 10-fold CV classification error of a linear SVM using the top-k ranked feats,
# for k = 1..length(rank). Used to choose the SVM-RFE panel size objectively.
svm_rfe_curve <- function(X, y, rank, cost = 1) {
  ks <- seq_along(rank)
  err <- sapply(ks, function(k) {
    set.seed(GLOBAL_SEED)
    acc <- svm(X[, rank[1:k], drop = FALSE], y, kernel = "linear",
               scale = TRUE, cost = cost, cross = 10)$tot.accuracy
    1 - acc / 100
  })
  list(k = ks, err = err, best = ks[which.min(err)], besterr = min(err))
}

# Tune the linear-SVM regularisation parameter cost (C) by 10-fold CV over a
# fixed grid; return the cost minimising cross-validated error. This is the
# hyperparameter tuned for the SVM-RFE selector.
SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)
tune_svm_cost <- function(X, y, grid = SVM_COST_GRID) {
  set.seed(GLOBAL_SEED)
  tc <- tune(svm, train.x = X, train.y = y, kernel = "linear", scale = TRUE,
             ranges = list(cost = grid),
             tunecontrol = tune.control(sampling = "cross", cross = 10))
  tc$best.parameters$cost
}

# Tune the random-forest mtry (number of variables tried per split) by 10-fold
# CV over a grid, at fixed ntree; return the best mtry (ties -> smaller mtry).
rf_mtry_grid <- function(p) sort(unique(pmin(p, c(1, 2, floor(sqrt(p)),
                                                  floor(p / 3), floor(p / 2), p))))
tune_rf_mtry <- function(X, y, ntree, grid) {
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE,
                       summaryFunction = twoClassSummary)
  set.seed(GLOBAL_SEED)
  fit <- train(x = X, y = y, method = "rf", metric = "ROC",
               trControl = ctrl, tuneGrid = expand.grid(mtry = grid),
               ntree = ntree, importance = TRUE)
  bt <- fit$bestTune$mtry
  list(mtry = bt, fit = fit)
}

# -----------------------------------------------------------------------------
# Per-sex feature selection.
# -----------------------------------------------------------------------------
run_fs <- function(sex_code) {
  sexlab <- if (sex_code == "F") "Female" else "Male"
  mr     <- if (sex_code == "F") mr_female else mr_male
  genes  <- mr$gene
  cols   <- meta$sample[meta$sex == sex_code]

  # design matrix: samples x candidate genes; response RA vs HC (HC reference)
  X <- t(expr[genes, cols, drop = FALSE])
  colnames(X) <- make.names(genes)                 # safe names for formulas
  lk   <- setNames(genes, make.names(genes))        # safe name -> gene symbol
  back <- function(v) unname(lk[v])
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))

  cat(sprintf("\n================  %s  ================\n", sexlab))
  cat(sprintf("candidate MR genes: %d | samples: %d (RA=%d, HC=%d)\n",
              ncol(X), nrow(X), sum(y == "RA"), sum(y == "HC")))

  ## (1) LASSO logistic regression -------------------------------------------
  set.seed(GLOBAL_SEED)
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10,
                  type.measure = "deviance")
  co_min <- coef(cv, s = "lambda.min")[-1, 1]
  co_1se <- coef(cv, s = "lambda.1se")[-1, 1]
  lasso_min <- back(names(co_min)[co_min != 0])
  lasso_1se <- back(names(co_1se)[co_1se != 0])
  lasso <- lasso_min                                # primary = lambda.min
  cat(sprintf("  LASSO : lambda.min=%.4f -> %d genes | lambda.1se=%.4f -> %d genes\n",
              cv$lambda.min, length(lasso_min), cv$lambda.1se, length(lasso_1se)))

  ## (2) Random Forest importance --------------------------------------------
  ## Hyperparameter: mtry (variables tried per split) is tuned by 10-fold CV
  ## (grid search, ROC metric) at a fixed ntree = 1000 (raised from the 500
  ## default for stable importance). The final 1000-tree forest is then fitted
  ## with the tuned mtry. A gene is RF-selected if its Mean Decrease in Gini
  ## exceeds the MEAN Gini across all candidate genes (above-average-importance
  ## cutoff), keeping only genes contributing more than average to node purity.
  mtry_grid <- rf_mtry_grid(ncol(X))
  rf_tune   <- tune_rf_mtry(X, y, ntree = 1000, grid = mtry_grid)
  best_mtry <- rf_tune$mtry
  set.seed(GLOBAL_SEED)
  rf <- randomForest(X, y, importance = TRUE, ntree = 1000, mtry = best_mtry)
  gini    <- sort(rf$importance[, "MeanDecreaseGini"], decreasing = TRUE)
  gini_thr <- mean(gini)
  rf_sel  <- back(names(gini)[gini > gini_thr])
  cat(sprintf("  RF    : ntree=1000, tuned mtry=%d (grid %s), Gini cutoff=%.3f -> %d genes\n",
              best_mtry, paste(mtry_grid, collapse = "/"), gini_thr, length(rf_sel)))

  ## (3) SVM-RFE -------------------------------------------------------------
  ## Hyperparameter: the linear-SVM cost C is tuned by 10-fold CV (grid search)
  ## before elimination. Using the tuned cost, features are eliminated one per
  ## round (finest granularity) to produce a full ranking, and the panel size
  ## is chosen as the number of top-ranked features minimising 10-fold CV error.
  best_cost <- tune_svm_cost(X, y)
  set.seed(GLOBAL_SEED)
  rank  <- svm_rfe_rank(X, y, cost = best_cost)
  curve <- svm_rfe_curve(X, y, rank, cost = best_cost)
  svm_sel <- back(rank[1:curve$best])
  cat(sprintf("  SVM-RFE: tuned cost=%s (grid %s), CV-optimal size=%d (err=%.4f) -> %d genes\n",
              best_cost, paste(SVM_COST_GRID, collapse = "/"),
              curve$best, curve$besterr, length(svm_sel)))

  ## consensus of the three independent selectors ----------------------------
  sets <- list(LASSO = lasso, RandomForest = rf_sel, SVM_RFE = svm_sel)
  consensus <- sort(Reduce(intersect, sets))
  union_all <- sort(Reduce(union, sets))
  cat(sprintf("  CONSENSUS (LASSO n RF n SVM-RFE): %d genes -> %s\n",
              length(consensus),
              ifelse(length(consensus), paste(consensus, collapse = ", "), "none")))

  ## NOTE: univariate ROC/AUC of these genes is deferred to a later script,
  ## on request; this script stops at the cross-method common-gene overlap.

  list(sexlab = sexlab, mr = mr, mr_genes = genes,
       sets = sets, consensus = consensus, union = union_all,
       lasso_min = lasso_min, lasso_1se = lasso_1se,
       cv = cv, rf = rf, gini = gini, gini_thr = gini_thr,
       svm_rank = rank, svm_curve = curve, samples = cols,
       tuned = list(lasso_lambda_min = cv$lambda.min,
                    lasso_lambda_1se = cv$lambda.1se,
                    rf_mtry = best_mtry, rf_mtry_grid = mtry_grid,
                    svm_cost = best_cost, svm_cost_grid = SVM_COST_GRID))
}

F <- run_fs("F")
M <- run_fs("M")

# ---- write per-sex tables ---------------------------------------------------
for (r in list(F, M)) {
  s <- tolower(r$sexlab)
  bym <- rbindlist(lapply(names(r$sets), function(n)
    data.table(method = n, gene = r$sets[[n]])))
  fwrite(bym, file.path(tab, sprintf("mr_fs_selected_bymethod_%s.csv", s)))
  fwrite(data.table(gene = r$consensus), file.path(tab, sprintf("mr_fs_consensus_%s.csv", s)))
}

# ---- cross-sex summary ------------------------------------------------------
shared <- intersect(F$consensus, M$consensus)
summary_tab <- data.table(
  sex             = c("Female", "Male"),
  n_samples       = c(length(F$samples), length(M$samples)),
  n_mr_candidates = c(length(F$mr_genes), length(M$mr_genes)),
  n_lasso         = c(length(F$sets$LASSO), length(M$sets$LASSO)),
  n_rf            = c(length(F$sets$RandomForest), length(M$sets$RandomForest)),
  n_svmrfe        = c(length(F$sets$SVM_RFE), length(M$sets$SVM_RFE)),
  n_consensus     = c(length(F$consensus), length(M$consensus)),
  tuned_lasso_lambda_min = c(round(F$tuned$lasso_lambda_min, 4), round(M$tuned$lasso_lambda_min, 4)),
  tuned_rf_mtry   = c(F$tuned$rf_mtry, M$tuned$rf_mtry),
  tuned_svm_cost  = c(F$tuned$svm_cost, M$tuned$svm_cost),
  consensus_genes = c(paste(F$consensus, collapse = "; "),
                      paste(M$consensus, collapse = "; ")))
fwrite(summary_tab, file.path(tab, "mr_fs_summary.csv"))

cat("\n=======================  SUMMARY  =======================\n")
print(summary_tab)
cat(sprintf("\nGenes shared between female and male consensus panels: %s\n",
            ifelse(length(shared), paste(shared, collapse = ", "), "none (no overlap)")))

# ---- save objects (drop-in structure for scripts 18/19) --------------------
saveRDS(list(female = F, male = M, expr = expr, meta = meta,
             seed = GLOBAL_SEED, built = "goal2_sex_stratified/12_feature_selection.R"),
        file.path(procN, "ml_features.rds"))

cat("\nSaved -> data/processed/new/ml_features.rds and results/tables/mr_fs_*.csv\n")
cat("\n---- sessionInfo (key packages) ----\n")
for (p in c("glmnet", "randomForest", "e1071", "caret", "data.table"))
  cat(sprintf("  %-14s %s\n", p, as.character(packageVersion(p))))
cat(sprintf("  R %s\n", paste0(R.version$major, ".", R.version$minor)))
