#!/usr/bin/env Rscript
# =============================================================================
# 18c_mr_nested_cv.R
# -----------------------------------------------------------------------------
# LEAKAGE-FREE nested cross-validation of the MR-anchored diagnostic pipeline,
# per sex. This is the HONEST replacement for the plain ("flat") training AUC in
# 18b: because the 3-method consensus panel in 18b was chosen using the whole
# training set and then scored by ordinary CV on that same set, its training and
# internal AUCs are inflated by feature-selection bias (Ambroise & McLachlan,
# PNAS 2002; Simon et al., JNCI 2003). Here EVERY step that looks at the outcome
# labels - the LASSO, Random Forest and SVM-RFE selection AND the 3-method
# consensus AND the logistic model - is redone INSIDE each outer training fold,
# and the held-out fold is never seen during selection or fitting.
#
# DESIGN
#   Candidate universe (fixed, legitimate): the within-sex MR-prioritised genes
#     (14 female / 40 male). These come from an EXTERNAL MR analysis on GWAS
#     summary data, not from this expression matrix, so keeping them fixed
#     across folds is not leakage.
#   Outer loop : repeated stratified k-fold.
#       Female : 10-fold x 5 repeats   (n = 145)
#       Male   : 5-fold  x 10 repeats  (n = 38; smaller k, more repeats to
#                                       stabilise the estimate at low n)
#   Inside each outer-TRAIN fold (never touching the outer-test fold):
#       (1) LASSO logistic (cv.glmnet, inner 5-fold CV tunes lambda) -> nonzero genes
#       (2) Random Forest (ntree = 500) -> genes with above-mean Gini importance
#       (3) SVM-RFE (linear, cost = 1) -> size by inner 5-fold CV error
#       consensus = LASSO n RF n SVM-RFE  (fallback: union, then LASSO, then the
#           full MR set, so the fold always yields a usable panel)
#       classifier = logistic regression on the consensus genes, standardised
#           with the FOLD-TRAIN mean/sd (frozen, applied to the test fold).
#   Predict the untouched outer-test fold -> pool out-of-fold probabilities.
#
# REPORTS
#   pooled out-of-fold AUC (+95% DeLong CI) and the per-repeat AUC distribution;
#   the in-fold consensus-gene re-selection frequency (stability of the panel).
# Females and males are handled completely separately throughout.
#
# Seeds: outer folds seeded 1000+repeat (reproducible). Outputs:
#   results/tables/mr_nested_cv_summary.csv
#   results/tables/mr_nested_cv_stability_{female,male}.csv
#   data/processed/mr_nested_objects.rds   (nested ROC coords for the figure)
# =============================================================================
suppressMessages({library(glmnet); library(randomForest); library(e1071)
                  library(pROC); library(caret); library(data.table)})
options(stringsAsFactors = FALSE)
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(procN, showWarnings = FALSE, recursive = TRUE); dir.create(tab, showWarnings = FALSE, recursive = TRUE)

o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr; meta <- as.data.table(o$meta)
mrF <- fread(file.path(tab, "FS_input_female.csv"))$gene
mrM <- fread(file.path(tab, "FS_input_male.csv"))$gene

# ---- in-fold selectors (fixed, modest hyperparameters; lambda self-tuned) ----
svm_rank <- function(X, y, cost = 1) {
  feats <- colnames(X); ranking <- character(0)
  while (length(feats) > 1) {
    m <- svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = TRUE, cost = cost)
    w2 <- ((t(m$coefs) %*% m$SV)[1, ])^2
    d  <- names(sort(w2))[1]; ranking <- c(d, ranking); feats <- setdiff(feats, d)
  }
  c(feats, ranking)
}
svm_size <- function(X, y, rank, cost = 1) {
  err <- sapply(seq_along(rank), function(k) {
    1 - svm(X[, rank[1:k], drop = FALSE], y, kernel = "linear", scale = TRUE,
            cost = cost, cross = 5)$tot.accuracy / 100 })
  rank[1:which.min(err)]
}
# return the in-fold 3-method consensus (make.names space), with fallbacks
select_infold <- function(X, y) {
  lasso <- character(0)
  cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5),
                 error = function(e) NULL)
  if (!is.null(cv)) { co <- coef(cv, s = "lambda.min")[-1, 1]; lasso <- colnames(X)[co != 0] }
  rf   <- randomForest(X, y, ntree = 500)
  gini <- rf$importance[, "MeanDecreaseGini"]; rfs <- names(gini)[gini > mean(gini)]
  svs  <- svm_size(X, y, svm_rank(X, y))
  cons <- Reduce(intersect, list(lasso, rfs, svs))
  uni  <- Reduce(union, list(lasso, rfs, svs))
  panel <- if (length(cons) >= 2) cons else if (length(uni) >= 1) uni
           else if (length(lasso) >= 1) lasso else colnames(X)
  list(panel = panel, consensus = cons)
}

# ---- nested CV for one sex --------------------------------------------------
run_nested <- function(sex_code, genes, kfold, repeats) {
  sexlab <- if (sex_code == "F") "Female" else "Male"
  cols <- meta$sample[meta$sex == sex_code]
  genes <- genes[genes %in% rownames(expr)]
  X0 <- t(expr[genes, cols, drop = FALSE]); colnames(X0) <- make.names(genes)
  lk <- setNames(genes, make.names(genes))
  y  <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))

  preds <- list(); cons_freq <- integer(0); rep_auc <- numeric(0); nfits <- 0L
  for (rp in seq_len(repeats)) {
    set.seed(1000 + rp)
    folds <- createFolds(y, k = kfold, returnTrain = FALSE)
    rp_pred <- data.table()
    for (fi in seq_along(folds)) {
      te <- folds[[fi]]; tr <- setdiff(seq_along(y), te)
      if (length(unique(y[tr])) < 2) next
      Xtr <- X0[tr, , drop = FALSE]; Xte <- X0[te, , drop = FALSE]; ytr <- y[tr]
      sel <- select_infold(Xtr, ytr); nfits <- nfits + 1L
      cons_freq <- c(cons_freq, sel$consensus)
      P <- sel$panel
      mu <- colMeans(Xtr[, P, drop = FALSE])
      sg <- apply(Xtr[, P, drop = FALSE], 2, sd); sg[sg == 0 | is.na(sg)] <- 1
      Ztr <- scale(Xtr[, P, drop = FALSE], center = mu, scale = sg)
      Zte <- scale(Xte[, P, drop = FALSE], center = mu, scale = sg)
      fit <- suppressWarnings(glm(ytr ~ ., data = data.frame(ytr, Ztr, check.names = FALSE),
                                  family = binomial))
      p <- as.numeric(predict(fit, newdata = data.frame(Zte, check.names = FALSE),
                              type = "response"))
      rp_pred <- rbind(rp_pred, data.table(sample = rownames(Xte), prob = p, obs = y[te]))
    }
    rr <- roc(rp_pred$obs, rp_pred$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
    rep_auc <- c(rep_auc, as.numeric(auc(rr))); preds[[rp]] <- rp_pred
  }
  allp <- rbindlist(preds)
  agg  <- allp[, .(prob = mean(prob), obs = obs[1]), by = sample]
  ro   <- roc(agg$obs, agg$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  ci   <- as.numeric(ci.auc(ro))
  freq <- sort(table(lk[cons_freq]), decreasing = TRUE)
  cat(sprintf("[%s] nested-CV pooled AUC = %.3f (95%% CI %.3f-%.3f) | per-repeat %.3f +/- %.3f | %d folds\n",
              sexlab, as.numeric(auc(ro)), ci[1], ci[3], mean(rep_auc), sd(rep_auc), nfits))
  list(sexlab = sexlab, auc = as.numeric(auc(ro)), ci = ci, rep_auc = rep_auc,
       roc = ro, agg = agg, cons_freq = freq, nfits = nfits)
}

cat("Nested CV (feature selection redone inside every fold; a few minutes)...\n")
F <- run_nested("F", mrF, kfold = 10, repeats = 5)
M <- run_nested("M", mrM, kfold = 5,  repeats = 10)

# ---- stability tables (how often each panel gene re-enters the consensus) ----
stab <- function(res, panel) {
  fr <- as.integer(res$cons_freq); names(fr) <- names(res$cons_freq)
  cnt <- fr[panel]; cnt[is.na(cnt)] <- 0L
  data.table(gene = panel, reselect_pct = round(100 * as.numeric(cnt) / res$nfits, 1))[order(-reselect_pct)]
}
ml <- readRDS(file.path(procN, "ml_features.rds"))
sfF <- stab(F, ml$female$consensus); sfM <- stab(M, ml$male$consensus)
fwrite(sfF, file.path(tab, "mr_nested_cv_stability_female.csv"))
fwrite(sfM, file.path(tab, "mr_nested_cv_stability_male.csv"))

summ <- data.table(
  sex = c("Female", "Male"),
  flat_train_CV_AUC = c(0.792, 0.966),          # optimistic numbers from 18b (for contrast)
  nested_CV_AUC = c(round(F$auc, 3), round(M$auc, 3)),
  nested_CI = c(sprintf("%.3f-%.3f", F$ci[1], F$ci[3]), sprintf("%.3f-%.3f", M$ci[1], M$ci[3])),
  per_repeat_mean = c(round(mean(F$rep_auc), 3), round(mean(M$rep_auc), 3)),
  per_repeat_sd   = c(round(sd(F$rep_auc), 3), round(sd(M$rep_auc), 3)))
fwrite(summ, file.path(tab, "mr_nested_cv_summary.csv"))

# ROC coordinates for the figure (labelled "Train (nested CV)")
roc_coords <- function(res) data.table(sex = res$sexlab, dataset = "Train (nested CV)",
  sens = res$roc$sensitivities, spec = res$roc$specificities, auc = res$auc)
saveRDS(list(female = F, male = M, summary = summ,
             roc = list(F = roc_coords(F), M = roc_coords(M))),
        file.path(procN, "mr_nested_objects.rds"))

cat("\n============== NESTED CV vs FLAT (leaky) TRAIN CV ==============\n"); print(summ)
cat("\nFemale consensus re-selection frequency:\n"); print(sfF)
cat("Male consensus re-selection frequency:\n");   print(sfM)
cat("\nSaved mr_nested_cv_summary.csv, stability tables, mr_nested_objects.rds\n")
