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
#   Candidate universe (fixed across folds). The within-sex MR-prioritised genes,
#     read from FS_input_{female,male}.csv. DO NOT quote a count here - it has
#     changed three times (14/40, then 74/55, then 32/25) and a header comment
#     is the wrong place to carry a result. Quote mr_fs_summary.csv.
#
#     SCOPE OF THIS DESIGN - stated precisely, because an earlier version of this
#     comment claimed the fixed universe "is not leakage", which is only half true.
#     The universe is an EXTERNALLY FILTERED INTERNAL LIST:
#       - genes SUBMITTED to MR = disease module (06) INTERSECT sex DEG (05),
#         both computed on the whole training partition USING THE LABELS;
#       - the FILTER applied to them = MR against external GWAS/eQTL summary
#         statistics, which used no expression data from this study.
#     Holding it fixed therefore means an outer-fold test sample helped define
#     the gene list its own prediction is built from. This loop corrects the
#     selection bias of the 3-selector stage and the classifier; it does NOT
#     correct the bias of the upstream DEG/WGCNA/MR stages. Treat the nested
#     number as an UPPER BOUND, and rest the diagnostic claim on the sealed
#     holdout and the external cohorts, which are unaffected. See thesis 2.9.2.
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

# -----------------------------------------------------------------------------
# THE FLAT ("LEAKY") COMPARISON, NOW COMPUTED RATHER THAN HARD-CODED.
# -----------------------------------------------------------------------------
# Earlier versions of this script wrote flat_train_CV_AUC = c(0.792, 0.966) as a
# LITERAL, carried over from a superseded run of a different script, and that
# literal was then written into mr_nested_cv_summary.csv as though it had been
# computed here. It had not been. Any reader comparing the two columns was being
# shown a number with no code behind it, and the value no longer corresponded to
# the current candidate sets (the cis filter and the FDR change in 10_MR.R both
# moved it). It is computed here instead.
#
# WHAT "FLAT" MEANS. The selection is done ONCE on the whole training set - the
# 3-method consensus from ml_features.rds - and only the logistic model is then
# cross-validated. Every fold therefore scores genes that were chosen using the
# held-out samples. That is precisely the feature-selection bias Ambroise &
# McLachlan (PNAS 2002) describe, and reproducing it here is the point: the gap
# between this column and the nested column IS the bias, measured on this data.
flat_cv_auc <- function(sex_code, panel_genes) {
  cols <- meta$sample[meta$sex == sex_code]
  g <- unique(panel_genes[panel_genes %in% rownames(expr)])
  if (length(g) < 2) return(NA_real_)
  X <- t(expr[g, cols, drop = FALSE]); colnames(X) <- make.names(g)
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))
  k <- if (sex_code == "F") 10 else 5
  reps <- if (sex_code == "F") 5 else 10
  probs <- data.table()
  for (rp in seq_len(reps)) {
    set.seed(1000 + rp)                       # same seed policy as the nested loop
    folds <- createFolds(y, k = k, returnTrain = FALSE)
    for (fi in seq_along(folds)) {
      te <- folds[[fi]]; tr <- setdiff(seq_along(y), te)
      if (length(unique(y[tr])) < 2) next
      mu <- colMeans(X[tr, , drop = FALSE])
      sg <- apply(X[tr, , drop = FALSE], 2, sd); sg[sg == 0 | is.na(sg)] <- 1
      Ztr <- scale(X[tr, , drop = FALSE], center = mu, scale = sg)
      Zte <- scale(X[te, , drop = FALSE], center = mu, scale = sg)
      fit <- suppressWarnings(glm(y[tr] ~ ., data = data.frame(y = y[tr], Ztr,
                                                               check.names = FALSE),
                                  family = binomial))
      p <- as.numeric(predict(fit, newdata = data.frame(Zte, check.names = FALSE),
                              type = "response"))
      probs <- rbind(probs, data.table(sample = rownames(X)[te], prob = p, obs = y[te]))
    }
  }
  agg <- probs[, .(prob = mean(prob), obs = obs[1]), by = sample]
  as.numeric(auc(roc(agg$obs, agg$prob, levels = c("HC", "RA"),
                     direction = "<", quiet = TRUE)))
}
# The APPARENT AUC - fixed panel, fitted and scored on the same samples, no
# resampling at all. This is the unambiguously optimistic number and is the
# correct upper anchor for a selection-bias comparison.
apparent_auc <- function(sex_code, panel_genes) {
  cols <- meta$sample[meta$sex == sex_code]
  g <- unique(panel_genes[panel_genes %in% rownames(expr)])
  if (length(g) < 2) return(NA_real_)
  X <- t(expr[g, cols, drop = FALSE]); colnames(X) <- make.names(g)
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))
  Z <- scale(X)
  fit <- suppressWarnings(glm(y ~ ., data = data.frame(y, Z, check.names = FALSE),
                              family = binomial))
  as.numeric(auc(roc(y, as.numeric(predict(fit, type = "response")),
                     levels = c("HC", "RA"), direction = "<", quiet = TRUE)))
}
ml_for_flat <- readRDS(file.path(procN, "ml_features.rds"))
flatF <- flat_cv_auc("F", ml_for_flat$female$consensus)
flatM <- flat_cv_auc("M", ml_for_flat$male$consensus)
appF  <- apparent_auc("F", ml_for_flat$female$consensus)
appM  <- apparent_auc("M", ml_for_flat$male$consensus)

cat(sprintf("apparent (no resampling)   : female %.3f | male %.3f\n", appF, appM))
cat(sprintf("flat CV (selection once)   : female %.3f | male %.3f\n", flatF, flatM))
cat(sprintf("nested CV (selection in-fold): female %.3f | male %.3f\n", F$auc, M$auc))
cat(sprintf("optimism (apparent - nested): female %+.3f | male %+.3f\n",
            appF - F$auc, appM - M$auc))

# ---------------------------------------------------------------------------
# A NOTE THE READER NEEDS, BECAUSE THE FLAT COLUMN DOES NOT BEHAVE AS EXPECTED.
# flat CV comes out BELOW nested CV here (female 0.801 vs 0.816; male 0.821 vs
# 0.896), which looks like negative selection bias and is not. Two effects run in
# opposite directions:
#   (+) flat CV IS inflated by having chosen the panel on all the data;
#   (-) nested CV pools out-of-fold probabilities across 5 (female) or 10 (male)
#       repeats in which EVERY FOLD SELECTS ITS OWN PANEL. Averaging predictions
#       over many different small panels is an ensemble, and ensembling raises
#       AUC. The pooled nested estimate therefore carries an ensemble bonus that
#       the fixed-panel flat estimate does not.
# The two are consequently NOT a clean bias decomposition, and the difference
# between them must not be reported as "the selection bias". The honest
# optimism estimate is APPARENT minus NESTED, which is positive in both sexes and
# is what the `optimism` column reports.
# ---------------------------------------------------------------------------

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
  procedure = "3-selector consensus re-derived in-fold (LASSO n RF n SVM-RFE)",
  apparent_AUC = c(round(appF, 3), round(appM, 3)),          # no resampling at all
  flat_CV_AUC  = c(round(flatF, 3), round(flatM, 3)),        # computed above, not literal
  nested_CV_AUC = c(round(F$auc, 3), round(M$auc, 3)),
  nested_CI = c(sprintf("%.3f-%.3f", F$ci[1], F$ci[3]), sprintf("%.3f-%.3f", M$ci[1], M$ci[3])),
  optimism_apparent_minus_nested = c(round(appF - F$auc, 3), round(appM - M$auc, 3)),
  flat_minus_nested = c(round(flatF - F$auc, 3), round(flatM - M$auc, 3)),
  flat_column_caveat = paste("flat CV is NOT a clean bias estimate: nested pooling",
                             "across repeats with fold-specific panels is an ensemble",
                             "and raises AUC. Use the optimism column."),
  per_repeat_mean = c(round(mean(F$rep_auc), 3), round(mean(M$rep_auc), 3)),
  per_repeat_sd   = c(round(sd(F$rep_auc), 3), round(sd(M$rep_auc), 3)),
  authoritative = "NO - see NESTED_CV_AUTHORITATIVE.csv (16d) for the reconciled figure")
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
