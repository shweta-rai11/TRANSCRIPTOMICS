#!/usr/bin/env Rscript
# =============================================================================
# 10c_nested_cv.R  —  NESTED cross-validation (leakage-free) diagnostic AUC,
# per sex. Fixes the selection-bias that inflated the earlier "flat-CV" AUCs
# (male = 1.000): ALL feature selection happens INSIDE each outer training fold
# and the held-out fold is never seen during selection or tuning.
#
#   Outer loop : repeated stratified k-fold (male 5-fold x 20; female 10-fold x 5)
#   Inside each outer-train fold, on TRAIN ONLY:
#       (1) variance prefilter -> top 2000 genes
#       (2) LASSO (cv.glmnet, inner 10-fold CV tunes lambda) selects features
#       (3) the fitted LASSO is the classifier
#   Predict the untouched outer-test fold -> pool out-of-fold probabilities.
# Reports pooled AUC (+95% CI), per-repeat AUC distribution, and gene-selection
# frequency (stability). Output: results/tables/nested_cv_*.csv + nested_cv.rds
# =============================================================================
suppressMessages({library(glmnet); library(pROC); library(caret); library(data.table); library(limma)})
proc <- "data/processed"; tab <- "results/tables"
o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr; meta <- as.data.table(o$meta)
ml <- readRDS(file.path(proc, "ml_features.rds"))
FC <- 0.5; PV <- 0.05; DEG_CAP <- 300
# WARNING (corrected 2026-07-27). The previous comment here claimed these were
# "the same thresholds as the pipeline (03_dge / 09)". They are NOT, and 03_dge.R
# belongs to the SYNOVIUM chapter, not this one. This chapter's DEG step
# (00_shared/05_dge.R lines 38, 56) uses |log2FC| > 0.1 and BH-FDR < 0.05.
#
# What this script actually benchmarks: limma -> LASSO, at the thresholds above.
# What the reported pipeline actually is: limma -> WGCNA disease modules ->
#   cis-eQTL MR -> LASSO n RF n SVM-RFE consensus.
# WGCNA, MR and the 3-way consensus are ABSENT from the fold loop below.
#
# => These numbers bound the magnitude of feature-selection bias. They are NOT a
#    leakage-free validation of the 7-gene / 4-gene panels, and must not be
#    reported as such. See THESIS_METHODOLOGY_AUDIT.md finding A2.

run_nested <- function(sex, kfold, repeats) {
  cols <- meta$sample[meta$sex == sex]
  E <- expr[, cols]                                      # genes x samples (for in-fold limma)
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))
  ds <- factor(meta$dataset[match(cols, meta$sample)])
  preds <- list(); sel_all <- character(0); rep_auc <- numeric(0)
  for (rp in seq_len(repeats)) {
    set.seed(1000 + rp)
    folds <- createFolds(y, k = kfold, returnTrain = FALSE)
    rp_pred <- data.table()
    for (fi in seq_along(folds)) {
      te <- folds[[fi]]; tr <- setdiff(seq_along(y), te)
      Etr <- E[, tr]; ytr <- y[tr]; dstr <- droplevels(ds[tr])
      # (1) in-fold DEG selection (limma RA vs HC, dataset-adjusted) on TRAIN ONLY
      design <- if (nlevels(dstr) > 1) model.matrix(~ ytr + dstr) else model.matrix(~ ytr)
      fit <- eBayes(lmFit(Etr, design))
      tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")   # coef 2 = ytrRA
      sig <- rownames(tt)[tt$P.Value < PV & abs(tt$logFC) > FC]
      if (length(sig) < 2) next
      if (sex == "M" && length(sig) > DEG_CAP)                         # match male cap
        sig <- rownames(tt[order(tt$adj.P.Val), ])[rownames(tt[order(tt$adj.P.Val), ]) %in% sig][1:DEG_CAP]
      Xtr <- t(Etr[sig, , drop = FALSE])
      # (2) LASSO selection + inner-CV lambda tuning on the in-fold DEGs (TRAIN ONLY)
      cv <- tryCatch(cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = 10),
                     error = function(e) NULL)
      if (is.null(cv)) next
      co <- coef(cv, s = "lambda.min")[-1, 1]; sel <- names(co)[co != 0]
      sel_all <- c(sel_all, sel)
      # (3) predict untouched TEST fold
      p <- as.numeric(predict(cv, t(E[sig, te, drop = FALSE]), s = "lambda.min", type = "response"))
      rp_pred <- rbind(rp_pred, data.table(sample = colnames(E)[te], prob = p,
                                           obs = y[te], rep = rp, nsel = length(sel)))
    }
    # per-repeat AUC (each sample once within a repeat)
    rr <- roc(rp_pred$obs, rp_pred$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
    rep_auc <- c(rep_auc, as.numeric(auc(rr)))
    preds[[rp]] <- rp_pred
  }
  allp <- rbindlist(preds)
  agg <- allp[, .(prob = mean(prob), obs = obs[1]), by = sample]     # avg prob per sample over repeats
  ro <- roc(agg$obs, agg$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(ro))
  cat(sprintf("\n[%s] NESTED-CV pooled AUC = %.3f (95%% CI %.3f-%.3f) | per-repeat %.3f +/- %.3f | median %d feats/fold\n",
      sex, as.numeric(auc(ro)), ci[1], ci[3], mean(rep_auc), sd(rep_auc), round(median(allp$nsel))))
  list(auc = as.numeric(auc(ro)), ci = ci, rep_auc = rep_auc, roc = ro,
       agg = agg, selfreq = sort(table(sel_all), decreasing = TRUE),
       nfits = kfold * repeats, sex = sex)
}

cat("Running nested CV (this recomputes feature selection inside every fold)...\n")
F <- run_nested("F", kfold = 10, repeats = 5)
M <- run_nested("M", kfold = 5,  repeats = 20)

# ---- how often were the ORIGINAL signature genes re-selected? ----------------
stab <- function(res, sig) {
  frv <- as.integer(res$selfreq); names(frv) <- names(res$selfreq)
  cnt <- frv[sig]; cnt[is.na(cnt)] <- 0L
  data.table(gene = sig, reselect_pct = round(100 * as.numeric(cnt) / res$nfits, 1))[order(-reselect_pct)]
}
sf <- stab(F, ml$female$signature); sm <- stab(M, ml$male$signature)
fwrite(sf, file.path(tab, "nested_cv_stability_female.csv"))
fwrite(sm, file.path(tab, "nested_cv_stability_male.csv"))

summ <- data.table(
  sex = c("Female", "Male"),
  flatCV_leaky_AUC = c(0.947, 1.000),                    # earlier best (with leakage)
  nested_AUC = c(round(F$auc, 3), round(M$auc, 3)),
  nested_CI = c(sprintf("%.3f-%.3f", F$ci[1], F$ci[3]), sprintf("%.3f-%.3f", M$ci[1], M$ci[3])),
  per_repeat_mean = c(round(mean(F$rep_auc), 3), round(mean(M$rep_auc), 3)))
fwrite(summ, file.path(tab, "nested_cv_summary.csv"))
cat("\n================ NESTED-CV vs LEAKY FLAT-CV ================\n"); print(summ)
cat("\nMale signature re-selection frequency (nested folds):\n"); print(sm)

saveRDS(list(female = F, male = M, summary = summ, stab_f = sf, stab_m = sm),
        file.path(proc, "nested_cv.rds"))
cat("\nSaved nested_cv.rds + nested_cv_*.csv\n")
