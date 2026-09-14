#!/usr/bin/env Rscript
# =============================================================================
# 38_testing_ml_algorithms_sametissue_crosstissue.R
# -----------------------------------------------------------------------------
# TESTING PHASE for the five-algorithm diagnostic-model comparison trained in
# 37_. Each of the five locked models (logistic regression, SVM-RBF, k-NN,
# random forest, ANN), fitted per sex on the blood training cohort, is applied
# ONCE to two independent test settings:
#
#   SAME-TISSUE  (blood -> blood): data/processed/internal_val_holdout_processed.rds,
#     the sealed 30% internal-validation holdout already used throughout this
#     chapter (frozen quantile + ComBat parameters; no holdout sample
#     influenced training). Female n=61 (36 RA/25 HC); Male n=13 (6 RA/7 HC).
#
#   CROSS-TISSUE (blood -> synovium): the synovial RNA-seq log2-CPM matrix
#     already computed and cached by 20_testing_synovium_external.R
#     (data/processed/new/val_synovium.rds, GSE89408, RA vs Normal). Female
#     n=120 (106 RA/14 Normal); Male n=60 (46 RA/14 Normal).
#
# Every model input is z-scored WITHIN the test dataset itself (own mean/sd),
# not by re-using the training mean/sd, for the same reason given in 37_: the
# synovial data is RNA-seq log2-CPM on a different platform and tissue from
# the microarray training data, so a frozen training scale would not transfer,
# and applying the identical unsupervised within-dataset standardisation to
# both test settings keeps same-tissue and cross-tissue comparable.
#
# "TRAIN" CURVE shown alongside each test curve is the pooled out-of-fold
# cross-validated prediction from 37_ (Train, resampled CV), NOT resubstitution
# on the full training set, so it is not a foregone conclusion that train beats
# test.
#
# EVIDENCE TIERS, carried over from 16_/18b_/20_: the male stratum (train
# n=38) is flagged EXPLORATORY throughout this comparison regardless of how
# strong an individual AUC looks, because at this sample size a handful of
# discordant pairs can drive an apparently perfect classifier (see 37_'s log:
# three of five male algorithms already reach CV ROC = 1.000 +/- 0.000 in
# training, a small-n separation artefact, not evidence of superiority).
# AUCs >= 0.999 are flagged SEPARATION in the table so they are never quoted
# as if they were precise performance estimates. 95% CIs use DeLong's method
# when n >= 20 and a stratified bootstrap (2000 resamples, seed 1234)
# otherwise (Carpenter & Bithell, 2000).
#
# Outputs:
#   results/tables/ML_performance_sametissue.csv
#   results/tables/ML_performance_crosstissue.csv
#   data/processed/new/ml_algo_roc.rds   (ROC coordinates for 39_figure_*.R)
# =============================================================================
suppressMessages({
  library(caret); library(data.table); library(pROC)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"

fit <- readRDS(file.path(procN, "ml_algo_models.rds"))
PANEL <- fit$panel

## ---- same-tissue test set: blood internal holdout -------------------------
h <- readRDS(file.path(proc, "internal_val_holdout_processed.rds"))
blood_test <- list(expr = h$expr,
                    group = factor(h$meta$group, levels = c("HC", "RA")),
                    sex = h$meta$sex, label = "Blood test (internal holdout)")

## ---- cross-tissue test set: synovium (cached, see header) -----------------
syn <- readRDS(file.path(procN, "val_synovium.rds"))
syn_group <- factor(as.character(syn$grp), levels = c("Normal", "RA"))
syn_test  <- list(expr = syn$logcpm, group = syn_group, sex = syn$sex,
                   label = "Synovium test (cross-tissue)")

cat(sprintf("Blood test:    n=%d (RA=%d/HC=%d) | F=%d M=%d\n",
            length(blood_test$group), sum(blood_test$group == "RA"), sum(blood_test$group == "HC"),
            sum(blood_test$sex == "F"), sum(blood_test$sex == "M")))
cat(sprintf("Synovium test: n=%d (RA=%d/Normal=%d) | F=%d M=%d\n",
            length(syn_test$group), sum(syn_test$group == "RA"), sum(syn_test$group == "Normal"),
            sum(syn_test$sex == "F"), sum(syn_test$sex == "M")))

## ---- helpers ---------------------------------------------------------------
zrows <- function(M) t(apply(M, 1, function(v) {
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
}))
auc_ci <- function(r) {
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < 20) { set.seed(GLOBAL_SEED)
    suppressWarnings(as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000)))
  } else as.numeric(ci.auc(r))
  c(auc = as.numeric(auc(r)), ci[c(1, 3)])
}
fmt <- function(a, n) {
  s <- sprintf("%.3f (%.3f-%.3f) [n=%d]", a[1], a[2], a[3], n)
  if (!is.na(a[1]) && a[1] >= 0.999) s <- paste0(s, " SEPARATION")
  s
}

## score a fitted caret model on a new (already sex-subset) expression block
score_model <- function(m, genes, colnms, X_genesXsamples, y) {
  present <- genes[genes %in% rownames(X_genesXsamples)]
  Z <- zrows(X_genesXsamples[present, , drop = FALSE])
  Zdf <- as.data.frame(t(Z), check.names = FALSE); colnames(Zdf) <- make.names(present)
  for (g in setdiff(colnms, colnames(Zdf))) Zdf[[g]] <- 0     # missing gene -> dataset mean (0 after z-score)
  Zdf <- Zdf[, colnms, drop = FALSE]
  p <- predict(m, newdata = Zdf, type = "prob")[, "RA"]
  r <- roc(y, p, direction = "<", levels = c("HC", "RA"), quiet = TRUE)
  list(roc = r, ci = auc_ci(r), n = length(y), n_RA = sum(y == "RA"), n_ctrl = sum(y == "HC"),
       missing = setdiff(genes, present))
}

## roc-coordinate row for the figure script
rc <- function(sexlab, algo, dataset, r, tier) {
  data.table(sex = sexlab, algorithm = algo, dataset = dataset,
             sens = r$sensitivities, spec = r$specificities,
             auc = as.numeric(auc(r)), evidence_tier = tier)
}

## -----------------------------------------------------------------------------
run_setting <- function(test_set, setting_label, fixed_group_levels) {
  perf_rows <- list(); roc_rows <- list()
  for (sx in c("F", "M")) {
    sexlab <- if (sx == "F") "Female" else "Male"
    res <- if (sx == "F") fit$female else fit$male
    n_train <- length(res$y)
    tier <- if (n_train < 50) "EXPLORATORY (underpowered)" else "primary"
    ti <- which(test_set$sex == sx)
    if (length(ti) < 3) { cat(sprintf("  [%s/%s] skipped: too few test samples\n", setting_label, sexlab)); next }
    y_te <- factor(test_set$group[ti], levels = fixed_group_levels)
    if (length(unique(y_te)) < 2) { cat(sprintf("  [%s/%s] skipped: single class in test\n", setting_label, sexlab)); next }
    # relabel to HC/RA for scoring (Normal == HC-equivalent "not RA" reference)
    y_te2 <- factor(ifelse(y_te == fixed_group_levels[2], "RA", "HC"), levels = c("HC", "RA"))
    Xte <- test_set$expr[, ti, drop = FALSE]

    for (a in names(res$models)) {
      mm <- res$models[[a]]
      spec_label <- ALGO_LABELS[[a]]
      # "Train (resampled CV)" honest curve, from 37_
      perf_rows[[length(perf_rows) + 1]] <- data.table(
        sex = sexlab, algorithm = spec_label, dataset = "Train (resampled CV)",
        setting = setting_label, evidence_tier = tier,
        n = n_train, n_RA = sum(res$y == "RA"), n_ctrl = sum(res$y == "HC"),
        AUC_CI = fmt(auc_ci(mm$oof), n_train), missing_genes = "")
      roc_rows[[length(roc_rows) + 1]] <- rc(sexlab, spec_label, "Train (resampled CV)", mm$oof, tier)

      te <- score_model(mm$model, mm$genes, mm$colnames, Xte, y_te2)
      perf_rows[[length(perf_rows) + 1]] <- data.table(
        sex = sexlab, algorithm = spec_label, dataset = setting_label,
        setting = setting_label, evidence_tier = tier,
        n = te$n, n_RA = te$n_RA, n_ctrl = te$n_ctrl,
        AUC_CI = fmt(te$ci, te$n), missing_genes = paste(te$missing, collapse = ";"))
      roc_rows[[length(roc_rows) + 1]] <- rc(sexlab, spec_label, setting_label, te$roc, tier)
    }
    cat(sprintf("  [%s/%s] done (train n=%d, test n=%d, tier=%s)\n",
                setting_label, sexlab, n_train, length(ti), tier))
  }
  list(perf = rbindlist(perf_rows), roc = rbindlist(roc_rows))
}

ALGO_LABELS <- list(logreg = "Logistic regression", svm = "SVM (RBF)", knn = "k-NN",
                     rf = "Random forest", ann = "ANN (nnet)")

cat("\n== SAME-TISSUE (blood -> blood) ==\n")
same <- run_setting(blood_test, "Blood test (internal holdout)", c("HC", "RA"))
cat("\n== CROSS-TISSUE (blood -> synovium) ==\n")
cross <- run_setting(syn_test, "Synovium test (cross-tissue)", c("Normal", "RA"))

fwrite(same$perf,  file.path(tab, "ML_performance_sametissue.csv"))
fwrite(cross$perf, file.path(tab, "ML_performance_crosstissue.csv"))

saveRDS(list(sametissue = same, crosstissue = cross, seed = GLOBAL_SEED,
             built = "goal2_sex_stratified/38_testing_ml_algorithms_sametissue_crosstissue.R"),
        file.path(procN, "ml_algo_roc.rds"))

cat("\n====================  SAME-TISSUE PERFORMANCE  ====================\n")
print(same$perf, width = 200)
cat("\n====================  CROSS-TISSUE PERFORMANCE  ====================\n")
print(cross$perf, width = 200)
cat("\nSaved ML_performance_sametissue.csv, ML_performance_crosstissue.csv, ml_algo_roc.rds\n")
