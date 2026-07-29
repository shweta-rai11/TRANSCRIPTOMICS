#!/usr/bin/env Rscript
# =============================================================================
# 16b_final_panel_noMHC.R  —  Head-to-head evaluation: primary vs MHC-free panel
#
# WHY THIS EXISTS
#   12b rebuilt the panels from the MHC-free candidate set. This script asks the
#   question that decides whether the rebuild costs anything: **how much
#   diagnostic performance was the MHC actually buying?**
#
#   That question has to be asked at two levels, because they answer different
#   things and only one of them is honest about selection bias:
#
#   (1) PROCEDURE level - NESTED CV. The whole pipeline (three selectors +
#       3-method consensus + logistic model) is re-run inside every outer
#       training fold, once starting from the MHC-retained candidate set and
#       once from the MHC-free set. The held-out fold is never seen during
#       selection. This is the leakage-free comparison and it is the one to
#       report. It compares two PROCEDURES, not two fixed gene lists.
#
#   (2) PANEL level - LOCKED MODEL. The fixed consensus panels from 12 and 12b
#       are trained once on the full training set and applied unchanged to the
#       internal hold-out and external blood. This is what a clinician would
#       actually deploy. Its training-set numbers are optimistic by construction
#       (the panel was chosen using those samples) and are reported as apparent,
#       never as validation.
#
# WHAT A NEGATIVE RESULT WOULD LOOK LIKE, AND WHY IT WOULD STILL BE FINE
#   If the MHC-free panel performs WORSE, that is not a failure of this analysis
#   - it quantifies how much of the original performance was attributable to
#   genes whose MR support came from LD with HLA-DRB1. A smaller, honest AUC
#   from a defensible gene set is worth more in a thesis than a larger one that
#   cannot survive its own sensitivity analysis. The result is reported either
#   way and no panel is declared the winner on the basis of a point estimate
#   alone; the DeLong test on paired predictions decides.
#
#   in : data/processed/combined_train.rds
#        data/processed/internal_val_holdout_processed.rds
#        data/raw/GSE15573_raw.rds
#        results/tables/FS_input_{female,male}{,_noMHC}.csv
#        data/processed/new/ml_features{,_noMHC}.rds
#   out: results/tables/PANEL_primary_vs_noMHC_performance.csv
#        results/tables/PANEL_primary_vs_noMHC_nestedcv.csv
#        results/tables/PANEL_primary_vs_noMHC_delong.csv
#        data/processed/new/panel_noMHC_objects.rds
#
#   Ambroise C, McLachlan GJ. PNAS 2002;99:6562-6566.    (selection bias)
#   DeLong ER, et al. Biometrics 1988;44:837-845.        (correlated ROC test)
#   Carpenter J, Bithell J. Stat Med 2000;19:1141-1164.  (bootstrap CIs)
# =============================================================================
suppressMessages({
  library(glmnet); library(randomForest); library(e1071)
  library(pROC); library(caret); library(Biobase); library(data.table)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
set.seed(GLOBAL_SEED)

proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
SMALL_N <- 20

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# =============================================================================
# STEP 1 — DATA, CANDIDATE SETS AND PANELS
# =============================================================================
hdr("STEP 1  INPUTS")
o <- readRDS(file.path(proc, "combined_train.rds"))
expr <- o$expr; meta <- as.data.table(o$meta)

cand <- list(
  primary = list(F = fread(file.path(tab, "FS_input_female.csv"))$gene,
                 M = fread(file.path(tab, "FS_input_male.csv"))$gene),
  noMHC   = list(F = fread(file.path(tab, "FS_input_female_noMHC.csv"))$gene,
                 M = fread(file.path(tab, "FS_input_male_noMHC.csv"))$gene))

mlP <- readRDS(file.path(procN, "ml_features.rds"))
mlN <- readRDS(file.path(procN, "ml_features_noMHC.rds"))
panel <- list(
  primary = list(F = mlP$female$consensus, M = mlP$male$consensus),
  noMHC   = list(F = mlN$female$consensus, M = mlN$male$consensus))

for (v in c("primary", "noMHC"))
  for (sx in c("F", "M"))
    say("%-8s %s : %2d candidates -> %d-gene panel (%s)", v, sx,
        length(cand[[v]][[sx]]), length(panel[[v]][[sx]]),
        paste(panel[[v]][[sx]], collapse = ", "))

load_internal <- function() {
  h <- readRDS(file.path(proc, "internal_val_holdout_processed.rds"))
  list(expr = h$expr, group = factor(h$meta$group, levels = c("HC", "RA")),
       sex = h$meta$sex, label = "Internal test")
}
load_blood <- function() {
  e <- readRDS("data/raw/GSE15573_raw.rds"); if (is.list(e)) e <- e[[1]]
  x <- exprs(e); if (max(x, na.rm = TRUE) > 50) x <- log2(x + 1)
  sym <- fData(e)[["Gene symbol"]]; keep <- !is.na(sym) & sym != ""
  x <- x[keep, ]; sym <- sym[keep]
  rmean <- rowMeans(x)
  best <- tapply(seq_along(sym), sym, function(ix) ix[which.max(rmean[ix])])
  xg <- x[unlist(best), ]; rownames(xg) <- names(best); p <- pData(e)
  grp <- ifelse(grepl("Rheumatoid|RA", p[["status:ch1"]], ignore.case = TRUE), "RA", "HC")
  sex <- ifelse(grepl("Female", p[["gender:ch1"]], ignore.case = TRUE), "F", "M")
  list(expr = xg, group = factor(grp, levels = c("HC", "RA")), sex = sex,
       label = "External blood")
}
internal <- load_internal(); blood <- load_blood()

# =============================================================================
# STEP 2 — HELPERS
# =============================================================================
zrows <- function(M) t(apply(M, 1, function(v) {
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s }))

auc_ci <- function(r) {
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < SMALL_N) { set.seed(GLOBAL_SEED)
    suppressWarnings(as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000)))
  } else suppressWarnings(as.numeric(ci.auc(r)))
  c(auc = as.numeric(auc(r)), lo = ci[1], hi = ci[3])
}
fmt <- function(a, n) {
  s <- sprintf("%.3f (%.3f-%.3f) [n=%d]", a[1], a[2], a[3], n)
  if (!is.na(a[1]) && a[1] >= 0.999) s <- paste0(s, " SEPARATION")
  s
}

# in-fold three-selector consensus, identical in spirit to 14_nested_cv.R
svm_rank <- function(X, y, cost = 1) {
  feats <- colnames(X); ranking <- character(0)
  while (length(feats) > 1) {
    m <- svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = TRUE, cost = cost)
    w2 <- ((t(m$coefs) %*% m$SV)[1, ])^2
    d <- names(sort(w2))[1]; ranking <- c(d, ranking); feats <- setdiff(feats, d)
  }
  c(feats, ranking)
}
svm_size <- function(X, y, rank, cost = 1) {
  err <- sapply(seq_along(rank), function(k)
    1 - svm(X[, rank[1:k], drop = FALSE], y, kernel = "linear", scale = TRUE,
            cost = cost, cross = 5)$tot.accuracy / 100)
  rank[1:which.min(err)]
}
select_infold <- function(X, y) {
  lasso <- character(0)
  cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5),
                 error = function(e) NULL)
  if (!is.null(cv)) { co <- coef(cv, s = "lambda.min")[-1, 1]; lasso <- colnames(X)[co != 0] }
  rf <- randomForest(X, y, ntree = 500)
  gini <- rf$importance[, "MeanDecreaseGini"]; rfs <- names(gini)[gini > mean(gini)]
  svs <- svm_size(X, y, svm_rank(X, y))
  cons <- Reduce(intersect, list(lasso, rfs, svs))
  uni  <- Reduce(union, list(lasso, rfs, svs))
  if (length(cons) >= 2) cons else if (length(uni) >= 1) uni
  else if (length(lasso) >= 1) lasso else colnames(X)
}

# =============================================================================
# STEP 3 — NESTED CV: PROCEDURE-LEVEL COMPARISON (the honest one)
# =============================================================================
hdr("STEP 3  NESTED CV - PROCEDURE-LEVEL COMPARISON")

run_nested <- function(genes, sx, kfold, repeats, label) {
  genes <- unique(genes[genes %in% rownames(expr)])
  cols <- meta$sample[meta$sex == sx]
  X0 <- t(expr[genes, cols, drop = FALSE]); colnames(X0) <- make.names(genes)
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))
  preds <- list(); nsel <- integer(0)
  for (rp in seq_len(repeats)) {
    set.seed(1000 + rp)
    folds <- createFolds(y, k = kfold, returnTrain = FALSE)
    rp_pred <- data.table()
    for (fi in seq_along(folds)) {
      te <- folds[[fi]]; tr <- setdiff(seq_along(y), te)
      if (length(unique(y[tr])) < 2) next
      Xtr <- X0[tr, , drop = FALSE]; Xte <- X0[te, , drop = FALSE]
      P <- select_infold(Xtr, y[tr]); nsel <- c(nsel, length(P))
      mu <- colMeans(Xtr[, P, drop = FALSE])
      sg <- apply(Xtr[, P, drop = FALSE], 2, sd); sg[sg == 0 | is.na(sg)] <- 1
      Ztr <- scale(Xtr[, P, drop = FALSE], center = mu, scale = sg)
      Zte <- scale(Xte[, P, drop = FALSE], center = mu, scale = sg)
      fit <- suppressWarnings(glm(y[tr] ~ ., data = data.frame(y = y[tr], Ztr,
                                                              check.names = FALSE),
                                  family = binomial))
      p <- as.numeric(predict(fit, newdata = data.frame(Zte, check.names = FALSE),
                              type = "response"))
      rp_pred <- rbind(rp_pred, data.table(sample = rownames(Xte), prob = p, obs = y[te]))
    }
    preds[[rp]] <- rp_pred
  }
  agg <- rbindlist(preds)[, .(prob = mean(prob), obs = obs[1]), by = sample]
  ro <- roc(agg$obs, agg$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  a <- auc_ci(ro)
  say("  [%-8s %s] nested-CV AUC = %s | median genes used %d",
      label, sx, fmt(a, length(y)), as.integer(median(nsel)))
  list(roc = ro, auc = a, agg = agg, n = length(y), med = as.integer(median(nsel)))
}

nested <- list()
for (v in c("primary", "noMHC")) {
  nested[[v]] <- list(
    F = run_nested(cand[[v]]$F, "F", 10, 5, v),
    M = run_nested(cand[[v]]$M, "M", 5, 10, v))
}

nest_tab <- rbindlist(lapply(c("F", "M"), function(sx) {
  sexlab <- if (sx == "F") "Female" else "Male"
  data.table(
    sex = sexlab,
    evidence_tier = if (sx == "M") "EXPLORATORY (underpowered)" else "primary",
    n = nested$primary[[sx]]$n,
    n_candidates_primary = length(cand$primary[[sx]]),
    n_candidates_noMHC   = length(cand$noMHC[[sx]]),
    nested_AUC_primary = fmt(nested$primary[[sx]]$auc, nested$primary[[sx]]$n),
    nested_AUC_noMHC   = fmt(nested$noMHC[[sx]]$auc,   nested$noMHC[[sx]]$n),
    delta_AUC = round(nested$noMHC[[sx]]$auc[1] - nested$primary[[sx]]$auc[1], 3),
    med_genes_primary = nested$primary[[sx]]$med,
    med_genes_noMHC   = nested$noMHC[[sx]]$med)
}))
fwrite(nest_tab, file.path(tab, "PANEL_primary_vs_noMHC_nestedcv.csv"))
print(nest_tab[, .(sex, n, nested_AUC_primary, nested_AUC_noMHC, delta_AUC)])

# DeLong on the paired pooled out-of-fold predictions
delong_rows <- lapply(c("F", "M"), function(sx) {
  a <- nested$primary[[sx]]$agg; b <- nested$noMHC[[sx]]$agg
  m <- merge(a, b, by = "sample", suffixes = c("_pri", "_no"))
  r1 <- roc(m$obs_pri, m$prob_pri, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  r2 <- roc(m$obs_pri, m$prob_no,  levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  p <- tryCatch(suppressWarnings(roc.test(r1, r2, method = "delong", paired = TRUE)$p.value),
                error = function(e) NA_real_)
  data.table(sex = if (sx == "F") "Female" else "Male",
             comparison = "nested CV: primary vs MHC-free",
             AUC_primary = round(as.numeric(auc(r1)), 3),
             AUC_noMHC = round(as.numeric(auc(r2)), 3),
             delong_p = signif(p, 3),
             verdict = fifelse(is.na(p), "indeterminate",
                        fifelse(p >= 0.05,
                                "NO significant difference - MHC removal costs nothing detectable",
                                fifelse(as.numeric(auc(r2)) < as.numeric(auc(r1)),
                                        "MHC-free panel significantly WORSE",
                                        "MHC-free panel significantly BETTER"))))
})
delong_tab <- rbindlist(delong_rows)
fwrite(delong_tab, file.path(tab, "PANEL_primary_vs_noMHC_delong.csv"))
say("")
print(delong_tab)

# =============================================================================
# STEP 4 — LOCKED PANELS APPLIED TO HELD-OUT DATA
# =============================================================================
hdr("STEP 4  LOCKED PANELS ON INTERNAL AND EXTERNAL DATA")

eval_locked <- function(genes, sx, variant) {
  genes <- unique(genes[genes %in% rownames(expr)])
  if (length(genes) < 2) return(NULL)
  tcols <- meta$sample[meta$sex == sx]
  ytr <- factor(meta$group[match(tcols, meta$sample)], levels = c("HC", "RA"))
  Ztr <- as.data.frame(t(zrows(expr[genes, tcols, drop = FALSE])))
  colnames(Ztr) <- make.names(genes)
  fit <- suppressWarnings(glm(ytr ~ ., data = cbind(ytr = ytr, Ztr), family = binomial))

  score <- function(ds) {
    sp <- which(ds$sex == sx); present <- genes[genes %in% rownames(ds$expr)]
    y <- factor(ds$group[sp], levels = c("HC", "RA"))
    if (length(sp) < 3 || length(unique(y)) < 2 || !length(present)) return(NULL)
    Z <- as.data.frame(t(zrows(ds$expr[present, sp, drop = FALSE])))
    colnames(Z) <- make.names(present)
    for (g in setdiff(make.names(genes), colnames(Z))) Z[[g]] <- 0
    Z <- Z[, make.names(genes), drop = FALSE]
    p <- as.numeric(predict(fit, newdata = Z, type = "response"))
    r <- roc(y, p, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
    list(a = auc_ci(r), n = length(y), miss = setdiff(genes, present))
  }
  ap <- roc(ytr, as.numeric(predict(fit, type = "response")),
            levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  list(train_apparent = list(a = auc_ci(ap), n = length(ytr), miss = character(0)),
       internal = score(internal), blood = score(blood),
       n_genes = length(genes), missing_ext = setdiff(genes, rownames(blood$expr)))
}

perf <- list()
rows <- list()
for (v in c("primary", "noMHC")) for (sx in c("F", "M")) {
  r <- eval_locked(panel[[v]][[sx]], sx, v)
  perf[[paste(v, sx)]] <- r
  if (is.null(r)) next
  sexlab <- if (sx == "F") "Female" else "Male"
  for (dn in c("train_apparent", "internal", "blood")) {
    x <- r[[dn]]; if (is.null(x)) next
    rows[[length(rows) + 1]] <- data.table(
      sex = sexlab, panel = v, n_panel_genes = r$n_genes,
      dataset = c(train_apparent = "Train (apparent)", internal = "Internal test",
                  blood = "External blood")[[dn]],
      n = x$n, AUC = round(x$a[1], 3), AUC_lo = round(x$a[2], 3), AUC_hi = round(x$a[3], 3),
      reported = fmt(x$a, x$n),
      genes_missing_in_dataset = paste(x$miss, collapse = ";"),
      evidence_tier = if (sx == "M" || x$n < SMALL_N) "EXPLORATORY (underpowered)" else "primary",
      note = if (dn == "train_apparent")
        "APPARENT - panel was selected using these samples; not a validation estimate" else "")
  }
}
perf_tab <- rbindlist(rows)
fwrite(perf_tab, file.path(tab, "PANEL_primary_vs_noMHC_performance.csv"))
print(perf_tab[, .(sex, panel, dataset, n, reported, evidence_tier)])

saveRDS(list(nested = nested, nested_tab = nest_tab, delong = delong_tab,
             performance = perf_tab, panels = panel, candidates = cand),
        file.path(procN, "panel_noMHC_objects.rds"))

# =============================================================================
# STEP 5 — VERDICT
# =============================================================================
hdr("STEP 5  VERDICT")
for (i in seq_len(nrow(delong_tab))) {
  d <- delong_tab[i]
  say("%s: nested-CV AUC %.3f (primary) vs %.3f (MHC-free), DeLong p = %s",
      d$sex, d$AUC_primary, d$AUC_noMHC, format.pval(d$delong_p, digits = 3))
  say("   -> %s", d$verdict)
}
say("")
say("READING RULE FOR THE THESIS")
say("  Report the NESTED-CV comparison, not the locked-panel training numbers.")
say("  The nested comparison is between two PROCEDURES and is leakage-free; the")
say("  locked 'Train (apparent)' rows are optimistic by construction and are")
say("  labelled as such in the output.")
say("  If the DeLong test shows no significant difference, the correct sentence")
say("  is: 'excluding the MHC did not measurably reduce diagnostic performance,")
say("  so the panel does not depend on HLA linkage disequilibrium'. That is a")
say("  STRENGTH of the MHC-free panel, and it is the version to carry forward.")
say("  Every Male row remains EXPLORATORY at n = 38/13/9 regardless of AUC.")
say("")
say("Wrote PANEL_primary_vs_noMHC_{performance,nestedcv,delong}.csv and")
say("panel_noMHC_objects.rds")
cat("\nDONE\n")
