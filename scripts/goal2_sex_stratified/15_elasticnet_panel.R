#!/usr/bin/env Rscript
# =============================================================================
# 18d_mr_elasticnet_panel.R
# -----------------------------------------------------------------------------
# LEGITIMATE attempt to raise the (honest) diagnostic AUC by giving the model
# more of the real signal that the strict 3-method consensus threw away.
#
# Rationale: the consensus panel is the intersection LASSO n RF n SVM-RFE, which
# yields 7 genes in females and 4 in males. Here we instead let an ELASTIC-NET
# logistic model use the FULL within-sex MR-screened gene set (74 female /
# 55 male) and pick the
# [header corrected 2026-07-27: previously described a superseded run in which
#  the female consensus collapsed to 2 genes (BNIP2, NMI) over a 14/40-gene MR
#  set. Current panels are 7 female (AZI2, CLSTN1, ESYT1, NMI, PCSK7, UBASH3A,
#  WDR46) and 4 male (FNDC3A, GABBR1, RPN2, SSRP1). See CODE_WALKTHROUGH.md.]
# sparsity itself (alpha + lambda tuned by inner CV). Elastic-net keeps
# correlated-but-informative MR genes that a hard intersection discards.
#
# Everything is validated the SAME honest way as 18b/18c:
#   * Candidate universe = the EXTERNAL MR causal genes (fixed; not chosen from
#     this expression matrix), so keeping them fixed across folds is not leakage.
#   * NESTED CV (leakage-free): inside every outer-train fold we tune alpha over
#     a grid and lambda by inner 5-fold CV, fit, and predict the untouched
#     outer-test fold. Female 10-fold x5, Male 5-fold x10 (as 18c).
#   * LOCKED model: one elastic-net fit on the FULL training set, applied ONCE
#     to the internal holdout and external blood. Per-dataset z-scoring (gene-
#     wise, unsupervised) for cross-platform transfer; direction fixed ("<");
#     missing external gene -> z=0 (mean). DeLong CI, or stratified BOOTSTRAP CI
#     when n<20 (male internal/blood), matching hardened 18b.
#
# This does NOT game the estimate: no auto-direction, no scoring train with the
# model that saw it, no seed-hunting. Whatever AUC comes out is reported, and it
# is placed side by side with the consensus-panel numbers from 18b/18c.
#
# Outputs:
#   results/tables/mr_elasticnet_summary.csv     (elastic-net vs consensus)
#   results/tables/mr_elasticnet_coefs_{female,male}.csv  (locked-model coefs)
#   data/processed/mr_elasticnet_objects.rds     (ROC coords for a figure)
# =============================================================================
suppressMessages({library(glmnet); library(pROC); library(caret)
                  library(Biobase); library(data.table)})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(procN,showWarnings=FALSE,recursive=TRUE); dir.create(tab,showWarnings=FALSE,recursive=TRUE)
ALPHAS <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

## ---- data -------------------------------------------------------------------
o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr
meta <- as.data.table(o$meta)
mrF <- fread(file.path(tab, "FS_input_female.csv"))$gene
mrM <- fread(file.path(tab, "FS_input_male.csv"))$gene

load_internal <- function() {
  h <- readRDS(file.path(proc, "internal_val_holdout_processed.rds"))
  list(expr = h$expr, group = factor(h$meta$group, levels = c("HC","RA")),
       sex = h$meta$sex, label = "Internal test")
}
load_blood <- function() {
  e <- readRDS("data/raw/GSE15573_raw.rds"); if (is.list(e)) e <- e[[1]]
  x <- exprs(e); if (max(x, na.rm = TRUE) > 50) x <- log2(x + 1)
  sym <- fData(e)[["Gene symbol"]]; keep <- !is.na(sym) & sym != ""
  x <- x[keep, ]; sym <- sym[keep]; rmean <- rowMeans(x)
  best <- tapply(seq_along(sym), sym, function(ix) ix[which.max(rmean[ix])])
  xg <- x[unlist(best), ]; rownames(xg) <- names(best); p <- pData(e)
  grp <- ifelse(grepl("Rheumatoid|RA", p[["status:ch1"]], ignore.case = TRUE), "RA","HC")
  sex <- ifelse(grepl("Female", p[["gender:ch1"]], ignore.case = TRUE), "F","M")
  list(expr = xg, group = factor(grp, levels = c("HC","RA")), sex = sex,
       label = "External blood")
}
internal <- load_internal(); blood <- load_blood()

zrows <- function(M) t(apply(M, 1, function(v) { s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s }))
auc_ci <- function(r) {                     # DeLong, or bootstrap when n<20
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < 20) { set.seed(GLOBAL_SEED)
    suppressWarnings(as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000)))
  } else as.numeric(ci.auc(r))
  c(auc = as.numeric(auc(r)), ci[c(1, 3)])
}

## ---- tuned elastic-net: pick (alpha,lambda) by inner CV, return fitted cv ----
fit_enet <- function(X, y, nfolds = 5) {
  best <- NULL; bestcvm <- Inf
  for (a in ALPHAS) {
    cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = a,
                             nfolds = nfolds, standardize = TRUE),
                   error = function(e) NULL)
    if (is.null(cv)) next
    if (min(cv$cvm) < bestcvm) { bestcvm <- min(cv$cvm); best <- list(cv = cv, alpha = a) }
  }
  best
}

## ---- nested CV for one sex (leakage-free) -----------------------------------
run_nested <- function(sex_code, genes, kfold, repeats) {
  sexlab <- if (sex_code == "F") "Female" else "Male"
  cols <- meta$sample[meta$sex == sex_code]
  genes <- genes[genes %in% rownames(expr)]
  X0 <- t(expr[genes, cols, drop = FALSE]); colnames(X0) <- make.names(genes)
  y  <- factor(meta$group[match(cols, meta$sample)], levels = c("HC","RA"))
  preds <- list(); rep_auc <- numeric(0); nz <- integer(0); nfits <- 0L
  for (rp in seq_len(repeats)) {
    set.seed(2000 + rp)
    folds <- createFolds(y, k = kfold, returnTrain = FALSE); rp_pred <- data.table()
    for (fi in seq_along(folds)) {
      te <- folds[[fi]]; tr <- setdiff(seq_along(y), te)
      if (length(unique(y[tr])) < 2) next
      b <- fit_enet(X0[tr, , drop = FALSE], y[tr]); if (is.null(b)) next
      nfits <- nfits + 1L
      p  <- as.numeric(predict(b$cv, newx = X0[te, , drop = FALSE],
                               s = "lambda.min", type = "response"))
      co <- as.numeric(coef(b$cv, s = "lambda.min"))[-1]; nz <- c(nz, sum(co != 0))
      rp_pred <- rbind(rp_pred, data.table(sample = rownames(X0)[te], prob = p, obs = y[te]))
    }
    rr <- roc(rp_pred$obs, rp_pred$prob, levels = c("HC","RA"), direction = "<", quiet = TRUE)
    rep_auc <- c(rep_auc, as.numeric(auc(rr))); preds[[rp]] <- rp_pred
  }
  allp <- rbindlist(preds); agg <- allp[, .(prob = mean(prob), obs = obs[1]), by = sample]
  ro <- roc(agg$obs, agg$prob, levels = c("HC","RA"), direction = "<", quiet = TRUE)
  ci <- auc_ci(ro)
  cat(sprintf("[%s] elastic-net NESTED AUC = %.3f (95%% CI %.3f-%.3f) | per-repeat %.3f +/- %.3f | median %d/%d genes\n",
              sexlab, ci[1], ci[2], ci[3], mean(rep_auc), sd(rep_auc),
              as.integer(median(nz)), length(genes)))
  list(sexlab = sexlab, auc = ci[1], ci = ci, rep_auc = rep_auc, roc = ro,
       nz_median = median(nz), ngene = length(genes))
}

## ---- locked model on full train -> internal + blood -------------------------
run_locked <- function(sex_code, genes) {
  sexlab <- if (sex_code == "F") "Female" else "Male"
  genes <- genes[genes %in% rownames(expr)]
  tcols <- meta$sample[meta$sex == sex_code]
  ytr <- factor(meta$group[match(tcols, meta$sample)], levels = c("HC","RA"))
  Ztr <- t(zrows(expr[genes, tcols, drop = FALSE])); colnames(Ztr) <- make.names(genes)
  b <- fit_enet(Ztr, ytr)                       # z-scored train, standardize handled
  co <- as.numeric(coef(b$cv, s = "lambda.min"))
  coefs <- data.table(sex = sexlab, term = c("(Intercept)", make.names(genes)),
                      coef = round(co, 4), alpha = b$alpha)[coef != 0 | term == "(Intercept)"]

  score_ds <- function(ds) {
    sp <- which(ds$sex == sex_code); present <- genes[genes %in% rownames(ds$expr)]
    y <- factor(ds$group[sp], levels = c("HC","RA"))
    if (length(sp) < 3 || length(unique(y)) < 2) return(NULL)
    Z <- zrows(ds$expr[present, sp, drop = FALSE]); Zdf <- as.data.frame(t(Z))
    colnames(Zdf) <- make.names(present)
    for (g in setdiff(make.names(genes), colnames(Zdf))) Zdf[[g]] <- 0
    Zmat <- as.matrix(Zdf[, make.names(genes), drop = FALSE])
    p <- as.numeric(predict(b$cv, newx = Zmat, s = "lambda.min", type = "response"))
    r <- roc(y, p, levels = c("HC","RA"), direction = "<", quiet = TRUE)
    list(label = ds$label, ci = auc_ci(r), roc = r, n = length(y),
         miss = setdiff(genes, present))
  }
  list(coefs = coefs, internal = score_ds(internal), blood = score_ds(blood), alpha = b$alpha)
}

cat("Elastic-net over FULL MR gene set (nested CV + locked external); a few minutes...\n\n")
nF <- run_nested("F", mrF, 10, 5);  nM <- run_nested("M", mrM, 5, 10)
lF <- run_locked("F", mrF);         lM <- run_locked("M", mrM)

fwrite(rbindlist(list(lF$coefs, lM$coefs)),
       file.path(tab, "mr_elasticnet_coefs_bysex.csv"))

## ---- side-by-side summary vs consensus (from 18b/18c) -----------------------
cons <- fread(file.path(tab, "mr_roc_panel_auc.csv"))     # consensus panel numbers
getc <- function(sx, ds) { r <- cons[sex == sx & dataset == ds]
  if (nrow(r)) sprintf("%.3f (%.3f-%.3f)", r$AUC, r$AUC_lo, r$AUC_hi) else NA }
fmt <- function(ci) sprintf("%.3f (%.3f-%.3f)", ci[1], ci[2], ci[3])
li <- function(x) if (is.null(x)) NA_character_ else fmt(x$ci)

summ <- data.table(
  sex = c("Female","Male"),
  consensus_nestedCV = c("0.767 (0.689-0.844)","0.952 (0.892-1.000)"),  # 18c
  enet_nestedCV      = c(fmt(nF$ci), fmt(nM$ci)),
  enet_genes_median  = c(sprintf("%d/%d", nF$nz_median, nF$ngene),
                         sprintf("%d/%d", nM$nz_median, nM$ngene)),
  consensus_internal = c(getc("Female","Internal test"), getc("Male","Internal test")),
  enet_internal      = c(li(lF$internal), li(lM$internal)),
  consensus_blood    = c(getc("Female","External blood"), getc("Male","External blood")),
  enet_blood         = c(li(lF$blood), li(lM$blood)))
fwrite(summ, file.path(tab, "mr_elasticnet_summary.csv"))

roc_coords <- function(res, ds) if (is.null(res)) NULL else
  data.table(sex = res$sexlab, dataset = ds, sens = res$roc$sensitivities,
             spec = res$roc$specificities, auc = res$auc)
saveRDS(list(nested = list(F = nF, M = nM), locked = list(F = lF, M = lM),
             summary = summ,
             roc = list(nestedF = data.table(sex="Female", dataset="Train (nested CV, enet)",
                          sens = nF$roc$sensitivities, spec = nF$roc$specificities, auc = nF$auc),
                        nestedM = data.table(sex="Male", dataset="Train (nested CV, enet)",
                          sens = nM$roc$sensitivities, spec = nM$roc$specificities, auc = nM$auc))),
        file.path(procN, "mr_elasticnet_objects.rds"))

cat("\n================  ELASTIC-NET (full MR set) vs CONSENSUS  ================\n")
print(summ, width = 200)
cat("\nFemale locked-model non-zero coefficients:\n"); print(lF$coefs)
cat("\nSaved mr_elasticnet_summary.csv, mr_elasticnet_coefs_bysex.csv, mr_elasticnet_objects.rds\n")
