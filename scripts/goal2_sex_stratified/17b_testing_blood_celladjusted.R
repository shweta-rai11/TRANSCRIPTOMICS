#!/usr/bin/env Rscript
# Tests whether the diagnostic panel adds signal beyond leukocyte composition (panel-only vs composition-only vs both vs residualised panel), by sex and dataset, with small-n exploratory flagging.
suppressMessages({
  library(pROC); library(data.table); library(Biobase)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
set.seed(GLOBAL_SEED)

proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)

CFG <- list(comp_pcs = 3, kfold = 10, boot_n = 2000, small_n = 20, min_frac = 0.001)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# STEP 1: load panels, expression and cell fractions
hdr("STEP 1  LOAD PANELS, EXPRESSION AND CELL FRACTIONS")
# Both panels (primary and MHC-free) are benchmarked
ml  <- readRDS(file.path(procN, "ml_features.rds"))
mlN_path <- file.path(procN, "ml_features_noMHC.rds")
mlN <- if (file.exists(mlN_path)) readRDS(mlN_path) else NULL

PANELS <- list(primary = list(F = ml$female$consensus, M = ml$male$consensus))
if (!is.null(mlN)) {
  PANELS$noMHC <- list(F = mlN$female$consensus, M = mlN$male$consensus)
} else {
  say("ml_features_noMHC.rds not found - only the primary panel benchmarked (run 12b).")
}

for (v in names(PANELS)) for (sx in c("F", "M"))
  say("%-8s %s: %d genes (%s)", v, sx, length(PANELS[[v]][[sx]]),
      paste(PANELS[[v]][[sx]], collapse = ", "))

tr <- readRDS(file.path(proc, "combined_train.rds"))
ho <- readRDS(file.path(proc, "internal_val_holdout_processed.rds"))

load_external <- function() {
  p <- "data/raw/GSE15573_raw.rds"
  if (!file.exists(p)) return(NULL)
  e <- readRDS(p); if (is.list(e) && !inherits(e, "ExpressionSet")) e <- e[[1]]
  ex <- exprs(e)
  if (stats::quantile(ex, 0.99, na.rm = TRUE) > 100) ex <- log2(ex + 1)
  fd <- fData(e)
  sc <- grep("^Symbol$|Gene.symbol|GENE_SYMBOL", colnames(fd), value = TRUE)[1]
  if (!is.na(sc)) {
    g <- as.character(fd[[sc]]); keep <- !is.na(g) & g != "" & !grepl("///", g)
    ex <- ex[keep, , drop = FALSE]; g <- g[keep]
    o <- order(g, -rowMeans(ex, na.rm = TRUE)); ex <- ex[o, , drop = FALSE]; g <- g[o]
    ex <- ex[!duplicated(g), , drop = FALSE]; rownames(ex) <- g[!duplicated(g)]
  }
  vm <- fread(file.path(tab, "GSE15573_verified_metadata.csv"))
  vm <- vm[sample %in% colnames(ex)]
  ex <- ex[, vm$sample, drop = FALSE]
  list(expr = ex, meta = data.table(sample = vm$sample, group = vm$group, sex = vm$sex))
}
ext <- load_external()

datasets <- list(
  train    = list(label = "Train",         expr = tr$expr,
                  meta = data.table(sample = colnames(tr$expr),
                                    group = tr$meta$group, sex = tr$meta$sex)),
  internal = list(label = "Internal test", expr = ho$expr,
                  meta = data.table(sample = colnames(ho$expr),
                                    group = ho$meta$group, sex = ho$meta$sex)))
if (!is.null(ext))
  datasets$external <- list(label = "External blood", expr = ext$expr, meta = ext$meta)

# CIBERSORT fractions written by 05c, one file per dataset
frac_file <- c(train = "CELL_fractions_train.csv",
               internal = "CELL_fractions_holdout.csv",
               external = "CELL_fractions_external.csv")
for (dn in names(datasets)) {
  fp <- file.path(tab, frac_file[[dn]])
  if (!file.exists(fp)) stop("Missing ", fp, " - run scripts/00_shared/05c_deconvolution.R first.")
  datasets[[dn]]$frac <- fread(fp)
}
say("datasets: %s", paste(sapply(datasets, function(d)
  sprintf("%s n=%d", d$label, nrow(d$meta))), collapse = " | "))

# STEP 2: composition PCs per dataset (CLR, unsupervised)
hdr("STEP 2  COMPOSITION PCs PER DATASET")
comp_pcs <- function(fr, samples, k) {
  cols <- grep("_CIBERSORT$", names(fr), value = TRUE)
  cols <- cols[!grepl("^P\\.value|^Correlation|^RMSE", cols)]
  M <- as.matrix(fr[match(samples, fr$sample), cols, with = FALSE])
  rownames(M) <- samples
  M <- M[, colMeans(M, na.rm = TRUE) >= CFG$min_frac, drop = FALSE]
  M[!is.finite(M)] <- 0
  pos <- M[M > 0]; M[M <= 0] <- if (length(pos)) min(pos) / 2 else 1e-6
  clr <- log(M) - rowMeans(log(M))
  p <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
  kk <- min(k, ncol(p$x))
  Z <- p$x[, seq_len(kk), drop = FALSE]
  colnames(Z) <- paste0("cPC", seq_len(kk))
  list(Z = Z, var_expl = sum((p$sdev^2 / sum(p$sdev^2))[seq_len(kk)]), n_subsets = ncol(M))
}
for (dn in names(datasets)) {
  d <- datasets[[dn]]
  cp <- comp_pcs(d$frac, d$meta$sample, CFG$comp_pcs)
  datasets[[dn]]$comp <- cp$Z
  say("%-15s %d LM22 subsets -> %d PCs, %.1f%% of composition variance",
      d$label, cp$n_subsets, ncol(cp$Z), 100 * cp$var_expl)
}

# STEP 3: helpers
auc_ci <- function(r) {
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < CFG$small_n) {
    set.seed(GLOBAL_SEED)
    suppressWarnings(as.numeric(ci.auc(r, method = "bootstrap", boot.n = CFG$boot_n)))
  } else suppressWarnings(as.numeric(ci.auc(r)))
  c(auc = as.numeric(auc(r)), lo = ci[1], hi = ci[3],
    ci_method = if (n < CFG$small_n) 1 else 0)
}

# standardise using FOLD-TRAIN moments only
freeze_scale <- function(Xtr, Xte) {
  mu <- colMeans(Xtr); sg <- apply(Xtr, 2, stats::sd)
  sg[!is.finite(sg) | sg == 0] <- 1
  list(tr = scale(Xtr, center = mu, scale = sg),
       te = scale(Xte, center = mu, scale = sg))
}

# residualise gene columns on composition, coefficients from FOLD-TRAIN only
freeze_resid <- function(Gtr, Gte, Ctr, Cte) {
  rtr <- Gtr; rte <- Gte
  for (j in seq_len(ncol(Gtr))) {
    fit <- stats::lm.fit(cbind(1, Ctr), Gtr[, j])
    rtr[, j] <- Gtr[, j] - cbind(1, Ctr) %*% fit$coefficients
    rte[, j] <- Gte[, j] - cbind(1, Cte) %*% fit$coefficients
  }
  list(tr = rtr, te = rte)
}

glm_oof <- function(X, y, k) {
  set.seed(GLOBAL_SEED)
  idx <- unlist(tapply(seq_along(y), y, function(ii) sample(ii)))
  fold <- integer(length(y)); fold[idx] <- rep_len(seq_len(k), length(idx))
  p <- rep(NA_real_, length(y))
  for (f in seq_len(k)) {
    te <- which(fold == f); trn <- which(fold != f)
    if (!length(te) || length(unique(y[trn])) < 2) next
    s <- freeze_scale(X[trn, , drop = FALSE], X[te, , drop = FALSE])
    dtr <- data.frame(y = y[trn], s$tr, check.names = FALSE)
    fit <- suppressWarnings(stats::glm(y ~ ., data = dtr, family = binomial))
    p[te] <- as.numeric(stats::predict(fit, newdata = data.frame(s$te, check.names = FALSE),
                                       type = "response"))
  }
  p
}

# out-of-fold probabilities for the residualised panel (model D)
glm_oof_resid <- function(G, C, y, k) {
  set.seed(GLOBAL_SEED)
  idx <- unlist(tapply(seq_along(y), y, function(ii) sample(ii)))
  fold <- integer(length(y)); fold[idx] <- rep_len(seq_len(k), length(idx))
  p <- rep(NA_real_, length(y))
  for (f in seq_len(k)) {
    te <- which(fold == f); trn <- which(fold != f)
    if (!length(te) || length(unique(y[trn])) < 2) next
    rs <- freeze_resid(G[trn, , drop = FALSE], G[te, , drop = FALSE],
                       C[trn, , drop = FALSE], C[te, , drop = FALSE])
    s <- freeze_scale(rs$tr, rs$te)
    dtr <- data.frame(y = y[trn], s$tr, check.names = FALSE)
    fit <- suppressWarnings(stats::glm(y ~ ., data = dtr, family = binomial))
    p[te] <- as.numeric(stats::predict(fit, newdata = data.frame(s$te, check.names = FALSE),
                                       type = "response"))
  }
  p
}

safe_roc <- function(y, p) {
  ok <- !is.na(p)
  if (sum(ok) < 4 || length(unique(y[ok])) < 2) return(NULL)
  suppressWarnings(roc(y[ok], p[ok], levels = c("HC", "RA"), direction = "<", quiet = TRUE))
}

# STEP 4: evaluate the four models
hdr("STEP 4  PANEL vs COMPOSITION")
rows <- list(); lrt_rows <- list(); rocs <- list()

for (variant in names(PANELS)) {
for (sx in c("F", "M")) {
  sexlab <- if (sx == "F") "Female" else "Male"
  genes  <- PANELS[[variant]][[sx]]
  cat(sprintf("\n---------------- %s panel [%s] ----------------\n", sexlab, variant))

  for (dn in names(datasets)) {
    d <- datasets[[dn]]
    sel <- which(d$meta$sex == sx)
    if (length(sel) < 6) { say("  [%s] skipped (n = %d)", d$label, length(sel)); next }
    y <- factor(d$meta$group[sel], levels = c("HC", "RA"))
    if (length(unique(y)) < 2) { say("  [%s] skipped (one class)", d$label); next }

    present <- genes[genes %in% rownames(d$expr)]
    if (length(present) < 2) { say("  [%s] skipped (<2 panel genes present)", d$label); next }
    G <- t(d$expr[present, d$meta$sample[sel], drop = FALSE])   # samples x genes
    C <- d$comp[sel, , drop = FALSE]                            # samples x cPCs
    k <- min(CFG$kfold, min(table(y)))

    models <- list(
      `A panel only`         = glm_oof(G,            y, k),
      `B composition only`   = glm_oof(C,            y, k),
      `C panel + composition`= glm_oof(cbind(G, C),  y, k),
      `D panel residualised` = glm_oof_resid(G, C,   y, k))

    for (mn in names(models)) {
      r <- safe_roc(y, models[[mn]])
      if (is.null(r)) next
      a <- auc_ci(r)
      rows[[length(rows) + 1]] <- data.table(
        sex = sexlab, panel = variant, dataset = d$label, model = mn,
        n = length(y), n_RA = sum(y == "RA"), n_HC = sum(y == "HC"),
        n_panel_genes = length(present),
        AUC = round(a["auc"], 3), AUC_lo = round(a["lo"], 3), AUC_hi = round(a["hi"], 3),
        CI_method = ifelse(a["ci_method"] == 1, "bootstrap (n<20)", "DeLong"),
        evidence_tier = ifelse(sx == "M" || length(y) < CFG$small_n,
                               "EXPLORATORY (underpowered)", "primary"))
      rocs[[paste(variant, sexlab, d$label, mn)]] <- r
    }

    # DeLong comparison: panel vs the composition-only benchmark
    rA <- safe_roc(y, models[["A panel only"]]); rB <- safe_roc(y, models[["B composition only"]])
    rC <- safe_roc(y, models[["C panel + composition"]])
    dl <- function(r1, r2) if (is.null(r1) || is.null(r2)) NA_real_ else
      tryCatch(suppressWarnings(roc.test(r1, r2, method = "delong",
                                         paired = TRUE)$p.value), error = function(e) NA_real_)

    # likelihood-ratio test: does the panel add signal beyond composition? (fitted in-sample)
    Zg <- scale(G); Zg[!is.finite(Zg)] <- 0
    Zc <- scale(C); Zc[!is.finite(Zc)] <- 0
    dfB <- data.frame(y = y, Zc, check.names = FALSE)
    dfC <- data.frame(y = y, Zg, Zc, check.names = FALSE)
    mB <- suppressWarnings(stats::glm(y ~ ., data = dfB, family = binomial))
    mC <- suppressWarnings(stats::glm(y ~ ., data = dfC, family = binomial))
    an <- suppressWarnings(stats::anova(mB, mC, test = "LRT"))
    lrt_p  <- an$`Pr(>Chi)`[2]; lrt_df <- an$Df[2]; lrt_dev <- an$Deviance[2]
    sep <- any(abs(stats::coef(mC)[-1]) > 10, na.rm = TRUE)   # quasi-separation guard

    lrt_rows[[length(lrt_rows) + 1]] <- data.table(
      sex = sexlab, panel = variant, dataset = d$label, n = length(y),
      AUC_panel = round(as.numeric(auc(rA)), 3),
      AUC_composition = round(as.numeric(auc(rB)), 3),
      AUC_both = if (is.null(rC)) NA_real_ else round(as.numeric(auc(rC)), 3),
      delta_panel_minus_composition = round(as.numeric(auc(rA)) - as.numeric(auc(rB)), 3),
      delong_p_panel_vs_composition = signif(dl(rA, rB), 3),
      LRT_df = lrt_df, LRT_deviance = round(lrt_dev, 2),
      LRT_p_panel_beyond_composition = signif(lrt_p, 3),
      separation_warning = sep,
      evidence_tier = ifelse(sx == "M" || length(y) < CFG$small_n,
                             "EXPLORATORY (underpowered)", "primary"))

    say("  [%-14s] n=%3d  A panel %.3f | B composition %.3f | C both %.3f | D residualised %.3f",
        d$label, length(y),
        as.numeric(auc(rA)), as.numeric(auc(rB)),
        if (is.null(rC)) NA_real_ else as.numeric(auc(rC)),
        {rD <- safe_roc(y, models[["D panel residualised"]])
         if (is.null(rD)) NA_real_ else as.numeric(auc(rD))})
    say("                    DeLong p (A vs B) = %s | LRT p (panel beyond composition) = %s%s",
        format.pval(dl(rA, rB), digits = 3), format.pval(lrt_p, digits = 3),
        ifelse(sep, "  [SEPARATION WARNING]", ""))
  }
}
}

auc_tab <- rbindlist(rows)
lrt_tab <- rbindlist(lrt_rows)
fwrite(auc_tab, file.path(tab, "PANEL_auc_celladjusted.csv"))
fwrite(lrt_tab, file.path(tab, "PANEL_incremental_value_LRT.csv"))

# STEP 5: verdict (a significant LRT and a materially better AUC are different claims; MIN_DELTA separates them)
hdr("STEP 5  VERDICT")
MIN_DELTA <- 0.02

verdict <- lrt_tab[, .(
  sex, panel, dataset, n, AUC_panel, AUC_composition, delta = delta_panel_minus_composition,
  LRT_p = LRT_p_panel_beyond_composition, evidence_tier,
  perfect_separation = AUC_panel >= 0.999,
  verdict = fifelse(evidence_tier != "primary",
                    "EXPLORATORY - not interpretable at this n",
             fifelse(is.na(LRT_p_panel_beyond_composition),
                    "INDETERMINATE",
             fifelse(LRT_p_panel_beyond_composition < 0.05 &
                     delta_panel_minus_composition >= MIN_DELTA,
                    "PANEL ADDS SIGNAL beyond composition",
             fifelse(LRT_p_panel_beyond_composition < 0.05 &
                     delta_panel_minus_composition > 0,
                    "MARGINAL - statistically additive but AUC gain < 0.02",
             fifelse(LRT_p_panel_beyond_composition < 0.05,
                    "panel adds signal but does not out-discriminate composition alone",
                    "NOT DISTINGUISHABLE from a differential white-cell count"))))))]
verdict[perfect_separation == TRUE,
        verdict := paste0(verdict, " [AUC = 1.000: perfect separation, treat as unstable]")]
fwrite(verdict, file.path(tab, "PANEL_auc_celladjusted_summary.csv"))
print(verdict)

cat("\n---- full AUC table ----\n")
print(auc_tab[, .(sex, panel, dataset, model, n, AUC, AUC_lo, AUC_hi, evidence_tier)])

saveRDS(list(auc = auc_tab, lrt = lrt_tab, verdict = verdict, rocs = rocs, config = CFG),
        file.path(procN, "panel_celladjusted_objects.rds"))

say("")
say("READING RULE FOR THE THESIS")
say("  Model A (panel only) may no longer be reported on its own. Every panel AUC")
say("  must be accompanied by model B, the composition-only benchmark, and by the")
say("  LRT of C against B. A panel that does not beat B is a white-cell count.")
say("  Every Male row is EXPLORATORY: n = 38/13/9 cannot support a diagnostic claim.")
say("")
say("Wrote PANEL_auc_celladjusted.csv, PANEL_incremental_value_LRT.csv,")
say("PANEL_auc_celladjusted_summary.csv and panel_celladjusted_objects.rds")
cat("\nDONE\n")
