#!/usr/bin/env Rscript
# Testing phase: apply the five locked models from 37_ to an independent external blood cohort (GSE15573), per sex, and merge into ml_algo_roc.rds.
suppressMessages({ library(caret); library(data.table); library(pROC); library(Biobase) })
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"

fit <- readRDS(file.path(procN, "ml_algo_models.rds"))

# external blood test set: GSE15573; fall back to sibling project's raw copy if data/raw symlink doesn't resolve
raw_candidates <- c(
  "data/raw/GSE15573_raw.rds",
  "/Users/swetarai/Library/CloudStorage/Dropbox/THESIS_SWETA_28_MAY/Thesis_chapters/Research_Q2_TRANSCRIPTOMICS_sexstratified/data/raw/GSE15573_raw.rds")
raw_path <- raw_candidates[file.exists(raw_candidates)][1]
if (is.na(raw_path)) stop("GSE15573_raw.rds not found in any known location.")

e <- readRDS(raw_path); if (is.list(e) && !is(e, "ExpressionSet")) e <- e[[1]]
x <- exprs(e); if (max(x, na.rm = TRUE) > 50) x <- log2(x + 1)
sym <- fData(e)[["Gene symbol"]]
keep <- !is.na(sym) & sym != ""; x <- x[keep, ]; sym <- sym[keep]
rmean <- rowMeans(x)
best <- tapply(seq_along(sym), sym, function(ix) ix[which.max(rmean[ix])])
xg <- x[unlist(best), ]; rownames(xg) <- names(best)
p <- pData(e)
grp <- ifelse(grepl("Rheumatoid|RA", p[["status:ch1"]], ignore.case = TRUE), "RA", "HC")
sex <- ifelse(grepl("Female", p[["gender:ch1"]], ignore.case = TRUE), "F", "M")
blood_ext <- list(expr = xg, group = factor(grp, levels = c("HC", "RA")), sex = sex,
                   label = "External blood test (GSE15573)")

cat(sprintf("External blood test: n=%d (RA=%d/HC=%d) | F=%d M=%d | source=%s\n",
            length(blood_ext$group), sum(blood_ext$group == "RA"), sum(blood_ext$group == "HC"),
            sum(blood_ext$sex == "F"), sum(blood_ext$sex == "M"), raw_path))

# helpers (identical to 38_)
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
score_model <- function(m, genes, colnms, X_genesXsamples, y) {
  present <- genes[genes %in% rownames(X_genesXsamples)]
  Z <- zrows(X_genesXsamples[present, , drop = FALSE])
  Zdf <- as.data.frame(t(Z), check.names = FALSE); colnames(Zdf) <- make.names(present)
  for (g in setdiff(colnms, colnames(Zdf))) Zdf[[g]] <- 0
  Zdf <- Zdf[, colnms, drop = FALSE]
  p <- predict(m, newdata = Zdf, type = "prob")[, "RA"]
  r <- roc(y, p, direction = "<", levels = c("HC", "RA"), quiet = TRUE)
  list(roc = r, ci = auc_ci(r), n = length(y), n_RA = sum(y == "RA"), n_ctrl = sum(y == "HC"),
       missing = setdiff(genes, present))
}
rc <- function(sexlab, algo, dataset, r, tier) {
  data.table(sex = sexlab, algorithm = algo, dataset = dataset,
             sens = r$sensitivities, spec = r$specificities,
             auc = as.numeric(auc(r)), evidence_tier = tier)
}
ALGO_LABELS <- list(logreg = "Logistic regression", svm = "SVM (RBF)", knn = "k-NN",
                     rf = "Random forest", ann = "ANN (nnet)")

perf_rows <- list(); roc_rows <- list()
for (sx in c("F", "M")) {
  sexlab <- if (sx == "F") "Female" else "Male"
  res <- if (sx == "F") fit$female else fit$male
  n_train <- length(res$y)
  tier <- if (n_train < 50) "EXPLORATORY (underpowered)" else "primary"
  ti <- which(blood_ext$sex == sx)
  if (length(ti) < 3) { cat(sprintf("  [%s] skipped: too few test samples\n", sexlab)); next }
  y_te <- factor(blood_ext$group[ti], levels = c("HC", "RA"))
  Xte <- blood_ext$expr[, ti, drop = FALSE]

  for (a in names(res$models)) {
    mm <- res$models[[a]]
    spec_label <- ALGO_LABELS[[a]]
    perf_rows[[length(perf_rows) + 1]] <- data.table(
      sex = sexlab, algorithm = spec_label, dataset = "Train (resampled CV)",
      setting = "External blood test (GSE15573)", evidence_tier = tier,
      n = n_train, n_RA = sum(res$y == "RA"), n_ctrl = sum(res$y == "HC"),
      AUC_CI = fmt(auc_ci(mm$oof), n_train), missing_genes = "")
    roc_rows[[length(roc_rows) + 1]] <- rc(sexlab, spec_label, "Train (resampled CV)", mm$oof, tier)

    te <- score_model(mm$model, mm$genes, mm$colnames, Xte, y_te)
    perf_rows[[length(perf_rows) + 1]] <- data.table(
      sex = sexlab, algorithm = spec_label, dataset = "External blood test (GSE15573)",
      setting = "External blood test (GSE15573)", evidence_tier = tier,
      n = te$n, n_RA = te$n_RA, n_ctrl = te$n_ctrl,
      AUC_CI = fmt(te$ci, te$n), missing_genes = paste(te$missing, collapse = ";"))
    roc_rows[[length(roc_rows) + 1]] <- rc(sexlab, spec_label, "External blood test (GSE15573)", te$roc, tier)
  }
  cat(sprintf("  [%s] done (train n=%d, external blood test n=%d, tier=%s)\n",
              sexlab, n_train, length(ti), tier))
}

externalblood <- list(perf = rbindlist(perf_rows), roc = rbindlist(roc_rows))
fwrite(externalblood$perf, file.path(tab, "ML_performance_externalblood.csv"))

# merge into the existing ml_algo_roc.rds (preserve 38_'s two settings)
prior <- readRDS(file.path(procN, "ml_algo_roc.rds"))
saveRDS(list(sametissue = prior$sametissue, crosstissue = prior$crosstissue,
             externalblood = externalblood, seed = GLOBAL_SEED,
             built = "08_crosstissue/38b_testing_ml_algorithms_externalblood.R"),
        file.path(procN, "ml_algo_roc.rds"))

cat("\n====================  EXTERNAL BLOOD PERFORMANCE  ====================\n")
print(externalblood$perf, width = 200)
cat("\nSaved ML_performance_externalblood.csv; updated ml_algo_roc.rds with $externalblood\n")
