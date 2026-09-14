#!/usr/bin/env Rscript
# Single authoritative nested-CV table: recomputes every candidate-set x selector x sex combination in one process under one seed policy, reconciled against the legacy per-script numbers.
suppressMessages({
  library(glmnet); library(randomForest); library(e1071)
  library(pROC); library(caret); library(data.table)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
set.seed(GLOBAL_SEED)

proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

o <- readRDS(file.path(proc, "combined_train.rds"))
expr <- o$expr; meta <- as.data.table(o$meta)

CAND <- list(
  primary = list(F = fread(file.path(tab, "FS_input_female.csv"))$gene,
                 M = fread(file.path(tab, "FS_input_male.csv"))$gene),
  noMHC   = list(F = fread(file.path(tab, "FS_input_female_noMHC.csv"))$gene,
                 M = fread(file.path(tab, "FS_input_male_noMHC.csv"))$gene))

RESAMPLE <- list(F = list(k = 10, reps = 5), M = list(k = 5, reps = 10))
ALPHAS <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

# Selectors: single definition, used by every variant
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
select_consensus <- function(X, y) {
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
fit_enet <- function(X, y) {
  best <- NULL; bcv <- Inf
  for (a in ALPHAS) {
    cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = a, nfolds = 5,
                             standardize = TRUE), error = function(e) NULL)
    if (!is.null(cv) && min(cv$cvm) < bcv) { bcv <- min(cv$cvm); best <- cv }
  }
  best
}

# One nested-CV engine; every variant goes through it, so reported differences trace to arguments alone
run_nested <- function(sex_code, genes, selector) {
  rs <- RESAMPLE[[sex_code]]
  genes <- unique(genes[genes %in% rownames(expr)])
  cols <- meta$sample[meta$sex == sex_code]
  X0 <- t(expr[genes, cols, drop = FALSE]); colnames(X0) <- make.names(genes)
  y  <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))

  preds <- list(); rep_auc <- numeric(0); nsel <- integer(0)
  for (rp in seq_len(rs$reps)) {
    set.seed(1000 + rp)                       # SINGLE seed policy for all variants
    folds <- createFolds(y, k = rs$k, returnTrain = FALSE)
    rp_pred <- data.table()
    for (fi in seq_along(folds)) {
      te <- folds[[fi]]; tr <- setdiff(seq_along(y), te)
      if (length(unique(y[tr])) < 2) next
      Xtr <- X0[tr, , drop = FALSE]; Xte <- X0[te, , drop = FALSE]; ytr <- y[tr]

      if (selector == "consensus") {
        P <- select_consensus(Xtr, ytr); nsel <- c(nsel, length(P))
        mu <- colMeans(Xtr[, P, drop = FALSE])
        sg <- apply(Xtr[, P, drop = FALSE], 2, sd); sg[sg == 0 | is.na(sg)] <- 1
        Ztr <- scale(Xtr[, P, drop = FALSE], center = mu, scale = sg)
        Zte <- scale(Xte[, P, drop = FALSE], center = mu, scale = sg)
        fit <- suppressWarnings(glm(ytr ~ ., data = data.frame(ytr, Ztr, check.names = FALSE),
                                    family = binomial))
        p <- as.numeric(predict(fit, newdata = data.frame(Zte, check.names = FALSE),
                                type = "response"))
      } else {
        b <- fit_enet(Xtr, ytr); if (is.null(b)) next
        co <- as.numeric(coef(b, s = "lambda.min"))[-1]; nsel <- c(nsel, sum(co != 0))
        p <- as.numeric(predict(b, newx = Xte, s = "lambda.min", type = "response"))
      }
      rp_pred <- rbind(rp_pred, data.table(sample = rownames(Xte), prob = p, obs = y[te]))
    }
    if (!nrow(rp_pred)) next
    rr <- roc(rp_pred$obs, rp_pred$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
    rep_auc <- c(rep_auc, as.numeric(auc(rr))); preds[[rp]] <- rp_pred
  }
  agg <- rbindlist(preds)[, .(prob = mean(prob), obs = obs[1]), by = sample]
  ro  <- roc(agg$obs, agg$prob, levels = c("HC", "RA"), direction = "<", quiet = TRUE)
  n <- length(y)
  ci <- if (n < 20) { set.seed(GLOBAL_SEED)
    suppressWarnings(as.numeric(ci.auc(ro, method = "bootstrap", boot.n = 2000)))
  } else suppressWarnings(as.numeric(ci.auc(ro)))
  list(auc = as.numeric(auc(ro)), lo = ci[1], hi = ci[3], n = n,
       rep_mean = mean(rep_auc), rep_sd = stats::sd(rep_auc),
       med_genes = as.integer(stats::median(nsel)), agg = agg, roc = ro)
}

# Run every candidate-set x selector x sex combination
hdr("NESTED CV - ALL VARIANTS, ONE PROCESS, ONE SEED POLICY")
grid <- CJ(candidate_set = c("primary", "noMHC"),
           selector = c("consensus", "elasticnet"),
           sex = c("F", "M"), sorted = FALSE)

res <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i]
  r <- run_nested(g$sex, CAND[[g$candidate_set]][[g$sex]], g$selector)
  key <- paste(g$candidate_set, g$selector, g$sex)
  res[[key]] <- r
  say("%-9s %-11s %s : AUC %.3f (%.3f-%.3f) | per-repeat %.3f +/- %.3f | median genes %d",
      g$candidate_set, g$selector, g$sex, r$auc, r$lo, r$hi, r$rep_mean, r$rep_sd, r$med_genes)
}

auth <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i]; r <- res[[paste(g$candidate_set, g$selector, g$sex)]]
  sexlab <- if (g$sex == "F") "Female" else "Male"
  data.table(
    sex = sexlab, candidate_set = g$candidate_set, selector = g$selector,
    n = r$n, resampling = sprintf("%d-fold x %d repeats",
                                  RESAMPLE[[g$sex]]$k, RESAMPLE[[g$sex]]$reps),
    nested_AUC = round(r$auc, 3),
    CI_lo = round(r$lo, 3), CI_hi = round(r$hi, 3),
    per_repeat_mean = round(r$rep_mean, 3), per_repeat_sd = round(r$rep_sd, 3),
    median_genes_used = r$med_genes,
    evidence_tier = if (g$sex == "M") "EXPLORATORY (underpowered)" else "primary",
    recommended = (g$candidate_set == "noMHC" & g$selector == "consensus"))
}))
setorder(auth, sex, candidate_set, selector)
fwrite(auth, file.path(tab, "NESTED_CV_AUTHORITATIVE.csv"))

hdr("AUTHORITATIVE TABLE")
print(auth[, .(sex, candidate_set, selector, n, nested_AUC, CI_lo, CI_hi,
               per_repeat_sd, recommended)])

# Reconciliation against the legacy per-script tables
hdr("RECONCILIATION WITH THE LEGACY NUMBERS")
legacy <- data.table(
  source = c("14_model_training_nested_cv.R", "14_model_training_nested_cv.R",
             "16_model_training_final_panel.R", "16_model_training_final_panel.R",
             "16b_model_training_final_panel_noMHC.R", "16b_model_training_final_panel_noMHC.R"),
  sex = c("Female", "Male", "Female", "Male", "Female", "Male"),
  legacy_AUC = c(0.816, 0.896, 0.809, 0.952, 0.805, 0.796),
  matches_variant = c("primary/consensus", "primary/consensus",
                      "primary/elasticnet", "primary/elasticnet",
                      "primary/consensus", "primary/consensus"))
legacy[, c("candidate_set", "selector") := tstrsplit(matches_variant, "/")]
rec <- merge(legacy, auth[, .(sex, candidate_set, selector, recomputed_AUC = nested_AUC)],
             by = c("sex", "candidate_set", "selector"), all.x = TRUE)
rec[, difference := round(recomputed_AUC - legacy_AUC, 3)]
rec[, note := fifelse(abs(difference) <= 0.02,
       "reproduces within 0.02 - Monte-Carlo variation in fold assignment",
       "DOES NOT reproduce - investigate before citing the legacy value")]
setcolorder(rec, c("source", "sex", "candidate_set", "selector",
                   "legacy_AUC", "recomputed_AUC", "difference", "note"))
setorder(rec, sex, source)
fwrite(rec, file.path(tab, "NESTED_CV_legacy_reconciliation.csv"))
print(rec[, .(source = substr(source, 1, 34), sex, legacy_AUC, recomputed_AUC, difference)])

saveRDS(list(authoritative = auth, reconciliation = rec, results = res),
        file.path(procN, "nested_cv_authoritative.rds"))

hdr("READING RULE")
say("  CITE NESTED_CV_AUTHORITATIVE.csv. Nothing else.")
say("  The recommended row is candidate_set = noMHC, selector = consensus:")
rec_row <- auth[recommended == TRUE & sex == "Female"]
if (nrow(rec_row))
  say("    Female %.3f (%.3f-%.3f), n = %d", rec_row$nested_AUC,
      rec_row$CI_lo, rec_row$CI_hi, rec_row$n)
say("  The other rows are NOT competing estimates of one quantity - they are")
say("  different procedures, and the columns say which. Report the recommended")
say("  row as the headline and the rest as a sensitivity grid.")
say("  Male rows remain EXPLORATORY at n = 38; their spread across procedures is")
say("  a statement about sample size, not about which procedure is better.")
say("")
say("Wrote NESTED_CV_AUTHORITATIVE.csv and NESTED_CV_legacy_reconciliation.csv")
cat("\nDONE\n")
