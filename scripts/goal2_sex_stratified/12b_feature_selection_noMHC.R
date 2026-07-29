#!/usr/bin/env Rscript
# =============================================================================
# 12b_feature_selection_noMHC.R
# -----------------------------------------------------------------------------
# REBUILD of the sex-stratified diagnostic panels from the MHC-FREE candidate
# set, using the identical three-selector consensus design as 12_feature_selection.R.
#
# WHY THIS EXISTS
#   10c_MR_mhc_sensitivity.R and 10d_coloc_panel_genes.R established two things
#   about the panels produced by 12_feature_selection.R:
#     * three of the six genes in each panel are MHC-affected (female GNL1 and
#       C6orf136, male VPS52 and HLA-DMA are untestable once MHC instruments are
#       removed; ESYT1 loses FDR in both sexes to BH re-ranking);
#     * six of the nine unique panel genes show PP.H3 >= 0.8 in colocalisation,
#       i.e. positive evidence that the eQTL and the RA association are driven by
#       DIFFERENT causal variants.
#
#   Diagnosing that and leaving the panel unchanged would be indefensible: it
#   would mean the sensitivity analysis was run but not acted on. This script
#   therefore re-derives the panels from the candidate set that survives MHC
#   exclusion, so the thesis can report a primary panel and an MHC-free panel
#   side by side and let the reader see exactly what the MHC was contributing.
#
#   NOTE ON WHAT THIS DOES AND DOES NOT REPAIR. Removing the MHC removes the
#   LD-confounding that is attributable to the HLA region. It does NOT make the
#   surviving genes colocalised - INPP5B (PP.H3 = 0.929) and ESYT1 (0.912) are
#   non-MHC and still fail colocalisation. The MHC-free panel is therefore a
#   CLEANER PRIORITISATION, not a causally validated one. No gene in either
#   panel may carry a causal claim. See 10d and results/RESULTS_ROBUSTNESS.md.
#
# DESIGN - IDENTICAL TO 12_feature_selection.R, DELIBERATELY
#   Same three selectors (LASSO lambda.min, Random Forest above-mean Gini,
#   SVM-RFE with CV-chosen size), same tuning grids, same 3-method intersection,
#   same seed (1234), same training matrix. The ONLY change is the input gene
#   list: FS_input_{sex}_noMHC.csv instead of FS_input_{sex}.csv. Holding the
#   procedure fixed is what makes the two panels comparable; if the design also
#   changed, any difference between the panels would be uninterpretable.
#
#   Candidate counts: 14 female / 14 male (from 32 / 25 with the MHC retained).
#   The candidate sets now overlap almost completely (13 of 14 genes shared),
#   which is expected: the MHC was supplying most of what distinguished them, and
#   it was supplying it through LD with a single dominant locus rather than
#   through sex-differential biology.
#
#   in : data/processed/combined_train.rds
#        results/tables/FS_input_{female,male}_noMHC.csv
#   out: data/processed/new/ml_features_noMHC.rds
#        results/tables/mr_fs_selected_bymethod_{female,male}_noMHC.csv
#        results/tables/mr_fs_consensus_{female,male}_noMHC.csv
#        results/tables/mr_fs_summary_noMHC.csv
#        results/tables/PANEL_primary_vs_noMHC_membership.csv
#
#   Friedman J, et al. J Stat Softw 2010;33:1-22.        (glmnet / LASSO)
#   Breiman L. Mach Learn 2001;45:5-32.                  (random forest)
#   Guyon I, et al. Mach Learn 2002;46:389-422.          (SVM-RFE)
#   Ambroise C, McLachlan GJ. PNAS 2002;99:6562-6566.    (selection bias)
# =============================================================================
suppressMessages({
  library(glmnet); library(randomForest); library(e1071); library(caret)
  library(data.table)
})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234
set.seed(GLOBAL_SEED)

proc  <- "data/processed"
procN <- "data/processed/new"
tab   <- "results/tables"
dir.create(procN, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# =============================================================================
# STEP 1 — INPUTS
# =============================================================================
hdr("STEP 1  MHC-FREE CANDIDATE SETS")
o    <- readRDS(file.path(proc, "combined_train.rds"))
expr <- o$expr
meta <- o$meta

mr_female <- fread(file.path(tab, "FS_input_female_noMHC.csv"))
mr_male   <- fread(file.path(tab, "FS_input_male_noMHC.csv"))
mr_female <- mr_female[gene %in% rownames(expr)]
mr_male   <- mr_male[gene   %in% rownames(expr)]
stopifnot(nrow(mr_female) >= 2L, nrow(mr_male) >= 2L)

say("MHC-free candidates present in the expression matrix: female %d, male %d",
    nrow(mr_female), nrow(mr_male))
say("  female: %s", paste(mr_female$gene, collapse = ", "))
say("  male  : %s", paste(mr_male$gene,   collapse = ", "))
say("shared between the two candidate sets: %d",
    length(intersect(mr_female$gene, mr_male$gene)))

# =============================================================================
# STEP 2 — THE THREE SELECTORS (verbatim from 12_feature_selection.R)
# =============================================================================
svm_rfe_rank <- function(X, y, cost = 1) {
  feats <- colnames(X); ranking <- character(0)
  while (length(feats) > 1) {
    m  <- svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = TRUE, cost = cost)
    w2 <- ((t(m$coefs) %*% m$SV)[1, ])^2
    drop <- names(sort(w2))[1]
    ranking <- c(drop, ranking); feats <- setdiff(feats, drop)
  }
  c(feats, ranking)
}
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
SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)
tune_svm_cost <- function(X, y, grid = SVM_COST_GRID) {
  set.seed(GLOBAL_SEED)
  tc <- tune(svm, train.x = X, train.y = y, kernel = "linear", scale = TRUE,
             ranges = list(cost = grid),
             tunecontrol = tune.control(sampling = "cross", cross = 10))
  tc$best.parameters$cost
}
rf_mtry_grid <- function(p) sort(unique(pmin(p, c(1, 2, floor(sqrt(p)),
                                                  floor(p / 3), floor(p / 2), p))))
tune_rf_mtry <- function(X, y, ntree, grid) {
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE,
                       summaryFunction = twoClassSummary)
  set.seed(GLOBAL_SEED)
  fit <- train(x = X, y = y, method = "rf", metric = "ROC", trControl = ctrl,
               tuneGrid = expand.grid(mtry = grid), ntree = ntree, importance = TRUE)
  list(mtry = fit$bestTune$mtry, fit = fit)
}

# =============================================================================
# STEP 3 — PER-SEX FEATURE SELECTION
# =============================================================================
hdr("STEP 3  THREE-METHOD CONSENSUS ON THE MHC-FREE SET")

run_fs <- function(sex_code) {
  sexlab <- if (sex_code == "F") "Female" else "Male"
  mr     <- if (sex_code == "F") mr_female else mr_male
  genes  <- mr$gene
  cols   <- meta$sample[meta$sex == sex_code]

  X <- t(expr[genes, cols, drop = FALSE])
  colnames(X) <- make.names(genes)
  lk   <- setNames(genes, make.names(genes))
  back <- function(v) unname(lk[v])
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))

  cat(sprintf("\n================  %s  ================\n", sexlab))
  say("candidates: %d | samples: %d (RA=%d, HC=%d)",
      ncol(X), nrow(X), sum(y == "RA"), sum(y == "HC"))

  set.seed(GLOBAL_SEED)
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10,
                  type.measure = "deviance")
  co_min <- coef(cv, s = "lambda.min")[-1, 1]
  co_1se <- coef(cv, s = "lambda.1se")[-1, 1]
  lasso_min <- back(names(co_min)[co_min != 0])
  lasso_1se <- back(names(co_1se)[co_1se != 0])
  say("  LASSO : lambda.min=%.4f -> %d genes | lambda.1se=%.4f -> %d genes",
      cv$lambda.min, length(lasso_min), cv$lambda.1se, length(lasso_1se))

  mtry_grid <- rf_mtry_grid(ncol(X))
  rf_tune   <- tune_rf_mtry(X, y, ntree = 1000, grid = mtry_grid)
  set.seed(GLOBAL_SEED)
  rf <- randomForest(X, y, importance = TRUE, ntree = 1000, mtry = rf_tune$mtry)
  gini     <- sort(rf$importance[, "MeanDecreaseGini"], decreasing = TRUE)
  gini_thr <- mean(gini)
  rf_sel   <- back(names(gini)[gini > gini_thr])
  say("  RF    : ntree=1000, tuned mtry=%d, Gini cutoff=%.3f -> %d genes",
      rf_tune$mtry, gini_thr, length(rf_sel))

  best_cost <- tune_svm_cost(X, y)
  set.seed(GLOBAL_SEED)
  rank  <- svm_rfe_rank(X, y, cost = best_cost)
  curve <- svm_rfe_curve(X, y, rank, cost = best_cost)
  svm_sel <- back(rank[1:curve$best])
  say("  SVM-RFE: tuned cost=%s, CV-optimal size=%d (err=%.4f) -> %d genes",
      best_cost, curve$best, curve$besterr, length(svm_sel))

  sets      <- list(LASSO = lasso_min, RandomForest = rf_sel, SVM_RFE = svm_sel)
  consensus <- sort(Reduce(intersect, sets))
  say("  CONSENSUS (LASSO n RF n SVM-RFE): %d genes -> %s",
      length(consensus),
      if (length(consensus)) paste(consensus, collapse = ", ") else "NONE")
  if (!length(consensus)) {
    say("  WARNING: empty consensus. With %d candidates the three selectors did", ncol(X))
    say("  not agree on any gene. Report as a NEGATIVE result, not as a panel.")
  }

  list(sexlab = sexlab, mr = mr, mr_genes = genes, sets = sets,
       consensus = consensus, union = sort(Reduce(union, sets)),
       lasso_min = lasso_min, lasso_1se = lasso_1se,
       cv = cv, rf = rf, gini = gini, gini_thr = gini_thr,
       svm_rank = rank, svm_curve = curve, samples = cols,
       tuned = list(lasso_lambda_min = cv$lambda.min,
                    lasso_lambda_1se = cv$lambda.1se,
                    rf_mtry = rf_tune$mtry, svm_cost = best_cost))
}

F <- run_fs("F")
M <- run_fs("M")

# =============================================================================
# STEP 4 — OUTPUT AND HEAD-TO-HEAD WITH THE PRIMARY PANELS
# =============================================================================
hdr("STEP 4  PRIMARY vs MHC-FREE PANEL MEMBERSHIP")

for (r in list(F, M)) {
  s <- tolower(r$sexlab)
  bym <- rbindlist(lapply(names(r$sets), function(n)
    data.table(method = n, gene = r$sets[[n]])))
  fwrite(bym, file.path(tab, sprintf("mr_fs_selected_bymethod_%s_noMHC.csv", s)))
  fwrite(data.table(gene = r$consensus),
         file.path(tab, sprintf("mr_fs_consensus_%s_noMHC.csv", s)))
}

summary_tab <- data.table(
  sex             = c("Female", "Male"),
  candidate_set   = "MHC-free (FS_input_*_noMHC)",
  n_samples       = c(length(F$samples), length(M$samples)),
  n_candidates    = c(length(F$mr_genes), length(M$mr_genes)),
  n_lasso         = c(length(F$sets$LASSO), length(M$sets$LASSO)),
  n_rf            = c(length(F$sets$RandomForest), length(M$sets$RandomForest)),
  n_svmrfe        = c(length(F$sets$SVM_RFE), length(M$sets$SVM_RFE)),
  n_consensus     = c(length(F$consensus), length(M$consensus)),
  consensus_genes = c(paste(F$consensus, collapse = "; "),
                      paste(M$consensus, collapse = "; ")))
fwrite(summary_tab, file.path(tab, "mr_fs_summary_noMHC.csv"))
print(summary_tab)

# membership comparison against the primary panels
prim <- list(Female = character(0), Male = character(0))
ml_path <- file.path(procN, "ml_features.rds")
if (file.exists(ml_path)) {
  ml <- readRDS(ml_path)
  prim$Female <- ml$female$consensus
  prim$Male   <- ml$male$consensus
}
newp <- list(Female = F$consensus, Male = M$consensus)

# colocalisation / MHC status per gene, so the comparison carries its provenance
cl <- if (file.exists(file.path(tab, "COLOC_results.csv")))
        fread(file.path(tab, "COLOC_results.csv")) else NULL
fa <- if (file.exists(file.path(tab, "MR_MHC_sensitivity_panel_fate.csv")))
        fread(file.path(tab, "MR_MHC_sensitivity_panel_fate.csv")) else NULL

memb <- rbindlist(lapply(c("Female", "Male"), function(sx) {
  gs <- union(prim[[sx]], newp[[sx]])
  if (!length(gs)) return(NULL)
  d <- data.table(sex = sx, gene = gs,
                  in_primary_panel = gs %in% prim[[sx]],
                  in_noMHC_panel   = gs %in% newp[[sx]])
  d[, status := fifelse(in_primary_panel & in_noMHC_panel, "retained in both",
                 fifelse(in_primary_panel, "DROPPED (was MHC-dependent)",
                                           "NEW (surfaced once MHC removed)"))]
  if (!is.null(cl)) d <- merge(d, unique(cl[, .(gene, MHC_gene, PP3 = round(PP3, 3),
                                                PP4 = round(PP4, 3), coloc_verdict)],
                                         by = "gene"), by = "gene", all.x = TRUE)
  if (!is.null(fa)) d <- merge(d, unique(fa[, .(gene, mhc_verdict = verdict)], by = "gene"),
                               by = "gene", all.x = TRUE)
  d
}), fill = TRUE)
setorder(memb, sex, -in_primary_panel, gene)
fwrite(memb, file.path(tab, "PANEL_primary_vs_noMHC_membership.csv"))
print(memb[, .(sex, gene, in_primary_panel, in_noMHC_panel, MHC_gene, PP3, status)])

saveRDS(list(female = F, male = M, expr = expr, meta = meta, seed = GLOBAL_SEED,
             candidate_set = "MHC-free",
             built = "goal2_sex_stratified/12b_feature_selection_noMHC.R"),
        file.path(procN, "ml_features_noMHC.rds"))

hdr("READING RULE")
say("  The MHC-free panel is a CLEANER PRIORITISATION, not a causally validated")
say("  set. Colocalisation still fails for its members (10d), so no gene here")
say("  carries a causal claim either. What the MHC-free panel does provide is a")
say("  diagnostic model that cannot be attributed to LD with HLA-DRB1.")
say("  Report both panels. Do not quietly replace one with the other.")
say("")
say("Wrote ml_features_noMHC.rds, mr_fs_*_noMHC.csv and")
say("PANEL_primary_vs_noMHC_membership.csv")
cat("\nDONE\n")
