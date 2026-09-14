#!/usr/bin/env Rscript
# =============================================================================
# 12c_feature_selection_diagnostics.R
# -----------------------------------------------------------------------------
# Algorithm-level diagnostic figures for the three feature-selection methods
# run in 12_feature_selection.R (LASSO, Random Forest, SVM-RFE), per sex.
#
# WHY THIS EXISTS
#   12_feature_selection.R fits and saves the three selector objects ($cv, $rf,
#   $svm_rank, $svm_curve) but never plots them; 13_feature_selection_venn.R
#   only draws the CONSENSUS overlap. Nothing in the project showed what each
#   individual algorithm's selection process looked like (the LASSO
#   regularisation path, the RF importance ranking, the SVM-RFE elimination
#   curve). This script draws those, straight from the objects 12 already
#   fitted and saved — no models are re-run, so figures are reproducible byte-
#   for-byte from ml_features.rds without re-fitting anything stochastic.
#
#   in : data/processed/new/ml_features.rds
#   out: results/figures/FIG_G2_05_lasso_cv_{female,male}.png/.pdf
#        results/figures/FIG_G2_06_rf_importance_{female,male}.png/.pdf
#        results/figures/FIG_G2_07_svmrfe_curve_{female,male}.png/.pdf
#        results/figures/FIG_G2_08_feature_selection_diagnostics_composite.png/.pdf
# =============================================================================
suppressMessages({
  library(glmnet); library(randomForest); library(ggplot2); library(patchwork)
})
options(stringsAsFactors = FALSE)

procN <- "data/processed/new"; fig <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

hdr("STEP 1  LOAD FITTED SELECTOR OBJECTS")
ml <- readRDS(file.path(procN, "ml_features.rds"))

pal <- list(female = c(light = "#F7E6EF", dark = "#C2185B"),
            male   = c(light = "#E3EEF7", dark = "#1565C0"))

# safe-name -> original gene symbol, same order used to build X in script 12
back_map <- function(r) setNames(r$mr_genes, make.names(r$mr_genes))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  say("  wrote %s.{png,pdf}", name)
}

# =============================================================================
# STEP 2  LASSO cross-validation (regularisation) curve
# =============================================================================
lasso_cv_plot <- function(r, sexlab, dark) {
  cv <- r$cv
  df <- data.frame(loglambda = log(cv$lambda), cvm = cv$cvm, cvlo = cv$cvlo, cvup = cv$cvup)
  ytop <- max(df$cvup)
  ggplot(df, aes(loglambda, cvm)) +
    geom_ribbon(aes(ymin = cvlo, ymax = cvup), fill = dark, alpha = 0.15) +
    geom_line(colour = dark, linewidth = 0.5) +
    geom_point(colour = dark, size = 1.1) +
    geom_vline(xintercept = log(cv$lambda.min), linetype = "dashed", colour = "grey30") +
    geom_vline(xintercept = log(cv$lambda.1se), linetype = "dotted", colour = "grey30") +
    annotate("text", x = log(cv$lambda.min), y = ytop, label = "min", vjust = -0.4,
             size = 3, colour = "grey20") +
    annotate("text", x = log(cv$lambda.1se), y = ytop, label = "1se", vjust = -0.4,
             size = 3, colour = "grey20") +
    coord_cartesian(clip = "off") +
    labs(subtitle = sprintf("lambda.min = %.4f -> %d genes | lambda.1se = %.4f -> %d genes",
                             cv$lambda.min, length(r$lasso_min), cv$lambda.1se, length(r$lasso_1se)),
         x = "log(lambda)", y = "Binomial deviance (10-fold CV)") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(t = 2, b = 10)),
          plot.margin = margin(16, 10, 6, 6))
}

# =============================================================================
# STEP 3  Random Forest importance ranking
# =============================================================================
rf_importance_plot <- function(r, sexlab, light, dark) {
  light <- unname(light); dark <- unname(dark)
  bm <- back_map(r)
  df <- data.frame(gene = unname(bm[names(r$gini)]), gini = as.numeric(r$gini))
  df <- df[order(-df$gini), ]
  df$gene <- factor(df$gene, levels = rev(df$gene))
  df$selected <- df$gini > r$gini_thr
  ggplot(df, aes(gene, gini, fill = selected)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = r$gini_thr, linetype = "dashed", colour = "grey30") +
    coord_flip() +
    scale_fill_manual(values = setNames(c(dark, light), c("TRUE", "FALSE")), guide = "none") +
    labs(subtitle = sprintf("ntree = 1000, tuned mtry = %d | selected: Gini > mean (%.3f) -> %d genes",
                             r$tuned$rf_mtry, r$gini_thr, sum(df$selected)),
         x = NULL, y = "Mean decrease in Gini") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(b = 10)),
          axis.text.y = element_text(size = 7.6))
}

# =============================================================================
# STEP 4  SVM-RFE elimination curve
# =============================================================================
svmrfe_curve_plot <- function(r, sexlab, dark) {
  cu <- r$svm_curve
  df <- data.frame(k = cu$k, err = cu$err)
  best_pt <- df[df$k == cu$best, ]
  ggplot(df, aes(k, err)) +
    geom_line(colour = dark, linewidth = 0.5) +
    geom_point(colour = dark, size = 1.2) +
    geom_vline(xintercept = cu$best, linetype = "dashed", colour = "grey30") +
    geom_point(data = best_pt, aes(k, err), shape = 21, colour = "black",
               fill = dark, size = 2.6) +
    labs(subtitle = sprintf("tuned cost = %s | CV-optimal size = %d genes (error = %.4f)",
                             format(r$tuned$svm_cost, trim = TRUE), cu$best, cu$besterr),
         x = "Number of top-ranked features (k)", y = "10-fold CV error rate") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(b = 10)))
}

# =============================================================================
# STEP 5  BUILD, SAVE
# =============================================================================
hdr("STEP 2  BUILD PLOTS")
plots <- list()
for (sx in c("female", "male")) {
  r <- ml[[sx]]; sexlab <- tools::toTitleCase(sx); cols <- pal[[sx]]
  plots[[paste0("lasso_", sx)]] <- lasso_cv_plot(r, sexlab, cols["dark"])
  plots[[paste0("rf_",    sx)]] <- rf_importance_plot(r, sexlab, cols["light"], cols["dark"])
  plots[[paste0("svm_",   sx)]] <- svmrfe_curve_plot(r, sexlab, cols["dark"])
  say("%s: LASSO %d genes | RF %d genes | SVM-RFE %d genes",
      sexlab, length(r$sets$LASSO), length(r$sets$RandomForest), length(r$sets$SVM_RFE))
}

hdr("STEP 3  SAVE FIGURES")
save_fig(plots$lasso_female, "FIG_G2_05_lasso_cv_female",        5.5, 4.6)
save_fig(plots$lasso_male,   "FIG_G2_05_lasso_cv_male",          5.5, 4.6)
save_fig(plots$rf_female,    "FIG_G2_06_rf_importance_female",   5.8, 6.6)
save_fig(plots$rf_male,      "FIG_G2_06_rf_importance_male",     5.8, 5.6)
save_fig(plots$svm_female,   "FIG_G2_07_svmrfe_curve_female",    5.5, 4.6)
save_fig(plots$svm_male,     "FIG_G2_07_svmrfe_curve_male",      5.5, 4.6)

composite <-
  (plots$lasso_female | plots$rf_female | plots$svm_female) /
  (plots$lasso_male   | plots$rf_male   | plots$svm_male) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13)))
save_fig(composite, "FIG_G2_08_feature_selection_diagnostics_composite", 15.5, 11)

cat("\nDONE\n")
