#!/usr/bin/env Rscript
# =============================================================================
# 12d_feature_selection_diagnostics_extra.R
# -----------------------------------------------------------------------------
# Additional algorithm-level diagnostic figures for the feature-selection step
# (12_feature_selection.R), one file per plot per sex (not composited), in the
# style used for conference/thesis figure decks: LASSO coefficient shrinkage
# path, LASSO coefficient direction at lambda.min, Random Forest OOB error vs
# ntree, and the SVM-RFE accuracy curve (companion to the error curve already
# drawn in 12c). All are built from the fitted objects already saved by 12 —
# no models are re-run.
#
# NOTE ON METHOD SCOPE
#   This pipeline's selectors are LASSO / Random Forest / SVM-RFE only; Boruta
#   is explicitly NOT used here (see 12_feature_selection.R header: the
#   candidate set is already small and causally pre-filtered by MR, so the
#   wrapper Boruta test was omitted in favour of the three-algorithm design).
#   No Boruta figure is produced — one would misrepresent what was run.
#
#   in : data/processed/new/ml_features.rds
#   out: results/figures/FIG_G2_09_lasso_coefpath_{female,male}.png/.pdf
#        results/figures/FIG_G2_10_lasso_coefdirection_{female,male}.png/.pdf
#        results/figures/FIG_G2_11_rf_oob_error_{female,male}.png/.pdf
#        results/figures/FIG_G2_12_svmrfe_accuracy_{female,male}.png/.pdf
# =============================================================================
suppressMessages({
  library(glmnet); library(randomForest); library(ggplot2)
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
dir_pal <- c(up = "#E8792A", down = "#3B6FA8")   # up in RA / down in RA

back_map <- function(r) setNames(r$mr_genes, make.names(r$mr_genes))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  say("  wrote %s.{png,pdf}", name)
}

# =============================================================================
# STEP 2  LASSO coefficient shrinkage path
# =============================================================================
lasso_coefpath_plot <- function(r, sexlab, dark) {
  bm   <- back_map(r)
  beta <- as.matrix(r$cv$glmnet.fit$beta)                 # genes x lambda steps
  lam  <- r$cv$glmnet.fit$lambda
  df <- do.call(rbind, lapply(seq_len(nrow(beta)), function(i) {
    data.frame(gene = unname(bm[rownames(beta)[i]]), loglambda = log(lam),
               coef = beta[i, ])
  }))
  df$selected <- df$gene %in% r$sets$LASSO
  sel_genes <- sort(unique(df$gene[df$selected]))
  cols <- setNames(colorRampPalette(c(dark, "#F2A65A", "#5E9C6E", "#3B6FA8",
                                       "#8E44AD", "#C0392B"))(length(sel_genes)),
                    sel_genes)
  ggplot(df, aes(loglambda, coef, group = gene)) +
    geom_line(data = df[!df$selected, ], colour = "grey80", linewidth = 0.35) +
    geom_line(data = df[df$selected, ], aes(colour = gene), linewidth = 0.8) +
    geom_vline(xintercept = log(r$cv$lambda.min), linetype = "dashed", colour = "grey30") +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(values = cols, name = "Selected gene\n(LASSO)") +
    labs(subtitle = "coloured = genes retained at lambda.min; grey = shrunk to zero",
         x = "log(lambda)", y = "Coefficient") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(b = 10)),
          legend.text = element_text(size = 7.5), legend.title = element_text(size = 8.5))
}

# =============================================================================
# STEP 3  LASSO coefficient direction at lambda.min
# =============================================================================
lasso_coefdirection_plot <- function(r, sexlab) {
  bm <- back_map(r)
  co <- coef(r$cv, s = "lambda.min")[-1, 1]
  co <- co[co != 0]
  df <- data.frame(gene = unname(bm[names(co)]), coef = as.numeric(co))
  df <- df[order(df$coef), ]
  df$gene <- factor(df$gene, levels = df$gene)
  df$direction <- ifelse(df$coef > 0, "up", "down")
  ggplot(df, aes(gene, coef, fill = direction)) +
    geom_col(width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = dir_pal, labels = c(up = "up in RA", down = "down in RA"),
                       name = NULL) +
    labs(subtitle = sprintf("%d genes retained at lambda.min = %.4f (control = reference)",
                             nrow(df), r$cv$lambda.min),
         x = NULL, y = "LASSO coefficient (log-odds scale)") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(b = 10)),
          legend.position = "top")
}

# =============================================================================
# STEP 4  Random Forest OOB error vs ntree
# =============================================================================
rf_oob_plot <- function(r, sexlab, dark) {
  dark <- unname(dark)
  er <- r$rf$err.rate
  df <- data.frame(ntree = seq_len(nrow(er)), OOB = er[, "OOB"],
                    HC = er[, "HC"], RA = er[, "RA"])
  dfl <- reshape(df, direction = "long", varying = c("OOB", "HC", "RA"),
                 v.names = "error", timevar = "series",
                 times = c("OOB", "HC", "RA"))
  series_cols <- setNames(c("grey25", "#3B6FA8", dark), c("OOB", "HC", "RA"))
  ggplot(dfl, aes(ntree, error, colour = series)) +
    geom_line(linewidth = 0.55) +
    scale_colour_manual(values = series_cols, name = NULL) +
    labs(subtitle = sprintf("ntree = %d, tuned mtry = %d | final OOB error = %.3f",
                             r$rf$ntree, r$tuned$rf_mtry, tail(er[, "OOB"], 1)),
         x = "Number of trees", y = "Error rate") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(b = 10)),
          legend.position = "top")
}

# =============================================================================
# STEP 5  SVM-RFE accuracy curve (companion to the 12c error curve)
# =============================================================================
svmrfe_accuracy_plot <- function(r, sexlab, dark) {
  cu <- r$svm_curve
  df <- data.frame(k = cu$k, acc = 1 - cu$err)
  best_pt <- df[df$k == cu$best, ]
  ggplot(df, aes(k, acc)) +
    geom_line(colour = dark, linewidth = 0.5) +
    geom_point(colour = dark, size = 1.2) +
    geom_vline(xintercept = cu$best, linetype = "dashed", colour = "grey30") +
    geom_point(data = best_pt, aes(k, acc), shape = 21, colour = "black",
               fill = dark, size = 2.6) +
    annotate("text", x = best_pt$k, y = best_pt$acc,
             label = sprintf("n=%d (%.3f)", best_pt$k, best_pt$acc),
             vjust = -0.9, size = 3, colour = "grey20") +
    labs(subtitle = sprintf("tuned cost = %s | CV-optimal size = %d genes (accuracy = %.3f)",
                             format(r$tuned$svm_cost, trim = TRUE), cu$best, 1 - cu$besterr),
         x = "Number of top-ranked features (k)", y = "10-fold CV accuracy") +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.2, margin = margin(t = 2, b = 10)))
}

# =============================================================================
# STEP 6  BUILD, SAVE
# =============================================================================
hdr("STEP 2  BUILD + SAVE FIGURES")
for (sx in c("female", "male")) {
  r <- ml[[sx]]; sexlab <- tools::toTitleCase(sx); cols <- pal[[sx]]

  save_fig(lasso_coefpath_plot(r, sexlab, cols["dark"]),
           sprintf("FIG_G2_09_lasso_coefpath_%s", sx), 6.4, 4.8)
  save_fig(lasso_coefdirection_plot(r, sexlab),
           sprintf("FIG_G2_10_lasso_coefdirection_%s", sx), 5.6, 4.2)
  save_fig(rf_oob_plot(r, sexlab, cols["dark"]),
           sprintf("FIG_G2_11_rf_oob_error_%s", sx), 5.6, 4.2)
  save_fig(svmrfe_accuracy_plot(r, sexlab, cols["dark"]),
           sprintf("FIG_G2_12_svmrfe_accuracy_%s", sx), 5.5, 4.6)
}

cat("\nDONE\n")
