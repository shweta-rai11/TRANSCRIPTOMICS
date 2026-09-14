#!/usr/bin/env Rscript
# ROC figures (Train vs Test, per sex x algorithm x tissue setting) for the five-algorithm comparison, plus a composite AUC dot plot.
suppressMessages({ library(data.table); library(ggplot2) })
procN <- "data/processed/new"; fig <- "results/figures/new"; tab <- "results/tables"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(procN, "ml_algo_roc.rds"))
perf_same  <- fread(file.path(tab, "ML_performance_sametissue.csv"))
perf_cross <- fread(file.path(tab, "ML_performance_crosstissue.csv"))

ALGO_TAG <- c("Logistic regression" = "logreg", "SVM (RBF)" = "svm", "k-NN" = "knn",
              "Random forest" = "rf", "ANN (nnet)" = "ann")
COLS <- c("Train" = "#C0392B", "Test" = "#1F3B99")   # train = red, test = blue (fig_syn_panel_roc.png palette)

# clean up the earlier faceted (superseded) figures if present
stale <- as.vector(outer(sprintf("fig_ml_%s_%s", rep(ALGO_TAG, each = 2), c("sametissue", "crosstissue")),
                          c("png", "pdf"), FUN = function(a, b) paste0(a, ".", b)))
unlink(file.path(fig, stale))

# single-panel Train(red)/Test(blue) ROC, one sex at a time
plot_one <- function(roc_dt, perf, algo, sexlab, test_label, group_word) {
  d <- roc_dt[algorithm == algo & sex == sexlab]
  d[, curve := ifelse(dataset == "Train (resampled CV)", "Train", "Test")]
  d[, curve := factor(curve, levels = c("Train", "Test"))]

  row_tr <- perf[sex == sexlab & algorithm == algo & dataset == "Train (resampled CV)"]
  row_te <- perf[sex == sexlab & algorithm == algo & dataset == test_label]
  sep_tr <- grepl("SEPARATION", row_tr$AUC_CI); sep_te <- grepl("SEPARATION", row_te$AUC_CI)
  auc_num <- function(s) as.numeric(sub("^([0-9.]+).*", "\\1", s))
  lab_tr <- sprintf("Train (resampled CV) (AUC=%.3f%s)", auc_num(row_tr$AUC_CI), if (sep_tr) ", SEPARATION" else "")
  lab_te <- sprintf("%s (AUC=%.3f%s)", test_label, auc_num(row_te$AUC_CI), if (sep_te) ", SEPARATION" else "")
  d[, lab := ifelse(curve == "Train", lab_tr, lab_te)]
  d[, lab := factor(lab, levels = c(lab_tr, lab_te))]

  cap <- sprintf("Train n=%d (%d RA/%d ctrl) | Test n=%d (%d RA/%d ctrl)%s",
                  row_tr$n, row_tr$n_RA, row_tr$n_ctrl, row_te$n, row_te$n_RA, row_te$n_ctrl,
                  if (sep_tr || sep_te) " -- SEPARATION = small-n artefact, not a precision estimate." else "")

  ggplot(d, aes(1 - spec, sens, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey75") +
    geom_path(linewidth = 1.1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_colour_manual(values = setNames(unname(COLS[c("Train","Test")]), c(lab_tr, lab_te)), name = NULL) +
    coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity",
         subtitle = sprintf("%s ROC, %s (%s)", algo, sexlab, group_word),
         caption = cap) +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 13),
          plot.caption = element_text(size = 7.5, colour = "grey40", hjust = 0))
}

save_pair <- function(g, path_noext) {
  ggsave(paste0(path_noext, ".png"), g, width = 6, height = 6, dpi = 300)
  ggsave(paste0(path_noext, ".pdf"), g, width = 6, height = 6)
}

for (algo in names(ALGO_TAG)) {
  tag <- ALGO_TAG[[algo]]
  for (sexlab in c("Female", "Male")) {
    g1 <- plot_one(obj$sametissue$roc, perf_same, algo, sexlab, "Blood test (internal holdout)", "blood")
    save_pair(g1, file.path(fig, sprintf("fig_ml_%s_sametissue_%s", tag, tolower(sexlab))))

    g2 <- plot_one(obj$crosstissue$roc, perf_cross, algo, sexlab, "Synovium test (cross-tissue)", "synovium")
    save_pair(g2, file.path(fig, sprintf("fig_ml_%s_crosstissue_%s", tag, tolower(sexlab))))
  }
  cat(sprintf("wrote fig_ml_%s_{sametissue,crosstissue}_{female,male}.{png,pdf}\n", tag))
}

# composite AUC comparison dot plot (train=red, test=blue)
parse_auc <- function(dt, setting_lab) {
  dt <- copy(dt)
  dt[, auc := as.numeric(sub("^([0-9.]+).*", "\\1", AUC_CI))]
  dt[, lo  := as.numeric(sub(".*\\(([0-9.]+)-.*", "\\1", AUC_CI))]
  dt[, hi  := as.numeric(sub(".*-([0-9.]+)\\).*", "\\1", AUC_CI))]
  dt[, sep := grepl("SEPARATION", AUC_CI)]
  dt[, setting := setting_lab]
  dt
}
comp <- rbindlist(list(
  parse_auc(perf_same,  "Same-tissue (blood to blood)"),
  parse_auc(perf_cross, "Cross-tissue (blood to synovium)")
))
comp[, dataset2 := ifelse(dataset == "Train (resampled CV)", "Train", "Test")]
comp[, algorithm := factor(algorithm, levels = rev(names(ALGO_TAG)))]
comp[, sex := factor(sex, levels = c("Female", "Male"))]

gcomp <- ggplot(comp, aes(x = auc, y = algorithm, colour = dataset2, shape = sep)) +
  geom_vline(xintercept = 0.5, linetype = 3, colour = "grey70") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15,
                position = position_dodge(width = 0.5), orientation = "y") +
  geom_point(aes(size = sep), position = position_dodge(width = 0.5)) +
  facet_grid(sex ~ setting) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 2.0, `FALSE` = 2.4), guide = "none") +
  labs(x = "AUC (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig, "fig_ml_algorithm_comparison_auc.png"), gcomp, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(fig, "fig_ml_algorithm_comparison_auc.pdf"), gcomp, width = 9, height = 5.5)
cat("wrote fig_ml_algorithm_comparison_auc.{png,pdf}\n")
