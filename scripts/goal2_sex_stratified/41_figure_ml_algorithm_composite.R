#!/usr/bin/env Rscript
# =============================================================================
# 41_figure_ml_algorithm_composite.R
# -----------------------------------------------------------------------------
# ONE composite figure per sex for the five-algorithm diagnostic comparison
# (37_/38_/38b_): a sub-figure per algorithm (logistic regression, SVM-RBF,
# k-NN, random forest, ANN), each showing all FOUR evaluation settings on the
# same axes -- Train (resampled CV) | Internal test (blood) | External test
# (blood, GSE15573) | Synovium (cross-tissue) -- so algorithm and tissue/
# validation-stage comparisons can both be read off one figure. Female and
# male are always separate files (never pooled).
#
# No plot title/subtitle text is drawn on either figure (captions are written
# separately, matching the convention already used in 28_figure_final_roc.R).
# The facet strip above each sub-panel names only the algorithm, which is
# structural (it is how the sub-figures are told apart), not a title.
#
# Outputs:
#   results/figures/new/fig_ml_composite_{female,male}.{png,pdf}
# =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
procN <- "data/processed/new"; fig <- "results/figures/new"; tab <- "results/tables"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(procN, "ml_algo_roc.rds"))

## ---- relabel each setting's dataset factor to a short, shared vocabulary --
lab_map <- c("Train (resampled CV)" = "Train",
             "Blood test (internal holdout)" = "Internal test (blood)",
             "External blood test (GSE15573)" = "External test (blood)",
             "Synovium test (cross-tissue)" = "Synovium")
ord <- c("Train", "Internal test (blood)", "External test (blood)", "Synovium")
pal <- c("Train" = "#1b6ca8", "Internal test (blood)" = "#e08214",
         "External test (blood)" = "#4d9221", "Synovium" = "#c0392b")

pick <- function(dt, keep_dataset) {
  d <- dt[dataset %in% keep_dataset]
  d[, dataset := lab_map[dataset]]
  d
}
roc_all <- rbindlist(list(
  pick(obj$sametissue$roc,    c("Train (resampled CV)", "Blood test (internal holdout)")),
  pick(obj$externalblood$roc, "External blood test (GSE15573)"),
  pick(obj$crosstissue$roc,   "Synovium test (cross-tissue)")
), use.names = TRUE)
roc_all[, dataset := factor(dataset, levels = ord)]

## per (sex, algorithm, dataset) AUC, for small in-panel annotations (not titles)
auc_dt <- unique(roc_all[, .(sex, algorithm, dataset, auc)])
auc_dt[, lab := sprintf("%s: %.3f", dataset, auc)]
setorder(auc_dt, sex, algorithm, dataset)
auc_dt[, y := 0.26 - 0.065 * (as.integer(dataset) - 1), by = .(sex, algorithm)]

ALGO_ORD <- c("Logistic regression", "SVM (RBF)", "k-NN", "Random forest", "ANN (nnet)")
roc_all[, algorithm := factor(algorithm, levels = ALGO_ORD)]
auc_dt[, algorithm  := factor(algorithm, levels = ALGO_ORD)]

plot_sex <- function(sexlab) {
  d <- roc_all[sex == sexlab]; a <- auc_dt[sex == sexlab]
  ggplot(d, aes(x = 1 - spec, y = sens, colour = dataset)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey75") +
    geom_path(linewidth = 0.9) +
    geom_text(data = a, aes(x = 0.98, y = y, label = lab, colour = dataset),
              hjust = 1, size = 2.7, show.legend = FALSE) +
    facet_wrap(~ algorithm, ncol = 3) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0), breaks = c(0, 0.5, 1)) +
    scale_colour_manual(values = pal, breaks = ord, name = NULL) +
    coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold"))
}

for (sexlab in c("Female", "Male")) {
  g <- plot_sex(sexlab); tag <- tolower(sexlab)
  ggsave(file.path(fig, sprintf("fig_ml_composite_%s.png", tag)), g, width = 9.5, height = 7, dpi = 300)
  ggsave(file.path(fig, sprintf("fig_ml_composite_%s.pdf", tag)), g, width = 9.5, height = 7)
  cat(sprintf("wrote fig_ml_composite_%s.{png,pdf}\n", tag))
}
