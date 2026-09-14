#!/usr/bin/env Rscript
# ROC curves for the final recommended MR-anchored elastic-net panels, per sex, across train (nested CV), internal test, and external blood.
suppressMessages({library(data.table); library(ggplot2)})
proc <- "data/processed"; fig <- "results/figures/new"
obj <- readRDS("data/processed/new/mr_final_objects.rds")
roc <- as.data.table(obj$roc); summ <- as.data.table(obj$summary)

ord <- c("Train (nested CV)", "Internal test", "External blood")
roc[, dataset := factor(dataset, levels = ord)]
pal <- c("Train (nested CV)" = "#1b6ca8", "Internal test" = "#e08214",
         "External blood" = "#4d9221")

# AUC label per dataset (from the saved summary, so CI matches the tables)
auc_lab <- function(sexlab) {
  s <- summ[sex == sexlab]
  strip_sep <- function(x) sub(" SEPARATION$", "", x)
  setNames(c(paste0("Train (nested CV): AUC ", sub(" \\(", "\n(", strip_sep(s$nested_CV))),
             paste0("Internal test: AUC ",     sub(" \\(", "\n(", strip_sep(s$internal_test))),
             paste0("External blood: AUC ",    sub(" \\(", "\n(", strip_sep(s$external_blood)))), ord)
}

plot_sex <- function(sexlab) {
  d <- roc[sex == sexlab]
  labs <- auc_lab(sexlab)
  d[, lab := labs[as.character(dataset)]]
  d[, lab := factor(lab, levels = labs[ord])]
  ggplot(d, aes(x = 1 - spec, y = sens, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
    geom_path(linewidth = 1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_colour_manual(values = setNames(unname(pal[ord]), labs[ord]), name = NULL) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = "grey80"),
          legend.margin = margin(4, 6, 4, 6),
          legend.key.height = unit(2, "lines"),
          legend.key.spacing.y = unit(4, "pt"),
          legend.text = element_text(size = 9, lineheight = 0.95),
          panel.grid.minor = element_blank())
}

for (sx in c("Female", "Male")) {
  g <- plot_sex(sx); tag <- tolower(sx)
  ggsave(file.path(fig, sprintf("fig_mr_final_roc_%s.png", tag)), g,
         width = 6, height = 6, dpi = 300)
  ggsave(file.path(fig, sprintf("fig_mr_final_roc_%s.pdf", tag)), g,
         width = 6, height = 6)
  cat(sprintf("wrote fig_mr_final_roc_%s.{png,pdf}\n", tag))
}
