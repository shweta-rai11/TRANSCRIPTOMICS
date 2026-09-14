#!/usr/bin/env Rscript
# Assembles the data-preprocessing summary figure (A-D) from panels saved by earlier 03_* scripts
suppressMessages({library(cowplot); library(magick)})
fig <- "results/figures"

P <- function(f) {
  p <- file.path(fig, f)
  if (!file.exists(p)) stop("missing panel: ", p)
  ggdraw() + draw_image(p)
}

# Pad a blank strip above each panel so the A)/B)/C)/D) label doesn't collide with panel titles
pad_top <- function(p, frac = 0.07) plot_grid(NULL, p, ncol = 1, rel_heights = c(frac, 1 - frac))

A <- pad_top(P("fig_combine_two_datasets.png"))
B <- pad_top(P("fig_combine_pca_combat.png"))
C <- pad_top(P("fig_dataset_gene_overlap_venn.png"))
D <- pad_top(P("fig_train_test_split_bar.png"))

composite <- plot_grid(A, B, C, D,
                        labels = c("A)", "B)", "C)", "D)"),
                        label_size = 20, label_fontface = "bold",
                        label_x = 0.01, label_y = 0.99,
                        ncol = 2, nrow = 2,
                        rel_heights = c(1.05, 0.95))

ggsave2(file.path(fig, "fig_data_preprocessing_composite.png"), composite,
        width = 10, height = 9.5, dpi = 300, bg = "white")
ggsave2(file.path(fig, "fig_data_preprocessing_composite.pdf"), composite,
        width = 10, height = 9.5, bg = "white")

cat("Wrote fig_data_preprocessing_composite.{png,pdf}\n")
