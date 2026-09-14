#!/usr/bin/env Rscript
# =============================================================================
# 03g_data_preprocessing_composite.R -- assemble the data-preprocessing summary
# figure (A-D) from panels already saved by earlier 03_* scripts:
#   A. fig_combine_two_datasets.png       <- 03_normalize_batch_figure.R
#   B. fig_combine_pca_combat.png         <- 03_normalize_batch_figure.R
#   C. fig_dataset_gene_overlap_venn.png  <- 03_dataset_gene_overlap_venn.R
#   D. fig_train_test_split_bar.png       <- 03f_train_test_split_barchart.R
# Panels are read as images and composited with cowplot/magick (B and C are
# grid/eulerr graphics, so image-compositing is the robust route).
#
# Output: results/figures/fig_data_preprocessing_composite.{png,pdf}
#         results/DataPreprocessing.png (copy, for direct sharing)
# Size  : 10 x 9.5 in @ 300 dpi -> 3000 x 2850 px
# =============================================================================
suppressMessages({library(cowplot); library(magick)})
fig <- "results/figures"

P <- function(f) {
  p <- file.path(fig, f)
  if (!file.exists(p)) stop("missing panel: ", p)
  ggdraw() + draw_image(p)
}

# Each source panel already prints its own text (e.g. "GSE93272") right at the
# top-left corner, which collides with a plot_grid panel label placed there.
# Pad a blank strip above each panel first, so the A)/B)/C)/D) label lands on
# blank space instead of on top of the panel's own title text.
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

file.copy(file.path(fig, "fig_data_preprocessing_composite.png"),
          "results/DataPreprocessing.png", overwrite = TRUE)

cat("Wrote fig_data_preprocessing_composite.{png,pdf} and results/DataPreprocessing.png\n")
