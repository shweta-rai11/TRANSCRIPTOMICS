#!/usr/bin/env Rscript
# =============================================================================
# 06b_wgcna_composite.R -- assemble the WGCNA network-construction summary
# figure from panels already saved by 06_WGCNA.R:
#   top-left    fig_wgcna_01_gene_variance.png     (gene variance + filter cutoff)
#   top-right   fig_wgcna_02_sample_clustering.png (sample clustering + outlier line)
#   bottom      fig_wgcna_03_soft_threshold.png    (scale independence + mean connectivity)
# Panels are read as images and composited with cowplot/magick.
#
# Output: results/figures/fig_wgcna_soft_threshold_composite.{png,pdf}
#         results/wgcna1.png (copy, for direct sharing)
# Size  : 10 x 9.5 in @ 300 dpi -> 3000 x 2850 px
# =============================================================================
suppressMessages({library(cowplot); library(magick)})
fig <- "results/figures"

P <- function(f) {
  p <- file.path(fig, f)
  if (!file.exists(p)) stop("missing panel: ", p)
  ggdraw() + draw_image(p)
}

A <- P("fig_wgcna_01_gene_variance.png")
B <- P("fig_wgcna_02_sample_clustering.png")
C <- P("fig_wgcna_03_soft_threshold.png")

# rel_widths/rel_heights chosen from each panel's native pixel aspect ratio
# (1100x800, 1400x650, 1300x620) so the composited row heights roughly match
# the image content instead of leaving large letterboxed gaps.
top <- plot_grid(A, B, ncol = 2, rel_widths = c(0.64, 1))
composite <- plot_grid(top, C, ncol = 1, rel_heights = c(0.6, 1))

ggsave2(file.path(fig, "fig_wgcna_soft_threshold_composite.png"), composite,
        width = 10, height = 9.5, dpi = 300, bg = "white")
ggsave2(file.path(fig, "fig_wgcna_soft_threshold_composite.pdf"), composite,
        width = 10, height = 9.5, bg = "white")

file.copy(file.path(fig, "fig_wgcna_soft_threshold_composite.png"),
          "results/wgcna1.png", overwrite = TRUE)

cat("Wrote fig_wgcna_soft_threshold_composite.{png,pdf} and results/wgcna1.png\n")
