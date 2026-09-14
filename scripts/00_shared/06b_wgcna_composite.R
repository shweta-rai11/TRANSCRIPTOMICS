#!/usr/bin/env Rscript
# Assemble the WGCNA network-construction summary figure from panels saved by 06_WGCNA.R
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

# rel_widths/rel_heights chosen to match each panel's native pixel aspect ratio
top <- plot_grid(A, B, ncol = 2, rel_widths = c(0.64, 1))
composite <- plot_grid(top, C, ncol = 1, rel_heights = c(0.6, 1))

ggsave2(file.path(fig, "fig_wgcna_soft_threshold_composite.png"), composite,
        width = 10, height = 9.5, dpi = 300, bg = "white")
ggsave2(file.path(fig, "fig_wgcna_soft_threshold_composite.pdf"), composite,
        width = 10, height = 9.5, bg = "white")

cat("Wrote fig_wgcna_soft_threshold_composite.{png,pdf}\n")
