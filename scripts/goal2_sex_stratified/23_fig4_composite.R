#!/usr/bin/env Rscript
# =============================================================================
# 23_fig4_composite.R  -- assemble the RA Fig 4 composite from the saved panels
# (A expression, B discovery logFC, C clustered heatmap) into one labelled
# figure. Panels are read as images (C is a ComplexHeatmap grid graphic, so
# image-compositing is the robust route).
#
# NOTE: former panels D and E (CIBERSORT stacked bar + boxplot) were REMOVED.
# They contrasted FEMALE RA vs MALE RA -- a between-sex (sex-SPECIFIC) test.
# This chapter reports SEX-STRATIFIED results only, so Fig 4 is now A-C.
#
# Inputs  (results/figures/new/):
#   fig_mr_gene_expression_groups.png     <- 19_fig4A_expression_groups.R
#   fig_mr_panelB_logfc_discovery.png     <- 20_fig4BC_logfc_heatmaps.R
#   fig_mr_panelC_clustered_heatmap.png   <- 22_fig4C_clustered_heatmap.R
# Output: results/figures/new/fig_mr_NEW_FIG4_composite.{png,pdf}
# =============================================================================
suppressMessages({library(cowplot); library(magick)})
fig <- "results/figures/new"
P <- function(f) {
  p <- file.path(fig, f)
  if (!file.exists(p)) stop("missing panel: ", p)
  ggdraw() + draw_image(p)
}

A <- P("fig_mr_gene_expression_groups.png")
B <- P("fig_mr_panelB_logfc_discovery.png")
C <- P("fig_mr_panelC_clustered_heatmap.png")

row1 <- plot_grid(A, labels = "A", label_size = 20)
row2 <- plot_grid(B, C, labels = c("B", "C"), label_size = 20,
                  rel_widths = c(1, 1.15), nrow = 1)

composite <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1.0, 0.95))

ggsave2(file.path(fig, "fig_mr_NEW_FIG4_composite.png"), composite,
        width = 13, height = 9, dpi = 220, bg = "white")
ggsave2(file.path(fig, "fig_mr_NEW_FIG4_composite.pdf"), composite,
        width = 13, height = 9, bg = "white")
cat("wrote fig_mr_NEW_FIG4_composite.{png,pdf} (panels A-C; D/E removed)\n")
