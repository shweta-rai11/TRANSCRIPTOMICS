#!/usr/bin/env Rscript
# =============================================================================
# 06c_wgcna_sex_dendro_moduletrait_composite.R -- assemble, for each sex, the
# gene dendrogram + module-trait heatmap into one summary figure, from panels
# already saved by 06_WGCNA.R and 08_module_trait_RA_control.R:
#   top    fig_wgcna_1{3,4}_dendro_{female,male}.png       (gene dendrogram)
#   bottom fig_wgcna_module_trait_RAvsControl_{female,male}.png (module-trait heatmap)
#
# The dendrogram (native 1300x600, wide-short) and the heatmap (native 640x950,
# narrow-tall) have very different aspect ratios. Forcing both into a shared
# row height (or stretching either to fill a box) distorts the panel -- text
# and dendrogram branches come out visibly elongated. Instead each panel is
# scaled UNIFORMLY (aspect preserved exactly) into a stacked layout: the
# dendrogram spans the full canvas width at its native aspect, and the
# heatmap is centered underneath at its native aspect, sized to exactly fill
# the remaining height. No raster distortion, and the canvas size is
# identical for both sexes.
#
# Output: results/figures/fig_wgcna_{sex}_dendro_moduletrait_composite.{png,pdf}
#         results/wgcna1_{sex}.png (copy, for direct sharing)
# Size  : 10 x 9.5 in @ 300 dpi -> 3000 x 2850 px (same for male and female)
# =============================================================================
suppressMessages({library(cowplot); library(magick)})
fig <- "results/figures"

w_px <- 3000; h_px <- 2850
dendro_ar  <- 1300 / 600   # native dendrogram aspect (width / height)
heatmap_ar <- 640 / 950    # native heatmap aspect (width / height)

h_top    <- w_px / dendro_ar          # full-width dendrogram row height
h_bottom <- h_px - h_top              # remaining row height for the heatmap
w_bottom <- h_bottom * heatmap_ar     # heatmap width at that height, native aspect

build_composite <- function(sex, dendro_file, heatmap_file) {
  Pd <- file.path(fig, dendro_file); Ph <- file.path(fig, heatmap_file)
  if (!file.exists(Pd)) stop("missing panel: ", Pd)
  if (!file.exists(Ph)) stop("missing panel: ", Ph)

  top    <- ggdraw() + draw_image(Pd)  # cell aspect == dendro_ar exactly: no distortion
  bottom <- ggdraw() + draw_image(Ph, width = w_bottom / w_px, height = 1,
                                   x = 0.5, hjust = 0.5)  # centered, native aspect

  composite <- plot_grid(top, bottom, ncol = 1, rel_heights = c(h_top, h_bottom))

  out_png <- sprintf("fig_wgcna_%s_dendro_moduletrait_composite.png", sex)
  out_pdf <- sprintf("fig_wgcna_%s_dendro_moduletrait_composite.pdf", sex)
  ggsave2(file.path(fig, out_png), composite, width = 10, height = 9.5, dpi = 300, bg = "white")
  ggsave2(file.path(fig, out_pdf), composite, width = 10, height = 9.5, bg = "white")

  file.copy(file.path(fig, out_png), sprintf("results/wgcna1_%s.png", sex), overwrite = TRUE)
  cat(sprintf("Wrote %s, %s and results/wgcna1_%s.png\n", out_png, out_pdf, sex))
}

build_composite("male",   "fig_wgcna_14_dendro_male.png",   "fig_wgcna_module_trait_RAvsControl_male.png")
build_composite("female", "fig_wgcna_13_dendro_female.png", "fig_wgcna_module_trait_RAvsControl_female.png")
