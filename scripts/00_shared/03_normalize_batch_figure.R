#!/usr/bin/env Rscript
# =============================================================================
# fig_two_datasets_boxplot.R
# Per-sample expression boxplots, BEFORE vs AFTER quantile normalization,
# one row per dataset (GSE93272 on top, GSE110169 below), all samples shown.
# Reads data/processed/combined_train.rds; writes PNG + PDF to results/figures/.
# =============================================================================
library(ggplot2)
library(reshape2)
library(patchwork)

# ---- I/O --------------------------------------------------------------------
proc <- "data/processed"
fig  <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

o    <- readRDS(file.path(proc, "combined_train.rds"))
meta <- o$meta
pre  <- o$expr_prenorm    # merged, before normalization
qn   <- o$expr_qnorm      # after quantile normalization (before ComBat)
cb   <- o$expr            # after quantile normalization + ComBat batch correction

# Per-dataset sample counts for the 70% TRAINING set (computed, never hardcoded,
# so subtitles stay correct if the split or cohort changes).
ds_tab   <- table(meta$dataset)
ds_label <- paste(sprintf("%s (n = %d)", names(ds_tab), as.integer(ds_tab)),
                  collapse = " and ")   # e.g. "GSE110169 (n = 84) and GSE93272 (n = 99)"

# save helper: writes both a 300-dpi PNG and a vector PDF (preferred for print)
save2 <- function(g, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), g, width = w, height = h, dpi = 300)
  ggsave(file.path(fig, paste0(name, ".pdf")), g, width = w, height = h, device = cairo_pdf)
}

# ---- build the long-format table (ALL samples) ------------------------------
s <- seq_len(ncol(pre))   # every training sample (no downsampling)

if (!is.matrix(qn) && !is.data.frame(qn)) {
  stop("Expected expr_qnorm to be a numeric matrix in combined_train.rds")
}

box_df <- function(mat, stage) {
  d <- melt(mat[, s]); names(d) <- c("gene", "sample", "value")
  d$stage   <- stage
  d$dataset <- meta$dataset[match(d$sample, meta$sample)]
  d
}

bn <- rbind(box_df(pre, "BEFORE normalization"),
            box_df(qn,  "AFTER quantile normalization"))
bn$stage <- factor(bn$stage,
                   levels = c("BEFORE normalization", "AFTER quantile normalization"))

# ---- colors: grey = before, teal = after (colorblind-safe) ------------------
stage_cols <- c("BEFORE normalization"         = "#BDBDBD",
                "AFTER quantile normalization" = "#00838F")

# ---- one before/after facet plot per dataset --------------------------------
make_ds_plot <- function(ds) {
  d <- bn[bn$dataset == ds, ]
  d$sample <- droplevels(factor(d$sample))   # drop the other dataset's ghost samples
  ggplot(d, aes(x = reorder(sample, value, median), y = value, fill = stage)) +
    geom_boxplot(outlier.size = .1, linewidth = .15) +
    facet_wrap(~ stage, ncol = 2) +          # before | after
    scale_fill_manual(values = stage_cols, guide = "none") +
    labs(title = ds, x = "Sample", y = expression(log[2]~expression)) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x      = element_blank(),    # 100+ samples: names unreadable, so blank them
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_text(size = 9),
      axis.title.y     = element_text(size = 9),
      plot.title       = element_text(face = "bold", size = 11),
      strip.text       = element_text(size = 9),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      panel.grid.minor = element_blank()
    )
}

# ---- stack the two datasets into one figure ("/" = vertical) ----------------
p_two <- make_ds_plot("GSE93272") / make_ds_plot("GSE110169") +
  plot_annotation(
    title    = "Per-sample expression: before vs after normalization",
    subtitle = paste0(ds_label, ", quantile normalization aligns per-sample scales"),
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9)))

save2(p_two, "fig_combine_two_datasets", 7.2, 10)   # double-column width, tall

cat("Wrote: fig_combine_two_datasets (PNG + PDF)\n")



##################################### yet to be done 





# ---- 1b. Density plots before vs after normalization -----------------------
density_df <- rbind(box_df(pre, "BEFORE normalization"), box_df(qn, "AFTER quantile normalization"))
p_density <- ggplot(density_df, aes(value, color = stage, linetype = dataset)) +
  geom_density(linewidth = 0.8, alpha = 0.2) +
  labs(title = "Distribution of expression values before vs after normalization",
       subtitle = "Quantile normalization makes the two datasets more comparable",
       x = "log2 expression", y = "density") +
  theme(legend.position = "top")
save2(p_density, "fig_combine_density_norm", 7, 4.5)



# ---- 4. PCA before vs after ComBat, WITH marginal densities -----------------
# Reference-style panels (Wang et al., "Managing batch effects"): a PC1/PC2
# scatter with a PC1 density strip on top and a PC2 density strip on the right,
# both split by study so the batch overlap is visible in the margins too.
# study colours (colourblind-safe): one hue per dataset
dcol <- c(GSE93272 = "#0072B2", GSE110169 = "#D55E00")

# publication theme: 7 pt base, thin lines, no clutter (journal-ready)
pub_theme <- theme_bw(base_size = 7) +
  theme(axis.title   = element_text(size = 7),
        axis.text    = element_text(size = 6, colour = "grey20"),
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        panel.border = element_rect(linewidth = .3, colour = "grey40"),
        axis.ticks   = element_line(linewidth = .3, colour = "grey40"))

# run PCA on one matrix; return the scores + the %variance carried by PC1/PC2.
pcp <- function(mat) {
  m <- mat[apply(mat, 1, sd) > 0, ]           # drop zero-variance genes (prcomp scaling)
  p <- prcomp(t(m), scale. = TRUE)
  v <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  list(df = data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2],
                       dataset = meta$dataset, group = meta$group),
       v1 = v[1], v2 = v[2])
}

# build one marginal-density panel (top density + scatter + right density).
# xlim / ylim keep the PC1 / PC2 axes identical across the before/after panels.
mk_panel <- function(pc, title, xlim, ylim) {
  d <- pc$df
  xlab <- sprintf("PC1: %.1f%% expl.var", pc$v1)
  ylab <- sprintf("PC2: %.1f%% expl.var", pc$v2)

  main <- ggplot(d, aes(PC1, PC2, color = dataset)) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = .25,
               colour = "grey45") +          # reference line through PC1 = 0
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = .25,
               colour = "grey45") +          # reference line through PC2 = 0
    stat_ellipse(aes(group = dataset), type = "norm", level = 0.9,
                 linewidth = .35, alpha = .6, show.legend = FALSE) +
    geom_point(aes(shape = group), size = 1.1, alpha = .8, stroke = .3) +
    scale_color_manual(values = dcol, name = "Study") +
    scale_shape_manual(values = c(HC = 16, RA = 17), name = "Group") +
    coord_cartesian(xlim = xlim, ylim = ylim) +   # shared PC1 + PC2 range
    labs(x = xlab, y = ylab) +
    pub_theme +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom")

  top <- ggplot(d, aes(PC1, fill = dataset, color = dataset)) +
    geom_density(alpha = .4, linewidth = .3, show.legend = FALSE) +
    scale_fill_manual(values = dcol) + scale_color_manual(values = dcol) +
    coord_cartesian(xlim = xlim) +          # match the scatter's PC1 range
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = .5, size = 8, face = "bold"),
          plot.margin = margin(2, 2, 0, 2))

  right <- ggplot(d, aes(PC2, fill = dataset, color = dataset)) +
    geom_density(alpha = .4, linewidth = .3, show.legend = FALSE) +
    scale_fill_manual(values = dcol) + scale_color_manual(values = dcol) +
    coord_flip(xlim = ylim) +               # match the scatter's PC2 range
    theme_void() + theme(plot.margin = margin(2, 2, 2, 0))

  # 2x2: top density | (corner)   //   scatter | right density
  top + plot_spacer() + main + right +
    plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
}

pc_before <- pcp(qn)
pc_after  <- pcp(cb)
# common PC1/PC2 ranges across both stages, padded 4%, so both axes are identical
pad <- function(r) r + c(-1, 1) * diff(r) * 0.04
xlim <- pad(range(c(pc_before$df$PC1, pc_after$df$PC1)))
ylim <- pad(range(c(pc_before$df$PC2, pc_after$df$PC2)))

panel_before <- mk_panel(pc_before, "Before ComBat (normalized)", xlim, ylim)
panel_after  <- mk_panel(pc_after,  "After ComBat",               xlim, ylim)

p_pca <- (panel_before | panel_after) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Batch correction by ComBat: PCA before vs after",
    subtitle = ds_label,
    theme    = theme(plot.title    = element_text(face = "bold", size = 9),
                     plot.subtitle = element_text(size = 7, colour = "grey30"))) &
  theme(legend.position = "bottom", legend.margin = margin(t = 0))

# double-column journal width (180 mm ~ 7.1 in); 600-dpi PNG + vector PDF
ggsave(file.path(fig, "fig_combine_pca_combat.png"), p_pca,
       width = 7.1, height = 3.9, dpi = 600)
ggsave(file.path(fig, "fig_combine_pca_combat.pdf"), p_pca,
       width = 7.1, height = 3.9, device = cairo_pdf)

cat("Wrote: fig_combine_density_norm, fig_combine_pca_combat (PNG+PDF)\n")
