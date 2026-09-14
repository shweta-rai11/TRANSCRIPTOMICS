#!/usr/bin/env Rscript
# =============================================================================
# 42_figure_crosstissue_summary.R -- summary figure for the blood/synovium
# consensus-panel overlap and concordance (companion to fig_crosstissue.png).
#  (A) Venn of female vs male consensus biomarker panels.
#  (B) Slopegraph of blood -> synovium log2FC per gene x sex, styled by
#      concordant/discordant direction.
# Uses data/processed/new/val_synovium.rds + ml_features.rds + dge_results.rds
# Output: results/figures/new/fig_crosstissue_summary.png/pdf
# =============================================================================
suppressMessages({
  library(data.table); library(ggplot2); library(ggforce); library(ggrepel)
  library(patchwork); library(eulerr)
})
proc <- "data/processed"; figN <- "results/figures/new"

v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
D  <- readRDS(file.path(proc, "dge_results.rds"))
panels <- list(Female = ml$female$consensus, Male = ml$male$consensus)

bloodFC <- function(g, s) { d <- as.data.table(D$res[[s]]); d$logFC[match(g, d$gene)] }
mk <- function(sexlab, tab) { g <- panels[[sexlab]]
  data.table(gene = g, sex = sexlab, blood = bloodFC(g, sexlab),
             syn = tab$syn_log2FC[match(g, tab$gene)]) }
cc <- rbindlist(list(mk("Female", v$sf), mk("Male", v$sm)))
cc[, concordant := factor(sign(blood) == sign(syn), levels = c(TRUE, FALSE),
                           labels = c("Concordant", "Discordant"))]

sex_cols <- c(Female = "#C0392B", Male = "#1F3B99")
base_theme <- theme_void(base_size = 12) +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9))

# ---- (A) Venn of the two consensus panels -----------------------------------
shared   <- intersect(panels$Female, panels$Male)
fem_only <- setdiff(panels$Female, shared)
male_only<- setdiff(panels$Male, shared)

# geometry solved by eulerr (exact, area-accurate circle placement) rather than
# hand-picked coordinates -- avoids labels landing on/inside the wrong circle.
fit <- euler(c(Female = length(fem_only), Male = length(male_only),
               "Female&Male" = length(shared)), shape = "circle")
ell   <- fit$ellipses
r_c   <- unname(ell["Female", "a"])
h_c   <- unname(abs(ell["Female", "h"]))
circ  <- data.table(x = c(-h_c, h_c), y = c(0, 0), sex = c("Female", "Male"))

label_offset <- 0.65 * r_c          # clear of the lens boundary (h_c - r_c) with margin
fem_label_x  <- -h_c - label_offset
male_label_x <-  h_c + label_offset

gA <- ggplot() +
  geom_circle(data = circ, aes(x0 = x, y0 = y, r = r_c, colour = sex),
              fill = NA, linewidth = 1) +
  annotate("text", x = -h_c - 0.35, y = r_c + 0.22, label = "Female",
           colour = sex_cols["Female"], fontface = "bold", size = 4, hjust = 0.5) +
  annotate("text", x = h_c + 0.35, y = r_c + 0.22, label = "Male",
           colour = sex_cols["Male"], fontface = "bold", size = 4, hjust = 0.5) +
  annotate("text", x = fem_label_x, y = 0, label = paste(fem_only, collapse = "\n"),
           fontface = "italic", size = 3.6, colour = sex_cols["Female"]) +
  annotate("text", x = 0, y = 0, label = paste(shared, collapse = "\n"),
           fontface = "bold.italic", size = 3.8, colour = "grey15") +
  annotate("text", x = male_label_x, y = 0, label = paste(male_only, collapse = "\n"),
           fontface = "italic", size = 3.6, colour = sex_cols["Male"]) +
  scale_colour_manual(values = sex_cols, guide = "none") +
  coord_fixed(clip = "off") +
  base_theme

# ---- (B) slopegraph: blood log2FC -> synovium log2FC ------------------------
long <- rbindlist(list(
  cc[, .(gene, sex, concordant, x = 1, y = blood)],
  cc[, .(gene, sex, concordant, x = 2, y = syn)]
))
long[, highlight := gene == "SMARCC2"]

gB <- ggplot(long, aes(x, y, group = interaction(gene, sex))) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_line(aes(colour = sex, linetype = concordant, linewidth = highlight)) +
  geom_point(aes(colour = sex), size = 1.8) +
  geom_text_repel(data = long[x == 2], aes(label = gene), fontface = "italic",
                   colour = "grey20", size = 3.2, direction = "y", hjust = 0,
                   nudge_x = 0.15, xlim = c(2.1, NA), segment.size = 0.3,
                   segment.colour = "grey60", box.padding = 0.15, seed = 1) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Blood", "Synovium"),
                      limits = c(0.8, 2.75), expand = c(0, 0)) +
  scale_colour_manual(values = sex_cols, name = "Sex") +
  scale_linetype_manual(values = c(Concordant = "solid", Discordant = "22"), name = "Direction") +
  scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.5), guide = "none") +
  labs(x = NULL, y = expression(log[2]*FC~(RA~vs~control))) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "grey40"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 10), legend.text = element_text(size = 9)
  )

g <- gA + gB + plot_layout(widths = c(0.85, 1))
ggsave(file.path(figN, "fig_crosstissue_summary.png"), g, width = 10.5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(figN, "fig_crosstissue_summary.pdf"), g, width = 10.5, height = 5, bg = "white")
cat("wrote fig_crosstissue_summary.png/pdf\n")
