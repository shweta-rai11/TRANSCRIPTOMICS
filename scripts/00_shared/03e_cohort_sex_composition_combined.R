#!/usr/bin/env Rscript
# =============================================================================
# 03e_cohort_sex_composition_combined.R  —  single 4-panel figure merging the
# dataset-composition bar chart (03c) and the sex-composition bar chart (03d):
#   a) Combined training set, by group (Control/RA)
#   b) Per-dataset composition, ordered by RA proportion
#   c) Combined training set, by sex (Female/Male)
#   d) RA/Control composition within each sex, ordered by RA proportion
# One shared Control/RA legend at the bottom (panels a/c need no legend --
# their bars are already labelled on the x-axis).
#
# Input  : results/tables/combined_cohort_summary.csv
# Outputs: results/figures/fig_cohort_sex_composition_bar.png
#          results/figures/fig_cohort_sex_composition_bar.pdf
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
tab <- "results/tables"; fig <- "results/figures"

d <- fread(file.path(tab, "combined_cohort_summary.csv"))
d[, group := factor(ifelse(group == "HC", "Control", "RA"), levels = c("Control", "RA"))]
d[, sex   := factor(ifelse(sex == "F", "Female", "Male"), levels = c("Female", "Male"))]

pal_grp <- c(Control = "black", RA = "grey65")
pal_sex <- c(Female = "black", Male = "grey65")
base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "plain", size = 13),
        axis.title = element_text(face = "bold"))

## ---- a) combined training set, by group ------------------------------------
comb <- d[, .(n = sum(Freq)), by = group]; comb[, N := sum(n)]; comb[, pct := 100 * n / N]
Ntot <- comb$N[1]
pa <- ggplot(comb, aes(group, pct, fill = group)) +
  geom_col(width = 0.6, colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", n, pct)), vjust = -0.3, size = 3.2, lineheight = 0.9) +
  scale_fill_manual(values = pal_grp, guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Group", y = "Percentage of samples",
       title = sprintf("a) Combined training set (N = %d)", Ntot)) +
  base_theme

## ---- b) per-dataset composition, ordered by RA proportion ------------------
byds <- d[, .(n = sum(Freq)), by = .(dataset, group)]
byds[, N := sum(n), by = dataset]; byds[, pct := 100 * n / N]
ordd <- byds[group == "RA"][order(pct), dataset]
byds[, dataset := factor(dataset, levels = ordd)]
byds[, xlab := sprintf("%s\n(N=%d)", dataset, N)]
byds[, xlab := factor(xlab, levels = byds[order(dataset), unique(xlab)])]
pb <- ggplot(byds, aes(xlab, pct, fill = group)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8), colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("%d\n%.1f%%", n, pct)), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = pal_grp, name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Dataset (ordered by RA proportion)", y = "Percentage of samples",
       title = "b) Composition by dataset") +
  base_theme

## ---- c) combined training set, by sex --------------------------------------
bysex_tot <- d[, .(n = sum(Freq)), by = sex]; bysex_tot[, N := sum(n)]; bysex_tot[, pct := 100 * n / N]
pc <- ggplot(bysex_tot, aes(sex, pct, fill = sex)) +
  geom_col(width = 0.6, colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", n, pct)), vjust = -0.3, size = 3.2, lineheight = 0.9) +
  scale_fill_manual(values = pal_sex, guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Sex", y = "Percentage of samples",
       title = sprintf("c) Combined training set (N = %d)", Ntot)) +
  base_theme

## ---- d) RA/Control composition within each sex, ordered by RA proportion ---
bysex <- d[, .(n = sum(Freq)), by = .(sex, group)]
bysex[, N := sum(n), by = sex]; bysex[, pct := 100 * n / N]
ords <- bysex[group == "RA"][order(pct), sex]
bysex[, sex := factor(sex, levels = ords)]
bysex[, xlab := sprintf("%s\n(N=%d)", sex, N)]
bysex[, xlab := factor(xlab, levels = bysex[order(sex), unique(xlab)])]
pd <- ggplot(bysex, aes(xlab, pct, fill = group)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8), colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("%d\n%.1f%%", n, pct)), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = pal_grp, name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Sex (ordered by RA proportion)", y = "Percentage of samples",
       title = "d) Composition by sex") +
  base_theme

g <- (pa + pb) / (pc + pd) +
  plot_layout(widths = c(1, 1.3), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig, "fig_cohort_sex_composition_bar.png"), g, width = 10, height = 9.5, dpi = 300, bg = "white")
ggsave(file.path(fig, "fig_cohort_sex_composition_bar.pdf"), g, width = 10, height = 9.5, bg = "white")
cat("Wrote fig_cohort_sex_composition_bar.png + .pdf\nDONE\n")
