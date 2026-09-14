#!/usr/bin/env Rscript
# =============================================================================
# 03d_sex_composition_barchart.R  —  two-panel sex composition bar chart for
# the training data: (a) sex distribution of the combined training set,
# (b) RA/Control composition within each sex, ordered by RA proportion.
# Same visual language as 03c_cohort_composition_barchart.R (black = Control,
# grey = RA, dashed line at 50%, n + % labelled above each bar).
#
# Input  : results/tables/combined_cohort_summary.csv
# Outputs: results/figures/fig_sex_composition_bar.png
#          results/figures/fig_sex_composition_bar.pdf
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
tab <- "results/tables"; fig <- "results/figures"

d <- fread(file.path(tab, "combined_cohort_summary.csv"))
d[, group := factor(ifelse(group == "HC", "Control", "RA"), levels = c("Control", "RA"))]
d[, sex := factor(ifelse(sex == "F", "Female", "Male"), levels = c("Female", "Male"))]

# ---- panel a: sex distribution of the whole training set -------------------
bysex_tot <- d[, .(n = sum(Freq)), by = sex]
bysex_tot[, N := sum(n)]; bysex_tot[, pct := 100 * n / N]
Ntot <- bysex_tot$N[1]

# ---- panel b: RA/Control composition within each sex, by RA proportion -----
bysex <- d[, .(n = sum(Freq)), by = .(sex, group)]
bysex[, N := sum(n), by = sex]
bysex[, pct := 100 * n / N]
ord <- bysex[group == "RA"][order(pct), sex]
bysex[, sex := factor(sex, levels = ord)]
bysex[, xlab := sprintf("%s\n(N=%d)", sex, N)]
bysex[, xlab := factor(xlab, levels = bysex[order(sex), unique(xlab)])]

pal_sex <- c(Female = "black", Male = "grey65")
pal_grp <- c(Control = "black", RA = "grey65")
base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "plain", size = 13),
        axis.title = element_text(face = "bold"))

pa <- ggplot(bysex_tot, aes(sex, pct, fill = sex)) +
  geom_col(width = 0.6, colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", n, pct)), vjust = -0.3, size = 3.4, lineheight = 0.9) +
  scale_fill_manual(values = pal_sex, guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Sex", y = "Percentage of samples",
       title = sprintf("a) Combined training set (N = %d)", Ntot)) +
  base_theme

pb <- ggplot(bysex, aes(xlab, pct, fill = group)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8), colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("%d\n%.1f%%", n, pct)), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.2, lineheight = 0.9) +
  scale_fill_manual(values = pal_grp, name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Sex (ordered by RA proportion)", y = "Percentage of samples", title = "b)") +
  base_theme + theme(legend.position = "bottom")

g <- pa + pb + plot_layout(widths = c(1, 1.3), guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(fig, "fig_sex_composition_bar.png"), g, width = 9.5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig, "fig_sex_composition_bar.pdf"), g, width = 9.5, height = 5, bg = "white")
cat("Wrote fig_sex_composition_bar.png + .pdf\nDONE\n")
