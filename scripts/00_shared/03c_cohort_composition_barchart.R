#!/usr/bin/env Rscript
# Two-panel RA/Control composition bar chart: (a) combined training set, (b) each platform, ordered by RA proportion
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
tab <- "results/tables"; fig <- "results/figures"

d <- fread(file.path(tab, "combined_cohort_summary.csv"))
d[, group := factor(ifelse(group == "HC", "Control", "RA"), levels = c("Control", "RA"))]

# ---- per-dataset composition, ordered by RA proportion ---------------------
byds <- d[, .(n = sum(Freq)), by = .(dataset, group)]
byds[, N := sum(n), by = dataset]
byds[, pct := 100 * n / N]
ord <- byds[group == "RA"][order(pct), dataset]
byds[, dataset := factor(dataset, levels = ord)]
byds[, xlab := sprintf("%s\n(N=%d)", dataset, N)]
byds[, xlab := factor(xlab, levels = byds[order(dataset), unique(xlab)])]

# ---- pooled/combined composition --------------------------------------------
comb <- d[, .(n = sum(Freq)), by = group]
comb[, N := sum(n)]; comb[, pct := 100 * n / N]
Ntot <- comb$N[1]

pal <- c(Control = "black", RA = "grey65")
base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "plain", size = 13),
        axis.title = element_text(face = "bold"))

pa <- ggplot(comb, aes(group, pct, fill = group)) +
  geom_col(width = 0.6, colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", n, pct)), vjust = -0.3, size = 3.4, lineheight = 0.9) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Group", y = "Percentage of samples",
       title = sprintf("a) Combined training set (N = %d)", Ntot)) +
  base_theme

pb <- ggplot(byds, aes(xlab, pct, fill = group)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8), colour = "black") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("%d\n%.1f%%", n, pct)), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.2, lineheight = 0.9) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = "Dataset (ordered by RA proportion)", y = "Percentage of samples", title = "b)") +
  base_theme + theme(legend.position = "bottom")

g <- pa + pb + plot_layout(widths = c(1, 1.3), guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(fig, "fig_cohort_composition_bar.png"), g, width = 9.5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig, "fig_cohort_composition_bar.pdf"), g, width = 9.5, height = 5, bg = "white")
cat("Wrote fig_cohort_composition_bar.png + .pdf\nDONE\n")
