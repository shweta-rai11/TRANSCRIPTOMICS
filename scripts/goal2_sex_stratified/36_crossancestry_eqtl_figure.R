#!/usr/bin/env Rscript
# =============================================================================
# 36_crossancestry_eqtl_figure.R  —  visual companion to 35. Shows whether the
# per-sex EUR MR-causal eQTL genes carry the same causal effect across three RA
# cohorts / two ancestries, using the transferability-aware ancestry_class from
# 35. Two panels per sex, log-log OR grammar (as in 13e):
#   A  Within-ancestry (EUR):  OR Okada 2014  vs  OR Stahl 2010
#   B  Cross-ancestry (EAS):   OR Okada 2014  vs  OR BBJ 2019 (East Asian)
# Colour = ancestry_class (shared EUR+EAS / EUR-replicated / EUR-only /
# untestable in EAS). Genes that transfer to East Asian are labelled.
# Plus a stacked-bar summary of the ancestry_class mix per sex.
#
# Inputs : results/tables/MR35_crossancestry_{female,male}.csv (from 35)
# Outputs: results/figures/fig_mr35_crossancestry_{female,male}.png
#          results/figures/fig_mr35_ancestry_class_summary.png
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(ggrepel); library(patchwork)})
tab <- "results/tables"; fig <- "results/figures"

CLASS <- c(`shared EUR+EAS`                    = "#C62828",
           `EUR-replicated, not EAS`           = "#1565C0",
           `EUR-discovery only`                = "grey65",
           `untestable in EAS (no instruments)`= "grey85")

plot_sex <- function(sx) {
  d <- fread(file.path(tab, sprintf("MR35_crossancestry_%s.csv", sx)))
  d[, ancestry_class := factor(ancestry_class, levels = names(CLASS))]
  rng <- range(c(d$OR_okada, d$OR_stahl, d$OR_bbj, 1), na.rm = TRUE)

  panel <- function(yvar, ytitle, subttl, need) {
    dd <- d[!is.na(OR_okada) & !is.na(get(need))]
    ggplot(dd, aes(OR_okada, .data[[yvar]], colour = ancestry_class)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
      geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey70") +
      geom_point(size = 2.4, alpha = 0.85) +
      ggrepel::geom_text_repel(data = dd[ancestry_class == "shared EUR+EAS"],
                               aes(label = gene), size = 2.9, max.overlaps = Inf,
                               min.segment.length = 0, segment.colour = "grey60",
                               show.legend = FALSE) +
      scale_x_log10() + scale_y_log10() +
      scale_colour_manual(values = CLASS, name = NULL, drop = FALSE) +
      coord_fixed(xlim = rng, ylim = rng) +
      labs(subtitle = subttl, x = "OR (Okada 2014, EUR)", y = ytitle) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            plot.subtitle = element_text(face = "bold", size = 10.5))
  }

  A <- panel("OR_stahl", "OR (Stahl 2010, EUR)",  "A  Within-ancestry (EUR): Okada vs Stahl", "OR_stahl")
  B <- panel("OR_bbj",   "OR (Biobank Japan 2019, East Asian)", "B  Cross-ancestry: European vs East Asian", "OR_bbj")

  nrep <- sum(d$replicated_EUR); ntr <- sum(d$transferable_EAS); nunt <- sum(!d$testable_EAS)
  CAVEAT <- paste(
    "MR design — EXPOSURE (all panels): European eQTL (eQTLGen), the SAME instruments everywhere.",
    "OUTCOME: RA GWAS per cohort — A: Okada & Stahl (European);  B: Biobank Japan (East Asian).",
    "Panel B is EXPLORATORY cross-ancestry: European eQTL exposure vs East Asian RA outcome (ancestry-mismatched);",
    "no East-Asian eQTL exists on OpenGWAS — instrument transferability reported separately.",
    sep = "\n")
  g <- (A + B) + plot_layout(guides = "collect") +
    plot_annotation(
      subtitle = sprintf("%d EUR-causal genes  |  %d replicate in EUR (Stahl)  |  %d transfer to East Asian (BBJ p<0.05, same direction)  |  %d untestable in EAS",
                         nrow(d), nrep, ntr, nunt),
      caption = CAVEAT,
      theme = theme(plot.subtitle = element_text(size = 9.5, colour = "grey35"),
                    plot.caption = element_text(hjust = 0, size = 8.6, colour = "grey20",
                                                face = "italic", lineheight = 1.15,
                                                margin = margin(t = 8)),
                    legend.position = "top")) &
    theme(legend.position = "top")
  ggsave(file.path(fig, sprintf("fig_mr35_crossancestry_%s.png", sx)), g,
         width = 11.5, height = 6.8, dpi = 300)
  cat(sprintf("  wrote fig_mr35_crossancestry_%s.png (%d genes; %d transfer to EAS)\n", sx, nrow(d), ntr))
  d[, sex := sx][, .(sex, gene, ancestry_class)]
}

alld <- rbindlist(list(plot_sex("female"), plot_sex("male")))

# ---- stacked-bar summary of ancestry_class per sex ---------------------------
alld[, ancestry_class := factor(ancestry_class, levels = rev(names(CLASS)))]
cnt <- alld[, .N, by = .(sex, ancestry_class)]
cnt[, sex := tools::toTitleCase(sex)]
gb <- ggplot(cnt, aes(sex, N, fill = ancestry_class)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(N > 0, N, "")), position = position_stack(vjust = 0.5),
            size = 3.1, colour = "white", fontface = "bold") +
  scale_fill_manual(values = CLASS, name = NULL, breaks = names(CLASS)) +
  labs(subtitle = "How many per-sex causal genes replicate in EUR vs transfer to East Asian",
       caption = paste0("Exposure = European eQTL (eQTLGen) for all genes; outcome = RA GWAS ",
                        "(EUR: Okada/Stahl; East Asian: Biobank Japan).\nEast-Asian arm is exploratory ",
                        "(ancestry-mismatched exposure — no East-Asian eQTL available)."),
       x = NULL, y = "Number of genes") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9.5),
        plot.caption = element_text(hjust = 0, size = 8.4, colour = "grey20",
                                    face = "italic", lineheight = 1.15, margin = margin(t = 8)),
        panel.grid.major.x = element_blank(), legend.position = "right")
ggsave(file.path(fig, "fig_mr35_ancestry_class_summary.png"), gb, width = 8.2, height = 4.6, dpi = 300)
cat("  wrote fig_mr35_ancestry_class_summary.png\n")
cat("DONE\n")
