#!/usr/bin/env Rscript
# Training vs testing (hold-out) split bar chart for the combined discovery cohort, showing the stratified 70:30 partition
suppressMessages({library(data.table); library(ggplot2)})
tab <- "results/tables"; fig <- "results/figures"

n_train <- fread(file.path(tab, "combined_cohort_summary.csv"))[, sum(Freq)]
n_hold  <- nrow(fread(file.path(tab, "internal_val_holdout_meta.csv")))
N <- n_train + n_hold

d <- data.table(split = factor(c("Training", "Testing (hold-out)"), levels = c("Training", "Testing (hold-out)")),
                 n = c(n_train, n_hold))
d[, pct := 100 * n / N]

pal <- c(Training = "black", "Testing (hold-out)" = "grey65")

p <- ggplot(d, aes(split, pct, fill = split)) +
  geom_col(width = 0.55, colour = "black") +
  geom_hline(yintercept = 70, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", n, pct)), vjust = -0.3, size = 3.6, lineheight = 0.9) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage of samples") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold"))

ggsave(file.path(fig, "fig_train_test_split_bar.png"), p, width = 5.5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig, "fig_train_test_split_bar.pdf"), p, width = 5.5, height = 5, bg = "white")
cat("Wrote fig_train_test_split_bar.png + .pdf\nDONE\n")
