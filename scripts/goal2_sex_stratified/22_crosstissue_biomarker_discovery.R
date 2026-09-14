#!/usr/bin/env Rscript
# =============================================================================
# new/23new_crosstissue.R  -- cross-tissue (blood -> synovium) summary.
#  (1) PANEL AUC in synovium (RA vs Normal), per sex: z-scored panel genes,
#      logistic model, apparent + 10-fold CV AUC (+ DeLong CI).
#  (2) Figure: (A) blood vs synovium log2FC concordance scatter per gene;
#              (B) per-gene synovium AUC bar (coloured by sex).
# Uses data/processed/new/val_synovium.rds (from 21) + blood DE + new panels.
# Output: results/tables/crosstissue_panel_auc.csv
#         results/figures/new/fig_crosstissue.png/pdf
# =============================================================================
suppressMessages({library(data.table); library(pROC); library(ggplot2); library(patchwork)})
proc <- "data/processed"; tabN <- "results/tables"; figN <- "results/figures/new"

v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
logcpm <- v$logcpm; grp <- v$grp; sex <- v$sex
panels <- list(Female = ml$female$consensus, Male = ml$male$consensus)
sc <- c(Female = "F", Male = "M")

# ---- (1) panel AUC in synovium: z-score genes, logistic, apparent + 10-fold CV
GLOBAL_SEED <- 1234
panel_auc <- function(sexlab) {
  genes <- panels[[sexlab]]; genes <- genes[genes %in% rownames(logcpm)]
  idx <- which(sex == sc[[sexlab]])
  y <- factor(grp[idx], levels = c("Normal", "RA"))
  X <- t(logcpm[genes, idx, drop = FALSE])
  Z <- scale(X); Z[is.na(Z)] <- 0; colnames(Z) <- make.names(genes)
  df <- data.frame(y = y, Z, check.names = FALSE)
  fit <- suppressWarnings(glm(y ~ ., df, family = binomial))
  pa <- as.numeric(predict(fit, type = "response"))
  ra <- roc(y, pa, levels = c("Normal","RA"), direction = "<", quiet = TRUE)
  # 10-fold CV
  set.seed(GLOBAL_SEED); fold <- sample(rep_len(1:10, length(y)))
  pcv <- rep(NA_real_, length(y))
  for (k in 1:10) { tr <- fold != k; te <- fold == k
    if (length(unique(y[tr])) < 2) next
    f <- suppressWarnings(glm(y ~ ., df[tr, , drop = FALSE], family = binomial))
    pcv[te] <- as.numeric(predict(f, df[te, , drop = FALSE], type = "response")) }
  rcv <- roc(y[!is.na(pcv)], pcv[!is.na(pcv)], levels = c("Normal","RA"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(rcv))
  data.table(sex = sexlab, n = length(y), n_RA = sum(y=="RA"), n_Normal = sum(y=="Normal"),
             genes = length(genes), apparent_AUC = round(as.numeric(auc(ra)), 3),
             CV_AUC = round(as.numeric(auc(rcv)), 3), CV_lo = round(ci[1],3), CV_hi = round(ci[3],3))
}
pa <- rbindlist(lapply(c("Female","Male"), panel_auc))
fwrite(pa, file.path(tabN, "crosstissue_panel_auc.csv"))
cat("=== PANEL AUC in synovium (RA vs Normal) ===\n"); print(pa)

# ---- (2A) blood vs synovium log2FC concordance ------------------------------
D <- readRDS(file.path(proc, "dge_results.rds"))
bloodFC <- function(g, s) { d <- as.data.table(D$res[[s]]); d$logFC[match(g, d$gene)] }
mk <- function(sexlab, tab) { g <- panels[[sexlab]]
  data.table(gene = g, sex = sexlab,
             blood = bloodFC(g, sexlab),
             syn = tab$syn_log2FC[match(g, tab$gene)],
             auc = tab$auc_all[match(g, tab$gene)]) }
cc <- rbindlist(list(mk("Female", v$sf), mk("Male", v$sm)))
cc[, concordant := sign(blood) == sign(syn)]

sex_cols <- c(Female = "#C0392B", Male = "#1F3B99")
base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.title = element_text(size = 11),
    axis.text = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    panel.border = element_rect(colour = "grey40", linewidth = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.9, "lines")
  )

gA <- ggplot(cc, aes(blood, syn, colour = sex, shape = concordant)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_point(size = 2.8, stroke = 1) +
  ggrepel::geom_text_repel(aes(label = gene), size = 3.2, colour = "grey20",
                            fontface = "italic", show.legend = FALSE,
                            max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3,
                            min.segment.length = 0, segment.colour = "grey50",
                            segment.size = 0.3, seed = GLOBAL_SEED) +
  scale_colour_manual(values = sex_cols, name = "Sex") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), name = "Concordant\ndirection") +
  labs(x = "Blood log2FC (RA vs HC)", y = "Synovium log2FC (RA vs Normal)",
       title = "A  Blood vs synovium concordance") +
  base_theme

gB <- ggplot(cc, aes(x = reorder(gene, auc), y = auc, colour = sex)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(xend = gene, y = 0.5, yend = auc), linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", auc)), hjust = -0.45, size = 3.2, colour = "grey20") +
  coord_flip(ylim = c(0.4, 1.0), clip = "off") +
  scale_y_continuous(breaks = seq(0.4, 1.0, 0.2), expand = expansion(mult = c(0.02, 0.12))) +
  scale_colour_manual(values = sex_cols, name = "Sex", guide = "none") +
  labs(x = NULL, y = "Synovium AUC (RA vs Normal)", title = "B  Per-gene synovium discrimination") +
  base_theme +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "italic"))

g <- gA + gB + plot_layout(widths = c(1.1, 1), guides = "collect") &
  theme(legend.position = "right")
ggsave(file.path(figN, "fig_crosstissue.png"), g, width = 12.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(figN, "fig_crosstissue.pdf"), g, width = 12.5, height = 5.5, bg = "white")
cat("\nwrote fig_crosstissue + crosstissue_panel_auc.csv\n")
