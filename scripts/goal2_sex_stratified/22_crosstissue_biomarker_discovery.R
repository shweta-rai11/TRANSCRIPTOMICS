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

gA <- ggplot(cc, aes(blood, syn, colour = sex, shape = concordant)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = gene), size = 3, show.legend = FALSE, max.overlaps = 20) +
  scale_colour_manual(values = c(Female = "#C0392B", Male = "#1F3B99")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), name = "Concordant") +
  labs(x = "Blood log2FC (RA vs HC)", y = "Synovium log2FC (RA vs Normal)",
       title = "A  Blood vs synovium concordance") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

gB <- ggplot(cc, aes(x = reorder(gene, auc), y = auc, fill = sex)) +
  geom_col() + geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_text(aes(label = sprintf("%.2f", auc)), hjust = -0.1, size = 3) +
  coord_flip(ylim = c(0.4, 1.0)) +
  scale_fill_manual(values = c(Female = "#C0392B", Male = "#1F3B99")) +
  labs(x = NULL, y = "Synovium AUC (RA vs Normal)", title = "B  Per-gene synovium discrimination") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

g <- gA + gB + plot_layout(widths = c(1.1, 1))
ggsave(file.path(figN, "fig_crosstissue.png"), g, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(figN, "fig_crosstissue.pdf"), g, width = 12, height = 5.5, bg = "white")
cat("\nwrote fig_crosstissue + crosstissue_panel_auc.csv\n")
