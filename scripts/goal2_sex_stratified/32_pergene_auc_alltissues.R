#!/usr/bin/env Rscript
# =============================================================================
# new/26new_pergene_auc_alltissues.R
# -----------------------------------------------------------------------------
# Per-gene AUC of every consensus panel gene across all four datasets:
#   Train | Internal test | External blood | External synovium.
# Train/Internal/Blood come from mr_roc_pergene_auc.csv (18b, train-fixed
# orientation). Synovium AUC put on the SAME train orientation (concordant ->
# auc_all; discordant -> 1-auc_all) so a direction reversal reads as AUC<0.5.
# Output: results/figures/new/fig_pergene_auc_alltissues.png/pdf
#         results/tables/pergene_auc_alltissues.csv
# =============================================================================
suppressMessages({library(data.table); library(ggplot2)})
proc <- "data/processed"; tabN <- "results/tables"; figN <- "results/figures/new"

pg <- fread(file.path(tabN, "mr_roc_pergene_auc.csv"))          # Train/Internal/Blood
v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
gF <- ml$female$consensus; gM <- ml$male$consensus
gene_panel <- setNames(c(rep("Female", length(gF)), rep("Male", length(gM))), c(gF, gM))

# blood/train/internal
bti <- pg[dataset %in% c("Train","Internal test","External blood"),
          .(gene, dataset, AUC)]
# synovium (train-fixed orientation)
syn <- rbind(as.data.table(v$sf), as.data.table(v$sm), fill = TRUE)[
         , .(gene, dataset = "External synovium",
             AUC = ifelse(concordant, auc_all, 1 - auc_all))]
d <- rbind(bti, syn)
d <- d[gene %in% names(gene_panel)]
d[, panel := gene_panel[gene]]
d[, dataset := factor(dataset, levels = c("Train","Internal test","External blood","External synovium"))]
d[, gene := factor(gene, levels = c(gF, gM))]
fwrite(dcast(d, gene ~ dataset, value.var = "AUC"), file.path(tabN, "pergene_auc_alltissues.csv"))

pal <- c("Train"="#1b6ca8", "Internal test"="#e08214",
         "External blood"="#4d9221", "External synovium"="#7B3294")
g <- ggplot(d, aes(gene, AUC, fill = dataset)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45") +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUC (RA vs control)",
       caption = "Per-gene AUC across datasets (train orientation; AUC < 0.50 = expression direction reversed vs blood training). Dashed line = chance (0.5).") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "top", strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.caption = element_text(size = 8, hjust = 0))

ggsave(file.path(figN, "fig_pergene_auc_alltissues.png"), g, width = 11, height = 6, dpi = 300, bg = "white")
ggsave(file.path(figN, "fig_pergene_auc_alltissues.pdf"), g, width = 11, height = 6, bg = "white")
cat("wrote fig_pergene_auc_alltissues + pergene_auc_alltissues.csv\n")
print(dcast(d, gene ~ dataset, value.var = "AUC"))
