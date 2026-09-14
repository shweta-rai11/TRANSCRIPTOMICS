#!/usr/bin/env Rscript
# Per-gene AUC of every consensus panel gene across all four datasets (train, internal test, external blood, external synovium), all on train-fixed orientation so reversal reads as AUC<0.5.
suppressMessages({library(data.table); library(ggplot2)})
proc <- "data/processed"; tabN <- "results/tables"; figN <- "results/figures/new"

pg <- fread(file.path(tabN, "mr_roc_pergene_auc.csv"))          # Train/Internal/Blood
v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
gF <- ml$female$consensus; gM <- ml$male$consensus
# Panels share genes (ESYT1, MED1, SMARCC2), so the key is (gene, panel), not gene alone
panel_map <- rbind(data.table(gene = gF, panel = "Female"),
                   data.table(gene = gM, panel = "Male"))

# blood/train/internal - keep sex; a shared gene has a DIFFERENT AUC in each panel
bti <- pg[dataset %in% c("Train","Internal test","External blood"),
          .(gene, panel = sex, dataset, AUC)]
# Synovium (train-fixed orientation, per-sex): uses auc_sex not auc_all since this is sex-stratified
syn <- rbind(cbind(as.data.table(v$sf), panel = "Female"),
             cbind(as.data.table(v$sm), panel = "Male"), fill = TRUE)[
         , .(gene, panel, dataset = "External synovium",
             AUC = ifelse(concordant, auc_sex, 1 - auc_sex))]
# inner join keeps only (gene, panel) pairs the gene actually belongs to
d <- merge(rbind(bti, syn), panel_map, by = c("gene", "panel"))
d[, dataset := factor(dataset, levels = c("Train","Internal test","External blood","External synovium"))]
d[, gene := factor(gene, levels = unique(c(gF, gM)))]
fwrite(dcast(d, gene + panel ~ dataset, value.var = "AUC"),
       file.path(tabN, "pergene_auc_alltissues.csv"))

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
print(dcast(d, gene + panel ~ dataset, value.var = "AUC"))
