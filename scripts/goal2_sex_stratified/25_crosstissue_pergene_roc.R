#!/usr/bin/env Rscript
# =============================================================================
# new/27new_pergene_roc_alltissues.R
# -----------------------------------------------------------------------------
# Per-gene ROC curves across all four datasets (Train | Internal test |
# External blood | External synovium), one facet per consensus gene, 4 coloured
# curves, AUC written in each facet. Train orientation kept throughout (a curve
# below the diagonal = expression direction reversed vs blood training).
# Output: results/figures/new/fig_pergene_roc_alltissues.png/pdf
# =============================================================================
suppressMessages({library(data.table); library(pROC); library(ggplot2)})
proc <- "data/processed"; figN <- "results/figures/new"

ro <- readRDS(file.path(proc, "new", "mr_roc_objects.rds"))$roc
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
D  <- readRDS(file.path(proc, "dge_results.rds"))
v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
gF <- ml$female$consensus; gM <- ml$male$consensus
gene_sex <- setNames(c(rep("Female", length(gF)), rep("Male", length(gM))), c(gF, gM))

# ---- Train / Internal / Blood gene ROC coords from saved store --------------
keys <- grep("^gene ", names(ro), value = TRUE)
bti <- rbindlist(lapply(keys, function(k) {
  toks <- strsplit(k, " ")[[1]]; sexlab <- toks[2]; gene <- toks[length(toks)]
  ds <- paste(toks[3:(length(toks)-1)], collapse = " ")
  if (!ds %in% c("Train","Internal test","External blood")) return(NULL)
  d <- as.data.table(ro[[k]]); d[, .(gene, dataset = ds, sens, spec, auc)]
}))

# ---- Synovium per-gene ROC, train-fixed orientation -------------------------
train_dir <- function(g, s) sign(as.data.table(D$res[[s]])$logFC[match(g, as.data.table(D$res[[s]])$gene)])
syn <- rbindlist(lapply(names(gene_sex), function(g) {
  s <- gene_sex[[g]]; if (!g %in% rownames(v$logcpm)) return(NULL)
  idx <- which(v$sex == ifelse(s=="Female","F","M")); y <- factor(v$grp[idx], levels=c("Normal","RA"))
  dir <- if (train_dir(g, s) > 0) "<" else ">"     # up-in-RA-blood -> higher=RA
  r <- roc(y, as.numeric(v$logcpm[g, idx]), direction = dir, levels = c("Normal","RA"), quiet = TRUE)
  data.table(gene = g, dataset = "External synovium", sens = r$sensitivities, spec = r$specificities,
             auc = as.numeric(auc(r)))
}))
d <- rbind(bti, syn)[gene %in% names(gene_sex)]
ORDER <- c("Train","Internal test","External blood","External synovium")
pal <- c("Train"="#1b6ca8","Internal test"="#e08214","External blood"="#4d9221","External synovium"="#7B3294")
d[, dataset := factor(dataset, levels = ORDER)]
d[, gene := factor(gene, levels = c(gF, gM))]

# AUC labels per gene/dataset (dataset name + value, stacked in each facet)
short <- c("Train"="Train", "Internal test"="Internal",
           "External blood"="Blood", "External synovium"="Synovium")
# AUC values go into each panel's STRIP header (outside the plot -> no overlap)
au <- dcast(unique(d[, .(gene, dataset, auc)]), gene ~ dataset, value.var = "auc")
f2 <- function(x) sub("^0", "", sprintf("%.2f", x))       # 0.78 -> .78 (saves width)
strip_lab <- setNames(
  sprintf("%s\nTrain %s   Internal %s\nBlood %s   Synovium %s",
          au$gene, f2(au$Train), f2(au$`Internal test`),
          f2(au$`External blood`), f2(au$`External synovium`)),
  au$gene)

g <- ggplot(d, aes(1 - spec, sens, colour = dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey80") +
  geom_path(linewidth = 1) +
  facet_wrap(~ gene, ncol = 3, labeller = as_labeller(strip_lab)) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       caption = "AUC per dataset in each panel header (train orientation; curve below diagonal = direction reversed vs blood training).") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 10.5, lineheight = 1.1),
        panel.grid.minor = element_blank(), plot.caption = element_text(size = 9, hjust = 0))

ggsave(file.path(figN, "fig_pergene_roc_alltissues.png"), g, width = 12, height = 15, dpi = 300, bg = "white")
ggsave(file.path(figN, "fig_pergene_roc_alltissues.pdf"), g, width = 12, height = 15, bg = "white")
cat("wrote fig_pergene_roc_alltissues\n")
