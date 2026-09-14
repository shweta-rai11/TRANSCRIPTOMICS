#!/usr/bin/env Rscript
# Per-gene ROC curves across all four datasets, one facet per consensus gene with AUC per facet, train orientation kept throughout; female and male plotted as separate figures.
suppressMessages({library(data.table); library(pROC); library(ggplot2)})
proc <- "data/processed"; figN <- "results/figures/new"

ro <- readRDS(file.path(proc, "new", "mr_roc_objects.rds"))$roc
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
D  <- readRDS(file.path(proc, "dge_results.rds"))
v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
gF <- ml$female$consensus; gM <- ml$male$consensus
# Panels share genes (ESYT1, MED1, SMARCC2), so the evaluation key is (gene, sex), not gene alone
panel_map <- rbind(data.table(gene = gF, sex = "Female"),
                   data.table(gene = gM, sex = "Male"))

# ---- Train / Internal / Blood gene ROC coords from saved store --------------
keys <- grep("^gene ", names(ro), value = TRUE)
bti <- rbindlist(lapply(keys, function(k) {
  toks <- strsplit(k, " ")[[1]]; sexlab <- toks[2]; gene <- toks[length(toks)]
  ds <- paste(toks[3:(length(toks)-1)], collapse = " ")
  if (!ds %in% c("Train","Internal test","External blood")) return(NULL)
  d <- as.data.table(ro[[k]]); d[, .(gene, sex = sexlab, dataset = ds, sens, spec, auc)]
}))

# ---- Synovium per-gene ROC, train-fixed orientation -------------------------
train_dir <- function(g, s) sign(as.data.table(D$res[[s]])$logFC[match(g, as.data.table(D$res[[s]])$gene)])
syn <- rbindlist(lapply(seq_len(nrow(panel_map)), function(i) {
  g <- panel_map$gene[i]; s <- panel_map$sex[i]
  if (!g %in% rownames(v$logcpm)) return(NULL)
  idx <- which(v$sex == ifelse(s=="Female","F","M")); y <- factor(v$grp[idx], levels=c("Normal","RA"))
  dir <- if (train_dir(g, s) > 0) "<" else ">"     # up-in-RA-blood -> higher=RA
  r <- roc(y, as.numeric(v$logcpm[g, idx]), direction = dir, levels = c("Normal","RA"), quiet = TRUE)
  data.table(gene = g, sex = s, dataset = "External synovium", sens = r$sensitivities,
             spec = r$specificities, auc = as.numeric(auc(r)))
}))
# inner join keeps only the (gene, sex) pairs the gene actually belongs to
d <- merge(rbind(bti, syn), panel_map, by = c("gene", "sex"))
ORDER <- c("Train","Internal test","External blood","External synovium")
pal <- c("Train"="#1b6ca8","Internal test"="#e08214","External blood"="#4d9221","External synovium"="#7B3294")
d[, dataset := factor(dataset, levels = ORDER)]
# one facet per (gene, sex): shared genes get a facet in EACH panel
FACETS <- sprintf("%s (%s)", panel_map$gene, panel_map$sex)
d[, gene := factor(sprintf("%s (%s)", gene, sex), levels = FACETS)]

# AUC labels per gene/dataset (dataset name + value, stacked in each facet)
short <- c("Train"="Train", "Internal test"="Internal",
           "External blood"="Blood", "External synovium"="Synovium")
# AUC values go into each panel's STRIP header (outside the plot -> no overlap)
au <- dcast(unique(d[, .(gene, dataset, auc)]), gene ~ dataset, value.var = "auc")
f2 <- function(x) sprintf("%.2f", x)                       # keep leading 0, e.g. 0.78
strip_lab <- setNames(
  sprintf("%s\nTrain %s   Internal %s\nBlood %s   Synovium %s",
          au$gene, f2(au$Train), f2(au$`Internal test`),
          f2(au$`External blood`), f2(au$`External synovium`)),
  au$gene)

plot_sex <- function(sexlab) {
  ds <- d[sex == sexlab]; ds[, gene := droplevels(gene)]
  ggplot(ds, aes(1 - spec, sens, colour = dataset)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey80") +
    geom_path(linewidth = 1) +
    facet_wrap(~ gene, ncol = 3, labeller = as_labeller(strip_lab)) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
    coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", legend.text = element_text(size = 12),
          strip.text = element_text(face = "bold", size = 10.5, lineheight = 1.1),
          panel.grid.minor = element_blank())
}

for (sexlab in c("Female", "Male")) {
  g <- plot_sex(sexlab); tag <- tolower(sexlab)
  ggsave(file.path(figN, sprintf("fig_pergene_roc_alltissues_%s.png", tag)), g,
         width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(figN, sprintf("fig_pergene_roc_alltissues_%s.pdf", tag)), g,
         width = 10, height = 7, bg = "white")
  cat(sprintf("wrote fig_pergene_roc_alltissues_%s\n", tag))
}
unlink(file.path(figN, c("fig_pergene_roc_alltissues.png", "fig_pergene_roc_alltissues.pdf")))
