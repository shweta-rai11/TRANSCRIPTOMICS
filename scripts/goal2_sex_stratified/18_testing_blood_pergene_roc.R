#!/usr/bin/env Rscript
# Individual-gene ROC overlay per sex: each MR-prioritised gene as its own univariate ROC curve, oriented to its better direction, legend sorted by AUC (within-train discrimination only, not a transfer result).
suppressMessages({library(pROC); library(data.table); library(ggplot2)})
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"; fig <- "results/figures/new"
dir.create(tab,showWarnings=FALSE,recursive=TRUE); dir.create(fig,showWarnings=FALSE,recursive=TRUE)

o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr
meta <- as.data.table(o$meta)
ml <- readRDS(file.path(procN, "ml_features.rds"))            # 3-method consensus panels
gsets <- list(Female = ml$female$consensus, Male = ml$male$consensus)
sexcode <- c(Female = "F", Male = "M")

pergene <- function(sexlab) {
  sx <- sexcode[[sexlab]]
  cols <- meta$sample[meta$sex == sx]
  y <- factor(meta$group[match(cols, meta$sample)], levels = c("HC", "RA"))
  genes <- gsets[[sexlab]]; genes <- genes[genes %in% rownames(expr)]
  crd <- list(); gene_auc <- numeric(0)
  for (g in genes) {
    v <- as.numeric(expr[g, cols])
    r <- roc(y, v, direction = "auto", levels = c("HC", "RA"), quiet = TRUE)  # AUC>=0.5
    a <- as.numeric(auc(r)); gene_auc[g] <- a
    crd[[g]] <- data.table(gene = g, sens = r$sensitivities, spec = r$specificities,
                           auc = a, lab = sprintf("%s (AUC=%.3f)", g, a))
  }
  d <- rbindlist(crd)
  ord <- names(sort(gene_auc, decreasing = TRUE))               # legend by AUC desc
  lab_levels <- sprintf("%s (AUC=%.3f)", ord, gene_auc[ord])
  d[, gene := factor(gene, levels = ord)]
  d[, lab := factor(lab, levels = lab_levels)]
  list(coords = d, auc = data.table(sex = sexlab, gene = ord, AUC = round(gene_auc[ord], 3)))
}

plot_overlay <- function(d, ncol_leg) {
  n <- nlevels(d$lab)
  pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  ggplot(d, aes(x = 1 - spec, y = sens, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey75") +
    geom_path(linewidth = 0.7) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    coord_equal() +
    guides(colour = guide_legend(ncol = ncol_leg)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          legend.text = element_text(size = if (n > 20) 6 else 8),
          legend.key.height = unit(if (n > 20) 0.65 else 0.9, "lines"),
          panel.grid.minor = element_blank())
}

auc_tabs <- list()
for (sexlab in c("Female", "Male")) {
  res <- pergene(sexlab); auc_tabs[[sexlab]] <- res$auc
  ncol_leg <- if (nlevels(res$coords$lab) > 20) 2 else 1
  wid <- if (ncol_leg == 2) 9 else 7.5
  g <- plot_overlay(res$coords, ncol_leg)
  tag <- tolower(sexlab)
  ggsave(file.path(fig, sprintf("fig_mr_pergene_roc_%s.png", tag)), g, width = wid, height = 6, dpi = 300)
  ggsave(file.path(fig, sprintf("fig_mr_pergene_roc_%s.pdf", tag)), g, width = wid, height = 6)
  cat(sprintf("wrote fig_mr_pergene_roc_%s.{png,pdf}  (%d genes)\n", tag, nlevels(res$coords$lab)))
}
fwrite(rbindlist(auc_tabs), file.path(tab, "mr_pergene_train_auc.csv"))
cat("wrote mr_pergene_train_auc.csv\n")
