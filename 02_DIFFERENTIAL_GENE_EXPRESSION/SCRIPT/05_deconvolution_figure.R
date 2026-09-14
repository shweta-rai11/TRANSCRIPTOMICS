#!/usr/bin/env Rscript
# Figure: diagnostic panel gene expression + LM22 immune-cell fraction differences (HC vs RA), by sex
suppressMessages({
  library(ggplot2); library(patchwork); library(data.table)
})

fig <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

tr  <- readRDS("data/processed/combined_train.rds")
cf  <- readRDS("data/processed/new/cell_fractions.rds")
ml  <- readRDS("data/processed/new/ml_features.rds")

train_expr <- tr$expr
train_meta <- as.data.table(tr$meta)

panels <- list(Female = ml$female$consensus, Male = ml$male$consensus)

make_expr_panel <- function(sx_code, sx_label) {
  genes <- panels[[sx_label]]
  md <- train_meta[sex == sx_code, .(sample, group)]
  E  <- train_expr[genes, md$sample, drop = FALSE]
  dt <- as.data.table(t(E)); dt[, sample := md$sample]
  dt <- melt(dt, id.vars = "sample", variable.name = "gene", value.name = "expr")
  dt <- merge(dt, md, by = "sample")
  dt[, group := factor(group, levels = c("HC", "RA"))]
  dt[, gene := factor(gene, levels = genes)]

  ggplot(dt, aes(group, expr, fill = group)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(HC = "#4C72B0", RA = "#C44E52")) +
    labs(x = NULL, y = "log2 expression", fill = NULL,
         title = sprintf("A. Diagnostic panel gene expression (%s, consensus panel n=%d)",
                          sx_label, length(genes))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.background = element_rect(fill = "grey90"))
}

make_immune_panel <- function(sx_label) {
  tst <- cf$tests[sex == sx_label][order(p)]
  sig <- tst[FDR < 0.05, cell]
  if (!length(sig)) sig <- head(tst$cell, 6)   # fall back to top-6 if none survive FDR

  fr <- cf$cibersort$train
  sx_code <- ifelse(sx_label == "Female", "F", "M")
  md <- train_meta[sex == sx_code, .(sample, group)]
  fr <- merge(md, fr, by = "sample")

  cols <- paste0(sig, "_CIBERSORT")
  cols <- intersect(cols, names(fr))
  dt <- melt(fr[, c("sample", "group", cols), with = FALSE],
             id.vars = c("sample", "group"), variable.name = "cell", value.name = "fraction")
  dt[, cell := sub("_CIBERSORT$", "", cell)]
  dt[, cell := gsub("_", " ", cell)]
  dt[, group := factor(group, levels = c("HC", "RA"))]
  ord <- gsub("_", " ", sig)
  dt[, cell := factor(cell, levels = ord)]

  ggplot(dt, aes(group, fraction, fill = group)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6) +
    facet_wrap(~ cell, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(HC = "#4C72B0", RA = "#C44E52")) +
    labs(x = NULL, y = "CIBERSORT LM22 fraction", fill = NULL,
         title = sprintf("B. Immune cell composition (%s, FDR<0.05 subsets)", sx_label)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.background = element_rect(fill = "grey90"))
}

for (sx in list(c("F", "Female"), c("M", "Male"))) {
  code <- sx[1]; label <- sx[2]
  A <- make_expr_panel(code, label)
  B <- make_immune_panel(label)
  composite <- A / B
  base <- file.path(fig, sprintf("fig_immune_deconvolution_panel_%s", tolower(label)))
  ggsave(paste0(base, ".png"), composite, width = 11, height = 8, dpi = 220, bg = "white")
  ggsave(paste0(base, ".pdf"), composite, width = 11, height = 8, bg = "white")
  cat("wrote", base, ".{png,pdf}\n")
}
