#!/usr/bin/env Rscript
# Fig 4C: clustered heatmap of per-sample z-scored expression of MR consensus genes, annotated by Group and Sex.
suppressMessages({library(ComplexHeatmap); library(circlize); library(data.table)})
proc <- "data/processed"; fig <- "results/figures/new"

o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr
meta <- as.data.table(o$meta)
ml <- readRDS("data/processed/new/ml_features.rds")
genesF <- ml$female$consensus; genesM <- ml$male$consensus
genes  <- c(genesF, genesM); genes <- genes[genes %in% rownames(expr)]

M  <- t(scale(t(expr[genes, ])))                       # z-score each gene (row)
M[is.na(M)] <- 0
lim <- 2.5; M[M >  lim] <-  lim; M[M < -lim] <- -lim    # cap for display

# Row split labels: sex-stratified panels, not a sex-specificity claim
row_split <- factor(ifelse(genes %in% genesF, "Female panel", "Male panel"),
                    levels = c("Female panel", "Male panel"))
col_ha <- HeatmapAnnotation(
  Group = factor(meta$group, levels = c("HC", "RA")),
  Sex   = factor(ifelse(meta$sex == "F", "Female", "Male")),
  col = list(Group = c(HC = "#4DAF4A", RA = "#984EA3"),
             Sex   = c(Female = "#C0392B", Male = "#2E86C1")),
  annotation_name_side = "left", simple_anno_size = unit(4, "mm"))

ht <- Heatmap(M, name = "z-score",
  col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
  top_annotation = col_ha, row_split = row_split,
  row_title_rot = 90, row_names_side = "left",
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_column_names = FALSE, show_column_dend = TRUE,
  column_title = "Samples (clustered)", row_gap = unit(2, "mm"),
  row_title_gp = gpar(fontface = "bold"),
  heatmap_legend_param = list(title = "Expression\n(z-score)"))

png(file.path(fig, "fig_mr_panelC_clustered_heatmap.png"), width = 2400, height = 1200, res = 300)
draw(ht, merge_legend = TRUE); dev.off()
pdf(file.path(fig, "fig_mr_panelC_clustered_heatmap.pdf"), width = 8, height = 4)
draw(ht, merge_legend = TRUE); dev.off()
cat("wrote fig_mr_panelC_clustered_heatmap.{png,pdf}\n")
