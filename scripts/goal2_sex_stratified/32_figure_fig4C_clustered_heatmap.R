#!/usr/bin/env Rscript
# =============================================================================
# 18m_panelC_clustered_heatmap.R  -- reference Fig 4C analog (clustered heatmap).
# Per-sample z-scored expression of the 6 MR consensus genes across the combined
# cohort, columns clustered (like the ADEx heatmap in the paper), annotated by
# Group (HC/RA) and Sex. Rows split Female-specific (red) / Male-specific (blue).
# This is the honest RA analog of their ADEx panel: RA has no ADEx database, so
# we show the study cohort's own per-sample expression, clustered.
# Output: results/figures/fig_mr_panelC_clustered_heatmap.png/pdf
# =============================================================================
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

# Labels corrected 2026-07-27: were "Female-specific" / "Male-specific" (see
# 20_fig4BC_logfc_heatmaps.R for the full rationale). This chapter makes no
# sex-specificity claim; these are sex-STRATIFIED panels.
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
