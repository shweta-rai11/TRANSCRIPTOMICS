#!/usr/bin/env Rscript
# =============================================================================
# 03_dataset_gene_overlap_venn.R  —  Venn of gene-symbol overlap between the two
# TRAINING datasets (GSE93272 [GPL570] and GSE110169 [GPL13667]). Shows how many
# genes are measured in each platform and how many are COMMON (the intersection
# used for the merged analysis).
#
# Inputs : data/raw/GSE93272_raw.rds, data/raw/GSE110169_raw.rds
# Outputs: results/figures/fig_dataset_gene_overlap_venn.png
#          results/tables/dataset_gene_overlap.csv
# ---- Reference: Langfelder & Horvath (WGCNA collapseRows). BMC Bioinf 2008;9:559.
# =============================================================================
suppressMessages({library(Biobase); library(data.table); library(ggVennDiagram); library(ggplot2)})
raw <- "data/raw"; fig <- "results/figures"; tab <- "results/tables"

# gene symbols measured on a dataset = valid, non-ambiguous fData symbols
genes_of <- function(rds) {
  e <- readRDS(file.path(raw, rds)); if (is(e, "list")) e <- e[[1]]
  fd <- fData(e)
  col <- grep("^gene[ ._]?symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  sym <- as.character(fd[[col]])
  sym <- sym[!is.na(sym) & sym != "" & !grepl("///", sym)]   # drop empty / ambiguous
  unique(sym)
}
g93  <- genes_of("GSE93272_raw.rds")
g110 <- genes_of("GSE110169_raw.rds")
common <- intersect(g93, g110)
cat(sprintf("GSE93272 genes: %d | GSE110169 genes: %d | COMMON: %d\n",
            length(g93), length(g110), length(common)))

# short set labels (the region counts already give per-dataset totals) so the
# labels don't run into the circles; extra left/right room via x expansion.
sets <- setNames(list(g93, g110), c("GSE93272", "GSE110169"))
p <- ggVennDiagram(sets, label = "count", label_alpha = 0, edge_size = 0.6, set_size = 4.5) +
  scale_fill_gradient(low = "#EAF3FB", high = "#4981BF") +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_cartesian(clip = "off") +
  labs(title = sprintf("Common genes between training datasets = %d", length(common)),
       subtitle = sprintf("GSE93272 = %d genes | GSE110169 = %d genes", length(g93), length(g110))) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        plot.margin = margin(10, 24, 10, 24))
ggsave(file.path(fig, "fig_dataset_gene_overlap_venn.png"), p, width = 7.2, height = 6.0, dpi = 300)

fwrite(data.table(dataset = c("GSE93272","GSE110169","common"),
                  n_genes = c(length(g93), length(g110), length(common))),
       file.path(tab, "dataset_gene_overlap.csv"))
cat("Wrote fig_dataset_gene_overlap_venn.png + dataset_gene_overlap.csv\nDONE\n")
