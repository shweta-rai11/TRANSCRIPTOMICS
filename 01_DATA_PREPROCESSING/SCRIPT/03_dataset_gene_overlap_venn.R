#!/usr/bin/env Rscript
# Venn/Euler diagram of gene-symbol overlap between the two training datasets (GSE93272, GSE110169)
suppressMessages({library(Biobase); library(data.table); library(eulerr); library(scales); library(grid)})
raw <- "data/raw"; fig <- "results/figures"; tab <- "results/tables"
for (d in c(fig, tab)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

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

# area-proportional Euler diagram, plain black-on-white, dataset names drawn outside the circles
only93 <- length(setdiff(g93, g110)); only110 <- length(setdiff(g110, g93))
tot93 <- only93 + length(common); tot110 <- only110 + length(common)
fit <- euler(c("GSE93272" = only93, "GSE110169" = only110,
               "GSE93272&GSE110169" = length(common)))
# label each circle with its platform's TOTAL gene count, not the exclusive-only count
p <- plot(fit,
  quantities = list(labels = c(comma(tot93), comma(tot110), comma(length(common))),
                     fontsize = 11, col = "black"),
  fills = list(fill = "white", alpha = 1),
  edges = list(col = "black", lwd = 1.6),
  labels = FALSE,
  legend = FALSE,
  main = NULL)
draw_venn <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.44, width = 0.8, height = 0.8))
  grid.draw(p)
  popViewport()
  grid.text("GSE93272", x = 0.16, y = 0.92, just = "left",
             gp = gpar(fontsize = 15, fontface = "bold"))
  grid.text("GSE110169", x = 0.84, y = 0.92, just = "right",
             gp = gpar(fontsize = 15, fontface = "bold"))
}
png(file.path(fig, "fig_dataset_gene_overlap_venn.png"), width = 1500, height = 1350, res = 260)
draw_venn(); dev.off()
cairo_pdf(file.path(fig, "fig_dataset_gene_overlap_venn.pdf"), width = 5.8, height = 5.2)
draw_venn(); dev.off()

fwrite(data.table(dataset = c("GSE93272","GSE110169","common"),
                  n_genes = c(length(g93), length(g110), length(common))),
       file.path(tab, "dataset_gene_overlap.csv"))
cat("Wrote fig_dataset_gene_overlap_venn.png + .pdf + dataset_gene_overlap.csv\nDONE\n")
