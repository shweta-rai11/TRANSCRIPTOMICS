#!/usr/bin/env Rscript
# =============================================================================
# new/19new_pathway_enrichment.R
# -----------------------------------------------------------------------------
# GO (biological process) + KEGG over-representation enrichment of the MR-prioritised
# gene sets, per sex (Female 74, Male 55). Background universe = expressed genes
# in the training cohort. Dotplot of the top terms per sex.
# Outputs (new/):
#   results/tables/enrich_GO_BP_{female,male}.csv
#   results/tables/enrich_KEGG_{female,male}.csv
#   results/figures/new/fig_enrich_{female,male}.png/pdf
# =============================================================================
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db)
                  library(ggplot2); library(data.table)})
proc <- "data/processed"; tabN <- "results/tables"; figN <- "results/figures/new"
dir.create(tabN, showWarnings = FALSE, recursive = TRUE)
dir.create(figN, showWarnings = FALSE, recursive = TRUE)

o <- readRDS(file.path(proc, "combined_train.rds"))
universe_sym <- rownames(o$expr)
uni <- suppressWarnings(bitr(universe_sym, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID

enrich_sex <- function(sx) {
  sexlab <- if (sx == "female") "Female" else "Male"
  genes <- fread(file.path(tabN, sprintf("FS_input_%s.csv", sx)))$gene
  eg <- suppressWarnings(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID
  cat(sprintf("%s: %d prioritised genes -> %d mapped to ENTREZ\n", sexlab, length(genes), length(eg)))

  go <- enrichGO(eg, org.Hs.eg.db, ont = "BP", universe = uni,
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  kg <- tryCatch(setReadable(enrichKEGG(eg, organism = "hsa", universe = uni,
                 pvalueCutoff = 0.05), org.Hs.eg.db, "ENTREZID"), error = function(e) NULL)

  if (!is.null(go) && nrow(as.data.frame(go)))
    fwrite(as.data.frame(go), file.path(tabN, sprintf("enrich_GO_BP_%s.csv", sx)))
  if (!is.null(kg) && nrow(as.data.frame(kg)))
    fwrite(as.data.frame(kg), file.path(tabN, sprintf("enrich_KEGG_%s.csv", sx)))

  # combine top terms for a dotplot (GO BP + KEGG)
  mk <- function(x, src) { if (is.null(x)) return(NULL); d <- as.data.frame(x)
    if (!nrow(d)) return(NULL); d <- head(d[order(d$p.adjust), ], 10)
    data.table(source = src, term = d$Description, count = d$Count,
               gene_ratio = sapply(strsplit(d$GeneRatio, "/"), function(v) as.numeric(v[1])/as.numeric(v[2])),
               padj = d$p.adjust) }
  tp <- rbindlist(list(mk(go, "GO:BP"), mk(kg, "KEGG")), fill = TRUE)
  if (is.null(tp) || !nrow(tp)) { cat(sprintf("  %s: no significant terms\n", sexlab)); return(invisible()) }
  tp[, term := factor(term, levels = rev(tp[order(source, gene_ratio)]$term))]

  g <- ggplot(tp, aes(x = gene_ratio, y = term, size = count, colour = padj)) +
    geom_point() +
    facet_grid(source ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#B2182B", high = "#2166AC", name = "p.adjust") +
    scale_size_continuous(name = "Count", range = c(2.5, 7)) +
    labs(x = "Gene ratio", y = NULL, subtitle = sprintf("%s MR-prioritised genes (n=%d)", sexlab, length(genes))) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8), strip.text.y = element_text(angle = 0, face = "bold"))
  ggsave(file.path(figN, sprintf("fig_enrich_%s.png", sx)), g, width = 8, height = 7, dpi = 300)
  ggsave(file.path(figN, sprintf("fig_enrich_%s.pdf", sx)), g, width = 8, height = 7)
  cat(sprintf("  wrote fig_enrich_%s + enrichment tables\n", sx))
}
enrich_sex("female")
enrich_sex("male")
cat("DONE\n")
