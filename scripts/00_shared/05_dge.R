#!/usr/bin/env Rscript
# Differential expression (limma) and functional enrichment (GO/KEGG/GSEA), RA vs HC, by sex.
suppressMessages({
  library(limma); library(data.table)
  library(clusterProfiler); library(org.Hs.eg.db)
})
options(stringsAsFactors = FALSE)
set.seed(1234)

proc <- "data/processed"; tab <- "results/tables"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)

FC <- 0.1; PVAL <- 0.05                       # DE thresholds: |log2FC| and (volcano) raw-p line
PADJ <- 0.05                                  # significance + enrichment BH (FDR) cutoff

o    <- readRDS(file.path(proc, "combined_train.rds"))
expr <- o$expr                                # genes x training samples (70% split)
meta <- o$meta
meta$group <- factor(meta$group, levels = c("HC","RA"))   # HC = reference

# Part 1: differential expression (limma), RA vs HC on a sample subset
run_de <- function(cols, label) {
  md <- meta[match(cols, meta$sample), ]
  E  <- expr[, cols, drop = FALSE]
  design <- model.matrix(~ group, data = md)
  # Empirical array quality weights: down-weight noisy arrays before eBayes moderation
  aw  <- limma::arrayWeights(E, design)
  fit <- eBayes(lmFit(E, design, weights = aw))
  tt  <- topTable(fit, coef = "groupRA", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt$sig  <- abs(tt$logFC) > FC & tt$adj.P.Val < PADJ        # |log2FC|>0.1 & FDR<0.05
  tt$dir  <- ifelse(!tt$sig, "ns", ifelse(tt$logFC > 0, "Up in RA", "Down in RA"))
  tt <- as.data.table(tt)[, .(gene, logFC, AveExpr, t, P.Value, adj.P.Val, sig, dir)]
  attr(tt, "n") <- c(RA = sum(md$group == "RA"), HC = sum(md$group == "HC"))
  attr(tt, "aw") <- data.table(comparison = label, sample = cols,
                               dataset = md$dataset, group = as.character(md$group),
                               sex = md$sex, rin = md$rin, weight = as.numeric(aw))
  cat(sprintf("  %-8s: n=%d (RA=%d, HC=%d) | sig=%d (up=%d, down=%d)\n",
              label, length(cols), sum(md$group=="RA"), sum(md$group=="HC"),
              sum(tt$sig), sum(tt$dir=="Up in RA"), sum(tt$dir=="Down in RA")))
  cat(sprintf("            arrayWeights: median=%.2f range=[%.2f, %.2f] | n(w<0.5)=%d n(w>2)=%d\n",
              median(aw), min(aw), max(aw), sum(aw < 0.5), sum(aw > 2)))
  tt
}

cat("DGE (limma), RA vs HC; sig = |log2FC|>", FC, "& FDR<", PADJ, "\n")
res <- list(
  All    = run_de(meta$sample,                      "All"),
  Female = run_de(meta$sample[meta$sex == "F"],     "Female"),
  Male   = run_de(meta$sample[meta$sex == "M"],     "Male"))

# Save DEG tables
for (nm in names(res)) {
  fwrite(res[[nm]],                    file.path(tab, sprintf("DEG_%s_full.csv", tolower(nm))))
  fwrite(res[[nm]][sig == TRUE][order(P.Value)],
                                       file.path(tab, sprintf("DEG_%s_significant.csv", tolower(nm))))
}

summ <- rbindlist(lapply(names(res), function(nm) {
  d <- res[[nm]]; n <- attr(d, "n")
  data.table(comparison = nm, RA = n["RA"], HC = n["HC"],
             genes_tested = nrow(d), significant = sum(d$sig),
             up_in_RA = sum(d$dir == "Up in RA"), down_in_RA = sum(d$dir == "Down in RA"))
}))
fwrite(summ, file.path(tab, "DEG_summary.csv"))
cat("\n===== DEG SUMMARY =====\n"); print(summ)

# Save array quality weights for the supplement
awt <- rbindlist(lapply(res, function(d) attr(d, "aw")))
fwrite(awt, file.path(tab, "DEG_array_weights.csv"))
cat("\n===== ARRAY QUALITY WEIGHTS (Ritchie et al. 2006) =====\n")
print(awt[, .(n = .N, median_w = round(median(weight), 3),
              min_w = round(min(weight), 3), max_w = round(max(weight), 3),
              n_low = sum(weight < 0.5)), by = .(comparison, dataset)])
cat("\nLowest-weighted arrays overall (candidate poor-quality):\n")
print(head(awt[comparison == "All"][order(weight)][
      , .(sample, dataset, group, sex, rin, weight = round(weight, 3))], 10))

# Sex-specific DEG lists removed: report sex-stratified DEG_female_*/DEG_male_*.csv instead
# Part 1b: TREAT sensitivity check (tests |logFC| against a threshold directly, McCarthy & Smyth 2009)
cat("\n===== TREAT SENSITIVITY (fold-change threshold inside the test) =====\n")
treat_tab <- rbindlist(lapply(c("All", "Female", "Male"), function(nm) {
  cols <- switch(nm, All = meta$sample,
                 Female = meta$sample[meta$sex == "F"],
                 Male   = meta$sample[meta$sex == "M"])
  md <- meta[match(cols, meta$sample), ]
  E  <- expr[, cols, drop = FALSE]
  d  <- model.matrix(~ group, data = md)
  aw <- limma::arrayWeights(E, d)
  f  <- lmFit(E, d, weights = aw)
  eb <- topTable(eBayes(f), coef = "groupRA", number = Inf, sort.by = "none")
  rbindlist(lapply(c(0.1, log2(1.2), 0.5), function(L) {
    tt <- topTreat(treat(f, lfc = L), coef = "groupRA", number = Inf, sort.by = "none")
    data.table(comparison = nm, lfc_threshold = round(L, 3),
               fold_change = sprintf("%.2fx", 2^L),
               posthoc_filter_n = sum(eb$adj.P.Val < PADJ & abs(eb$logFC) > L),
               treat_n = sum(tt$adj.P.Val < PADJ))
  }))
}))

fwrite(treat_tab, file.path(tab, "DEG_treat_sensitivity.csv"))
print(treat_tab)
cat("treat_n << posthoc_filter_n is EXPECTED: treat() demands positive evidence\n")
cat("that the effect EXCEEDS the threshold, not merely that it is non-zero.\n")

sigF <- res$Female[sig == TRUE]$gene
sigM <- res$Male[sig == TRUE]$gene
cat(sprintf("\nSignificant DEGs -- female: %d | male: %d | overlapping: %d\n",
            length(sigF), length(sigM), length(intersect(sigF, sigM))))

saveRDS(list(res = res, thr = list(FC = FC, PVAL = PVAL)),
        file.path(proc, "dge_results.rds"))
cat("Saved DEG tables + data/processed/dge_results.rds\n")

# Part 2: functional enrichment of each DEG set (GO + KEGG)
cat("\nEnrichment (GO BP/CC/MF + KEGG, BH adj.P <", PADJ, ")\n")
enrich_group <- function(comp) {
  genes <- res[[comp]][sig == TRUE]$gene
  eg <- suppressWarnings(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID
  eg <- unique(eg)
  go <- tryCatch(as.data.frame(enrichGO(eg, org.Hs.eg.db, ont = "ALL",
            pAdjustMethod = "BH", pvalueCutoff = PADJ, qvalueCutoff = 0.2,
            readable = TRUE)), error = function(e) NULL)
  kg <- tryCatch(as.data.frame(enrichKEGG(eg, organism = "hsa",
            pAdjustMethod = "BH", pvalueCutoff = PADJ, qvalueCutoff = 0.2)),
            error = function(e) NULL)
  terms <- c(if (!is.null(go) && nrow(go)) go$ID, if (!is.null(kg) && nrow(kg)) kg$ID)
  tab_df <- rbindlist(list(
    if (!is.null(go) && nrow(go)) data.table(group = comp, source = go$ONTOLOGY,
        ID = go$ID, Description = go$Description, GeneRatio = go$GeneRatio,
        Count = go$Count, p.adjust = go$p.adjust, geneID = go$geneID),
    if (!is.null(kg) && nrow(kg)) data.table(group = comp, source = "KEGG",
        ID = kg$ID, Description = kg$Description, GeneRatio = kg$GeneRatio,
        Count = kg$Count, p.adjust = kg$p.adjust, geneID = kg$geneID)),
    fill = TRUE)
  cat(sprintf("  %-6s: %d DEGs -> %d enriched terms (GO %d, KEGG %d)\n", comp,
      length(genes), length(unique(terms)),
      if (is.null(go)) 0 else nrow(go), if (is.null(kg)) 0 else nrow(kg)))
  list(terms = unique(terms), table = tab_df)
}

enr <- lapply(c("Female", "Male", "All"), enrich_group)
names(enr) <- c("Female", "Male", "All")
sets <- lapply(enr, `[[`, "terms")

for (g in names(enr))
  fwrite(enr[[g]]$table, file.path(tab, sprintf("dge_enriched_terms_%s.csv", tolower(g))))

ov <- data.table(
  Female_only = length(setdiff(sets$Female, union(sets$Male, sets$All))),
  Male_only   = length(setdiff(sets$Male,   union(sets$Female, sets$All))),
  All_only    = length(setdiff(sets$All,    union(sets$Female, sets$Male))),
  Female_Male = length(setdiff(intersect(sets$Female, sets$Male), sets$All)),
  Female_All  = length(setdiff(intersect(sets$Female, sets$All),  sets$Male)),
  Male_All    = length(setdiff(intersect(sets$Male,   sets$All),  sets$Female)),
  all_three   = length(Reduce(intersect, sets)))
cat("\nEnriched-term overlap (Female/Male/All):\n"); print(t(ov))

saveRDS(list(sets = sets, res = enr), file.path(proc, "dge_enrich.rds"))
cat("\nSaved dge_enrich.rds + dge_enriched_terms_*.csv\n")

# Part 3: KEGG GSEA, genes ranked by the limma t-statistic (RA vs HC), per comparison
cat("\nKEGG GSEA (genes ranked by limma t-statistic)\n")
gsea_kegg <- function(comp) {
  d  <- as.data.frame(res[[comp]])
  eg <- suppressWarnings(bitr(d$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db))
  d  <- merge(d, eg, by.x = "gene", by.y = "SYMBOL")
  d  <- d[order(-abs(d$t)), ]; d <- d[!duplicated(d$ENTREZID), ]   # 1 rank per gene
  ranks <- sort(setNames(d$t, d$ENTREZID), decreasing = TRUE)
  gk <- tryCatch(gseKEGG(ranks, organism = "hsa", pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, verbose = FALSE, seed = TRUE),
                 error = function(e) NULL)
  gdf <- if (!is.null(gk)) as.data.frame(gk) else NULL
  if (is.null(gdf) || !nrow(gdf)) { cat(sprintf("  %-6s: 0 GSEA pathways\n", comp)); return(NULL) }
  cat(sprintf("  %-6s: %d GSEA pathways (%d activated / %d suppressed in RA)\n",
              comp, nrow(gdf), sum(gdf$NES > 0), sum(gdf$NES < 0)))
  tab <- as.data.table(gdf)[, .(group = comp, ID, Description, setSize, NES,
                         pvalue, p.adjust, core_enrichment)]
  list(obj = gk, table = tab)                    # keep the gseaResult for running-score plots
}
gres <- lapply(c("Female", "Male", "All"), gsea_kegg)
names(gres) <- c("Female", "Male", "All")
gres <- gres[!vapply(gres, is.null, logical(1))]
gsea      <- rbindlist(lapply(gres, `[[`, "table"), fill = TRUE)
gsea_objs <- lapply(gres, `[[`, "obj")
fwrite(gsea, file.path(tab, "gsea_kegg.csv"))
saveRDS(gsea, file.path(proc, "dge_gsea.rds"))
saveRDS(gsea_objs, file.path(proc, "dge_gsea_obj.rds"))   # for enrichplot::gseaplot2
cat("\nSaved dge_gsea.rds + dge_gsea_obj.rds + gsea_kegg.csv\nDONE\n")

