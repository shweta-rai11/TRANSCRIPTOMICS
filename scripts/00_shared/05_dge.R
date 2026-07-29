#!/usr/bin/env Rscript
# =============================================================================
# Differential expression AND functional enrichment (ANALYSIS).
#
# ---- References (methods) --------------------------------------------------
# Differential expression (limma linear models + empirical-Bayes moderation):
#   # Smyth GK. Linear models and empirical Bayes methods for assessing differential
#   #   expression in microarray experiments. Stat Appl Genet Mol Biol 2004;3:Article 3.
#   # Ritchie ME, et al. limma powers differential expression analyses for RNA-seq
#   #   and microarray studies. Nucleic Acids Res 2015;43(7):e47.
# Multiple-testing correction (Benjamini-Hochberg FDR):
#   # Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
#   #   powerful approach to multiple testing. J R Stat Soc B 1995;57(1):289-300.
# Over-representation & GSEA framework (clusterProfiler):
#   # Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological
#   #   themes among gene clusters. OMICS 2012;16(5):284-287.
#   # Wu T, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting
#   #   omics data. Innovation (Camb) 2021;2(3):100141.
# Gene-set / pathway annotations:
#   # Ashburner M, et al. Gene Ontology: tool for the unification of biology.
#   #   Nat Genet 2000;25(1):25-29.  (GO)
#   # Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes.
#   #   Nucleic Acids Res 2000;28(1):27-30.  (KEGG)
# Gene Set Enrichment Analysis (ranked test):
#   # Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for
#   #   interpreting genome-wide expression profiles. PNAS 2005;102(43):15545-15550.
# =============================================================================
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

# =============================================================================
# PART 1 — DIFFERENTIAL EXPRESSION (limma), RA vs HC on a sample subset
# =============================================================================
run_de <- function(cols, label) {
  md <- meta[match(cols, meta$sample), ]
  E  <- expr[, cols, drop = FALSE]
  design <- model.matrix(~ group, data = md)
  # ---- empirical array quality weights --------------------------------------
  # REML per-array weights: arrays whose residuals are consistently large across
  # ALL genes are down-weighted, rather than excluded by an arbitrary QC rule.
  # Motivated here by (i) the two-platform merge (GPL570 + GPL13667), (ii) variable
  # RNA integrity (RIN 7.0-9.0 in GSE93272) and (iii) the small male stratum
  # (n=38), where a single poor-quality array has high leverage on logFC.
  # Weights act PER ARRAY and are orthogonal to the eBayes variance moderation,
  # which borrows strength PER GENE across the transcriptome.
  #   # Ritchie ME, Diyagama D, Neilson J, et al. Empirical array quality weights
  #   #   in the analysis of microarray data. BMC Bioinformatics 2006;7:261.
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

# ---- save DEG tables --------------------------------------------------------
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

# ---- array quality weights: save for the supplement -------------------------
awt <- rbindlist(lapply(res, function(d) attr(d, "aw")))
fwrite(awt, file.path(tab, "DEG_array_weights.csv"))
cat("\n===== ARRAY QUALITY WEIGHTS (Ritchie et al. 2006) =====\n")
print(awt[, .(n = .N, median_w = round(median(weight), 3),
              min_w = round(min(weight), 3), max_w = round(max(weight), 3),
              n_low = sum(weight < 0.5)), by = .(comparison, dataset)])
cat("\nLowest-weighted arrays overall (candidate poor-quality):\n")
print(head(awt[comparison == "All"][order(weight)][
      , .(sample, dataset, group, sex, rin, weight = round(weight, 3))], 10))

# ---- REMOVED (sex-specific) -------------------------------------------------
# Was: "female-only" / "male-only" significant DEG lists built as a set
# difference (setdiff(sigF, sigM)) -> DEG_female_specific.csv / DEG_male_specific.csv.
# A set difference of two separately-thresholded lists is a SEX-SPECIFICITY
# claim and is not a valid test of one (a gene can miss significance in men
# purely from the smaller male n). This chapter reports SEX-STRATIFIED DEGs
# only: DEG_female_*.csv and DEG_male_*.csv, each RA vs Control within that sex.
# The overlap is still reported descriptively in the counts below.
# =============================================================================
# TREAT SENSITIVITY ANALYSIS (this block makes the McCarthy & Smyth citation
# below TRUE. It previously said "evaluated as sensitivity check" with no
# treat() call anywhere in the script.)
#
# WHY: our significance rule filters on |logFC| AFTER testing against H0: lfc=0.
# treat() instead tests H0: |true logFC| <= lfc, i.e. the fold-change threshold
# is inside the hypothesis, which is the statistically correct way to demand a
# minimum effect size. It is deliberately CONSERVATIVE.
#   # McCarthy DJ, Smyth GK. Bioinformatics 2009;25(6):765-771.
# =============================================================================
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

# =============================================================================
# PART 2 — FUNCTIONAL ENRICHMENT of each DEG set (GO + KEGG)
# =============================================================================
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

# =============================================================================
# PART 3 — KEGG GSEA (ranked whole-genome test, per comparison)
#   Genes ranked by the limma moderated t-statistic (RA vs HC); gseKEGG returns
#   a signed NES per pathway (NES > 0 = up in RA, NES < 0 = down in RA).
# =============================================================================
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

# =============================================================================
# FULL PIPELINE METHODS REFERENCES (01 -> 05)
# -----------------------------------------------------------------------------
# DATA ACQUISITION
#   Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus. Nucleic Acids Res 2002;30(1):207-210.
#   Barrett T, et al. NCBI GEO: archive for functional genomics data sets - update. Nucleic Acids Res 2013;41(D1):D991-D995.
#   Davis S, Meltzer PS. GEOquery: a bridge between GEO and BioConductor. Bioinformatics 2007;23(14):1846-1847.
#   Huber W, et al. Orchestrating high-throughput genomic analysis with Bioconductor. Nat Methods 2015;12(2):115-121.
# PREPROCESSING
#   Miller JA, et al. Strategies for aggregating gene expression data: the collapseRows R function. BMC Bioinformatics 2011;12:322.
#   Bolstad BM, et al. A comparison of normalization methods for high density oligonucleotide array data. Bioinformatics 2003;19(2):185-193.
#   Ritchie ME, et al. limma powers differential expression analyses for RNA-seq and microarray studies. Nucleic Acids Res 2015;43(7):e47.
#   Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 2007;8(1):118-127.
#   Leek JT, et al. The sva package for removing batch effects. Bioinformatics 2012;28(6):882-883.
#   Nygaard V, Rodland EA, Hovig E. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics 2016;17(1):29-39.
# LEAKAGE-SAFE 70:30 SPLIT
#   Ambroise C, McLachlan GJ. Selection bias in gene extraction on the basis of microarray gene-expression data. PNAS 2002;99(10):6562-6566.
#   Simon R, et al. Pitfalls in the use of DNA microarray data for classification. J Natl Cancer Inst 2003;95(1):14-18.
#   Kaufman S, et al. Leakage in data mining: formulation, detection, and avoidance. ACM TKDD 2012;6(4):Article 15.
#   Kuhn M. Building predictive models in R using the caret package. J Stat Softw 2008;28(5):1-26.
# DIFFERENTIAL EXPRESSION
#   Smyth GK. Linear models and empirical Bayes methods for assessing differential expression. Stat Appl Genet Mol Biol 2004;3:Article 3.
#   Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc B 1995;57(1):289-300.   <- significance rule: |log2FC|>0.1 & FDR<0.05
#   McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. Bioinformatics 2009;25(6):765-771. (evaluated as sensitivity check) doi:10.1093/bioinformatics/btp053
# FUNCTIONAL ENRICHMENT
#   Yu G, et al. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 2012;16(5):284-287.
#   Wu T, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation (Camb) 2021;2(3):100141.
#   Ashburner M, et al. Gene Ontology: tool for the unification of biology. Nat Genet 2000;25(1):25-29.
#   Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res 2000;28(1):27-30.
#   Subramanian A, et al. Gene set enrichment analysis. PNAS 2005;102(43):15545-15550.
# WGCNA + SEX-STRATIFIED ANALYSIS
#   Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008;9:559.
#   Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. Stat Appl Genet Mol Biol 2005;4:Article 17.
#   Langfelder P, et al. Is my network module preserved and reproducible? PLoS Comput Biol 2011;7(1):e1001057.
#   Gelman A, Stern H. The difference between "significant" and "not significant" is not itself statistically significant. Am Stat 2006;60(4):328-331.
# =============================================================================
