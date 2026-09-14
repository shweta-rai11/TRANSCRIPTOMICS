#!/usr/bin/env Rscript
# WGCNA: weighted gene co-expression network analysis, disease-module detection, and figures, in one script

options(stringsAsFactors = FALSE)
suppressMessages({
  library(WGCNA); library(data.table); library(ggplot2); library(pheatmap)
  library(igraph); library(ggrepel); library(reshape2); library(patchwork)
  library(scales); library(clusterProfiler); library(org.Hs.eg.db)
})

# Step 0: configuration - every tunable lives here
CFG <- list(
  seed              = 1234,   # global RNG seed
  # no variance filter: any gene removed pre-network can never enter a disease module
  var_quantile      = 0,
  drop_outliers     = TRUE,   # actually remove sample outliers (old script only printed them)
  outlier_sd        = 3,      # outlier cut = mean(merge height) + k*sd
  network_type      = "signed",
  tom_type          = "signed",
  # pearson, with outlier robustness handled at the sample level (outliers dropped in STEP 2) rather than via bicor
  cor_type          = "pearson",
  max_p_outliers    = 0.05,   # retained; used only when cor_type = "bicor"
  # soft power fixed at 12 for well-conditioned connectivity (see WGCNA_01_soft_threshold.csv for the R^2/connectivity tradeoff)
  power_mode        = 12,
  rsq_cut           = 0.85,   # scale-free R^2 target reported for the auto estimate
  max_block_size    = 20000,  # keep all genes in ONE block (default 5000 splits them)
  min_module_size   = 30,
  merge_cut_height  = 0.25,   # merge MEs correlated > 0.75
  disease_trait     = "RA",
  # --- disease-module rule: DATA-DRIVEN, not colour-based (rule R1) ----------
  dm_min_abs_cor    = 0.50,   # |cor(ME, RA)| threshold
  dm_max_p          = 1e-8,   # module-trait p threshold
  hub_kME           = 0.80,   # hub gene: |kME| >
  hub_GS            = 0.20,   # hub gene: |GS|  >
  pres_permutations = 200,    # modulePreservation permutations (>=200 for publication)
  do_preservation   = TRUE,
  use_cache         = TRUE    # cache is keyed to a parameter hash (rule R3)
)
set.seed(CFG$seed)

proc <- "data/processed"; fig <- "results/figures"; tab <- "results/tables"
for (d in c(proc, fig, tab)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# parameter + data fingerprint -> cache key, so a stale network is never silently reused
.tmp_in <- readRDS(file.path("data/processed", "combined_train.rds"))
data_fp <- paste(dim(.tmp_in$expr)[1], dim(.tmp_in$expr)[2],
                 sprintf("%.10g", sum(.tmp_in$expr)),
                 sprintf("%.10g", sum(abs(.tmp_in$expr[, 1]))), sep = ":")
rm(.tmp_in)
net_key <- substr(paste(
  CFG$seed, CFG$var_quantile, CFG$drop_outliers, CFG$outlier_sd,
  CFG$network_type, CFG$tom_type, CFG$cor_type, CFG$max_p_outliers,
  CFG$power_mode, CFG$rsq_cut,
  CFG$max_block_size, CFG$min_module_size, CFG$merge_cut_height,
  data_fp,                                   # <- input data fingerprint
  sep = "|"), 1, 300)
key_hash <- sprintf("%08x", sum(utf8ToInt(net_key) * seq_len(nchar(net_key))))
cache_of <- function(what) file.path(proc, sprintf("wgcna_%s_%s.rds", what, key_hash))

# map corType once to the function-name form used by pickSoftThreshold/adjacency
COR_FNC  <- if (identical(CFG$cor_type, "bicor")) "bicor" else "cor"
COR_OPTS <- if (identical(CFG$cor_type, "bicor"))
              list(maxPOutliers = CFG$max_p_outliers) else list(use = "p")

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

hdr("STEP 0  CONFIGURATION")
say("cache key            : %s", key_hash)
say("network / TOM type   : %s / %s", CFG$network_type, CFG$tom_type)
say("variance filter      : %s", if (CFG$var_quantile <= 0) "DISABLED (all genes kept)"
    else sprintf("drop below the %.0f%% quantile", 100 * CFG$var_quantile))
say("correlation          : %s (maxPOutliers = %.2f)", CFG$cor_type, CFG$max_p_outliers)
say("disease-module rule  : |cor(ME, %s)| >= %.2f  AND  p < %g",
    CFG$disease_trait, CFG$dm_min_abs_cor, CFG$dm_max_p)

# Step 1: load the 70% training split (leakage-safe; holdout never seen here)
hdr("STEP 1  LOAD DATA")
o <- readRDS(file.path(proc, "combined_train.rds"))
stopifnot(o$role == "combined_training_70pct")
expr_all <- as.matrix(o$expr)                       # genes x samples
meta_all <- o$meta
rownames(meta_all) <- meta_all$sample
stopifnot(identical(colnames(expr_all), meta_all$sample))
say("expression : %d genes x %d samples", nrow(expr_all), ncol(expr_all))
say("groups     : RA=%d  HC=%d", sum(meta_all$group == "RA"), sum(meta_all$group == "HC"))
say("sex        : F=%d   M=%d",  sum(meta_all$sex == "F"),    sum(meta_all$sex == "M"))

# Step 2: QC - goodSamplesGenes, then sample-outlier detection and removal
hdr("STEP 2  QUALITY CONTROL")
datExpr0 <- t(expr_all)                              # WGCNA wants samples x genes
gsg <- goodSamplesGenes(datExpr0, verbose = 0)
say("goodSamplesGenes: %d/%d genes OK, %d/%d samples OK",
    sum(gsg$goodGenes), length(gsg$goodGenes), sum(gsg$goodSamples), length(gsg$goodSamples))
if (!gsg$allOK) datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

sample_tree <- hclust(dist(datExpr0), method = "average")
cut_height  <- mean(sample_tree$height) + CFG$outlier_sd * sd(sample_tree$height)
outlier_cl  <- cutreeStatic(sample_tree, cutHeight = cut_height, minSize = 1)
keep_samp   <- outlier_cl == names(sort(table(outlier_cl), decreasing = TRUE))[1]
n_out       <- sum(!keep_samp)
say("outlier cut height : %.2f  (mean + %g*sd of merge heights)", cut_height, CFG$outlier_sd)
say("sample outliers    : %d  -> %s", n_out,
    if (CFG$drop_outliers && n_out > 0) "REMOVED" else "retained")
outlier_ids <- rownames(datExpr0)[!keep_samp]
if (n_out > 0) say("  %s", paste(outlier_ids, collapse = ", "))
if (CFG$drop_outliers && n_out > 0) datExpr0 <- datExpr0[keep_samp, , drop = FALSE]

# Step 3: gene variance filter (STEP 14 quantifies sex-DEGs lost to it)
hdr("STEP 3  GENE FILTER")
gene_variance <- apply(datExpr0, 2, var)
if (CFG$var_quantile <= 0) {
  var_thresh    <- -Inf
  datExpr       <- datExpr0
  genes_dropped <- character(0)
  say("variance filter : DISABLED (var_quantile = 0)")
  say("  Rationale: filtering before the network permanently excludes those genes")
  say("  from every disease module, hence from `disease-module n sex-DEG`.")
} else {
  var_thresh    <- quantile(gene_variance, CFG$var_quantile)
  datExpr       <- datExpr0[, gene_variance > var_thresh, drop = FALSE]
  genes_dropped <- setdiff(colnames(datExpr0), colnames(datExpr))
  say("variance threshold (%.0fth pct) : %.5f", 100 * CFG$var_quantile, var_thresh)
}
say("genes retained : %d of %d  (dropped %d)",
    ncol(datExpr), ncol(datExpr0), length(genes_dropped))
say("final matrix   : %d samples x %d genes", nrow(datExpr), ncol(datExpr))

meta <- meta_all[rownames(datExpr), , drop = FALSE]
stopifnot(identical(rownames(datExpr), meta$sample))

# Step 4: trait matrix (binary indicators aligned to datExpr rows)
hdr("STEP 4  TRAIT MATRIX")
traits <- data.frame(
  RA        = as.numeric(meta$group == "RA"),
  HC        = as.numeric(meta$group == "HC"),
  Female    = as.numeric(meta$sex   == "F"),
  Male      = as.numeric(meta$sex   == "M"),
  RA_Female = as.numeric(meta$group == "RA" & meta$sex == "F"),
  RA_Male   = as.numeric(meta$group == "RA" & meta$sex == "M"),
  HC_Female = as.numeric(meta$group == "HC" & meta$sex == "F"),
  HC_Male   = as.numeric(meta$group == "HC" & meta$sex == "M"),
  Age       = as.numeric(meta$age),
  row.names = rownames(datExpr))
stopifnot(identical(rownames(traits), rownames(datExpr)))
print(colSums(traits[, setdiff(names(traits), "Age")]))

# Step 5: soft-thresholding power (data-driven estimate by default; FAQ table reported alongside)
hdr("STEP 5  SOFT-THRESHOLDING POWER")
enableWGCNAThreads()
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = CFG$network_type,
                         corFnc = COR_FNC, corOptions = COR_OPTS,
                         RsquaredCut = CFG$rsq_cut, verbose = 0)
fitI <- sft$fitIndices
sfR2 <- -sign(fitI[, 3]) * fitI[, 2]
faq_power <- if (nrow(datExpr) < 20) 18 else if (nrow(datExpr) <= 30) 16 else
             if (nrow(datExpr) <= 40) 14 else 12
auto_power <- sft$powerEstimate
fit_fails  <- is.na(auto_power)

soft_power <- if (is.numeric(CFG$power_mode)) CFG$power_mode else
              if (CFG$power_mode == "faq") faq_power else
              if (!fit_fails) auto_power else faq_power

say("scale-free fit reaches R^2 >= %.2f : %s", CFG$rsq_cut, if (fit_fails) "NO" else "YES")
say("pickSoftThreshold estimate        : %s", ifelse(fit_fails, "NA", auto_power))
say("WGCNA FAQ table (n=%d, signed)     : %d   [applies only if fit FAILS]",
    nrow(datExpr), faq_power)
say("POWER USED                        : %d   (mode = %s)", soft_power, CFG$power_mode)
say("  scale-free R^2 at chosen power   : %.3f", sfR2[fitI[, 1] == soft_power])
say("  mean connectivity at that power  : %.1f", fitI[fitI[, 1] == soft_power, 5])
print(data.frame(power = fitI[, 1], signed_R2 = round(sfR2, 3),
                 slope = round(fitI[, 3], 2), mean_k = round(fitI[, 5], 1)))
fwrite(data.table(power = fitI[, 1], signed_R2 = sfR2, slope = fitI[, 3],
                  mean_k = fitI[, 5], median_k = fitI[, 6], max_k = fitI[, 7]),
       file.path(tab, "WGCNA_01_soft_threshold.csv"))

# Step 6: network construction + module detection (WGCNA::cor masking restored via on.exit)
hdr("STEP 6  NETWORK CONSTRUCTION")
build_net <- function(X, pw) {
  cor <- WGCNA::cor
  on.exit({ cor <- stats::cor }, add = TRUE)
  blockwiseModules(
    X, power = pw,
    networkType = CFG$network_type, TOMType = CFG$tom_type,
    corType = CFG$cor_type, maxPOutliers = CFG$max_p_outliers,
    maxBlockSize = CFG$max_block_size, minModuleSize = CFG$min_module_size,
    reassignThreshold = 0, mergeCutHeight = CFG$merge_cut_height,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = FALSE, randomSeed = CFG$seed, verbose = 0)
}
net_file <- cache_of("net")
if (CFG$use_cache && file.exists(net_file)) {
  net <- readRDS(net_file); say("loaded cached network (key %s)", key_hash)
} else {
  say("building network (%d genes, power %d) ...", ncol(datExpr), soft_power)
  net <- build_net(datExpr, soft_power)
  saveRDS(net, net_file)
}
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
mod_sizes <- sort(table(moduleColors), decreasing = TRUE)
n_modules <- length(setdiff(names(mod_sizes), "grey"))
say("modules detected: %d (plus grey = %d unassigned genes)",
    n_modules, if ("grey" %in% names(mod_sizes)) mod_sizes[["grey"]] else 0)
print(mod_sizes)

# Step 7: module eigengenes
hdr("STEP 7  MODULE EIGENGENES")
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
say("eigengenes: %d modules x %d samples", ncol(MEs), nrow(MEs))

# Step 8: module-trait correlation
hdr("STEP 8  MODULE-TRAIT CORRELATION")
nSamples <- nrow(datExpr)
mtCor <- stats::cor(MEs, traits, use = "pairwise.complete.obs")
mtP   <- corPvalueStudent(mtCor, nSamples)

mt_tab <- data.table(
  module    = sub("^ME", "", rownames(mtCor)),
  size      = as.integer(mod_sizes[sub("^ME", "", rownames(mtCor))]),
  cor_RA    = round(mtCor[, CFG$disease_trait], 3),
  p_RA      = signif(mtP[, CFG$disease_trait], 3),
  cor_Male  = round(mtCor[, "Male"], 3),
  p_Male    = signif(mtP[, "Male"], 3),
  direction = ifelse(mtCor[, CFG$disease_trait] > 0, "UP in RA", "DOWN in RA"))
setorder(mt_tab, -cor_RA)
fwrite(mt_tab, file.path(tab, "WGCNA_02_module_trait.csv"))
print(mt_tab)

# Step 9: disease module selection, data-driven by correlation with RA (not by colour name)
hdr("STEP 9  DISEASE MODULE SELECTION (data-driven)")
dm <- mt_tab[module != "grey" &
             abs(cor_RA) >= CFG$dm_min_abs_cor & p_RA < CFG$dm_max_p]
if (!nrow(dm)) {
  say("WARNING: no module meets |cor|>=%.2f & p<%g; falling back to the single",
      CFG$dm_min_abs_cor, CFG$dm_max_p)
  say("         most RA-correlated non-grey module.")
  dm <- mt_tab[module != "grey"][which.max(abs(cor_RA))]
}
disease_modules <- dm$module
disease_genes   <- names(moduleColors)[moduleColors %in% disease_modules]
say("rule            : |cor(ME,%s)| >= %.2f AND p < %g",
    CFG$disease_trait, CFG$dm_min_abs_cor, CFG$dm_max_p)
say("DISEASE MODULES : %s", paste(sprintf("%s (r=%+.2f, n=%d)",
    dm$module, dm$cor_RA, dm$size), collapse = " | "))
say("disease genes   : %d", length(disease_genes))
fwrite(dm, file.path(tab, "WGCNA_03_disease_modules.csv"))
fwrite(data.table(gene = disease_genes,
                  module = moduleColors[disease_genes]),
       file.path(tab, "WGCNA_04_disease_module_genes.csv"))

# Step 10: kME (module membership) and GS (gene significance for RA)
hdr("STEP 10  kME AND GENE SIGNIFICANCE")
kME   <- signedKME(datExpr, MEs)
GS_RA <- as.numeric(stats::cor(datExpr, traits[[CFG$disease_trait]],
                               use = "pairwise.complete.obs"))
names(GS_RA) <- colnames(datExpr)
gene_tab <- data.table(
  gene   = colnames(datExpr),
  module = moduleColors[colnames(datExpr)],
  GS_RA  = round(GS_RA, 3),
  kME_own = round(vapply(colnames(datExpr), function(g)
              kME[g, paste0("kME", moduleColors[[g]])], numeric(1)), 3),
  is_disease_module = moduleColors[colnames(datExpr)] %in% disease_modules)
fwrite(gene_tab, file.path(tab, "WGCNA_05_gene_module_assignment.csv"))
say("gene table: %d genes | in a disease module: %d",
    nrow(gene_tab), sum(gene_tab$is_disease_module))

# Step 11: hub genes within each disease module
hdr("STEP 11  HUB GENES PER DISEASE MODULE")
hub_list <- rbindlist(lapply(disease_modules, function(mc) {
  g <- names(moduleColors)[moduleColors == mc]
  a <- adjacency(datExpr[, g, drop = FALSE], power = soft_power,
                 type = CFG$network_type, corFnc = COR_FNC, corOptions = COR_OPTS)
  diag(a) <- 0
  d <- data.table(module = mc, gene = g,
                  kME = round(kME[g, paste0("kME", mc)], 3),
                  GS_RA = round(GS_RA[g], 3),
                  connectivity = round(rowSums(a), 2))
  setorder(d, -connectivity)
  d[, is_hub := abs(kME) > CFG$hub_kME & abs(GS_RA) > CFG$hub_GS]
  say("  %-11s %5d genes | hubs (|kME|>%.2f & |GS|>%.2f): %d",
      mc, nrow(d), CFG$hub_kME, CFG$hub_GS, sum(d$is_hub))
  d }))
fwrite(hub_list, file.path(tab, "WGCNA_06_disease_module_hubs.csv"))
fwrite(hub_list[is_hub == TRUE][order(module, -kME)],
       file.path(tab, "WGCNA_07_hub_genes_only.csv"))

# Step 12: functional enrichment of the disease-module genes (universe = all network genes)
hdr("STEP 12  ENRICHMENT OF DISEASE-MODULE GENES")
uni <- suppressWarnings(bitr(colnames(datExpr), "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID
dg  <- suppressWarnings(bitr(disease_genes,     "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID
ego <- tryCatch(enrichGO(dg, OrgDb = org.Hs.eg.db, ont = "ALL", universe = uni,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, readable = TRUE),
                error = function(e) NULL)
ekg <- tryCatch(setReadable(enrichKEGG(dg, organism = "hsa", universe = uni,
                                       pAdjustMethod = "BH", pvalueCutoff = 0.05),
                            org.Hs.eg.db, "ENTREZID"),
                error = function(e) { message("KEGG unavailable: ", conditionMessage(e)); NULL })
go_dt <- if (!is.null(ego) && nrow(ego)) as.data.table(ego@result)[p.adjust < 0.05][order(p.adjust)] else data.table()
kg_dt <- if (!is.null(ekg) && nrow(ekg)) as.data.table(ekg@result)[p.adjust < 0.05][order(p.adjust)] else data.table()
if (nrow(go_dt)) fwrite(go_dt, file.path(tab, "WGCNA_08_disease_GO.csv"))
if (nrow(kg_dt)) fwrite(kg_dt, file.path(tab, "WGCNA_09_disease_KEGG.csv"))
say("GO terms (p.adj<0.05)   : %d%s", nrow(go_dt),
    if (nrow(go_dt)) sprintf("  top: %s", go_dt$Description[1]) else "")
say("KEGG pathways (p.adj<.05): %d%s", nrow(kg_dt),
    if (nrow(kg_dt)) sprintf("  top: %s", kg_dt$Description[1]) else "")

# Step 13: sex-stratified networks + module preservation (Female -> Male)
hdr("STEP 13  SEX NETWORKS AND MODULE PRESERVATION")
fs <- rownames(datExpr)[meta$sex == "F"]; ms <- rownames(datExpr)[meta$sex == "M"]
say("female n=%d | male n=%d  (WGCNA FAQ minimum 15, preferred >=20)",
    length(fs), length(ms))
eF <- datExpr[fs, , drop = FALSE]; eM <- datExpr[ms, , drop = FALSE]

sex_file <- cache_of("sexnet")
if (CFG$use_cache && file.exists(sex_file)) {
  sn <- readRDS(sex_file); say("loaded cached sex networks")
} else {
  sn <- list(female = build_net(eF, soft_power), male = build_net(eM, soft_power))
  saveRDS(sn, sex_file)
}
colF <- labels2colors(sn$female$colors); colM <- labels2colors(sn$male$colors)
say("female modules: %d | male modules: %d",
    length(setdiff(unique(colF), "grey")), length(setdiff(unique(colM), "grey")))

pres_dt <- data.table()
if (CFG$do_preservation) {
  gF <- goodSamplesGenes(eF, verbose = 0); gM <- goodSamplesGenes(eM, verbose = 0)
  kp <- gF$goodGenes & gM$goodGenes
  say("genes passing QC in BOTH sexes: %d / %d", sum(kp), length(kp))
  pres_file <- cache_of(sprintf("pres%d", CFG$pres_permutations))
  if (CFG$use_cache && file.exists(pres_file)) {
    pres <- readRDS(pres_file); say("loaded cached preservation")
  } else {
    say("modulePreservation: %d permutations (slow) ...", CFG$pres_permutations)
    # reference colours are the combined-network modules, not the female-network modules
    pres <- modulePreservation(
      list(Female = list(data = eF[, kp, drop = FALSE]),
           Male   = list(data = eM[, kp, drop = FALSE])),
      list(Female = moduleColors[kp]), referenceNetworks = 1,
      nPermutations = CFG$pres_permutations, randomSeed = CFG$seed, verbose = 0)
    saveRDS(pres, pres_file)
  }
  Z <- pres$preservation$Z$ref.Female$inColumnsAlsoPresentIn.Male
  pres_dt <- data.table(module = rownames(Z), moduleSize = Z$moduleSize,
                        Zsummary = round(Z$Zsummary.pres, 2))
  pres_dt[, status := cut(Zsummary, c(-Inf, 2, 10, Inf),
                          labels = c("Not preserved", "Moderate", "Strong"))]
  setorder(pres_dt, Zsummary)
  fwrite(pres_dt, file.path(tab, "WGCNA_10_module_preservation.csv"))
  print(pres_dt)
}

# Step 14: the deliverable - disease-module genes n sex-stratified DEGs
hdr("STEP 14  DISEASE-MODULE  n  SEX-DEG   (the candidate step)")
deg_file <- file.path(proc, "dge_results.rds")
cand <- list()
if (!file.exists(deg_file)) {
  say("dge_results.rds not found - run 05_dge.R first. Skipping STEP 14.")
} else {
  D <- readRDS(deg_file)
  say("DEG cutoff in force: |log2FC| > %.2f & FDR < 0.05", D$thr$FC)
  ov <- rbindlist(lapply(c("Female", "Male"), function(sx) {
    sig  <- as.data.table(D$res[[sx]])[sig == TRUE]
    hit  <- intersect(sig$gene, disease_genes)
    # sex-DEGs lost before the network, to the variance filter
    lost <- intersect(sig$gene, genes_dropped)
    cand[[sx]] <<- merge(sig[gene %in% hit, .(gene, logFC, adj.P.Val)],
                         gene_tab[, .(gene, module, GS_RA, kME_own)], by = "gene")
    # i-order (not setorder) since it must evaluate abs()
    cand[[sx]] <<- cand[[sx]][order(-abs(logFC))]
    fwrite(cand[[sx]], file.path(tab, sprintf("WGCNA_11_candidates_%s.csv", tolower(sx))))
    data.table(sex = sx, n_sigDEG = nrow(sig), n_disease_genes = length(disease_genes),
               n_candidates = length(hit),
               DEGs_lost_to_variance_filter = length(lost),
               pct_DEG_lost = round(100 * length(lost) / nrow(sig), 1)) }))
  fwrite(ov, file.path(tab, "WGCNA_12_candidate_summary.csv"))
  print(ov)
  say("")
  say("Candidates written: WGCNA_11_candidates_female.csv / _male.csv")
  say("Female n=%d | Male n=%d | shared=%d",
      nrow(cand$Female), nrow(cand$Male),
      length(intersect(cand$Female$gene, cand$Male$gene)))
}

# Step 15: explicit, named export
hdr("STEP 15  EXPORT")
saveRDS(list(
  config = CFG, cache_key = key_hash, soft_power = soft_power,
  power_auto = auto_power, power_faq = faq_power, sft_fit = fitI,
  datExpr = datExpr, meta = meta, traits = traits,
  moduleColors = moduleColors, MEs = MEs, kME = kME, GS_RA = GS_RA,
  gene_tab = gene_tab, moduleTraitCor = mtCor, moduleTraitP = mtP,
  module_trait_table = mt_tab,
  disease_modules = disease_modules, disease_genes = disease_genes,
  hub_table = hub_list, preservation = pres_dt,
  outliers_removed = outlier_ids, genes_dropped = genes_dropped),
  file.path(proc, "wgcna_results.rds"))
say("wrote data/processed/wgcna_results.rds")

# Figures (all inline, no separate figure script)
hdr("FIGURES")
pngf <- function(n, w = 1100, h = 800) png(file.path(fig, paste0(n, ".png")),
                                           width = w, height = h, res = 110)
gg <- function(p, n, w, h) { ggsave(file.path(fig, paste0(n, ".png")), p,
                                    width = w, height = h, dpi = 300)
                             ggsave(file.path(fig, paste0(n, ".pdf")), p,
                                    width = w, height = h) }

# F01 gene variance + filter cutoff
pngf("fig_wgcna_01_gene_variance"); hist(log10(gene_variance), breaks = 60,
  main = "Gene variance and filter cutoff", xlab = "log10(variance)",
  col = "lightblue", border = "white")
abline(v = log10(var_thresh), col = "red", lty = 2, lwd = 2)
legend("topright", sprintf("%.0fth pct cutoff", 100 * CFG$var_quantile),
       col = "red", lty = 2, lwd = 2, bty = "n"); dev.off()

# F02 sample clustering + outlier line
pngf("fig_wgcna_02_sample_clustering", 1400, 650); par(mar = c(2, 5, 3, 2))
plot(sample_tree, main = sprintf("Sample clustering (outliers removed: %d)", n_out),
     sub = "", xlab = "", cex = 0.6)
abline(h = cut_height, col = "red", lty = 2, lwd = 2); dev.off()

# F03 soft threshold
pngf("fig_wgcna_03_soft_threshold", 1300, 620); par(mfrow = c(1, 2))
plot(fitI[, 1], sfR2, type = "n", xlab = "Soft threshold (power)",
     ylab = "Scale-free topology fit, signed R^2", main = "Scale independence")
text(fitI[, 1], sfR2, labels = powers, col = "red", cex = 0.9)
abline(h = CFG$rsq_cut, col = "red", lty = 2)
abline(v = soft_power, col = "blue", lty = 3, lwd = 2)
plot(fitI[, 1], fitI[, 5], type = "n", xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", main = "Mean connectivity")
text(fitI[, 1], fitI[, 5], labels = powers, col = "red", cex = 0.9)
abline(v = soft_power, col = "blue", lty = 3, lwd = 2)
par(mfrow = c(1, 1)); dev.off()

# F04 dendrogram + module colours
pngf("fig_wgcna_04_dendrogram", 1300, 800)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
  "Module colours", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
  guideHang = 0.05, main = sprintf("Gene dendrogram, %d modules (power %d)",
                                   n_modules, soft_power)); dev.off()

# F05 module-trait heatmap (disease + sex), disease modules starred
Cm <- cbind(`Disease\n(RA vs HC)` = mtCor[, "RA"], `Sex\n(M vs F)` = mtCor[, "Male"])
Pm <- cbind(mtP[, "RA"], mtP[, "Male"])
oo <- order(-Cm[, 1]); Cm <- Cm[oo, , drop = FALSE]; Pm <- Pm[oo, , drop = FALSE]
txt <- paste0(signif(Cm, 2), "\n(", signif(Pm, 1), ")"); dim(txt) <- dim(Cm)
ylab <- sub("^ME", "", rownames(Cm))
ylab <- ifelse(ylab %in% disease_modules, paste0("* ", ylab), ylab)
pngf("fig_wgcna_05_module_trait", 700, 1000); par(mar = c(6, 10, 4, 4))
labeledHeatmap(Cm, xLabels = colnames(Cm), yLabels = rownames(Cm), ySymbols = ylab,
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt,
  setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1, 1),
  main = "Module-trait relationships\n(* = disease module)"); dev.off()

# F06 module eigengenes across samples
ac <- traits; bn <- setdiff(names(ac), "Age"); ac[bn] <- lapply(ac[bn], factor)
pngf("fig_wgcna_06_module_eigengenes", 1300, 1000)
print(pheatmap(t(MEs), scale = "row", annotation_col = ac, show_colnames = FALSE,
  clustering_method = "average", fontsize = 8,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Module eigengenes across samples")); dev.off()

# F07 ME correlation heatmap
mec <- stats::cor(MEs, use = "pairwise.complete.obs")
rownames(mec) <- colnames(mec) <- sub("^ME", "", colnames(MEs))
pngf("fig_wgcna_07_ME_correlation", 1100, 1000)
print(pheatmap(mec, color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-1, 1, length.out = 51), fontsize = 8,
  main = "Module eigengene correlation")); dev.off()

# F08 kME vs GS for each disease module
p8 <- ggplot(hub_list, aes(kME, GS_RA, colour = is_hub)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_vline(xintercept = c(-CFG$hub_kME, CFG$hub_kME), linetype = 2, colour = "grey50") +
  geom_hline(yintercept = c(-CFG$hub_GS, CFG$hub_GS), linetype = 2, colour = "grey50") +
  facet_wrap(~ module) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"), name = "Hub") +
  labs(title = "Module membership vs gene significance for RA",
       x = "kME (module membership)", y = "GS (correlation with RA)") +
  theme_bw(base_size = 11)
gg(p8, "fig_wgcna_08_kME_vs_GS", 9, 4.5)

# F09 top hub genes by connectivity
tops <- hub_list[, head(.SD[order(-connectivity)], 15), by = module]
tops[, gene_module := interaction(gene, module, drop = TRUE)]
p9 <- ggplot(tops, aes(reorder(gene_module, connectivity), connectivity, fill = module)) +
  geom_col(width = 0.6) + coord_flip() + facet_wrap(~ module, scales = "free") +
  scale_x_discrete(labels = function(x) sub("\\..*$", "", x)) +
  scale_fill_identity() +
  labs(x = NULL, y = "Intramodular connectivity") + theme_bw(base_size = 10)
gg(p9, "fig_wgcna_09_hub_genes", 10, 5)

# F10 disease-module expression heatmap (top 80 by |kME|)
sel <- hub_list[order(-abs(kME))][1:min(80, .N)]$gene
sa  <- data.frame(RA = factor(meta$group, c("HC", "RA")),
                  Sex = factor(meta$sex, c("F", "M")), row.names = rownames(datExpr))
pngf("fig_wgcna_10_disease_heatmap", 1300, 1000)
print(pheatmap(t(datExpr[, sel]), scale = "row", annotation_col = sa,
  show_rownames = FALSE, show_colnames = FALSE,
  clustering_distance_rows = "correlation",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_colors = list(RA = c(HC = "white", RA = "darkred"),
                           Sex = c(F = "darkgreen", M = "darkblue")),
  main = sprintf("Disease modules (%s): top %d genes by |kME|",
                 paste(disease_modules, collapse = "+"), length(sel)))); dev.off()

# F11 disease-module eigengene by group x sex
med <- melt(data.table(sample = rownames(MEs), group = meta$group, sex = meta$sex,
                       MEs[, paste0("ME", disease_modules), drop = FALSE]),
            id.vars = c("sample", "group", "sex"),
            variable.name = "module", value.name = "ME")
p11 <- ggplot(med, aes(group, ME, fill = group)) +
  geom_boxplot(outlier.size = .6, alpha = .85) + facet_grid(module ~ sex) +
  scale_fill_manual(values = c(HC = "#1565C0", RA = "#6A1B9A")) +
  labs(x = NULL, y = "Module eigengene") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
gg(p11, "fig_wgcna_11_ME_by_group_sex", 7, 6)

# F12 enrichment bubble
if (nrow(go_dt) || nrow(kg_dt)) {
  mk <- function(d, lab, n = 12) { if (!nrow(d)) return(NULL)
    d <- head(d, n)
    data.table(cat = lab, Description = substr(d$Description, 1, 46),
      GeneRatio = sapply(strsplit(d$GeneRatio, "/"),
                         function(x) as.numeric(x[1]) / as.numeric(x[2])),
      p.adjust = d$p.adjust, Count = d$Count) }
  pdf12 <- rbindlist(list(
    if (nrow(go_dt)) mk(go_dt[ONTOLOGY == "BP"], "GO: BP"),
    if (nrow(go_dt)) mk(go_dt[ONTOLOGY == "CC"], "GO: CC"),
    if (nrow(go_dt)) mk(go_dt[ONTOLOGY == "MF"], "GO: MF"),
    if (nrow(kg_dt)) mk(kg_dt, "KEGG")), use.names = TRUE)
  if (nrow(pdf12)) {
    p12 <- ggplot(pdf12, aes(GeneRatio, reorder(Description, GeneRatio))) +
      geom_point(aes(size = Count, colour = p.adjust)) +
      facet_wrap(~ cat, scales = "free", ncol = 2) +
      scale_colour_gradient(low = "#B2182B", high = "#2166AC", name = "adj. p") +
      labs(title = "Disease-module functional enrichment", x = "Gene ratio", y = NULL) +
      theme_bw(base_size = 9)
    gg(p12, "fig_wgcna_12_enrichment", 13, 9)
  }
}

# F13/F14 sex dendrograms
pngf("fig_wgcna_13_dendro_female", 1300, 600); par(mar = c(2, 5, 3, 2))
plotDendroAndColors(sn$female$dendrograms[[1]], colF[sn$female$blockGenes[[1]]],
  "Female modules", dendroLabels = FALSE, hang = .03, addGuide = TRUE,
  main = sprintf("Gene dendrogram - FEMALE (n=%d)", length(fs))); dev.off()
pngf("fig_wgcna_14_dendro_male", 1300, 600); par(mar = c(2, 5, 3, 2))
plotDendroAndColors(sn$male$dendrograms[[1]], colM[sn$male$blockGenes[[1]]],
  "Male modules", dendroLabels = FALSE, hang = .03, addGuide = TRUE,
  main = sprintf("Gene dendrogram - MALE (n=%d, underpowered)", length(ms))); dev.off()

# F15 module preservation
if (nrow(pres_dt)) {
  p15 <- ggplot(pres_dt, aes(moduleSize, Zsummary, colour = status, label = module)) +
    geom_hline(yintercept = c(2, 10), linetype = 2, colour = "grey50") +
    geom_point(size = 4, alpha = .85) + geom_text_repel(size = 3, max.overlaps = 20) +
    scale_colour_manual(values = c(`Not preserved` = "#d62728",
                                   Moderate = "#ff7f0e", Strong = "#2ca02c")) +
    scale_x_log10() +
    labs(x = "Module size (genes, log scale)", y = "Preservation Zsummary",
         caption = "Zsummary > 10 strong | 2-10 moderate | < 2 not preserved") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  gg(p15, "fig_wgcna_15_preservation", 10, 7)
}

# F16 the deliverable: candidates per sex
if (length(cand)) {
  cdf <- rbindlist(lapply(names(cand), function(s)
           data.table(sex = s, cand[[s]])), use.names = TRUE)
  p16 <- ggplot(cdf, aes(logFC, GS_RA, colour = module)) +
    geom_point(alpha = .7, size = 1.6) + facet_wrap(~ sex) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    labs(title = "Candidates: disease-module genes that are also sex-stratified DEGs",
         x = "log2 fold change (RA vs HC, within sex)",
         y = "Gene significance for RA (WGCNA)") +
    theme_bw(base_size = 11)
  gg(p16, "fig_wgcna_16_candidates", 10, 5)
}

