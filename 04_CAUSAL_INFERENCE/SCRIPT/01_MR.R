#!/usr/bin/env Rscript
# Two-sample Mendelian randomisation, end to end: cis-eQTL instruments (eQTLGen) -> Okada 2014 RA GWAS outcome, per-sex FDR, FS_input. Supersedes the old 10/10b/11 script split.
suppressMessages({
  library(TwoSampleMR); library(dplyr); library(data.table)
  library(org.Hs.eg.db); library(AnnotationDbi); library(writexl)
  library(EnsDb.Hsapiens.v75); library(ensembldb)   # GRCh37 gene coords for the cis filter (D8)
})
set.seed(2024)

# Step 0: configuration
CFG <- list(
  outcome        = "ieu-a-832",   # Okada 2014 EUR RA
  eqtl_p         = 5e-8,          # instrument selection threshold
  min_F          = 10,            # weak-instrument filter (NON-BINDING, see note)
  cis_window     = 1e6,           # D8: cis = same chr AND within 1Mb of the gene body
  mhc            = c(25e6, 34e6), # extended MHC on chr6 (GRCh37), flagged not dropped
  exclude_mhc    = FALSE,         # TRUE = drop MHC instruments outright
  chunk_size     = 250,           # SNPs per outcome query
  retries        = 3,
  fdr_cut        = 0.05,
  fresh_extract  = FALSE          # TRUE = ignore cache, re-extract every gene
)
# min_F is a non-binding verification: p<5e-8 already implies F>=29.7, so it removes nothing

proc <- "data/processed/new"; tab <- "results/tables"
dir.create(proc, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(proc, "MR_instruments.rds")

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

hdr("STEP 0  CONFIGURATION")
say("outcome GWAS      : %s", CFG$outcome)
say("instrument p      : %g   |  F >= %d", CFG$eqtl_p, CFG$min_F)
say("outcome chunking  : %d SNPs/query, proxies on, %d retries",
    CFG$chunk_size, CFG$retries)
say("fresh extraction  : %s", CFG$fresh_extract)

# Step 1: load candidate genes
hdr("STEP 1  CANDIDATE GENES")
fem <- fread(file.path(tab, "candidates_female_disease.csv"))$gene
mal <- fread(file.path(tab, "candidates_male_disease.csv"))$gene
allg <- union(fem, mal)
say("female %d | male %d | union %d (MR is run ONCE on the union)",
    length(fem), length(mal), length(allg))

# Step 2: instrument extraction (cis-eQTL exposure), cached and resumable per gene
hdr("STEP 2  INSTRUMENT EXTRACTION")
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = allg,
        keytype = "SYMBOL", columns = "ENSEMBL"))
setDT(map)

extract_one <- function(g) {
  ensgs <- unique(map[SYMBOL == g & !is.na(ENSEMBL), ENSEMBL])
  for (e in ensgs) {
    for (att in seq_len(CFG$retries)) {
      x <- tryCatch(extract_instruments(outcomes = paste0("eqtl-a-", e),
                                        p1 = CFG$eqtl_p, clump = TRUE),
                    error = function(err) NULL)
      if (!is.null(x) && nrow(x)) { x$gene <- g; return(x) }
      if (is.null(x)) Sys.sleep(1) else break   # NULL = error -> retry
    }
  }
  NULL
}

cache_obj <- if (!CFG$fresh_extract && file.exists(CACHE)) readRDS(CACHE) else NULL
cached_inst <- if (!is.null(cache_obj)) as.data.table(cache_obj$inst) else data.table()
cached_none <- if (!is.null(cache_obj)) cache_obj$no_instrument else character(0)

# D4 fix: never index an empty data.table by a column that does not exist.
have <- if (nrow(cached_inst)) intersect(allg, unique(cached_inst$gene)) else character(0)
skip <- intersect(allg, cached_none)
need <- setdiff(allg, union(have, skip))
say("cache: %d genes with instruments reused, %d known-empty skipped, %d to query",
    length(have), length(skip), length(need))

inst_list <- list(); found <- character(0); empty <- character(0)
for (i in seq_along(need)) {
  x <- extract_one(need[i])
  if (!is.null(x)) { inst_list[[need[i]]] <- x; found <- c(found, need[i]) }
  else empty <- c(empty, need[i])
  if (i %% 100 == 0)
    say("  ...%d/%d queried, %d with instruments", i, length(need), length(found))
}

# cache accumulates across all genes ever queried; only the analysis object is subset to the current candidate set
cache_all <- rbindlist(c(if (nrow(cached_inst)) list(cached_inst) else NULL,
                         inst_list), fill = TRUE)
if (nrow(cache_all)) cache_all <- unique(cache_all, by = c("gene", "SNP"))
saveRDS(list(inst = cache_all, no_instrument = unique(c(cached_none, empty))), CACHE)
say("cache now holds %d genes (%d SNPs) for reuse",
    uniqueN(cache_all$gene), nrow(cache_all))

inst <- copy(cache_all)
stopifnot(nrow(inst) > 0)
inst[, Fstat := (beta.exposure / se.exposure)^2]
inst <- inst[Fstat >= CFG$min_F][gene %in% allg]
say("instruments before cis filter: %d SNPs / %d genes", nrow(inst), uniqueN(inst$gene))

# cis filter: restrict instruments to same-chromosome, within cis_window of the gene body (GRCh37, via EnsDb.Hsapiens.v75), excluding trans associations
inst[, chr.exposure := as.character(chr.exposure)]
gr <- ensembldb::genes(EnsDb.Hsapiens.v75,
                       filter = AnnotationFilter::GeneNameFilter(unique(inst$gene)))
gcoord <- unique(data.table(gene    = gr$gene_name,
                            g_chr   = as.character(GenomeInfoDb::seqnames(gr)),
                            g_start = BiocGenerics::start(gr),
                            g_end   = BiocGenerics::end(gr)), by = "gene")
inst <- merge(inst, gcoord, by = "gene", all.x = TRUE)

n_nocoord <- uniqueN(inst[is.na(g_chr)]$gene)
inst[, cis := !is.na(g_chr) & chr.exposure == g_chr &
              pos.exposure >= (g_start - CFG$cis_window) &
              pos.exposure <= (g_end   + CFG$cis_window)]
inst[, MHC := chr.exposure == "6" &
              pos.exposure > CFG$mhc[1] & pos.exposure < CFG$mhc[2]]

say("  cis filter (+/-%.0f kb of gene body, GRCh37):", CFG$cis_window / 1e3)
say("    kept (cis)          : %d SNPs", sum(inst$cis))
say("    dropped: other chr  : %d SNPs", sum(!is.na(inst$g_chr) & inst$chr.exposure != inst$g_chr))
say("    dropped: same chr >1Mb : %d SNPs",
    sum(!is.na(inst$g_chr) & inst$chr.exposure == inst$g_chr & !inst$cis))
say("    genes with no GRCh37 coordinate (dropped): %d", n_nocoord)
inst <- inst[cis == TRUE]

if (isTRUE(CFG$exclude_mhc)) {
  n_mhc <- sum(inst$MHC)
  inst <- inst[MHC == FALSE]
  say("    MHC excluded (exclude_mhc=TRUE): %d SNPs removed", n_mhc)
} else {
  say("    MHC instruments RETAINED and flagged: %d SNPs in %d genes",
      sum(inst$MHC), uniqueN(inst[MHC == TRUE]$gene))
  say("    NOTE: MHC LD is extensive; single-SNP estimates there cannot be")
  say("          attributed to individual genes independently of HLA-DRB1.")
}

say("instruments after cis filter: %d SNPs / %d genes (min F = %.1f)",
    nrow(inst), uniqueN(inst$gene), min(inst$Fstat))
say("  %d of %d candidate genes have no usable cis instrument",
    length(allg) - uniqueN(inst$gene), length(allg))
stopifnot(nrow(inst) > 0)

# Step 3: outcome (RA GWAS) extraction, chunked with proxies and retries
hdr("STEP 3  OUTCOME EXTRACTION")
snps <- unique(inst$SNP)
chunks <- split(snps, ceiling(seq_along(snps) / CFG$chunk_size))
say("fetching %s for %d SNPs in %d chunks", CFG$outcome, length(snps), length(chunks))
out_list <- list()
for (i in seq_along(chunks)) {
  o <- NULL
  for (att in seq_len(CFG$retries)) {
    o <- tryCatch(extract_outcome_data(snps = chunks[[i]], outcomes = CFG$outcome,
                                       proxies = TRUE),
                  error = function(e) {
                    say("   chunk %d attempt %d failed: %s", i, att, conditionMessage(e))
                    NULL })
    if (!is.null(o)) break
  }
  if (!is.null(o)) out_list[[i]] <- o
  say("  chunk %d/%d done (%d SNPs)", i, length(chunks), if (is.null(o)) 0L else nrow(o))
}
out <- as.data.frame(rbindlist(out_list, fill = TRUE))
say("outcome SNPs retrieved: %d of %d", nrow(out), length(snps))

# Step 4: harmonise (action=2 infers strand for palindromic SNPs, drops ambiguous ones)
hdr("STEP 4  HARMONISATION")
dat <- harmonise_data(as.data.frame(inst), out, action = 2)

# match on (SNP, id.exposure) pair, not SNP alone, so a multi-gene SNP isn't mislabelled
dat$gene <- inst$gene[match(paste(dat$SNP, dat$id.exposure),
                            paste(inst$SNP, inst$id.exposure))]
stopifnot(!any(is.na(dat$gene)))
say("harmonised rows: %d | genes: %d | mr_keep: %d",
    nrow(dat), length(unique(dat$gene)), sum(dat$mr_keep))

# Step 5: MR per gene (STEP 6 sensitivity computed in the same pass)
hdr("STEP 5  MR ESTIMATION")
res <- list(); het <- list(); pleio <- list()
genes <- unique(dat$gene)
for (k in seq_along(genes)) {
  g <- genes[k]
  d <- dat[dat$gene == g & dat$mr_keep, , drop = FALSE]
  if (nrow(d) < 1) next
  ms <- if (nrow(d) >= 3) c("mr_ivw", "mr_egger_regression", "mr_weighted_median") else
        if (nrow(d) == 2) "mr_ivw" else "mr_wald_ratio"
  r <- tryCatch(mr(d, method_list = ms), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) next
  r$gene <- g; r$nSNP <- nrow(d); res[[g]] <- r
  if (nrow(d) >= 3) {
    het[[g]]   <- tryCatch(cbind(gene = g, mr_heterogeneity(d)), error = function(e) NULL)
    pleio[[g]] <- tryCatch(cbind(gene = g, mr_pleiotropy_test(d)), error = function(e) NULL)
  }
}
res <- bind_rows(res)
stopifnot(nrow(res) > 0)
res_or <- as.data.table(generate_odds_ratios(res))
say("MR completed for %d genes", uniqueN(res$gene))

# every dplyr verb namespace-qualified: IRanges::slice masks dplyr::slice once org.Hs.eg.db is attached
primary <- res %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(factor(method, levels = c("Inverse variance weighted",
                                           "Wald ratio",
                                           "Weighted median",
                                           "MR Egger")), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  as.data.table()
primary[, OR := exp(b)][, OR_lo := exp(b - 1.96 * se)][, OR_hi := exp(b + 1.96 * se)]
primary[, FDR_pooled := p.adjust(pval, "BH")]     # legacy, NOT reported (D5)
primary[, risk := ifelse(pval >= 0.05, "ns",
                  ifelse(OR > 1, "risk (OR>1)", "protective (OR<1)"))]

# carry the MHC flag and instrument count through to every result table
gene_flags <- inst[, .(MHC_gene = any(MHC), instrument_chr = chr.exposure[1]), by = gene]
primary <- merge(primary, gene_flags, by = "gene", all.x = TRUE)
primary[, sensitivity_testable := nSNP >= 3]
say("primary table flagged: %d genes in MHC, %d genes with >=3 instruments",
    sum(primary$MHC_gene, na.rm = TRUE), sum(primary$sensitivity_testable))

het <- bind_rows(het); pleio <- bind_rows(pleio)

hdr("STEP 6  SENSITIVITY")
say("heterogeneity tests : %d genes (>=3 instruments)",
    if (nrow(het)) uniqueN(het$gene) else 0L)
say("pleiotropy tests    : %d genes (>=3 instruments)",
    if (nrow(pleio)) uniqueN(pleio$gene) else 0L)
say("NOTE: genes with 1-2 instruments carry NO pleiotropy assessment; for those")
say("      the exclusion-restriction assumption is untested, not satisfied.")

saveRDS(list(primary = primary, res_or = res_or, het = het, pleio = pleio,
             inst = inst, dat = as.data.frame(dat)),
        file.path(proc, "MR_primary_objects.rds"))

# Step 7: per-sex output tables, within-stratum FDR, FS_input = FDR-surviving set
hdr("STEP 7  PER-SEX TABLES")
emit <- function(sx, sex_genes) {
  pr  <- copy(primary[gene %in% sex_genes])
  ro  <- res_or[gene %in% sex_genes]
  ins <- inst[gene %in% sex_genes]
  ht  <- if (nrow(het))   as.data.table(het)[gene %in% sex_genes]   else data.table()
  pl  <- if (nrow(pleio)) as.data.table(pleio)[gene %in% sex_genes] else data.table()

  # correction denominator = genes tested in this stratum
  pr[, FDR_stratum := p.adjust(pval, "BH")]
  pr[, risk_fdr := ifelse(FDR_stratum >= CFG$fdr_cut, "ns",
                   ifelse(OR > 1, "risk (OR>1)", "protective (OR<1)"))]

  fwrite(ins, file.path(tab, sprintf("MR_%s_TABLE1_instruments.csv", sx)))
  fwrite(ro,  file.path(tab, sprintf("MR_%s_TABLE2_results_allmethods.csv", sx)))
  fwrite(pl,  file.path(tab, sprintf("MR_%s_TABLE3_pleiotropy.csv", sx)))
  fwrite(ht,  file.path(tab, sprintf("MR_%s_TABLE4_heterogeneity.csv", sx)))
  fwrite(pr[order(pval)], file.path(tab, sprintf("MR_%s_primary_okada.csv", sx)))

  sheets <- list(`1_instruments`   = as.data.frame(ins),
                 `2_MR_allmethods` = as.data.frame(ro),
                 `2b_MR_primary`   = as.data.frame(pr),
                 `3_pleiotropy`    = if (nrow(pl)) as.data.frame(pl),
                 `4_heterogeneity` = if (nrow(ht)) as.data.frame(ht))
  sheets <- sheets[!vapply(sheets, is.null, logical(1))]
  write_xlsx(sheets, file.path(tab, sprintf("MR_%s_all_tables.xlsx", sx)))

  # FS_input is the FDR-surviving set, not the nominal screen
  fs <- pr[risk_fdr != "ns", .(gene, direction = risk_fdr,
                               MR_OR = round(OR, 3),
                               MR_pval = signif(pval, 3),
                               MR_FDR_stratum = signif(FDR_stratum, 3),
                               nSNP, method)]
  setorder(fs, MR_pval)
  fwrite(fs, file.path(tab, sprintf("FS_input_%s.csv", sx)))
  fwrite(pr[FDR_stratum < CFG$fdr_cut][order(pval)],
         file.path(tab, sprintf("MR_causal_FDR_%s.csv", sx)))

  n_nom <- nrow(pr[risk != "ns"])
  say("== %s ==", toupper(sx))
  say("  candidates %d | MR-tested %d | nominal p<0.05 %d | FDR<%.2f %d",
      length(sex_genes), nrow(pr), n_nom, CFG$fdr_cut, nrow(fs))
  say("  expected false positives had the nominal set been used: %.1f of %d",
      0.05 * nrow(pr), n_nom)
  say("  FS_input (%d genes): %s", nrow(fs), paste(fs$gene, collapse = ", "))
  invisible(fs)
}
fs_f <- emit("female", fem)
fs_m <- emit("male",   mal)

hdr("SUMMARY")
# "list only" not "specific": a set difference between two FDR-thresholded lists is not a sex-specificity test
say("shared prioritised genes : %d", length(intersect(fs_f$gene, fs_m$gene)))
say("in female list only : %s", paste(setdiff(fs_f$gene, fs_m$gene), collapse = ", "))
say("in male list only   : %s", paste(setdiff(fs_m$gene, fs_f$gene), collapse = ", "))
say("")
say("REMINDER: shared genes carry IDENTICAL estimates in both files. Both GWAS")
say("inputs are sex-combined, so the MR contributes no sex-specificity; that")
say("comes solely from the upstream within-sex differential expression.")
cat("\nDONE\n")
