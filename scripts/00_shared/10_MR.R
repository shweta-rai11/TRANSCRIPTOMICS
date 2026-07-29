#!/usr/bin/env Rscript
# =============================================================================
# 10_MR.R  —  Two-sample Mendelian randomisation, END TO END in ONE script.
#
# SUPERSEDES scripts 10_mr_extract_instruments.R, 10b_mr_freshextract_primary.R
# and 11_mr_primary_and_FSinput.R. Those three were one workflow that had been
# broken into pieces by two successive failures:
#
#   10   extracted instruments (reusing a cache from ANOTHER project directory)
#        then fetched the whole outcome in ONE query -> that query timed out at
#        2,773 SNPs, so the script never reached the MR.
#   11   was written as a RESUME for 10: it skips extraction entirely, assumes
#        the instruments are already cached, and replaces the single outcome
#        query with 250-SNP chunks + LD proxies + retries. It cannot bootstrap
#        from nothing.
#   10b  was a clean rewrite doing both halves with no cache dependency, and
#        added the within-stratum FDR. It crashed on an unqualified slice().
#
# THIS SCRIPT is all of it, in order, with every one of those defects fixed.
#
# ---- DEFECTS CARRIED BY THE OLD SCRIPTS, AND THE FIX HERE -------------------
#  D1  NAMESPACE COLLISION (fatal in 10 and 10b).
#      org.Hs.eg.db -> AnnotationDbi -> IRanges is attached AFTER dplyr and
#      masks dplyr::slice. Unqualified slice() dispatched to IRanges::slice and
#      died with "Rle of type 'list' is not supported" AFTER all the MR had
#      been computed. Every dplyr verb here is namespace-qualified.
#  D2  CROSS-PROJECT CACHE (10). Script 10 read MR29_candidate_instruments.rds
#      from a different project tree, mixing provenance. The cache here is
#      written by THIS script, into THIS project, keyed by gene.
#  D3  UNCHUNKED OUTCOME FETCH (10). Replaced by 250-SNP chunks with proxies
#      and 3 retries.
#  D4  EMPTY-CACHE CRASH. `cached[gene %in% have]` errors on an empty
#      data.table ("Object 'gene' not found amongst []"). Guarded here.
#  D5  POOLED FDR (10 and 11). Both corrected across the union of both sexes.
#      In a sex-stratified design the denominator must be the genes tested IN
#      THAT STRATUM. Computed per stratum here; the pooled value is retained
#      only as a legacy column.
#  D6  FS_input BUILT FROM THE NOMINAL SCREEN. p<0.05 returned 113 female /
#      115 male genes against ~74.5 / ~74.7 expected false positives - about
#      two thirds noise - and that list was feeding feature selection.
#      FS_input is the FDR-surviving set.
#
# ---- DESIGN ----------------------------------------------------------------
#   Exposure : eQTLGen whole-blood eQTL, eqtl-a-<ENSG>, p<5e-8, LD-clumped
#              (r2<0.001, 10,000kb - TwoSampleMR defaults), then RESTRICTED TO CIS
#              (same chromosome, within 1Mb of the gene body, GRCh37) - see D8.
#              F = (beta/se)^2 >= 10 is a verification, not a filter: it is
#              implied by p<5e-8 (min F = 29.7) and removes nothing.
#   Outcome  : Okada 2014 European RA GWAS, ieu-a-832.
#   Estimator: >=3 SNPs -> IVW + MR-Egger + weighted median
#              2 SNPs   -> IVW
#              1 SNP    -> Wald ratio
#              primary ordered IVW > Wald > weighted median > MR-Egger.
#   Sensitivity: Cochran's Q and MR-Egger intercept where >=3 instruments.
#
#   NOTE ON SEX. Both eQTLGen and Okada are SEX-COMBINED, and a survey of all
#   37 RA datasets in OpenGWAS confirms no sex-stratified RA GWAS exists. The
#   MR estimate for a gene is therefore IDENTICAL in the female and male
#   outputs; the two differ only in which genes appear and in the FDR
#   denominator. This is sex-stratified DISCOVERY followed by sex-combined
#   CAUSAL VALIDATION. No sex-specific causal claim is supported.
#
#   in : results/tables/candidates_{female,male}_disease.csv
#   out: data/processed/new/MR_instruments.rds
#        data/processed/new/MR_primary_objects.rds
#        results/tables/FS_input_{female,male}.csv          (FDR-surviving)
#        results/tables/MR_causal_FDR_{female,male}.csv
#        results/tables/MR_{sex}_TABLE1-4.csv
#        results/tables/MR_{sex}_primary_okada.csv
#        results/tables/MR_{sex}_all_tables.xlsx
#
# ---- References -------------------------------------------------------------
#   Vosa U, et al. Nat Genet 2021;53:1300-1310.           (eQTLGen)
#   Okada Y, et al. Nature 2014;506:376-381.              (RA GWAS)
#   Burgess S, et al. Genet Epidemiol 2013;37:658-665.    (IVW)
#   Bowden J, et al. Int J Epidemiol 2015;44:512-525.     (MR-Egger)
#   Bowden J, et al. Genet Epidemiol 2016;40:304-314.     (weighted median)
#   Benjamini Y, Hochberg Y. J R Stat Soc B 1995;57:289-300.
# =============================================================================
suppressMessages({
  library(TwoSampleMR); library(dplyr); library(data.table)
  library(org.Hs.eg.db); library(AnnotationDbi); library(writexl)
  library(EnsDb.Hsapiens.v75); library(ensembldb)   # GRCh37 gene coords for the cis filter (D8)
})
set.seed(2024)

# =============================================================================
# STEP 0 — CONFIGURATION
# =============================================================================
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
# NOTE ON min_F. F = (beta/se)^2 = Z^2. At p < 5e-8, |Z| >= 5.4513, so F >= 29.72
# for EVERY instrument by construction. This filter removed 0 SNPs and cannot
# remove any. It is retained as a VERIFICATION that instrument strength exceeds
# the conventional threshold, not as a selection step. Report it as such.

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

# =============================================================================
# STEP 1 — LOAD CANDIDATE GENES
# =============================================================================
hdr("STEP 1  CANDIDATE GENES")
fem <- fread(file.path(tab, "candidates_female_disease.csv"))$gene
mal <- fread(file.path(tab, "candidates_male_disease.csv"))$gene
allg <- union(fem, mal)
say("female %d | male %d | union %d (MR is run ONCE on the union)",
    length(fem), length(mal), length(allg))

# =============================================================================
# STEP 2 — INSTRUMENTS (cis-eQTL exposure), with resume
#   Instruments are a property of the GENE, not of the candidate list, so a
#   cache written by this script is a legitimate reuse (cf. D2). Genes that
#   return nothing are recorded too, so they are not re-queried every run.
# =============================================================================
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

# D7 fix (2026-07-28). The cache must ACCUMULATE, never shrink. Previously the
# saved object was filtered to `allg` first, so running the script with a
# smaller candidate list silently discarded instruments for every gene outside
# it - a later run with the original list then had to re-query them (observed:
# the cache fell from 1,980 genes to 1,400). Instruments are a property of the
# GENE, so anything ever extracted is kept; only the ANALYSIS object is
# subset to the current candidate set.
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

# ---------------------------------------------------------------------------
# D8  CIS FILTER (added 2026-07-28).
#   The design has always claimed "cis-eQTL" instruments, but NOTHING enforced
#   it: extract_instruments() returns whatever OpenGWAS holds for eqtl-a-<ENSG>,
#   and eQTLGen's datasets include TRANS associations. 996 of 4,932 instruments
#   (20%) were on a DIFFERENT CHROMOSOME from the gene they instrumented.
#
#   The two largest effect estimates in the whole analysis came from trans
#   instruments sitting on the two strongest RA loci in the genome:
#     HNRNPM (chr19) instrumented by chr6:32,431,962  -> MHC class II, OR 287.3
#     FOXP3  (chrX)  instrumented by chr1:114,303,808 -> PTPN22 locus,  OR 263.9
#   Both are horizontal pleiotropy by construction: the variant reaches RA
#   through HLA-DRB1 / PTPN22, not through the gene's expression. The exclusion
#   restriction is VIOLATED, not merely untested.
#
#   Gene coordinates come from EnsDb.Hsapiens.v75 (Ensembl 75 = GRCh37), which
#   MATCHES the build of the eQTLGen/OpenGWAS SNP positions. Do not substitute a
#   GRCh38 annotation here - the coordinates would not be comparable.
# ---------------------------------------------------------------------------
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

# =============================================================================
# STEP 3 — OUTCOME (RA GWAS), chunked with proxies and retries  (D3)
# =============================================================================
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

# =============================================================================
# STEP 4 — HARMONISE
#   action = 2 infers strand for palindromic SNPs from allele frequency and
#   drops those that remain ambiguous.
# =============================================================================
hdr("STEP 4  HARMONISATION")
dat <- harmonise_data(as.data.frame(inst), out, action = 2)

# D9 fix (2026-07-28). Gene labels were previously re-attached with
#   dat$gene <- inst$gene[match(dat$SNP, inst$SNP)]
# match() returns the FIRST hit, so any SNP instrumenting more than one gene had
# ALL of its harmonised rows relabelled with whichever gene appeared first in
# `inst` - silently moving instruments between genes. Before the cis filter this
# affected 163 SNPs (one instrumented 83 genes); after it, 24 SNPs / 50 rows.
# The exposure dataset id (eqtl-a-<ENSG>) is what actually identifies the gene,
# so match on the (SNP, id.exposure) PAIR.
dat$gene <- inst$gene[match(paste(dat$SNP, dat$id.exposure),
                            paste(inst$SNP, inst$id.exposure))]
stopifnot(!any(is.na(dat$gene)))
say("harmonised rows: %d | genes: %d | mr_keep: %d",
    nrow(dat), length(unique(dat$gene)), sum(dat$mr_keep))

# =============================================================================
# STEP 5 — MR PER GENE  (+ STEP 6 sensitivity, computed in the same pass)
# =============================================================================
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

# D1 fix: EVERY dplyr verb namespace-qualified. IRanges::slice masks
# dplyr::slice once org.Hs.eg.db is attached, and the failure is silent until
# it kills the run after all the MR is done.
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

# D8: carry the MHC flag and the instrument count through to every result table,
# so no estimate is ever read without knowing (a) whether it sits in the MHC and
# (b) whether any pleiotropy assessment was possible for it at all.
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

# =============================================================================
# STEP 7 — PER-SEX OUTPUT, within-stratum FDR (D5), FS_input = FDR set (D6)
# =============================================================================
hdr("STEP 7  PER-SEX TABLES")
emit <- function(sx, sex_genes) {
  pr  <- copy(primary[gene %in% sex_genes])
  ro  <- res_or[gene %in% sex_genes]
  ins <- inst[gene %in% sex_genes]
  ht  <- if (nrow(het))   as.data.table(het)[gene %in% sex_genes]   else data.table()
  pl  <- if (nrow(pleio)) as.data.table(pleio)[gene %in% sex_genes] else data.table()

  # D5: correction denominator = genes tested IN THIS STRATUM
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

  # D6: FS_input is the FDR-SURVIVING set, not the nominal screen.
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
# Labels deliberately say "list only", NOT "specific". A set difference of two
# separately-FDR-thresholded lists is not a test of sex-specificity (same reason
# the setdiff sex-specific DEG lists were removed from 05_dge.R). Here it is
# weaker still: shared genes carry IDENTICAL estimates because both GWAS are
# sex-combined, so these differences reflect only which genes were eligible
# upstream, plus near-identical FDR denominators (1,477 vs 1,478 tested).
say("shared prioritised genes : %d", length(intersect(fs_f$gene, fs_m$gene)))
say("in female list only : %s", paste(setdiff(fs_f$gene, fs_m$gene), collapse = ", "))
say("in male list only   : %s", paste(setdiff(fs_m$gene, fs_f$gene), collapse = ", "))
say("")
say("REMINDER: shared genes carry IDENTICAL estimates in both files. Both GWAS")
say("inputs are sex-combined, so the MR contributes no sex-specificity; that")
say("comes solely from the upstream within-sex differential expression.")
cat("\nDONE\n")
