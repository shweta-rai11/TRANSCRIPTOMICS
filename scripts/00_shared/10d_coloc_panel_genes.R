#!/usr/bin/env Rscript
# =============================================================================
# 10d_coloc_panel_genes.R  —  Bayesian colocalisation for every MR-prioritised gene
#
# WHY THIS EXISTS — THE HOLE IT CLOSES
#   A cis-MR estimate answers "is the eQTL for gene G associated with RA?". It
#   does NOT answer "is the eQTL for gene G and the RA risk signal the SAME
#   underlying causal variant?". Those two questions come apart whenever the
#   eQTL variant is merely in LD with a distinct disease-causing variant, and in
#   that situation cis-MR returns a confident, highly significant, and entirely
#   spurious causal estimate. This is the single most-cited weakness of cis-eQTL
#   MR (Zhu 2016; Wallace 2020), and it is precisely why the MHC results in
#   10_MR.R cannot be taken at face value: the extended MHC has the longest-range
#   LD in the genome and contains the dominant RA locus, HLA-DRB1.
#
#   10c_MR_mhc_sensitivity.R answers the question by BLUNT EXCISION (drop the
#   MHC and see what survives). This script answers it DIRECTLY, for every gene
#   including the non-MHC ones, by testing the two competing hypotheses formally.
#
# WHAT coloc.abf TESTS
#   Over a cis window it evaluates five mutually exclusive hypotheses:
#     H0  no causal variant for either trait
#     H1  causal variant for expression only
#     H2  causal variant for RA only
#     H3  BOTH traits have a causal variant, but they are DIFFERENT variants
#     H4  BOTH traits share ONE causal variant  (true colocalisation)
#   H3 is the failure mode that invalidates a cis-MR estimate. H4 is the state of
#   the world the MR estimate silently assumes. Reporting PP.H4 alongside every
#   MR estimate therefore converts an untested assumption into a measured one.
#
# INTERPRETATION THRESHOLDS ADOPTED (conventional; Giambartolomei 2014)
#   PP.H4 >= 0.80                     strong support for a shared causal variant
#   0.50 <= PP.H4 < 0.80              suggestive
#   PP.H3 >= 0.80                     evidence AGAINST colocalisation: distinct
#                                     causal variants, MR estimate LD-confounded
#   PP.H4/(PP.H3+PP.H4) >= 0.80       conditional support given both traits are
#                                     associated (robust when power is low, i.e.
#                                     when PP.H0/H1/H2 absorb most of the mass)
#
# PRIOR SENSITIVITY
#   The p12 prior (probability a variant is causal for BOTH traits) is the main
#   subjective input and the usual point of attack. Every gene is therefore run
#   at the default p12 = 1e-5 AND at the conservative p12 = 1e-6, and both are
#   reported. A colocalisation that only appears at the permissive prior is not
#   reported as support.
#
# **THE SINGLE-CAUSAL-VARIANT ASSUMPTION, AND WHERE IT BREAKS**
#   coloc.abf assumes AT MOST ONE causal variant per trait in the region. That
#   assumption is reasonable for a typical cis-eQTL window. It is FALSE in the
#   MHC, where rheumatoid arthritis has multiple independent, well-mapped signals
#   (HLA-DRB1 positions 11/71/74, HLA-B position 9, HLA-DPB1 position 9;
#   Raychaudhuri 2012), and where the regions queried here contain 1,026-5,193
#   SNPs in extended LD.
#
#   When the assumption is violated, PP.H3 is inflated: two traits each driven by
#   several variants will look like "different causal variants" even when they
#   share one. **A high PP.H3 inside the MHC therefore CANNOT be read as clean
#   evidence of distinct causal variants.** The correct and weaker conclusion is:
#
#       In the MHC, cis-MR and coloc.abf are BOTH unreliable, so no causal claim
#       can be supported there in either direction.
#
#   That conclusion is sufficient for this project's purposes - it removes the
#   causal claim either way - but it must not be overstated into a positive
#   finding. Outside the MHC the assumption is defensible and PP.H3 there
#   (INPP5B 0.929, ESYT1 0.912, CDC37 1.000, NCOA5 1.000) does carry its usual
#   meaning. The MHC and non-MHC verdicts are therefore flagged separately in the
#   output via the `assumption_valid` column.
#
#   The proper fix for the MHC is coloc.susie, which permits multiple causal
#   variants per region. It requires an LD reference matrix for the region and is
#   recorded as outstanding work in results/RESULTS_ROBUSTNESS.md.
#
# DATA
#   Exposure : eQTLGen whole-blood cis-eQTL, full regional summary statistics
#              (NOT the clumped instruments - coloc requires all SNPs in the
#              window, which is why this cannot be done from the cached
#              instrument file and must query OpenGWAS).
#   Outcome  : Okada 2014 European RA GWAS, ieu-a-832
#              (14,361 cases / 43,923 controls; case proportion s = 0.246).
#   Window   : gene body +/- 250 kb, GRCh37 coordinates from EnsDb.Hsapiens.v75
#              (the same coordinate source 10_MR.R uses for its cis filter, so
#              the two analyses cannot disagree about where a gene is).
#
#   NETWORK. This script queries OpenGWAS and is therefore the one step here that
#   is not offline. Every regional extract is CACHED to
#   data/processed/new/coloc_regions.rds, so re-runs are deterministic and the
#   reported numbers can be regenerated without network access.
#
#   in : data/processed/new/MR_primary_objects.rds
#        results/tables/FS_input_{female,male}.csv
#        data/processed/new/ml_features.rds
#   out: results/tables/COLOC_results.csv          per gene, both priors
#        results/tables/COLOC_panel_genes.csv      panel genes, MR + MHC + coloc
#        results/tables/COLOC_summary.csv
#        data/processed/new/coloc_regions.rds      cached regional stats
#        data/processed/new/coloc_objects.rds
#
#   Giambartolomei C, et al. PLoS Genet 2014;10:e1004383.   (coloc.abf)
#   Wallace C. PLoS Genet 2020;16:e1008720.                 (priors, p12)
#   Zhu Z, et al. Nat Genet 2016;48:481-487.                (SMR/HEIDI rationale)
#   Vosa U, et al. Nat Genet 2021;53:1300-1310.             (eQTLGen)
#   Okada Y, et al. Nature 2014;506:376-381.                (RA GWAS)
# =============================================================================
suppressMessages({
  library(coloc); library(ieugwasr); library(data.table)
  library(EnsDb.Hsapiens.v75); library(ensembldb)
})
options(stringsAsFactors = FALSE)
set.seed(2024)

proc <- "data/processed/new"; tab <- "results/tables"
dir.create(proc, showWarnings = FALSE, recursive = TRUE)

CFG <- list(
  outcome    = "ieu-a-832",
  window     = 250e3,      # cis window each side of the gene body
  n_cases    = 14361,      # Okada 2014 European
  n_controls = 43923,
  p12_main   = 1e-5,       # coloc default
  p12_cons   = 1e-6,       # conservative sensitivity prior
  min_snps   = 50,         # below this the posterior is not interpretable
  pp4_strong = 0.80,
  pp3_strong = 0.80,
  refresh    = FALSE       # TRUE = ignore the region cache and re-query
)
CFG$s <- CFG$n_cases / (CFG$n_cases + CFG$n_controls)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

hdr("STEP 0  CONFIGURATION")
say("outcome GWAS  : %s (%d cases / %d controls, s = %.3f)",
    CFG$outcome, CFG$n_cases, CFG$n_controls, CFG$s)
say("cis window    : +/- %.0f kb around the gene body (GRCh37)", CFG$window / 1e3)
say("priors        : p12 = %g (main) and %g (conservative)", CFG$p12_main, CFG$p12_cons)

# =============================================================================
# STEP 1 — GENE LIST: every FDR-surviving prioritised gene, plus the panels
# =============================================================================
hdr("STEP 1  GENE LIST")
o       <- readRDS(file.path(proc, "MR_primary_objects.rds"))
primary <- as.data.table(o$primary)

fsF <- fread(file.path(tab, "FS_input_female.csv"))
fsM <- fread(file.path(tab, "FS_input_male.csv"))

panels <- list(female = character(0), male = character(0))
ml_path <- file.path(proc, "ml_features.rds")
if (file.exists(ml_path)) {
  ml <- readRDS(ml_path)
  panels$female <- ml$female$consensus
  panels$male   <- ml$male$consensus
}
panel_all <- union(panels$female, panels$male)

genes <- sort(unique(c(fsF$gene, fsM$gene, panel_all)))
say("FDR-surviving prioritised genes : female %d, male %d", nrow(fsF), nrow(fsM))
say("final panel genes          : female %d, male %d (%d unique)",
    length(panels$female), length(panels$male), length(panel_all))
say("genes to colocalise        : %d unique", length(genes))

# ENSG id per gene, taken from the exposure id the MR actually used, so the
# coloc window and the MR estimate always refer to the same eQTLGen dataset.
ens_map <- primary[gene %in% genes, .(gene, ensg = sub("^eqtl-a-", "", id.exposure))]
ens_map <- unique(ens_map, by = "gene")

# GRCh37 gene bodies (same source as the 10_MR.R cis filter)
gr <- suppressWarnings(
  ensembldb::genes(EnsDb.Hsapiens.v75,
                   filter = GeneIdFilter(ens_map$ensg),
                   columns = c("gene_id", "seq_name", "gene_seq_start", "gene_seq_end")))
coords <- data.table(ensg   = gr$gene_id,
                     g_chr  = as.character(GenomicRanges::seqnames(gr)),
                     g_start = GenomicRanges::start(gr),
                     g_end   = GenomicRanges::end(gr))
gl <- merge(ens_map, coords, by = "ensg", all.x = TRUE)
say("gene coordinates resolved  : %d of %d", sum(!is.na(gl$g_chr)), nrow(gl))
if (any(is.na(gl$g_chr)))
  say("  no GRCh37 coordinate (skipped): %s",
      paste(gl$gene[is.na(gl$g_chr)], collapse = ", "))
gl <- gl[!is.na(g_chr)]

# =============================================================================
# STEP 2 — REGIONAL SUMMARY STATISTICS (cached)
# =============================================================================
hdr("STEP 2  FETCH REGIONAL SUMMARY STATISTICS")
CACHE <- file.path(proc, "coloc_regions.rds")
regions <- if (file.exists(CACHE) && !CFG$refresh) readRDS(CACHE) else list()
say("cache: %d regions already stored", length(regions))

fetch_region <- function(g, ensg, chr, start, end) {
  rng <- sprintf("%s:%.0f-%.0f", chr, max(1, start - CFG$window), end + CFG$window)
  eq <- tryCatch(as.data.table(associations(variants = rng,
                                            id = paste0("eqtl-a-", ensg))),
                 error = function(e) NULL)
  ra <- tryCatch(as.data.table(associations(variants = rng, id = CFG$outcome)),
                 error = function(e) NULL)
  list(gene = g, ensg = ensg, range = rng, eqtl = eq, gwas = ra)
}

todo <- gl[!(gene %in% names(regions))]
if (nrow(todo)) {
  say("fetching %d regions from OpenGWAS ...", nrow(todo))
  for (i in seq_len(nrow(todo))) {
    r <- todo[i]
    reg <- suppressMessages(fetch_region(r$gene, r$ensg, r$g_chr, r$g_start, r$g_end))
    regions[[r$gene]] <- reg
    say("  [%3d/%3d] %-10s %-24s eQTL %5s SNPs | GWAS %5s SNPs",
        i, nrow(todo), r$gene, reg$range,
        if (is.null(reg$eqtl)) "ERR" else nrow(reg$eqtl),
        if (is.null(reg$gwas)) "ERR" else nrow(reg$gwas))
    if (i %% 10 == 0) saveRDS(regions, CACHE)   # checkpoint against a dropped connection
  }
  saveRDS(regions, CACHE)
} else {
  say("all regions served from cache - no network access required")
}

# =============================================================================
# STEP 3 — HARMONISE AND COLOCALISE
# -----------------------------------------------------------------------------
# Harmonisation rules, all conservative:
#   * merge on rsid; keep only SNPs present in BOTH datasets
#   * align the GWAS effect allele to the eQTL effect allele, flipping the sign
#     of beta where the alleles are swapped
#   * DROP strand-ambiguous SNPs (A/T, C/G). The Okada extract carries no allele
#     frequency, so strand cannot be resolved by MAF and a wrong flip would be
#     silent. Dropping them costs power and buys correctness.
#   * DROP SNPs whose allele pairs do not match after flipping (multi-allelic or
#     annotation mismatch)
#   * MAF is taken from the eQTL dataset (both studies are European)
# =============================================================================
hdr("STEP 3  HARMONISE AND RUN coloc.abf")

AMBIG <- c("AT", "TA", "CG", "GC")

harmonise_region <- function(reg) {
  eq <- reg$eqtl; ra <- reg$gwas
  if (is.null(eq) || is.null(ra) || !nrow(eq) || !nrow(ra)) return(NULL)
  eq <- eq[!is.na(beta) & !is.na(se) & se > 0 & !is.na(eaf) & eaf > 0 & eaf < 1]
  ra <- ra[!is.na(beta) & !is.na(se) & se > 0]
  eq <- unique(eq, by = "rsid"); ra <- unique(ra, by = "rsid")
  m <- merge(eq[, .(rsid, ea_e = toupper(ea), nea_e = toupper(nea),
                    b_e = beta, se_e = se, maf = pmin(eaf, 1 - eaf), n_e = n)],
             ra[, .(rsid, ea_g = toupper(ea), nea_g = toupper(nea),
                    b_g = beta, se_g = se)],
             by = "rsid")
  if (!nrow(m)) return(NULL)
  m <- m[!paste0(ea_e, nea_e) %in% AMBIG]                       # strand-ambiguous
  m[, same := ea_e == ea_g & nea_e == nea_g]
  m[, flip := ea_e == nea_g & nea_e == ea_g]
  m <- m[same | flip]                                            # resolvable only
  m[flip == TRUE, b_g := -b_g]
  m <- m[maf > 0 & maf < 1 & se_e > 0 & se_g > 0]
  unique(m, by = "rsid")
}

run_coloc <- function(m, p12) {
  d_eqtl <- list(snp = m$rsid, beta = m$b_e, varbeta = m$se_e^2,
                 MAF = m$maf, N = as.integer(stats::median(m$n_e, na.rm = TRUE)),
                 type = "quant")
  d_gwas <- list(snp = m$rsid, beta = m$b_g, varbeta = m$se_g^2,
                 MAF = m$maf, N = CFG$n_cases + CFG$n_controls,
                 type = "cc", s = CFG$s)
  r <- suppressMessages(coloc.abf(dataset1 = d_eqtl, dataset2 = d_gwas, p12 = p12))
  as.list(r$summary)
}

out <- list()
for (g in gl$gene) {
  reg <- regions[[g]]
  m <- harmonise_region(reg)
  nsnp <- if (is.null(m)) 0L else nrow(m)
  if (nsnp < CFG$min_snps) {
    out[[g]] <- data.table(gene = g, nsnp_coloc = nsnp,
                           PP0 = NA_real_, PP1 = NA_real_, PP2 = NA_real_,
                           PP3 = NA_real_, PP4 = NA_real_,
                           PP4_cons = NA_real_,
                           coloc_status = sprintf("NOT RUN (%d shared SNPs < %d)",
                                                  nsnp, CFG$min_snps))
    say("  %-10s SKIPPED - only %d usable shared SNPs", g, nsnp)
    next
  }
  s1 <- run_coloc(m, CFG$p12_main)
  s2 <- run_coloc(m, CFG$p12_cons)
  out[[g]] <- data.table(
    gene = g, nsnp_coloc = nsnp,
    PP0 = s1$PP.H0.abf, PP1 = s1$PP.H1.abf, PP2 = s1$PP.H2.abf,
    PP3 = s1$PP.H3.abf, PP4 = s1$PP.H4.abf,
    PP4_cons = s2$PP.H4.abf,
    coloc_status = NA_character_)
  say("  %-10s %5d SNPs | PP3 %.3f  PP4 %.3f  (PP4 at p12=1e-6: %.3f)",
      g, nsnp, s1$PP.H3.abf, s1$PP.H4.abf, s2$PP.H4.abf)
}
res <- rbindlist(out)

# conditional support: PP4 given that both traits are associated in the region.
# Robust when the region is underpowered and PP0/PP1/PP2 hold most of the mass.
res[, PP4_conditional := PP4 / (PP3 + PP4)]

res[, coloc_verdict := fifelse(
      !is.na(coloc_status), coloc_status,
 fifelse(PP4 >= CFG$pp4_strong & PP4_cons >= CFG$pp4_strong,
                                  "COLOCALISED (shared causal variant, both priors)",
 fifelse(PP4 >= CFG$pp4_strong,   "COLOCALISED at default prior only (not robust to p12)",
 fifelse(PP3 >= CFG$pp3_strong,   "DISTINCT VARIANTS (LD-confounded: MR estimate not causal)",
 fifelse(PP4 >= 0.5,              "SUGGESTIVE colocalisation",
 fifelse(PP4_conditional >= 0.8,  "UNDERPOWERED but conditionally consistent",
                                  "INCONCLUSIVE (no strong evidence either way)"))))))]

# =============================================================================
# STEP 4 — JOIN TO THE MR AND MHC RESULTS
# =============================================================================
hdr("STEP 4  JOIN TO MR AND MHC SENSITIVITY")
mr_cols <- primary[, .(gene, MHC_gene, instrument_chr, nSNP, method,
                       OR = round(OR, 3), MR_pval = signif(pval, 3))]
res <- merge(res, mr_cols, by = "gene", all.x = TRUE)

# Flag where coloc.abf's single-causal-variant assumption is tenable. Inside the
# MHC it is not (multiple independent RA signals + extended LD), so a high PP.H3
# there is inflated and must be reported as "unreliable in both directions"
# rather than as positive evidence of distinct causal variants.
res[, assumption_valid := fifelse(MHC_gene == TRUE,
      "NO - MHC: multiple independent RA signals violate the single-causal-variant assumption",
      "yes - single-causal-variant assumption tenable")]
res[MHC_gene == TRUE & grepl("^DISTINCT", coloc_verdict),
    coloc_verdict := "MHC - UNRELIABLE BOTH WAYS (PP.H3 inflated by assumption violation; no causal claim either direction)"]
say("MHC genes reclassified from DISTINCT VARIANTS to UNRELIABLE: %d",
    sum(res$MHC_gene == TRUE & grepl("^MHC - UNRELIABLE", res$coloc_verdict), na.rm = TRUE))

fate_path <- file.path(tab, "MR_MHC_sensitivity_panel_fate.csv")
if (file.exists(fate_path)) {
  fate <- fread(fate_path)
  vf <- unique(fate[, .(gene, mhc_verdict = verdict)], by = "gene")
  res <- merge(res, vf, by = "gene", all.x = TRUE)
} else {
  res[, mhc_verdict := NA_character_]
  say("MR_MHC_sensitivity_panel_fate.csv not found - run 10c first for the join.")
}

res[, in_panel_female := gene %in% panels$female]
res[, in_panel_male   := gene %in% panels$male]
res[, in_FS_female    := gene %in% fsF$gene]
res[, in_FS_male      := gene %in% fsM$gene]
setorder(res, -PP4)
fwrite(res, file.path(tab, "COLOC_results.csv"))

pan <- res[in_panel_female | in_panel_male]
setorder(pan, -PP4)
fwrite(pan, file.path(tab, "COLOC_panel_genes.csv"))

say("")
say("---- FINAL PANEL GENES ----")
print(pan[, .(gene, F = in_panel_female, M = in_panel_male, MHC = MHC_gene,
              nSNP_MR = nSNP, OR, nsnp_coloc, PP3 = round(PP3, 3),
              PP4 = round(PP4, 3), PP4_1e6 = round(PP4_cons, 3), coloc_verdict)])

# =============================================================================
# STEP 5 — SUMMARY
# =============================================================================
hdr("STEP 5  SUMMARY")
tally <- function(d) data.table(
  n_genes            = nrow(d),
  n_run              = sum(is.na(d$coloc_status)),
  colocalised        = sum(grepl("^COLOCALISED \\(", d$coloc_verdict)),
  coloc_prior_fragile= sum(grepl("default prior only", d$coloc_verdict)),
  suggestive         = sum(grepl("^SUGGESTIVE", d$coloc_verdict)),
  distinct_variants  = sum(grepl("^DISTINCT", d$coloc_verdict)),
  mhc_unreliable     = sum(grepl("^MHC - UNRELIABLE", d$coloc_verdict)),
  inconclusive       = sum(grepl("^INCONCLUSIVE|^UNDERPOWERED", d$coloc_verdict)),
  not_run            = sum(!is.na(d$coloc_status)))

summ <- rbindlist(list(
  cbind(set = "All prioritised genes",       tally(res)),
  cbind(set = "MHC genes",              tally(res[MHC_gene == TRUE])),
  cbind(set = "non-MHC genes",          tally(res[MHC_gene == FALSE])),
  cbind(set = "Final panels (F union M)", tally(pan))))
fwrite(summ, file.path(tab, "COLOC_summary.csv"))
print(summ)

saveRDS(list(results = res, panel = pan, summary = summ, config = CFG),
        file.path(proc, "coloc_objects.rds"))

say("")
say("READING RULE FOR THE THESIS")
say("  A causal claim for a gene requires BOTH:")
say("    (1) ROBUST in the MHC sensitivity analysis (10c), AND")
say("    (2) PP.H4 >= %.2f at both priors here.", CFG$pp4_strong)
say("  PP.H3 >= %.2f is positive evidence that the MR estimate is LD-confounded,",
    CFG$pp3_strong)
say("  not merely absence of evidence, and must be reported as such.")
say("  Genes that are INCONCLUSIVE are underpowered for coloc - Okada 2014 is a")
say("  2014-vintage imputation and the regional SNP density is the limiting")
say("  factor, not the method. Report them as untested, never as negative.")
say("")
say("Wrote COLOC_results.csv, COLOC_panel_genes.csv, COLOC_summary.csv,")
say("coloc_regions.rds (cache) and coloc_objects.rds")
cat("\nDONE\n")
