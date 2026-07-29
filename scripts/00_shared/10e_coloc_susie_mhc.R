#!/usr/bin/env Rscript
# =============================================================================
# 10e_coloc_susie_mhc.R  —  Multiple-causal-variant colocalisation for the MHC
#
# WHY THIS EXISTS
#   10d runs coloc.abf, which assumes AT MOST ONE causal variant per trait per
#   region. That assumption is defensible for a typical cis-eQTL window and is
#   FALSE in the MHC, where rheumatoid arthritis has several independent,
#   well-mapped signals (HLA-DRB1 positions 11/71/74, HLA-B position 9, HLA-DPB1
#   position 9; Raychaudhuri 2012) and where the regions queried carry 1,026-5,193
#   SNPs in extended LD.
#
#   The consequence was recorded honestly in 10d: inside the MHC a high PP.H3 is
#   INFLATED by the assumption violation, so it cannot be read as evidence of
#   distinct causal variants, and those genes were labelled "UNRELIABLE BOTH
#   WAYS". That is an accurate statement of ignorance, but it is still ignorance.
#   This script removes it.
#
# WHAT coloc.susie DOES DIFFERENTLY
#   SuSiE (Sum of Single Effects; Wang et al. 2020) decomposes each trait's
#   regional association into MULTIPLE credible sets, each representing one
#   independent causal signal. coloc.susie then tests colocalisation between
#   every pair of credible sets - eQTL signal i against GWAS signal j - and
#   returns a posterior per pair. A gene can therefore colocalise with ONE of
#   several RA signals while being independent of the others, which is exactly
#   the situation coloc.abf cannot represent and exactly the situation the MHC
#   presents.
#
# WHAT THIS NEEDS THAT coloc.abf DID NOT
#   An LD matrix for the region, in the SAME allele orientation as the summary
#   statistics. This is fetched from the OpenGWAS 1000G EUR reference panel via
#   ieugwasr::ld_matrix(). Two consequences must be stated rather than buried:
#
#     * REFERENCE-PANEL MISMATCH. The LD matrix comes from 1000G EUR (n ~ 500),
#       not from the eQTLGen or Okada samples. In the MHC, where haplotype
#       structure is extreme, an out-of-sample LD reference is a genuine source
#       of error and can make SuSiE report spurious or missing credible sets.
#       This is a real limitation of doing susie from summary statistics without
#       in-sample LD, and it is why the result below is reported as INDICATIVE.
#     * SNP CAP. ld_matrix is capped server-side, so each region is reduced to
#       the CFG$max_snps most significant shared SNPs (by the smaller of the two
#       traits' p-values). Truncating a region can drop a causal variant.
#
#   Both mean this analysis IMPROVES on coloc.abf in the MHC without settling it
#   definitively. In-sample LD would be required for that, and neither eQTLGen
#   nor Okada release it.
#
# HOW TO READ THE OUTPUT
#   Per gene, the best-supported eQTL-signal x GWAS-signal pair is reported with
#   its PP.H4. The verdict vocabulary matches 10d so the two can sit in one table:
#     PP.H4 >= 0.80 for some pair    COLOCALISED (with that RA signal)
#     PP.H3 >= 0.80 for every pair   DISTINCT VARIANTS (now a valid reading,
#                                    because multiple signals were modelled)
#     SuSiE found < 1 credible set   NOT RESOLVABLE (region underpowered or LD
#                                    reference inadequate) - NOT a negative result
#
#   in : data/processed/new/coloc_regions.rds   (cached by 10d; no re-download)
#        results/tables/COLOC_results.csv
#   out: results/tables/COLOC_SUSIE_mhc.csv
#        results/tables/COLOC_combined_abf_susie.csv
#        data/processed/new/coloc_susie_objects.rds
#
#   Wang G, et al. J R Stat Soc B 2020;82:1273-1300.        (SuSiE)
#   Wallace C. PLoS Genet 2021;17:e1009440.                 (coloc.susie)
#   Zou Y, et al. PLoS Genet 2022;18:e1010299.              (LD mismatch in susie)
#   Raychaudhuri S, et al. Nat Genet 2012;44:291-296.       (MHC fine-mapping, RA)
# =============================================================================
suppressMessages({
  library(coloc); library(susieR); library(ieugwasr); library(data.table)
})
options(stringsAsFactors = FALSE)
set.seed(2024)

proc <- "data/processed/new"; tab <- "results/tables"

CFG <- list(
  max_snps   = 450,     # server-side ld_matrix cap, with headroom
  min_snps   = 60,      # below this SuSiE is not worth running
  L          = 10,      # maximum credible sets per trait
  p12        = 1e-5,
  n_cases    = 14361,
  n_controls = 43923,
  coverage   = 0.95,
  pp4_strong = 0.80
)
CFG$s <- CFG$n_cases / (CFG$n_cases + CFG$n_controls)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# =============================================================================
# STEP 1 — TARGET GENES: everything 10d could not resolve
# =============================================================================
hdr("STEP 1  TARGETS")
regions <- readRDS(file.path(proc, "coloc_regions.rds"))
abf <- fread(file.path(tab, "COLOC_results.csv"))

targets <- abf[MHC_gene == TRUE]$gene
say("MHC genes unresolved by coloc.abf : %d", length(targets))
say("  %s", paste(targets, collapse = ", "))
panel_genes <- abf[in_panel_female == TRUE | in_panel_male == TRUE]$gene
say("of which in a final panel          : %s",
    paste(intersect(targets, panel_genes), collapse = ", "))

# =============================================================================
# STEP 2 — HARMONISE, FETCH LD, RUN SuSiE ON EACH TRAIT
# =============================================================================
hdr("STEP 2  SuSiE PER TRAIT, THEN coloc.susie")

AMBIG <- c("AT", "TA", "CG", "GC")

harmonise <- function(reg) {
  eq <- reg$eqtl; ra <- reg$gwas
  if (is.null(eq) || is.null(ra) || !nrow(eq) || !nrow(ra)) return(NULL)
  eq <- unique(as.data.table(eq)[!is.na(beta) & !is.na(se) & se > 0 &
                                 !is.na(eaf) & eaf > 0 & eaf < 1], by = "rsid")
  ra <- unique(as.data.table(ra)[!is.na(beta) & !is.na(se) & se > 0], by = "rsid")
  m <- merge(eq[, .(rsid, ea_e = toupper(ea), nea_e = toupper(nea),
                    b_e = beta, se_e = se, p_e = p,
                    maf = pmin(eaf, 1 - eaf), n_e = n)],
             ra[, .(rsid, ea_g = toupper(ea), nea_g = toupper(nea),
                    b_g = beta, se_g = se, p_g = p)], by = "rsid")
  if (!nrow(m)) return(NULL)
  m <- m[!paste0(ea_e, nea_e) %in% AMBIG]
  m[, same := ea_e == ea_g & nea_e == nea_g]
  m[, flip := ea_e == nea_g & nea_e == ea_g]
  m <- m[same | flip]
  m[flip == TRUE, b_g := -b_g]
  m[, best_p := pmin(p_e, p_g)]
  setorder(m, best_p)
  head(m, CFG$max_snps)
}

run_gene <- function(g) {
  reg <- regions[[g]]
  if (is.null(reg)) return(data.table(gene = g, status = "no cached region"))
  m <- harmonise(reg)
  if (is.null(m) || nrow(m) < CFG$min_snps)
    return(data.table(gene = g, n_snps = if (is.null(m)) 0L else nrow(m),
                      status = "too few shared SNPs"))

  ld <- tryCatch(suppressMessages(ld_matrix(m$rsid, with_alleles = TRUE, pop = "EUR")),
                 error = function(e) NULL)
  if (is.null(ld) || !nrow(ld))
    return(data.table(gene = g, n_snps = nrow(m), status = "LD matrix unavailable"))

  # ld_matrix returns rownames "rsid_A1_A2"; align orientation to our effect allele
  parts <- do.call(rbind, strsplit(rownames(ld), "_", fixed = TRUE))
  ldsnp <- parts[, 1]; ld_a1 <- toupper(parts[, 2]); ld_a2 <- toupper(parts[, 3])
  keep <- match(ldsnp, m$rsid)
  ok <- !is.na(keep)
  ld <- ld[ok, ok, drop = FALSE]; ldsnp <- ldsnp[ok]
  ld_a1 <- ld_a1[ok]; ld_a2 <- ld_a2[ok]
  mm <- m[match(ldsnp, rsid)]

  # flip the sign of any SNP whose LD-panel effect allele is our other allele
  needs_flip <- mm$ea_e == ld_a2 & mm$nea_e == ld_a1
  drop <- !((mm$ea_e == ld_a1 & mm$nea_e == ld_a2) | needs_flip)
  if (any(drop)) {
    ld <- ld[!drop, !drop, drop = FALSE]; mm <- mm[!drop]; needs_flip <- needs_flip[!drop]
  }
  if (nrow(mm) < CFG$min_snps)
    return(data.table(gene = g, n_snps = nrow(mm), status = "too few after LD alignment"))
  sgn <- ifelse(needs_flip, -1, 1)
  ld <- ld * outer(sgn, sgn)                       # re-orient LD to our alleles
  rownames(ld) <- colnames(ld) <- mm$rsid

  D_eqtl <- list(snp = mm$rsid, beta = mm$b_e, varbeta = mm$se_e^2, MAF = mm$maf,
                 N = as.integer(stats::median(mm$n_e, na.rm = TRUE)),
                 type = "quant", LD = ld, position = seq_len(nrow(mm)))
  D_gwas <- list(snp = mm$rsid, beta = mm$b_g, varbeta = mm$se_g^2, MAF = mm$maf,
                 N = CFG$n_cases + CFG$n_controls, type = "cc", s = CFG$s,
                 LD = ld, position = seq_len(nrow(mm)))

  s1 <- tryCatch(suppressWarnings(suppressMessages(
          runsusie(D_eqtl, L = CFG$L, coverage = CFG$coverage, repeat_until_convergence = FALSE))),
        error = function(e) NULL)
  s2 <- tryCatch(suppressWarnings(suppressMessages(
          runsusie(D_gwas, L = CFG$L, coverage = CFG$coverage, repeat_until_convergence = FALSE))),
        error = function(e) NULL)
  n_cs <- function(s) if (is.null(s) || is.null(s$sets$cs)) 0L else length(s$sets$cs)
  if (n_cs(s1) < 1 || n_cs(s2) < 1)
    return(data.table(gene = g, n_snps = nrow(mm),
                      n_cs_eqtl = n_cs(s1), n_cs_gwas = n_cs(s2),
                      status = "SuSiE found no credible set in one or both traits"))

  cs <- tryCatch(suppressWarnings(suppressMessages(
          coloc.susie(s1, s2, p12 = CFG$p12))), error = function(e) NULL)
  if (is.null(cs) || is.null(cs$summary) || !nrow(cs$summary))
    return(data.table(gene = g, n_snps = nrow(mm),
                      n_cs_eqtl = n_cs(s1), n_cs_gwas = n_cs(s2),
                      status = "coloc.susie returned no pair"))

  su <- as.data.table(cs$summary)
  setorder(su, -PP.H4.abf)
  b <- su[1]
  data.table(gene = g, n_snps = nrow(mm),
             n_cs_eqtl = n_cs(s1), n_cs_gwas = n_cs(s2),
             n_pairs = nrow(su),
             best_PP_H3 = round(b$PP.H3.abf, 3),
             best_PP_H4 = round(b$PP.H4.abf, 3),
             max_PP_H4_any_pair = round(max(su$PP.H4.abf), 3),
             status = "ok")
}

out <- list()
for (g in targets) {
  r <- run_gene(g)
  out[[g]] <- r
  say("  %-10s %s", g,
      if (identical(r$status, "ok"))
        sprintf("%d SNPs | eQTL sets %d, GWAS sets %d | best PP.H4 %.3f (PP.H3 %.3f)",
                r$n_snps, r$n_cs_eqtl, r$n_cs_gwas, r$best_PP_H4, r$best_PP_H3)
      else r$status)
}
res <- rbindlist(out, fill = TRUE)

res[, susie_verdict := fifelse(
      status != "ok",                       paste0("NOT RESOLVABLE - ", status),
 fifelse(max_PP_H4_any_pair >= CFG$pp4_strong,
                                            "COLOCALISED with an RA signal (multi-variant model)",
 fifelse(best_PP_H3 >= 0.80,                "DISTINCT VARIANTS (valid: multiple signals modelled)",
                                            "INCONCLUSIVE under the multi-variant model")))]
res[, in_panel := gene %in% panel_genes]
setorder(res, -max_PP_H4_any_pair, na.last = TRUE)
fwrite(res, file.path(tab, "COLOC_SUSIE_mhc.csv"))

hdr("SuSiE RESULTS")
print(res[, .(gene, in_panel, n_snps, n_cs_eqtl, n_cs_gwas,
              best_PP_H4 = max_PP_H4_any_pair, susie_verdict)])

# =============================================================================
# STEP 3 — COMBINED abf + susie VIEW
# =============================================================================
hdr("STEP 3  COMBINED VIEW")
comb <- merge(abf[, .(gene, MHC_gene, in_panel_female, in_panel_male,
                      abf_PP3 = round(PP3, 3), abf_PP4 = round(PP4, 3),
                      abf_verdict = coloc_verdict)],
              res[, .(gene, susie_n_cs_eqtl = n_cs_eqtl, susie_n_cs_gwas = n_cs_gwas,
                      susie_PP4 = max_PP_H4_any_pair, susie_verdict)],
              by = "gene", all.x = TRUE)
comb[, final_verdict := fifelse(
       MHC_gene == FALSE, abf_verdict,
  fifelse(is.na(susie_verdict), "MHC - unresolved (susie not run)",
  fifelse(grepl("^COLOCALISED", susie_verdict), susie_verdict,
  fifelse(grepl("^DISTINCT",    susie_verdict), susie_verdict,
          paste0("MHC - ", susie_verdict)))))]
setorder(comb, -abf_PP4)
fwrite(comb, file.path(tab, "COLOC_combined_abf_susie.csv"))
print(comb[MHC_gene == TRUE, .(gene, abf_PP4, susie_PP4, final_verdict)])

saveRDS(list(susie = res, combined = comb, config = CFG),
        file.path(proc, "coloc_susie_objects.rds"))

hdr("READING RULE")
say("  This SUPERSEDES the coloc.abf verdict for MHC genes only. Outside the MHC,")
say("  coloc.abf's single-causal-variant assumption is tenable and 10d stands.")
say("  Report susie results as INDICATIVE, not definitive: the LD matrix is an")
say("  out-of-sample 1000G EUR reference, and in the MHC an LD mismatch can create")
say("  or destroy credible sets (Zou 2022). In-sample LD would be needed to settle")
say("  the region, and neither eQTLGen nor Okada release it.")
say("  A 'NOT RESOLVABLE' row is a statement about power and LD reference quality,")
say("  never evidence of absence.")
say("")
say("Wrote COLOC_SUSIE_mhc.csv and COLOC_combined_abf_susie.csv")
cat("\nDONE\n")
