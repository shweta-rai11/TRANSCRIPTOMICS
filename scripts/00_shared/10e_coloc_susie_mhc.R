#!/usr/bin/env Rscript
# Multiple-causal-variant colocalisation (coloc.susie) for MHC genes, where coloc.abf's single-causal-variant assumption fails. Uses SuSiE credible sets per trait + an out-of-sample 1000G EUR LD reference; results are indicative, not definitive, and supersede 10d's coloc.abf verdict for MHC genes only.
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

# Step 1: target genes - everything 10d could not resolve
hdr("STEP 1  TARGETS")
regions <- readRDS(file.path(proc, "coloc_regions.rds"))
abf <- fread(file.path(tab, "COLOC_results.csv"))

targets <- abf[MHC_gene == TRUE]$gene
say("MHC genes unresolved by coloc.abf : %d", length(targets))
say("  %s", paste(targets, collapse = ", "))
panel_genes <- abf[in_panel_female == TRUE | in_panel_male == TRUE]$gene
say("of which in a final panel          : %s",
    paste(intersect(targets, panel_genes), collapse = ", "))

# Step 2: harmonise, fetch LD, run SuSiE on each trait, then coloc.susie
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

  # align LD-matrix orientation to our effect allele (rownames "rsid_A1_A2")
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

# Step 3: combined abf + susie view
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
