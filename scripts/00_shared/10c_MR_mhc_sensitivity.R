#!/usr/bin/env Rscript
# =============================================================================
# 10c_MR_mhc_sensitivity.R  —  MHC-EXCLUDED sensitivity analysis for 10_MR.R
#
# WHY THIS EXISTS
#   10_MR.R runs with CFG$exclude_mhc = FALSE: instruments inside the extended
#   MHC (chr6:25-34 Mb, GRCh37) are FLAGGED but RETAINED. That is the correct
#   default for a primary analysis, but it leaves the strongest objection to the
#   causal claim untested, because:
#
#     (a) HLA-DRB1 is by a wide margin the dominant RA susceptibility locus, and
#     (b) the extended MHC carries the most extensive long-range LD in the human
#         genome, so a cis-eQTL for ANY gene in the region is correlated with the
#         HLA-DRB1 risk haplotype whether or not that gene is causal.
#
#   The consequence is that a cis-MR estimate for an MHC gene cannot be
#   distinguished from the HLA-DRB1 signal read through a proxy. In the primary
#   run this is not a marginal issue: 14 of 32 female (44%) and 10 of 25 male
#   (40%) FDR-surviving "causal" genes are MHC-flagged, and 2 of 6 female
#   (C6orf136, GNL1) and 3 of 6 male (HLA-DMA, VPS52, plus MHC-adjacent) final
#   panel genes sit in the region - most on a SINGLE-SNP Wald ratio, for which
#   no pleiotropy test is even possible.
#
#   This script therefore re-runs the ENTIRE MR with the MHC instruments removed
#   and reports the result as a FULL PARALLEL COLUMN against the primary run, so
#   that every reported causal gene is read together with its MHC-excluded fate.
#
# WHAT IS AND IS NOT RE-DONE
#   Re-done : instrument set (MHC SNPs dropped), estimator choice (the hierarchy
#             is re-applied because dropping SNPs can demote a gene from IVW to
#             Wald ratio), heterogeneity/pleiotropy tests, per-stratum BH FDR
#             (the denominator legitimately shrinks - genes with no non-MHC
#             instrument are not "tested and null", they are UNTESTABLE).
#   NOT re-done : instrument extraction, LD clumping, cis filtering and outcome
#             harmonisation. All of that is read from the cached objects written
#             by 10_MR.R, so this script is FULLY OFFLINE and deterministic. It
#             cannot drift from the primary run, because it starts from the
#             primary run's own harmonised data.
#
# INTERPRETATION RULE ADOPTED FOR THE THESIS
#   A gene is reported as ROBUST only if it survives FDR < 0.05 in BOTH the
#   primary and the MHC-excluded analysis. A gene that survives only in the
#   primary analysis is reported as MHC-DEPENDENT and is NOT eligible to carry a
#   causal claim. Genes are not silently dropped: the fate of every one of them
#   is tabulated.
#
#   in : data/processed/new/MR_primary_objects.rds   (from 10_MR.R)
#        results/tables/candidates_{female,male}_disease.csv
#        results/tables/FS_input_{female,male}.csv
#        data/processed/new/ml_features.rds          (final panels, if present)
#   out: results/tables/MR_MHC_sensitivity_{female,male}.csv   parallel columns
#        results/tables/MR_MHC_sensitivity_summary.csv
#        results/tables/MR_MHC_sensitivity_panel_fate.csv
#        results/tables/FS_input_{female,male}_noMHC.csv
#        data/processed/new/MR_mhc_sensitivity_objects.rds
#
#   Vosa U, et al. Nat Genet 2021;53:1300-1310.        (eQTLGen)
#   Okada Y, et al. Nature 2014;506:376-381.           (RA GWAS; HLA-DRB1)
#   Raychaudhuri S, et al. Nat Genet 2012;44:291-296.  (MHC fine-mapping in RA)
#   Trynka G, et al. Nat Genet 2011;43:1193-1201.      (long-range MHC LD)
# =============================================================================
suppressMessages({
  library(TwoSampleMR); library(dplyr); library(data.table)
})
options(stringsAsFactors = FALSE)
set.seed(2024)

proc <- "data/processed/new"; tab <- "results/tables"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)
FDR_CUT <- 0.05

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# =============================================================================
# STEP 0 — LOAD THE CACHED PRIMARY RUN
# =============================================================================
hdr("STEP 0  LOAD CACHED PRIMARY MR OBJECTS")
obj_path <- file.path(proc, "MR_primary_objects.rds")
if (!file.exists(obj_path))
  stop("Missing ", obj_path, " - run scripts/00_shared/10_MR.R first.")
o <- readRDS(obj_path)

primary <- as.data.table(o$primary)      # one row per gene, MHC RETAINED
inst    <- as.data.table(o$inst)         # cis instruments, with $MHC flag
dat     <- as.data.frame(o$dat)          # harmonised exposure-outcome rows

fem <- fread(file.path(tab, "candidates_female_disease.csv"))$gene
mal <- fread(file.path(tab, "candidates_male_disease.csv"))$gene

say("cached harmonised rows : %d over %d genes", nrow(dat), uniqueN(dat$gene))
say("cached cis instruments : %d SNPs (%d MHC-flagged in %d genes)",
    nrow(inst), sum(inst$MHC), uniqueN(inst[MHC == TRUE]$gene))
say("candidate genes        : female %d | male %d", length(fem), length(mal))

# =============================================================================
# STEP 1 — DROP MHC INSTRUMENTS FROM THE HARMONISED DATA
# -----------------------------------------------------------------------------
# The MHC flag lives on `inst`, keyed by the (SNP, id.exposure) PAIR - the same
# key 10_MR.R uses to attach gene labels, and for the same reason: a SNP can
# instrument more than one gene, so matching on SNP alone would mislabel rows.
# =============================================================================
hdr("STEP 1  REMOVE MHC INSTRUMENTS")
key_dat  <- paste(dat$SNP, dat$id.exposure)
key_inst <- paste(inst$SNP, inst$id.exposure)
dat$MHC  <- inst$MHC[match(key_dat, key_inst)]
stopifnot(!any(is.na(dat$MHC)))

n_mhc_rows <- sum(dat$MHC)
dat_no <- dat[!dat$MHC, , drop = FALSE]

genes_before <- unique(dat$gene[dat$mr_keep])
genes_after  <- unique(dat_no$gene[dat_no$mr_keep])
lost_genes   <- setdiff(genes_before, genes_after)

say("harmonised rows removed        : %d of %d", n_mhc_rows, nrow(dat))
say("genes with >=1 usable instrument: %d -> %d", length(genes_before), length(genes_after))
say("genes left with NO instrument   : %d  (untestable after MHC exclusion)", length(lost_genes))

# =============================================================================
# STEP 2 — RE-RUN MR ON THE MHC-FREE INSTRUMENT SET
# -----------------------------------------------------------------------------
# The estimator hierarchy is re-applied rather than reused. Dropping a SNP can
# take a gene from 3 instruments (IVW + Egger + weighted median) to 2 (IVW only)
# or to 1 (Wald ratio), and the primary estimate must follow the instrument
# count that actually remains.
# =============================================================================
hdr("STEP 2  RE-ESTIMATE MR WITHOUT MHC INSTRUMENTS")
res <- list(); het <- list(); pleio <- list()
gs <- unique(dat_no$gene)
for (g in gs) {
  d <- dat_no[dat_no$gene == g & dat_no$mr_keep, , drop = FALSE]
  if (nrow(d) < 1) next
  ms <- if (nrow(d) >= 3) c("mr_ivw", "mr_egger_regression", "mr_weighted_median") else
        if (nrow(d) == 2) "mr_ivw" else "mr_wald_ratio"
  r <- tryCatch(mr(d, method_list = ms), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) next
  r$gene <- g; r$nSNP <- nrow(d); res[[g]] <- r
  if (nrow(d) >= 3) {
    het[[g]]   <- tryCatch(cbind(gene = g, mr_heterogeneity(d)),   error = function(e) NULL)
    pleio[[g]] <- tryCatch(cbind(gene = g, mr_pleiotropy_test(d)), error = function(e) NULL)
  }
}
res <- bind_rows(res)
stopifnot(nrow(res) > 0)
say("MR re-estimated for %d genes", uniqueN(res$gene))

# identical primary-estimator hierarchy to 10_MR.R (all dplyr verbs qualified: D1)
prim_no <- res %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(factor(method, levels = c("Inverse variance weighted",
                                           "Wald ratio",
                                           "Weighted median",
                                           "MR Egger")), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  as.data.table()
prim_no[, OR := exp(b)][, OR_lo := exp(b - 1.96 * se)][, OR_hi := exp(b + 1.96 * se)]
prim_no[, sensitivity_testable := nSNP >= 3]
het <- bind_rows(het); pleio <- bind_rows(pleio)

say("estimator mix after MHC exclusion:")
print(prim_no[, .N, by = method][order(-N)])

# =============================================================================
# STEP 3 — PER-STRATUM FDR AND THE PARALLEL COMPARISON TABLE
# -----------------------------------------------------------------------------
# NOTE ON THE FDR DENOMINATOR. It shrinks, and that is correct. A gene whose only
# instrument was an MHC SNP has not been "tested and found null"; it has become
# UNTESTABLE without the MHC. Keeping it in the BH denominator would penalise the
# genes that ARE testable for the absence of evidence about genes that are not.
# The shrinkage is reported explicitly so the change is never silent.
# =============================================================================
hdr("STEP 3  PER-STRATUM FDR AND PARALLEL COLUMNS")

compare_stratum <- function(sx, sex_genes) {
  pr_p <- copy(primary[gene %in% sex_genes])          # primary: MHC retained
  pr_p[, FDR_stratum := p.adjust(pval, "BH")]

  pr_n <- copy(prim_no[gene %in% sex_genes])          # sensitivity: MHC excluded
  pr_n[, FDR_stratum := p.adjust(pval, "BH")]

  cmp <- merge(
    pr_p[, .(gene, MHC_gene, instrument_chr,
             nSNP_primary = nSNP, method_primary = method,
             OR_primary = round(OR, 3),
             OR_lo_primary = round(OR_lo, 3), OR_hi_primary = round(OR_hi, 3),
             p_primary = signif(pval, 3), FDR_primary = signif(FDR_stratum, 3))],
    pr_n[, .(gene, nSNP_noMHC = nSNP, method_noMHC = method,
             OR_noMHC = round(OR, 3),
             OR_lo_noMHC = round(OR_lo, 3), OR_hi_noMHC = round(OR_hi, 3),
             p_noMHC = signif(pval, 3), FDR_noMHC = signif(FDR_stratum, 3))],
    by = "gene", all.x = TRUE)

  sig_p <- !is.na(cmp$FDR_primary) & cmp$FDR_primary < FDR_CUT
  sig_n <- !is.na(cmp$FDR_noMHC)   & cmp$FDR_noMHC   < FDR_CUT
  testable_n <- !is.na(cmp$FDR_noMHC)

  # A gene can lose FDR significance for TWO completely different reasons, and
  # conflating them would overstate the damage:
  #   (i)  its own instrument set changed - MHC SNPs were dropped from it, so the
  #        estimate itself is different. This is genuine MHC dependence.
  #   (ii) its instruments were all non-MHC and its estimate is bit-identical,
  #        but Benjamini-Hochberg is a RANK procedure: removing the very strong
  #        MHC genes above it shifts its rank i, and FDR = p * n / i rises even
  #        though nothing about the gene's own evidence changed.
  # Case (ii) is a multiplicity-bookkeeping effect, not evidence of confounding,
  # and is labelled separately.
  cmp[, estimate_unchanged := !is.na(p_noMHC) &
        abs(p_primary - p_noMHC) < 1e-12 & nSNP_primary == nSNP_noMHC]

  # Verdict vocabulary is deliberately blunt: an examiner should be able to read
  # a single column and know whether a causal claim is defensible.
  cmp[, verdict := fifelse(
        sig_p & sig_n,                "ROBUST (significant with and without MHC)",
   fifelse(sig_p & !testable_n,       "UNTESTABLE without MHC (no non-MHC instrument)",
   fifelse(sig_p & !sig_n & estimate_unchanged,
                                      "FDR-RANK ONLY (estimate identical; lost to BH re-ranking)",
   fifelse(sig_p & !sig_n,            "MHC-DEPENDENT (lost when MHC excluded)",
   fifelse(!sig_p & sig_n,            "GAINED (significant only without MHC)",
                                      "ns in both")))))]
  # direction agreement, only meaningful where both estimates exist
  cmp[, direction_consistent := fifelse(
        is.na(OR_noMHC), NA,
        (OR_primary > 1) == (OR_noMHC > 1))]

  setorder(cmp, p_primary)
  fwrite(cmp, file.path(tab, sprintf("MR_MHC_sensitivity_%s.csv", sx)))

  # the MHC-free replacement candidate list for feature selection
  fs_no <- pr_n[FDR_stratum < FDR_CUT,
                .(gene, direction = fifelse(OR > 1, "risk (OR>1)", "protective (OR<1)"),
                  MR_OR = round(OR, 3), MR_pval = signif(pval, 3),
                  MR_FDR_stratum = signif(FDR_stratum, 3), nSNP, method)]
  setorder(fs_no, MR_pval)
  fwrite(fs_no, file.path(tab, sprintf("FS_input_%s_noMHC.csv", sx)))

  say("== %s ==", toupper(sx))
  say("  genes tested       : %d primary -> %d after MHC exclusion (denominator -%d)",
      nrow(pr_p), nrow(pr_n), nrow(pr_p) - nrow(pr_n))
  say("  FDR<%.2f survivors  : %d primary -> %d after MHC exclusion",
      FDR_CUT, sum(sig_p), sum(sig_n))
  print(cmp[verdict != "ns in both", .N, by = verdict][order(-N)])
  say("  direction agreement where both estimable: %d of %d",
      sum(cmp$direction_consistent, na.rm = TRUE), sum(!is.na(cmp$direction_consistent)))
  say("  MHC-free causal set (%d): %s", nrow(fs_no),
      if (nrow(fs_no)) paste(fs_no$gene, collapse = ", ") else "none")

  list(sex = sx, cmp = cmp, fs_no = fs_no,
       n_tested_primary = nrow(pr_p), n_tested_noMHC = nrow(pr_n),
       n_sig_primary = sum(sig_p), n_sig_noMHC = sum(sig_n))
}

F <- compare_stratum("female", fem)
M <- compare_stratum("male",   mal)

# =============================================================================
# STEP 4 — FATE OF THE FDR-SURVIVING SET AND OF THE FINAL PANELS
# -----------------------------------------------------------------------------
# This is the table the thesis actually needs: for every gene that carries a
# causal claim in the primary analysis, and for every gene in the final
# diagnostic panels, what happens when the MHC is removed.
# =============================================================================
hdr("STEP 4  FATE OF CAUSAL GENES AND FINAL PANELS")

panels <- list(female = character(0), male = character(0))
ml_path <- file.path(proc, "ml_features.rds")
if (file.exists(ml_path)) {
  ml <- readRDS(ml_path)
  panels$female <- ml$female$consensus
  panels$male   <- ml$male$consensus
  say("final panels read from ml_features.rds: female %d, male %d genes",
      length(panels$female), length(panels$male))
} else {
  say("ml_features.rds not found - panel fate limited to the FS_input sets.")
}

fate <- rbindlist(lapply(list(F, M), function(r) {
  fs_primary <- fread(file.path(tab, sprintf("FS_input_%s.csv", r$sex)))$gene
  pn <- panels[[r$sex]]
  x <- r$cmp[gene %in% union(fs_primary, pn)]
  x[, sex := r$sex]
  x[, in_FS_input_primary := gene %in% fs_primary]
  x[, in_final_panel := gene %in% pn]
  x[, .(sex, gene, in_FS_input_primary, in_final_panel, MHC_gene, instrument_chr,
        nSNP_primary, method_primary, OR_primary, p_primary, FDR_primary,
        nSNP_noMHC, method_noMHC, OR_noMHC, p_noMHC, FDR_noMHC,
        direction_consistent, estimate_unchanged, verdict)]
}))
setorder(fate, sex, -in_final_panel, FDR_primary)
fwrite(fate, file.path(tab, "MR_MHC_sensitivity_panel_fate.csv"))

for (sx in c("female", "male")) {
  say("")
  say("-- %s: fate of the %d FDR-surviving causal genes --", toupper(sx),
      sum(fate$sex == sx & fate$in_FS_input_primary))
  print(fate[sex == sx & in_FS_input_primary == TRUE,
             .(gene, chr = instrument_chr, MHC = MHC_gene, panel = in_final_panel,
               OR_pri = OR_primary, FDR_pri = FDR_primary,
               OR_no = OR_noMHC, FDR_no = FDR_noMHC, verdict)])
}

# =============================================================================
# STEP 5 — HEADLINE SUMMARY
# =============================================================================
hdr("STEP 5  SUMMARY")
panel_fate_counts <- function(sx) {
  x <- fate[sex == sx & in_final_panel == TRUE]
  if (!nrow(x)) return(c(n = 0L, robust = 0L, untestable = 0L,
                         mhc_dep = 0L, rank_only = 0L))
  c(n          = nrow(x),
    robust     = sum(grepl("^ROBUST",        x$verdict)),
    untestable = sum(grepl("^UNTESTABLE",    x$verdict)),
    mhc_dep    = sum(grepl("^MHC-DEPENDENT", x$verdict)),
    rank_only  = sum(grepl("^FDR-RANK ONLY", x$verdict)))
}
pfF <- panel_fate_counts("female"); pfM <- panel_fate_counts("male")
cnt <- function(r, pat) sum(grepl(pat, r$cmp$verdict))

summ <- data.table(
  sex                      = c("Female", "Male"),
  genes_tested_primary     = c(F$n_tested_primary, M$n_tested_primary),
  genes_tested_noMHC       = c(F$n_tested_noMHC,   M$n_tested_noMHC),
  causal_FDR05_primary     = c(F$n_sig_primary,    M$n_sig_primary),
  causal_FDR05_noMHC       = c(F$n_sig_noMHC,      M$n_sig_noMHC),
  causal_robust            = c(cnt(F, "^ROBUST"),        cnt(M, "^ROBUST")),
  causal_untestable_noMHC  = c(cnt(F, "^UNTESTABLE"),    cnt(M, "^UNTESTABLE")),
  causal_MHC_dependent     = c(cnt(F, "^MHC-DEPENDENT"), cnt(M, "^MHC-DEPENDENT")),
  causal_FDR_rank_only     = c(cnt(F, "^FDR-RANK ONLY"), cnt(M, "^FDR-RANK ONLY")),
  panel_n                  = c(pfF["n"],          pfM["n"]),
  panel_robust             = c(pfF["robust"],     pfM["robust"]),
  panel_untestable         = c(pfF["untestable"], pfM["untestable"]),
  panel_MHC_dependent      = c(pfF["mhc_dep"],    pfM["mhc_dep"]),
  panel_FDR_rank_only      = c(pfF["rank_only"],  pfM["rank_only"]))
fwrite(summ, file.path(tab, "MR_MHC_sensitivity_summary.csv"))
print(summ)

saveRDS(list(summary = summ, female = F, male = M, fate = fate,
             primary_noMHC = prim_no, het_noMHC = het, pleio_noMHC = pleio,
             n_mhc_rows_removed = n_mhc_rows, lost_genes = lost_genes),
        file.path(proc, "MR_mhc_sensitivity_objects.rds"))

say("")
say("REPORTING RULE: a gene may carry a causal claim in the thesis only if its")
say("verdict is ROBUST. MHC-DEPENDENT and UNTESTABLE genes are reported as")
say("MHC-confounded and are described as associated, not causal.")
say("")
say("Wrote MR_MHC_sensitivity_{female,male}.csv, _summary.csv, _panel_fate.csv,")
say("FS_input_{female,male}_noMHC.csv and MR_mhc_sensitivity_objects.rds")
cat("\nDONE\n")
