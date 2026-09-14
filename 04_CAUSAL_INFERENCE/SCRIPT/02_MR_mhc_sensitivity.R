#!/usr/bin/env Rscript
# MHC-excluded sensitivity analysis for 10_MR.R: re-runs MR offline from cached harmonised data with MHC instruments dropped, and reports a parallel column against the primary run so causal claims can be checked for HLA-DRB1 confounding.
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

# Step 0: load the cached primary MR run
hdr("STEP 0  LOAD CACHED PRIMARY MR OBJECTS")
obj_path <- file.path(proc, "MR_primary_objects.rds")
if (!file.exists(obj_path))
  stop("Missing ", obj_path, " - run scripts/04_causal/10_MR.R first.")
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

# Step 1: drop MHC instruments from the harmonised data, keyed on (SNP, id.exposure)
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

# Step 2: re-run MR on the MHC-free instrument set (estimator hierarchy re-applied since instrument count can change)
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

# Step 3: per-stratum FDR (denominator shrinks - genes untestable without MHC are excluded) and the parallel comparison table
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

  # distinguish genuine MHC-dependence (instrument set changed) from pure BH rank shift (estimate unchanged)
  cmp[, estimate_unchanged := !is.na(p_noMHC) &
        abs(p_primary - p_noMHC) < 1e-12 & nSNP_primary == nSNP_noMHC]

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
  say("  MHC-free prioritised set (%d): %s", nrow(fs_no),
      if (nrow(fs_no)) paste(fs_no$gene, collapse = ", ") else "none")

  list(sex = sx, cmp = cmp, fs_no = fs_no,
       n_tested_primary = nrow(pr_p), n_tested_noMHC = nrow(pr_n),
       n_sig_primary = sum(sig_p), n_sig_noMHC = sum(sig_n))
}

F <- compare_stratum("female", fem)
M <- compare_stratum("male",   mal)

# Step 4: fate of the FDR-surviving set and of the final panels when MHC is removed
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
  say("-- %s: fate of the %d FDR-surviving prioritised genes --", toupper(sx),
      sum(fate$sex == sx & fate$in_FS_input_primary))
  print(fate[sex == sx & in_FS_input_primary == TRUE,
             .(gene, chr = instrument_chr, MHC = MHC_gene, panel = in_final_panel,
               OR_pri = OR_primary, FDR_pri = FDR_primary,
               OR_no = OR_noMHC, FDR_no = FDR_noMHC, verdict)])
}

# Step 5: headline summary
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
