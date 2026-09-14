#!/usr/bin/env Rscript
# Comparative MR of the per-sex EUR MR-prioritised genes against RA GWAS from three cohorts (Okada EUR discovery, Stahl EUR replication, BBJ EAS), with instrument-transferability diagnostics for the ancestry-mismatched EAS arm.
suppressMessages({library(TwoSampleMR); library(dplyr); library(data.table)})
set.seed(2024)
proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"

OKADA <- "ieu-a-832"   # EUR discovery
STAHL <- "ieu-a-834"   # EUR replication
BBJ   <- "bbj-a-151"   # EAS cross-ancestry

# ---- 1. EUR-prioritised gene set (per sex) + European eQTLGen instruments ----------
fsF <- fread(file.path(tab, "FS_input_female.csv"))
fsM <- fread(file.path(tab, "FS_input_male.csv"))
causal <- unique(c(fsF$gene, fsM$gene))
cat(sprintf("EUR-prioritised genes: %d female + %d male -> %d unique\n",
            nrow(fsF), nrow(fsM), length(causal)))

o    <- readRDS(file.path(procN, "MR_primary_objects.rds"))
inst <- as.data.table(o$inst)[gene %in% causal]        # F>10 European eQTLGen SNPs
cat(sprintf("Instruments for prioritised genes: %d SNPs / %d genes\n",
            nrow(inst), length(unique(inst$gene))))

# ---- 2. helpers --------------------------------------------------------------
# per-gene primary MR against one outcome (pipeline-identical method ladder)
mr_primary_vs <- function(dat) {
  res <- list()
  for (g in unique(dat$gene)) {
    d <- dat[dat$gene == g & dat$mr_keep, , drop = FALSE]; if (nrow(d) < 1) next
    ms <- if (nrow(d) >= 3) c("mr_ivw","mr_egger_regression","mr_weighted_median") else
          if (nrow(d) == 2) "mr_ivw" else "mr_wald_ratio"
    r <- tryCatch(mr(d, method_list = ms), error = function(e) NULL)
    if (!is.null(r) && nrow(r)) { r$gene <- g; r$nSNP_used <- nrow(d); res[[g]] <- r }
  }
  res <- bind_rows(res)
  if (!nrow(res)) return(data.table())
  p <- res %>% group_by(gene) %>%
    arrange(factor(method, levels = c("Inverse variance weighted","Wald ratio",
                                      "Weighted median","MR Egger")), .by_group = TRUE) %>%
    slice(1) %>% ungroup() %>% as.data.table()
  p[, OR := exp(b)][, `:=`(OR_lo = exp(b - 1.96*se), OR_hi = exp(b + 1.96*se))]
  p[, risk := ifelse(pval >= 0.05, "ns", ifelse(OR > 1, "risk (OR>1)", "protective (OR<1)"))]
  p[, .(gene, method, nSNP_used, b, se, pval, OR, OR_lo, OR_hi, risk)]
}

# fetch outcome, harmonise (action=2, pipeline-consistent), return dat + per-SNP eaf
fetch_dat <- function(outcome_id, label) {
  cat(sprintf("  fetching %s (%s) for %d SNPs ...\n", outcome_id, label, length(unique(inst$SNP))))
  out <- extract_outcome_data(snps = unique(inst$SNP), outcomes = outcome_id)
  if (is.null(out) || !nrow(out)) stop(sprintf("no outcome data for %s", outcome_id))
  dat <- harmonise_data(as.data.frame(inst), out, action = 2)
  dat$gene <- inst$gene[match(dat$SNP, inst$SNP)]
  as.data.table(dat)
}

# ---- 3. Okada (EUR) - reuse pipeline results; no refetch ---------------------
okada <- as.data.table(o$primary)[gene %in% causal,
           .(gene, OR_okada = OR, p_okada = pval, risk_okada = risk)]
dat_ok <- as.data.table(o$dat)[gene %in% causal]        # harmonised vs Okada (for SNP availability)

# ---- 4. Stahl (EUR) + BBJ (EAS) - fetch & MR ---------------------------------
dat_st <- fetch_dat(STAHL, "Stahl 2010 EUR")
dat_bj <- fetch_dat(BBJ,   "BBJ 2019 EAS")

st <- mr_primary_vs(dat_st)[, .(gene, OR_stahl = OR, p_stahl = pval, risk_stahl = risk)]
bj <- mr_primary_vs(dat_bj)[, .(gene, OR_bbj = OR, p_bbj = pval, risk_bbj = risk,
                                b_bbj = b, se_bbj = se)]

# ---- 5. instrument transferability across cohorts ----------------------------
# SNPs available (mr_keep) per gene in each outcome, + EUR-vs-EAS allele-freq gap
avail <- function(d, nm) d[mr_keep == TRUE, .(n = .N,
            eaf_gap = mean(abs(eaf.exposure - eaf.outcome), na.rm = TRUE)), by = gene] |>
            setNames(c("gene", paste0("nSNP_", nm), paste0("eafgap_", nm)))
tr <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
             list(data.table(gene = unique(inst$gene),
                             nSNP_instrument = inst[, .N, by = gene][match(unique(inst$gene), gene), N]),
                  avail(dat_ok, "okada"), avail(dat_st, "stahl"), avail(dat_bj, "bbj")))
tr[is.na(tr)] <- 0

# ---- 6. assemble comparison + per-sex classification -------------------------
comp <- Reduce(function(a, b) merge(a, b, by = "gene", all.x = TRUE),
               list(okada, st, bj))
comp <- merge(comp, tr[, .(gene, nSNP_instrument, nSNP_okada, nSNP_stahl, nSNP_bbj,
                           eafgap_bbj)], by = "gene", all.x = TRUE)

sd1 <- function(x) sign(x - 1)                          # OR direction vs the null
comp[, `:=`(
  dir_okada_eq_stahl = !is.na(OR_stahl) & sd1(OR_okada) == sd1(OR_stahl),
  dir_okada_eq_bbj   = !is.na(OR_bbj)   & sd1(OR_okada) == sd1(OR_bbj))]
comp[, `:=`(
  replicated_EUR  = dir_okada_eq_stahl & !is.na(p_stahl) & p_stahl < 0.05,
  transferable_EAS= dir_okada_eq_bbj   & !is.na(p_bbj)   & p_bbj   < 0.05,
  testable_EAS    = nSNP_bbj >= 1)]
comp[, ancestry_class := fifelse(!testable_EAS,          "untestable in EAS (no instruments)",
                        fifelse(transferable_EAS,        "shared EUR+EAS",
                        fifelse(replicated_EUR,          "EUR-replicated, not EAS",
                                                         "EUR-discovery only")))]

# ---- 7. emit per sex ---------------------------------------------------------
emit <- function(sx, genes) {
  d <- comp[gene %in% genes][order(p_okada)]
  keep <- d[, .(gene,
      OR_okada = round(OR_okada,3), p_okada = signif(p_okada,3), risk_okada,
      OR_stahl = round(OR_stahl,3), p_stahl = signif(p_stahl,3),
      OR_bbj   = round(OR_bbj,3),   p_bbj   = signif(p_bbj,3),
      dir_okada_eq_stahl, dir_okada_eq_bbj,
      replicated_EUR, transferable_EAS, testable_EAS, ancestry_class,
      nSNP_okada, nSNP_stahl, nSNP_bbj, eafgap_bbj = round(eafgap_bbj,3))]
  fwrite(keep, file.path(tab, sprintf("MR35_crossancestry_%s.csv", sx)))
  fwrite(d[, .(gene, nSNP_instrument, nSNP_okada, nSNP_stahl, nSNP_bbj,
               eafgap_bbj = round(eafgap_bbj,3), testable_EAS)],
         file.path(tab, sprintf("MR35_instrument_transferability_%s.csv", sx)))
  cat(sprintf("\n== %s (%d prioritised genes) ==\n", toupper(sx), nrow(d)))
  cl <- d[, .N, by = ancestry_class][order(-N)]; print(cl)
  cat(sprintf("  replicated in EUR (Stahl p<0.05, same dir): %d\n", sum(d$replicated_EUR)))
  cat(sprintf("  transferable to EAS (BBJ p<0.05, same dir): %d\n", sum(d$transferable_EAS)))
  cat(sprintf("  untestable in EAS (no surviving instruments): %d\n", sum(!d$testable_EAS)))
  data.table(sex = sx, n_causal = nrow(d),
             n_replicated_EUR = sum(d$replicated_EUR),
             n_transferable_EAS = sum(d$transferable_EAS),
             n_untestable_EAS = sum(!d$testable_EAS),
             n_EUR_only = sum(d$ancestry_class == "EUR-discovery only"))
}
summ <- rbindlist(list(emit("female", fsF$gene), emit("male", fsM$gene)))
fwrite(summ, file.path(tab, "MR35_crossancestry_summary.csv"))

saveRDS(list(comp = comp, transfer = tr, okada = okada, stahl = st, bbj = bj,
             dat_stahl = dat_st, dat_bbj = dat_bj, summary = summ),
        file.path(proc, "MR35_crossancestry_objects.rds"))

cat("\n=== SUMMARY (per sex) ===\n"); print(summ)
cat("\nDONE. Wrote MR35_crossancestry_{female,male}.csv, _transferability_*, _summary.csv\n")
