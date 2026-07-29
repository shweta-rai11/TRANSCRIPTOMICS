#!/usr/bin/env Rscript
# =============================================================================
# 05b_dge_sensitivity.R  —  Two checks the thesis needs and 05_dge.R does not do.
#
# CHECK 1  ComBat-then-DE vs batch-in-the-model
#   05_dge.R runs limma on the ComBat-corrected matrix with design ~ group.
#   The model therefore treats ComBat's output as raw data and never debits the
#   degrees of freedom ComBat spent, which inflates significance - worst when the
#   design is unbalanced across batches (Nygaard, Rodland & Hovig, Biostatistics
#   2016;17:29-39, already cited in 05_dge.R). This re-runs DE on the PRE-ComBat
#   quantile-normalised matrix with batch as a model term, and compares.
#   Motivating anomaly: the male stratum (n=38) returned MORE DEGs (5,820) than
#   the female stratum (n=145, 5,131). If that inverts here, the ComBat->DE path
#   was the cause.
#
# CHECK 2  group x sex INTERACTION
#   The thesis claims sex-stratified biology but never tests a sex interaction.
#   ComBat protected only main effects (mod = ~group + sex), so the interaction
#   is tested here on the PRE-ComBat matrix with batch in the model, where it was
#   never at risk of being absorbed.
# =============================================================================
suppressMessages({library(limma); library(data.table)})
set.seed(1234)
tab <- "results/tables"
o <- readRDS("data/processed/combined_train.rds")
qn <- o$expr_qnorm; cb <- o$expr; me <- as.data.table(o$meta)
me[, group := factor(group, levels = c("HC","RA"))]
me[, sex   := factor(sex,   levels = c("F","M"))]
me[, batch_full := factor(batch_full)]
FDR <- 0.05; LFC <- 0.1
say <- function(...) cat(sprintf(...), "\n", sep = "")

fit_de <- function(E, md, form, coefname) {
  d <- model.matrix(form, data = md)
  qrd <- qr(d)
  if (qrd$rank < ncol(d)) {                    # drop aliased columns, report it
    keep <- qrd$pivot[seq_len(qrd$rank)]
    say("      rank-deficient: %d of %d columns aliased -> dropped",
        ncol(d) - qrd$rank, ncol(d))
    d <- d[, sort(keep), drop = FALSE]
  }
  if (!coefname %in% colnames(d)) return(NULL)
  aw <- limma::arrayWeights(E, d)
  tt <- topTable(eBayes(lmFit(E, d, weights = aw)), coef = coefname,
                 number = Inf, sort.by = "none")
  list(n_sig = sum(tt$adj.P.Val < FDR & abs(tt$logFC) > LFC),
       n_fdr = sum(tt$adj.P.Val < FDR), tt = tt, resid_df = ncol(E) - ncol(d))
}

cat("\n===== CHECK 1: ComBat->DE  vs  batch-in-model (pre-ComBat) =====\n")
res <- rbindlist(lapply(list(
  list(lab="All",    idx = rep(TRUE, nrow(me))),
  list(lab="Female", idx = me$sex == "F"),
  list(lab="Male",   idx = me$sex == "M")), function(s) {
    md <- droplevels(me[s$idx]); cols <- md$sample
    say("  [%s] n=%d (RA=%d HC=%d), batches=%d", s$lab, length(cols),
        sum(md$group=="RA"), sum(md$group=="HC"), nlevels(md$batch_full))
    a <- fit_de(cb[, cols], md, ~ group, "groupRA")                  # as published
    b <- fit_de(qn[, cols], md, ~ group + batch_full, "groupRA")     # honest
    data.table(comparison = s$lab, n = length(cols),
               combat_then_DE = a$n_sig, batch_in_model = b$n_sig,
               ratio = round(b$n_sig / max(a$n_sig, 1), 3),
               resid_df_combat = a$resid_df, resid_df_batch = b$resid_df)
  }))
print(res); fwrite(res, file.path(tab, "DEG_sensitivity_combat_vs_batch.csv"))

cat("\n===== CHECK 2: group x sex INTERACTION (pre-ComBat, batch in model) =====\n")
i_qn <- fit_de(qn, me, ~ group * sex + batch_full, "groupRA:sexM")
i_cb <- fit_de(cb, me, ~ group * sex,              "groupRA:sexM")
say("  pre-ComBat + batch : %d genes with FDR<0.05 interaction", i_qn$n_fdr)
say("  ComBat matrix      : %d genes with FDR<0.05 interaction", i_cb$n_fdr)
int <- data.table(model = c("qnorm + batch_full", "ComBat"),
                  n_interaction_FDR05 = c(i_qn$n_fdr, i_cb$n_fdr))
print(int); fwrite(int, file.path(tab, "DEG_interaction_summary.csv"))
if (i_qn$n_fdr > 0) {
  top <- as.data.table(i_qn$tt, keep.rownames = "gene")[order(adj.P.Val)][1:min(25,.N),
          .(gene, interaction_logFC = round(logFC,3), adj.P.Val = signif(adj.P.Val,3))]
  print(top); fwrite(top, file.path(tab, "DEG_interaction_top.csv"))
}
cat("\nDONE -> DEG_sensitivity_combat_vs_batch.csv, DEG_interaction_summary.csv\n")
