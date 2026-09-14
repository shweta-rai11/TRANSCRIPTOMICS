#!/usr/bin/env Rscript
# Full report of the diagnosis x sex interaction: primary (pre-ComBat+batch) vs ComBat-sensitivity models, per-sex effect patterns, and enrichment
suppressMessages({
  library(limma); library(data.table)
})
options(stringsAsFactors = FALSE)
set.seed(1234)

proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)

CFG <- list(fdr = 0.05, lfc = 0.1, comp_pcs = 3, top_n = 25)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# Step 1: load data
hdr("STEP 1  DATA")
o  <- readRDS(file.path(proc, "combined_train.rds"))
qn <- o$expr_qnorm          # pre-ComBat, quantile-normalised -> primary model
cb <- o$expr                # ComBat-corrected -> sensitivity only
me <- as.data.table(o$meta)
me[, group := factor(group, levels = c("HC", "RA"))]
me[, sex   := factor(sex,   levels = c("F", "M"))]
me[, batch_full := factor(batch_full)]

say("samples %d | female %d (RA %d) | male %d (RA %d) | batches %d",
    nrow(me), sum(me$sex == "F"), sum(me$sex == "F" & me$group == "RA"),
    sum(me$sex == "M"), sum(me$sex == "M" & me$group == "RA"),
    nlevels(me$batch_full))
say("genes %d", nrow(qn))
say("NOTE: interaction power is ~1/4 of main-effect power at the same n, and the")
say("      male stratum contributes only %d RA cases. Treat the count as a FLOOR.",
    sum(me$sex == "M" & me$group == "RA"))

# Step 2: fit the diagnosis x sex interaction model
hdr("STEP 2  DIAGNOSIS x SEX INTERACTION")

fit_int <- function(E, md, form, coefname, label) {
  d <- stats::model.matrix(form, data = md)
  qrd <- qr(d)
  if (qrd$rank < ncol(d)) {
    keep <- sort(qrd$pivot[seq_len(qrd$rank)])
    say("  [%s] rank-deficient: %d of %d columns aliased -> dropped",
        label, ncol(d) - qrd$rank, ncol(d))
    d <- d[, keep, drop = FALSE]
  }
  if (!coefname %in% colnames(d)) stop("coef ", coefname, " absent in ", label)
  aw  <- limma::arrayWeights(E, d)
  fit <- limma::eBayes(limma::lmFit(E, d, weights = aw))
  tt  <- limma::topTable(fit, coef = coefname, number = Inf, sort.by = "none")
  say("  [%-28s] FDR<%.2f: %5d genes | with |logFC|>%.1f: %5d | resid df %d",
      label, CFG$fdr, sum(tt$adj.P.Val < CFG$fdr), CFG$lfc,
      sum(tt$adj.P.Val < CFG$fdr & abs(tt$logFC) > CFG$lfc), ncol(E) - ncol(d))
  list(tt = tt, resid_df = ncol(E) - ncol(d), n_fdr = sum(tt$adj.P.Val < CFG$fdr))
}

primary <- fit_int(qn, me, ~ group * sex + batch_full, "groupRA:sexM",
                   "PRIMARY qnorm + batch")
combat  <- fit_int(cb, me, ~ group * sex, "groupRA:sexM",
                   "sensitivity ComBat")

# composition-adjusted interaction (05c output, if available)
adjusted <- NULL
cf_path <- file.path(procN, "cell_fractions.rds")
if (file.exists(cf_path)) {
  cf <- readRDS(cf_path)
  cpc <- as.data.table(cf$comp_pcs)
  md2 <- merge(me, cpc, by = "sample", sort = FALSE)
  E2  <- qn[, md2$sample, drop = FALSE]
  frm <- stats::as.formula(paste("~ group * sex + batch_full +",
                                 paste0("cPC", seq_len(CFG$comp_pcs), collapse = " + ")))
  adjusted <- fit_int(E2, md2, frm, "groupRA:sexM", "composition-adjusted")
} else {
  say("  cell_fractions.rds not found - composition-adjusted interaction skipped")
}

mc <- data.table(
  model = c("PRIMARY: qnorm + batch_full (reported)",
            "sensitivity: ComBat matrix, no batch term",
            "sensitivity: qnorm + batch + composition PCs"),
  matrix_used = c("pre-ComBat quantile-normalised", "ComBat-corrected",
                  "pre-ComBat quantile-normalised"),
  n_interaction_FDR05 = c(primary$n_fdr, combat$n_fdr,
                          if (is.null(adjusted)) NA_integer_ else adjusted$n_fdr),
  resid_df = c(primary$resid_df, combat$resid_df,
               if (is.null(adjusted)) NA_integer_ else adjusted$resid_df),
  note = c("conservative; batch modelled explicitly, df correctly debited",
           "INFLATED: ComBat mod = ~group+sex protects main effects only, not the interaction",
           "tests whether the interaction survives leukocyte-composition adjustment"))
fwrite(mc, file.path(tab, "DEG_interaction_model_comparison.csv"))
print(mc)

# Step 3: within-sex effects behind each interaction, classified into patterns
hdr("STEP 3  WITHIN-SEX EFFECTS AND PATTERN CLASSIFICATION")

fit_within <- function(sx) {
  md <- droplevels(me[sex == sx]); cols <- md$sample
  d <- stats::model.matrix(~ group + batch_full, data = md)
  qrd <- qr(d)
  if (qrd$rank < ncol(d)) d <- d[, sort(qrd$pivot[seq_len(qrd$rank)]), drop = FALSE]
  E <- qn[, cols, drop = FALSE]
  aw <- limma::arrayWeights(E, d)
  tt <- limma::topTable(limma::eBayes(limma::lmFit(E, d, weights = aw)),
                        coef = "groupRA", number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, P = tt$P.Value, FDR = tt$adj.P.Val)
}
wf <- fit_within("F"); wm <- fit_within("M")
say("within-sex models fitted on the SAME pre-ComBat matrix and batch design as")
say("the interaction, so the three estimates are mutually consistent.")

it <- as.data.table(primary$tt, keep.rownames = "gene")
full <- Reduce(function(a, b) merge(a, b, by = "gene"), list(
  it[, .(gene, interaction_logFC = logFC, interaction_P = P.Value,
         interaction_FDR = adj.P.Val)],
  wf[, .(gene, female_logFC = logFC, female_FDR = FDR)],
  wm[, .(gene, male_logFC = logFC, male_FDR = FDR)]))

if (!is.null(adjusted)) {
  ad <- as.data.table(adjusted$tt, keep.rownames = "gene")
  full <- merge(full, ad[, .(gene, interaction_logFC_celladj = logFC,
                             interaction_FDR_celladj = adj.P.Val)],
                by = "gene", all.x = TRUE)
}

f_sig <- full$female_FDR < CFG$fdr
m_sig <- full$male_FDR   < CFG$fdr
full[, pattern := fifelse(
      f_sig & m_sig & sign(female_logFC) != sign(male_logFC), "OPPOSITE direction",
 fifelse(f_sig & m_sig,                                       "MAGNITUDE difference",
 fifelse(f_sig & !m_sig,                                      "FEMALE-restricted",
 fifelse(!f_sig & m_sig,                                      "MALE-restricted",
                                                              "neither sex significant"))))]
setorder(full, interaction_P)
fwrite(full, file.path(tab, "DEG_interaction_full.csv"))

sig <- full[interaction_FDR < CFG$fdr]
fwrite(sig, file.path(tab, "DEG_interaction_significant.csv"))

say("")
say("interaction genes at FDR < %.2f : %d", CFG$fdr, nrow(sig))
pat <- sig[, .N, by = pattern][order(-N)]
fwrite(pat, file.path(tab, "DEG_interaction_patterns.csv"))
print(pat)

if (!is.null(adjusted) && "interaction_FDR_celladj" %in% names(sig)) {
  surv <- sum(sig$interaction_FDR_celladj < CFG$fdr, na.rm = TRUE)
  say("of these, %d of %d survive leukocyte-composition adjustment (%.0f%%)",
      surv, nrow(sig), 100 * surv / max(1, nrow(sig)))
}

say("")
say("top %d interaction genes:", CFG$top_n)
print(head(sig[, .(gene,
                   int_logFC = round(interaction_logFC, 3),
                   int_FDR = signif(interaction_FDR, 3),
                   F_logFC = round(female_logFC, 3), F_FDR = signif(female_FDR, 2),
                   M_logFC = round(male_logFC, 3),   M_FDR = signif(male_FDR, 2),
                   pattern)], CFG$top_n))

# panel cross-reference: are any panel genes sex-differential?
ml_path <- file.path(procN, "ml_features.rds")
if (file.exists(ml_path)) {
  ml <- readRDS(ml_path)
  pg <- union(ml$female$consensus, ml$male$consensus)
  ov <- sig[gene %in% pg]
  say("")
  say("panel genes among the interaction set: %d of %d", nrow(ov), length(pg))
  if (nrow(ov)) print(ov[, .(gene, interaction_FDR = signif(interaction_FDR, 3), pattern)])
  else say("  none - the panels are sex-STRATIFIED, not sex-SPECIFIC. This is the")
  say("  expected result and must be stated plainly in the thesis: the panel genes")
  say("  are not shown to behave differently between the sexes.")
}

# Step 4: functional enrichment of the interaction genes
hdr("STEP 4  ENRICHMENT OF THE INTERACTION GENES")
enr <- NULL
ok <- requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (!ok) {
  say("clusterProfiler / org.Hs.eg.db unavailable - enrichment skipped")
} else if (nrow(sig) < 10) {
  say("only %d interaction genes - too few for over-representation analysis", nrow(sig))
} else {
  suppressMessages({library(clusterProfiler); library(org.Hs.eg.db)})
  univ <- full$gene
  e2 <- suppressMessages(bitr(sig$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db))
  eu <- suppressMessages(bitr(univ,     "SYMBOL", "ENTREZID", org.Hs.eg.db))
  say("mapped %d of %d interaction genes to Entrez (universe %d)",
      nrow(e2), nrow(sig), nrow(eu))
  go <- tryCatch(as.data.table(suppressMessages(
          enrichGO(e2$ENTREZID, OrgDb = org.Hs.eg.db, universe = eu$ENTREZID,
                   ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2,
                   readable = TRUE))), error = function(e) NULL)
  kg <- tryCatch(as.data.table(suppressMessages(
          enrichKEGG(e2$ENTREZID, organism = "hsa", universe = eu$ENTREZID,
                     pAdjustMethod = "BH", qvalueCutoff = 0.2))),
          error = function(e) NULL)
  parts <- list()
  if (!is.null(go) && nrow(go)) parts$GO   <- cbind(source = "GO:BP", go)
  if (!is.null(kg) && nrow(kg)) parts$KEGG <- cbind(source = "KEGG",  kg)
  if (length(parts)) {
    enr <- rbindlist(parts, fill = TRUE)
    fwrite(enr, file.path(tab, "DEG_interaction_enrichment.csv"))
    say("enriched terms: %d", nrow(enr))
    print(head(enr[, .(source, Description, GeneRatio, p.adjust = signif(p.adjust, 3))], 15))
  } else {
    say("no terms enriched at q < 0.2 - report this as a NEGATIVE result, not as")
    say("an absence of analysis. With %d genes the ORA is underpowered.", nrow(sig))
    fwrite(data.table(source = character(), Description = character(),
                      note = "no term enriched at q<0.2"),
           file.path(tab, "DEG_interaction_enrichment.csv"))
  }
}

saveRDS(list(primary = primary, combat = combat, adjusted = adjusted,
             full = full, significant = sig, patterns = pat,
             model_comparison = mc, enrichment = enr,
             within_female = wf, within_male = wm, config = CFG),
        file.path(procN, "interaction_objects.rds"))

hdr("READING RULE FOR THE THESIS")
say("  REPORT the %d-gene interaction set from the PRIMARY model. It is the only", nrow(sig))
say("  direct evidence in this project that the RA blood transcriptome differs")
say("  between the sexes, and it must appear in the Results.")
say("  DO NOT report the ComBat count (%d) as the headline - it is inflated.", combat$n_fdr)
say("  DO NOT claim the panels are sex-specific on the strength of it: the panels")
say("  and the interaction set are separate analyses answering separate questions.")
say("  DO state that interaction power is ~1/4 of main-effect power, so this is a")
say("  floor on the number of sex-differential genes, not an estimate of it.")
say("")
say("Wrote DEG_interaction_full.csv, DEG_interaction_significant.csv,")
say("DEG_interaction_patterns.csv, DEG_interaction_model_comparison.csv,")
say("DEG_interaction_enrichment.csv and interaction_objects.rds")
cat("\nDONE\n")
