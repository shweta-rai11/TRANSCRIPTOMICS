#!/usr/bin/env Rscript
# Estimates blood cell composition (CIBERSORT/LM22, checked against MCP-counter) and re-fits DE with composition PCs (CLR) as covariates, to separate cell-intrinsic from compositional signal, within sex
suppressMessages({
  library(IOBR); library(limma); library(data.table)
})
options(stringsAsFactors = FALSE)
set.seed(1234)

proc <- "data/processed"; procN <- "data/processed/new"; tab <- "results/tables"
dir.create(procN, showWarnings = FALSE, recursive = TRUE)
dir.create(tab,   showWarnings = FALSE, recursive = TRUE)

CFG <- list(
  perm      = 100,     # CIBERSORT permutations for the deconvolution p-value
  comp_pcs  = 3,       # composition PCs entering the DE design
  fdr       = 0.05,
  lfc       = 0.1,     # matches 05_dge.R
  min_frac  = 0.001    # subsets below this mean fraction are dropped from CLR
)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# STEP 1: load the three blood datasets
hdr("STEP 1  LOAD DATASETS")
tr <- readRDS(file.path(proc, "combined_train.rds"))
ho <- readRDS(file.path(proc, "internal_val_holdout_processed.rds"))

train_expr <- tr$expr; train_meta <- as.data.table(tr$meta)
hold_expr  <- ho$expr; hold_meta  <- as.data.table(ho$meta)
say("train   : %d genes x %d samples", nrow(train_expr), ncol(train_expr))
say("holdout : %d genes x %d samples", nrow(hold_expr),  ncol(hold_expr))

# external blood (GSE15573), loaded exactly as script 13 does
ext_expr <- NULL; ext_meta <- NULL
ext_path <- "data/raw/GSE15573_raw.rds"
if (file.exists(ext_path)) {
  e <- readRDS(ext_path); if (is.list(e) && !inherits(e, "ExpressionSet")) e <- e[[1]]
  ex <- Biobase::exprs(e); ph <- Biobase::pData(e)
  if (stats::quantile(ex, 0.99, na.rm = TRUE) > 100) ex <- log2(ex + 1)
  sym <- Biobase::fData(e)
  symcol <- grep("^Symbol$|Gene.symbol|GENE_SYMBOL", colnames(sym), value = TRUE)[1]
  if (!is.na(symcol)) {
    g <- as.character(sym[[symcol]])
    keep <- !is.na(g) & g != "" & !grepl("///", g)
    ex <- ex[keep, , drop = FALSE]; g <- g[keep]
    o <- order(g, -rowMeans(ex, na.rm = TRUE))
    ex <- ex[o, , drop = FALSE]; g <- g[o]
    ex <- ex[!duplicated(g), , drop = FALSE]; rownames(ex) <- g[!duplicated(g)]
  }
  ttl <- as.character(ph$title)
  grp <- ifelse(grepl("control|healthy|normal", ttl, ignore.case = TRUE), "HC", "RA")
  ext_expr <- ex
  ext_meta <- data.table(sample = colnames(ex), group = grp)
  vm <- file.path(tab, "GSE15573_verified_metadata.csv")
  if (file.exists(vm)) {                       # prefer the curated sex/group calls
    v <- fread(vm)
    idc <- intersect(c("sample", "geo_accession", "gsm"), names(v))[1]
    if (!is.na(idc)) {
      ext_meta <- merge(ext_meta[, .(sample)], v, by.x = "sample", by.y = idc, all.x = TRUE)
      setnames(ext_meta, old = names(ext_meta), new = names(ext_meta))
    }
  }
  say("external: %d genes x %d samples", nrow(ext_expr), ncol(ext_expr))
} else {
  say("external GSE15573 not found at %s - skipped", ext_path)
}

# STEP 2: deconvolution (CIBERSORT LM22 + MCP-counter), fitted independently per dataset (holdout stays sealed)
hdr("STEP 2  CIBERSORT (LM22) AND MCP-COUNTER")

run_cibersort <- function(expr, label) {
  lin <- 2^expr                       # CIBERSORT expects linear-scale intensities
  lin[!is.finite(lin)] <- 0
  t0 <- Sys.time()
  r <- suppressMessages(suppressWarnings(
    deconvo_tme(eset = lin, method = "cibersort", arrays = TRUE, perm = CFG$perm)))
  r <- as.data.table(r)
  setnames(r, "ID", "sample")
  say("  CIBERSORT %-9s %d samples x %d columns  (%.1f s)",
      label, nrow(r), ncol(r) - 1L,
      as.numeric(difftime(Sys.time(), t0, units = "secs")))
  r
}

run_mcp <- function(expr, label) {
  r <- suppressMessages(suppressWarnings(
    deconvo_tme(eset = expr, method = "mcpcounter")))
  r <- as.data.table(r); setnames(r, "ID", "sample")
  say("  MCP-counter %-7s %d samples x %d scores", label, nrow(r), ncol(r) - 1L)
  r
}

cib <- list(train = run_cibersort(train_expr, "train"),
            holdout = run_cibersort(hold_expr, "holdout"))
mcp <- list(train = run_mcp(train_expr, "train"))
if (!is.null(ext_expr)) {
  cib$external <- run_cibersort(ext_expr, "external")
  mcp$external <- run_mcp(ext_expr, "external")
}

# LM22 fraction columns only (IOBR suffixes them "_CIBERSORT"; the three
# diagnostics P-value / Correlation / RMSE are kept aside, not modelled).
frac_cols <- function(d) setdiff(grep("_CIBERSORT$", names(d), value = TRUE),
                                 grep("^(P.value|Correlation|RMSE)",
                                      grep("_CIBERSORT$", names(d), value = TRUE),
                                      value = TRUE))
fc <- frac_cols(cib$train)
fc <- fc[!grepl("^P\\.value|^Correlation|^RMSE", fc)]
say("LM22 fraction columns modelled: %d", length(fc))

for (nm in names(cib)) fwrite(cib[[nm]], file.path(tab, sprintf("CELL_fractions_%s.csv", nm)))
for (nm in names(mcp)) fwrite(mcp[[nm]], file.path(tab, sprintf("CELL_mcpcounter_%s.csv", nm)))

# deconvolution quality on the training set
if ("P.value_CIBERSORT" %in% names(cib$train)) {
  pv <- cib$train$P.value_CIBERSORT
  say("CIBERSORT deconvolution p < 0.05 in %d of %d training samples (%.0f%%)",
      sum(pv < 0.05, na.rm = TRUE), length(pv), 100 * mean(pv < 0.05, na.rm = TRUE))
}

# STEP 3: does cell composition differ between RA and control, within sex
hdr("STEP 3  COMPOSITION VS DISEASE STATUS, WITHIN SEX")
ct <- merge(train_meta[, .(sample, group, sex)], cib$train[, c("sample", fc), with = FALSE],
            by = "sample")

tests <- rbindlist(lapply(c("F", "M"), function(sx) {
  d <- ct[sex == sx]
  rbindlist(lapply(fc, function(k) {
    x <- d[[k]]; g <- factor(d$group, levels = c("HC", "RA"))
    if (length(unique(g)) < 2 || stats::sd(x, na.rm = TRUE) == 0) return(NULL)
    w <- suppressWarnings(stats::wilcox.test(x ~ g))
    data.table(sex = ifelse(sx == "F", "Female", "Male"),
               cell = sub("_CIBERSORT$", "", k),
               mean_HC = round(mean(x[g == "HC"], na.rm = TRUE), 4),
               mean_RA = round(mean(x[g == "RA"], na.rm = TRUE), 4),
               diff    = round(mean(x[g == "RA"], na.rm = TRUE) -
                               mean(x[g == "HC"], na.rm = TRUE), 4),
               p = w$p.value)
  }))
}))
tests[, FDR := p.adjust(p, "BH"), by = sex]
setorder(tests, sex, p)
fwrite(tests, file.path(tab, "CELL_fraction_group_tests.csv"))

for (sx in c("Female", "Male")) {
  n_sig <- sum(tests$sex == sx & tests$FDR < CFG$fdr)
  say("%s: %d of %d LM22 subsets differ by disease status at FDR < %.2f",
      sx, n_sig, length(fc), CFG$fdr)
  print(head(tests[sex == sx, .(cell, mean_HC, mean_RA, diff, p = signif(p, 3),
                                FDR = signif(FDR, 3))], 8))
}

# STEP 4: composition principal components (CLR space)
hdr("STEP 4  COMPOSITION PCs (CLR)")
Fm <- as.matrix(ct[, fc, with = FALSE]); rownames(Fm) <- ct$sample
keep <- colMeans(Fm, na.rm = TRUE) >= CFG$min_frac
say("LM22 subsets retained for CLR (mean fraction >= %.3f): %d of %d",
    CFG$min_frac, sum(keep), length(keep))
Fk <- Fm[, keep, drop = FALSE]
Fk[Fk <= 0] <- min(Fk[Fk > 0], na.rm = TRUE) / 2      # zero replacement before log
clr <- log(Fk) - rowMeans(log(Fk))                     # centred log-ratio
pca <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
ve  <- pca$sdev^2 / sum(pca$sdev^2)
say("variance explained by PC1-PC%d: %s (cumulative %.1f%%)",
    CFG$comp_pcs, paste(sprintf("%.1f%%", 100 * ve[seq_len(CFG$comp_pcs)]), collapse = ", "),
    100 * sum(ve[seq_len(CFG$comp_pcs)]))

cpc <- as.data.table(pca$x[, seq_len(CFG$comp_pcs), drop = FALSE])
setnames(cpc, paste0("cPC", seq_len(CFG$comp_pcs)))
cpc[, sample := rownames(pca$x)]
fwrite(cbind(data.table(sample = rownames(pca$x)),
             as.data.table(pca$x[, seq_len(min(5, ncol(pca$x))), drop = FALSE]),
             data.table(var_expl_PC1 = ve[1])),
       file.path(tab, "CELL_composition_pca.csv"))

# STEP 5: composition-adjusted differential expression, within sex
hdr("STEP 5  COMPOSITION-ADJUSTED DIFFERENTIAL EXPRESSION")
md_all <- merge(train_meta[, .(sample, group, sex)], cpc, by = "sample")
setkey(md_all, sample)

fit_sex <- function(sx, adjust) {
  d <- md_all[sex == sx]
  cols <- d$sample
  E <- train_expr[, cols, drop = FALSE]
  g <- factor(d$group, levels = c("HC", "RA"))
  design <- if (adjust)
    stats::model.matrix(stats::as.formula(
      paste("~ g +", paste0("cPC", seq_len(CFG$comp_pcs), collapse = " + "))), data = d)
  else stats::model.matrix(~ g)
  aw  <- limma::arrayWeights(E, design)
  fit <- limma::eBayes(limma::lmFit(E, design, weights = aw))
  tt  <- limma::topTable(fit, coef = "gRA", number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, P = tt$P.Value, FDR = tt$adj.P.Val,
             resid_df = fit$df.residual[1])
}

res <- list()
for (sx in c("F", "M")) {
  lab <- ifelse(sx == "F", "Female", "Male")
  un <- fit_sex(sx, adjust = FALSE)
  ad <- fit_sex(sx, adjust = TRUE)
  m <- merge(un[, .(gene, logFC_unadj = logFC, P_unadj = P, FDR_unadj = FDR)],
             ad[, .(gene, logFC_adj = logFC, P_adj = P, FDR_adj = FDR)], by = "gene")
  sig_u <- m$FDR_unadj < CFG$fdr & abs(m$logFC_unadj) > CFG$lfc
  sig_a <- m$FDR_adj   < CFG$fdr & abs(m$logFC_adj)   > CFG$lfc
  m[, sig_unadjusted := sig_u][, sig_adjusted := sig_a]
  m[, status := fifelse(sig_u & sig_a, "retained",
                fifelse(sig_u & !sig_a, "lost to composition",
                fifelse(!sig_u & sig_a, "gained after adjustment", "ns in both")))]
  setorder(m, P_adj)
  fwrite(m, file.path(tab, sprintf("DEG_celladjusted_%s.csv", tolower(lab))))
  res[[lab]] <- list(m = m, un = un, ad = ad, lab = lab,
                     n_u = sum(sig_u), n_a = sum(sig_a),
                     retained = sum(sig_u & sig_a),
                     lost = sum(sig_u & !sig_a),
                     rho = stats::cor(m$logFC_unadj, m$logFC_adj, method = "spearman"),
                     df_u = un$resid_df[1], df_a = ad$resid_df[1])
  say("%s: DEGs %d (unadjusted) -> %d (composition-adjusted) | retained %d, lost %d",
      lab, sum(sig_u), sum(sig_a), sum(sig_u & sig_a), sum(sig_u & !sig_a))
  say("      logFC Spearman rho unadjusted vs adjusted = %.3f | resid df %d -> %d",
      res[[lab]]$rho, un$resid_df[1], ad$resid_df[1])
}

summ <- rbindlist(lapply(res, function(r) data.table(
  sex = r$lab, n_DEG_unadjusted = r$n_u, n_DEG_adjusted = r$n_a,
  pct_retained = round(100 * r$retained / max(1, r$n_u), 1),
  n_lost_to_composition = r$lost,
  logFC_spearman = round(r$rho, 3),
  resid_df_unadjusted = r$df_u, resid_df_adjusted = r$df_a,
  comp_PCs = CFG$comp_pcs)))
fwrite(summ, file.path(tab, "DEG_celladjusted_summary.csv"))
print(summ)

# STEP 6: what happens to the panel genes under composition adjustment
hdr("STEP 6  PANEL GENES UNDER COMPOSITION ADJUSTMENT")
panels <- list(Female = character(0), Male = character(0))
ml_path <- file.path(procN, "ml_features.rds")
if (file.exists(ml_path)) {
  ml <- readRDS(ml_path)
  panels$Female <- ml$female$consensus
  panels$Male   <- ml$male$consensus
}
pg <- rbindlist(lapply(names(panels), function(sx) {
  gs <- panels[[sx]]
  if (!length(gs)) return(NULL)
  m <- res[[sx]]$m[gene %in% gs]
  m[, sex := sx]
  m[, .(sex, gene, logFC_unadj = round(logFC_unadj, 3), FDR_unadj = signif(FDR_unadj, 3),
        logFC_adj = round(logFC_adj, 3), FDR_adj = signif(FDR_adj, 3),
        pct_logFC_retained = round(100 * logFC_adj / logFC_unadj, 1), status)]
}))
if (!is.null(pg) && nrow(pg)) {
  fwrite(pg, file.path(tab, "CELL_panel_gene_adjustment.csv"))
  print(pg)
} else say("no panel genes available (ml_features.rds missing)")

saveRDS(list(cibersort = cib, mcpcounter = mcp, frac_cols = fc,
             comp_pca = pca, comp_pcs = cpc, tests = tests,
             de = res, summary = summ, panel = pg, config = CFG),
        file.path(procN, "cell_fractions.rds"))

say("")
say("READING RULE FOR THE THESIS")
say("  Unadjusted model = TOTAL RA signal -> the correct basis for a DIAGNOSTIC claim.")
say("  Adjusted model   = signal NOT explained by leukocyte composition -> the")
say("                     only basis on which a CELL-INTRINSIC or mechanistic")
say("                     claim may be made for a gene.")
say("  Composition is plausibly a MEDIATOR (RA -> neutrophilia -> transcriptome),")
say("  so the adjusted estimate is deliberately conservative and is NOT the")
say("  'unconfounded' estimate. Report both; never silently substitute one.")
say("")
say("Wrote CELL_fractions_*.csv, CELL_fraction_group_tests.csv,")
say("CELL_composition_pca.csv, DEG_celladjusted_*.csv,")
say("CELL_panel_gene_adjustment.csv and cell_fractions.rds")
cat("\nDONE\n")
