#!/usr/bin/env Rscript
# =============================================================================
# 20_testing_synovium_external.R  —  CROSS-TISSUE validation of the blood sex-stratified
# ML signatures in RA SYNOVIUM (GSE89408 RNA-seq, n=218). RA (152) vs Normal
# (28); sex assigned by position from the GSM pData (146 F / 72 M). Counts are
# TMM-normalized (edgeR) and log2-CPM transformed. For each signature gene we
# compute synovium RA-vs-Normal log2FC (limma-voom, sex-adjusted), significance,
# and AUC — overall and sex-stratified — plus direction concordance vs the blood
# training DE. Output: results/tables/val_synovium_*.csv + val_synovium.rds
# =============================================================================
suppressMessages({library(edgeR); library(limma); library(Biobase); library(pROC); library(data.table)})
options(stringsAsFactors = FALSE)
proc <- "data/processed"; tab <- "results/tables"; dir.create(tab,showWarnings=FALSE,recursive=TRUE)

ml <- readRDS("data/processed/new/ml_features.rds")
D  <- readRDS(file.path(proc, "dge_results.rds"))
fsig <- ml$female$consensus; msig <- ml$male$consensus; allsig <- union(fsig, msig)

# ---- load counts + map disease (col prefix) & sex (pData position) ----------
cnt <- as.matrix(read.delim(gzfile("data/raw/GSE89408_counts.txt.gz"), row.names = 1, check.names = FALSE))
pref <- gsub("_[0-9]+$", "", colnames(cnt))
dis  <- c(normal_tissue="Normal", RA_tissue="RA", OA_tissue="OA",
          AG_tissue="AG", undiff_tissue="UA")[pref]
e <- readRDS("data/raw/GSE89408_raw.rds"); sx <- pData(e)[["Sex:ch1"]]   # position-aligned (verified 218/218)
stopifnot(length(sx) == ncol(cnt))

keep_s <- dis %in% c("RA", "Normal")
cnt <- cnt[, keep_s]; grp <- factor(dis[keep_s], levels = c("Normal","RA")); sex <- sx[keep_s]
cat(sprintf("GSE89408 synovium: RA=%d Normal=%d | F=%d M=%d\n",
            sum(grp=="RA"), sum(grp=="Normal"), sum(sex=="F"), sum(sex=="M")))
cat(sprintf("  RA:  F=%d M=%d | Normal: F=%d M=%d\n",
            sum(grp=="RA"&sex=="F"), sum(grp=="RA"&sex=="M"),
            sum(grp=="Normal"&sex=="F"), sum(grp=="Normal"&sex=="M")))

# ---- TMM + voom, sex-adjusted RA-vs-Normal DE ------------------------------
dge <- DGEList(cnt); keepg <- filterByExpr(dge, group = grp)
dge <- dge[keepg, , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
design <- model.matrix(~ sex + grp)                       # adjust for sex; coef grpRA = RA effect
v <- voom(dge, design); fit <- eBayes(lmFit(v, design))
tt <- topTable(fit, coef = "grpRA", number = Inf, sort.by = "none")
tt$gene <- rownames(tt); setDT(tt)

train_dir <- function(s){ d <- D$res[[ifelse(s=="F","Female","Male")]]; setNames(sign(d$logFC), d$gene) }

# -----------------------------------------------------------------------------
# PER-GENE SYNOVIUM AUC — TWO ORIENTATION CONVENTIONS, BOTH EMITTED.
#
# This function previously returned ONLY the best-direction AUC (it flipped the
# ROC direction whenever AUC < 0.5), so every stored value was >= 0.5 by
# construction. That silently contradicted the methods text, which promised that
# a gene reversing direction out of sample would show as AUC < 0.5. The table is
# now explicit and carries both, with self-describing column names:
#
#   *_bestdir   BEST-DIRECTION. Orientation chosen inside the synovial data.
#               >= 0.5 by construction. Answers "how much information does this
#               gene carry in synovium, irrespective of direction?".
#               NOT a transfer result. Never place beside blood AUCs.
#
#   *_trainorient  TRAIN-FIXED. Orientation fixed by the gene's blood training
#               direction for that sex. AUC < 0.5 means the association REVERSED
#               between blood and synovium. This is the ONLY convention valid for
#               cross-dataset comparison, and equals bestdir when concordant,
#               1 - bestdir when discordant (the transform applied in 24_).
#
# See thesis 2.9 "Orientation conventions" and 2.11.
# -----------------------------------------------------------------------------
gene_stat <- function(g, samples, y) {
  if (!g %in% rownames(logcpm)) return(NA_real_)
  vv <- as.numeric(logcpm[g, samples])
  r <- roc(y, vv, direction="<", levels=c("Normal","RA"), quiet=TRUE)
  a <- as.numeric(auc(r))
  if (a < 0.5) a <- 1 - a          # best-direction value
  a
}
val_syn <- function(sig, s) {
  td <- train_dir(s)
  ov <- seq_along(grp); sxi <- which(sex==s)
  rbindlist(lapply(sig, function(g) {
    row <- tt[gene==g]
    if (!nrow(row)) return(data.table(gene=g, present=FALSE))
    conc <- sign(row$logFC) == unname(td[g])
    a_all <- gene_stat(g, ov,  grp)
    a_sex <- gene_stat(g, sxi, droplevels(grp[sxi]))
    data.table(gene=g, present=TRUE, syn_log2FC=row$logFC, syn_adjP=row$adj.P.Val,
               syn_dir=sign(row$logFC), train_dir=unname(td[g]),
               concordant=conc,
               auc_all_bestdir = a_all,
               auc_sex_bestdir = a_sex,
               auc_all_trainorient = ifelse(conc, a_all, 1 - a_all),
               auc_sex_trainorient = ifelse(conc, a_sex, 1 - a_sex),
               # legacy aliases so downstream scripts keep working unchanged
               auc_all = a_all, auc_sex = a_sex)
  }), fill=TRUE)
}
sf <- val_syn(fsig, "F"); sm <- val_syn(msig, "M")

for (s in c("female","male")) fwrite(if(s=="female") sf else sm,
     file.path(tab, sprintf("val_synovium_pergene_%s.csv", s)))

summ <- function(x, lab) cat(sprintf(
  "== %s synovium: sig&concordant %d/%d | concordant %d/%d | median AUC(all)=%.3f\n",
  lab, sum(x$concordant & x$syn_adjP<0.05, na.rm=TRUE), nrow(x[present==TRUE]),
  sum(x$concordant, na.rm=TRUE), nrow(x[present==TRUE]), median(x$auc_all, na.rm=TRUE)))
summ(sf,"FEMALE"); summ(sm,"MALE")
cat("\nFEMALE (synovium RA vs Normal):\n"); print(sf[present==TRUE][order(-auc_all),
   .(gene, syn_log2FC=round(syn_log2FC,2), syn_adjP=signif(syn_adjP,2), concordant, auc_all=round(auc_all,3))])
cat("\nMALE (synovium RA vs Normal):\n"); print(sm[present==TRUE][order(-auc_all),
   .(gene, syn_log2FC=round(syn_log2FC,2), syn_adjP=signif(syn_adjP,2), concordant, auc_all=round(auc_all,3))])

saveRDS(list(logcpm=logcpm, grp=grp, sex=sex, tt=tt, sf=sf, sm=sm,
             fsig=fsig, msig=msig), file.path("data/processed/new", "val_synovium.rds"))
cat("\nSaved val_synovium.rds + per-gene tables\n")
