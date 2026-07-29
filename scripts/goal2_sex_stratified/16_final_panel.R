#!/usr/bin/env Rscript
# =============================================================================
# 18f_final_panel.R
# -----------------------------------------------------------------------------
# FINAL recommended MR-anchored diagnostic models, evaluated honestly across all
# three datasets, and saved for the figure (18f_figure.R).
#
#   Female : elastic-net over the MR-CAUSAL BH-FDR<0.05 set (32 genes,
#            FS_input_female). The p<0.10 widening considered in 18e was
#            REJECTED to preserve stringency.
#   Male   : elastic-net over the 25 MR-causal FDR<0.05 genes (FS_input_male).
# [header corrected 2026-07-28. Prior versions said 14/40 then 74/55, and called
#  these "MR-screened p<0.05". BOTH the counts and the label are now wrong:
#    (a) 10_MR.R fix D6 made FS_input the BH-FDR<0.05 SURVIVING set, so these are
#        FDR-adjusted, not a nominal screen. "MR-causal" is now correct.
#    (b) 10_MR.R fix D8 added the missing CIS filter (996/4,932 instruments were
#        trans; HNRNPM OR=287 and FOXP3 OR=264 were pleiotropic by construction),
#        which re-ran the MR and changed the counts again.
#  Current: 32 female / 25 male, from 1,477 / 1,478 genes tested.]
#
# For each sex: (1) leakage-free NESTED CV on train (10x5 F / 5x10 M);
#               (2) LOCKED elastic-net on full train applied ONCE to internal
#                   holdout and external blood (per-dataset z-score transfer,
#                   fixed direction "<", missing gene -> z=0). DeLong CI, or
#                   bootstrap CI when n<20.
# No estimate gaming (no auto-direction, no scoring train with its own model,
# no seed-hunting).
#
# ---- EVIDENCE TIERS (added 2026-07-28) --------------------------------------
# The male arm is NOT a co-equal result and this table must not present it as
# one. Its sample sizes are:
#       train n = 38 (17 RA / 21 HC) | internal n = 13 (6/7) | external n = 9 (4/5)
# At n = 13 there are 42 case-control pairs, so a perfect AUC is an unremarkable
# event rather than evidence of a perfect classifier, and the bootstrap interval
# collapses to the degenerate "1.000 (1.000-1.000)" that earlier versions of this
# table reported as if it were a validated performance estimate. Two things are
# therefore added here:
#   * an `evidence_tier` column - "primary" or "EXPLORATORY (underpowered)" -
#     assigned from n, not from how good the number looks;
#   * an explicit `separation_flag` wherever AUC >= 0.999, so a degenerate
#     interval can never again be read as precision.
# The male panel is reported as a POWER-LIMITED EXPLORATORY analysis. It is not
# a validated diagnostic panel and no diagnostic claim may rest on it.
#
# Outputs:
#   results/tables/mr_final_panel_summary.csv
#   results/tables/mr_final_coefs_bysex.csv
#   data/processed/mr_final_objects.rds
#
#   Simon R, et al. J Natl Cancer Inst 2003;95:14-18.     (small-n validation)
#   Carpenter J, Bithell J. Stat Med 2000;19:1141-1164.   (bootstrap CIs)
# =============================================================================
suppressMessages({library(glmnet); library(pROC); library(caret)
                  library(Biobase); library(data.table)})
options(stringsAsFactors = FALSE)
GLOBAL_SEED <- 1234; ALPHAS <- c(0.1,0.3,0.5,0.7,0.9,1.0)
proc <- "data/processed"; tab <- "results/tables"

o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr
meta <- as.data.table(o$meta)

## candidate sets ------------------------------------------------------------
genesF <- fread(file.path(tab, "FS_input_female.csv"))$gene   # 32 MR-causal, BH-FDR<0.05
genesM <- fread(file.path(tab, "FS_input_male.csv"))$gene     # 25 MR-causal, BH-FDR<0.05
CAND <- list(F = genesF, M = genesM)

## external datasets ---------------------------------------------------------
load_internal <- function() { h <- readRDS(file.path(proc,"internal_val_holdout_processed.rds"))
  list(expr=h$expr, group=factor(h$meta$group,levels=c("HC","RA")), sex=h$meta$sex, label="Internal test") }
load_blood <- function() {
  e <- readRDS("data/raw/GSE15573_raw.rds"); if (is.list(e)) e <- e[[1]]
  x <- exprs(e); if (max(x,na.rm=TRUE)>50) x <- log2(x+1)
  sym <- fData(e)[["Gene symbol"]]; keep <- !is.na(sym)&sym!=""; x <- x[keep,]; sym <- sym[keep]
  rmean <- rowMeans(x); best <- tapply(seq_along(sym), sym, function(ix) ix[which.max(rmean[ix])])
  xg <- x[unlist(best),]; rownames(xg) <- names(best); p <- pData(e)
  grp <- ifelse(grepl("Rheumatoid|RA",p[["status:ch1"]],ignore.case=TRUE),"RA","HC")
  sex <- ifelse(grepl("Female",p[["gender:ch1"]],ignore.case=TRUE),"F","M")
  list(expr=xg, group=factor(grp,levels=c("HC","RA")), sex=sex, label="External blood") }
internal <- load_internal(); blood <- load_blood()

zrows <- function(M) t(apply(M,1,function(v){ s<-sd(v,na.rm=TRUE)
  if(is.na(s)||s==0) rep(0,length(v)) else (v-mean(v,na.rm=TRUE))/s }))
auc_ci <- function(r){ n <- length(r$cases)+length(r$controls)
  ci <- if(n<20){ set.seed(GLOBAL_SEED); suppressWarnings(as.numeric(ci.auc(r,method="bootstrap",boot.n=2000))) }
        else as.numeric(ci.auc(r)); c(auc=as.numeric(auc(r)), ci[c(1,3)]) }
fit_enet <- function(X,y,nfolds=5){ best<-NULL; bcv<-Inf
  for(a in ALPHAS){ cv<-tryCatch(cv.glmnet(X,y,family="binomial",alpha=a,nfolds=nfolds,standardize=TRUE),error=function(e)NULL)
    if(!is.null(cv)&&min(cv$cvm)<bcv){bcv<-min(cv$cvm); best<-list(cv=cv,alpha=a)} }; best }

## nested CV -----------------------------------------------------------------
run_nested <- function(sx, genes, kfold, repeats){
  sexlab <- if(sx=="F") "Female" else "Male"
  genes <- unique(genes[genes %in% rownames(expr)]); cols <- meta$sample[meta$sex==sx]
  X0 <- t(expr[genes,cols,drop=FALSE]); colnames(X0) <- make.names(genes)
  y <- factor(meta$group[match(cols,meta$sample)],levels=c("HC","RA"))
  preds<-list(); rep_auc<-numeric(0); nz<-integer(0)
  for(rp in seq_len(repeats)){ set.seed(2000+rp); folds<-createFolds(y,k=kfold,returnTrain=FALSE); rpp<-data.table()
    for(fi in seq_along(folds)){ te<-folds[[fi]]; tr<-setdiff(seq_along(y),te)
      if(length(unique(y[tr]))<2) next
      b<-fit_enet(X0[tr,,drop=FALSE],y[tr]); if(is.null(b)) next
      p<-as.numeric(predict(b$cv,newx=X0[te,,drop=FALSE],s="lambda.min",type="response"))
      co<-as.numeric(coef(b$cv,s="lambda.min"))[-1]; nz<-c(nz,sum(co!=0))
      rpp<-rbind(rpp,data.table(sample=rownames(X0)[te],prob=p,obs=y[te])) }
    rr<-roc(rpp$obs,rpp$prob,levels=c("HC","RA"),direction="<",quiet=TRUE); rep_auc<-c(rep_auc,as.numeric(auc(rr))); preds[[rp]]<-rpp }
  allp<-rbindlist(preds); agg<-allp[,.(prob=mean(prob),obs=obs[1]),by=sample]
  ro<-roc(agg$obs,agg$prob,levels=c("HC","RA"),direction="<",quiet=TRUE); a<-auc_ci(ro)
  list(sexlab=sexlab, auc=a[1], ci=a, rep_mean=mean(rep_auc), rep_sd=sd(rep_auc),
       roc=ro, med_used=as.integer(median(nz)), n_cand=length(genes)) }

## locked model on full train -> internal + blood ----------------------------
run_locked <- function(sx, genes){
  sexlab <- if(sx=="F") "Female" else "Male"
  genes <- unique(genes[genes %in% rownames(expr)]); tcols <- meta$sample[meta$sex==sx]
  ytr <- factor(meta$group[match(tcols,meta$sample)],levels=c("HC","RA"))
  Ztr <- t(zrows(expr[genes,tcols,drop=FALSE])); colnames(Ztr) <- make.names(genes)
  b <- fit_enet(Ztr,ytr); co <- as.numeric(coef(b$cv,s="lambda.min"))
  coefs <- data.table(sex=sexlab, term=c("(Intercept)",make.names(genes)), coef=round(co,4), alpha=b$alpha)[coef!=0|term=="(Intercept)"]
  score <- function(ds){ sp<-which(ds$sex==sx); present<-genes[genes %in% rownames(ds$expr)]
    y<-factor(ds$group[sp],levels=c("HC","RA")); if(length(sp)<3||length(unique(y))<2) return(NULL)
    Z<-zrows(ds$expr[present,sp,drop=FALSE]); Zdf<-as.data.frame(t(Z)); colnames(Zdf)<-make.names(present)
    for(g in setdiff(make.names(genes),colnames(Zdf))) Zdf[[g]]<-0
    Zmat<-as.matrix(Zdf[,make.names(genes),drop=FALSE])
    p<-as.numeric(predict(b$cv,newx=Zmat,s="lambda.min",type="response"))
    r<-roc(y,p,levels=c("HC","RA"),direction="<",quiet=TRUE)
    list(label=ds$label, ci=auc_ci(r), roc=r, n=length(y), miss=setdiff(genes,present)) }
  list(coefs=coefs, internal=score(internal), blood=score(blood), alpha=b$alpha) }

cat("Final recommended models (nested CV + locked external)...\n\n")
nF<-run_nested("F",CAND$F,10,5); nM<-run_nested("M",CAND$M,5,10)
lF<-run_locked("F",CAND$F);      lM<-run_locked("M",CAND$M)
fwrite(rbindlist(list(lF$coefs,lM$coefs)), file.path(tab,"mr_final_coefs_bysex.csv"))

# An AUC of 1.000 at n = 13 is a statement about how few discordant pairs exist,
# not about classifier quality. Annotate rather than silently print it.
SMALL_N <- 20
fmt <- function(ci, n = NA_integer_) {
  s <- sprintf("%.3f (%.3f-%.3f)", ci[1], ci[2], ci[3])
  if (!is.na(n)) s <- sprintf("%s [n=%d]", s, n)
  if (!is.na(ci[1]) && ci[1] >= 0.999) s <- paste0(s, " SEPARATION")
  s
}
li  <- function(x) if (is.null(x)) NA_character_ else fmt(x$ci, x$n)
lin <- function(x) if (is.null(x)) NA_integer_   else x$n

n_train  <- c(sum(meta$sex == "F"), sum(meta$sex == "M"))
n_int    <- c(lin(lF$internal), lin(lM$internal))
n_ext    <- c(lin(lF$blood),    lin(lM$blood))
# tier from sample size alone - never from how strong the estimate looks
tier <- ifelse(n_train < 50 | n_int < SMALL_N | n_ext < SMALL_N,
               "EXPLORATORY (underpowered)", "primary")

summ <- data.table(
  sex=c("Female","Male"),
  evidence_tier=tier,
  candidate=c(sprintf("MR-causal BH-FDR<0.05 (%d)",nF$n_cand),
              sprintf("MR-causal BH-FDR<0.05 (%d)",nM$n_cand)),
  n_train=n_train, n_internal=n_int, n_external=n_ext,
  nested_CV=c(fmt(nF$ci, n_train[1]), fmt(nM$ci, n_train[2])),
  nested_per_repeat=c(sprintf("%.3f +/- %.3f",nF$rep_mean,nF$rep_sd), sprintf("%.3f +/- %.3f",nM$rep_mean,nM$rep_sd)),
  internal_test=c(li(lF$internal), li(lM$internal)),
  external_blood=c(li(lF$blood), li(lM$blood)),
  median_genes_used=c(nF$med_used, nM$med_used),
  interpretation=ifelse(tier == "primary",
    "Reportable as a diagnostic performance estimate.",
    paste("POWER-LIMITED. Not a validated panel. Any AUC marked SEPARATION is a",
          "perfect-separation artefact of small n, not a performance estimate.")))
fwrite(summ, file.path(tab,"mr_final_panel_summary.csv"))

rc <- function(sexlab,dataset,r,a) data.table(sex=sexlab,dataset=dataset,sens=r$sensitivities,spec=r$specificities,auc=a)
roc_coords <- rbindlist(list(
  rc("Female","Train (nested CV)",nF$roc,nF$auc), rc("Female","Internal test",lF$internal$roc,lF$internal$ci[1]),
  rc("Female","External blood",lF$blood$roc,lF$blood$ci[1]),
  rc("Male","Train (nested CV)",nM$roc,nM$auc), rc("Male","Internal test",lM$internal$roc,lM$internal$ci[1]),
  rc("Male","External blood",lM$blood$roc,lM$blood$ci[1])))
saveRDS(list(summary=summ, roc=roc_coords, nested=list(F=nF,M=nM), locked=list(F=lF,M=lM)),
        file.path("data/processed/new","mr_final_objects.rds"))

cat("\n==================  FINAL RECOMMENDED PANEL  ==================\n"); print(summ, width=200)
cat("\nFemale locked coefficients:\n"); print(lF$coefs)

cat("\n------------------  EVIDENCE TIER WARNINGS  ------------------\n")
for (i in seq_len(nrow(summ))) {
  if (summ$evidence_tier[i] == "primary") next
  cat(sprintf("  %s: EXPLORATORY - train n=%d, internal n=%s, external n=%s.\n",
              summ$sex[i], summ$n_train[i],
              ifelse(is.na(summ$n_internal[i]), "NA", summ$n_internal[i]),
              ifelse(is.na(summ$n_external[i]), "NA", summ$n_external[i])))
  cat("     Report as a power-limited exploratory analysis. Do NOT present its\n")
  cat("     AUCs alongside the female estimates as if they were comparable, and\n")
  cat("     do NOT quote any interval marked SEPARATION as a performance figure.\n")
}
cat("\nSaved mr_final_panel_summary.csv, mr_final_coefs_bysex.csv, mr_final_objects.rds\n")
