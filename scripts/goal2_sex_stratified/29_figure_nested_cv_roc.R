#!/usr/bin/env Rscript
# Nested-CV figures: (A) leaky flat-CV vs leakage-free nested AUC, (B) nested-CV pooled ROC, (C) male signature-gene re-selection stability across folds.
suppressMessages({library(ggplot2); library(data.table); library(pROC)})
theme_set(theme_bw(base_size = 12))
proc <- "data/processed"; fig <- "results/figures"
N <- readRDS(file.path(proc, "nested_cv.rds"))
save2 <- function(g,n,w,h){ ggsave(file.path(fig,paste0(n,".png")),g,width=w,height=h,dpi=300)
  ggsave(file.path(fig,paste0(n,".pdf")),g,width=w,height=h,device=cairo_pdf) }
sexcol <- c(Female = "#C0392B", Male = "#2C5F8A")

# ---- A: leaky vs nested AUC -------------------------------------------------
au <- rbind(
  data.table(sex=c("Female","Male"), type="Flat CV (leaky)",
             auc=c(0.947,1.000), lo=NA, hi=NA),
  data.table(sex=c("Female","Male"), type="Nested CV (leakage-free)",
             auc=c(N$female$auc, N$male$auc),
             lo=c(N$female$ci[1], N$male$ci[1]), hi=c(N$female$ci[3], N$male$ci[3])))
au[, type := factor(type, levels=c("Flat CV (leaky)","Nested CV (leakage-free)"))]
pA <- ggplot(au, aes(type, auc, fill=sex)) +
  geom_col(position=position_dodge(0.8), width=0.7, alpha=0.9) +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(0.8), width=0.2, na.rm=TRUE) +
  geom_text(aes(label=sprintf("%.3f", auc)), position=position_dodge(0.8), vjust=-0.5, size=3.4) +
  facet_wrap(~sex) +
  scale_fill_manual(values=sexcol, guide="none") +
  coord_cartesian(ylim=c(0.5,1.05)) +
  labs(title="Leakage-free nested CV corrects the inflated flat-CV AUC",
       subtitle="Feature selection (in-fold limma DEG -> LASSO) redone inside every outer fold; nested bars show 95% CI",
       x=NULL, y="AUC") +
  theme(plot.subtitle=element_text(size=9), axis.text.x=element_text(size=9))
save2(pA, "fig_nested_auc_compare", 9, 5.5)

# ---- B: nested pooled ROC ---------------------------------------------------
roc_dt <- function(res, lab) data.table(sex=lab, fpr=1-res$roc$specificities,
                                        sens=res$roc$sensitivities, auc=res$auc)
rc <- rbind(roc_dt(N$female,"Female"), roc_dt(N$male,"Male"))
setorder(rc, sex, fpr, sens)
labs <- rc[, .(auc=auc[1]), by=sex][, lab:=sprintf("%s  AUC %.3f", sex, auc)]
rc <- merge(rc, labs[,.(sex,lab)], by="sex")
pB <- ggplot(rc, aes(fpr, sens, color=sex)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey65") +
  geom_step(direction="hv", linewidth=1) +
  scale_color_manual(values=sexcol, name=NULL,
                     labels=setNames(labs$lab, labs$sex)) +
  coord_equal() +
  labs(title="Nested-CV pooled ROC (leakage-free)", x="1 - specificity", y="Sensitivity") +
  theme(legend.position=c(0.98,0.02), legend.justification=c(1,0))
save2(pB, "fig_nested_roc", 6.5, 6)

# ---- C: male signature stability -------------------------------------------
sm <- as.data.table(N$stab_m); setorder(sm, reselect_pct)
sm[, gene := factor(gene, levels=gene)]
pC <- ggplot(sm, aes(reselect_pct, gene)) +
  geom_col(fill="#2C5F8A", width=0.7) +
  geom_text(aes(label=sprintf("%.0f%%", reselect_pct)), hjust=-0.15, size=3.4) +
  coord_cartesian(xlim=c(0,100)) +
  labs(title="Male signature-gene re-selection frequency (nested folds)",
       subtitle="How often each original male signature gene is re-picked when selection is redone leakage-free",
       x="Re-selected in % of outer folds (n=100)", y=NULL) +
  theme(plot.subtitle=element_text(size=9))
save2(pC, "fig_nested_male_stability", 8, 5)

cat("Wrote: fig_nested_auc_compare, fig_nested_roc, fig_nested_male_stability\n")
