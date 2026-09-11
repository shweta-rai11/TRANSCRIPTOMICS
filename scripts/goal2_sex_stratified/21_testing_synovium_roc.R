#!/usr/bin/env Rscript
# =============================================================================
# new/24new_synovium_roc.R  -- ROC curves in the SYNOVIUM external test
# (GSE89408, RA vs Normal), per sex, for the consensus panels.
#   (1) per-gene univariate ROC overlay. ORIENTATION = BEST-DIRECTION
#       (direction="auto", so AUC>=0.5 by construction). Direction reversal vs
#       blood is NOT visible here; it is carried by the `concordant` flag from
#       20_. The train-oriented version is built in 24_. Do not read this
#       overlay as evidence of cross-tissue transfer.
#   (2) panel ROC: APPARENT (resubstitution) ONLY. The model is refitted on the
#       synovial samples and scored on those same samples. There is NO 10-fold
#       CV in this script (an earlier header claimed one) and no CI. Only the
#       GENE SET transfers from blood, never the coefficients. See thesis 2.11.
# x-axis = 1 - Specificity (standard). Uses data/processed/new/val_synovium.rds.
# Output: results/figures/new/fig_syn_pergene_roc_{female,male}.png/pdf
#         results/figures/new/fig_syn_panel_roc.png/pdf
# =============================================================================
suppressMessages({library(data.table); library(pROC); library(ggplot2)})
proc <- "data/processed"; figN <- "results/figures/new"
v  <- readRDS(file.path(proc, "new", "val_synovium.rds"))
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
logcpm <- v$logcpm; grp <- v$grp; sex <- v$sex
panels <- list(Female = ml$female$consensus, Male = ml$male$consensus)
sc <- c(Female = "F", Male = "M")

# ---- (1) per-gene ROC overlay -----------------------------------------------
pergene <- function(sexlab) {
  genes <- panels[[sexlab]]; genes <- genes[genes %in% rownames(logcpm)]
  idx <- which(sex == sc[[sexlab]]); y <- factor(grp[idx], levels = c("Normal","RA"))
  crd <- list(); au <- numeric(0)
  for (g in genes) {
    r <- roc(y, as.numeric(logcpm[g, idx]), direction = "auto", levels = c("Normal","RA"), quiet = TRUE)
    au[g] <- as.numeric(auc(r))
    crd[[g]] <- data.table(gene = g, sens = r$sensitivities, spec = r$specificities, auc = au[g],
                           lab = sprintf("%s (AUC=%.3f)", g, au[g]))
  }
  d <- rbindlist(crd); ord <- names(sort(au, decreasing = TRUE))
  d[, lab := factor(lab, levels = sprintf("%s (AUC=%.3f)", ord, au[ord]))]
  n <- length(genes); pal <- grDevices::hcl.colors(n, "Dark 3")
  g <- ggplot(d, aes(1 - spec, sens, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey75") +
    geom_path(linewidth = 0.8) +
    scale_x_continuous(limits = c(0,1), expand = c(0.01,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0.01,0)) +
    scale_colour_manual(values = pal, name = NULL) + coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity",
         subtitle = sprintf("%s panel in synovium (RA vs Normal, n=%d)", sexlab, length(idx))) +
    theme_bw(base_size = 12) + theme(legend.position = "right",
      legend.text = element_text(size = 8), panel.grid.minor = element_blank())
  ggsave(file.path(figN, sprintf("fig_syn_pergene_roc_%s.png", tolower(sexlab))), g, width = 7.5, height = 5.5, dpi = 300)
  ggsave(file.path(figN, sprintf("fig_syn_pergene_roc_%s.pdf", tolower(sexlab))), g, width = 7.5, height = 5.5)
  invisible(au)
}
for (s in c("Female","Male")) pergene(s)

# ---- (2) panel ROC (combined logistic) --------------------------------------
panel_roc <- function(sexlab) {
  genes <- panels[[sexlab]]; genes <- genes[genes %in% rownames(logcpm)]
  idx <- which(sex == sc[[sexlab]]); y <- factor(grp[idx], levels = c("Normal","RA"))
  Z <- scale(t(logcpm[genes, idx, drop=FALSE])); Z[is.na(Z)] <- 0; colnames(Z) <- make.names(genes)
  fit <- suppressWarnings(glm(y ~ ., data.frame(y=y, Z, check.names=FALSE), family=binomial))
  r <- roc(y, as.numeric(predict(fit, type="response")), levels=c("Normal","RA"), direction="<", quiet=TRUE)
  data.table(sex = sexlab, sens = r$sensitivities, spec = r$specificities,
             lab = sprintf("%s panel (AUC=%.3f)", sexlab, as.numeric(auc(r))))
}
pr <- rbindlist(lapply(c("Female","Male"), panel_roc))
gp <- ggplot(pr, aes(1 - spec, sens, colour = lab)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey75") +
  geom_path(linewidth=1.1) +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0.01,0)) +
  scale_colour_manual(values=c("#C0392B","#1F3B99"), name=NULL) + coord_equal() +
  labs(x="1 - Specificity", y="Sensitivity", subtitle="Panel ROC in synovium (RA vs Normal), apparent") +
  theme_bw(base_size=12) + theme(legend.position=c(0.98,0.02), legend.justification=c(1,0),
    panel.grid.minor=element_blank())
ggsave(file.path(figN, "fig_syn_panel_roc.png"), gp, width=6, height=6, dpi=300)
ggsave(file.path(figN, "fig_syn_panel_roc.pdf"), gp, width=6, height=6)
cat("wrote fig_syn_pergene_roc_{female,male} + fig_syn_panel_roc\n")
