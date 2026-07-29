#!/usr/bin/env Rscript
# =============================================================================
# 18k_panelBC_logfc_heatmaps.R  -- reference Fig 4B/4C analogs for RA data.
#   Panel B: logFC (RA - HC, log2) of the 6 MR consensus genes in each DISCOVERY
#            dataset (GSE93272, GSE110169), split by sex. (= their per-GEO 4B)
#   Panel C: the same logFC in the INDEPENDENT VALIDATION cohorts (internal
#            holdout, external blood GSE15573), split by sex. This is the honest
#            RA analog of their ADEx external-database panel 4C -- RA has no ADEx
#            equivalent, so we use the independent cohorts and label it as such.
# Genes annotated female-specific (red) / male-specific (blue) as in the paper.
# Outputs: fig_mr_panelB_logfc_discovery.{png,pdf}, fig_mr_panelC_logfc_validation.{png,pdf}
#          results/tables/mr_panelBC_logfc.csv
# =============================================================================
suppressMessages({library(Biobase); library(data.table); library(ggplot2)})
proc <- "data/processed"; tab <- "results/tables"; fig <- "results/figures/new"

ml <- readRDS("data/processed/new/ml_features.rds")
genesF <- ml$female$consensus; genesM <- ml$male$consensus
genes  <- c(genesF, genesM)
# Labels corrected 2026-07-27: were "Female-specific" / "Male-specific", which is
# a sex-SPECIFICITY claim this chapter explicitly does not make (see
# results/README_GOALS.md, "Terminology"). A gene appears in only one panel
# because it passed selection in that stratum, not because it is specific to it.
# The formal interaction test (Q4 04_parent_panel_sexspecificity.R) finds 0 of 11
# panel genes sex-specific, minimum FDR 0.238.
gene_type <- setNames(c(rep("Female panel", length(genesF)),
                        rep("Male panel", length(genesM))), genes)

# logFC (RA - HC) per gene, within a given expr matrix + meta subset -----------
logfc <- function(expr, samples, grp, sx, cohort) {
  rbindlist(lapply(c("F", "M"), function(s) {
    idx <- which(sx == s); if (!length(idx)) return(NULL)
    g <- grp[idx]; sm <- samples[idx]
    ra <- sm[g == "RA"]; hc <- sm[g == "HC"]
    if (length(ra) < 2 || length(hc) < 2) return(NULL)
    present <- genes[genes %in% rownames(expr)]
    lfc <- rowMeans(expr[present, ra, drop = FALSE]) - rowMeans(expr[present, hc, drop = FALSE])
    data.table(gene = present, logFC = as.numeric(lfc),
               cohort = cohort, sex = ifelse(s == "F", "Female", "Male"),
               n_RA = length(ra), n_HC = length(hc))
  }))
}

## discovery: split combined_train by dataset --------------------------------
o <- readRDS(file.path(proc, "combined_train.rds")); me <- as.data.table(o$meta)
disc <- rbindlist(lapply(unique(me$dataset), function(ds) {
  ii <- me$dataset == ds
  logfc(o$expr, me$sample[ii], me$group[ii], me$sex[ii], ds)
}))
disc[, panel := "B: Discovery datasets"]

## validation: internal holdout + external blood -----------------------------
h <- readRDS(file.path(proc, "internal_val_holdout_processed.rds")); mh <- as.data.table(h$meta)
val_int <- logfc(h$expr, mh$sample, mh$group, mh$sex, "Internal test")

e <- readRDS("data/raw/GSE15573_raw.rds"); if (is.list(e)) e <- e[[1]]
x <- exprs(e); if (max(x, na.rm = TRUE) > 50) x <- log2(x + 1)
sym <- fData(e)[["Gene symbol"]]; keep <- !is.na(sym) & sym != ""; x <- x[keep, ]; sym <- sym[keep]
rmean <- rowMeans(x); best <- tapply(seq_along(sym), sym, function(ix) ix[which.max(rmean[ix])])
xg <- x[unlist(best), ]; rownames(xg) <- names(best); p <- pData(e)
bgrp <- ifelse(grepl("Rheumatoid|RA", p[["status:ch1"]], ignore.case = TRUE), "RA", "HC")
bsex <- ifelse(grepl("Female", p[["gender:ch1"]], ignore.case = TRUE), "F", "M")
val_bld <- logfc(xg, colnames(xg), bgrp, bsex, "External blood")
val <- rbindlist(list(val_int, val_bld)); val[, panel := "C: Validation cohorts"]

all <- rbindlist(list(disc, val), fill = TRUE)
all[, gene := factor(gene, levels = rev(genes))]
all[, col := paste(cohort, sex, sep = "\n")]
all[, type := gene_type[as.character(gene)]]
fwrite(all, file.path(tab, "mr_panelBC_logfc.csv"))

lim <- max(abs(all$logFC), na.rm = TRUE)
heat <- function(d, subtitle) {
  ggplot(d, aes(x = col, y = gene, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", logFC)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim), name = "logFC\n(RA - HC)") +
    facet_grid(type ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL, subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 8), panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0, face = "bold"),
          axis.text.y = element_text(face = "bold"))
}
gB <- heat(all[panel == "B: Discovery datasets"], "Discovery datasets (per GEO series)")
gC <- heat(all[panel == "C: Validation cohorts"], "Independent validation cohorts")
ggsave(file.path(fig, "fig_mr_panelB_logfc_discovery.png"), gB, width = 6.5, height = 5, dpi = 300)
ggsave(file.path(fig, "fig_mr_panelB_logfc_discovery.pdf"), gB, width = 6.5, height = 5)
ggsave(file.path(fig, "fig_mr_panelC_logfc_validation.png"), gC, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig, "fig_mr_panelC_logfc_validation.pdf"), gC, width = 6, height = 5)
cat("wrote panel B and C heatmaps\n"); print(all[, .(cohort, sex, gene, logFC = round(logFC, 2))])
