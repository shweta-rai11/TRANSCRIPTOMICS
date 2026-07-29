#!/usr/bin/env Rscript
# =============================================================================
# eda_figure.R  —  EDA FIGURES (plotting only) for GSE93272 + GSE110169.
# Reads data/processed/eda_results.rds (from eda.R) and writes PNG + PDF to
# results/figures/. Focus: sex x disease per dataset, metadata distributions,
# and missing-value audit (all as bar plots), plus expression QC.
# =============================================================================
suppressMessages({library(ggplot2); library(data.table)})
theme_set(theme_bw(base_size = 12))
proc <- "data/processed"; fig <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
res <- readRDS(file.path(proc, "eda_results.rds"))

save2 <- function(g, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), g, width = w, height = h, dpi = 300)
  ggsave(file.path(fig, paste0(name, ".pdf")), g, width = w, height = h, device = cairo_pdf)
}
sexcol <- c(F = "#D81B60", M = "#1E88E5")

# ---- 1. SEX x DISEASE per dataset (RA vs control) --------------------------
sx <- rbindlist(lapply(res, `[[`, "sxdis"))
sx <- sx[group %in% c("RA","HC") & sex %in% c("F","M")]
sx$group <- factor(sx$group, levels = c("HC","RA"), labels = c("Healthy control","RA"))
p1 <- ggplot(sx, aes(group, Freq, fill = sex)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_text(aes(label = Freq), position = position_dodge(.8), vjust = -0.3, size = 4, fontface = "bold") +
  facet_wrap(~dataset, scales = "free_x") +
  scale_fill_manual(values = sexcol, name = "Sex", labels = c("Female","Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "Sex distribution in RA vs healthy control, per dataset",
       subtitle = "GSE93272 = unique patients (longitudinal repeats collapsed); GSE110169 = subjects (SLE excluded)",
       x = NULL, y = "Number of subjects") +
  theme(legend.position = "top")
save2(p1, "eda_sex_by_disease", 9, 5)

# ---- 2. full disease distribution (incl. SLE, shown as excluded) -----------
dd <- rbindlist(lapply(names(res), function(nm){ d <- as.data.table(res[[nm]]$disease_tab); d$dataset <- nm; d }))
dd$group <- factor(dd$group, levels = c("HC","RA","SLE","other"))
dd$use <- ifelse(dd$group %in% c("RA","HC"), "used (RA vs control)", "excluded")
p2 <- ggplot(dd, aes(group, Freq, fill = use)) +
  geom_col(width = .7) + geom_text(aes(label = Freq), vjust = -0.3, size = 4) +
  facet_wrap(~dataset, scales = "free") +
  scale_fill_manual(values = c(`used (RA vs control)` = "#2E7D32", excluded = "#BDBDBD"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "Disease-group composition per dataset (SLE excluded from analysis)",
       x = NULL, y = "Number of subjects") + theme(legend.position = "top")
save2(p2, "eda_disease_distribution", 9, 5)

# ---- 3. missing-value audit per dataset ------------------------------------
miss <- rbindlist(lapply(res, `[[`, "miss"))
miss_plot <- miss[pct_missing > 0]
if (nrow(miss_plot)) {
  p3 <- ggplot(miss_plot, aes(reorder(field, pct_missing), pct_missing, fill = dataset)) +
    geom_col() + coord_flip() +
    geom_text(aes(label = paste0(pct_missing, "%")), hjust = -0.1, size = 3) +
    facet_wrap(~dataset, scales = "free_y") +
    scale_fill_manual(values = c(GSE93272 = "#5E35B1", GSE110169 = "#00838F"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = "Missing values per metadata field (fields with >0% missing)",
         x = NULL, y = "% missing") + theme(strip.text = element_text(face = "bold"))
  save2(p3, "eda_missing_values", 11, 7)
} else {
  save2(ggplot() + annotate("text", 0, 0, label = "No missing values detected") + theme_void(),
        "eda_missing_values", 6, 3)
}
# complete-vs-missing overview (all fields)
mc <- copy(miss); mc$complete <- 100 - mc$pct_missing
mcl <- melt(mc[, .(dataset, field, complete, pct_missing)], id.vars = c("dataset","field"),
            variable.name = "status", value.name = "pct")
mcl$status <- factor(mcl$status, levels = c("pct_missing","complete"), labels = c("missing","present"))

# ---- 4. metadata distributions per dataset (categorical bars) --------------
plot_cats <- function(nm) {
  cats <- res[[nm]]$cats
  if (is.null(cats)) return(invisible())
  cats <- as.data.table(cats)
  cats[, value := as.character(value)]
  cats[, value := ifelse(nchar(value) > 22, paste0(substr(value,1,20),".."), value)]
  g <- ggplot(cats, aes(value, Freq)) +
    geom_col(fill = "#3949AB") +
    geom_text(aes(label = Freq), vjust = -0.2, size = 2.6) +
    facet_wrap(~field, scales = "free", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = sprintf("%s — categorical metadata distributions", nm), x = NULL, y = "n samples") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 9))
  n_fields <- length(unique(cats$field))
  save2(g, paste0("eda_metadata_", nm), 12, max(4, 2.2 * ceiling(n_fields/3)))
}
for (nm in names(res)) plot_cats(nm)

# ---- 5. expression QC: per-sample distribution + PCA -----------------------
ps <- rbindlist(lapply(res, function(r) r$qc$persample))
p5 <- ggplot(ps, aes(reorder(sample, median), median)) +
  geom_linerange(aes(ymin = q25, ymax = q75), color = "#90A4AE", linewidth = .3) +
  geom_point(size = .5, color = "#1565C0") +
  facet_wrap(~dataset, scales = "free", ncol = 1) +
  labs(title = "Expression QC — per-sample median with IQR (aligned = normalized)",
       x = "Sample (ranked)", y = "expression (median, IQR)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save2(p5, "eda_expression_distribution", 11, 6)

# RA vs HC only — PCA RECOMPUTED with SLE removed (all-groups version dropped as redundant)
pca_rh <- rbindlist(lapply(res, function(r) r$qc$pca_rh$scores))
pca_rh <- pca_rh[group %in% c("RA","HC")]
p6b <- ggplot(pca_rh, aes(PC1, PC2, color = group)) +
  geom_point(size = 1.8, alpha = .85) + facet_wrap(~dataset, scales = "free") +
  scale_color_manual(values = c(HC = "#43A047", RA = "#FB8C00"),
                     labels = c("Healthy control","RA")) +
  labs(title = "Expression QC — PCA per dataset, RA vs healthy control (SLE removed)",
       subtitle = "PCA recomputed on RA + HC samples only (top-2000 variable probes)",
       x = "PC1", y = "PC2", color = NULL) + theme(legend.position = "top")
save2(p6b, "eda_pca_ra_vs_hc", 10, 5)

cat("EDA figures written to results/figures/ (eda_*.png/.pdf)\n")
