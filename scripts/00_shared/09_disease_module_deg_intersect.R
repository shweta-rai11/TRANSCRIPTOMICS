#!/usr/bin/env Rscript
# =============================================================================
# 09_disease_module_deg_intersect.R      *** THE CANDIDATE STEP ***
# -----------------------------------------------------------------------------
#        candidates  =  disease-module genes   n   sex-stratified DEGs
#
# LOGIC (symmetric, data-driven, no sex-arbitrary module assignment)
#   1. Take the DISEASE modules selected in 06_WGCNA.R from the OVERALL
#      RA-vs-Control module-trait correlation on ALL samples:
#           |cor(ME, RA)| >= dm_min_abs_cor  AND  p < dm_max_p
#      The SAME module set is used for BOTH sexes, so the candidate pool is not
#      biased by which module a sex happens to be assigned to.
#   2. Intersect that disease-gene background with each sex's limma DEGs
#      SEPARATELY:
#           Female candidates = disease modules  n  Female DEGs
#           Male   candidates = disease modules  n  Male   DEGs
#      Sex-stratification therefore enters ONLY through the DEG contrast, on a
#      common, disease-relevant gene background.
#
# WHAT CHANGED (2026-07-27)
#   - Reads data/processed/wgcna_results.rds (06_WGCNA.R). The old input
#     wgcna_new_results.rds came from the deleted 06_WGCNA_NEW.R.
#   - Disease modules are read from the analysis object, never re-derived here
#     and never hardcoded. The previous version printed per-module counts using
#     literal "green"/"brown" tests, which silently returned 0 whenever the
#     network produced different colour names (colours track module SIZE RANK).
#   - Reports how many sex-DEGs were LOST to the WGCNA variance filter before
#     the network was built. Those genes could never become candidates, so the
#     loss is stated rather than hidden.
#
# INTERPRETATION LIMIT (sex-stratified, not sex-specific)
#   A gene in the Female list but not the Male list is NOT shown to behave
#   differently in women. With ~4x fewer males, absence from the male list is
#   expected from power alone. No interaction test is performed anywhere here.
#     # Gelman A, Stern H. Am Stat 2006;60(4):328-331.
#
# Inputs : data/processed/wgcna_results.rds   (gene_tab, disease_modules, kME)
#          data/processed/dge_results.rds     (res$Female / res$Male; sig, dir)
# Outputs: results/tables/candidates_female_disease.csv
#          results/tables/candidates_male_disease.csv
#          results/tables/disease_module_selection.csv
#          results/tables/candidate_summary.csv
#          results/figures/fig_venn_{female,male}_disease_candidates.png
#
# ---- References -------------------------------------------------------------
#   Langfelder P, Horvath S. WGCNA. BMC Bioinformatics 2008;9:559.
#   Ritchie ME, et al. limma. Nucleic Acids Res 2015;43(7):e47.
#   Gelman A, Stern H. Am Stat 2006;60(4):328-331.
# =============================================================================
suppressMessages({
  library(data.table); library(ggVennDiagram); library(ggplot2)
})
options(stringsAsFactors = FALSE)

proc <- "data/processed"; tab <- "results/tables"; fig <- "results/figures"
for (d in c(tab, fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---- load -------------------------------------------------------------------
wf <- file.path(proc, "wgcna_results.rds")
df <- file.path(proc, "dge_results.rds")
if (!file.exists(wf)) stop("Missing ", wf, " - run 06_WGCNA.R first.")
if (!file.exists(df)) stop("Missing ", df, " - run 05_dge.R first.")
w <- readRDS(wf); D <- readRDS(df)

gt       <- as.data.table(w$gene_tab)          # gene, module, GS_RA, kME_own
kME      <- as.data.frame(w$kME)
dis_mods <- w$disease_modules                  # data-driven (06_WGCNA.R STEP 9)
dropped  <- w$genes_dropped                    # removed by the variance filter

disease_genes <- gt[module %in% dis_mods, gene]

cat("\n===== INPUTS =====\n")
say("WGCNA power              : %d", w$soft_power)
say("genes in network         : %d  (variance filter dropped %d)",
    nrow(gt), length(dropped))
say("disease modules          : %s", paste(dis_mods, collapse = " + "))
say("  selection rule         : |cor(ME,RA)| >= %.2f & p < %g",
    w$config$dm_min_abs_cor, w$config$dm_max_p)
say("disease-module genes     : %d", length(disease_genes))
say("DEG cutoff in force      : |log2FC| > %.2f & FDR < 0.05", D$thr$FC)

# per-module detail of the selection, for the supplement
sel_tab <- as.data.table(w$module_trait_table)[module %in% dis_mods]
fwrite(sel_tab, file.path(tab, "disease_module_selection.csv"))
print(sel_tab)

# ---- intersect with each sex's DEGs -----------------------------------------
build <- function(sx) {
  d    <- as.data.table(D$res[[sx]])
  sig  <- d[sig == TRUE, gene]
  hit  <- intersect(disease_genes, sig)
  lost <- intersect(sig, dropped)              # DEGs the variance filter removed

  mod_of <- gt$module[match(hit, gt$gene)]
  kme_of <- vapply(seq_along(hit), function(i) {
    kc <- paste0("kME", mod_of[i])
    if (kc %in% colnames(kME) && hit[i] %in% rownames(kME)) kME[hit[i], kc] else NA_real_
  }, numeric(1))

  out <- data.table(
    sex     = sx,
    gene    = hit,
    module  = mod_of,
    logFC   = d$logFC[match(hit, d$gene)],
    DEG_FDR = d$adj.P.Val[match(hit, d$gene)],
    dir     = d$dir[match(hit, d$gene)],
    kME     = round(kme_of, 3),
    GS_RA   = gt$GS_RA[match(hit, gt$gene)])

  # NB: setorder() accepts column NAMES only - it cannot evaluate abs().
  # Ordering in i does evaluate expressions, so use that instead.
  out <- out[order(-abs(logFC))]
  fwrite(out, file.path(tab, sprintf("candidates_%s_disease.csv", tolower(sx))))

  list(out = out, sig = sig, n_lost = length(lost))
}

cat("\n===== INTERSECTION =====\n")
R <- setNames(lapply(c("Female", "Male"), build), c("Female", "Male"))

# ---- summary, with the per-module breakdown computed dynamically ------------
summ <- rbindlist(lapply(names(R), function(sx) {
  o <- R[[sx]]$out
  bymod <- table(factor(o$module, levels = dis_mods))
  data.table(
    sex                          = sx,
    n_sig_DEG                    = length(R[[sx]]$sig),
    n_disease_module_genes       = length(disease_genes),
    n_candidates                 = nrow(o),
    n_up_in_RA                   = sum(o$dir == "Up in RA"),
    n_down_in_RA                 = sum(o$dir == "Down in RA"),
    DEGs_lost_to_variance_filter = R[[sx]]$n_lost,
    pct_DEG_lost                 = round(100 * R[[sx]]$n_lost / length(R[[sx]]$sig), 1),
    per_module                   = paste(sprintf("%s=%d", names(bymod), as.integer(bymod)),
                                         collapse = " "))
}))
fwrite(summ, file.path(tab, "candidate_summary.csv"))
print(summ)

shared <- intersect(R$Female$out$gene, R$Male$out$gene)
say("")
say("Female candidates : %d", nrow(R$Female$out))
say("Male candidates   : %d", nrow(R$Male$out))
say("shared by both    : %d", length(shared))
say("union             : %d", length(union(R$Female$out$gene, R$Male$out$gene)))

# ---- 2-set Venn per sex ------------------------------------------------------
save_venn <- function(sets, file, title, hi) {
  p <- ggVennDiagram(sets, label = "count", label_alpha = 0,
                     edge_size = 0.5, set_size = 3.6) +
    scale_fill_gradient(low = "#F7FBFF", high = hi) +
    scale_x_continuous(expand = expansion(mult = 0.28)) +
    coord_cartesian(clip = "off") + labs(title = title) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11),
          plot.margin = margin(10, 26, 10, 26))
  ggsave(file.path(fig, file), p, width = 6.5, height = 5.6, dpi = 300)
  say("  wrote %s", file)
}

mlab <- sprintf("Disease modules\n(%s)", paste(dis_mods, collapse = " + "))
for (sx in c("Female", "Male")) {
  save_venn(setNames(list(disease_genes, R[[sx]]$sig), c(mlab, paste(sx, "DEGs"))),
            sprintf("fig_venn_%s_disease_candidates.png", tolower(sx)),
            sprintf("%s: disease modules n %s DEGs = %d", sx, sx, nrow(R[[sx]]$out)),
            if (sx == "Female") "#2E7D32" else "#795548")
}

# =============================================================================
# COMPOSITE FIGURE: the two candidate Venns beside the per-stratum module-trait
# heatmap. One panel that shows WHICH modules were selected, HOW strongly each
# stratum supports them, and HOW MANY candidates each sex yields.
# Requires 08_module_trait_RA_control.R to have run (writes the per-stratum CSVs).
# =============================================================================
mt_files <- file.path(tab, sprintf("module_trait_RAvsControl_%s.csv",
                                   c("all", "female", "male")))
if (!all(file.exists(mt_files))) {
  say("")
  say("Composite figure SKIPPED - run 08_module_trait_RA_control.R first.")
} else {
  suppressMessages({library(patchwork); library(ggplotify)})
  labs3 <- c(sprintf("All\n(n=%d)",    fread(mt_files[1])$n_samples[1]),
             sprintf("Female\n(n=%d)", fread(mt_files[2])$n_samples[1]),
             sprintf("Male\n(n=%d)",   fread(mt_files[3])$n_samples[1]))
  hm <- rbindlist(lapply(seq_along(mt_files), function(i)
          fread(mt_files[i])[, .(module, size, cor = cor_RA, p = p_RA, set = labs3[i])]))
  ordm <- fread(mt_files[1])[order(-cor_RA)]$module
  hm[, module := factor(module, levels = rev(ordm))]
  hm[, set    := factor(set, levels = labs3)]
  star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**",
                       ifelse(p < 5e-2, "*", "")))
  hm[, lab := sprintf("%.2f%s", cor, star(p))]
  hm[, selected := as.character(module) %in% dis_mods]
  ylabs <- setNames(sprintf("%s (%d)", ordm,
             fread(mt_files[1])$size[match(ordm, fread(mt_files[1])$module)]), ordm)

  p_hm <- ggplot(hm, aes(set, module, fill = cor)) +
    geom_tile(colour = "grey88", linewidth = 0.4) +
    geom_tile(data = hm[selected == TRUE], colour = "black",
              linewidth = 1.2, fill = NA) +
    geom_text(aes(label = lab), size = 3.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "cor\n(RA)") +
    scale_y_discrete(labels = ylabs) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "WGCNA module (size)",
         title = sprintf(paste("Boxed = selected disease modules",
                               "(|cor_all| >= %.1f & p < %g).  *p<.05  **p<.01  ***p<.001"),
                         w$config$dm_min_abs_cor, w$config$dm_max_p)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = ifelse(rev(ordm) %in% dis_mods,
                                                   "bold", "plain")),
          plot.title = element_text(size = 8.6, hjust = 0.5))

  venn_panel <- function(sx, hi) {
    ggVennDiagram(setNames(list(R[[sx]]$sig, disease_genes),
                           c(paste(sx, "DEGs"), mlab)),
                  label = "count", label_alpha = 0, edge_size = 0.5,
                  set_size = 3.1, label_size = 3.6) +
      scale_fill_gradient(low = "#FFFFFF", high = hi) +
      scale_x_continuous(expand = expansion(mult = 0.30)) +
      scale_y_continuous(expand = expansion(mult = 0.14)) +
      coord_cartesian(clip = "off") +
      labs(title = sprintf("%s: disease modules n %s DEGs = %d",
                           sx, sx, nrow(R[[sx]]$out))) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 9.5),
            plot.margin = margin(6, 20, 6, 20))
  }

  comp <- (venn_panel("Male", "#795548") / venn_panel("Female", "#2E7D32")) |
          p_hm
  comp <- comp + plot_layout(widths = c(1, 1.55))
  ggsave(file.path(fig, "fig_disease_module_selection_composite.png"), comp,
         width = 15, height = 8.5, dpi = 300)
  ggsave(file.path(fig, "fig_disease_module_selection_composite.pdf"), comp,
         width = 15, height = 8.5)
  say("")
  say("  wrote fig_disease_module_selection_composite.{png,pdf}")
}

say("")
say("Wrote candidates_{female,male}_disease.csv -> consumed by the MR step (10/11)")
cat("DONE\n")
