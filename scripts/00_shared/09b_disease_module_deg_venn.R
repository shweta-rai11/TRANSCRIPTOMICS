#!/usr/bin/env Rscript
# =============================================================================
# 09b_disease_module_deg_venn.R
# -----------------------------------------------------------------------------
# Directional Venn diagrams: for EACH disease module, its overlap with the
# UP-in-RA and DOWN-in-RA DEGs, within each sex.
#
# ***  RENAMED FROM 09_greenmodule_deg_venn.R  ***
#   The old script was written around `gt[module == "green"]` with the header
#   "the GREEN module (cor with RA = +0.57)". That is not a valid identifier.
#   WGCNA assigns colour names by module SIZE RANK (labels2colors): the largest
#   module becomes turquoise, the 2nd blue, 3rd brown, 4th yellow, 5th green.
#   Change the gene filter, the soft power or the sample set and the SAME
#   biology is renamed. That is exactly why the project README says the disease
#   modules are "green + brown" while the previous script hardcoded
#   "yellow + blue" - both were snapshots of different runs, and a script keyed
#   to a colour string silently returns ZERO genes when the colour moves.
#
#   This version takes the disease modules from the analysis object and loops
#   over whatever they turn out to be. No colour literal appears anywhere.
#
# DIRECTIONAL LOGIC
#   A module with cor(ME, RA) > 0 is UP in RA, so its biologically consistent
#   overlap is with the UP-in-RA DEGs; the DOWN overlap is the inconsistent one
#   and should be small. For a module with cor(ME, RA) < 0 the expectation
#   reverses. The script labels which overlap is the CONSISTENT one for each
#   module rather than assuming a direction.
#
# INTERPRETATION LIMIT
#   Female and Male panels are separate within-sex contrasts. A larger overlap
#   in one sex is not evidence of a sex difference - see Gelman & Stern (2006).
#
# Inputs : data/processed/wgcna_results.rds  (moduleColors, disease_modules, kME)
#          data/processed/dge_results.rds    (res$Female / res$Male; sig, dir)
# Outputs: results/figures/fig_diseasemod_venn_{module}_{sex}_{up,down}.png
#          results/tables/diseasemod_DEG_intersection_{sex}.csv
#          results/tables/diseasemod_DEG_direction_summary.csv
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

wf <- file.path(proc, "wgcna_results.rds")
df <- file.path(proc, "dge_results.rds")
if (!file.exists(wf)) stop("Missing ", wf, " - run 06_WGCNA.R first.")
if (!file.exists(df)) stop("Missing ", df, " - run 05_dge.R first.")
w <- readRDS(wf); D <- readRDS(df)

gt       <- as.data.table(w$gene_tab)
kME      <- as.data.frame(w$kME)
dis_mods <- w$disease_modules
mt       <- as.data.table(w$module_trait_table)

say("disease modules: %s", paste(dis_mods, collapse = " + "))
for (m in dis_mods)
  say("  %-11s n=%-5d cor(ME,RA)=%+.3f  ->  UP-in-RA overlap is the CONSISTENT one: %s",
      m, gt[module == m, .N], mt[module == m, cor_RA],
      ifelse(mt[module == m, cor_RA] > 0, "YES", "NO (module is DOWN in RA)"))

save_venn <- function(sets, file, title, hi) {
  p <- ggVennDiagram(sets, label = "count", label_alpha = 0,
                     edge_size = 0.5, set_size = 3.8) +
    scale_fill_gradient(low = "#F7FBFF", high = hi) +
    scale_x_continuous(expand = expansion(mult = 0.28)) +
    scale_y_continuous(expand = expansion(mult = 0.10)) +
    coord_cartesian(clip = "off") + labs(title = title) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 10.5),
          plot.margin = margin(8, 22, 8, 22))
  ggsave(file.path(fig, file), p, width = 6.6, height = 5.4, dpi = 300)
}

all_rows <- list(); summary_rows <- list()

for (sx in c("Female", "Male")) {
  d    <- as.data.table(D$res[[sx]])
  up   <- d[sig == TRUE & dir == "Up in RA",   gene]
  down <- d[sig == TRUE & dir == "Down in RA", gene]
  cat(sprintf("\n===== %s : %d up-in-RA / %d down-in-RA DEGs =====\n",
              sx, length(up), length(down)))

  for (m in dis_mods) {
    mg      <- gt[module == m, gene]
    r       <- mt[module == m, cor_RA]
    consist <- if (r > 0) "up" else "down"          # direction expected for THIS module
    int_up   <- intersect(mg, up)
    int_down <- intersect(mg, down)

    save_venn(setNames(list(mg, up), c(sprintf("%s module\n(cor=%+.2f)", m, r),
                                       paste0(sx, " DEG_up"))),
              sprintf("fig_diseasemod_venn_%s_%s_up.png", m, tolower(sx)),
              sprintf("%s | %s module n DEG_up = %d%s",
                      sx, m, length(int_up),
                      if (consist == "up") "   [direction-consistent]" else ""),
              "#C62828")
    save_venn(setNames(list(mg, down), c(sprintf("%s module\n(cor=%+.2f)", m, r),
                                         paste0(sx, " DEG_down"))),
              sprintf("fig_diseasemod_venn_%s_%s_down.png", m, tolower(sx)),
              sprintf("%s | %s module n DEG_down = %d%s",
                      sx, m, length(int_down),
                      if (consist == "down") "   [direction-consistent]" else ""),
              "#1565C0")

    kcol <- paste0("kME", m)
    ann <- function(genes, dirlab) if (!length(genes)) NULL else
      data.table(sex = sx, module = m, cor_ME_RA = round(r, 3),
                 direction = dirlab,
                 consistent = (substr(dirlab, 1, 2) == consist),
                 gene = genes,
                 logFC = d$logFC[match(genes, d$gene)],
                 DEG_FDR = d$adj.P.Val[match(genes, d$gene)],
                 kME = if (kcol %in% colnames(kME))
                         round(kME[genes, kcol], 3) else NA_real_)
    all_rows[[paste(sx, m)]] <- rbind(ann(int_up, "up in RA"),
                                      ann(int_down, "down in RA"), fill = TRUE)

    summary_rows[[paste(sx, m)]] <- data.table(
      sex = sx, module = m, module_size = length(mg), cor_ME_RA = round(r, 3),
      module_direction = ifelse(r > 0, "UP in RA", "DOWN in RA"),
      n_overlap_up = length(int_up), n_overlap_down = length(int_down),
      n_consistent = if (consist == "up") length(int_up) else length(int_down),
      n_inconsistent = if (consist == "up") length(int_down) else length(int_up))

    say("  %-11s cor=%+.2f | n DEG_up=%-5d n DEG_down=%-5d  (consistent = %s)",
        m, r, length(int_up), length(int_down), consist)
  }

  rows <- rbindlist(all_rows[grep(paste0("^", sx, " "), names(all_rows))], fill = TRUE)
  if (nrow(rows))
    fwrite(rows, file.path(tab, sprintf("diseasemod_DEG_intersection_%s.csv", tolower(sx))))
}

sm <- rbindlist(summary_rows)
fwrite(sm, file.path(tab, "diseasemod_DEG_direction_summary.csv"))
cat("\n===== DIRECTION SUMMARY =====\n"); print(sm)
cat("\nA well-behaved disease module has n_consistent >> n_inconsistent.\n")
cat("DONE\n")
