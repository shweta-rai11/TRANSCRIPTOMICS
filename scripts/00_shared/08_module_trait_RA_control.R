#!/usr/bin/env Rscript
# =============================================================================
# 08_module_trait_RA_control.R
# -----------------------------------------------------------------------------
# Module-trait association, RA vs Control, computed SEPARATELY within each sex
# stratum (All / Female / Male), plus a combined selection figure showing which
# modules were chosen as DISEASE modules and why.
#
# WHY BOTH HALVES ARE IN ONE FILE
#   The previous version of this file contained TWO complete scripts
#   concatenated (a second `#!/usr/bin/env Rscript` began at line 59). Only the
#   first half ran under Rscript in some invocations, and the second half
#   silently depended on tables the first half had just written. They are one
#   logical step and are now one script.
#
# WHAT CHANGED (2026-07-27)
#   - Reads data/processed/wgcna_results.rds (written by 06_WGCNA.R).
#     The old file wgcna_new_results.rds came from the deleted 06_WGCNA_NEW.R.
#   - NO hardcoded module colours. The previous version hardcoded
#     `sel <- c("green","brown")`. WGCNA colour names are assigned by module
#     SIZE RANK, so they are not stable identifiers; the disease modules are
#     read from the analysis object instead.
#   - Sample counts in figure labels are COMPUTED, not hardcoded ("n=183",
#     "n=145", "n=38" were stale after outlier removal).
#
# IMPORTANT INTERPRETATION NOTE (sex-stratified, not sex-specific)
#   The Female and Male panels are each an RA-vs-Control contrast computed
#   WITHIN that sex. A module correlating more strongly with RA in males than in
#   females is NOT evidence of a sex difference - the male stratum is far
#   smaller, so its correlations are noisier and its extremes larger. Comparing
#   the two columns is invalid without a formal interaction test.
#     # Gelman A, Stern H. The difference between "significant" and "not
#     #   significant" is not itself statistically significant.
#     #   Am Stat 2006;60(4):328-331.
#
# Inputs : data/processed/wgcna_results.rds
# Outputs: results/tables/module_trait_RAvsControl_{all,female,male}.csv
#          results/figures/fig_wgcna_module_trait_RAvsControl_{all,female,male}.png
#          results/figures/fig_module_trait_disease_selection.{png,pdf}
#
# ---- References -------------------------------------------------------------
#   Langfelder P, Horvath S. WGCNA: an R package for weighted correlation
#     network analysis. BMC Bioinformatics 2008;9:559.
#   Langfelder P, Horvath S. Eigengene networks for studying the relationships
#     between co-expression modules. BMC Syst Biol 2007;1:54.
# =============================================================================
suppressMessages({
  library(WGCNA); library(data.table); library(ggplot2)
})
options(stringsAsFactors = FALSE)

proc <- "data/processed"; tab <- "results/tables"; fig <- "results/figures"
for (d in c(tab, fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---- load -------------------------------------------------------------------
wf <- file.path(proc, "wgcna_results.rds")
if (!file.exists(wf)) stop("Missing ", wf, " - run scripts/00_shared/06_WGCNA.R first.")
w        <- readRDS(wf)
MEs      <- w$MEs
meta     <- w$meta
sizes    <- table(w$moduleColors)
dis_mods <- w$disease_modules          # data-driven, from 06_WGCNA.R STEP 9
stopifnot(identical(rownames(MEs), meta$sample))

say("loaded wgcna_results.rds : %d modules x %d samples", ncol(MEs), nrow(MEs))
say("disease modules (from 06): %s", paste(dis_mods, collapse = " + "))
say("selection rule           : |cor(ME,RA)| >= %.2f & p < %g",
    w$config$dm_min_abs_cor, w$config$dm_max_p)

# =============================================================================
# PART A — per-stratum module-trait tables + labeledHeatmap figures
# =============================================================================
make_stratum <- function(keep, tag, label) {
  me  <- MEs[keep, , drop = FALSE]
  grp <- meta$group[keep]
  n   <- sum(keep)
  if (length(unique(grp)) < 2) { say("  %s: only one group present - skipped", label); return(NULL) }

  trait <- cbind(Control = as.numeric(grp == "HC"),
                 RA      = as.numeric(grp == "RA"))
  C <- stats::cor(me, trait, use = "pairwise.complete.obs")
  P <- corPvalueStudent(C, n)
  ord <- order(-C[, "RA"])
  C <- C[ord, , drop = FALSE]; P <- P[ord, , drop = FALSE]

  # star the disease modules on the y-axis so the figure is self-explaining
  msym <- sub("^ME", "", rownames(C))
  ysym <- ifelse(msym %in% dis_mods, paste0("* ", msym), msym)

  txt <- paste0(signif(C, 2), "\n(", signif(P, 1), ")"); dim(txt) <- dim(C)
  png(file.path(fig, sprintf("fig_wgcna_module_trait_RAvsControl_%s.png", tag)),
      width = 640, height = 950, res = 120)
  par(mar = c(4, 10, 4, 3))
  labeledHeatmap(Matrix = C, xLabels = colnames(C), yLabels = rownames(C),
    ySymbols = ysym, colorLabels = FALSE, colors = blueWhiteRed(50),
    textMatrix = txt, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1, 1),
    main = sprintf("Module-trait: Control vs RA\n%s (n = %d)   * = disease module",
                   label, n))
  dev.off()

  out <- data.table(module      = msym,
                    size        = as.integer(sizes[msym]),
                    n_samples   = n,
                    cor_Control = round(C[, "Control"], 3),
                    p_Control   = signif(P[, "Control"], 3),
                    cor_RA      = round(C[, "RA"], 3),
                    p_RA        = signif(P[, "RA"], 3),
                    is_disease_module = msym %in% dis_mods)
  fwrite(out, file.path(tab, sprintf("module_trait_RAvsControl_%s.csv", tag)))
  say("  %-8s n=%-4d  wrote table + heatmap", label, n)
  out[, stratum := label][]
}

cat("\n--- PART A: module-trait per stratum ---\n")
res <- list(
  all    = make_stratum(rep(TRUE, nrow(meta)),  "all",    "All samples"),
  female = make_stratum(meta$sex == "F",        "female", "Female"),
  male   = make_stratum(meta$sex == "M",        "male",   "Male"))
res <- Filter(Negate(is.null), res)

# =============================================================================
# PART B — combined selection figure (All / Female / Male side by side)
# =============================================================================
cat("\n--- PART B: combined disease-module selection figure ---\n")
d <- rbindlist(lapply(res, function(x)
       x[, .(module, cor = cor_RA, p = p_RA, n = n_samples, stratum)]))
# stratum labels carry their COMPUTED n (previously hardcoded and stale)
d[, set := sprintf("%s\n(n=%d)", stratum, n)]
lev <- unique(d[order(match(stratum, c("All samples", "Female", "Male")))]$set)
d[, set := factor(set, levels = lev)]

ord <- res$all[order(-cor_RA)]$module              # order rows by overall cor
d[, module := factor(module, levels = rev(ord))]
star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 5e-2, "*", "")))
d[, lab := sprintf("%.2f%s", cor, star(p))]
d[, selected := as.character(module) %in% dis_mods]

ylab <- setNames(sprintf("%s (%d)", ord, as.integer(sizes[ord])), ord)

g <- ggplot(d, aes(x = set, y = module, fill = cor)) +
  geom_tile(colour = "grey90", linewidth = 0.4) +
  geom_tile(data = d[selected == TRUE], colour = "black", linewidth = 1.1, fill = NA) +
  geom_text(aes(label = lab), size = 3.1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "cor\n(RA)") +
  scale_y_discrete(labels = ylab) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "WGCNA module (size)",
       subtitle = sprintf(paste("Boxed = disease modules (|cor_all| >= %.1f & p < %g),",
                                "selected on ALL samples and applied to BOTH sexes.",
                                "\n*p<.05 **p<.01 ***p<.001.",
                                "Female and Male columns are SEX-STRATIFIED contrasts;",
                                "\nthey must NOT be compared to each other without an interaction test."),
                          w$config$dm_min_abs_cor, w$config$dm_max_p)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = ifelse(rev(ord) %in% dis_mods, "bold", "plain")),
        plot.subtitle = element_text(size = 7.6))

ggsave(file.path(fig, "fig_module_trait_disease_selection.png"), g,
       width = 7.2, height = 6.4, dpi = 300)
ggsave(file.path(fig, "fig_module_trait_disease_selection.pdf"), g,
       width = 7.2, height = 6.4)
say("wrote fig_module_trait_disease_selection.{png,pdf}")

# combined long table for the supplement
fwrite(rbindlist(res), file.path(tab, "module_trait_RAvsControl_ALLSTRATA.csv"))

cat("\n=== SUMMARY: cor(ME, RA) for the disease modules, per stratum ===\n")
print(dcast(d[selected == TRUE], module ~ set, value.var = "lab"))
cat("\nDONE\n")
