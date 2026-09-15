#!/usr/bin/env Rscript
# Venn diagrams of the feature-selection step: per-sex selector overlap (primary and MHC-free candidate sets) and female-vs-male consensus panel membership.
suppressMessages({
  library(ggVennDiagram); library(ggplot2); library(patchwork); library(data.table)
})
options(stringsAsFactors = FALSE)

procN <- "data/processed/new"; tab <- "results/tables"; fig <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# Step 1: load the selector sets
hdr("STEP 1  SELECTOR SETS")
mlP <- readRDS(file.path(procN, "ml_features.rds"))
mlN_path <- file.path(procN, "ml_features_noMHC.rds")
mlN <- if (file.exists(mlN_path)) readRDS(mlN_path) else NULL
if (is.null(mlN))
  say("ml_features_noMHC.rds not found - MHC-free panels will be skipped (run 12b).")

get_sets <- function(ml, sx) {
  r <- ml[[sx]]
  list(LASSO = r$sets$LASSO, `Random Forest` = r$sets$RandomForest,
       `SVM-RFE` = r$sets$SVM_RFE)
}
sets <- list(
  primary_female = get_sets(mlP, "female"),
  primary_male   = get_sets(mlP, "male"))
if (!is.null(mlN)) {
  sets$noMHC_female <- get_sets(mlN, "female")
  sets$noMHC_male   <- get_sets(mlN, "male")
}
cons <- list(
  primary_female = mlP$female$consensus, primary_male = mlP$male$consensus,
  noMHC_female = if (!is.null(mlN)) mlN$female$consensus else character(0),
  noMHC_male   = if (!is.null(mlN)) mlN$male$consensus   else character(0))

for (nm in names(sets))
  say("%-15s LASSO %2d | RF %2d | SVM-RFE %2d -> consensus %d (%s)",
      nm, length(sets[[nm]]$LASSO), length(sets[[nm]]$`Random Forest`),
      length(sets[[nm]]$`SVM-RFE`), length(cons[[nm]]),
      if (length(cons[[nm]])) paste(cons[[nm]], collapse = ", ") else "none")

# Step 2: draw Venn diagrams (consensus genes printed as a caption below each, not inside it)
hdr("STEP 2  DRAW VENN DIAGRAMS")

venn3 <- function(s, title, subtitle, fill_low, fill_high) {
  ggVennDiagram(s, label = "count", label_alpha = 0, edge_size = 0.4,
                set_size = 3.4, label_size = 3.6) +
    scale_fill_gradient(low = fill_low, high = fill_high, guide = "none") +
    scale_colour_manual(values = rep("grey25", 3), guide = "none") +
    labs(title = title, subtitle = subtitle) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, hjust = 0.5,
                                       margin = margin(t = 2, b = 6)),
          plot.margin   = margin(10, 16, 10, 16))
}

wrap_genes <- function(g, width = 46) {
  if (!length(g)) return("consensus: none")
  paste0("consensus (", length(g), "): ",
         paste(strwrap(paste(g, collapse = ", "), width = width), collapse = "\n"))
}

pal <- list(female = c("#F7E6EF", "#C2185B"), male = c("#E3EEF7", "#1565C0"))

plots <- list()
for (sx in c("female", "male")) {
  k <- paste0("primary_", sx)
  plots[[k]] <- venn3(sets[[k]],
    sprintf("%s panel - three-selector consensus", tools::toTitleCase(sx)),
    wrap_genes(cons[[k]]), pal[[sx]][1], pal[[sx]][2])
  if (!is.null(mlN)) {
    k2 <- paste0("noMHC_", sx)
    plots[[k2]] <- venn3(sets[[k2]],
      sprintf("%s panel - MHC-free candidate set", tools::toTitleCase(sx)),
      wrap_genes(cons[[k2]]), pal[[sx]][1], pal[[sx]][2])
  }
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  say("  wrote %s.{png,pdf}", name)
}

save_fig(plots$primary_female, "FIG_G2_01_panel_venn_female", 5.2, 4.6)
save_fig(plots$primary_male,   "FIG_G2_01_panel_venn_male",   5.2, 4.6)
if (!is.null(mlN)) {
  save_fig(plots$noMHC_female, "FIG_G2_02_panel_venn_noMHC_female", 5.2, 4.6)
  save_fig(plots$noMHC_male,   "FIG_G2_02_panel_venn_noMHC_male",   5.2, 4.6)
}

# ---- female vs male panel membership ----------------------------------------
fm_sets <- list(Female = cons$primary_female, Male = cons$primary_male)
shared  <- intersect(cons$primary_female, cons$primary_male)
# Widen x range so two-set category labels (drawn outside the circles) aren't clipped
p_fm <- ggVennDiagram(fm_sets, label = "count", label_alpha = 0, edge_size = 0.4,
                      set_size = 3.6, label_size = 4) +
  scale_fill_gradient(low = "#F3F0FA", high = "#5E35B1", guide = "none") +
  scale_colour_manual(values = rep("grey25", 2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.30)) +
  scale_y_continuous(expand = expansion(mult = 0.16)) +
  labs(title = "Female vs male consensus panel",
       subtitle = paste0(
         "shared (", length(shared), "): ",
         if (length(shared)) paste(shared, collapse = ", ") else "none", "\n",
         paste(strwrap(paste("shared genes carry IDENTICAL MR estimates -",
                             "both GWAS inputs are sex-combined"),
                       width = 52), collapse = "\n"))) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(t = 2, b = 6)),
        plot.margin = margin(10, 26, 10, 26))
save_fig(p_fm, "FIG_G2_03_panel_venn_female_vs_male", 5.8, 4.8)

# ---- composite ---------------------------------------------------------------
# Panel E spans the full width (its labels clip if shared with a row at half width)
comp <- if (!is.null(mlN))
  (plots$primary_female | plots$primary_male) /
  (plots$noMHC_female   | plots$noMHC_male)   /
  p_fm +
  plot_layout(heights = c(1, 1, 1.05)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13))) else
  (plots$primary_female | plots$primary_male) / p_fm +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13)))
save_fig(comp, "FIG_G2_04_feature_selection_venn_composite", 10.4,
         if (!is.null(mlN)) 14 else 9.4)

# Step 3: membership table
hdr("STEP 3  MEMBERSHIP TABLE")
memb <- rbindlist(lapply(names(sets), function(nm) {
  s <- sets[[nm]]
  g <- sort(Reduce(union, s))
  data.table(candidate_set = sub("_.*$", "", nm),
             sex = sub("^.*_", "", nm),
             gene = g,
             LASSO = g %in% s$LASSO,
             RandomForest = g %in% s$`Random Forest`,
             SVM_RFE = g %in% s$`SVM-RFE`,
             n_methods = (g %in% s$LASSO) + (g %in% s$`Random Forest`) +
                         (g %in% s$`SVM-RFE`),
             in_consensus = g %in% cons[[nm]])
}))
setorder(memb, candidate_set, sex, -n_methods, gene)
fwrite(memb, file.path(tab, "FS_venn_membership.csv"))
say("wrote results/tables/FS_venn_membership.csv (%d rows)", nrow(memb))
print(memb[in_consensus == TRUE, .(candidate_set, sex, gene, n_methods)])

say("")
say("Every panel figure now has a script behind it (FIG_G2_01_* included).")
cat("\nDONE\n")
