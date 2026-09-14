#!/usr/bin/env Rscript
# All figures for the DGE step (plotting only): volcano, DEG counts, Venn overlaps, heatmaps, enrichment dot/GSEA plots
suppressMessages({
  library(ggplot2); library(data.table); library(VennDiagram)
  library(EnhancedVolcano); library(patchwork); library(grid); library(ggVennDiagram)
  library(ComplexHeatmap); library(circlize); library(enrichplot)
})
theme_set(theme_bw(base_size = 11))
set.seed(1234)
proc <- "data/processed"; fig <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

D   <- readRDS(file.path(proc, "dge_results.rds"))
res <- D$res; FC <- D$thr$FC; PVAL <- D$thr$PVAL
o   <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr; meta <- o$meta

save2 <- function(g, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), g, width = w, height = h, dpi = 300)
  ggsave(file.path(fig, paste0(name, ".pdf")), g, width = w, height = h, device = cairo_pdf)
}
dcols <- c(`Up in RA` = "#C62828", `Down in RA` = "#1565C0", ns = "grey78")

# 1. Volcano: three comparisons faceted (overview)
volc <- rbindlist(lapply(names(res), function(nm) { d <- copy(res[[nm]]); d$comparison <- nm; d }))
volc$comparison <- factor(volc$comparison, levels = c("All","Female","Male"))
volc$y <- -log10(pmax(volc$P.Value, 1e-300))
lab <- volc[sig == TRUE][order(P.Value)][, head(.SD, 8), by = comparison]
pv <- ggplot(volc, aes(logFC, y, color = dir)) +
  geom_point(size = .7, alpha = .55) +
  geom_vline(xintercept = c(-FC, FC), linetype = 2, colour = "grey50", linewidth = .3) +
  geom_hline(yintercept = -log10(PVAL), linetype = 2, colour = "grey50", linewidth = .3) +
  geom_text(data = lab, aes(label = gene), size = 2.3, vjust = -0.5, show.legend = FALSE, check_overlap = TRUE) +
  facet_wrap(~comparison, scales = "free") +
  scale_color_manual(values = dcols, name = NULL) +
  labs(x = "log2 fold-change (RA vs HC)", y = "-log10(p-value)") +
  theme(legend.position = "top")
save2(pv, "fig_dge_volcano", 9.5, 5)

# 2. Volcano (EnhancedVolcano), Female | Male side by side, locked axes; byFC vs byP label selection
resF <- as.data.frame(res$Female); resM <- as.data.frame(res$Male)

# locked axis limits shared across both sexes
allFC <- c(resF$logFC, resM$logFC)
allP  <- pmax(c(resF$P.Value, resM$P.Value), 1e-300)
xmax  <- ceiling(max(abs(allFC)) * 1.05 * 10) / 10
ymax  <- ceiling(max(-log10(allP)) * 1.05)
XLIM  <- c(-xmax, xmax); YLIM <- c(0, ymax)
cat(sprintf("Sex volcano - locked axes x=[%.2f,%.2f] y=[0,%.1f]\n", -xmax, xmax, ymax))

# top-10 up + top-10 down labels (significant genes only)
pick_labs <- function(df, by = c("fc","pval"), n = 10) {
  by <- match.arg(by); s <- df[df$sig, , drop = FALSE]
  up <- s[s$logFC > 0, , drop = FALSE]; dn <- s[s$logFC < 0, , drop = FALSE]
  if (by == "fc") { up <- up[order(-up$logFC), ]; dn <- dn[order(dn$logFC), ] }
  else            { up <- up[order(up$P.Value), ]; dn <- dn[order(dn$P.Value), ] }
  c(head(up$gene, n), head(dn$gene, n))
}

# one EnhancedVolcano panel (its own legend suppressed; shared legend added later)
ev_panel <- function(df, sextitle, labs) {
  sub <- sprintf("%d significant genes  (%d up / %d down)",
                 sum(df$sig), sum(df$sig & df$logFC > 0), sum(df$sig & df$logFC < 0))
  EnhancedVolcano(df, lab = df$gene, x = "logFC", y = "P.Value",
    selectLab = labs, xlim = XLIM, ylim = YLIM,
    pCutoff = PVAL, FCcutoff = FC,
    title = sextitle, subtitle = sub, caption = NULL,
    titleLabSize = 13, subtitleLabSize = 9,
    pointSize = 1.4, labSize = 3.2, colAlpha = 0.6,
    drawConnectors = TRUE, widthConnectors = 0.4, arrowheads = FALSE,
    max.overlaps = Inf, boxedLabels = FALSE, borderWidth = 0.6,
    gridlines.major = FALSE, gridlines.minor = FALSE,
    legendPosition = "none") +
    theme(plot.title = element_text(face = "bold"))
}

# one shared 4-category legend (built from a helper plot where every category
# has a real point, so colours are correct); extracted as a grob and attached once
shared_legend <- function() {
  d <- data.frame(x = 1:4, y = 1,
        Sig = factor(c("NS","FC","P","FC_P"), levels = c("NS","FC","P","FC_P")))
  lp <- ggplot(d, aes(x, y, colour = Sig)) + geom_point(size = 3.5) +
    scale_colour_manual(
      values = c(NS = "grey30", FC = "forestgreen", P = "royalblue", FC_P = "red2"),
      labels = c("NS", "Log2 FC", "p-value", "p-value and Log2 FC"), name = NULL) +
    guides(colour = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 11))
  g  <- ggplotGrob(lp)
  gb <- g$grobs[grepl("guide-box", vapply(g$grobs, `[[`, character(1), "name"))]
  gb <- gb[vapply(gb, function(x) !inherits(x, "zeroGrob"), logical(1))]
  gb[[1]]
}
LEG <- shared_legend()

make_sex_volcano <- function(by, tag, ranklabel) {
  pF <- ev_panel(resF, "Female", pick_labs(resF, by))
  pM <- ev_panel(resM, "Male",   pick_labs(resM, by))
  g <- ((pF | pM) / wrap_elements(full = LEG)) +
    plot_layout(heights = c(1, 0.07))
  save2(g, sprintf("fig_dge_volcano_sex_%s", tag), 9.5, 5)
  cat(sprintf("Wrote fig_dge_volcano_sex_%s.{png,pdf}\n", tag))

  # same panels, saved as two SEPARATE single-sex figures (each with its own
  # copy of the shared legend attached) so Female/Male can be used on their own
  gF <- (pF / wrap_elements(full = LEG)) + plot_layout(heights = c(1, 0.09))
  gM <- (pM / wrap_elements(full = LEG)) + plot_layout(heights = c(1, 0.09))
  save2(gF, sprintf("fig_dge_volcano_sex_%s_female", tag), 9.5, 5)
  save2(gM, sprintf("fig_dge_volcano_sex_%s_male",   tag), 9.5, 5)
  cat(sprintf("Wrote fig_dge_volcano_sex_%s_{female,male}.{png,pdf}\n", tag))
}
make_sex_volcano("fc",   "byFC", "log2 fold-change magnitude")
make_sex_volcano("pval", "byP",  "statistical significance (p-value)")

# 3. DEG counts barplot
cnt <- rbindlist(lapply(names(res), function(nm) {
  d <- res[[nm]]
  data.table(comparison = nm, `Up in RA` = sum(d$dir == "Up in RA"),
             `Down in RA` = sum(d$dir == "Down in RA")) }))
cntl <- melt(cnt, id.vars = "comparison", variable.name = "direction", value.name = "n")
cntl$comparison <- factor(cntl$comparison, levels = c("All","Female","Male"))
pc <- ggplot(cntl, aes(comparison, n, fill = direction)) +
  geom_col(position = position_dodge(.55), width = .45, colour = "grey25", linewidth = .3) +
  geom_text(aes(label = n), position = position_dodge(.55), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c(`Up in RA` = "grey35", `Down in RA` = "grey78"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = sprintf("Significant DEGs per comparison (|log2FC| > %.1f, p < %.2f)", FC, PVAL),
       x = NULL, y = "n significant genes") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", axis.text = element_text(colour = "black"))
save2(pc, "fig_dge_counts", 9.5, 5)

# 4. Venn of significant DEGs (All / Female / Male), region fill = gene count (sqrt scale)
vsets <- list(All = res$All[sig==TRUE]$gene, Female = res$Female[sig==TRUE]$gene,
              Male = res$Male[sig==TRUE]$gene)
names(vsets) <- c(sprintf("All (n = %d)",    length(vsets[[1]])),
                  sprintf("Female (n = %d)", length(vsets[[2]])),
                  sprintf("Male (n = %d)",   length(vsets[[3]])))
vcol <- c("#2E7D32", "#C2185B", "#1565C0")            # All / Female / Male

vdat <- process_data(Venn(vsets))
vpoly <- venn_regionedge(vdat); vedg <- venn_setedge(vdat)
vsl <- venn_setlabel(vdat);     vrl <- venn_regionlabel(vdat)
vtot <- sum(vrl$count)
vrl$lab <- sprintf("%s\n(%.1f%%)", format(vrl$count, big.mark = ","), 100 * vrl$count / vtot)
vrl$txt <- ifelse(vrl$count > 500, "white", "grey12")           # white on dark fills
vsl$col <- vcol[seq_len(nrow(vsl))]
vedg$col <- vcol[as.integer(factor(vedg$id))]
vcx <- mean(range(vpoly$X)); vcy <- mean(range(vpoly$Y))        # push set labels outward
vsl$X <- vcx + (vsl$X - vcx) * 1.10; vsl$Y <- vcy + (vsl$Y - vcy) * 1.18

pvenn <- ggplot() +
  geom_polygon(data = vpoly, aes(X, Y, group = id, fill = count), alpha = 0.9) +
  geom_path(data = vedg, aes(X, Y, group = id, color = col), linewidth = 1.3, show.legend = FALSE) +
  geom_text(data = vsl, aes(X, Y, label = name, color = col), fontface = "bold", size = 5.2, show.legend = FALSE) +
  geom_text(data = vrl, aes(X, Y, label = lab), size = 4, fontface = "bold", color = vrl$txt, lineheight = 0.9) +
  scale_fill_gradientn(colours = c("#F7FBFF","#C6DBEF","#6BAED6","#2171B5","#08306B"),
                       trans = "sqrt", name = "genes\nper region") +
  scale_color_identity() + coord_equal(clip = "off") +
  theme_void(base_size = 13) +
  theme(plot.margin = margin(34, 26, 34, 26), legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9))
save2(pvenn, "fig_dge_venn", 9.5, 5)

# 5. Heatmaps of significant DEGs (Female, Male, All): rows = top-N up/down DEGs, columns = samples split Control|RA, z-scored
NTOP     <- 25                                      # top-N up + top-N down per panel
grp_col   <- c(Control = "#2E7D32", RA = "#C62828")
sex_col   <- c(F = "#D81B60", M = "#1565C0")
study_col <- c(GSE93272 = "#5E35B1", GSE110169 = "#00838F")
z_col     <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
mdt <- as.data.table(meta)

build_hm <- function(comp, samples, add_sex) {
  r <- res[[comp]]; sg <- r[sig == TRUE]
  up <- head(sg[dir == "Up in RA"][order(adj.P.Val)]$gene,   NTOP)
  dn <- head(sg[dir == "Down in RA"][order(adj.P.Val)]$gene, NTOP)
  genes  <- c(up, dn)
  rsplit <- factor(c(rep("Up in RA", length(up)), rep("Down in RA", length(dn))),
                   levels = c("Up in RA", "Down in RA"))
  md <- mdt[match(samples, sample)]
  z  <- t(scale(t(expr[genes, samples, drop = FALSE]))); z[is.na(z)] <- 0
  z[z > 2] <- 2; z[z < -2] <- -2
  csplit <- factor(ifelse(md$group == "HC", "Control", "RA"), levels = c("Control", "RA"))

  ann <- list(Group = ifelse(md$group == "HC", "Control", "RA"), Study = md$dataset)
  acol <- list(Group = grp_col, Study = study_col)
  if (add_sex) { ann$Sex <- md$sex; acol$Sex <- sex_col }
  top_anno <- ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(ann), col = acol,
    annotation_name_gp = grid::gpar(fontsize = 8), simple_anno_size = grid::unit(3.2, "mm"))

  ComplexHeatmap::Heatmap(z, name = "row z-score", col = z_col, top_annotation = top_anno,
    row_split = rsplit, column_split = csplit,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    show_column_names = FALSE, show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6.5, fontface = "italic"),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    row_gap = grid::unit(2, "mm"), column_gap = grid::unit(1.5, "mm"), border = TRUE,
    heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2), legend_height = grid::unit(3, "cm")))
}

draw_hm <- function(comp, samples, add_sex) {
  ht  <- build_hm(comp, samples, add_sex)
  ttl <- sprintf("%s : top differentially expressed genes, RA vs Control", comp)
  out <- file.path(fig, sprintf("fig_dge_heatmap_%s", tolower(comp)))
  for (dev_open in list(
        function() png(paste0(out, ".png"), width = 9.5, height = 5, units = "in", res = 300),
        function() cairo_pdf(paste0(out, ".pdf"), width = 9.5, height = 5))) {
    dev_open()
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
      merge_legend = TRUE, column_title = ttl,
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"))
    dev.off()
  }
  cat(sprintf("  heatmap %s: %d genes, %d samples\n", comp,
              nrow(ht@matrix), length(samples)))
}
draw_hm("Female", meta$sample[meta$sex == "F"], add_sex = FALSE)
draw_hm("Male",   meta$sample[meta$sex == "M"], add_sex = FALSE)
draw_hm("All",    meta$sample,                  add_sex = TRUE)

# 5b. Composite: single-sex volcano stacked on its DEG heatmap, one figure per sex
composite_deg <- function(sex_label, resSex, samples) {
  pv <- ev_panel(resSex, sex_label, pick_labs(resSex, "fc"))
  gv <- wrap_elements(full = (pv / wrap_elements(full = LEG)) + plot_layout(heights = c(1, 0.09)))
  ht  <- build_hm(sex_label, samples, add_sex = FALSE)
  ttl <- sprintf("%s : top differentially expressed genes, RA vs Control", sex_label)
  hgrob <- grid.grabExpr(ComplexHeatmap::draw(ht,
    heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
    column_title = ttl, column_title_gp = grid::gpar(fontsize = 13, fontface = "bold")))
  (gv / wrap_elements(full = hgrob)) + plot_layout(heights = c(1, 1.1))
}
gDEGm <- composite_deg("Male", resM, meta$sample[meta$sex == "M"])
ggsave(file.path(fig, "fig_dge_composite_male.png"), gDEGm, width = 10, height = 9.5, dpi = 300)
ggsave(file.path(fig, "fig_dge_composite_male.pdf"), gDEGm, width = 10, height = 9.5, device = cairo_pdf)
cat("Wrote fig_dge_composite_male.{png,pdf}\n")

gDEGf <- composite_deg("Female", resF, meta$sample[meta$sex == "F"])
ggsave(file.path(fig, "fig_dge_composite_female.png"), gDEGf, width = 10, height = 9.5, dpi = 300)
ggsave(file.path(fig, "fig_dge_composite_female.pdf"), gDEGf, width = 10, height = 9.5, device = cairo_pdf)
cat("Wrote fig_dge_composite_female.{png,pdf}\n")

# 5c. Composite: 3-panel volcano (All/Female/Male) stacked on the "All" DEG heatmap
ht_all    <- build_hm("All", meta$sample, add_sex = TRUE)
ttl_all   <- "All : top differentially expressed genes, RA vs Control"
hgrob_all <- grid.grabExpr(ComplexHeatmap::draw(ht_all,
  heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
  column_title = ttl_all, column_title_gp = grid::gpar(fontsize = 13, fontface = "bold")))
gDEGall <- (wrap_elements(full = pv) / wrap_elements(full = hgrob_all)) + plot_layout(heights = c(1, 1.1))
ggsave(file.path(fig, "fig_dge_composite_all.png"), gDEGall, width = 10, height = 9.5, dpi = 300)
ggsave(file.path(fig, "fig_dge_composite_all.pdf"), gDEGall, width = 10, height = 9.5, device = cairo_pdf)
cat("Wrote fig_dge_composite_all.{png,pdf}\n")

# 6. Venn of enriched GO/KEGG terms (Female / Male / All) [needs dge_enrich.rds]
enrich_path <- file.path(proc, "dge_enrich.rds")
if (file.exists(enrich_path)) {
  es <- readRDS(enrich_path)$sets
  esets <- list(All = es$All, Female = es$Female, Male = es$Male)     # same style as fig_dge_venn
  names(esets) <- c(sprintf("All (n = %d)",    length(esets[[1]])),
                    sprintf("Female (n = %d)", length(esets[[2]])),
                    sprintf("Male (n = %d)",   length(esets[[3]])))
  ed <- process_data(Venn(esets))
  epoly <- venn_regionedge(ed); eedg <- venn_setedge(ed)
  esl <- venn_setlabel(ed);     erl <- venn_regionlabel(ed)
  etot <- sum(erl$count)
  erl$lab <- sprintf("%s\n(%.1f%%)", format(erl$count, big.mark = ","), 100 * erl$count / etot)
  erl$txt <- ifelse(erl$count > 0.45 * max(erl$count), "white", "grey12")
  esl$col <- vcol[seq_len(nrow(esl))]
  eedg$col <- vcol[as.integer(factor(eedg$id))]
  ecx <- mean(range(epoly$X)); ecy <- mean(range(epoly$Y))
  esl$X <- ecx + (esl$X - ecx) * 1.10; esl$Y <- ecy + (esl$Y - ecy) * 1.18
  pev <- ggplot() +
    geom_polygon(data = epoly, aes(X, Y, group = id, fill = count), alpha = 0.9) +
    geom_path(data = eedg, aes(X, Y, group = id, color = col), linewidth = 1.3, show.legend = FALSE) +
    geom_text(data = esl, aes(X, Y, label = name, color = col), fontface = "bold", size = 5.2, show.legend = FALSE) +
    geom_text(data = erl, aes(X, Y, label = lab), size = 4, fontface = "bold", color = erl$txt, lineheight = 0.9) +
    scale_fill_gradientn(colours = c("#F7FBFF","#C6DBEF","#6BAED6","#2171B5","#08306B"),
                         trans = "sqrt", name = "enriched\nterms / region") +
    scale_color_identity() + coord_equal(clip = "off") +
    theme_void(base_size = 13) +
    theme(plot.margin = margin(34, 26, 34, 26), legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9))
  save2(pev, "fig_dge_enriched_terms_venn", 9.5, 5)
  cat(sprintf("  enriched-terms venn: Female=%d Male=%d All=%d\n",
              length(es$Female), length(es$Male), length(es$All)))
} else {
  cat("  [skip] fig_dge_enriched_terms_venn - run 03_dge.R (Part 2) first to make dge_enrich.rds\n")
}

# 7. Dotplots of functional enrichment (Female / Male / All), faceted by group [needs dge_enrich.rds]
if (file.exists(enrich_path)) {
  er   <- readRDS(enrich_path)$res
  etab <- rbindlist(lapply(names(er), function(g) er[[g]]$table), fill = TRUE)
  etab[, ratio := vapply(strsplit(GeneRatio, "/"),
                         function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))]
  etab[, group := factor(group, levels = c("Female", "Male", "All"))]

  dotplot_enrich <- function(dat, topn, ttl, fname, w, h) {
    d <- dat[order(p.adjust)][, head(.SD, topn), by = group]
    if (!nrow(d)) { cat(sprintf("  [skip] %s - no terms\n", fname)); return(invisible()) }
    ord <- d[, .(mp = min(p.adjust)), by = Description][order(-mp), Description]
    d[, Description := factor(Description, levels = ord)]
    p <- ggplot(d, aes(ratio, Description, size = Count, colour = p.adjust)) +
      geom_point() +
      facet_grid(group ~ ., scales = "free_y", space = "free_y") +
      scale_colour_gradient(low = "#B2182B", high = "#4575B4", name = "adj. p",
                            labels = function(x) format(x, scientific = FALSE)) +
      scale_size_continuous(range = c(2, 7), name = "gene count") +
      labs(title = ttl, x = "Gene ratio", y = NULL) +
      theme_bw(base_size = 11) +
      theme(strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8))
    save2(p, fname, w, h)
    cat(sprintf("  %s: %d terms across groups\n", fname, nrow(d)))
  }

  dotplot_enrich(etab, 8, "Functional enrichment of DEGs (top GO/KEGG terms per group)",
                 "fig_dge_enrich_dotplot", 9.5, 5)
  dotplot_enrich(etab[source == "KEGG"], 10, "KEGG pathway enrichment of DEGs",
                 "fig_dge_kegg_dotplot", 9.5, 5)
} else {
  cat("  [skip] enrichment dotplots - run 03_dge.R (Part 2) first to make dge_enrich.rds\n")
}

# 8. KEGG GSEA dotplot (Female / Male / All): NES > 0 activated, NES < 0 suppressed in RA [needs dge_gsea.rds]
gsea_path <- file.path(proc, "dge_gsea.rds")
if (file.exists(gsea_path) && nrow(as.data.table(readRDS(gsea_path)))) {
  gs <- as.data.table(readRDS(gsea_path))
  gs[, group := factor(group, levels = c("Female", "Male", "All"))]
  gtop <- gs[order(p.adjust)][, head(.SD, 12), by = group]
  ord  <- gtop[, .(mn = mean(NES)), by = Description][order(mn), Description]
  gtop[, Description := factor(Description, levels = ord)]
  pg <- ggplot(gtop, aes(NES, Description, size = setSize, colour = p.adjust)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_point() +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#B2182B", high = "#4575B4", name = "adj. p",
                          labels = function(x) format(x, scientific = FALSE)) +
    scale_size_continuous(range = c(2, 7), name = "set size") +
    labs(title = "KEGG GSEA: RA vs control (genes ranked by t-statistic)",
         subtitle = "NES > 0 = activated in RA, NES < 0 = suppressed in RA",
         x = "Normalised enrichment score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8))
  save2(pg, "fig_dge_kegg_gsea", 9.5, 5)
  cat(sprintf("  fig_dge_kegg_gsea: %d pathways across groups\n", nrow(gtop)))
} else {
  cat("  [skip] fig_dge_kegg_gsea - run 03_dge.R (Part 3) first to make dge_gsea.rds\n")
}

# 9. Classic GSEA running-score plots, top pathways per sex (enrichplot::gseaplot2) [dge_gsea_obj.rds]
gobj_path <- file.path(proc, "dge_gsea_obj.rds")
if (file.exists(gobj_path)) {
  gobjs <- readRDS(gobj_path)
  for (comp in c("Female", "Male")) {
    gk <- gobjs[[comp]]
    if (is.null(gk) || !nrow(as.data.frame(gk))) {
      cat(sprintf("  [skip] gsea running plot %s - no pathways\n", comp)); next }
    d   <- as.data.frame(gk)
    top <- d$ID[order(d$p.adjust)][seq_len(min(5, nrow(d)))]    # top 5 pathways by adj.P
    gk2 <- gk                                                   # show adj.P in the legend labels
    gk2@result$Description <- sprintf("%s  (p.adj = %s)", gk@result$Description,
                                      format(signif(gk@result$p.adjust, 3), scientific = FALSE, trim = TRUE))
    pr  <- gseaplot2(gk2, geneSetID = top, pvalue_table = FALSE, # built-in table dropped: overlapped curves
                     base_size = 11, ES_geom = "line",
                     title = sprintf("%s: top 5 KEGG GSEA pathways (RA vs control)", comp))
    pr[[1]] <- pr[[1]] +                                        # single legend (with adj.P) on the right
      theme(legend.position = "right", legend.title = element_blank(),
            legend.text = element_text(size = 8)) +
      guides(colour = guide_legend(ncol = 1))
    pr[[2]] <- pr[[2]] + theme(legend.position = "none")        # drop duplicate tick/metric legends
    pr[[3]] <- pr[[3]] + theme(legend.position = "none")
    prg <- wrap_plots(pr, ncol = 1, heights = c(1.5, 0.4, 0.9)) +
      plot_layout(guides = "collect")
    fn <- sprintf("fig_dge_gsea_running_%s", tolower(comp))
    ggsave(file.path(fig, paste0(fn, ".png")), prg, width = 9.5, height = 5, dpi = 300)
    ggsave(file.path(fig, paste0(fn, ".pdf")), prg, width = 9.5, height = 5, device = cairo_pdf)
    cat(sprintf("  fig_dge_gsea_running_%s: %s\n", tolower(comp),
                paste(d$Description[match(top, d$ID)], collapse = " | ")))
  }
} else {
  cat("  [skip] gsea running plots - run 03_dge.R (Part 3) to make dge_gsea_obj.rds\n")
}

cat("\nWrote: fig_dge_volcano, fig_dge_volcano_sex_{byFC,byP}, fig_dge_counts, fig_dge_venn,\n",
    "       fig_dge_heatmap_{female,male,all}, DEG_{female,male,all}, fig_dge_enriched_terms_venn,\n",
    "       fig_dge_enrich_dotplot, fig_dge_kegg_dotplot, fig_dge_kegg_gsea,\n",
    "       fig_dge_gsea_running_{female,male}\nDONE\n")
