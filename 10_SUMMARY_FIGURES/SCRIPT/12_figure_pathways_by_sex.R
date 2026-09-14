#!/usr/bin/env Rscript
# Composite figure: established RA/immune KEGG pathways compared between female and male RA (DEG ORA, GSEA, MR-prioritised layers).
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
tabN <- "results/tables"; figN <- "results/figures/new"
dir.create(figN, showWarnings = FALSE, recursive = TRUE)

sexcol <- c(Female = "#C0392B", Male = "#1F5FA8")
FDR    <- 0.05

# 1. read the three layers
ora <- rbindlist(list(
  fread(file.path(tabN, "9_DEG_enrichment_female.csv")),
  fread(file.path(tabN, "10_DEG_enrichment_male.csv"))))[source == "KEGG"]
ora <- ora[, .(ID, Description, sex = group, n = Count, fdr = p.adjust)]

gsea <- fread(file.path(tabN, "gsea_kegg.csv"))[group %in% c("Female", "Male"),
  .(ID, Description, sex = group, NES, fdr = p.adjust)]

mr <- rbindlist(list(
  fread(file.path(tabN, "enrich_KEGG_female.csv"))[, sex := "Female"],
  fread(file.path(tabN, "enrich_KEGG_male.csv"))[, sex := "Male"]))
mr <- mr[, .(ID, Description, sex, ratio = GeneRatio, n = Count, fdr = p.adjust, genes = geneID)]

# 2. full cross-tabulation (all KEGG terms, all layers)
w <- function(d, lab, cols) {
  x <- dcast(d, ID + Description ~ sex, value.var = cols)
  setnames(x, setdiff(names(x), c("ID", "Description")),
           paste0(lab, "_", setdiff(names(x), c("ID", "Description"))))
  x }
all3 <- Reduce(function(a, b) merge(a, b, by = c("ID", "Description"), all = TRUE),
               list(w(ora, "DEG_ORA", c("n", "fdr")), w(gsea, "GSEA", c("NES", "fdr")),
                    w(mr, "MRgenes", c("n", "fdr", "genes"))))
pat <- function(f, m, lab) fifelse(!is.na(f) & is.na(m), paste0(lab, ":female_only"),
                fifelse(is.na(f) & !is.na(m), paste0(lab, ":male_only"),
                fifelse(!is.na(f) & !is.na(m), paste0(lab, ":both"), NA_character_)))
all3[, sex_pattern := {
  v <- c(pat(DEG_ORA_fdr_Female, DEG_ORA_fdr_Male, "DEG_ORA"),
         pat(GSEA_fdr_Female,    GSEA_fdr_Male,    "GSEA"),
         pat(MRgenes_fdr_Female, MRgenes_fdr_Male, "MR_genes"))
  paste(v[!is.na(v)], collapse = "; ") }, by = ID]
setorder(all3, ID)
fwrite(all3, file.path(tabN, "pathways_by_sex_all_layers.csv"))

# 3. curated "established RA/immune" KEGG set for the figure
cur <- fread(text = "
ID,label,axis
hsa04662,B cell receptor signalling,Adaptive
hsa04660,T cell receptor signalling,Adaptive
hsa04659,Th17 cell differentiation,Adaptive
hsa04658,Th1 / Th2 cell differentiation,Adaptive
hsa05235,PD-L1 / PD-1 checkpoint,Adaptive
hsa05340,Primary immunodeficiency,Adaptive
hsa04064,NF-kB signalling,Cytokine / signalling
hsa04668,TNF signalling,Cytokine / signalling
hsa04062,Chemokine signalling,Cytokine / signalling
hsa04010,MAPK signalling,Cytokine / signalling
hsa04068,FoxO signalling,Cytokine / signalling
hsa04115,p53 signalling,Cytokine / signalling
hsa04330,Notch signalling,Cytokine / signalling
hsa04370,VEGF signalling,Cytokine / signalling
hsa04933,AGE-RAGE signalling,Cytokine / signalling
hsa04621,NOD-like receptor signalling,Innate
hsa04623,Cytosolic DNA sensing,Innate
hsa04620,Toll-like receptor signalling,Innate
hsa04625,C-type lectin receptor signalling,Innate
hsa04666,Fc gamma R-mediated phagosome formation,Innate
hsa04664,Fc epsilon RI signalling,Innate
hsa04145,Phagocytosis,Innate
hsa04611,Platelet activation,Innate
hsa04380,Osteoclast differentiation,Bone
hsa04140,Autophagy,Cell fate / metabolism
hsa04210,Apoptosis,Cell fate / metabolism
hsa04216,Ferroptosis,Cell fate / metabolism
hsa04066,HIF-1 signalling,Cell fate / metabolism
hsa04150,mTOR signalling,Cell fate / metabolism
hsa00190,Oxidative phosphorylation,Cell fate / metabolism
hsa00020,Citrate (TCA) cycle,Cell fate / metabolism
hsa00010,Glycolysis / gluconeogenesis,Cell fate / metabolism
hsa04142,Lysosome,Cell fate / metabolism
hsa04141,Protein processing in ER,Cell fate / metabolism
hsa00061,Fatty acid biosynthesis,Cell fate / metabolism
hsa04935,Growth hormone synthesis / action,Hormonal
hsa04914,Progesterone-mediated oocyte maturation,Hormonal
hsa01522,Endocrine resistance,Hormonal
hsa04919,Thyroid hormone signalling,Hormonal
hsa05323,Rheumatoid arthritis,Adaptive
hsa04612,Antigen processing and presentation,Adaptive
hsa04640,Haematopoietic cell lineage,Adaptive
hsa04514,Cell adhesion molecules,Adaptive
hsa04672,Intestinal immune network for IgA,Adaptive
")
cur[, axis := factor(axis, c("Adaptive", "Cytokine / signalling", "Innate", "Bone",
                             "Cell fate / metabolism", "Hormonal"))]

# 4. panel data
# A: ORA, both sexes
oA <- ora[ID %in% cur$ID]
both_ids <- intersect(oA[sex == "Female", ID], oA[sex == "Male", ID])
pA <- merge(oA[ID %in% both_ids], cur, by = "ID")
pA[, nl := -log10(fdr)]
ordA <- pA[, .(m = mean(nl)), by = .(label, axis)][order(axis, m)]
pA[, label := factor(label, ordA$label)]

# B: ORA, exactly one sex
one_ids <- setdiff(oA$ID, both_ids)
pB <- merge(oA[ID %in% one_ids], cur, by = "ID")
pB[, nl := -log10(fdr) * fifelse(sex == "Female", -1, 1)]
# annotate GSEA where the same pathway is also displaced in that sex
pB <- merge(pB, gsea[, .(ID, sex, NES)], by = c("ID", "sex"), all.x = TRUE)
pB[, nes_lab := fifelse(is.na(NES), "", sprintf("NES %.2f", NES))]
ordB <- pB[order(axis, sex, -abs(nl))]
pB[, label := factor(label, rev(unique(ordB$label)))]

# C: GSEA, exactly one sex, restricted to curated set
gC <- gsea[ID %in% cur$ID]
gC_one <- gC[, .N, by = ID][N == 1, ID]
pC <- merge(gC[ID %in% gC_one], cur, by = "ID")
pC[, label := factor(label, pC[order(sex, NES)]$label)]

# D: MR-prioritised genes, curated RA-relevant terms, both sexes as a grid
dD <- c("hsa05323", "hsa04659", "hsa04658", "hsa04612", "hsa04640", "hsa04145",
        "hsa04514", "hsa04672")
pD <- CJ(ID = dD, sex = c("Female", "Male"))
pD <- merge(pD, mr[, .(ID, sex, n, ratio, fdr, genes)], by = c("ID", "sex"), all.x = TRUE)
pD <- merge(pD, cur[, .(ID, label)], by = "ID")
pD[, nl := -log10(fdr)]
pD[, txt := fifelse(is.na(fdr), "ns", sprintf("%s\nFDR %.1e", ratio, fdr))]
pD[, label := factor(label, rev(dD |> (\(i) cur[match(i, ID), label])()))]

fig_dat <- rbindlist(list(
  pA[, .(panel = "A_shared_ORA", ID, label, axis, sex, n, fdr, NES = NA_real_, genes = NA_character_)],
  pB[, .(panel = "B_sex_restricted_ORA", ID, label, axis, sex, n, fdr, NES, genes = NA_character_)],
  pC[, .(panel = "C_GSEA_one_sex", ID, label, axis, sex, n = NA_integer_, fdr, NES, genes = NA_character_)],
  pD[, .(panel = "D_MR_genes", ID, label, axis = NA, sex, n, fdr, NES = NA_real_, genes)]), fill = TRUE)
fwrite(fig_dat, file.path(tabN, "pathways_by_sex_figure_data.csv"))

# 5. plots
base <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "top", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        axis.text.y = element_text(size = 8.5))

gA <- ggplot(pA, aes(x = nl, y = label)) +
  geom_line(aes(group = label), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(colour = sex), size = 3) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = sexcol) +
  geom_vline(xintercept = -log10(FDR), linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  labs(x = expression(-log[10]~FDR~"(DEG over-representation)"), y = NULL,
       title = "A  Shared core: enriched in both sexes",
       subtitle = "KEGG pathways significant in female and male RA DEG sets") +
  base + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))

lim <- max(abs(pB$nl)) + 2.2
gB <- ggplot(pB, aes(x = nl, y = label, fill = sex)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = nes_lab, x = nl + sign(nl) * 0.15,
                hjust = fifelse(sex == "Female", 1, 0)), size = 2.6, colour = "grey20") +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = sexcol) +
  scale_x_continuous(limits = c(-lim, lim), labels = function(v) abs(v),
                     breaks = function(l) { b <- seq(0, floor(l[2]), by = 2); c(-rev(b[-1]), b) }) +
  labs(x = expression(-log[10]~FDR~"(DEG over-representation)"), y = NULL,
       title = "B  Sex-restricted: enriched in one sex only",
       subtitle = "Female-only left, male-only right; NES where GSEA also significant") +
  base + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))

gC <- ggplot(pC, aes(x = NES, y = label, colour = sex)) +
  geom_segment(aes(x = 0, xend = NES, yend = label), linewidth = 1.2) +
  geom_point(size = 3.2) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  facet_grid(sex ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = sexcol, guide = "none") +
  labs(x = "Normalised enrichment score (RA vs HC; < 0 = down in RA)", y = NULL,
       title = "C  Ranked GSEA: displaced in one sex only",
       subtitle = "Curated pathways with FDR < 0.05 in exactly one sex") +
  base

gD <- ggplot(pD, aes(x = sex, y = label)) +
  geom_tile(aes(fill = nl), colour = "white", linewidth = 1.5) +
  geom_text(aes(label = txt), size = 2.6, colour = "grey15", lineheight = 0.9) +
  scale_fill_gradient(low = "#DDEFEF", high = "#4FA3A8", na.value = "grey93",
                      name = expression(-log[10]~FDR), limits = c(0, NA)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "D  MR-prioritised causal genes",
       subtitle = "KEGG enrichment of each sex's candidate set (expressed-gene background)") +
  base + theme(legend.position = "right", legend.title = element_text(size = 8),
               panel.grid = element_blank(),
               axis.text.x = element_text(face = "bold", size = 9.5))

fig <- (gA | gB) / (gC | gD) + plot_layout(heights = c(1.35, 1), widths = c(1, 1.1))
ggsave(file.path(figN, "fig_pathways_by_sex.png"), fig, width = 15, height = 13, dpi = 300, bg = "white")
ggsave(file.path(figN, "fig_pathways_by_sex.pdf"), fig, width = 15, height = 13, bg = "white")

cat(sprintf("A: %d shared pathways | B: %d sex-restricted (F %d, M %d) | C: %d GSEA one-sex (F %d, M %d) | D: %d MR terms\n",
            uniqueN(pA$ID), uniqueN(pB$ID), pB[sex == "Female", .N], pB[sex == "Male", .N],
            uniqueN(pC$ID), pC[sex == "Female", .N], pC[sex == "Male", .N], uniqueN(pD$ID)))
cat("DONE\n")
