#!/usr/bin/env Rscript
# Fig 4A: per-gene expression by group (sex x RA/HC) for the MR consensus genes, with bootstrap CI and Wilcoxon significance stars.
suppressMessages({library(data.table); library(ggplot2)})
set.seed(1234)
proc <- "data/processed"; tab <- "results/tables"; fig <- "results/figures/new"

o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr
meta <- as.data.table(o$meta)
ml <- readRDS("data/processed/new/ml_features.rds")
genesF <- ml$female$consensus; genesM <- ml$male$consensus
genes  <- c(genesF, genesM)
gene_order <- genes[genes %in% rownames(expr)]

# long table --------------------------------------------------------------
L <- rbindlist(lapply(gene_order, function(g) data.table(
  gene = g, sample = colnames(expr), value = as.numeric(expr[g, ]),
  sex = meta$sex[match(colnames(expr), meta$sample)],
  group = meta$group[match(colnames(expr), meta$sample)])))
L <- L[sex %in% c("F", "M") & group %in% c("HC", "RA")]
L[, sexf := ifelse(sex == "F", "Female", "Male")]
L[, grp := factor(paste(sexf, ifelse(group == "HC", "Control", "RA")),
                  levels = c("Female Control", "Female RA", "Male Control", "Male RA"))]
L[, gene := factor(gene, levels = gene_order)]

# median + bootstrap 95% CI per (gene, grp) -------------------------------
boot_ci <- function(x, B = 2000) {
  x <- x[!is.na(x)]; if (length(x) < 2) return(c(median(x), NA, NA))
  bm <- replicate(B, median(sample(x, replace = TRUE)))
  c(median(x), quantile(bm, c(0.025, 0.975)))
}
summ <- L[, { ci <- boot_ci(value)
  .(med = ci[1], lo = ci[2], hi = ci[3], n = .N) }, by = .(gene, grp)]

# Wilcoxon RA vs HC within each sex, per gene -> stars ---------------------
star <- function(p) ifelse(is.na(p), "ns", ifelse(p <= 1e-4, "****",
  ifelse(p <= 1e-3, "***", ifelse(p <= 1e-2, "**", ifelse(p <= 5e-2, "*", "ns")))))
stat <- L[, {
  p <- tryCatch(wilcox.test(value[group == "RA"], value[group == "HC"])$p.value,
                error = function(e) NA_real_)
  .(p = p, star = star(p),
    med_HC = median(value[group == "HC"]), med_RA = median(value[group == "RA"]),
    ymax = max(value)) }, by = .(gene, sexf)]
fwrite(stat[, .(gene, sex = sexf, p = signif(p, 3), star,
                med_HC = round(med_HC, 3), med_RA = round(med_RA, 3))],
       file.path(tab, "mr_gene_group_stats.csv"))

# annotation positions (star centred over each sex's HC/RA pair) -----------
xpos <- c(Female = 1.5, Male = 3.5)
gy <- L[, .(top = max(value), bot = min(value)), by = gene]
ann <- merge(stat[, .(gene, sexf, star)], gy, by = "gene")
ann[, x := xpos[sexf]][, y := top + 0.08 * (top - bot)]

pal <- c("Female Control" = "#F5A623", "Female RA" = "#C0392B",
         "Male Control" = "#5DADE2", "Male RA" = "#1E8449")

g <- ggplot(L, aes(x = grp, y = value)) +
  geom_jitter(aes(colour = grp), width = 0.18, height = 0, size = 0.9, alpha = 0.55) +
  geom_errorbar(data = summ, aes(x = grp, ymin = lo, ymax = hi), width = 0.28,
                linewidth = 0.6, colour = "grey20", inherit.aes = FALSE) +
  geom_point(data = summ, aes(x = grp, y = med), inherit.aes = FALSE,
             size = 1.8, colour = "black") +
  geom_text(data = ann, aes(x = x, y = y, label = star), inherit.aes = FALSE,
            size = 3.6, vjust = 0) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "mRNA expression level") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig, "fig_mr_gene_expression_groups.png"), g, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig, "fig_mr_gene_expression_groups.pdf"), g, width = 9, height = 6)
cat("wrote fig_mr_gene_expression_groups.{png,pdf}\n")
print(stat[, .(gene, sex = sexf, p = signif(p, 3), star)])
