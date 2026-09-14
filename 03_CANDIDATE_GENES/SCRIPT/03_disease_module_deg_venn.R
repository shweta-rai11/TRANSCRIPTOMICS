#!/usr/bin/env Rscript
# Directional Venn diagrams: each disease module's overlap with UP-in-RA and DOWN-in-RA DEGs, per sex. Disease modules read from the analysis object, never hardcoded by colour.
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

# nudge long two-line set-name labels away from the anchor so they don't cross into the ellipse
fix_label_overlap <- function(p) {
  is_setname_layer <- vapply(p$layers, function(l)
    inherits(l$geom, "GeomText") && "name" %in% names(l$data), logical(1))
  i <- which(is_setname_layer)[1]
  if (!is.na(i)) {
    d <- p$layers[[i]]$data
    d$hjust <- ifelse(d$X <= mean(range(d$X)), 1, 0)
    d$X <- d$X + ifelse(d$hjust == 1, -0.3, 0.3)
    p$layers[[i]]$data <- d
    p$layers[[i]]$mapping <- modifyList(p$layers[[i]]$mapping,
                                        aes(hjust = .data$hjust))
  }
  p
}

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
  p <- fix_label_overlap(p)
  ggsave(file.path(fig, file), p, width = 6.6, height = 5.4, dpi = 300)
}

# fill colour per module: matches the module's WGCNA colour, dark=up/light=down
venn_hi <- function(m, dirn) {
  if (m == "yellow") {
    if (dirn == "up") "#F9A825" else "#FFF59D"
  } else if (m == "brown") {
    if (dirn == "up") "#6D4C41" else "#D7CCC8"
  } else {
    if (dirn == "up") "#E57373" else "#64B5F6"
  }
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
              venn_hi(m, "up"))
    save_venn(setNames(list(mg, down), c(sprintf("%s module\n(cor=%+.2f)", m, r),
                                         paste0(sx, " DEG_down"))),
              sprintf("fig_diseasemod_venn_%s_%s_down.png", m, tolower(sx)),
              sprintf("%s | %s module n DEG_down = %d%s",
                      sx, m, length(int_down),
                      if (consist == "down") "   [direction-consistent]" else ""),
              venn_hi(m, "down"))

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
