#!/usr/bin/env Rscript
# =============================================================================
# 29_figure_final_roc_synovium.R -- ROC curves for the FINAL recommended
# MR-anchored elastic-net panels (16_model_training_final_panel.R), per sex,
# with the External-blood test swapped for a synovium test (GSE89408, RA vs
# Normal): Train (nested CV, honest) | Internal test | Test synovium.
#
# The LOCKED coefficients saved in mr_final_objects.rds are applied, UNREFIT,
# to per-dataset z-scored synovium expression -- the same transfer rule used
# for internal/external blood in 16_ (missing gene -> z=0, fixed direction
# "<"). Only the coefficients transfer; the model is never refit on synovium.
#
# Outputs:
#   results/tables/mr_final_panel_synovium_auc.csv
#   results/figures/new/fig_mr_final_roc_{female,male}.png
# =============================================================================
suppressMessages({library(data.table); library(pROC); library(ggplot2)})
proc <- "data/processed"; tab <- "results/tables"; fig <- "results/figures/new"

obj <- readRDS(file.path(proc, "new", "mr_final_objects.rds"))
v   <- readRDS(file.path(proc, "new", "val_synovium.rds"))
genesF <- fread(file.path(tab, "FS_input_female.csv"))$gene   # same MR-prioritised set as 16_
genesM <- fread(file.path(tab, "FS_input_male.csv"))$gene
CAND <- list(F = genesF, M = genesM)

zrows <- function(M) t(apply(M, 1, function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}))
auc_ci <- function(r) {
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < 20) { set.seed(1234); suppressWarnings(as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000))) }
        else as.numeric(ci.auc(r))
  c(auc = as.numeric(auc(r)), ci[c(1, 3)])
}
fmt <- function(ci, n) {
  s <- sprintf("%.3f (%.3f-%.3f) [n=%d]", ci[1], ci[2], ci[3], n)
  if (ci[1] >= 0.999) s <- paste0(s, " SEPARATION")
  s
}

# ---- score the LOCKED blood-trained model on synovium (coefficients only) --
score_synovium <- function(sx) {
  sexlab <- if (sx == "F") "Female" else "Male"
  genes  <- unique(CAND[[sx]])
  coefs  <- obj$locked[[sx]]$coefs
  idx <- which(v$sex == sx); y <- factor(v$grp[idx], levels = c("Normal", "RA"))
  present <- genes[genes %in% rownames(v$logcpm)]
  Z <- zrows(v$logcpm[present, idx, drop = FALSE])
  Zdf <- as.data.frame(t(Z)); colnames(Zdf) <- make.names(present)
  for (g in setdiff(make.names(genes), colnames(Zdf))) Zdf[[g]] <- 0   # missing gene -> z=0
  Zdf <- Zdf[, make.names(genes), drop = FALSE]

  b0 <- coefs$coef[coefs$term == "(Intercept)"]
  bterm <- setNames(coefs$coef[coefs$term != "(Intercept)"], coefs$term[coefs$term != "(Intercept)"])
  bterm <- bterm[make.names(genes)]; bterm[is.na(bterm)] <- 0

  lp <- b0 + as.matrix(Zdf) %*% bterm
  p  <- as.numeric(1 / (1 + exp(-lp)))
  r  <- roc(y, p, levels = c("Normal", "RA"), direction = "<", quiet = TRUE)
  list(sexlab = sexlab, roc = r, ci = auc_ci(r), n = length(y))
}

resF <- score_synovium("F"); resM <- score_synovium("M")

out <- data.table(sex = c("Female", "Male"), n_synovium = c(resF$n, resM$n),
                   synovium_test = c(fmt(resF$ci, resF$n), fmt(resM$ci, resM$n)))
fwrite(out, file.path(tab, "mr_final_panel_synovium_auc.csv"))
cat("wrote mr_final_panel_synovium_auc.csv\n"); print(out)

# ---- figure: Train (nested CV) + Internal test + Test synovium -------------
roc0 <- as.data.table(obj$roc)
keep <- roc0[dataset %in% c("Train (nested CV)", "Internal test")]
syn_rc <- rbindlist(list(
  data.table(sex = "Female", dataset = "Test synovium", sens = resF$roc$sensitivities,
             spec = resF$roc$specificities, auc = resF$ci[1]),
  data.table(sex = "Male", dataset = "Test synovium", sens = resM$roc$sensitivities,
             spec = resM$roc$specificities, auc = resM$ci[1])))
roc_all <- rbindlist(list(keep, syn_rc), use.names = TRUE)

ord <- c("Train (nested CV)", "Internal test", "Test synovium")
roc_all[, dataset := factor(dataset, levels = ord)]
pal <- c("Train (nested CV)" = "#1b6ca8", "Internal test" = "#e08214", "Test synovium" = "#c0392b")

summ <- as.data.table(obj$summary)
auc_lab <- function(sexlab) {
  s <- summ[sex == sexlab]; syn <- out[sex == sexlab]$synovium_test
  strip_sep <- function(x) sub(" SEPARATION$", "", x)
  setNames(c(paste0("Train (nested CV): AUC ", sub(" \\(", "\n(", strip_sep(s$nested_CV))),
             paste0("Internal test: AUC ",     sub(" \\(", "\n(", strip_sep(s$internal_test))),
             paste0("Test synovium: AUC ",     sub(" \\(", "\n(", strip_sep(syn)))), ord)
}

plot_sex <- function(sexlab) {
  d <- roc_all[sex == sexlab]
  labs <- auc_lab(sexlab)
  d[, lab := labs[as.character(dataset)]]
  d[, lab := factor(lab, levels = labs[ord])]
  ggplot(d, aes(x = 1 - spec, y = sens, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
    geom_path(linewidth = 1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_colour_manual(values = setNames(unname(pal[ord]), labs[ord]), name = NULL) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = "grey80"),
          legend.margin = margin(4, 6, 4, 6),
          legend.key.height = unit(2, "lines"),
          legend.key.spacing.y = unit(4, "pt"),
          legend.text = element_text(size = 9, lineheight = 0.95),
          panel.grid.minor = element_blank())
}

for (sx in c("Female", "Male")) {
  g <- plot_sex(sx); tag <- tolower(sx)
  ggsave(file.path(fig, sprintf("fig_mr_final_roc_%s.png", tag)), g, width = 6, height = 6, dpi = 300)
  cat(sprintf("wrote fig_mr_final_roc_%s.png\n", tag))
}
