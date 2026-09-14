#!/usr/bin/env Rscript
# Clinical diagnostic-model validation of the sex-stratified MR consensus panels: nomogram, calibration, decision curve, and clinical impact, per sex.
suppressMessages({library(rms); library(data.table); library(magick)})
proc <- "data/processed"; figN <- "results/figures/new"; tabN <- "results/tables"
dir.create(figN, showWarnings = FALSE, recursive = TRUE)
GLOBAL_SEED <- 1234  # matches GLOBAL_SEED convention used throughout the goal2 pipeline scripts

o <- readRDS(file.path(proc, "combined_train.rds")); expr <- o$expr
meta <- as.data.table(o$meta)
ml <- readRDS("data/processed/new/ml_features.rds")
panels <- list(Female = ml$female$consensus, Male = ml$male$consensus)
sexcode <- c(Female = "F", Male = "M")

# net benefit (Vickers DCA): NB = TP/n - FP/n * pt/(1-pt) ---------------------
dca <- function(y, p, th) {
  n <- length(y); ev <- mean(y == 1)
  model <- sapply(th, function(pt) { pos <- p >= pt
    sum(pos & y == 1)/n - (sum(pos & y == 0)/n) * (pt/(1 - pt)) })
  all <- sapply(th, function(pt) ev - (1 - ev) * (pt/(1 - pt)))
  list(model = model, all = all)
}

for (sx in names(panels)) {
  genes <- panels[[sx]]; genes <- genes[genes %in% rownames(expr)]
  cols  <- meta$sample[meta$sex == sexcode[[sx]]]
  df <- as.data.frame(t(expr[genes, cols, drop = FALSE])); colnames(df) <- make.names(genes)
  df$RA <- as.integer(meta$group[match(cols, meta$sample)] == "RA")
  dd <- datadist(df); options(datadist = "dd")
  form <- as.formula(paste("RA ~", paste(make.names(genes), collapse = " + ")))
  # mild penalty stabilises the fit (male panel separates perfectly at n=38)
  pen <- if (sx == "Male") 5 else 0
  fit <- lrm(form, data = df, x = TRUE, y = TRUE, penalty = pen)
  p <- plogis(predict(fit)); ev <- mean(df$RA == 1)
  th <- seq(0.01, 0.99, 0.01); d <- dca(df$RA, p, th)
  fwrite(data.table(threshold = th, NB_panel = d$model, NB_all = d$all, NB_none = 0),
         file.path(tabN, sprintf("diag_dca_%s.csv", tolower(sx))))

  tag <- tolower(sx); tmp <- function(x) file.path(figN, sprintf(".tmp_%s_%s.png", tag, x))

  ## (A) nomogram -- its own device (does not compose inside layout)
  png(tmp("A"), width = 1500, height = 1250, res = 200); par(mar = c(2,1,3,1))
  plot(nomogram(fit, fun = plogis, funlabel = "Risk of RA",
                fun.at = c(0.05,0.1,0.3,0.5,0.7,0.9,0.99)))
  title(main = sprintf("A   %s nomogram (%d genes)", sx, length(genes)), adj = 0, cex.main = 1.3)
  dev.off()

  ## (B) calibration
  png(tmp("B"), width = 1300, height = 1250, res = 200); par(mar = c(4,4,3,1))
  set.seed(GLOBAL_SEED)  # calibrate() bootstraps internally and was previously unseeded -> non-reproducible MAE/MSE
  cal <- tryCatch(calibrate(fit, B = 200), error = function(e) NULL)
  if (!is.null(cal)) plot(cal, xlab = "Predicted probability", ylab = "Actual probability", subtitles = FALSE, main = "")
  else plot.new()
  title(main = "B   Calibration", adj = 0, cex.main = 1.3); dev.off()

  # secondary Cost:Benefit axis (cost:benefit = pt/(1-pt)), rmda convention
  cb_at  <- c(1/101, 0.2, 0.4, 0.6, 0.8, 100/101)
  cb_lab <- c("1:100", "1:4", "2:3", "3:2", "4:1", "100:1")
  add_cb <- function() { axis(1, at = cb_at, labels = cb_lab, line = 3.2, cex.axis = 0.85)
    mtext("Cost:Benefit Ratio", side = 1, line = 5.4) }

  ## (C) decision curve -- clip y so treat-all's dive doesn't squash the panel
  png(tmp("C"), width = 1300, height = 1350, res = 200); par(mar = c(7.5,4,3,1))
  ymax <- max(c(d$model, ev)) * 1.1
  plot(th, d$model, type = "l", col = "#C0392B", lwd = 2.5, ylim = c(-0.05, ymax),
       xlab = "", ylab = "Net benefit", main = "")
  mtext("Threshold probability", side = 1, line = 2.2)
  lines(th, d$all, col = "grey55", lwd = 1.5); abline(h = 0, col = "black", lwd = 1.2)
  add_cb()
  legend("topright", c("Panel genes", "Treat all", "Treat none"),
         col = c("#C0392B", "grey55", "black"), lwd = c(2.5, 1.5, 1.2), bty = "n")
  title(main = "C   Decision curve", adj = 0, cex.main = 1.3); dev.off()

  ## (D) clinical impact -- point estimates + bootstrap 95% confidence bands
  N <- 1000
  nhigh  <- sapply(th, function(pt) mean(p >= pt) * N)
  nevent <- sapply(th, function(pt) mean(p >= pt & df$RA == 1) * N)
  B <- 500; bh <- matrix(NA_real_, length(th), B); be <- matrix(NA_real_, length(th), B)
  for (b in seq_len(B)) {
    set.seed(b); idx <- sample(nrow(df), replace = TRUE)
    fb <- tryCatch(lrm(form, data = df[idx, ], penalty = pen), error = function(e) NULL)
    if (is.null(fb)) next
    pb <- plogis(predict(fb, newdata = df[idx, ])); yb <- df$RA[idx]
    bh[, b] <- sapply(th, function(pt) mean(pb >= pt) * N)
    be[, b] <- sapply(th, function(pt) mean(pb >= pt & yb == 1) * N)
  }
  q <- function(M, pr) apply(M, 1, quantile, pr, na.rm = TRUE)
  png(tmp("D"), width = 1300, height = 1350, res = 200); par(mar = c(7.5,4,3,1))
  plot(th, nhigh, type = "l", col = "#C0392B", lwd = 2.5, ylim = c(0, N),
       xlab = "", ylab = "Number high risk (per 1000)", main = "")
  mtext("High-risk threshold", side = 1, line = 2.2)
  lines(th, q(bh, .025), col = "#C0392B", lwd = 0.8); lines(th, q(bh, .975), col = "#C0392B", lwd = 0.8)
  lines(th, nevent, col = "#2E86C1", lwd = 2.5, lty = 2)
  lines(th, q(be, .025), col = "#2E86C1", lwd = 0.8, lty = 3); lines(th, q(be, .975), col = "#2E86C1", lwd = 0.8, lty = 3)
  add_cb()
  legend("topright", c("Number high risk", "Number high risk with RA", "95% CI (bootstrap)"),
         col = c("#C0392B", "#2E86C1", "grey40"), lwd = c(2.5, 2.5, 0.8), lty = c(1, 2, 1), bty = "n")
  title(main = "D   Clinical impact", adj = 0, cex.main = 1.3); dev.off()

  ## compose 2x2 (A|B over C|D) with magick
  ga <- image_read(tmp("A")); gb <- image_read(tmp("B"))
  gc <- image_read(tmp("C")); gd <- image_read(tmp("D"))
  h <- 1250; sc <- function(x) image_resize(x, sprintf("x%d", h))
  row1 <- image_append(c(sc(ga), sc(gb))); row2 <- image_append(c(sc(gc), sc(gd)))
  comp <- image_append(c(row1, row2), stack = TRUE)
  image_write(comp, file.path(figN, sprintf("fig_diag_validation_%s.png", tag)))
  image_write(image_convert(comp, format = "pdf"), file.path(figN, sprintf("fig_diag_validation_%s.pdf", tag)))
  file.remove(tmp("A"), tmp("B"), tmp("C"), tmp("D"))
  cat(sprintf("wrote fig_diag_validation_%s.png  (panel %d genes, n=%d, RA=%d)\n",
              tag, length(genes), nrow(df), sum(df$RA)))
}
