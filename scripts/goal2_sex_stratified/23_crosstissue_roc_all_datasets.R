#!/usr/bin/env Rscript
# One figure per sex overlaying panel ROC across all four datasets: train (10-fold CV), internal test, external blood, and external synovium.
suppressMessages({library(data.table); library(pROC); library(ggplot2)})
proc <- "data/processed"; tabN <- "results/tables"; figN <- "results/figures/new"

ro <- readRDS(file.path(proc, "new", "mr_roc_objects.rds"))$roc   # named list of panel/gene ROC coords
ml <- readRDS(file.path(proc, "new", "ml_features.rds"))
panels <- list(Female = ml$female$consensus, Male = ml$male$consensus)

# ---- blood/train/internal panel ROC coords from saved store -----------------
want <- c("Train (10-fold CV)", "Internal test", "External blood")
grab <- function(sexlab) {
  keys <- names(ro)
  rbindlist(lapply(want, function(ds) {
    k <- keys[keys == sprintf("panel %s %s", sexlab, ds)]
    if (!length(k)) return(NULL)
    d <- as.data.table(ro[[k]]); d[, dataset := ds]; d[, .(dataset, sens, spec, auc)]
  }))
}

# ---- synovium panel ROC (compute here) --------------------------------------
v <- readRDS(file.path(proc, "new", "val_synovium.rds"))
syn_roc <- function(sexlab) {   # 10-fold CV out-of-fold ROC (honest, not apparent)
  genes <- panels[[sexlab]]; genes <- genes[genes %in% rownames(v$logcpm)]
  idx <- which(v$sex == ifelse(sexlab=="Female","F","M"))
  y <- factor(v$grp[idx], levels = c("Normal","RA"))
  Z <- scale(t(v$logcpm[genes, idx, drop=FALSE])); Z[is.na(Z)] <- 0; colnames(Z) <- make.names(genes)
  df <- data.frame(y=y, Z, check.names=FALSE)
  set.seed(1234); fold <- sample(rep_len(1:10, length(y))); p <- rep(NA_real_, length(y))
  for (k in 1:10) { tr <- fold!=k; te <- fold==k
    if (length(unique(y[tr]))<2) next
    f <- suppressWarnings(glm(y ~ ., df[tr,,drop=FALSE], family=binomial))
    p[te] <- as.numeric(predict(f, df[te,,drop=FALSE], type="response")) }
  ok <- !is.na(p)
  r <- roc(y[ok], p[ok], levels=c("Normal","RA"), direction="<", quiet=TRUE)
  data.table(dataset = "External synovium", sens = r$sensitivities, spec = r$specificities, auc = as.numeric(auc(r)))
}

ORDER <- c("Train (10-fold CV)", "Internal test", "External blood", "External synovium")
pal <- c("Train (10-fold CV)"="#1b6ca8", "Internal test"="#e08214",
         "External blood"="#4d9221", "External synovium"="#7B3294")

build <- function(sexlab) {
  d <- rbindlist(list(grab(sexlab), syn_roc(sexlab)))
  d[, sex := sexlab]
  au <- d[, .(auc = auc[1]), by = dataset]
  d[, lab := factor(sprintf("%s (AUC=%.2f)", dataset, au$auc[match(dataset, au$dataset)]),
                    levels = sprintf("%s (AUC=%.2f)", ORDER, au$auc[match(ORDER, au$dataset)]))]
  d
}
dd <- rbindlist(lapply(c("Female","Male"), build))
dd[, sex := factor(sex, levels = c("Female","Male"))]

# AUC summary table
fwrite(unique(dd[, .(sex, dataset, AUC = round(auc,3))])[order(sex, factor(dataset, levels=ORDER))],
       file.path(tabN, "allvalidation_panel_auc.csv"))

# colour per dataset but labels carry AUC (per facet) -> map by dataset colour
dd[, dcol := pal[dataset]]
g <- ggplot(dd, aes(1 - spec, sens, group = lab, colour = dataset)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey75") +
  geom_path(linewidth=1) +
  facet_wrap(~ sex) +
  scale_colour_manual(values = pal, guide = "none") +
  geom_text(data = unique(dd[, .(sex, dataset, auc)])[, .SD[1], by=.(sex,dataset)][
              , .(sex, dataset, auc, lab = sprintf("%s: AUC %.2f", dataset, auc))][order(sex, factor(dataset,levels=ORDER))][
              , y := 0.30 - 0.07*(match(dataset,ORDER)-1), by=sex],
            aes(x = 0.55, y = y, label = lab, colour = dataset), hjust = 0, size = 3.2, inherit.aes = FALSE) +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0.01,0)) +
  coord_equal() +
  labs(x="1 - Specificity", y="Sensitivity") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(figN, "fig_allvalidation_roc.png"), g, width=11, height=6, dpi=300, bg="white")
ggsave(file.path(figN, "fig_allvalidation_roc.pdf"), g, width=11, height=6, bg="white")
cat("wrote fig_allvalidation_roc + allvalidation_panel_auc.csv\n")
print(unique(dd[, .(sex, dataset, AUC=round(auc,3))])[order(sex, factor(dataset, levels=ORDER))])
