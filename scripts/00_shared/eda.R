#!/usr/bin/env Rscript
# =============================================================================
# eda.R  —  TWO-DATASET EXPLORATORY DATA ANALYSIS (computation only)
#           GSE93272 (GPL570) + GSE110169 (GPL13667), whole blood.
# Goal: sex-stratified biomarkers in RHEUMATOID ARTHRITIS (SLE is NOT a target;
#       SLE samples are flagged but excluded from the RA-vs-control counts).
#
# Produces results/tables/eda_*.csv and data/processed/eda_results.rds
# (bundle consumed by eda_figure.R). No plotting here.
# Covers: sex x disease per dataset, disease/sex/other-metadata distributions,
#         MISSING-VALUE audit, and expression-level QC (scale + PCA).
# =============================================================================
suppressMessages({library(Biobase); library(data.table)})
options(stringsAsFactors = FALSE); set.seed(1234)
raw <- "data/raw"; proc <- "data/processed"; tab <- "results/tables"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)

MISSING <- c(NA, "", "NA", "N/A", "n/a", "null", "NULL", "--", "unknown")

# ---- harmonize disease group + sex from a raw GEO eset ----------------------
harmonize <- function(eset, name) {
  pd <- pData(eset); clin <- grep(":ch1$", colnames(pd), value = TRUE)
  discol <- grep("disease|status", clin, value = TRUE, ignore.case = TRUE)[1]
  sexcol <- grep("gender|sex",     clin, value = TRUE, ignore.case = TRUE)[1]
  idcol  <- grep("individual|patient|subject", clin, value = TRUE, ignore.case = TRUE)
  dis <- as.character(pd[[discol]])
  grp <- ifelse(grepl("sle|lupus", dis, ignore.case = TRUE), "SLE",
          ifelse(grepl("normal|healthy|control", dis, ignore.case = TRUE), "HC",
          ifelse(grepl("rheumatoid|(^|[^a-z])ra([^a-z]|$)", dis, ignore.case = TRUE), "RA", "other")))
  sx <- toupper(substr(gsub("[^A-Za-z]", "", as.character(pd[[sexcol]])), 1, 1))
  sx[!sx %in% c("F","M")] <- NA
  meta <- data.frame(sample = rownames(pd), dataset = name,
                     disease_raw = dis, group = grp, sex = sx,
                     id = if (length(idcol)) as.character(pd[[idcol]]) else rownames(pd))
  list(meta = meta, pd = pd, clin = clin, discol = discol, sexcol = sexcol, has_id = length(idcol) > 0)
}

# ---- missing-value audit over all clinical fields ---------------------------
miss_audit <- function(pd, clin, name) {
  rbindlist(lapply(clin, function(cn) {
    v <- as.character(pd[[cn]]); m <- v %in% MISSING | is.na(v)
    data.table(dataset = name, field = sub(":ch1$", "", cn),
               n_missing = sum(m), pct_missing = round(100 * mean(m), 1),
               n_levels = length(unique(v[!m])))
  }))
}

# ---- categorical metadata distributions (<=15 levels) -----------------------
cat_dists <- function(pd, clin, name) {
  out <- list()
  for (cn in clin) {
    v <- as.character(pd[[cn]]); v[v %in% MISSING] <- NA
    u <- unique(v[!is.na(v)])
    if (length(u) >= 2 && length(u) <= 15) {
      tb <- as.data.frame(table(value = v, useNA = "no"))
      tb$dataset <- name; tb$field <- sub(":ch1$", "", cn); out[[cn]] <- tb
    }
  }
  if (length(out)) rbindlist(out) else NULL
}

# ---- expression QC (scale + PCA on top-variance probes) ---------------------
# Two PCAs: (1) all samples, (2) RA + HC only, RECOMPUTED with SLE removed so the
# axes reflect RA-vs-control variation (not SLE-driven variance).
expr_qc <- function(eset, name, grp) {
  ex <- exprs(eset); ex <- ex[rowSums(is.na(ex)) == 0, , drop = FALSE]
  mx <- max(ex); logged <- mx < 100
  persample <- data.table(dataset = name, sample = colnames(ex),
                          median = apply(ex, 2, median), q25 = apply(ex, 2, quantile, .25),
                          q75 = apply(ex, 2, quantile, .75))
  vt <- apply(ex, 1, var); top <- names(sort(vt, decreasing = TRUE))[1:min(2000, nrow(ex))]
  g <- grp[colnames(ex)]
  do_pca <- function(cols) {
    m <- ex[top, cols, drop = FALSE]; m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
    p <- prcomp(t(m), scale. = TRUE); v <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
    list(scores = data.table(dataset = name, sample = cols, PC1 = p$x[,1], PC2 = p$x[,2],
                             group = g[cols]), var = v[1:2])
  }
  list(logged = logged, max = round(mx, 1), persample = persample,
       pca_all = do_pca(colnames(ex)),
       pca_rh  = do_pca(colnames(ex)[g %in% c("RA","HC")]))   # SLE removed
}

# ---- run for both datasets --------------------------------------------------
datasets <- list(GSE93272 = "GSE93272_raw.rds", GSE110169 = "GSE110169_raw.rds")
res <- list()
for (nm in names(datasets)) {
  message("== ", nm, " ==")
  eset <- readRDS(file.path(raw, datasets[[nm]]))
  h <- harmonize(eset, nm)
  meta <- h$meta
  # subject-level (dedup by id — collapses GSE93272 longitudinal repeats)
  meta_subj <- meta[!duplicated(meta$id), ]
  sxdis <- as.data.frame(table(group = meta_subj$group, sex = meta_subj$sex, useNA = "ifany"))
  sxdis$dataset <- nm
  res[[nm]] <- list(
    meta = meta, meta_subj = meta_subj,
    n_samples = nrow(meta), n_subjects = nrow(meta_subj),
    sxdis = sxdis,
    disease_tab = as.data.frame(table(group = meta_subj$group)),
    sex_tab = as.data.frame(table(sex = meta_subj$sex, useNA = "ifany")),
    miss = miss_audit(h$pd, h$clin, nm),
    cats = cat_dists(h$pd, h$clin, nm),
    qc = expr_qc(eset, nm, setNames(meta$group, meta$sample)))
  cat(sprintf("  samples=%d, subjects=%d, clinical fields=%d\n",
              nrow(meta), nrow(meta_subj), length(h$clin)))
  cat("  sex x disease (subjects):\n"); print(table(meta_subj$group, meta_subj$sex, useNA="ifany"))
}

# ---- write CSV summaries ----------------------------------------------------
sxdis_all <- rbindlist(lapply(res, `[[`, "sxdis"))
fwrite(sxdis_all, file.path(tab, "eda_sex_by_disease.csv"))
miss_all <- rbindlist(lapply(res, `[[`, "miss"))
fwrite(miss_all[order(dataset, -pct_missing)], file.path(tab, "eda_missing_values.csv"))
overview <- rbindlist(lapply(names(res), function(nm) data.table(
  dataset = nm, samples = res[[nm]]$n_samples, subjects = res[[nm]]$n_subjects,
  RA = sum(res[[nm]]$meta_subj$group=="RA"), HC = sum(res[[nm]]$meta_subj$group=="HC"),
  SLE = sum(res[[nm]]$meta_subj$group=="SLE"),
  logged = res[[nm]]$qc$logged, max_expr = res[[nm]]$qc$max)))
fwrite(overview, file.path(tab, "eda_overview.csv"))
cat("\n== OVERVIEW ==\n"); print(overview)

saveRDS(res, file.path(proc, "eda_results.rds"))
cat("\nEDA computation done -> data/processed/eda_results.rds + results/tables/eda_*.csv\n")
