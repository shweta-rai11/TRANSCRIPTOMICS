#!/usr/bin/env Rscript
# R version 4.4.2
#
# 03_apply_holdout.R
# ---------------------------------------------------------------------------
# PURPOSE
#   Project the FROZEN 30% internal-validation holdout onto the TRAINING scale,
#   using normalisation + batch-correction parameters ESTIMATED ON TRAINING ONLY.
#   No holdout sample contributes to any estimated parameter -> leakage-free.
#
# WHEN TO RUN THIS  (read before running)
#   Run ONCE, at the very end, AFTER the biomarker panel and ML model are fully
#   LOCKED on the 70% training data:
#     - differential expression / WGCNA / feature discovery  -> done on train
#     - the panel (which genes)                              -> chosen on train
#     - the model, hyper-parameters and decision threshold   -> fixed via CV WITHIN train
#   Until then the holdout stays sealed. Process it here, evaluate the locked model
#   ONCE, and report that number. Do NOT revisit the model after seeing the holdout.
#
# WHY frozen parameters (and not a fresh joint normalisation/ComBat)
#   Quantile normalisation learns a shared reference distribution across samples, and
#   ComBat estimates batch + covariate effects (it even uses the outcome labels via `mod`).
#   Re-estimating either on train+holdout together lets the holdout leak into the model's
#   input representation -> optimistic bias. We therefore FIX the training parameters and
#   only APPLY them to the holdout (an "addon"/frozen adjustment).
#
# REFERENCES
#   # Ambroise C, McLachlan GJ. Selection bias in gene extraction... PNAS 2002;99(10):6562-6566.
#   # Simon R, et al. Pitfalls in the use of DNA microarray data for classification. J Natl Cancer Inst 2003;95(1):14-18.
#   # Kaufman S, et al. Leakage in data mining. ACM TKDD 2012;6(4):Article 15.
#   # Bolstad BM, et al. A comparison of normalization methods... Bioinformatics 2003;19(2):185-193.  (quantile norm)
#   # Ritchie ME, et al. limma powers differential expression analyses. Nucleic Acids Res 2015;43(7):e47.
#   # Johnson WE, Li C, Rabinovic A. Adjusting batch effects... empirical Bayes. Biostatistics 2007;8(1):118-127.  (ComBat)
#   # Leek JT, et al. The sva package for removing batch effects. Bioinformatics 2012;28(6):882-883.
#   # Nygaard V, Rodland EA, Hovig E. Methodological variations in ...ComBat and their impact. Brief Bioinform 2016;17(1):29-39.  (caution: fit batch params on train, keep design balanced)
# ---------------------------------------------------------------------------

library(preprocessCore)   # normalize.quantiles.use.target : frozen quantile normalisation to a fixed target
library(data.table)

proc <- "data/processed"
tab  <- "results/tables"

train_obj   <- readRDS(file.path(proc, "combined_train.rds"))          # written by 02_normalize_batch.R (70% train)
holdout_obj <- readRDS(file.path(proc, "internal_val_holdout.rds"))    # 30% holdout, pre-norm, pre-ComBat

stopifnot(train_obj$role == "combined_training_70pct")

# Align holdout to the training gene order (both use the same `common` gene set).
genes <- rownames(train_obj$expr)
train_prenorm <- train_obj$expr_prenorm[genes, , drop = FALSE]
train_qnorm   <- train_obj$expr_qnorm[genes, , drop = FALSE]
train_combat  <- train_obj$expr[genes, , drop = FALSE]
train_meta    <- train_obj$meta

hold_prenorm  <- holdout_obj$expr_prenorm[genes, , drop = FALSE]
hold_meta     <- holdout_obj$meta

# ---- STEP 1: frozen quantile normalisation --------------------------------
# The training quantile reference = the common (averaged, sorted) distribution that
# limma::normalizeBetweenArrays produced on the training set. After quantile norm every
# training column shares this sorted distribution, so we recover the target as the
# row-mean of the column-sorted training matrix and map each HOLDOUT sample onto it.
# Skip if quantile normalisation was NOT applied in training (mirror the training decision).
if (isTRUE(train_obj$normalization_applied)) {
  target <- rowMeans(apply(train_qnorm, 2, sort))                       # frozen training quantile target
  hold_qnorm <- preprocessCore::normalize.quantiles.use.target(as.matrix(hold_prenorm), target)
  dimnames(hold_qnorm) <- dimnames(hold_prenorm)
  message("STEP 1: frozen quantile normalisation applied to holdout (target from training).")
} else {
  hold_qnorm <- as.matrix(hold_prenorm)
  message("STEP 1: quantile normalisation was NOT applied in training -> holdout left unchanged (consistent).")
}

# ---- parametric ComBat: fit on TRAIN, apply to HOLDOUT ---------------------
# Faithful re-implementation of the parametric empirical-Bayes ComBat estimation
# (Johnson et al. 2007) so the batch parameters can be FROZEN and applied out-of-sample.
# Priors / iterative posterior solution follow the sva reference implementation.
aprior <- function(d) { m <- mean(d); s2 <- var(d); (2 * s2 + m^2) / s2 }
bprior <- function(d) { m <- mean(d); s2 <- var(d); (m * s2 + m^3) / s2 }
postmean <- function(g.hat, g.bar, n, d.star, t2) (t2 * n * g.hat + d.star * g.bar) / (t2 * n + d.star)
postvar  <- function(sum2, n, a, b) (0.5 * sum2 + b) / (n / 2 + a - 1)
it.sol <- function(sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 1e-4) {
  n <- rowSums(!is.na(sdat)); g.old <- g.hat; d.old <- d.hat; change <- 1
  while (change > conv) {
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2  <- rowSums((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, na.rm = TRUE)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old, na.rm = TRUE)
    g.old <- g.new; d.old <- d.new
  }
  list(g.star = g.old, d.star = d.old)
}

# Fit standardisation + EB batch parameters on TRAINING data only.
fit_combat <- function(dat, batch, mod) {
  batch <- factor(batch)
  lev <- levels(batch)
  batches <- lapply(lev, function(b) which(batch == b))
  n.batch <- length(lev); n.batches <- sapply(batches, length); n.array <- sum(n.batches)
  batchmod <- model.matrix(~ -1 + batch)
  design <- cbind(batchmod, mod[, -1, drop = FALSE])                    # drop intercept of mod (batch dummies carry it)
  B.hat <- solve(crossprod(design), t(design) %*% t(dat))              # (params x genes)
  grand.mean <- as.vector(crossprod(n.batches / n.array, B.hat[1:n.batch, , drop = FALSE]))
  var.pooled <- as.vector(((dat - t(design %*% B.hat))^2) %*% rep(1 / n.array, n.array))
  cov.coef <- B.hat[(n.batch + 1):nrow(B.hat), , drop = FALSE]         # covariate (group/sex) coefficients
  # standardise
  stand.mean <- grand.mean + t(mod[, -1, drop = FALSE] %*% cov.coef)
  s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, n.array)))
  # batch EB parameters
  gamma.hat <- t(sapply(batches, function(ix) rowMeans(s.data[, ix, drop = FALSE])))
  delta.hat <- t(sapply(batches, function(ix) apply(s.data[, ix, drop = FALSE], 1, var)))
  gamma.bar <- rowMeans(gamma.hat); t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior); b.prior <- apply(delta.hat, 1, bprior)
  gamma.star <- delta.star <- matrix(NA, n.batch, nrow(dat))
  for (i in 1:n.batch) {
    tmp <- it.sol(s.data[, batches[[i]], drop = FALSE], gamma.hat[i, ], delta.hat[i, ],
                  gamma.bar[i], t2[i], a.prior[i], b.prior[i])
    gamma.star[i, ] <- tmp$g.star; delta.star[i, ] <- tmp$d.star
  }
  list(levels = lev, grand.mean = grand.mean, var.pooled = var.pooled,
       cov.coef = cov.coef, gamma.star = gamma.star, delta.star = delta.star)
}

# Apply FROZEN training parameters to new samples (holdout).
apply_combat <- function(par, dat, batch, mod) {
  batch <- as.character(batch)
  if (!all(batch %in% par$levels))
    stop("Holdout contains batch(es) absent from training: ",
         paste(setdiff(unique(batch), par$levels), collapse = ", "),
         " -> cannot freeze batch parameters for them.")
  stand.mean <- par$grand.mean + t(mod[, -1, drop = FALSE] %*% par$cov.coef)
  s.data <- (dat - stand.mean) / (sqrt(par$var.pooled) %*% t(rep(1, ncol(dat))))
  out <- s.data
  for (j in seq_len(ncol(dat))) {
    bi <- match(batch[j], par$levels)
    out[, j] <- (s.data[, j] - par$gamma.star[bi, ]) / sqrt(par$delta.star[bi, ])
  }
  out * (sqrt(par$var.pooled) %*% t(rep(1, ncol(dat)))) + stand.mean
}

# Build model matrices consistently for train and holdout (protect group + sex).
mk_mod <- function(m) {
  m$group <- factor(m$group, levels = c("HC", "RA"))
  m$sex   <- factor(m$sex,   levels = c("F", "M"))
  model.matrix(~ group + sex, data = m)
}
mod_train <- mk_mod(train_meta)
mod_hold  <- mk_mod(hold_meta)

# The training script uses batch = study+internal (batch_full) with a fallback to
# study-level batch. Pick whichever batch definition BEST reproduces the saved training
# ComBat output, then use that SAME definition to freeze-and-apply to the holdout.
recon_err <- function(batch_train) {
  par <- fit_combat(train_qnorm, batch_train, mod_train)
  rec <- apply_combat(par, train_qnorm, batch_train, mod_train)        # re-apply to train == should match saved combat
  list(par = par, err = max(abs(rec - train_combat)))
}
fit_full <- tryCatch(recon_err(train_meta$batch_full), error = function(e) NULL)
fit_study <- recon_err(train_meta$dataset)
use_full <- !is.null(fit_full) && fit_full$err <= fit_study$err
par     <- if (use_full) fit_full$par else fit_study$par
batch_h <- if (use_full) hold_meta$batch_full else hold_meta$dataset
message(sprintf("STEP 2: ComBat batch definition = %s (train reconstruction max|err| = %.3g).",
                if (use_full) "study+internal (batch_full)" else "study-level (dataset)",
                if (use_full) fit_full$err else fit_study$err))
if ((if (use_full) fit_full$err else fit_study$err) > 1e-3)
  warning("Train reconstruction error is larger than expected; verify the batch definition matches 02_normalize_batch.R.")

hold_combat <- apply_combat(par, hold_qnorm, batch_h, mod_hold)
dimnames(hold_combat) <- dimnames(hold_qnorm)

# ---- SAVE ------------------------------------------------------------------
saveRDS(list(expr = hold_combat,                       # holdout on the TRAINING scale (ready for the locked model)
             expr_qnorm = hold_qnorm,
             expr_prenorm = as.matrix(hold_prenorm),
             meta = hold_meta, genes = genes,
             role = "internal_validation_processed_30pct",
             provenance = paste("Frozen from training:",
                                "quantile target = rowMeans(sort(train_qnorm));",
                                if (use_full) "ComBat batch = study+internal;" else "ComBat batch = study-level;",
                                "parameters estimated on 70% training only (leakage-free)."),
             batch_definition = if (use_full) "batch_full" else "dataset"),
        file.path(proc, "internal_val_holdout_processed.rds"))
fwrite(data.table(gene = rownames(hold_combat), as.data.table(hold_combat)),
       file.path(proc, "internal_val_holdout_processed.csv.gz"))

cat("\n===== HOLDOUT PROCESSED (frozen from training) =====\n")
cat(sprintf("Holdout: %d samples x %d genes | RA=%d HC=%d | F=%d M=%d\n",
            ncol(hold_combat), nrow(hold_combat),
            sum(hold_meta$group == "RA"), sum(hold_meta$group == "HC"),
            sum(hold_meta$sex == "F"), sum(hold_meta$sex == "M")))
cat("Saved: data/processed/internal_val_holdout_processed.rds (+ .csv.gz)\n")
cat("\nNEXT: load the LOCKED panel/model and evaluate on $expr ONCE. Do not modify the model afterwards.\n")
