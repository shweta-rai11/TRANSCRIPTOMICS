#!/usr/bin/env Rscript
# R version 4.4.2
#
# Projects the frozen 30% internal-validation holdout onto the training scale, using normalization + ComBat parameters estimated on training only (leakage-free). Run once, after the panel/model are locked on the 70% training data.

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
# Recover the training quantile target (row-mean of column-sorted training matrix) and map holdout onto it; skip if training skipped quantile norm.
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
# Re-implementation of parametric empirical-Bayes ComBat so batch parameters can be frozen and applied out-of-sample
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

# Pick whichever batch definition (batch_full vs dataset) best reproduces the saved training ComBat output, then use it for the holdout
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
