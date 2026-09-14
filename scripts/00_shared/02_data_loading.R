#!/usr/bin/env Rscript
# R version 4.4.2
# Retrieves GEO datasets (GSE93272, GSE110169, GSE15573; TRAINING/VALIDATION blood), writes raw ExpressionSets to data/raw/*.rds, and prints phenotype summaries

suppressMessages({
  library(GEOquery)
  library(Biobase)
})
options(timeout = 3600)                       # large series matrices
Sys.setenv(VROOM_CONNECTION_SIZE = 5e6) # 5,000,000 bytes 

# loads the data
raw_dir <- "data/raw"
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

gses <- c("GSE93272", "GSE110169", "GSE15573")

get_one <- function(gse) {
  rds <- file.path(raw_dir, paste0(gse, "_raw.rds"))
  if (file.exists(rds)) {
    message("...", gse, "loading")
    return(readRDS(rds))
  }
  message("... downloading ", gse)
  i <- getGEO(gse, destdir = raw_dir, GSEMatrix = TRUE, getGPL = TRUE,
              AnnotGPL = TRUE)
  eset <- if (length(i) == 1) i[[1]] else i     # this is a list of ExpressionSets if multi-platform
  saveRDS(eset, rds)
  eset
}

esets <- setNames(lapply(gses, get_one), gses)
names(esets) <- gses

# metadata informations

metadata <- function(gse, x) {cat("\n", gse, "\n")
  
  if (is(x, "list")) {
    cat("multi-platform: ", length(x), " ExpressionSets\n")
    for (i in seq_along(x)) metadata(paste0(gse, "_part", i), x[[i]])
    return(invisible())
  }
  cat("Platform :", annotation(x), "\n")
  cat("Features :", nrow(x), "  Samples:", ncol(x), "\n")
  ex <- exprs(x)
  cat("Expr range: [", round(min(ex, na.rm = TRUE), 2), ",",
      round(max(ex, na.rm = TRUE), 2), "]  (log2 if max>100)\n")
  pd <- pData(x)
  
# show phenodata columns 
  vary <- names(pd)[sapply(pd, function(c) length(unique(c)) > 1 &
                                            length(unique(c)) <= 25)]
  cat("Variable pheno columns:\n")
  for (v in vary) {
    vals <- table(pd[[v]])
    cat("   -", v, "::", paste(names(vals), vals, sep="=", collapse=" | "), "\n")
  }
  cat("Feature-data columns:", paste(head(names(fData(x)), 20), collapse=", "), "\n")
}

for (g in gses) metadata(g, esets[[g]])

cat("\n\n gene-symbol column \n")
for (g in gses) {
  x <- esets[[g]]
  if (is(x, "list")) x <- x[[1]]
  fd <- fData(x)
  cat("\n", g, "(", annotation(x), ") fData head \n")
  print(utils::head(fd, 3))
}
cat("\ncompleted_retreiving_data\n")
