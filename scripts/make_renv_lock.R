#!/usr/bin/env Rscript
# =============================================================================
# make_renv_lock.R  —  build a valid renv.lock recording the EXACT package
# versions used to produce the reported results.
#
# Deliberately does NOT call renv::init(): that would relocate the project
# library and install an auto-activating .Rprofile, which risks breaking a
# working Bioconductor setup. This script only WRITES the lockfile; your
# existing library is untouched. `renv::restore()` can consume it later.
#
# Dependency closure is computed from LOCALLY installed metadata, so no network
# access is required.
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

top <- c("limma","sva","WGCNA","edgeR","glmnet","randomForest","e1071","caret","pROC",
         "TwoSampleMR","clusterProfiler","org.Hs.eg.db","data.table","ggplot2","rms",
         "preprocessCore","Biobase","IOBR","ggVennDiagram","ggrepel","patchwork",
         "AnnotationDbi","writexl","dplyr","igraph","pheatmap","magick","ggraph",
         "reshape2","scales","jsonlite")
top <- top[vapply(top, function(p) requireNamespace(p, quietly = TRUE), logical(1))]

ip <- installed.packages()

getdeps <- function(p) {
  if (!(p %in% rownames(ip))) return(character(0))
  f <- paste(ip[p, c("Depends", "Imports", "LinkingTo")], collapse = ",")
  f <- gsub("[(][^)]*[)]", "", f)          # strip version constraints
  f <- trimws(unlist(strsplit(f, ",")))
  f[nzchar(f) & f != "NA" & f != "R"]
}

closure <- top; frontier <- top
while (length(frontier)) {
  nx <- setdiff(unique(unlist(lapply(frontier, getdeps))), closure)
  nx <- nx[nx %in% rownames(ip)]
  closure <- c(closure, nx); frontier <- nx
}
base    <- rownames(installed.packages(priority = "base"))
closure <- sort(setdiff(unique(closure), base))

entry <- function(p) {
  d    <- packageDescription(p)
  bioc <- !is.null(d$biocViews)
  sprintf(paste0('    "%s": {\n',
                 '      "Package": "%s",\n',
                 '      "Version": "%s",\n',
                 '      "Source": "%s",\n',
                 '      "Repository": "%s"\n',
                 '    }'),
          p, p, as.character(packageVersion(p)),
          if (bioc) "Bioconductor" else "Repository",
          if (bioc) "Bioconductor" else "CRAN")
}

bv <- tryCatch(as.character(BiocManager::version()), error = function(e) "3.20")

lock <- sprintf(paste0(
  '{\n  "R": {\n    "Version": "%s",\n    "Repositories": [\n',
  '      { "Name": "CRAN", "URL": "https://cloud.r-project.org" },\n',
  '      { "Name": "BioCsoft", "URL": "https://bioconductor.org/packages/%s/bioc" },\n',
  '      { "Name": "BioCann", "URL": "https://bioconductor.org/packages/%s/data/annotation" }\n',
  '    ]\n  },\n  "Bioconductor": { "Version": "%s" },\n  "Packages": {\n%s\n  }\n}\n'),
  paste0(R.version$major, ".", R.version$minor), bv, bv, bv,
  paste(vapply(closure, entry, character(1)), collapse = ",\n"))

writeLines(lock, "renv.lock")
cat(sprintf("top-level packages: %d | full dependency closure: %d\n", length(top), length(closure)))
