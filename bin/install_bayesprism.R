#!/usr/bin/env Rscript
# Install BayesPrism 2.2.3 (Danko-Lab commit 19052052) and its dependencies into a task-local library.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1]]; ncpus <- as.integer(args[[2]])
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = BiocManager::repositories())
BiocManager::install(c("snowfall", "NMF", "gplots", "scran", "BiocParallel"), lib = lib, update = FALSE, ask = FALSE, Ncpus = ncpus)
commit <- "19052052a6f30833b27e2148294459b0ba2f923e"
tgz <- file.path(tempdir(), "bayesprism.tar.gz")
download.file(sprintf("https://github.com/Danko-Lab/BayesPrism/archive/%s.tar.gz", commit), tgz, mode = "wb")
untar(tgz, exdir = tempdir())
install.packages(file.path(tempdir(), paste0("BayesPrism-", commit), "BayesPrism"), repos = NULL, type = "source", lib = lib)
if (as.character(packageVersion("BayesPrism", lib.loc = lib)) != "2.2.3") stop("Unexpected BayesPrism version")
cat("BayesPrism", as.character(packageVersion("BayesPrism", lib.loc = lib)), "commit", commit, "\n")
