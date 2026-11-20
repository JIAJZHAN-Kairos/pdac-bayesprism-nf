#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BayesPrism))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: run_bayesprism_source.R SOURCE MIXTURE_RDS REFERENCE_RDS OUTPUT_DIR N_CORES")
}

source_name <- args[[1]]
mixture_path <- normalizePath(args[[2]], mustWork = TRUE)
reference_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
n_cores <- as.integer(args[[5]])
if (is.na(n_cores) || n_cores < 1L) stop("N_CORES must be positive")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_output_files <- c(
  "initial_state_fraction.tsv", "initial_type_fraction.tsv", "final_type_fraction.tsv",
  "initial_state_cv.tsv", "final_type_cv.tsv", "bayesprism_full_result.rds", "run_receipt.tsv"
)
unlink(file.path(output_dir, expected_output_files))

write_matrix <- function(x, path) {
  output <- data.frame(sample_id = rownames(x), x, check.names = FALSE)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

start_time <- Sys.time()
reference <- readRDS(reference_path)
mixture <- readRDS(mixture_path)

if (!identical(reference$key, "Malignant")) stop("Reference malignant key is not Malignant")
if (is.null(rownames(mixture)) || anyDuplicated(rownames(mixture))) {
  stop("Mixture sample IDs are missing or duplicated")
}
if (any(!is.finite(mixture)) || any(mixture < 0)) stop("Mixture contains invalid values")

prism <- new.prism(
  reference = reference$sc.dat,
  input.type = "count.matrix",
  cell.type.labels = reference$cell.type.labels,
  cell.state.labels = reference$cell.state.labels,
  key = reference$key,
  mixture = mixture,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

result <- run.prism(prism = prism, n.cores = n_cores)
saveRDS(result, file.path(output_dir, "bayesprism_full_result.rds"), compress = "gzip")

initial_state <- get.fraction(result, which.theta = "first", state.or.type = "state")
initial_type <- get.fraction(result, which.theta = "first", state.or.type = "type")
final_type <- get.fraction(result, which.theta = "final", state.or.type = "type")
initial_state_cv <- result@posterior.initial.cellState@theta.cv
initial_type_cv <- result@posterior.initial.cellType@theta.cv
final_type_cv <- result@posterior.theta_f@theta.cv
initial_type_cv_available <- identical(dim(initial_type_cv), dim(initial_type)) &&
  identical(dimnames(initial_type_cv), dimnames(initial_type))

outputs <- list(
  initial_state_fraction = initial_state,
  initial_type_fraction = initial_type,
  final_type_fraction = final_type,
  initial_state_cv = initial_state_cv,
  final_type_cv = final_type_cv
)

nonfinite_cv_counts <- integer(0)
for (output_name in names(outputs)) {
  x <- outputs[[output_name]]
  if (!identical(rownames(x), rownames(mixture))) {
    stop(output_name, " sample IDs or order do not match the input mixture")
  }
  if (grepl("fraction$", output_name)) {
    if (any(!is.finite(x))) stop(output_name, " contains non-finite values")
    if (any(x < 0)) stop(output_name, " contains negative fractions")
    if (any(abs(rowSums(x) - 1) > 1e-6)) stop(output_name, " rows do not sum to one")
  } else {
    nonfinite_cv_counts[[output_name]] <- sum(!is.finite(x))
    x[!is.finite(x)] <- NA_real_
    if (any(x < 0, na.rm = TRUE)) stop(output_name, " contains negative finite CV values")
  }
  write_matrix(x, file.path(output_dir, paste0(output_name, ".tsv")))
}

end_time <- Sys.time()
receipt <- data.frame(
  metric = c(
    "rna_source", "bayesprism_version", "reference_path", "mixture_path",
    "samples", "genes_in_input_mixture", "genes_used_by_prism", "cell_types",
    "cell_states", "n_cores", "start_time", "end_time", "elapsed_minutes",
    "maximum_final_fraction_sum_error", "initial_type_cv_available",
    "nonfinite_initial_state_cv", "nonfinite_final_type_cv"
  ),
  value = c(
    source_name, as.character(packageVersion("BayesPrism")), reference_path, mixture_path,
    nrow(mixture), ncol(mixture), ncol(result@prism@mixture), ncol(final_type),
    ncol(initial_state), n_cores,
    format(start_time, tz = "UTC", usetz = TRUE),
    format(end_time, tz = "UTC", usetz = TRUE),
    as.numeric(difftime(end_time, start_time, units = "mins")),
    max(abs(rowSums(final_type) - 1)), initial_type_cv_available,
    nonfinite_cv_counts[["initial_state_cv"]],
    nonfinite_cv_counts[["final_type_cv"]]
  ),
  stringsAsFactors = FALSE
)
write.table(
  receipt,
  file.path(output_dir, "run_receipt.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Completed BayesPrism source:", source_name, "\n")
print(receipt)
