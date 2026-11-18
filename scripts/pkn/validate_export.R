# Validate all assembled PKN exports against the installed cosmosR cleaner.
# Run from repository root; --validate-only checks saved CSVs without writing.
# No omic input vectors or MOON scoring are constructed here.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/pkn/validate_export.R [--output=PKN_DIRECTORY] [--validate-only]\n")
  quit(save = "no")
}
if (any(!grepl("^--output=.+$|^--validate-only$", args)) ||
    anyDuplicated(sub("=.*$", "", args))) stop("Unknown or repeated option; use --help")
output_arg <- args[startsWith(args, "--output=")]
output_dir <- if (length(output_arg)) sub("^--output=", "", output_arg) else
  "results/causal_modelling/cosmos_mouse_pkn"
validate_only <- "--validate-only" %in% args
stopifnot(requireNamespace("cosmosR", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))
files <- sort(list.files(output_dir, pattern = "^pkn_.*[.]csv$", full.names = TRUE))
stopifnot(length(files) >= 3L)
checks <- lapply(files, function(path) {
  pkn <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(identical(names(pkn), c("source", "interaction", "target")))
  stopifnot(nrow(pkn) > 0L, !anyNA(pkn), all(pkn$interaction %in% c(-1, 1)))
  stopifnot(all(nzchar(pkn$source)), all(nzchar(pkn$target)))
  stopifnot(!any(pkn$source == pkn$target))
  stopifnot(!anyDuplicated(pkn[c("source", "target")]))
  clean <- cosmosR::meta_network_cleanup(pkn)
  # The exported graph must already satisfy the cleaner's exact contract.
  stopifnot(nrow(clean) == nrow(pkn))
  canonical <- function(x) {
    x <- x[order(x$source, x$target, x$interaction), c("source", "interaction", "target")]
    rownames(x) <- NULL
    x$interaction <- as.numeric(x$interaction)
    x
  }
  stopifnot(identical(canonical(clean), canonical(pkn)))
  if (!validate_only) saveRDS(pkn, sub("[.]csv$", ".rds", path), compress = "xz")
  cat(basename(path), ":", nrow(pkn), "edges; cleaner unchanged\n")
  list(file = path, edges = nrow(pkn), nodes = length(unique(c(pkn$source, pkn$target))),
       valid_signed_graph = TRUE, cosmos_cleaner_unchanged = TRUE)
})
if (!validate_only) {
jsonlite::write_json(list(
  validation = checks,
  cosmosR_version = as.character(utils::packageVersion("cosmosR")),
  R_version = R.version.string,
  no_omics_scoring_performed = TRUE
), file.path(output_dir, "r_validation.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
}
if (validate_only) cat("Validated", length(checks), "saved PKN exports; no files written.\n")
