# Validate the separately versioned direction-safe STITCH exploratory overlay.
# Run from the repository root after build_stitch_directional.py.
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else "results/causal_modelling/cosmos_mouse_pkn_stitch_directional"
stopifnot(requireNamespace("cosmosR", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))
path <- file.path(output_dir, "pkn_exploratory_stitch_directional.csv")
stopifnot(file.exists(path))
pkn <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(names(pkn), c("source", "interaction", "target")))
stopifnot(nrow(pkn) > 0L, !anyNA(pkn), all(pkn$interaction %in% c(-1, 1)))
stopifnot(all(nzchar(pkn$source)), all(nzchar(pkn$target)))
stopifnot(!any(pkn$source == pkn$target))
stopifnot(!anyDuplicated(pkn[c("source", "target")]))
clean <- cosmosR::meta_network_cleanup(pkn)
canonical <- function(x) {
  x <- x[order(x$source, x$target, x$interaction), c("source", "interaction", "target")]
  rownames(x) <- NULL
  x$interaction <- as.numeric(x$interaction)
  x
}
stopifnot(nrow(clean) == nrow(pkn), identical(canonical(clean), canonical(pkn)))
saveRDS(pkn, sub("[.]csv$", ".rds", path), compress = "xz")
jsonlite::write_json(list(
  validation = list(file = path, edges = nrow(pkn), nodes = length(unique(c(pkn$source, pkn$target))),
                    valid_signed_graph = TRUE, cosmos_cleaner_unchanged = TRUE),
  cosmosR_version = as.character(utils::packageVersion("cosmosR")),
  R_version = R.version.string,
  no_omics_scoring_performed = TRUE
), file.path(output_dir, "r_validation.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat(basename(path), ":", nrow(pkn), "edges; cosmosR cleaner unchanged\n")
