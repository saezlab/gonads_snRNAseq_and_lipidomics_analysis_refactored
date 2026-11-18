# Validate and export the separate measured-species PKNs; no MOON scoring.
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else
  "results/causal_modelling/cosmos_mouse_pkn_expanded_species"
stopifnot(requireNamespace("cosmosR", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE))
canonical <- function(x) {
  x <- x[order(x$source, x$target, x$interaction),
         c("source", "interaction", "target")]
  rownames(x) <- NULL
  x$interaction <- as.numeric(x$interaction)
  x
}
validation <- list()
for (name in c("pkn_expanded_chain", "pkn_expanded_glycan_conditional")) {
  path <- file.path(output_dir, paste0(name, ".csv"))
  pkn <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(identical(names(pkn), c("source", "interaction", "target")),
            !anyNA(pkn), all(pkn$interaction %in% c(-1, 1)),
            all(nzchar(pkn$source)), all(nzchar(pkn$target)),
            !any(pkn$source == pkn$target),
            !anyDuplicated(pkn[c("source", "target")]))
  clean <- cosmosR::meta_network_cleanup(pkn)
  stopifnot(identical(canonical(clean), canonical(pkn)))
  saveRDS(pkn, file.path(output_dir, paste0(name, ".rds")), compress = "xz")
  validation[[name]] <- list(
    edges = nrow(pkn), nodes = length(unique(c(pkn$source, pkn$target))),
    cosmos_cleaner_unchanged = TRUE, valid_signed_projection = TRUE)
  cat(name, ":", nrow(pkn), "edges; cleaner unchanged\n")
}
jsonlite::write_json(list(
  validation = validation,
  cosmosR_version = as.character(utils::packageVersion("cosmosR")),
  R_version = R.version.string, no_omics_scoring_performed = TRUE
), file.path(output_dir, "r_validation.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(output_dir, "R_sessionInfo.txt"))
