# Verify the frozen imported bytes without requiring the sibling source archive.
# Run from the repository root; the original archive remains a provenance record.
source("scripts/support_functions.R")
source("scripts/config/analysis_release.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
if (!requireNamespace("digest", quietly = TRUE)) stop("Missing package: digest")
manifest <- utils::read.csv("data/wp1_wp2_input_migration_manifest.csv",
                            check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 6L,
          !anyDuplicated(manifest$destination_relative_to_repository),
          all(toupper(as.character(manifest$source_archive_preserved)) == "TRUE"),
          all(manifest$source_sha256 == manifest$destination_sha256),
          all(file.exists(analysis_resolve_recorded_path(manifest$destination_relative_to_repository))))
actual_paths <- analysis_resolve_recorded_path(manifest$destination_relative_to_repository,mustWork=TRUE)
actual_sha <- vapply(actual_paths,
                     digest::digest, character(1L), file = TRUE, algo = "sha256")
stopifnot(identical(unname(actual_sha), manifest$destination_sha256),
          identical(as.numeric(file.info(actual_paths)$size),
                    as.numeric(manifest$bytes)))
stopifnot(sum(manifest$bytes) == 30965755,
          all(grepl("^(data/(RNA|metabolomics)/|docs/reference_figures/)",
                    manifest$destination_relative_to_repository)))
cat("Six WP1/WP2 repository input copies match frozen SHA256 and byte counts.\n")
