#!/usr/bin/env Rscript
# Export only the seven declared sample-level pseudobulk inputs; no nuclei.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript from this repository")
setwd(normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE))
source("scripts/support_functions.R")
source("scripts/core/pseudobulk_csv.R")
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/maintenance/export_pseudobulk_csv.R [--validate-only] [--dataset=ROLE]\n",
      "Exports seven declared pseudobulks to data/RNA/pseudobulk_csv; originals and approved results are preserved.\n",
      "Existing manifests are verified/reused. Local source RDS is required for exact round-trip validation.\n")
  quit(save = "no")
}
if (any(!grepl("^--validate-only$|^--dataset=.+$", args)) || anyDuplicated(sub("=.*$", "", args)))
  stop("Unknown or repeated option")
sources <- c(
  current_fcg_zonation = "data/RNA/zonation/current/fcg_zonation_pseudobulks_score.rds",
  current_gdx_zonation = "data/RNA/zonation/current/gdx_zonation_pseudobulks_score.rds",
  current_wdsd_zonation = "data/RNA/zonation/current/wdsd_zonation_pseudobulks_score.rds",
  wdsd_hepatocytes = "data/RNA/wdsd/wdsd_hep_pseudobulks.rds",
  historical_wdsd_zonation = "miscellaneous/older_inputs/two_marker_zonation/RNA/wdsd/wdsd_zonation_pseudobulks.rds",
  wd_cell_types = "data/RNA/wd/wd_pseudobulks.rds",
  wd_zonation = "data/RNA/wd/wd_zonation_pseudobulks.rds")
selected <- args[startsWith(args, "--dataset=")]
if (length(selected)) {
  name <- sub("^--dataset=", "", selected)
  if (!name %in% names(sources)) stop("Unknown pseudobulk role: ", name)
  sources <- sources[name]
}
results <- list()
for (role in names(sources)) {
  source <- sources[[role]]
  stopifnot(file.exists(source))
  output <- file.path("data/RNA/pseudobulk_csv", role)
  manifest <- file.path(output, "manifest.json")
  before <- pb_sha256(source)
  if (!file.exists(manifest)) {
    if ("--validate-only" %in% args) stop("Transport missing: ", manifest)
    pb_csv_write(source, output)
  }
  original <- pb_components(source)
  imported <- pb_read(manifest)
  stopifnot(identical(original$features, imported$features),
    identical(original$cells, imported$cells),
    identical(original$meta.data, imported$meta.data),
    identical(original$identities, imported$identities),
    identical(original$identity_levels, imported$identity_levels),
    identical(original$identity_class, imported$identity_class),
    identical(original$recorded_commands, imported$recorded_commands),
    identical(original$object_version, imported$object_version),
    identical(original$layer_classes, imported$layer_classes),
    identical(original$layers, imported$layers),
    identical(before, pb_sha256(source)), identical(before, imported$source$sha256))
  files <- pb_source_files(manifest)
  csvs <- files[grepl("[.]csv$", files)]
  stopifnot(all(file.info(csvs)$size < 50000000))
  record <- list(role = role, source_path = source, source_sha256 = before,
    manifest_path = manifest, manifest_sha256 = pb_sha256(manifest),
    n_features = length(imported$features), n_profiles = length(imported$profile_ids),
    layers = names(imported$layers), counts_nonnegative_integer = TRUE,
    layer_values_bitwise_identical = TRUE, metadata_identical = TRUE,
    identities_and_levels_identical = TRUE, source_preserved = TRUE,
    maximum_absolute_error = 0, csv_files = length(csvs),
    maximum_csv_bytes = max(file.info(csvs)$size),
    payloads = imported$manifest$files)
  results[[role]] <- record
  cat(role, ":", record$n_features, "genes x", record$n_profiles,
      "profiles; exact three-layer/metadata round trip; largest CSV", record$maximum_csv_bytes, "bytes\n")
  rm(original, imported); invisible(gc(FALSE))
}
if (!length(selected) && !"--validate-only" %in% args) {
  dir.create("docs/provenance", showWarnings = FALSE)
  jsonlite::write_json(list(format = "pseudobulk_csv_export_audit", format_version = 1L,
    maximum_csv_bytes_exclusive = 50000000L, individual_nuclei_exported = FALSE,
    scientific_models_run = FALSE, source_rds_preserved = TRUE, datasets = results),
    "docs/provenance/pseudobulk_csv_roundtrip.json", pretty = TRUE,
    auto_unbox = TRUE, digits = NA, null = "null")
}
cat("Pseudobulk transport verification passed; no scientific models were run.\n")
