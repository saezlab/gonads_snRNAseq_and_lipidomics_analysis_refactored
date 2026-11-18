# CF459 intake: locale parsing, real mouse pairing, stored-layer fidelity,
# unavailable QC, immutable imported bytes and explicit metadata failures.
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/wd_intake.R")
source("scripts/core/pseudobulk_csv.R")
source("scripts/config/analysis_release.R")
source("scripts/core/wdsd_intake.R")
source("scripts/core/cyp2c38_reference.R")

expect_error <- function(expression, pattern) {
  error <- tryCatch({force(expression); NULL}, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected error matching: ", pattern)
  }
}

manifest <- utils::read.csv("data/wp3_input_migration_manifest.csv", stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 3L, !anyDuplicated(manifest$destination_relative_to_repository),
  all(manifest$source_sha256 == manifest$destination_sha256),
  all(toupper(manifest$source_archive_preserved) == "TRUE"), sum(manifest$bytes) == 13118772)
if (!requireNamespace("digest", quietly = TRUE)) stop("Missing digest package.")
source_verification<-analysis_verify_inputs(manifest$destination_relative_to_repository,manifest$destination_sha256)
stopifnot(all(source_verification$verified))

mouse_ids <- paste0("cf459-", c(1L:4L, 13L:16L))
meta <- wd_lipid_sample_metadata(rev(mouse_ids))
stopifnot(identical(meta$mouse_id, mouse_ids), identical(meta$lipid_sample, mouse_ids),
          all(meta$diet == "WD"), identical(meta$sex, rep(c("female", "male"), each = 4L)))
expect_error(wd_lipid_sample_metadata(c(mouse_ids[-1L], mouse_ids[2L])), "eight unique")
expect_error(wd_lipid_sample_metadata(sub("cf459-16", "cf459-17", mouse_ids)), "eight unique")

input <- wd_load_inputs()
stopifnot(identical(input$metadata$mouse_id, mouse_ids),
  identical(dim(input$lipid), c(1395L, 8L)), identical(colnames(input$lipid), mouse_ids),
  nrow(input$lipid_ledger) == 1463L, sum(input$lipid_ledger$keep_primary) == 1395L,
  sum(input$lipid_ledger$exclusion_reason == "repository_excluded_class") == 68L,
  !any(input$lipid_ledger$internal_standard), !any(input$lipid_ledger$isotope_labelled),
  sum(duplicated(input$lipid_ledger[["Alignment ID"]])) == 244L,
  sum(input$lipid_ledger$pipe_delimited_annotation) == 1103L,
  !anyDuplicated(input$lipid_ledger$feature_id), !anyDuplicated(input$lipid_ledger$canonical_lipid),
  all(input$lipid_ledger$canonical_lipid == wp1_canonical_lipid(input$lipid_ledger[["Metabolite name"]])),
  all(is.na(input$lipid_ledger$qc_cv_percent)), all(is.na(input$lipid_ledger$keep_qc30)),
  !any(input$lipid_ledger$analytical_qc_available), !any(input$lipid_ledger$adduct_available),
  !any(input$lipid_ledger$fill_percentage_available), nrow(input$analytical_qc) == 0L,
  max(abs(input$lipid - log2(input$lipid_supplied))) == 0,
  input$lipid_supplied["WD_0001", "cf459-1"] == 0.000111221,
  all(is.na(input$metadata$nuclei_contributing)))
expected_paths <- c(lipidomics = "data/metabolomics/wd/lipidomics_wd.csv",
  cell_types = analysis_path("wd_cell_types_transport"), zonation = analysis_path("wd_zonation_transport"))
stopifnot(identical(input$input_paths[names(expected_paths)], vapply(expected_paths, normalizePath, character(1L), mustWork = TRUE)),
  nrow(input$rna_profile_metadata) == 56L, nrow(input$rna_profile_qc) == 56L,
  all(is.na(input$rna_profile_qc$nuclei_contributing)), nrow(input$layer_inventory) == 6L,
  identical(as.integer(input$layer_inventory$n_profiles), c(40L, 40L, 40L, 16L, 16L, 16L)),
  all(input$layer_inventory$nonfinite_values == 0),
  identical(input$recorded_commands$cell_types, "ScaleData.RNA"),
  identical(input$recorded_commands$zonation, "ScaleData.RNA"),
  all(input$object_versions == "5.3.0"), all(input$zonation_marker_diagnostic$expected_marker_order),
  all(input$zonation_marker_diagnostic$aggregate_data_threshold_rule),
  identical(dim(input$rna_hepatocytes), c(20123L, 8L)),
  identical(colnames(input$rna_hepatocytes), mouse_ids),
  max(abs(input$rna_hepatocytes - input$rna_hepatocytes_data_scaled)) < 1e-10)
for (zone in c("Periportal", "Pericentral")) {
  stopifnot(identical(dim(input$rna[[zone]]), c(20123L, 8L)),
    identical(colnames(input$rna[[zone]]), mouse_ids),
    identical(colnames(input$rna_data[[zone]]), mouse_ids),
    identical(colnames(input$rna_counts[[zone]]), mouse_ids),
    all(input$rna_counts[[zone]] >= 0), all(input$rna_counts[[zone]] == round(input$rna_counts[[zone]])),
    max(abs(input$rna[[zone]] - input$rna_data_scaled[[zone]])) < 1e-10)
}

# Independent original-layer lookup and deliberately different matrix/metadata
# order must preserve exact animal pairing, including cf459-13 versus cf459-1.
for (kind in c("cell_types", "zonation")) {
  original <- pb_components(input$input_paths[[kind]])
  m <- original$meta.data
  validated <- wd_validate_rna_metadata(m, original$cells, meta, kind)
  bad <- m
  bad$group[1L] <- if (bad$group[1L] == "f-ctrl") "m-ctrl" else "f-ctrl"
  expect_error(wd_validate_rna_metadata(bad, original$cells, meta, kind), "sex assignments disagree")
  bad <- m
  bad$sample_name[1L] <- "cf459-99"
  expect_error(wd_validate_rna_metadata(bad, original$cells, meta, kind), "identities or sex")
  bad <- m
  if (kind == "zonation") bad$zonation_type[1L] <- bad$zonation_type[2L] else bad$cell_type[1L] <- bad$cell_type[2L]
  expect_error(wd_validate_rna_metadata(bad, original$cells, meta, kind), "Missing or repeated")
  expect_error(wd_validate_rna_metadata(m, rev(original$cells), meta, kind), "Inconsistent")
  original_scale <- original$layers$scale.data
  if (kind == "zonation") {
    set.seed(20261014)
    scrambled <- original_scale[, sample(ncol(original_scale)), drop = FALSE]
    shuffled <- validated[sample(nrow(validated)), , drop = FALSE]
    for (zone in c("Periportal", "Pericentral")) {
      filtered <- m[m$zonation_type == zone, , drop = FALSE]
      profiles <- rownames(filtered)[match(mouse_ids, filtered$sample_name)]
      expected <- as.matrix(original_scale[, profiles, drop = FALSE])
      aligned <- wp2_align_rna_zone(scrambled, shuffled, mouse_ids, zone)
      stopifnot(identical(unname(input$rna[[zone]]), unname(expected)),
                identical(unname(aligned), unname(expected)))
    }
    duplicated <- validated
    duplicated$mouse_id[2L] <- duplicated$mouse_id[4L]
    expect_error(wp2_align_rna_zone(original_scale, duplicated, mouse_ids, "Periportal"), "exactly one")
    expect_error(wp2_align_rna_zone(original_scale[, -1L], validated, mouse_ids, "Periportal"), "unmatched")
  } else {
    filtered <- m[m$cell_type == "Hepatocytes", , drop = FALSE]
    profiles <- rownames(filtered)[match(mouse_ids, filtered$sample_name)]
    stopifnot(identical(unname(input$rna_hepatocytes), unname(as.matrix(original_scale[, profiles, drop = FALSE]))))
  }
}
cat("CF459 input hashes, locale, pairing/permutations, metadata failures, exclusions and layer fidelity passed.\n")
