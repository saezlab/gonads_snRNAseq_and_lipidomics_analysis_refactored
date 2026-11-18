# Historical WP2 two-marker fixture; current science checks are qa_current_zonation_release.R.
# Meaningful WDSD intake checks: mouse pairing, stable identity, exclusion/QC
# policy and stored-layer fidelity. Run from the repository root.
source("scripts/support_functions.R")
source("scripts/core/wdsd_intake.R")
source("scripts/core/pseudobulk_csv.R")
source("scripts/core/cyp2c38_reference.R")

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected error matching: ", pattern)
  }
}
message("Testing synthetic alignment and integrity failures...")
mouse_ids <- paste0("cf453-", seq_len(28L))
lipid_columns <- paste0(toupper(mouse_ids), "_", rep(c("F-SD", "F-WD", "M-SD", "M-WD"), each = 7L))
meta <- wp2_lipid_sample_metadata(rev(lipid_columns))
stopifnot(identical(meta$mouse_id, mouse_ids),
          identical(meta$lipid_sample, lipid_columns),
          identical(as.integer(table(meta$group)), rep(7L, 4L)))
expect_error(wp2_lipid_sample_metadata(c(lipid_columns[-1L], lipid_columns[2L])), "uniquely")
wrong_labels <- lipid_columns
wrong_labels[1L] <- "CF453-1_M-SD"
expect_error(wp2_lipid_sample_metadata(wrong_labels), "sex/diet")

profile <- data.frame(profile_id = c("a_PP", "b_PC", "a_PC", "b_PP"),
                      mouse_id = c("a", "b", "a", "b"),
                      zonation_type = c("Periportal", "Pericentral", "Pericentral", "Periportal"),
                      stringsAsFactors = FALSE)
matrix <- matrix(c(11, 101, 22, 202, 33, 303, 44, 404), nrow = 2L,
                 dimnames = list(c("GeneA", "GeneB"), profile$profile_id))
scrambled <- matrix[, c(4L, 3L, 1L, 2L)]
pp <- wp2_align_rna_zone(scrambled, profile[c(3L, 1L, 4L, 2L), ], c("b", "a"), "Periportal")
pc <- wp2_align_rna_zone(scrambled, profile, c("b", "a"), "Pericentral")
stopifnot(identical(unname(pp), matrix[, c(4L, 1L), drop = FALSE] |> unname()),
          identical(unname(pc), matrix[, c(2L, 3L), drop = FALSE] |> unname()),
          identical(colnames(pp), c("b", "a")), identical(colnames(pc), c("b", "a")))
duplicate_profile <- profile
duplicate_profile$mouse_id[4L] <- "a"
expect_error(wp2_align_rna_zone(scrambled, duplicate_profile, c("b", "a"), "Periportal"), "exactly one")
expect_error(wp2_align_rna_zone(scrambled[, -1L], profile, c("a", "b"), "Periportal"), "unmatched")

synthetic <- rbind(varying = c(1, 2, 3, 5), constant = rep(7, 4))
scaled <- wp2_row_scale(synthetic)
stopifnot(max(abs(scaled[1L, ] - as.numeric(scale(synthetic[1L, ])))) < 1e-14,
          all(scaled[2L, ] == 0), all(is.finite(scaled)),
          identical(dimnames(scaled), dimnames(synthetic)))
bad <- synthetic
bad[1L, 1L] <- NA_real_
expect_error(wp2_row_scale(bad), "finite numeric")

message("Loading and checking supplied WDSD inputs...")
input <- wp2_load_inputs(input_zonation=analysis_path("historical_wdsd_zonation_transport"))
expected_local_inputs <- c(
  lipidomics = "data/metabolomics/wdsd/lipidomics_cf453.csv",
  hepatocytes = analysis_path("wdsd_hepatocytes_transport"),
  zonation = analysis_path("historical_wdsd_zonation_transport")
)
stopifnot(identical(unname(input$input_paths[names(expected_local_inputs)]),
                    unname(vapply(expected_local_inputs, normalizePath, character(1L), mustWork = TRUE))))
stopifnot(identical(input$metadata$mouse_id, mouse_ids),
          identical(dim(input$lipid), c(1177L, 28L)),
          identical(colnames(input$lipid), mouse_ids),
          nrow(input$lipid_ledger) == 1240L,
          sum(input$lipid_ledger$keep_primary) == 1177L,
          sum(input$lipid_ledger$keep_qc30) == 1129L,
          sum(input$lipid_ledger$exclusion_reason == "internal_standard") == 14L,
          sum(input$lipid_ledger$exclusion_reason == "repository_excluded_class") == 49L,
          !anyDuplicated(input$lipid_ledger$feature_id),
          !anyDuplicated(input$lipid_ledger$canonical_lipid),
          sum(duplicated(input$lipid_ledger[["Alignment ID"]])) == 67L,
          all(input$lipid_ledger$alignment_adduct_multiplicity == 1L),
          max(abs(input$lipid - log2(input$lipid_supplied))) == 0,
          all(input$lipid_ledger$keep_qc30 <= input$lipid_ledger$keep_primary),
          all(input$lipid_ledger$canonical_lipid == wp1_canonical_lipid(input$lipid_ledger[["Metabolite name"]])))
for (zone in c("Periportal", "Pericentral")) {
  stopifnot(identical(dim(input$rna[[zone]]), c(20123L, 28L)),
            identical(colnames(input$rna[[zone]]), mouse_ids),
            identical(colnames(input$rna_data[[zone]]), mouse_ids),
            identical(colnames(input$rna_counts[[zone]]), mouse_ids),
            all(input$rna_counts[[zone]] >= 0),
            all(input$rna_counts[[zone]] == round(input$rna_counts[[zone]])),
            max(abs(input$rna[[zone]] - input$rna_data_scaled[[zone]])) < 1e-10)
}
stopifnot(all(input$zonation_marker_diagnostic$expected_marker_order),
          all(input$zonation_marker_diagnostic$aggregate_data_threshold_rule),
          nrow(input$rna_profile_qc) == 84L,
          nrow(input$layer_inventory) == 6L,
          all(input$layer_inventory$nonfinite_values == 0L),
          identical(input$recorded_commands$hepatocytes, "ScaleData.RNA"),
          identical(input$recorded_commands$zonation, "ScaleData.RNA"))

# Independently access the original Assay5 layers and RNA metadata. Tests must
# catch accidental lexicographic or exported-lipid order pairing mistakes.
original <- pb_components(input$input_paths[["zonation"]])
original_metadata <- original$meta.data
original_scale <- original$layers$scale.data
for (zone in c("Periportal", "Pericentral")) {
  # Build the roster from filtered metadata before applying mouse ordering.
  filtered <- original_metadata[original_metadata$zonation_type == zone, , drop = FALSE]
  names <- rownames(filtered)[match(mouse_ids, filtered$cf453_id)]
  values <- as.matrix(original_scale[, names, drop = FALSE])
  stopifnot(identical(unname(input$rna[[zone]]), unname(values)))
}
message("All WDSD intake checks passed.")
