source("scripts/config/analysis_release.R")
# Historical two-marker normalization audit; current scientific QA is separate.
#!/usr/bin/env Rscript
# Read-only input audit; writes only this analysis's QA/provenance directory.
# Run from the repository root. No differential analysis or plotting occurs here.

options(stringsAsFactors = FALSE, width = 160)
root <- normalizePath(".", mustWork = TRUE)
local_library <- file.path(root, ".local", "R-library")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
output <- file.path(root, "miscellaneous/local/repository_cleanup/historical_input_provenance")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name) {
  utils::write.csv(x, file.path(output, name), row.names = FALSE, na = "NA")
}

packages <- c("limma", "edgeR", "statmod", "ggplot2", "svglite", "jsonlite", "digest", "SeuratObject", "Matrix")
package_audit <- do.call(rbind, lapply(packages, function(package) {
  available <- requireNamespace(package, quietly = TRUE)
  data.frame(package = package, available = available,
             version = if (available) as.character(utils::packageVersion(package)) else NA_character_,
             path = if (available) find.package(package) else NA_character_)
}))
write_csv(package_audit, "package_availability.csv")
print(package_audit[, c("package", "available", "version")])
stopifnot(all(package_audit$available[package_audit$package %in% c("SeuratObject", "Matrix")]))
source("scripts/support_functions.R")
source("scripts/core/wdsd_intake.R")
intake <- wp2_load_inputs(root,input_zonation=analysis_path("historical_wdsd_zonation_rds",root))
write_csv(intake$metadata, "mouse_pairing.csv")
write_csv(intake$layer_inventory, "RNA_layer_inventory.csv")
write_csv(intake$rna_profile_qc, "RNA_profile_QC.csv")
write_csv(intake$zonation_marker_diagnostic, "zonation_marker_diagnostics.csv")
write_csv(intake$marker_by_profile, "marker_expression.csv")
write_csv(intake$qc_summary, "input_QC_summary.csv")
writeLines(intake$limitations, file.path(output, "inherited_input_limitations.txt"))

relationships <- list()
column_relationships <- list()
objects <- list()
object_paths <- intake$input_paths[c("hepatocytes", "zonation")]
history_text <- character()
for (kind in names(object_paths)) {
  object <- readRDS(object_paths[[kind]])
  objects[[kind]] <- object
  counts <- as.matrix(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"))
  data <- as.matrix(SeuratObject::LayerData(object, assay = "RNA", layer = "data"))
  scaled <- as.matrix(SeuratObject::LayerData(object, assay = "RNA", layer = "scale.data"))
  stopifnot(identical(dimnames(counts), dimnames(data)), identical(dimnames(data), dimnames(scaled)),
            all(counts >= 0), all(counts == round(counts)), all(is.finite(counts)))
  library_size <- colSums(counts)
  relative <- sweep(counts, 2L, library_size, "/")
  candidates <- list(
    log1p_counts = log1p(counts),
    log1p_library_normalized_1e4 = log1p(relative * 1e4),
    log2p1_library_normalized_1e4 = log2(1 + relative * 1e4),
    log10p1_library_normalized_1e4 = log10(1 + relative * 1e4),
    log1p_library_normalized_1e6 = log1p(relative * 1e6))
  relationships[[kind]] <- do.call(rbind, lapply(names(candidates), function(name) {
    difference <- candidates[[name]] - data
    data.frame(object_kind = kind, proposed_data_transform = name,
               max_abs_error = max(abs(difference)), rms_error = sqrt(mean(difference^2)),
               pearson_all_entries = stats::cor(as.vector(candidates[[name]]), as.vector(data)))
  }))
  inferred_scales <- vapply(seq_len(ncol(counts)), function(j) {
    positive <- counts[, j] > 0
    stats::median(expm1(data[positive, j]) / counts[positive, j]) * library_size[j]
  }, numeric(1L))
  inferred_prediction <- log1p(sweep(relative, 2L, inferred_scales, "*"))
  column_relationships[[kind]] <- do.call(rbind, lapply(seq_len(ncol(counts)), function(j) {
    positive <- counts[, j] > 0
    k <- expm1(data[positive, j]) / counts[positive, j]
    data.frame(object_kind = kind, profile_id = colnames(counts)[j],
               supplied_count_total = library_size[j], count_detected_genes = sum(positive),
               positive_data_at_zero_count = sum(data[, j] > 0 & !positive),
               zero_data_at_positive_count = sum(data[, j] == 0 & positive),
               inferred_log1p_normalization_scale = inferred_scales[j],
               inferred_count_multiplier_min = min(k), inferred_count_multiplier_median = stats::median(k),
               inferred_count_multiplier_max = max(k),
               inferred_transform_max_abs_error = max(abs(inferred_prediction[, j] - data[, j])),
               log1p_norm1e4_max_abs_error = max(abs(candidates$log1p_library_normalized_1e4[, j] - data[, j])),
               data_count_spearman = stats::cor(data[, j], counts[, j], method = "spearman"))
  }))
  relationships[[kind]] <- rbind(relationships[[kind]],
    data.frame(object_kind = kind, proposed_data_transform = "log1p_library_normalized_fitted_per_profile_scale",
               max_abs_error = max(abs(inferred_prediction - data)),
               rms_error = sqrt(mean((inferred_prediction - data)^2)),
               pearson_all_entries = stats::cor(as.vector(inferred_prediction), as.vector(data))),
    data.frame(object_kind = kind, proposed_data_transform = "row_zscore_stored_data_to_scale.data",
               max_abs_error = max(abs(wp2_row_scale(data) - scaled)),
               rms_error = sqrt(mean((wp2_row_scale(data) - scaled)^2)),
               pearson_all_entries = stats::cor(as.vector(wp2_row_scale(data)), as.vector(scaled))))
  history_text <- c(history_text, paste0("OBJECT: ", kind),
    paste0("Path: ", object_paths[[kind]]),
    paste0("Object version: ", as.character(methods::slot(object, "version"))),
    "Metadata columns:", names(object[[]]),
    "Recorded commands and parameters:",
    capture.output(str(methods::slot(object, "commands"), max.level = 4L)),
    "Object misc:", capture.output(str(methods::slot(object, "misc"), max.level = 3L)),
    "Assay misc:", capture.output(str(methods::slot(object[["RNA"]], "misc"), max.level = 3L)), "")
  cat(kind, ": ", nrow(counts), " genes x ", ncol(counts), " profiles; count library range ",
      min(library_size), "--", max(library_size), "\n", sep = "")
  print(relationships[[kind]][, c("proposed_data_transform", "max_abs_error", "rms_error")])
}
relationship_table <- do.call(rbind, relationships)
column_table <- do.call(rbind, column_relationships)
write_csv(relationship_table, "RNA_transform_relationships.csv")
write_csv(column_table, "RNA_profile_transform_relationships.csv")
writeLines(history_text, file.path(output, "RNA_object_metadata_and_history.txt"))

# A diagnostic of aggregate count nesting, not proof of the original nuclear gate.
hepa_counts <- as.matrix(SeuratObject::LayerData(objects$hepatocytes, assay = "RNA", layer = "counts"))
hepa_metadata <- objects$hepatocytes[[]]
hepa_counts <- hepa_counts[, match(intake$metadata$mouse_id, as.character(hepa_metadata$cf453_id)), drop = FALSE]
colnames(hepa_counts) <- intake$metadata$mouse_id
zone_sum <- intake$rna_counts$Periportal + intake$rna_counts$Pericentral
stopifnot(identical(dimnames(hepa_counts), dimnames(zone_sum)))
count_nesting <- data.frame(mouse_id = intake$metadata$mouse_id,
  whole_hepatocyte_count_total = colSums(hepa_counts),
  two_zone_count_total = colSums(zone_sum),
  two_zone_to_hepatocyte_library_ratio = colSums(zone_sum) / colSums(hepa_counts),
  genes_two_zones_exceed_hepatocyte = colSums(zone_sum > hepa_counts),
  maximum_two_zone_excess_count = apply(pmax(zone_sum - hepa_counts, 0), 2L, max))
write_csv(count_nesting, "RNA_zone_count_nesting.csv")

documents <- c("data/WP1_WP2_INPUTS.md", "data/wp1_wp2_input_migration_manifest.csv",
  "docs/research_plans/2026-10-05_Western_diet_and_Cyp2c38_analysis_plan.md",
  "docs/workflow_plans/05_zonation_analysis.md", "docs/workflow_plans/12_wdsd_network_comparison.md",
  "scripts/core/wdsd_intake.R", "scripts/workflows/00_snrnaseq_preprocessing.R")
evidence <- do.call(rbind, lapply(documents, function(path) {
  lines <- readLines(path, warn = FALSE)
  hit <- which(grepl("cf453|wdsd|AggregateExpression|normaliz|imput|pseudobulk|ScaleData", lines, ignore.case = TRUE))
  data.frame(path = path, line = hit, text = lines[hit])
}))
write_csv(evidence, "source_documentation_evidence.csv")
archive_paths <- c("../original_datasets/wdsd_snrnaseq", "../original_datasets/wdsd_lipidomics")
archive_inventory <- do.call(rbind, lapply(archive_paths, function(path) {
  files <- list.files(path, full.names = TRUE, recursive = TRUE)
  data.frame(path = files, bytes = file.info(files)$size)
}))
write_csv(archive_inventory, "source_archive_file_inventory.csv")
if (requireNamespace("digest", quietly = TRUE)) {
  paths <- c(unname(intake$input_paths), file.path(root, "tests/qa_cf453_provenance.R"))
  write_csv(data.frame(path = paths, bytes = file.info(paths)$size,
    sha256 = vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1L))),
    "input_and_audit_script_SHA256.csv")
}

exact_ln <- relationship_table$max_abs_error[
  relationship_table$proposed_data_transform == "log1p_library_normalized_1e4"]
exact_scale <- relationship_table$max_abs_error[
  relationship_table$proposed_data_transform == "row_zscore_stored_data_to_scale.data"]
stopifnot(length(exact_ln) == 2L, all(exact_ln < 1e-10),
          length(exact_scale) == 2L, all(exact_scale < 1e-10),
          all(intake$zonation_marker_diagnostic$expected_marker_order),
          nrow(intake$metadata) == 28L, all(table(intake$metadata$group) == 7L))
summary_lines <- c(
  "CF453 supplied RNA provenance and layer audit",
  paste0("Run date: ", format(Sys.Date(), "%Y-%m-%d")),
  "Scope: input provenance and aggregate QC only; no differential expression, networks or figures.",
  "28 paired biological mice, seven per sex/diet group; 28 whole-hepatocyte and 56 paired-zone profiles.",
  "All supplied counts are finite nonnegative integers; counts/data/scale.data have identical named gene/profile rosters.",
  paste0("Maximum errors against log1p(count / profile_count_total * 10000), hepatocytes/zonation: ",
         paste(format(exact_ln, digits = 12), collapse = ", ")),
  paste0("Maximum errors for per-gene z-scoring stored data, hepatocytes/zonation: ",
         paste(format(exact_scale, digits = 12), collapse = ", ")),
  paste0("Expected aggregate marker order: ", sum(intake$zonation_marker_diagnostic$expected_marker_order), "/56 profiles."),
  paste0("Mouse/gene entries where two-zone count sum exceeds whole-hepatocyte count: ",
         sum(count_nesting$genes_two_zones_exceed_hepatocyte), "; two-zone/whole-hepatocyte library ratio range: ",
         paste(format(range(count_nesting$two_zone_to_hepatocyte_library_ratio), digits = 8), collapse = "--")),
  "Integer counts and numeric transformation replay establish layer properties, not how counts were originally aggregated.",
  "Nonnegative integer counts, library-normalization replay and nested zone counts are consistent with aggregate molecular counts, but do not certify that unnormalized counts were summed from the intended nuclear population.",
  "Search scope: repository scripts/docs, byte-import manifest and original_datasets/wdsd_snrnaseq plus wdsd_lipidomics archive inventories.",
  "No CF453 upstream count-summing/aggregation script or source methods document was identified in those locations.",
  "The two source RNA archive files are processed Seurat pseudobulks; command histories and misc content are captured separately.",
  "Per-profile nucleus counts, upstream nuclear QC, original gates, batch, litter and randomization metadata remain unestablished.",
  "Do not treat count totals as numbers of nuclei.",
  "DE scale recommendation: supplied unscaled data is a numerically verified log1p-library-normalized expression scale, suitable for a transparently described primary limma analysis with paired-zone mouse blocking and an appropriate abundance filter.",
  "Report coefficients as differences in natural-log normalized expression. They are not exact log2 fold changes because of the pseudocount and natural-log base.",
  "For direct Figure 5 analytical comparability, supplied scale.data is also available: its coefficients are standardized expression differences, not log2FC. This is a different cross-gene modeling scale; do not silently switch between the two.",
  "Count-based TMM/voom analysis requires an explicit count-origin assumption until source confirmation of summed unnormalized molecular counts; label it as a sensitivity under that assumption.",
  "Marker diagnostics are aggregate checks and cannot validate original nucleus-level gating.")
writeLines(summary_lines, file.path(output, "audit_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
cat("Saved provenance audit to ", output, "\n", sep = "")
cat(paste(summary_lines, collapse = "\n"), "\n")
