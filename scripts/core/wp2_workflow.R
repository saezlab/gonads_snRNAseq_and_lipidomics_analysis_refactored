source("scripts/config/analysis_release.R")
# Workflow 12 phase orchestration. Core statistics and feature universes are
# unchanged. Source scripts/support_functions.R first.

wp2_pooled_run <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(!grepl("^--(permutations|bootstraps|paired-bootstraps|seed|output)=|^--validate-only$", args))) {
    stop("Unknown pooled follow-up option: ", paste(args, collapse = " "))
  }
  option <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (length(hit) > 1L) stop("Duplicate option: ", name)
    if (length(hit)) substring(hit, nchar(name) + 4L) else default
  }
  P <- as.integer(option("permutations", "19999"))
  B <- as.integer(option("bootstraps", "4999"))
  paired_B <- as.integer(option("paired-bootstraps", "999"))
  seed <- as.integer(option("seed", "20261007"))
  stopifnot(all(is.finite(c(P, B, paired_B, seed))), min(P, B, paired_B) >= 19L)
  out <- option("output", file.path(analysis_path("historical_wf12"),"followup"))
  for (file in c("cyp2c38_reference.R", "wdsd_intake.R", "wdsd_pooled_reference.R",
                 "wdsd_lipid_overlap_audit.R", "wdsd_network_statistics.R",
                 "wdsd_pooled_statistics.R", "wdsd_pooled_setup.R")) source(file.path("scripts/core", file))
  data <- wp2_load_inputs(input_zonation=analysis_path("historical_wdsd_zonation_transport"))
  reference <- wp2_pooled_reference()
  audit <- wp2_load_lipid_overlap_audit(data = data)
  strict <- wp2_prepare_pooled_context(data, reference, audit$safe_map)
  conditional <- wp2_prepare_pooled_context(data, reference, audit$conditional_compatibility_map,
                                             "conditional_compatibility")
  coverage <- rbind(strict$coverage, conditional$coverage)
  print(coverage)
  if ("--validate-only" %in% args) {
    cat("Original pooled reference replay, molecular audit, 14-mouse SD pairing and mapping checks passed.\n")
    return(invisible(list(data = data, reference = reference, audit = audit, strict = strict, conditional = conditional)))
  }
  if (file.exists(file.path(out, "pooled_analysis.rds")) || file.exists(file.path(out, "README.md"))) {
    stop("Pooled follow-up already exists; choose a fresh --output destination.")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
  audit_dir <- file.path(out, "lipid_overlap_audit")
  if (dir.exists(audit_dir)) {
    staged <- readRDS(file.path(audit_dir, "lipid_overlap_audit.rds"))
    for (field in c("safe_map", "conditional_compatibility_map", "metrics", "candidate_pairs")) {
      stopifnot(isTRUE(all.equal(staged[[field]], audit[[field]], check.attributes = FALSE)))
    }
  } else wp2_write_lipid_overlap_audit(audit, audit_dir)
  saveRDS(audit, file.path(out, "audit.rds"))
  write_csv(data$metadata, "metadata.csv")
  write_csv(strict$context$metadata_SD, "SD_mouse_manifest.csv")
  write_csv(reference$metadata, "FCG_mouse_manifest.csv")
  write_csv(reference$replay_metrics, "pooled_reference_replay.csv")
  write_csv(coverage, "pooled_coverage.csv")
  write_csv(strict$ledger, "pooled_edge_eligibility_ledger.csv")
  write_csv(conditional$ledger, "conditional_edge_eligibility_ledger.csv")
  write_csv(data$lipid_ledger, "lipid_annotation_and_QC_ledger.csv")
  write_csv(data$qc_summary, "input_qc_summary.csv")
  saveRDS(strict$context, file.path(out, "paired_network_inputs.rds"))
  saveRDS(reference$provenance, file.path(out, "reference_provenance.rds"))
  writeLines(c(data$limitations, audit$limitations,
    "Pooled standard diet compares 14 wild-type WDSD mice with 43 original FCG mice; age, genetic background and cohort differences cannot be separated.",
    "Low supplied Fill values may indicate upstream missingness despite positive exported intensities; peak-level absence/imputation provenance unavailable."),
    file.path(out, "input_limitations.txt"))
  analysis <- wp2_pooled_analysis(strict$context, permutations = P, B = B, paired_B = paired_B, seed = seed)
  analysis$edges <- wp2_pooled_annotate_edges(analysis$edges, data$lipid_ledger)
  saveRDS(analysis, file.path(out, "pooled_analysis.rds"))
  tables <- c(edges = "matched_pooled_edges.csv", inference = "pooled_inference.csv",
    permutation_withinsex_draws = "withinsex_pairing_permutation_draws.csv",
    permutation_unrestricted_draws = "unrestricted_pairing_permutation_draws.csv",
    bootstrap_withinsex_draws = "withinsex_null_bootstrap_draws.csv",
    bootstrap_unrestricted_draws = "unrestricted_null_bootstrap_draws.csv",
    paired_stability_draws = "paired_SD_stability_draws.csv",
    paired_reference_stability_draws = "paired_FCG_SD_stability_draws.csv",
    resampling_indices = "resampling_mouse_indices.csv")
  for (nm in names(tables)) write_csv(analysis[[nm]], tables[[nm]])
  contexts <- list(strict = strict$context, conditional_compatibility = conditional$context)
  mapping_analysis <- wp2_pooled_analysis(conditional$context, permutations = P, B = B,
                                           paired_B = paired_B, seed = seed)
  saveRDS(mapping_analysis, file.path(out, "conditional_compatibility_analysis.rds"))
  mapping_inference <- rbind(cbind(mapping = "strict", analysis$inference),
    cbind(mapping = "conditional_compatibility", mapping_analysis$inference))
  write_csv(mapping_inference, "mapping_sensitivity_inference.csv")
  restrictions <- list(QC_CV_le30 = data$lipid_ledger$feature_id[data$lipid_ledger$keep_qc30],
    supplied_fill_ge50 = data$lipid_ledger$feature_id[data$lipid_ledger$keep_primary &
      is.finite(data$lipid_ledger[["Fill %"]]) & data$lipid_ledger[["Fill %"]] >= 0.5])
  sensitivity <- list()
  for (nm in names(restrictions)) {
    context <- strict$context
    context$edges <- context$edges[context$edges$feature_id %in% restrictions[[nm]], , drop = FALSE]
    contexts[[nm]] <- context
    z <- wp2_pooled_analysis(context, permutations = P, B = B, paired_B = paired_B, seed = seed)
    sensitivity[[nm]] <- cbind(sensitivity = nm, z$inference)
    saveRDS(z, file.path(out, paste0(nm, "_analysis.rds")))
  }
  sensitivity <- do.call(rbind, sensitivity)
  write_csv(sensitivity, "pooled_QC_sensitivity.csv")
  saveRDS(contexts, file.path(out, "pooled_contexts.rds"))
  # The stored RNA data layer is an affine equivalent, not independent evidence.
  alt <- data$rna_hepatocytes_data[rownames(strict$context$rna_SD),
                                  colnames(strict$context$rna_SD), drop = FALSE]
  alt_r <- wp2_cor_edges(alt, strict$context$lipid_SD, strict$context$edges)
  stopifnot(max(abs(alt_r - analysis$edges$r_SD)) < 1e-10)
  write_csv(data.frame(sensitivity = "stored_RNA_data_affine_replay",
    max_abs_r_difference = max(abs(alt_r - analysis$edges$r_SD))), "RNA_data_layer_sensitivity.csv")
  dir.create(file.path(out, "networks"), showWarnings = FALSE)
  for (universe in c("primary", "full")) {
    e <- analysis$edges[universe == "full" | analysis$edges$primary_edge, ]
    write_csv(e[, c("gene", "reference_lipid", "feature_id", "r_reference", "r_SD", "primary_edge")],
              paste0("networks/pooled_SD_", universe, "_SIF.csv"))
  }
  source("scripts/core/wdsd_pooled_presentation.R")
  wp2_pooled_make_plots(out, data, strict$context, analysis, audit, sensitivity)
  wp2_pooled_make_report(out)
  paths <- unique(c(unname(data$input_paths), reference$source_paths, unname(audit$input_paths),
    "scripts/support_functions.R", "scripts/core/cyp2c38_reference.R",
    paste0("scripts/core/", c("wdsd_intake.R", "wdsd_pooled_reference.R", "wdsd_lipid_overlap_audit.R",
      "wdsd_network_statistics.R", "wdsd_pooled_statistics.R", "wdsd_pooled_setup.R", "wdsd_pooled_presentation.R")),
    "scripts/workflows/12_wdsd_network_comparison.R", "scripts/core/wp2_workflow.R"))
  hash <- function(path) {
    value <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
    if (length(value) != 1L || !grepl("^[a-f0-9]{64} ", value)) stop("SHA256 failed for ", path)
    substr(value, 1L, 64L)
  }
  write_csv(data.frame(path = normalizePath(paths), SHA256 = vapply(paths, hash, character(1L))), "input_and_code_SHA256.csv")
  jsonlite::write_json(list(date_UTC = format(Sys.time(), tz = "UTC"), args = args,
    settings = analysis$metadata, coverage = coverage,
    interpretation = "strict species annotation primary; conditional alias compatibility sensitivity; pooled replication and beyond-sex checks distinct"),
    file.path(out, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("Completed lipid audit and pooled Figure 3 follow-up: ", normalizePath(out), "\n", sep = "")
  invisible(list(analysis = analysis, sensitivities = sensitivity, contexts = contexts))
}

wp2_run <- function(args = commandArgs(trailingOnly = TRUE), allow_pooled = FALSE) {
  if (any(!grepl("^--(bootstraps|paired-bootstraps|seed|output)=|^--validate-only$", args))) {
    stop("Unknown WP2 option: ", paste(args, collapse = " "))
  }
  option <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (length(hit) > 1L) stop("Duplicate option: ", name)
    if (length(hit)) substring(hit, nchar(name) + 4L) else default
  }
  B <- as.integer(option("bootstraps", "4999"))
  paired_B <- as.integer(option("paired-bootstraps", "999"))
  seed <- as.integer(option("seed", "20261007"))
  stopifnot(!is.na(B), !is.na(paired_B), !is.na(seed), B >= 19L, paired_B >= 19L)
  out <- option("output", analysis_path("historical_wf12"))
  source("scripts/core/cyp2c38_reference.R")
  source("scripts/core/wdsd_intake.R")
  source("scripts/core/wdsd_reference.R")
  source("scripts/core/wdsd_network_statistics.R")
  source("scripts/core/wdsd_comparison_setup.R")
  data <- wp2_load_inputs(input_zonation=analysis_path("historical_wdsd_zonation_transport"))
  reference <- wp2_load_reference()
  setup <- wp2_prepare_contexts(data, reference)
  contexts <- setup$contexts
  print(setup$coverage)
  if ("--validate-only" %in% args) {
    cat("WP2 intake, archived reference replay, strict molecular matching and seven-mouse pairing checks passed.\n")
    return(invisible(list(data = data, reference = reference, setup = setup)))
  }
  existing <- if (dir.exists(out)) list.files(out, all.files = TRUE, no.. = TRUE) else character()
  permitted <- if (allow_pooled) "followup" else character()
  if (length(setdiff(existing, permitted))) {
    stop("Output directory contains existing zonation results; choose a fresh --output destination.")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
  write_csv(data$metadata, "metadata.csv")
  write_csv(data$lipid_ledger, "lipid_annotation_and_QC_ledger.csv")
  write_csv(data$qc_summary, "input_qc_summary.csv")
  write_csv(data$rna_profile_metadata, "RNA_profile_metadata.csv")
  write_csv(data$rna_profile_qc, "RNA_profile_QC.csv")
  write_csv(data$layer_inventory, "RNA_layer_inventory.csv")
  write_csv(data$marker_by_profile, "zonation_markers_by_profile.csv")
  write_csv(data$marker_summary, "zonation_marker_summary.csv")
  write_csv(data$zonation_marker_diagnostic, "zonation_marker_diagnostics.csv")
  write_csv(attr(reference, "summary"), "reference_replay_summary.csv")
  write_csv(do.call(rbind, lapply(names(reference), function(name) {
    cbind(context = name, reference[[name]]$metadata)
  })), "reference_mouse_manifest.csv")
  write_csv(setup$coverage, "network_coverage.csv")
  write_csv(setup$lipid_crosswalk, "lipid_identity_crosswalk.csv")
  write_csv(setup$edge_eligibility_ledger, "edge_eligibility_ledger.csv")
  utils::write.table(cbind(feature_id = rownames(data$lipid), data$lipid),
                     file.path(out, "log2_lipid_analysis_matrix.csv"), sep = ",",
                     row.names = FALSE, col.names = TRUE, quote = TRUE)
  writeLines(c(data$limitations,
    "Some supplied Fill % values are below 0.5 despite positive exported intensities; upstream missingness and imputation provenance are unavailable.",
    "Cross-cohort comparisons involve older WDSD wild-type animals and younger original FCG animals; age, background and genotype differences cannot be separated here."),
    file.path(out, "input_limitations.txt"))
  saveRDS(contexts, file.path(out, "paired_network_inputs.rds"))
  saveRDS(list(preprocessing = data$preprocessing, commands = data$recorded_commands,
               object_versions = data$object_versions, reference = attr(reference, "provenance")),
          file.path(out, "input_provenance.rds"))

  focus <- c("Mtmr11", "Sult3a1", "Abcb1a", "Cyp2a4", "Mindy3", "Esco2", "Kitl")
  analysis <- wp2_network_analysis(contexts, B = B, paired_B = paired_B,
                                    seed = seed, focus_genes = focus)
  analysis$edges$canonical_lipid <- data$lipid_ledger$canonical_lipid[
    match(analysis$edges$feature_id, data$lipid_ledger$feature_id)]
  analysis$edges$lipid_class <- data$lipid_ledger$lipid_class[
    match(analysis$edges$feature_id, data$lipid_ledger$feature_id)]
  table_names <- c(edges = "network_edges.csv", inference = "baseline_inference.csv",
    comparison_summary = "diet_comparison_summary.csv", permutation_draws = "pairing_permutation_draws.csv",
    bootstrap_null_draws = "independent_marginal_null_bootstrap_draws.csv",
    paired_stability_draws = "paired_mouse_stability_draws.csv",
    paired_reference_stability_draws = "paired_reference_stability_draws.csv",
    resampling_indices = "resampling_mouse_indices.csv",
    pairing_permutation_indices = "exact_pairing_permutation_indices.csv",
    loco_omission_ledger = "LOCO_mouse_omission_ledger.csv",
    loco_aggregate_summary = "loco_aggregate_summary.csv")
  for (name in names(table_names)) write_csv(analysis[[name]], table_names[[name]])
  saveRDS(analysis$loco_matrices, file.path(out, "LOCO_edge_draws.rds"))
  saveRDS(analysis, file.path(out, "network_analysis.rds"))
  write_csv(analysis$edges[analysis$edges$gene == "Mtmr11", ], "Mtmr11_edges.csv")

  # QC-only restrictions use annotation information, never observed correlation.
  # They retain the corresponding subset of the same originally frozen pairs.
  restrictions <- list(QC_CV_le30 = data$lipid_ledger$feature_id[data$lipid_ledger$keep_qc30],
    supplied_fill_ge50 = data$lipid_ledger$feature_id[data$lipid_ledger$keep_primary &
      is.finite(data$lipid_ledger[["Fill %"]]) & data$lipid_ledger[["Fill %"]] >= 0.5])
  sensitivities <- list()
  for (name in names(restrictions)) {
    restricted <- lapply(contexts, function(x) {
      x$edges <- x$edges[x$edges$feature_id %in% restrictions[[name]], , drop = FALSE]
      x
    })
    z <- wp2_network_analysis(restricted, B = B, paired_B = paired_B, seed = seed,
                              focus_genes = "Mtmr11")
    sensitivities[[name]] <- cbind(sensitivity = name, z$inference)
    write_csv(z$comparison_summary, paste0(name, "_diet_comparison.csv"))
    saveRDS(z, file.path(out, paste0(name, "_analysis.rds")))
  }
  # Stored data is an affine version of scale.data: demonstrate rather than
  # portray this as an independent normalization validation.
  data_sensitivity <- do.call(rbind, lapply(names(contexts), function(name) {
    x <- contexts[[name]]; alt <- data$rna_data[[x$zone]]
    r_sd <- wp2_cor_edges(alt[, colnames(x$rna_SD), drop = FALSE], x$lipid_SD, x$edges)
    r_wd <- wp2_cor_edges(alt[, colnames(x$rna_WD), drop = FALSE], x$lipid_WD, x$edges)
    observed <- analysis$edges[analysis$edges$context == name, ]
    stopifnot(identical(x$edges$edge_id, observed$edge_id),
              max(abs(r_sd - observed$r_SD)) < 1e-10, max(abs(r_wd - observed$r_WD)) < 1e-10)
    data.frame(context = name, sensitivity = "stored_RNA_data_affine_replay",
      max_abs_SD_r_difference = max(abs(r_sd - observed$r_SD)),
      max_abs_WD_r_difference = max(abs(r_wd - observed$r_WD)))
  }))
  write_csv(data_sensitivity, "RNA_data_layer_sensitivity.csv")
  sensitivity_table <- do.call(rbind, sensitivities)
  write_csv(sensitivity_table, "sensitivities.csv")

  change_summary <- function(column) {
    a <- analysis$edges
    groups <- split(seq_len(nrow(a)), interaction(a$context, a[[column]], drop = TRUE))
    do.call(rbind, lapply(groups, function(ix) {
      z <- a[ix, ]; data.frame(context = z$context[1], category = z[[column]][1],
        n_edges = nrow(z), n_primary_edges = sum(z$primary_edge),
        median_abs_delta_r = stats::median(z$abs_delta_r),
        sign_change_fraction = mean(z$sign_change),
        all_49_same_delta_direction_fraction = mean(z$loco_all_49_same_delta_direction),
        all_49_sign_reversal_fraction = mean(z$loco_all_49_sign_reversal))
    }))
  }
  write_csv(change_summary("gene"), "gene_diet_change_summary.csv")
  write_csv(change_summary("lipid_class"), "lipid_class_diet_change_summary.csv")
  write_csv(analysis$edges[order(-analysis$edges$abs_delta_r), ], "descriptive_ranked_diet_changes.csv")
  dir.create(file.path(out, "networks"), showWarnings = FALSE)
  for (name in names(contexts)) for (diet in c("SD", "WD")) for (universe in c("full", "primary")) {
    a <- analysis$edges[analysis$edges$context == name &
                         (universe == "full" | analysis$edges$primary_edge), ]
    sif <- data.frame(source = a$gene, target = a$reference_lipid, corr = a[[paste0("r_", diet)]],
      feature_id = a$feature_id, n_mice = 7L, context = name, diet = diet,
      original_displayed_edge = a$primary_edge)
    write_csv(sif, paste0("networks/", name, "_", diet, "_", universe, "_SIF.csv"))
  }
  source("scripts/core/wdsd_presentation.R")
  wp2_make_plots(out, data, contexts, analysis, sensitivity_table)
  wp2_make_report(out)
  paths <- unique(c(unname(data$input_paths), attr(reference, "provenance")$inputs,
    "scripts/support_functions.R", "scripts/core/cyp2c38_reference.R",
    "scripts/core/wdsd_intake.R", "scripts/core/wdsd_reference.R",
    "scripts/core/wdsd_comparison_setup.R", "scripts/core/wdsd_network_statistics.R",
    "scripts/core/wdsd_presentation.R", "scripts/workflows/12_wdsd_network_comparison.R", "scripts/core/wp2_workflow.R",
    "scripts/reports/build_wdsd_wp2_pdf.py"))
  hash <- function(path) {
    value <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
    if (length(value) != 1L || !grepl("^[a-f0-9]{64} ", value)) stop("SHA256 failed for ", path)
    substr(value, 1L, 64L)
  }
  write_csv(data.frame(path = normalizePath(paths), SHA256 = vapply(paths, hash, character(1L))),
            "input_and_code_SHA256.csv")
  jsonlite::write_json(list(date_UTC = format(Sys.time(), tz = "UTC"), args = args,
    settings = analysis$metadata, focus_genes = focus,
    primary_interpretation = "frozen displayed edges; full universe complementary; four context comparisons",
    fill_sensitivity = "QC sensitivity added after input inspection and before association analysis",
    biological_units = "28 mice; 7 per sex/diet; shared bulk lipid profile across paired RNA zones",
    interpretation = "Native ULM p descriptive; exact pairing and marginal bootstrap impose no RNA-lipid association; diet changes descriptive"),
    file.path(out, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("Completed WP2; results: ", normalizePath(out), "\n", sep = "")
  invisible(list(analysis = analysis, sensitivities = sensitivity_table, setup = setup))
}

# Complete work package: pooled Figure 3 first, then sex/zonation and diet.
wp2_package_run <- function(args = commandArgs(trailingOnly = TRUE)) {
  allowed <- "^--(permutations|bootstraps|paired-bootstraps|seed|output)=|^--(validate-only|report-only)$"
  if (any(!grepl(allowed, args))) stop("Unknown Workflow 12 option")
  if (all(c("--validate-only", "--report-only") %in% args)) stop("Choose validation or report refresh")
  option <- function(key, default) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (length(hit) > 1L) stop("Duplicate option: ", key)
    if (length(hit)) substring(hit, nchar(key) + 4L) else default
  }
  out <- option("output", analysis_path("historical_wf12"))
  make_report <- function() {
    source("scripts/core/wp2_report.R")
    wp2_write_package_report(out)
  }
  if ("--report-only" %in% args) {
    cat("\n2. REFRESH POOLED-FIRST WORK PACKAGE REPORT\n")
    make_report()
    return(invisible(list(output = out, mode = "report-only")))
  }
  common <- c(paste0("--bootstraps=", option("bootstraps", "4999")),
    paste0("--paired-bootstraps=", option("paired-bootstraps", "999")),
    paste0("--seed=", option("seed", "20261007")))
  pooled_args <- c(common, paste0("--permutations=", option("permutations", "19999")),
    paste0("--output=", file.path(out, "followup")))
  zone_args <- c(common, paste0("--output=", out))
  if ("--validate-only" %in% args) {
    cat("\n2. VALIDATE POOLED FIGURE 3 AND LIPID ANNOTATION AUDIT\n")
    pooled <- wp2_pooled_run(c(pooled_args, "--validate-only"))
    cat("\n3. VALIDATE ZONATION REFERENCES AND DIET PAIRING\n")
    zonation <- wp2_run(c(zone_args, "--validate-only"))
    cat("Workflow 12 validation passed without reading saved WP2 outputs.\n")
    return(invisible(list(pooled = pooled, zonation = zonation)))
  }
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; choose a fresh --output destination or --report-only.")
  }
  cat("\n2. POOLED FIGURE 3 CONSERVATION AND LIPID ANNOTATION AUDIT\n")
  pooled <- wp2_pooled_run(pooled_args)
  cat("\n3. SEX/ZONATION CONSERVATION AND WESTERN-DIET COMPARISON\n")
  zonation <- wp2_run(zone_args, allow_pooled = TRUE)
  cat("\n4. WRITE POOLED-FIRST WORK PACKAGE REPORT\n")
  make_report()
  jsonlite::write_json(list(workflow = "12_wdsd_network_comparison",
    phases = c("pooled Figure 3 and lipid overlap", "sex/zonation and diet", "combined report"),
    args = args, pooled_args = pooled_args, zonation_args = zone_args,
    finished_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(out, "workflow_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  cat("Workflow 12 complete: ", normalizePath(out), "\n", sep = "")
  invisible(list(pooled = pooled, zonation = zonation))
}
