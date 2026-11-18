source("scripts/config/analysis_release.R")
# Complete WP3 orchestration; statistical modules define functions only.
wp3_check_native_ULM <- function(contexts, inference = NULL) {
  records <- list()
  for (name in names(contexts)) {
    ct <- contexts[[name]]; e <- ct$edges
    a <- wp2_cor_edges(ct$rna_CF453, ct$lipid_CF453, e)
    b <- wp2_cor_edges(ct$rna_CF459, ct$lipid_CF459, e)
    for (universe in c("full", "displayed")) for (scope in c("network", "Mtmr11")) {
      use <- which((universe == "full" | e$primary_edge) & (scope == "network" | e$gene == scope))
      if (length(use) < 3L) next
      native <- wp2_vector_stats(a[use], b[use])
      net <- data.frame(source = "CF453_WD", target = e$edge_id[use], mor = a[use])
      mat <- matrix(b[use], ncol = 1L, dimnames = list(e$edge_id[use], "CF459_WD"))
      ulm <- decoupleR::run_ulm(mat, net, minsize = 3L)
      if (nrow(ulm) != 1L || max(abs(c(native["ULM_score"] - ulm$score,
               native["ULM_p_value"] - ulm$p_value))) > 1e-10) stop("Native ULM disagrees with decoupleR: ", name)
      spec <- paste(name, universe, if (scope == "network") "network" else paste0("gene::", scope), sep = "__")
      if (!is.null(inference)) {
        row <- inference[inference$specification == spec, , drop = FALSE]
        if (nrow(row) != 1L || max(abs(c(row$ULM_score - ulm$score, row$ULM_p_value - ulm$p_value))) > 1e-10) stop("Exported ULM does not replay: ", spec)
      }
      records[[spec]] <- data.frame(specification = spec, n_edges = length(use),
        native_ULM = unname(native["ULM_score"]), decoupleR_ULM = ulm$score,
        native_p = unname(native["ULM_p_value"]), decoupleR_p = ulm$p_value,
        max_abs_difference = max(abs(c(native["ULM_score"] - ulm$score, native["ULM_p_value"] - ulm$p_value))))
    }
  }
  do.call(rbind, records)
}

wp3_run <- function(args = commandArgs(trailingOnly = TRUE)) {
  allowed <- "^--(output|comparison-output|seed|bootstraps|paired-bootstraps)=|^--(validate-only|report-only|fcg-comparison-only)$"
  if (any(!grepl(allowed, args)) || anyDuplicated(args)) stop("Unknown or repeated Workflow 13 option.")
  option <- function(key, fallback) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (length(hit) > 1L) stop("Duplicate option: ", key)
    if (length(hit)) substring(hit, nchar(key) + 4L) else fallback
  }
  modes <- intersect(args, c("--validate-only", "--report-only", "--fcg-comparison-only"))
  if (length(modes) > 1L) stop("Choose validation, report-only or fcg-comparison-only mode.")
  if (any(startsWith(args, "--comparison-output=")) && !"--fcg-comparison-only" %in% args) stop("--comparison-output requires --fcg-comparison-only.")
  out <- option("output", analysis_path("historical_wf13"))
  if (!nzchar(out)) stop("Output directory must not be empty.")
  if ("--report-only" %in% args) {
    source("scripts/core/wd_crosscohort_presentation.R")
    wp3_write_report(out)
    return(invisible(list(output = out, mode = "report-only")))
  }
  for (f in c("cyp2c38_reference.R", "wdsd_intake.R", "wd_intake.R", "wdsd_reference.R",
              "wdsd_comparison_setup.R", "wdsd_network_statistics.R", "wd_crosscohort_setup.R",
              "wd_crosscohort_statistics.R", "wd_fcg_comparison.R")) source(file.path("scripts/core", f))
  B <- wp3_count(as.numeric(option("bootstraps", "4999")), "bootstraps", 19L)
  paired_B <- wp3_count(as.numeric(option("paired-bootstraps", "999")), "paired-bootstraps", 19L)
  seed <- wp3_count(as.numeric(option("seed", "20261014")), "seed", 0L)
  for (p in c("SeuratObject", "Matrix", "decoupleR", "jsonlite", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  if ("--fcg-comparison-only" %in% args) {
    if (any(grepl("^--(bootstraps|paired-bootstraps)=", args))) stop("Extension reuses saved draw banks; do not supply bootstrap counts.")
    required <- file.path(out, c("contexts.rds", "analysis.rds"))
    if (any(!file.exists(required))) stop("Extension requires saved original WP3 contexts.rds and analysis.rds.")
    extension_out <- option("comparison-output", file.path(out, "fcg_comparison"))
    wp3_fcg_output_guard(extension_out)
    contexts <- readRDS(required[1L]); analysis <- readRDS(required[2L])
    reference <- wp2_load_reference()
    comparison <- wp3_fcg_comparison(contexts, reference, analysis, seed = seed)
    wp3_save_fcg_comparison(comparison, reference, extension_out, args, required)
    return(invisible(list(output = extension_out, mode = "fcg-comparison-only", inference = comparison$inference)))
  }
  if (!"--validate-only" %in% args && file.exists(out) &&
      (!dir.exists(out) || length(list.files(out, all.files = TRUE, no.. = TRUE)))) stop("Output directory is not empty; select a fresh --output or explicit --report-only.")
  cat("\n2. INTAKE, FIXED REFERENCES AND STRICT MATCHING\n")
  cf453 <- wp2_load_inputs(input_zonation=analysis_path("historical_wdsd_zonation_transport")); cf459 <- wd_load_inputs(); reference <- wp2_load_reference()
  setup <- wp3_prepare_contexts(cf453, cf459, reference)
  check <- wp3_check_native_ULM(setup$contexts)
  if ("--validate-only" %in% args) {
    cat("Workflow 13 validation passed; no results written.\n")
    return(invisible(list(setup = setup, native_ULM = check)))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (d in c("QA", "plots", "networks")) dir.create(file.path(out, d), showWarnings = FALSE)
  write_csv <- function(x, file) utils::write.csv(x, file.path(out, file), row.names = FALSE, na = "")
  write_csv(setup$coverage, "coverage.csv")
  write_csv(setup$lipid_crosswalk, "lipid_crosswalk.csv")
  write_csv(setup$edge_eligibility_ledger, "edge_eligibility_ledger.csv")
  write_csv(setup$WP2_replay, "WP2_CF453_WD_replay.csv")
  write_csv(attr(reference, "summary"), "original_FCG_reference_replay.csv")
  saveRDS(setup$contexts, file.path(out, "contexts.rds"))
  for (cohort in c("CF453", "CF459")) {
    x <- if (cohort == "CF453") cf453 else cf459
    for (field in c("metadata", "lipid_ledger", "qc_summary", "rna_profile_qc", "layer_inventory", "marker_by_profile", "zonation_marker_diagnostic")) {
      write_csv(x[[field]], paste0(cohort, "_", field, ".csv"))
    }
  }
  cat("\n3. CROSS-COHORT NETWORKS, EXACT NULLS AND STABILITY\n")
  analysis <- wp3_crosscohort_analysis(setup$contexts, B = B, paired_B = paired_B, seed = seed)
  check <- wp3_check_native_ULM(setup$contexts, analysis$inference)
  write_csv(check, "QA/native_ULM_replay.csv")
  for (field in names(analysis)) if (is.data.frame(analysis[[field]])) {
    file <- switch(field, edges = "network_edges.csv", inference = "crosscohort_inference.csv", paste0(field, ".csv"))
    write_csv(analysis[[field]], file)
  }
  saveRDS(analysis, file.path(out, "analysis.rds"))
  write_csv(wp3_feature_support(setup$contexts, analysis$resampling_index_matrices), "feature_bootstrap_degeneracy.csv")
  for (name in names(setup$contexts)) for (cohort in c("CF453", "CF459")) for (universe in c("full", "displayed")) {
    a <- analysis$edges[analysis$edges$context == name & (universe == "full" | analysis$edges$primary_edge), , drop = FALSE]
    sif <- data.frame(source = a$gene, target = a$reference_lipid,
      corr = a[[paste0("r_", cohort, "_WD")]], n_mice = if (cohort == "CF453") 7L else 4L,
      original_displayed_edge = a$primary_edge, feature_id_CF453 = a$feature_id_CF453, feature_id_CF459 = a$feature_id_CF459)
    write_csv(sif, paste0("networks/", name, "_", cohort, "_WD_", universe, "_SIF.csv"))
  }
  cat("\n4. PRESPECIFIED ONE-COHORT QC AND STORED RNA-LAYER SENSITIVITIES\n")
  sensitivity_fits <- list()
  for (type in c("CF453_QC30", "CF453_Fill50", "RNA_data")) {
    cat("Sensitivity: ", type, "\n", sep = "")
    ct <- wp3_context_sensitivity(setup$contexts, cf453, cf459, type)
    fit <- wp3_crosscohort_analysis(ct, B = B, paired_B = paired_B, seed = seed)
    fit$inference$sensitivity <- type
    sensitivity_fits[[type]] <- fit
  }
  sensitivities <- do.call(rbind, lapply(sensitivity_fits, `[[`, "inference"))
  write_csv(sensitivities, "sensitivities.csv")
  saveRDS(sensitivity_fits, file.path(out, "sensitivity_analyses.rds"))
  cat("\n5. ORIGINAL MALE FCG SD ANCHOR AND DIRECT CLOSENESS\n")
  fcg_comparison <- wp3_fcg_comparison(setup$contexts, reference, analysis, seed = seed)
  wp3_save_fcg_comparison(fcg_comparison, reference, file.path(out, "fcg_comparison"), args,
                          file.path(out, c("contexts.rds", "analysis.rds")))
  paths <- unique(c(unname(cf453$input_paths), unname(cf459$input_paths),
    attr(reference, "provenance")$inputs, setup$WP2_replay_path,
    "data/wp3_input_migration_manifest.csv", "docs/workflow_plans/13_wd_crosscohort_robustness.md",
    "scripts/support_functions.R", "scripts/core/cyp2c38_reference.R", "scripts/core/wdsd_intake.R",
    "scripts/core/wdsd_reference.R", "scripts/core/wdsd_comparison_setup.R", "scripts/core/wdsd_network_statistics.R",
    "scripts/core/wd_intake.R", "scripts/core/wd_crosscohort_setup.R", "scripts/core/wd_crosscohort_statistics.R",
    "scripts/core/wd_crosscohort_workflow.R", "scripts/core/wd_fcg_comparison.R", "scripts/core/wd_crosscohort_presentation.R",
    "scripts/workflows/13_wd_crosscohort_robustness.R", "scripts/reports/build_wp3_report_pdf.py"))
  hash <- function(p) {
    value <- system2("shasum", c("-a", "256", shQuote(p)), stdout = TRUE)
    if (length(value) != 1L || !grepl("^[a-f0-9]{64} ", value)) stop("SHA256 failed: ", p)
    substr(value, 1L, 64L)
  }
  write_csv(data.frame(path = normalizePath(paths, mustWork = TRUE), SHA256 = vapply(paths, hash, character(1L))), "input_and_code_SHA256.csv")
  jsonlite::write_json(list(workflow = "13_wd_crosscohort_robustness", args = args,
    completed_analysis_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), settings = analysis$metadata,
    primary_contexts = names(setup$contexts), sensitivities = names(sensitivity_fits),
    predefined_gene = "Mtmr11", imported_inputs = "data/wp3_input_migration_manifest.csv",
    CF459_WD_provenance = "Supplied WD-only cohort/meeting designation; f-ctrl/m-ctrl metadata is not SD",
    cohort_ages_months_user_20261014 = list(CF459 = 3L, original_FCG = 3L, CF453 = 10L),
    cohort_metadata_unknown = c( "strain_background", "diet_formulation_and_exposure", "batch_litter", "upstream_processing"),
    limitations = unique(c(cf453$limitations, cf459$limitations))), file.path(out, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\n6. SCIENTIFIC FIGURES AND SELF-CONTAINED REPORT\n")
  source("scripts/core/wd_crosscohort_presentation.R")
  wp3_make_plots(out, cf453, cf459, setup, analysis, sensitivities, fcg_comparison = fcg_comparison)
  wp3_write_report(out)
  cat("Workflow 13 complete: ", normalizePath(out), "\n", sep = "")
  print(analysis$inference[analysis$inference$scope == "network", c("context", "universe", "coefficient_vector_r", "ULM_score", "ULM_p_value", "pairing_CF459_given_CF453_p", "pairing_CF453_given_CF459_p", "null_bootstrap_p")])
  invisible(list(output = out, inference = analysis$inference, sensitivities = sensitivities, fcg_comparison = fcg_comparison$inference))
}

# Diagnose support/constant-feature behavior without changing the frozen universe.
# These counts explain undefined draws; they are neither an exclusion rule nor
# an additional hypothesis test.
wp3_feature_support <- function(contexts, indices) {
  records <- list()
  for (name in names(contexts)) for (cohort in c("CF453", "CF459")) for (modality in c("rna", "lipid")) {
    x <- contexts[[name]][[paste0(modality, "_", cohort)]]
    null_index <- if (cohort == "CF459") indices[[paste0("CF459_null_", modality)]] else NULL
    paired_index <- indices[[paste0(cohort, "_paired")]]
    constant_counts <- function(bank) {
      if (is.null(bank)) return(rep(NA_integer_, nrow(x)))
      result <- integer(nrow(x))
      for (j in seq_len(nrow(bank))) result <- result + !wp2_unit_rows(x[, bank[j, ], drop = FALSE])$valid
      result
    }
    records[[paste(name, cohort, modality)]] <- data.frame(context = name, cohort = cohort,
      modality = modality, feature = rownames(x), n_observed_mice = ncol(x),
      observed_unique_values = apply(x, 1L, function(z) length(unique(z))),
      null_draws = if (is.null(null_index)) 0L else nrow(null_index),
      null_constant_draws = constant_counts(null_index), paired_draws = nrow(paired_index),
      paired_constant_draws = constant_counts(paired_index), stringsAsFactors = FALSE)
  }
  do.call(rbind, records)
}


# A saved extension never overwrites either an existing extension or base outputs.
wp3_fcg_output_guard <- function(out) {
  if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out)) stop("Comparison output must be a nonempty path.")
  if (file.exists(out) && (!dir.exists(out) || length(list.files(out, all.files = TRUE, no.. = TRUE)))) stop("Comparison output directory is not empty; choose a fresh --comparison-output.")
  invisible(TRUE)
}

wp3_save_fcg_comparison <- function(comparison, reference, out, args = character(), base_files = character()) {
  wp3_fcg_output_guard(out)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(comparison, file.path(out, "analysis.rds"))
  for (field in names(comparison)) if (is.data.frame(comparison[[field]])) {
    utils::write.csv(comparison[[field]], file.path(out, paste0(field, ".csv")), row.names = FALSE, na = "")
  }
  saveRDS(comparison$FCG_paired_index_matrix, file.path(out, "FCG_paired_index_matrix.rds"))
  saveRDS(reference[c("male_periportal", "male_pericentral")], file.path(out, "reference_contexts.rds"))
  paths <- unique(c(base_files, attr(reference, "provenance")$inputs,
    "scripts/core/wd_fcg_comparison.R", "scripts/core/wd_crosscohort_workflow.R",
    "scripts/core/wd_crosscohort_statistics.R", "scripts/core/wdsd_network_statistics.R",
    "scripts/core/wdsd_reference.R", "scripts/workflows/13_wd_crosscohort_robustness.R",
    "docs/workflow_plans/13_wd_crosscohort_robustness.md"))
  hash <- function(p) {
    value <- system2("shasum", c("-a", "256", shQuote(p)), stdout = TRUE)
    if (length(value) != 1L || !grepl("^[a-f0-9]{64} ", value)) stop("SHA256 failed: ", p)
    substr(value, 1L, 64L)
  }
  utils::write.csv(data.frame(path = normalizePath(paths, mustWork = TRUE),
    SHA256 = vapply(paths, hash, character(1L))), file.path(out, "input_and_code_SHA256.csv"), row.names = FALSE)
  jsonlite::write_json(list(workflow = "13_wd_crosscohort_robustness_FCG_extension",
    completed_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), args = args,
    base_files_preserved = normalizePath(base_files, mustWork = TRUE), settings = comparison$metadata),
    file.path(out, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  invisible(out)
}
