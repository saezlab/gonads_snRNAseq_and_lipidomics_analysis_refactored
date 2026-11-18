#!/usr/bin/env Rscript
# Parameter search for relaxed Factor 1 lipid-to-TF candidates.
# Fits are score-first and use only strict consecutive score-layer extraction.
options(stringsAsFactors = FALSE, digits = 17)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run this entry point with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../.."), mustWork = TRUE))

parse_cli <- function(args) {
  value_keys <- c("config", "output", "pruning-steps", "layers", "tf-panels", "primary", "secondary")
  flag_keys <- c("help", "validate-only")
  parsed <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) stop("Expected --option=value, received: ", arg)
    key <- sub("=.*$", "", sub("^--", "", arg))
    if (!key %in% c(value_keys, flag_keys)) stop("Unknown option: --", key)
    if (key %in% names(parsed)) stop("Repeated option: --", key)
    if (key %in% flag_keys) {
      if (grepl("=", arg, fixed = TRUE)) stop("Flag takes no value: --", key)
      parsed[[key]] <- TRUE
    } else {
      if (!grepl("=", arg, fixed = TRUE)) stop("Use --", key, "=VALUE")
      value <- sub("^[^=]*=", "", arg)
      if (!nzchar(value)) stop("Empty value for --", key)
      parsed[[key]] <- value
    }
  }
  parsed
}

parse_numbers <- function(value, name, integer = FALSE, strictly_positive = FALSE) {
  if (grepl("(^,|,$|,,)", value)) stop("Empty numeric entry in --", name)
  items <- strsplit(value, ",", fixed = TRUE)[[1L]]
  values <- suppressWarnings(as.numeric(items))
  if (!length(values) || any(!is.finite(values)) || any(values < 0) ||
      (strictly_positive && any(values <= 0)) ||
      (integer && any(values != floor(values) | values > .Machine$integer.max))) {
    stop("Invalid --", name)
  }
  sort(unique(if (integer) as.integer(values) else values))
}

parse_panels <- function(value) {
  panels <- sort(unique(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(panels) || any(!panels %in% c("all_tfs", "strong_tfs"))) {
    stop("--tf-panels must contain one or both of all_tfs,strong_tfs")
  }
  panels
}

frontier <- function(x, objectives) {
  if (!nrow(x)) return(logical())
  values <- as.matrix(x[, objectives, drop = FALSE])
  !vapply(seq_len(nrow(x)), function(i) {
    at_least <- rowSums(sweep(values, 2L, values[i, ], `>=`)) == ncol(values)
    strictly <- rowSums(sweep(values, 2L, values[i, ], `>`)) > 0L
    any(at_least & strictly)
  }, logical(1))
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(cli$help)) {
  cat(paste(
    "Usage: Rscript scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R [OPTIONS]",
    "  --config=PATH             Default: scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json",
    "  --output=NEW_DIRECTORY    New/empty result directory",
    "  --pruning-steps=4,6,8     Override configured directed pruning distances",
    "  --layers=2,4,6            Override configured recursive ULM depths",
    "  --tf-panels=all_tfs,strong_tfs  Override configured observation panels",
    "  --primary=0.5,1           Override primary score thresholds",
    "  --secondary=0,0.5         Override secondary score thresholds",
    "  --validate-only            Validate grids and frozen dependencies without fitting",
    "  --help                     Show this message",
    "Each topology/depth/panel fit is calculated once. Every valid threshold pair is recorded.",
    "Only strict consecutive score-layer paths can be materialized; ordinary sign-consistent connectivity is excluded from this search.",
    "Compact networks require <500 edges, a nonzero secondary threshold, at least one independent direct-edge fingerprint group, and at least two observed TFs.",
    "Identical SIF topologies are materialized once and summarized in score_layer_stability.csv. No automatic optimum is selected.",
    sep = "\n"), "\n")
  quit(save = "no")
}

source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
source("scripts/moon/lib/relaxed_lipid_tf.R")

search_config_path <- if (is.null(cli$config)) {
  "scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json"
} else cli$config
search_cfg <- relaxed_read_search_config(search_config_path)
base_cfg <- relaxed_read_config(search_cfg$base_config)
get_value <- function(key, default) if (is.null(cli[[key]])) default else cli[[key]]
pruning <- parse_numbers(get_value("pruning-steps", paste(search_cfg$pruning_steps, collapse = ",")),
                         "pruning-steps", integer = TRUE, strictly_positive = TRUE)
layers <- parse_numbers(get_value("layers", paste(search_cfg$n_layers, collapse = ",")),
                        "layers", integer = TRUE, strictly_positive = TRUE)
panels <- parse_panels(get_value("tf-panels", paste(search_cfg$tf_panels, collapse = ",")))
primary <- parse_numbers(get_value("primary", paste(search_cfg$primary_thresholds, collapse = ",")),
                         "primary", strictly_positive = TRUE)
secondary <- parse_numbers(get_value("secondary", paste(search_cfg$secondary_thresholds, collapse = ",")),
                           "secondary")
thresholds <- expand.grid(primary = primary, secondary = secondary, KEEP.OUT.ATTRS = FALSE)
thresholds$threshold_id <- sprintf("T%04d", seq_len(nrow(thresholds)))
thresholds$valid_pair <- thresholds$primary >= thresholds$secondary
if (!any(thresholds$valid_pair)) stop("No valid threshold pairs: primary must be at least secondary")
fit_grid <- expand.grid(tf_panel = panels, pruning_steps = pruning, n_layers = layers,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
fit_grid$run_id <- paste("direct", fit_grid$tf_panel, paste0("pruning", fit_grid$pruning_steps),
                         paste0("layers", fit_grid$n_layers), sep = "_")
output_path <- get_value("output", search_cfg$output)
data <- relaxed_load_data(base_cfg)
if (isTRUE(cli[["validate-only"]])) {
  cat("Validated ", nrow(fit_grid), " relaxed fits and ", sum(thresholds$valid_pair),
      " valid strict threshold pairs per fit (", nrow(fit_grid) * sum(thresholds$valid_pair),
      " total). No files written.\n", sep = "")
  quit(save = "no")
}
if (dir.exists(output_path) && length(list.files(output_path, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; choose a new --output= directory: ", output_path)
}

policy_cfg <- base_cfg
policy_cfg$maximum_solution_edges_exclusive <- search_cfg$maximum_solution_edges_exclusive
policy_cfg$minimum_independent_edge_fingerprint_groups <- search_cfg$minimum_independent_edge_fingerprint_groups
policy_cfg$minimum_downstream_tfs <- search_cfg$minimum_downstream_tfs
out <- moon_new_output(output_path)
relaxed_write_common_inputs(data, base_cfg, search_cfg$base_config, out)
invisible(file.copy(search_config_path, file.path(out, "parameter_search_config.json"), overwrite = FALSE))
moon_write_csv(fit_grid, file.path(out, "fit_grid.csv"))
moon_write_csv(thresholds, file.path(out, "threshold_grid.csv"))
dependency_paths <- unique(c(
  search_config_path, search_cfg$base_config, names(base_cfg$input_md5),
  "scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R",
  "scripts/moon/lib/relaxed_lipid_tf.R", "scripts/moon/lib/workflow.R", "scripts/moon/moon_safe_helpers.R"
))
moon_write_csv(data.frame(path = dependency_paths, md5 = unname(tools::md5sum(dependency_paths))),
               file.path(out, "dependency_manifest.csv"))
relaxed_write_json(list(
  analysis = search_cfg$analysis, base_analysis = base_cfg$analysis,
  pruning_steps = pruning, n_layers = layers, tf_panels = panels,
  primary_thresholds = primary, secondary_thresholds = secondary,
  threshold_policy = "primary >= secondary >= 0; strict absolute cutoffs use >",
  extraction_policy = "Strict score-layer extraction only: source level = target level + 1, signed-edge consistency, source-to-observed-TF walk, and copied-measurement sign agreement.",
  eligibility_policy = "0 < edges < 500; secondary > 0; at least one independent direct-edge fingerprint group; at least two observed TFs.",
  duplicate_policy = "Identical signed SIFs are materialized once and all threshold/run occurrences are retained in score_layer_stability.csv.",
  selection_policy = search_cfg$selection_policy,
  command_args = commandArgs(trailingOnly = TRUE), started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(out, "search_config.json"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

fit_rows <- grid_rows <- status_rows <- node_rows <- list()
materialized <- new.env(parent = emptyenv())
persist <- function() {
  if (length(fit_rows)) moon_write_csv(do.call(rbind, fit_rows), file.path(out, "fit_results.csv"))
  if (length(grid_rows)) moon_write_csv(do.call(rbind, grid_rows), file.path(out, "grid_results.csv"))
  if (length(status_rows)) moon_write_csv(do.call(rbind, status_rows), file.path(out, "candidate_assignment_status_all_runs.csv"))
  if (length(node_rows)) moon_write_csv(do.call(rbind, node_rows), file.path(out, "candidate_node_status_all_runs.csv"))
}

for (i in seq_len(nrow(fit_grid))) {
  current <- fit_grid[i, , drop = FALSE]
  run_dir <- file.path(out, "runs", current$run_id)
  dir.create(run_dir, recursive = TRUE)
  panel <- relaxed_panel_data(data, current$tf_panel)
  started <- proc.time()[["elapsed"]]
  message("Starting relaxed fit: ", current$run_id)
  fit <- tryCatch(relaxed_fit(panel, current$pruning_steps, current$n_layers, base_cfg$max_iterations), error = identity)
  failed <- inherits(fit, "error")
  fit_error <- if (failed) conditionMessage(fit) else ""
  fit_status <- if (failed && grepl("did not converge", fit_error, fixed = TRUE)) "nonconvergent" else if (failed) "failed" else "converged"
  status <- node_status <- seeds <- NULL
  if (!failed) {
    moon_save_fit(fit, run_dir)
    moon_write_csv(fit$scores, file.path(run_dir, "moon_scores_unconstrained.csv"))
    status <- relaxed_candidate_status(data$assignments, fit)
    node_status <- relaxed_candidate_node_summary(status)
    seeds <- relaxed_sign_consensus_seeds(status)
    moon_write_csv(status, file.path(run_dir, "candidate_assignment_status.csv"))
    moon_write_csv(node_status, file.path(run_dir, "candidate_node_status.csv"))
    moon_write_csv(data.frame(source = names(seeds), copied_measurement_sign = unname(seeds)),
                   file.path(run_dir, "candidate_sign_consensus_seeds.csv"))
    status_rows[[length(status_rows) + 1L]] <- cbind(run_id = current$run_id, status)
    node_rows[[length(node_rows) + 1L]] <- cbind(run_id = current$run_id, node_status)
  }
  fit_rows[[length(fit_rows) + 1L]] <- data.frame(
    run_id = current$run_id, tf_panel = current$tf_panel,
    pruning_steps = current$pruning_steps, n_layers = current$n_layers,
    status = fit_status, converged = !failed, validated = !failed,
    iterations = if (failed) NA_integer_ else nrow(fit$coherence_ledger),
    final_scoring_edges = if (failed) NA_integer_ else nrow(fit$pkn_scoring),
    final_uncompressed_edges = if (failed) NA_integer_ else nrow(fit$pkn_final),
    scored_nodes = if (failed) NA_integer_ else nrow(fit$scores),
    sign_consensus_candidate_nodes = if (failed) NA_integer_ else length(seeds),
    elapsed_seconds = proc.time()[["elapsed"]] - started, error = fit_error,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(nrow(thresholds))) {
    cutoff <- thresholds[j, , drop = FALSE]
    row <- data.frame(
      run_id = current$run_id, tf_panel = current$tf_panel,
      pruning_steps = current$pruning_steps, n_layers = current$n_layers,
      threshold_id = cutoff$threshold_id, primary = cutoff$primary, secondary = cutoff$secondary,
      status = "pending", converged = !failed, validated = FALSE,
      score_layer_edges = NA_integer_, score_layer_nodes = NA_integer_,
      raw_upstream_candidate_nodes = NA_integer_, independent_edge_fingerprint_groups = NA_integer_,
      observed_tfs = NA_integer_, measurement_features = "", upstream_source_nodes = "",
      edge_fingerprint_groups = "", topology_hash = "", materialized = FALSE,
      duplicate_of_solution_ref = "", solution_ref = "", eligible = FALSE, error = "",
      stringsAsFactors = FALSE
    )
    if (!cutoff$valid_pair) {
      row$status <- "invalid_threshold_pair"
      row$error <- "primary must be at least secondary"
    } else if (failed) {
      row$status <- paste0("fit_", fit_status)
      row$error <- fit_error
    } else {
      solution <- tryCatch({
        value <- relaxed_extract_score_layer_paths(
          fit$scores, fit$pkn_final, seeds, cutoff$primary, cutoff$secondary, data$rna,
          pkn_already_checked = TRUE
        )
        relaxed_validate_score_layer_solution(value, fit, seeds, cutoff$primary, cutoff$secondary)
        value
      }, error = identity)
      if (inherits(solution, "error")) {
        row$status <- "extraction_failed"
        row$error <- conditionMessage(solution)
      } else {
        metadata <- relaxed_solution_metadata(solution, status)
        eligibility <- relaxed_solution_is_eligible(solution, status, policy_cfg)
        row$status <- if (nrow(solution$SIF)) "ok" else "empty"
        row$validated <- TRUE
        row$score_layer_edges <- nrow(solution$SIF)
        row$score_layer_nodes <- nrow(solution$ATT)
        row$raw_upstream_candidate_nodes <- metadata$raw_upstream_candidate_nodes
        row$independent_edge_fingerprint_groups <- metadata$independent_edge_fingerprint_groups
        row$observed_tfs <- metadata$retained_observed_tfs
        row$measurement_features <- relaxed_collapse(metadata$measurement_features)
        row$upstream_source_nodes <- relaxed_collapse(metadata$source_nodes)
        row$edge_fingerprint_groups <- relaxed_collapse(metadata$edge_fingerprint_groups)
        row$eligible <- isTRUE(eligibility$eligible)
        if (row$eligible) {
          row$topology_hash <- relaxed_sif_hash(solution$SIF)
          if (exists(row$topology_hash, envir = materialized, inherits = FALSE)) {
            row$duplicate_of_solution_ref <- get(row$topology_hash, envir = materialized, inherits = FALSE)
            row$solution_ref <- row$duplicate_of_solution_ref
          } else {
            row$solution_ref <- file.path("runs", current$run_id, "solutions", cutoff$threshold_id)
            solution_dir <- file.path(out, row$solution_ref)
            dir.create(solution_dir, recursive = TRUE)
            solution$ATT <- relaxed_annotate_solution(solution, status)
            saveRDS(solution, file.path(solution_dir, "solution.rds"))
            moon_write_csv(solution$SIF, file.path(solution_dir, "network_SIF.csv"))
            moon_write_csv(solution$ATT, file.path(solution_dir, "network_ATT.csv"))
            relaxed_write_json(list(
              extraction = solution$audit,
              candidate_metadata = metadata,
              eligibility = eligibility,
              run_id = current$run_id,
              tf_panel = current$tf_panel,
              pruning_steps = current$pruning_steps,
              n_layers = current$n_layers,
              topology_hash = row$topology_hash,
              interpretation = "Strict layer-consecutive score-path sensitivity; not a calibrated causal confidence estimate."
            ), file.path(solution_dir, "extraction_audit.json"))
            assign(row$topology_hash, row$solution_ref, envir = materialized)
            row$materialized <- TRUE
          }
        }
      }
    }
    grid_rows[[length(grid_rows) + 1L]] <- row
  }
  persist()
  rm(fit)
  invisible(gc())
}

grid <- do.call(rbind, grid_rows)
eligible <- grid[grid$eligible, , drop = FALSE]
if (nrow(eligible)) {
  eligible$coverage_threshold_frontier <- frontier(eligible,
    c("primary", "secondary", "observed_tfs", "independent_edge_fingerprint_groups"))
  eligible$threshold_frontier <- frontier(eligible, c("primary", "secondary"))
  eligible <- eligible[order(-eligible$secondary, -eligible$primary, -eligible$observed_tfs,
                             -eligible$independent_edge_fingerprint_groups, eligible$score_layer_edges), , drop = FALSE]
} else {
  eligible$coverage_threshold_frontier <- logical()
  eligible$threshold_frontier <- logical()
}
moon_write_csv(eligible, file.path(out, "eligible_under_500.csv"))
moon_write_csv(eligible[eligible$coverage_threshold_frontier, , drop = FALSE],
               file.path(out, "coverage_threshold_frontier.csv"))
moon_write_csv(eligible[eligible$threshold_frontier, , drop = FALSE],
               file.path(out, "threshold_frontier.csv"))

if (nrow(eligible)) {
  grouped <- split(seq_len(nrow(eligible)), eligible$topology_hash)
  stability <- do.call(rbind, lapply(names(grouped), function(hash) {
    x <- eligible[grouped[[hash]], , drop = FALSE]
    data.frame(
      topology_hash = hash,
      n_grid_rows = nrow(x),
      materialized_solution_ref = x$solution_ref[which(x$materialized)[1L]],
      tf_panels = relaxed_collapse(x$tf_panel),
      pruning_steps = relaxed_collapse(x$pruning_steps),
      n_layers = relaxed_collapse(x$n_layers),
      threshold_pairs = relaxed_collapse(paste(x$primary, x$secondary, sep = "/")),
      measurement_features = relaxed_collapse(x$measurement_features),
      upstream_source_nodes = relaxed_collapse(x$upstream_source_nodes),
      observed_tfs = relaxed_collapse(x$observed_tfs),
      edges = x$score_layer_edges[1L], nodes = x$score_layer_nodes[1L],
      stringsAsFactors = FALSE
    )
  }))
} else {
  stability <- data.frame(topology_hash = character(), n_grid_rows = integer(), materialized_solution_ref = character(),
                          tf_panels = character(), pruning_steps = character(), n_layers = character(),
                          threshold_pairs = character(), measurement_features = character(),
                          upstream_source_nodes = character(), observed_tfs = character(),
                          edges = integer(), nodes = integer())
}
moon_write_csv(stability, file.path(out, "score_layer_stability.csv"))
summary <- list(
  analysis = search_cfg$analysis, fits_requested = nrow(fit_grid),
  fits_converged = sum(vapply(fit_rows, function(x) x$converged, logical(1))),
  valid_threshold_pairs_per_fit = sum(thresholds$valid_pair),
  strict_extractions_evaluated = nrow(fit_grid) * sum(thresholds$valid_pair),
  extraction_status = as.list(table(grid$status)),
  eligible_grid_rows = nrow(eligible),
  unique_materialized_topologies = nrow(stability),
  interpretation = "The grid describes score-layer and coverage sensitivity within this exploratory PKN. It does not select a global optimum or establish causal regulation."
)
relaxed_write_json(summary, file.path(out, "search_summary.json"))
writeLines(c(
  "# Relaxed Factor 1 lipid-to-TF parameter search",
  "",
  "This ancillary workflow uses the same score-first candidate semantics and strict consecutive score-layer extraction as the selected LPS refit.",
  "Candidate measurements are copied only as post-score sign-consistency evidence; they are not numerical ULM inputs.",
  "",
  "`fit_results.csv` records every requested compressed fit. Each successful `runs/<run_id>/` directory contains its fit state, compression map, coherence ledger, scores and candidate-status tables.",
  "`grid_results.csv` records every threshold pair, including invalid, failed and empty extractions. Ordinary sign-consistent connectivity is intentionally not generated as a solution output.",
  "`eligible_under_500.csv` requires a nonempty strict solution with fewer than 500 edges, a nonzero secondary cutoff, at least one independent direct-edge fingerprint group and at least two observed TFs.",
  "`score_layer_stability.csv` groups identical signed SIFs across parameter settings. Only the first occurrence of each topology is materialized; later rows point to the same `solution_ref`.",
  "Threshold and coverage frontiers describe tradeoffs only and do not define an optimal or high-confidence network."
), file.path(out, "README.md"))
outputs <- list.files(out, recursive = TRUE, full.names = TRUE)
moon_write_csv(data.frame(path = substring(outputs, nchar(out) + 2L), md5 = unname(tools::md5sum(outputs))),
               file.path(out, "output_manifest.csv"))
cat("Parameter search saved to ", out, ": ", nrow(eligible), " eligible grid rows, ",
    nrow(stability), " unique strict topologies.\n", sep = "")
