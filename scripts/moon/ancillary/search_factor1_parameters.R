#!/usr/bin/env Rscript
# Compressed-only sensitivity runner. Run from the repository root.
# Frozen input semantics are shared with the selected TF-to-lipid workflow.

search_cli <- function(args) {
  value_keys <- c("config", "output", "pruning-steps", "layers", "primary", "secondary")
  flag_keys <- c("help", "validate-only")
  result <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) stop("Expected --option=value, received: ", arg)
    name <- sub("=.*$", "", sub("^--", "", arg))
    if (!name %in% c(value_keys, flag_keys)) stop("Unknown option: --", name)
    if (name %in% names(result)) stop("Duplicate option: --", name)
    if (name %in% flag_keys) {
      if (grepl("=", arg, fixed = TRUE)) stop("Flag takes no value: --", name)
      result[[name]] <- TRUE
    } else {
      if (!grepl("=", arg, fixed = TRUE)) stop("Use --", name, "=VALUE")
      value <- sub("^[^=]*=", "", arg)
      if (!nzchar(value)) stop("Empty value for --", name)
      result[[name]] <- value
    }
  }
  result
}

search_numbers <- function(value, name, integer = FALSE, strictly_positive = FALSE) {
  if (grepl("(^,|,$|,,)", value)) stop("Empty numeric entry in --", name)
  parts <- strsplit(value, ",", fixed = TRUE)[[1L]]
  numbers <- suppressWarnings(as.numeric(parts))
  if (!length(numbers) || any(!is.finite(numbers)) ||
      any(numbers < 0) || (strictly_positive && any(numbers <= 0)) ||
      (integer && any(numbers != floor(numbers) | numbers > .Machine$integer.max))) {
    stop("Invalid --", name, ": expected comma-separated ",
         if (integer) "positive integers" else if (strictly_positive) "positive finite numbers" else "nonnegative finite numbers")
  }
  sort(unique(if (integer) as.integer(numbers) else numbers))
}

search_write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = 16)
}

search_frontier <- function(x, objectives) {
  if (!nrow(x)) return(logical())
  values <- as.matrix(x[, objectives, drop = FALSE])
  !vapply(seq_len(nrow(x)), function(i) {
    at_least <- rowSums(sweep(values, 2L, values[i, ], `>=`)) == ncol(values)
    strictly <- rowSums(sweep(values, 2L, values[i, ], `>`)) > 0L
    any(at_least & strictly)
  }, logical(1))
}

main <- function() {
  options(stringsAsFactors = FALSE)
  cli <- search_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) {
    cat(paste(
      "Usage: Rscript scripts/moon/ancillary/search_factor1_parameters.R [OPTIONS]",
      "Run from the repository root. Each pruning/layer pair is scored once with compression.",
      "  --config=PATH       Default: scripts/moon/config/factor1_tf_to_lipid.json",
      "  --output=PATH       New/empty directory; default is a timestamped causal_modelling directory",
      "  --pruning-steps=6,12  Default: 6; forward/backward pruning steps",
      "  --layers=4,8        Default: 4,5,6,8; MOON scoring layers",
      "  --primary=1.85,2.15 Default: 0.5 through 3 by 0.1",
      "  --secondary=0.55,1.075 Default: 0 through 2 by 0.1",
      "  --validate-only     Validate options, frozen inputs and packages; do not score or write",
      "  --help              Show this message",
      "Primary must be at least secondary; equal cutoffs are allowed. Invalid pairs are recorded but not extracted.",
      "All grid rows and fit states are saved; eligible solutions (<500 edges, >=2 TFs and >=2 lipids)",
      "are materialized. No cache is reused and no universal optimum is selected.", sep = "\n"), "\n")
    return(invisible(NULL))
  }
  get_value <- function(name, default) if (is.null(cli[[name]])) default else cli[[name]]
  pruning <- search_numbers(get_value("pruning-steps", "6"), "pruning-steps", TRUE, TRUE)
  layers <- search_numbers(get_value("layers", "4,5,6,8"), "layers", TRUE, TRUE)
  primary <- search_numbers(get_value("primary", paste(seq(0.5, 3, 0.1), collapse = ",")), "primary", FALSE, TRUE)
  secondary <- search_numbers(get_value("secondary", paste(seq(0, 2, 0.1), collapse = ",")), "secondary")
  thresholds <- expand.grid(primary = primary, secondary = secondary, KEEP.OUT.ATTRS = FALSE)
  thresholds$threshold_id <- sprintf("T%04d", seq_len(nrow(thresholds)))
  thresholds$valid_pair <- thresholds$primary >= thresholds$secondary
  if (!any(thresholds$valid_pair)) stop("No valid threshold pairs: primary must be at least secondary")
  config_path <- get_value("config", "scripts/moon/config/factor1_tf_to_lipid.json")
  output_path <- get_value("output", file.path("results/causal_modelling",
    paste0("factor1_tf_to_lipid_parameter_search_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())))
  source("scripts/moon/moon_safe_helpers.R")
  source("scripts/moon/lib/workflow.R")
  cfg <- moon_read_config(config_path)
  data <- moon_load_factor1(cfg)
  fit_grid <- expand.grid(pruning_steps = pruning, n_layers = layers, KEEP.OUT.ATTRS = FALSE)
  fit_grid$run_id <- paste0("pruning", fit_grid$pruning_steps, "_layers", fit_grid$n_layers, "_compressed")
  if (dir.exists(output_path) && length(list.files(output_path, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; choose a new --output= directory: ", output_path)
  }
  if (isTRUE(cli[["validate-only"]])) {
    cat("Validated", nrow(fit_grid), "compressed fits and", sum(thresholds$valid_pair),
        "valid threshold pairs per fit;", sum(!thresholds$valid_pair), "invalid pairs will be recorded.\n")
    return(invisible(NULL))
  }
  output <- moon_new_output(output_path)
  dependency_paths <- unique(c(config_path, "scripts/moon/ancillary/search_factor1_parameters.R",
    "scripts/moon/lib/workflow.R", "scripts/moon/moon_safe_helpers.R", names(cfg$input_md5)))
  hashes <- tools::md5sum(dependency_paths)
  moon_write_csv(data.frame(path = names(hashes), md5 = unname(hashes)), file.path(output, "dependency_manifest.csv"))
  file.copy(config_path, file.path(output, "analysis_config.json"), overwrite = FALSE)
  packages <- stats::setNames(vapply(names(cfg$packages), function(pkg)
    as.character(utils::packageVersion(pkg)), character(1)), names(cfg$packages))
  search_write_json(list(analysis = cfg$analysis, compression = "safe_same_children",
    compression_protection = "upstream and downstream inputs plus TF-regulon sources and targets",
    pruning_steps = pruning, n_layers = layers, primary = primary, secondary = secondary,
    max_iterations = cfg$max_iterations, tf_cutoff = cfg$tf_cutoff,
    threshold_policy = "primary >= secondary >= 0, with primary > 0; absolute inferred-score cutoffs remain strict >",
    eligible_policy = "validated, converged, 0 < edges < 500, at least 2 upstream TFs and 2 measured downstream lipids",
    cache_policy = "No cache reuse; every topology/layer fit is recomputed once in this new output directory",
    run_reference_policy = "run_ref and solution_ref are relative to this output directory",
    packages = packages, cosmos_commit = utils::packageDescription("cosmosR")$RemoteSha,
    command_args = commandArgs(trailingOnly = TRUE), started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(output, "search_config.json"))
  utils::write.table(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  moon_write_csv(thresholds, file.path(output, "threshold_grid.csv"))
  moon_write_csv(fit_grid, file.path(output, "fit_grid.csv"))
  rows <- runs <- list()
  persist <- function() {
    if (length(rows)) moon_write_csv(do.call(rbind, rows), file.path(output, "grid_results.csv"))
    if (length(runs)) moon_write_csv(do.call(rbind, runs), file.path(output, "fit_results.csv"))
  }
  for (i in seq_len(nrow(fit_grid))) {
    current <- fit_grid[i, ]
    run_ref <- file.path("runs", current$run_id)
    run_dir <- file.path(output, run_ref)
    dir.create(run_dir, recursive = TRUE)
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    message("Starting compressed fit: ", current$run_id)
    fit <- tryCatch(withCallingHandlers({
      value <- moon_fit(data, current$pruning_steps, current$n_layers, cfg$max_iterations)
      stopifnot(isTRUE(value$converged), nrow(value$coherence_ledger) >= 1L,
                utils::tail(value$coherence_ledger$edges_removed, 1L) == 0L)
      moon_save_fit(value, run_dir)
      value
    }, warning = function(w) warnings <<- unique(c(warnings, conditionMessage(w)))), error = identity)
    failed <- inherits(fit, "error")
    fit_error <- if (failed) conditionMessage(fit) else ""
    fit_status <- if (failed && grepl("did not converge", fit_error, fixed = TRUE)) "nonconvergent" else if (failed) "failed" else "converged"
    runs[[length(runs) + 1L]] <- data.frame(run_id = current$run_id, run_ref = run_ref,
      pruning_steps = current$pruning_steps, n_layers = current$n_layers,
      status = fit_status, converged = !failed, validated = !failed,
      iterations = if (failed) NA_integer_ else nrow(fit$coherence_ledger),
      final_scoring_edges = if (failed) NA_integer_ else nrow(fit$pkn_scoring),
      final_uncompressed_edges = if (failed) NA_integer_ else nrow(fit$pkn_final),
      scored_nodes = if (failed) NA_integer_ else nrow(fit$scores),
      realized_max_level = if (failed) NA_real_ else max(fit$scores$level),
      elapsed_seconds = proc.time()[["elapsed"]] - started, error = fit_error,
      warnings = paste(warnings, collapse = " | "))
    search_write_json(as.list(runs[[length(runs)]]), file.path(run_dir, "fit_audit.json"))
    for (j in seq_len(nrow(thresholds))) {
      cutoffs <- thresholds[j, ]
      row <- data.frame(run_id = current$run_id, run_ref = run_ref,
        pruning_steps = current$pruning_steps, n_layers = current$n_layers,
        threshold_id = cutoffs$threshold_id, primary = cutoffs$primary, secondary = cutoffs$secondary,
        status = "pending", converged = !failed, validated = FALSE,
        nodes = NA_integer_, edges = NA_integer_, upstream_tfs = NA_integer_, downstream_lipids = NA_integer_,
        candidate_edges = NA_integer_, consistent_edges = NA_integer_, eligible = FALSE,
        solution_ref = "", error = "")
      if (!cutoffs$valid_pair) {
        row$status <- "invalid_threshold_pair"
        row$error <- "primary must be at least secondary"
      } else if (failed) {
        row$status <- paste0("fit_", fit_status)
        row$error <- fit_error
      } else {
        solution <- tryCatch({
          value <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
            primary_thresh = cutoffs$primary, secondary_thresh = cutoffs$secondary, RNA_input = data$rna)
          moon_validate_solution(value, fit, cutoffs$primary, cutoffs$secondary)
          value
        }, error = identity)
        if (inherits(solution, "error")) {
          row$status <- "extraction_failed"
          row$error <- conditionMessage(solution)
        } else {
          row$nodes <- nrow(solution$ATT)
          row$edges <- nrow(solution$SIF)
          row$upstream_tfs <- solution$audit$retained_upstream
          row$downstream_lipids <- solution$audit$retained_observed
          row$candidate_edges <- solution$audit$candidate_edges
          row$consistent_edges <- solution$audit$consistent_edges
          row$validated <- TRUE
          row$status <- if (row$edges == 0L) "empty" else "ok"
          row$eligible <- row$edges > 0L && row$edges < 500L && row$upstream_tfs >= 2L && row$downstream_lipids >= 2L
          if (row$eligible) {
            row$solution_ref <- file.path(run_ref, "solutions", cutoffs$threshold_id)
            solution_dir <- file.path(output, row$solution_ref)
            dir.create(solution_dir, recursive = TRUE)
            solution$ATT$upstream_activity <- unname(fit$prepared$upstream[solution$ATT$source])
            saveRDS(solution, file.path(solution_dir, "solution.rds"))
            moon_write_csv(solution$SIF, file.path(solution_dir, "network_SIF.csv"))
            moon_write_csv(solution$ATT, file.path(solution_dir, "network_ATT.csv"))
            search_write_json(c(list(run_id = current$run_id, pruning_steps = current$pruning_steps,
              n_layers = current$n_layers, threshold_id = cutoffs$threshold_id), solution$audit),
              file.path(solution_dir, "extraction_audit.json"))
          }
        }
      }
      rows[[length(rows) + 1L]] <- row
      if (j %% 50L == 0L) {
        persist()
        message(current$run_id, ": extracted/recorded ", j, "/", nrow(thresholds), " threshold pairs")
      }
    }
    persist()
    rm(fit)
    invisible(gc())
  }
  grid <- do.call(rbind, rows)
  eligible <- grid[grid$eligible, , drop = FALSE]
  eligible$coverage_threshold_frontier <- search_frontier(eligible,
    c("primary", "secondary", "upstream_tfs", "downstream_lipids"))
  eligible$threshold_frontier <- search_frontier(eligible, c("primary", "secondary"))
  eligible <- eligible[order(-eligible$secondary, -eligible$primary,
    -eligible$downstream_lipids, -eligible$upstream_tfs, eligible$edges), , drop = FALSE]
  moon_write_csv(eligible, file.path(output, "eligible_under_500.csv"))
  moon_write_csv(eligible[eligible$coverage_threshold_frontier, , drop = FALSE], file.path(output, "coverage_threshold_frontier.csv"))
  moon_write_csv(eligible[eligible$threshold_frontier, , drop = FALSE], file.path(output, "threshold_frontier.csv"))
  summary <- list(fits_requested = nrow(fit_grid), fits_converged = sum(vapply(runs, function(x) x$converged, logical(1))),
    grid_rows = nrow(grid), extraction_status = as.list(table(grid$status)), eligible_solutions = nrow(eligible),
    threshold_frontier_solutions = sum(eligible$threshold_frontier),
    interpretation = "Frontiers describe threshold/coverage tradeoffs only within this tested grid. Scores and thresholds are not calibrated confidence measures across different graph backgrounds. No universal optimum or new manuscript selection is asserted.",
    empty_solution_policy = "Empty, validated extractions remain in grid_results.csv; they are not eligible.",
    failed_run_policy = "Failed and nonconvergent fits are recorded for every requested threshold pair; no solution is extracted from them.",
    finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  search_write_json(summary, file.path(output, "search_summary.json"))
  writeLines(c("# Compressed Factor 1 TF-to-lipid parameter search", "",
    "This ancillary run uses the frozen input mapping and scoring semantics of the main workflow.",
    "`analysis_config.json`, `search_config.json`, `dependency_manifest.csv` and `session_info.txt` record inputs and software.",
    "`fit_results.csv` records convergence/errors. Each `run_ref` contains its validated scoring state, compression map and coherence ledger when fitting succeeded.",
    "`grid_results.csv` retains all requested threshold combinations, including invalid pairs, failed fits, failed extractions and empty networks.",
    "`eligible_under_500.csv` requires fewer than 500 edges, at least two TFs and two measured lipid nodes, convergence and validation.",
    "Eligible networks are saved at `solution_ref`; other networks can be re-extracted from their run state and recorded thresholds.",
    "Both frontier tables describe threshold/coverage tradeoffs within the requested grid. Scores and thresholds are not calibrated confidence measures across different graph backgrounds. The frontiers do not establish a universal optimum or replace the selected manuscript analysis.",
    "All run and solution references are relative to this directory. No existing archive or fit cache was used."), file.path(output, "README.md"))
  cat("Parameter search saved to", output, "with", nrow(eligible), "eligible solutions.\n")
  if (!any(grid$validated)) stop("No validated extractions; inspect the saved fit_results.csv and grid_results.csv")
  invisible(summary)
}

main()
