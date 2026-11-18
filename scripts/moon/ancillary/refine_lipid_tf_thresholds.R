#!/usr/bin/env Rscript
# Targeted dual-threshold refinement from frozen lipid-to-TF fits; never refits MOON.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run using Rscript")
setwd(normalizePath(file.path(dirname(script), "../../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R [OPTIONS]\n",
      "--results=DIR  Frozen main run; default results/causal_modelling/factor1_lipid_to_tf\n",
      "--output=DIR   New/empty refinement output; default RESULTS/refinement\n",
      "--step=0.025   Fine threshold spacing\n",
      "--workers=2    Independent extraction workers (POSIX fork)\n",
      "--select-only  Materialize/export the agreed 1.65/1.10 all_tfs candidate from an existing refinement.\n",
      "Default runs refinement then selects that candidate. Existing selected output is never overwritten.\n",
      "Reproduce: Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R --output=NEW_EMPTY_DIRECTORY\n",
      "Only runs with nonempty coarse solutions are refined. No scores or PKN edges change.\n",
      "The primary range extends from each run's highest nonempty coarse primary to the next coarse primary.\n",
      "The secondary range extends from the highest nonempty coarse secondary to just above the first-edge bottleneck.\n",
      "Outputs retain all extractions and separate all_tfs from strong_tfs observation panels.\n")
  quit(save = "no")
}
stopifnot(all(grepl("^--(results|output|step|workers)=.+$|^--select-only$", args)))
if (sum(args == "--select-only") > 1L) stop("Repeated --select-only")
option <- function(name, default) {
  values <- args[startsWith(args, paste0("--", name, "="))]
  if (length(values) > 1L) stop("Duplicate option: ", name)
  if (length(values)) sub("^[^=]*=", "", values) else default
}
results <- option("results", "results/causal_modelling/factor1_lipid_to_tf")
out_path <- option("output", file.path(results, "refinement"))
step <- as.numeric(option("step", "0.025"))
workers <- as.integer(option("workers", "2"))
stopifnot(is.finite(step), step > 0, step <= 0.25, is.finite(workers), workers >= 1L)
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
materialize_all_tfs <- function(results, output, inputs, primary = 1.65, secondary = 1.10) {
  id <- "ph73_lipid_precursor_all_tfs_pruning12_layers12"
  grid_path <- file.path(output, "runs", id, "threshold_grid.csv")
  grid <- utils::read.csv(grid_path)
  selected_row <- grid[abs(grid$primary - primary) < 1e-12 & abs(grid$secondary - secondary) < 1e-12, , drop = FALSE]
  stopifnot(nrow(selected_row) == 1L, isTRUE(selected_row$eligible), isTRUE(selected_row$validated))
  selected <- moon_new_output(file.path(output, "selected_all_tfs"))
  fit_path <- file.path(results, "runs", id, "run_state.rds")
  fit <- readRDS(fit_path)
  stopifnot(isTRUE(fit$converged))
  solution <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
    primary_thresh = primary, secondary_thresh = secondary, RNA_input = inputs$hepatocyte_gene_weight)
  moon_validate_solution(solution, fit, primary, secondary)
  stopifnot(nrow(solution$SIF) == selected_row$edges, nrow(solution$ATT) == selected_row$nodes,
    solution$audit$retained_upstream == selected_row$upstream_lipids,
    solution$audit$retained_observed == selected_row$downstream_tfs)
  solution$ATT$upstream_lipid_weight <- unname(fit$prepared$upstream[solution$ATT$source])
  solution$ATT$downstream_tf_activity <- ifelse(solution$ATT$level == 0,
    unname(fit$prepared$downstream[solution$ATT$source]), NA_real_)
  solution$ATT$label <- sub("^Gene__", "", solution$ATT$source)
  saveRDS(solution, file.path(selected, "solution.rds"))
  moon_write_csv(solution$SIF, file.path(selected, "network_SIF.csv"))
  moon_write_csv(solution$ATT, file.path(selected, "network_ATT.csv"))
  jsonlite::write_json(solution$audit, file.path(selected, "extraction_audit.json"),
    pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
  neighborhood <- grid[grid$primary >= 1.60 - 1e-12 & grid$primary <= 1.70 + 1e-12 &
    grid$secondary >= 1.05 - 1e-12 & grid$secondary <= 1.175 + 1e-12, , drop = FALSE]
  moon_write_csv(neighborhood, file.path(selected, "threshold_neighborhood.csv"))
  pkn <- "results/causal_modelling/cosmos_mouse_pkn"
  mapping <- "miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv"
  exporter <- "scripts/moon/ancillary/export_lipid_tf_network.py"
  export_args <- c(exporter, "--results", selected, "--pkn", pkn,
    "--input-mapping", mapping, "--prior-csv", file.path(pkn, "pkn_lipid_curated_ph73.csv"),
    "--node-catalog", file.path(pkn, "nodes_ph73.jsonl"), "--node-aliases", file.path(pkn, "chemical_node_aliases.csv"))
  status <- system2("python3", vapply(export_args, shQuote, character(1)),
    stdout = file.path(selected, "export.log"), stderr = file.path(selected, "export.log"))
  if (status != 0L) stop("Exporter failed; inspect ", file.path(selected, "export.log"))
  dependencies <- c(fit_path, grid_path, exporter, file.path(pkn, "pkn_lipid_curated_ph73.csv"),
    file.path(pkn, "nodes_ph73.jsonl"), file.path(pkn, "chemical_node_aliases.csv"), mapping)
  jsonlite::write_json(list(run_id = id, tf_panel = "all_tfs", exploratory_policy = "ph73_lipid_precursor",
    pruning_steps = fit$pruning_steps, n_layers = fit$n_layers, primary = primary, secondary = secondary,
    edges = nrow(solution$SIF), nodes = nrow(solution$ATT), upstream_lipids = solution$audit$retained_upstream,
    downstream_tfs = solution$audit$retained_observed, refitted = FALSE, manually_trimmed = FALSE,
    choice = "Representative exploratory all_tfs network below 500 edges; rounded secondary 1.10 avoids selecting immediately below the upstream score boundary 1.177648. Not a statistical or biological optimum.",
    limitations = "Only the existing acid/base-normalized sensitivity supports this network. Reaction participation encodes metabolic availability hypotheses, not demonstrated flux or direct measured-lipid signaling. All nonzero TF observations remain the primary panel; strong_tfs is separate.",
    reproduction = paste("Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R --results=", results,
      " --output=NEW_EMPTY_DIRECTORY", sep = ""),
    materialize_existing_refinement = paste("Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R --select-only --results=", results,
      " --output=", output, sep = ""),
    export_command = paste("python3", paste(vapply(export_args, shQuote, character(1)), collapse = " ")),
    input_md5 = as.list(tools::md5sum(dependencies)), selection_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(selected, "selection.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
  message("Selected exploratory all_tfs network: ", nrow(solution$SIF), " edges, ", nrow(solution$ATT),
    " nodes; primary ", primary, ", secondary ", secondary)
  invisible(selected)
}

select_only <- "--select-only" %in% args
out <- if (select_only) normalizePath(out_path, mustWork = TRUE) else moon_new_output(out_path)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
coarse <- utils::read.csv(file.path(results, "threshold_grid.csv"))
manifest <- jsonlite::read_json(file.path(results, "run_manifest.json"), simplifyVector = TRUE)
inputs <- readRDS(manifest$config$inputs)
if (select_only) {
  materialize_all_tfs(results, out, inputs)
  quit(save = "no")
}
ids <- unique(coarse$run_id[coarse$edges > 0])
stopifnot(length(ids) > 0L)
run_paths <- file.path(results, "runs", ids, "run_state.rds")
stopifnot(all(file.exists(run_paths)))
dependencies <- c("scripts/moon/ancillary/refine_lipid_tf_thresholds.R", "scripts/moon/moon_safe_helpers.R",
  "scripts/moon/lib/workflow.R", "scripts/moon/ancillary/export_lipid_tf_network.py", file.path(results, "threshold_grid.csv"),
  file.path(results, "run_manifest.json"), manifest$config$inputs, run_paths)
moon_write_csv(data.frame(path = dependencies, md5 = unname(tools::md5sum(dependencies))),
               file.path(out, "dependency_manifest.csv"))
refine_run <- function(id) {
  fit_path <- file.path(results, "runs", id, "run_state.rds")
  fit <- readRDS(fit_path)
  stopifnot(isTRUE(fit$converged))
  original_scores <- fit$scores
  original_prior <- fit$pkn_final
  previous <- coarse[coarse$run_id == id, , drop = FALSE]
  nonempty <- previous[previous$edges > 0, , drop = FALSE]
  sm <- stats::setNames(fit$scores$score, fit$scores$source)
  levels <- stats::setNames(fit$scores$level, fit$scores$source)
  seed <- fit$scores[fit$scores$source %in% names(fit$prepared$upstream), , drop = FALSE]
  first <- fit$pkn_final[fit$pkn_final$source %in% seed$source, , drop = FALSE]
  first$source_score <- unname(sm[first$source])
  first$target_score <- unname(sm[first$target])
  first <- first[is.finite(first$source_score) & is.finite(first$target_score) &
                   sign(first$source_score * first$target_score) == first$interaction, , drop = FALSE]
  stopifnot(nrow(first) > 0L)
  first$secondary_bound <- pmin(abs(first$source_score),
    ifelse(levels[first$target] == 0, Inf, abs(first$target_score)))
  secondary_bound <- max(first$secondary_bound)
  primary_start <- max(nonempty$primary)
  upper_options <- unique(previous$primary[previous$primary > primary_start])
  primary_end <- if (length(upper_options)) min(upper_options) else primary_start + 0.5
  secondary_start <- max(nonempty$secondary)
  secondary_end <- ceiling(secondary_bound / step) * step
  primary <- unique(c(seq(primary_start, primary_end, by = step), primary_end))
  secondary <- unique(c(seq(secondary_start, secondary_end, by = step), secondary_end))
  grid <- expand.grid(primary = primary, secondary = secondary)
  grid <- grid[grid$primary >= grid$secondary, , drop = FALSE]
  run_dir <- file.path(out, "runs", id)
  dir.create(run_dir, recursive = TRUE)
  moon_write_csv(first, file.path(run_dir, "upstream_first_edge_bottlenecks.csv"))
  write_json(list(run_id = id, tf_panel = previous$tf_panel[1],
    primary_start = primary_start, primary_end = primary_end,
    secondary_start = secondary_start, secondary_end = secondary_end,
    upstream_first_edge_secondary_bound = secondary_bound, strict_greater_than = TRUE,
    step = step, threshold_pairs = nrow(grid), fit_path = fit_path), file.path(run_dir, "refinement_bounds.json"))
  rows <- vector("list", nrow(grid))
  for (j in seq_len(nrow(grid))) {
    p <- grid$primary[j]; s <- grid$secondary[j]
    value <- tryCatch({
      value <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
        primary_thresh = p, secondary_thresh = s, RNA_input = inputs$hepatocyte_gene_weight)
      moon_validate_solution(value, fit, p, s)
      value
    }, error = identity)
    failed <- inherits(value, "error")
    rows[[j]] <- data.frame(run_id = id, tf_panel = previous$tf_panel[1],
      policy = previous$policy[1], pruning_steps = fit$pruning_steps, n_layers = fit$n_layers,
      primary = p, secondary = s, validated = !failed,
      edges = if (failed) NA_integer_ else nrow(value$SIF),
      nodes = if (failed) NA_integer_ else nrow(value$ATT),
      upstream_lipids = if (failed) NA_integer_ else value$audit$retained_upstream,
      downstream_tfs = if (failed) NA_integer_ else value$audit$retained_observed,
      eligible = !failed && nrow(value$SIF) > 0L && nrow(value$SIF) < 500L &&
        value$audit$retained_upstream >= 1L && value$audit$retained_observed >= 2L,
      error = if (failed) conditionMessage(value) else "")
  }
  stopifnot(identical(original_scores, fit$scores), identical(original_prior, fit$pkn_final))
  table <- do.call(rbind, rows)
  moon_write_csv(table, file.path(run_dir, "threshold_grid.csv"))
  message("Refined ", id, ": ", nrow(table), " pairs; ", sum(table$eligible), " eligible")
  table
}
tables <- parallel::mclapply(ids, refine_run, mc.cores = min(workers, length(ids)), mc.preschedule = FALSE)
if (any(vapply(tables, inherits, logical(1), "try-error"))) stop("A refinement worker failed; inspect partial run outputs")
grid <- do.call(rbind, tables)
grid <- grid[order(grid$tf_panel, grid$run_id, grid$primary, grid$secondary), , drop = FALSE]
moon_write_csv(grid, file.path(out, "threshold_grid.csv"))
for (panel in unique(grid$tf_panel)) {
  tab <- grid[grid$tf_panel == panel & grid$eligible, , drop = FALSE]
  tab <- tab[order(-tab$secondary, -tab$primary, -tab$downstream_tfs, tab$edges), , drop = FALSE]
  moon_write_csv(tab, file.path(out, paste0("eligible_", panel, ".csv")))
}
write_json(list(analysis = "targeted_lipid_to_tf_threshold_refinement", frozen_runs = length(ids),
  threshold_evaluations = nrow(grid), validated = sum(grid$validated), eligible = sum(grid$eligible),
  eligible_by_panel = as.list(table(grid$tf_panel[grid$eligible])),
  no_refitting = TRUE, no_edge_trimming = TRUE, no_observation_or_mapping_changes = TRUE,
  selection = "No universal optimum asserted. all_tfs is the primary observation panel; strong_tfs is a separately reported sensitivity.",
  scope = "Targeted fine grids around the greatest coarse primary/secondary values yielding nonempty solutions. Not an exhaustive continuous search.",
  threshold_semantics = "Both cutoffs are strict absolute inferred-MOON-score thresholds. Observed TF activities are terminal and exempt. Thresholds are not calibrated confidence probabilities.",
  command_args = args, finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
  file.path(out, "refinement_summary.json"))
capture.output(utils::sessionInfo(), file = file.path(out, "sessionInfo.txt"))
writeLines(c("# Lipid-to-TF threshold refinement", "",
  "This targeted refinement reads the frozen fits listed in dependency_manifest.csv; no MOON scores, observations, mappings or PKN edges are changed.",
  "Only runs with nonempty coarse solutions are explored. Bounds and first-edge bottlenecks are explicit in each runs/<id>/ directory.",
  "All extractions use moon_extract_paths and moon_validate_solution. Eligible means 0 < edges < 500, >=1 measured upstream lipid and >=2 observed downstream TFs.",
  "all_tfs is the primary observation panel. strong_tfs is an observation-set sensitivity and is not pooled with it for biological ranking.",
  "A compact output is descriptive; it does not establish direct lipid signaling, metabolic flux or a causal relationship."), file.path(out, "README.md"))
cat("Refinement complete:", nrow(grid), "threshold evaluations;", sum(grid$eligible), "eligible.\n")

materialize_all_tfs(results, out, inputs)
