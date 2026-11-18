#!/usr/bin/env Rscript
# Separate compressed sensitivity: expanded chain-compatible PKN, frozen or extended observations.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R --mode=expanded|frozen_control --output=NEW_DIR [--mapping=CSV] [--currency-policy=full|no_generic_substrate] [--currency-edges=CSV] [--layer-consistent] [--pruning-steps=6] [--layers=4] [--validate-only]\n")
  quit(save = "no")
}
if (any(!grepl("^--(mode|output|mapping|pruning-steps|layers|currency-policy|currency-edges)=.+$|^--(validate-only|layer-consistent)$", args))) stop("Unknown option")
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
cfg_path <- "scripts/moon/config/factor1_tf_to_lipid_expanded_species.json"
cfg <- jsonlite::read_json(cfg_path, simplifyVector = TRUE)
# Check frozen measurements/software using the baseline checker; load the graph separately.
baseline_cfg <- moon_read_config(cfg$baseline_config)
mode <- option("--mode", "expanded")
stopifnot(mode %in% c("expanded", "frozen_control"))
steps <- as.numeric(option("--pruning-steps", as.character(cfg$pruning_steps)))
layers <- as.numeric(option("--layers", as.character(cfg$n_layers)))
stopifnot(length(steps) == 1L, is.finite(steps), steps >= 1, steps == floor(steps),
          length(layers) == 1L, is.finite(layers), layers >= 1, layers == floor(layers))
for (path in names(cfg$input_md5)) {
  stopifnot(file.exists(path), unname(tools::md5sum(path)) == cfg$input_md5[[path]])
}
data <- moon_load_factor1(baseline_cfg)
data$pkn <- moon_check_pkn(readRDS(cfg$pkn))
full_pkn_edges <- nrow(data$pkn)
currency_policy <- option("--currency-policy", "full")
stopifnot(currency_policy %in% c("full", "no_generic_substrate"))
policy_path <- option("--currency-edges", "support/cosmos_mouse/moon_currency_substrate_edges.csv")
policy_edges_removed <- data.frame(source = character(), target = character(), interaction = numeric())
if (currency_policy == "no_generic_substrate") {
  exclusions <- utils::read.csv(policy_path, check.names = FALSE)
  stopifnot(all(exclusions$interaction == 1), all(grepl("^Enzyme__", exclusions$target)),
    all(moon_edge_keys(exclusions) %in% moon_edge_keys(data$pkn)))
  policy_edges_removed <- data$pkn[moon_edge_keys(data$pkn) %in% moon_edge_keys(exclusions), , drop = FALSE]
  data$pkn <- data$pkn[!moon_edge_keys(data$pkn) %in% moon_edge_keys(exclusions), , drop = FALSE]
  stopifnot(nrow(policy_edges_removed) == nrow(exclusions))
}
frozen <- data$downstream
mapping_path <- option("--mapping", cfg$input_mapping)
if (mode == "expanded") {
  mapping <- utils::read.csv(mapping_path, check.names = FALSE)
  stopifnot(nrow(mapping) == cfg$expected_mapped_observations,
    !anyDuplicated(mapping$node_id), !anyDuplicated(mapping$feature_id),
    all(mapping$node_id %in% moon_nodes(data$pkn)),
    all(grepl("^Metab__", mapping$node_id)),
    !any(grepl("CandidateMeasuredLipid", mapping$node_id, fixed = TRUE)))
  input <- readRDS(baseline_cfg$inputs)
  exact_weights <- unname(input$original_lipid_weight[mapping$feature_id]) * input$orientation$multiplier
  stopifnot(all(is.finite(exact_weights)), max(abs(exact_weights - mapping$weight)) < 1e-13)
  data$downstream <- stats::setNames(exact_weights, mapping$node_id)
  # Preserve the original R doubles, rather than passing through rounded CSV text.
  data$downstream[names(frozen)] <- frozen
  moon_vector_check(data$downstream)
  stopifnot(all(names(frozen) %in% names(data$downstream)),
    identical(unname(frozen), unname(data$downstream[names(frozen)])))
}
stopifnot(!length(intersect(names(data$upstream), names(data$downstream))))
cat("Mode:", mode, "; PKN:", nrow(data$pkn), "edges;", length(moon_nodes(data$pkn)), "nodes\n",
    "Inputs:", length(data$upstream), "TFs and", length(data$downstream), "lipids; pruning", steps,
    "; layers", layers, "; orientation multiplier", data$orientation$multiplier,
    "; currency policy", currency_policy, "; removed", nrow(policy_edges_removed), "edges\n")
if ("--validate-only" %in% args) {
  cat("PASS: packages, source hashes, unique inputs, preserved frozen values and IDs. No output written.\n")
  quit(save = "no")
}
out <- moon_new_output(option("--output", file.path(cfg$output, mode)))
file.copy(script, file.path(out, "runner_source_snapshot.R"))
moon_write_csv(policy_edges_removed, file.path(out, "currency_policy_removed_edges.csv"))
write_json <- function(x, name) jsonlite::write_json(x, file.path(out, name),
  pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
write_json(list(config = cfg, mode = mode, pruning_steps = steps, n_layers = layers,
  orientation = data$orientation, command_args = args, currency_policy = currency_policy,
  currency_edges_path = if (currency_policy != "full") policy_path else NULL,
  currency_policy_removed_edges = nrow(policy_edges_removed), layer_consistent_extraction = "--layer-consistent" %in% args,
  input_policy = if (mode == "expanded") "Preserve frozen155 identities and values; add268 one-to-one chain-compatible measured product features; exclude glycan alternatives" else "Frozen155 observations verbatim; new graph only",
  compression = "Identical signed children; protect all input and TF-regulon IDs",
  terminal_observation_policy = "Remove all outgoing edges from observed lipids before preparation/scoring",
  RNA_policy = "Same frozen RNA TF-target sign coherence; no new expression filter",
  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), "run_manifest.json")
saveRDS(data, file.path(out, "analysis_inputs.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
fit <- moon_fit(data, steps, layers, cfg$max_iterations)
moon_save_fit(fit, out)
moon_write_csv(fit$pkn_final, file.path(out, "pkn_final_uncompressed.csv"))
moon_write_csv(data.frame(source = names(fit$prepared$upstream), score = unname(fit$prepared$upstream)),
  file.path(out, "upstream_input.csv"))
moon_write_csv(data.frame(source = names(fit$prepared$downstream), score = unname(fit$prepared$downstream)),
  file.path(out, "downstream_input.csv"))
lost <- data.frame(source = names(data$downstream), weight = unname(data$downstream),
  retained_after_preparation = names(data$downstream) %in% names(fit$prepared$downstream))
moon_write_csv(lost, file.path(out, "lipid_input_retention.csv"))
fixed <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
  cfg$primary, cfg$secondary, data$rna)
moon_validate_solution(fixed, fit, cfg$primary, cfg$secondary)
dir.create(file.path(out, "fixed_thresholds"))
moon_write_csv(fixed$SIF, file.path(out, "fixed_thresholds/network_SIF.csv"))
moon_write_csv(fixed$ATT, file.path(out, "fixed_thresholds/network_ATT.csv"))
jsonlite::write_json(fixed$audit, file.path(out, "fixed_thresholds/extraction_audit.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
extraction_pkn <- fit$pkn_final
if ("--layer-consistent" %in% args) {
  level <- stats::setNames(fit$scores$level, fit$scores$source)
  keep <- level[extraction_pkn$source] == level[extraction_pkn$target] + 1L
  extraction_pkn <- extraction_pkn[!is.na(keep) & keep, , drop = FALSE]
}
moon_write_csv(extraction_pkn, file.path(out, "pkn_extraction.csv"))
grid <- expand.grid(primary = cfg$primary_grid, secondary = cfg$secondary_grid, KEEP.OUT.ATTRS = FALSE)
grid <- grid[grid$primary >= grid$secondary, , drop = FALSE]
rows <- solutions <- list()
for (i in seq_len(nrow(grid))) {
  p <- grid$primary[i]; s <- grid$secondary[i]
  ans <- moon_extract_paths(fit$scores, extraction_pkn, fit$prepared$upstream, p, s, data$rna)
  moon_validate_solution(ans, fit, p, s)
  rows[[i]] <- data.frame(primary = p, secondary = s, edges = nrow(ans$SIF), nodes = nrow(ans$ATT),
    upstream_tfs = ans$audit$retained_upstream, observed_lipids = ans$audit$retained_observed,
    validated = TRUE, eligible = nrow(ans$SIF) > 0L && nrow(ans$SIF) < 500L &&
      ans$audit$retained_upstream >= 2L && ans$audit$retained_observed >= 2L)
  solutions[[i]] <- ans
}
grid_result <- do.call(rbind, rows)
moon_write_csv(grid_result, file.path(out, "threshold_grid_results.csv"))
eligible <- which(grid_result$eligible)
if (length(eligible)) {
  # Explicit representative policy: higher secondary then primary, followed by coverage.
  ranked <- eligible[order(-grid_result$secondary[eligible], -grid_result$primary[eligible],
    -grid_result$observed_lipids[eligible], -grid_result$upstream_tfs[eligible], grid_result$edges[eligible])]
  chosen <- ranked[1L]
  solution <- solutions[[chosen]]
  primary <- grid_result$primary[chosen]; secondary <- grid_result$secondary[chosen]
  selection <- "Grid representative: <500 edges, >=2TFs and >=2lipids; descending secondary, primary, lipid/TF coverage"
} else {
  solution <- moon_extract_paths(fit$scores, extraction_pkn, fit$prepared$upstream,
    cfg$primary, cfg$secondary, data$rna)
  moon_validate_solution(solution, fit, cfg$primary, cfg$secondary)
  primary <- cfg$primary; secondary <- cfg$secondary
  selection <- "No grid solution met compact coverage criteria; extract original cutoffs under the declared layer policy without claiming compact selection"
}
all_scores <- moon_diagnose(fit, data, solution)
focus <- data.frame(source = cfg$focus_tfs,
  in_initial_upstream = cfg$focus_tfs %in% names(data$upstream),
  in_pruned_pkn = cfg$focus_tfs %in% moon_nodes(fit$prepared$net),
  unchanged_by_compression = fit$compression$node_map$compressed[
    match(cfg$focus_tfs, fit$compression$node_map$original)] == cfg$focus_tfs,
  in_final_scoring_pkn = cfg$focus_tfs %in% moon_nodes(fit$pkn_scoring))
focus <- merge(focus, all_scores, by = "source", all.x = TRUE, sort = FALSE)
moon_write_csv(focus, file.path(out, "focus_tf_retention.csv"))
moon_write_csv(all_scores, file.path(out, "moon_scores_all.csv"))
moon_write_csv(solution$SIF, file.path(out, "network_SIF.csv"))
moon_write_csv(solution$ATT, file.path(out, "network_ATT.csv"))
saveRDS(solution, file.path(out, "solution.rds"))
write_json(c(solution$audit, list(selection = selection, pruning_steps = steps,
  n_layers = layers, compression = "safe_same_children")), "extraction_audit.json")
summary <- list(mode = mode, PKN_edges = nrow(data$pkn), PKN_nodes = length(moon_nodes(data$pkn)),
  full_resource_edges = full_pkn_edges, currency_policy = currency_policy,
  currency_edges_path = if (currency_policy != "full") policy_path else NULL,
  currency_policy_removed_edges = nrow(policy_edges_removed), layer_consistent_extraction = "--layer-consistent" %in% args,
  initial_TFs = length(data$upstream), initial_lipids = length(data$downstream),
  retained_TFs = length(fit$prepared$upstream), retained_lipids = length(fit$prepared$downstream),
  observed_outgoing_edges_removed = sum(data$pkn$source %in% names(data$downstream)),
  original_prepared_edges = nrow(fit$prepared$net), compressed_prepared_edges = nrow(fit$compression$pkn_scoring),
  final_scoring_edges = nrow(fit$pkn_scoring), final_original_edges = nrow(fit$pkn_final),
  compressed_returned_scores = nrow(fit$scores_compressed), decompressed_returned_scores = nrow(fit$scores),
  coherence_iterations = nrow(fit$coherence_ledger), converged = fit$converged,
  pruning_steps = steps, n_layers = layers,
  fixed_thresholds = fixed$audit, selected_primary = primary, selected_secondary = secondary,
  selected_edges = nrow(solution$SIF), selected_nodes = nrow(solution$ATT),
  selected_TFs = solution$ATT$source[solution$ATT$type == "upstream_input"],
  selected_lipids = solution$ATT$source[solution$ATT$level == 0],
  grid_combinations = nrow(grid_result), eligible_combinations = length(eligible), selection = selection,
  score_interpretation = "Prior-guided consistency, not calibrated significance, flux or causal validation")
write_json(summary, "results_summary.json")
dependencies <- unique(c(cfg_path, cfg$baseline_config, names(cfg$input_md5),
  names(baseline_cfg$input_md5), if (mode == "expanded") mapping_path,
  if (currency_policy != "full") policy_path,
  "scripts/moon/lib/workflow.R", "scripts/moon/moon_safe_helpers.R", script))
moon_write_csv(data.frame(path = dependencies, md5 = unname(tools::md5sum(dependencies))),
  file.path(out, "dependency_manifest.csv"))
write_json(list(completed = TRUE, mode = mode, config = cfg, orientation = data$orientation,
  command_args = args, selected_primary = primary, selected_secondary = secondary,
  currency_policy = currency_policy, layer_consistent_extraction = "--layer-consistent" %in% args,
  currency_edges_path = if (currency_policy != "full") policy_path else NULL,
  pruning_steps = steps, n_layers = layers, packages = baseline_cfg$packages,
  cosmos_commit = baseline_cfg$cosmos_commit, finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "completed_manifest.json")
cat("PASS:", mode, ";", nrow(solution$SIF), "edges;", nrow(solution$ATT), "nodes;",
  solution$audit$retained_upstream, "TFs;", solution$audit$retained_observed, "lipids; cutoffs",
  primary, secondary, ". Results:", out, "\n")
