#!/usr/bin/env Rscript
# Diagnose saved final MOON scoring graphs. Does not refit coherence or extract networks.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/audit_lipid_tf_scores.R [--analysis=ANALYSIS_DIRECTORY] [--output=NEW_DIRECTORY]\n",
      "Recomputes unconstrained scores on saved final scoring graphs only; no refit or extraction.\n")
  quit(save = "no")
}
stopifnot(all(grepl("^--(analysis|output)=.+$", args)))
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
analysis_dir <- normalizePath(option("--analysis", "results/causal_modelling/factor1_lipid_to_tf"))
manifest_path <- file.path(analysis_dir, "run_manifest.json")
manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
cfg <- manifest$config
stopifnot(identical(cfg$analysis, "factor1_lipid_to_tf"), isTRUE(manifest$completed))
out <- moon_new_output(option("--output", file.path(analysis_dir, "score_audit")))
input <- readRDS(cfg$inputs)
moon_vector_check(input$lipid_weight)
moon_vector_check(input$tf_activity)
mapping_path <- file.path(cfg$input_directory, "factor1_lipid_node_mapping.csv")
mapping <- utils::read.csv(mapping_path, check.names = FALSE)
mapping <- mapping[mapping$included, , drop = FALSE]
stopifnot(!anyDuplicated(mapping$node_id), setequal(mapping$node_id, names(input$lipid_weight)))
catalog_path <- file.path(cfg$pkn_directory, "nodes.csv")
catalog <- utils::read.csv(catalog_path, check.names = FALSE)
node_labels <- stats::setNames(catalog$label, catalog$node_id)
runs_path <- file.path(analysis_dir, "run_summary.csv")
runs <- utils::read.csv(runs_path)
scored <- runs[runs$status == "scored", , drop = FALSE]
stopifnot(nrow(scored) > 0L)
policy_files <- c(lipid_precursor = "pkn_lipid_precursor.csv",
                  ph73_lipid_precursor = "pkn_ph73_lipid_precursor.csv")
policy_nodes <- lapply(policy_files, function(file) {
  path <- file.path(analysis_dir, "coverage", file)
  if (file.exists(path)) moon_nodes(utils::read.csv(path)) else character()
})
tf_focus <- paste0("Gene__", c("Esr1", "Nfic", "Dbp", "Rxra"))
chemical_focus <- paste0("Metab__CHEBI:", c("72998", "73858", "74963", "77016", "32395", "57756", "60119"), "_unknown")
reaction_focus <- catalog$node_id[catalog$reaction_id %in% c("RHEA:37476", "RHEA:41232")]
focus <- unique(c(tf_focus, chemical_focus, reaction_focus))
lipid_rows <- tf_rows <- focal_rows <- support_rows <- summary_rows <- list()
source_paths <- c(script, "scripts/moon/lib/workflow.R", "scripts/moon/moon_safe_helpers.R",
                  manifest_path, cfg$inputs, mapping_path, catalog_path, runs_path)

for (run_index in seq_len(nrow(scored))) {
  run <- scored[run_index, , drop = FALSE]
  message("Diagnosing ", run$run_id)
  state_path <- file.path(analysis_dir, "runs", run$run_id, "run_state.rds")
  source_paths <- c(source_paths, state_path)
  fit <- readRDS(state_path)
  stopifnot(isTRUE(fit$converged), fit$n_layers == run$n_layers,
            fit$pruning_steps == run$pruning_steps)
  # NULL removes the package's final upstream sign filter only. No iteration,
  # input mapping, coherence filtering, or solution extraction is rerun here.
  raw_compressed <- moon_score(fit$pkn_scoring, NULL, fit$prepared$downstream, fit$n_layers)
  raw <- moon_safe_decompress_scores(raw_compressed, fit$compression$node_map)
  si <- match(fit$scores$source, raw$source)
  ci <- match(fit$scores_compressed$source, raw_compressed$source)
  stopifnot(!anyNA(si), !anyNA(ci),
    all(fit$scores$level == raw$level[si]),
    all(fit$scores_compressed$level == raw_compressed$level[ci]),
    max(abs(fit$scores$score - raw$score[si])) < 1e-10,
    max(abs(fit$scores_compressed$score - raw_compressed$score[ci])) < 1e-10)
  excluded <- setdiff(raw$source, fit$scores$source)
  stopifnot(all(excluded %in% names(fit$prepared$upstream)))
  if (length(excluded)) stopifnot(all(sign(raw$score[match(excluded, raw$source)]) !=
    sign(fit$prepared$upstream[excluded])))
  raw$returned_with_upstream_sign_filter <- raw$source %in% fit$scores$source
  raw$measured_lipid_weight <- unname(input$lipid_weight[raw$source])
  raw$measured_tf_activity <- unname(input$tf_activity[raw$source])
  moon_write_csv(raw, file.path(out, paste0(run$run_id, "_unfiltered_scores.csv")))
  node_map <- stats::setNames(fit$compression$node_map$compressed, fit$compression$node_map$original)
  score_map <- stats::setNames(raw_compressed$score, raw_compressed$source)
  level_map <- stats::setNames(raw_compressed$level, raw_compressed$source)
  original_final_nodes <- moon_nodes(fit$pkn_final)
  final_graph <- igraph::graph_from_data_frame(fit$pkn_scoring[, c("source", "target")], directed = TRUE)
  final_ids <- igraph::V(final_graph)$name
  target_ids <- intersect(names(fit$prepared$downstream), final_ids)
  query_ids <- unique(na.omit(unname(node_map[unique(c(mapping$node_id, focus))])))
  query_ids <- intersect(query_ids, final_ids)
  distances <- if (length(query_ids) && length(target_ids))
    igraph::distances(final_graph, v = query_ids, to = target_ids, mode = "out") else
      matrix(numeric(), nrow = 0L, ncol = 0L)
  shortest <- if (nrow(distances)) apply(distances, 1L, min) else numeric()

  describe <- function(node) {
    mapped <- unname(node_map[node])
    if (is.na(mapped)) mapped <- node
    index <- match(node, raw$source)
    has_score <- !is.na(index)
    level <- if (has_score) raw$level[index] else NA_integer_
    score <- if (has_score) raw$score[index] else NA_real_
    retained <- node %in% fit$scores$source
    in_prepared <- node %in% moon_nodes(fit$prepared$net)
    in_final <- node %in% original_final_nodes
    distance <- unname(shortest[mapped])
    if (!length(distance)) distance <- NA_real_
    status <- if (has_score && level == 0L) "observed_level0" else
      if (has_score && retained) "inferred_returned" else
      if (has_score && !retained) "inferred_upstream_sign_conflict" else
      if (!in_prepared) "not_retained_after_topology_pruning" else
      if (!in_final) "lost_from_final_coherence_filtered_graph" else
      if (is.finite(distance) && distance > fit$n_layers) "shortest_endpoint_distance_exceeds_requested_layers" else
        "not_scored_despite_topology_support_within_requested_layers"
    previous_ids <- if (has_score && level > 0L) names(level_map)[level_map == level - 1L] else character()
    children <- fit$pkn_scoring[fit$pkn_scoring$source == mapped &
      fit$pkn_scoring$target %in% previous_ids, , drop = FALSE]
    children$target_score <- unname(score_map[children$target])
    children$target_level <- unname(level_map[children$target])
    children$signed_target_score <- children$interaction * children$target_score
    children$agrees_with_parent_score <- if (has_score) sign(children$signed_target_score) == sign(score) else logical(nrow(children))
    list(node = node, mapped = mapped, score = score, level = level, retained = retained,
      in_prepared = in_prepared, in_final = in_final, distance = distance, status = status,
      previous_universe = length(previous_ids), children = children)
  }
  for (j in seq_len(nrow(mapping))) {
    m <- mapping[j, , drop = FALSE]
    d <- describe(m$node_id)
    in_policy <- if (run$policy %in% names(policy_nodes)) m$node_id %in% policy_nodes[[run$policy]] else NA
    status <- if (identical(in_policy, FALSE)) "absent_from_policy_graph" else d$status
    lipid_rows[[length(lipid_rows) + 1L]] <- data.frame(run_id = run$run_id,
      policy = run$policy, tf_panel = run$tf_panel, pruning_steps = fit$pruning_steps,
      n_layers = fit$n_layers, feature = m$feature, node_id = m$node_id,
      measured_weight = unname(input$lipid_weight[m$node_id]),
      present_in_policy_graph = in_policy, retained_prepared_input = m$node_id %in% names(fit$prepared$upstream),
      present_in_final_scoring_graph = d$in_final, minimum_final_graph_distance_to_observation = d$distance,
      inferred_score_unfiltered = d$score, level = d$level,
      returned_by_moon = d$retained,
      sign_agrees = if (is.finite(d$score)) sign(d$score) == sign(input$lipid_weight[m$node_id]) else NA,
      status = status, previous_layer_compressed_universe = d$previous_universe,
      signed_children_at_previous_layer = nrow(d$children),
      supporting_children = sum(d$children$agrees_with_parent_score),
      opposing_children = sum(!d$children$agrees_with_parent_score),
      interpretation = "MOON score is prior-guided signed alignment; no statistical confidence or direct ligand activity is inferred.")
  }
  for (node in focus) {
    d <- describe(node)
    row <- data.frame(run_id = run$run_id, policy = run$policy, tf_panel = run$tf_panel,
      pruning_steps = fit$pruning_steps, n_layers = fit$n_layers,
      node_id = node, label = unname(node_labels[node]), scoring_node = d$mapped,
      measured_tf_activity = unname(input$tf_activity[node]),
      measured_lipid_weight = unname(input$lipid_weight[node]),
      selected_as_downstream_observation = node %in% names(fit$prepared$downstream),
      present_in_prepared_graph = d$in_prepared, present_in_final_graph = d$in_final,
      score = d$score, level = d$level, returned_by_moon = d$retained, status = d$status,
      minimum_final_graph_distance_to_observation = d$distance,
      previous_layer_compressed_universe = d$previous_universe,
      signed_children_at_previous_layer = nrow(d$children),
      supporting_children = sum(d$children$agrees_with_parent_score),
      opposing_children = sum(!d$children$agrees_with_parent_score))
    focal_rows[[length(focal_rows) + 1L]] <- row
    if (node %in% tf_focus) tf_rows[[length(tf_rows) + 1L]] <- row
    if (nrow(d$children)) {
      ch <- d$children
      ch$run_id <- run$run_id
      ch$original_parent <- node
      ch$parent_score <- d$score
      ch$parent_level <- d$level
      ch$previous_layer_compressed_universe <- d$previous_universe
      support_rows[[length(support_rows) + 1L]] <- ch
    }
  }
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(run_id = run$run_id,
    tf_panel = run$tf_panel, pruning_steps = fit$pruning_steps, n_layers = fit$n_layers,
    final_compressed_edges = nrow(fit$pkn_scoring),
    final_original_edges = nrow(fit$pkn_final), prepared_lipid_inputs = length(fit$prepared$upstream),
    prepared_tf_observations = length(fit$prepared$downstream),
    inferred_and_observed_compressed_nodes = nrow(raw_compressed),
    inferred_and_observed_original_nodes = nrow(raw),
    input_lipids_with_any_inferred_score = sum(names(input$lipid_weight) %in% raw$source),
    input_lipids_returned = sum(names(input$lipid_weight) %in% fit$scores$source),
    input_lipids_sign_excluded = length(excluded), realized_max_level = max(raw$level),
    max_saved_score_difference = max(abs(fit$scores$score - raw$score[si])),
    max_saved_compressed_score_difference = max(abs(fit$scores_compressed$score - raw_compressed$score[ci])),
    excluded_node_ids = paste(excluded, collapse = "|"), retained_score_comparison_passed = TRUE)
}
lipids <- do.call(rbind, lipid_rows)
focus_tfs <- do.call(rbind, tf_rows)
focals <- do.call(rbind, focal_rows)
summaries <- do.call(rbind, summary_rows)
moon_write_csv(lipids, file.path(out, "lipid_input_audit.csv"))
moon_write_csv(focus_tfs, file.path(out, "focus_tf_audit.csv"))
moon_write_csv(focals, file.path(out, "focal_node_scores.csv"))
moon_write_csv(summaries, file.path(out, "run_score_validation.csv"))
if (length(support_rows)) moon_write_csv(do.call(rbind, support_rows), file.path(out, "scoring_support_edges.csv"))
moon_write_csv(lipids[lipids$node_id %in% c("Metab__CHEBI:72998_unknown", "Metab__CHEBI:73858_unknown"), ],
               file.path(out, "two_lpc_score_summary.csv"))
source_paths <- unique(c(source_paths, file.path(analysis_dir, "coverage", unname(policy_files))))
source_paths <- source_paths[file.exists(source_paths)]
moon_write_csv(data.frame(path = source_paths, md5 = unname(tools::md5sum(source_paths))),
               file.path(out, "dependency_manifest.csv"))
jsonlite::write_json(list(analysis = "factor1_lipid_to_tf_saved_score_diagnostic",
  final_graph_rescoring_only = TRUE, coherence_refitting = FALSE, extraction_rerun = FALSE,
  input_mapping_changed = FALSE, saved_runs = nrow(scored), frozen_input_lipids = nrow(mapping),
  all_saved_scores_reproduced = all(summaries$retained_score_comparison_passed),
  maximum_score_difference = max(summaries$max_saved_score_difference),
  statuses = as.list(table(lipids$status)),
  score_semantics = "Signed ULM enrichment/activity scores, not p-values, probabilities, confidence bounds, reaction flux, or direct ligand activity.",
  support_semantics = "Signed outgoing children in the final compressed PKN at the immediately previous returned MOON level, with that compressed layer as background. This describes structural score support, not independent experiments.",
  pruning_semantics = "Distances use the final compressed scoring graph and do not replace the recorded iterative MOON pruning ledger.",
  normalization_caveat = "All scored runs belong to the optional existing pH7.3-normalized precursor sensitivity; primary graph hypotheses remain absent."),
  file.path(out, "summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines(c("# Saved lipid-to-TF score audit", "",
  "Scores were recomputed only on each saved final compressed scoring graph, with upstream_input=NULL to expose the package's final sign-filter exclusions. Every originally returned score and level was compared with the recomputation. No coherence loop, input mapping, parameter fitting or network extraction was rerun.", "",
  "`lipid_input_audit.csv` diagnoses all 155 frozen accepted lipid inputs in every saved scored run. `two_lpc_score_summary.csv` isolates the two LPCs that support normalized precursor topology. Absent policy nodes, topology loss, final-graph loss, missing layer support, and upstream sign exclusion are distinct statuses.", "",
  "`focus_tf_audit.csv` distinguishes Esr1, Nfic, Dbp and Rxra supplied as level-0 observations from intermediates with inferred scores. An observed activity is not an independently predicted TF activity.", "",
  "`focal_node_scores.csv` and `scoring_support_edges.csv` report the relevant LPC/PC/DHA intermediates, reaction instances, and TFs, using final compressed score-support edges. Counts refer to signed children in the previous scoring layer and its compressed universe. They are not independent evidence counts or statistical confidence intervals.", "",
  "Low or conflicting inferred lipid scores cannot be repaired by lowering extraction cutoffs alone. Precursor paths can depend on unmeasured cosubstrates; the LPC16/PC/DHA acyl-remodeling route regenerates LPC and does not establish net DHA production."),
  file.path(out, "README.md"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Validated", nrow(scored), "saved scoring runs and", nrow(lipids), "lipid/run statuses. Outputs:", out, "\n")
