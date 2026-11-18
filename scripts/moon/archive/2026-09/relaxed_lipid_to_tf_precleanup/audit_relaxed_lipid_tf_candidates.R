#!/usr/bin/env Rscript
# Audit the relaxed Factor 1 lipid -> TF candidate screen without refitting MOON.
# It collapses only identical direct-edge fingerprints for reporting, while retaining
# every original feature-to-chemical mapping row as provenance.
options(stringsAsFactors = FALSE, digits = 17)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript scripts/moon/ancillary/audit_relaxed_lipid_tf_candidates.R ",
    "[--analysis=PATH] [--output=PATH]\n",
    "Reads a completed relaxed lipid-to-TF MOON run, checks its dependency hashes, ",
    "and writes fingerprint-aware candidate and route audits without refitting.\n"
  ))
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

analysis_dir <- normalizePath(option(
  "--analysis", "results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory"
))
output_dir <- option("--output", file.path(analysis_dir, "candidate_audit"))
if (!dir.exists(analysis_dir)) stop("Analysis directory does not exist: ", analysis_dir)
out <- moon_new_output(output_dir)

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}
write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", null = "null", digits = NA)
}
collapse_text <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  paste(x, collapse = ";")
}
edge_key <- function(x) paste(x$source, x$target, x$interaction, sep = "\r")

# The fit must be audited against the dependencies it declared at scoring time.
deps <- utils::read.csv(file.path(analysis_dir, "dependency_manifest.csv"), check.names = FALSE)
require_columns(deps, c("path", "md5"), "Dependency manifest")
for (i in seq_len(nrow(deps))) {
  path <- deps$path[i]
  if (!file.exists(path)) stop("Declared dependency is absent: ", path)
  actual <- unname(tools::md5sum(path))
  if (!identical(actual, as.character(deps$md5[i]))) {
    stop("Dependency changed since the relaxed MOON run: ", path)
  }
}

manifest <- jsonlite::read_json(file.path(analysis_dir, "run_manifest.json"), simplifyVector = FALSE)
if (is.null(manifest$config$inputs) || is.null(manifest$config$mapping)) {
  stop("Run manifest lacks the input/mapping paths required for this audit")
}
input <- readRDS(manifest$config$inputs)
mapping_source_path <- manifest$config$mapping
mapping_source_dir <- dirname(mapping_source_path)

mapping <- utils::read.csv(file.path(analysis_dir, "relaxed_candidate_mapping.csv"), check.names = FALSE)
require_columns(mapping, c(
  "assignment_id", "measurement_feature", "candidate_measurement_value",
  "pkn_source_node", "pkn_source_name", "edge_fingerprint_group",
  "assignment_role", "primary_mapping_policy"
), "Relaxed mapping copy")
if (anyDuplicated(mapping$assignment_id)) stop("Mapping has duplicate assignment IDs")

grid <- utils::read.csv(file.path(analysis_dir, "threshold_grid.csv"), check.names = FALSE)
require_columns(grid, c(
  "run_id", "tf_panel", "pruning_steps", "n_layers", "seed_mode",
  "primary", "secondary", "edges", "nodes", "upstream_candidates",
  "downstream_tfs", "eligible"
), "Threshold grid")
choice <- jsonlite::read_json(file.path(analysis_dir, "selected_sign_consensus", "selection.json"),
                              simplifyVector = TRUE)
choice_cols <- c("run_id", "tf_panel", "pruning_steps", "n_layers", "seed_mode",
                 "primary", "secondary", "edges", "nodes", "upstream_candidates",
                 "downstream_tfs")
# Do not turn the mixed-type data frame into a character matrix: `apply()` pads
# numeric columns with whitespace, which makes an otherwise identical JSON row
# fail to match its CSV counterpart.
selected_idx <- which(
  grid$run_id == as.character(choice$run_id) &
    grid$tf_panel == as.character(choice$tf_panel) &
    grid$pruning_steps == as.integer(choice$pruning_steps) &
    grid$n_layers == as.integer(choice$n_layers) &
    grid$seed_mode == as.character(choice$seed_mode) &
    abs(grid$primary - as.numeric(choice$primary)) < .Machine$double.eps &
    abs(grid$secondary - as.numeric(choice$secondary)) < .Machine$double.eps &
    grid$edges == as.integer(choice$edges) &
    grid$nodes == as.integer(choice$nodes) &
    grid$upstream_candidates == as.integer(choice$upstream_candidates) &
    grid$downstream_tfs == as.integer(choice$downstream_tfs)
)
if (length(selected_idx) != 1L) stop("Selected-sign-consensus JSON does not resolve to exactly one threshold-grid row")
selected_row <- grid[selected_idx, , drop = FALSE]
if (!isTRUE(selected_row$eligible) || selected_row$seed_mode != "sign_consensus") {
  stop("Selected-sign-consensus row is not an eligible sign-consensus extraction")
}

node_summary <- function(status) {
  require_columns(status, c(
    "assignment_id", "pkn_source_node", "pkn_source_chebi", "measurement_feature",
    "candidate_measurement_value", "candidate_in_prepared_topology", "candidate_scored",
    "moon_score", "moon_level", "assignment_sign_agrees"
  ), "Candidate-assignment status")
  grouped <- split(seq_len(nrow(status)), status$pkn_source_node)
  rows <- lapply(names(grouped), function(node) {
    x <- status[grouped[[node]], , drop = FALSE]
    measured_sign <- sort(unique(sign(as.numeric(x$candidate_measurement_value))))
    present <- as.logical(x$candidate_scored)
    data.frame(
      pkn_source_node = node,
      candidate_scored = all(present),
      moon_score = if (all(present)) as.numeric(x$moon_score[1]) else NA_real_,
      moon_level = if (all(present)) as.numeric(x$moon_level[1]) else NA_real_,
      all_assignments_sign_concordant = all(present) && all(as.logical(x$assignment_sign_agrees)),
      measurement_signs = paste(measured_sign, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
sign_consensus_seeds <- function(status) {
  summary <- node_summary(status)
  keep <- summary$candidate_scored & summary$all_assignments_sign_concordant &
    summary$moon_level > 0 & summary$moon_score != 0
  stats::setNames(sign(summary$moon_score[keep]), summary$pkn_source_node[keep])
}
read_run <- function(run_id) {
  run_dir <- file.path(analysis_dir, "runs", run_id)
  fit <- readRDS(file.path(run_dir, "run_state.rds"))
  status <- utils::read.csv(file.path(run_dir, "candidate_assignment_status.csv"), check.names = FALSE)
  list(fit = fit, status = status)
}
extract_checked <- function(run_id, primary, secondary) {
  run <- read_run(run_id)
  seeds <- sign_consensus_seeds(run$status)
  solution <- moon_extract_paths(run$fit$scores, run$fit$pkn_final, seeds, primary, secondary,
                                 input$hepatocyte_gene_weight)
  moon_validate_solution(solution, run$fit, primary, secondary)
  list(fit = run$fit, status = run$status, seeds = seeds, solution = solution)
}

# Independently recreate the materialized selection before changing only its reporting unit.
selected <- extract_checked(selected_row$run_id, selected_row$primary, selected_row$secondary)
raw_sif <- utils::read.csv(file.path(analysis_dir, "selected_sign_consensus", "network_SIF.csv"), check.names = FALSE)
raw_att <- utils::read.csv(file.path(analysis_dir, "selected_sign_consensus", "network_ATT.csv"), check.names = FALSE)
require_columns(raw_sif, c("source", "target", "interaction", "consistency"), "Materialized selected SIF")
require_columns(raw_att, c("source", "score", "level", "type"), "Materialized selected ATT")
if (!setequal(edge_key(raw_sif), edge_key(selected$solution$SIF)) ||
    !setequal(raw_att$source, selected$solution$ATT$source)) {
  stop("Materialized selected network does not reproduce from saved fit and threshold row")
}

raw_sources <- sort(raw_att$source[raw_att$type == "upstream_input"])
if (!length(raw_sources)) stop("Selected solution has no upstream candidate source")
status <- selected$status
status$in_selected_raw_solution <- status$pkn_source_node %in% raw_sources
status$selected_raw_source_count <- sum(status$in_selected_raw_solution)

# Each fingerprint is an observational equivalence group within this direct PKN:
# one representative is retained for an illustrative network, but every candidate
# mapping remains in candidate_hypothesis_summary.csv.
groups <- split(seq_len(nrow(status)), status$edge_fingerprint_group)
fingerprint_rows <- lapply(names(groups), function(group) {
  x <- status[groups[[group]], , drop = FALSE]
  nodes <- sort(unique(x$pkn_source_node))
  selected_nodes <- intersect(nodes, raw_sources)
  representative <- if (length(selected_nodes)) sort(selected_nodes)[1] else ""
  data.frame(
    edge_fingerprint_group = group,
    n_assignment_rows = nrow(x),
    n_candidate_source_nodes = length(nodes),
    candidate_source_nodes = collapse_text(nodes),
    candidate_chemical_names = collapse_text(x$pkn_source_name),
    n_measurement_features = length(unique(x$measurement_feature)),
    measurement_features = collapse_text(x$measurement_feature),
    measurement_values = collapse_text(format(as.numeric(x$candidate_measurement_value), scientific = FALSE, trim = TRUE)),
    mapping_policies = collapse_text(x$primary_mapping_policy),
    n_scored_sign_concordant_assignments = sum(x$assignment_status == "scored_sign_concordant"),
    n_selected_raw_source_nodes = length(selected_nodes),
    selected_raw_source_nodes = collapse_text(selected_nodes),
    fingerprint_representative_source_node = representative,
    fingerprint_representative_chemical_name = if (nzchar(representative)) {
      x$pkn_source_name[match(representative, x$pkn_source_node)]
    } else "",
    selected_in_fingerprint_representative_solution = length(selected_nodes) > 0L,
    stringsAsFactors = FALSE
  )
})
fingerprint_summary <- do.call(rbind, fingerprint_rows)
fingerprint_summary <- fingerprint_summary[order(!fingerprint_summary$selected_in_fingerprint_representative_solution,
                                                 fingerprint_summary$edge_fingerprint_group), , drop = FALSE]
selected_fingerprint_summary <- fingerprint_summary[
  fingerprint_summary$selected_in_fingerprint_representative_solution, , drop = FALSE
]
if (!nrow(selected_fingerprint_summary)) stop("No selected fingerprint group was recovered")

# Select one deterministic representative per selected direct-edge fingerprint.
representatives <- selected_fingerprint_summary$fingerprint_representative_source_node
node_status <- node_summary(status)
score_idx <- match(representatives, node_status$pkn_source_node)
if (anyNA(score_idx) || any(!is.finite(node_status$moon_score[score_idx]))) {
  stop("A selected fingerprint representative has no finite MOON score")
}
representative_seeds <- stats::setNames(sign(node_status$moon_score[score_idx]), representatives)
if (any(representative_seeds == 0)) stop("A selected fingerprint representative has zero MOON score")
representative_solution <- moon_extract_paths(
  selected$fit$scores, selected$fit$pkn_final, representative_seeds,
  selected_row$primary, selected_row$secondary, input$hepatocyte_gene_weight
)
moon_validate_solution(representative_solution, selected$fit,
                       selected_row$primary, selected_row$secondary)

status$in_fingerprint_representative_solution <- status$pkn_source_node %in% representatives
status$candidate_hypothesis_disposition <- ifelse(
  status$in_fingerprint_representative_solution,
  "retained_as_fingerprint_representative",
  ifelse(status$in_selected_raw_solution,
         "same_edge_fingerprint_as_retained_representative",
         ifelse(status$assignment_status == "scored_sign_concordant",
                "sign_concordant_but_not_in_selected_threshold_solution",
                status$assignment_status)))
status <- status[order(status$measurement_feature, status$edge_fingerprint_group,
                       status$pkn_source_node), , drop = FALSE]
utils::write.csv(status, file.path(out, "candidate_hypothesis_summary.csv"), row.names = FALSE)
utils::write.csv(fingerprint_summary, file.path(out, "candidate_edge_fingerprint_summary.csv"), row.names = FALSE)

# Materialize the fingerprint-representative version of the selected graph.
rep_dir <- file.path(out, "fingerprint_representative_solution")
dir.create(rep_dir, recursive = TRUE)
utils::write.csv(representative_solution$SIF, file.path(rep_dir, "network_SIF.csv"), row.names = FALSE)
rep_att <- representative_solution$ATT
att_idx <- match(rep_att$source, status$pkn_source_node)
annotation_cols <- intersect(c(
  "assignment_id", "measurement_feature", "candidate_measurement_value", "pkn_source_chebi",
  "pkn_source_name", "primary_mapping_policy", "edge_fingerprint_group",
  "assignment_sign_agrees", "candidate_hypothesis_disposition"
), names(status))
for (column in annotation_cols) rep_att[[column]] <- status[[column]][att_idx]
utils::write.csv(rep_att, file.path(rep_dir, "network_ATT.csv"), row.names = FALSE)
write_json(list(
  derived_from_selection = as.list(selected_row[, c(
    "run_id", "tf_panel", "pruning_steps", "n_layers", "seed_mode", "primary", "secondary"
  ), drop = FALSE]),
  direct_edge_fingerprint_policy = "One lexical source-node representative per selected identical signed outgoing-edge fingerprint; all alternate feature-to-chemical mappings retained in candidate_hypothesis_summary.csv.",
  raw_selected_source_nodes = raw_sources,
  representative_source_nodes = representatives,
  raw_selected_edges = nrow(raw_sif),
  raw_selected_nodes = nrow(raw_att),
  representative_edges = nrow(representative_solution$SIF),
  representative_nodes = nrow(representative_solution$ATT),
  representative_upstream_fingerprint_groups = nrow(selected_fingerprint_summary),
  representative_observed_tfs = representative_solution$audit$retained_observed
), file.path(rep_dir, "selection.json"))

# Retain evidence only for direct candidate -> protein transitions that survived
# dual-threshold path extraction.
evidence_path <- file.path(mapping_source_dir, "candidate_assignment_regulatory_edges.csv")
if (!file.exists(evidence_path)) stop("Candidate direct-evidence file is absent: ", evidence_path)
evidence <- utils::read.csv(evidence_path, check.names = FALSE)
require_columns(evidence, c("pkn_source_node", "target_node", "interaction", "resources", "source_record_ids"),
                "Candidate direct-evidence table")
first_edges <- representative_solution$SIF[
  representative_solution$SIF$source %in% representatives, c("source", "target", "interaction", "consistency"), drop = FALSE
]
first_edges$.join_key <- paste(first_edges$source, first_edges$target, first_edges$interaction, sep = "\r")
evidence$.join_key <- paste(evidence$pkn_source_node, evidence$target_node, evidence$interaction, sep = "\r")
idx <- match(first_edges$.join_key, evidence$.join_key)
if (anyNA(idx)) stop("A selected direct regulatory edge lacks mapping-ledger evidence")
first_evidence <- cbind(first_edges[, setdiff(names(first_edges), ".join_key"), drop = FALSE],
                        evidence[idx, setdiff(names(evidence), ".join_key"), drop = FALSE])
utils::write.csv(first_evidence, file.path(rep_dir, "selected_direct_regulatory_edges.csv"), row.names = FALSE)

# `moon_extract_paths()` is a sign-consistent topology extractor. Its routes can
# connect nodes assigned nonconsecutive recursive-ULM levels, so it must not be
# used as an explanation of a candidate's ULM score. The package's dedicated
# helper returns the levelwise subnetwork from which each candidate was scored.
# Generate those small, candidate-specific score explanations separately.
score_dir <- file.path(out, "candidate_score_explanations")
dir.create(score_dir, recursive = TRUE)
node_status <- node_summary(status)
score_rows <- list()
for (source_node in sort(unique(status$pkn_source_node))) {
  source_rows <- status[status$pkn_source_node == source_node, , drop = FALSE]
  source_info <- node_status[match(source_node, node_status$pkn_source_node), , drop = FALSE]
  score_subdir <- file.path(score_dir, gsub("[^A-Za-z0-9_.-]", "_", source_node))
  base_row <- data.frame(
    pkn_source_node = source_node,
    pkn_source_chebi = source_rows$pkn_source_chebi[1],
    pkn_source_name = source_rows$pkn_source_name[1],
    edge_fingerprint_group = source_rows$edge_fingerprint_group[1],
    measurement_features = collapse_text(source_rows$measurement_feature),
    measurement_values = collapse_text(format(as.numeric(source_rows$candidate_measurement_value), scientific = FALSE, trim = TRUE)),
    mapping_policies = collapse_text(source_rows$primary_mapping_policy),
    surviving_resources = collapse_text(source_rows$surviving_resources),
    candidate_scored = source_info$candidate_scored,
    candidate_moon_score = source_info$moon_score,
    candidate_moon_level = source_info$moon_level,
    all_assignment_signs_concordant = source_info$all_assignments_sign_concordant,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(source_info$candidate_scored)) {
    score_rows[[length(score_rows) + 1L]] <- cbind(
      base_row,
      score_explanation_status = "candidate_not_scored",
      score_explanation_edges = 0L, score_explanation_nodes = 0L,
      score_explanation_observed_tfs = 0L, observed_tf_nodes = "",
      all_edges_descend_one_score_level = NA, all_edges_sign_consistent = NA,
      all_nodes_on_candidate_to_observed_walk = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  # Candidate inputs were explicitly protected during compression, so their IDs
  # must be present in the compressed score space used by cosmosR's helper.
  if (!source_node %in% selected$fit$scores_compressed$source) {
    stop("Scored candidate was unexpectedly absent from compressed score space: ", source_node)
  }
  score_network <- NULL
  invisible(capture.output(
    score_network <- cosmosR::get_moon_scoring_network(
      source_node, selected$fit$pkn_scoring, selected$fit$scores_compressed
    )
  ))
  require_columns(score_network$SIF, c("source", "target", "interaction"),
                  "Candidate score-explanation SIF")
  require_columns(score_network$ATT, c("source", "score", "level"),
                  "Candidate score-explanation ATT")
  if (!source_node %in% score_network$ATT$source) {
    stop("Candidate score explanation omitted its input node: ", source_node)
  }
  score_map <- stats::setNames(score_network$ATT$score, score_network$ATT$source)
  level_map <- stats::setNames(score_network$ATT$level, score_network$ATT$source)
  if (!all(score_network$SIF$source %in% names(level_map)) ||
      !all(score_network$SIF$target %in% names(level_map))) {
    stop("Candidate score explanation has an edge endpoint missing from ATT: ", source_node)
  }
  score_network$SIF$source_level <- unname(level_map[score_network$SIF$source])
  score_network$SIF$target_level <- unname(level_map[score_network$SIF$target])
  score_network$SIF$edge_sign_consistent <- sign(
    unname(score_map[score_network$SIF$source]) * unname(score_map[score_network$SIF$target])
  ) == score_network$SIF$interaction
  descending <- all(score_network$SIF$source_level == score_network$SIF$target_level + 1L)
  observed_nodes <- score_network$ATT$source[score_network$ATT$level == 0L]
  graph_scoring <- igraph::graph_from_data_frame(score_network$SIF[, c("source", "target")], directed = TRUE)
  root_reachable <- igraph::subcomponent(graph_scoring, source_node, mode = "out")
  reverse_graph <- igraph::reverse_edges(graph_scoring, igraph::E(graph_scoring))
  observable <- unique(unlist(lapply(observed_nodes, function(node) {
    igraph::as_ids(igraph::subcomponent(reverse_graph, node, mode = "out"))
  })))
  all_on_walk <- setequal(score_network$ATT$source, intersect(igraph::as_ids(root_reachable), observable))
  if (!descending || !all_on_walk || !length(observed_nodes)) {
    stop("cosmosR score-explanation helper failed its levelwise candidate-to-observation contract for ", source_node)
  }
  dir.create(score_subdir, recursive = TRUE)
  utils::write.csv(score_network$SIF, file.path(score_subdir, "network_SIF_compressed.csv"), row.names = FALSE)
  utils::write.csv(score_network$ATT, file.path(score_subdir, "network_ATT_compressed.csv"), row.names = FALSE)
  direct <- score_network$SIF[score_network$SIF$source == source_node,
                              c("source", "target", "interaction", "source_level", "target_level", "edge_sign_consistent"), drop = FALSE]
  direct$.join_key <- paste(direct$source, direct$target, direct$interaction, sep = "\r")
  evidence_idx <- match(direct$.join_key, evidence$.join_key)
  if (anyNA(evidence_idx)) stop("Score-explanation direct edge lacks mapping-ledger evidence: ", source_node)
  direct_evidence <- cbind(direct[, setdiff(names(direct), ".join_key"), drop = FALSE],
                           evidence[evidence_idx, setdiff(names(evidence), ".join_key"), drop = FALSE])
  utils::write.csv(direct_evidence, file.path(score_subdir, "direct_regulatory_evidence.csv"), row.names = FALSE)
  score_status <- if (isTRUE(source_info$all_assignments_sign_concordant) &&
                      all(score_network$SIF$edge_sign_consistent)) {
    "layerwise_and_sign_consistent"
  } else if (isTRUE(source_info$all_assignments_sign_concordant)) {
    "layerwise_with_mixed_edge_sign_evidence"
  } else {
    "measurement_sign_discordant"
  }
  score_rows[[length(score_rows) + 1L]] <- cbind(
    base_row,
    score_explanation_status = score_status,
    score_explanation_edges = nrow(score_network$SIF),
    score_explanation_nodes = nrow(score_network$ATT),
    score_explanation_observed_tfs = length(observed_nodes),
    observed_tf_nodes = collapse_text(observed_nodes),
    all_edges_descend_one_score_level = descending,
    all_edges_sign_consistent = all(score_network$SIF$edge_sign_consistent),
    all_nodes_on_candidate_to_observed_walk = all_on_walk,
    stringsAsFactors = FALSE
  )
}
score_summary <- do.call(rbind, score_rows)
score_summary <- score_summary[order(score_summary$score_explanation_status,
                                     score_summary$measurement_features,
                                     score_summary$pkn_source_node), , drop = FALSE]
utils::write.csv(score_summary, file.path(out, "candidate_score_explanation_summary.csv"), row.names = FALSE)
summary_idx <- match(status$pkn_source_node, score_summary$pkn_source_node)
for (column in c("score_explanation_status", "score_explanation_edges", "score_explanation_nodes",
                 "score_explanation_observed_tfs", "observed_tf_nodes",
                 "all_edges_descend_one_score_level", "all_edges_sign_consistent",
                 "all_nodes_on_candidate_to_observed_walk")) {
  status[[column]] <- score_summary[[column]][summary_idx]
}
utils::write.csv(status, file.path(out, "candidate_hypothesis_summary.csv"), row.names = FALSE)

# One directed shortest path per representative-to-observed-TF pair. This is a
# readable route index, not a claim that the shortest route is uniquely operative.
graph <- igraph::graph_from_data_frame(representative_solution$SIF[, c("source", "target")], directed = TRUE)
observed <- representative_solution$ATT$source[representative_solution$ATT$type == "level0"]
interaction_map <- stats::setNames(
  as.character(representative_solution$SIF$interaction),
  paste(representative_solution$SIF$source, representative_solution$SIF$target, sep = "\r")
)
att_score <- stats::setNames(representative_solution$ATT$score, representative_solution$ATT$source)
path_rows <- list()
for (source_node in representatives) {
  for (target_node in sort(observed)) {
    path <- igraph::shortest_paths(graph, from = source_node, to = target_node,
                                   mode = "out", output = "vpath")$vpath[[1]]
    nodes <- igraph::as_ids(path)
    if (!length(nodes)) next
    pairs <- if (length(nodes) > 1L) paste(head(nodes, -1L), tail(nodes, -1L), sep = "\r") else character()
    path_rows[[length(path_rows) + 1L]] <- data.frame(
      source_node = source_node,
      target_observed_tf = target_node,
      path_length_edges = length(nodes) - 1L,
      path_nodes = paste(nodes, collapse = " -> "),
      path_interactions = paste(interaction_map[pairs], collapse = ";"),
      source_moon_score = unname(att_score[source_node]),
      target_tf_activity = unname(att_score[target_node]),
      stringsAsFactors = FALSE
    )
  }
}
paths <- if (length(path_rows)) do.call(rbind, path_rows) else data.frame(
  source_node = character(), target_observed_tf = character(), path_length_edges = integer(),
  path_nodes = character(), path_interactions = character(), source_moon_score = numeric(),
  target_tf_activity = numeric()
)
utils::write.csv(paths, file.path(rep_dir, "candidate_to_observed_tf_shortest_paths.csv"), row.names = FALSE)

# Recalculate every eligible sign-consensus extraction and report the number of
# independent direct-edge fingerprints rather than raw duplicated source nodes.
eligible <- grid[grid$eligible & grid$seed_mode == "sign_consensus", , drop = FALSE]
eligible_rows <- lapply(seq_len(nrow(eligible)), function(i) {
  x <- eligible[i, , drop = FALSE]
  reconstructed <- extract_checked(x$run_id, x$primary, x$secondary)
  sources <- reconstructed$solution$ATT$source[reconstructed$solution$ATT$type == "upstream_input"]
  source_rows <- reconstructed$status[reconstructed$status$pkn_source_node %in% sources, , drop = FALSE]
  groups <- sort(unique(source_rows$edge_fingerprint_group))
  data.frame(
    x,
    independent_edge_fingerprint_groups = length(groups),
    edge_fingerprint_groups = collapse_text(groups),
    unique_measurement_features = length(unique(source_rows$measurement_feature)),
    measurement_features = collapse_text(source_rows$measurement_feature),
    source_nodes = collapse_text(sources),
    stringsAsFactors = FALSE
  )
})
eligible_fingerprint <- do.call(rbind, eligible_rows)
eligible_fingerprint <- eligible_fingerprint[order(-eligible_fingerprint$primary,
                                                   -eligible_fingerprint$secondary,
                                                   eligible_fingerprint$edges,
                                                   eligible_fingerprint$run_id), , drop = FALSE]
utils::write.csv(eligible_fingerprint, file.path(out, "eligible_sign_consensus_fingerprint_summary.csv"), row.names = FALSE)

utils::write.csv(data.frame(path = c(
  "scripts/moon/ancillary/audit_relaxed_lipid_tf_candidates.R",
  file.path(analysis_dir, "dependency_manifest.csv"),
  file.path(analysis_dir, "threshold_grid.csv"),
  file.path(analysis_dir, "selected_sign_consensus", "selection.json"),
  file.path(analysis_dir, "selected_sign_consensus", "network_SIF.csv"),
  file.path(analysis_dir, "selected_sign_consensus", "network_ATT.csv"),
  evidence_path
), md5 = unname(tools::md5sum(c(
  "scripts/moon/ancillary/audit_relaxed_lipid_tf_candidates.R",
  file.path(analysis_dir, "dependency_manifest.csv"),
  file.path(analysis_dir, "threshold_grid.csv"),
  file.path(analysis_dir, "selected_sign_consensus", "selection.json"),
  file.path(analysis_dir, "selected_sign_consensus", "network_SIF.csv"),
  file.path(analysis_dir, "selected_sign_consensus", "network_ATT.csv"),
  evidence_path
))), stringsAsFactors = FALSE), file.path(out, "dependency_manifest.csv"), row.names = FALSE)

writeLines(c(
  "# Fingerprint-aware relaxed lipid-to-TF candidate audit",
  "",
  "This audit reads a completed score-first relaxed MOON run and does not refit its PKN, pruning, compression, RNA-coherence loop, or node scores.",
  "",
  "Candidate chemicals are topology seeds in the parent analysis. Factor 1 lipid weights are assignment-level annotations used for post-score sign agreement; they are not ULM inputs. A feature may map to multiple candidates, and a candidate may share an identical signed direct-outgoing-edge fingerprint with another candidate. The audit retains every assignment in `candidate_hypothesis_summary.csv` and collapses only selected identical fingerprint groups in `fingerprint_representative_solution/` for a non-duplicated illustrative connectivity network.",
  "",
  "`eligible_sign_consensus_fingerprint_summary.csv` recalculates the compact sign-consensus candidates and reports both raw source count and the number of independent direct-edge fingerprint groups. `candidate_edge_fingerprint_summary.csv` identifies the counterparts represented by each retained fingerprint. `candidate_to_observed_tf_shortest_paths.csv` provides one directed shortest route per retained representative and observed TF; it is a route index, not proof that a shortest path is biologically active. Those connectivity paths are not necessarily consecutive recursive-MOON score layers and must not be used as a ULM-score explanation.",
  "",
  "`candidate_score_explanation_summary.csv` and `candidate_score_explanations/` are the candidate-local, layerwise score traces returned by `cosmosR::get_moon_scoring_network()`. They identify the observed TFs that numerically support each candidate score, flag mixed edge-sign evidence, and are the appropriate files for score interpretation. The direct regulatory-edge provenance for the illustrative connectivity network is in `fingerprint_representative_solution/selected_direct_regulatory_edges.csv`. In the current run, all selected direct evidence remains labelled exploratory, including historical STITCH v5 source data and relaxed chemical identity mappings."
), file.path(out, "README.md"))
write_json(list(
  analysis = "factor1_lipid_to_tf_relaxed_candidate_audit",
  source_analysis = analysis_dir,
  selected_raw_source_nodes = raw_sources,
  selected_raw_source_count = length(raw_sources),
  selected_independent_fingerprint_groups = nrow(selected_fingerprint_summary),
  representative_source_nodes = representatives,
  representative_network_edges = nrow(representative_solution$SIF),
  representative_network_nodes = nrow(representative_solution$ATT),
  observed_tfs = representative_solution$audit$retained_observed,
  shortest_paths_written = nrow(paths),
  candidate_nodes_scored = sum(score_summary$candidate_scored),
  layerwise_sign_consistent_candidate_nodes = sum(score_summary$score_explanation_status == "layerwise_and_sign_consistent"),
  interpretation = "Fingerprint collapse is a reporting deduplication of identical signed direct outgoing PKN profiles. The representative connectivity network is not a ULM score trace; candidate_score_explanations contains the layerwise score support. Neither output resolves chemical identity or establishes ligand activity."
), file.path(out, "audit_summary.json"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed fingerprint-aware relaxed lipid-to-TF candidate audit: ",
    nrow(status), " mapping rows; ", nrow(selected_fingerprint_summary),
    " selected independent direct-edge fingerprint group(s); ", nrow(paths),
    " indexed connectivity paths; ", sum(score_summary$candidate_scored),
    " candidate-local score traces. Results: ", out, "\n", sep = "")
