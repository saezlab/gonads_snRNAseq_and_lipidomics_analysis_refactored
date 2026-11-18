#!/usr/bin/env Rscript
# Candidate-local score audit for a completed selected relaxed lipid-to-TF refit.
# This script reads saved state only; it never repeats fitting or coherence.
options(stringsAsFactors = FALSE, digits = 17)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run this entry point with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../.."), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R ",
    "[--analysis=PATH] [--output=NEW_DIRECTORY]\n",
    "Audits every relaxed mapping candidate against one completed selected refit. ",
    "Writes candidate-local cosmosR score traces, mapping/sign status and direct-edge evidence without refitting.\n"
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
source("scripts/moon/lib/relaxed_lipid_tf.R")

analysis_dir <- normalizePath(option("--analysis", "results/causal_modelling/factor1_lipid_to_tf_relaxed"))
out <- moon_new_output(option("--output", file.path(analysis_dir, "candidate_audit")))
required <- c("run_state.rds", "candidate_assignment_status.csv", "network_SIF.csv", "network_ATT.csv",
              "dependency_manifest.csv", "analysis_config.json")
missing <- required[!file.exists(file.path(analysis_dir, required))]
if (length(missing)) stop("Completed selected relaxed analysis is missing: ", paste(missing, collapse = ", "))

deps <- utils::read.csv(file.path(analysis_dir, "dependency_manifest.csv"), check.names = FALSE)
relaxed_required_columns(deps, c("path", "md5"), "Selected relaxed dependency manifest")
for (i in seq_len(nrow(deps))) {
  if (!file.exists(deps$path[i]) ||
      !identical(unname(tools::md5sum(deps$path[i])), as.character(deps$md5[i]))) {
    stop("Selected relaxed dependency changed since fitting: ", deps$path[i])
  }
}
cfg <- jsonlite::read_json(file.path(analysis_dir, "analysis_config.json"), simplifyVector = TRUE)
fit <- readRDS(file.path(analysis_dir, "run_state.rds"))
status <- utils::read.csv(file.path(analysis_dir, "candidate_assignment_status.csv"), check.names = FALSE)
relaxed_required_columns(status, c(
  "assignment_id", "measurement_feature", "candidate_measurement_value", "pkn_source_node",
  "pkn_source_chebi", "pkn_source_name", "primary_mapping_policy", "edge_fingerprint_group",
  "candidate_scored", "moon_score", "moon_level", "assignment_sign_agrees", "assignment_status"
), "Selected relaxed candidate status")
evidence <- utils::read.csv(cfg$direct_evidence, check.names = FALSE)
relaxed_required_columns(evidence, c("pkn_source_node", "target_node", "interaction"),
                         "Relaxed direct-evidence ledger")
evidence$.key <- paste(evidence$pkn_source_node, evidence$target_node, evidence$interaction, sep = "\r")

node_summary <- relaxed_candidate_node_summary(status)
fingerprint_rows <- lapply(split(seq_len(nrow(status)), status$edge_fingerprint_group), function(index) {
  x <- status[index, , drop = FALSE]
  data.frame(
    edge_fingerprint_group = x$edge_fingerprint_group[1],
    n_assignment_rows = nrow(x),
    n_candidate_source_nodes = length(unique(x$pkn_source_node)),
    candidate_source_nodes = relaxed_collapse(x$pkn_source_node),
    candidate_chemical_names = relaxed_collapse(x$pkn_source_name),
    measurement_features = relaxed_collapse(x$measurement_feature),
    mapping_policies = relaxed_collapse(x$primary_mapping_policy),
    scored_sign_concordant_assignments = sum(as.logical(x$assignment_sign_agrees)),
    stringsAsFactors = FALSE
  )
})
fingerprints <- do.call(rbind, fingerprint_rows)
fingerprints <- fingerprints[order(fingerprints$measurement_features, fingerprints$edge_fingerprint_group), , drop = FALSE]
utils::write.csv(fingerprints, file.path(out, "candidate_edge_fingerprint_summary.csv"), row.names = FALSE)

score_rows <- list()
score_dir <- file.path(out, "candidate_score_explanations")
dir.create(score_dir, recursive = TRUE)
for (source_node in sort(unique(status$pkn_source_node))) {
  source_rows <- status[status$pkn_source_node == source_node, , drop = FALSE]
  source_info <- node_summary[match(source_node, node_summary$pkn_source_node), , drop = FALSE]
  base <- data.frame(
    pkn_source_node = source_node,
    pkn_source_chebi = source_rows$pkn_source_chebi[1],
    pkn_source_name = source_rows$pkn_source_name[1],
    edge_fingerprint_group = source_rows$edge_fingerprint_group[1],
    measurement_features = relaxed_collapse(source_rows$measurement_feature),
    measurement_values = relaxed_collapse(format(as.numeric(source_rows$candidate_measurement_value), scientific = FALSE, trim = TRUE)),
    mapping_policies = relaxed_collapse(source_rows$primary_mapping_policy),
    candidate_scored = source_info$candidate_scored,
    candidate_moon_score = source_info$moon_score,
    candidate_moon_level = source_info$moon_level,
    all_assignment_signs_concordant = source_info$all_assignments_sign_concordant,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(source_info$candidate_scored)) {
    score_rows[[length(score_rows) + 1L]] <- cbind(
      base, score_explanation_status = "candidate_not_scored", score_explanation_edges = 0L,
      score_explanation_nodes = 0L, score_explanation_observed_tfs = 0L, observed_tf_nodes = "",
      all_edges_descend_one_score_level = NA, all_edges_sign_consistent = NA,
      all_nodes_on_candidate_to_observed_walk = NA, stringsAsFactors = FALSE
    )
    next
  }
  if (!source_node %in% fit$scores_compressed$source) {
    stop("A scored candidate was absent from the protected compressed score space: ", source_node)
  }
  score_network <- NULL
  invisible(capture.output(
    score_network <- cosmosR::get_moon_scoring_network(source_node, fit$pkn_scoring, fit$scores_compressed)
  ))
  relaxed_required_columns(score_network$SIF, c("source", "target", "interaction"),
                           "Candidate score-explanation SIF")
  relaxed_required_columns(score_network$ATT, c("source", "score", "level"),
                           "Candidate score-explanation ATT")
  score_map <- stats::setNames(score_network$ATT$score, score_network$ATT$source)
  level_map <- stats::setNames(score_network$ATT$level, score_network$ATT$source)
  if (!source_node %in% names(score_map) ||
      !all(score_network$SIF$source %in% names(score_map)) ||
      !all(score_network$SIF$target %in% names(score_map))) {
    stop("Candidate score explanation has inconsistent endpoints: ", source_node)
  }
  score_network$SIF$source_level <- unname(level_map[score_network$SIF$source])
  score_network$SIF$target_level <- unname(level_map[score_network$SIF$target])
  score_network$SIF$edge_sign_consistent <- sign(
    unname(score_map[score_network$SIF$source]) * unname(score_map[score_network$SIF$target])
  ) == score_network$SIF$interaction
  graph <- igraph::graph_from_data_frame(score_network$SIF[, c("source", "target")], directed = TRUE)
  observed <- score_network$ATT$source[score_network$ATT$level == 0L]
  forward <- igraph::as_ids(igraph::subcomponent(graph, source_node, mode = "out"))
  reversed <- igraph::reverse_edges(graph, igraph::E(graph))
  backward <- unique(unlist(lapply(observed, function(node) {
    igraph::as_ids(igraph::subcomponent(reversed, node, mode = "out"))
  })))
  all_on_walk <- setequal(score_network$ATT$source, intersect(forward, backward))
  descending <- all(score_network$SIF$source_level == score_network$SIF$target_level + 1L)
  source_dir <- file.path(score_dir, gsub("[^A-Za-z0-9_.-]", "_", source_node))
  dir.create(source_dir, recursive = TRUE)
  utils::write.csv(score_network$SIF, file.path(source_dir, "network_SIF_compressed.csv"), row.names = FALSE)
  utils::write.csv(score_network$ATT, file.path(source_dir, "network_ATT_compressed.csv"), row.names = FALSE)
  direct <- score_network$SIF[score_network$SIF$source == source_node,
                              c("source", "target", "interaction", "source_level", "target_level", "edge_sign_consistent"), drop = FALSE]
  direct$.key <- paste(direct$source, direct$target, direct$interaction, sep = "\r")
  evidence_index <- match(direct$.key, evidence$.key)
  if (length(direct) && anyNA(evidence_index)) {
    stop("Candidate score explanation direct edge lacks mapping evidence: ", source_node)
  }
  direct_evidence <- if (nrow(direct)) {
    cbind(direct[, setdiff(names(direct), ".key"), drop = FALSE],
          evidence[evidence_index, setdiff(names(evidence), ".key"), drop = FALSE])
  } else {
    data.frame(source = character(), target = character(), interaction = numeric(), source_level = numeric(),
               target_level = numeric(), edge_sign_consistent = logical(), stringsAsFactors = FALSE)
  }
  utils::write.csv(direct_evidence, file.path(source_dir, "direct_regulatory_evidence.csv"), row.names = FALSE)
  score_status <- if (!descending || !all_on_walk || !length(observed)) {
    "score_trace_requires_topology_review"
  } else if (isTRUE(source_info$all_assignments_sign_concordant) &&
                      all(score_network$SIF$edge_sign_consistent)) {
    "layerwise_and_sign_consistent"
  } else if (isTRUE(source_info$all_assignments_sign_concordant)) {
    "layerwise_with_mixed_edge_sign_evidence"
  } else {
    "measurement_sign_discordant"
  }
  score_rows[[length(score_rows) + 1L]] <- cbind(
    base, score_explanation_status = score_status,
    score_explanation_edges = nrow(score_network$SIF), score_explanation_nodes = nrow(score_network$ATT),
    score_explanation_observed_tfs = length(observed), observed_tf_nodes = relaxed_collapse(observed),
    all_edges_descend_one_score_level = descending,
    all_edges_sign_consistent = all(score_network$SIF$edge_sign_consistent),
    all_nodes_on_candidate_to_observed_walk = all_on_walk, stringsAsFactors = FALSE
  )
}
score_summary <- do.call(rbind, score_rows)
score_summary <- score_summary[order(score_summary$score_explanation_status,
                                     score_summary$measurement_features, score_summary$pkn_source_node), , drop = FALSE]
utils::write.csv(score_summary, file.path(out, "candidate_score_explanation_summary.csv"), row.names = FALSE)

summary_index <- match(status$pkn_source_node, score_summary$pkn_source_node)
for (field in setdiff(names(score_summary), c("pkn_source_node", "pkn_source_chebi", "pkn_source_name",
                                               "edge_fingerprint_group", "measurement_features", "measurement_values",
                                               "mapping_policies", "candidate_scored", "candidate_moon_score",
                                               "candidate_moon_level", "all_assignment_signs_concordant"))) {
  status[[field]] <- score_summary[[field]][summary_index]
}
utils::write.csv(status, file.path(out, "candidate_hypothesis_summary.csv"), row.names = FALSE)
relaxed_write_json(list(
  analysis = "factor1_lipid_to_tf_relaxed_candidate_audit", source_analysis = analysis_dir,
  mapping_rows = nrow(status), candidate_nodes = nrow(node_summary),
  scored_candidates = sum(score_summary$candidate_scored),
  layerwise_sign_consistent_candidates = sum(score_summary$score_explanation_status == "layerwise_and_sign_consistent"),
  interpretation = "Candidate-local score traces explain recursive ULM support within the saved selected fit. Mapping remains conditional, copied lipid values are post-score annotations, and direct PKN evidence does not establish liver exposure or causal regulation."
), file.path(out, "audit_summary.json"))
writeLines(c(
  "# Relaxed lipid-to-TF candidate audit",
  "",
  "This ancillary audit reads the completed selected refit and does not repeat PKN preparation, pruning, compression, RNA coherence or ULM scoring.",
  "`candidate_hypothesis_summary.csv` retains every feature-to-chemical mapping, its copied Factor 1 loading, candidate score/status and score-explanation flags.",
  "`candidate_edge_fingerprint_summary.csv` prevents chemically unresolved candidates with identical signed direct PKN neighborhoods from being counted as independent evidence.",
  "`candidate_score_explanations/` contains the compressed, layerwise `get_moon_scoring_network()` trace and direct mapping evidence for each scored candidate.",
  "These traces support conditional hypotheses only; they do not resolve lipid identity or demonstrate signaling in hepatocytes."
), file.path(out, "README.md"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed relaxed candidate audit: ", nrow(status), " mapping rows and ",
    sum(score_summary$candidate_scored), " scored candidate nodes. Results: ", out, "\n", sep = "")
