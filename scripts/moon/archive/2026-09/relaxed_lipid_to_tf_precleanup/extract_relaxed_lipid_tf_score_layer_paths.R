#!/usr/bin/env Rscript
# Re-extract saved relaxed lipid -> TF MOON fits as layer-consecutive score paths.
# This script never refits pruning, compression, RNA coherence, or ULM scores.
options(stringsAsFactors = FALSE, digits = 17)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript scripts/moon/ancillary/extract_relaxed_lipid_tf_score_layer_paths.R ",
    "[--analysis=PATH] [--output=PATH] [--self-test]\n",
    "Re-extracts a completed relaxed lipid-to-TF MOON analysis using only ",
    "edges that descend exactly one recursive ULM score layer.\n"
  ))
  quit(save = "no")
}
stopifnot(all(grepl("^--(analysis|output)=.+$|^--self-test$", args)))
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}

source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}
collapse_text <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  paste(x, collapse = ";")
}
edge_key <- function(x) paste(x$source, x$target, x$interaction, sep = "\r")
write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", null = "null", digits = NA)
}

# This differs deliberately from moon_extract_paths(): every SIF edge is an
# actual recursive scoring transition, source level = target level + 1.
extract_score_layer_paths <- function(scores, pkn, upstream_input,
                                      primary, secondary, RNA_input = NULL,
                                      pkn_already_checked = FALSE) {
  require_columns(scores, c("source", "score", "level"), "MOON scores")
  stopifnot(!anyDuplicated(scores$source), all(is.finite(scores$score)),
            all(is.finite(scores$level)), primary >= secondary, secondary >= 0,
            is.numeric(upstream_input), !is.null(names(upstream_input)),
            !anyDuplicated(names(upstream_input)), all(is.finite(upstream_input)))
  if (!pkn_already_checked) pkn <- moon_check_pkn(pkn)
  scores <- scores[, c("source", "score", "level"), drop = FALSE]
  scores$source <- as.character(scores$source)
  observed <- scores$source[scores$level == 0L & scores$score != 0]
  primary_nodes <- scores$source[scores$level > 0L & abs(scores$score) > primary]
  secondary_nodes <- union(observed, scores$source[scores$level > 0L & abs(scores$score) > secondary])
  score_map <- stats::setNames(scores$score, scores$source)
  level_map <- stats::setNames(scores$level, scores$source)
  input_index <- match(names(upstream_input), scores$source)
  input_ok <- !is.na(input_index)
  input_ok[input_ok] <- scores$level[input_index[input_ok]] > 0L &
    sign(scores$score[input_index[input_ok]]) == sign(upstream_input[input_ok])
  seeds <- intersect(names(upstream_input)[input_ok], secondary_nodes)
  failed_inputs <- setdiff(intersect(names(upstream_input), scores$source), seeds)
  secondary_nodes <- setdiff(secondary_nodes, failed_inputs)
  primary_nodes <- intersect(primary_nodes, secondary_nodes)
  net <- pkn[pkn$source %in% secondary_nodes & pkn$target %in% secondary_nodes &
               (pkn$source %in% primary_nodes | pkn$target %in% primary_nodes) &
               !pkn$source %in% observed, , drop = FALSE]
  candidate_edges <- nrow(net)
  net <- net[level_map[net$source] == level_map[net$target] + 1L, , drop = FALSE]
  layer_edges <- nrow(net)
  net$consistency <- sign(score_map[net$source] * score_map[net$target]) == net$interaction
  net <- net[net$consistency, , drop = FALSE]
  consistent_edges <- nrow(net)
  empty <- function() list(
    SIF = data.frame(source = character(), target = character(), interaction = numeric(), consistency = logical()),
    ATT = data.frame(source = character(), score = numeric(), level = numeric(), type = character(), RNA_input = numeric()),
    audit = list(primary = primary, secondary = secondary, score_layer_constrained = TRUE,
                 candidate_edges = candidate_edges, score_layer_edges = layer_edges,
                 consistent_edges = consistent_edges, final_edges = 0L, final_nodes = 0L,
                 retained_upstream = 0L, retained_observed = 0L,
                 all_edges_descend_one_score_level = TRUE,
                 all_nodes_on_seed_to_observation_walk = TRUE)
  )
  if (!nrow(net) || !length(seeds) || !length(observed)) return(empty())
  graph <- igraph::graph_from_data_frame(net[, c("source", "target")], directed = TRUE)
  present <- igraph::V(graph)$name
  seeds <- intersect(seeds, present)
  observed <- intersect(observed, present)
  if (!length(seeds) || !length(observed)) return(empty())
  reachable <- function(graph, start, reverse = FALSE) {
    root <- "__MOON_SCORE_LAYER_ROOT__"
    while (root %in% igraph::V(graph)$name) root <- paste0(root, "_")
    if (reverse) graph <- igraph::reverse_edges(graph, igraph::E(graph))
    graph <- igraph::add_vertices(graph, 1L, name = root)
    graph <- igraph::add_edges(graph, as.vector(rbind(rep(root, length(start)), start)))
    setdiff(igraph::as_ids(igraph::subcomponent(graph, root, mode = "out")), root)
  }
  keep <- intersect(reachable(graph, seeds), reachable(graph, observed, reverse = TRUE))
  net <- net[net$source %in% keep & net$target %in% keep, , drop = FALSE]
  if (!nrow(net)) return(empty())
  nodes <- sort(unique(c(net$source, net$target)))
  att <- scores[match(nodes, scores$source), , drop = FALSE]
  att$type <- ifelse(att$source %in% seeds, "upstream_input",
                     ifelse(att$level == 0L, "level0", "other"))
  if (!is.null(RNA_input)) att$RNA_input <- unname(RNA_input[att$source]) else att$RNA_input <- NA_real_
  final_graph <- igraph::graph_from_data_frame(net[, c("source", "target")], directed = TRUE)
  retained_seeds <- intersect(seeds, nodes)
  retained_observed <- intersect(observed, nodes)
  verified <- intersect(reachable(final_graph, retained_seeds),
                        reachable(final_graph, retained_observed, reverse = TRUE))
  stopifnot(setequal(nodes, verified), all(net$consistency),
            all(level_map[net$source] == level_map[net$target] + 1L),
            setequal(nodes, att$source))
  rownames(net) <- rownames(att) <- NULL
  list(
    SIF = net,
    ATT = att,
    audit = list(primary = primary, secondary = secondary, score_layer_constrained = TRUE,
                 candidate_edges = candidate_edges, score_layer_edges = layer_edges,
                 consistent_edges = consistent_edges, final_edges = nrow(net),
                 final_nodes = nrow(att), retained_upstream = length(retained_seeds),
                 retained_observed = length(retained_observed),
                 all_edges_descend_one_score_level = TRUE,
                 all_nodes_on_seed_to_observation_walk = TRUE)
  )
}
validate_score_layer_solution <- function(solution, fit, primary, secondary) {
  if (!nrow(solution$SIF)) return(invisible(TRUE))
  score_map <- stats::setNames(fit$scores$score, fit$scores$source)
  level_map <- stats::setNames(fit$scores$level, fit$scores$source)
  stopifnot(
    all(edge_key(solution$SIF) %in% edge_key(fit$pkn_final)),
    all(solution$SIF$source %in% names(level_map)),
    all(solution$SIF$target %in% names(level_map)),
    all(level_map[solution$SIF$source] == level_map[solution$SIF$target] + 1L),
    all(sign(score_map[solution$SIF$source] * score_map[solution$SIF$target]) == solution$SIF$interaction),
    all(abs(solution$ATT$score[solution$ATT$level > 0L]) > secondary),
    solution$audit$all_nodes_on_seed_to_observation_walk
  )
  invisible(TRUE)
}
self_test <- function() {
  net <- data.frame(source = c("U", "A", "U"), target = c("A", "T", "T"), interaction = 1)
  scores <- data.frame(source = c("U", "A", "T"), score = c(2, 2, 1), level = c(2, 1, 0))
  out <- extract_score_layer_paths(scores, net, c(U = 1), primary = 1.5, secondary = 1)
  stopifnot(nrow(out$SIF) == 2L, all(out$SIF$source %in% c("U", "A")),
            out$audit$all_edges_descend_one_score_level)
  wrong <- extract_score_layer_paths(scores, net, c(U = -1), primary = 1.5, secondary = 1)
  stopifnot(nrow(wrong$SIF) == 0L)
  cat("PASS: strict score-layer extraction rejects a level-skipping shortcut and wrong input sign.\n")
}
if ("--self-test" %in% args) {
  self_test()
  quit(save = "no")
}

analysis_dir <- normalizePath(option("--analysis", "results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory"))
output_dir <- option("--output", file.path(analysis_dir, "score_layer_constrained"))
if (!dir.exists(analysis_dir)) stop("Analysis directory does not exist: ", analysis_dir)
out <- moon_new_output(output_dir)

# Validate the original scored-fit dependency contract before derived extraction.
deps <- utils::read.csv(file.path(analysis_dir, "dependency_manifest.csv"), check.names = FALSE)
require_columns(deps, c("path", "md5"), "Scored-fit dependency manifest")
for (i in seq_len(nrow(deps))) {
  if (!file.exists(deps$path[i])) stop("Declared scored-fit dependency is missing: ", deps$path[i])
  if (!identical(unname(tools::md5sum(deps$path[i])), as.character(deps$md5[i]))) {
    stop("Declared scored-fit dependency changed: ", deps$path[i])
  }
}
manifest <- jsonlite::read_json(file.path(analysis_dir, "run_manifest.json"), simplifyVector = FALSE)
input <- readRDS(manifest$config$inputs)
mapping_dir <- dirname(manifest$config$mapping)
evidence_path <- file.path(mapping_dir, "candidate_assignment_regulatory_edges.csv")
if (!file.exists(evidence_path)) stop("Missing mapping-ledger direct-evidence table: ", evidence_path)
evidence <- utils::read.csv(evidence_path, check.names = FALSE)
require_columns(evidence, c("pkn_source_node", "target_node", "interaction"), "Mapping-ledger direct evidence")
evidence$.key <- paste(evidence$pkn_source_node, evidence$target_node, evidence$interaction, sep = "\r")

grid <- utils::read.csv(file.path(analysis_dir, "threshold_grid.csv"), check.names = FALSE)
require_columns(grid, c("run_id", "tf_panel", "pruning_steps", "n_layers", "seed_mode", "primary", "secondary"),
                "Threshold grid")
# The topology mode uses inferred score signs as seed values. This report instead
# tests the measured assignment sign, so it only evaluates the pre-existing
# sign-consensus candidate panel.
grid <- grid[grid$seed_mode == "sign_consensus", , drop = FALSE]

node_summary <- function(status) {
  require_columns(status, c("pkn_source_node", "measurement_feature", "candidate_measurement_value",
                            "candidate_scored", "moon_score", "moon_level", "assignment_sign_agrees",
                            "edge_fingerprint_group", "pkn_source_name"), "Candidate status")
  groups <- split(seq_len(nrow(status)), status$pkn_source_node)
  rows <- lapply(names(groups), function(node) {
    x <- status[groups[[node]], , drop = FALSE]
    values <- sort(unique(sign(as.numeric(x$candidate_measurement_value))))
    data.frame(
      pkn_source_node = node,
      candidate_score = if (all(x$candidate_scored)) as.numeric(x$moon_score[1]) else NA_real_,
      candidate_level = if (all(x$candidate_scored)) as.integer(x$moon_level[1]) else NA_integer_,
      candidate_scored = all(x$candidate_scored),
      all_assignment_signs_concordant = all(x$candidate_scored) && all(x$assignment_sign_agrees),
      measurement_sign = if (length(values) == 1L) values else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
seeds_from_status <- function(status) {
  summary <- node_summary(status)
  keep <- summary$candidate_scored & summary$all_assignment_signs_concordant &
    summary$candidate_level > 0L & summary$candidate_score != 0
  stats::setNames(summary$measurement_sign[keep], summary$pkn_source_node[keep])
}

rows <- list()
materialized <- list()
run_cache <- new.env(parent = emptyenv())
read_cached_run <- function(run_id) {
  if (!exists(run_id, envir = run_cache, inherits = FALSE)) {
    run_dir <- file.path(analysis_dir, "runs", run_id)
    fit <- readRDS(file.path(run_dir, "run_state.rds"))
    # The saved uncompressed PKN is already the canonical final scoring graph.
    # Check it once per fit rather than repeating full canonical cleanup for all
    # threshold pairs from the same scores.
    fit$pkn_final <- moon_check_pkn(fit$pkn_final)
    status <- utils::read.csv(file.path(run_dir, "candidate_assignment_status.csv"), check.names = FALSE)
    assign(run_id, list(fit = fit, status = status, seeds = seeds_from_status(status)), envir = run_cache)
  }
  get(run_id, envir = run_cache, inherits = FALSE)
}
for (i in seq_len(nrow(grid))) {
  config <- grid[i, , drop = FALSE]
  run <- read_cached_run(config$run_id)
  fit <- run$fit
  status <- run$status
  solution <- extract_score_layer_paths(fit$scores, fit$pkn_final, run$seeds,
                                        config$primary, config$secondary,
                                        input$hepatocyte_gene_weight,
                                        pkn_already_checked = TRUE)
  validate_score_layer_solution(solution, fit, config$primary, config$secondary)
  sources <- solution$ATT$source[solution$ATT$type == "upstream_input"]
  source_rows <- status[status$pkn_source_node %in% sources, , drop = FALSE]
  fingerprints <- sort(unique(source_rows$edge_fingerprint_group))
  features <- sort(unique(source_rows$measurement_feature))
  # Candidate scores must clear a nonzero secondary threshold for a compact
  # result to be called score-supported. This prevents the prior p=2/s=0
  # connectivity screen from being counted as an eligible candidate result.
  eligible <- nrow(solution$SIF) > 0L && nrow(solution$SIF) < 500L &&
    length(fingerprints) >= 1L && solution$audit$retained_observed >= 2L && config$secondary > 0
  rows[[length(rows) + 1L]] <- data.frame(
    config,
    score_layer_edges = nrow(solution$SIF),
    score_layer_nodes = nrow(solution$ATT),
    raw_upstream_candidate_nodes = length(sources),
    independent_edge_fingerprint_groups = length(fingerprints),
    edge_fingerprint_groups = collapse_text(fingerprints),
    unique_measurement_features = length(features),
    measurement_features = collapse_text(features),
    upstream_source_nodes = collapse_text(sources),
    observed_tfs = collapse_text(solution$ATT$source[solution$ATT$type == "level0"]),
    score_layer_eligible = eligible,
    stringsAsFactors = FALSE
  )
  if (eligible) {
    key <- paste(config$run_id, paste0("p", format(config$primary, trim = TRUE, scientific = FALSE)),
                 paste0("s", format(config$secondary, trim = TRUE, scientific = FALSE)), sep = "__")
    materialized[[key]] <- list(config = config, solution = solution, status = status, fingerprints = fingerprints)
  }
}
summary <- do.call(rbind, rows)
summary <- summary[order(summary$tf_panel, summary$pruning_steps, summary$n_layers,
                         summary$primary, summary$secondary), , drop = FALSE]
utils::write.csv(summary, file.path(out, "score_layer_threshold_grid.csv"), row.names = FALSE)
eligible <- summary[summary$score_layer_eligible, , drop = FALSE]
eligible <- eligible[order(eligible$tf_panel, eligible$score_layer_edges,
                           eligible$pruning_steps, eligible$n_layers,
                           -eligible$primary, -eligible$secondary), , drop = FALSE]
utils::write.csv(eligible, file.path(out, "eligible_score_layer_solutions.csv"), row.names = FALSE)

# Write each compact, score-supported solution, including direct regulatory
# provenance from the relaxed mapping ledger. No one row is labelled "optimal".
solutions_dir <- file.path(out, "solutions")
dir.create(solutions_dir, recursive = TRUE)
for (key in names(materialized)) {
  item <- materialized[[key]]
  dir <- file.path(solutions_dir, key)
  dir.create(dir, recursive = TRUE)
  utils::write.csv(item$solution$SIF, file.path(dir, "network_SIF.csv"), row.names = FALSE)
  att <- item$solution$ATT
  lookup <- item$status[match(att$source, item$status$pkn_source_node), , drop = FALSE]
  for (column in intersect(c("assignment_id", "measurement_feature", "candidate_measurement_value",
                             "factor1_abs_weight", "factor1_abs_rank_desc", "factor1_abs_percentile",
                             "pkn_source_chebi", "pkn_source_name", "primary_mapping_policy",
                             "edge_fingerprint_group", "assignment_sign_agrees"), names(lookup))) {
    att[[column]] <- lookup[[column]]
  }
  utils::write.csv(att, file.path(dir, "network_ATT.csv"), row.names = FALSE)
  direct <- item$solution$SIF[item$solution$SIF$source %in%
                                item$solution$ATT$source[item$solution$ATT$type == "upstream_input"],
                              c("source", "target", "interaction", "consistency"), drop = FALSE]
  direct$.key <- paste(direct$source, direct$target, direct$interaction, sep = "\r")
  idx <- match(direct$.key, evidence$.key)
  if (anyNA(idx)) stop("A materialized score-layer direct edge lacks source evidence")
  direct <- cbind(direct[, setdiff(names(direct), ".key"), drop = FALSE],
                  evidence[idx, setdiff(names(evidence), ".key"), drop = FALSE])
  utils::write.csv(direct, file.path(dir, "direct_regulatory_evidence.csv"), row.names = FALSE)
  write_json(list(
    run_id = item$config$run_id,
    tf_panel = item$config$tf_panel,
    pruning_steps = item$config$pruning_steps,
    n_layers = item$config$n_layers,
    primary = item$config$primary,
    secondary = item$config$secondary,
    score_layer_constrained = TRUE,
    raw_upstream_candidate_nodes = item$solution$audit$retained_upstream,
    independent_edge_fingerprint_groups = length(item$fingerprints),
    edge_fingerprint_groups = item$fingerprints,
    observed_tfs = item$solution$audit$retained_observed,
    edges = nrow(item$solution$SIF),
    nodes = nrow(item$solution$ATT),
    interpretation = "Layer-consecutive, sign-consistent recursive-MOON score paths. Candidate chemical identity and mapping remain conditional; this extraction does not establish ligand exposure or causality."
  ), file.path(dir, "selection.json"))
}

utils::write.csv(data.frame(path = c(
  "scripts/moon/ancillary/extract_relaxed_lipid_tf_score_layer_paths.R",
  file.path(analysis_dir, "dependency_manifest.csv"),
  file.path(analysis_dir, "threshold_grid.csv"), evidence_path
), md5 = unname(tools::md5sum(c(
  "scripts/moon/ancillary/extract_relaxed_lipid_tf_score_layer_paths.R",
  file.path(analysis_dir, "dependency_manifest.csv"),
  file.path(analysis_dir, "threshold_grid.csv"), evidence_path
))), stringsAsFactors = FALSE), file.path(out, "dependency_manifest.csv"), row.names = FALSE)
writeLines(c(
  "# Layer-constrained relaxed lipid-to-TF MOON extractions",
  "",
  "This is a derived extraction from the completed score-first relaxed lipid-to-TF fits. It does not refit the PKN, pruning, compression, RNA coherence, or ULM scores. Before extracting routes, it checks the original dependency manifest.",
  "",
  "Only sign-concordant feature-to-chemical assignments are eligible as source seeds. Each retained edge must descend exactly one recursive ULM scoring level and agree with the signed PKN and node scores. This distinguishes score-supporting routes from the broader sign-consistent connectivity screen in the parent analysis.",
  "",
  "`score_layer_threshold_grid.csv` retains every threshold result. A compact result is listed in `eligible_score_layer_solutions.csv` only if it has fewer than 500 edges, at least one unique direct-edge fingerprint group, at least two observed TFs, and a nonzero secondary threshold. It is not a statistical model selection or proof of mechanism. Each materialized solution retains raw candidate-node count and fingerprint-group count so alternative mappings with identical direct PKN profiles are not counted as independent regulators.",
  "",
  "The source measurement may be copied to several chemically plausible nodes. Mapping tier, Factor 1 loading/rank, and direct-edge provenance are saved with each materialized network. Bulk abundance, ligand exposure, receptor engagement, temporal ordering and causal regulation remain unmeasured."
), file.path(out, "README.md"))
write_json(list(
  analysis = "factor1_lipid_to_tf_relaxed_score_layer_extraction",
  source_analysis = analysis_dir,
  threshold_evaluations = nrow(summary),
  nonempty_score_layer_solutions = sum(summary$score_layer_edges > 0L),
  compact_score_layer_solutions = nrow(eligible),
  materialized_solutions = length(materialized),
  score_layer_contract = "Every retained edge descends exactly one recursive ULM score level, is sign-consistent, and lies on a directed candidate-to-observed-TF walk.",
  eligibility_contract = "<500 edges; >=1 unique direct-edge fingerprint; >=2 observed TFs; secondary threshold >0."
), file.path(out, "analysis_summary.json"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed strict score-layer extraction: ", nrow(summary), " threshold evaluations; ",
    nrow(eligible), " compact materialized solution(s). Results: ", out, "\n", sep = "")
