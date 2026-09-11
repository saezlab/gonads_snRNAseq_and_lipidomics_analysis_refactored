#!/usr/bin/env Rscript
# Read-only scoring audit: independently validate saved candidate SIF/ATT exports.
source("scripts/support_functions.R")
base <- "results/cosmos_factor1_parameter_search"
grid <- utils::read.csv(file.path(base, "parameter_grid.csv"))
write_audit <- function(x, name) utils::write.csv(x, file.path(base, name), row.names = FALSE, na = "")
edge_keys <- function(e) paste(e$source, e$target, e$interaction, sep = "\t")
nodes <- function(e) unique(c(e$source, e$target))
walk <- function(e, seeds, reverse = FALSE) {
  src <- if (reverse) e$target else e$source
  trg <- if (reverse) e$source else e$target
  seen <- intersect(seeds, unique(c(src, trg))); frontier <- seen
  while (length(frontier)) {
    new <- setdiff(unique(trg[src %in% frontier]), seen)
    seen <- union(seen, new); frontier <- new
  }
  seen
}
# This reproduces the declared set operation independently, without invoking
# moon_extract_paths or any MOON scoring/coherence/compression function.
expected_edges <- function(st, primary, secondary) {
  s <- st$scores; e <- st$pkn_final
  sm <- stats::setNames(s$score, s$source)
  observed <- s$source[s$level == 0 & s$score != 0]
  prim <- s$source[s$level > 0 & abs(s$score) > primary]
  eligible <- union(observed, s$source[s$level > 0 & abs(s$score) > secondary])
  seeds <- intersect(names(st$tf), eligible)
  seeds <- seeds[sign(sm[seeds]) == sign(st$tf[seeds])]
  eligible <- setdiff(eligible, setdiff(intersect(names(st$tf), s$source), seeds))
  prim <- intersect(prim, eligible)
  e <- e[e$source %in% eligible & e$target %in% eligible &
           (e$source %in% prim | e$target %in% prim) & !e$source %in% observed, , drop = FALSE]
  e <- e[sign(sm[e$source] * sm[e$target]) == e$interaction, , drop = FALSE]
  keep <- intersect(walk(e, seeds), walk(e, observed, reverse = TRUE))
  e[e$source %in% keep & e$target %in% keep, , drop = FALSE]
}
files <- list.files(file.path(base, "runs"), pattern = "^solution_.*_SIF\\.csv$",
                    recursive = TRUE, full.names = TRUE)
files <- files[grepl("/radius6_", files)]
stopifnot(length(files) > 0)
checks <- list(); routes <- list(); endpoints <- list()
for (path in files) {
  run_dir <- dirname(path); st <- readRDS(file.path(run_dir, "run_state.rds"))
  stem <- sub("_SIF\\.csv$", "", basename(path))
  bits <- regmatches(stem, regexec("^solution_p([0-9]+_[0-9]+)_s([0-9]+_[0-9]+)$", stem))[[1]]
  primary <- as.numeric(sub("_", ".", bits[2])); secondary <- as.numeric(sub("_", ".", bits[3]))
  e <- utils::read.csv(path); a <- utils::read.csv(sub("_SIF\\.csv$", "_ATT.csv", path))
  sm <- stats::setNames(st$scores$score, st$scores$source)
  lm <- stats::setNames(st$scores$level, st$scores$source)
  measured <- a$source[a$level == 0]; seeds <- intersect(a$source, names(st$tf))
  inferred <- a$source[a$level > 0]
  strong <- inferred[abs(sm[inferred]) > primary]
  exp <- expected_edges(st, primary, secondary)
  manifest <- jsonlite::read_json(file.path(run_dir, "run_manifest.json"), simplifyVector = TRUE)
  ledger <- utils::read.csv(file.path(run_dir, "coherence_ledger.csv"))
  gi <- grid$run_id == st$run_id & abs(grid$primary - primary) < 1e-9 &
    abs(grid$secondary - secondary) < 1e-9
  gg <- igraph::graph_from_data_frame(e[, c("source", "target")], directed = TRUE)
  dd <- igraph::distances(gg, v = seeds, to = measured, mode = "out")
  rec <- data.frame(run_id = st$run_id, solution = stem, depth = st$n_layers,
    compression = st$compression, primary = primary, secondary = secondary,
    nodes = nrow(a), edges = nrow(e), upstream_tfs = length(seeds),
    measured_lipids = length(measured), strict_under_500 = nrow(e) > 0 && nrow(e) < 500,
    unique_edges_and_nodes = !anyDuplicated(edge_keys(e)) && !anyDuplicated(a$source),
    nodes_equal_edge_endpoints = setequal(a$source, nodes(e)),
    edges_in_final_pkn = all(edge_keys(e) %in% edge_keys(st$pkn_final)),
    scores_match_saved_state = all(a$source %in% st$scores$source) &&
      max(abs(a$score - sm[a$source])) < 1e-12 && all(a$level == lm[a$source]),
    all_edges_sign_consistent = all(sign(sm[e$source] * sm[e$target]) == e$interaction),
    all_inferred_pass_secondary = all(abs(sm[inferred]) > secondary),
    all_edges_have_primary_inferred_endpoint = all(e$source %in% strong | e$target %in% strong),
    observed_weights_match_inputs = all(measured %in% names(st$lipid)) &&
      max(abs(sm[measured] - st$lipid[measured])) < 1e-12,
    measured_nodes_are_sinks = !any(e$source %in% measured),
    intended_tf_signs_match = length(seeds) > 0 && all(sign(sm[seeds]) == sign(st$tf[seeds])),
    every_node_reached_from_tf = all(nodes(e) %in% walk(e, seeds)),
    every_node_can_reach_observation = all(nodes(e) %in% walk(e, measured, reverse = TRUE)),
    complete_declared_extraction = setequal(edge_keys(e), edge_keys(exp)),
    convergence_confirmed = isTRUE(st$converged) && isTRUE(manifest$converged) &&
      tail(ledger$edges_removed, 1) == 0,
    table_counts_match = any(gi) && all(grid$edges[gi] == nrow(e)) &&
      all(grid$nodes[gi] == nrow(a)) && all(grid$upstream_tfs[gi] == length(seeds)) &&
      all(grid$downstream_lipids[gi] == length(measured)),
    reachable_tf_lipid_pairs = sum(is.finite(dd)), possible_tf_lipid_pairs = length(dd),
    min_shortest_route_edges = min(dd[is.finite(dd)]),
    max_shortest_route_edges = max(dd[is.finite(dd)]),
    min_inferred_secondary_margin = min(abs(sm[inferred]) - secondary),
    min_tf_secondary_margin = min(abs(sm[seeds]) - secondary))
  logical_cols <- vapply(rec, is.logical, logical(1))
  rec$all_checks_pass <- all(unlist(rec[, logical_cols, drop = FALSE]))
  checks[[length(checks) + 1L]] <- rec
  for (s in seeds) for (m in measured) routes[[length(routes) + 1L]] <- data.frame(
    run_id = st$run_id, solution = stem, tf = s, lipid = m,
    shortest_path_edges = unname(dd[s, m]), tf_score = unname(sm[s]),
    tf_input_activity = unname(st$tf[s]), tf_secondary_margin = abs(sm[s]) - secondary)
  ep <- e[e$target %in% measured, , drop = FALSE]
  ep$run_id <- st$run_id; ep$solution <- stem
  ep$source_score <- unname(sm[ep$source]); ep$lipid_weight <- unname(st$lipid[ep$target])
  endpoints[[length(endpoints) + 1L]] <- ep
}
validation <- do.call(rbind, checks)
write_audit(validation, "refined_candidate_validation.csv")
write_audit(do.call(rbind, routes), "refined_candidate_routes.csv")
write_audit(do.call(rbind, endpoints), "refined_candidate_endpoint_edges.csv")
stopifnot(all(validation$all_checks_pass))

# Compare the two refined uncompressed representatives, then remove the
# unequal-primary-cutoff confounding with a matched-threshold comparison.
r4 <- readRDS(file.path(base, "runs/radius6_depth4_uncompressed/run_state.rds"))
r5 <- readRDS(file.path(base, "runs/radius6_depth5_uncompressed/run_state.rds"))
stopifnot(identical(r4$tf, r5$tf), identical(r4$lipid, r5$lipid),
          r4$topology_radius == r5$topology_radius)
shared_scores <- merge(r4$scores[, c("source", "score", "level")],
                       r5$scores[, c("source", "score", "level")],
                       by = "source", suffixes = c("_4", "_5"))
write_audit(data.frame(depth4_scored_nodes = nrow(r4$scores),
                      depth5_scored_nodes = nrow(r5$scores),
                      shared_scored_nodes = nrow(shared_scores),
                      shared_score_changes = sum(shared_scores$score_4 != shared_scores$score_5),
                      max_shared_score_difference = max(abs(shared_scores$score_4 - shared_scores$score_5)),
                      additional_level5_nodes = sum(r5$scores$level == 5)),
             "depth4_depth5_score_stability.csv")
comparisons <- list()
for (pair in list(c(1.9, 1.6, 2.05, 1.6), c(1.9, 1.6, 1.9, 1.6), c(2.05, 1.6, 2.05, 1.6))) {
  a <- expected_edges(r4, pair[1], pair[2]); b <- expected_edges(r5, pair[3], pair[4])
  ak <- edge_keys(a); bk <- edge_keys(b)
  comparisons[[length(comparisons) + 1L]] <- data.frame(
    depth4_primary = pair[1], depth4_secondary = pair[2],
    depth5_primary = pair[3], depth5_secondary = pair[4],
    depth4_edges = nrow(a), depth5_edges = nrow(b),
    shared_edges = length(intersect(ak, bk)), edge_union = length(union(ak, bk)),
    edge_jaccard = length(intersect(ak, bk)) / length(union(ak, bk)),
    depth4_tfs = paste(sort(intersect(names(r4$tf), nodes(a))), collapse = ";"),
    depth5_tfs = paste(sort(intersect(names(r5$tf), nodes(b))), collapse = ";"),
    depth4_lipids = paste(sort(intersect(names(r4$lipid), nodes(a))), collapse = ";"),
    depth5_lipids = paste(sort(intersect(names(r5$lipid), nodes(b))), collapse = ";"))
}
write_audit(do.call(rbind, comparisons), "depth4_depth5_refined_comparison.csv")
cat("All", nrow(validation), "saved candidates passed independent validation.\n")
print(validation[, c("run_id", "solution", "nodes", "edges", "upstream_tfs", "measured_lipids",
                       "max_shortest_route_edges", "min_tf_secondary_margin", "all_checks_pass")])
print(do.call(rbind, comparisons))
