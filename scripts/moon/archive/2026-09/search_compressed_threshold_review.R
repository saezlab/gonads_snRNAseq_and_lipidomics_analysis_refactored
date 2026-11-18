#!/usr/bin/env Rscript
# Extraction-only sensitivity search: all scoring states are immutable inputs.
options(stringsAsFactors = FALSE)
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
stopifnot(requireNamespace("igraph", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
# The extractor cleans the same immutable final PKN for every threshold pair.
# Memoize only that exact-input cleanup; extraction and graph audits are unchanged.
moon_check_pkn_uncached <- moon_check_pkn
pkn_cleanup_cache <- new.env(parent = emptyenv())
moon_check_pkn <- function(pkn) {
  if (exists("input", envir = pkn_cleanup_cache, inherits = FALSE) &&
      identical(pkn, pkn_cleanup_cache$input)) return(pkn_cleanup_cache$output)
  checked <- moon_check_pkn_uncached(pkn)
  pkn_cleanup_cache$input <- pkn
  pkn_cleanup_cache$output <- checked
  checked
}
sensitivity <- "--sensitivity" %in% commandArgs(trailingOnly = TRUE)
out <- if (sensitivity) "results/cosmos_factor1_compressed_review/sensitivity_threshold_search" else
  "results/cosmos_factor1_compressed_review/threshold_search"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, file) utils::write.csv(x, file.path(out, file), row.names = FALSE, na = "")
edge_key <- function(e) paste(e$source, e$target, e$interaction, sep = "\t")
grid_key <- function(x) paste(x$primary, x$secondary, sep = ":")
depths <- c(4L, 5L, 6L, 7L, 8L, 10L, 12L)
paths <- if (sensitivity) unlist(lapply(c("inputs_only_sensitivity", "topology_sensitivity", "prescore_sensitivity"),
  function(z) list.files(file.path("results/cosmos_factor1_compressed_review", z, "runs"),
    pattern = "^run_state\\.rds$", full.names = TRUE, recursive = TRUE)), use.names = FALSE) else
  sprintf("results/cosmos_factor1_parameter_search/runs/radius6_depth%d_compressed/run_state.rds", depths)
states <- lapply(paths, readRDS)
names(states) <- vapply(states, `[[`, character(1), "run_id")
depths <- vapply(states, `[[`, integer(1), "n_layers")
stopifnot(all(vapply(states, function(s) isTRUE(s$converged) && s$compression == "compressed", logical(1))))
same_state <- function(a, b) {
  identical(a$scores, b$scores) && identical(a$pkn_final, b$pkn_final) &&
    identical(a$tf, b$tf) && identical(a$lipid, b$lipid)
}
aliases <- character(length(states))
for (i in seq_along(states)) {
  previous <- if (i == 1L) integer() else which(vapply(states[seq_len(i - 1L)], same_state, logical(1), b = states[[i]]))
  aliases[i] <- if (length(previous)) aliases[previous[1L]] else names(states)[i]
}
alias_table <- data.frame(run_id = names(states), depth = depths, canonical_run = aliases,
  realized_max_level = vapply(states, `[[`, numeric(1), "realized_max_level"))
write_csv(alias_table, "score_state_aliases.csv")
active <- states[unique(aliases)]
catalog <- utils::read.csv("results/cosmos_mouse_pkn/nodes.csv", check.names = FALSE)
labels <- stats::setNames(catalog$label, catalog$node_id)
label <- function(x) { z <- unname(labels[x]); z[is.na(z) | !nzchar(z)] <- x[is.na(z) | !nzchar(z)]; z }
input <- readRDS("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
validation <- list()
for (s in active) {
  mapping_path <- if (sensitivity) file.path(s$run_dir, "compression.rds") else
    "results/cosmos_factor1_parameter_search/topology_radius_6/compressed/compression.rds"
  mapping <- readRDS(mapping_path)$node_map
  map <- stats::setNames(mapping$compressed, mapping$original)
  original <- s$pkn_final
  projected <- original
  projected$source <- unname(map[projected$source])
  projected$target <- unname(map[projected$target])
  scoring <- readRDS(file.path(s$run_dir, "pkn_final_scoring.rds"))
  pass <- !anyNA(projected) && setequal(edge_key(unique(projected)), edge_key(scoring)) &&
    !anyDuplicated(s$scores$source) && all(is.finite(s$scores$score)) &&
    all(s$scores$compressed_source == unname(map[s$scores$source])) &&
    all(vapply(split(s$scores$score, s$scores$compressed_source), function(v) length(unique(v)) == 1L, logical(1)))
  stopifnot(pass)
  validation[[s$run_id]] <- data.frame(run_id = s$run_id, original_final_edges = nrow(original),
    scoring_final_edges = nrow(scoring), merged_original_nodes = sum(mapping$merged),
    same_compressed_group_same_score = TRUE, final_original_projects_to_scored_pkn = TRUE)
}
write_csv(do.call(rbind, validation), "state_validation.csv")

extract_grid <- function(s, thresholds, phase) {
  net_keys <- edge_key(s$pkn_final)
  score <- stats::setNames(s$scores$score, s$scores$source)
  rows <- vector("list", nrow(thresholds))
  cat(s$run_id, phase, nrow(thresholds), "threshold pairs\n")
  for (j in seq_len(nrow(thresholds))) {
    p <- thresholds$primary[j]; z <- thresholds$secondary[j]
    sol <- moon_extract_paths(s$scores, s$pkn_final, s$tf,
      primary_thresh = p, secondary_thresh = z)
    e <- sol$SIF; a <- sol$ATT
    valid <- all(edge_key(e) %in% net_keys) &&
      all(sign(score[e$source] * score[e$target]) == e$interaction) &&
      all(a$source %in% s$scores$source) && !any(e$source %in% names(s$lipid)) &&
      isTRUE(sol$audit$all_nodes_on_seed_to_observation_walk)
    stopifnot(valid)
    tf <- a$source[a$type == "upstream_input"]
    lipid <- a$source[a$level == 0]
    route_pairs <- 0L; min_route <- NA_real_; max_route <- NA_real_
    if (nrow(e) > 0 && nrow(e) < 700) {
      g <- igraph::graph_from_data_frame(e[, c("source", "target")], directed = TRUE)
      distance <- igraph::distances(g, v = tf, to = lipid, mode = "out")
      route_pairs <- sum(is.finite(distance))
      min_route <- min(distance[is.finite(distance)])
      max_route <- max(distance[is.finite(distance)])
    }
    rows[[j]] <- data.frame(run_id = s$run_id, depth = s$n_layers, phase = phase,
      primary = p, secondary = z, nodes = nrow(a), edges = nrow(e),
      upstream_tfs = length(tf), downstream_lipids = length(lipid),
      tf_names = paste(sort(label(tf)), collapse = ";"),
      lipid_names = paste(sort(label(lipid)), collapse = ";"),
      route_pairs = route_pairs, shortest_route_min = min_route, shortest_route_max = max_route,
      converged = s$converged, valid = valid,
      candidate_edges = sol$audit$candidate_edges, consistent_edges = sol$audit$consistent_edges)
    if (j %% 100L == 0L) cat(" ", j, "/", nrow(thresholds), "\n")
  }
  do.call(rbind, rows)
}
pareto <- function(x, objectives) {
  if (!nrow(x)) return(logical())
  values <- as.matrix(x[, objectives, drop = FALSE])
  !vapply(seq_len(nrow(values)), function(i) {
    at_least <- rowSums(sweep(values, 2, values[i, ], `>=`)) == ncol(values)
    strictly <- rowSums(sweep(values, 2, values[i, ], `>`)) > 0
    any(at_least & strictly)
  }, logical(1))
}
threshold_grid <- function(p, s) {
  x <- expand.grid(primary = round(p, 6), secondary = round(s, 6), KEEP.OUT.ATTRS = FALSE)
  x <- unique(x[x$secondary <= x$primary & x$primary > 0 & x$secondary >= 0, ])
  x[order(-x$secondary, -x$primary), ]
}
coarse <- threshold_grid(seq(0.5, 3, 0.1), seq(0, 2, 0.1))
all_rows <- list()
for (s in active) {
  cache <- file.path(out, paste0(s$run_id, "_coarse.csv"))
  rows <- if (file.exists(cache)) utils::read.csv(cache, check.names = FALSE) else extract_grid(s, coarse, "coarse")
  utils::write.csv(rows, cache, row.names = FALSE, na = "")
  all_rows[[s$run_id]] <- rows
}
grid <- do.call(rbind, all_rows)
write_csv(grid, "grid_coarse.csv")
for (s in active) {
  current <- grid[grid$run_id == s$run_id, ]
  eligible <- current[current$edges > 0 & current$edges < 500 & current$upstream_tfs >= 2 & current$downstream_lipids >= 2, ]
  if (!nrow(eligible)) eligible <- current[current$edges > 0 & current$edges < 500, ]
  frontier <- if (nrow(eligible)) eligible[pareto(eligible, c("primary", "secondary", "upstream_tfs", "downstream_lipids")), ] else eligible
  # Include near-cap infeasible points so coarse jumps across 500 are refined.
  near <- current[current$edges >= 500 & current$edges < 1000 & current$upstream_tfs >= 2 & current$downstream_lipids >= 2, ]
  if (nrow(near)) near <- near[pareto(near, c("primary", "secondary")), ]
  anchors <- unique(rbind(frontier, near))
  if (!nrow(anchors)) next
  local <- do.call(rbind, lapply(seq_len(nrow(anchors)), function(i) {
    threshold_grid(seq(max(0.1, anchors$primary[i] - 0.1), anchors$primary[i] + 0.1, by = 0.025),
      seq(max(0, anchors$secondary[i] - 0.1), anchors$secondary[i] + 0.1, by = 0.025))
  }))
  local <- unique(local)
  local <- local[!grid_key(local) %in% grid_key(current), ]
  cache <- file.path(out, paste0(s$run_id, "_refined.csv"))
  rows <- if (file.exists(cache)) utils::read.csv(cache, check.names = FALSE) else extract_grid(s, local, "refined")
  utils::write.csv(rows, cache, row.names = FALSE, na = "")
  grid <- rbind(grid, rows)
  write_csv(grid, "grid_complete.csv")
}
grid$eligible <- grid$valid & grid$converged & grid$edges > 0 & grid$edges < 500
grid$multiple_tfs_and_lipids <- grid$eligible & grid$upstream_tfs >= 2 & grid$downstream_lipids >= 2
write_csv(grid, "grid_complete.csv")
eligible <- grid[grid$eligible, ]
eligible$coverage_threshold_pareto <- pareto(eligible, c("primary", "secondary", "upstream_tfs", "downstream_lipids"))
eligible$threshold_frontier <- pareto(eligible, c("primary", "secondary"))
eligible <- eligible[order(-eligible$secondary, -eligible$primary, -eligible$downstream_lipids, -eligible$upstream_tfs, eligible$edges), ]
write_csv(eligible, "eligible_under_500.csv")
write_csv(eligible[eligible$coverage_threshold_pareto, ], "coverage_threshold_pareto.csv")
write_csv(eligible[eligible$threshold_frontier, ], "threshold_frontier.csv")
interpretable <- eligible[eligible$multiple_tfs_and_lipids, ]
if (nrow(interpretable)) {
  picks <- list()
  for (run in unique(interpretable$run_id)) {
    q <- interpretable[interpretable$run_id == run, ]
    strict <- q[order(-q$secondary, -q$primary, -q$downstream_lipids, -q$upstream_tfs, q$edges), ][1, ]
    broad <- q[order(-q$downstream_lipids, -q$upstream_tfs, -q$secondary, -q$primary, q$edges), ][1, ]
    picks[[paste0(run, "_strict")]] <- strict
    picks[[paste0(run, "_coverage")]] <- broad
    if (!sensitivity) {
      balanced <- q[abs(q$primary - 1.85) < 1e-8 & abs(q$secondary - 0.55) < 1e-8, ]
      if (nrow(balanced)) picks[[paste0(run, "_balanced")]] <- balanced[1, ]
    }
  }
  selected <- unique(do.call(rbind, picks))
  selected$candidate_id <- sprintf("C%02d", seq_len(nrow(selected)))
  selected$solution_dir <- NA_character_
  endpoint_rows <- list(); path_rows <- list(); checks <- list()
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]; s <- active[[row$run_id]]
    # Independently materialize every exported network with the unmodified,
    # uncached helper contract, and verify its sizes against the cached grid.
    moon_check_pkn <- moon_check_pkn_uncached
    sol <- moon_extract_paths(s$scores, s$pkn_final, s$tf,
      primary_thresh = row$primary, secondary_thresh = row$secondary,
      RNA_input = input$hepatocyte_gene_weight)
    dir <- file.path(out, row$candidate_id)
    dir.create(dir, showWarnings = FALSE)
    selected$solution_dir[i] <- dir
    a <- sol$ATT; a$label <- label(a$source)
    a$compressed_source <- s$scores$compressed_source[match(a$source, s$scores$source)]
    a$upstream_activity <- unname(s$tf[a$source])
    a$secondary_margin <- ifelse(a$level > 0, abs(a$score) - row$secondary, NA_real_)
    a$above_primary <- a$level > 0 & abs(a$score) > row$primary
    e <- sol$SIF; e$source_label <- label(e$source); e$target_label <- label(e$target)
    utils::write.csv(e, file.path(dir, "network_SIF.csv"), row.names = FALSE, na = "")
    utils::write.csv(a, file.path(dir, "network_ATT.csv"), row.names = FALSE, na = "")
    saveRDS(sol, file.path(dir, "solution.rds"))
    jsonlite::write_json(c(list(run_id = s$run_id, candidate_id = row$candidate_id), sol$audit),
      file.path(dir, "extraction_audit.json"), pretty = TRUE, auto_unbox = TRUE)
    vertex <- a; names(vertex)[names(vertex) == "source"] <- "name"
    g <- igraph::graph_from_data_frame(e, directed = TRUE, vertices = vertex)
    igraph::write_graph(g, file.path(dir, "network.graphml"), format = "graphml")
    tf <- a$source[a$type == "upstream_input"]; lipid <- a$source[a$level == 0]
    endpoints <- a[a$source %in% c(tf, lipid), ]
    endpoints$candidate_id <- row$candidate_id
    endpoint_rows[[i]] <- endpoints
    routes <- list()
    for (t in tf) for (l in lipid) {
      if (!is.finite(igraph::distances(g, v = t, to = l, mode = "out")[1])) next
      ids <- igraph::as_ids(igraph::shortest_paths(g, from = t, to = l, mode = "out")$vpath[[1]])
      routes[[length(routes) + 1L]] <- data.frame(candidate_id = row$candidate_id, tf = label(t),
        lipid = label(l), edge_count = length(ids) - 1L,
        node_ids = paste(ids, collapse = " -> "), labels = paste(label(ids), collapse = " -> "))
    }
    path_rows[[i]] <- do.call(rbind, routes)
    checks[[i]] <- data.frame(candidate_id = row$candidate_id, edges = nrow(e), nodes = nrow(a),
      matches_grid = nrow(e) == row$edges && nrow(a) == row$nodes,
      edges_in_final_pkn = all(edge_key(e) %in% edge_key(s$pkn_final)),
      score_exactly_matches_state = identical(unname(a$score), unname(s$scores$score[match(a$source, s$scores$source)])),
      all_inferred_above_secondary = all(a$secondary_margin[a$level > 0] > 0),
      all_edges_incident_to_primary = all(e$source %in% a$source[a$above_primary] | e$target %in% a$source[a$above_primary]),
      all_edges_sign_consistent = all(e$consistency),
      all_nodes_on_seed_observation_walk = sol$audit$all_nodes_on_seed_to_observation_walk)
  }
  write_csv(selected, "selected_candidates.csv")
  write_csv(do.call(rbind, endpoint_rows), "selected_endpoints_and_margins.csv")
  write_csv(do.call(rbind, path_rows), "selected_shortest_routes.csv")
  check <- do.call(rbind, checks)
  stopifnot(all(as.matrix(check[, -(1:3)])))
  write_csv(check, "selected_validation.csv")
}
jsonlite::write_json(list(created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  mode = if (sensitivity) "Extraction only; new compressed topology, pre-score RNA, and input-only protection sensitivity states" else
    "Extraction only; unchanged converged compressed states from the previous search",
  cleanup_cache = "Exact-input memoization of repeated PKN cleanup; materialized candidates re-extracted without the cache",
  unique_scoring_states = length(active), aliases = alias_table,
  coarse_primary = seq(0.5, 3, 0.1), coarse_secondary = seq(0, 2, 0.1),
  refined = "0.025 increments within +/-0.1 of threshold/coverage Pareto and near-cap 500-999-edge coarse anchors",
  retention = "All grid rows, including empty; strict edges <500; at least 2 TF and 2 lipid observations preferred",
  level0_policy = "Original measured weights are exempt from inferred-node score cutoffs and cannot serve as primary nodes",
  md5 = as.list(tools::md5sum(c(paths, "scripts/moon/search_compressed_threshold_review.R", "scripts/moon/moon_safe_helpers.R")))),
  file.path(out, "search_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed", nrow(grid), "unique-state threshold pairs;", nrow(eligible), "eligible <500 edges;", nrow(interpretable), "with >=2 TFs and >=2 lipids.\n")
