#!/usr/bin/env Rscript
# Topology-only depth sensitivity. Does not run MOON, coherence, or compression.
# Run from repository root: Rscript scripts/moon/audit_factor1_depth_topology.R
source("scripts/support_functions.R")
stopifnot(requireNamespace("igraph", quietly = TRUE), requireNamespace("cosmosR", quietly = TRUE))
out <- "results/cosmos_factor1_parameter_search"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input_file <- "results/cosmos_factor1_pilot/inputs/factor1_inputs.rds"
pkn_file <- "results/cosmos_mouse_pkn/pkn_lipid_rules.rds"
x <- readRDS(input_file)
pkn <- cosmosR::meta_network_cleanup(readRDS(pkn_file))[, c("source", "interaction", "target")]
node_ids <- function(e) unique(c(e$source, e$target))
key <- function(e) paste(e$source, e$interaction, e$target, sep = "\t")
tf <- names(x$tf_activity)[abs(x$tf_activity) > 2 & names(x$tf_activity) %in% node_ids(pkn)]
lipid <- names(x$lipid_weight)[x$lipid_weight != 0 & names(x$lipid_weight) %in% node_ids(pkn)]
terminal <- pkn[!pkn$source %in% lipid, ]
reach <- function(e, seeds, mode = "out", depth) {
  if (mode == "in") e <- data.frame(source = e$target, target = e$source)
  seen <- intersect(seeds, node_ids(e)); frontier <- seen; step <- 0L
  while (length(frontier) && step < depth) {
    new <- setdiff(unique(e$target[e$source %in% frontier]), seen)
    seen <- c(seen, new); frontier <- new; step <- step + 1L
  }
  seen
}
graph_of <- function(e) igraph::graph_from_data_frame(e[, c("source", "target")], directed = TRUE)
distances <- function(g, from, to = igraph::V(g)$name, mode = "out") {
  igraph::distances(g, v = from, to = to, mode = mode, algorithm = "unweighted")
}
full_g <- graph_of(terminal)
present <- igraph::V(full_g)$name
full_pair <- matrix(Inf, nrow = length(tf), ncol = length(lipid), dimnames = list(tf, lipid))
full_pair[intersect(tf, present), intersect(lipid, present)] <- distances(full_g, intersect(tf, present), intersect(lipid, present))
full_reachable_pairs <- sum(is.finite(full_pair))
summary <- list(); histograms <- list(); pair_rows <- list(); endpoint_rows <- list()
previous_keys <- character(); previous_pair <- NULL
for (depth in 4:14) {
  net <- terminal; retained_tf <- tf; retained_lipid <- lipid; iteration <- 0L
  repeat {
    old <- nrow(net); iteration <- iteration + 1L
    forward <- reach(net, retained_tf, depth = depth)
    net <- net[net$source %in% forward & net$target %in% forward, ]
    backward <- reach(net, retained_lipid, mode = "in", depth = depth)
    net <- net[net$source %in% backward & net$target %in% backward, ]
    retained_tf <- intersect(retained_tf, node_ids(net))
    retained_lipid <- intersect(retained_lipid, node_ids(net))
    if (nrow(net) == old) break
    stopifnot(iteration < 50)
  }
  stopifnot(nrow(net) > 0, !any(net$source %in% lipid))
  g <- graph_of(net); vertices <- igraph::V(g)$name
  from_tf <- distances(g, retained_tf)
  to_lipid <- distances(g, retained_lipid, mode = "in")
  nearest_from_tf <- apply(from_tf, 2, min)
  nearest_to_lipid <- apply(to_lipid, 2, min)
  pair <- matrix(Inf, nrow = length(tf), ncol = length(lipid), dimnames = list(tf, lipid))
  pair[retained_tf, retained_lipid] <- from_tf[, retained_lipid, drop = FALSE]
  stopifnot(all(pair >= full_pair))
  finite <- as.numeric(pair[is.finite(pair)])
  # These unions capture the shortest route from any seed to each lipid, plus
  # the shortest route from each TF to any lipid. They are not all pairwise paths.
  endpoint_minimal_edge <- rep(FALSE, nrow(net))
  for (j in seq_along(retained_lipid)) {
    shortest <- min(pair[, retained_lipid[j]])
    endpoint_minimal_edge <- endpoint_minimal_edge |
      (nearest_from_tf[net$source] + 1 + to_lipid[j, net$target] == shortest)
  }
  for (j in seq_along(retained_tf)) {
    shortest <- min(pair[retained_tf[j], ])
    endpoint_minimal_edge <- endpoint_minimal_edge |
      (from_tf[j, net$source] + 1 + nearest_to_lipid[net$target] == shortest)
  }
  components <- igraph::components(g, mode = "strong")
  membership <- components$membership
  cyclic_node <- components$csize[membership] > 1
  cyclic_edge <- membership[net$source] == membership[net$target]
  added_edge <- !key(net) %in% previous_keys
  if (!is.null(previous_pair)) {
    stopifnot(all(pair <= previous_pair))
    new_pairs <- sum(is.finite(pair) & !is.finite(previous_pair))
    shorter_pairs <- sum(is.finite(previous_pair) & pair < previous_pair)
  } else {
    new_pairs <- NA_integer_; shorter_pairs <- NA_integer_
  }
  summary[[length(summary) + 1L]] <- data.frame(depth = depth,
    pruning_iterations = iteration, retained_edges = nrow(net), retained_nodes = length(vertices),
    candidate_tfs = length(tf), retained_tfs = length(retained_tf), candidate_lipids = length(lipid),
    retained_lipids = length(retained_lipid), reachable_tf_lipid_pairs = length(finite),
    full_terminal_reachable_pairs = full_reachable_pairs,
    pair_coverage_of_full_terminal = length(finite) / full_reachable_pairs,
    pairs_shortest_length_le_depth = sum(finite <= depth),
    pairs_shortest_length_gt_depth = sum(finite > depth),
    shortest_pair_min = min(finite), shortest_pair_median = median(finite),
    shortest_pair_p90 = as.numeric(quantile(finite, 0.9)), shortest_pair_max = max(finite),
    new_reachable_pairs_vs_previous_depth = new_pairs,
    shortened_pairs_vs_previous_depth = shorter_pairs,
    endpoint_minimal_route_edges = sum(endpoint_minimal_edge),
    edges_on_directed_cycles = sum(cyclic_edge), nodes_on_directed_cycles = sum(cyclic_node),
    largest_strong_component = max(components$csize),
    added_edges_vs_previous_depth = if (depth == 4) NA_integer_ else sum(added_edge),
    added_edges_on_cycles = if (depth == 4) NA_integer_ else sum(added_edge & cyclic_edge),
    added_endpoint_minimal_route_edges = if (depth == 4) NA_integer_ else sum(added_edge & endpoint_minimal_edge))
  hist <- as.data.frame(table(factor(as.vector(pair), levels = c(1:28, Inf))))
  names(hist) <- c("shortest_path_edges", "tf_lipid_pairs")
  hist$depth <- depth
  histograms[[length(histograms) + 1L]] <- hist[, c("depth", "shortest_path_edges", "tf_lipid_pairs")]
  endpoints <- rbind(data.frame(depth = depth, node_id = tf, type = "TF",
                retained = tf %in% retained_tf, shortest_to_opposite_input = apply(pair, 1, min)),
             data.frame(depth = depth, node_id = lipid, type = "lipid",
                retained = lipid %in% retained_lipid, shortest_to_opposite_input = apply(pair, 2, min)))
  endpoint_rows[[length(endpoint_rows) + 1L]] <- endpoints
  pairs <- as.data.frame(as.table(pair)); names(pairs) <- c("tf_node", "lipid_node", "shortest_path_edges")
  pairs$depth <- depth
  pair_rows[[length(pair_rows) + 1L]] <- pairs[, c("depth", "tf_node", "lipid_node", "shortest_path_edges")]
  cat("Depth", depth, ":", nrow(net), "edges,", length(vertices), "nodes,",
      length(retained_tf), "TFs,", length(retained_lipid), "lipids,", length(finite), "reachable pairs;",
      "max shortest path", max(finite), "; cycle edges", sum(cyclic_edge), "\n")
  previous_keys <- key(net); previous_pair <- pair
  # Validate the existing depth-six topology exactly, without writing over it.
  if (depth == 6 && file.exists("results/cosmos_factor1_pilot/pkn_pruned_uncompressed.rds")) {
    old_net <- readRDS("results/cosmos_factor1_pilot/pkn_pruned_uncompressed.rds")
    stopifnot(setequal(key(net), key(old_net)))
  }
}
write.csv(do.call(rbind, summary), file.path(out, "depth_topology_audit.csv"), row.names = FALSE)
write.csv(do.call(rbind, histograms), file.path(out, "depth_shortest_path_distribution.csv"), row.names = FALSE)
write.csv(do.call(rbind, endpoint_rows), file.path(out, "depth_input_coverage.csv"), row.names = FALSE)
write.csv(do.call(rbind, pair_rows), file.path(out, "depth_tf_lipid_shortest_paths.csv"), row.names = FALSE)
sources <- c(input_file, pkn_file, "scripts/moon/audit_factor1_depth_topology.R")
write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, "depth_topology_source_manifest.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "depth_topology_session_info.txt"))
print(do.call(rbind, summary), row.names = FALSE)
