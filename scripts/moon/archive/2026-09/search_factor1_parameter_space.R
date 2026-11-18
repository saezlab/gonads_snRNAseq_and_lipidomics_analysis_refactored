#!/usr/bin/env Rscript
# Systematic depth and dual-threshold search for the Factor 1 MOON pilot.
#
# The input mapping, TF cutoff and PKN are intentionally fixed.  This script
# varies only MOON propagation depth, safe same-child compression, and the
# post-score dual thresholds.  Topology pruning radius is kept separate from
# n_layers so that depth is not confounded with a changing input background.
options(stringsAsFactors = FALSE)

source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")

required <- c("cosmosR", "decoupleR", "igraph", "jsonlite")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(x) x %in% args
option_value <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", key, "="), "", hit[[length(hit)]])
}
parse_numbers <- function(text_value, as_integer = FALSE) {
  pieces <- trimws(unlist(strsplit(as.character(text_value[[1L]]), ",", fixed = TRUE), use.names = FALSE))
  value <- suppressWarnings(as.numeric(pieces))
  if (!length(value) || any(!is.finite(value)) || any(value <= 0)) {
    stop("Expected a nonempty comma-separated positive numeric option; got: ", text_value[[1L]])
  }
  if (as_integer) {
    stopifnot(all(value == as.integer(value)))
    value <- as.integer(value)
  }
  unique(value)
}

phase <- option_value("--phase", "coarse")
stopifnot(phase %in% c("coarse", "refine"))
depths <- parse_numbers(option_value("--depths", "4,5,6,7,8,10,12"), as_integer = TRUE)
topology_radius <- as.integer(option_value("--topology-radius", "6"))
stopifnot(is.finite(topology_radius), topology_radius >= 1L)
compression_modes <- strsplit(option_value("--compression", "uncompressed,compressed"), ",", fixed = TRUE)[[1]]
stopifnot(all(compression_modes %in% c("uncompressed", "compressed")))
force <- has_flag("--force")
out_root <- "results/cosmos_factor1_parameter_search"

if (has_flag("--validate-setup")) {
  cat("Parameter-search runner parsed. Packages:", paste(required, collapse = ", "), "\n")
  cat("Default depth grid:", paste(depths, collapse = ", "),
      "; fixed topology radius:", topology_radius,
      "; compression:", paste(compression_modes, collapse = ", "), "\n")
  quit(save = "no")
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "runs"), recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
node_ids <- function(e) unique(c(e$source, e$target))
edge_key <- function(e) paste(e$source, e$interaction, e$target, sep = "\t")

reach <- function(e, seeds, direction = "out", depth = Inf) {
  if (direction == "in") e <- data.frame(source = e$target, target = e$source)
  seen <- intersect(seeds, node_ids(e))
  frontier <- seen
  step <- 0L
  while (length(frontier) && step < depth) {
    new <- setdiff(unique(e$target[e$source %in% frontier]), seen)
    seen <- c(seen, new)
    frontier <- new
    step <- step + 1L
  }
  seen
}

validate_vector <- function(x) {
  stopifnot(is.numeric(x), length(x) > 0, !is.null(names(x)),
            !anyDuplicated(names(x)), all(nzchar(names(x))), all(is.finite(x)))
}

input_path <- "results/cosmos_factor1_pilot/inputs/factor1_inputs.rds"
pkn_path <- "results/cosmos_mouse_pkn/pkn_lipid_rules.rds"
input <- readRDS(input_path)
needed <- c("tf_activity", "lipid_weight", "hepatocyte_gene_weight", "tf_regulon")
stopifnot(all(needed %in% names(input)))
invisible(lapply(input[needed[1:3]], validate_vector))

pkn_raw <- cosmosR::meta_network_cleanup(readRDS(pkn_path))
pkn_raw <- pkn_raw[, c("source", "interaction", "target")]
stopifnot(!anyDuplicated(pkn_raw[, c("source", "target")]),
          all(pkn_raw$interaction %in% c(-1, 1)), all(pkn_raw$source != pkn_raw$target))

tf_all <- input$tf_activity
tf_seed <- tf_all[abs(tf_all) > 2 & names(tf_all) %in% node_ids(pkn_raw)]
lipid_seed <- input$lipid_weight[input$lipid_weight != 0 & names(input$lipid_weight) %in% node_ids(pkn_raw)]
rna <- input$hepatocyte_gene_weight
reg <- input$tf_regulon[, c("source", "target", "mor")]
stopifnot(length(tf_seed) > 0, length(lipid_seed) > 2,
          all(grepl("^Gene__", names(tf_seed))),
          all(grepl("^Metab__", names(lipid_seed))),
          all(reg$mor %in% c(-1, 1)), all(is.finite(reg$mor)))

# Keep observed lipid nodes terminal, exactly as moon() will do.  The topology
# radius is fixed within a sweep so the same PKN background is scored.  It is
# deliberately allowed to be lower than a tested n_layers value: this tests
# additional propagation in a frozen, auditable PKN rather than silently
# expanding the observation background with depth.
prepare_topology <- function(radius) {
  net <- pkn_raw[!pkn_raw$source %in% names(lipid_seed), ]
  tf <- tf_seed
  lipid <- lipid_seed
  ledger <- list()
  repeat {
    old <- nrow(net)
    forward <- reach(net, names(tf), depth = radius)
    net <- net[net$source %in% forward & net$target %in% forward, , drop = FALSE]
    backward <- reach(net, names(lipid), direction = "in", depth = radius)
    net <- net[net$source %in% backward & net$target %in% backward, , drop = FALSE]
    present <- node_ids(net)
    tf <- tf[names(tf) %in% present]
    lipid <- lipid[names(lipid) %in% present]
    ledger[[length(ledger) + 1L]] <- data.frame(
      iteration = length(ledger) + 1L, edges = nrow(net), nodes = length(present),
      upstream_tfs = length(tf), downstream_lipids = length(lipid)
    )
    if (old == nrow(net)) break
  }
  stopifnot(nrow(net) > 0, length(tf) > 0, length(lipid) > 2)
  list(net = moon_check_pkn(net), tf = tf, lipid = lipid, ledger = do.call(rbind, ledger))
}

prepared <- prepare_topology(topology_radius)
prepared_dir <- file.path(out_root, paste0("topology_radius_", topology_radius))
dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(prepared, file.path(prepared_dir, "prepared_topology.rds"))
write_csv(prepared$ledger, file.path(prepared_dir, "pruning_ledger.csv"))
write_csv(data.frame(source = names(prepared$tf), score = unname(prepared$tf)),
          file.path(prepared_dir, "upstream_input.csv"))
write_csv(data.frame(source = names(prepared$lipid), score = unname(prepared$lipid)),
          file.path(prepared_dir, "downstream_input.csv"))

protected <- unique(c(names(prepared$tf), names(prepared$lipid), reg$source, reg$target))
compression_cache <- list(
  compressed = moon_safe_compress(prepared$net, protected),
  uncompressed = moon_safe_compress(prepared$net, node_ids(prepared$net))
)
for (mode in names(compression_cache)) {
  cache_dir <- file.path(prepared_dir, mode)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(compression_cache[[mode]], file.path(cache_dir, "compression.rds"))
  write_csv(compression_cache[[mode]]$node_map, file.path(cache_dir, "compression_mapping.csv"))
}

score_once <- function(score_net, n_layers) {
  ans <- cosmosR::moon(
    upstream_input = prepared$tf,
    downstream_input = prepared$lipid,
    meta_network = score_net,
    n_layers = n_layers,
    statistic = "ulm"
  )
  ans <- as.data.frame(ans[, c("source", "score", "level")])
  stopifnot(!anyDuplicated(ans$source), all(is.finite(ans$score)),
            all(names(prepared$lipid) %in% ans$source[ans$level == 0]))
  ans
}

annotate_scores <- function(scores) {
  catalog <- utils::read.csv("results/cosmos_mouse_pkn/nodes.csv", check.names = FALSE)
  x <- scores
  x$type <- ifelse(x$source %in% names(prepared$tf), "upstream_input",
                   ifelse(x$level == 0, "level0", "other"))
  x$RNA_input <- unname(rna[x$source])
  x$upstream_activity <- unname(tf_all[x$source])
  merge(x, catalog, by.x = "source", by.y = "node_id", all.x = TRUE, sort = FALSE)
}

run_id_for <- function(mode, depth) {
  paste0("radius", topology_radius, "_depth", depth, "_", mode)
}

run_one <- function(mode, depth) {
  run_id <- run_id_for(mode, depth)
  run_dir <- file.path(out_root, "runs", run_id)
  state_path <- file.path(run_dir, "run_state.rds")
  manifest_path <- file.path(run_dir, "run_manifest.json")
  if (!force && file.exists(state_path) && file.exists(manifest_path)) {
    cached <- readRDS(state_path)
    if (identical(cached$run_id, run_id) && identical(cached$n_layers, depth) &&
        identical(cached$compression, mode) && identical(cached$topology_radius, topology_radius)) {
      cat("Reusing completed scoring run:", run_id, "\n")
      return(cached)
    }
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- compression_cache[[mode]]
  score_net <- cache$pkn_scoring
  mapping <- cache$node_map
  ledger <- list()
  removed <- list()
  converged <- FALSE
  for (iteration in seq_len(10L)) {
    used <- score_net
    cat("[", run_id, "] coherence iteration ", iteration, ": ", nrow(used), " edges\n", sep = "")
    scored_compressed <- score_once(used, depth)
    filtered <- cosmosR::filter_incohrent_TF_target(scored_compressed, reg, used, rna)
    delta <- used[!edge_key(used) %in% edge_key(filtered), , drop = FALSE]
    if (nrow(delta)) {
      delta$iteration <- iteration
      removed[[length(removed) + 1L]] <- delta
    }
    ledger[[iteration]] <- data.frame(
      iteration = iteration, edges_scored = nrow(used), edges_removed = nrow(delta),
      scored_nodes = nrow(scored_compressed),
      upstream_scored = sum(names(prepared$tf) %in% scored_compressed$source),
      max_level = max(scored_compressed$level)
    )
    score_net <- filtered
    if (identical(edge_key(used), edge_key(filtered))) {
      converged <- TRUE
      break
    }
  }
  if (!converged) scored_compressed <- score_once(score_net, depth)
  scores <- moon_safe_decompress_scores(scored_compressed, mapping)
  scores <- scores[, c("source", "score", "level", "compressed_source")]
  stopifnot(!anyDuplicated(scores$source), all(is.finite(scores$score)))
  map_vector <- stats::setNames(mapping$compressed, mapping$original)
  projected <- prepared$net
  projected$source <- unname(map_vector[projected$source])
  projected$target <- unname(map_vector[projected$target])
  final <- prepared$net[edge_key(projected) %in% edge_key(score_net), , drop = FALSE]
  stopifnot(!any(final$source %in% names(prepared$lipid)))
  write_csv(do.call(rbind, ledger), file.path(run_dir, "coherence_ledger.csv"))
  write_csv(if (length(removed)) do.call(rbind, removed) else
    data.frame(source = character(), interaction = numeric(), target = character(), iteration = integer()),
    file.path(run_dir, "coherence_removed_edges.csv"))
  write_csv(scores, file.path(run_dir, "moon_scores.csv"))
  write_csv(annotate_scores(scores), file.path(run_dir, "moon_scores_annotated.csv"))
  saveRDS(score_net, file.path(run_dir, "pkn_final_scoring.rds"))
  saveRDS(final, file.path(run_dir, "pkn_final_uncompressed.rds"))
  state <- list(run_id = run_id, run_dir = run_dir, compression = mode,
    n_layers = depth, topology_radius = topology_radius, converged = converged,
    iterations = length(ledger), tf = prepared$tf, lipid = prepared$lipid,
    scores = scores, pkn_final = final,
    realized_max_level = max(scores$level), scored_nodes = nrow(scores),
    final_scoring_edges = nrow(score_net), final_original_edges = nrow(final),
    topology_reference = normalizePath(file.path(prepared_dir, "prepared_topology.rds")))
  saveRDS(state, state_path)
  manifest <- list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    command = "Rscript scripts/moon/search_factor1_parameter_space.R --phase=coarse",
    run_id = run_id, compression = mode, n_layers = depth, topology_radius = topology_radius,
    topology_radius_policy = "fixed across n_layers to isolate MOON propagation depth",
    pkn = pkn_path, input = input_path, statistic = "ulm", tf_cutoff = "abs(activity) > 2",
    expression_filter = "none: Factor 1 RNA weights are not expression-presence measurements",
    downstream_scaling = "none; observed level-0 weights are exempt from inferred-score thresholds",
    converged = converged, iterations = length(ledger),
    upstream_inputs = length(prepared$tf), downstream_inputs = length(prepared$lipid),
    scored_nodes = nrow(scores), realized_max_level = max(scores$level),
    final_scoring_edges = nrow(score_net), final_original_edges = nrow(final),
    orientation = input$orientation,
    input_and_script_md5 = as.list(tools::md5sum(c(
      pkn_path, input_path, "scripts/moon/search_factor1_parameter_space.R",
      "scripts/moon/moon_safe_helpers.R"))),
    template_url = "https://github.com/saezlab/COSMOS_basic/blob/7108c6b849c712b1282440dc122c9d783f563da8/scripts/net_compr_MOON.R"
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  state
}

coarse_thresholds <- function() {
  primary <- c(1.7, seq(2, 5, by = 0.25))
  secondary <- c(0.5, 0.75, seq(1, 3.5, by = 0.25))
  grid <- expand.grid(primary_threshold = primary, secondary_threshold = secondary,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[grid$secondary_threshold <= grid$primary_threshold, , drop = FALSE]
  grid[order(-grid$secondary_threshold, -grid$primary_threshold), , drop = FALSE]
}

extract_threshold_grid <- function(state, thresholds, phase_name) {
  rows <- list()
  for (i in seq_len(nrow(thresholds))) {
    primary <- thresholds$primary_threshold[[i]]
    secondary <- thresholds$secondary_threshold[[i]]
    solution <- moon_extract_paths(state$scores, state$pkn_final, state$tf,
                                   primary_thresh = primary, secondary_thresh = secondary)
    att <- solution$ATT
    sif <- solution$SIF
    score_map <- stats::setNames(state$scores$score, state$scores$source)
    valid <- all(edge_key(sif) %in% edge_key(state$pkn_final)) &&
      all(sign(score_map[sif$source] * score_map[sif$target]) == sif$interaction) &&
      !any(sif$source %in% names(state$lipid)) &&
      isTRUE(solution$audit$all_edges_consistent) &&
      isTRUE(solution$audit$all_nodes_on_seed_to_observation_walk)
    rows[[i]] <- data.frame(
      run_id = state$run_id, run_dir = state$run_dir,
      compression = state$compression, depth = state$n_layers,
      n_layers = state$n_layers, topology_radius = state$topology_radius,
      phase = phase_name, primary = primary, secondary = secondary,
      primary_threshold = primary, secondary_threshold = secondary,
      status = "complete", converged = state$converged, audit_passed = valid, valid = valid,
      nodes = nrow(att), edges = nrow(sif), upstream_tfs = sum(att$type == "upstream_input"),
      downstream_lipids = sum(att$level == 0), downstream_nodes = sum(att$level == 0),
      realized_max_level = state$realized_max_level, scored_nodes = state$scored_nodes,
      candidate_edges = solution$audit$candidate_edges,
      consistent_edges = solution$audit$consistent_edges,
      route_pairs = NA_integer_, short_route_pairs = NA_integer_,
      solution_stem = NA_character_, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

grid_path <- file.path(out_root, "parameter_grid.csv")
update_grid <- function(new_rows) {
  old <- if (file.exists(grid_path)) utils::read.csv(grid_path, check.names = FALSE, stringsAsFactors = FALSE) else new_rows[0, ]
  if (nrow(old)) {
    same <- old$run_id %in% new_rows$run_id &
      old$primary_threshold %in% new_rows$primary_threshold &
      old$secondary_threshold %in% new_rows$secondary_threshold
    old <- old[!same, , drop = FALSE]
  }
  all_rows <- rbind(old, new_rows)
  all_rows <- all_rows[order(all_rows$compression, all_rows$depth,
                             all_rows$secondary_threshold, all_rows$primary_threshold), , drop = FALSE]
  write_csv(all_rows, grid_path)
  all_rows
}

materialize_solution <- function(state, row) {
  stem <- paste0("solution_p", formatC(row$primary_threshold, format = "f", digits = 2),
                 "_s", formatC(row$secondary_threshold, format = "f", digits = 2))
  stem <- gsub("\\.", "_", stem)
  solution <- moon_extract_paths(state$scores, state$pkn_final, state$tf,
                                 primary_thresh = row$primary_threshold,
                                 secondary_thresh = row$secondary_threshold)
  output <- state$run_dir
  annotated <- annotate_scores(solution$ATT)
  write_csv(solution$SIF, file.path(output, paste0(stem, "_SIF.csv")))
  write_csv(annotated, file.path(output, paste0(stem, "_ATT.csv")))
  stem
}

if (phase == "coarse") {
  all_rows <- list()
  run_rows <- list()
  for (depth in depths) {
    for (mode in compression_modes) {
      run_id <- run_id_for(mode, depth)
      state <- tryCatch(run_one(mode, depth), error = function(e) {
        run_dir <- file.path(out_root, "runs", run_id)
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines(conditionMessage(e), file.path(run_dir, "ERROR.txt"))
        NULL
      })
      if (is.null(state)) {
        run_rows[[length(run_rows) + 1L]] <- data.frame(
          run_id = run_id, run_dir = file.path(out_root, "runs", run_id), compression = mode,
          depth = depth, n_layers = depth, topology_radius = topology_radius,
          status = "failed", converged = FALSE, audit_passed = FALSE, valid = FALSE,
          nodes = NA_integer_, edges = NA_integer_, upstream_tfs = NA_integer_,
          downstream_lipids = NA_integer_, downstream_nodes = NA_integer_,
          realized_max_level = NA_integer_, scored_nodes = NA_integer_,
          final_scoring_edges = NA_integer_, final_original_edges = NA_integer_, stringsAsFactors = FALSE)
      } else {
        run_rows[[length(run_rows) + 1L]] <- data.frame(
          run_id = state$run_id, run_dir = state$run_dir, compression = state$compression,
          depth = state$n_layers, n_layers = state$n_layers, topology_radius = state$topology_radius,
          status = "complete", converged = state$converged, audit_passed = TRUE, valid = TRUE,
          nodes = NA_integer_, edges = NA_integer_, upstream_tfs = length(state$tf),
          downstream_lipids = length(state$lipid), downstream_nodes = length(state$lipid),
          realized_max_level = state$realized_max_level, scored_nodes = state$scored_nodes,
          final_scoring_edges = state$final_scoring_edges, final_original_edges = state$final_original_edges,
          stringsAsFactors = FALSE)
        all_rows[[length(all_rows) + 1L]] <- extract_threshold_grid(state, coarse_thresholds(), "coarse")
      }
    }
  }
  write_csv(do.call(rbind, run_rows), file.path(out_root, "run_summary.csv"))
  combined <- update_grid(do.call(rbind, all_rows))
  cat("Coarse threshold search complete:", nrow(combined), "rows across", length(run_rows), "runs.\n")
} else {
  stopifnot(file.exists(grid_path))
  existing <- utils::read.csv(grid_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_grid <- c("run_id", "run_dir", "depth", "primary_threshold", "secondary_threshold",
                     "edges", "upstream_tfs", "downstream_lipids", "converged", "audit_passed")
  stopifnot(all(required_grid %in% names(existing)))
  in_requested_scope <- existing$topology_radius == topology_radius &
    existing$compression %in% compression_modes & existing$depth %in% depths
  eligible <- existing[in_requested_scope & existing$status == "complete" & existing$converged & existing$audit_passed &
                         existing$edges > 0 & existing$edges < 500 &
                         existing$upstream_tfs > 0 & existing$downstream_lipids > 0, , drop = FALSE]
  if (!nrow(eligible)) stop("No completed <500-edge candidate exists in the requested topology, compression, and depth scope to refine.")
  # Pick one high-threshold candidate per independently scored run.  Requiring
  # two TFs and two lipids when possible avoids using a trivial single path to
  # define the local search window.
  interpretable <- eligible[eligible$upstream_tfs >= 2 & eligible$downstream_lipids >= 2, , drop = FALSE]
  pool <- if (nrow(interpretable)) interpretable else eligible
  selected <- do.call(rbind, lapply(split(pool, pool$run_id), function(x) {
    x <- x[order(-x$secondary_threshold, -x$primary_threshold,
                 -x$downstream_lipids, -x$upstream_tfs, x$edges), , drop = FALSE]
    x[1, , drop = FALSE]
  }))
  all_rows <- list()
  for (i in seq_len(nrow(selected))) {
    pick <- selected[i, , drop = FALSE]
    state <- readRDS(file.path(pick$run_dir, "run_state.rds"))
    p <- pick$primary_threshold
    s <- pick$secondary_threshold
    pvals <- round(seq(max(1.7, p - 0.25), min(5, p + 0.25), by = 0.05), 2)
    svals <- round(seq(max(0.25, s - 0.25), min(4, s + 0.25), by = 0.05), 2)
    grid <- expand.grid(primary_threshold = pvals, secondary_threshold = svals,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid <- grid[grid$secondary_threshold <= grid$primary_threshold, , drop = FALSE]
    all_rows[[length(all_rows) + 1L]] <- extract_threshold_grid(state, grid, "refined")
  }
  combined <- update_grid(do.call(rbind, all_rows))
  cat("Refined threshold search complete:", nrow(combined), "rows.\n")
}

# The independent ranking script enforces convergence/audit/edge constraints
# and writes the full-grid and Pareto summaries.  It does not alter scoring.
source("scripts/moon/rank_parameter_candidates.R")
grid <- utils::read.csv(grid_path, check.names = FALSE, stringsAsFactors = FALSE)
ranked <- moon_rank_parameter_candidates(grid, max_edges = 500L)
eligible_rank <- ranked[ranked$eligible, , drop = FALSE]
frontier <- ranked[ranked$eligible & ranked$threshold_frontier_within_run, , drop = FALSE]
interpretable_rank <- ranked[ranked$eligible & ranked$multiple_tfs_and_lipids, , drop = FALSE]
write_csv(ranked, file.path(out_root, "parameter_grid_ranked.csv"))
write_csv(eligible_rank, file.path(out_root, "eligible_under_500_edges.csv"))
write_csv(frontier, file.path(out_root, "within_run_threshold_frontier.csv"))
write_csv(interpretable_rank, file.path(out_root, "interpretable_candidates_under_500_edges.csv"))

# Materialize only the current parameter-ranked leader for each scoring run.
# Complete SIF/ATT tables are kept separate from the exhaustive grid so they
# can be independently provenance-exported after review.
leaders <- interpretable_rank
if (!nrow(leaders)) leaders <- eligible_rank
if (nrow(leaders)) {
  leaders <- do.call(rbind, lapply(split(leaders, leaders$run_id), function(x) x[1, , drop = FALSE]))
  leaders$solution_stem <- vapply(seq_len(nrow(leaders)), function(i) {
    state <- readRDS(file.path(leaders$run_dir[[i]], "run_state.rds"))
    materialize_solution(state, leaders[i, , drop = FALSE])
  }, character(1))
  write_csv(leaders, file.path(out_root, "run_leaders_under_500_edges.csv"))
}

plan <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  phase = phase, depths = depths, topology_radius = topology_radius,
  compression_modes = compression_modes,
  score_settings = list(statistic = "ulm", tf_cutoff = "abs(activity) > 2",
                        max_coherence_iterations = 10),
  threshold_policy = list(
    coarse_primary = c(1.7, seq(2, 5, by = 0.25)),
    coarse_secondary = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5),
    refined = "0.05 local grid around one best <=500-edge candidate per independently scored run",
    ranking = "secondary descending, primary descending, coverage descending, edges ascending; all rows retained"),
  mapping = list(upstream = "Factor 1 hepatocyte TF footprint activities", downstream = "exact-mapped Factor 1 lipid weights", rna = "Factor 1 hepatocyte RNA weights for coherence"),
  input_and_script_md5 = as.list(tools::md5sum(c(input_path, pkn_path,
    "scripts/moon/search_factor1_parameter_space.R", "scripts/moon/moon_safe_helpers.R",
    "scripts/moon/rank_parameter_candidates.R")))
)
jsonlite::write_json(plan, file.path(out_root, "search_plan.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
cat("Wrote parameter search artifacts to", out_root, "\n")
