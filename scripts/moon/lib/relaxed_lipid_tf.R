# Shared operations for the relaxed Factor 1 lipid-to-TF screen.
# Source after moon_safe_helpers.R and workflow.R. This module deliberately
# keeps the candidate-chemical semantics separate from ordinary MOON inputs.

relaxed_required_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

relaxed_truth <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes")
relaxed_protein_node <- function(x) grepl("^(Gene|UniProt)__", x)
relaxed_edge_key <- function(x) paste(x$source, x$interaction, x$target, sep = "\r")
relaxed_write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", null = "null", digits = NA)
}
relaxed_collapse <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  paste(x, collapse = ";")
}

relaxed_read_config <- function(path) {
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)
  required <- c("analysis", "inputs", "pkn", "mapping", "direct_evidence", "output", "selected",
                "max_iterations", "strong_tf_cutoff", "maximum_solution_edges_exclusive",
                "minimum_independent_edge_fingerprint_groups", "minimum_downstream_tfs",
                "packages", "cosmos_commit", "input_md5")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("Relaxed selected config missing: ", paste(missing, collapse = ", "))
  if (!identical(as.character(cfg$analysis), "factor1_lipid_to_tf_relaxed")) {
    stop("Unexpected relaxed analysis identifier: ", cfg$analysis)
  }
  selected_required <- c("tf_panel", "pruning_steps", "n_layers", "primary", "secondary",
                         "expected_measurement_feature", "expected_source_node", "expected_observed_tfs",
                         "reference_fit", "reference_solution", "reference_fit_md5",
                         "reference_solution_sif_md5")
  missing_selected <- setdiff(selected_required, names(cfg$selected))
  if (length(missing_selected)) stop("Selected relaxed config missing: ", paste(missing_selected, collapse = ", "))
  stopifnot(cfg$max_iterations >= 1L, cfg$strong_tf_cutoff > 0,
            cfg$maximum_solution_edges_exclusive > 0L,
            cfg$minimum_independent_edge_fingerprint_groups >= 1L,
            cfg$minimum_downstream_tfs >= 2L,
            cfg$selected$tf_panel %in% c("all_tfs", "strong_tfs"),
            cfg$selected$pruning_steps >= 1L, cfg$selected$n_layers >= 1L,
            cfg$selected$primary >= cfg$selected$secondary,
            cfg$selected$secondary > 0)
  for (pkg in names(cfg$packages)) {
    if (!requireNamespace(pkg, quietly = TRUE) ||
        as.character(utils::packageVersion(pkg)) != as.character(cfg$packages[[pkg]])) {
      stop("Relaxed selected replay requires ", pkg, " ", cfg$packages[[pkg]])
    }
  }
  sha <- utils::packageDescription("cosmosR")$RemoteSha
  if (!identical(as.character(sha), as.character(cfg$cosmos_commit))) {
    stop("cosmosR commit differs from the selected relaxed analysis")
  }
  for (dependency in names(cfg$input_md5)) {
    if (!file.exists(dependency) ||
        !identical(unname(tools::md5sum(dependency)), as.character(cfg$input_md5[[dependency]]))) {
      stop("Missing or changed relaxed dependency: ", dependency)
    }
  }
  for (dependency in c(cfg$selected$reference_fit, file.path(cfg$selected$reference_solution, "network_SIF.csv"))) {
    if (!file.exists(dependency)) stop("Selected relaxed reference is missing: ", dependency)
  }
  if (!identical(unname(tools::md5sum(cfg$selected$reference_fit)),
                 as.character(cfg$selected$reference_fit_md5))) {
    stop("Selected relaxed reference fit changed: ", cfg$selected$reference_fit)
  }
  reference_sif <- file.path(cfg$selected$reference_solution, "network_SIF.csv")
  if (!identical(unname(tools::md5sum(reference_sif)),
                 as.character(cfg$selected$reference_solution_sif_md5))) {
    stop("Selected relaxed reference SIF changed: ", reference_sif)
  }
  attr(cfg, "config_path") <- path
  cfg
}

relaxed_read_search_config <- function(path) {
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)
  required <- c("analysis", "base_config", "output", "pruning_steps", "n_layers", "tf_panels",
                "primary_thresholds", "secondary_thresholds", "maximum_solution_edges_exclusive",
                "minimum_independent_edge_fingerprint_groups", "minimum_downstream_tfs")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("Relaxed parameter-search config missing: ", paste(missing, collapse = ", "))
  if (!identical(as.character(cfg$analysis), "factor1_lipid_to_tf_relaxed_parameter_search")) {
    stop("Unexpected relaxed search identifier: ", cfg$analysis)
  }
  cfg$pruning_steps <- as.integer(cfg$pruning_steps)
  cfg$n_layers <- as.integer(cfg$n_layers)
  cfg$primary_thresholds <- as.numeric(cfg$primary_thresholds)
  cfg$secondary_thresholds <- as.numeric(cfg$secondary_thresholds)
  stopifnot(length(cfg$pruning_steps) > 0L, all(is.finite(cfg$pruning_steps)), all(cfg$pruning_steps > 0),
            length(cfg$n_layers) > 0L, all(is.finite(cfg$n_layers)), all(cfg$n_layers > 0),
            length(cfg$primary_thresholds) > 0L, all(is.finite(cfg$primary_thresholds)), all(cfg$primary_thresholds > 0),
            length(cfg$secondary_thresholds) > 0L, all(is.finite(cfg$secondary_thresholds)), all(cfg$secondary_thresholds >= 0),
            length(cfg$tf_panels) > 0L, all(cfg$tf_panels %in% c("all_tfs", "strong_tfs")),
            cfg$maximum_solution_edges_exclusive > 0L,
            cfg$minimum_independent_edge_fingerprint_groups >= 1L,
            cfg$minimum_downstream_tfs >= 2L,
            file.exists(cfg$base_config))
  attr(cfg, "config_path") <- path
  cfg
}

relaxed_load_data <- function(cfg) {
  input <- readRDS(cfg$inputs)
  relaxed_required_columns(input, c("tf_activity", "hepatocyte_gene_weight", "tf_regulon", "orientation"),
                           "Factor 1 input bundle")
  moon_vector_check(input$tf_activity)
  moon_vector_check(input$hepatocyte_gene_weight)
  overlay <- moon_check_pkn(readRDS(cfg$pkn))
  assignments <- utils::read.csv(cfg$mapping, check.names = FALSE)
  relaxed_required_columns(assignments,
    c("assignment_id", "measurement_feature", "candidate_measurement_value", "pkn_source_node",
      "pkn_source_chebi", "pkn_source_name", "primary_mapping_policy", "edge_fingerprint_group",
      "default_relaxed_candidate"),
    "Relaxed mapping ledger")
  assignments <- assignments[relaxed_truth(assignments$default_relaxed_candidate), , drop = FALSE]
  if (!nrow(assignments)) stop("No default relaxed assignments in mapping ledger")
  assignments$candidate_measurement_value <- as.numeric(assignments$candidate_measurement_value)
  if (any(!is.finite(assignments$candidate_measurement_value)) ||
      any(assignments$candidate_measurement_value == 0)) {
    stop("Relaxed mapping ledger has non-finite or zero Factor 1 measurements")
  }
  if (anyDuplicated(assignments$assignment_id)) stop("Relaxed assignment IDs must be unique")
  if (!all(grepl("^Metab__", assignments$pkn_source_node))) {
    stop("Relaxed candidates must be Metab__ PKN source nodes")
  }
  candidate_nodes <- sort(unique(assignments$pkn_source_node))
  if (!all(candidate_nodes %in% moon_nodes(overlay))) {
    stop("Mapping ledger contains candidate node(s) absent from the directional overlay")
  }
  direct_pkn <- overlay[
    (grepl("^Metab__", overlay$source) & relaxed_protein_node(overlay$target)) |
      (relaxed_protein_node(overlay$source) & relaxed_protein_node(overlay$target)),
    c("source", "target", "interaction"), drop = FALSE
  ]
  direct_pkn <- moon_check_pkn(direct_pkn)
  if (!all(candidate_nodes %in% moon_nodes(direct_pkn))) {
    stop("Candidate node(s) lack a direct regulatory representation in the direct signaling graph")
  }
  topology_seed <- stats::setNames(rep(1, length(candidate_nodes)), candidate_nodes)
  moon_vector_check(topology_seed)
  tf <- input$tf_activity[input$tf_activity != 0 & names(input$tf_activity) %in% moon_nodes(direct_pkn)]
  if (!length(tf)) stop("No Factor 1 TF activities are present in the direct signaling PKN")
  panels <- list(all_tfs = tf, strong_tfs = tf[abs(tf) > cfg$strong_tf_cutoff])
  if (any(vapply(panels, length, integer(1)) == 0L)) stop("A configured TF panel is empty")
  regulon <- input$tf_regulon[, c("source", "target", "mor"), drop = FALSE]
  relaxed_required_columns(regulon, c("source", "target", "mor"), "TF regulon")
  if (anyNA(regulon) || !all(regulon$mor %in% c(-1, 1))) stop("TF regulon must have signed complete rows")
  list(input = input, overlay = overlay, assignments = assignments, candidate_nodes = candidate_nodes,
       topology_seed = topology_seed, direct_pkn = direct_pkn, panels = panels,
       regulon = regulon, rna = input$hepatocyte_gene_weight)
}

relaxed_panel_data <- function(data, panel_name) {
  if (!panel_name %in% names(data$panels)) stop("Unknown relaxed TF panel: ", panel_name)
  list(pkn = data$direct_pkn, topology_seed = data$topology_seed, downstream = data$panels[[panel_name]],
       regulon = data$regulon, rna = data$rna)
}

# Candidate numerical values do not enter ULM score construction. The seed keeps
# candidate IDs in the directed corridor; values are checked after scoring.
relaxed_fit <- function(data, pruning_steps, n_layers, max_iterations = 10L) {
  prepared <- moon_prepare(data$pkn, data$topology_seed, data$downstream, pruning_steps)
  protected <- unique(c(names(prepared$upstream), names(prepared$downstream),
                         data$regulon$source, data$regulon$target))
  compression <- moon_safe_compress(prepared$net, protected)
  scoring <- compression$pkn_scoring
  ledger <- removed <- list()
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    message("MOON iteration ", iteration, ": ", nrow(scoring), " compressed edges")
    scores <- moon_score(scoring, NULL, prepared$downstream, n_layers)
    filtered <- cosmosR::filter_incohrent_TF_target(scores, data$regulon, scoring, data$rna)
    stopifnot(all(moon_edge_keys(filtered) %in% moon_edge_keys(scoring)))
    delta <- scoring[!moon_edge_keys(scoring) %in% moon_edge_keys(filtered), , drop = FALSE]
    ledger[[iteration]] <- data.frame(
      iteration = iteration, edges_scored = nrow(scoring), edges_removed = nrow(delta),
      scored_nodes = nrow(scores), topology_seed_nodes = length(prepared$upstream),
      max_level = max(scores$level), stringsAsFactors = FALSE
    )
    if (nrow(delta)) {
      delta$iteration <- iteration
      removed[[length(removed) + 1L]] <- delta
    }
    if (identical(moon_edge_keys(scoring), moon_edge_keys(filtered))) {
      converged <- TRUE
      break
    }
    scoring <- filtered
  }
  if (!converged) stop("RNA coherence did not converge in ", max_iterations, " iterations")
  node_map <- stats::setNames(compression$node_map$compressed, compression$node_map$original)
  projected <- prepared$net
  projected$source <- unname(node_map[projected$source])
  projected$target <- unname(node_map[projected$target])
  final <- prepared$net[moon_edge_keys(projected) %in% moon_edge_keys(scoring), , drop = FALSE]
  stopifnot(setequal(
    moon_edge_keys(projected[moon_edge_keys(projected) %in% moon_edge_keys(scoring), , drop = FALSE]),
    moon_edge_keys(scoring)
  ))
  list(
    prepared = prepared, compression = compression, pkn_scoring = scoring,
    pkn_final = final, scores_compressed = scores,
    scores = moon_safe_decompress_scores(scores, compression$node_map),
    coherence_ledger = do.call(rbind, ledger),
    removed_edges = if (length(removed)) do.call(rbind, removed) else data.frame(
      source = character(), target = character(), interaction = numeric(), iteration = integer()),
    pruning_steps = pruning_steps, n_layers = n_layers, converged = converged,
    scoring_mode = "unconstrained_candidate_scores"
  )
}

relaxed_candidate_status <- function(assignments, fit) {
  score_map <- stats::setNames(fit$scores$score, fit$scores$source)
  level_map <- stats::setNames(fit$scores$level, fit$scores$source)
  retained_nodes <- moon_nodes(fit$prepared$net)
  out <- assignments
  out$candidate_in_prepared_topology <- out$pkn_source_node %in% retained_nodes
  out$candidate_scored <- out$pkn_source_node %in% names(score_map)
  out$moon_score <- unname(score_map[out$pkn_source_node])
  out$moon_level <- unname(level_map[out$pkn_source_node])
  out$assignment_sign_agrees <- out$candidate_scored &
    sign(out$moon_score) == sign(out$candidate_measurement_value)
  out$assignment_status <- ifelse(!out$candidate_in_prepared_topology, "pruned_no_candidate_to_TF_corridor",
    ifelse(!out$candidate_scored, "not_scored_within_realized_layers",
      ifelse(out$assignment_sign_agrees, "scored_sign_concordant", "scored_sign_discordant")))
  out
}

relaxed_candidate_node_summary <- function(status) {
  groups <- split(seq_len(nrow(status)), status$pkn_source_node)
  rows <- lapply(names(groups), function(node) {
    x <- status[groups[[node]], , drop = FALSE]
    present <- as.logical(x$candidate_scored)
    values <- sort(unique(sign(as.numeric(x$candidate_measurement_value))))
    data.frame(
      pkn_source_node = node,
      pkn_source_chebi = x$pkn_source_chebi[1],
      pkn_source_name = x$pkn_source_name[1],
      n_measurement_assignments = nrow(x),
      measurement_features = relaxed_collapse(x$measurement_feature),
      assigned_loading_signs = relaxed_collapse(values),
      assignment_signs_are_consistent = length(values) == 1L,
      candidate_in_prepared_topology = all(as.logical(x$candidate_in_prepared_topology)),
      candidate_scored = all(present),
      moon_score = if (all(present)) x$moon_score[1] else NA_real_,
      moon_level = if (all(present)) x$moon_level[1] else NA_real_,
      assignments_sign_concordant = sum(as.logical(x$assignment_sign_agrees)),
      assignments_sign_discordant = sum(present & !as.logical(x$assignment_sign_agrees)),
      all_assignments_sign_concordant = all(present) && all(as.logical(x$assignment_sign_agrees)),
      edge_fingerprint_groups = relaxed_collapse(x$edge_fingerprint_group),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

relaxed_sign_consensus_seeds <- function(status) {
  summary <- relaxed_candidate_node_summary(status)
  keep <- summary$candidate_scored & summary$all_assignments_sign_concordant &
    summary$moon_level > 0 & summary$moon_score != 0
  signs <- vapply(summary$pkn_source_node[keep], function(node) {
    values <- unique(sign(as.numeric(status$candidate_measurement_value[status$pkn_source_node == node])))
    if (length(values) != 1L) stop("Conflicting measured signs for candidate source: ", node)
    values
  }, numeric(1))
  stats::setNames(signs, summary$pkn_source_node[keep])
}

# Unlike moon_extract_paths(), this keeps only actual recursive-score transitions.
relaxed_extract_score_layer_paths <- function(scores, pkn, upstream_input,
                                              primary, secondary, RNA_input = NULL,
                                              pkn_already_checked = FALSE) {
  relaxed_required_columns(scores, c("source", "score", "level"), "MOON scores")
  stopifnot(!anyDuplicated(scores$source), all(is.finite(scores$score)), all(is.finite(scores$level)),
            primary >= secondary, secondary >= 0, is.numeric(upstream_input),
            !is.null(names(upstream_input)), !anyDuplicated(names(upstream_input)),
            all(is.finite(upstream_input)))
  if (!pkn_already_checked) pkn <- moon_check_pkn(pkn)
  scores <- scores[, c("source", "score", "level"), drop = FALSE]
  scores$source <- as.character(scores$source)
  observed <- scores$source[scores$level == 0L & scores$score != 0]
  primary_nodes <- scores$source[scores$level > 0L & abs(scores$score) > primary]
  secondary_nodes <- union(observed, scores$source[scores$level > 0L & abs(scores$score) > secondary])
  score_map <- stats::setNames(scores$score, scores$source)
  level_map <- stats::setNames(scores$level, scores$source)
  score_index <- match(names(upstream_input), scores$source)
  input_ok <- !is.na(score_index)
  input_ok[input_ok] <- scores$level[score_index[input_ok]] > 0L &
    sign(scores$score[score_index[input_ok]]) == sign(upstream_input[input_ok])
  seeds <- intersect(names(upstream_input)[input_ok], secondary_nodes)
  failed_inputs <- setdiff(intersect(names(upstream_input), scores$source), seeds)
  secondary_nodes <- setdiff(secondary_nodes, failed_inputs)
  primary_nodes <- intersect(primary_nodes, secondary_nodes)
  net <- pkn[pkn$source %in% secondary_nodes & pkn$target %in% secondary_nodes &
               (pkn$source %in% primary_nodes | pkn$target %in% primary_nodes) &
               !pkn$source %in% observed, , drop = FALSE]
  candidate_edges <- nrow(net)
  net <- net[level_map[net$source] == level_map[net$target] + 1L, , drop = FALSE]
  score_layer_edges <- nrow(net)
  net$consistency <- sign(score_map[net$source] * score_map[net$target]) == net$interaction
  net <- net[net$consistency, , drop = FALSE]
  consistent_edges <- nrow(net)
  empty <- function() list(
    SIF = data.frame(source = character(), target = character(), interaction = numeric(), consistency = logical()),
    ATT = data.frame(source = character(), score = numeric(), level = numeric(), type = character(), RNA_input = numeric()),
    audit = list(primary = primary, secondary = secondary, score_layer_constrained = TRUE,
                 candidate_edges = candidate_edges, score_layer_edges = score_layer_edges,
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
    root <- "__RELAXED_MOON_SCORE_LAYER_ROOT__"
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
  att$RNA_input <- if (is.null(RNA_input)) NA_real_ else unname(RNA_input[att$source])
  final_graph <- igraph::graph_from_data_frame(net[, c("source", "target")], directed = TRUE)
  retained_seeds <- intersect(seeds, nodes)
  retained_observed <- intersect(observed, nodes)
  verified <- intersect(reachable(final_graph, retained_seeds),
                        reachable(final_graph, retained_observed, reverse = TRUE))
  stopifnot(setequal(nodes, verified), all(net$consistency),
            all(level_map[net$source] == level_map[net$target] + 1L), setequal(nodes, att$source))
  rownames(net) <- rownames(att) <- NULL
  list(
    SIF = net,
    ATT = att,
    audit = list(primary = primary, secondary = secondary, score_layer_constrained = TRUE,
                 candidate_edges = candidate_edges, score_layer_edges = score_layer_edges,
                 consistent_edges = consistent_edges, final_edges = nrow(net), final_nodes = nrow(att),
                 retained_upstream = length(retained_seeds), retained_observed = length(retained_observed),
                 all_edges_descend_one_score_level = TRUE,
                 all_nodes_on_seed_to_observation_walk = TRUE)
  )
}

relaxed_validate_score_layer_solution <- function(solution, fit, upstream_input, primary, secondary) {
  if (!nrow(solution$SIF)) return(invisible(TRUE))
  score_map <- stats::setNames(fit$scores$score, fit$scores$source)
  level_map <- stats::setNames(fit$scores$level, fit$scores$source)
  source_nodes <- solution$ATT$source[solution$ATT$type == "upstream_input"]
  stopifnot(
    all(relaxed_edge_key(solution$SIF) %in% relaxed_edge_key(fit$pkn_final)),
    all(sign(score_map[solution$SIF$source] * score_map[solution$SIF$target]) == solution$SIF$interaction),
    all(level_map[solution$SIF$source] == level_map[solution$SIF$target] + 1L),
    all(abs(solution$ATT$score[solution$ATT$level > 0L]) > secondary),
    all(abs(solution$ATT$score[solution$ATT$level > 0L & solution$ATT$source %in% source_nodes]) > secondary),
    all(sign(score_map[source_nodes]) == sign(upstream_input[source_nodes])),
    all(solution$SIF$source %in% solution$ATT$source[solution$ATT$level > 0L & abs(solution$ATT$score) > primary] |
        solution$SIF$target %in% solution$ATT$source[solution$ATT$level > 0L & abs(solution$ATT$score) > primary]),
    isTRUE(solution$audit$all_nodes_on_seed_to_observation_walk)
  )
  invisible(TRUE)
}

relaxed_solution_metadata <- function(solution, status) {
  sources <- solution$ATT$source[solution$ATT$type == "upstream_input"]
  source_rows <- status[status$pkn_source_node %in% sources, , drop = FALSE]
  list(
    source_nodes = sort(unique(sources)),
    measurement_features = sort(unique(source_rows$measurement_feature)),
    edge_fingerprint_groups = sort(unique(source_rows$edge_fingerprint_group)),
    raw_upstream_candidate_nodes = length(unique(sources)),
    independent_edge_fingerprint_groups = length(unique(source_rows$edge_fingerprint_group)),
    observed_tfs = sort(solution$ATT$source[solution$ATT$level == 0L]),
    retained_observed_tfs = solution$audit$retained_observed
  )
}

relaxed_solution_is_eligible <- function(solution, status, cfg) {
  metadata <- relaxed_solution_metadata(solution, status)
  eligible <- nrow(solution$SIF) > 0L && nrow(solution$SIF) < cfg$maximum_solution_edges_exclusive &&
    metadata$independent_edge_fingerprint_groups >= cfg$minimum_independent_edge_fingerprint_groups &&
    metadata$retained_observed_tfs >= cfg$minimum_downstream_tfs && solution$audit$secondary > 0
  c(list(eligible = eligible, edges = nrow(solution$SIF), nodes = nrow(solution$ATT)), metadata)
}

relaxed_annotate_solution <- function(solution, status) {
  att <- solution$ATT
  sources <- sort(unique(status$pkn_source_node))
  source_summary <- lapply(sources, function(source) {
    x <- status[status$pkn_source_node == source, , drop = FALSE]
    data.frame(
      source = source,
      candidate_assignment_ids = relaxed_collapse(x$assignment_id),
      candidate_measurement_features = relaxed_collapse(x$measurement_feature),
      candidate_measurement_values = relaxed_collapse(format(as.numeric(x$candidate_measurement_value), scientific = FALSE, trim = TRUE)),
      candidate_factor1_abs_ranks = relaxed_collapse(x$factor1_abs_rank_desc),
      candidate_pkn_source_chebi = relaxed_collapse(x$pkn_source_chebi),
      candidate_pkn_source_names = relaxed_collapse(x$pkn_source_name),
      candidate_mapping_policies = relaxed_collapse(x$primary_mapping_policy),
      candidate_edge_fingerprint_groups = relaxed_collapse(x$edge_fingerprint_group),
      candidate_assignment_signs_agree = relaxed_collapse(x$assignment_sign_agrees),
      candidate_assignment_statuses = relaxed_collapse(x$assignment_status),
      stringsAsFactors = FALSE
    )
  })
  source_summary <- do.call(rbind, source_summary)
  idx <- match(att$source, source_summary$source)
  for (field in setdiff(names(source_summary), "source")) att[[field]] <- source_summary[[field]][idx]
  att
}

relaxed_write_common_inputs <- function(data, cfg, config_path, out) {
  moon_write_csv(data$assignments, file.path(out, "relaxed_candidate_mapping.csv"))
  moon_write_csv(data.frame(pkn_source_node = names(data$topology_seed), topology_seed_value = unname(data$topology_seed)),
                 file.path(out, "candidate_topology_seeds.csv"))
  moon_write_csv(data$direct_pkn, file.path(out, "pkn_direct_signaling.csv"))
  file.copy(config_path, file.path(out, "analysis_config.json"), overwrite = FALSE)
  relaxed_write_json(list(
    user_direction = "lipids -> TFs",
    factor = "Factor1",
    orientation = data$input$orientation,
    mapping = "Relaxed feature-to-chemical assignments; collisions remain separate rows and are never summed or averaged.",
    upstream = "Candidate chemical node IDs are topology seeds only. Their numeric seed values are not MOON measurements.",
    downstream = cfg$downstream_semantics,
    score_construction = "ULM scores are calculated with upstream_input=NULL. Factor 1 lipid loadings are joined afterwards for assignment-level sign concordance.",
    pkn_policy = cfg$pkn_semantics,
    timing = "Cross-sectional signatures; no temporal ordering or direct hepatocyte ligand exposure is inferred.",
    limitations = cfg$limitations
  ), file.path(out, "input_decision.json"))
}

relaxed_compare_reference <- function(fit, solution, cfg) {
  old_fit <- readRDS(cfg$selected$reference_fit)
  reference_sif <- utils::read.csv(file.path(cfg$selected$reference_solution, "network_SIF.csv"), check.names = FALSE)
  reference_att <- utils::read.csv(file.path(cfg$selected$reference_solution, "network_ATT.csv"), check.names = FALSE)
  score_index <- match(old_fit$scores$source, fit$scores$source)
  att_index <- match(reference_att$source, solution$ATT$source)
  stopifnot(
    setequal(old_fit$scores$source, fit$scores$source), !anyNA(score_index),
    setequal(moon_edge_keys(old_fit$pkn_final), moon_edge_keys(fit$pkn_final)),
    isTRUE(all.equal(old_fit$compression$node_map, fit$compression$node_map, check.attributes = FALSE)),
    all(old_fit$scores$level == fit$scores$level[score_index]),
    max(abs(old_fit$scores$score - fit$scores$score[score_index])) < 1e-10,
    setequal(relaxed_edge_key(reference_sif), relaxed_edge_key(solution$SIF)),
    !anyNA(att_index), all(reference_att$level == solution$ATT$level[att_index]),
    max(abs(reference_att$score - solution$ATT$score[att_index])) < 1e-10,
    nrow(solution$SIF) == 5L, nrow(solution$ATT) == 6L
  )
  list(
    passed = TRUE,
    exact_signed_edges = TRUE,
    exact_score_node_ids_and_levels = TRUE,
    exact_compression_map = TRUE,
    tolerance = 1e-10,
    max_score_difference = max(abs(old_fit$scores$score - fit$scores$score[score_index])),
    max_selected_att_score_difference = max(abs(reference_att$score - solution$ATT$score[att_index]))
  )
}

relaxed_sif_hash <- function(sif) {
  path <- tempfile("relaxed_sif_", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(sort(relaxed_edge_key(sif)), path)
  unname(tools::md5sum(path))
}

relaxed_write_graphml <- function(solution, path) {
  graph <- igraph::graph_from_data_frame(solution$SIF[, c("source", "target", "interaction", "consistency"), drop = FALSE],
                                         directed = TRUE, vertices = solution$ATT)
  igraph::write_graph(graph, path, format = "graphml")
}

relaxed_self_test <- function() {
  scores <- data.frame(
    source = c("Metab__candidate", "Gene__middle", "Gene__tf"),
    score = c(-1, -1, -1), level = c(2, 1, 0), stringsAsFactors = FALSE
  )
  pkn <- data.frame(
    source = c("Metab__candidate", "Gene__middle", "Metab__candidate"),
    target = c("Gene__middle", "Gene__tf", "Gene__tf"),
    interaction = c(1, 1, 1), stringsAsFactors = FALSE
  )
  good <- relaxed_extract_score_layer_paths(scores, pkn, c(Metab__candidate = -1), 0.5, 0.5)
  stopifnot(nrow(good$SIF) == 2L,
            all(good$SIF$source != "Metab__candidate" | good$SIF$target != "Gene__tf"))
  bad <- relaxed_extract_score_layer_paths(scores, pkn, c(Metab__candidate = 1), 0.5, 0.5)
  stopifnot(!nrow(bad$SIF))
  cat("PASS: relaxed strict score-layer extraction rejects a level-skipping shortcut and wrong input sign.\n")
  invisible(TRUE)
}
