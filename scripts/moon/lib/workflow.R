# Shared, direction-neutral graph/scoring operations. Source after moon_safe_helpers.R.
# Input semantics and thresholds belong in the calling analysis configuration.

moon_nodes <- function(net) unique(c(net$source, net$target))
moon_edge_keys <- function(net) paste(net$source, net$interaction, net$target, sep = "\t")
moon_write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")

moon_vector_check <- function(x) {
  stopifnot(is.numeric(x), length(x) > 0L, !is.null(names(x)),
            !anyDuplicated(names(x)), all(nzchar(names(x))), all(is.finite(x)))
  invisible(x)
}

moon_read_config <- function(path) {
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)
  stopifnot(identical(cfg$analysis, "factor1_tf_to_lipid"),
            cfg$pruning_steps >= 1L, cfg$n_layers >= 1L,
            cfg$max_iterations >= 1L, cfg$primary > cfg$secondary,
            cfg$secondary >= 0, cfg$tf_cutoff > 0)
  for (pkg in names(cfg$packages)) {
    if (!requireNamespace(pkg, quietly = TRUE) ||
        as.character(utils::packageVersion(pkg)) != cfg$packages[[pkg]]) {
      stop("Selected-result replay requires ", pkg, " ", cfg$packages[[pkg]])
    }
  }
  sha <- utils::packageDescription("cosmosR")$RemoteSha
  if (!identical(sha, cfg$cosmos_commit)) stop("cosmosR commit differs from the selected analysis")
  for (name in names(cfg$input_md5)) {
    if (!file.exists(name) || unname(tools::md5sum(name)) != cfg$input_md5[[name]]) {
      stop("Missing or changed frozen dependency: ", name)
    }
  }
  cfg
}

moon_load_factor1 <- function(cfg) {
  input <- readRDS(cfg$inputs)
  invisible(lapply(input[c("tf_activity", "lipid_weight", "hepatocyte_gene_weight")], moon_vector_check))
  pkn <- moon_check_pkn(readRDS(cfg$pkn))
  upstream <- input$tf_activity[abs(input$tf_activity) > cfg$tf_cutoff &
                                  names(input$tf_activity) %in% moon_nodes(pkn)]
  downstream <- input$lipid_weight[input$lipid_weight != 0 &
                                    names(input$lipid_weight) %in% moon_nodes(pkn)]
  moon_vector_check(upstream)
  moon_vector_check(downstream)
  stopifnot(all(grepl("^Gene__", names(upstream))),
            all(grepl("^Metab__", names(downstream))),
            !length(intersect(names(upstream), names(downstream))))
  regulon <- input$tf_regulon[, c("source", "target", "mor")]
  stopifnot(!anyNA(regulon), all(regulon$mor %in% c(-1, 1)))
  list(pkn = pkn, upstream = upstream, downstream = downstream,
       rna = input$hepatocyte_gene_weight, regulon = regulon,
       orientation = input$orientation)
}

moon_reach <- function(net, seeds, steps = Inf, mode = "out") {
  if (mode == "in") net <- data.frame(source = net$target, target = net$source)
  seen <- intersect(seeds, moon_nodes(net))
  frontier <- seen
  step <- 0L
  while (length(frontier) && step < steps) {
    frontier <- setdiff(unique(net$target[net$source %in% frontier]), seen)
    seen <- c(seen, frontier)
    step <- step + 1L
  }
  seen
}

moon_prepare <- function(pkn, upstream, downstream, pruning_steps) {
  moon_vector_check(upstream)
  moon_vector_check(downstream)
  stopifnot(!length(intersect(names(upstream), names(downstream))))
  net <- pkn[!pkn$source %in% names(downstream), , drop = FALSE]
  ledger <- list()
  repeat {
    previous <- moon_edge_keys(net)
    forward <- moon_reach(net, names(upstream), pruning_steps)
    net <- net[net$source %in% forward & net$target %in% forward, , drop = FALSE]
    backward <- moon_reach(net, names(downstream), pruning_steps, "in")
    net <- net[net$source %in% backward & net$target %in% backward, , drop = FALSE]
    ids <- moon_nodes(net)
    upstream <- upstream[names(upstream) %in% ids]
    downstream <- downstream[names(downstream) %in% ids]
    if (!length(upstream) || !length(downstream) || !nrow(net)) {
      stop("No upstream-to-observation topology remains at pruning_steps=", pruning_steps)
    }
    ledger[[length(ledger) + 1L]] <- data.frame(iteration = length(ledger) + 1L,
      edges = nrow(net), nodes = length(ids), upstream_inputs = length(upstream),
      downstream_inputs = length(downstream))
    if (identical(previous, moon_edge_keys(net))) break
  }
  list(net = moon_check_pkn(net), upstream = upstream, downstream = downstream,
       ledger = do.call(rbind, ledger))
}

moon_score <- function(net, upstream, downstream, n_layers) {
  ans <- as.data.frame(cosmosR::moon(upstream_input = upstream,
    downstream_input = downstream, meta_network = net, n_layers = n_layers,
    statistic = "ulm"))[, c("source", "score", "level")]
  stopifnot(!anyDuplicated(ans$source), all(is.finite(ans$score)),
            setequal(names(downstream), ans$source[ans$level == 0]))
  ans
}

moon_fit <- function(data, pruning_steps, n_layers, max_iterations = 10L) {
  prepared <- moon_prepare(data$pkn, data$upstream, data$downstream, pruning_steps)
  protected <- unique(c(names(prepared$upstream), names(prepared$downstream),
                         data$regulon$source, data$regulon$target))
  compression <- moon_safe_compress(prepared$net, protected)
  scoring <- compression$pkn_scoring
  ledger <- removed <- list()
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    message("MOON iteration ", iteration, ": ", nrow(scoring), " compressed edges")
    scores <- moon_score(scoring, prepared$upstream, prepared$downstream, n_layers)
    filtered <- cosmosR::filter_incohrent_TF_target(scores, data$regulon, scoring, data$rna)
    stopifnot(all(moon_edge_keys(filtered) %in% moon_edge_keys(scoring)))
    delta <- scoring[!moon_edge_keys(scoring) %in% moon_edge_keys(filtered), , drop = FALSE]
    ledger[[iteration]] <- data.frame(iteration = iteration, edges_scored = nrow(scoring),
      edges_removed = nrow(delta), scored_nodes = nrow(scores),
      upstream_scored = sum(names(prepared$upstream) %in% scores$source),
      max_level = max(scores$level))
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
  if (!converged) stop("RNA coherence did not converge in ", max_iterations,
                       " iterations; no solution was selected")
  mapping <- stats::setNames(compression$node_map$compressed, compression$node_map$original)
  projected <- prepared$net
  projected$source <- unname(mapping[projected$source])
  projected$target <- unname(mapping[projected$target])
  final <- prepared$net[moon_edge_keys(projected) %in% moon_edge_keys(scoring), , drop = FALSE]
  stopifnot(setequal(moon_edge_keys(projected[moon_edge_keys(projected) %in%
                                               moon_edge_keys(scoring), ]), moon_edge_keys(scoring)))
  list(prepared = prepared, compression = compression, pkn_scoring = scoring,
       pkn_final = final, scores_compressed = scores,
       scores = moon_safe_decompress_scores(scores, compression$node_map),
       coherence_ledger = do.call(rbind, ledger),
       removed_edges = if (length(removed)) do.call(rbind, removed) else
         data.frame(source = character(), target = character(), interaction = numeric(), iteration = integer()),
       pruning_steps = pruning_steps, n_layers = n_layers, converged = converged)
}

moon_diagnose <- function(fit, data, solution) {
  # NULL removes only the package's final sign filter, never the scoring graph.
  raw <- moon_score(fit$pkn_scoring, NULL, fit$prepared$downstream, fit$n_layers)
  raw <- moon_safe_decompress_scores(raw, fit$compression$node_map)
  idx <- match(fit$scores$source, raw$source)
  stopifnot(!anyNA(idx), max(abs(fit$scores$score - raw$score[idx])) < 1e-10,
            all(fit$scores$level == raw$level[idx]))
  raw$upstream_activity <- unname(fit$prepared$upstream[raw$source])
  raw$input_tf <- raw$source %in% names(fit$prepared$upstream)
  raw$returned_by_moon <- raw$source %in% fit$scores$source
  raw$sign_agrees <- ifelse(raw$input_tf, sign(raw$score) == sign(raw$upstream_activity), NA)
  raw$retention_reason <- ifelse(raw$returned_by_moon, "returned",
    ifelse(raw$input_tf & !raw$sign_agrees, "upstream_sign_conflict", "other"))
  stopifnot(!any(raw$retention_reason == "other"))
  raw$in_solution <- raw$source %in% solution$ATT$source
  raw
}

moon_validate_solution <- function(solution, fit, primary, secondary) {
  a <- solution$ATT
  e <- solution$SIF
  sm <- stats::setNames(fit$scores$score, fit$scores$source)
  primary_nodes <- a$source[a$level > 0 & abs(a$score) > primary]
  stopifnot(all(moon_edge_keys(e) %in% moon_edge_keys(fit$pkn_final)),
            all(sign(sm[e$source] * sm[e$target]) == e$interaction),
            all(abs(a$score[a$level > 0]) > secondary),
            all(e$source %in% primary_nodes | e$target %in% primary_nodes),
            !any(e$source %in% names(fit$prepared$downstream)),
            solution$audit$all_nodes_on_seed_to_observation_walk)
  invisible(TRUE)
}

moon_compare_reference <- function(fit, solution, cfg, all_scores) {
  old <- readRDS(file.path(cfg$reference_run, "run_state.rds"))
  old_compression <- readRDS(cfg$reference_compression)
  old_sif <- utils::read.csv(file.path(cfg$reference_solution, "network_SIF.csv"))
  old_all <- utils::read.csv(cfg$reference_all_scores)
  si <- match(old$scores$source, fit$scores$source)
  ai <- match(old_all$source, all_scores$source)
  stopifnot(setequal(old$scores$source, fit$scores$source), !anyNA(si),
    setequal(moon_edge_keys(old$pkn_final), moon_edge_keys(fit$pkn_final)),
    setequal(moon_edge_keys(old_sif), moon_edge_keys(solution$SIF)),
    identical(old_compression$node_map, fit$compression$node_map),
    all(old$scores$level == fit$scores$level[si]),
    max(abs(old$scores$score - fit$scores$score[si])) < 1e-10,
    setequal(old_all$source, all_scores$source), !anyNA(ai),
    max(abs(old_all$score - all_scores$score[ai])) < 1e-10,
    nrow(solution$SIF) == 388L, nrow(solution$ATT) == 153L,
    nrow(fit$coherence_ledger) == 3L,
    solution$audit$retained_upstream == 4L, solution$audit$retained_observed == 3L)
  list(passed = TRUE, exact_signed_edges = TRUE, exact_score_node_ids_and_levels = TRUE,
       exact_compression_map = TRUE, tolerance = 1e-10,
       max_returned_score_difference = max(abs(old$scores$score - fit$scores$score[si])),
       max_diagnostic_score_difference = max(abs(old_all$score - all_scores$score[ai])))
}

moon_save_fit <- function(fit, out) {
  saveRDS(fit, file.path(out, "run_state.rds"))
  saveRDS(fit$prepared, file.path(out, "prepared_topology.rds"))
  saveRDS(fit$compression, file.path(out, "compression.rds"))
  saveRDS(fit$pkn_scoring, file.path(out, "pkn_final_scoring.rds"))
  saveRDS(fit$pkn_final, file.path(out, "pkn_final_uncompressed.rds"))
  moon_write_csv(fit$scores, file.path(out, "moon_scores.csv"))
  moon_write_csv(fit$prepared$ledger, file.path(out, "pruning_ledger.csv"))
  moon_write_csv(fit$coherence_ledger, file.path(out, "coherence_ledger.csv"))
  moon_write_csv(fit$removed_edges, file.path(out, "coherence_removed_edges.csv"))
  moon_write_csv(fit$compression$node_map, file.path(out, "compression_mapping.csv"))
}

moon_new_output <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; choose a new --output= directory: ", path)
  }
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(path)) {
    stop("Cannot create output directory: ", path)
  }
  normalizePath(path)
}
