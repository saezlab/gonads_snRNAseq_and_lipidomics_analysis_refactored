# Compatibility helpers for cosmosR 1.21.2. No omic files are read here.
# Run `Rscript scripts/moon/moon_safe_helpers.R` for synthetic contract tests.

moon_check_pkn <- function(pkn) {
  stopifnot(is.data.frame(pkn), all(c("source", "target", "interaction") %in% names(pkn)))
  pkn <- as.data.frame(pkn[, c("source", "target", "interaction"), drop = FALSE])
  pkn$source <- as.character(pkn$source)
  pkn$target <- as.character(pkn$target)
  stopifnot(is.numeric(pkn$interaction), !anyNA(pkn),
            all(nzchar(pkn$source)), all(nzchar(pkn$target)),
            all(pkn$interaction %in% c(-1, 1)))
  pkn <- suppressMessages(cosmosR::meta_network_cleanup(pkn))
  pkn <- pkn[order(pkn$source, pkn$target, pkn$interaction),
             c("source", "target", "interaction"), drop = FALSE]
  rownames(pkn) <- NULL
  pkn
}

moon_safe_compress <- function(pkn, protected_nodes = character()) {
  pkn <- moon_check_pkn(pkn)
  nodes <- sort(unique(c(pkn$source, pkn$target)))
  protected_nodes <- intersect(unique(as.character(protected_nodes)), nodes)
  # Sorting canonicalizes signed children before calling the positional helper.
  # Protect an ENTIRE signature group if any of its parents is protected.
  by_parent <- split(seq_len(nrow(pkn)), pkn$source)
  signatures <- vapply(by_parent, function(i) {
    paste0("parent_of_", paste0(pkn$target[i], pkn$interaction[i], collapse = "_____"))
  }, character(1))
  protected_signatures <- unique(signatures[intersect(names(signatures), protected_nodes)])
  expanded_protected <- union(protected_nodes,
                              names(signatures)[signatures %in% protected_signatures])
  duplicated_group <- duplicated(signatures) | duplicated(signatures, fromLast = TRUE)
  eligible <- duplicated_group & !names(signatures) %in% expanded_protected
  mapping <- data.frame(original = nodes, compressed = nodes, merged = FALSE,
                        stringsAsFactors = FALSE)
  compression <- NULL
  scoring <- pkn
  if (any(eligible)) {
    protect_vector <- stats::setNames(rep(0, length(expanded_protected)), expanded_protected)
    compression <- cosmosR::compress_same_children(
      pkn, sig_input = protect_vector, metab_input = numeric())
    replacements <- compression$duplicated_signatures
    stopifnot(length(replacements) > 0,
              !any(names(replacements) %in% expanded_protected))
    # A virtual signature must not coincide with any existing original ID.
    stopifnot(!any(unname(replacements) %in% nodes))
    idx <- match(names(replacements), mapping$original)
    stopifnot(!anyNA(idx))
    mapping$compressed[idx] <- unname(replacements)
    mapping$merged[idx] <- TRUE
    scoring <- moon_check_pkn(compression$compressed_network)
  }
  mapped_protected <- mapping$compressed[match(protected_nodes, mapping$original)]
  scoring_nodes <- unique(c(scoring$source, scoring$target))
  stopifnot(identical(unname(mapped_protected), unname(protected_nodes)),
            all(protected_nodes %in% scoring_nodes),
            setequal(mapping$compressed, scoring_nodes))
  # Verify the complete original->virtual map reproduces the exported graph.
  map_vector <- stats::setNames(mapping$compressed, mapping$original)
  projected <- pkn
  projected$source <- unname(map_vector[projected$source])
  projected$target <- unname(map_vector[projected$target])
  stopifnot(identical(moon_check_pkn(projected), scoring))
  message("MOON compression: ", nrow(pkn), " -> ", nrow(scoring), " edges; ",
          sum(mapping$merged), " original nodes in ",
          length(unique(mapping$compressed[mapping$merged])), " merged groups; ",
          length(protected_nodes), " protected IDs verified.")
  list(pkn_scoring = scoring, pkn_uncompressed = pkn, node_map = mapping,
       compression = compression, protected_nodes = protected_nodes,
       expanded_protected_nodes = expanded_protected,
       n_merged_nodes = sum(mapping$merged),
       n_groups = length(unique(mapping$compressed[mapping$merged])))
}

moon_safe_decompress_scores <- function(scores, node_map) {
  stopifnot(all(c("source", "score", "level") %in% names(scores)),
            all(c("original", "compressed") %in% names(node_map)),
            !anyDuplicated(scores$source), !anyDuplicated(node_map$original),
            all(is.finite(scores$score)), all(is.finite(scores$level)),
            all(scores$source %in% node_map$compressed))
  # A complete map avoids the package helper's zero-compression/leaf-only bug.
  mapping <- node_map[node_map$compressed %in% scores$source,
                      c("original", "compressed"), drop = FALSE]
  idx <- match(mapping$compressed, scores$source)
  out <- data.frame(source = mapping$original, score = scores$score[idx],
                    level = scores$level[idx], compressed_source = mapping$compressed,
                    stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(out$source),
            setequal(out$compressed_source, scores$source))
  out[order(out$source), , drop = FALSE]
}

moon_extract_paths <- function(scores, pkn, upstream_input,
                               primary_thresh = 1.7, secondary_thresh = 1,
                               RNA_input = NULL) {
  stopifnot(all(c("source", "score", "level") %in% names(scores)),
            !anyDuplicated(scores$source), all(is.finite(scores$score)),
            all(is.finite(scores$level)), length(primary_thresh) == 1,
            length(secondary_thresh) == 1, is.finite(primary_thresh),
            is.finite(secondary_thresh), primary_thresh >= secondary_thresh,
            secondary_thresh >= 0, is.numeric(upstream_input),
            !is.null(names(upstream_input)), !anyDuplicated(names(upstream_input)),
            all(is.finite(upstream_input)))
  pkn <- moon_check_pkn(pkn)
  scores <- as.data.frame(scores[, c("source", "score", "level"), drop = FALSE])
  scores$source <- as.character(scores$source)
  observed <- scores$source[scores$level == 0 & scores$score != 0]
  primary <- scores$source[scores$level > 0 & abs(scores$score) > primary_thresh]
  secondary <- union(observed,
                     scores$source[scores$level > 0 & abs(scores$score) > secondary_thresh])
  score_map <- stats::setNames(scores$score, scores$source)
  # Reapply upstream sign agreement even when the caller supplies unfiltered scores.
  seed_idx <- match(names(upstream_input), scores$source)
  seed_ok <- !is.na(seed_idx)
  seed_ok[seed_ok] <- scores$level[seed_idx[seed_ok]] > 0 &
    sign(scores$score[seed_idx[seed_ok]]) == sign(upstream_input[seed_ok])
  seeds <- intersect(names(upstream_input)[seed_ok], secondary)
  wrong_upstream <- setdiff(intersect(names(upstream_input), scores$source), seeds)
  secondary <- setdiff(secondary, wrong_upstream)
  primary <- intersect(primary, secondary)
  net <- pkn[pkn$source %in% secondary & pkn$target %in% secondary &
               (pkn$source %in% primary | pkn$target %in% primary) &
               !pkn$source %in% observed, , drop = FALSE]
  initial_edges <- nrow(net)
  net$consistency <- sign(score_map[net$source] * score_map[net$target]) == net$interaction
  net <- net[net$consistency, , drop = FALSE]
  consistent_edges <- nrow(net)
  empty_result <- function() {
    list(SIF = data.frame(source = character(), target = character(),
                          interaction = numeric(), consistency = logical()),
         ATT = data.frame(source = character(), score = numeric(), level = numeric(),
                          type = character(), RNA_input = numeric()),
         audit = list(primary_thresh = primary_thresh, secondary_thresh = secondary_thresh,
                      observed_exempt_from_secondary = TRUE,
                      observed_eligible_as_primary = FALSE,
                      candidate_edges = initial_edges, consistent_edges = consistent_edges,
                      final_edges = 0L, final_nodes = 0L, retained_upstream = 0L,
                      retained_observed = 0L, all_edges_consistent = TRUE,
                      all_nodes_on_seed_to_observation_walk = TRUE,
                      observed_outgoing_edges = 0L))
  }
  if (!nrow(net) || !length(seeds) || !length(observed)) return(empty_result())
  g <- igraph::graph_from_data_frame(net[, c("source", "target")], directed = TRUE)
  present <- igraph::V(g)$name
  seeds <- intersect(seeds, present)
  sinks <- intersect(observed, present)
  if (!length(seeds) || !length(sinks)) return(empty_result())
  # A multi-seed BFS via a synthetic super-root avoids all-pairs distance matrices.
  reachable <- function(graph, seed_names, mode) {
    root <- "__MOON_TEMP_REACHABILITY_ROOT__"
    while (root %in% igraph::V(graph)$name) root <- paste0(root, "_")
    if (mode == "in") graph <- igraph::reverse_edges(graph, igraph::E(graph))
    graph <- igraph::add_vertices(graph, 1, name = root)
    edges <- as.vector(rbind(rep(root, length(seed_names)), seed_names))
    graph <- igraph::add_edges(graph, edges)
    reached <- igraph::subcomponent(graph, root, mode = "out")
    setdiff(igraph::as_ids(reached), root)
  }
  keep <- intersect(reachable(g, seeds, "out"), reachable(g, sinks, "in"))
  net <- net[net$source %in% keep & net$target %in% keep, , drop = FALSE]
  if (!nrow(net)) return(empty_result())
  nodes <- sort(unique(c(net$source, net$target)))
  att <- scores[match(nodes, scores$source), , drop = FALSE]
  att$type <- ifelse(att$source %in% seeds, "upstream_input",
                    ifelse(att$level == 0, "level0", "other"))
  if (!is.null(RNA_input)) {
    stopifnot(is.numeric(RNA_input), !is.null(names(RNA_input)),
              !anyDuplicated(names(RNA_input)))
    att$RNA_input <- unname(RNA_input[att$source])
  } else att$RNA_input <- NA_real_
  final_g <- igraph::graph_from_data_frame(net[, c("source", "target")], directed = TRUE)
  retained_seeds <- intersect(seeds, nodes)
  retained_observed <- intersect(observed, nodes)
  verified_keep <- intersect(reachable(final_g, retained_seeds, "out"),
                             reachable(final_g, retained_observed, "in"))
  stopifnot(setequal(nodes, verified_keep), all(net$consistency),
            !any(net$source %in% observed),
            all(sign(score_map[net$source] * score_map[net$target]) == net$interaction))
  rownames(net) <- rownames(att) <- NULL
  list(SIF = net, ATT = att,
       audit = list(primary_thresh = primary_thresh, secondary_thresh = secondary_thresh,
                    observed_exempt_from_secondary = TRUE, observed_eligible_as_primary = FALSE,
                    candidate_edges = initial_edges, consistent_edges = consistent_edges,
                    final_edges = nrow(net), final_nodes = nrow(att),
                    retained_upstream = length(retained_seeds),
                    retained_observed = length(retained_observed),
                    all_edges_consistent = TRUE,
                    all_nodes_on_seed_to_observation_walk = TRUE,
                    observed_outgoing_edges = 0L))
}

moon_safe_helpers_self_test <- function() {
  message("Running synthetic MOON compatibility tests...")
  # Correct handling of original column order and a protected parent LAST among
  # three duplicate signatures, while another unprotected group really merges.
  net <- data.frame(source = c("A", "B", "P", "X", "Y", "U", "U", "U", "C"),
                    interaction = c(1, 1, 1, 1, 1, 1, 1, 1, -1),
                    target = c("L", "L", "L", "Z", "Z", "P", "A", "X", "L"))
  z <- moon_safe_compress(net, c("U", "P", "L", "Z"))
  m <- stats::setNames(z$node_map$compressed, z$node_map$original)
  stopifnot(all(m[c("A", "B", "P", "C")] == c("A", "B", "P", "C")),
            m["X"] == m["Y"], m["X"] != "X", z$n_merged_nodes == 2L,
            any(z$pkn_scoring$source == "U" & z$pkn_scoring$target == m["X"]))
  c_nodes <- unique(z$node_map$compressed)
  s <- data.frame(source = c_nodes, score = seq_along(c_nodes), level = 1)
  d <- moon_safe_decompress_scores(s, z$node_map)
  stopifnot(setequal(d$source, z$node_map$original),
            d$score[d$source == "X"] == d$score[d$source == "Y"])
  # Canonical order yields the same mapping when input edges are permuted.
  zz <- moon_safe_compress(net[rev(seq_len(nrow(net))), ], c("U", "P", "L", "Z"))
  stopifnot(identical(z$node_map, zz$node_map))
  plain <- data.frame(source = c("A", "B"), target = c("B", "L"), interaction = 1)
  zp <- moon_safe_compress(plain, c("A", "L"))
  dp <- moon_safe_decompress_scores(data.frame(source = c("A", "B", "L"),
                                               score = 1, level = 2:0), zp$node_map)
  stopifnot(zp$n_merged_nodes == 0, setequal(dp$source, c("A", "B", "L")))
  # Low raw observations are retained; high raw observations cannot create primary
  # nodes. Dead-end cycles, disconnected coherent edges and measured outgoing
  # edges are removed independently of whether any edge has a sign conflict.
  ex <- data.frame(source = c("U", "A", "U", "B", "C", "U", "Q", "L",
                              "U", "D", "E", "X", "U", "I"),
                   target = c("A", "L", "B", "C", "L", "Q", "M", "K",
                              "D", "E", "D", "Y", "I", "L"),
                   interaction = 1)
  sc <- data.frame(source = c("U", "A", "B", "C", "L", "Q", "M", "K",
                              "D", "E", "X", "Y", "I"),
                   score = c(2.5, 2, 1.2, 2, 0.02, 1.2, 10, 3, 3, 3, 3, 3, -3),
                   level = c(3, 1, 2, 1, 0, 1, 0, 1, 2, 2, 2, 1, 1))
  out <- moon_extract_paths(sc, ex, c(U = 1))
  stopifnot(nrow(out$SIF) == 5L,
            setequal(out$ATT$source, c("U", "A", "B", "C", "L")),
            out$ATT$score[out$ATT$source == "L"] == 0.02,
            !any(c("Q", "M", "K", "D", "E", "X", "Y", "I") %in% out$ATT$source),
            out$audit$all_nodes_on_seed_to_observation_walk,
            out$audit$observed_outgoing_edges == 0L)
  coherent <- moon_extract_paths(sc[sc$source != "I", ],
                                  ex[!ex$source %in% "I" & !ex$target %in% "I", ], c(U = 1))
  stopifnot(identical(out$SIF, coherent$SIF))
  wrong_sign <- moon_extract_paths(sc, ex, c(U = -1))
  stopifnot(nrow(wrong_sign$SIF) == 0L)
  message("PASS: protected groups, canonical columns/order, complete decompression, ",
          "mixed score thresholds, sign coherence, measurement sinks, cycle/disconnection ",
          "pruning and empty upstream case.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) moon_safe_helpers_self_test()
