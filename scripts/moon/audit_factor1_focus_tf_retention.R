#!/usr/bin/env Rscript
# Read-only diagnostic of where focus TFs are lost. No production state changes.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
options(digits = 17, stringsAsFactors = FALSE)
root <- "results/cosmos_factor1_parameter_search"
out <- "results/cosmos_factor1_compressed_review/tf_audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
ids <- c("Gene__Esr1", "Gene__Nfic", "Gene__Dbp")
input <- readRDS("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
raw <- readRDS("results/cosmos_mouse_pkn/pkn_lipid_rules.rds")
catalog <- utils::read.csv("results/cosmos_mouse_pkn/nodes.csv")
nodes <- function(e) unique(c(e$source, e$target))
ekey <- function(e) paste(e$source, e$target, e$interaction, sep = "\t")
score_unfiltered <- function(pkn, state) {
  x <- as.data.frame(cosmosR::moon(upstream_input = NULL, downstream_input = state$lipid,
       meta_network = pkn, n_layers = state$n_layers, statistic = "ulm"))
  stopifnot(!anyDuplicated(x$source), all(is.finite(x$score)))
  x
}
paths <- list.files(file.path(root, "runs"), pattern = "^run_state\\.rds$", recursive = TRUE, full.names = TRUE)
stopifnot(length(paths) > 0)
rows <- list(); support <- list(); background <- list(); edges <- list(); all_tf <- list()
for (path in paths) {
  state <- readRDS(path)
  if (!((state$topology_radius == 6L && state$compression == "compressed") ||
        (state$topology_radius %in% c(6L, 7L, 12L) && state$n_layers %in% c(4L, 5L, 6L) && state$compression == "uncompressed"))) next
  cat("Audit:", state$run_id, "\n")
  prepared <- readRDS(state$topology_reference)
  cache <- readRDS(file.path(dirname(state$topology_reference), state$compression, "compression.rds"))
  final_pkn <- readRDS(file.path(dirname(path), "pkn_final_scoring.rds"))
  initial <- score_unfiltered(cache$pkn_scoring, state)
  final <- score_unfiltered(final_pkn, state)
  final_decompressed <- moon_safe_decompress_scores(final, cache$node_map)
  saved_match <- match(state$scores$source, final_decompressed$source)
  stopifnot(!anyNA(saved_match), max(abs(state$scores$score - final_decompressed$score[saved_match])) < 1e-10)
  saveRDS(list(initial = initial, final = final, final_decompressed = final_decompressed),
          file.path(out, paste0(state$run_id, "_diagnostic_scores.rds")))
  for (id in names(state$tf)) {
    fi <- match(id, final$source); si <- match(id, state$scores$source)
    all_tf[[length(all_tf) + 1L]] <- data.frame(run_id = state$run_id, source = id,
      activity = unname(state$tf[id]), raw_score = final$score[fi], level = final$level[fi],
      sign_compatible = is.finite(final$score[fi]) && sign(final$score[fi]) == sign(state$tf[id]),
      present_in_saved_scores = !is.na(si))
  }
  removed <- utils::read.csv(file.path(dirname(path), "coherence_removed_edges.csv"))
  for (id in ids) {
    fi <- match(id, final$source); ii <- match(id, initial$source)
    mapped <- cache$node_map$compressed[match(id, cache$node_map$original)]
    r <- data.frame(run_id = state$run_id, compression = state$compression,
      topology_radius = state$topology_radius, requested_layers = state$n_layers,
      source = id, activity = unname(input$tf_activity[id]),
      in_raw_pkn = id %in% nodes(raw), raw_outdegree = sum(raw$source == id),
      passes_activity_cutoff = abs(input$tf_activity[id]) > 2,
      in_pruned_pkn = id %in% nodes(prepared$net), pruned_outdegree = sum(prepared$net$source == id),
      in_upstream_input = id %in% names(state$tf),
      compression_protected = id %in% cache$protected_nodes,
      compressed_id_unchanged = !is.na(mapped) && mapped == id,
      in_final_scoring_pkn = id %in% nodes(final_pkn),
      final_outdegree = sum(final_pkn$source == id),
      outgoing_edges_removed_by_coherence = sum(removed$source == id),
      incoming_edges_removed_by_coherence = sum(removed$target == id),
      initial_score = initial$score[ii], initial_level = initial$level[ii],
      raw_final_score = final$score[fi], final_level = final$level[fi],
      initial_sign_compatible = is.finite(initial$score[ii]) && sign(initial$score[ii]) == sign(input$tf_activity[id]),
      final_sign_compatible = is.finite(final$score[fi]) && sign(final$score[fi]) == sign(input$tf_activity[id]),
      in_saved_scores = id %in% state$scores$source)
    r$exclusion_reason <- if (!r$in_upstream_input) "input_or_topology_filter" else if (!r$in_final_scoring_pkn) "coherence_edge_pruning" else if (!is.finite(r$raw_final_score)) "not_scored_within_realized_layers" else if (!r$final_sign_compatible) "post_scoring_upstream_sign_conflict" else "retained"
    rows[[length(rows) + 1L]] <- r
    if (!is.finite(r$raw_final_score)) next
    prev <- if (r$final_level == 1) data.frame(source = names(state$lipid), score = unname(state$lipid), level = 0) else final[final$level == r$final_level - 1, ]
    active <- final_pkn[final_pkn$source == id & final_pkn$target %in% prev$source, , drop = FALSE]
    coef <- stats::setNames(rep(0, nrow(prev)), prev$source)
    coef[active$target] <- active$interaction
    corr <- stats::cor(coef, prev$score)
    reconstructed <- corr * sqrt((nrow(prev) - 2) / ((1-corr+1e-20)*(1+corr+1e-20)))
    stopifnot(abs(reconstructed - r$raw_final_score) < 1e-10)
    cov_terms <- (coef - mean(coef)) * (prev$score - mean(prev$score))
    background[[length(background) + 1L]] <- data.frame(run_id = state$run_id, source = id,
      score = r$raw_final_score, level = r$final_level,
      layer_background_nodes = nrow(prev), scored_children = nrow(active),
      positive_children = sum(active$interaction > 0), negative_children = sum(active$interaction < 0),
      background_mean_score = mean(prev$score), target_mean_score = mean(prev$score[coef != 0]),
      nontarget_mean_score = mean(prev$score[coef == 0]), reconstructed_score = reconstructed,
      target_covariance_contribution = sum(cov_terms[coef != 0]),
      nontarget_covariance_contribution = sum(cov_terms[coef == 0]))
    for (j in seq_len(nrow(prev))) {
      id2 <- prev$source[j]
      # Leave-one-layer-feature-out influence is descriptive, no graph rescoring.
      cr <- if (nrow(prev) > 3) stats::cor(coef[-j], prev$score[-j]) else NA_real_
      loo <- cr * sqrt((nrow(prev) - 3) / ((1-cr+1e-20)*(1+cr+1e-20)))
      support[[length(support) + 1L]] <- data.frame(run_id = state$run_id, source = id,
        tf_score = r$raw_final_score, target = id2, target_label = catalog$label[match(id2, catalog$node_id)],
        target_type = catalog$node_type[match(id2, catalog$node_id)],
        target_score = prev$score[j], target_level = prev$level[j], coefficient = unname(coef[j]),
        is_scoring_child = coef[j] != 0, covariance_contribution = unname(cov_terms[j]),
        score_without_feature = loo, leave_one_feature_out_score_change = loo - r$raw_final_score)
    }
    for (stage in c("raw", "pruned", "final")) {
      net <- switch(stage, raw = raw, pruned = prepared$net, final = state$pkn_final)
      ee <- net[net$source == id, , drop = FALSE]
      if (!nrow(ee)) next
      ee$run_id <- state$run_id; ee$stage <- stage
      ee$target_label <- catalog$label[match(ee$target, catalog$node_id)]
      ee$target_type <- catalog$node_type[match(ee$target, catalog$node_id)]
      ee$target_RNA_weight <- unname(input$hepatocyte_gene_weight[ee$target])
      ee$tf_activity <- unname(input$tf_activity[id])
      ee$prescore_coherent_with_edge_sign <- sign(ee$interaction * ee$target_RNA_weight * ee$tf_activity) >= 0
      reg_key <- paste(input$tf_regulon$source, input$tf_regulon$target, input$tf_regulon$mor, sep = "\t")
      ee$matches_regulon_signed_pair <- ekey(ee) %in% reg_key
      edges[[length(edges) + 1L]] <- ee
    }
  }
}
audit <- do.call(rbind, rows)
write_csv(audit, "focus_tf_stage_audit.csv")
write_csv(do.call(rbind, all_tf), "all_upstream_tf_sign_audit.csv")
write_csv(do.call(rbind, background), "focus_tf_score_background.csv")
write_csv(do.call(rbind, support), "focus_tf_layer_influence.csv")
write_csv(do.call(rbind, edges), "focus_tf_outgoing_edges_by_stage.csv")
saveRDS(list(stages = audit, all_tf = do.call(rbind, all_tf), backgrounds = do.call(rbind, background),
             influences = do.call(rbind, support)), file.path(out, "audit_full_precision.rds"))
capture.output(cosmosR::moon, file = file.path(out, "installed_moon_source.txt"))
capture.output(cosmosR::filter_incohrent_TF_target, file = file.path(out, "installed_coherence_source.txt"))
capture.output(utils::sessionInfo(), file = file.path(out, "session_info.txt"))
print(audit[, c("run_id", "source", "initial_score", "raw_final_score", "final_level", "exclusion_reason")], row.names = FALSE)
cat("TF audit completed:", nrow(audit), "TF/run comparisons.\n")
