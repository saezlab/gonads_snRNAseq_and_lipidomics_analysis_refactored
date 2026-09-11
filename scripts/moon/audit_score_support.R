#!/usr/bin/env Rscript
# Diagnostic scoring only. Main pilot files are read, never overwritten.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
out <- "results/cosmos_factor1_pilot"
write_audit <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
pkn <- readRDS(file.path(out, "pkn_final_uncompressed.rds"))
compressed <- readRDS(file.path(out, "pkn_final_compressed.rds"))
mapping <- utils::read.csv(file.path(out, "compression_mapping.csv"))
main <- utils::read.csv(file.path(out, "moon_scores.csv"))
tf <- utils::read.csv(file.path(out, "upstream_input.csv"))
lipid <- utils::read.csv(file.path(out, "downstream_input.csv"))
catalog <- utils::read.csv("results/cosmos_mouse_pkn/nodes.csv")
# CSV scores have limited printed precision; use exact RDS values for the repeat
# ULM call so floating-point perturbations of near-zero scores do not change signs.
exact_input <- readRDS(file.path(out, "inputs/factor1_inputs.rds"))
tf$score <- unname(exact_input$tf_activity[tf$source])
lipid$score <- unname(exact_input$lipid_weight[lipid$source])
tf_input <- stats::setNames(tf$score, tf$source)
lipid_input <- stats::setNames(lipid$score, lipid$source)
score_map <- stats::setNames(main$score, main$source)
level_map <- stats::setNames(main$level, main$source)

# Endpoint support is exact at level 1: all input lipids form the ULM background.
ends <- pkn[pkn$target %in% lipid$source, , drop = FALSE]
by_parent <- split(seq_len(nrow(ends)), ends$source)
support <- do.call(rbind, lapply(names(by_parent), function(parent) {
  e <- ends[by_parent[[parent]], , drop = FALSE]
  cat_row <- catalog[match(parent, catalog$node_id), , drop = FALSE]
  coefs <- stats::setNames(rep(0, nrow(lipid)), lipid$source)
  coefs[e$target] <- e$interaction
  r <- stats::cor(coefs, lipid_input)
  reconstructed <- r * sqrt((nrow(lipid) - 2) / ((1-r+1e-20)*(1+r+1e-20)))
  data.frame(source = parent, compressed_source = mapping$compressed[match(parent, mapping$original)],
    label = cat_row$label, resource = cat_row$resource, catalysts = cat_row$catalysts,
    reaction_id = cat_row$reaction_id, observed_target_count = length(unique(e$target)),
    observed_targets = paste(sort(unique(e$target)), collapse = ";"),
    observed_coefficients = paste(e$interaction, collapse = ";"),
    observation_background_n = nrow(lipid), all_pkn_children = sum(pkn$source == parent),
    score = unname(score_map[parent]), level = unname(level_map[parent]),
    reconstructed_ulm_score = unname(reconstructed),
    sign_consistent_observed_edges = sum(sign(score_map[parent] * score_map[e$target]) == e$interaction),
    passes_primary_1_7 = abs(score_map[parent]) > 1.7,
    passes_secondary_1 = abs(score_map[parent]) > 1,
    stringsAsFactors = FALSE)
}))
stopifnot(all(support$level == 1),
          max(abs(support$score - support$reconstructed_ulm_score)) < 1e-10)
write_audit(support, "score_support_audit.csv")

# Identify immediate incoming dependencies of strong consistent endpoint enzymes.
strong_ends <- support$source[support$passes_primary_1_7 & support$sign_consistent_observed_edges > 0]
deps <- pkn[pkn$target %in% strong_ends, , drop = FALSE]
deps$source_score <- unname(score_map[deps$source])
deps$source_level <- unname(level_map[deps$source])
deps$target_score <- unname(score_map[deps$target])
deps$consistent <- sign(deps$source_score * deps$target_score) == deps$interaction
deps$source_passes_secondary_1 <- abs(deps$source_score) > 1
deps$source_label <- catalog$label[match(deps$source, catalog$node_id)]
deps$target_label <- catalog$label[match(deps$target, catalog$node_id)]
write_audit(deps, "strong_endpoint_predecessors.csv")

# Unconstrained scores expose TFs discarded solely for sign disagreement.
# The sole change in this diagnostic is whether same-child parents are merged.
# The frozen final coherence-filtered graphs are not filtered or refitted again.
cat("Scoring frozen final compressed graph without upstream sign filtering...\n")
score_c <- cosmosR::moon(upstream_input = NULL, downstream_input = lipid_input,
  meta_network = compressed, n_layers = 6, statistic = "ulm")
stopifnot(!anyDuplicated(score_c$source), all(is.finite(score_c$score)))
score_c <- moon_safe_decompress_scores(score_c, mapping)
cat("Scoring frozen final UNCOMPRESSED graph without upstream sign filtering...\n")
score_u <- cosmosR::moon(upstream_input = NULL, downstream_input = lipid_input,
  meta_network = pkn, n_layers = 6, statistic = "ulm")
stopifnot(!anyDuplicated(score_u$source), all(is.finite(score_u$score)))
write_audit(score_u, "uncompressed_sensitivity_scores.csv")
write_audit(score_c, "compressed_unconstrained_scores.csv")
shared_main <- match(main$source, score_c$source)
stopifnot(!anyNA(shared_main), max(abs(main$score - score_c$score[shared_main])) < 1e-10)
both <- merge(score_c[, c("source", "score", "level")],
              score_u[, c("source", "score", "level")], by = "source", all = TRUE,
              suffixes = c("_compressed", "_uncompressed"))
both$score_change <- both$score_uncompressed - both$score_compressed
both$same_sign <- sign(both$score_compressed) == sign(both$score_uncompressed)
both$same_level <- both$level_compressed == both$level_uncompressed
write_audit(both, "compression_score_comparison.csv")
tf_audit <- merge(tf, both, by = "source", all.x = TRUE, sort = FALSE)
names(tf_audit)[names(tf_audit) == "score"] <- "input_activity"
tf_audit$compressed_inferred <- is.finite(tf_audit$score_compressed)
tf_audit$uncompressed_inferred <- is.finite(tf_audit$score_uncompressed)
tf_audit$compressed_sign_compatible <- tf_audit$compressed_inferred &
  sign(tf_audit$input_activity) == sign(tf_audit$score_compressed)
tf_audit$uncompressed_sign_compatible <- tf_audit$uncompressed_inferred &
  sign(tf_audit$input_activity) == sign(tf_audit$score_uncompressed)
tf_audit <- tf_audit[order(tf_audit$source), ]
write_audit(tf_audit, "tf_unconstrained_compression_audit.csv")
level_summary <- do.call(rbind, lapply(sort(unique(both$level_compressed)), function(lev) {
  b <- both[both$level_compressed == lev & is.finite(both$score_uncompressed), , drop = FALSE]
  data.frame(level = lev, shared_nodes = nrow(b),
             score_correlation = stats::cor(b$score_compressed, b$score_uncompressed),
             sign_changes = sum(!b$same_sign), level_changes = sum(!b$same_level),
             median_abs_score_change = stats::median(abs(b$score_change)),
             max_abs_score_change = max(abs(b$score_change)))
}))
write_audit(level_summary, "compression_level_comparison.csv")
diag_extract <- list()
for (kind in c("compressed", "uncompressed")) {
  s <- if (kind == "compressed") score_c else score_u
  for (j in 1:2) {
    prim <- c(1, 1.7)[j]; sec <- c(0, 1)[j]
    sol <- moon_extract_paths(s, pkn, tf_input, prim, sec)
    diag_extract[[length(diag_extract) + 1L]] <- data.frame(scoring = kind,
      primary = prim, secondary = sec, nodes = nrow(sol$ATT), edges = nrow(sol$SIF),
      tfs = sum(sol$ATT$type == "upstream_input"), lipids = sum(sol$ATT$level == 0))
  }
}
write_audit(do.call(rbind, diag_extract), "compression_extraction_comparison.csv")
cat("Endpoint support count distribution:\n"); print(table(support$observed_target_count))
cat("Strong consistent endpoint predecessors:\n"); print(deps)
cat("Compression level comparison:\n"); print(level_summary)
cat("TFs of interest:\n"); print(tf_audit[tf_audit$source %in% c("Gene__Esr1", "Gene__Dbp", "Gene__Nfic"), ])
cat("TF sign compatibility:", sum(tf_audit$compressed_sign_compatible), "compressed;",
    sum(tf_audit$uncompressed_sign_compatible), "uncompressed; overlap",
    sum(tf_audit$compressed_sign_compatible & tf_audit$uncompressed_sign_compatible), "\n")
print(do.call(rbind, diag_extract))
