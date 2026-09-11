#!/usr/bin/env Rscript
source("scripts/moon/moon_safe_helpers.R")
s <- readRDS("results/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/run_state.rds")
base_dir <- "results/cosmos_factor1_compressed_review/threshold_search"
reference <- readRDS(file.path(base_dir, "C03", "solution.rds"))
key <- function(e) paste(e$source, e$target, e$interaction, sep = "\t")
reference_keys <- key(reference$SIF)
grid <- expand.grid(primary = c(1.8, 1.85, 1.875), secondary = c(0.55, 0.575, 0.60))
rows <- lapply(seq_len(nrow(grid)), function(i) {
  p <- grid$primary[i]; z <- grid$secondary[i]
  sol <- moon_extract_paths(s$scores, s$pkn_final, s$tf, p, z)
  a <- sol$ATT; e <- sol$SIF
  tf <- a[a$type == "upstream_input", ]
  data.frame(primary = p, secondary = z, nodes = nrow(a), edges = nrow(e),
    upstream_tfs = nrow(tf), downstream_lipids = sum(a$level == 0),
    hnf4a_secondary_margin = abs(s$scores$score[s$scores$source == "Gene__Hnf4a"]) - z,
    min_inferred_node_secondary_margin = min(abs(a$score[a$level > 0]) - z),
    min_primary_node_margin = min(abs(a$score[a$level > 0 & abs(a$score) > p]) - p),
    identical_edge_set_to_C03 = setequal(reference_keys, key(e)),
    edge_jaccard_to_C03 = length(intersect(reference_keys, key(e))) / length(union(reference_keys, key(e))))
})
utils::write.csv(do.call(rbind, rows), file.path(base_dir, "C03_local_threshold_neighborhood.csv"), row.names = FALSE)
print(do.call(rbind, rows))

# Compare the same two illustrative cutoffs across supplied compressed states.
# This does not rescore any state or alter the threshold-search selections.
paths <- c("results/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/run_state.rds",
  unlist(lapply(c("inputs_only_sensitivity", "topology_sensitivity", "prescore_sensitivity"),
    function(z) list.files(file.path("results/cosmos_factor1_compressed_review", z, "runs"),
      pattern = "^run_state\\.rds$", full.names = TRUE, recursive = TRUE)), use.names = FALSE))
matched <- list()
for (path in paths) {
  state <- readRDS(path)
  for (setting in list(c(1.85, 0.55), c(2, 0.85))) {
    sol <- moon_extract_paths(state$scores, state$pkn_final, state$tf, setting[1], setting[2])
    a <- sol$ATT
    matched[[length(matched) + 1L]] <- data.frame(run_id = state$run_id,
      primary = setting[1], secondary = setting[2], nodes = nrow(a), edges = nrow(sol$SIF),
      upstream_tfs = sum(a$type == "upstream_input"), downstream_lipids = sum(a$level == 0),
      exact_edge_set_equals_baseline_C03 = setequal(reference_keys, key(sol$SIF)),
      edge_jaccard_to_baseline_C03 = length(intersect(reference_keys, key(sol$SIF))) /
        length(union(reference_keys, key(sol$SIF))),
      tf_ids = paste(sort(a$source[a$type == "upstream_input"]), collapse = ";"),
      lipid_ids = paste(sort(a$source[a$level == 0]), collapse = ";"))
  }
}
utils::write.csv(do.call(rbind, matched),
  "results/cosmos_factor1_compressed_review/sensitivity_threshold_search/matched_threshold_sensitivity.csv",
  row.names = FALSE)
