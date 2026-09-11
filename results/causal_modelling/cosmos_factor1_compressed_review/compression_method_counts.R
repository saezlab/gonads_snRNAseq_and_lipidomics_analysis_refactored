#!/usr/bin/env Rscript
# Read-only diagnostic of existing scoring graphs; no scoring or PKN edits.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
out <- "results/cosmos_factor1_compressed_review"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base <- "results/cosmos_factor1_parameter_search/topology_radius_6"
p <- readRDS(file.path(base, "prepared_topology.rds"))
input <- readRDS("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
full <- readRDS(file.path(base, "compressed/compression.rds"))
none <- readRDS(file.path(base, "uncompressed/compression.rds"))
minimal <- moon_safe_compress(p$net, unique(c(names(p$tf), names(p$lipid))))
signatures <- function(e) {
  e <- e[order(e$source, e$target, e$interaction), ]
  vapply(split(seq_len(nrow(e)), e$source), function(i) {
    paste0(e$target[i], "\t", e$interaction[i], collapse = "\n")
  }, character(1))
}
duplicate_stats <- function(e, protected) {
  s <- signatures(e)
  groups <- split(names(s), s)
  groups <- groups[lengths(groups) > 1L]
  has_protected <- vapply(groups, function(g) any(g %in% protected), logical(1))
  c(duplicate_groups = length(groups), duplicate_parents = sum(lengths(groups)),
    groups_with_protected_parent = sum(has_protected),
    groups_without_protected_parent = sum(!has_protected))
}
objects <- list(uncompressed = none, existing_compressed = full,
                inputs_only_compression_diagnostic = minimal)
counts <- do.call(rbind, lapply(names(objects), function(n) {
  x <- objects[[n]]
  data.frame(mode = n, edges = nrow(x$pkn_scoring),
    nodes = length(unique(c(x$pkn_scoring$source, x$pkn_scoring$target))),
    merged_original_nodes = x$n_merged_nodes, merged_groups = x$n_groups,
    protected = length(x$protected_nodes),
    expanded_protected = length(x$expanded_protected_nodes),
    as.list(duplicate_stats(x$pkn_scoring, x$expanded_protected_nodes)))
}))
utils::write.csv(counts, file.path(out, "compression_method_counts.csv"), row.names = FALSE)
print(counts)
additional <- minimal$node_map[minimal$node_map$merged &
  !minimal$node_map$original %in% full$node_map$original[full$node_map$merged], ]
additional$is_regulon_source <- additional$original %in% input$tf_regulon$source
additional$is_regulon_target <- additional$original %in% input$tf_regulon$target
utils::write.csv(additional, file.path(out, "additional_nodes_merged_with_inputs_only_protection.csv"), row.names = FALSE)
cat("Additional merged:", nrow(additional), "; regulon sources:", sum(additional$is_regulon_source),
    "; regulon targets:", sum(additional$is_regulon_target), "\n")
final <- do.call(rbind, lapply(c(4, 5, 6, 7, 8, 10, 12), function(d) {
  path <- paste0("results/cosmos_factor1_parameter_search/runs/radius6_depth", d, "_compressed/pkn_final_scoring.rds")
  e <- readRDS(path)
  data.frame(depth = d, edges = nrow(e), nodes = length(unique(c(e$source, e$target))),
             as.list(duplicate_stats(e, full$expanded_protected_nodes)))
}))
utils::write.csv(final, file.path(out, "compression_final_graph_redundancy.csv"), row.names = FALSE)
print(final)
stopifnot(none$n_groups == 0L, identical(none$pkn_scoring, none$pkn_uncompressed),
  full$n_groups > 0L, all(full$protected_nodes %in% full$node_map$original[!full$node_map$merged]),
  nrow(minimal$pkn_scoring) < nrow(full$pkn_scoring))
