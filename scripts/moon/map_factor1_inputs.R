#!/usr/bin/env Rscript
# Run after prepare_factor1_inputs.R and audit_factor1_model.py.
source("scripts/support_functions.R")
out_dir <- "results/cosmos_factor1_pilot/inputs"
pkn_dir <- "results/cosmos_mouse_pkn"
raw <- readRDS(file.path(out_dir, "factor1_inputs_raw.rds"))
coords <- read.csv(file.path(out_dir, "factor1_sample_coordinates.csv"), stringsAsFactors = FALSE)
group_means <- tapply(coords$factor1_coordinate_saved, coords$gonadal_type, mean)
multiplier <- if (group_means[["Testes"]] > group_means[["Ovaries"]]) 1 else -1
nodes <- read.csv(file.path(pkn_dir, "nodes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
pkn <- readRDS(file.path(pkn_dir, "pkn_lipid_rules.rds"))
present <- unique(c(pkn$source, pkn$target))
nodes <- nodes[nodes$node_id %in% present, ]
crosswalk <- read.csv(file.path(pkn_dir, "mapping_features.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(names(raw$lipid_weight), crosswalk$feature_id))

mapping <- do.call(rbind, lapply(names(raw$lipid_weight), function(feature) {
  info <- crosswalk[crosswalk$feature_id == feature, ]
  measured_node <- paste0("Metab__MeasuredLipid:", feature, "_unknown")
  ids <- if (!is.na(info$exact_chebi_ids) && nzchar(info$exact_chebi_ids)) unique(strsplit(info$exact_chebi_ids, "\\|")[[1]]) else character()
  chosen <- character()
  method <- "excluded"
  reason <- "no_exact_measured_identity_in_pkn"
  if (measured_node %in% present) {
    chosen <- measured_node
    method <- "measured_species_rule_node"
    reason <- ""
  } else if (length(ids) == 1) {
    chosen <- nodes$node_id[nodes$node_type == "metabolite" & nodes$chemical_id == ids]
    if (length(chosen)) {
      method <- "unambiguous_exact_reported_resolution_chebi"
      reason <- ""
    }
  } else if (length(ids) > 1) {
    reason <- "multiple_exact_chebi_identities"
  }
  if (!length(chosen)) chosen <- NA_character_
  idx <- match(chosen, nodes$node_id)
  data.frame(feature = feature, node_id = chosen, mapping_method = method,
             included = !is.na(chosen), exclusion_reason = reason,
             raw_weight = unname(raw$lipid_weight[[feature]]),
             oriented_weight = unname(raw$lipid_weight[[feature]]) * multiplier,
             in_previous_sparse_network_64 = raw$lipid_table$in_previous_sparse_network_64[match(feature, raw$lipid_table$feature)],
             chemical_id = nodes$chemical_id[idx], compartment = nodes$compartment[idx],
             compartment_is_measured = FALSE, nodes_per_measured_feature = sum(!is.na(chosen)),
             compartment_policy = if (sum(!is.na(chosen)) > 1) "same_unlocalized_measurement_copied_to_alternative_compartment_nodes" else "unlocalized_measurement",
             exact_chebi_alternatives = info$exact_chebi_ids,
             resolution = info$nomenclature_level,
             prior_mapping_ambiguity = info$ambiguity, stringsAsFactors = FALSE)
}))
# No feature averaging: distinct observations that point to the same chemical
# node are left unresolved instead of silently replacing one another.
by_node <- split(mapping$feature[mapping$included], mapping$node_id[mapping$included])
duplicate_nodes <- names(by_node)[vapply(by_node, function(x) length(unique(x)) > 1, logical(1))]
conflict <- mapping$included & mapping$node_id %in% duplicate_nodes
mapping$included[conflict] <- FALSE
mapping$exclusion_reason[conflict] <- "multiple_measured_features_for_same_node"
write.csv(mapping, file.path(out_dir, "factor1_lipid_node_mapping.csv"), row.names = FALSE)
included <- mapping[mapping$included, ]
stopifnot(!anyDuplicated(included$node_id), all(is.finite(included$oriented_weight)))
tf_regulon <- read.csv(file.path(out_dir, "hepatocyte_tf_inference_regulons.csv"), stringsAsFactors = FALSE)
tf_regulon$source <- paste0("Gene__", tf_regulon$source)
tf_regulon$target <- paste0("Gene__", tf_regulon$target)
tf_activity <- setNames(raw$tf_activity * multiplier, paste0("Gene__", names(raw$tf_activity)))
hepatocyte_gene_weight <- setNames(raw$hepatocyte_gene_weight * multiplier, paste0("Gene__", names(raw$hepatocyte_gene_weight)))
lipid_weight <- setNames(included$oriented_weight, included$node_id)
tf_statistics <- raw$tf_statistics
tf_statistics$node_id <- paste0("Gene__", tf_statistics$source)
tf_statistics$oriented_activity <- tf_statistics$saved_factor1_activity * multiplier
tf_statistics$in_pkn <- tf_statistics$node_id %in% present
write.csv(tf_statistics, file.path(out_dir, "factor1_hepatocyte_tf_activities_mapped.csv"), row.names = FALSE)
write.csv(tf_regulon, file.path(out_dir, "hepatocyte_tf_regulon_pkn_ids.csv"), row.names = FALSE)
orientation <- list(label = "testes_positive_Factor1_coordinate", multiplier = multiplier,
                    saved_coordinate_group_means = group_means,
                    scope = "same_multiplier_for_TF_activities_RNA_loadings_lipid_loadings_and_sample_coordinates",
                    interpretation = "positive_agrees_with_testes_associated_Factor1_not_a_direct_differential_effect")
saveRDS(list(tf_activity = tf_activity, lipid_weight = lipid_weight,
             hepatocyte_gene_weight = hepatocyte_gene_weight, tf_regulon = tf_regulon,
             mapping = mapping, orientation = orientation, tf_statistics = tf_statistics,
             original_tf_activity = raw$tf_activity, original_lipid_weight = raw$lipid_weight,
             original_hepatocyte_gene_weight = raw$hepatocyte_gene_weight,
             source_manifest = raw$source_manifest), file.path(out_dir, "factor1_inputs.rds"))
cat("Mapped lipid nodes:", length(lipid_weight), "; unique measured lipid features:", length(unique(included$feature)), "\n")
print(table(included$mapping_method))
cat("TF activities:", length(tf_activity), "; TFs present in PKN:", sum(names(tf_activity) %in% present), "\n")
cat("Orientation multiplier:", multiplier, "; saved coordinate group means:\n")
print(group_means)
cat("Multiple-node compartment copies:\n")
print(included[included$nodes_per_measured_feature > 1, c("feature", "node_id", "raw_weight")])
