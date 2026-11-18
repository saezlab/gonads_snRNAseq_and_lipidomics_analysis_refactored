#!/usr/bin/env Rscript
# Run from the repository root. This extracts observed inputs only: no PKN node
# mapping, threshold selection, sign orientation, or MOON inference is performed.
source("scripts/support_functions.R")

out_dir <- "results/cosmos_factor1_pilot/inputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
factor_name <- "Factor1"
weights_file <- "results/mofa/main_model_RNA_lipid_weight.RData"
weights_env <- new.env(parent = emptyenv())
load(weights_file, envir = weights_env)
weights <- weights_env$weights
stopifnot(is.list(weights), all(c("Hepatocytes", "lipidomic", "Hepatocytes_TFs") %in% names(weights)))
stopifnot(all(vapply(weights, function(x) factor_name %in% colnames(x), logical(1))))
cat("Input class:", class(weights), "\n")
print(data.frame(view = names(weights), features = vapply(weights, nrow, integer(1)),
                 factors = vapply(weights, ncol, integer(1))))

all_weights <- do.call(rbind, lapply(names(weights), function(view) {
  values <- weights[[view]][, factor_name]
  data.frame(feature_id = names(values), view = view,
             feature = sub(paste0("_", view, "$"), "", names(values)),
             factor = factor_name, raw_value = as.numeric(values),
             measurement_type = if (view == "Hepatocytes_TFs") "TF_ULM_activity_score" else "MOFA_loading",
             stringsAsFactors = FALSE)
}))
stopifnot(all(is.finite(all_weights$raw_value)))
rownames(all_weights) <- NULL
write.csv(all_weights, file.path(out_dir, "factor1_all_views_raw.csv"), row.names = FALSE)
lipids <- all_weights[all_weights$view == "lipidomic", ]
names(lipids)[names(lipids) == "raw_value"] <- "factor1_weight"
cor_env <- new.env(parent = emptyenv())
load("results/mofa/networks/gene_lipid_correlations_factor1.RData", envir = cor_env)
prior_selected_lipids <- colnames(cor_env$hepatocyte_cor_mat_RNA_lip)
lipids$in_previous_sparse_network_64 <- lipids$feature %in% prior_selected_lipids
lipids$absolute_weight_rank <- rank(-abs(lipids$factor1_weight), ties.method = "min")
lipids$loading_p_value <- NA_real_ # MOFA loadings have no saved feature-level p-value.
write.csv(lipids, file.path(out_dir, "factor1_lipid_weights_all.csv"), row.names = FALSE)

# Recreate the exact input regulons used by workflow 04, preserving full ULM
# statistics instead of only the legacy score matrix. These are target-based TF
# activity estimates, never TF transcript expression.
stopifnot(requireNamespace("decoupleR", quietly = TRUE))
reg_env <- new.env(parent = emptyenv())
load("support/TF_regulons.RData", envir = reg_env)
regulons <- as.data.frame(reg_env$TF_regulons)
regulons <- regulons[complete.cases(regulons), c("source", "target", "mor")]
hormone <- read.csv("support/Gonadal_Hormone_Genes.csv", stringsAsFactors = FALSE)
hormone$source <- gsub("_spec", "", hormone$source)
hormone$source <- gsub("ESR1", "Esr1", hormone$source)
hormone$source <- gsub("AR", "Ar", hormone$source)
regulons <- rbind(regulons[!regulons$source %in% c("Ar", "Esr1"), ], hormone[, c("source", "target", "mor")])
rna <- weights$Hepatocytes[, factor_name, drop = FALSE]
rownames(rna) <- sub("_Hepatocytes$", "", rownames(rna))
cat("Hepatocyte loading matrix:", nrow(rna), "x", ncol(rna), "\n")
cat("Regulon rows:", nrow(regulons), "\n")
tf_statistics <- as.data.frame(decoupleR::run_ulm(mat = rna, network = regulons, minsize = 20))
tf_statistics$p_adj_BH <- p.adjust(tf_statistics$p_value, method = "BH")
tf_statistics$n_measured_targets <- vapply(tf_statistics$source, function(tf) {
  length(unique(regulons$target[regulons$source == tf & regulons$target %in% rownames(rna)]))
}, integer(1))
saved_scores <- weights$Hepatocytes_TFs[, factor_name]
names(saved_scores) <- sub("_Hepatocytes_TFs$", "", names(saved_scores))
tf_statistics$saved_factor1_activity <- unname(saved_scores[tf_statistics$source])
tf_statistics$score_difference_from_saved <- tf_statistics$score - tf_statistics$saved_factor1_activity
stopifnot(setequal(names(saved_scores), tf_statistics$source))
stopifnot(max(abs(tf_statistics$score_difference_from_saved)) < 1e-7)
tf_statistics <- tf_statistics[order(-abs(tf_statistics$saved_factor1_activity)), ]
write.csv(tf_statistics, file.path(out_dir, "factor1_hepatocyte_tf_activities.csv"), row.names = FALSE)
write.csv(regulons, file.path(out_dir, "hepatocyte_tf_inference_regulons.csv"), row.names = FALSE)

enrichment_env <- new.env(parent = emptyenv())
load("results/mofa/factor_analysis/factor_enrichment_gonadal.RData", envir = enrichment_env)
write.csv(enrichment_env$factor_enrichment, file.path(out_dir, "saved_factor_gonadal_enrichment.csv"), row.names = FALSE)
input_files <- c(weights_file, "support/TF_regulons.RData", "support/Gonadal_Hormone_Genes.csv",
                 "results/mofa/RNA/mofa_res_12factor.hdf5",
                 "results/mofa/networks/gene_lipid_correlations_factor1.RData")
manifest <- data.frame(path = input_files, md5 = unname(tools::md5sum(input_files)),
                       bytes = file.info(input_files)$size)
write.csv(manifest, file.path(out_dir, "input_source_manifest.csv"), row.names = FALSE)
saveRDS(list(factor = factor_name, tf_activity = saved_scores,
             lipid_weight = setNames(lipids$factor1_weight, lipids$feature),
             hepatocyte_gene_weight = setNames(as.numeric(rna[, 1]), rownames(rna)),
             tf_statistics = tf_statistics, lipid_table = lipids,
             orientation = "unchanged_saved_Factor1_coordinate",
             source_manifest = manifest), file.path(out_dir, "factor1_inputs_raw.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "preparation_session_info.txt"))
cat("All lipid weights:", nrow(lipids), "; previous sparse-network lipids:", sum(lipids$in_previous_sparse_network_64), "\n")
cat("TF activities:", nrow(tf_statistics), "; max reproduced score error:", max(abs(tf_statistics$score_difference_from_saved)), "\n")
print(head(tf_statistics[, c("source", "saved_factor1_activity", "p_value", "p_adj_BH", "n_measured_targets")], 12))
cat("Raw input extraction completed. Factor signs and thresholds remain unmodified.\n")
