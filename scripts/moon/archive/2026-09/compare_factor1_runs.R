#!/usr/bin/env Rscript
# Compare independently converged runs; write diagnostic artifacts only.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
out <- "results/cosmos_factor1_pilot"
uncompressed_dir <- file.path(out, "uncompressed_sensitivity")
exact <- readRDS(file.path(out, "inputs/factor1_inputs.rds"))
read_table <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)
write_table <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
tf_ids <- read_table(file.path(out, "upstream_input.csv"))$source
tf_values <- exact$tf_activity[tf_ids]
stopifnot(!anyNA(tf_values), !anyDuplicated(tf_ids),
          setequal(tf_ids, read_table(file.path(uncompressed_dir, "upstream_input.csv"))$source))

score_unconstrained <- function(run_dir, diagnostic_name) {
  net <- readRDS(file.path(run_dir, "pkn_final_compressed.rds"))
  map <- read_table(file.path(run_dir, "compression_mapping.csv"))
  observed <- read_table(file.path(run_dir, "downstream_input.csv"))$source
  values <- exact$lipid_weight[observed]
  stopifnot(!anyNA(values), !anyDuplicated(observed), all(is.finite(values)))
  cat("Unconstrained final-graph diagnostic:", run_dir, "\n")
  s <- cosmosR::moon(upstream_input = NULL, downstream_input = values,
                    meta_network = net, n_layers = 6, statistic = "ulm")
  stopifnot(!anyDuplicated(s$source), all(is.finite(s$score)))
  s <- moon_safe_decompress_scores(s, map)
  main <- read_table(file.path(run_dir, "moon_scores.csv"))
  idx <- match(main$source, s$source)
  stopifnot(!anyNA(idx), max(abs(main$score - s$score[idx])) < 1e-10,
            all(main$level == s$level[idx]))
  write_table(s, file.path(run_dir, diagnostic_name))
  s
}

# These calls expose TF values omitted by moon's post-scoring upstream sign filter.
# Neither PKN nor any main scoring/extraction output is changed or refiltered.
sc <- score_unconstrained(out, "fullcomparison_compressed_unconstrained.csv")
su <- score_unconstrained(uncompressed_dir, "moon_scores_unconstrained_audit.csv")
tc <- sc$score[match(tf_ids, sc$source)]
tu <- su$score[match(tf_ids, su$source)]
tf <- data.frame(source = tf_ids, input_activity = unname(tf_values),
                 compressed_score = tc, uncompressed_score = tu,
                 compressed_inferred = is.finite(tc), uncompressed_inferred = is.finite(tu),
                 compressed_sign_compatible = is.finite(tc) & sign(tc) == sign(tf_values),
                 uncompressed_sign_compatible = is.finite(tu) & sign(tu) == sign(tf_values),
                 compressed_near_zero = is.finite(tc) & abs(tc) < 1e-10,
                 uncompressed_near_zero = is.finite(tu) & abs(tu) < 1e-10,
                 row.names = NULL)
tf <- tf[order(tf$source), ]
write_table(tf, file.path(out, "fullcomparison_tf.csv"))
ct <- tf$source[tf$compressed_sign_compatible]
ut <- tf$source[tf$uncompressed_sign_compatible]
stopifnot(setequal(ct, intersect(tf_ids, read_table(file.path(out, "moon_scores.csv"))$source)),
          setequal(ut, intersect(tf_ids, read_table(file.path(uncompressed_dir, "moon_scores.csv"))$source)))
overlap <- function(a, b) {
  a <- unique(a); b <- unique(b)
  total <- length(union(a, b))
  list(compressed_count = length(a), uncompressed_count = length(b),
       shared_count = length(intersect(a, b)), union_count = total,
       jaccard = if (total) length(intersect(a, b)) / total else NA_real_)
}
network_sets <- function(run_dir, setting) {
  e <- read_table(file.path(run_dir, paste0("solution_", setting, "_SIF.csv")))
  a <- read_table(file.path(run_dir, paste0("solution_", setting, "_ATT.csv")))
  list(edges = paste(e$source, e$target, e$interaction, sep = "\t"),
       nodes = a$source, upstream_tfs = a$source[a$type == "upstream_input"],
       lipids = a$source[a$level == 0])
}
records <- list()
for (setting in c("relaxed", "primary", "strict")) {
  a <- network_sets(out, setting); b <- network_sets(uncompressed_dir, setting)
  for (entity in names(a)) records[[length(records) + 1L]] <- data.frame(
    setting = setting, entity = entity, as.data.frame(overlap(a[[entity]], b[[entity]])))
}
network_overlap <- do.call(rbind, records)
write_table(network_overlap, file.path(out, "fullcomparison_overlap.csv"))
paths <- c(file.path(out, "pkn_final_uncompressed.rds"),
           file.path(uncompressed_dir, "pkn_final_uncompressed.rds"))
summary <- list(
  comparison = paste("Two independently converged ULM runs from the same initial pruned PKN:",
                     "protected same-child compression versus no compression.",
                     "Final TF diagnostic scores are before sign filtering."),
  command = "Rscript scripts/moon/compare_factor1_runs.R",
  exact_input_source = file.path(out, "inputs/factor1_inputs.rds"),
  compressed_iterations = nrow(read_table(file.path(out, "coherence_ledger.csv"))),
  uncompressed_iterations = nrow(read_table(file.path(uncompressed_dir, "coherence_ledger.csv"))),
  tf_inferred = overlap(tf$source[tf$compressed_inferred], tf$source[tf$uncompressed_inferred]),
  tf_sign_compatible = overlap(ct, ut),
  shared_sign_compatible_tfs = sort(intersect(ct, ut)),
  compressed_only_sign_compatible_tfs = sort(setdiff(ct, ut)),
  uncompressed_only_sign_compatible_tfs = sort(setdiff(ut, ct)),
  key_tfs = tf[tf$source %in% c("Gene__Esr1", "Gene__Dbp", "Gene__Nfic"), ],
  network_overlap = network_overlap,
  final_graph_md5 = as.list(tools::md5sum(paths)),
  comparison_script_md5 = unname(tools::md5sum("scripts/moon/compare_factor1_runs.R")))
jsonlite::write_json(summary, file.path(out, "fullcomparison.json"),
                     pretty = TRUE, auto_unbox = TRUE, na = "null")
print(summary$tf_sign_compatible)
print(summary$key_tfs)
print(network_overlap)
cat("Full comparison diagnostics written; main run outputs unchanged.\n")
