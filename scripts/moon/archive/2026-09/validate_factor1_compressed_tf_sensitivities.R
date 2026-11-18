#!/usr/bin/env Rscript
# Independent saved-state validation; does not refit or change the sensitivities.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
base <- "results/cosmos_factor1_compressed_review"
out <- file.path(base, "tf_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input <- readRDS("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
nodes <- function(e) unique(c(e$source, e$target))
key <- function(e) paste(e$source, e$interaction, e$target, sep = "\t")
reg <- input$tf_regulon[, c("source", "target", "mor")]
reg_signed <- reg; names(reg_signed)[3] <- "interaction"
focus <- paste0("Gene__", c("Esr1", "Nfic", "Dbp"))
dirs <- file.path(base, c("topology_sensitivity", "prescore_sensitivity", "inputs_only_sensitivity"))
paths <- unlist(lapply(dirs, function(x) list.files(x, pattern = "^run_state\\.rds$", recursive = TRUE, full.names = TRUE)), use.names = FALSE)
stopifnot(length(paths) == 8L)
summary <- list(); tf_summary <- list()
for (path in paths) {
  state <- readRDS(path); run <- dirname(path)
  comp <- readRDS(file.path(run, "compression.rds"))
  score_net <- readRDS(file.path(run, "pkn_final_scoring.rds"))
  final <- readRDS(file.path(run, "pkn_final_uncompressed.rds"))
  raw_csv <- utils::read.csv(file.path(run, "moon_scores_before_upstream_sign_filter.csv"))
  ledger <- utils::read.csv(file.path(run, "coherence_ledger.csv"))
  mv <- stats::setNames(comp$node_map$compressed, comp$node_map$original)
  projected <- final
  projected$source <- unname(mv[projected$source]); projected$target <- unname(mv[projected$target])
  stopifnot(state$converged, tail(ledger$edges_removed, 1) == 0,
            identical(moon_check_pkn(projected), moon_check_pkn(score_net)),
            identical(moon_check_pkn(final), moon_check_pkn(state$pkn_final)),
            all(mv[names(state$tf)] == names(state$tf)),
            all(mv[names(state$lipid)] == names(state$lipid)),
            comp$n_merged_nodes > 0L,
            !any(final$source %in% names(state$lipid)),
            !anyDuplicated(state$scores$source))
  raw_idx <- match(state$scores$source, raw_csv$source)
  stopifnot(!anyNA(raw_idx), max(abs(state$scores$score - raw_csv$score[raw_idx])) < 1e-10)
  # An independent check that the claimed converged graph is a fixed point.
  coherent <- cosmosR::filter_incohrent_TF_target(state$scores, reg, final, input$hepatocyte_gene_weight)
  stopifnot(identical(key(final), key(coherent)))
  if (grepl("inputs_only", state$run_id)) {
    checkpoint <- readRDS(file.path(run, paste0("iteration_", state$iterations, ".rds")))
    from_checkpoint <- moon_safe_decompress_scores(checkpoint$scores_compressed, checkpoint$compression$node_map)
    stopifnot(identical(from_checkpoint, state$scores),
              identical(moon_check_pkn(checkpoint$pkn_original), moon_check_pkn(comp$pkn_uncompressed)),
              identical(moon_check_pkn(checkpoint$compression$pkn_scoring), moon_check_pkn(score_net)),
              setequal(comp$protected_nodes, c(names(state$tf), names(state$lipid))))
    checkpoints <- list.files(run, pattern = "^iteration_[0-9]+\\.rds$", full.names = TRUE)
    stopifnot(length(checkpoints) == state$iterations)
    for (cp in checkpoints) {
      ck <- readRDS(cp)
      proj <- ck$pkn_original
      map <- stats::setNames(ck$compression$node_map$compressed, ck$compression$node_map$original)
      proj$source <- unname(map[proj$source]); proj$target <- unname(map[proj$target])
      stopifnot(identical(moon_check_pkn(proj), moon_check_pkn(ck$compression$pkn_scoring)))
    }
  }
  for (tf in focus) {
    rr <- raw_csv[match(tf, raw_csv$source), ]
    stopifnot(tf %in% names(state$tf), tf %in% nodes(final), rr$score > 0,
              input$tf_activity[tf] < 0, !tf %in% state$scores$source)
    prev <- raw_csv[raw_csv$level == rr$level - 1, ]
    e <- score_net[score_net$source == tf, ]
    # Raw CSV is decompressed, so restore the unique virtual previous layer.
    prev <- prev[!duplicated(prev$compressed_source), ]
    coef <- stats::setNames(rep(0, nrow(prev)), prev$compressed_source)
    keep <- e$target %in% names(coef); coef[e$target[keep]] <- e$interaction[keep]
    cr <- stats::cor(coef, prev$score)
    score_reconstructed <- cr * sqrt((nrow(prev) - 2) / ((1-cr+1e-20)*(1+cr+1e-20)))
    stopifnot(abs(score_reconstructed - rr$score) < 1e-10)
    tf_summary[[length(tf_summary) + 1L]] <- data.frame(run_id = state$run_id,
      source = tf, input_activity = unname(input$tf_activity[tf]), raw_score = rr$score,
      level = rr$level, scored_children = sum(coef != 0), background_n = length(coef),
      retained_before_scoring = TRUE, retained_after_sign_filter = FALSE,
      exclusion_reason = "post_scoring_upstream_sign_conflict")
  }
  summary[[length(summary) + 1L]] <- data.frame(run_id = state$run_id,
    n_tfs = length(state$tf), n_lipids = length(state$lipid),
    scoring_edges = nrow(score_net), original_edges = nrow(final),
    merged_nodes = comp$n_merged_nodes, merged_groups = comp$n_groups,
    iterations = state$iterations, convergence_validated = TRUE,
    projection_validated = TRUE, score_copy_validated = TRUE,
    all_focus_tf_positive_raw_scores = TRUE)
}

# Reconstruct the initial-only prescore filter independently from exact evidence.
full_pkn <- moon_check_pkn(readRDS("results/cosmos_mouse_pkn/pkn_lipid_rules.rds"))
tf_input <- input$tf_activity[abs(input$tf_activity) > 2 & names(input$tf_activity) %in% nodes(full_pkn)]
lipid_input <- input$lipid_weight[input$lipid_weight != 0 & names(input$lipid_weight) %in% nodes(full_pkn)]
terminal <- full_pkn[!full_pkn$source %in% names(lipid_input), ]
evidence <- utils::read.csv(file.path(out, "exact_transcription_evidence.csv"))
stopifnot(all(evidence$mechanism == "transcription"), all(evidence$resource == "OmniPath:grn"))
allowed <- intersect(key(evidence), key(reg_signed))
ta <- unname(tf_input[terminal$source]); wa <- unname(input$hepatocyte_gene_weight[terminal$target])
bad <- key(terminal) %in% allowed & is.finite(ta) & is.finite(wa) & sign(ta * wa * terminal$interaction) < 0
expected <- terminal[bad, ]
actual <- utils::read.csv(file.path(base, "prescore_sensitivity/prescore_removed_edges.csv"))
stopifnot(setequal(key(expected), key(actual)), nrow(expected) == nrow(actual), length(tf_input) == 49L)
prescore_summary <- data.frame(all_initial_measured_tf_inputs = length(tf_input),
  removed_edges = nrow(actual), affected_tfs = length(unique(actual$source)),
  all_removed_triples_in_grn_and_input_regulon = all(key(actual) %in% allowed),
  exact_expected_removal_set = TRUE)
utils::write.csv(do.call(rbind, summary), file.path(out, "additional_sensitivity_validation.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, tf_summary), file.path(out, "additional_sensitivity_tf_audit.csv"), row.names = FALSE)
utils::write.csv(prescore_summary, file.path(out, "prescore_filter_validation.csv"), row.names = FALSE)
print(do.call(rbind, summary), row.names = FALSE)
print(do.call(rbind, tf_summary), row.names = FALSE)
print(prescore_summary, row.names = FALSE)
cat("All eight additional sensitivity states validated without refitting.\n")
