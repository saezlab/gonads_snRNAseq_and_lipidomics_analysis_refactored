#!/usr/bin/env Rscript
# Extract a bounded threshold grid from a completed expanded-PKN fit; never rescore.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(script), "../../.."), mustWork = TRUE)
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R --fit=COMPLETED_FIT_DIR --output=NEW_DIR [--primary=CSV] [--secondary=CSV] [--layer-consistent] [--validate-only]\n",
      "Reuses saved scores/graphs without any MOON scoring. Retains every valid grid row.\n",
      "Default grid combines original cutoffs with lower cutoffs: primary .25 through10, secondary0 through3 (explicit finite grids).\n",
      "--layer-consistent restricts extraction to source_level = target_level + 1; scoring remains unchanged.\n",
      "Eligible: 0 < edges < 500, at least two TFs and two observed lipids.\n",
      "Representative: descending secondary, primary, observed-lipid count, TF count, then fewer edges.\n")
  quit(save = "no")
}
if (any(!grepl("^--(fit|output|primary|secondary)=.+$|^--(validate-only|layer-consistent)$", args))) stop("Unknown or empty option; use --help")
option <- function(key, default = NULL) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (!length(hit)) return(default)
  substring(hit, nchar(key) + 2L)
}
numbers <- function(value, positive = FALSE) {
  if (grepl("(^,|,$|,,)", value)) stop("Empty numeric grid entry")
  x <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  stopifnot(length(x) > 0L, all(is.finite(x)), all(x >= 0), !positive || all(x > 0))
  sort(unique(x))
}
fit_path <- option("--fit")
out_path <- option("--output")
if (is.null(fit_path) || is.null(out_path)) stop("--fit and --output are required")
fit_path <- normalizePath(fit_path, mustWork = TRUE)
primary <- numbers(option("--primary", ".25,.5,.75,1,1.25,1.5,1.85,2.15,2.5,3,4,5,6,8,10"), TRUE)
secondary <- numbers(option("--secondary", "0,.05,.1,.25,.4,.55,.75,1,1.075,1.25,1.5,2,2.5,3"))
layer_consistent <- "--layer-consistent" %in% args
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
required <- file.path(fit_path, c("run_state.rds", "analysis_inputs.rds", "completed_manifest.json", "moon_scores_all.csv"))
if (any(!file.exists(required))) stop("Missing completed-fit dependencies: ", paste(required[!file.exists(required)], collapse = ", "))
manifest <- jsonlite::read_json(file.path(fit_path, "completed_manifest.json"), simplifyVector = TRUE)
stopifnot(isTRUE(manifest$completed), identical(utils::packageDescription("cosmosR")$RemoteSha, manifest$cosmos_commit))
fit <- readRDS(file.path(fit_path, "run_state.rds"))
data <- readRDS(file.path(fit_path, "analysis_inputs.rds"))
all_scores <- utils::read.csv(file.path(fit_path, "moon_scores_all.csv"), check.names = FALSE)
stopifnot(isTRUE(fit$converged), utils::tail(fit$coherence_ledger$edges_removed, 1L) == 0L,
  !anyDuplicated(fit$scores$source), all(is.finite(fit$scores$score)),
  setequal(fit$scores$source[fit$scores$level == 0], names(fit$prepared$downstream)),
  !anyDuplicated(all_scores$source))
extraction_net <- fit$pkn_final
if (layer_consistent) {
  level <- stats::setNames(fit$scores$level, fit$scores$source)
  a <- unname(level[extraction_net$source]); b <- unname(level[extraction_net$target])
  extraction_net <- extraction_net[which(!is.na(a) & !is.na(b) & a == b + 1), , drop = FALSE]
}
grid <- expand.grid(primary = primary, secondary = secondary, KEEP.OUT.ATTRS = FALSE)
grid$valid_pair <- grid$primary >= grid$secondary
if (!any(grid$valid_pair)) stop("No valid primary >= secondary grid pair")
if ("--validate-only" %in% args) {
  cat("PASS: completed fit and source commit;", sum(grid$valid_pair), "valid threshold pairs. No scoring or writes.\n")
  quit(save = "no")
}
out <- moon_new_output(out_path)
file.copy(script, file.path(out, "reextraction_script_snapshot.R"), overwrite = FALSE)
moon_write_csv(data.frame(source = names(fit$prepared$upstream), score = unname(fit$prepared$upstream)),
  file.path(out, "upstream_input.csv"))
moon_write_csv(data.frame(source = names(fit$prepared$downstream), score = unname(fit$prepared$downstream)),
  file.path(out, "downstream_input.csv"))
moon_write_csv(data.frame(source = names(data$rna), score = unname(data$rna)),
  file.path(out, "rna_coherence_input.csv"))
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, na = "null", null = "null", digits = NA)
policy <- "Within-grid representative: 0<edges<500, >=2TFs and >=2observations; descending secondary then primary then observation/TF coverage then fewer edges. Not a calibrated optimum."
write_json(list(fit_reference = fit_path, no_rescoring = TRUE, primary = primary, secondary = secondary,
  selection_policy = policy, source_fit_mode = manifest$mode, command_args = args,
  layer_consistent_extraction = layer_consistent, extraction_edges = nrow(extraction_net), scoring_graph_unchanged = TRUE,
  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, "search_manifest.json"))
rows <- solutions <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  p <- grid$primary[i]; s <- grid$secondary[i]
  row <- data.frame(primary = p, secondary = s, valid_pair = grid$valid_pair[i],
    status = "invalid_threshold_pair", edges = NA_integer_, nodes = NA_integer_,
    upstream_tfs = NA_integer_, observed_lipids = NA_integer_, validated = FALSE, eligible = FALSE)
  if (grid$valid_pair[i]) {
    ans <- moon_extract_paths(fit$scores, extraction_net, fit$prepared$upstream, p, s, data$rna)
    moon_validate_solution(ans, fit, p, s)
    row$status <- if (nrow(ans$SIF)) "ok" else "empty"
    row$edges <- nrow(ans$SIF); row$nodes <- nrow(ans$ATT)
    row$upstream_tfs <- ans$audit$retained_upstream
    row$observed_lipids <- ans$audit$retained_observed
    row$validated <- TRUE
    row$eligible <- row$edges > 0L && row$edges < 500L && row$upstream_tfs >= 2L && row$observed_lipids >= 2L
    solutions[[i]] <- ans
  }
  rows[[i]] <- row
  if (i %% 20L == 0L) {
    moon_write_csv(do.call(rbind, rows[seq_len(i)]), file.path(out, "threshold_grid_results.csv"))
    cat("Recorded", i, "/", nrow(grid), "threshold rows; no rescoring.\n")
  }
}
results <- do.call(rbind, rows)
moon_write_csv(results, file.path(out, "threshold_grid_results.csv"))
eligible <- which(results$eligible)
if (length(eligible)) {
  ranked <- eligible[order(-results$secondary[eligible], -results$primary[eligible],
    -results$observed_lipids[eligible], -results$upstream_tfs[eligible], results$edges[eligible])]
  chosen <- ranked[1L]
  solution <- solutions[[chosen]]
  selected <- as.list(results[chosen, ])
  audit <- solution$audit
  saveRDS(solution, file.path(out, "solution.rds"))
} else {
  solution <- list(SIF = data.frame(source = character(), target = character(), interaction = numeric(), consistency = logical()),
                   ATT = data.frame(source = character(), score = numeric(), level = numeric(), type = character(), RNA_input = numeric()))
  selected <- NULL
  audit <- list(selected_solution_available = FALSE, reason = "No threshold pair met compact coverage criteria")
}
moon_write_csv(solution$SIF, file.path(out, "network_SIF.csv"))
moon_write_csv(solution$ATT, file.path(out, "network_ATT.csv"))
write_json(c(audit, list(fit_reference = fit_path, no_rescoring = TRUE, selection_policy = policy, layer_consistent_extraction = layer_consistent)),
  file.path(out, "extraction_audit.json"))
# Diagnostics use the saved pre-upstream-sign-filter score table, never moon().
ids <- sort(unique(c(moon_nodes(fit$prepared$net), all_scores$source, names(data$upstream), names(data$downstream))))
diagnostic <- data.frame(source = ids, input_TF = ids %in% names(data$upstream),
  input_lipid = ids %in% names(data$downstream), retained_after_preparation = ids %in% moon_nodes(fit$prepared$net),
  returned_by_fit = ids %in% fit$scores$source, in_solution = ids %in% solution$ATT$source,
  score_before_upstream_filter = all_scores$score[match(ids, all_scores$source)],
  level = all_scores$level[match(ids, all_scores$source)],
  upstream_activity = unname(data$upstream[ids]))
diagnostic$upstream_sign_agrees <- ifelse(diagnostic$input_TF & !is.na(diagnostic$score_before_upstream_filter),
  sign(diagnostic$score_before_upstream_filter) == sign(diagnostic$upstream_activity), NA)
diagnostic$retention_reason <- ifelse(!diagnostic$retained_after_preparation, "topology_loss",
  ifelse(is.na(diagnostic$score_before_upstream_filter), "not_scored_within_layers",
    ifelse(!diagnostic$returned_by_fit & diagnostic$input_TF & !diagnostic$upstream_sign_agrees,
      "upstream_sign_conflict", ifelse(diagnostic$in_solution, "selected",
        ifelse(diagnostic$returned_by_fit, "below_threshold_or_no_complete_consistent_path", "not_returned")))))
moon_write_csv(diagnostic, file.path(out, "node_retention_diagnostics.csv"))
summary <- list(fit_reference = fit_path, no_rescoring = TRUE, source_fit_mode = manifest$mode,
  layer_consistent_extraction = layer_consistent, extraction_edges = nrow(extraction_net), scoring_graph_unchanged = TRUE,
  pruning_steps = fit$pruning_steps, n_layers = fit$n_layers, compression = "same_saved_safe_compressed_fit",
  grid_rows = nrow(results), valid_pairs = sum(results$valid_pair), validated_pairs = sum(results$validated),
  eligible_combinations = length(eligible), selected_solution_available = length(eligible) > 0L,
  selected = selected, selection_policy = policy,
  selected_TFs = solution$ATT$source[solution$ATT$type == "upstream_input"],
  selected_lipids = solution$ATT$source[solution$ATT$level == 0],
  interpretation = "Threshold sensitivity of saved fit; lower cutoffs indicate weaker numerical evidence, not calibrated significance.",
  finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
write_json(summary, file.path(out, "results_summary.json"))
dependencies <- c(required, "scripts/moon/moon_safe_helpers.R", "scripts/moon/lib/workflow.R", script, file.path(out, "reextraction_script_snapshot.R"))
moon_write_csv(data.frame(path = dependencies, md5 = unname(tools::md5sum(dependencies))), file.path(out, "dependency_manifest.csv"))
writeLines(c("# Expanded Factor 1 threshold re-extraction", "", "This ancillary analysis reuses the completed compressed fit named in results_summary.json. No MOON scoring is performed.",
  "All valid, invalid and empty threshold rows are retained. If no pair meets the compact coverage criteria, network files are empty and no solution.rds is written.",
  "node_retention_diagnostics.csv uses the saved pre-upstream-filter scores. Thresholds are exploratory numerical settings, not calibrated significance measures.", policy), file.path(out, "README.md"))
cat("PASS: re-extraction;", length(eligible), "eligible combinations; selected", nrow(solution$SIF),
    "edges,", nrow(solution$ATT), "nodes. Output:", out, "\n")
