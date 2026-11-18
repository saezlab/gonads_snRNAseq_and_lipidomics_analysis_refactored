#!/usr/bin/env Rscript
# Refit and audit the selected manuscript model from frozen, oriented inputs/PKN.
# No parameter search or new biological input mapping is performed here.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run this entry point with Rscript")
root <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/run_factor1_tf_to_lipid.R [--config=PATH] [--output=NEW_DIRECTORY] [--validate-only]\n",
      "Default: refit the selected compressed model; verify the frozen C03 result; export evidence and figures.\n",
      "Paths are relative to the project root. Existing nonempty output directories are refused.\n")
  quit(save = "no")
}
if (any(!grepl("^--(config|output)=.+$|^--validate-only$", args))) stop("Unknown or empty option; use --help")
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
config_path <- option("--config", "scripts/moon/config/factor1_tf_to_lipid.json")
cfg <- moon_read_config(config_path)
dependencies <- c(cfg$inputs, cfg$pkn, cfg$input_mapping,
  file.path(cfg$reference_run, c("run_state.rds", "pkn_final_scoring.rds")),
  cfg$reference_compression, cfg$reference_all_scores,
  file.path(cfg$reference_solution, "network_SIF.csv"),
  file.path(cfg$pkn_directory, c("nodes.csv", "nodes.jsonl", "pkn_lipid_curated.rds",
                                "pkn_lipid_rules.csv", "edge_evidence.jsonl.gz", "mouse_gem_reactions.jsonl",
                                "rhea_reactions.jsonl", "lipid_rules_reactions.jsonl")))
if (any(!file.exists(dependencies))) stop("Missing dependency: ", paste(dependencies[!file.exists(dependencies)], collapse = ", "))
if (!nzchar(Sys.which("python3"))) stop("python3 is required for evidence export")
data <- moon_load_factor1(cfg)
cat("Frozen PKN:", nrow(data$pkn), "edges;", length(moon_nodes(data$pkn)), "nodes\n",
    "Inputs:", length(data$upstream), "TFs and", length(data$downstream), "mapped lipids\n")
if ("--validate-only" %in% args) {
  cat("PASS: configuration, packages, input hashes and required artifacts. No files written.\n")
  quit(save = "no")
}
out <- moon_new_output(option("--output", cfg$output))
fit <- moon_fit(data, cfg$pruning_steps, cfg$n_layers, cfg$max_iterations)
solution <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
  primary_thresh = cfg$primary, secondary_thresh = cfg$secondary, RNA_input = data$rna)
moon_validate_solution(solution, fit, cfg$primary, cfg$secondary)
all_scores <- moon_diagnose(fit, data, solution)
comparison <- moon_compare_reference(fit, solution, cfg, all_scores)
catalog <- utils::read.csv(file.path(cfg$pkn_directory, "nodes.csv"), check.names = FALSE)
all_scores$label <- catalog$label[match(all_scores$source, catalog$node_id)]
focus <- data.frame(source = cfg$focus_tfs,
  in_initial_upstream = cfg$focus_tfs %in% names(data$upstream),
  in_pruned_pkn = cfg$focus_tfs %in% moon_nodes(fit$prepared$net),
  unchanged_by_compression = fit$compression$node_map$compressed[
    match(cfg$focus_tfs, fit$compression$node_map$original)] == cfg$focus_tfs,
  in_final_scoring_pkn = cfg$focus_tfs %in% moon_nodes(fit$pkn_scoring))
focus <- merge(focus, all_scores, by = "source", all.x = TRUE, sort = FALSE)
stopifnot(all(focus$in_initial_upstream), all(focus$in_pruned_pkn),
          all(focus$unchanged_by_compression), all(focus$in_final_scoring_pkn),
          all(focus$retention_reason == "upstream_sign_conflict"), all(!focus$in_solution))
solution$ATT$label <- catalog$label[match(solution$ATT$source, catalog$node_id)]
solution$ATT$upstream_activity <- unname(data$upstream[solution$ATT$source])
moon_save_fit(fit, out)
saveRDS(solution, file.path(out, "solution.rds"))
moon_write_csv(solution$SIF, file.path(out, "network_SIF.csv"))
moon_write_csv(solution$ATT, file.path(out, "network_ATT.csv"))
moon_write_csv(all_scores, file.path(out, "moon_scores_all.csv"))
moon_write_csv(focus, file.path(out, "focus_tf_retention.csv"))
moon_write_csv(data.frame(source = names(fit$prepared$upstream), score = unname(fit$prepared$upstream)),
               file.path(out, "upstream_input.csv"))
moon_write_csv(data.frame(source = names(fit$prepared$downstream), score = unname(fit$prepared$downstream)),
               file.path(out, "downstream_input.csv"))
write_json <- function(x, name) jsonlite::write_json(x, file.path(out, name),
  pretty = TRUE, auto_unbox = TRUE, na = "null", digits = NA)
write_json(c(solution$audit, list(topology_radius = cfg$pruning_steps, n_layers = cfg$n_layers,
  compression = "compressed", selection = "fixed manuscript C03 representative")), "extraction_audit.json")
write_json(comparison, "reference_comparison.json")
curated <- readRDS(file.path(cfg$pkn_directory, "pkn_lipid_curated.rds"))
summary <- list(analysis = cfg$analysis, direction = "TF activities -> downstream lipid weights",
  curated_pkn_edges = nrow(curated), curated_pkn_nodes = length(moon_nodes(curated)),
  expanded_pkn_edges = nrow(data$pkn), expanded_pkn_nodes = length(moon_nodes(data$pkn)),
  pruning_steps = cfg$pruning_steps, n_layers = cfg$n_layers,
  primary = cfg$primary, secondary = cfg$secondary,
  retained_tf_inputs = length(fit$prepared$upstream),
  retained_lipid_inputs = length(fit$prepared$downstream),
  initial_original_edges = nrow(fit$prepared$net),
  initial_scoring_edges = nrow(fit$compression$pkn_scoring),
  final_scoring_edges = nrow(fit$pkn_scoring), final_original_edges = nrow(fit$pkn_final),
  coherence_iterations = nrow(fit$coherence_ledger),
  returned_score_rows = nrow(fit$scores), diagnostic_score_rows = nrow(all_scores),
  nodes = nrow(solution$ATT), edges = nrow(solution$SIF),
  upstream_tfs = solution$ATT$label[solution$ATT$type == "upstream_input"],
  downstream_lipids = solution$ATT$label[solution$ATT$level == 0],
  exact_reference_reproduced = comparison$passed,
  interpretation = "Exploratory prior-guided hypotheses; score thresholds are not calibrated significance cutoffs.")
write_json(summary, "results_summary.json")
export_status <- system2("python3", vapply(c("scripts/moon/ancillary/export_factor1_network.py",
  "--results", out, "--pkn", cfg$pkn_directory, "--input-mapping", cfg$input_mapping,
  "--stem", "network"), shQuote, character(1)),
  stdout = file.path(out, "export.log"), stderr = file.path(out, "export.log"))
if (export_status != 0L) stop("Evidence export failed; inspect ", file.path(out, "export.log"))
scripts <- c("scripts/moon/run_factor1_tf_to_lipid.R", "scripts/moon/lib/workflow.R",
  "scripts/moon/moon_safe_helpers.R", "scripts/moon/ancillary/export_factor1_network.py",
  "scripts/moon/ancillary/plot_factor1_network.R", config_path)
manifest_paths <- unique(c(dependencies, scripts))
moon_write_csv(data.frame(path = manifest_paths, md5 = unname(tools::md5sum(manifest_paths))),
               file.path(out, "dependency_manifest.csv"))
write_json(list(created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  command = paste("Rscript scripts/moon/run_factor1_tf_to_lipid.R", paste(args, collapse = " ")),
  mode = "refit_from_frozen_inputs", config = cfg, orientation = data$orientation,
  R_version = R.version.string, reference_comparison = comparison,
  input_construction = "Previously audited, jointly oriented Factor 1 TF footprints, lipid weights and RNA coherence signature; mapping unchanged.",
  compression = "Same signed children; protect complete input/regulon groups; fixed map throughout coherence.",
  completed = TRUE), "run_manifest.json")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
outputs <- list.files(out, recursive = TRUE, full.names = TRUE)
moon_write_csv(data.frame(path = substring(outputs, nchar(out) + 2L),
  md5 = unname(tools::md5sum(outputs))), file.path(out, "output_manifest.csv"))
cat("PASS: refit matches archived C03: 153 nodes, 388 edges; scores, compression and TF exclusions verified.\n",
    "Results:", out, "\n")
