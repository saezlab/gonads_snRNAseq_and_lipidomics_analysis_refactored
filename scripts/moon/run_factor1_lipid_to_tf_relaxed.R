#!/usr/bin/env Rscript
# Refit the fixed, strict score-layer LPS -> TF sensitivity from frozen inputs.
# Parameter exploration belongs in ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R.
options(stringsAsFactors = FALSE, digits = 17)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run this entry point with Rscript")
setwd(normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R ",
    "[--config=PATH] [--output=NEW_DIRECTORY] [--validate-only]\n",
    "Default: refit the fixed low-threshold LPS(18:0) candidate sensitivity, ",
    "validate it against the saved strict score-layer reference, and export provenance.\n",
    "Existing nonempty output directories are refused. Factor 1 lipid weights remain ",
    "post-score sign-consistency annotations, not ULM inputs.\n"
  ))
  quit(save = "no")
}
stopifnot(all(grepl("^--(config|output)=.+$|^--validate-only$", args)))
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}

source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
source("scripts/moon/lib/relaxed_lipid_tf.R")

config_path <- option("--config", "scripts/moon/config/factor1_lipid_to_tf_relaxed.json")
cfg <- relaxed_read_config(config_path)
exporter <- "scripts/moon/ancillary/export_factor1_lipid_to_tf_relaxed_network.py"
if (!file.exists(exporter) || !nzchar(Sys.which("python3"))) {
  stop("The relaxed network exporter and python3 are required")
}
data <- relaxed_load_data(cfg)
panel <- relaxed_panel_data(data, cfg$selected$tf_panel)
if ("--validate-only" %in% args) {
  cat("Validated selected relaxed LPS sensitivity: ", nrow(data$assignments),
      " mapping rows, ", length(data$candidate_nodes), " candidate chemicals, ",
      nrow(data$direct_pkn), " direct signaling edges, ", length(panel$downstream),
      " terminal TF observations. No files written.\n", sep = "")
  quit(save = "no")
}

out <- moon_new_output(option("--output", cfg$output))
relaxed_write_common_inputs(data, cfg, config_path, out)
fit <- relaxed_fit(panel, cfg$selected$pruning_steps, cfg$selected$n_layers, cfg$max_iterations)
status <- relaxed_candidate_status(data$assignments, fit)
node_status <- relaxed_candidate_node_summary(status)
seeds <- relaxed_sign_consensus_seeds(status)
solution <- relaxed_extract_score_layer_paths(
  fit$scores, fit$pkn_final, seeds, cfg$selected$primary, cfg$selected$secondary,
  data$rna
)
relaxed_validate_score_layer_solution(solution, fit, seeds, cfg$selected$primary, cfg$selected$secondary)
metadata <- relaxed_solution_metadata(solution, status)
eligibility <- relaxed_solution_is_eligible(solution, status, cfg)
if (!isTRUE(eligibility$eligible)) stop("Fixed selected relaxed route did not meet its predeclared strict extraction contract")
stopifnot(
  identical(metadata$source_nodes, as.character(cfg$selected$expected_source_node)),
  identical(metadata$measurement_features, as.character(cfg$selected$expected_measurement_feature)),
  identical(metadata$observed_tfs, sort(as.character(cfg$selected$expected_observed_tfs)))
)
comparison <- relaxed_compare_reference(fit, solution, cfg)

moon_save_fit(fit, out)
moon_write_csv(fit$scores, file.path(out, "moon_scores_unconstrained.csv"))
moon_write_csv(status, file.path(out, "candidate_assignment_status.csv"))
moon_write_csv(node_status, file.path(out, "candidate_node_status.csv"))
moon_write_csv(data.frame(source = names(panel$downstream), activity = unname(panel$downstream)),
               file.path(out, "downstream_tf_observations.csv"))
moon_write_csv(data.frame(source = names(seeds), copied_measurement_sign = unname(seeds)),
               file.path(out, "candidate_sign_consensus_seeds.csv"))

solution$ATT <- relaxed_annotate_solution(solution, status)
saveRDS(solution, file.path(out, "solution.rds"))
moon_write_csv(solution$SIF, file.path(out, "network_SIF.csv"))
moon_write_csv(solution$ATT, file.path(out, "network_ATT.csv"))
relaxed_write_graphml(solution, file.path(out, "network.graphml"))
relaxed_write_json(list(
  extraction = solution$audit,
  candidate_metadata = metadata,
  eligibility = eligibility,
  topology_radius = cfg$selected$pruning_steps,
  n_layers = cfg$selected$n_layers,
  compression = "safe_same_children",
  extraction_policy = "strict_consecutive_score_layers",
  selection = "fixed exploratory LPS sensitivity"
), file.path(out, "extraction_audit.json"))
relaxed_write_json(comparison, file.path(out, "reference_comparison.json"))

# The package helper explains the selected candidate's recursive ULM score using
# the same compressed scoring PKN as the fit. It is saved separately from the
# displayed strict solution because the two are related but distinct audits.
score_dir <- file.path(out, "score_explanation")
dir.create(score_dir, recursive = TRUE)
selected_source <- metadata$source_nodes[[1L]]
if (!selected_source %in% fit$scores_compressed$source) {
  stop("Selected candidate was unexpectedly absent from the compressed score space")
}
score_explanation <- NULL
invisible(capture.output(
  score_explanation <- cosmosR::get_moon_scoring_network(
    selected_source, fit$pkn_scoring, fit$scores_compressed
  )
))
relaxed_required_columns(score_explanation$SIF, c("source", "target", "interaction"),
                         "Selected candidate score-explanation SIF")
relaxed_required_columns(score_explanation$ATT, c("source", "score", "level"),
                         "Selected candidate score-explanation ATT")
saveRDS(score_explanation, file.path(score_dir, "score_explanation.rds"))
moon_write_csv(score_explanation$SIF, file.path(score_dir, "network_SIF.csv"))
moon_write_csv(score_explanation$ATT, file.path(score_dir, "network_ATT.csv"))

summary <- list(
  analysis = cfg$analysis,
  direction = "conditional relaxed lipid candidate -> TF footprint activity",
  tf_panel = cfg$selected$tf_panel,
  pruning_steps = cfg$selected$pruning_steps,
  n_layers = cfg$selected$n_layers,
  primary = cfg$selected$primary,
  secondary = cfg$selected$secondary,
  mapping_rows = nrow(data$assignments),
  unique_candidate_chemical_nodes = length(data$candidate_nodes),
  direct_signaling_edges = nrow(data$direct_pkn),
  candidate_measurement_features = metadata$measurement_features,
  candidate_source_nodes = metadata$source_nodes,
  independent_edge_fingerprint_groups = metadata$independent_edge_fingerprint_groups,
  observed_tfs = metadata$observed_tfs,
  nodes = nrow(solution$ATT),
  edges = nrow(solution$SIF),
  score_explanation_edges = nrow(score_explanation$SIF),
  reference_reproduced = comparison$passed,
  interpretation = "Fixed low-threshold exploratory sensitivity. The copied Factor 1 lipid loading is a post-score sign-consistency annotation; Ar and Esr1 are supplied terminal TF observations, not independent validation targets."
)
relaxed_write_json(summary, file.path(out, "results_summary.json"))

export_status <- system2("python3", c(
  shQuote(exporter), "--results", shQuote(out), "--pkn", shQuote(cfg$pkn),
  "--mapping-evidence", shQuote(cfg$direct_evidence)
), stdout = file.path(out, "export.log"), stderr = file.path(out, "export.log"))
if (export_status != 0L) stop("Relaxed evidence export failed; inspect ", file.path(out, "export.log"))

dependency_paths <- unique(c(
  config_path, cfg$inputs, cfg$pkn, cfg$mapping, cfg$direct_evidence,
  cfg$selected$reference_fit, file.path(cfg$selected$reference_solution, c("network_SIF.csv", "network_ATT.csv")),
  "scripts/moon/run_factor1_lipid_to_tf_relaxed.R", "scripts/moon/lib/relaxed_lipid_tf.R",
  "scripts/moon/lib/workflow.R", "scripts/moon/moon_safe_helpers.R", exporter
))
moon_write_csv(data.frame(path = dependency_paths, md5 = unname(tools::md5sum(dependency_paths))),
               file.path(out, "dependency_manifest.csv"))
relaxed_write_json(list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  command = paste("Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R", paste(args, collapse = " ")),
  config = cfg, orientation = data$input$orientation, R_version = R.version.string,
  completed = TRUE, reference_comparison = comparison,
  score_construction = "Candidate chemical IDs seed directed topology; ULM scoring uses terminal TF observations with upstream_input=NULL; copied lipid measurements are checked after scoring.",
  compression = "Same signed children with complete candidate/TF/regulon protection; a fixed compression map is retained through RNA coherence.",
  selection = "Fixed strong-TF, pruning-4, layer-4, primary/secondary-0.5/0.5 strict score-layer sensitivity."
), file.path(out, "run_manifest.json"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
outputs <- list.files(out, recursive = TRUE, full.names = TRUE)
moon_write_csv(data.frame(path = substring(outputs, nchar(out) + 2L), md5 = unname(tools::md5sum(outputs))),
               file.path(out, "output_manifest.csv"))
cat("PASS: selected relaxed LPS sensitivity reproduced: ", nrow(solution$ATT), " nodes, ",
    nrow(solution$SIF), " edges; reference scores and compression verified. Results: ", out, "\n", sep = "")
