#!/usr/bin/env Rscript
# Factor 1 lipid -> TF analysis with explicit topology/measurement feasibility gates.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/run_factor1_lipid_to_tf.R [--output=NEW_DIRECTORY] [--config=PATH] [--include-ph73] [--validate-only]\n",
      "Default: exact-mapped direct signalling and lipid-only precursor policies.\n",
      "--include-ph73 additionally scores the existing acid/base-normalized prior as an exploratory sensitivity.\n",
      "Empty topology or insufficient ULM observations is reported, never replaced by fabricated scores.\n")
  quit(save = "no")
}
stopifnot(all(grepl("^--(output|config)=.+$|^--(include-ph73|validate-only)$", args)))
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
config_path <- option("--config", "scripts/moon/config/factor1_lipid_to_tf.json")
cfg <- jsonlite::read_json(config_path, simplifyVector = TRUE)
stopifnot(identical(cfg$analysis, "factor1_lipid_to_tf"))
# Reuse the validated frozen input/software contract, not the forward model's direction.
contract <- moon_read_config(cfg$dependency_contract)
stopifnot(identical(cfg$inputs, contract$inputs))
input <- readRDS(cfg$inputs)
raw <- moon_check_pkn(readRDS(contract$pkn))
moon_vector_check(input$lipid_weight)
moon_vector_check(input$tf_activity)
tf <- input$tf_activity[input$tf_activity != 0 & names(input$tf_activity) %in% moon_nodes(raw)]
panels <- list(all_tfs = tf, strong_tfs = tf[abs(tf) > cfg$strong_tf_cutoff])
if ("--validate-only" %in% args) {
  cat("Validated", length(input$lipid_weight), "exact-mapped lipids;", length(tf),
      "TF observations; optional normalized sensitivity:", "--include-ph73" %in% args, "\n")
  quit(save = "no")
}
out <- moon_new_output(option("--output", cfg$output))
write_json <- function(x, name) jsonlite::write_json(x, file.path(out, name),
  pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
write_json(list(user_direction = "lipids -> TFs", factor = "Factor1",
  orientation = input$orientation, upstream = cfg$upstream_semantics,
  downstream = cfg$downstream_semantics, rna = cfg$rna_semantics,
  ph73_sensitivity_enabled = "--include-ph73" %in% args,
  inherited_mapping = "Frozen exact species mapping; no hierarchy proxy or new chemical identity assignment.",
  factor_view_variance_percent = list(hepatocyte_RNA = 8.906, lipid = 29.958),
  stability = "Previously audited descriptive consistency across fitted factor-number choices; no resampling validation.",
  timing = "Cross-sectional bulk-liver lipid and hepatocyte RNA signatures; no temporal ordering inferred.",
  limitations = cfg$limitations), "input_decision.json")
coverage_dir <- file.path(out, "coverage")
status <- system2("python3", vapply(c("scripts/moon/ancillary/audit_lipid_tf_coverage.py",
  "--output", coverage_dir, "--pkn", cfg$pkn_directory, "--inputs", cfg$input_directory), shQuote, character(1)),
  stdout = file.path(out, "coverage.log"), stderr = file.path(out, "coverage.log"))
if (status != 0L) stop("Coverage audit failed; inspect ", file.path(out, "coverage.log"))
protein <- function(x) grepl("^(Gene|UniProt)__", x)
direct <- raw[(protein(raw$source) & protein(raw$target)) |
                (grepl("^Metab__", raw$source) & protein(raw$target)), , drop = FALSE]
graphs <- list(direct_regulation = direct,
  lipid_precursor = utils::read.csv(file.path(coverage_dir, "pkn_lipid_precursor.csv")))
if ("--include-ph73" %in% args) {
  graphs$ph73_lipid_precursor <- utils::read.csv(file.path(coverage_dir, "pkn_ph73_lipid_precursor.csv"))
}
grid <- expand.grid(primary = cfg$primary_thresholds, secondary = cfg$secondary_thresholds)
grid <- grid[grid$primary >= grid$secondary, , drop = FALSE]
run_rows <- grid_rows <- list()
catalog <- utils::read.csv(file.path(cfg$pkn_directory, "nodes.csv"), check.names = FALSE)
for (policy in names(graphs)) {
  pkn <- moon_check_pkn(graphs[[policy]])
  for (panel in names(panels)) for (steps in cfg$pruning_steps) {
    upstream <- input$lipid_weight[names(input$lipid_weight) %in% moon_nodes(pkn)]
    downstream <- panels[[panel]][names(panels[[panel]]) %in% moon_nodes(pkn)]
    prepared <- tryCatch(moon_prepare(pkn, upstream, downstream, steps), error = identity)
    reason <- if (!length(upstream)) "no_mapped_upstream_nodes" else
      if (inherits(prepared, "error")) "no_directed_input_to_observation_topology" else
        if (length(prepared$downstream) < 3L) "insufficient_observations_for_ULM" else "ready"
    for (layers in cfg$n_layers) {
      id <- paste(policy, panel, paste0("pruning", steps), paste0("layers", layers), sep = "_")
      row <- data.frame(run_id = id, policy = policy, tf_panel = panel,
        pruning_steps = steps, n_layers = layers, status = reason,
        upstream_lipids = if (inherits(prepared, "error")) 0L else length(prepared$upstream),
        downstream_tfs = if (inherits(prepared, "error")) 0L else length(prepared$downstream),
        converged = FALSE, iterations = NA_integer_, scoring_edges = NA_integer_, error = "")
      fit <- NULL
      if (reason == "ready") {
        message("Fitting ", id)
        data <- list(pkn = pkn, upstream = upstream, downstream = downstream,
          regulon = input$tf_regulon[, c("source", "target", "mor")], rna = input$hepatocyte_gene_weight)
        fit <- tryCatch(moon_fit(data, steps, layers, cfg$max_iterations), error = identity)
        if (inherits(fit, "error")) {
          row$status <- "scoring_failed"
          row$error <- conditionMessage(fit)
        } else {
          row$status <- "scored"
          row$converged <- fit$converged
          row$iterations <- nrow(fit$coherence_ledger)
          row$scoring_edges <- nrow(fit$pkn_scoring)
          run_dir <- file.path(out, "runs", id)
          dir.create(run_dir, recursive = TRUE)
          moon_save_fit(fit, run_dir)
          for (j in seq_len(nrow(grid))) {
            s <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
              primary_thresh = grid$primary[j], secondary_thresh = grid$secondary[j],
              RNA_input = data$rna)
            moon_validate_solution(s, fit, grid$primary[j], grid$secondary[j])
            eligible <- nrow(s$SIF) > 0L && nrow(s$SIF) < cfg$maximum_solution_edges_exclusive &&
              s$audit$retained_upstream >= cfg$minimum_upstream_lipids &&
              s$audit$retained_observed >= cfg$minimum_downstream_tfs
            grid_rows[[length(grid_rows) + 1L]] <- data.frame(run_id = id, policy = policy,
              tf_panel = panel, pruning_steps = steps, n_layers = layers,
              primary = grid$primary[j], secondary = grid$secondary[j],
              edges = nrow(s$SIF), nodes = nrow(s$ATT), upstream_lipids = s$audit$retained_upstream,
              downstream_tfs = s$audit$retained_observed, eligible = eligible)
          }
        }
      }
      run_rows[[length(run_rows) + 1L]] <- row
      moon_write_csv(do.call(rbind, run_rows), file.path(out, "run_summary.csv"))
      if (length(grid_rows)) moon_write_csv(do.call(rbind, grid_rows), file.path(out, "threshold_grid.csv"))
    }
  }
}
all_runs <- do.call(rbind, run_rows)
thresholds <- if (length(grid_rows)) do.call(rbind, grid_rows) else data.frame(
  run_id = character(), policy = character(), tf_panel = character(), pruning_steps = integer(),
  n_layers = integer(), primary = numeric(), secondary = numeric(), edges = integer(), nodes = integer(),
  upstream_lipids = integer(), downstream_tfs = integer(), eligible = logical())
moon_write_csv(thresholds, file.path(out, "threshold_grid.csv"))
eligible <- thresholds[thresholds$eligible, , drop = FALSE]
moon_write_csv(eligible, file.path(out, "eligible_candidates.csv"))
# Candidate comparison is explicit; no automated cross-background biological ranking.
write_json(list(analysis = cfg$analysis, ph73_sensitivity_enabled = "--include-ph73" %in% args,
  policies = names(graphs), runs = nrow(all_runs), runs_scored = sum(all_runs$status == "scored"),
  status_counts = as.list(table(all_runs$status)), threshold_evaluations = nrow(thresholds),
  eligible_candidates = nrow(eligible),
  strict_direct_signalling_supported = any(eligible$policy == "direct_regulation"),
  primary_lipid_precursor_supported = any(eligible$policy == "lipid_precursor"),
  interpretation = "No direct measured-lipid signaling claim without an eligible regulatory path. Normalized-policy candidates, if present, remain separate metabolic-availability hypotheses requiring co-substrate inspection."), "analysis_summary.json")
sources <- c(config_path, cfg$dependency_contract, cfg$inputs,
  "scripts/moon/run_factor1_lipid_to_tf.R", "scripts/moon/lib/workflow.R", "scripts/moon/moon_safe_helpers.R",
  "scripts/moon/ancillary/audit_lipid_tf_coverage.py")
moon_write_csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))), file.path(out, "dependency_manifest.csv"))
write_json(list(created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  config = cfg, command_args = args, R_version = R.version.string,
  packages = contract$packages, cosmos_commit = contract$cosmos_commit,
  completed = TRUE), "run_manifest.json")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed lipid-to-TF feasibility and scoring workflow:", nrow(eligible), "eligible candidates.\n",
    "Results:", out, "\n")
