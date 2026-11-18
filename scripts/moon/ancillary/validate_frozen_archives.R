#!/usr/bin/env Rscript
# Replay saved MOON extraction/reference contracts without fitting or scoring.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../../.."), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/validate_frozen_archives.R [--output=NEW_JSON]\n",
      "Checks frozen dependency hashes, saved scores/compression and selected extraction; no model fit or new scores.\n",
      "Software version differences are reported, not waived for a full fit. Full-fit --validate-only gates remain unchanged.\n")
  quit(save = "no")
}
if (any(!grepl("^--output=.+$", args)) || length(args) > 1L) stop("Unknown or repeated option")
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
source("scripts/moon/lib/relaxed_lipid_tf.R")
check_md5 <- function(cfg) {
  paths <- names(cfg$input_md5)
  stopifnot(all(file.exists(paths)))
  observed <- unname(tools::md5sum(paths))
  stopifnot(identical(observed, as.character(unlist(cfg$input_md5, use.names = FALSE))))
  length(paths)
}
cfg <- jsonlite::read_json("scripts/moon/config/factor1_tf_to_lipid.json", simplifyVector = TRUE)
forward_hashes <- check_md5(cfg)
data <- moon_load_factor1(cfg)
fit <- readRDS(file.path(cfg$output, "run_state.rds"))
saved <- readRDS(file.path(cfg$output, "solution.rds"))
all_scores <- utils::read.csv(file.path(cfg$output, "moon_scores_all.csv"), check.names = FALSE)
forward <- moon_compare_reference(fit, saved, cfg, all_scores)
replayed <- moon_extract_paths(fit$scores, fit$pkn_final, fit$prepared$upstream,
  primary_thresh = cfg$primary, secondary_thresh = cfg$secondary, RNA_input = data$rna)
moon_validate_solution(replayed, fit, cfg$primary, cfg$secondary)
stopifnot(setequal(moon_edge_keys(replayed$SIF), moon_edge_keys(saved$SIF)))
forward$saved_extraction_identical <- TRUE
forward$edges <- nrow(saved$SIF)
forward$nodes <- nrow(saved$ATT)
forward$dependency_hashes <- forward_hashes

rcfg <- jsonlite::read_json("scripts/moon/config/factor1_lipid_to_tf_relaxed.json", simplifyVector = TRUE)
relaxed_hashes <- check_md5(rcfg)
stopifnot(unname(tools::md5sum(rcfg$selected$reference_fit)) == rcfg$selected$reference_fit_md5,
  unname(tools::md5sum(file.path(rcfg$selected$reference_solution, "network_SIF.csv"))) ==
    rcfg$selected$reference_solution_sif_md5)
rfit <- readRDS(file.path(rcfg$output, "run_state.rds"))
rsaved <- readRDS(file.path(rcfg$output, "solution.rds"))
relaxed <- relaxed_compare_reference(rfit, rsaved, rcfg)
seed_table <- utils::read.csv(file.path(rcfg$output, "candidate_sign_consensus_seeds.csv"))
seeds <- stats::setNames(seed_table$copied_measurement_sign, seed_table$source)
rreplayed <- relaxed_extract_score_layer_paths(rfit$scores, rfit$pkn_final, seeds,
  rcfg$selected$primary, rcfg$selected$secondary, data$rna)
relaxed_validate_score_layer_solution(rreplayed, rfit, seeds,
  rcfg$selected$primary, rcfg$selected$secondary)
stopifnot(setequal(relaxed_edge_key(rreplayed$SIF), relaxed_edge_key(rsaved$SIF)))
relaxed$saved_extraction_identical <- TRUE
relaxed$edges <- nrow(rsaved$SIF)
relaxed$nodes <- nrow(rsaved$ATT)
relaxed$dependency_hashes <- relaxed_hashes

reverse_cfg <- jsonlite::read_json("scripts/moon/config/factor1_lipid_to_tf.json", simplifyVector = TRUE)
stopifnot(identical(reverse_cfg$inputs, cfg$inputs), dir.exists(reverse_cfg$input_directory),
  dir.exists(reverse_cfg$pkn_directory))
software <- data.frame(package = names(cfg$packages), required = unname(unlist(cfg$packages)),
  installed = vapply(names(cfg$packages), function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1)))
software$matches <- !is.na(software$installed) & software$required == software$installed
commit <- utils::packageDescription("cosmosR")$RemoteSha
report <- list(passed = TRUE, scope = "Saved-artifact hashes, reference comparisons and extraction replay only",
  model_refitted = FALSE, new_scores_computed = FALSE, frozen_artifacts_written = FALSE,
  forward = forward, relaxed = relaxed, reverse_dependency_contract = TRUE,
  software = software, cosmos_commit_matches = identical(commit, cfg$cosmos_commit),
  full_fit_environment_ready = all(software$matches) && identical(commit, cfg$cosmos_commit),
  full_fit_contract = "Unchanged: full-fit validate-only still requires the pinned software versions")
if (length(args)) {
  output <- sub("^--output=", "", args)
  if (file.exists(output)) stop("Refuse to overwrite audit; choose a new --output")
  jsonlite::write_json(report, output, pretty = TRUE, auto_unbox = TRUE, digits = NA)
}
cat(jsonlite::toJSON(report, pretty = TRUE, auto_unbox = TRUE, digits = NA), "\n")
