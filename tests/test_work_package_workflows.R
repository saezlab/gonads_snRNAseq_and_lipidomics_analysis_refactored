source("scripts/config/analysis_release.R")
# Orchestration contract checks: phase order, fresh intake, CLI routing and guards.
# Run from the repository root: Rscript tests/test_work_package_workflows.R
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
expect_error <- function(expr, text) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, "error"), grepl(text, conditionMessage(err), fixed = TRUE))
}
get_option <- function(args, key) {
  hit <- args[startsWith(args, paste0("--", key, "="))]
  stopifnot(length(hit) == 1L)
  substring(hit, nchar(key) + 4L)
}

wp1 <- new.env(parent = globalenv())
sys.source("scripts/core/wp1_workflow.R", envir = wp1)
calls <- list()
raw_intake <- list(y = matrix(1, 1L, 1L), metadata = "raw", ledger = "raw")
wp1$wp1_run <- function(args) {
  calls[[length(calls) + 1L]] <<- list(phase = "intake", args = args)
  stopifnot("--validate-only" %in% args)
  raw_intake
}
wp1$wp1_ulm_run <- function(args, previous, intake) {
  calls[[length(calls) + 1L]] <<- list(phase = "ULM", args = args, previous = previous)
  stopifnot("--validate-only" %in% args, identical(intake, raw_intake))
  "validated"
}
wp1$wp1_package_run(c("--validate-only", "--output=never_created_wp1",
  "--bootstraps=29", "--base-bootstraps=19", "--permutations=39", "--workers=1"))
stopifnot(identical(vapply(calls, `[[`, character(1L), "phase"), c("intake", "ULM")),
  get_option(calls[[1L]]$args, "bootstraps") == "19",
  get_option(calls[[2L]]$args, "bootstraps") == "29",
  get_option(calls[[1L]]$args, "permutations") == "39",
  get_option(calls[[2L]]$args, "output") == "never_created_wp1/ulm_followup",
  calls[[2L]]$previous == "never_created_wp1", !dir.exists("never_created_wp1"))
expect_error(wp1$wp1_package_run("--unexpected"), "Unknown Workflow 11 option")
expect_error(wp1$wp1_package_run(c("--output=a", "--output=b")), "Duplicate option")

wp2 <- new.env(parent = globalenv())
sys.source("scripts/core/wp2_workflow.R", envir = wp2)
calls <- list()
wp2$wp2_pooled_run <- function(args) {
  calls[[length(calls) + 1L]] <<- list(phase = "pooled", args = args)
  stopifnot("--validate-only" %in% args)
  "pooled"
}
wp2$wp2_run <- function(args) {
  calls[[length(calls) + 1L]] <<- list(phase = "zonation", args = args)
  stopifnot("--validate-only" %in% args)
  "zonation"
}
wp2$wp2_package_run(c("--validate-only", "--output=never_created_wp2",
  "--permutations=29", "--bootstraps=39", "--paired-bootstraps=19"))
stopifnot(identical(vapply(calls, `[[`, character(1L), "phase"), c("pooled", "zonation")),
  get_option(calls[[1L]]$args, "output") == "never_created_wp2/followup",
  get_option(calls[[2L]]$args, "output") == "never_created_wp2",
  get_option(calls[[1L]]$args, "permutations") == "29",
  get_option(calls[[2L]]$args, "bootstraps") == "39",
  !any(startsWith(calls[[2L]]$args, "--permutations=")), !dir.exists("never_created_wp2"))
expect_error(wp2$wp2_package_run(c("--validate-only", "--report-only")), "Choose validation or report refresh")
expect_error(wp2$wp2_package_run(c("--output=a", "--output=b")), "Duplicate option")
existing <- tempfile("existing_wp2_")
dir.create(existing)
writeLines("preserve", file.path(existing, "existing_result.txt"))
before <- length(calls)
expect_error(wp2$wp2_package_run(paste0("--output=", existing)), "Output directory is not empty")
stopifnot(length(calls) == before, readLines(file.path(existing, "existing_result.txt")) == "preserve")
unlink(existing, recursive = TRUE)
cat("Work package orchestration contracts passed.\n")

# Optional actual-data replay checks after a fresh workflow run. WP1 compares
# observed estimates even when its smoke run uses fewer Monte Carlo draws;
# WP2 inference replay requires the original default draw counts and seed.
# Example: --wp1-output=/private/tmp/wp1_new --wp2-output=/private/tmp/wp2_new
args <- commandArgs(trailingOnly = TRUE)
if (any(!grepl("^--wp[12]-output=", args))) stop("Unknown replay-check option")
replay_path <- function(key) {
  hit <- args[startsWith(args, paste0("--", key, "="))]
  if (length(hit) > 1L) stop("Duplicate replay-check option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 4L) else NULL
}
compare_tables <- function(reference, regenerated, files) {
  for (name in files) {
    a <- utils::read.csv(file.path(reference, name), check.names = FALSE)
    b <- utils::read.csv(file.path(regenerated, name), check.names = FALSE)
    if (!isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))) {
      stop("Observed/output replay differs: ", name)
    }
    cat("Actual output replay passed: ", name, "\n", sep = "")
  }
}
wp1_out <- replay_path("wp1-output")
if (!is.null(wp1_out)) {
  old <- "results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo"
  compare_tables(old, wp1_out, c("differential_primary_all.csv", "t_statistic_matrix.csv",
    "feature_annotation_and_QC_ledger.csv", "fixed_reference_and_matching_ledger.csv"))
  columns <- c("reference_universe", "contrast", "reference_set", "n_matched_lipids",
    "ULM_score", "ULM_p_value")
  a <- utils::read.csv(file.path(old, "ulm_followup/ULM_summary.csv"))[, columns]
  b <- utils::read.csv(file.path(wp1_out, "ulm_followup/ULM_summary.csv"))[, columns]
  stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE)))
  cat("Actual WP1 native ULM replay passed.\n")
}
wp2_out <- replay_path("wp2-output")
if (!is.null(wp2_out)) {
  old <- analysis_path("historical_wf12")
  compare_tables(old, wp2_out, c("network_coverage.csv", "network_edges.csv",
    "baseline_inference.csv", "diet_comparison_summary.csv", "sensitivities.csv",
    "gene_diet_change_summary.csv", "lipid_class_diet_change_summary.csv"))
  compare_tables(file.path(old, "followup"), file.path(wp2_out, "followup"),
    c("pooled_coverage.csv", "matched_pooled_edges.csv", "pooled_inference.csv",
      "mapping_sensitivity_inference.csv", "pooled_QC_sensitivity.csv"))
}
