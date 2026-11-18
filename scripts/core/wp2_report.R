# Package-wide WP2 reporting. The existing component generators read completed
# CSV/RDS outputs; this composer preserves their numerical and methodological
# content while placing pooled Figure 3 first and zonation/diet second.
source("scripts/config/analysis_release.R")

wp2_report_section <- function(lines, heading) {
  start <- match(paste0("## ", heading), lines)
  if (is.na(start)) stop("Missing component report section: ", heading)
  next_heading <- which(startsWith(lines, "## ") & seq_along(lines) > start)
  end <- if (length(next_heading)) next_heading[1L] - 1L else length(lines)
  body <- lines[seq.int(start + 1L, end)]
  while (length(body) && !nzchar(body[1L])) body <- body[-1L]
  while (length(body) && !nzchar(body[length(body)])) body <- body[-length(body)]
  body
}

wp2_report_prefix_links <- function(lines, prefix) {
  gsub("(\\]\\()(?!(?:[a-z]+:|/|#))", paste0("\\1", prefix), lines, perl = TRUE)
}

wp2_report_subsection <- function(heading, body) {
  # Promote original subsection titles below the new package-level heading.
  body <- sub("^### ", "#### ", body)
  c(paste0("### ", heading), "", body, "")
}

wp2_report_stage <- function(source_dir, destination) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  paths <- list.files(source_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  paths <- paths[basename(paths) != "README.md"]
  ok <- file.symlink(paths, file.path(destination, basename(paths)))
  if (!all(ok)) stop("Cannot stage completed report inputs: ", paste(paths[!ok], collapse = ", "))
}

wp2_report_relative_input <- function(path, output_root) {
  source <- strsplit(normalizePath(path, mustWork = TRUE), "/", fixed = TRUE)[[1L]]
  destination <- strsplit(normalizePath(output_root, mustWork = TRUE), "/", fixed = TRUE)[[1L]]
  common <- 0L
  for (i in seq_len(min(length(source), length(destination)))) {
    if (!identical(source[i], destination[i])) break
    common <- i
  }
  paste(c(rep("..", length(destination) - common), source[seq.int(common + 1L, length(source))]), collapse = "/")
}

#' Refresh the single pooled-first WP2 summary without rerunning statistics.
#' All completed component outputs must exist below the supplied output root.
#' The original root zonation report is archived once before replacement.
wp2_write_package_report <- function(output_root = analysis_path("historical_wf12"),
                                     report_date = format(Sys.Date(), "%d %B %Y")) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  pooled_root <- file.path(output_root, "followup")
  if (!dir.exists(pooled_root)) stop("Completed pooled outputs are required in ", pooled_root)
  required <- c("baseline_inference.csv", "network_analysis.rds",
                "followup/pooled_inference.csv", "followup/pooled_analysis.rds")
  if (any(!file.exists(file.path(output_root, required)))) {
    stop("Both WP2 analysis phases must be complete before composing the package report.")
  }
  stage <- tempfile("wp2_package_report_")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  zone_stage <- file.path(stage, "zonation")
  pool_stage <- file.path(stage, "pooled")
  wp2_report_stage(output_root, zone_stage)
  wp2_report_stage(pooled_root, pool_stage)
  generators <- new.env(parent = globalenv())
  sys.source("scripts/core/wdsd_presentation.R", envir = generators)
  sys.source("scripts/core/wdsd_pooled_presentation.R", envir = generators)
  generators$wp2_make_report(zone_stage)
  generators$wp2_pooled_make_report(pool_stage)
  zonal <- readLines(file.path(zone_stage, "README.md"), warn = FALSE)
  pooled <- wp2_report_prefix_links(readLines(file.path(pool_stage, "README.md"), warn = FALSE), "followup/")
  zsec <- function(name) wp2_report_section(zonal, name)
  psec <- function(name) wp2_report_section(pooled, name)
  split_paragraphs <- function(x) {
    pieces <- split(x, cumsum(!nzchar(x)))
    unname(lapply(pieces, function(z) z[nzchar(z)]))
  }
  p_findings <- split_paragraphs(psec("Main findings"))
  z_findings <- split_paragraphs(zsec("Main findings"))
  if (length(p_findings) < 2L || length(z_findings) < 2L) stop("Component findings contract changed.")
  pooled_interpretation <- psec("Interpretation and reproduction")
  zonal_interpretation <- zsec("Interpretation limits and manuscript use")
  before_limitations <- function(x) {
    end <- match("Recorded input limitations:", x)
    if (is.na(end)) stop("Missing recorded input limitations in a component report.")
    head(x, end - 1L)
  }
  limitations <- unique(c(
    readLines(file.path(output_root, "input_limitations.txt"), warn = FALSE),
    readLines(file.path(pooled_root, "input_limitations.txt"), warn = FALSE)))
  settings_zone <- readRDS(file.path(output_root, "network_analysis.rds"))$metadata
  settings_pool <- readRDS(file.path(pooled_root, "pooled_analysis.rds"))$metadata
  settings_equal <- identical(as.integer(settings_zone$seed), as.integer(settings_pool$seed)) &&
    identical(as.integer(settings_zone$B), as.integer(settings_pool$null_bootstrap_draws)) &&
    identical(as.integer(settings_zone$paired_B), as.integer(settings_pool$paired_stability_draws))
  if (!settings_equal) stop("A single package reproduction command requires shared seed/bootstrap settings.")
  reproduction_command <- paste0(
    "Rscript scripts/workflows/12_wdsd_network_comparison.R --permutations=", settings_pool$random_permutations,
    " --bootstraps=", settings_pool$null_bootstrap_draws,
    " --paired-bootstraps=", settings_pool$paired_stability_draws,
    " --seed=", settings_pool$seed, " --output=results/wdsd_network_comparison_reproduction")
  input_link <- function(label, path) paste0("[", label, "](", wp2_report_relative_input(path, output_root), ")")
  input_manifest_note <- paste0("Repository input provenance: ",
    input_link("input guide", "data/WP1_WP2_INPUTS.md"), " and ",
    input_link("source/destination SHA256 manifest", "data/wp1_wp2_input_migration_manifest.csv"),
    " document the byte-preserving copies and original archive locations.")
  input_files_note <- paste0("WP2 local inputs: ",
    input_link("whole-hepatocyte Seurat pseudobulks", "data/RNA/wdsd/wdsd_hep_pseudobulks.rds"), "; ",
    input_link("frozen two-marker paired zonation Seurat pseudobulks", analysis_path("historical_wdsd_zonation_transport")), "; ",
    input_link("CF453 SD/WD lipidomics export", "data/metabolomics/wdsd/lipidomics_cf453.csv"), "; ",
    input_link("raw FCG lipid annotation export", "data/metabolomics/fcg/lipidomics_fcg.csv"), "; and ",
    input_link("cleaned FCG lipid matrix", "data/metabolomics/Lipid_data_file_JPC_cleaned.csv"),
    ". The raw FCG annotation table supplements the fixed cleaned matrix rather than replacing it.")
  report <- c(
    "# WP2: pooled and zonation network conservation, followed by diet-associated changes", "",
    paste0("Report reorganized: ", report_date, ". The numerical analyses use the recorded 7 October 2026 run settings; this report combines both completed WP2 analyses without changing their results."), "",
    "## Overview", "", p_findings[[2L]], "", z_findings[[1L]], "", z_findings[[2L]], "",
    "The pooled whole-hepatocyte Figure 3 analysis is presented first, followed by the four sex-specific zonation networks from Figure 5 and their standard-to-Western-diet changes. These comparisons answer different questions: the pooled pattern includes between-sex covariance, while each zonation comparison uses one sex and seven SD mice. The lipid annotation audit and shared quality controls follow both network analyses. Native ULM scores and nominal p-values are retained beside explicitly defined mouse-level null tests; diet changes remain descriptive.", "",
    "## 1. Pooled whole-hepatocyte network comparison: Figure 3", "",
    wp2_report_subsection("Frozen network, original selection and cohort coverage", psec("2. Which Figure 3 network was compared?")),
    wp2_report_subsection("Pooled coefficient conservation and animal-level null tests", psec("Pooled coefficient conservation")),
    wp2_report_subsection("Removing the between-sex means", psec("Removing the between-sex means")),
    wp2_report_subsection("Four predefined Cyp genes", psec("Four predefined Cyp genes")),
    "## 2. Sex-specific zonation networks and diet comparison: Figure 5", "",
    wp2_report_subsection("Frozen Figure 5 reference and molecular coverage", zsec("Frozen Figure 5 reference and molecular coverage")),
    wp2_report_subsection("Standard-diet conservation within each sex and zone", zsec("Stage A: standard-diet conservation")),
    wp2_report_subsection("Western versus standard diet within CF453", zsec("Stage B: Western versus standard diet")),
    wp2_report_subsection("Predefined Mtmr11 focus", zsec("Predefined Mtmr11 focus")),
    "## 3. Lipid annotation overlap and unresolved molecular identity", "", p_findings[[1L]], "",
    psec("1. Why does the lipid overlap appear limited?"), "",
    "## 4. Shared cohort, processing, QC and sensitivity analyses", "",
    wp2_report_subsection("Cohort, processing and QC", zsec("Cohort, processing and QC")),
    wp2_report_subsection("Pooled molecular mapping and analytical-quality sensitivities", psec("Molecular mapping and analytical-quality sensitivities")),
    wp2_report_subsection("Zonation quality-field and RNA-layer sensitivities", zsec("Prespecified sensitivities")),
    "## 5. Interpretation and remaining input limitations", "",
    before_limitations(pooled_interpretation), "", before_limitations(zonal_interpretation), "",
    "Recorded input limitations:", "", paste0("- ", limitations), "",
    "## 6. Artifacts, local inputs and reproduction", "",
    "Workflow 12 is the single WP2 entry point. It executes pooled Figure 3 / annotation analysis first, then the four sex/zone baseline and diet comparisons, and finally writes this package report. Processed inputs are stored in the repository data tree; the data manifest records source locations and integrity checks. Shared helpers retain the original calculations and frozen historical feature membership.", "",
    input_manifest_note, "", input_files_note, "",
    "- Pooled: [inference](followup/pooled_inference.csv), [complete matched coefficients](followup/matched_pooled_edges.csv), [coverage](followup/pooled_coverage.csv), [reference replay](followup/pooled_reference_replay.csv), [mapping sensitivities](followup/mapping_sensitivity_inference.csv) and [QC sensitivities](followup/pooled_QC_sensitivity.csv).",
    "- Zonation/diet: [baseline inference](baseline_inference.csv), [diet summaries](diet_comparison_summary.csv), [complete species coefficients](network_edges.csv), [joint omission summaries](loco_aggregate_summary.csv), [sensitivity results](sensitivities.csv) and [Figure 6 draft](plots/Figure6_draft.pdf).",
    "- Pairing and provenance: [all 28 mice](metadata.csv), [14 pooled SD mice](followup/SD_mouse_manifest.csv), [43 pooled FCG mice](followup/FCG_mouse_manifest.csv), [zonal reference mice](reference_mouse_manifest.csv), [strict lipid crosswalk](lipid_identity_crosswalk.csv) and [all molecular candidates](followup/lipid_overlap_audit/candidate_pairs.csv).",
    "- Figures: [pooled figure manifest](followup/plot_manifest.csv) and [zonation/diet figure manifest](plot_manifest.csv). The combined PDF preserves all scientific figure pages and presents pooled comparisons before zonation/diet comparisons, followed by annotation and shared QC diagnostics.", "",
    paste0("Recorded seed: ", settings_pool$seed, "; pooled random pairings per null: ", settings_pool$random_permutations,
           "; null-bootstrap draws: ", settings_pool$null_bootstrap_draws,
           "; paired stability draws: ", settings_pool$paired_stability_draws,
           "; exhaustive zonation pairings per sex: 5040. Identity and shared mouse indices remain in each phase's audit files."), "",
    "Run from the repository root. A full rerun requires a fresh output directory:", "", "```sh",
    "Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only",
    reproduction_command, "```", "",
    "Refresh the combined report from completed outputs without recomputing analyses:", "", "```sh",
    "Rscript scripts/workflows/12_wdsd_network_comparison.R --report-only",
    "python3 scripts/reports/build_wp2_report_pdf.py", "```", "",
    "Each analysis phase retains its input/code SHA256, settings, session versions, shared resampling indices and null draws. The original FCG results and supplied data are fixed inputs. The preceding zonation-only README is archived once as README_zonation_original_2026-10-07.md; the pooled component report remains in followup/README.md.")
  audit_files <- c("run_manifest.json", "input_and_code_SHA256.csv", "sessionInfo.txt")
  if (all(file.exists(file.path(output_root, audit_files))) && all(file.exists(file.path(pooled_root, audit_files)))) {
    report <- c(report, "", "Run audits: [zonal settings](run_manifest.json), [zonal input/code SHA256](input_and_code_SHA256.csv), [zonal session](sessionInfo.txt), [pooled settings](followup/run_manifest.json), [pooled input/code SHA256](followup/input_and_code_SHA256.csv) and [pooled session](followup/sessionInfo.txt).")
  }
  links <- regmatches(report, gregexpr("\\]\\(([^)]+)\\)", report, perl = TRUE))
  targets <- unlist(lapply(links, function(z) sub("^\\]\\((.*)\\)$", "\\1", z)), use.names = FALSE)
  local <- targets[!grepl("^[a-z]+://", targets)]
  absent <- local[!file.exists(file.path(output_root, local))]
  if (length(absent)) stop("Combined WP2 report references absent artifacts: ", paste(unique(absent), collapse = ", "))
  destination <- file.path(output_root, "README.md")
  archive <- file.path(output_root, "README_zonation_original_2026-10-07.md")
  if (file.exists(destination) && !file.exists(archive)) {
    old_title <- readLines(destination, n = 1L, warn = FALSE)
    if (!identical(old_title, report[1L]) && !file.copy(destination, archive, overwrite = FALSE)) {
      stop("Could not preserve the preceding zonation report.")
    }
  }
  writeLines(report, destination, useBytes = TRUE)
  cat("Combined pooled-first WP2 report saved: ", destination, "\n", sep = "")
  invisible(destination)
}
