source("scripts/config/analysis_release.R")
# Build the Workflow 07b scientific report from complete saved results.
# This helper does not refit models or alter result tables. The caller provides
# its completed output directory and makes jsonlite available in .libPaths().
wp1_write_report <- function(out) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required.")
  out <- normalizePath(out, mustWork = TRUE)
  read_result <- function(name) {
    path <- file.path(out, name)
    if (!file.exists(path)) stop("Missing report input: ", path)
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  provenance_file <- file.path(out, "provenance.json")
  if (!file.exists(provenance_file)) stop("Missing report input: ", provenance_file)
  provenance <- jsonlite::fromJSON(provenance_file)
  concordance <- read_result("concordance_primary_and_secondary.csv")
  primary <- read_result("differential_primary_all.csv")
  sensitivities <- read_result("differential_sensitivities_all.csv")
  discovery <- read_result("differential_discovery_summary.csv")
  samples <- read_result("sample_manifest.csv")
  features <- read_result("feature_annotation_and_QC_ledger.csv")
  reference <- read_result("fixed_reference_and_matching_ledger.csv")
  reference_mice <- read_result("FCG_paired_reference_mouse_manifest.csv")
  loo_ko <- read_result("leave_one_KO_mouse_out.csv")
  loo_fcg <- read_result("leave_one_FCG_mouse_out.csv")
  ulm <- read_result("signed_ULM_descriptive_scores.csv")
  sens_concordance <- read_result("concordance_sensitivities.csv")
  fmt <- function(x, digits = 3L) formatC(as.numeric(x), format = "f", digits = digits)
  sci <- function(x) formatC(as.numeric(x), format = "e", digits = 2L)
  label <- c(pooled_KO_vs_WT = "Pooled KO − WT", female_KO_vs_WT = "Female KO − WT",
             male_KO_vs_WT = "Male KO − WT", interaction_M_minus_F = "Male − female KO effect")
  main_contrasts <- c("pooled_KO_vs_WT", "female_KO_vs_WT", "male_KO_vs_WT")
  main <- concordance[concordance$primary_test, , drop = FALSE]
  main <- main[match(main_contrasts, main$contrast), , drop = FALSE]
  if (nrow(main) != 3L || anyNA(main$contrast) ||
      !identical(main$reference_set, c("pooled", "female", "male")) ||
      any(!is.finite(main$permutation_fdr)) ||
      any(abs(main$concordance_score + main$pearson_r) > 1e-10)) {
    stop("The report requires the complete three-test primary concordance family.")
  }
  # Verify discovery statements directly from complete tables, not prose.
  actual_discoveries <- vapply(discovery$contrast, function(ct) {
    sum(primary$adj.P.Val[primary$contrast == ct] <= 0.05)
  }, numeric(1L))
  if (!identical(as.numeric(discovery$FDR_le_005), unname(actual_discoveries))) {
    stop("Discovery summary disagrees with the complete primary result tables.")
  }
  sensitivity_models <- unique(sensitivities$model)
  sensitivity_counts <- lapply(sensitivity_models, function(model) {
    subset <- sensitivities[sensitivities$model == model, , drop = FALSE]
    counts <- vapply(unique(subset$contrast), function(ct) {
      sum(subset$adj.P.Val[subset$contrast == ct] <= 0.05)
    }, numeric(1L))
    list(model = model, features = length(unique(subset$feature_id)),
         contrasts = length(counts), maximum_discoveries = max(counts))
  })
  sample_counts <- table(factor(samples$group, levels = c("F_WT", "F_KO", "M_WT", "M_KO")))
  exclusion_counts <- table(features$exclusion_reason)
  exclusions <- function(reason) {
    if (reason %in% names(exclusion_counts)) as.integer(exclusion_counts[[reason]]) else 0L
  }
  if (nrow(samples) != 32L || !identical(as.integer(sample_counts), c(8L, 8L, 7L, 9L)) ||
      sum(features$keep_primary) != provenance$primary_features ||
      sum(features$keep_qc30) != provenance$QC30_features ||
      sum(reference$match_status == "unique_exact") != provenance$shared_reference_features) {
    stop("Sample, feature or reference manifests disagree with the completed-run provenance.")
  }
  main_rows <- vapply(seq_len(nrow(main)), function(i) {
    x <- main[i, ]
    paste0("| ", label[[x$contrast]], " / ", x$reference_set, " reference | ",
           fmt(x$pearson_r), " | ", fmt(x$concordance_score), " | [",
           fmt(x$bootstrap_ci_low), ", ", fmt(x$bootstrap_ci_high), "] | ",
           fmt(x$permutation_p_two_sided, 4L), " | ", fmt(x$permutation_fdr, 4L), " |")
  }, character(1L))
  discovery_rows <- vapply(seq_len(nrow(discovery)), function(i) {
    x <- discovery[i, ]
    paste0("| ", label[[x$contrast]], " | ", x$tested, " | ", x$FDR_le_005,
           " | ", x$KO_higher, " | ", x$KO_lower, " |")
  }, character(1L))
  model_labels <- c(unweighted_log2 = "Unweighted joint log2", qc_CV_le30_log2 = "QC CV ≤30%, weighted log2",
                    supplied_scale_sex_weighted = "Supplied-scale weighted joint", separate_sex_log2 = "Separate-sex log2")
  sensitivity_rows <- vapply(sensitivity_counts, function(x) {
    paste0("| ", if (x$model %in% names(model_labels)) model_labels[[x$model]] else x$model,
           " | ", x$features, " | ", x$contrasts, " | ", x$maximum_discoveries, " |")
  }, character(1L))
  loo_rows <- vapply(seq_len(nrow(main)), function(i) {
    key <- paste(main$contrast[i], main$reference_set[i], sep = "__")
    if (!key %in% names(loo_ko) || !key %in% names(loo_fcg)) stop("Incomplete influence tables.")
    a <- range(loo_ko[[key]]); b <- range(loo_fcg[[key]])
    paste0("| ", label[[main$contrast[i]]], " | [", fmt(a[1]), ", ", fmt(a[2]),
           "] | [", fmt(b[1]), ", ", fmt(b[2]), "] |")
  }, character(1L))
  influential <- loo_ko[loo_ko$omitted_mouse == "Cyp-15_M-Ctrl", , drop = FALSE]
  if (nrow(influential) != 1L) stop("Expected the retained Cyp-15_M-Ctrl influence record.")
  ulm_male <- ulm[ulm$condition == "male_KO_vs_WT" & ulm$reference_set == "male", , drop = FALSE]
  if (nrow(ulm_male) != 1L) stop("Expected one primary male signed-ULM row.")
  secondary <- concordance[!concordance$primary_test, , drop = FALSE]
  secondary_rows <- vapply(seq_len(nrow(secondary)), function(i) {
    x <- secondary[i, ]
    paste0(label[[x$contrast]], ": r = ", fmt(x$pearson_r), ", score = ", fmt(x$concordance_score),
           ", two-sided animal p = ", fmt(x$permutation_p_two_sided, 4L))
  }, character(1L))
  historical_file <- file.path(out, "..", "global_network_validation_scores.csv")
  if (!file.exists(historical_file)) historical_file <- analysis_resolve_recorded_path("results/validation/global_network_validation_scores.csv")
  if (!file.exists(historical_file)) stop("Missing existing Workflow 07 siRNA validation context.")
  historical <- utils::read.csv(historical_file, check.names = FALSE, stringsAsFactors = FALSE)
  historical <- historical[historical$source == "Cyp2c38" & historical$condition == "cyp2c38", , drop = FALSE]
  if (nrow(historical) != 1L) stop("Expected one saved on-target Cyp2c38 siRNA ULM result.")
  relative_path <- function(from, target) {
    a <- strsplit(normalizePath(from, mustWork = TRUE), "/", fixed = TRUE)[[1L]]
    b <- strsplit(normalizePath(target, mustWork = TRUE), "/", fixed = TRUE)[[1L]]
    common <- 0L
    while (common < min(length(a), length(b)) && a[common + 1L] == b[common + 1L]) common <- common + 1L
    paste(c(rep("..", length(a) - common), b[seq.int(common + 1L, length(b))]), collapse = "/")
  }
  historical_link <- relative_path(out, historical_file)
  qc_shared <- sum(reference$feature_id %in% features$feature_id[features$keep_qc30])
  intervals_include_zero <- main$bootstrap_ci_low <= 0 & main$bootstrap_ci_high >= 0
  interval_statement <- if (all(intervals_include_zero)) "All three intervals include zero in this run." else
    paste0(sum(intervals_include_zero), " of the three intervals include zero in this run.")
  no_concordance <- all(main$permutation_fdr > 0.05)
  no_discoveries <- all(discovery$FDR_le_005 == 0L)
  no_sensitivity_discoveries <- all(vapply(sensitivity_counts, function(x) x$maximum_discoveries == 0L, logical(1L)))
  result_statement <- if (no_concordance && no_discoveries) {
    "**This dataset does not establish robust in vivo validation of the fixed Cyp2c38–lipid predictions.** None of the three primary animal-calibrated concordance tests passes BH FDR 0.05, and no individual primary lipid contrast produces an FDR discovery. This is inconclusive evidence, not proof that Cyp2c38 has no lipid effect."
  } else {
    "The complete primary discovery and animal-calibrated concordance results are tabulated below. Interpret lipid discoveries and predictive validation separately; statistical association does not establish direct substrates or flux."
  }
  female_weight <- unique(samples$precision_weight[samples$sex == "F"])
  male_weight <- unique(samples$precision_weight[samples$sex == "M"])
  if (length(female_weight) != 1L || length(male_weight) != 1L) stop("Expected one precision weight per sex.")
  invalid_perm <- provenance$permutations - main$permutation_valid_draws
  invalid_boot <- provenance$bootstraps - main$bootstrap_valid_draws
  rerun <- paste0("Rscript scripts/workflows/07b_cyp2c38_in_vivo_validation.R \\\n",
                  "  --permutations=", provenance$permutations, " --bootstraps=", provenance$bootstraps,
                  " --seed=", provenance$seed, " --workers=", provenance$workers, " \\\n",
                  "  --output=results/validation/cyp2c38_in_vivo_replay")
  lines <- c(
    "# Cyp2c38 in vivo heterozygous knockout: lipidomics and predictive validation",
    "", paste0("Completed run: ", provenance$completed_at, ". All numerical statements below are read from the saved result tables."),
    "", result_statement,
    "", "## Primary results", "",
    "The validation statistic compares the original gene–lipid Pearson coefficient with the new KO-minus-WT moderated lipid t-statistic. Negative Pearson r has the predicted perturbation direction. The signed concordance score is **−r**, so positive scores support agreement. Bootstrap intervals below apply to that score, not directly to the displayed r.",
    "", "| KO contrast / original reference | Pearson r | Signed score | Combined-bootstrap 95% score interval | Animal permutation p, two-sided | BH FDR, three tests |",
    "| --- | ---: | ---: | --- | ---: | ---: |", main_rows,
    "", paste0("Each primary comparison uses ", provenance$shared_reference_features, " fixed shared lipid features. ", interval_statement, " These intervals independently resample animals in both the KO and original paired reference cohorts; they remain conditional on the established feature selection."),
    "", "| Primary contrast | Features tested | FDR ≤0.05 | Higher after KO | Lower after KO |",
    "| --- | ---: | ---: | ---: | ---: |", discovery_rows,
    "", "FDR is BH-adjusted across all tested lipid features separately within each differential contrast. The interaction tests whether the male KO effect differs from the female KO effect; it is exploratory and is not a fourth predictive-validation test. Nominally small feature P values are not counted as discoveries.",
    "", "## Sensitivities and influence", "",
    "| Differential sensitivity | Features | Contrasts available | Maximum FDR discoveries in any contrast |",
    "| --- | ---: | ---: | ---: |", sensitivity_rows,
    "", if (no_sensitivity_discoveries) paste0("Every sensitivity model also has zero lipid discoveries at FDR ≤0.05. The QC restriction leaves ", qc_shared, " of the ", provenance$shared_reference_features, " shared reference features. Sensitivity concordance is descriptive: the animal-calibrated three-test family is the primary weighted-log2 analysis, and favorable secondary scores are not substitutes for its outcome.") else "Sensitivity discoveries are shown explicitly above; their results do not replace the prespecified primary family.",
    "", "The supplied-scale sensitivity uses the original positive numerical intensities without log2. Although limma names its coefficient column `logFC`, its values in that sensitivity are **linear mean-intensity differences**, not log fold changes. These units must not be pooled with or exponentiated as log2 effects.",
    "", "| Main signed score | Range omitting one KO-cohort animal | Range omitting one paired FCG animal |",
    "| --- | --- | --- |", loo_rows,
    "", paste0("Omitting male control **Cyp-15_M-Ctrl** changes the pooled score to ", fmt(influential$pooled_KO_vs_WT__pooled),
                " and the male sex-matched score to ", fmt(influential$male_KO_vs_WT__male),
                ", reversing their observed positive signs. This is an observational influence finding, not evidence that the animal is technically invalid. It remains in all primary analyses. Female reference leave-one-out scores also cross zero. Sample identities and all omission results are retained for review."),
    "", paste0("The male sex-matched ULM score is ", fmt(ulm_male$score), " with edge-level p = ", sci(ulm_male$p_value),
                ". That apparent significance treats lipid targets as independent and is **descriptive**. The whole-animal permutation p for the corresponding main comparison is ", fmt(main$permutation_p_two_sided[3], 4L), "; it governs the inferential conclusion."),
    "", paste0("Secondary sex-specific effects versus the pooled reference are also reported: ", paste(secondary_rows, collapse = "; "),
                ". These comparisons have no BH value in the primary three-test family. Their direction differs between sexes; reporting both prevents selecting only the more favorable pooled-reference comparison."),
    "", "## Samples, matching and preprocessing", "",
    paste0("There are **32 biological mice**: female WT 8, female KO 8, male WT 7 and male KO 9. The supplied control labels are interpreted as WT as described in the meeting. Analytical QC injections are not biological replicates. From ", nrow(features), " supplied measurement rows, exclude ",
           exclusions("internal_standard"), " explicit internal standards, ", exclusions("isotope_labelled_feature"), " additional isotope-labelled feature and ", exclusions("repository_excluded_class"),
           " rows in the existing Workflow 07 excluded classes. This leaves **", provenance$primary_features, " primary features**; **", provenance$QC30_features, "** meet QC CV ≤30% across nine analytical QC injections."),
    "", "Primary models use log2 supplied positive intensities. No new normalization, gap filling or imputation is applied. The prior normalization and processing provenance is unavailable; finite positive exported values do not prove that they were never imputed upstream. Retain complete original annotations, QC CV, feature IDs and the biological sample manifest.",
    "", paste0("The fixed prediction universe contains **", nrow(reference), " Cyp2c38 edges** from the original global full MOFA export. Of these, **", provenance$shared_reference_features, "** have unique exact KO annotation matches and ", nrow(reference) - provenance$shared_reference_features,
                " remain unmatched. Matching only takes the last pipe-delimited annotation and converts class/space formatting to parentheses. Chain order, slash versus underscore, ether prefixes, oxygenation and sum-versus-resolved molecular structure are preserved. Two cardiolipin reorder candidates remain excluded; no ambiguous structure is forced into a match."),
    "", paste0("Original reference measurements are already scaled hepatocyte RNA paired with log2 lipidomics from **", provenance$reference_mice, " FCG mice: 21 ovary-bearing and 22 testes-bearing**. Pooled coefficients reproduce all archived edges with maximum discrepancy **", sci(provenance$checks$reference_replay_max_difference),
                "**, below the required 1e−12 tolerance. Within-gonad coefficients are recomputed for the same universe. These estimate a different association and are generally weaker than pooled coefficients, which can include between-gonad covariance. There is no KO RNA dataset and no MOFA refit."),
    "", "## Models and animal-level uncertainty", "",
    paste0("A joint four-group limma design estimates F_WT, F_KO, M_WT and M_KO means. Female effect is F_KO−F_WT; male effect is M_KO−M_WT; pooled effect is **0.5×(female effect+male effect)**; interaction is male effect−female effect. The pooled estimand gives the sexes equal weight despite unequal male genotype counts. It is not an average of t-statistics. Sex-group precision weights from `arrayWeights(var.group=sex)` are ", fmt(female_weight, 4L), " for female and ", fmt(male_weight, 4L), " for male samples. These weights affect precision, not the pooled estimand. Weighted/unweighted group-mean coefficients agree to ", sci(provenance$checks$weighted_unweighted_logFC_max_difference), "."),
    "", "Models use robust, intensity-trend empirical Bayes moderation. Full tables contain effects, confidence intervals, moderated t, P values and BH FDR. Separate-sex sensitivity fits moderate each sex independently. Per-lipid signed consistency is `−reference_r × KO_t`; it is an agreement score, not a direct causal or enzyme-substrate estimate.",
    "", paste0("For inference, **", provenance$permutations, " genotype permutations within sex** preserve the original counts and refit precision weights and moderation. Two-sided empirical p values use the +1 numerator/denominator correction; BH adjusts only the three primary tests. **", provenance$bootstraps,
                " combined whole-animal bootstraps** independently resample KO mice within four sex/genotype groups and paired FCG mice within XXO=10, XYO=11, XXT=11 and XYT=11. Reference RNA and lipids remain paired. Established feature selection is held fixed. Percentile 2.5%/97.5% score limits are conditional intervals."),
    "", paste0("For the three main comparisons, invalid permutation draws are ", paste(invalid_perm, collapse = "/"), " and invalid bootstrap draws are ", paste(invalid_boot, collapse = "/"),
                " (pooled/female/male). Seed ", provenance$seed, " and ", provenance$workers, " workers are recorded; each draw has a saved seed. Whole-animal leave-one-out checks refit KO models or recompute reference correlations as appropriate."),
    "", "Permutation inference assumes genotype labels are exchangeable within sex. Unavailable litter/batch metadata cannot be controlled. Small cohorts, heterozygous protein reduction, upstream processing uncertainty and dependent lipid species limit interpretation. A null result here neither establishes biological equivalence nor rules out smaller, context-dependent or indirect effects. This in vivo heterozygous experiment is distinct from Workflow 07’s earlier siRNA knockdown in primary hepatocytes. Neither experiment alone identifies direct substrates, metabolic flux or the full impact of complete deletion.",
    "", paste0("For context, the [existing Workflow 07 siRNA scores](", historical_link, ") contain an on-target Cyp2c38 row with ULM score ", fmt(historical$score, 5L), " and p = ", sci(historical$p_value), ". This is also an edge-based descriptive p value. Different interventions and experimental settings mean comparing these p values does not quantify relative validation strength."),
    "", "## Artifacts", "",
    "- Inputs and QC: [sample manifest](sample_manifest.csv), [complete annotation/QC ledger](feature_annotation_and_QC_ledger.csv), [log2 analysis matrix](log2_analysis_matrix.csv), [reference matching ledger](fixed_reference_and_matching_ledger.csv), [paired FCG mouse manifest](FCG_paired_reference_mouse_manifest.csv), [contrast matrix](contrast_matrix.csv).",
    "- Differential results: [complete primary table](differential_primary_all.csv), [t-statistic matrix](t_statistic_matrix.csv), [discovery summary](differential_discovery_summary.csv), [pooled](differential_pooled_KO_vs_WT.csv), [female](differential_female_KO_vs_WT.csv), [male](differential_male_KO_vs_WT.csv), [interaction](differential_interaction_M_minus_F.csv), [all sensitivities](differential_sensitivities_all.csv), [saved model fits](limma_fits.rds).",
    "- Validation: [main/secondary concordance](concordance_primary_and_secondary.csv), [sensitivity concordance](concordance_sensitivities.csv), [per-lipid predictive consistency](per_lipid_predictive_consistency.csv), [sensitivity consistency](per_lipid_consistency_sensitivities.csv), [descriptive signed ULM](signed_ULM_descriptive_scores.csv), [descriptive lipid-class summary](lipid_class_descriptive_summary.csv).",
    "- Uncertainty: [all permutation draws](permutation_draws.csv), [all combined-bootstrap draws](combined_bootstrap_draws.csv), [KO omission checks](leave_one_KO_mouse_out.csv), [FCG omission checks](leave_one_FCG_mouse_out.csv).",
    "- Figures: [QC/PCA and model diagnostics](plots/QC_PCA.pdf), [differential analysis](plots/differential.pdf), [predictive validation](plots/concordance.pdf).",
    "- Provenance: [run settings/checks/package versions](provenance.json), [input/code SHA256 manifest](input_and_code_SHA256.csv), [R session information](sessionInfo.txt).",
    "", "## Reproduction", "",
    "Run from the analysis repository root with the recorded packages available. Use a fresh replay destination to preserve this completed run:",
    "", "```sh", rerun, "```", "",
    "Run the meaningful model/sign/invariant checks separately:", "", "```sh",
    "Rscript tests/test_cyp2c38_in_vivo.R", "```", "",
    "The workflow’s `--validate-only` option checks intake and original-reference reconstruction without fitting KO effects. The saved SHA256 manifest and provenance file identify the exact inputs, code, software versions and resampling settings. Retain matching/reporting conventions when interpreting a replay."
  )
  # Do not emit broken artifact links when invoked on an incomplete run.
  links <- regmatches(paste(lines, collapse = "\n"),
                      gregexpr("\\]\\(([^)]+)\\)", paste(lines, collapse = "\n")))[[1L]]
  targets <- sub("^\\]\\((.*)\\)$", "\\1", links)
  missing_artifacts <- targets[!file.exists(file.path(out, targets))]
  if (length(missing_artifacts)) stop("Missing linked report artifact(s): ", paste(missing_artifacts, collapse = ", "))
  destination <- file.path(out, "README.md")
  writeLines(lines, destination, useBytes = TRUE)
  cat("Wrote Cyp2c38 scientific report: ", destination, "\n", sep = "")
  invisible(destination)
}
