# ULM-focused figures and report for Workflow 07c.
# Sourcing defines functions without reading outputs, fitting models or plotting.

wp1_ulm_format_p <- function(x) {
  if (length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stop("A displayed p-value must be one finite probability.")
  }
  if (x < 0.001) formatC(x, format = "e", digits = 2L) else
    formatC(x, format = "f", digits = 4L)
}

wp1_ulm_summary_check <- function(summary) {
  columns <- c("reference_universe", "contrast", "reference_set",
               "n_matched_lipids", "ULM_score", "ULM_p_value",
               "bootstrap_p_KO_only", "bootstrap_p_KO_and_reference",
               "bootstrap_p_restricted_residuals", "permutation_p",
               "primary_comparison")
  if (!is.data.frame(summary) || !all(columns %in% names(summary))) {
    stop("The complete Workflow 07c ULM summary schema is required.")
  }
  keys <- paste(summary$reference_universe, summary$contrast,
                summary$reference_set, sep = "__")
  if (anyNA(keys) || anyDuplicated(keys) ||
      !all(summary$reference_universe %in% c("full", "fixed828")) ||
      !all(summary$reference_set %in% c("pooled", "female", "male")) ||
      anyNA(summary$primary_comparison)) {
    stop("ULM comparison keys and primary-comparison flags must be valid and unique.")
  }
  if (any(!is.finite(summary$ULM_score)) ||
      any(!is.finite(summary$n_matched_lipids)) || any(summary$n_matched_lipids < 3L)) {
    stop("All ULM scores and matched-lipid counts must be finite.")
  }
  p_columns <- c("ULM_p_value", "bootstrap_p_KO_only",
                 "bootstrap_p_KO_and_reference", "bootstrap_p_restricted_residuals",
                 "permutation_p")
  for (column in p_columns) {
    p <- summary[[column]]
    if (!is.numeric(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) {
      stop("Invalid probability in ULM summary column ", column, ".")
    }
  }
  expected <- expand.grid(reference_universe = c("full", "fixed828"),
                          stratum = c("pooled", "female", "male"),
                          stringsAsFactors = FALSE)
  expected_keys <- paste(expected$reference_universe,
                         paste0(expected$stratum, "_KO_vs_WT"),
                         expected$stratum, sep = "__")
  main_keys <- keys[summary$primary_comparison %in% TRUE]
  if (!setequal(main_keys, expected_keys)) {
    stop("Require the six complete pooled and sex-matched primary ULM comparisons.")
  }
  invisible(keys)
}

#' Plot original Cyp2c38 coefficients against KO-minus-WT moderated lipid t-values.
#'
#' Figures display the native ULM statistic and the supplied nominal and
#' whole-mouse bootstrap p-values. No statistic is sign-reversed for presentation.
#' @param out Completed workflow output root.
#' @param edges All matched lipid edges, keyed by universe, contrast and reference.
#' @param summary Complete observed ULM and resampling summary.
#' @return Absolute path to a multipage scientific figure PDF.
wp1_ulm_followup_plots <- function(out, edges, summary) {
  keys <- wp1_ulm_summary_check(summary)
  edge_columns <- c("reference_universe", "contrast", "reference_set",
                    "reference_r", "t", "canonical_lipid", "feature_id")
  if (!is.data.frame(edges) || !all(edge_columns %in% names(edges)) ||
      any(!is.finite(edges$reference_r)) || any(!is.finite(edges$t))) {
    stop("Figures require complete finite matched-edge annotations and statistics.")
  }
  edge_keys <- paste(edges$reference_universe, edges$contrast,
                     edges$reference_set, sep = "__")
  if (any(!edge_keys %in% keys) ||
      anyDuplicated(paste(edge_keys, edges$feature_id, sep = "__"))) {
    stop("Every figure edge must belong to a unique exported ULM comparison.")
  }
  counts <- tabulate(match(edge_keys, keys), nbins = length(keys))
  if (!identical(as.numeric(counts), as.numeric(summary$n_matched_lipids))) {
    stop("Figure edge counts must agree with observed ULM summary counts.")
  }
  for (i in seq_len(nrow(summary))) {
    data <- edges[edge_keys == keys[i], , drop = FALSE]
    slope <- summary(stats::lm(t ~ reference_r, data = data))$coefficients[2L, ]
    if (!is.finite(slope[3L]) || !is.finite(slope[4L]) ||
        abs(summary$ULM_score[i] - slope[3L]) > 1e-8 * (1 + abs(slope[3L])) ||
        abs(summary$ULM_p_value[i] - slope[4L]) > 1e-10) {
      stop("Figure edges must reproduce their exported native ULM score and p-value.")
    }
  }
  plot_dir <- file.path(out, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(plot_dir, "ULM_concordance.pdf")
  grDevices::pdf(destination, width = 11.7, height = 7.3, onefile = TRUE,
                 family = "Helvetica", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  titles <- c(pooled_KO_vs_WT = "Pooled KO - WT",
               female_KO_vs_WT = "Female KO - WT",
               male_KO_vs_WT = "Male KO - WT")
  reference_titles <- c(pooled = "Pooled FCG reference",
                         female = "Ovarian FCG reference",
                         male = "Testicular FCG reference")
  colors <- c(pooled = "#526b7e", female = "#ca4747", male = "#3c6daf")
  draw_page <- function(rows, title, subtitle, ncol) {
    graphics::par(mfrow = c(ceiling(length(rows) / ncol), ncol),
                  mar = c(5.3, 4.3, 4, 1), oma = c(3.1, 0.5, 3.2, 0.5),
                  las = 1)
    for (row in rows) {
      s <- summary[row, , drop = FALSE]
      data <- edges[edge_keys == keys[row], , drop = FALSE]
      # A single y scale per page makes differences between sex groups visible.
      page_keys <- keys[rows]
      page_y <- edges$t[edge_keys %in% page_keys]
      ymax <- max(abs(page_y), 1)
      yrange <- c(-ymax * 1.08, ymax * 1.45)
      color <- colors[[s$reference_set]]
      graphics::plot(data$reference_r, data$t,
                     pch = 16L, cex = 0.55,
                     col = grDevices::adjustcolor(color, alpha.f = 0.43),
                     xlim = c(-1, 1), ylim = yrange,
                     xlab = "Original Cyp2c38-lipid\nPearson coefficient",
                     ylab = "KO - WT moderated lipid t-statistic",
                     main = paste(titles[[s$contrast]],
                                  reference_titles[[s$reference_set]], sep = "\n"),
                     cex.main = 1.05, cex.lab = 0.92, cex.axis = 0.9)
      graphics::abline(h = 0, v = 0, col = "#dce3e8", lty = 3L)
      graphics::abline(stats::lm(t ~ reference_r, data = data),
                       col = color, lwd = 2L)
      annotation <- c(sprintf("ULM score = %.3f", s$ULM_score),
                       paste0("ULM p = ", wp1_ulm_format_p(s$ULM_p_value)),
                       paste0("Bootstrap KO p = ", wp1_ulm_format_p(s$bootstrap_p_KO_only)),
                       paste0("Bootstrap KO + FCG p = ",
                              wp1_ulm_format_p(s$bootstrap_p_KO_and_reference)),
                       paste0("Matched lipids: ", s$n_matched_lipids))
      graphics::legend("topleft", legend = annotation, bty = "o",
                       bg = "white", box.col = "#cad5df", text.col = "#263d50",
                       cex = 0.84, inset = c(0.015, 0.01), x.intersp = 0.1,
                       y.intersp = 1.08)
    }
    graphics::mtext(title, outer = TRUE, side = 3, line = 1.3,
                    font = 2L, cex = 1.1)
    graphics::mtext(subtitle, outer = TRUE, side = 3, line = 0.1, cex = 0.84)
    graphics::mtext("Negative ULM scores support the expected response to reduced Cyp2c38. All displayed p-values are two-sided.",
                    outer = TRUE, side = 1, line = 0.45, cex = 0.83)
    graphics::mtext("Bootstrap KO fixes the FCG reference; bootstrap KO + FCG additionally resamples paired reference mice.",
                    outer = TRUE, side = 1, line = 1.35, cex = 0.8)
  }
  primary_contrasts <- paste0(c("pooled", "female", "male"), "_KO_vs_WT")
  for (universe in c("full", "fixed828")) {
    rows <- which(summary$reference_universe == universe &
                    summary$primary_comparison %in% TRUE)
    rows <- rows[match(primary_contrasts, summary$contrast[rows])]
    if (anyNA(rows) || length(rows) != 3L) stop("Incomplete primary plot page.")
    full <- universe == "full"
    title <- if (full) "Cyp2c38 in vivo validation: regenerated full FCG networks" else
      "Cyp2c38 in vivo validation: historical feature-universe sensitivity"
    subtitle <- paste0(if (full) "All 1,254 original lipidomics features" else
                         "Historical 828-species MOFA network", "; ",
                       unique(summary$n_matched_lipids[rows]),
                       " strict molecular matches; ovarian reference n = 21, testicular reference n = 22")
    draw_page(rows, title, subtitle, ncol = 3L)
  }
  secondary <- which(!summary$primary_comparison &
                       summary$reference_set == "pooled" &
                       summary$contrast %in% c("female_KO_vs_WT", "male_KO_vs_WT"))
  if (length(secondary) == 4L) {
    order_key <- paste(summary$reference_universe[secondary],
                       summary$contrast[secondary], sep = "__")
    wanted <- c("full__female_KO_vs_WT", "full__male_KO_vs_WT",
                "fixed828__female_KO_vs_WT", "fixed828__male_KO_vs_WT")
    secondary <- secondary[match(wanted, order_key)]
    if (anyNA(secondary)) stop("Incomplete secondary figure page.")
    draw_page(secondary, "Secondary comparisons: sex-specific KO effects against the pooled FCG reference",
              "Full measured-species universe in the top row; historical feature universe in the bottom row",
              ncol = 2L)
  }
  cat("Wrote ULM concordance figure: ", destination, "\n", sep = "")
  invisible(normalizePath(destination, mustWork = TRUE))
}

#' Write a self-contained ULM-focused report from a completed 07c run.
#' Missing result artifacts or incomplete calibration inputs stop explicitly.
wp1_ulm_followup_report <- function(out) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required.")
  out <- normalizePath(out, mustWork = TRUE)
  read_result <- function(name) {
    path <- file.path(out, name)
    if (!file.exists(path)) stop("Missing completed-run report input: ", path)
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  provenance_file <- file.path(out, "provenance.json")
  if (!file.exists(provenance_file)) stop("Missing completed-run provenance: ", provenance_file)
  provenance <- jsonlite::fromJSON(provenance_file)
  settings <- c("completed_at", "seed", "bootstraps", "permutations", "workers")
  if (!all(settings %in% names(provenance))) stop("Report provenance settings are incomplete.")
  summary <- read_result("ULM_summary.csv")
  wp1_ulm_summary_check(summary)
  network_manifest <- read_result("network_manifest.csv")
  manifest_columns <- c("stratum", "gonadal_type", "universe", "n_edges", "n_mice", "file")
  if (!all(manifest_columns %in% names(network_manifest)) || nrow(network_manifest) != 6L) {
    stop("Report requires the complete six-network export manifest.")
  }
  null_diagnostics <- read_result("null_distribution_diagnostics.csv")
  null_columns <- c("reference_universe", "contrast", "reference_set",
                     "primary_comparison", "draws", "ULM_null_SD",
                     "nominal_ULM_p_le005_fraction", "method")
  if (!all(null_columns %in% names(null_diagnostics))) {
    stop("Report requires the completed fitted-null dependence diagnostics.")
  }
  calibration <- read_result("diagnostics/calibration_null_rejection_summary.csv")
  calibration_columns <- c("reference_universe", "contrast", "reference_set", "n_lipids",
                           "simulated_datasets", "inner_bootstrap_draws", "rejections_p_le005",
                           "rejection_fraction", "monte_carlo_se", "median_p")
  if (!all(calibration_columns %in% names(calibration))) {
    stop("Report requires complete null-calibration counts, rates and Monte Carlo errors.")
  }
  power <- read_result("diagnostics/aligned_effect_injection_summary.csv")
  # Normalize earlier contract names inside this reader, never rewrite outputs.
  power_aliases <- c(lambda = "effect_SD_multiplier", n_datasets = "simulated_datasets",
                     threshold = "fixed_null_abs95",
                     detected_fraction = "fraction_detected_two_sided",
                     median_abs_injected_log2FC = "median_abs_injected_log2_FC")
  for (name in names(power_aliases)) {
    alternative <- power_aliases[[name]]
    if (!name %in% names(power) && alternative %in% names(power)) {
      power[[name]] <- power[[alternative]]
    }
  }
  power_columns <- c("reference_universe", "contrast", "reference_set", "lambda", "n_datasets",
                     "median_ULM", "mean_ULM", "threshold", "detected_fraction",
                     "median_abs_injected_log2FC")
  if (!all(power_columns %in% names(power))) {
    stop("Report requires complete aligned-effect injection summaries.")
  }
  primary <- summary[summary$primary_comparison %in% TRUE, , drop = FALSE]
  desired_keys <- paste(rep(c("full", "fixed828"), each = 3L),
                       rep(paste0(c("pooled", "female", "male"), "_KO_vs_WT"), 2L),
                       rep(c("pooled", "female", "male"), 2L), sep = "__")
  key <- function(x) paste(x$reference_universe, x$contrast, x$reference_set, sep = "__")
  primary <- primary[match(desired_keys, key(primary)), , drop = FALSE]
  calibration <- calibration[match(desired_keys, key(calibration)), , drop = FALSE]
  if (anyNA(primary$contrast) || anyNA(calibration$contrast) ||
      nrow(calibration) != 6L || any(!is.finite(calibration$rejection_fraction)) ||
      any(!is.finite(calibration$monte_carlo_se)) ||
      any(calibration$simulated_datasets < 1L) ||
      any(calibration$inner_bootstrap_draws < 1L) ||
      any(calibration$rejection_fraction < 0 | calibration$rejection_fraction > 1) ||
      any(abs(calibration$rejection_fraction - calibration$rejections_p_le005 /
                calibration$simulated_datasets) > 1e-10)) {
    stop("Report calibration must match all six primary comparison keys and its integer counts.")
  }
  labels <- c(pooled_KO_vs_WT = "Pooled KO − WT / pooled FCG",
              female_KO_vs_WT = "Female KO − WT / ovarian FCG",
              male_KO_vs_WT = "Male KO − WT / testicular FCG")
  fmt <- function(x, digits = 3L) formatC(as.numeric(x), digits = digits, format = "f")
  result_table <- function(universe) {
    selected <- primary[primary$reference_universe == universe, , drop = FALSE]
    vapply(seq_len(nrow(selected)), function(i) {
      x <- selected[i, , drop = FALSE]
      paste0("| ", labels[[x$contrast]], " | ", x$n_matched_lipids, " | ",
             fmt(x$ULM_score), " | ", wp1_ulm_format_p(x$ULM_p_value), " | ",
             wp1_ulm_format_p(x$bootstrap_p_KO_only), " | ",
             wp1_ulm_format_p(x$bootstrap_p_KO_and_reference), " |")
    }, character(1L))
  }
  full <- primary[primary$reference_universe == "full", , drop = FALSE]
  historical <- primary[primary$reference_universe == "fixed828", , drop = FALSE]
  if (length(unique(full$n_matched_lipids)) != 1L ||
      length(unique(historical$n_matched_lipids)) != 1L) {
    stop("All contrasts must retain their universe's complete shared species roster.")
  }
  dependence <- null_diagnostics[null_diagnostics$reference_universe == "full" &
                                 null_diagnostics$primary_comparison %in% TRUE &
                                 null_diagnostics$method == "contrast_null_wild", , drop = FALSE]
  dependence <- dependence[match(key(full), key(dependence)), , drop = FALSE]
  if (nrow(dependence) != 3L || anyNA(dependence$contrast) ||
      any(!is.finite(dependence$ULM_null_SD)) || any(dependence$ULM_null_SD <= 0) ||
      any(!is.finite(dependence$nominal_ULM_p_le005_fraction)) ||
      any(dependence$nominal_ULM_p_le005_fraction < 0 |
            dependence$nominal_ULM_p_le005_fraction > 1) ||
      any(dependence$draws != provenance$bootstraps)) {
    stop("The three full-network dependence diagnostic rows must match completed null-draw counts.")
  }
  dependence_rows <- vapply(seq_len(nrow(dependence)), function(i) {
    x <- dependence[i, , drop = FALSE]
    paste0("| ", sub(" KO − WT.*$", "", labels[[x$contrast]]), " | ",
           x$draws, " | ", fmt(x$ULM_null_SD), " | ",
           fmt(100 * x$nominal_ULM_p_le005_fraction, 1L), "% |")
  }, character(1L))
  agreement_statement <- paste0(
    "The KO-only bootstrap p-values are ",
    paste(vapply(full$bootstrap_p_KO_only, wp1_ulm_format_p, character(1L)), collapse = "/"),
    ", compared with genotype-permutation p-values ",
    paste(vapply(full$permutation_p, wp1_ulm_format_p, character(1L)), collapse = "/"),
    " (pooled/female/male). Adding paired FCG uncertainty gives ",
    paste(vapply(full$bootstrap_p_KO_and_reference, wp1_ulm_format_p, character(1L)), collapse = "/"),
    ". Their agreement is evidence that the discrepancy with nominal ULM p-values is not caused solely by adding reference-estimation uncertainty."
  )
  female <- full[full$contrast == "female_KO_vs_WT", , drop = FALSE]
  male <- full[full$contrast == "male_KO_vs_WT", , drop = FALSE]
  male_statement <- paste0("The regenerated testicular-network comparison gives male ULM = ",
                           fmt(male$ULM_score), " and nominal ULM p = ",
                           wp1_ulm_format_p(male$ULM_p_value),
                           "; its KO-only bootstrap p is ",
                           wp1_ulm_format_p(male$bootstrap_p_KO_only),
                           " and KO-plus-reference bootstrap p is ",
                           wp1_ulm_format_p(male$bootstrap_p_KO_and_reference), ". ")
  if (male$ULM_score < 0 && male$ULM_p_value < 0.05 && male$bootstrap_p_KO_only >= 0.05) {
    male_statement <- paste0(male_statement,
                             "The male association has the predicted direction and a small nominal p-value, but the whole-animal test does not establish robust evidence under its contrast-null model.")
  } else if (male$ULM_score < 0 && male$bootstrap_p_KO_only < 0.05) {
    male_statement <- paste0(male_statement,
                             "The male association has the predicted direction and passes the KO-only whole-animal test under its stated contrast-null model; reference-estimation uncertainty is reported separately.")
  } else {
    male_statement <- paste0(male_statement,
                             "Direction and animal-level evidence should be read together under the stated contrast-null model.")
  }
  female_statement <- if (female$ULM_score > 0) {
    "The regenerated ovarian-network comparison has a positive ULM score, so its observed direction does not support the predicted inverse response to reduced Cyp2c38."
  } else {
    "The regenerated ovarian-network comparison has the predicted negative ULM direction; its strength is quantified by the p-values in the table."
  }
  exact_calibration_p <- vapply(seq_len(nrow(calibration)), function(i) {
    stats::binom.test(calibration$rejections_p_le005[i],
                      calibration$simulated_datasets[i], p = 0.05,
                      alternative = "two.sided")$p.value
  }, numeric(1L))
  calibration_rows <- vapply(seq_len(nrow(calibration)), function(i) {
    x <- calibration[i, , drop = FALSE]
    paste0("| ", if (x$reference_universe == "full") "Full" else "Historical", " / ",
           sub(" KO − WT.*$", "", labels[[x$contrast]]), " | ",
           x$rejections_p_le005, " / ", x$simulated_datasets, " | ",
           fmt(100 * x$rejection_fraction, 1L), "% | ",
           fmt(100 * x$monte_carlo_se, 1L), " | ",
           fmt(x$median_p), " | ", wp1_ulm_format_p(exact_calibration_p[i]), " |")
  }, character(1L))
  below <- calibration$rejection_fraction < 0.05 & exact_calibration_p < 0.05
  above <- calibration$rejection_fraction > 0.05 & exact_calibration_p < 0.05
  calibration_statement <- if (any(below)) {
    paste0("Some simulated comparisons reject below the nominal 5% rate and have an exact binomial diagnostic p<0.05 against that target. This is a fitted-scenario indication of conservatism for those comparisons, so a blanket claim that the procedure is not conservative would be too strong. ",
           if (any(above)) "Other comparisons reject above the target rate; calibration is not uniform across contrasts. " else "",
           "These are unadjusted diagnostic checks across six comparisons; their interpretation remains conditional on the fitted sign-symmetric residual scenario.")
  } else if (any(above)) {
    "These fitted-scenario simulations do not demonstrate excessive conservatism by the exact binomial check. Some comparisons reject above the 5% target with diagnostic p<0.05. These unadjusted checks across six comparisons do not establish universal calibration, particularly with a small number of simulated datasets."
  } else {
    "None of the six exact binomial diagnostic p-values is below 0.05 against the 5% target. This fitted-scenario check does not demonstrate excessive conservatism; with the limited number of simulated datasets, it also does not establish exact or universal calibration. The binomial checks are unadjusted diagnostic comparisons, conditional on the fitted sign-symmetric residual scenario."
  }
  full_power <- power[power$reference_universe == "full", , drop = FALSE]
  strengths <- c(0, 0.1, 0.2, 0.4, 0.8)
  contrast_order <- paste0(c("pooled", "female", "male"), "_KO_vs_WT")
  power_keys <- paste(full_power$contrast, full_power$lambda, sep = "__")
  expected_power_keys <- as.vector(outer(contrast_order, strengths,
                                        function(x, y) paste(x, y, sep = "__")))
  if (nrow(full_power) != 15L || anyDuplicated(power_keys) ||
      !setequal(power_keys, expected_power_keys) ||
      any(!is.finite(full_power$detected_fraction)) ||
      any(full_power$detected_fraction < 0 | full_power$detected_fraction > 1) ||
      any(!is.finite(full_power$median_abs_injected_log2FC)) ||
      any(!is.finite(full_power$n_datasets)) || any(full_power$n_datasets < 1L)) {
    stop("Report requires all five controlled strengths in each full-network contrast.")
  }
  strongest <- power[power$lambda == max(power$lambda) & power$reference_universe == "full", , drop = FALSE]
  strongest <- strongest[match(c("pooled", "female", "male"), strongest$reference_set), , drop = FALSE]
  if (nrow(strongest) != 3L || anyNA(strongest$contrast) ||
      any(!is.finite(strongest$detected_fraction)) ||
      any(strongest$detected_fraction < 0 | strongest$detected_fraction > 1) ||
      any(!is.finite(strongest$median_abs_injected_log2FC))) {
    stop("Report requires all three complete full-network aligned-effect scenarios.")
  }
  power_rows <- vapply(contrast_order, function(contrast) {
    x <- full_power[full_power$contrast == contrast, , drop = FALSE]
    x <- x[match(strengths, x$lambda), , drop = FALSE]
    if (length(unique(x$n_datasets)) != 1L) {
      stop("The controlled-strength scenarios must use the same replicate count.")
    }
    rates <- paste0(fmt(100 * x$detected_fraction, 1L), "%")
    paste0("| ", sub(" KO − WT.*$", "", labels[[contrast]]), " | ",
           unique(x$n_datasets), " | ", paste(rates, collapse = " | "), " |")
  }, character(1L))
  physical_effect_statement <- paste0(
    "At λ=0.8, the median absolute injected log2 effects are ",
    paste0(c("pooled ", "female ", "male "),
           fmt(strongest$median_abs_injected_log2FC), collapse = "; "),
    ". Every setting's physical effect size is retained in the linked complete injection table."
  )
  sensitivity_rows <- vapply(seq_len(nrow(primary)), function(i) {
    x <- primary[i, , drop = FALSE]
    paste0("| ", if (x$reference_universe == "full") "Full" else "Historical", " / ",
           sub(" KO − WT.*$", "", labels[[x$contrast]]), " | ",
           wp1_ulm_format_p(x$bootstrap_p_restricted_residuals), " | ",
           wp1_ulm_format_p(x$permutation_p), " |")
  }, character(1L))
  relative_path <- function(target) {
    source_parts <- strsplit(out, "/", fixed = TRUE)[[1L]]
    target_parts <- strsplit(normalizePath(target, mustWork = TRUE), "/", fixed = TRUE)[[1L]]
    common <- 0L
    while (common < min(length(source_parts), length(target_parts)) &&
           source_parts[common + 1L] == target_parts[common + 1L]) common <- common + 1L
    suffix <- if (common < length(target_parts)) target_parts[seq.int(common + 1L, length(target_parts))] else character()
    paste(c(rep("..", length(source_parts) - common), suffix), collapse = "/")
  }
  network_links <- vapply(seq_len(nrow(network_manifest)), function(i) {
    x <- network_manifest[i, , drop = FALSE]
    target <- if (file.exists(x$file)) x$file else file.path(out, x$file)
    if (!file.exists(target)) stop("Missing exported network: ", x$file)
    paste0("- [", x$gonadal_type, ", ", x$universe, " network](", relative_path(target),
           "): ", x$n_edges, " molecular edges, ", x$n_mice, " paired FCG mice.")
  }, character(1L))
  rerun <- paste0("Rscript scripts/workflows/07c_cyp2c38_ulm_resampling.R \\\n",
                  "  --bootstraps=", provenance$bootstraps, " --permutations=", provenance$permutations,
                  " --seed=", provenance$seed, " --workers=", provenance$workers,
                  " \\\n  --output=results/validation/cyp2c38_in_vivo/ulm_followup_replay")
  lines <- c(
    "# Cyp2c38 in vivo lipidomics: ULM scores, bootstrap tests and regenerated gonadal networks",
    "", paste0("Completed run: ", provenance$completed_at, ". Every result table below comes from the saved outputs."),
    "", male_statement, "", female_statement,
    "", "## Regenerated full-network results", "",
    "The native ULM score is the slope t-statistic from a linear model with an intercept: KO-minus-WT moderated lipid t-statistics are regressed on original Cyp2c38–lipid Pearson coefficients. A **negative ULM score** supports the expected perturbation direction. ULM p-values are nominal target-level p-values; bootstrap p-values resample whole mice and retain dependence among lipid species. All reported tests are two-sided.",
    "", "| KO contrast / FCG reference | Matched lipids | ULM score | ULM p-value | Bootstrap p, KO only | Bootstrap p, KO + reference |",
    "| --- | ---: | ---: | ---: | ---: | ---: |", result_table("full"),
    "", paste0("All **1,254 original lipidomics features** were regenerated without MOFA selection or a correlation/significance threshold. There are **", unique(full$n_matched_lipids),
                " unique strict molecular matches** to the retained KO feature roster. The ovarian and testicular networks use 21 and 22 original paired FCG mice, respectively; pooled coefficients use all 43 mice. The 32-mouse KO cohort contains female WT=8, female KO=8, male WT=7 and male KO=9. The primary model retains 1,237 curated endogenous KO features."),
    "", "## Historical feature-universe sensitivity", "",
    "| KO contrast / FCG reference | Matched lipids | ULM score | ULM p-value | Bootstrap p, KO only | Bootstrap p, KO + reference |",
    "| --- | ---: | ---: | ---: | ---: | ---: |", result_table("fixed828"),
    "", paste0("The historical 828-edge MOFA-derived universe retains **", unique(historical$n_matched_lipids),
                " shared species**. Its pooled and within-gonad coefficients numerically reproduce the original WP1 reference. Reporting this subset alongside the complete regenerated network makes the expanded feature universe explicit rather than attributing changed ULM scores to a different KO model."),
    "", "## What the bootstrap p-values test", "",
    "The earlier ordinary bootstrap resampled the observed cohorts to estimate uncertainty. Its draws were not a null distribution; the earlier reported animal-level p-values came from genotype permutations. This follow-up creates an explicit null distribution for the **ULM statistic**, with each tested KO-minus-WT mean contrast constrained to zero for every lipid while the other group-mean dimensions remain unrestricted and are estimated under that constraint. This is stronger than a null of zero ULM alignment with arbitrary nonaligned KO effects still present; the bootstrap p-values must be interpreted under this stated contrast-mean null.",
    "", "For the main KO-only wild bootstrap, fit the weighted four-group model, impose the tested mean-contrast restriction, and form null means. Use unrestricted four-group residuals divided by the square root of one minus their leverage (HC2). One independent Rademacher multiplier, −1 or +1 with equal probability, multiplies each whole-mouse residual lipid vector. This retains covariance among lipids and allows differing residual scales. Refit sex precision weights, robust/trend limma moderation, the KO t-statistics and ULM in every draw. This construction follows the null-imposed wild-bootstrap framework described by [Davidson and MacKinnon, Bootstrap Methods in Econometrics](https://www.russell-davidson.research.mcgill.ca/articles/bsurvey.pdf); finite-sample validity remains an assumption assessed below.",
    "", paste0("The reference network stays fixed in the KO-only column. The KO-plus-reference sensitivity additionally resamples paired original FCG RNA/lipid mice within XXO=10, XYO=11, XXT=11 and XYT=11 and recomputes each reference network. It therefore adds network-estimation uncertainty and answers a broader question. The run uses **", provenance$bootstraps,
                " null bootstrap draws**, **", provenance$permutations, " within-sex genotype permutations**, seed **", provenance$seed,
                "**, and ", provenance$workers, " workers. Empirical p-values use `(1 + number of null |ULM| values at least as large as observed)/(1 + valid draws)`. Invalid-draw counts and reproducible draw seeds are retained in the full outputs."),
    "", "A null-restricted-residual bootstrap is a sensitivity because those residuals can absorb a real KO signal and change the null distribution. The original within-sex genotype permutation is also retained as an independent comparison; it assumes exchangeability of genotype labels within sex. Neither bootstrap nor permutation treats molecular edges as independent experimental animals.",
    "", "| Universe / contrast | Bootstrap p, restricted residuals | Genotype permutation p |",
    "| --- | ---: | ---: |", sensitivity_rows,
    "", "## Why nominal ULM and whole-mouse p-values differ", "",
    "The fixed fitted-null bank gives a quantitative dependence diagnostic. Under the tested contrast's zero-mean-effect scenario, the whole-mouse ULM distribution is wide and nominal independent-target ULM p-values are frequently small:",
    "", "| Full-network contrast | Fitted-null draws | Null ULM SD | Fraction with nominal ULM p≤0.05 |",
    "| --- | ---: | ---: | ---: |", dependence_rows,
    "", "These fractions describe draws from this fitted contrast-null bank; they are not empirical population-wide error rates or the nested end-to-end calibration results below. They show that many nominal ULM p-values can be small even under the imposed mean-null when shared whole-mouse lipid variation is retained.",
    "", agreement_statement,
    "", "## Is the bootstrap overly conservative because the KO contrast is underpowered?", "",
    "Low power and conservatism are different properties. Low power is a low probability of detecting a real effect; conservatism is a rejection probability below the nominated level when the null is true. A weak individual-lipid contrast can still contribute collectively to ULM, and its lack of significance does not automatically make the bootstrap conservative. Conversely, a large nominal ULM magnitude can arise from correlated lipids sharing noise from a few mice. Whole-mouse resampling includes that dependence.",
    "", "The null-calibration audit generates globally null datasets with no female or male KO effects, using the original HC2-corrected whole-mouse residual vectors and independent Rademacher signs. It preserves the actual group sizes and estimated sex-dependent lipid covariance. Each dataset constructs its own contrast-specific null, refits the complete pipeline, and receives a fresh inner wild-bootstrap test. The nominal rejection threshold is p≤0.05. This is a calibration check within the fitted sign-symmetric residual scenario, using the same resampling family as the test. The results below report counts, rejection rates and Monte Carlo standard errors rather than assuming ideal calibration.",
    "", "| Universe / contrast | Null rejections / datasets | Null rejection rate | Monte Carlo SE, percentage points | Median bootstrap p | Exact binomial p vs 5% |",
    "| --- | ---: | ---: | ---: | ---: | ---: |", calibration_rows,
    "", calibration_statement,
    "", paste0("With only ", unique(calibration$simulated_datasets),
                " outer null datasets per comparison, moderate conservatism cannot be ruled out regardless of the exact binomial diagnostic. Non-rejection of the 5% target is not evidence that the true rejection rate equals exactly 5%."),
    "", paste0("Each null dataset uses ", unique(calibration$inner_bootstrap_draws),
                " inner bootstrap draws. These simulations condition on the fitted residual scenario and observed feature selection. They do not rule out miscalibration from unrecorded batch/litter structure, asymmetric errors or an inaccurate fitted covariance model."),
    "", "Known-effect injection curves add deliberately concordant, network-shaped KO effects at controlled strengths while retaining the actual group sizes and fitted noise. Detection fractions across the entire full-network curve are:",
    "", "| Full-network contrast | Datasets per setting | λ=0 | λ=0.1 | λ=0.2 | λ=0.4 | λ=0.8 |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: |", power_rows,
    "", physical_effect_statement,
    "", "Injection strength is a controlled scenario parameter: shifts oppose the centered, SD-standardized full-network correlation pattern and are scaled by the original within-sex residual SD. Detection compares each simulated absolute ULM score with the observed-data null bank's 95th percentile, held fixed; it does not run a fresh inner test for every injected dataset. These are fitted-scenario sensitivity calculations, not an estimate of post-hoc power for the unknown true biological effect. Poor detection of small injected effects supports a sample-size/noise limitation; it does not by itself diagnose an overly conservative test. These diagnostics do not justify interpreting a large bootstrap p-value as proof that Cyp2c38 has no lipid effect.",
    "", "## Pairing, feature identity and interpretation", "",
    "Original scaled hepatocyte Cyp2c38 expression is paired with log2 original lipid intensities; reference mice retain their recorded gonadal and chromosome groups. The separate networks describe ovarian or testicular FCG associations, while pooled associations can also contain between-gonad covariance. There is no new KO RNA measurement or MOFA refit. Molecular matching only standardizes the last pipe-delimited annotation and class-parenthesis formatting: chain order, slash versus underscore, ether prefixes, oxygenation and structural resolution remain intact. Internal standards, isotope-labelled features and repository-excluded lipid classes retain the WP1 exclusions. No animal or lipid is removed because its outcome is unfavorable.",
    "", "KO models use log2 supplied positive intensities, four group means, sex-specific precision weights and robust intensity-trend empirical Bayes. The pooled effect is half the female plus male KO-minus-WT effects; t-statistics are computed from that contrast rather than averaged. Upstream normalization/gap-filling, litter and acquisition-batch provenance remain unavailable. The perturbation is heterozygous in vivo Cyp2c38 reduction; the earlier primary-hepatocyte siRNA experiment is a distinct intervention. Correlation and concordance do not establish direct substrates, flux or the consequences of complete deletion.",
    "", "## Saved artifacts", "",
    "- [Complete ULM and p-value summary](ULM_summary.csv), [network export manifest](network_manifest.csv), [ULM concordance figures](plots/ULM_concordance.pdf), [run provenance](provenance.json).",
    "- [Fitted-null dependence diagnostics](null_distribution_diagnostics.csv).",
    "- [Null-calibration summary](diagnostics/calibration_null_rejection_summary.csv) and [complete aligned-effect injection curves](diagnostics/aligned_effect_injection_summary.csv).",
    network_links,
    "", "## Reproduce", "",
    "Run from the repository root, preserving the completed output with a fresh destination:",
    "", "```sh", rerun, "```", "",
    "The saved provenance, input hashes, R session details and all resampling draws identify the numerical run. Keep the same null, contrast orientation, molecular matching and feature universes when comparing replays."
  )
  text <- paste(lines, collapse = "\n")
  matches <- regmatches(text, gregexpr("\\]\\(([^)]+)\\)", text))[[1L]]
  targets <- sub("^\\]\\((.*)\\)$", "\\1", matches)
  local_targets <- targets[!grepl("^https?://", targets)]
  missing <- local_targets[!file.exists(file.path(out, local_targets))]
  if (length(missing)) stop("Missing linked report artifact(s): ", paste(missing, collapse = ", "))
  destination <- file.path(out, "README.md")
  writeLines(lines, destination, useBytes = TRUE)
  cat("Wrote ULM-focused report: ", destination, "\n", sep = "")
  invisible(destination)
}
