# Conditional calibration and known-effect sensitivity for Workflow 07c.
# Source support_functions.R, cyp2c38_models.R and cyp2c38_ulm_bootstrap.R first.
# Sourcing this file defines one function; no models or simulations run at source.

#' Audit contrast-specific whole-mouse bootstrap at the actual WP1 sample sizes.
#'
#' Calibration data are generated under a global no-KO-effect null using the
#' original full-model HC2 residual vectors and independent mouse Rademacher
#' multipliers. Each synthetic dataset then constructs its own contrast-specific
#' constrained null, preserving its estimated nuisance contrasts, and refits
#' weights and moderation for every inner bootstrap draw. Both reference feature
#' universes share the same synthetic data and model fits.
#'
#' Known-effect injections are idealized sensitivity diagnostics: shifts align
#' with every matched full-network coefficient, and detection uses fixed critical
#' values from the observed-data null bank. They are neither prospective power
#' estimates nor estimates of power calculated from the observed KO effects.
#'
#' @param out Follow-up output directory containing null_wild_draws.csv.
#' @param y Finite log2 feature-by-mouse matrix, ordered as meta.
#' @param meta Ordered KO metadata accepted by wp1_fit().
#' @param reference Paired original FCG reference; checked for feature coverage.
#' @param map Strict full-network matches with feature_id, reference_lipid,
#'   fixed_reference, r_pooled, r_female and r_male columns.
#' @param specs Comparison table with reference_universe (full/fixed828),
#'   contrast, reference_set (pooled/female/male) and primary_comparison.
#' @param seed Integer RNG seed. Per-dataset seeds are saved before parallel work.
#' @param workers Number of independent simulated datasets to fit concurrently.
#' @param n_datasets Number of calibration datasets (default 120).
#' @param n_inner Null-bootstrap draws per calibration dataset (default 199).
#' @param n_injection Synthetic datasets per injection setting (default 399).
#' @return Invisibly a list of calibration/injection summary and draw tables.
wp1_calibration_audit <- function(out, y, meta, reference, map, specs, seed,
                                  workers, n_datasets = 120L, n_inner = 199L,
                                  n_injection = 399L) {
  required_functions <- c("wp1_fit", "wp1_null_components", "wp1_wild_null",
                          "wp1_ulm_stats")
  if (any(!vapply(required_functions, exists, logical(1L), mode = "function"))) {
    stop("Source the WP1 model and ULM bootstrap helpers before calibration.")
  }
  controls <- c(workers = workers, n_datasets = n_datasets,
                n_inner = n_inner, n_injection = n_injection)
  if (length(controls) != 4L || !is.numeric(controls) || anyNA(controls) ||
      any(controls != as.integer(controls)) || any(controls < 1L) ||
      n_inner < 19L || length(seed) != 1L || is.na(seed) ||
      !is.numeric(seed) || seed != as.integer(seed)) {
    stop("Supply positive integer controls, at least 19 inner draws and one integer seed.")
  }
  if (!is.character(out) || length(out) != 1L || is.na(out) ||
      !dir.exists(out) || !file.exists(file.path(out, "null_wild_draws.csv"))) {
    stop("out must contain the main null_wild_draws.csv before diagnostics run.")
  }
  needed_map <- c("feature_id", "reference_lipid", "fixed_reference",
                  "r_pooled", "r_female", "r_male")
  needed_specs <- c("reference_universe", "contrast", "reference_set",
                    "primary_comparison")
  if (!is.data.frame(map) || !all(needed_map %in% names(map)) ||
      !is.data.frame(specs) || !all(needed_specs %in% names(specs)) ||
      !is.logical(map$fixed_reference) || anyNA(map$fixed_reference) ||
      !is.logical(specs$primary_comparison) || anyNA(specs$primary_comparison) ||
      anyNA(map$feature_id) || anyDuplicated(map$feature_id) ||
      anyDuplicated(map$reference_lipid) ||
      !all(map$feature_id %in% rownames(y)) ||
      !is.list(reference) || !is.matrix(reference$lipids) ||
      !all(map$reference_lipid %in% rownames(reference$lipids))) {
    stop("Supply validated, unique strict matches and the paired full reference.")
  }
  cases <- specs[specs$primary_comparison, needed_specs, drop = FALSE]
  rownames(cases) <- NULL
  if (nrow(cases) != 6L ||
      !setequal(cases$reference_universe, c("full", "fixed828")) ||
      !setequal(cases$contrast, c("pooled_KO_vs_WT", "female_KO_vs_WT",
                                "male_KO_vs_WT")) ||
      anyDuplicated(cases[, c("reference_universe", "contrast")])) {
    stop("Expected three primary contrasts in each of the two feature universes.")
  }
  labels <- paste(cases$reference_universe, cases$contrast,
                  cases$reference_set, sep = "__")
  selected <- lapply(cases$reference_universe, function(u) {
    if (u == "full") seq_len(nrow(map)) else which(map$fixed_reference)
  })
  predictions <- lapply(seq_len(nrow(cases)), function(j) {
    r <- map[[paste0("r_", cases$reference_set[j])]][selected[[j]]]
    if (!is.numeric(r) || any(!is.finite(r)) || length(r) < 10L ||
        stats::sd(r) == 0) stop("Invalid reference predictor for ", labels[j])
    r
  })
  cases$n_lipids <- vapply(selected, length, integer(1L))
  null_bank <- utils::read.csv(file.path(out, "null_wild_draws.csv"),
                               check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(labels %in% names(null_bank))) {
    stop("Main null bank lacks the labeled primary comparisons.")
  }
  critical <- vapply(labels, function(lab) {
    z <- null_bank[[lab]]
    if (!is.numeric(z) || mean(is.finite(z)) < 0.95) {
      stop("Insufficient finite main null draws for ", lab)
    }
    unname(stats::quantile(abs(z), 0.95, na.rm = TRUE))
  }, numeric(1L))
  diagnostic_dir <- file.path(out, "diagnostics")
  dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, name) {
    utils::write.csv(x, file.path(diagnostic_dir, name), row.names = FALSE, na = "")
  }
  score_fit <- function(fit, indices = seq_len(nrow(cases))) {
    vapply(indices, function(j) {
      tv <- fit$fit$t[map$feature_id[selected[[j]]], cases$contrast[j]]
      unname(wp1_ulm_stats(predictions[[j]], tv)["score"])
    }, numeric(1L))
  }
  original_fit <- wp1_fit(y, meta)
  generating_null <- wp1_null_components(y, meta, restricted = FALSE,
                                         contrast = NULL, fit = original_fit)
  contrast_names <- unique(cases$contrast)
  set.seed(seed)
  seeds_calibration <- sample.int(.Machine$integer.max, n_datasets)
  seeds_injection <- sample.int(.Machine$integer.max, n_injection)
  run_chunks <- function(draw_seeds, fun, label) {
    output <- vector("list", length(draw_seeds))
    chunk_size <- max(1L, 2L * as.integer(workers))
    for (start in seq.int(1L, length(draw_seeds), by = chunk_size)) {
      ii <- start:min(start + chunk_size - 1L, length(draw_seeds))
      output[ii] <- parallel::mclapply(draw_seeds[ii], function(s) {
        set.seed(s)
        fun(s)
      }, mc.cores = as.integer(workers), mc.set.seed = FALSE)
      failed <- vapply(output[ii], inherits, logical(1L), what = "try-error")
      if (any(failed)) stop("A ", label, " simulation failed: ", output[ii][[which(failed)[1L]]])
      if (!startsWith(label, "Aligned injection") || max(ii) == length(draw_seeds)) {
        cat(label, max(ii), "/", length(draw_seeds), "\n")
      }
    }
    output
  }
  cat("Calibration: global no-KO-effect data; contrast-specific HC2 wild tests; ",
      n_datasets, " datasets x ", n_inner, " inner draws.\n", sep = "")
  calibration_one <- function(draw_seed) {
    synthetic <- wp1_wild_null(generating_null)
    fitted <- wp1_fit(synthetic, meta)
    observed <- score_fit(fitted)
    components <- setNames(lapply(contrast_names, function(ct) {
      wp1_null_components(synthetic, meta, restricted = FALSE,
                           contrast = ct, fit = fitted)
    }), contrast_names)
    inner <- matrix(NA_real_, nrow = n_inner, ncol = nrow(cases))
    for (b in seq_len(n_inner)) {
      signs <- sample(c(-1, 1), ncol(y), replace = TRUE)
      for (ct in contrast_names) {
        jj <- which(cases$contrast == ct)
        null_y <- wp1_wild_null(components[[ct]], signs)
        inner[b, jj] <- score_fit(wp1_fit(null_y, meta), jj)
      }
    }
    valid <- colSums(is.finite(inner))
    if (any(!is.finite(observed)) || any(valid < 0.95 * n_inner)) {
      stop("Insufficient finite calibration scores at seed ", draw_seed)
    }
    p <- vapply(seq_len(nrow(cases)), function(j) {
      (1 + sum(abs(inner[, j]) >= abs(observed[j]), na.rm = TRUE)) /
        (1 + valid[j])
    }, numeric(1L))
    cbind(cases, data.frame(replicate_seed = draw_seed,
                            observed_ULM = observed, bootstrap_p = p,
                            valid_null_draws = valid,
                            null_ULM_SD = apply(inner, 2L, stats::sd, na.rm = TRUE),
                            stringsAsFactors = FALSE))
  }
  cal <- do.call(rbind, run_chunks(seeds_calibration, calibration_one,
                                   "Conditional null calibration"))
  cal_summary <- do.call(rbind, lapply(seq_len(nrow(cases)), function(j) {
    ii <- cal$reference_universe == cases$reference_universe[j] &
      cal$contrast == cases$contrast[j]
    p <- cal$bootstrap_p[ii]
    rejected <- sum(p <= 0.05)
    fraction <- rejected / length(p)
    interval <- stats::binom.test(rejected, length(p))$conf.int
    cbind(cases[j, , drop = FALSE],
          data.frame(simulated_datasets = length(p),
                     inner_bootstrap_draws = n_inner,
                     rejections_p_le005 = rejected,
                     rejection_fraction = fraction,
                     monte_carlo_se = sqrt(fraction * (1 - fraction) / length(p)),
                     binomial_low = interval[1L], binomial_high = interval[2L],
                     median_p = stats::median(p), stringsAsFactors = FALSE))
  }))
  write_csv(cal, "calibration_null_dataset_pvalues.csv")
  write_csv(cal_summary, "calibration_null_rejection_summary.csv")

  # The DGP HC2 residuals estimate the original within-sex feature dispersion.
  # Inject one known reference-aligned pattern into full-network matched species,
  # then score both universes from each common refitted KO contrast.
  sex_sd <- vapply(c("F", "M"), function(sx) {
    sqrt(rowMeans(generating_null$residuals_hc2[, meta$sex == sx, drop = FALSE]^2))
  }, numeric(nrow(y)))
  rownames(sex_sd) <- rownames(y)
  lambda <- c(0, 0.1, 0.2, 0.4, 0.8)
  injection_output <- vector("list", length(contrast_names) * length(lambda))
  setting <- 0L
  for (ct in contrast_names) {
    jj <- which(cases$contrast == ct)
    full_case <- jj[cases$reference_universe[jj] == "full"]
    if (length(full_case) != 1L) stop("Each injection needs one full-network case.")
    r <- predictions[[full_case]]
    predictor <- (r - mean(r)) / stats::sd(r)
    sexes <- if (ct == "pooled_KO_vs_WT") c("F", "M") else {
      if (ct == "female_KO_vs_WT") "F" else "M"
    }
    for (strength in lambda) {
      effect <- matrix(0, nrow = nrow(y), ncol = ncol(y), dimnames = dimnames(y))
      shifted_columns <- which(meta$genotype == "KO" & meta$sex %in% sexes)
      for (sx in sexes) {
        kk <- which(meta$sex == sx & meta$genotype == "KO")
        delta <- -strength * predictor * sex_sd[map$feature_id, sx]
        effect[map$feature_id, kk] <- delta
      }
      injection_one <- function(draw_seed) {
        synthetic <- wp1_wild_null(generating_null) + effect
        score_fit(wp1_fit(synthetic, meta), jj)
      }
      scores <- do.call(rbind, run_chunks(seeds_injection, injection_one,
                         paste("Aligned injection", ct, "lambda", strength)))
      if (any(colMeans(is.finite(scores)) < 0.95)) {
        stop("Insufficient finite injection scores for ", ct, " at ", strength)
      }
      setting <- setting + 1L
      injection_output[[setting]] <- do.call(rbind, lapply(seq_along(jj), function(k) {
        j <- jj[k]
        z <- scores[, k]
        keep <- is.finite(z)
        detected <- sum(abs(z[keep]) >= critical[j])
        fraction <- detected / sum(keep)
        ii <- map$feature_id[selected[[j]]]
        cbind(cases[j, , drop = FALSE],
          data.frame(effect_SD_multiplier = strength,
                     simulated_datasets = sum(keep),
                     mean_ULM = mean(z[keep]), median_ULM = stats::median(z[keep]),
                     fixed_null_abs95 = unname(critical[j]),
                     detected_two_sided = detected,
                     fraction_detected_two_sided = fraction,
                     monte_carlo_se = sqrt(fraction * (1 - fraction) / sum(keep)),
                     median_abs_injected_log2_FC =
                       stats::median(abs(effect[ii, shifted_columns, drop = FALSE])),
                     injected_feature_universe = "full", stringsAsFactors = FALSE))
      }))
      draw_table <- data.frame(replicate_seed = seeds_injection, scores,
                               check.names = FALSE)
      names(draw_table)[-1L] <- labels[jj]
      filename <- paste0("injection_draws_", ct, "_lambda_",
                          gsub("\\.", "p", as.character(strength)), ".csv")
      write_csv(draw_table, filename)
    }
  }
  injection_summary <- do.call(rbind, injection_output)
  rownames(cal_summary) <- rownames(injection_summary) <- NULL
  write_csv(injection_summary, "aligned_effect_injection_summary.csv")
  write_csv(data.frame(kind = "calibration", replicate = seq_along(seeds_calibration),
                       seed = seeds_calibration), "calibration_seeds.csv")
  write_csv(data.frame(kind = "injection", replicate = seq_along(seeds_injection),
                       seed = seeds_injection), "injection_seeds.csv")
  provenance <- list(
    seed = seed, workers = workers, calibration_datasets = n_datasets,
    inner_bootstrap_draws = n_inner, injection_datasets_per_setting = n_injection,
    lambda = lambda, matched_lipids_full = nrow(map),
    matched_lipids_fixed = sum(map$fixed_reference),
    data_generating_null = "No KO mean effects in either sex; original sex means and full four-group HC2 mouse residual covariance",
    tested_null = "Each tested contrast mean is zero for every lipid; estimated nuisance contrasts retained using the production weighted projection",
    shared_multipliers = "One independent Rademacher sign per whole mouse; shared across all lipid species and, within an inner draw, across three tested-contrast simulations",
    fitted_model = "Production sex-weighted limma; weights, trend and robust moderation refitted for every simulated dataset",
    references = "Both full and historical networks fixed; original feature selection conditioned upon",
    injection = "Known negative mean effects proportional to centered SD-standardized full-network r times original within-sex residual SD; unmatched features unshifted",
    injection_detection = "Fixed observed-data null bank 95th absolute ULM percentile; not a refitted prospective power calculation",
    conditioning = "Estimated residual covariance, independent/sign-symmetric mouse errors, actual sample counts, measured feature universe; unknown batch/litter effects unmodeled",
    interpretation = "Finite Monte Carlo calibration under a fitted scenario; does not establish universal calibration or power for the unknown true biological effect"
  )
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to save calibration provenance.")
  }
  jsonlite::write_json(provenance,
                        file.path(diagnostic_dir, "calibration_provenance.json"),
                        auto_unbox = TRUE, pretty = TRUE)
  cat("Conditional ULM calibration and known-effect sensitivity complete.\n")
  invisible(list(calibration = cal_summary, calibration_draws = cal,
                  injection = injection_summary, provenance = provenance))
}
