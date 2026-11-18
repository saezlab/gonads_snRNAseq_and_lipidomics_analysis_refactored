# Shared statistical helpers for Workflow 07c.
# Sourcing defines functions without fitting models or changing the RNG state.

#' ULM slope t-statistic and its nominal two-sided p-value.
#'
#' This is the intercept-containing OLS model effect ~ prediction used by
#' decoupleR::run_ulm. The score retains its original sign: a negative score
#' supports concordance of gene-lipid predictions with a KO-minus-WT effect.
#' A whole invalid vector remains invalid; missing/nonfinite edges are never
#' silently dropped. Constant or numerically near-constant inputs return NA.
wp1_ulm_stats <- function(prediction, effect) {
  if (!is.numeric(prediction) || !is.numeric(effect) ||
      !is.null(dim(prediction)) || !is.null(dim(effect)) ||
      length(prediction) != length(effect)) {
    stop("prediction and effect must be equal-length numeric vectors.")
  }
  n <- length(prediction)
  unavailable <- c(score = NA_real_, p_value = NA_real_, n_lipids = n,
                   slope = NA_real_)
  if (n < 3L || any(!is.finite(prediction)) || any(!is.finite(effect))) {
    return(unavailable)
  }
  px <- max(abs(prediction))
  py <- max(abs(effect))
  if (px == 0 || py == 0) return(unavailable)
  x <- prediction / px
  y <- effect / py
  x <- x - mean(x)
  y <- y - mean(y)
  xx <- sum(x * x)
  yy <- sum(y * y)
  # Relative tolerance avoids classifying an otherwise variable small-scale
  # predictor as constant merely because all its values are close to zero.
  threshold <- n * .Machine$double.eps
  if (xx <= threshold || yy <= threshold) return(unavailable)
  xy <- sum(x * y)
  correlation <- max(-1, min(1, xy / sqrt(xx * yy)))
  score <- correlation * sqrt((n - 2) / (1 - correlation^2))
  c(score = score,
    p_value = 2 * stats::pt(abs(score), df = n - 2, lower.tail = FALSE),
    n_lipids = n, slope = (xy / xx) * (py / px))
}

#' Prepare a contrast-null whole-mouse wild bootstrap on the supplied scale.
#'
#' @param y Numeric finite feature-by-mouse matrix, in metadata column order.
#' @param meta Ordered metadata with sample_id/sample, sex and group.
#' @param restricted FALSE uses full four-group residuals with HC2 leverage;
#'   TRUE uses residuals and leverage from the null-restricted model.
#' @param contrast NULL imposes zero female AND male KO effects (global null).
#'   Otherwise a named wp1_contrasts() column or a named four-group vector.
#'   A tested-contrast null permits effects in the remaining nuisance contrasts.
#' @param fit A wp1_fit() result supplying constant-within-group precision
#'   weights. If a contrast is supplied and fit is NULL, fit the original model.
#'
#' The null restricts the chosen genotype mean contrast to zero FOR EVERY
#' lipid. This is stronger than a zero ULM slope while arbitrary nonaligned KO
#' effects remain. It calibrates the observed ULM statistic under that explicit
#' mean-null. Reference coefficients are handled separately by the caller.
#' HC2 divides each residual mouse vector by sqrt(1 - its model leverage).
wp1_null_components <- function(y, meta, restricted = FALSE,
                                contrast = NULL, fit = NULL) {
  if (!is.matrix(y) || !is.numeric(y) || any(!is.finite(y)) ||
      !is.data.frame(meta) || nrow(meta) != ncol(y) ||
      !all(c("sex", "group") %in% names(meta)) ||
      length(restricted) != 1L || is.na(restricted) ||
      !is.logical(restricted)) {
    stop("Supply a finite feature-by-mouse matrix, ordered metadata and logical restricted.")
  }
  sample_fields <- intersect(c("sample_id", "sample"), names(meta))
  for (field in sample_fields) {
    if (!identical(as.character(meta[[field]]), colnames(y))) {
      stop("Metadata ", field, " must match y columns in the same order.")
    }
  }
  group_levels <- c("F_WT", "F_KO", "M_WT", "M_KO")
  group <- factor(as.character(meta$group), levels = group_levels)
  if (anyNA(group)) stop("Unknown or missing genotype/sex group.")
  counts <- tabulate(as.integer(group), nbins = 4L)
  names(counts) <- group_levels
  if (any(counts < 2L)) stop("Every genotype/sex group requires at least two mice.")
  sex <- toupper(substr(as.character(meta$sex), 1L, 1L))
  if (anyNA(sex) || !identical(sex, substr(as.character(group), 1L, 1L))) {
    stop("Metadata sex must agree with its genotype/sex group.")
  }
  coefficients <- vapply(group_levels, function(g) {
    rowMeans(y[, group == g, drop = FALSE])
  }, numeric(nrow(y)))
  if (nrow(y) == 1L) coefficients <- matrix(coefficients, nrow = 1L)
  dimnames(coefficients) <- list(rownames(y), group_levels)
  full_mean <- coefficients[, as.integer(group), drop = FALSE]
  dimnames(full_mean) <- dimnames(y)
  weights <- rep(1, ncol(y))
  if (!is.null(contrast) && is.null(fit)) {
    if (!exists("wp1_fit", mode = "function")) {
      stop("Source cyp2c38_models.R before requesting a fitted contrast null.")
    }
    fit <- wp1_fit(y, meta)
  }
  if (!is.null(fit)) {
    if (!is.list(fit) || is.null(fit$weights) ||
        !is.numeric(fit$weights) || length(fit$weights) != ncol(y) ||
        any(!is.finite(fit$weights)) || any(fit$weights <= 0)) {
      stop("fit must be a wp1_fit() result with positive mouse precision weights.")
    }
    if (!is.null(names(fit$weights)) &&
        !identical(names(fit$weights), colnames(y))) {
      stop("Fitted precision weights must retain y column order.")
    }
    weights <- as.numeric(fit$weights)
  }
  group_weights <- vapply(group_levels, function(g) {
    w <- weights[group == g]
    if (diff(range(w)) > 1e-10 * max(w)) {
      stop("Contrast projection requires constant precision weights within each group.")
    }
    mean(w)
  }, numeric(1L))
  null_coefficients <- coefficients
  full_leverage <- 1 / counts[as.integer(group)]
  if (is.null(contrast)) {
    null_leverage <- rep(NA_real_, ncol(y))
    for (sx in c("F", "M")) {
      columns <- startsWith(group_levels, paste0(sx, "_"))
      selected <- sex == sx
      null_coefficients[, columns] <-
        rowSums(sweep(y[, selected, drop = FALSE], 2L, weights[selected], "*")) /
        sum(weights[selected])
      null_leverage[selected] <- weights[selected] / sum(weights[selected])
    }
    hypothesis <- "Zero female and male KO-minus-WT mean effects for every lipid"
    contrast_vector <- NULL
  } else {
    if (is.character(contrast)) {
      if (length(contrast) != 1L ||
          !exists("wp1_contrasts", mode = "function") ||
          !contrast %in% colnames(wp1_contrasts())) {
        stop("contrast must name one column of wp1_contrasts().")
      }
      contrast_vector <- wp1_contrasts()[, contrast]
      contrast_label <- contrast
    } else {
      if (!is.numeric(contrast) || length(contrast) != 4L ||
          any(!is.finite(contrast)) || all(contrast == 0) ||
          is.null(names(contrast)) ||
          !setequal(names(contrast), group_levels)) {
        stop("A numeric contrast requires four uniquely named finite group coefficients.")
      }
      contrast_vector <- contrast[group_levels]
      contrast_label <- "supplied four-group contrast"
    }
    v <- 1 / (counts * group_weights)
    vc <- v * contrast_vector
    denominator <- sum(contrast_vector * vc)
    effect <- as.numeric(coefficients %*% contrast_vector)
    null_coefficients <- coefficients - outer(effect, vc / denominator)
    null_leverage <- (1 / counts - group_weights * vc^2 / denominator)[as.integer(group)]
    hypothesis <- paste("Zero", contrast_label, "mean effect for every lipid; other contrasts retained")
  }
  null_mean <- null_coefficients[, as.integer(group), drop = FALSE]
  dimnames(null_mean) <- dimnames(y)
  leverage <- if (restricted) null_leverage else full_leverage
  leverage <- as.numeric(leverage)
  if (any(!is.finite(leverage)) || any(leverage < -1e-12) || any(leverage >= 1)) {
    stop("Invalid residual leverage for HC2 correction.")
  }
  leverage <- pmax(leverage, 0)
  residuals <- y - if (restricted) null_mean else full_mean
  residuals_hc2 <- sweep(residuals, 2L, sqrt(1 - leverage), "/")
  list(null_mean = null_mean, residuals_hc2 = residuals_hc2,
       group_counts = counts, groups = counts, leverage = leverage,
       group_coefficients = coefficients, null_group_coefficients = null_coefficients,
       contrast = contrast_vector, precision_weights = weights,
       restricted = restricted, null_hypothesis = hypothesis)
}

#' Generate one null draw with a single shared sign per whole mouse vector.
#'
#' Independent Rademacher signs are drawn only if multipliers is omitted. A
#' caller controls seeds; this function never creates or resets a random seed.
wp1_wild_null <- function(components, multipliers = NULL) {
  if (!is.list(components) || !is.matrix(components$null_mean) ||
      !is.matrix(components$residuals_hc2) ||
      !identical(dim(components$null_mean), dim(components$residuals_hc2)) ||
      !identical(dimnames(components$null_mean), dimnames(components$residuals_hc2))) {
    stop("components must come from wp1_null_components().")
  }
  if (is.null(multipliers)) {
    multipliers <- sample(c(-1, 1), ncol(components$null_mean), replace = TRUE)
  }
  if (!is.numeric(multipliers) || length(multipliers) != ncol(components$null_mean) ||
      anyNA(multipliers) || any(!multipliers %in% c(-1, 1))) {
    stop("Supply exactly one Rademacher (-1 or +1) multiplier per mouse.")
  }
  components$null_mean + sweep(components$residuals_hc2, 2L, multipliers, "*")
}

#' Recompute pooled, ovarian and testicular references from paired mouse draws.
#'
#' indices selects whole paired RNA/lipid mice; duplicates are allowed for a
#' bootstrap. Original lipid row order is retained. Undefined correlations,
#' including a draw with constant Cyp2c38 expression, remain NA explicitly.
wp1_reference_recompute <- function(reference, indices) {
  if (!is.list(reference) || !is.numeric(reference$gene) ||
      !is.matrix(reference$lipids) || !is.numeric(reference$lipids) ||
      !is.data.frame(reference$metadata) ||
      !"gonadal_type" %in% names(reference$metadata) ||
      length(reference$gene) != ncol(reference$lipids) ||
      nrow(reference$metadata) != ncol(reference$lipids) ||
      any(!is.finite(reference$gene)) || any(!is.finite(reference$lipids)) ||
      !is.numeric(indices) || !length(indices) || any(!is.finite(indices)) ||
      any(indices != as.integer(indices)) ||
      any(indices < 1L | indices > length(reference$gene))) {
    stop("Supply a paired reference and valid integer mouse indices.")
  }
  if (!is.null(names(reference$gene)) && !is.null(colnames(reference$lipids)) &&
      !identical(names(reference$gene), colnames(reference$lipids))) {
    stop("Reference gene names must match lipid mouse columns.")
  }
  if ("sample_id" %in% names(reference$metadata) &&
      !identical(as.character(reference$metadata$sample_id), colnames(reference$lipids))) {
    stop("Reference metadata must match lipid mouse columns.")
  }
  gonad <- reference$metadata$gonadal_type[indices]
  if (anyNA(gonad) || !all(gonad %in% c("Ovaries", "Testes"))) {
    stop("Reference gonadal types must be Ovaries or Testes.")
  }
  gene <- reference$gene[indices]
  lipid <- reference$lipids[, indices, drop = FALSE]
  output <- data.frame(reference_lipid = rownames(reference$lipids),
                       stringsAsFactors = FALSE)
  selected <- list(pooled = seq_along(indices),
                   female = which(gonad == "Ovaries"),
                   male = which(gonad == "Testes"))
  for (sx in names(selected)) {
    ii <- selected[[sx]]
    values <- rep(NA_real_, nrow(lipid))
    if (length(ii) >= 3L && stats::sd(gene[ii]) > 0) {
      variable <- apply(lipid[, ii, drop = FALSE], 1L, stats::sd) > 0
      if (any(variable)) {
        values[variable] <- as.numeric(stats::cor(gene[ii],
                                  t(lipid[variable, ii, drop = FALSE])))
      }
    }
    output[[paste0("r_", sx)]] <- values
  }
  output
}
