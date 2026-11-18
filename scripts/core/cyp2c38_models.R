# Statistical models for Workflow 07b: in vivo Cyp2c38 lipidomics.
# Input y is a finite numeric feature-by-mouse matrix on the log2 scale.
# Sample metadata must already be ordered to match the matrix columns.
# The pooled effect gives female and male KO-minus-WT effects equal weight.

wp1_group_levels <- function() {
  c("F_WT", "F_KO", "M_WT", "M_KO")
}

wp1_design <- function(meta) {
  if (!is.data.frame(meta) || !"group" %in% names(meta)) {
    stop("Metadata must be a data frame with a group column.")
  }
  levels <- wp1_group_levels()
  group <- factor(as.character(meta$group), levels = levels)
  if (anyNA(group)) {
    stop("Metadata group values must be F_WT, F_KO, M_WT or M_KO.")
  }
  if (any(tabulate(as.integer(group), nbins = length(levels)) < 2L)) {
    stop("Each sex/genotype group needs at least two biological mice.")
  }
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels
  if (qr(design)$rank != ncol(design)) {
    stop("The four-group design is not full rank.")
  }
  design
}

wp1_contrasts <- function() {
  matrix(
    c(-0.5, 0.5, -0.5, 0.5,
      -1, 1, 0, 0,
      0, 0, -1, 1,
      1, -1, -1, 1),
    nrow = 4L,
    dimnames = list(
      wp1_group_levels(),
      c("pooled_KO_vs_WT", "female_KO_vs_WT", "male_KO_vs_WT",
        "interaction_M_minus_F")
    )
  )
}

wp1_validate_model_input <- function(y, meta) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package limma is required for Cyp2c38 models.")
  }
  if (!is.matrix(y) || !is.numeric(y) || !all(is.finite(y))) {
    stop("y must be a finite numeric feature-by-mouse matrix on the log2 scale.")
  }
  if (!is.data.frame(meta) || nrow(meta) != ncol(y)) {
    stop("The number of metadata rows must equal the number of mouse columns.")
  }
  if (is.null(rownames(y)) || anyNA(rownames(y)) ||
      any(!nzchar(rownames(y))) || anyDuplicated(rownames(y))) {
    stop("y must have nonempty, unique feature row names.")
  }
  if (is.null(colnames(y)) || anyNA(colnames(y)) ||
      any(!nzchar(colnames(y))) || anyDuplicated(colnames(y))) {
    stop("y must have nonempty, unique mouse column names.")
  }
  sample_columns <- intersect(c("sample", "sample_id"), names(meta))
  for (column in sample_columns) {
    if (!identical(as.character(meta[[column]]), colnames(y))) {
      stop("Metadata ", column, " must match y column names in the same order.")
    }
  }
  if (!all(c("group", "sex") %in% names(meta))) {
    stop("Metadata must contain group and sex columns.")
  }
  sex <- tolower(as.character(meta$sex))
  sex[sex == "female"] <- "f"
  sex[sex == "male"] <- "m"
  if (anyNA(sex) || !all(sex %in% c("f", "m")) ||
      !identical(sex, tolower(substr(as.character(meta$group), 1L, 1L)))) {
    stop("Metadata sex must be F/M or female/male and agree with group.")
  }
  invisible(TRUE)
}

wp1_fit <- function(y, meta, weighted = TRUE, trend = TRUE, robust = TRUE) {
  wp1_validate_model_input(y, meta)
  design <- wp1_design(meta)
  contrasts <- wp1_contrasts()
  if (!identical(rownames(contrasts), colnames(design))) {
    stop("Contrast rows must match the design columns.")
  }
  weights <- rep(1, ncol(y))
  if (weighted) {
    # Equal weight within each sex; precision differs between sexes.
    # This is distinct from the equal-sex weighting of the pooled estimand.
    weights <- limma::arrayWeights(
      y, design = design, var.group = factor(meta$sex)
    )
  }
  if (length(weights) != ncol(y) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Sample precision weights must be positive and finite.")
  }
  names(weights) <- colnames(y)
  fit <- limma::lmFit(y, design = design, weights = weights)
  fit <- limma::contrasts.fit(fit, contrasts = contrasts)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  list(fit = fit, weights = weights, design = design, contrasts = contrasts)
}

wp1_fit_separate <- function(y, meta, trend = TRUE, robust = TRUE) {
  wp1_validate_model_input(y, meta)
  # Validate the full intake before fitting two independently moderated models.
  wp1_design(meta)
  groups <- as.character(meta$group)
  output <- list()
  for (sex_label in c("female", "male")) {
    prefix <- if (sex_label == "female") "F" else "M"
    selected <- substr(groups, 1L, 1L) == prefix
    group <- factor(groups[selected],
                    levels = paste0(prefix, c("_WT", "_KO")))
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    contrast_name <- paste0(sex_label, "_KO_vs_WT")
    contrast <- matrix(
      c(-1, 1), ncol = 1L,
      dimnames = list(colnames(design), contrast_name)
    )
    fit <- limma::lmFit(y[, selected, drop = FALSE], design = design)
    fit <- limma::contrasts.fit(fit, contrasts = contrast)
    output[[sex_label]] <- limma::eBayes(fit, trend = trend, robust = robust)
  }
  output
}

wp1_table <- function(fit, contrast) {
  if (!inherits(fit, "MArrayLM") ||
      !is.character(contrast) || length(contrast) != 1L ||
      !contrast %in% colnames(fit$coefficients)) {
    stop("Specify a moderated limma fit and one of its named contrasts.")
  }
  result <- limma::topTable(
    fit, coef = contrast, number = Inf, sort.by = "none", adjust.method = "BH", confint = TRUE
  )
  required <- c("t", "logFC", "CI.L", "CI.R", "AveExpr", "P.Value", "adj.P.Val", "B")
  if (!all(required %in% names(result)) ||
      !identical(rownames(result), rownames(fit$coefficients))) {
    stop("The full result table must retain all features in their input order.")
  }
  data.frame(
    feature_id = rownames(result), result[, required, drop = FALSE],
    row.names = NULL, check.names = FALSE
  )
}
