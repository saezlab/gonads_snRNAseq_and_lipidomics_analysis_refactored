# Independent actual-data replay for direct CF453 diet interaction models.
# Run: Rscript tests/qa_cf453_diet_impact.R ANALYSIS_RDS [QA_DIR]
# No workflow/statistical/plotting helper is sourced by this verifier.

qa_cf453_diet_impact <- function(analysis_file,
    qa_dir = "tests/outputs/cf453_diet_impact") {
  if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
  a <- readRDS(analysis_file)
  checks <- list()
  record <- function(label, difference = 0, detail = "") {
    if (!is.finite(difference) || difference > 1e-8) {
      stop("Independent direct-diet QA failed: ", label, "; difference=", difference)
    }
    checks[[length(checks) + 1L]] <<- data.frame(check = label,
      max_abs_difference = difference, detail = detail, stringsAsFactors = FALSE)
  }
  same <- function(actual, expected, label) {
    actual <- as.numeric(actual); expected <- as.numeric(expected)
    stopifnot(length(actual) == length(expected),
              identical(is.na(actual), is.na(expected)))
    valid <- is.finite(actual) & is.finite(expected)
    stopifnot(all(valid | (is.na(actual) & is.na(expected)) |
                    (is.infinite(actual) & actual == expected)))
    record(label, if (any(valid)) max(abs(actual[valid] - expected[valid])) else 0)
  }
  required <- c("modeled_rna", "profile_metadata", "metadata", "design", "rho",
    "base_fit", "primary_contrasts", "conditional_contrasts", "fit_primary",
    "fit_conditional", "overview", "conditional", "omnibus", "heatmap_selection",
    "heatmap", "example_selection", "example_points", "example_means")
  stopifnot(all(required %in% names(a)))
  y <- as.matrix(a$modeled_rna); x <- as.matrix(a$design)
  pm <- a$profile_metadata; md <- a$metadata
  stopifnot(nrow(pm) == 56L, nrow(md) == 28L, !anyDuplicated(md$mouse_id),
    !anyDuplicated(pm$profile_id), identical(colnames(y), pm$profile_id),
    identical(rownames(x), pm$profile_id), identical(dim(x), c(56L, 8L)),
    all(table(md$sex, md$diet) == 7L), all(table(pm$mouse_id, pm$zone) == 1L),
    all(pm$sex == md$sex[match(pm$mouse_id, md$mouse_id)]),
    all(pm$diet == md$diet[match(pm$mouse_id, md$mouse_id)]),
    setequal(pm$zone, c("PP", "PC")), all(x %in% c(0, 1)),
    all(rowSums(x) == 1L), all(colSums(x) == 7L))
  record("28 independent mice and 56 correctly paired RNA profiles; balanced eight-cell design")
  cell_metadata <- do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    row <- unique(pm[x[, j] == 1L, c("sex", "diet", "zone"), drop = FALSE])
    stopifnot(nrow(row) == 1L)
    row
  }))
  cell <- function(sex, diet, zone) {
    as.numeric(cell_metadata$sex == sex & cell_metadata$diet == diet & cell_metadata$zone == zone)
  }
  ds <- function(zone) {
    cell("male", "WD", zone) - cell("female", "WD", zone) -
      cell("male", "SD", zone) + cell("female", "SD", zone)
  }
  dz <- function(sex) {
    cell(sex, "WD", "PP") - cell(sex, "WD", "PC") -
      cell(sex, "SD", "PP") + cell(sex, "SD", "PC")
  }
  expected_primary <- cbind(diet_sex = (ds("PC") + ds("PP")) / 2,
    diet_zone = (dz("female") + dz("male")) / 2,
    diet_sex_zone = ds("PP") - ds("PC"))
  pc <- as.matrix(a$primary_contrasts); sc <- as.matrix(a$conditional_contrasts)
  stopifnot(identical(colnames(pc), colnames(expected_primary)),
            identical(rownames(pc), colnames(x)), identical(rownames(sc), colnames(x)),
            qr(pc)$rank == 3L, ncol(sc) == 4L)
  same(pc, expected_primary, "Every primary contrast sign, marginal half-weight and rank-three estimand")
  expected_conditional <- cbind(ds_PC = ds("PC"), ds_PP = ds("PP"),
    dz_female = dz("female"), dz_male = dz("male"))
  map <- vapply(seq_len(ncol(sc)), function(j) {
    match <- which(vapply(seq_len(ncol(expected_conditional)), function(k)
      identical(as.numeric(sc[, j]), expected_conditional[, k]), logical(1L)))
    stopifnot(length(match) == 1L)
    match
  }, integer(1L))
  stopifnot(setequal(map, 1:4))
  same(sc, expected_conditional[, map, drop = FALSE], "All four secondary conditional interaction vectors")
  same(expected_conditional[, 1L], pc[, "diet_sex"] - pc[, "diet_sex_zone"] / 2,
       "Diet-sex change in PC decomposes as marginal minus half three-way")
  same(expected_conditional[, 2L], pc[, "diet_sex"] + pc[, "diet_sex_zone"] / 2,
       "Diet-sex change in PP decomposes as marginal plus half three-way")
  same(expected_conditional[, 3L], pc[, "diet_zone"] - pc[, "diet_sex_zone"] / 2,
       "Female diet-zonation change decomposes as marginal minus half three-way")
  same(expected_conditional[, 4L], pc[, "diet_zone"] + pc[, "diet_sex_zone"] / 2,
       "Male diet-zonation change decomposes as marginal plus half three-way")

  # Base GLS uses the shared mouse covariance; no limma fitting helper is used.
  rho <- as.numeric(a$rho); stopifnot(length(rho) == 1L, abs(rho) < 1)
  v <- matrix(0, 56L, 56L); v[outer(pm$mouse_id, pm$mouse_id, "==")] <- rho; diag(v) <- 1
  vi <- solve(v); covariance <- solve(crossprod(x, vi %*% x))
  beta <- t(covariance %*% crossprod(x, vi %*% t(y)))
  residual <- y - beta %*% t(x)
  residual_df <- 56L - qr(x)$rank
  sigma2 <- rowSums((residual %*% vi) * residual) / residual_df
  same(a$base_fit$coefficients, beta, "All eight-cell coefficients independently reconstructed by base GLS")
  same(a$base_fit$cov.coefficients, covariance, "Full paired eight-cell coefficient covariance")
  same(a$base_fit$sigma^2, sigma2, "Every RNA residual variance from independent GLS")

  replay <- function(fit, contrast, result, family, label) {
    effect <- beta %*% contrast
    contrast_covariance <- crossprod(contrast, covariance %*% contrast)
    unscaled <- sqrt(diag(contrast_covariance))
    prior_df <- rep_len(fit$df.prior, nrow(y)); prior_var <- rep_len(fit$s2.prior, nrow(y))
    posterior <- (prior_df * prior_var + residual_df * sigma2) / (prior_df + residual_df)
    posterior[is.infinite(prior_df)] <- prior_var[is.infinite(prior_df)]
    se <- sqrt(posterior) %o% unscaled
    t <- effect / se
    df <- rep_len(fit$df.total, nrow(y))
    p <- matrix(2 * stats::pt(-abs(as.vector(t)), rep(df, ncol(t))),
                nrow(y), dimnames = dimnames(effect))
    same(fit$coefficients, effect, paste(label, "all coefficients"))
    same(fit$cov.coefficients, contrast_covariance, paste(label, "full contrast covariance"))
    same(fit$s2.post, posterior, paste(label, "moderated residual variances"))
    same(fit$df.total, pmin(prior_df + residual_df, residual_df * nrow(y)),
         paste(label, "moderated degrees of freedom"))
    same(fit$t, t, paste(label, "all moderated statistics"))
    same(fit$p.value, p, paste(label, "all moderated probabilities"))
    stopifnot(nrow(result) == nrow(y) * ncol(contrast),
              !anyDuplicated(paste(result$gene, result$contrast)),
              all(result$family == family))
    ii <- cbind(match(result$gene, rownames(y)), match(result$contrast, colnames(contrast)))
    stopifnot(!anyNA(ii))
    same(result$effect, effect[ii], paste(label, "exported effect estimates"))
    same(result$se, se[ii], paste(label, "exported standard errors"))
    same(result$t, t[ii], paste(label, "exported t values"))
    same(result$p, p[ii], paste(label, "exported probabilities"))
    same(result$q, stats::p.adjust(result$p, "BH"), paste(label, "entire declared BH family"))
    if (all(c("low", "high") %in% names(result))) {
      critical <- stats::qt(.975, df[ii[, 1L]])
      same(result$low, result$effect - critical * result$se, paste(label, "nominal interval lower limits"))
      same(result$high, result$effect + critical * result$se, paste(label, "nominal interval upper limits"))
    }
    list(effect = effect, covariance = contrast_covariance, posterior = posterior)
  }
  primary <- replay(a$fit_primary, pc, a$overview, "primary_direct_diet", "Primary three-contrast fit")
  secondary <- replay(a$fit_conditional, sc, a$conditional, "conditional_direct_diet", "Secondary four-contrast fit")
  same(primary$posterior, secondary$posterior, "Primary and conditional fits preserve identical empirical Bayes moderation")

  # Independent covariance-aware joint test, restricted to the three primary
  # columns. A general diet main effect is outside this hypothesis subspace.
  quadratic <- rowSums((primary$effect %*% solve(primary$covariance)) * primary$effect)
  f <- quadratic / (3 * primary$posterior)
  o <- a$omnibus; oi <- match(o$gene, rownames(y))
  stopifnot(nrow(o) == nrow(y), !anyDuplicated(o$gene), !anyNA(oi), all(o$df1 == 3L))
  same(o$F, f[oi], "Every covariance-aware three-df moderated omnibus F statistic")
  same(o$df2, rep_len(a$fit_primary$df.prior, nrow(y))[oi] + residual_df,
       "Omnibus denominator degrees of freedom")
  same(o$p, stats::pf(f[oi], df1 = 3, df2 = o$df2, lower.tail = FALSE),
       "Every omnibus F probability from recorded denominator degrees of freedom")
  same(o$q, stats::p.adjust(o$p, "BH"), "Omnibus BH family across all eligible genes")

  # All selection is based on direct interactions or the explicit omnibus test.
  expected_heatmap <- o[is.finite(o$q) & o$q <= .05, , drop = FALSE]
  expected_heatmap <- head(expected_heatmap[order(expected_heatmap$q, expected_heatmap$p,
                                                expected_heatmap$gene), , drop = FALSE], 20L)
  hs <- a$heatmap_selection
  stopifnot(identical(hs$gene, expected_heatmap$gene), identical(as.integer(hs$rank), seq_len(nrow(hs))))
  for (field in c("F", "p", "q")) same(hs[[field]], expected_heatmap[[field]], paste("Heatmap ranking", field))
  hm <- a$heatmap
  stopifnot(nrow(hm) == 4L * nrow(hs), setequal(hm$gene, hs$gene),
            !anyDuplicated(paste(hm$gene, hm$contrast)))
  ci <- match(paste(hm$gene, hm$contrast), paste(a$conditional$gene, a$conditional$contrast))
  stopifnot(!anyNA(ci))
  for (field in c("effect", "p", "q")) same(hm[[field]], a$conditional[[field]][ci], paste("Heatmap source", field))
  same(hm$joint_q, o$q[match(hm$gene, o$gene)], "Heatmap direct omnibus adjusted probabilities")
  threeway <- a$overview[a$overview$contrast == "diet_sex_zone", , drop = FALSE]
  same(hm$threeway_q, threeway$q[match(hm$gene, threeway$gene)], "Heatmap new-family three-way adjusted probabilities")

  slots <- c("diet_sex", "diet_zone", "diet_sex_zone", "diet_sex_zone")
  chosen <- character(); selections <- list()
  for (j in seq_along(slots)) {
    d <- a$overview[a$overview$contrast == slots[j] & a$overview$q <= .05 &
                      !a$overview$gene %in% chosen, , drop = FALSE]
    d <- d[order(-abs(d$effect), d$q, d$p, d$gene), , drop = FALSE]
    if (!nrow(d)) next
    d <- d[1L, , drop = FALSE]; chosen <- c(chosen, d$gene)
    selections[[length(selections) + 1L]] <- data.frame(gene = d$gene,
      example_rank = j, selected_contrast = slots[j], q_primary = d$q,
      stringsAsFactors = FALSE)
  }
  es <- a$example_selection
  expected_selection <- if (length(selections)) do.call(rbind, selections) else
    data.frame(gene = character(), example_rank = integer(), selected_contrast = character(), q_primary = numeric())
  stopifnot(identical(es$gene, expected_selection$gene),
            identical(es$selected_contrast, expected_selection$selected_contrast))
  same(es$example_rank, expected_selection$example_rank, "Distinct direct-effect example rank slots")
  same(es$q_primary, expected_selection$q_primary, "Examples selected solely by direct primary BH evidence")
  for (j in seq_len(nrow(es))) {
    d <- es[j, ]; point <- a$example_points[a$example_points$gene == d$gene, , drop = FALSE]
    stopifnot(nrow(point) == 28L, !anyDuplicated(point$mouse_id),
      setequal(point$mouse_id, md$mouse_id), all(table(point$sex, point$diet) == 7L))
    pp <- pm$profile_id[match(paste(point$mouse_id, "PP"), paste(pm$mouse_id, pm$zone))]
    pc_profile <- pm$profile_id[match(paste(point$mouse_id, "PC"), paste(pm$mouse_id, pm$zone))]
    values <- if (d$selected_contrast == "diet_sex")
      (y[d$gene, pp] + y[d$gene, pc_profile]) / 2 else y[d$gene, pp] - y[d$gene, pc_profile]
    same(point$value, values, paste(d$gene, "all independent-mouse averages or paired-zone differences"))
    stopifnot(all(point$sex == md$sex[match(point$mouse_id, md$mouse_id)]),
              all(point$diet == md$diet[match(point$mouse_id, md$mouse_id)]))
    means <- a$example_means[a$example_means$gene == d$gene, , drop = FALSE]
    stopifnot(nrow(means) == 4L, !anyDuplicated(paste(means$sex, means$diet)))
    g <- match(d$gene, rownames(y))
    for (k in seq_len(nrow(means))) {
      sex <- means$sex[k]; diet <- means$diet[k]
      w <- if (d$selected_contrast == "diet_sex")
        (cell(sex, diet, "PP") + cell(sex, diet, "PC")) / 2 else
          cell(sex, diet, "PP") - cell(sex, diet, "PC")
      estimate <- as.numeric(beta[g, ] %*% w)
      standard_error <- sqrt(as.numeric(crossprod(w, covariance %*% w)) * primary$posterior[g])
      critical <- stats::qt(.975, a$fit_primary$df.total[g])
      label <- paste(d$gene, sex, diet)
      same(means$estimate[k], estimate, paste(label, "model mean"))
      same(means$se[k], standard_error, paste(label, "paired-covariance posterior standard error"))
      same(means$low[k], estimate - critical * standard_error, paste(label, "nominal lower interval"))
      same(means$high[k], estimate + critical * standard_error, paste(label, "nominal upper interval"))
      same(means$estimate[k], mean(point$value[point$sex == sex & point$diet == diet]),
           paste(label, "actual animal summary matches model mean"))
    }
  }

  # The prior fit is a validated same-CF453 input, not another cohort or a
  # network-selected feature roster. Raw three-way evidence must be unchanged.
  source_paths <- unlist(a$source_paths, use.names = TRUE)
  parent_path <- source_paths[basename(source_paths) == "models.rds"]
  stopifnot(length(parent_path) == 1L, file.exists(parent_path),
            !any(grepl("CF459|scaled_RNA|results/zonation/|results/mofa/", source_paths)))
  parent <- readRDS(parent_path)
  same(y, parent$modeled_rna, "All modeled RNA values unchanged from validated same-cohort parent")
  same(x, parent$design, "Eight-cell design unchanged from validated same-cohort parent")
  same(a$rho, parent$rho, "Paired covariance unchanged from validated same-cohort parent")
  old <- parent$rna_results[parent$rna_results$contrast == "sex_diet_zone", , drop = FALSE]
  old <- old[match(threeway$gene, old$gene), , drop = FALSE]
  for (field in c("effect", "se", "t", "p")) same(threeway[[field]], old[[field]],
    paste("Raw three-way parent replay; new multiplicity family only", field))

  prior_manifest <- utils::read.csv(file.path(dirname(parent_path), "provenance_sha256.csv"),
                                    stringsAsFactors = FALSE)
  raw_paths <- unname(parent$input_paths)
  hashes <- vapply(raw_paths, digest::digest, character(1L), file = TRUE, algo = "sha256")
  expected_hashes <- prior_manifest$sha256[match(raw_paths, prior_manifest$path)]
  stopifnot(!anyNA(expected_hashes), identical(unname(hashes), unname(expected_hashes)),
            all(raw_paths %in% source_paths))
  record("All three original CF453 input SHA256 values match validated parent provenance")

  # Reconstruct the whole-liver lipid model with 28 independent observations;
  # sex-by-diet is a difference of two diet effects and has no zone replication.
  lipid <- as.matrix(a$lipid_matrix); lx <- as.matrix(a$lipid_design)
  lc <- as.matrix(a$lipid_contrasts); lf <- a$lipid_fit
  stopifnot(identical(colnames(lipid), md$mouse_id), ncol(lipid) == 28L,
            identical(dim(lx), c(28L, 4L)), all(rowSums(lx) == 1L),
            all(colSums(lx) == 7L))
  same(lipid, parent$lipid_matrix, "All primary whole-liver lipid inputs unchanged")
  lipid_cell <- function(sex, diet) {
    result <- numeric(ncol(lx))
    for (j in seq_len(ncol(lx))) {
      rows <- md[lx[, j] == 1L, , drop = FALSE]
      stopifnot(length(unique(rows$sex)) == 1L, length(unique(rows$diet)) == 1L)
      result[j] <- rows$sex[1L] == sex && rows$diet[1L] == diet
    }
    stopifnot(sum(result) == 1L)
    result
  }
  lipid_expected <- list(
    sex_SD = lipid_cell("male", "SD") - lipid_cell("female", "SD"),
    sex_WD = lipid_cell("male", "WD") - lipid_cell("female", "WD"),
    diet_female = lipid_cell("female", "WD") - lipid_cell("female", "SD"),
    diet_male = lipid_cell("male", "WD") - lipid_cell("male", "SD"))
  lipid_expected$sex_diet <- lipid_expected$diet_male - lipid_expected$diet_female
  for (name in names(lipid_expected)) same(lc[, name], lipid_expected[[name]],
    paste("Independent lipid contrast", name))
  lipid_covariance <- solve(crossprod(lx))
  lipid_beta <- t(lipid_covariance %*% crossprod(lx, t(lipid)))
  lipid_effect <- lipid_beta %*% lc
  lipid_resid <- lipid - lipid_beta %*% t(lx)
  lipid_sigma2 <- rowSums(lipid_resid^2) / 24
  lipid_prior_df <- rep_len(lf$df.prior, nrow(lipid))
  lipid_prior_var <- rep_len(lf$s2.prior, nrow(lipid))
  lipid_post <- (lipid_prior_df * lipid_prior_var + 24 * lipid_sigma2) / (lipid_prior_df + 24)
  lipid_post[is.infinite(lipid_prior_df)] <- lipid_prior_var[is.infinite(lipid_prior_df)]
  lipid_se <- sqrt(lipid_post) %o% sqrt(diag(crossprod(lc, lipid_covariance %*% lc)))
  lipid_t <- lipid_effect / lipid_se
  lipid_df <- rep_len(lf$df.total, nrow(lipid))
  lipid_p <- matrix(2 * stats::pt(-abs(as.vector(lipid_t)), rep(lipid_df, ncol(lipid_t))),
                    nrow(lipid), dimnames = dimnames(lipid_effect))
  same(lf$coefficients, lipid_effect, "All whole-liver lipid coefficients from independent OLS")
  same(lf$sigma^2, lipid_sigma2, "All whole-liver lipid residual variances")
  same(lf$s2.post, lipid_post, "All whole-liver lipid posterior variances")
  same(lf$t, lipid_t, "All whole-liver lipid moderated statistics")
  same(lf$p.value, lipid_p, "All whole-liver lipid moderated probabilities")
  lr <- a$lipid_results
  for (family in unique(lr$family)) {
    rows <- which(lr$family == family)
    same(lr$q[rows], stats::p.adjust(lr$p[rows], "BH"), paste("Retained lipid BH family", family))
  }
  ls <- a$lipid_scatter
  stopifnot(nrow(ls) == nrow(lipid), !anyDuplicated(ls$lipid),
            setequal(ls$lipid, rownames(lipid)))
  li <- match(ls$lipid, rownames(lipid))
  same(ls$effect_female, lipid_effect[li, "diet_female"], "Lipid scatter female WD-minus-SD effects")
  same(ls$effect_male, lipid_effect[li, "diet_male"], "Lipid scatter male WD-minus-SD effects")
  same(ls$effect_interaction, lipid_effect[li, "sex_diet"], "Lipid scatter direct sex-by-diet effects")
  same(ls$effect_interaction, ls$effect_male - ls$effect_female,
       "Every lipid interaction equals male-minus-female diet change")
  same(ls$effect_interaction, ls$sex_WD - ls$sex_SD,
       "Every lipid interaction equals WD-minus-SD sex change")
  same(ls$p_interaction, lipid_p[li, "sex_diet"], "Lipid scatter direct interaction probabilities")
  same(ls$q_interaction, stats::p.adjust(ls$p_interaction, "BH"),
       "Lipid scatter BH family includes every primary interaction species")
  critical <- stats::qt(.975, lipid_df[li])
  for (name in c("female", "male", "interaction")) {
    contrast <- if (name == "interaction") "sex_diet" else paste0("diet_", name)
    estimate <- lipid_effect[li, contrast]; standard_error <- lipid_se[li, contrast]
    same(ls[[paste0("low_", name)]], estimate - critical * standard_error,
         paste("Lipid scatter nominal lower interval", name))
    same(ls[[paste0("high_", name)]], estimate + critical * standard_error,
         paste("Lipid scatter nominal upper interval", name))
  }

  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- do.call(rbind, checks)
  utils::write.csv(summary, file.path(qa_dir, "independent_checks.csv"), row.names = FALSE)
  saveRDS(list(checks = summary, source_analysis = normalizePath(analysis_file),
               primary_contrast_covariance = primary$covariance,
               expected_omnibus_F = f, session = utils::sessionInfo()),
          file.path(qa_dir, "independent_qa.rds"))
  cat("Independent CF453 direct-diet QA passed:", nrow(summary), "check groups.\n")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || identical(args, "--help")) {
    cat("Usage: Rscript tests/qa_cf453_diet_impact.R ANALYSIS_RDS [QA_DIR]\n")
  } else {
    qa_cf453_diet_impact(args[1L], if (length(args) >= 2L) args[2L] else
      "tests/outputs/cf453_diet_impact")
  }
}
