# Whole-mouse correlation-network assessment for the WDSD cohort.
# Sourcing this file only defines functions; it does not change the RNG state.
# The elementary sampling unit is a mouse, including both RNA zones and its
# whole-liver lipid vector. Correlation-network edges are not independent mice.

wp2_vector_stats <- function(reference_r, observed_r) {
  if (!is.numeric(reference_r) || !is.numeric(observed_r) ||
      !is.null(dim(reference_r)) || !is.null(dim(observed_r)) ||
      length(reference_r) != length(observed_r)) {
    stop("Supply two equal-length numeric coefficient vectors.")
  }
  n <- length(reference_r)
  result <- c(coefficient_vector_r = NA_real_, ULM_score = NA_real_,
              ULM_p_value = NA_real_, slope = NA_real_, n_edges = n)
  # No pairwise deletion: an undefined edge makes this frozen vector undefined.
  if (n < 3L || any(!is.finite(reference_r)) || any(!is.finite(observed_r))) {
    return(result)
  }
  xscale <- max(abs(reference_r)); yscale <- max(abs(observed_r))
  if (xscale == 0 || yscale == 0) return(result)
  x <- reference_r / xscale; y <- observed_r / yscale
  x <- x - mean(x); y <- y - mean(y)
  xx <- sum(x * x); yy <- sum(y * y); xy <- sum(x * y)
  if (xx <= n * .Machine$double.eps || yy <= n * .Machine$double.eps) return(result)
  rho <- max(-1, min(1, xy / sqrt(xx * yy)))
  score <- rho * sqrt((n - 2) / (1 - rho * rho))
  c(coefficient_vector_r = rho, ULM_score = score,
    ULM_p_value = 2 * stats::pt(abs(score), df = n - 2, lower.tail = FALSE),
    slope = (xy / xx) * (yscale / xscale), n_edges = n)
}

wp2_validate_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) ||
      is.null(colnames(x)) || anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(x)) || anyNA(rownames(x)) || anyNA(colnames(x))) {
    stop(label, " must be a numeric feature-by-mouse matrix with unique names.")
  }
  invisible(TRUE)
}

wp2_unit_rows <- function(x) {
  # Scale first to make constant-feature detection relative to the feature's
  # numerical scale; a variable but very small abundance is not constant.
  # pmax operates over whole feature vectors in native code. At seven mice,
  # this avoids invoking an R function once per feature on every bootstrap draw
  # while preserving the same scale and nonfinite-feature behavior.
  scale <- do.call(pmax, lapply(seq_len(ncol(x)), function(j) abs(x[, j])))
  finite <- rowSums(!is.finite(x)) == 0L & is.finite(scale) & scale > 0
  scale[!finite] <- 1
  z <- x / scale
  z <- z - rowMeans(z)
  ss <- rowSums(z * z)
  valid <- finite & is.finite(ss) & ss > ncol(x) * .Machine$double.eps
  z[!valid, ] <- 0
  divisor <- sqrt(ss); divisor[!valid] <- 1
  list(unit = z / divisor, valid = valid)
}

wp2_prepare_edges <- function(rna, lipids, edges) {
  wp2_validate_matrix(rna, "RNA"); wp2_validate_matrix(lipids, "Lipids")
  if (!identical(colnames(rna), colnames(lipids))) {
    stop("RNA and lipid mouse IDs must be identical and in the same order.")
  }
  if (!all(c("gene", "feature_id") %in% names(edges))) {
    stop("edges requires gene and feature_id (or use wp2_edge_correlations with lipid).")
  }
  genes <- unique(as.character(edges$gene))
  features <- unique(as.character(edges$feature_id))
  # Missing input rows remain undefined; they are never silently fabricated.
  gm <- match(genes, rownames(rna)); lm <- match(features, rownames(lipids))
  g <- matrix(NA_real_, length(genes), ncol(rna), dimnames = list(genes, colnames(rna)))
  l <- matrix(NA_real_, length(features), ncol(lipids), dimnames = list(features, colnames(lipids)))
  if (any(!is.na(gm))) g[!is.na(gm), ] <- rna[gm[!is.na(gm)], , drop = FALSE]
  if (any(!is.na(lm))) l[!is.na(lm), ] <- lipids[lm[!is.na(lm)], , drop = FALSE]
  gu <- wp2_unit_rows(g); lu <- wp2_unit_rows(l)
  list(rna = g, lipids = l, rna_unit = gu, lipid_unit = lu,
       gene_index = match(as.character(edges$gene), genes),
       lipid_index = match(as.character(edges$feature_id), features),
       mouse_ids = colnames(rna), n_edges = nrow(edges))
}

wp2_prepared_cor <- function(prepared, rna_indices = NULL, lipid_indices = rna_indices,
                             fixed_margins = FALSE) {
  n <- ncol(prepared$rna)
  if (is.null(rna_indices)) rna_indices <- seq_len(n)
  if (is.null(lipid_indices)) lipid_indices <- seq_len(n)
  if (length(rna_indices) != length(lipid_indices) || length(rna_indices) < 3L ||
      anyNA(rna_indices) || anyNA(lipid_indices) ||
      any(rna_indices < 1L | rna_indices > n) ||
      any(lipid_indices < 1L | lipid_indices > n)) {
    stop("Resampling indices must be equal length, valid mouse positions, with at least three entries.")
  }
  if (fixed_margins) {
    # Only a complete permutation has unchanged centering and normalization.
    if (!identical(sort(as.integer(rna_indices)), seq_len(n)) ||
        !identical(sort(as.integer(lipid_indices)), seq_len(n))) {
      stop("fixed_margins requires complete permutations without replacement.")
    }
    gu <- prepared$rna_unit; lu <- prepared$lipid_unit
    gu$unit <- gu$unit[, rna_indices, drop = FALSE]
    lu$unit <- lu$unit[, lipid_indices, drop = FALSE]
  } else {
    gu <- wp2_unit_rows(prepared$rna[, rna_indices, drop = FALSE])
    lu <- wp2_unit_rows(prepared$lipids[, lipid_indices, drop = FALSE])
  }
  cross <- tcrossprod(gu$unit, lu$unit)
  answer <- cross[cbind(prepared$gene_index, prepared$lipid_index)]
  bad <- !gu$valid[prepared$gene_index] | !lu$valid[prepared$lipid_index]
  answer[bad] <- NA_real_
  answer[!bad] <- pmax(-1, pmin(1, answer[!bad]))
  answer
}

wp2_cor_edges <- function(rna, lipids, edges) {
  wp2_prepared_cor(wp2_prepare_edges(rna, lipids, edges))
}

wp2_edge_correlations <- function(rna, lipid, edges, rna_indices = NULL,
                                  lipid_indices = rna_indices) {
  if (!"feature_id" %in% names(edges) && "lipid" %in% names(edges)) {
    edges$feature_id <- as.character(edges$lipid)
  }
  wp2_prepared_cor(wp2_prepare_edges(rna, lipid, edges), rna_indices, lipid_indices)
}

wp2_all_permutations <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L || n > 8L) {
    stop("Exact pairing permutation is supported for one to eight mice.")
  }
  recurse <- function(values) {
    if (length(values) == 1L) return(matrix(values, nrow = 1L))
    do.call(rbind, lapply(values, function(first) {
      cbind(first, recurse(values[values != first]))
    }))
  }
  unname(recurse(seq_len(n)))
}

wp2_normalize_context <- function(x, name, global_edges = NULL) {
  aliases <- c(RNA_SD = "rna_SD", RNA_WD = "rna_WD", lipidSD = "lipid_SD", lipidWD = "lipid_WD")
  for (field in names(aliases)) if (is.null(x[[aliases[[field]]]]) && !is.null(x[[field]])) {
    x[[aliases[[field]]]] <- x[[field]]
  }
  if (!is.null(global_edges)) x$edges <- global_edges[global_edges$context == name, , drop = FALSE]
  required <- c("rna_SD", "rna_WD", "lipid_SD", "lipid_WD", "edges")
  if (!all(required %in% names(x))) stop("Context ", name, " requires ", paste(required, collapse = ", "), ".")
  e <- x$edges
  if (!"feature_id" %in% names(e) && "lipid" %in% names(e)) e$feature_id <- as.character(e$lipid)
  if (!"r_reference" %in% names(e) && "reference_r" %in% names(e)) e$r_reference <- e$reference_r
  if (!all(c("gene", "feature_id", "r_reference", "primary_edge") %in% names(e))) {
    stop("Context edges require gene, feature_id/lipid, r_reference/reference_r and primary_edge.")
  }
  if (!is.logical(e$primary_edge) || anyNA(e$primary_edge)) stop("primary_edge must be a complete logical vector.")
  if (!"edge_id" %in% names(e)) e$edge_id <- paste(e$gene, e$feature_id, sep = "::")
  if (anyDuplicated(e$edge_id) || anyDuplicated(paste(e$gene, e$feature_id, sep = "\r"))) {
    stop("Duplicate edge identities in context ", name, ".")
  }
  x$edges <- e
  x$SD <- wp2_prepare_edges(x$rna_SD, x$lipid_SD, e)
  x$WD <- wp2_prepare_edges(x$rna_WD, x$lipid_WD, e)
  if (length(intersect(x$SD$mouse_ids, x$WD$mouse_ids))) stop("SD and WD must contain different mice.")
  if (ncol(x$rna_SD) != 7L || ncol(x$rna_WD) != 7L) {
    stop("This WDSD analysis contract requires seven SD and seven WD mice per sex.")
  }
  e$r_SD <- wp2_prepared_cor(x$SD)
  e$r_WD <- wp2_prepared_cor(x$WD)
  e$eligible <- is.finite(e$r_reference) & is.finite(e$r_SD) & is.finite(e$r_WD)
  e$exclusion_reason <- ifelse(!is.finite(e$r_reference), "undefined_reference",
                              ifelse(!is.finite(e$r_SD), "undefined_SD",
                                     ifelse(!is.finite(e$r_WD), "undefined_WD", "")))
  e$context <- name
  e$sex <- if (!is.null(x$sex)) x$sex else sub("_.*$", "", name)
  e$zone <- if (!is.null(x$zone)) x$zone else sub("^[^_]+_", "", name)
  e$delta_r <- e$r_WD - e$r_SD
  e$abs_delta_r <- abs(e$delta_r)
  e$sign_change <- sign(e$r_WD) != sign(e$r_SD)
  e$n_SD <- ncol(x$rna_SD); e$n_WD <- ncol(x$rna_WD)
  x$edges <- e
  # Freeze the observed eligible set before resampling, identically for SD/WD.
  x$eligible_index <- which(e$eligible)
  if (!is.null(x$reference_rna) && !is.null(x$reference_lipid)) {
    x$REF <- wp2_prepare_edges(x$reference_rna, x$reference_lipid, e)
    m <- x$reference_metadata
    if (is.null(m) || !all(c("mouse_id", "group") %in% names(m)) ||
        !identical(as.character(m$mouse_id), x$REF$mouse_ids) || anyNA(m$group)) {
      stop("reference_metadata must give ordered reference mouse_id/group.")
    }
    replay <- wp2_prepared_cor(x$REF)
    finite <- which(e$eligible)
    if (length(finite) && (any(!is.finite(replay[finite])) ||
                          max(abs(replay[finite] - e$r_reference[finite])) > 1e-8)) {
      stop("Reference input matrices do not replay frozen coefficients in ", name, ".")
    }
  }
  x
}

wp2_specifications <- function(contexts, focus_genes) {
  specs <- list()
  for (name in names(contexts)) {
    e <- contexts[[name]]$edges
    for (universe in c("full", "primary")) {
      index <- which(e$eligible & (universe == "full" | e$primary_edge))
      if (!length(index)) next
      scopes <- list(network = index)
      for (gene in intersect(focus_genes, unique(e$gene[index]))) {
        scopes[[paste0("gene::", gene)]] <- index[e$gene[index] == gene]
      }
      for (scope_name in names(scopes)) {
        ix <- scopes[[scope_name]]
        if (length(ix) < 3L) next
        key <- paste(name, universe, scope_name, sep = "__")
        specs[[key]] <- list(context = name, universe = universe,
                             scope = if (scope_name == "network") "network" else "gene",
                             gene = if (scope_name == "network") NA_character_ else sub("^gene::", "", scope_name),
                             index = ix, reference = e$r_reference[ix])
      }
    }
  }
  specs
}

wp2_p_from_null <- function(observed, null, exact = FALSE) {
  valid <- is.finite(null)
  nvalid <- sum(valid)
  if (!is.finite(observed) || !nvalid) return(NA_real_)
  # Floating tolerance retains the identity permutation and exact ties.
  exceed <- sum(abs(null[valid]) >= abs(observed) - 1e-12)
  if (exact) exceed / nvalid else (1 + exceed) / (1 + nvalid)
}

wp2_edge_change_stability <- function(contexts, edges = NULL) {
  contexts <- lapply(names(contexts), function(name) wp2_normalize_context(contexts[[name]], name, edges)) |>
    stats::setNames(names(contexts))
  matrices <- list(); ledgers <- list(); summaries <- list(); all_edges <- list()
  for (name in names(contexts)) {
    x <- contexts[[name]]; e <- x$edges
    # All 49 pairs drop one independent SD mouse and one independent WD mouse.
    # The same mouse-position omissions are applied to both zones by the caller
    # contract checked in wp2_network_analysis.
    delta <- matrix(NA_real_, nrow(e), 49L)
    reversals <- matrix(NA, nrow(e), 49L)
    ledger <- expand.grid(omit_SD = seq_len(7L), omit_WD = seq_len(7L))
    for (j in seq_len(nrow(ledger))) {
      sd <- wp2_prepared_cor(x$SD, setdiff(seq_len(7L), ledger$omit_SD[j]))
      wd <- wp2_prepared_cor(x$WD, setdiff(seq_len(7L), ledger$omit_WD[j]))
      delta[, j] <- wd - sd
      reversals[, j] <- sign(wd) != sign(sd)
    }
    rownames(delta) <- rownames(reversals) <- e$edge_id
    colnames(delta) <- colnames(reversals) <- paste0("SD", ledger$omit_SD, "_WD", ledger$omit_WD)
    valid <- is.finite(delta)
    e$loco_valid_combinations <- rowSums(valid)
    e$loco_delta_min <- apply(delta, 1L, function(z) if (any(is.finite(z))) min(z[is.finite(z)]) else NA_real_)
    e$loco_delta_max <- apply(delta, 1L, function(z) if (any(is.finite(z))) max(z[is.finite(z)]) else NA_real_)
    e$loco_min_abs_delta <- apply(delta, 1L, function(z) if (any(is.finite(z))) min(abs(z[is.finite(z)])) else NA_real_)
    e$loco_same_delta_direction_fraction <- rowSums(sign(delta) == sign(e$delta_r), na.rm = TRUE) / e$loco_valid_combinations
    e$loco_sign_reversal_fraction <- rowSums(reversals, na.rm = TRUE) / e$loco_valid_combinations
    e$loco_delta_crosses_zero <- e$loco_delta_min <= 0 & e$loco_delta_max >= 0
    e$loco_all_49_same_delta_direction <- e$loco_valid_combinations == 49L & e$loco_same_delta_direction_fraction == 1
    e$loco_all_49_sign_reversal <- e$loco_valid_combinations == 49L & e$loco_sign_reversal_fraction == 1
    ledger$context <- name; ledger$sex <- unique(e$sex)
    ledger$omitted_SD_mouse <- x$SD$mouse_ids[ledger$omit_SD]
    ledger$omitted_WD_mouse <- x$WD$mouse_ids[ledger$omit_WD]
    matrices[[name]] <- list(delta_r = delta, sign_reversal = reversals)
    ledgers[[name]] <- ledger; all_edges[[name]] <- e
    for (universe in c("full", "primary")) {
      ix <- which(e$eligible & (universe == "full" | e$primary_edge))
      if (!length(ix)) next
      summaries[[paste(name, universe)]] <- data.frame(
        context = name, sex = unique(e$sex), zone = unique(e$zone), universe = universe,
        n_edges = length(ix), median_abs_delta_r = stats::median(e$abs_delta_r[ix]),
        mean_abs_delta_r = mean(e$abs_delta_r[ix]), sign_change_fraction = mean(e$sign_change[ix]),
        all_49_defined_fraction = mean(e$loco_valid_combinations[ix] == 49L),
        all_49_same_delta_direction_fraction = mean(e$loco_all_49_same_delta_direction[ix]),
        all_49_sign_reversal_fraction = mean(e$loco_all_49_sign_reversal[ix]),
        stringsAsFactors = FALSE)
    }
  }
  list(edges = do.call(rbind, all_edges), loco_matrices = matrices,
       loco_omission_ledger = do.call(rbind, ledgers),
       loco_aggregate_summary = do.call(rbind, summaries))
}

#' Network conservation with explicit mouse-level nulls and sensitivity checks.
#'
#' Exact pairing permutations condition on complete observed marginal vectors
#' and require exchangeable mice within sex x diet. They test absence of RNA-
#' lipid pairing information, using reference alignment as the chosen statistic.
#' The independent empirical-marginal bootstrap samples RNA mouse vectors and
#' lipid mouse vectors independently with replacement, imposing independence on
#' the empirical margins. This is an approximate small-sample null bootstrap.
#' Neither null is the weaker composite hypothesis of zero alignment with FCG
#' while arbitrary nonaligned RNA-lipid associations remain.
#'
#' Native ULM p-values assume independent edge residuals and are descriptive.
#' Ordinary paired bootstrap draws retain associations and provide stability
#' fractions only; they do not generate null p-values. All features and both RNA
#' zones share the mouse-level draw, and the frozen edge set never shrinks within
#' a draw. Undefined edges invalidate the affected whole summary statistic.
wp2_network_analysis <- function(contexts, B = 4999L, paired_B = 999L,
                                 seed = 20261007L, focus_genes = c("Mtmr11"),
                                 workers = 1L, edges = NULL, verbose = TRUE) {
  if (is.null(names(contexts)) || anyDuplicated(names(contexts)) || !length(contexts)) {
    stop("contexts must be a nonempty, uniquely named list.")
  }
  if (B < 1L || paired_B < 1L || workers < 1L) stop("Draw counts and workers must be positive.")
  # workers is recorded for compatibility; calculations remain deterministic
  # in a single process and share exact draw banks across the two zones.
  contexts <- lapply(names(contexts), function(name) wp2_normalize_context(contexts[[name]], name, edges)) |>
    stats::setNames(names(contexts))
  sexes <- unique(vapply(contexts, function(x) unique(x$edges$sex), character(1L)))
  for (sex in sexes) {
    cx <- contexts[vapply(contexts, function(x) unique(x$edges$sex) == sex, logical(1L))]
    for (diet in c("SD", "WD")) {
      first <- cx[[1L]][[diet]]$mouse_ids
      if (any(!vapply(cx, function(x) identical(first, x[[diet]]$mouse_ids), logical(1L)))) {
        stop("Both zones must have the same ordered mouse IDs within sex/diet.")
      }
    }
    refs <- cx[vapply(cx, function(x) !is.null(x$REF), logical(1L))]
    if (length(refs) && (length(refs) != length(cx) ||
                         any(!vapply(refs, function(x) identical(refs[[1L]]$REF$mouse_ids, x$REF$mouse_ids) &&
                                                    identical(refs[[1L]]$reference_metadata$group, x$reference_metadata$group), logical(1L))))) {
      stop("Both zones must share the same ordered reference mice and group strata.")
    }
  }
  specs <- wp2_specifications(contexts, unique(as.character(focus_genes)))
  if (!length(specs)) stop("No eligible network vector with at least three edges.")
  labels <- names(specs)
  observed <- lapply(specs, function(s) wp2_vector_stats(s$reference, contexts[[s$context]]$edges$r_SD[s$index]))
  permutations <- wp2_all_permutations(7L)
  perm_bank <- matrix(NA_real_, nrow(permutations), length(specs), dimnames = list(NULL, labels))
  null_bank <- matrix(NA_real_, B, length(specs), dimnames = list(NULL, labels))
  paired_bank <- combined_bank <- matrix(NA_real_, paired_B, length(specs), dimnames = list(NULL, labels))
  set.seed(seed)
  draws <- list(); draw_ledgers <- list()
  for (sex in sexes) {
    sd0 <- contexts[[which(vapply(contexts, function(x) unique(x$edges$sex) == sex, logical(1L)))[1L]]]$SD
    bank <- list(null_rna = matrix(sample.int(7L, 7L * B, replace = TRUE), B, 7L),
                 null_lipid = matrix(sample.int(7L, 7L * B, replace = TRUE), B, 7L),
                 paired_SD = matrix(sample.int(7L, 7L * paired_B, replace = TRUE), paired_B, 7L))
    cx <- contexts[vapply(contexts, function(x) unique(x$edges$sex) == sex, logical(1L))]
    if (!is.null(cx[[1L]]$REF)) {
      strata <- split(seq_len(ncol(cx[[1L]]$REF$rna)), as.character(cx[[1L]]$reference_metadata$group))
      bank$paired_reference <- t(vapply(seq_len(paired_B), function(j) {
        index <- seq_len(ncol(cx[[1L]]$REF$rna))
        for (stratum in strata) index[stratum] <- stratum[sample.int(length(stratum), length(stratum), replace = TRUE)]
        index
      }, integer(ncol(cx[[1L]]$REF$rna))))
    }
    draws[[sex]] <- bank
    for (method in names(bank)) {
      z <- bank[[method]]
      ids <- if (method == "paired_reference") cx[[1L]]$REF$mouse_ids else sd0$mouse_ids
      draw_ledgers[[paste(sex, method)]] <- data.frame(
        sex = sex, method = method, draw = rep(seq_len(nrow(z)), each = ncol(z)),
        position = rep(seq_len(ncol(z)), times = nrow(z)),
        source_position = as.vector(t(z)), source_mouse = ids[as.vector(t(z))], stringsAsFactors = FALSE)
    }
  }
  fill_stats <- function(bank, row, context, values, reference_values = NULL) {
    hits <- which(vapply(specs, function(s) s$context == context, logical(1L)))
    for (k in hits) {
      s <- specs[[k]]
      ref <- if (is.null(reference_values)) s$reference else reference_values[s$index]
      bank[row, k] <- wp2_vector_stats(ref, values[s$index])["coefficient_vector_r"]
    }
    bank
  }
  for (name in names(contexts)) {
    x <- contexts[[name]]; sex <- unique(x$edges$sex); draw <- draws[[sex]]
    hits <- which(vapply(specs, function(s) s$context == name, logical(1L)))
    calculate <- function(values, reference_values = NULL) {
      vapply(hits, function(k) {
        s <- specs[[k]]
        ref <- if (is.null(reference_values)) s$reference else reference_values[s$index]
        wp2_vector_stats(ref, values[s$index])["coefficient_vector_r"]
      }, numeric(1L))
    }
    if (verbose) message("WDSD resampling ", name, ": ", nrow(x$edges), " measured edge candidates; ", length(x$eligible_index), " eligible.")
    for (j in seq_len(nrow(permutations))) {
      val <- wp2_prepared_cor(x$SD, seq_len(7L), permutations[j, ], fixed_margins = TRUE)
      perm_bank[j, hits] <- calculate(val)
    }
    for (j in seq_len(B)) {
      val <- wp2_prepared_cor(x$SD, draw$null_rna[j, ], draw$null_lipid[j, ])
      null_bank[j, hits] <- calculate(val)
    }
    for (j in seq_len(paired_B)) {
      val <- wp2_prepared_cor(x$SD, draw$paired_SD[j, ])
      paired_bank[j, hits] <- calculate(val)
      if (!is.null(x$REF)) {
        ref <- wp2_prepared_cor(x$REF, draw$paired_reference[j, ])
        combined_bank[j, hits] <- calculate(val, ref)
      }
    }
  }
  inference <- do.call(rbind, lapply(seq_along(specs), function(k) {
    s <- specs[[k]]; x <- contexts[[s$context]]; z <- observed[[k]]
    nv <- sum(is.finite(null_bank[, k])); pv <- sum(is.finite(perm_bank[, k]))
    paired_valid <- sum(is.finite(paired_bank[, k])); combined_valid <- sum(is.finite(combined_bank[, k]))
    rho <- unname(z["coefficient_vector_r"])
    data.frame(
      specification = labels[k], context = s$context, sex = unique(x$edges$sex), zone = unique(x$edges$zone),
      universe = s$universe, scope = s$scope, gene = s$gene, n_edges = unname(z["n_edges"]),
      n_SD = 7L, n_WD = 7L, coefficient_vector_r = rho,
      ULM_score = unname(z["ULM_score"]), ULM_p_value = unname(z["ULM_p_value"]), slope = unname(z["slope"]),
      pairing_permutation_p = wp2_p_from_null(rho, perm_bank[, k], exact = TRUE),
      null_bootstrap_p = wp2_p_from_null(rho, null_bank[, k]),
      exact_pairing_permutations = nrow(permutations), permutation_valid_draws = pv,
      permutation_invalid_draws = nrow(permutations) - pv,
      null_bootstrap_draws = B, bootstrap_valid_draws = nv, bootstrap_invalid_draws = B - nv,
      paired_stability_draws = paired_B, paired_valid_draws = paired_valid,
      paired_positive_fraction = if (paired_valid) mean(paired_bank[is.finite(paired_bank[, k]), k] > 0) else NA_real_,
      paired_same_direction_fraction = if (paired_valid && is.finite(rho)) mean(sign(paired_bank[is.finite(paired_bank[, k]), k]) == sign(rho)) else NA_real_,
      paired_reference_valid_draws = combined_valid,
      paired_reference_positive_fraction = if (combined_valid) mean(combined_bank[is.finite(combined_bank[, k]), k] > 0) else NA_real_,
      paired_reference_same_direction_fraction = if (combined_valid && is.finite(rho)) mean(sign(combined_bank[is.finite(combined_bank[, k]), k]) == sign(rho)) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  comparison <- do.call(rbind, lapply(specs, function(s) {
    x <- contexts[[s$context]]; e <- x$edges[s$index, , drop = FALSE]
    z <- wp2_vector_stats(e$r_SD, e$r_WD)
    data.frame(context = s$context, sex = unique(e$sex), zone = unique(e$zone),
               universe = s$universe, scope = s$scope, gene = s$gene,
               n_edges = nrow(e), n_SD = 7L, n_WD = 7L,
               coefficient_vector_r_SD_WD = unname(z["coefficient_vector_r"]),
               ULM_score_SD_WD = unname(z["ULM_score"]), ULM_p_value_SD_WD = unname(z["ULM_p_value"]),
               mean_abs_delta_r = mean(e$abs_delta_r), median_abs_delta_r = stats::median(e$abs_delta_r),
               sign_change_fraction = mean(e$sign_change), stringsAsFactors = FALSE)
  }))
  stability <- wp2_edge_change_stability(contexts)
  pairing_ledger <- data.frame(draw = rep(seq_len(nrow(permutations)), each = 7L),
                              position = rep(seq_len(7L), times = nrow(permutations)),
                              lipid_source_position = as.vector(t(permutations)))
  list(edges = stability$edges, inference = inference, comparison_summary = comparison,
       permutation_draws = data.frame(draw = seq_len(nrow(perm_bank)), perm_bank, check.names = FALSE),
       bootstrap_null_draws = data.frame(draw = seq_len(B), null_bank, check.names = FALSE),
       paired_stability_draws = data.frame(draw = seq_len(paired_B), paired_bank, check.names = FALSE),
       paired_reference_stability_draws = data.frame(draw = seq_len(paired_B), combined_bank, check.names = FALSE),
       resampling_indices = do.call(rbind, draw_ledgers), pairing_permutation_indices = pairing_ledger,
       loco_matrices = stability$loco_matrices, loco_omission_ledger = stability$loco_omission_ledger,
       loco_aggregate_summary = stability$loco_aggregate_summary,
       metadata = list(seed = seed, B = B, paired_B = paired_B, workers_requested = workers,
                       workers_used = 1L, exact_permutations = factorial(7L),
                       statistic = "Pearson association between frozen edge-coefficient vectors",
                       null = "RNA-vector/lipid-vector independence within sex and standard diet",
                       nominal_ULM = "intercept OLS slope t-statistic; dependent edge residuals, descriptive only",
                       bootstrap = "independent empirical-marginal whole-mouse null bootstrap",
                       paired_bootstrap = "paired whole mice within cohort; reference groups retained; stability only",
                       diet_inference = "descriptive changes with 49 shared-zone joint mouse omissions; no diet p-value"))
}

wp2_baseline_assessment <- function(contexts, edges, bootstrap_draws = 4999L,
                                     seed = 20261007L, permutations = "exact",
                                     paired_draws = 999L, focus_genes = c("Mtmr11"),
                                     workers = 1L, verbose = TRUE) {
  if (!identical(permutations, "exact")) stop("Only the exact 7! pairing bank is supported.")
  wp2_network_analysis(contexts, B = bootstrap_draws, paired_B = paired_draws,
                       seed = seed, focus_genes = focus_genes, workers = workers,
                       edges = edges, verbose = verbose)
}
