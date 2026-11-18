# WP3: frozen, male Western-diet networks in independent CF453/CF459 mice.
# Source wdsd_network_statistics.R first. Sourcing this helper performs no
# analysis, changes no RNG state, and reads no data or results.

wp3_require_statistics <- function() {
  required <- c("wp2_validate_matrix", "wp2_prepare_edges", "wp2_prepared_cor",
                "wp2_vector_stats", "wp2_all_permutations")
  missing <- required[!vapply(required, exists, logical(1L), mode = "function", inherits = TRUE)]
  if (length(missing)) stop("Source wdsd_network_statistics.R first; missing: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

wp3_count <- function(x, label, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x > .Machine$integer.max || x < minimum || x != trunc(x)) {
    stop(label, " must be one integer >= ", minimum, ".")
  }
  as.integer(x)
}

wp3_normalize_context <- function(x, name) {
  wp3_require_statistics()
  allowed <- c("male_periportal", "male_pericentral")
  if (!name %in% allowed) stop("Unknown WP3 context: ", name, "; use ", paste(allowed, collapse = ", "), ".")
  required <- c("rna_CF453", "lipid_CF453", "rna_CF459", "lipid_CF459", "edges")
  if (!is.list(x) || !all(required %in% names(x))) stop("Context ", name, " requires ", paste(required, collapse = ", "), ".")
  for (field in required[1:4]) {
    wp2_validate_matrix(x[[field]], paste(name, field))
    if (any(!nzchar(rownames(x[[field]]))) || any(!nzchar(colnames(x[[field]])))) {
      stop("Feature and mouse IDs must not be empty in ", name, ".")
    }
  }
  if (ncol(x$rna_CF453) != 7L || ncol(x$lipid_CF453) != 7L ||
      ncol(x$rna_CF459) != 4L || ncol(x$lipid_CF459) != 4L) {
    stop("WP3 requires seven CF453 WD and four CF459 WD male mice, with matched modalities.")
  }
  if (length(intersect(colnames(x$rna_CF453), colnames(x$rna_CF459)))) {
    stop("CF453 and CF459 must contain independent, differently named mice.")
  }
  expected_zone <- sub("^male_", "", name)
  if ((!is.null(x$sex) && !identical(as.character(x$sex), "male")) ||
      (!is.null(x$zone) && !identical(tolower(as.character(x$zone)), expected_zone))) {
    stop("Context name and sex/zone annotations disagree in ", name, ".")
  }
  e <- x$edges
  if (!is.data.frame(e) || !nrow(e)) stop("edges must be a nonempty data frame in ", name, ".")
  if (!"feature_id" %in% names(e) && "lipid" %in% names(e)) e$feature_id <- as.character(e$lipid)
  if (!"r_reference" %in% names(e) && "reference_r" %in% names(e)) e$r_reference <- e$reference_r
  if (!all(c("gene", "feature_id", "r_reference", "primary_edge") %in% names(e))) {
    stop("Edges require gene, feature_id/lipid, r_reference, and primary_edge (displayed FCG membership).")
  }
  e$gene <- as.character(e$gene); e$feature_id <- as.character(e$feature_id)
  if (anyNA(e$gene) || anyNA(e$feature_id) || any(!nzchar(e$gene)) || any(!nzchar(e$feature_id))) {
    stop("Edge gene/lipid identities must be complete and nonempty.")
  }
  if (!is.numeric(e$r_reference) || !is.logical(e$primary_edge) || anyNA(e$primary_edge)) {
    stop("r_reference must be numeric and primary_edge a complete logical vector.")
  }
  if (!"edge_id" %in% names(e)) e$edge_id <- paste(e$gene, e$feature_id, sep = "::")
  e$edge_id <- as.character(e$edge_id)
  if (anyNA(e$edge_id) || any(!nzchar(e$edge_id)) || anyDuplicated(e$edge_id) ||
      anyDuplicated(paste(e$gene, e$feature_id, sep = "\r"))) stop("Duplicate or missing edge identities in ", name, ".")
  x$CF453 <- wp2_prepare_edges(x$rna_CF453, x$lipid_CF453, e)
  x$CF459 <- wp2_prepare_edges(x$rna_CF459, x$lipid_CF459, e)
  e$r_CF453_WD <- wp2_prepared_cor(x$CF453)
  e$r_CF459_WD <- wp2_prepared_cor(x$CF459)
  # The original FCG coefficient supplies provenance, not an observed outcome
  # selection rule. Freeze only the predeclared, observed-measurable common set.
  e$eligible <- is.finite(e$r_CF453_WD) & is.finite(e$r_CF459_WD)
  e$exclusion_reason <- ifelse(!is.finite(e$r_CF453_WD), "undefined_CF453_WD",
                                ifelse(!is.finite(e$r_CF459_WD), "undefined_CF459_WD", ""))
  e$context <- name; e$sex <- "male"; e$zone <- expected_zone
  e$n_CF453 <- 7L; e$n_CF459 <- 4L
  e$sign_agreement <- sign(e$r_CF453_WD) == sign(e$r_CF459_WD)
  x$edges <- e; x$sex <- "male"; x$zone <- expected_zone
  x
}

wp3_specifications <- function(contexts, focus_genes = "Mtmr11") {
  specs <- list(); skipped <- list()
  for (name in names(contexts)) {
    e <- contexts[[name]]$edges
    for (universe in c("full", "displayed")) {
      ix <- which(e$eligible & (universe == "full" | e$primary_edge))
      scopes <- list(network = ix)
      for (gene in unique(focus_genes)) scopes[[paste0("gene::", gene)]] <- ix[e$gene[ix] == gene]
      for (scope_name in names(scopes)) {
        key <- paste(name, universe, scope_name, sep = "__")
        s <- list(context = name, universe = universe,
                  scope = if (scope_name == "network") "network" else "gene",
                  gene = if (scope_name == "network") NA_character_ else sub("^gene::", "", scope_name),
                  index = scopes[[scope_name]])
        if (length(s$index) < 3L) {
          skipped[[key]] <- data.frame(specification = key, context = name, universe = universe,
                                       scope = s$scope, gene = s$gene, n_edges = length(s$index),
                                       reason = "fewer_than_three_observed_eligible_edges", stringsAsFactors = FALSE)
        } else specs[[key]] <- s
      }
    }
  }
  list(specs = specs, skipped = if (length(skipped)) do.call(rbind, skipped) else data.frame())
}

wp3_null_tail <- function(observed_r, null_r, exact = FALSE) {
  # With a fixed vector length, absolute ULM is monotone in absolute Pearson r.
  # Use r for numerical stability near perfect alignment. Undefined whole-vector
  # draws are recorded, never rescued by deleting individual edges.
  valid <- is.finite(null_r); n <- length(null_r); nv <- sum(valid)
  out <- c(p_value = NA_real_, valid_draws = nv, invalid_draws = n - nv,
           exceedances = NA_real_, completion_lower = NA_real_, completion_upper = NA_real_)
  if (!is.finite(observed_r)) return(out)
  exceed <- sum(abs(null_r[valid]) >= abs(observed_r) - 1e-12)
  plus <- if (exact) 0 else 1
  out[c("p_value", "exceedances", "completion_lower", "completion_upper")] <-
    c(if (nv) (plus + exceed) / (plus + nv) else NA_real_, exceed,
      (plus + exceed) / (plus + n), (plus + exceed + n - nv) / (plus + n))
  out
}

wp3_index_ledger <- function(index, cohort, method, ids) {
  data.frame(cohort = cohort, method = method,
             draw = rep(seq_len(nrow(index)), each = ncol(index)),
             position = rep(seq_len(ncol(index)), times = nrow(index)),
             source_position = as.vector(t(index)), source_mouse = ids[as.vector(t(index))],
             stringsAsFactors = FALSE)
}

wp3_loco <- function(contexts, specs) {
  # A global omission bank guarantees that both RNA zones use the same mice.
  joint <- expand.grid(omit_CF453 = seq_len(7L), omit_CF459 = seq_len(4L))
  joint$method <- "joint_omission"
  a <- data.frame(omit_CF453 = seq_len(7L), omit_CF459 = NA_integer_, method = "CF453_only")
  b <- data.frame(omit_CF453 = NA_integer_, omit_CF459 = seq_len(4L), method = "CF459_only")
  ledger <- rbind(joint, a, b); ledger$draw <- seq_len(nrow(ledger))
  ledger$omitted_CF453_mouse <- contexts[[1L]]$CF453$mouse_ids[ledger$omit_CF453]
  ledger$omitted_CF459_mouse <- contexts[[1L]]$CF459$mouse_ids[ledger$omit_CF459]
  index_ledger <- list(); matrices <- list(); rows <- list(); summaries <- list()
  kept <- lapply(seq_len(nrow(ledger)), function(j) {
    ia <- if (is.na(ledger$omit_CF453[j])) seq_len(7L) else setdiff(seq_len(7L), ledger$omit_CF453[j])
    ib <- if (is.na(ledger$omit_CF459[j])) seq_len(4L) else setdiff(seq_len(4L), ledger$omit_CF459[j])
    for (cohort in c("CF453", "CF459")) {
      ix <- if (cohort == "CF453") ia else ib
      index_ledger[[paste(j, cohort)]] <<- data.frame(
        draw = j, method = ledger$method[j], cohort = cohort, position = seq_along(ix),
        source_position = ix, source_mouse = contexts[[1L]][[cohort]]$mouse_ids[ix], stringsAsFactors = FALSE)
    }
    list(CF453 = ia, CF459 = ib)
  })
  for (name in names(contexts)) {
    x <- contexts[[name]]; e <- x$edges
    ca <- cb <- matrix(NA_real_, nrow(e), nrow(ledger), dimnames = list(e$edge_id, paste0("loco_", ledger$draw)))
    hits <- which(vapply(specs, function(s) s$context == name, logical(1L)))
    for (j in seq_len(nrow(ledger))) {
      ca[, j] <- wp2_prepared_cor(x$CF453, kept[[j]]$CF453)
      cb[, j] <- wp2_prepared_cor(x$CF459, kept[[j]]$CF459)
      for (k in hits) {
        s <- specs[[k]]; z <- wp2_vector_stats(ca[s$index, j], cb[s$index, j])
        rows[[paste(name, j, k)]] <- data.frame(
          specification = names(specs)[k], context = name, universe = s$universe, scope = s$scope,
          gene = s$gene, draw = j, method = ledger$method[j],
          n_CF453 = length(kept[[j]]$CF453), n_CF459 = length(kept[[j]]$CF459),
          coefficient_vector_r = unname(z["coefficient_vector_r"]), ULM_score = unname(z["ULM_score"]),
          n_edges = length(s$index), stringsAsFactors = FALSE)
      }
    }
    matrices[[name]] <- list(r_CF453_WD = ca, r_CF459_WD = cb)
    ij <- which(ledger$method == "joint_omission")
    defined <- is.finite(ca[, ij, drop = FALSE]) & is.finite(cb[, ij, drop = FALSE])
    sign_match <- sign(ca[, ij, drop = FALSE]) == sign(cb[, ij, drop = FALSE])
    e$loco_joint_valid_combinations <- rowSums(defined)
    e$loco_joint_sign_agreement_fraction <- rowSums(sign_match & defined, na.rm = TRUE) / rowSums(defined)
    e$loco_joint_sign_agreement_fraction[!rowSums(defined)] <- NA_real_
    contexts[[name]]$edges <- e
  }
  values <- do.call(rbind, rows); rownames(values) <- NULL
  for (key in names(specs)) {
    s <- specs[[key]]; e <- contexts[[s$context]]$edges
    observed <- wp2_vector_stats(e$r_CF453_WD[s$index], e$r_CF459_WD[s$index])["coefficient_vector_r"]
    for (method in unique(ledger$method)) {
      v <- values$coefficient_vector_r[values$specification == key & values$method == method]
      valid <- is.finite(v)
      summaries[[paste(key, method)]] <- data.frame(
        specification = key, context = s$context, universe = s$universe, scope = s$scope,
        gene = s$gene, method = method, draws = length(v), valid_draws = sum(valid),
        invalid_draws = sum(!valid), minimum_r = if (any(valid)) min(v[valid]) else NA_real_,
        maximum_r = if (any(valid)) max(v[valid]) else NA_real_,
        positive_fraction = if (any(valid)) mean(v[valid] > 0) else NA_real_,
        same_direction_fraction = if (any(valid) && is.finite(observed)) mean(sign(v[valid]) == sign(observed)) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  list(edges = do.call(rbind, lapply(contexts, `[[`, "edges")), matrices = matrices,
       omission_ledger = ledger, indices = do.call(rbind, index_ledger),
       draws = values, summary = do.call(rbind, summaries))
}

#' Cross-cohort conservation of frozen Western-diet correlation networks.
#'
#' CF453 and CF459 are independent cohorts, each with paired RNA/lipid mice.
#' The primary exact test conditions on observed CF453 coefficients and permutes
#' the four complete CF459 lipid mouse vectors (24 possibilities). The reciprocal
#' 7! test conditions on CF459 and is a complementary perspective, not a second
#' replication. These are pairing-independence nulls, not tests of equality of
#' cohort correlations or of zero alignment allowing arbitrary other coupling.
#' Nominal intercept-ULM p-values treat dependent edges as independent residuals
#' and are descriptive. Paired bootstraps and omissions measure stability only.
wp3_crosscohort_analysis <- function(contexts, B = 4999L, paired_B = 999L,
                                    seed = 20261014L, focus_genes = "Mtmr11",
                                    verbose = TRUE) {
  wp3_require_statistics()
  B <- wp3_count(B, "B"); paired_B <- wp3_count(paired_B, "paired_B")
  seed <- wp3_count(seed, "seed", minimum = 0L)
  if (!is.list(contexts) || !length(contexts) || is.null(names(contexts)) ||
      anyNA(names(contexts)) || any(!nzchar(names(contexts))) || anyDuplicated(names(contexts))) {
    stop("contexts must be a nonempty, uniquely named list.")
  }
  if (!is.character(focus_genes) || anyNA(focus_genes) || any(!nzchar(focus_genes))) {
    stop("focus_genes must be complete nonempty gene names, or character(0).")
  }
  contexts <- stats::setNames(lapply(names(contexts), function(name) wp3_normalize_context(contexts[[name]], name)), names(contexts))
  for (cohort in c("CF453", "CF459")) {
    ids <- contexts[[1L]][[cohort]]$mouse_ids
    if (any(!vapply(contexts, function(x) identical(ids, x[[cohort]]$mouse_ids), logical(1L)))) {
      stop("Both zones must have the same ordered mouse IDs within ", cohort, ".")
    }
  }
  specification <- wp3_specifications(contexts, unique(focus_genes)); specs <- specification$specs
  if (!length(specs)) stop("No frozen eligible vector contains at least three edges.")
  labels <- names(specs)
  new_bank <- function(n) matrix(NA_real_, n, length(specs), dimnames = list(NULL, labels))
  perm459 <- wp2_all_permutations(4L); perm453 <- wp2_all_permutations(7L)
  banks <- list(permutation_CF459 = new_bank(24L), permutation_CF453 = new_bank(5040L),
                bootstrap_null = new_bank(B), paired_stability = new_bank(paired_B))
  score_banks <- banks
  set.seed(seed)
  indices <- list(
    CF459_null_rna = matrix(sample.int(4L, 4L * B, replace = TRUE), B, 4L),
    CF459_null_lipid = matrix(sample.int(4L, 4L * B, replace = TRUE), B, 4L),
    CF453_paired = matrix(sample.int(7L, 7L * paired_B, replace = TRUE), paired_B, 7L),
    CF459_paired = matrix(sample.int(4L, 4L * paired_B, replace = TRUE), paired_B, 4L))
  index_ledger <- do.call(rbind, lapply(names(indices), function(method) {
    cohort <- sub("_.*$", "", method)
    wp3_index_ledger(indices[[method]], cohort, method, contexts[[1L]][[cohort]]$mouse_ids)
  }))
  draw_summary <- do.call(rbind, lapply(names(indices), function(method) data.frame(
    method = method, cohort = sub("_.*$", "", method), draw = seq_len(nrow(indices[[method]])),
    n_positions = ncol(indices[[method]]),
    n_unique_mice = apply(indices[[method]], 1L, function(z) length(unique(z))), stringsAsFactors = FALSE)))
  observed <- lapply(specs, function(s) {
    e <- contexts[[s$context]]$edges
    wp2_vector_stats(e$r_CF453_WD[s$index], e$r_CF459_WD[s$index])
  })
  for (name in names(contexts)) {
    x <- contexts[[name]]
    hits <- which(vapply(specs, function(s) s$context == name, logical(1L)))
    if (verbose) message("WP3 resampling ", name, ": ", sum(x$edges$eligible), "/", nrow(x$edges), " frozen eligible edges; CF453 n=7, CF459 n=4.")
    calculate <- function(a, b) {
      vapply(hits, function(k) {
        ix <- specs[[k]]$index
        wp2_vector_stats(a[ix], b[ix])[c("coefficient_vector_r", "ULM_score")]
      }, numeric(2L))
    }
    for (j in seq_len(24L)) {
      z <- calculate(x$edges$r_CF453_WD, wp2_prepared_cor(x$CF459, seq_len(4L), perm459[j, ], fixed_margins = TRUE))
      banks$permutation_CF459[j, hits] <- z[1L, ]; score_banks$permutation_CF459[j, hits] <- z[2L, ]
    }
    for (j in seq_len(5040L)) {
      z <- calculate(wp2_prepared_cor(x$CF453, seq_len(7L), perm453[j, ], fixed_margins = TRUE), x$edges$r_CF459_WD)
      banks$permutation_CF453[j, hits] <- z[1L, ]; score_banks$permutation_CF453[j, hits] <- z[2L, ]
    }
    for (j in seq_len(B)) {
      z <- calculate(x$edges$r_CF453_WD, wp2_prepared_cor(x$CF459, indices$CF459_null_rna[j, ], indices$CF459_null_lipid[j, ]))
      banks$bootstrap_null[j, hits] <- z[1L, ]; score_banks$bootstrap_null[j, hits] <- z[2L, ]
    }
    for (j in seq_len(paired_B)) {
      z <- calculate(wp2_prepared_cor(x$CF453, indices$CF453_paired[j, ]), wp2_prepared_cor(x$CF459, indices$CF459_paired[j, ]))
      banks$paired_stability[j, hits] <- z[1L, ]; score_banks$paired_stability[j, hits] <- z[2L, ]
    }
  }
  inference <- do.call(rbind, lapply(seq_along(specs), function(k) {
    s <- specs[[k]]; z <- observed[[k]]; rho <- unname(z["coefficient_vector_r"])
    p459 <- wp3_null_tail(rho, banks$permutation_CF459[, k], exact = TRUE)
    p453 <- wp3_null_tail(rho, banks$permutation_CF453[, k], exact = TRUE)
    nb <- wp3_null_tail(rho, banks$bootstrap_null[, k])
    pv <- is.finite(banks$paired_stability[, k]); e <- contexts[[s$context]]$edges[s$index, , drop = FALSE]
    data.frame(specification = labels[k], context = s$context, sex = "male", zone = unique(e$zone),
               universe = s$universe, analysis_role = if (s$universe == "full") "primary" else "secondary_displayed",
               scope = s$scope, gene = s$gene, n_edges = length(s$index), n_genes = length(unique(e$gene)),
               n_lipids = length(unique(e$feature_id)), n_CF453 = 7L, n_CF459 = 4L,
               coefficient_vector_r = rho, ULM_score = unname(z["ULM_score"]),
               ULM_p_value = unname(z["ULM_p_value"]), slope_CF459_on_CF453 = unname(z["slope"]),
               sign_agreement_fraction = mean(e$sign_agreement),
               pairing_CF459_given_CF453_p = unname(p459["p_value"]),
               pairing_CF453_given_CF459_p = unname(p453["p_value"]),
               null_bootstrap_p = unname(nb["p_value"]),
               permutation_CF459_draws = 24L, permutation_CF459_valid_draws = unname(p459["valid_draws"]),
               permutation_CF459_invalid_draws = unname(p459["invalid_draws"]),
               permutation_CF453_draws = 5040L, permutation_CF453_valid_draws = unname(p453["valid_draws"]),
               permutation_CF453_invalid_draws = unname(p453["invalid_draws"]),
               null_bootstrap_draws = B, bootstrap_valid_draws = unname(nb["valid_draws"]),
               bootstrap_invalid_draws = unname(nb["invalid_draws"]),
               bootstrap_completion_lower = unname(nb["completion_lower"]),
               bootstrap_completion_upper = unname(nb["completion_upper"]),
               paired_stability_draws = paired_B, paired_valid_draws = sum(pv), paired_invalid_draws = sum(!pv),
               paired_positive_fraction = if (any(pv)) mean(banks$paired_stability[pv, k] > 0) else NA_real_,
               paired_same_direction_fraction = if (any(pv) && is.finite(rho)) mean(sign(banks$paired_stability[pv, k]) == sign(rho)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  loco <- wp3_loco(contexts, specs)
  as_draw_frame <- function(x) data.frame(draw = seq_len(nrow(x)), x, check.names = FALSE)
  result <- list(edges = loco$edges, inference = inference, skipped_specifications = specification$skipped,
                 resampling_indices = index_ledger, resampling_draw_summary = draw_summary,
                 resampling_index_matrices = indices,
                 pairing_permutation_indices = rbind(
                   wp3_index_ledger(perm459, "CF459", "exact_pairing", contexts[[1L]]$CF459$mouse_ids),
                   wp3_index_ledger(perm453, "CF453", "exact_pairing", contexts[[1L]]$CF453$mouse_ids)),
                 loco_matrices = loco$matrices, loco_omission_ledger = loco$omission_ledger,
                 loco_indices = loco$indices, loco_draws = loco$draws, loco_summary = loco$summary)
  for (method in names(banks)) {
    result[[paste0(method, "_draws")]] <- as_draw_frame(banks[[method]])
    result[[paste0(method, "_ULM_draws")]] <- as_draw_frame(score_banks[[method]])
  }
  result$metadata <- list(seed = seed, B = B, paired_B = paired_B,
    statistic = "Pearson association between frozen CF453 WD and CF459 WD edge-coefficient vectors",
    primary_universe = "all observed-measurable common frozen FCG-context edges; no CF453/CF459 coefficient threshold",
    secondary_universe = "displayed FCG edges within the same observed-measurable common set",
    exact_primary = "24 complete CF459 lipid-vector pairing permutations, conditional on observed CF453 coefficients",
    exact_complementary = "5040 complete CF453 lipid-vector pairing permutations, conditional on observed CF459 coefficients",
    exact_identity_included = TRUE, null_statistic = "absolute coefficient-vector r, equivalent to absolute ULM at fixed edge count",
    null_hypothesis = "RNA-vector/lipid-vector independence in the permuted or null-resampled cohort; not coefficient equality",
    bootstrap = "independent empirical-marginal CF459 RNA/lipid whole-mouse vectors; fixed observed CF453 anchor; approximate null",
    p_value_convention = "exact exceedances/valid enumeration; Monte Carlo (1+exceedances)/(1+valid draws), including ties",
    undefined_draws = "undefined edge invalidates the entire frozen vector; no draw-specific edge deletion; P conditional on defined draws",
    completion_bounds = "allocate all undefined null draws to either nonexceedance or exceedance; mathematical bounds, not confidence intervals or calibrated corrections",
    nominal_ULM = "intercept OLS slope t-statistic; dependent edge residuals make nominal P descriptive",
    paired_bootstrap = "paired whole mice sampled independently within each cohort; both zones share indices; stability only",
    loco = "28 joint omissions plus seven CF453-only and four CF459-only omissions; same mice in both zones; descriptive stability",
    limitations = "CF459 n=4; neither cohort age/background equivalence nor within-CF459 SD versus WD effects are identified")
  result
}
