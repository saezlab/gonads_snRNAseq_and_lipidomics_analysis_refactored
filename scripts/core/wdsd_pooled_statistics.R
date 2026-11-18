# Pooled whole-hepatocyte Figure 3 conservation, using whole mice as units.
# Source scripts/core/wdsd_network_statistics.R before calling these helpers.
# Sourcing this file defines functions only and does not change the RNG state.

wp2_pooled_center_sex <- function(x, sex) {
  wp2_validate_matrix(x, "Pooled centering matrix")
  if (length(sex) != ncol(x) || anyNA(sex) || any(!nzchar(as.character(sex)))) {
    stop("Provide one complete recorded sex label for every ordered mouse.")
  }
  answer <- x
  for (ix in split(seq_len(ncol(x)), as.character(sex))) {
    answer[, ix] <- x[, ix, drop = FALSE] - rowMeans(x[, ix, drop = FALSE])
  }
  answer
}

wp2_pooled_null_tails <- function(observed, null, tolerance = 1e-12) {
  valid <- is.finite(null)
  n <- sum(valid)
  result <- c(absolute_p = NA_real_, upper_p = NA_real_, lower_p = NA_real_,
              equal_tail_p = NA_real_, valid_draws = n,
              invalid_draws = length(null) - n, null_median = NA_real_)
  if (!is.finite(observed) || n == 0L) return(result)
  z <- null[valid]
  # The observed statistic supplies the +1 term, independently of the saved
  # identity audit row. Random null draws do not contain that extra audit row.
  upper <- (1 + sum(z >= observed - tolerance)) / (1 + n)
  lower <- (1 + sum(z <= observed + tolerance)) / (1 + n)
  c(absolute_p = (1 + sum(abs(z) >= abs(observed) - tolerance)) / (1 + n),
    upper_p = upper, lower_p = lower, equal_tail_p = min(1, 2 * min(upper, lower)),
    valid_draws = n, invalid_draws = length(null) - n,
    null_median = stats::median(z))
}

wp2_pooled_index_bank <- function(n, strata, draws, replace = FALSE) {
  if (length(strata) != n || anyNA(strata) || draws < 1L) stop("Invalid draw-bank strata/count.")
  groups <- split(seq_len(n), as.character(strata))
  t(vapply(seq_len(draws), function(j) {
    index <- seq_len(n)
    for (ix in groups) index[ix] <- ix[sample.int(length(ix), length(ix), replace = replace)]
    index
  }, integer(n)))
}

wp2_pooled_context <- function(context) {
  required <- c("rna_SD", "lipid_SD", "metadata_SD", "edges", "reference_rna",
                "reference_lipid", "reference_metadata")
  if (!all(required %in% names(context))) stop("Pooled context requires ", paste(required, collapse = ", "), ".")
  e <- context$edges
  if (!all(c("gene", "feature_id", "r_reference", "primary_edge") %in% names(e)) ||
      !is.logical(e$primary_edge) || anyNA(e$primary_edge)) {
    stop("Pooled edges require gene, feature_id, r_reference and complete logical primary_edge.")
  }
  if (!"edge_id" %in% names(e)) e$edge_id <- paste(e$gene, e$feature_id, sep = "::")
  if (anyDuplicated(e$edge_id) || anyDuplicated(paste(e$gene, e$feature_id, sep = "\r"))) {
    stop("Pooled edge identities must be unique.")
  }
  context$SD <- wp2_prepare_edges(context$rna_SD, context$lipid_SD, e)
  context$REF <- wp2_prepare_edges(context$reference_rna, context$reference_lipid, e)
  for (field in c("SD", "REF")) {
    m <- context[[if (field == "SD") "metadata_SD" else "reference_metadata"]]
    if (!all(c("mouse_id", "sex") %in% names(m)) ||
        !identical(as.character(m$mouse_id), context[[field]]$mouse_ids) ||
        anyNA(m$sex) || !all(as.character(m$sex) %in% c("female", "male"))) {
      stop(field, " metadata must give the identical ordered mouse_id and female/male labels.")
    }
    if (anyDuplicated(m$mouse_id) || length(unique(as.character(m$sex))) != 2L) {
      stop(field, " must have distinct mice in both recorded sexes.")
    }
  }
  sd_sex <- as.character(context$metadata_SD$sex)
  if (ncol(context$rna_SD) != 14L || !identical(as.integer(table(factor(sd_sex, c("female", "male")))), c(7L, 7L))) {
    stop("Pooled WDSD SD contract requires 14 mice, seven female and seven male.")
  }
  ref_meta <- context$reference_metadata
  if ("group" %in% names(ref_meta)) {
    if (anyNA(ref_meta$group) || any(!nzchar(as.character(ref_meta$group))) ||
        any(vapply(split(as.character(ref_meta$sex), as.character(ref_meta$group)),
                   function(z) length(unique(z)) != 1L, logical(1L)))) {
      stop("FCG reference groups must be complete and nested within recorded sex.")
    }
  }
  e$r_SD <- wp2_prepared_cor(context$SD)
  replay <- wp2_prepared_cor(context$REF)
  checked <- is.finite(e$r_reference)
  if (any(!is.finite(replay[checked])) ||
      (any(checked) && max(abs(replay[checked] - e$r_reference[checked])) > 1e-8)) {
    stop("Supplied reference matrices do not replay frozen Figure 3 coefficients.")
  }
  e$eligible <- is.finite(e$r_reference) & is.finite(e$r_SD)
  e$exclusion_reason <- ifelse(!is.finite(e$r_reference), "undefined_reference",
                              ifelse(!is.finite(e$r_SD), "undefined_SD", ""))
  rna_SD_centered <- wp2_pooled_center_sex(context$rna_SD, sd_sex)
  lipid_SD_centered <- wp2_pooled_center_sex(context$lipid_SD, sd_sex)
  rna_REF_centered <- wp2_pooled_center_sex(context$reference_rna, as.character(ref_meta$sex))
  lipid_REF_centered <- wp2_pooled_center_sex(context$reference_lipid, as.character(ref_meta$sex))
  e$r_SD_withinsex <- wp2_cor_edges(rna_SD_centered, lipid_SD_centered, e)
  e$r_reference_withinsex <- wp2_cor_edges(rna_REF_centered, lipid_REF_centered, e)
  # Residual-invariant features stay visible as undefined on the same frozen
  # edge identities. They do not create a smaller, newly selected universe.
  e$withinsex_defined <- is.finite(e$r_SD_withinsex) & is.finite(e$r_reference_withinsex)
  context$edges <- e
  context
}

wp2_pooled_specifications <- function(context, focus_genes) {
  e <- context$edges
  specifications <- list()
  for (universe in c("full", "primary")) {
    ix <- which(e$eligible & (universe == "full" | e$primary_edge))
    if (length(ix) < 3L) next
    scopes <- list(network = ix)
    for (gene in intersect(unique(as.character(focus_genes)), unique(as.character(e$gene[ix])))) {
      gi <- ix[e$gene[ix] == gene]
      if (length(gi) >= 3L) scopes[[paste0("gene::", gene)]] <- gi
    }
    for (scope in names(scopes)) {
      key <- paste(universe, scope, sep = "__")
      specifications[[key]] <- list(universe = universe, scope = if (scope == "network") "network" else "gene",
                                    gene = if (scope == "network") NA_character_ else sub("^gene::", "", scope),
                                    index = scopes[[scope]])
    }
  }
  specifications
}

#' Pooled Figure 3 replication with distinct overall and conditional nulls.
#'
#' The broad null independently rearranges or resamples RNA and lipid mouse
#' vectors across all 14 SD animals. Rejection may result from shared gonadal
#' sex structure; it does not isolate within-sex correspondence. The conditional
#' null independently rearranges/resamples vectors within sex, retaining the
#' between-sex margins. Its prespecified upper tail asks whether positive FCG
#' alignment exceeds that produced by the shared sex structure alone. Each draw
#' is shared by every feature and edge. Monte Carlo p-values use +1 correction.
#'
#' Native ULM p-values treat dependent edge residuals as independent and are
#' descriptive. Paired within-sex/FCG-group bootstraps describe stability and
#' do not supply null p-values. Undefined edges invalidate their entire frozen
#' summary vector; probabilities explicitly condition on valid draws. Bootstrap
#' calibration at n=7 per sex is not established. Sex-centered correlations are
#' a different, descriptive network sensitivity, not archival Figure 3 replay.
wp2_pooled_analysis <- function(context, permutations = 19999L, B = 4999L,
                                paired_B = 999L, seed = 20261007L,
                                focus_genes = c("Cyp2b9", "Cyp2b13", "Cyp2c38", "Cyp2d9"),
                                unrestricted = TRUE, verbose = TRUE) {
  for (count in list(permutations, B, paired_B)) {
    if (length(count) != 1L || is.na(count) || !is.finite(count) || count < 1L || count != as.integer(count)) {
      stop("Permutation, null-bootstrap and paired-bootstrap counts must be positive integers.")
    }
  }
  context <- wp2_pooled_context(context)
  e <- context$edges
  specs <- wp2_pooled_specifications(context, focus_genes)
  if (!length(specs)) stop("No eligible pooled fixed coefficient vector with at least three edges.")
  labels <- names(specs)
  observed <- lapply(specs, function(s) wp2_vector_stats(e$r_reference[s$index], e$r_SD[s$index]))
  centered <- lapply(specs, function(s) wp2_vector_stats(e$r_reference_withinsex[s$index], e$r_SD_withinsex[s$index]))
  set.seed(seed)
  sd_sex <- as.character(context$metadata_SD$sex)
  reference_strata <- if ("group" %in% names(context$reference_metadata)) {
    as.character(context$reference_metadata$group)
  } else as.character(context$reference_metadata$sex)
  n_reference <- ncol(context$reference_rna)
  banks <- list(
    permutation_withinsex = wp2_pooled_index_bank(14L, sd_sex, permutations),
    bootstrap_withinsex_rna = wp2_pooled_index_bank(14L, sd_sex, B, replace = TRUE),
    bootstrap_withinsex_lipid = wp2_pooled_index_bank(14L, sd_sex, B, replace = TRUE),
    paired_SD = wp2_pooled_index_bank(14L, sd_sex, paired_B, replace = TRUE),
    paired_reference = wp2_pooled_index_bank(n_reference, reference_strata, paired_B, replace = TRUE))
  if (isTRUE(unrestricted)) {
    banks$permutation_unrestricted <- wp2_pooled_index_bank(14L, rep("all", 14L), permutations)
    banks$bootstrap_unrestricted_rna <- wp2_pooled_index_bank(14L, rep("all", 14L), B, replace = TRUE)
    banks$bootstrap_unrestricted_lipid <- wp2_pooled_index_bank(14L, rep("all", 14L), B, replace = TRUE)
  }
  calculate <- function(values, reference = e$r_reference) {
    vapply(specs, function(s) unname(wp2_vector_stats(reference[s$index], values[s$index])["coefficient_vector_r"]), numeric(1L))
  }
  # Identity is recorded separately from the random draws, avoiding double
  # counting its +1 contribution. Both permutation banks have an audit row 0.
  identity <- calculate(wp2_prepared_cor(context$SD, seq_len(14L), seq_len(14L), fixed_margins = TRUE))
  if (!isTRUE(all.equal(unname(identity), vapply(observed, function(z) unname(z["coefficient_vector_r"]), numeric(1L)),
                        tolerance = 1e-12, check.attributes = FALSE))) stop("Identity pairing failed to reproduce pooled observation.")
  within_perm <- within_boot <- pooled_perm <- pooled_boot <- paired <- paired_reference <- NULL
  empty_bank <- function(n) matrix(NA_real_, n, length(specs), dimnames = list(NULL, labels))
  within_perm <- empty_bank(permutations); within_boot <- empty_bank(B)
  paired <- paired_reference <- empty_bank(paired_B)
  if (isTRUE(unrestricted)) { pooled_perm <- empty_bank(permutations); pooled_boot <- empty_bank(B) }
  if (verbose) message("Pooled Figure 3: ", nrow(e), " candidates; ", sum(e$eligible),
                       " eligible; ", length(specs), " network/gene summaries; 14 SD / ", n_reference, " FCG mice.")
  for (j in seq_len(permutations)) {
    within_perm[j, ] <- calculate(wp2_prepared_cor(context$SD, seq_len(14L), banks$permutation_withinsex[j, ], fixed_margins = TRUE))
    if (isTRUE(unrestricted)) pooled_perm[j, ] <- calculate(wp2_prepared_cor(context$SD, seq_len(14L), banks$permutation_unrestricted[j, ], fixed_margins = TRUE))
  }
  if (verbose) message("Pooled whole-mouse permutation banks completed.")
  for (j in seq_len(B)) {
    within_boot[j, ] <- calculate(wp2_prepared_cor(context$SD, banks$bootstrap_withinsex_rna[j, ], banks$bootstrap_withinsex_lipid[j, ]))
    if (isTRUE(unrestricted)) pooled_boot[j, ] <- calculate(wp2_prepared_cor(context$SD, banks$bootstrap_unrestricted_rna[j, ], banks$bootstrap_unrestricted_lipid[j, ]))
  }
  for (j in seq_len(paired_B)) {
    val <- wp2_prepared_cor(context$SD, banks$paired_SD[j, ])
    paired[j, ] <- calculate(val)
    reference <- wp2_prepared_cor(context$REF, banks$paired_reference[j, ])
    paired_reference[j, ] <- calculate(val, reference)
  }
  make_inference <- function(k, method = "pooled") {
    s <- specs[[k]]
    z <- if (method == "pooled") observed[[k]] else centered[[k]]
    answer <- data.frame(specification = labels[k], analysis = method, universe = s$universe,
                         scope = s$scope, gene = s$gene, n_edges = length(s$index),
                         n_SD = 14L, n_SD_female = 7L, n_SD_male = 7L,
                         n_reference = n_reference,
                         n_reference_female = sum(context$reference_metadata$sex == "female"),
                         n_reference_male = sum(context$reference_metadata$sex == "male"),
                         coefficient_vector_r = unname(z["coefficient_vector_r"]),
                         ULM_score = unname(z["ULM_score"]), ULM_p_value = unname(z["ULM_p_value"]),
                         slope = unname(z["slope"]),
                         observed_defined_edges = if (method == "pooled") length(s$index) else sum(e$withinsex_defined[s$index]),
                         stringsAsFactors = FALSE)
    for (bank in c("withinsex_permutation", "withinsex_bootstrap", "unrestricted_permutation", "unrestricted_bootstrap")) {
      null <- switch(bank, withinsex_permutation = within_perm[, k], withinsex_bootstrap = within_boot[, k],
                     unrestricted_permutation = if (isTRUE(unrestricted)) pooled_perm[, k] else numeric(),
                     unrestricted_bootstrap = if (isTRUE(unrestricted)) pooled_boot[, k] else numeric())
      tails <- if (method == "pooled") wp2_pooled_null_tails(answer$coefficient_vector_r, null) else {
        stats::setNames(rep(NA_real_, 7L), c("absolute_p", "upper_p", "lower_p", "equal_tail_p", "valid_draws", "invalid_draws", "null_median"))
      }
      for (field in names(tails)) answer[[paste(bank, field, sep = "_")]] <- unname(tails[field])
    }
    for (bank in c("paired", "paired_reference")) {
      null <- if (bank == "paired") paired[, k] else paired_reference[, k]
      valid <- is.finite(null)
      answer[[paste0(bank, "_valid_draws")]] <- if (method == "pooled") sum(valid) else NA_integer_
      answer[[paste0(bank, "_invalid_draws")]] <- if (method == "pooled") paired_B - sum(valid) else NA_integer_
      answer[[paste0(bank, "_positive_fraction")]] <- if (method == "pooled" && any(valid)) mean(null[valid] > 0) else NA_real_
      answer[[paste0(bank, "_same_direction_fraction")]] <- if (method == "pooled" && any(valid) && is.finite(answer$coefficient_vector_r)) {
        mean(sign(null[valid]) == sign(answer$coefficient_vector_r))
      } else NA_real_
    }
    answer
  }
  inference <- do.call(rbind, c(lapply(seq_along(specs), make_inference),
                                lapply(seq_along(specs), make_inference, method = "withinsex_centered")))
  draw_frame <- function(z, identity_row = FALSE) {
    if (is.null(z)) return(NULL)
    if (identity_row) z <- rbind(identity, z)
    data.frame(draw = if (identity_row) 0:permutations else seq_len(nrow(z)), z, check.names = FALSE)
  }
  draw_ledger <- lapply(names(banks), function(method) {
    ix <- banks[[method]]
    source_ids <- if (method == "paired_reference") context$REF$mouse_ids else context$SD$mouse_ids
    source_sex <- if (method == "paired_reference") as.character(context$reference_metadata$sex) else sd_sex
    position <- rep(seq_len(ncol(ix)), times = nrow(ix))
    source_position <- as.vector(t(ix))
    data.frame(method = method, draw = rep(seq_len(nrow(ix)), each = ncol(ix)),
               position = position, source_position = source_position,
               position_mouse = source_ids[position], source_mouse = source_ids[source_position],
               position_sex = source_sex[position], source_sex = source_sex[source_position], stringsAsFactors = FALSE)
  })
  list(edges = e, inference = inference,
       permutation_withinsex_draws = draw_frame(within_perm, identity_row = TRUE),
       permutation_unrestricted_draws = draw_frame(pooled_perm, identity_row = TRUE),
       bootstrap_withinsex_draws = draw_frame(within_boot),
       bootstrap_unrestricted_draws = draw_frame(pooled_boot),
       paired_stability_draws = draw_frame(paired),
       paired_reference_stability_draws = draw_frame(paired_reference),
       resampling_indices = do.call(rbind, draw_ledger), index_banks = banks,
       metadata = list(seed = seed, random_permutations = permutations,
                       identity_audit_draw = 0L, identity_in_null_draw_count = FALSE,
                       null_bootstrap_draws = B, paired_stability_draws = paired_B,
                       broad_null = "RNA-vector/lipid-vector independence across all 14 SD mice; sex structure may account for rejection",
                       broad_report_p = "two-sided absolute coefficient-vector association tail with +1 correction",
                       conditional_null = "RNA-vector/lipid-vector independence within each SD sex, retaining between-sex margins",
                       conditional_report_p = "one-sided upper tail for positive reference alignment beyond shared sex structure, with +1 correction",
                       paired_bootstrap = "paired SD mice within sex; paired reference mice within FCG group when available; stability only",
                       fixed_universe = "observed eligible edges fixed before resampling; one undefined edge invalidates its summary vector",
                       inference_limit = "conditional on fixed FCG reference and exchangeable mice; invalid-draw probabilities condition on valid draws; bootstrap size/power not calibrated",
                       centered_sensitivity = "featurewise subtraction of each cohort's sex means, then Pearson; descriptive r/native ULM only on unchanged edge identities",
                       nominal_ULM = "intercept OLS slope t and nominal p; dependent edges are not independent animals"))
}
