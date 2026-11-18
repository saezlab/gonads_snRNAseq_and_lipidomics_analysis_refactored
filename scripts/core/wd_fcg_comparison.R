# WP3 extension: CF459 WD closeness to original male FCG SD versus CF453 WD.
# Frozen edges and original WP3 index banks are retained. Source loading defines
# functions only; the three-cohort bootstrap measures stability, not a null P.

wp3_fcg_indices <- function(metadata, paired_B, seed) {
  if (!is.data.frame(metadata) || !all(c("mouse_id", "group") %in% names(metadata)) ||
      nrow(metadata) != 22L || anyNA(metadata[, c("mouse_id", "group")]) ||
      anyDuplicated(metadata$mouse_id) || !setequal(metadata$group, c("XXT", "XYT")) ||
      any(table(metadata$group) != 11L)) stop("Original male FCG metadata must contain 11 XXT and 11 XYT paired mice.")
  paired_B <- wp3_count(paired_B, "paired_B"); seed <- wp3_count(seed, "seed", 0L)
  # Local RNG scope leaves the caller and original WP3 banks unchanged.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  FCG_seed <- as.integer((as.double(seed) + 104729) %% .Machine$integer.max)
  set.seed(FCG_seed)
  ans <- matrix(NA_integer_, paired_B, 22L)
  for (group in c("XXT", "XYT")) {
    positions <- which(metadata$group == group)
    ans[, positions] <- matrix(positions[sample.int(length(positions), length(positions) * paired_B, replace = TRUE)], paired_B, length(positions))
  }
  ans
}

# A source-identity contingency with exact empirical independence certifies
# zero covariance for every pair of feature functions. This integer identity
# detects balanced n=4 cancellations without thresholding legitimate weak r.
wp3_fcg_exact_zero_pairing <- function(rna_indices, lipid_indices) {
  if (length(rna_indices) != length(lipid_indices) || !length(rna_indices) ||
      anyNA(rna_indices) || anyNA(lipid_indices)) stop("Complete equal-length mouse indices are required.")
  joint <- table(factor(rna_indices), factor(lipid_indices))
  all(length(rna_indices) * joint == outer(rowSums(joint), colSums(joint)))
}

wp3_fcg_triplet_stats <- function(fcg, cf453, cf459) {
  a <- wp2_vector_stats(fcg, cf459)
  b <- wp2_vector_stats(cf453, cf459)
  c <- wp2_vector_stats(fcg, cf453)
  c(r_FCG_CF459 = unname(a["coefficient_vector_r"]),
    r_CF453_CF459 = unname(b["coefficient_vector_r"]),
    r_FCG_CF453 = unname(c["coefficient_vector_r"]),
    delta_r = unname(a["coefficient_vector_r"] - b["coefficient_vector_r"]),
    ULM_FCG_CF459 = unname(a["ULM_score"]), ULM_p_FCG_CF459 = unname(a["ULM_p_value"]))
}

#' Compare CF459 with both age-labelled anchors on precisely the same edges.
#' @param contexts Original WP3 male contexts, unchanged.
#' @param reference Fresh wp2_load_reference() output with exact archived replay.
#' @param analysis Existing wp3_crosscohort_analysis() result (banks are reused).
#' @return Original-aligned edge coefficients, inference, null/paired banks,
#'   stratified FCG indices, omission diagnostics and explicit interpretations.
wp3_fcg_comparison <- function(contexts, reference, analysis, seed = 20261014L,
                               focus_genes = "Mtmr11", verbose = TRUE) {
  wp3_require_statistics()
  expected <- c("male_periportal", "male_pericentral")
  if (!is.list(contexts) || !identical(names(contexts), expected) ||
      !is.list(reference) || !all(expected %in% names(reference))) stop("Both ordered male WP3 contexts and original references are required.")
  required <- c("edges", "inference", "resampling_index_matrices", "permutation_CF459_draws", "bootstrap_null_draws", "paired_stability_draws", "loco_omission_ledger")
  if (!is.list(analysis) || !all(required %in% names(analysis))) stop("Supply the complete original WP3 analysis and saved index banks.")
  contexts <- stats::setNames(lapply(expected, function(name) wp3_normalize_context(contexts[[name]], name)), expected)
  specset <- wp3_specifications(contexts, focus_genes); specs <- specset$specs; labels <- names(specs)
  if (!length(specs) || !all(labels %in% analysis$inference$specification)) stop("Frozen specifications do not match the original WP3 analysis.")
  indices <- analysis$resampling_index_matrices
  if (!all(c("CF459_null_rna", "CF459_null_lipid", "CF453_paired", "CF459_paired") %in% names(indices))) stop("Missing original WP3 mouse banks.")
  validate_index <- function(x, n, label) {
    if (!is.matrix(x) || !is.numeric(x) || ncol(x) != n || !nrow(x) || anyNA(x) || any(x != trunc(x)) || any(x < 1L | x > n)) stop("Invalid original mouse-index bank: ", label)
  }
  for (field in names(indices)) validate_index(indices[[field]], if (grepl("^CF459", field)) 4L else 7L, field)
  B <- nrow(indices$CF459_null_rna); paired_B <- nrow(indices$CF459_paired)
  if (nrow(indices$CF459_null_lipid) != B || nrow(indices$CF453_paired) != paired_B ||
      nrow(analysis$bootstrap_null_draws) != B || nrow(analysis$paired_stability_draws) != paired_B ||
      nrow(analysis$permutation_CF459_draws) != 24L) stop("Original WP3 bank dimensions disagree.")
  first_meta <- reference[[expected[1L]]]$metadata
  FCG_paired <- wp3_fcg_indices(first_meta, paired_B, seed)
  index_ledger <- wp3_index_ledger(FCG_paired, "FCG", "FCG_stratified_paired", first_meta$mouse_id)
  index_ledger$source_group <- first_meta$group[index_ledger$source_position]
  index_ledger$position_group <- first_meta$group[index_ledger$position]
  if (any(index_ledger$source_group != index_ledger$position_group)) stop("FCG resampling must preserve chromosome strata.")
  unique_counts <- data.frame(draw = seq_len(paired_B),
    n_unique_FCG = apply(FCG_paired, 1L, function(x) length(unique(x))),
    n_unique_FCG_XXT = apply(FCG_paired[, first_meta$group == "XXT", drop = FALSE], 1L, function(x) length(unique(x))),
    n_unique_FCG_XYT = apply(FCG_paired[, first_meta$group == "XYT", drop = FALSE], 1L, function(x) length(unique(x))),
    n_unique_CF453 = apply(indices$CF453_paired, 1L, function(x) length(unique(x))),
    n_unique_CF459 = apply(indices$CF459_paired, 1L, function(x) length(unique(x))))
  certified_zero_null <- vapply(seq_len(B), function(j) wp3_fcg_exact_zero_pairing(indices$CF459_null_rna[j, ], indices$CF459_null_lipid[j, ]), logical(1L))
  degeneracy <- list()
  perms <- wp2_all_permutations(4L)
  newbank <- function(n) matrix(NA_real_, n, length(specs), dimnames = list(NULL, labels))
  nulls <- list(permutation_CF459 = newbank(24L), bootstrap_null = newbank(B))
  historical_null_corrected <- newbank(B)
  scores <- nulls; paired <- list(); edge_tables <- list(); observed <- list(); prepared_FCG <- list()
  replay <- list()
  for (name in expected) {
    x <- contexts[[name]]; e <- x$edges; ref <- reference[[name]]
    if (!identical(ref$metadata$mouse_id, first_meta$mouse_id) || !identical(ref$metadata$group, first_meta$group)) stop("Original male zones must share ordered mice and genotype strata.")
    if (!identical(colnames(ref$rna), first_meta$mouse_id) || !identical(colnames(ref$lipids), first_meta$mouse_id)) stop("Original FCG RNA/lipid/metadata order disagrees.")
    if (length(intersect(first_meta$mouse_id, c(x$CF453$mouse_ids, x$CF459$mouse_ids)))) stop("The three cohorts must be independent mice.")
    if (!"reference_lipid" %in% names(e)) stop("Original measured lipid identity is required for the FCG join.")
    key <- paste(e$gene, e$reference_lipid, sep = "\r")
    refkey <- paste(ref$edges$gene, ref$edges$reference_lipid, sep = "\r")
    if (anyDuplicated(refkey)) stop("Ambiguous original reference edges.")
    ix <- match(key, refkey)
    if (anyNA(ix) || !identical(e$primary_edge, ref$edges$primary_edge[ix])) stop("Frozen full/displayed FCG membership changed.")
    fe <- data.frame(gene = e$gene, feature_id = e$reference_lipid)
    f <- wp2_prepare_edges(ref$rna, ref$lipids, fe); fr <- wp2_prepared_cor(f)
    if (any(!is.finite(fr)) || max(abs(fr - ref$edges$r_reference[ix])) > 1e-10 || max(abs(fr - e$r_reference)) > 1e-10) stop("Original FCG matrices fail exact frozen-edge replay.")
    prior <- analysis$edges[analysis$edges$context == name, , drop = FALSE]
    pi <- match(e$edge_id, prior$edge_id)
    if (anyNA(pi) || max(abs(e$r_CF453_WD - prior$r_CF453_WD[pi])) > 1e-10 || max(abs(e$r_CF459_WD - prior$r_CF459_WD[pi])) > 1e-10) stop("The two WD cohorts no longer replay the original WP3 inputs.")
    replay[[name]] <- data.frame(context = name, n_edges = nrow(e),
      max_abs_FCG_archive_difference = max(abs(fr - ref$edges$r_reference[ix])),
      max_abs_CF453_difference = max(abs(e$r_CF453_WD - prior$r_CF453_WD[pi])),
      max_abs_CF459_difference = max(abs(e$r_CF459_WD - prior$r_CF459_WD[pi])))
    e$r_FCG_SD <- fr; e$n_FCG <- 22L; edge_tables[[name]] <- e; prepared_FCG[[name]] <- f
    hits <- which(vapply(specs, function(s) s$context == name, logical(1L)))
    if (verbose) message("WP3 original-cohort comparison ", name, ": ", nrow(e), " unchanged edges; FCG n22 / CF453 n7 / CF459 n4.")
    for (k in hits) {
      s <- specs[[k]]; observed[[labels[k]]] <- wp3_fcg_triplet_stats(fr[s$index], e$r_CF453_WD[s$index], e$r_CF459_WD[s$index])
      old <- analysis$inference[match(labels[k], analysis$inference$specification), ]
      if (abs(observed[[labels[k]]]["r_CF453_CF459"] - old$coefficient_vector_r) > 1e-10) stop("Original CF453-CF459 observed summary changed.")
    }
    for (method in names(nulls)) {
      nr <- nrow(nulls[[method]])
      for (j in seq_len(nr)) {
        cr <- if (method == "permutation_CF459") wp2_prepared_cor(x$CF459, seq_len(4L), perms[j, ], fixed_margins = TRUE) else wp2_prepared_cor(x$CF459, indices$CF459_null_rna[j, ], indices$CF459_null_lipid[j, ])
        for (k in hits) {
          s <- specs[[k]]; raw <- cr[s$index]
          raw_new <- wp2_vector_stats(fr[s$index], raw)
          # First replay the historical calculation before applying the explicit
          # extension-only correction. Original result bytes stay preserved.
          check <- wp2_vector_stats(e$r_CF453_WD[s$index], raw)["coefficient_vector_r"]
          old <- analysis[[paste0(method, "_draws")]][[labels[k]]][j]
          if (!identical(unname(is.finite(check)), unname(is.finite(old))) || (is.finite(check) && abs(check - old) > 1e-10)) stop("Original conditional null-bank replay failed.")
          corrected <- raw
          if (method == "bootstrap_null" && certified_zero_null[j]) {
            corrected[is.finite(corrected)] <- 0
            degeneracy[[paste(name, j, k)]] <- data.frame(specification = labels[k], context = name,
              universe = s$universe, scope = s$scope, gene = s$gene, draw = j,
              n_edges = length(s$index), n_undefined_edges_before = sum(!is.finite(raw)),
              max_abs_finite_raw_coefficient = if (any(is.finite(raw))) max(abs(raw[is.finite(raw)])) else NA_real_,
              original_CF453_anchor_r = unname(old), raw_FCG_anchor_r = unname(raw_new["coefficient_vector_r"]),
              original_CF453_anchor_defined = is.finite(old), raw_FCG_anchor_defined = is.finite(raw_new["coefficient_vector_r"]),
              corrected_FCG_anchor_defined = FALSE,
              correction_reason = "exact_integer_source_contingency_independence_certifies_zero_covariance", stringsAsFactors = FALSE)
          }
          if (method == "bootstrap_null") historical_null_corrected[j, k] <- wp2_vector_stats(e$r_CF453_WD[s$index], corrected)["coefficient_vector_r"]
          z <- wp2_vector_stats(fr[s$index], corrected)
          nulls[[method]][j, k] <- z["coefficient_vector_r"]; scores[[method]][j, k] <- z["ULM_score"]
        }
      }
    }
    for (j in seq_len(paired_B)) {
      a <- wp2_prepared_cor(f, FCG_paired[j, ])
      b <- wp2_prepared_cor(x$CF453, indices$CF453_paired[j, ])
      c <- wp2_prepared_cor(x$CF459, indices$CF459_paired[j, ])
      for (k in hits) {
        s <- specs[[k]]; z <- wp3_fcg_triplet_stats(a[s$index], b[s$index], c[s$index])
        old <- analysis$paired_stability_draws[[labels[k]]][j]
        if (!identical(unname(is.finite(z["r_CF453_CF459"])), unname(is.finite(old))) || (is.finite(old) && abs(z["r_CF453_CF459"] - old) > 1e-10)) stop("Original paired-bank replay failed.")
        paired[[paste(name, j, k)]] <- data.frame(specification = labels[k], context = name,
          universe = s$universe, scope = s$scope, gene = s$gene, draw = j,
          n_edges = length(s$index), as.list(z), stringsAsFactors = FALSE)
      }
    }
  }
  paired <- do.call(rbind, paired); rownames(paired) <- NULL
  inference <- paired_summary <- historical_correction <- list()
  for (k in seq_along(specs)) {
    s <- specs[[k]]; z <- observed[[labels[k]]]
    exact <- wp3_null_tail(z["r_FCG_CF459"], nulls$permutation_CF459[, k], exact = TRUE)
    nb <- wp3_null_tail(z["r_FCG_CF459"], nulls$bootstrap_null[, k])
    old <- analysis$inference[match(labels[k], analysis$inference$specification), ]
    corrected_old <- wp3_null_tail(old$coefficient_vector_r, historical_null_corrected[, k])
    historical_correction[[k]] <- data.frame(specification = labels[k], context = s$context,
      universe = s$universe, scope = s$scope, gene = s$gene,
      observed_CF453_CF459_r = old$coefficient_vector_r,
      original_null_bootstrap_p = old$null_bootstrap_p,
      corrected_null_bootstrap_p = unname(corrected_old["p_value"]),
      original_valid_draws = old$bootstrap_valid_draws,
      original_invalid_draws = old$bootstrap_invalid_draws,
      corrected_valid_draws = unname(corrected_old["valid_draws"]),
      corrected_invalid_draws = unname(corrected_old["invalid_draws"]),
      newly_invalid_draws = old$bootstrap_valid_draws - unname(corrected_old["valid_draws"]),
      corrected_completion_lower = unname(corrected_old["completion_lower"]),
      corrected_completion_upper = unname(corrected_old["completion_upper"]),
      correction = "exact source-contingency zero correction; original WP3 file preserved", stringsAsFactors = FALSE)
    pv <- paired[paired$specification == labels[k], ]; vd <- is.finite(pv$delta_r)
    inference[[k]] <- data.frame(specification = labels[k], context = s$context, universe = s$universe,
      analysis_role = if (s$universe == "full") "primary" else "secondary_displayed",
      scope = s$scope, gene = s$gene, n_edges = length(s$index), n_FCG = 22L, n_CF453 = 7L, n_CF459 = 4L,
      as.list(z), pairing_CF459_given_FCG_p = unname(exact["p_value"]),
      permutation_CF459_draws = 24L, permutation_CF459_valid_draws = unname(exact["valid_draws"]),
      null_bootstrap_p = unname(nb["p_value"]), null_bootstrap_draws = B,
      bootstrap_valid_draws = unname(nb["valid_draws"]), bootstrap_invalid_draws = unname(nb["invalid_draws"]),
      bootstrap_completion_lower = unname(nb["completion_lower"]), bootstrap_completion_upper = unname(nb["completion_upper"]),
      paired_draws = paired_B, paired_delta_valid_draws = sum(vd), paired_delta_invalid_draws = sum(!vd),
      paired_delta_positive_fraction = if (any(vd)) mean(pv$delta_r[vd] > 0) else NA_real_, stringsAsFactors = FALSE)
    for (field in c("r_FCG_CF459", "r_CF453_CF459", "r_FCG_CF453", "delta_r")) {
      v <- pv[[field]]; valid <- is.finite(v)
      paired_summary[[paste(k, field)]] <- data.frame(specification = labels[k], context = s$context, universe = s$universe,
        scope = s$scope, gene = s$gene, statistic = field, draws = paired_B,
        valid_draws = sum(valid), invalid_draws = sum(!valid), positive_fraction = if (any(valid)) mean(v[valid] > 0) else NA_real_,
        minimum = if (any(valid)) min(v[valid]) else NA_real_, maximum = if (any(valid)) max(v[valid]) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  # Original 28 joint + 7 CF453-only + 4 CF459-only omissions, followed by
  # 22 FCG-only omissions. No edge is removed when an omitted vector is constant.
  omission <- analysis$loco_omission_ledger
  omission$omit_FCG <- NA_integer_; omission$omitted_FCG_mouse <- NA_character_
  extra <- omission[rep(NA_integer_, 22L), , drop = FALSE]
  extra$omit_CF453 <- extra$omit_CF459 <- NA_integer_
  extra$omit_FCG <- seq_len(22L); extra$omitted_FCG_mouse <- first_meta$mouse_id
  extra$method <- "FCG_only"; extra$draw <- max(omission$draw) + seq_len(22L)
  omission <- rbind(omission, extra); rownames(omission) <- NULL
  omdraw <- list(); omindex <- list()
  for (j in seq_len(nrow(omission))) {
    keep <- list(FCG = if (is.na(omission$omit_FCG[j])) seq_len(22L) else setdiff(seq_len(22L), omission$omit_FCG[j]),
                 CF453 = if (is.na(omission$omit_CF453[j])) seq_len(7L) else setdiff(seq_len(7L), omission$omit_CF453[j]),
                 CF459 = if (is.na(omission$omit_CF459[j])) seq_len(4L) else setdiff(seq_len(4L), omission$omit_CF459[j]))
    for (cohort in names(keep)) {
      ix <- keep[[cohort]]; ids <- if (cohort == "FCG") first_meta$mouse_id else contexts[[1L]][[cohort]]$mouse_ids
      omindex[[paste(j, cohort)]] <- data.frame(draw = omission$draw[j], method = omission$method[j], cohort = cohort, position = seq_along(ix), source_position = ix, source_mouse = ids[ix])
    }
    for (name in expected) {
      x <- contexts[[name]]; a <- wp2_prepared_cor(prepared_FCG[[name]], keep$FCG)
      b <- wp2_prepared_cor(x$CF453, keep$CF453); c <- wp2_prepared_cor(x$CF459, keep$CF459)
      hits <- which(vapply(specs, function(s) s$context == name, logical(1L)))
      for (k in hits) {
        s <- specs[[k]]; z <- wp3_fcg_triplet_stats(a[s$index], b[s$index], c[s$index])
        omdraw[[paste(j, k)]] <- data.frame(specification = labels[k], context = name, universe = s$universe,
          scope = s$scope, gene = s$gene, draw = omission$draw[j], method = omission$method[j],
          n_edges = length(s$index), n_FCG = length(keep$FCG), n_CF453 = length(keep$CF453), n_CF459 = length(keep$CF459),
          as.list(z), stringsAsFactors = FALSE)
      }
    }
  }
  omdraw <- do.call(rbind, omdraw); rownames(omdraw) <- NULL
  omsummary <- list()
  for (key in labels) for (method in unique(omission$method)) {
    rows <- omdraw[omdraw$specification == key & omdraw$method == method, ]
    v <- rows$delta_r; valid <- is.finite(v)
    omsummary[[paste(key, method)]] <- data.frame(specification = key, context = rows$context[1L], universe = rows$universe[1L],
      scope = rows$scope[1L], gene = rows$gene[1L], method = method, draws = nrow(rows), valid_draws = sum(valid), invalid_draws = sum(!valid),
      minimum_delta_r = if (any(valid)) min(v[valid]) else NA_real_, maximum_delta_r = if (any(valid)) max(v[valid]) else NA_real_,
      delta_positive_fraction = if (any(valid)) mean(v[valid] > 0) else NA_real_, stringsAsFactors = FALSE)
  }
  asbank <- function(x) data.frame(draw = seq_len(nrow(x)), x, check.names = FALSE)
  result <- list(edges = do.call(rbind, edge_tables), inference = do.call(rbind, inference),
    skipped_specifications = specset$skipped, replay = do.call(rbind, replay),
    historical_null_correction_summary = do.call(rbind, historical_correction),
    historical_null_corrected_draws = asbank(historical_null_corrected),
    null_degeneracy_ledger = if (length(degeneracy)) do.call(rbind, degeneracy) else data.frame(),
    exact_zero_null_draws = data.frame(draw = seq_len(B), exact_zero_source_contingency = certified_zero_null),
    paired_draws = paired, paired_summary = do.call(rbind, paired_summary),
    FCG_paired_index_matrix = FCG_paired, FCG_paired_indices = index_ledger, paired_draw_summary = unique_counts,
    omission_ledger = omission, omission_indices = do.call(rbind, omindex), omission_draws = omdraw,
    omission_summary = do.call(rbind, omsummary),
    original_WP3_index_matrices = indices,
    exact_CF459_indices = wp3_index_ledger(perms, "CF459", "exact_pairing", contexts[[1L]]$CF459$mouse_ids))
  for (method in names(nulls)) {
    result[[paste0(method, "_draws")]] <- asbank(nulls[[method]])
    result[[paste0(method, "_ULM_draws")]] <- asbank(scores[[method]])
  }
  result$metadata <- list(seed = seed, FCG_paired_seed = as.integer((as.double(seed) + 104729) %% .Machine$integer.max), B = B, paired_B = paired_B,
    contrast = "delta_r = cor(CF459 WD, original FCG SD) - cor(CF459 WD, CF453 WD); positive means closer to original FCG SD on the same fixed edges",
    cohort_age_labels = "User-corrected 14 October 2026: CF459 3-month WD, original FCG 3-month SD, CF453 10-month WD. Age at sampling is user-supplied cohort provenance; no diet-duration or strain inference.",
    age_report_note = "Ages clarified by the user on 14 October 2026: CF459 male WD and original male FCG SD mice were 3 months old; CF453 male WD mice were 10 months old. These are ages at sampling, not established diet durations.",
    age_figure_note = "User-confirmed ages: CF459 WD 3 months; original FCG SD 3 months; CF453 WD 10 months.",
    fixed_universe = "Identical original WP3 strict shared 610-lipid male full/displayed edge universe; no new feature selection",
    FCG_replay = "Fresh original male RNA/log2-lipid matrices replay archived Figure 5 coefficients before matching by gene and exact reference_lipid",
    null = "24 CF459 pairing permutations and existing 4999 empirical-marginal CF459 null draws conditional on fixed original FCG coefficients; two-sided absolute ULM equivalent to absolute r at fixed edge count",
    paired = "Original independent CF453/CF459 paired mouse banks reused; one additional independent FCG bank paired across RNA/lipids within 11 XXT and 11 XYT; the SAME CF459 draw enters both anchor comparisons and both zones",
    direct_comparison = "Observed delta and paired/omission delta stability are descriptive. Delta-positive bootstrap fraction is conditional on defined draws, not a null P or posterior probability.",
    undefined = "Any undefined frozen edge invalidates its entire coefficient vector; delta requires both anchor comparisons defined; no dynamic edge deletion",
    exact_zero_null_correction = "Integer source-identity contingency independence certifies all feature covariances mathematically zero; finite roundoff coefficients are set exactly zero, making coefficient-vector r undefined. No absolute coefficient threshold is used. Original WP3 null outputs are preserved and historical discrepancies recorded in null_degeneracy_ledger.",
    omissions = "Original 39 WD omission configurations plus all 22 original-FCG single omissions; both zones share indices",
    nominal_ULM = "Native intercept slope t; nominal edge-level P descriptive because edges share genes/lipids/mice",
    limitations = "FC459 n4, original selected FCG network, genotype/sample-size and cohort processing differences; age and diet are confounded across cohorts. No causal separation of age, diet, strain, diet duration or technical batch.")
  result
}
