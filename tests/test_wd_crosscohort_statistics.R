# Run from the repository root: Rscript tests/test_wd_crosscohort_statistics.R
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/wdsd_network_statistics.R")
set.seed(117L); rng_before <- .Random.seed
source("scripts/core/wd_crosscohort_statistics.R")
stopifnot(identical(rng_before, .Random.seed))

fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
same <- function(a, b, tol = 1e-10) stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol)))

# Independently fitted OLS and the actual decoupleR ULM reproduce native scores.
x <- stats::rnorm(67L); y <- 0.45 * x + stats::rnorm(67L)
native <- wp2_vector_stats(x, y)
ols <- summary(stats::lm(y ~ x))$coefficients[2L, ]
net <- data.frame(source = "frozen", target = paste0("edge_", seq_along(x)), mor = x)
ulm <- decoupleR::run_ulm(matrix(y, ncol = 1L, dimnames = list(net$target, "CF459")), net, minsize = 3L)
same(native["ULM_score"], ols["t value"]); same(native["ULM_p_value"], ols["Pr(>|t|)"])
same(native["ULM_score"], ulm$score); same(native["ULM_p_value"], ulm$p_value)
same(native["coefficient_vector_r"], stats::cor(x, y)); same(native["slope"], ols["Estimate"])
stopifnot(is.na(wp2_vector_stats(x, replace(y, 1L, NA_real_))["ULM_score"]))

g453 <- matrix(stats::rnorm(4L * 7L), 4L, 7L,
               dimnames = list(c("Mtmr11", "GeneB", "GeneC", "Rare"), paste0("CF453_", seq_len(7L))))
g459 <- matrix(stats::rnorm(4L * 4L), 4L, 4L,
               dimnames = list(rownames(g453), paste0("CF459_", seq_len(4L))))
g453["Rare", ] <- c(1, rep(0, 6L)); g459["Rare", ] <- c(1, 0, 0, 0)
l453 <- matrix(stats::rnorm(4L * 7L), 4L, 7L, dimnames = list(paste0("L", 1:4), colnames(g453)))
l459 <- matrix(stats::rnorm(4L * 4L), 4L, 4L, dimnames = list(rownames(l453), colnames(g459)))
edge <- expand.grid(gene = rownames(g453), feature_id = rownames(l453), stringsAsFactors = FALSE)
edge$primary_edge <- edge$gene != "Rare" & edge$feature_id != "L4"
edge$r_reference <- seq(-0.8, 0.8, length.out = nrow(edge))
edge$r_reference[1L] <- NA_real_ # FCG coefficient is provenance, not an outcome-selection gate.
context <- list(rna_CF453 = g453, lipid_CF453 = l453, rna_CF459 = g459,
                lipid_CF459 = l459, edges = edge)
central <- context
central$rna_CF453[c("Mtmr11", "GeneB"), ] <- central$rna_CF453[c("GeneB", "Mtmr11"), ]
central$rna_CF459["GeneC", ] <- -central$rna_CF459["GeneC", ]
contexts <- list(male_periportal = context, male_pericentral = central)
fit <- wp3_crosscohort_analysis(contexts, B = 53L, paired_B = 31L, focus_genes = c("Mtmr11", "absent"), verbose = FALSE)
stopifnot(nrow(fit$inference) == 8L, nrow(fit$skipped_specifications) == 4L,
          all(fit$inference$n_CF453 == 7L), all(fit$inference$n_CF459 == 4L),
          all(fit$inference$permutation_CF459_valid_draws == 24L),
          all(fit$inference$permutation_CF453_valid_draws == 5040L),
          nrow(fit$edges) == 32L, all(fit$edges$eligible),
          all(fit$inference$bootstrap_valid_draws + fit$inference$bootstrap_invalid_draws == 53L),
          any(fit$inference$bootstrap_invalid_draws > 0L))

# Independent Pearson edges, identity, exact tails, and reciprocal enumeration.
for (name in names(contexts)) {
  cx <- contexts[[name]]
  actual <- fit$edges[fit$edges$context == name, ]
  manual453 <- vapply(seq_len(nrow(edge)), function(i) stats::cor(cx$rna_CF453[edge$gene[i], ], cx$lipid_CF453[edge$feature_id[i], ]), numeric(1L))
  manual459 <- vapply(seq_len(nrow(edge)), function(i) stats::cor(cx$rna_CF459[edge$gene[i], ], cx$lipid_CF459[edge$feature_id[i], ]), numeric(1L))
  same(actual$r_CF453_WD, manual453); same(actual$r_CF459_WD, manual459)
  key <- paste(name, "full", "network", sep = "__")
  rec <- fit$inference[fit$inference$specification == key, ]
  expected <- summary(stats::lm(manual459 ~ manual453))$coefficients[2L, ]
  same(rec$ULM_score, expected["t value"]); same(rec$ULM_p_value, expected["Pr(>|t|)"])
  for (cohort in c("CF459", "CF453")) {
    n <- if (cohort == "CF459") 4L else 7L
    perm <- wp2_all_permutations(n)
    stopifnot(nrow(perm) == factorial(n), identical(as.integer(perm[1L, ]), seq_len(n)),
              !anyDuplicated(apply(perm, 1L, paste, collapse = ",")),
              all(apply(perm, 1L, function(ix) identical(sort(as.integer(ix)), seq_len(n)))))
    bank <- fit[[paste0("permutation_", cohort, "_draws")]][[key]]
    same(bank[1L], rec$coefficient_vector_r)
    pcol <- paste0("pairing_", cohort, "_given_", if (cohort == "CF459") "CF453" else "CF459", "_p")
    same(rec[[pcol]], mean(abs(bank) >= abs(rec$coefficient_vector_r) - 1e-12))
    stopifnot(rec[[pcol]] >= 1 / factorial(n))
    # Explicit Pearson reconstruction at several nonidentity pairings detects
    # routing errors between the unequal-sized cohort banks.
    for (j in unique(c(1L, 2L, floor(nrow(perm) / 2L), nrow(perm)))) {
      cxg <- cx[[paste0("rna_", cohort)]]; cxl <- cx[[paste0("lipid_", cohort)]]
      rr <- vapply(seq_len(nrow(edge)), function(i) stats::cor(cxg[edge$gene[i], ], cxl[edge$feature_id[i], perm[j, ]]), numeric(1L))
      expected_r <- if (cohort == "CF459") stats::cor(manual453, rr) else stats::cor(rr, manual459)
      same(bank[j], expected_r)
    }
  }
}

# Raw r and ULM draw statistics share the expected fixed-n monotone relation.
for (method in c("permutation_CF459", "permutation_CF453", "bootstrap_null", "paired_stability")) {
  for (k in seq_len(nrow(fit$inference))) {
    rec <- fit$inference[k, ]; key <- rec$specification
    r <- fit[[paste0(method, "_draws")]][[key]]
    score <- fit[[paste0(method, "_ULM_draws")]][[key]]
    valid <- is.finite(r)
    same(score[valid], r[valid] * sqrt((rec$n_edges - 2) / (1 - r[valid]^2)), tol = 1e-8)
    stopifnot(identical(is.na(r), is.na(score)))
  }
}

# Every feature and both zones use one common whole-mouse draw, not edge-level
# pseudo-replicates. Null modalities are independent; paired modalities share
# their cohort indices. Unique-mouse counts quantify effective n=4 degeneracy.
ix <- fit$resampling_index_matrices
stopifnot(identical(dim(ix$CF459_null_rna), c(53L, 4L)),
          identical(dim(ix$CF453_paired), c(31L, 7L)),
          !identical(ix$CF459_null_rna, ix$CF459_null_lipid),
          all(fit$resampling_draw_summary$n_unique_mice <= fit$resampling_draw_summary$n_positions),
          any(fit$resampling_draw_summary$n_unique_mice < 3L))
for (method in names(ix)) {
  ledger <- fit$resampling_indices[fit$resampling_indices$method == method, ]
  same(ledger$source_position, as.vector(t(ix[[method]])))
  cohort <- sub("_.*$", "", method)
  ids <- colnames(context[[paste0("rna_", cohort)]])
  stopifnot(identical(ledger$source_mouse, ids[ledger$source_position]))
}
for (name in names(contexts)) {
  cx <- contexts[[name]]; e <- fit$edges[fit$edges$context == name, ]
  for (j in c(1L, 5L, 19L, 31L)) {
    a <- wp2_edge_correlations(cx$rna_CF453, cx$lipid_CF453, edge, ix$CF453_paired[j, ])
    b <- wp2_edge_correlations(cx$rna_CF459, cx$lipid_CF459, edge, ix$CF459_paired[j, ])
    expected <- wp2_vector_stats(a, b)["coefficient_vector_r"]
    same(fit$paired_stability_draws[[paste(name, "full", "network", sep = "__")]][j], expected)
  }
  for (j in c(1L, 3L, 17L, 53L)) {
    b <- wp2_edge_correlations(cx$rna_CF459, cx$lipid_CF459, edge, ix$CF459_null_rna[j, ], ix$CF459_null_lipid[j, ])
    expected <- wp2_vector_stats(e$r_CF453_WD, b)["coefficient_vector_r"]
    same(fit$bootstrap_null_draws[[paste(name, "full", "network", sep = "__")]][j], expected)
  }
}
key <- "male_periportal__full__network"
drop_rare <- which(apply(ix$CF459_null_rna, 1L, function(z) !1L %in% z))
stopifnot(length(drop_rare) > 0L, all(is.na(fit$bootstrap_null_draws[[key]][drop_rare])))
stopifnot(all(is.na(wp2_edge_correlations(g459, l459, edge, rep(1L, 4L), rep(2L, 4L)))))

# Bootstrap P is conditional on the defined frozen-vector draws, with a +1
# estimate. Completion bounds expose undefined null mass without calling them
# confidence intervals or claiming calibration under this composite setting.
for (k in seq_len(nrow(fit$inference))) {
  rec <- fit$inference[k, ]; values <- fit$bootstrap_null_draws[[rec$specification]]
  valid <- is.finite(values); exceed <- sum(abs(values[valid]) >= abs(rec$coefficient_vector_r) - 1e-12)
  same(rec$null_bootstrap_p, (1 + exceed) / (1 + sum(valid)))
  same(rec$bootstrap_completion_lower, (1 + exceed) / (1 + length(values)))
  same(rec$bootstrap_completion_upper, (1 + exceed + sum(!valid)) / (1 + length(values)))
  stopifnot(rec$bootstrap_completion_lower <= rec$null_bootstrap_p,
            rec$null_bootstrap_p <= rec$bootstrap_completion_upper)
}
undefined_tail <- wp3_null_tail(0.5, rep(NA_real_, 5L))
stopifnot(is.na(undefined_tail["p_value"]), undefined_tail["completion_lower"] == 1 / 6,
          undefined_tail["completion_upper"] == 1)

# Exactly 28 joint and 7/4 single-cohort omissions retain both zones' mice and
# never alter the observed eligible edge universe. One explicit replay checks
# the n=6/n=3 correlation calculation and mouse identity ledger.
lc <- fit$loco_omission_ledger
stopifnot(nrow(lc) == 39L, sum(lc$method == "joint_omission") == 28L,
          sum(lc$method == "CF453_only") == 7L, sum(lc$method == "CF459_only") == 4L,
          !anyDuplicated(paste(lc$omit_CF453[1:28], lc$omit_CF459[1:28])),
          identical(lc$omitted_CF453_mouse[1:7], colnames(g453)))
for (name in names(contexts)) {
  cx <- contexts[[name]]
  a <- wp2_edge_correlations(cx$rna_CF453, cx$lipid_CF453, edge, 2:7)
  b <- wp2_edge_correlations(cx$rna_CF459, cx$lipid_CF459, edge, 2:4)
  same(fit$loco_matrices[[name]]$r_CF453_WD[, 1L], a)
  same(fit$loco_matrices[[name]]$r_CF459_WD[, 1L], b)
  stopifnot(all(fit$loco_draws$n_edges[fit$loco_draws$context == name & fit$loco_draws$universe == "full" & fit$loco_draws$scope == "network"] == 16L))
}

# Reproducibility and feature-row invariance preserve all banks exactly.
again <- wp3_crosscohort_analysis(contexts, B = 53L, paired_B = 31L, focus_genes = c("Mtmr11", "absent"), verbose = FALSE)
stopifnot(identical(fit, again))
swapped <- lapply(contexts, function(cx) {
  for (field in c("rna_CF453", "lipid_CF453", "rna_CF459", "lipid_CF459")) cx[[field]] <- cx[[field]][rev(seq_len(nrow(cx[[field]]))), , drop = FALSE]
  cx
})
row_replay <- wp3_crosscohort_analysis(swapped, B = 53L, paired_B = 31L, focus_genes = c("Mtmr11", "absent"), verbose = FALSE)
stopifnot(identical(fit$inference, row_replay$inference), identical(fit$bootstrap_null_draws, row_replay$bootstrap_null_draws))

# Simultaneous RNA/lipid order changes preserve the observed networks and exact
# conditional tails, although finite Monte Carlo draws route to different IDs.
column_order <- lapply(contexts, function(cx) {
  for (cohort in c("CF453", "CF459")) for (prefix in c("rna_", "lipid_")) {
    field <- paste0(prefix, cohort); cx[[field]] <- cx[[field]][, rev(seq_len(ncol(cx[[field]]))), drop = FALSE]
  }
  cx
})
ordered <- wp3_crosscohort_analysis(column_order, B = 2L, paired_B = 2L, focus_genes = c("Mtmr11", "absent"), verbose = FALSE)
for (field in c("coefficient_vector_r", "ULM_score", "ULM_p_value", "pairing_CF459_given_CF453_p", "pairing_CF453_given_CF459_p")) same(fit$inference[[field]], ordered$inference[[field]])

# Input contracts reject ambiguous IDs, orientations, contexts and sample size.
fails(wp3_crosscohort_analysis(list(unknown = context), B = 1L, paired_B = 1L, verbose = FALSE))
fails(wp3_crosscohort_analysis(contexts, B = 0L)); fails(wp3_crosscohort_analysis(contexts, B = 1.5))
fails(wp3_crosscohort_analysis(contexts, B = 1e20))
fails(wp3_crosscohort_analysis(contexts, seed = NA_real_))
bad <- context; bad$rna_CF459 <- t(bad$rna_CF459)
fails(wp3_normalize_context(bad, "male_periportal"))
bad <- context; bad$lipid_CF459 <- bad$lipid_CF459[, 4:1]
fails(wp3_normalize_context(bad, "male_periportal"))
bad <- context; colnames(bad$rna_CF453)[2L] <- colnames(bad$rna_CF453)[1L]
fails(wp3_normalize_context(bad, "male_periportal"))
bad <- context; bad$rna_CF459 <- bad$rna_CF459[, 1:3]; bad$lipid_CF459 <- bad$lipid_CF459[, 1:3]
fails(wp3_normalize_context(bad, "male_periportal"))
bad <- context; bad$edges <- rbind(bad$edges, bad$edges[1L, ])
fails(wp3_normalize_context(bad, "male_periportal"))
bad <- central
bad$rna_CF459 <- bad$rna_CF459[, 4:1]; bad$lipid_CF459 <- bad$lipid_CF459[, 4:1]
fails(wp3_crosscohort_analysis(list(male_periportal = context, male_pericentral = bad), B = 1L, paired_B = 1L, verbose = FALSE))
missing <- context; missing$edges$gene[1L] <- "not_measured"
z <- wp3_normalize_context(missing, "male_periportal")
stopifnot(!z$edges$eligible[1L], z$edges$exclusion_reason[1L] == "undefined_CF453_WD")
cat("WP3 cross-cohort statistics: native ULM, exact 4!/7! nulls, fixed vectors, shared-zone indices, bootstrap and 28+7+4 omissions passed.\n")
