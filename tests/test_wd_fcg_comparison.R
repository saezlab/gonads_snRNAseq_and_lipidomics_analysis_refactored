# Rscript tests/test_wd_fcg_comparison.R (repository root).
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/wdsd_network_statistics.R")
source("scripts/core/wd_crosscohort_statistics.R")
set.seed(941L); state <- .Random.seed
source("scripts/core/wd_fcg_comparison.R")
stopifnot(identical(state, .Random.seed))
same <- function(a, b, tol = 1e-10) stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol)))
fails <- function(x) stopifnot(inherits(try(force(x), silent = TRUE), "try-error"))
mat <- function(genes, n, prefix) matrix(stats::rnorm(length(genes) * n), length(genes), n, dimnames = list(genes, paste0(prefix, seq_len(n))))
genes <- c("Mtmr11", "GeneB", "Rare"); lipids <- paste0("L", 1:4)
base <- list(rna_CF453 = mat(genes, 7L, "CF453_"), rna_CF459 = mat(genes, 4L, "CF459_"),
             lipid_CF453 = mat(lipids, 7L, "CF453_"), lipid_CF459 = mat(lipids, 4L, "CF459_"))
base$rna_CF459["Rare", ] <- c(1, 0, 0, 0)
base$rna_CF453["Rare", ] <- c(1, rep(0, 6L))
ref <- list(rna = mat(genes, 22L, "FCG_"), lipids = mat(lipids, 22L, "FCG_"))
ref$metadata <- data.frame(mouse_id = colnames(ref$rna), group = rep(c("XXT", "XYT"), 11L))
e <- expand.grid(gene = genes, reference_lipid = lipids, stringsAsFactors = FALSE)
e$feature_id <- e$reference_lipid; e$primary_edge <- e$gene != "Rare"
e$r_reference <- wp2_cor_edges(ref$rna, ref$lipids, e)
base$edges <- e; ref$edges <- e[, c("gene", "reference_lipid", "r_reference", "primary_edge")]
central <- base; central$rna_CF459["GeneB", ] <- -central$rna_CF459["GeneB", ]
cref <- ref; cref$rna["GeneB", ] <- -cref$rna["GeneB", ]; cref$edges$r_reference <- wp2_cor_edges(cref$rna, cref$lipids, e)
central$edges$r_reference <- cref$edges$r_reference
contexts <- list(male_periportal = base, male_pericentral = central)
reference <- list(male_periportal = ref, male_pericentral = cref)
legacy <- wp3_crosscohort_analysis(contexts, B = 19L, paired_B = 31L, verbose = FALSE)
legacy_copy <- legacy; context_copy <- contexts; reference_copy <- reference
state <- .Random.seed
fit <- wp3_fcg_comparison(contexts, reference, legacy, verbose = FALSE)
stopifnot(identical(state, .Random.seed), identical(legacy, legacy_copy), identical(contexts, context_copy), identical(reference, reference_copy),
          nrow(fit$inference) == 8L, nrow(fit$edges) == 24L, nrow(fit$omission_ledger) == 61L,
          sum(fit$omission_ledger$method == "FCG_only") == 22L,
          identical(fit$original_WP3_index_matrices, legacy$resampling_index_matrices))
# FCG sampling is paired across modalities, stratified 11+11 and uses a distinct
# seed stream. Even reduced original null banks cannot overlap its seed stream.
ix <- fit$FCG_paired_index_matrix; seed <- fit$metadata$FCG_paired_seed
stopifnot(seed != legacy$metadata$seed, identical(dim(ix), c(31L, 22L)),
          all(fit$FCG_paired_indices$source_group == fit$FCG_paired_indices$position_group))
set.seed(seed); independent <- matrix(NA_integer_, 31L, 22L)
for (g in c("XXT", "XYT")) {
  positions <- which(ref$metadata$group == g)
  independent[, positions] <- matrix(positions[sample.int(11L, 11L * 31L, replace = TRUE)], 31L, 11L)
}
stopifnot(identical(ix, independent))
# Independent per-specification base Pearson, OLS and decoupleR agreement.
for (name in names(contexts)) {
  cx <- contexts[[name]]; rf <- reference[[name]]
  manual <- function(g, l, ri = seq_len(ncol(g)), li = ri) {
    suppressWarnings(vapply(seq_len(nrow(e)), function(j) stats::cor(g[e$gene[j], ri], l[e$reference_lipid[j], li]), numeric(1L)))
  }
  a <- manual(rf$rna, rf$lipids); b <- manual(cx$rna_CF453, cx$lipid_CF453); c <- manual(cx$rna_CF459, cx$lipid_CF459)
  edge <- fit$edges[fit$edges$context == name, ]; same(edge$r_FCG_SD, a); same(edge$r_CF453_WD, b); same(edge$r_CF459_WD, c)
  for (k in which(fit$inference$context == name)) {
    rec <- fit$inference[k, ]; use <- if (rec$universe == "full") seq_len(nrow(e)) else which(e$primary_edge)
    if (rec$scope == "gene") use <- use[e$gene[use] == rec$gene]
    same(rec$r_FCG_CF459, stats::cor(a[use], c[use])); same(rec$r_CF453_CF459, stats::cor(b[use], c[use])); same(rec$r_FCG_CF453, stats::cor(a[use], b[use]))
    same(rec$delta_r, stats::cor(a[use], c[use]) - stats::cor(b[use], c[use]))
    ols <- summary(stats::lm(c[use] ~ a[use]))$coefficients[2L, ]
    same(rec$ULM_FCG_CF459, ols["t value"]); same(rec$ULM_p_FCG_CF459, ols["Pr(>|t|)"])
    network <- data.frame(source = "FCG", target = paste0("edge", use), mor = a[use])
    dec <- decoupleR::run_ulm(matrix(c[use], ncol = 1L, dimnames = list(network$target, "CF459")), network, minsize = 3L)
    same(rec$ULM_FCG_CF459, dec$score); same(rec$ULM_p_FCG_CF459, dec$p_value)
    # Every exact pairing is independently reconstructed; identity is first.
    pp <- wp2_all_permutations(4L)
    exact <- vapply(seq_len(24L), function(j) stats::cor(a[use], manual(cx$rna_CF459, cx$lipid_CF459, 1:4, pp[j, ])[use]), numeric(1L))
    same(fit$permutation_CF459_draws[[rec$specification]], exact)
    same(rec$pairing_CF459_given_FCG_p, mean(abs(exact) >= abs(rec$r_FCG_CF459) - 1e-12))
    # Shared CF459 draw simultaneously enters BOTH comparisons. Undefined
    # vectors retain all edges and invalidate delta; no dynamic edge deletion.
    for (j in c(1L, 5L, 19L, 31L)) {
      aa <- manual(rf$rna, rf$lipids, ix[j, ])
      bb <- manual(cx$rna_CF453, cx$lipid_CF453, legacy$resampling_index_matrices$CF453_paired[j, ])
      cc <- manual(cx$rna_CF459, cx$lipid_CF459, legacy$resampling_index_matrices$CF459_paired[j, ])
      z <- wp3_fcg_triplet_stats(aa[use], bb[use], cc[use])
      row <- fit$paired_draws[fit$paired_draws$specification == rec$specification & fit$paired_draws$draw == j, ]
      for (field in c("r_FCG_CF459", "r_CF453_CF459", "r_FCG_CF453", "delta_r")) same(row[[field]], z[field])
    }
    null <- fit$bootstrap_null_draws[[rec$specification]]; valid <- is.finite(null)
    same(rec$null_bootstrap_p, (1 + sum(abs(null[valid]) >= abs(rec$r_FCG_CF459) - 1e-12)) / (1 + sum(valid)))
    stopifnot(rec$bootstrap_valid_draws == sum(valid), rec$bootstrap_invalid_draws == sum(!valid))
    rows <- fit$paired_draws[fit$paired_draws$specification == rec$specification, ]; vd <- is.finite(rows$delta_r)
    stopifnot(rec$paired_delta_valid_draws == sum(vd), rec$paired_delta_invalid_draws == sum(!vd),
              all(is.na(rows$delta_r) == (!is.finite(rows$r_FCG_CF459) | !is.finite(rows$r_CF453_CF459))))
    same(rec$paired_delta_positive_fraction, mean(rows$delta_r[vd] > 0))
    # Explicit one-omission replay for all three mouse cohorts.
    for (cohort in c("CF459", "CF453", "FCG")) {
      ledger <- fit$omission_ledger[fit$omission_ledger$method == paste0(cohort, "_only"), ][1L, ]
      aa <- if (cohort == "FCG") manual(rf$rna, rf$lipids, 2:22) else a
      bb <- if (cohort == "CF453") manual(cx$rna_CF453, cx$lipid_CF453, 2:7) else b
      cc <- if (cohort == "CF459") manual(cx$rna_CF459, cx$lipid_CF459, 2:4) else c
      z <- wp3_fcg_triplet_stats(aa[use], bb[use], cc[use])
      row <- fit$omission_draws[fit$omission_draws$specification == rec$specification & fit$omission_draws$draw == ledger$draw, ]
      same(row$delta_r, z["delta_r"])
    }
  }
}
for (j in seq_len(nrow(fit$historical_null_correction_summary))) {
  row <- fit$historical_null_correction_summary[j, ]; old <- legacy$inference[match(row$specification, legacy$inference$specification), ]
  bank <- fit$historical_null_corrected_draws[[row$specification]]; valid <- is.finite(bank)
  same(row$original_null_bootstrap_p, old$null_bootstrap_p)
  same(row$corrected_null_bootstrap_p, (1 + sum(abs(bank[valid]) >= abs(old$coefficient_vector_r) - 1e-12)) / (1 + sum(valid)))
  stopifnot(row$corrected_valid_draws == sum(valid), row$corrected_invalid_draws == sum(!valid),
            row$newly_invalid_draws == old$bootstrap_valid_draws - sum(valid))
}
stopifnot(any(fit$inference$paired_delta_invalid_draws > 0L))
again <- wp3_fcg_comparison(contexts, reference, legacy, verbose = FALSE)
stopifnot(identical(fit, again))
# Align by names, not row order. Mismatched mice/genotypes/archive are rejected.
swapped <- lapply(reference, function(x) { x$rna <- x$rna[3:1, , drop = FALSE]; x$lipids <- x$lipids[4:1, , drop = FALSE]; x$edges <- x$edges[rev(seq_len(nrow(x$edges))), ]; x })
row_fit <- wp3_fcg_comparison(contexts, swapped, legacy, verbose = FALSE)
stopifnot(identical(fit, row_fit))
bad <- reference; bad$male_periportal$metadata$group[1L] <- "XYT"; fails(wp3_fcg_comparison(contexts, bad, legacy, verbose = FALSE))
bad <- reference; bad$male_periportal$lipids <- bad$male_periportal$lipids[, 22:1]; fails(wp3_fcg_comparison(contexts, bad, legacy, verbose = FALSE))
bad <- reference; bad$male_periportal$edges$r_reference[1L] <- 0.2; fails(wp3_fcg_comparison(contexts, bad, legacy, verbose = FALSE))
bad <- legacy; bad$resampling_index_matrices$CF459_paired[1L, 1L] <- 5L; fails(wp3_fcg_comparison(contexts, reference, bad, verbose = FALSE))
bad <- legacy; bad$bootstrap_null_draws[1L, 2L] <- 0.999; fails(wp3_fcg_comparison(contexts, reference, bad, verbose = FALSE))
# Certify mathematical zero from the integer source-identity contingency;
# retain weak varying coefficient vectors rather than imposing an amplitude cut.
stopifnot(wp3_fcg_exact_zero_pairing(c(2, 2, 4, 4), c(4, 2, 2, 4)),
          !wp3_fcg_exact_zero_pairing(c(2, 2, 4, 4), c(2, 2, 4, 4)),
          is.finite(wp2_vector_stats(c(1, 2, 4, 8), c(-2, -1, 1, 2) * 1e-14)["coefficient_vector_r"]))
forced <- legacy
forced$resampling_index_matrices$CF459_null_rna[1L, ] <- c(2, 2, 4, 4)
forced$resampling_index_matrices$CF459_null_lipid[1L, ] <- c(4, 2, 2, 4)
for (name in names(contexts)) {
  cx <- wp3_normalize_context(contexts[[name]], name)
  raw <- wp2_prepared_cor(cx$CF459, c(2, 2, 4, 4), c(4, 2, 2, 4))
  specs <- wp3_specifications(stats::setNames(list(cx), name))$specs
  for (key in names(specs)) {
    use <- specs[[key]]$index; z <- wp2_vector_stats(cx$edges$r_CF453_WD[use], raw[use])
    forced$bootstrap_null_draws[[key]][1L] <- z["coefficient_vector_r"]
    forced$bootstrap_null_ULM_draws[[key]][1L] <- z["ULM_score"]
  }
}
fixed <- wp3_fcg_comparison(contexts, reference, forced, verbose = FALSE)
stopifnot(all(is.na(as.numeric(fixed$bootstrap_null_draws[1L, -1L]))),
          fixed$exact_zero_null_draws$exact_zero_source_contingency[1L],
          sum(fixed$null_degeneracy_ledger$draw == 1L) == 8L,
          all(!fixed$null_degeneracy_ledger$corrected_FCG_anchor_defined),
          identical(fixed$paired_draws, fit$paired_draws), identical(fixed$omission_draws, fit$omission_draws))
cat("WP3 original FCG comparison: native ULM, all 24 exact pairings, same-edge delta, independent stratified FCG stream, shared CF459 pairing, exact zero-contingency guard without amplitude threshold, undefined vectors and 61 omissions passed.\n")
