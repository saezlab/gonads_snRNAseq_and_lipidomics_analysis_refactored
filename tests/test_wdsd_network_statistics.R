source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
set.seed(20261007L)
rng_before <- .Random.seed
source("scripts/core/wdsd_network_statistics.R")
stopifnot(identical(rng_before, .Random.seed))

# The vectorized row-maximum optimization is exactly equivalent to the
# pre-optimization normalization, including extremes and undefined rows.
legacy_unit_rows <- function(x) {
  scale <- apply(abs(x), 1L, max)
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
special <- rbind(
  ordinary = stats::rnorm(7L), huge = seq(-3, 3) * 1e307,
  tiny = seq(-3, 3) * 1e-307, zero = rep(0, 7L),
  constant = rep(10, 7L), near_constant = 1 + seq_len(7L) * 1e-12,
  missing = c(NA_real_, stats::rnorm(6L)), nan = c(NaN, stats::rnorm(6L)),
  infinite = c(Inf, stats::rnorm(6L)), all_infinite = rep(-Inf, 7L),
  mixed_undefined = c(NaN, NA_real_, Inf, -Inf, 1, 2, 3))
stopifnot(identical(wp2_unit_rows(special), legacy_unit_rows(special)),
          identical(wp2_unit_rows(special["ordinary", , drop = FALSE]),
                    legacy_unit_rows(special["ordinary", , drop = FALSE])))

# Native intercept ULM agrees with independently fitted OLS and decoupleR.
x <- stats::rnorm(120L); y <- 0.6 * x + stats::rnorm(120L)
native <- wp2_vector_stats(x, y)
ols <- summary(stats::lm(y ~ x))$coefficients[2L, ]
network <- data.frame(source = "example", target = paste0("edge_", seq_along(x)), mor = x)
observations <- matrix(y, ncol = 1L, dimnames = list(network$target, "observed"))
ulm <- decoupleR::run_ulm(observations, network, minsize = 3L)
stopifnot(abs(native["ULM_score"] - ols["t value"]) < 1e-10,
          abs(native["ULM_p_value"] - ols["Pr(>|t|)"]) < 1e-10,
          abs(native["slope"] - ols["Estimate"]) < 1e-10,
          abs(native["coefficient_vector_r"] - stats::cor(x, y)) < 1e-12,
          abs(native["ULM_score"] - ulm$score) < 1e-10,
          abs(native["ULM_p_value"] - ulm$p_value) < 1e-10,
          abs(wp2_vector_stats(x, -y)["ULM_score"] + native["ULM_score"]) < 1e-10,
          is.na(wp2_vector_stats(rep(1, 120L), y)["ULM_score"]),
          is.na(wp2_vector_stats(x, replace(y, 1L, NA_real_))["ULM_score"]))

# Matrix-based edge correlations agree with independent Pearson calculations,
# including row permutations, resampled repeated mice and explicit constants.
g <- matrix(stats::rnorm(4L * 7L), 4L, 7L,
            dimnames = list(c("Mtmr11", "GeneB", "GeneC", "Rare"), paste0("SD_", seq_len(7L))))
g["Rare", ] <- c(1, rep(0, 6L))
l <- matrix(stats::rnorm(4L * 7L), 4L, 7L,
            dimnames = list(paste0("lipid_", seq_len(4L)), colnames(g)))
e <- expand.grid(gene = rownames(g), feature_id = rownames(l), stringsAsFactors = FALSE)
e$edge_id <- paste(e$gene, e$feature_id, sep = "::")
e$primary_edge <- e$gene != "Rare"
answer <- wp2_cor_edges(g, l, e)
manual <- vapply(seq_len(nrow(e)), function(j) stats::cor(g[e$gene[j], ], l[e$feature_id[j], ]), numeric(1L))
stopifnot(max(abs(answer - manual)) < 1e-12,
          identical(answer, wp2_cor_edges(g[c(3, 1, 4, 2), ], l[c(4, 2, 1, 3), ], e)))
ix <- c(1L, 2L, 2L, 4L, 5L, 6L, 7L)
lx <- c(7L, 7L, 4L, 3L, 2L, 2L, 1L)
repeated <- wp2_edge_correlations(g, l, e, ix, lx)
manual <- vapply(seq_len(nrow(e)), function(j) stats::cor(g[e$gene[j], ix], l[e$feature_id[j], lx]), numeric(1L))
stopifnot(max(abs(repeated - manual)) < 1e-12)
constant <- wp2_edge_correlations(g, l, e, rep(2L, 7L), rep(3L, 7L))
stopifnot(all(is.na(constant)),
          inherits(try(wp2_cor_edges(g, l[, rev(seq_len(7L))], e), silent = TRUE), "try-error"))
missing_edge <- e[1L, ]; missing_edge$gene <- "unmeasured"
stopifnot(is.na(wp2_cor_edges(g, l, missing_edge)))

# Exact 7! enumeration includes identity, has no duplicate pairing and each row
# uses every animal once. Enumeration creates no pseudo-replicated lipid draw.
permutations <- wp2_all_permutations(7L)
stopifnot(nrow(permutations) == 5040L, ncol(permutations) == 7L,
          identical(as.integer(permutations[1L, ]), seq_len(7L)),
          !anyDuplicated(apply(permutations, 1L, paste, collapse = ",")),
          all(apply(permutations, 1L, function(z) identical(sort(as.integer(z)), seq_len(7L)))))

wd_g <- matrix(stats::rnorm(4L * 7L), 4L, 7L,
               dimnames = list(rownames(g), paste0("WD_", seq_len(7L))))
wd_l <- matrix(stats::rnorm(4L * 7L), 4L, 7L,
               dimnames = list(rownames(l), colnames(wd_g)))
ref_g <- matrix(stats::rnorm(4L * 11L), 4L, 11L,
                dimnames = list(rownames(g), paste0("FCG_", seq_len(11L))))
ref_l <- matrix(stats::rnorm(4L * 11L), 4L, 11L,
                dimnames = list(rownames(l), colnames(ref_g)))
ref_meta <- data.frame(mouse_id = colnames(ref_g), group = rep(c("XXT", "XYT"), c(6L, 5L)))
e$r_reference <- wp2_cor_edges(ref_g, ref_l, e)
context <- list(rna_SD = g, lipid_SD = l, rna_WD = wd_g, lipid_WD = wd_l,
                reference_rna = ref_g, reference_lipid = ref_l,
                reference_metadata = ref_meta, edges = e)
central <- context
central$rna_SD[c("Mtmr11", "GeneB"), ] <- central$rna_SD[c("GeneB", "Mtmr11"), ]
central$reference_rna[c("Mtmr11", "GeneB"), ] <- central$reference_rna[c("GeneB", "Mtmr11"), ]
central$edges$r_reference <- wp2_cor_edges(central$reference_rna, central$reference_lipid, e)
contexts <- list(male_periportal = context, male_pericentral = central)
fit <- wp2_network_analysis(contexts, B = 49L, paired_B = 29L, focus_genes = c("Mtmr11", "Rare"), verbose = FALSE)
fit2 <- wp2_network_analysis(contexts, B = 49L, paired_B = 29L, focus_genes = c("Mtmr11", "Rare"), workers = 3L, verbose = FALSE)
stopifnot(identical(fit$inference, fit2$inference),
          identical(fit$bootstrap_null_draws, fit2$bootstrap_null_draws),
          identical(fit$resampling_indices, fit2$resampling_indices),
          all(fit$inference$permutation_valid_draws == 5040L),
          any(fit$inference$bootstrap_invalid_draws > 0L),
          all(fit$inference$bootstrap_valid_draws + fit$inference$bootstrap_invalid_draws == 49L))

# Identity must reproduce the observed coefficient-vector statistic exactly;
# exact p-value is the full-enumeration tail, not a Monte Carlo +1 estimator.
for (j in seq_len(nrow(fit$inference))) {
  summary <- fit$inference[j, ]
  column <- summary$specification
  vals <- fit$permutation_draws[[column]]
  stopifnot(abs(vals[1L] - summary$coefficient_vector_r) < 1e-12,
            abs(summary$pairing_permutation_p - mean(abs(vals) >= abs(summary$coefficient_vector_r) - 1e-12)) < 1e-12)
  null <- fit$bootstrap_null_draws[[column]]
  valid <- is.finite(null)
  expected <- (1 + sum(abs(null[valid]) >= abs(summary$coefficient_vector_r) - 1e-12)) / (1 + sum(valid))
  stopifnot(abs(expected - summary$null_bootstrap_p) < 1e-12)
}

# Each modality's entire mouse vector is resampled once, not per edge. The same
# sex-level index bank is used for both zones; reference group strata are kept.
ledger <- fit$resampling_indices
nr <- ledger[ledger$method == "null_rna", ]
nl <- ledger[ledger$method == "null_lipid", ]
stopifnot(!identical(nr$source_position, nl$source_position),
          nrow(nr) == 49L * 7L, nrow(nl) == 49L * 7L)
rr <- ledger[ledger$method == "paired_reference", ]
stopifnot(all(ref_meta$group[rr$position] == ref_meta$group[rr$source_position]))
# Reconstruct a null draw from its saved bank, then verify the network statistic.
draw <- 7L
ri <- nr$source_position[nr$draw == draw]; li <- nl$source_position[nl$draw == draw]
for (name in names(contexts)) {
  z <- contexts[[name]]
  val <- wp2_edge_correlations(z$rna_SD, z$lipid_SD, z$edges, ri, li)
  column <- paste(name, "primary", "network", sep = "__")
  expected <- wp2_vector_stats(z$edges$r_reference[z$edges$primary_edge], val[z$edges$primary_edge])["coefficient_vector_r"]
  stopifnot(isTRUE(all.equal(unname(expected), fit$bootstrap_null_draws[[column]][draw], tolerance = 1e-12)))
}

# Joint LOCO computes all 49 independent group omissions, keeps shared zones
# synchronized, and flags undefined rare-feature comparisons instead of hiding
# them. Every defined result is independently checked against direct Pearson r.
rare <- fit$edges$gene == "Rare"
stopifnot(all(fit$edges$loco_valid_combinations[rare] == 42L),
          all(!fit$edges$loco_all_49_same_delta_direction[rare]),
          all(!fit$edges$loco_all_49_sign_reversal[rare]))
for (name in names(contexts)) {
  z <- contexts[[name]]
  omissions <- fit$loco_omission_ledger[fit$loco_omission_ledger$context == name, ]
  bank <- fit$loco_matrices[[name]]$delta_r
  stopifnot(nrow(omissions) == 49L, ncol(bank) == 49L,
            identical(as.character(omissions$omitted_SD_mouse), colnames(g)[omissions$omit_SD]),
            identical(as.character(omissions$omitted_WD_mouse), colnames(wd_g)[omissions$omit_WD]))
  edge_index <- which(z$edges$gene == "Mtmr11")[1L]
  for (j in seq_len(49L)) {
    keep_SD <- setdiff(seq_len(7L), omissions$omit_SD[j])
    keep_WD <- setdiff(seq_len(7L), omissions$omit_WD[j])
    edge <- z$edges[edge_index, ]
    expected <- stats::cor(z$rna_WD[edge$gene, keep_WD], z$lipid_WD[edge$feature_id, keep_WD]) -
      stats::cor(z$rna_SD[edge$gene, keep_SD], z$lipid_SD[edge$feature_id, keep_SD])
    stopifnot(abs(bank[edge_index, j] - expected) < 1e-12)
  }
}

# Invalid cross-zone mouse order and a mismatched historical replay are errors.
bad <- contexts; bad$male_pericentral$rna_SD <- bad$male_pericentral$rna_SD[, rev(seq_len(7L))]
bad$male_pericentral$lipid_SD <- bad$male_pericentral$lipid_SD[, rev(seq_len(7L))]
stopifnot(inherits(try(wp2_network_analysis(bad, B = 1L, paired_B = 1L, verbose = FALSE), silent = TRUE), "try-error"))
bad <- contexts; bad$male_pericentral$edges$r_reference[1L] <- 0.999
stopifnot(inherits(try(wp2_network_analysis(bad, B = 1L, paired_B = 1L, verbose = FALSE), silent = TRUE), "try-error"))

cat("WDSD network-statistics tests passed: exact Pearson/ULM equivalence, whole-mouse nulls, frozen edges, shared zones, reference strata, and joint LOCO.\n")
