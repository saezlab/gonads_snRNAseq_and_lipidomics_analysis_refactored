source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/cyp2c38_models.R")
source("scripts/core/cyp2c38_reference.R")
set.seed(20261007)
rng_before <- .Random.seed
source("scripts/core/cyp2c38_ulm_bootstrap.R")
stopifnot(identical(rng_before, .Random.seed))

# Intercept ULM equivalence and original KO direction on nontrivial synthetic data.
prediction <- stats::rnorm(120L)
effect <- -0.2 * prediction + stats::rnorm(120L)
u <- wp1_ulm_stats(prediction, effect)
fit_ols <- summary(stats::lm(effect ~ prediction))$coefficients[2L, ]
network <- data.frame(source = "Cyp2c38", target = paste0("lipid_", seq_along(effect)), mor = prediction)
mat <- matrix(effect, ncol = 1L, dimnames = list(network$target, "KO_minus_WT"))
d <- decoupleR::run_ulm(mat, network, minsize = 10L)
stopifnot(abs(u["score"] - fit_ols["t value"]) < 1e-10,
          abs(u["slope"] - fit_ols["Estimate"]) < 1e-10,
          abs(u["p_value"] - fit_ols["Pr(>|t|)"]) < 1e-10,
          abs(u["score"] - d$score) < 1e-10,
          abs(u["p_value"] - d$p_value) < 1e-10,
          abs(wp1_ulm_stats(prediction, -effect)["score"] + u["score"]) < 1e-10,
          abs(wp1_ulm_stats(prediction * 1e-20, effect)["score"] - u["score"]) < 1e-10,
          is.na(wp1_ulm_stats(rep(0, 120L), effect)["score"]),
          is.na(wp1_ulm_stats(rep(1, 120L), effect)["score"]),
          is.na(wp1_ulm_stats(1 + 1e-12 * prediction, effect)["score"]),
          is.na(wp1_ulm_stats(prediction, rep(0, 120L))["score"]),
          is.na(wp1_ulm_stats(prediction, replace(effect, 1L, NA_real_))["score"]))

# Replay all historical nominal ULM scores/p-values directly from persisted input.
out <- "results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo"
old <- utils::read.csv(file.path(out, "signed_ULM_descriptive_scores.csv"))
ledger <- utils::read.csv(file.path(out, "fixed_reference_and_matching_ledger.csv"))
ledger <- ledger[ledger$match_status == "unique_exact", ]
tm <- utils::read.csv(file.path(out, "t_statistic_matrix.csv"), check.names = FALSE)
for (j in seq_len(nrow(old))) {
  values <- tm[[old$condition[j]]][match(ledger$feature_id, tm$feature_id)]
  replay <- wp1_ulm_stats(ledger[[paste0("r_", old$reference_set[j])]], values)
  stopifnot(abs(replay["score"] - old$score[j]) < 1e-10,
            abs(replay["p_value"] - old$p_value[j]) < 1e-10)
}

meta <- data.frame(sample_id = paste0("mouse_", seq_len(32L)),
                   group = rep(c("F_WT", "F_KO", "M_WT", "M_KO"), c(8L, 8L, 7L, 9L)))
meta$sex <- substr(meta$group, 1L, 1L)
y <- matrix(stats::rnorm(120L * 32L), 120L, 32L)
y[, meta$sex == "M"] <- y[, meta$sex == "M"] * 1.8 + 3
y[, meta$group == "F_KO"] <- y[, meta$group == "F_KO"] + prediction
y[, meta$group == "M_KO"] <- y[, meta$group == "M_KO"] - 0.6 * prediction
dimnames(y) <- list(paste0("lipid_", seq_len(120L)), meta$sample_id)
primary <- wp1_fit(y, meta)
X <- primary$design
w <- primary$weights
for (contrast in colnames(wp1_contrasts())) {
  a <- wp1_null_components(y, meta, contrast = contrast, fit = primary)
  b <- wp1_null_components(y, meta, restricted = TRUE, contrast = contrast, fit = primary)
  cvec <- wp1_contrasts()[, contrast]
  stopifnot(max(abs(a$null_group_coefficients %*% cvec)) < 1e-10,
            max(abs(b$null_group_coefficients %*% cvec)) < 1e-10,
            identical(a$null_mean, b$null_mean),
            identical(dimnames(y), dimnames(a$residuals_hc2)),
            all(a$leverage > b$leverage - 1e-12),
            all(a$leverage >= 0 & a$leverage < 1),
            all(b$leverage >= 0 & b$leverage < 1))
  # Independently fit a three-dimensional kernel parameterization of the
  # constrained weighted least-squares model and reconstruct four group means.
  Z <- qr.Q(qr(matrix(cvec, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  X0 <- X %*% Z
  constrained <- limma::lmFit(y, X0, weights = w)$coefficients %*% t(Z)
  H0 <- X0 %*% solve(crossprod(X0, X0 * w)) %*% t(X0 * w)
  stopifnot(max(abs(constrained - a$null_group_coefficients)) < 1e-10,
            max(abs(diag(H0) - b$leverage)) < 1e-10)
}
female <- wp1_null_components(y, meta, contrast = "female_KO_vs_WT", fit = primary)
male <- wp1_null_components(y, meta, contrast = "male_KO_vs_WT", fit = primary)
pooled <- wp1_null_components(y, meta, contrast = "pooled_KO_vs_WT", fit = primary)
stopifnot(identical(female$group_coefficients[, c("M_WT", "M_KO")],
                    female$null_group_coefficients[, c("M_WT", "M_KO")]),
          identical(male$group_coefficients[, c("F_WT", "F_KO")],
                    male$null_group_coefficients[, c("F_WT", "F_KO")]),
          max(abs(pooled$null_group_coefficients[, "F_KO"] -
                    pooled$null_group_coefficients[, "F_WT"])) > 0.1)
global <- wp1_null_components(y, meta)
stopifnot(max(abs(global$null_group_coefficients[, "F_WT"] -
                    global$null_group_coefficients[, "F_KO"])) < 1e-12,
          max(abs(global$null_group_coefficients[, "M_WT"] -
                    global$null_group_coefficients[, "M_KO"])) < 1e-12,
          identical(as.integer(global$group_counts), c(8L, 8L, 7L, 9L)))
signs <- rep(c(-1, 1), 16L)
draw <- wp1_wild_null(female, signs)
stopifnot(identical(dimnames(draw), dimnames(y)),
          max(abs((draw - female$null_mean) -
                    sweep(female$residuals_hc2, 2L, signs, "*"))) < 1e-12)
set.seed(418)
draw1 <- wp1_wild_null(female)
set.seed(418)
draw2 <- wp1_wild_null(female)
stopifnot(identical(draw1, draw2),
          inherits(try(wp1_wild_null(female, rep(0, 32L)), silent = TRUE), "try-error"),
          inherits(try(wp1_null_components(y, meta[c(2L, 1L, 3:32), ]), silent = TRUE), "try-error"))

# Reference bootstrap maintains paired RNA/lipid mouse indexing and feature order.
reference <- wp1_load_reference()
replayed <- wp1_reference_recompute(reference, seq_along(reference$gene))
stopifnot(identical(replayed$reference_lipid, reference$coefficients$reference_lipid),
          max(abs(replayed$r_pooled - reference$coefficients$r_pooled)) < 1e-12,
          max(abs(replayed$r_female - reference$coefficients$r_female)) < 1e-12,
          max(abs(replayed$r_male - reference$coefficients$r_male)) < 1e-12)
indices <- unlist(lapply(split(seq_along(reference$gene), reference$metadata$group),
                        function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
rebuilt <- wp1_reference_recompute(reference, indices)
for (sex in c("female", "male")) {
  gonad <- if (sex == "female") "Ovaries" else "Testes"
  ii <- indices[reference$metadata$gonadal_type[indices] == gonad]
  direct <- vapply(seq_len(nrow(reference$lipids)), function(j) {
    stats::cor(reference$gene[ii], reference$lipids[j, ii])
  }, numeric(1L))
  stopifnot(max(abs(rebuilt[[paste0("r_", sex)]] - direct)) < 1e-12)
}
cat("ULM equivalence, historical replay, contrast-null projection, nuisance means, HC2 leverage, shared whole-mouse signs, RNG and paired-reference tests passed.\n")
