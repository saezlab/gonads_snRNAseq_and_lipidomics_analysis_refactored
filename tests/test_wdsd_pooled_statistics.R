source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/wdsd_network_statistics.R")
set.seed(20261007L)
rng_before <- .Random.seed
source("scripts/core/wdsd_pooled_statistics.R")
stopifnot(identical(rng_before, .Random.seed))

# Independent synthetic FCG/SD cohorts, with explicit between-sex mean effects
# and a rare observed gene that makes some fixed-universe bootstrap draws
# undefined. No resample may quietly discard that gene's edges.
genes <- c("Cyp2b9", "Cyp2b13", "Cyp2c38", "Cyp2d9", "Rare")
features <- paste0("lipid_", seq_len(5L))
sd_ids <- paste0("SD_", seq_len(14L)); ref_ids <- paste0("FCG_", seq_len(43L))
sd_sex <- rep(c("female", "male"), each = 7L)
ref_sex <- rep(c("female", "male"), c(21L, 22L))
gene_effect <- c(2, -3, 1, -1.5, 0)
lipid_effect <- c(1, -2, 3, -1, 0.5)
sd_gene <- matrix(stats::rnorm(5L * 14L), 5L, 14L, dimnames = list(genes, sd_ids)) +
  gene_effect %o% ifelse(sd_sex == "male", 1, -1)
sd_lipid <- matrix(stats::rnorm(5L * 14L), 5L, 14L, dimnames = list(features, sd_ids)) +
  lipid_effect %o% ifelse(sd_sex == "male", 1, -1)
sd_gene["Rare", ] <- c(1, rep(0, 13L))
ref_gene <- matrix(stats::rnorm(5L * 43L), 5L, 43L, dimnames = list(genes, ref_ids)) +
  gene_effect %o% ifelse(ref_sex == "male", 1, -1)
ref_lipid <- matrix(stats::rnorm(5L * 43L), 5L, 43L, dimnames = list(features, ref_ids)) +
  lipid_effect %o% ifelse(ref_sex == "male", 1, -1)
edges <- expand.grid(gene = genes, feature_id = features, stringsAsFactors = FALSE)
edges$primary_edge <- edges$gene != "Rare"
edges$edge_id <- paste(edges$gene, edges$feature_id, sep = "::")
edges$r_reference <- vapply(seq_len(nrow(edges)), function(i) {
  stats::cor(ref_gene[edges$gene[i], ], ref_lipid[edges$feature_id[i], ])
}, numeric(1L))
context <- list(rna_SD = sd_gene, lipid_SD = sd_lipid,
                metadata_SD = data.frame(mouse_id = sd_ids, sex = sd_sex),
                reference_rna = ref_gene, reference_lipid = ref_lipid,
                reference_metadata = data.frame(mouse_id = ref_ids, sex = ref_sex,
                  group = c(rep("XXO", 11L), rep("XYO", 10L), rep("XXT", 11L), rep("XYT", 11L))),
                edges = edges)
fit <- wp2_pooled_analysis(context, permutations = 129L, B = 79L, paired_B = 39L, verbose = FALSE)
repeat_fit <- wp2_pooled_analysis(context, permutations = 129L, B = 79L, paired_B = 39L, verbose = FALSE)
stopifnot(identical(fit$inference, repeat_fit$inference),
          identical(fit$resampling_indices, repeat_fit$resampling_indices),
          nrow(fit$edges) == 25L, all(fit$edges$eligible),
          fit$metadata$identity_audit_draw == 0L,
          !fit$metadata$identity_in_null_draw_count)

# Native ULM is independently checked against intercept OLS and decoupleR,
# while its p remains nominal because these synthetic edge residuals share
# the same mice/features just as the real network does.
primary <- fit$edges$primary_edge
ols <- summary(stats::lm(fit$edges$r_SD[primary] ~ fit$edges$r_reference[primary]))$coefficients[2L, ]
summary_primary <- fit$inference[fit$inference$analysis == "pooled" &
                                  fit$inference$scope == "network" & fit$inference$universe == "primary", ]
network <- data.frame(source = "Figure3", target = fit$edges$edge_id[primary], mor = fit$edges$r_reference[primary])
observations <- matrix(fit$edges$r_SD[primary], ncol = 1L, dimnames = list(network$target, "SD"))
ulm <- decoupleR::run_ulm(observations, network, minsize = 3L)
stopifnot(abs(summary_primary$ULM_score - ols["t value"]) < 1e-10,
          abs(summary_primary$ULM_p_value - ols["Pr(>|t|)"]) < 1e-10,
          abs(summary_primary$ULM_score - ulm$score) < 1e-10,
          abs(summary_primary$ULM_p_value - ulm$p_value) < 1e-10)

# Both identity audit rows are the observed network and do not enter the
# Monte Carlo random-null bank. Upper/absolute/equal tails are independently
# recalculated, preserving their different questions when the null has a
# positive alignment caused by the shared sex structure.
pooled <- fit$inference[fit$inference$analysis == "pooled", ]
bank_fields <- c(withinsex_permutation = "permutation_withinsex_draws",
                 unrestricted_permutation = "permutation_unrestricted_draws",
                 withinsex_bootstrap = "bootstrap_withinsex_draws",
                 unrestricted_bootstrap = "bootstrap_unrestricted_draws")
for (i in seq_len(nrow(pooled))) {
  row <- pooled[i, ]; label <- row$specification; observed <- row$coefficient_vector_r
  for (method in names(bank_fields)) {
    bank <- fit[[bank_fields[[method]]]]
    if (grepl("permutation", method)) {
      stopifnot(bank$draw[1L] == 0L, nrow(bank) == 130L,
                abs(bank[[label]][1L] - observed) < 1e-12)
      null <- bank[[label]][bank$draw > 0L]
    } else null <- bank[[label]]
    null <- null[is.finite(null)]
    expected_absolute <- (1 + sum(abs(null) >= abs(observed) - 1e-12)) / (1 + length(null))
    expected_upper <- (1 + sum(null >= observed - 1e-12)) / (1 + length(null))
    expected_lower <- (1 + sum(null <= observed + 1e-12)) / (1 + length(null))
    stopifnot(abs(row[[paste0(method, "_absolute_p")]] - expected_absolute) < 1e-12,
              abs(row[[paste0(method, "_upper_p")]] - expected_upper) < 1e-12,
              abs(row[[paste0(method, "_lower_p")]] - expected_lower) < 1e-12,
              abs(row[[paste0(method, "_equal_tail_p")]] - min(1, 2 * min(expected_upper, expected_lower))) < 1e-12)
  }
}
stopifnot(any(pooled$withinsex_bootstrap_invalid_draws > 0L),
          all(pooled$withinsex_bootstrap_valid_draws + pooled$withinsex_bootstrap_invalid_draws == 79L),
          all(pooled$unrestricted_bootstrap_valid_draws + pooled$unrestricted_bootstrap_invalid_draws == 79L))
ties <- wp2_pooled_null_tails(0.4, rep(0.4, 20L))
stopifnot(all(ties[c("absolute_p", "upper_p", "lower_p", "equal_tail_p")] == 1),
          is.na(wp2_pooled_null_tails(NA_real_, c(0.1, 0.2))["absolute_p"]),
          is.na(wp2_pooled_null_tails(0.2, c(NA_real_, NaN))["upper_p"]))
asymmetric <- wp2_pooled_null_tails(-0.8, c(0.2, 0.3, 0.4))
stopifnot(asymmetric["upper_p"] == 1, asymmetric["absolute_p"] == 0.25,
          asymmetric["equal_tail_p"] == 0.5)

# Every conditional permutation/sample stays within its recorded sex; FCG
# paired draws additionally preserve the four genotype-group allocations.
# Independence bootstrap banks sample the modalities separately, while each
# modality's entire feature vector shares one mouse draw.
ledger <- fit$resampling_indices
conditional <- ledger$method %in% c("permutation_withinsex", "bootstrap_withinsex_rna", "bootstrap_withinsex_lipid", "paired_SD", "paired_reference")
stopifnot(all(ledger$position_sex[conditional] == ledger$source_sex[conditional]),
          any(ledger$position_sex[ledger$method == "permutation_unrestricted"] !=
                ledger$source_sex[ledger$method == "permutation_unrestricted"]),
          !identical(fit$index_banks$bootstrap_withinsex_rna, fit$index_banks$bootstrap_withinsex_lipid))
reference_ledger <- ledger[ledger$method == "paired_reference", ]
stopifnot(all(context$reference_metadata$group[reference_ledger$position] ==
                context$reference_metadata$group[reference_ledger$source_position]))
for (i in seq_len(nrow(fit$index_banks$permutation_withinsex))) {
  stopifnot(identical(sort(as.integer(fit$index_banks$permutation_withinsex[i, ])), seq_len(14L)))
}
for (i in seq_len(nrow(fit$index_banks$permutation_unrestricted))) {
  stopifnot(identical(sort(as.integer(fit$index_banks$permutation_unrestricted[i, ])), seq_len(14L)))
}
draw <- 4L
for (kind in c("withinsex", "unrestricted")) {
  ri <- fit$index_banks[[paste0("bootstrap_", kind, "_rna")]][draw, ]
  li <- fit$index_banks[[paste0("bootstrap_", kind, "_lipid")]][draw, ]
  direct <- vapply(which(primary), function(i) {
    stats::cor(sd_gene[edges$gene[i], ri], sd_lipid[edges$feature_id[i], li])
  }, numeric(1L))
  expected <- stats::cor(edges$r_reference[primary], direct)
  stopifnot(abs(fit[[paste0("bootstrap_", kind, "_draws")]][["primary__network"]][draw] - expected) < 1e-12)
}

# Within-sex centered sensitivity independently removes each cohort's gene
# and lipid sex means. Affine offsets/feature scales do not affect the result.
independent_center <- function(x, sex) {
  t(apply(x, 1L, function(z) as.numeric(stats::residuals(stats::lm(z ~ factor(sex))))))
}
cg <- independent_center(sd_gene, sd_sex); cl <- independent_center(sd_lipid, sd_sex)
rg <- independent_center(ref_gene, ref_sex); rl <- independent_center(ref_lipid, ref_sex)
for (i in seq_len(nrow(edges))) {
  stopifnot(abs(fit$edges$r_SD_withinsex[i] - stats::cor(cg[edges$gene[i], ], cl[edges$feature_id[i], ])) < 1e-12,
            abs(fit$edges$r_reference_withinsex[i] - stats::cor(rg[edges$gene[i], ], rl[edges$feature_id[i], ])) < 1e-12)
}
centered <- fit$inference[fit$inference$analysis == "withinsex_centered", ]
stopifnot(all(is.na(centered$withinsex_permutation_upper_p)),
          all(is.na(centered$unrestricted_bootstrap_absolute_p)),
          all(centered$n_edges == pooled$n_edges))

# Pure sex-mediated associations yield apparently perfect pooled replication
# but conditional p=1, and undefined centered networks. This checks that the
# main conditional null does not call shared sex differences evidence of
# within-sex correspondence or invent residual correlations.
pure <- context
pure$rna_SD <- gene_effect %o% ifelse(sd_sex == "male", 1, -1)
pure$lipid_SD <- lipid_effect %o% ifelse(sd_sex == "male", 1, -1)
pure$reference_rna <- gene_effect %o% ifelse(ref_sex == "male", 1, -1)
pure$reference_lipid <- lipid_effect %o% ifelse(ref_sex == "male", 1, -1)
dimnames(pure$rna_SD) <- list(genes, sd_ids); dimnames(pure$lipid_SD) <- list(features, sd_ids)
dimnames(pure$reference_rna) <- list(genes, ref_ids); dimnames(pure$reference_lipid) <- list(features, ref_ids)
pure$edges <- edges[edges$primary_edge, ]
pure$edges$r_reference <- wp2_cor_edges(pure$reference_rna, pure$reference_lipid, pure$edges)
sex_only <- wp2_pooled_analysis(pure, permutations = 29L, B = 19L, paired_B = 9L, unrestricted = FALSE, verbose = FALSE)
sex_network <- sex_only$inference[sex_only$inference$analysis == "pooled" & sex_only$inference$scope == "network", ]
stopifnot(all(abs(sex_network$coefficient_vector_r - 1) < 1e-12),
          all(sex_network$withinsex_permutation_upper_p == 1),
          all(sex_network$withinsex_bootstrap_upper_p == 1),
          all(!sex_only$edges$withinsex_defined),
          all(is.na(sex_only$inference$ULM_score[sex_only$inference$analysis == "withinsex_centered"])))

# Explicit failures prevent accidental mouse reordering, incomplete sex
# groups, incompatible FCG strata, duplicate edges, or altered reference r.
bad <- context; bad$metadata_SD <- bad$metadata_SD[14:1, ]
stopifnot(inherits(try(wp2_pooled_context(bad), silent = TRUE), "try-error"))
bad <- context; bad$lipid_SD <- bad$lipid_SD[, 14:1]
stopifnot(inherits(try(wp2_pooled_context(bad), silent = TRUE), "try-error"))
bad <- context; bad$metadata_SD$sex[1L] <- "male"
stopifnot(inherits(try(wp2_pooled_context(bad), silent = TRUE), "try-error"))
bad <- context; bad$reference_metadata$group <- "merged_across_sex"
stopifnot(inherits(try(wp2_pooled_context(bad), silent = TRUE), "try-error"))
bad <- context; bad$edges$r_reference[1L] <- 0.99
stopifnot(inherits(try(wp2_pooled_context(bad), silent = TRUE), "try-error"))
bad <- context; bad$edges <- rbind(edges, edges[1L, ])
stopifnot(inherits(try(wp2_pooled_context(bad), silent = TRUE), "try-error"))

cat("Pooled Figure 3 statistics passed: ordered replay, native ULM, complete-mouse nulls, identity/+1 tails, sex/group strata, undefined frozen draws, and centered sensitivity.\n")
