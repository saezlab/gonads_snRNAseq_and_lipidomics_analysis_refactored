# Workflow 11 phase orchestration. Statistical implementations remain in their
# existing cyp2c38_* helpers. Source scripts/support_functions.R first.

wp1_run <- function(args = commandArgs(trailingOnly = TRUE)) {
  source("scripts/core/cyp2c38_reference.R")
  source("scripts/core/cyp2c38_models.R")
  option <- function(key, default) {
    a <- args[startsWith(args, paste0("--", key, "="))]
    if (length(a) > 1L) stop("Repeated option: ", key)
    if (length(a)) substring(a, nchar(key) + 4L) else default
  }
  allowed <- "^--(permutations|bootstraps|seed|workers|output)=|^--validate-only$"
  if (any(!grepl(allowed, args))) stop("Unknown CLI option")
  B <- as.integer(option("permutations", "1999"))
  K <- as.integer(option("bootstraps", "999"))
  seed <- as.integer(option("seed", "20261006"))
  workers <- as.integer(option("workers", "3"))
  stopifnot(B >= 19, K >= 19, workers >= 1, !is.na(seed))
  out <- option("output", "results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo")
  required <- c("limma", "statmod", "digest", "jsonlite", "decoupleR")
  for (p in required) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")

  input <- "data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv"
  raw <- utils::read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
  samples <- grep("^Cyp-[0-9]+_[FM]-(Ctrl|KO)$", names(raw), value = TRUE)
  qc <- grep("^QC[0-9]+$", names(raw), value = TRUE)
  meta <- data.frame(sample_id = samples,
                     mouse_id = sub("_.*$", "", samples),
                     sex = sub("^.*_([FM])-.*$", "\\1", samples),
                     genotype = ifelse(grepl("-Ctrl$", samples), "WT", "KO"),
                     stringsAsFactors = FALSE)
  meta$group <- paste(meta$sex, meta$genotype, sep = "_")
  stopifnot(!anyDuplicated(meta$sample_id), !anyDuplicated(meta$mouse_id),
            identical(as.integer(table(factor(meta$group, levels = c("F_WT", "F_KO", "M_WT", "M_KO")))), c(8L, 8L, 7L, 9L)),
            length(qc) == 9L)
  ledger <- raw[, setdiff(names(raw), c(samples, qc)), drop = FALSE]
  ledger$feature_id <- sprintf("KO_%04d", seq_len(nrow(raw)))
  ledger$canonical_lipid <- wp1_canonical_lipid(raw[["Metabolite name"]])
  ledger$lipid_class <- sub("\\(.*$", "", ledger$canonical_lipid)
  ledger$internal_standard <- grepl("^IS:", raw[["Metabolite name"]]) | raw$Ontology == "IS"
  ledger$isotope_labelled <- grepl("\\([dD][0-9]+\\)|\\b[0-9]+[Cc]\\b", raw[["Metabolite name"]])
  excluded_classes <- c("BA", "DGGA", "HBMP", "LPE-N", "NATau", "PEtOH", "PMeOH", "SL", "CAR", "SPB", "VAE")
  ledger$excluded_class <- ledger$lipid_class %in% excluded_classes
  qcmat <- as.matrix(raw[, qc, drop = FALSE])
  ledger$qc_cv_percent <- 100 * apply(qcmat, 1, stats::sd) / rowMeans(qcmat)
  ledger$exclusion_reason <- ifelse(ledger$internal_standard, "internal_standard",
                                  ifelse(ledger$isotope_labelled, "isotope_labelled_feature",
                                         ifelse(ledger$excluded_class, "repository_excluded_class", "")))
  ledger$keep_primary <- ledger$exclusion_reason == ""
  ledger$keep_qc30 <- ledger$keep_primary & ledger$qc_cv_percent <= 30
  x <- as.matrix(raw[ledger$keep_primary, samples, drop = FALSE])
  storage.mode(x) <- "double"
  rownames(x) <- ledger$feature_id[ledger$keep_primary]
  stopifnot(all(is.finite(x)), all(x > 0), !anyDuplicated(ledger$canonical_lipid[ledger$keep_primary]),
            sum(ledger$keep_primary) == 1237L, sum(ledger$keep_qc30) == 1215L)
  y <- log2(x)
  reference <- wp1_load_reference()
  refs <- reference$coefficients
  refs$feature_id <- ledger$feature_id[match(refs$reference_lipid, ledger$canonical_lipid)]
  refs$feature_id[!refs$feature_id %in% rownames(y)] <- NA_character_
  refs$match_status <- ifelse(is.na(refs$feature_id), "unmatched_exact_molecular_annotation", "unique_exact")
  # CL chain reordering is candidate annotation equivalence, never a primary join.
  cl_key <- function(v) vapply(v, function(s) {
    if (!grepl("^CL\\([^()]+\\)$", s)) return(s)
    paste0("CL(", paste(sort(strsplit(sub("^CL\\((.*)\\)$", "\\1", s), "_", fixed = TRUE)[[1]]), collapse = "_"), ")")
  }, character(1))
  refs$unordered_CL_candidate_annotations <- ""
  refs$unordered_CL_candidate_alignment_ids <- ""
  refs$unordered_CL_reference_multiplicity <- 0L
  rk <- cl_key(refs$reference_lipid)
  lk <- cl_key(ledger$canonical_lipid)
  for (j in which(is.na(refs$feature_id) & startsWith(refs$reference_lipid, "CL("))) {
    ii <- which(ledger$keep_primary & lk == rk[j])
    refs$unordered_CL_candidate_annotations[j] <- paste(ledger$canonical_lipid[ii], collapse = ";")
    refs$unordered_CL_candidate_alignment_ids[j] <- paste(ledger[["Alignment ID"]][ii], collapse = ";")
    refs$unordered_CL_reference_multiplicity[j] <- sum(rk == rk[j])
  }
  stopifnot(sum(!is.na(refs$feature_id)) == 552L,
            !anyDuplicated(stats::na.omit(refs$feature_id)), identical(colnames(y), meta$sample_id))
  map <- refs[!is.na(refs$feature_id), , drop = FALSE]
  cat("WP1 intake:", nrow(y), "endogenous features;", ncol(y), "mice;", nrow(map), "fixed reference matches\n")
  if ("--validate-only" %in% args) {
    cat("Input, sample counts, one-to-one mapping and archived reference replay passed.\n")
    return(invisible(list(metadata = meta, ledger = ledger, reference = reference, map = map, y = y, input = input)))
  }
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; choose a fresh --output destination.")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "plots"), showWarnings = FALSE)
  annotation <- ledger[, c("feature_id", "canonical_lipid", "lipid_class", "qc_cv_percent", "Fill %")]
  append_annotation <- function(d) {
    a <- annotation[match(d$feature_id, annotation$feature_id), , drop = FALSE]
    cbind(d, a[, setdiff(names(a), "feature_id"), drop = FALSE])
  }
  tabulate_fit <- function(result, model) {
    ans <- lapply(colnames(result$fit$coefficients), function(ct) {
      d <- wp1_table(result$fit, ct)
      d$contrast <- ct; d$model <- model
      d$effect_units <- if (model == "supplied_scale_sex_weighted") "supplied_intensity_mean_difference" else "log2_fold_change"
      append_annotation(d)
    })
    do.call(rbind, ans)
  }
  cat("Fitting fixed primary weighted model and prespecified sensitivities...\n")
  primary <- wp1_fit(y, meta)
  unweighted <- wp1_fit(y, meta, weighted = FALSE)
  qc30 <- wp1_fit(y[ledger$feature_id[ledger$keep_qc30], , drop = FALSE], meta)
  supplied_scale <- wp1_fit(x, meta)
  separate <- lapply(wp1_fit_separate(y, meta), function(f) list(fit = f))
  differential <- tabulate_fit(primary, "primary_log2_sex_weighted")
  sensitivity <- rbind(tabulate_fit(unweighted, "unweighted_log2"),
                       tabulate_fit(qc30, "qc_CV_le30_log2"),
                       tabulate_fit(supplied_scale, "supplied_scale_sex_weighted"),
                       do.call(rbind, lapply(separate, function(z) tabulate_fit(z, "separate_sex_log2"))))
  invariant <- c(reference_replay_max_difference = reference$max_archive_difference,
                 weighted_unweighted_logFC_max_difference = max(abs(primary$fit$coefficients - unweighted$fit$coefficients)),
                 pooled_identity_max_difference = max(abs(primary$fit$coefficients[, 1] - rowMeans(primary$fit$coefficients[, 2:3]))),
                 interaction_identity_max_difference = max(abs(primary$fit$coefficients[, 4] - (primary$fit$coefficients[, 3] - primary$fit$coefficients[, 2]))))
  stopifnot(all(invariant < 1e-10))
  meta$precision_weight <- primary$weights
  meta$total_supplied_intensity <- colSums(as.matrix(raw[, samples, drop = FALSE]))
  meta$total_primary_intensity <- colSums(x)
  write_csv(meta, "sample_manifest.csv")
  write_csv(ledger, "feature_annotation_and_QC_ledger.csv")
  write_csv(refs, "fixed_reference_and_matching_ledger.csv")
  write_csv(reference$metadata, "FCG_paired_reference_mouse_manifest.csv")
  write_csv(data.frame(feature_id = rownames(y), y, check.names = FALSE), "log2_analysis_matrix.csv")
  write_csv(data.frame(group = rownames(primary$contrasts), primary$contrasts, check.names = FALSE), "contrast_matrix.csv")
  write_csv(differential, "differential_primary_all.csv")
  write_csv(sensitivity, "differential_sensitivities_all.csv")
  write_csv(append_annotation(data.frame(feature_id = rownames(primary$fit$t), primary$fit$t, check.names = FALSE)), "t_statistic_matrix.csv")
  for (ct in unique(differential$contrast)) write_csv(differential[differential$contrast == ct, ], paste0("differential_", ct, ".csv"))

  specs <- data.frame(contrast = c("pooled_KO_vs_WT", "female_KO_vs_WT", "male_KO_vs_WT", "female_KO_vs_WT", "male_KO_vs_WT"),
                      reference_set = c("pooled", "female", "male", "pooled", "pooled"),
                      primary_test = c(TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
  score_vector <- function(tmat, coeff = refs) {
    vapply(seq_len(nrow(specs)), function(i) {
      ok <- !is.na(coeff$feature_id) & coeff$feature_id %in% rownames(tmat)
      r <- coeff[[paste0("r_", specs$reference_set[i])]][ok]
      tval <- tmat[coeff$feature_id[ok], specs$contrast[i]]
      if (any(!is.finite(r)) || any(!is.finite(tval)) || stats::sd(r) == 0 || stats::sd(tval) == 0) return(NA_real_)
      -stats::cor(r, tval)
    }, numeric(1))
  }
  consistency_for <- function(result, model) {
    rows <- lapply(seq_len(nrow(specs)), function(i) {
      m <- map[map$feature_id %in% rownames(result$fit$t), , drop = FALSE]
      if (!specs$contrast[i] %in% colnames(result$fit$t)) return(NULL)
      r <- m[[paste0("r_", specs$reference_set[i])]]
      d <- data.frame(feature_id = m$feature_id, contrast = specs$contrast[i], reference_set = specs$reference_set[i],
                      reference_r = r, t = result$fit$t[m$feature_id, specs$contrast[i]],
                      logFC = result$fit$coefficients[m$feature_id, specs$contrast[i]], model = model,
                      stringsAsFactors = FALSE)
      d$consistency <- -d$reference_r * d$t
      d$effect_units <- if (model == "supplied_scale_sex_weighted") "supplied_intensity_mean_difference" else "log2_fold_change"
      append_annotation(d)
    })
    do.call(rbind, rows)
  }
  consistency <- consistency_for(primary, "primary_log2_sex_weighted")
  sens_consistency <- rbind(consistency_for(unweighted, "unweighted_log2"),
                           consistency_for(qc30, "qc_CV_le30_log2"),
                           consistency_for(supplied_scale, "supplied_scale_sex_weighted"),
                           do.call(rbind, lapply(separate, function(z) consistency_for(z, "separate_sex_log2"))))
  summarize_consistency <- function(d) {
    keys <- unique(d[, c("contrast", "reference_set", "model")])
    do.call(rbind, lapply(seq_len(nrow(keys)), function(j) {
      a <- d[d$contrast == keys$contrast[j] & d$reference_set == keys$reference_set[j] & d$model == keys$model[j], ]
      cbind(keys[j, ], data.frame(n_lipids = nrow(a), pearson_r = stats::cor(a$reference_r, a$t),
                                       concordance_score = -stats::cor(a$reference_r, a$t),
                                       spearman_r = stats::cor(a$reference_r, a$t, method = "spearman"),
                                       slope_t_per_r = unname(stats::coef(stats::lm(t ~ reference_r, a))[2]),
                                       mean_consistency = mean(a$consistency), fraction_positive_consistency = mean(a$consistency > 0)))
    }))
  }
  concordance <- summarize_consistency(consistency)
  concordance$primary_test <- specs$primary_test
  write_csv(consistency, "per_lipid_predictive_consistency.csv")
  write_csv(sens_consistency, "per_lipid_consistency_sensitivities.csv")
  write_csv(summarize_consistency(sens_consistency), "concordance_sensitivities.csv")
  observed <- score_vector(primary$fit$t)

  # Seeds are generated in advance so results are identical for 1 or multiple workers.
  set.seed(seed)
  seeds_perm <- sample.int(.Machine$integer.max, B)
  seeds_boot <- sample.int(.Machine$integer.max, K)
  run_draws <- function(draw_seeds, fun, label) {
    ans <- vector("list", length(draw_seeds))
    for (start in seq.int(1L, length(draw_seeds), by = 100L)) {
      ii <- start:min(start + 99L, length(draw_seeds))
      ans[ii] <- parallel::mclapply(draw_seeds[ii], function(s) {
        set.seed(s)
        tryCatch(fun(), error = function(e) rep(NA_real_, nrow(specs)))
      }, mc.cores = workers, mc.set.seed = FALSE)
      cat(label, max(ii), "/", length(draw_seeds), "\n")
    }
    do.call(rbind, ans)
  }
  perm <- run_draws(seeds_perm, function() {
    mp <- meta
    for (sx in c("F", "M")) {
      ii <- which(mp$sex == sx)
      mp$genotype[ii] <- sample(mp$genotype[ii])
    }
    mp$group <- paste(mp$sex, mp$genotype, sep = "_")
    score_vector(wp1_fit(y, mp)$fit$t)
  }, "Within-sex KO permutations")
  # Preserve reference FCG chromosome/gonad groups when these metadata are supplied.
  strata_ref <- if ("group" %in% names(reference$metadata)) reference$metadata$group else reference$metadata$gonadal_type
  stopifnot(!anyNA(strata_ref))
  boot <- run_draws(seeds_boot, function() {
    ik <- unlist(lapply(split(seq_len(ncol(y)), meta$group), function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
    ib <- unlist(lapply(split(seq_len(ncol(reference$lipids)), strata_ref), function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
    yc <- y[, ik, drop = FALSE]
    mb <- meta[ik, , drop = FALSE]
    colnames(yc) <- mb$sample_id <- paste0("bootstrap_", seq_along(ik))
    rb <- refs
    gb <- reference$gene[ib]
    lb <- reference$lipids[, ib, drop = FALSE]
    for (sx in c("pooled", "female", "male")) {
      ii <- if (sx == "pooled") seq_along(ib) else which(reference$metadata$gonadal_type[ib] == if (sx == "female") "Ovaries" else "Testes")
      rb[[paste0("r_", sx)]] <- as.numeric(stats::cor(gb[ii], t(lb[, ii, drop = FALSE])))
    }
    score_vector(wp1_fit(yc, mb)$fit$t, rb)
  }, "Combined KO and paired FCG bootstrap")
  for (z in list(perm, boot)) stopifnot(all(colMeans(is.finite(z)) >= 0.95))
  labels <- paste(specs$contrast, specs$reference_set, sep = "__")
  colnames(perm) <- colnames(boot) <- labels
  write_csv(data.frame(draw = seq_len(B), seed = seeds_perm, perm, check.names = FALSE), "permutation_draws.csv")
  write_csv(data.frame(draw = seq_len(K), seed = seeds_boot, boot, check.names = FALSE), "combined_bootstrap_draws.csv")
  concordance$permutation_valid_draws <- colSums(is.finite(perm))
  concordance$bootstrap_valid_draws <- colSums(is.finite(boot))
  concordance$permutation_p_two_sided <- vapply(seq_len(nrow(specs)), function(i) (1 + sum(abs(perm[, i]) >= abs(observed[i]), na.rm = TRUE)) / (1 + sum(is.finite(perm[, i]))), numeric(1))
  concordance$permutation_p_directional <- vapply(seq_len(nrow(specs)), function(i) (1 + sum(perm[, i] >= observed[i], na.rm = TRUE)) / (1 + sum(is.finite(perm[, i]))), numeric(1))
  concordance$permutation_fdr <- NA_real_
  concordance$permutation_fdr[specs$primary_test] <- stats::p.adjust(concordance$permutation_p_two_sided[specs$primary_test], method = "BH")
  concordance$bootstrap_ci_low <- apply(boot, 2, stats::quantile, 0.025, na.rm = TRUE)
  concordance$bootstrap_ci_high <- apply(boot, 2, stats::quantile, 0.975, na.rm = TRUE)
  write_csv(concordance, "concordance_primary_and_secondary.csv")

  cat("Whole-animal influence checks and signed ULM...\n")
  loo_ko <- t(vapply(seq_len(nrow(meta)), function(j) score_vector(wp1_fit(y[, -j, drop = FALSE], meta[-j, ])$fit$t), numeric(nrow(specs))))
  loo_ref <- t(vapply(seq_len(nrow(reference$metadata)), function(j) {
    rb <- refs
    for (sx in c("pooled", "female", "male")) {
      ii <- setdiff(if (sx == "pooled") seq_along(reference$gene) else which(reference$metadata$gonadal_type == if (sx == "female") "Ovaries" else "Testes"), j)
      rb[[paste0("r_", sx)]] <- as.numeric(stats::cor(reference$gene[ii], t(reference$lipids[, ii, drop = FALSE])))
    }
    score_vector(primary$fit$t, rb)
  }, numeric(nrow(specs))))
  colnames(loo_ko) <- colnames(loo_ref) <- labels
  write_csv(data.frame(omitted_mouse = meta$sample_id, loo_ko, check.names = FALSE), "leave_one_KO_mouse_out.csv")
  write_csv(data.frame(omitted_mouse = reference$metadata$sample_id, loo_ref, check.names = FALSE), "leave_one_FCG_mouse_out.csv")
  ulm <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    mat <- matrix(primary$fit$t[map$feature_id, specs$contrast[i]], ncol = 1, dimnames = list(map$reference_lipid, specs$contrast[i]))
    network <- data.frame(source = "Cyp2c38", target = map$reference_lipid, mor = map[[paste0("r_", specs$reference_set[i])]])
    d <- as.data.frame(decoupleR::run_ulm(mat, network, minsize = 10L))
    d$reference_set <- specs$reference_set[i]
    d$inference_note <- "Edge-independence p-value is descriptive; negative signed score supports concordance"
    d
  }))
  write_csv(ulm, "signed_ULM_descriptive_scores.csv")
  discovery <- do.call(rbind, lapply(split(differential, differential$contrast), function(d) data.frame(contrast = d$contrast[1], tested = nrow(d), FDR_le_005 = sum(d$adj.P.Val <= 0.05), KO_higher = sum(d$adj.P.Val <= 0.05 & d$logFC > 0), KO_lower = sum(d$adj.P.Val <= 0.05 & d$logFC < 0))))
  write_csv(discovery, "differential_discovery_summary.csv")
  classes <- do.call(rbind, lapply(split(differential, interaction(differential$contrast, differential$lipid_class, drop = TRUE)), function(d) data.frame(contrast = d$contrast[1], lipid_class = d$lipid_class[1], tested = nrow(d), FDR_le_005 = sum(d$adj.P.Val <= 0.05), median_logFC = stats::median(d$logFC))))
  write_csv(classes, "lipid_class_descriptive_summary.csv")
  saveRDS(list(primary = primary, unweighted = unweighted, qc30 = qc30, supplied_scale = supplied_scale, separate_sex = separate), file.path(out, "limma_fits.rds"))
  capture.output(utils::sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  inputs <- c(input, paste0("data/RNA/scaled_RNA_", 1:3, ".csv"), "data/RNA/meta_data_df.csv", "data/metabolomics/Lipid_data_file_JPC_cleaned.csv", "data/Suppl_Table_1_OverviewData.xlsx", "results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv", "results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv", "scripts/workflows/11_cyp2c38_in_vivo_validation.R", "scripts/core/wp1_workflow.R", "scripts/core/cyp2c38_reference.R", "scripts/core/cyp2c38_models.R", "scripts/core/cyp2c38_plots.R", "scripts/core/cyp2c38_report.R", "tests/test_cyp2c38_in_vivo.R")
  write_csv(data.frame(path = inputs, sha256 = vapply(inputs, digest::digest, character(1), algo = "sha256", file = TRUE)), "input_and_code_SHA256.csv")
  provenance <- list(completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE), seed = seed, permutations = B, bootstraps = K, workers = workers,
                     primary_features = nrow(y), QC30_features = sum(ledger$keep_qc30), shared_reference_features = nrow(map),
                     reference_mice = nrow(reference$metadata), reference_bootstrap_strata = as.list(setNames(as.integer(table(strata_ref)), names(table(strata_ref)))),
                     checks = as.list(invariant), package_versions = setNames(lapply(required, function(p) as.character(utils::packageVersion(p))), required),
                     primary_model = "log2 supplied intensities; sex precision weights; robust/trend empirical Bayes",
                     normalization_provenance = "Not available; no additional normalization/imputation applied",
                     inferential_family = "Three primary concordance tests, two-sided stratified-permutation p with BH adjustment",
                     exchangeability = "Genotype labels exchangeable within sex; no batch/litter metadata available",
                     bootstrap_scope = "Whole animals; paired reference modalities; fixed feature selection; conditional percentile 95% intervals",
                     limitations = "Relative/scaled abundance; heterozygous KO; unknown processing/litter/batch provenance; indirect associations do not identify substrates or flux")
  jsonlite::write_json(provenance, file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
  source("scripts/core/cyp2c38_plots.R")
  wp1_make_plots(out, y, meta, ledger, primary$fit, differential, consistency, concordance)
  source("scripts/core/cyp2c38_report.R")
  wp1_write_report(out)
  cat("WP1 complete. Results:", normalizePath(out), "\n")
  print(discovery)
  print(concordance)
  invisible(list(out = out, differential = differential, concordance = concordance, discovery = discovery))
}


wp1_ulm_run <- function(args = commandArgs(trailingOnly = TRUE),
                        previous = "results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo",
                        intake = NULL) {
  for (file in c("cyp2c38_reference.R", "cyp2c38_sex_networks.R",
                 "cyp2c38_models.R", "cyp2c38_ulm_bootstrap.R")) {
    source(file.path("scripts/core", file))
  }
  get_option <- function(key, fallback) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (length(hit) > 1L) stop("Duplicate option: ", key)
    if (length(hit)) substring(hit, nchar(key) + 4L) else fallback
  }
  stopifnot(all(grepl("^--(bootstraps|permutations|seed|workers|output|calibration-datasets|calibration-inner|injection-datasets)=|^--validate-only$", args)))
  B <- as.integer(get_option("bootstraps", "4999"))
  P <- as.integer(get_option("permutations", "1999"))
  seed <- as.integer(get_option("seed", "20261007"))
  workers <- as.integer(get_option("workers", "3"))
  ncal <- as.integer(get_option("calibration-datasets", "120"))
  ninner <- as.integer(get_option("calibration-inner", "199"))
  ninjection <- as.integer(get_option("injection-datasets", "399"))
  stopifnot(B >= 19, P >= 19, workers >= 1, !is.na(seed), ncal >= 2, ninner >= 19, ninjection >= 3)
  out <- get_option("output", "results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/ulm_followup")
  if (is.null(intake)) {
    m <- utils::read.csv(file.path(previous, "log2_analysis_matrix.csv"), check.names = FALSE)
    y <- as.matrix(m[, -1, drop = FALSE]); rownames(y) <- m$feature_id
    meta <- utils::read.csv(file.path(previous, "sample_manifest.csv"), check.names = FALSE)
    ledger <- utils::read.csv(file.path(previous, "feature_annotation_and_QC_ledger.csv"), check.names = FALSE)
  } else {
    # Validation can start from raw intake without any completed WP1 outputs.
    y <- intake$y; meta <- intake$metadata; ledger <- intake$ledger
  }
  stopifnot(identical(colnames(y), meta$sample_id), all(is.finite(y)), nrow(y) == 1237L)
  reference <- wp1_full_reference()
  coefficients <- reference$coefficients
  coefficients$feature_id <- ledger$feature_id[match(coefficients$reference_lipid, ledger$canonical_lipid)]
  coefficients$feature_id[!coefficients$feature_id %in% rownames(y)] <- NA_character_
  coefficients$match_status <- ifelse(is.na(coefficients$feature_id), "unmatched_exact_molecular_annotation", "unique_exact")
  map <- coefficients[!is.na(coefficients$feature_id), , drop = FALSE]
  stopifnot(nrow(map) == 856L, sum(map$fixed_reference) == 552L, !anyDuplicated(map$feature_id))
  primary <- wp1_fit(y, meta)
  replay_difference <- NA_real_
  if (is.null(intake)) {
    saved_fit <- readRDS(file.path(previous, "limma_fits.rds"))$primary$fit
    replay_difference <- max(abs(primary$fit$t - saved_fit$t))
    stopifnot(replay_difference < 1e-10)
  }
  stopifnot(all(is.finite(primary$fit$t)))
  specs <- do.call(rbind, lapply(c("full", "fixed828"), function(universe) {
    data.frame(reference_universe = universe,
                contrast = c("pooled_KO_vs_WT", "female_KO_vs_WT", "male_KO_vs_WT", "female_KO_vs_WT", "male_KO_vs_WT"),
                reference_set = c("pooled", "female", "male", "pooled", "pooled"),
                primary_comparison = c(TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
  }))
  labels <- paste(specs$reference_universe, specs$contrast, specs$reference_set, sep = "__")
  if ("--validate-only" %in% args) {
    cat("Original/full sex-network replay, 856/552 one-to-one matches and finite KO model checks passed.\n")
    return(invisible(list(reference = reference, map = map, specs = specs)))
  }
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
    stop("ULM output directory is not empty; choose a fresh --output destination.")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(d, name) utils::write.csv(d, file.path(out, name), row.names = FALSE, na = "")
  write_csv(wp1_export_sex_networks(out, reference), "network_manifest.csv")
  write_csv(coefficients, "full_reference_and_identity_ledger.csv")
  write_csv(reference$metadata, "FCG_reference_mouse_manifest.csv")
  write_csv(meta, "KO_mouse_manifest.csv")
  evaluate <- function(tmat, ref = reference$coefficients) {
    vapply(seq_len(nrow(specs)), function(i) {
      selected <- if (specs$reference_universe[i] == "full") rep(TRUE, nrow(map)) else map$fixed_reference
      r <- ref[[paste0("r_", specs$reference_set[i])]][match(map$reference_lipid[selected], coefficients$reference_lipid)]
      wp1_ulm_stats(r, tmat[map$feature_id[selected], specs$contrast[i]])[["score"]]
    }, numeric(1L))
  }
  edges <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    a <- map[if (specs$reference_universe[i] == "full") rep(TRUE, nrow(map)) else map$fixed_reference, , drop = FALSE]
    data.frame(reference_universe = specs$reference_universe[i], contrast = specs$contrast[i], reference_set = specs$reference_set[i],
                feature_id = a$feature_id, canonical_lipid = a$reference_lipid,
                reference_r = a[[paste0("r_", specs$reference_set[i])]],
                t = primary$fit$t[a$feature_id, specs$contrast[i]], stringsAsFactors = FALSE)
  }))
  summary <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    a <- edges[edges$reference_universe == specs$reference_universe[i] & edges$contrast == specs$contrast[i] & edges$reference_set == specs$reference_set[i], ]
    u <- wp1_ulm_stats(a$reference_r, a$t)
    d <- as.data.frame(decoupleR::run_ulm(matrix(a$t, ncol = 1L, dimnames = list(a$canonical_lipid, specs$contrast[i])),
                                         data.frame(source = "Cyp2c38", target = a$canonical_lipid, mor = a$reference_r), minsize = 10L))
    stopifnot(abs(u[["score"]] - d$score) < 1e-10, abs(u[["p_value"]] - d$p_value) < 1e-10)
    cbind(specs[i, ], data.frame(n_matched_lipids = nrow(a), ULM_score = u[["score"]], ULM_p_value = u[["p_value"]],
                               reference_mice = if (specs$reference_set[i] == "pooled") 43L else if (specs$reference_set[i] == "female") 21L else 22L,
                               KO_mice = if (specs$contrast[i] == "pooled_KO_vs_WT") 32L else 16L))
  }))
  rownames(summary) <- NULL
  write_csv(edges, "ULM_matched_lipid_edges.csv")
  contrasts <- colnames(primary$fit$t)[1:3]
  components <- setNames(lapply(contrasts, function(ct) wp1_null_components(y, meta, contrast = ct, fit = primary)), contrasts)
  restricted <- setNames(lapply(contrasts, function(ct) wp1_null_components(y, meta, contrast = ct, fit = primary, restricted = TRUE)), contrasts)
  strata <- split(seq_along(reference$gene), reference$metadata$group)
  set.seed(seed)
  seeds_wild <- sample.int(.Machine$integer.max, B)
  seeds_restricted <- sample.int(.Machine$integer.max, min(B, 1999L))
  seeds_perm <- sample.int(.Machine$integer.max, P)
  run_draws <- function(seeds, fun, name, nstats = nrow(specs)) {
    chunks <- vector("list", length(seeds))
    for (start in seq.int(1L, length(seeds), by = 100L)) {
      ii <- start:min(start + 99L, length(seeds))
      chunks[ii] <- parallel::mclapply(seeds[ii], function(s) {
        set.seed(s)
        tryCatch(fun(), error = function(e) { warning(conditionMessage(e)); rep(NA_real_, nstats) })
      }, mc.cores = workers, mc.set.seed = FALSE)
      cat(name, max(ii), "/", length(seeds), "\n")
    }
    ans <- do.call(rbind, chunks)
    stopifnot(all(colMeans(is.finite(ans)) >= 0.99))
    ans
  }
  fit_null <- function(comp, z) {
    tt <- primary$fit$t
    for (ct in contrasts) tt[, ct] <- wp1_fit(wp1_wild_null(comp[[ct]], z), meta)$fit$t[, ct]
    tt
  }
  # Each tested contrast is null-restricted separately; other effects remain nuisance means.
  bank <- run_draws(seeds_wild, function() {
    z <- sample(c(-1, 1), ncol(y), replace = TRUE)
    tt <- fit_null(components, z)
    ib <- unlist(lapply(strata, function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
    rb <- wp1_reference_recompute(reference, ib)
    c(evaluate(tt), evaluate(tt, rb))
  }, "Contrast-null wild bootstrap", nstats = 2L * nrow(specs))
  wild <- bank[, seq_len(nrow(specs)), drop = FALSE]
  wild_reference <- bank[, nrow(specs) + seq_len(nrow(specs)), drop = FALSE]
  wild_restricted <- run_draws(seeds_restricted, function() {
    evaluate(fit_null(restricted, sample(c(-1, 1), ncol(y), replace = TRUE)))
  }, "Null-restricted residual sensitivity")
  permutation <- run_draws(seeds_perm, function() {
    mp <- meta
    for (sx in c("F", "M")) {
      ii <- which(mp$sex == sx)
      mp$genotype[ii] <- sample(mp$genotype[ii])
    }
    mp$group <- paste(mp$sex, mp$genotype, sep = "_")
    evaluate(wp1_fit(y, mp)$fit$t)
  }, "Within-sex KO-label permutation")
  write_bank <- function(z, seeds, file) {
    colnames(z) <- labels
    write_csv(data.frame(draw = seq_len(nrow(z)), seed = seeds, z, check.names = FALSE), file)
  }
  write_bank(wild, seeds_wild, "null_wild_draws.csv")
  write_bank(wild_reference, seeds_wild, "null_wild_reference_resampled_draws.csv")
  write_bank(wild_restricted, seeds_restricted, "null_wild_restricted_residual_draws.csv")
  write_bank(permutation, seeds_perm, "null_permutation_draws.csv")
  p_empirical <- function(z) vapply(seq_len(nrow(specs)), function(i) {
    (1 + sum(abs(z[, i]) >= abs(summary$ULM_score[i]), na.rm = TRUE)) / (1 + sum(is.finite(z[, i])))
  }, numeric(1L))
  summary$bootstrap_p_KO_only <- p_empirical(wild)
  summary$bootstrap_p_KO_and_reference <- p_empirical(wild_reference)
  summary$bootstrap_p_restricted_residuals <- p_empirical(wild_restricted)
  summary$permutation_p <- p_empirical(permutation)
  summary$bootstrap_valid_draws <- colSums(is.finite(wild))
  summary$combined_valid_draws <- colSums(is.finite(wild_reference))
  summary$bootstrap_MCSE <- sqrt(summary$bootstrap_p_KO_only * (1 - summary$bootstrap_p_KO_only) / summary$bootstrap_valid_draws)
  summary$bootstrap_null <- "Tested KO contrast equals zero for every lipid; retain nuisance contrasts"
  write_csv(summary, "ULM_summary.csv")
  diagnostics <- do.call(rbind, lapply(list(contrast_null_wild = wild, wild_plus_reference = wild_reference,
                                         restricted_residual_wild = wild_restricted, genotype_permutation = permutation), function(z) {
    cbind(specs, data.frame(draws = nrow(z), ULM_null_mean = colMeans(z, na.rm = TRUE), ULM_null_SD = apply(z, 2L, stats::sd, na.rm = TRUE),
                            nominal_ULM_p_le005_fraction = vapply(seq_len(nrow(specs)), function(i) mean(2 * stats::pt(-abs(z[, i]), summary$n_matched_lipids[i] - 2L) <= 0.05, na.rm = TRUE), numeric(1L))))
  }))
  diagnostics$method <- rep(c("contrast_null_wild", "wild_plus_reference", "restricted_residual_wild", "genotype_permutation"), each = nrow(specs))
  write_csv(diagnostics, "null_distribution_diagnostics.csv")
  old_boot <- utils::read.csv(file.path(previous, "combined_bootstrap_draws.csv"), check.names = FALSE)
  old_perm <- utils::read.csv(file.path(previous, "permutation_draws.csv"), check.names = FALSE)
  historical <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    if (specs$reference_universe[i] != "fixed828") return(NULL)
    key <- paste(specs$contrast[i], specs$reference_set[i], sep = "__")
    convert <- function(v) -v * sqrt((summary$n_matched_lipids[i] - 2L) / (1 - v^2))
    data.frame(contrast = specs$contrast[i], reference_set = specs$reference_set[i],
                previous_ordinary_combined_bootstrap_ULM_SD = stats::sd(convert(old_boot[[key]])),
                previous_permutation_null_ULM_SD = stats::sd(convert(old_perm[[key]])))
  }))
  write_csv(historical, "previous_resampling_variance_audit.csv")

  source("scripts/core/cyp2c38_ulm_calibration.R")
  wp1_calibration_audit(out, y, meta, reference, map, specs, seed + 1000L, workers,
                       n_datasets = ncal, n_inner = ninner, n_injection = ninjection)
  provenance <- list(completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE), seed = seed, bootstraps = B, permutations = P, workers = workers,
                     calibration_datasets = ncal, calibration_inner = ninner, injection_datasets = ninjection,
                     KO_model_replay_max_difference = replay_difference, network_archive_max_difference = reference$max_archive_difference,
                     full_network_features = 1254L, fixed_network_features = 828L, full_matches = 856L, fixed_matches = 552L,
                     model = "Unchanged 07b four-group sex-precision-weighted robust/trend limma on log2 supplied abundances",
                     null = "Tested KO contrast equals zero for every lipid, retaining other group-effect nuisance parameters",
                     wild_bootstrap = "Null-restricted means plus HC2-adjusted full-model residuals, shared Rademacher multiplier per mouse; refit weights and moderation",
                     reference_sensitivity = "Same KO null draws; independent paired FCG bootstrap within XXO/XYO/XXT/XYT",
                     nominal_ULM = "OLS-with-intercept slope t; nominal p treats lipid residual observations as independent",
                     scope = "No signed correlation score, confidence intervals or FDR in follow-up concordance reporting; previous outputs retained")
  jsonlite::write_json(provenance, file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
  inputs <- c(file.path(previous, c("log2_analysis_matrix.csv", "sample_manifest.csv", "feature_annotation_and_QC_ledger.csv", "limma_fits.rds", "combined_bootstrap_draws.csv", "permutation_draws.csv")),
              "data/metabolomics/Lipid_data_file_JPC_cleaned.csv", paste0("data/RNA/scaled_RNA_", 1:3, ".csv"), "data/RNA/meta_data_df.csv", "results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv",
              "results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv", "scripts/workflows/11_cyp2c38_in_vivo_validation.R", "scripts/core/wp1_workflow.R",
              file.path("scripts/core", c("cyp2c38_reference.R", "cyp2c38_models.R", "cyp2c38_sex_networks.R", "cyp2c38_ulm_bootstrap.R", "cyp2c38_ulm_calibration.R", "cyp2c38_ulm_presentation.R")))
  write_csv(data.frame(path = inputs, sha256 = vapply(inputs, digest::digest, character(1), algo = "sha256", file = TRUE)), "input_and_code_SHA256.csv")
  capture.output(utils::sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  source("scripts/core/cyp2c38_ulm_presentation.R")
  wp1_ulm_followup_plots(out, edges, summary)
  wp1_ulm_followup_report(out)
  cat("ULM follow-up complete:", normalizePath(out), "\n")
  print(summary[summary$primary_comparison, c("reference_universe", "contrast", "ULM_score", "ULM_p_value", "bootstrap_p_KO_only", "bootstrap_p_KO_and_reference")])
  invisible(summary)
}

# Complete work package: KO intake/differential analysis, then ULM and calibration.
wp1_package_run <- function(args = commandArgs(trailingOnly = TRUE)) {
  allowed <- "^--(permutations|bootstraps|base-permutations|base-bootstraps|base-seed|seed|workers|output|calibration-datasets|calibration-inner|injection-datasets)=|^--validate-only$"
  if (any(!grepl(allowed, args))) stop("Unknown Workflow 11 option")
  option <- function(key, default) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (length(hit) > 1L) stop("Duplicate option: ", key)
    if (length(hit)) substring(hit, nchar(key) + 4L) else default
  }
  out <- option("output", "results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo")
  base_args <- c(paste0("--output=", out),
    paste0("--permutations=", option("base-permutations", option("permutations", "1999"))),
    paste0("--bootstraps=", option("base-bootstraps", option("bootstraps", "999"))),
    paste0("--seed=", option("base-seed", "20261006")),
    paste0("--workers=", option("workers", "3")))
  follow_args <- c(paste0("--output=", file.path(out, "ulm_followup")),
    paste0("--permutations=", option("permutations", "1999")),
    paste0("--bootstraps=", option("bootstraps", "4999")),
    paste0("--seed=", option("seed", "20261007")),
    paste0("--workers=", option("workers", "3")),
    paste0("--calibration-datasets=", option("calibration-datasets", "120")),
    paste0("--calibration-inner=", option("calibration-inner", "199")),
    paste0("--injection-datasets=", option("injection-datasets", "399")))
  if ("--validate-only" %in% args) {
    cat("\n2. VALIDATE KO INTAKE AND FIXED FCG REFERENCE\n")
    intake <- wp1_run(c(base_args, "--validate-only"))
    cat("\n3. VALIDATE FULL/SEX-SPECIFIC NETWORKS AND ULM MODEL\n")
    follow <- wp1_ulm_run(c(follow_args, "--validate-only"), previous = out, intake = intake)
    cat("Workflow 11 validation passed without reading saved WP1 outputs.\n")
    return(invisible(list(intake = intake, ULM = follow)))
  }
  cat("\n2. KO INTAKE, DIFFERENTIAL ANALYSIS AND FIXED-REFERENCE CHECKS\n")
  base <- wp1_run(base_args)
  cat("\n3. FULL/SEX-SPECIFIC NETWORKS, ULM AND NULL-BOOTSTRAP CALIBRATION\n")
  follow <- wp1_ulm_run(follow_args, previous = out)
  jsonlite::write_json(list(workflow = "11_cyp2c38_in_vivo_validation",
    phases = c("KO intake and differential analysis", "native ULM and null-bootstrap calibration"),
    args = args, base_args = base_args, ULM_args = follow_args,
    current_concordance_report = "ulm_followup/README.md",
    finished_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(out, "workflow_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  cat("Workflow 11 complete: ", normalizePath(out), "\n", sep = "")
  invisible(list(base = base, ULM = follow))
}
