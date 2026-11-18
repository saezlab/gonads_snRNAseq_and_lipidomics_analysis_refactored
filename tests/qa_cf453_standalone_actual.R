# Independent actual-data checks for the standalone CF453 analysis.
# Run: Rscript tests/qa_cf453_standalone_actual.R ANALYSIS_RDS [QA_DIR]
# Focused table replay: Rscript tests/qa_cf453_standalone_actual.R
#                      --class-density-only ANALYSIS_RDS [QA_DIR]
# This verifier does not source the statistical or plotting implementation.

qa_cf453_standalone_actual <- function(analysis_file,
    qa_dir = "tests/outputs/cf453_standalone_actual", root = ".") {
  if (!file.exists(analysis_file)) stop("Analysis RDS is missing: ", analysis_file)
  if (dir.exists(file.path(root, ".local/R-library"))) {
    .libPaths(c(normalizePath(file.path(root, ".local/R-library")), .libPaths()))
  }
  a <- readRDS(analysis_file)
  required <- c("metadata", "modeled_rna", "profile_metadata", "design",
                "contrasts", "rho", "fit", "rna_results", "gene_eligibility",
                "candidates", "lipid_matrix", "lipid_ledger",
                "settings", "input_paths")
  stopifnot(all(required %in% names(a)))
  checks <- list()
  record <- function(label, difference = 0, detail = "") {
    if (!is.finite(difference) || difference > 1e-8) {
      stop("Independent CF453 QA failed: ", label, "; difference=", difference)
    }
    checks[[length(checks) + 1L]] <<- data.frame(
      check = label, max_abs_difference = difference, detail = detail,
      stringsAsFactors = FALSE)
  }
  same <- function(observed, expected, label, tolerance = 1e-8) {
    observed <- as.numeric(observed); expected <- as.numeric(expected)
    if (length(observed) != length(expected) ||
        !identical(is.na(observed), is.na(expected))) {
      stop("Independent CF453 QA length/missingness mismatch: ", label)
    }
    good <- is.finite(observed) & is.finite(expected)
    stopifnot(all(good | (is.na(observed) & is.na(expected)) |
                    (is.infinite(observed) & observed == expected)))
    difference <- if (any(good)) max(abs(observed[good] - expected[good])) else 0
    if (difference > tolerance) stop("Independent CF453 QA failed: ", label,
                                    "; difference=", difference)
    record(label, difference)
  }
  require_columns <- function(x, fields, label) {
    if (!is.data.frame(x) || !all(fields %in% names(x))) {
      stop("Missing independent-QA fields in ", label, ": ",
           paste(setdiff(fields, names(x)), collapse = ", "))
    }
  }
  finish <- function(scope, observed_r = NULL, omission_r = NULL) {
    dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
    summary <- do.call(rbind, checks)
    utils::write.csv(summary, file.path(qa_dir, "independent_checks.csv"), row.names = FALSE)
    saveRDS(list(checks = summary, scope = scope,
                 source_analysis = normalizePath(analysis_file),
                 observed_r = observed_r, omission_r = omission_r,
                 session = utils::sessionInfo()), file.path(qa_dir, "independent_qa.rds"))
    cat("Independent CF453 actual-data QA passed:", nrow(summary),
        "check groups; scope:", scope, "\n")
    invisible(summary)
  }
  m <- a$metadata; pm <- a$profile_metadata
  require_columns(m, c("mouse_id", "sex", "diet"), "metadata")
  require_columns(pm, c("profile_id", "mouse_id", "sex", "diet", "zone"),
                  "profile metadata")
  y <- as.matrix(a$modeled_rna); x <- as.matrix(a$design)
  cc <- as.matrix(a$contrasts); lipid <- as.matrix(a$lipid_matrix)
  stopifnot(nrow(m) == 28L, nrow(pm) == 56L,
            !anyDuplicated(m$mouse_id), !anyDuplicated(pm$profile_id),
            identical(colnames(y), pm$profile_id),
            identical(colnames(lipid), m$mouse_id),
            setequal(m$mouse_id, paste0("cf453-", 1:28)),
            all(table(m$sex, m$diet) == 7L),
            all(table(pm$mouse_id, pm$zone) == 1L),
            setequal(pm$zone, c("PP", "PC")),
            all(is.finite(y)), all(is.finite(lipid)),
            all(pm$sex == m$sex[match(pm$mouse_id, m$mouse_id)]),
            all(pm$diet == m$diet[match(pm$mouse_id, m$mouse_id)]))
  expected_sex <- rep(c("female", "male"), each = 14L)
  expected_diet <- rep(rep(c("SD", "WD"), each = 7L), 2L)
  roster <- m[match(paste0("cf453-", 1:28), m$mouse_id), , drop = FALSE]
  stopifnot(identical(roster$sex, expected_sex), identical(roster$diet, expected_diet))
  record("28 independent CF453 mice, seven per sex/diet; paired zones and lipid alignment")

  # Reconstruct the eight cell memberships without using model helpers.
  stopifnot(identical(dim(x), c(56L, 8L)), all(x %in% c(0, 1)),
            all(rowSums(x) == 1), all(colSums(x) == 7L),
            identical(rownames(cc), colnames(x)))
  cells <- do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    z <- unique(pm[x[, j] == 1, c("sex", "diet", "zone"), drop = FALSE])
    stopifnot(nrow(z) == 1L)
    z
  }))
  rownames(cells) <- colnames(x)
  cell <- function(sex, diet, zone) {
    z <- as.numeric(cells$sex == sex & cells$diet == diet & cells$zone == zone)
    stopifnot(sum(z) == 1L)
    z
  }
  zone_effect <- function(sex, diet) {
    cell(sex, diet, "PP") - cell(sex, diet, "PC")
  }
  interaction <- function(diet) zone_effect("male", diet) - zone_effect("female", diet)
  expected_primary <- cbind(interaction("SD"), interaction("WD"),
                            interaction("WD") - interaction("SD"))
  rr <- a$rna_results
  require_columns(rr, c("gene", "contrast", "effect", "se", "t", "p", "q", "family"),
                  "RNA results")
  primary_names <- unique(rr$contrast[grepl("primary", rr$family, ignore.case = TRUE)])
  stopifnot(length(primary_names) == 3L, all(primary_names %in% colnames(cc)))
  observed_vectors <- cc[, primary_names, drop = FALSE]
  vector_keys <- function(v) apply(v, 2L, paste, collapse = ",")
  stopifnot(setequal(vector_keys(observed_vectors), vector_keys(expected_primary)))
  record("Primary contrasts independently reconstruct signed sex-by-zone and three-way effects")
  expected_other <- list()
  for (sex in c("female", "male")) for (diet in c("SD", "WD")) {
    expected_other[[paste("zone", sex, diet, sep = "_")]] <- zone_effect(sex, diet)
  }
  for (sex in c("female", "male")) for (zone in c("PC", "PP")) {
    expected_other[[paste("diet", sex, zone, sep = "_")]] <-
      cell(sex, "WD", zone) - cell(sex, "SD", zone)
  }
  for (sex in c("female", "male")) {
    expected_other[[paste0("diet_zone_", sex)]] <- zone_effect(sex, "WD") - zone_effect(sex, "SD")
  }
  for (diet in c("SD", "WD")) for (zone in c("PC", "PP")) {
    expected_other[[paste("sex", diet, zone, sep = "_")]] <-
      cell("male", diet, zone) - cell("female", diet, zone)
  }
  for (contrast in names(expected_other)) {
    same(cc[, contrast], expected_other[[contrast]], paste("Secondary contrast sign", contrast))
  }
  for (family in unique(rr$family)) {
    j <- which(rr$family == family)
    same(rr$q[j], stats::p.adjust(rr$p[j], method = "BH"),
         paste("Complete RNA BH family", family))
    stopifnot(!anyDuplicated(paste(rr$gene[j], rr$contrast[j])))
  }

  # Base generalized least squares replay: duplicate-correlation covariance is
  # shared by genes, with no fictitious third observation from bulk lipids.
  rho <- as.numeric(a$rho)
  stopifnot(length(rho) == 1L, is.finite(rho), abs(rho) < 1,
            !is.null(a$fit$s2.post), !is.null(a$fit$df.total))
  v <- matrix(0, nrow(pm), nrow(pm))
  v[outer(pm$mouse_id, pm$mouse_id, "==")] <- rho
  diag(v) <- 1
  vi <- solve(v); info <- crossprod(x, vi %*% x); info_inv <- solve(info)
  beta <- t(info_inv %*% crossprod(x, vi %*% t(y)))
  residual <- y - beta %*% t(x)
  residual_df <- ncol(y) - qr(x)$rank
  sigma2 <- rowSums((residual %*% vi) * residual) / residual_df
  effect <- beta %*% cc
  unscaled <- sqrt(diag(t(cc) %*% info_inv %*% cc))
  se <- sqrt(a$fit$s2.post) %o% unscaled
  tt <- effect / se
  df_total <- rep_len(a$fit$df.total, nrow(y))
  same(df_total, pmin(rep_len(a$fit$df.prior, nrow(y)) + residual_df,
                      residual_df * nrow(y)), "Moderated RNA degrees of freedom and pooled bound")
  pp <- matrix(2 * stats::pt(-abs(as.vector(tt)),
                             df = rep(df_total, ncol(tt))),
               nrow(y), ncol(tt), dimnames = dimnames(effect))
  same(a$fit$coefficients, effect, "All modeled coefficients from independent base GLS")
  same(a$fit$sigma^2, sigma2, "All residual variances from independent base GLS")
  same(a$fit$stdev.unscaled, matrix(rep(unscaled, each = nrow(y)), nrow(y)),
       "All contrast standard errors before empirical Bayes")
  # Verify the moderation identity using independently reconstructed residuals.
  prior_df <- rep_len(a$fit$df.prior, nrow(y))
  prior_var <- rep_len(a$fit$s2.prior, nrow(y))
  posterior <- (prior_df * prior_var + residual_df * sigma2) / (prior_df + residual_df)
  posterior[is.infinite(prior_df)] <- prior_var[is.infinite(prior_df)]
  same(a$fit$s2.post, posterior, "Empirical Bayes posterior variances")
  same(a$fit$t, tt, "All moderated t statistics from independent GLS coefficients")
  same(a$fit$p.value, pp, "All moderated probabilities from independent statistics")
  row_index <- match(rr$gene, rownames(y)); col_index <- match(rr$contrast, colnames(cc))
  stopifnot(!anyNA(row_index), !anyNA(col_index),
            nrow(rr) == nrow(y) * ncol(cc))
  index <- cbind(row_index, col_index)
  same(rr$effect, effect[index], "RNA exported effects")
  same(rr$se, se[index], "RNA exported posterior standard errors")
  same(rr$t, tt[index], "RNA exported moderated t statistics")
  same(rr$p, pp[index], "RNA exported probabilities")
  if (all(c("low", "high") %in% names(rr))) {
    critical <- stats::qt(0.975, df_total[row_index])
    same(rr$low, rr$effect - critical * rr$se, "RNA effect interval lower limits")
    same(rr$high, rr$effect + critical * rr$se, "RNA effect interval upper limits")
  }

  # Independently reconstruct candidate memberships from the direct interaction
  # and corresponding simple zonation contrasts, retaining legitimate overlaps.
  candidate_expected <- list()
  for (sex in c("female", "male")) for (diet in c("SD", "WD")) {
    interaction_column <- which(vapply(seq_len(ncol(cc)), function(k)
      identical(as.numeric(cc[, k]), interaction(diet)), logical(1L)))
    zone_column <- which(vapply(seq_len(ncol(cc)), function(k)
      identical(as.numeric(cc[, k]), zone_effect(sex, diet)), logical(1L)))
    stopifnot(length(interaction_column) == 1L, length(zone_column) == 1L)
    ir <- rr[rr$contrast == colnames(cc)[interaction_column], , drop = FALSE]
    zr <- rr[rr$contrast == colnames(cc)[zone_column], , drop = FALSE]
    zr <- zr[match(ir$gene, zr$gene), , drop = FALSE]
    keep <- is.finite(ir$q) & is.finite(zr$q) & ir$q <= 0.05 & zr$q <= 0.05
    if (any(keep)) candidate_expected[[length(candidate_expected) + 1L]] <- data.frame(
      gene = ir$gene[keep], sex = sex, diet = diet,
      zone = ifelse(zr$effect[keep] > 0, "PP", "PC"),
      effect = zr$effect[keep], q_interaction = ir$q[keep], q_zonation = zr$q[keep],
      stringsAsFactors = FALSE)
  }
  ce <- if (length(candidate_expected)) do.call(rbind, candidate_expected) else
    data.frame(gene = character(), sex = character(), diet = character(),
               zone = character(), effect = numeric(), q_interaction = numeric(),
               q_zonation = numeric())
  ca <- a$candidates
  key <- function(d) paste(d$gene, d$sex, d$diet, d$zone, sep = "|")
  stopifnot(!anyDuplicated(key(ca)), setequal(key(ca), key(ce)))
  ca <- ca[match(key(ce), key(ca)), , drop = FALSE]
  for (field in c("effect", "q_interaction", "q_zonation")) {
    same(ca[[field]], ce[[field]], paste("Candidate membership", field))
  }
  record("All CF453-only RNA candidate memberships, signs and shared memberships")

  # Direct count-layer and lipid-source checks; no workflow/reference loader.
  paths <- unlist(a$input_paths, use.names = TRUE)
  stopifnot(length(paths) >= 2L, all(file.exists(paths)),
            !any(grepl("CF459|scaled_RNA|results/zonation|results/mofa|Lipid_data_file_JPC", paths)))
  zonation_path <- paths[grepl("zonation", names(paths))]
  lipid_path <- paths[grepl("lipid", names(paths))]
  stopifnot(length(zonation_path) == 1L, length(lipid_path) == 1L,
            requireNamespace("SeuratObject", quietly = TRUE))
  object <- readRDS(zonation_path)
  raw_counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  raw_data <- as.matrix(SeuratObject::LayerData(object, assay = "RNA", layer = "data"))
  ge <- a$gene_eligibility
  require_columns(ge, c("gene", "detected_profiles", "variable", "keep"), "gene eligibility")
  raw_counts <- raw_counts[ge$gene, pm$profile_id, drop = FALSE]
  raw_data <- raw_data[ge$gene, pm$profile_id, drop = FALSE]
  detected <- as.numeric(Matrix::rowSums(raw_counts > 0))
  variable <- apply(raw_data, 1L, stats::sd) > 1e-8
  same(ge$detected_profiles, detected, "All count-layer detection totals")
  stopifnot(identical(as.logical(ge$variable), unname(variable)),
            identical(as.logical(ge$keep), unname(detected >= 7L & variable)),
            identical(rownames(y), ge$gene[ge$keep]))
  same(y, raw_data[rownames(y), , drop = FALSE], "Modeled log-normalized RNA from supplied data layer")
  count_matrix <- as.matrix(raw_counts)
  expected_data <- log1p(sweep(count_matrix, 2L, colSums(count_matrix), "/") * 10000)
  same(raw_data, expected_data, "Independently reconstructed ln(1+10000*counts/library) transform")
  record("Gene eligibility is fixed before lipid outcomes; supplied log-normalized model scale")
  raw <- utils::read.csv(lipid_path, check.names = FALSE, stringsAsFactors = FALSE)
  ll <- a$lipid_ledger
  require_columns(ll, c("feature_id", "lipid_class", "keep_primary", "keep_broad", "keep_qc30",
                       "qc_cv_percent"), "lipid ledger")
  stopifnot(nrow(ll) == nrow(raw), identical(ll$feature_id, sprintf("WDSD_%04d", seq_len(nrow(raw)))))
  qc_columns <- grep("^QC[0-9]+$", names(raw), value = TRUE)
  biological <- raw[, m$lipid_sample, drop = FALSE]
  qc <- as.matrix(raw[, qc_columns, drop = FALSE])
  cv <- 100 * apply(qc, 1L, stats::sd) / rowMeans(qc)
  same(ll$qc_cv_percent, cv, "All independently calculated seven-injection QC CVs")
  stopifnot(identical(as.logical(ll$keep_qc30), ll$keep_broad & cv <= 30),
            identical(as.logical(ll$keep_primary), as.logical(ll$keep_qc30)),
            identical(rownames(lipid), ll$feature_id[ll$keep_primary]))
  expected_lipid <- log2(as.matrix(biological[ll$keep_primary, , drop = FALSE]))
  same(lipid, expected_lipid, "All primary mouse lipid log2 values and QC eligibility")
  broad_lipid <- as.matrix(a$lipid_broad)
  stopifnot(identical(rownames(broad_lipid), ll$feature_id[ll$keep_broad]),
            identical(colnames(broad_lipid), m$mouse_id))
  same(broad_lipid, log2(as.matrix(biological[ll$keep_broad, , drop = FALSE])),
       "All broad lipid sensitivity values")
  record("Input manifest contains local CF453 inputs without original-cohort statistical references")

  # Independent ordinary least squares and moderation replay for the 28-mouse
  # lipid model. Zone labels cannot enter this design.
  ld <- as.matrix(a$lipid_design); lc <- as.matrix(a$lipid_contrasts); lf <- a$lipid_fit
  stopifnot(identical(dim(ld), c(28L, 4L)), all(ld %in% c(0, 1)),
            all(rowSums(ld) == 1L), all(colSums(ld) == 7L))
  lcell <- function(sex, diet) {
    z <- numeric(ncol(ld))
    for (j in seq_len(ncol(ld))) {
      rows <- m[ld[, j] == 1L, , drop = FALSE]
      stopifnot(length(unique(rows$sex)) == 1L, length(unique(rows$diet)) == 1L)
      z[j] <- rows$sex[1L] == sex && rows$diet[1L] == diet
    }
    stopifnot(sum(z) == 1L)
    z
  }
  lipid_expected <- list(
    sex_SD = lcell("male", "SD") - lcell("female", "SD"),
    sex_WD = lcell("male", "WD") - lcell("female", "WD"),
    diet_female = lcell("female", "WD") - lcell("female", "SD"),
    diet_male = lcell("male", "WD") - lcell("male", "SD"))
  lipid_expected$sex_diet <- lipid_expected$diet_male - lipid_expected$diet_female
  for (contrast in names(lipid_expected)) {
    same(lc[, contrast], lipid_expected[[contrast]], paste("Lipid contrast sign", contrast))
  }
  linfo <- solve(crossprod(ld))
  lbeta <- t(linfo %*% crossprod(ld, t(lipid)))
  leffect <- lbeta %*% lc
  lresid <- lipid - lbeta %*% t(ld)
  lsigma2 <- rowSums(lresid^2) / 24
  lunscaled <- sqrt(diag(t(lc) %*% linfo %*% lc))
  lpost <- (lf$df.prior * lf$s2.prior + 24 * lsigma2) / (lf$df.prior + 24)
  if (length(lf$df.prior) == 1L && is.infinite(lf$df.prior)) lpost[] <- lf$s2.prior
  lse <- sqrt(lpost) %o% lunscaled
  ltt <- leffect / lse
  lpp <- matrix(2 * stats::pt(-abs(as.vector(ltt)),
    df = rep(rep_len(lf$df.total, nrow(lipid)), ncol(ltt))), nrow(lipid),
    dimnames = dimnames(leffect))
  same(lf$coefficients, leffect, "All 28-mouse lipid effects from independent OLS")
  same(lf$sigma^2, lsigma2, "All independent lipid residual variances")
  same(lf$s2.post, lpost, "All independent lipid moderated variances")
  lr <- a$lipid_results
  lindex <- cbind(match(lr$lipid, rownames(lipid)), match(lr$contrast, colnames(lc)))
  same(lr$effect, leffect[lindex], "All exported lipid effect sizes")
  same(lr$se, lse[lindex], "All exported lipid posterior standard errors")
  same(lr$t, ltt[lindex], "All exported lipid moderated statistics")
  same(lr$p, lpp[lindex], "All exported lipid moderated probabilities")
  for (family in unique(lr$family)) {
    j <- which(lr$family == family)
    same(lr$q[j], stats::p.adjust(lr$p[j], "BH"), paste("Lipid BH family", family))
  }
  if (is.null(a$edges)) return(finish("models and actual input provenance; network analysis not yet supplied"))

  # Every observed association and every one-mouse omission is recalculated by
  # base cor on explicitly matched animals. No molecular edge is a mouse.
  edges <- a$edges
  require_columns(edges, c("gene", "lipid", "lipid_class", "sex", "zone", "diet",
                          "r", "p", "q", "n", "loo_min", "loo_max", "loo_sign_fraction"),
                  "full association edge table")
  context_keys <- paste(edges$sex, edges$zone, edges$diet, sep = "|")
  edge_key <- function(d) paste(d$gene, d$lipid, d$sex, d$zone, d$diet, sep = "|")
  stopifnot(!anyDuplicated(edge_key(edges)), all(edges$n == 7L))
  expected_r <- expected_p <- rep(NA_real_, nrow(edges))
  expected_loo <- matrix(NA_real_, nrow(edges), 7L)
  for (sex in c("female", "male")) for (zone in c("PP", "PC")) {
    genes <- unique(ce$gene[ce$sex == sex & ce$zone == zone])
    lipids <- rownames(broad_lipid)
    valid_gene <- rep(TRUE, length(genes)); valid_lipid <- rep(TRUE, length(lipids))
    for (diet in c("SD", "WD")) {
      mice <- m$mouse_id[m$sex == sex & m$diet == diet]
      profile <- pm$profile_id[match(paste(mice, zone), paste(pm$mouse_id, pm$zone))]
      if (length(genes)) valid_gene <- valid_gene &
        apply(y[genes, profile, drop = FALSE], 1L, stats::sd) > 1e-8
      valid_lipid <- valid_lipid &
        apply(broad_lipid[lipids, mice, drop = FALSE], 1L, stats::sd) > 1e-8
    }
    genes <- sort(genes[valid_gene]); lipids <- lipids[valid_lipid]
    for (diet in c("SD", "WD")) {
      context <- paste(sex, zone, diet, sep = "|")
      j <- which(context_keys == context)
      mice <- m$mouse_id[m$sex == sex & m$diet == diet]
      profile <- pm$profile_id[match(paste(mice, zone), paste(pm$mouse_id, pm$zone))]
      rna <- y[genes, profile, drop = FALSE]
      lp <- broad_lipid[lipids, mice, drop = FALSE]
      if (!length(genes)) {
        stopifnot(!length(j))
        next
      }
      expected_pairs <- expand.grid(gene = genes, lipid = lipids, stringsAsFactors = FALSE)
      stopifnot(setequal(paste(edges$gene[j], edges$lipid[j]),
                        paste(expected_pairs$gene, expected_pairs$lipid)),
                nrow(expected_pairs) == length(j))
      observed <- suppressWarnings(stats::cor(t(rna), t(lp)))
      ij <- cbind(match(edges$gene[j], rownames(observed)), match(edges$lipid[j], colnames(observed)))
      expected_r[j] <- observed[ij]
      expected_p[j] <- 2 * stats::pt(-abs(expected_r[j] * sqrt(5 / (1 - expected_r[j]^2))), df = 5)
      for (omit in seq_len(7L)) {
        z <- suppressWarnings(stats::cor(t(rna[, -omit, drop = FALSE]),
                                         t(lp[, -omit, drop = FALSE])))
        expected_loo[j, omit] <- z[ij]
      }
      stopifnot(all(edges$lipid_class[j] == ll$lipid_class[match(edges$lipid[j], ll$feature_id)]))
      record(paste(context, "all mouse-paired coefficients and seven complete-mouse omissions"),
             detail = paste(length(j), "eligible molecular pairs from seven mice"))
    }
  }
  same(edges$r, expected_r, "All eight-context observed Pearson coefficients")
  same(edges$p, expected_p, "All edge nominal Pearson probabilities")
  expected_q <- rep(NA_real_, nrow(edges)); primary <- edges$primary_lipid
  expected_q[primary] <- stats::p.adjust(expected_p[primary], method = "BH")
  same(edges$q, expected_q, "Joint primary-lipid edge BH family across all evaluated contexts")
  safe_min <- function(z) if (any(is.finite(z))) min(z, na.rm = TRUE) else NA_real_
  safe_max <- function(z) if (any(is.finite(z))) max(z, na.rm = TRUE) else NA_real_
  same(edges$loo_min, apply(expected_loo, 1L, safe_min), "Every edge LOO minimum")
  same(edges$loo_max, apply(expected_loo, 1L, safe_max), "Every edge LOO maximum")
  valid_loo <- rowSums(is.finite(expected_loo))
  sign_fraction <- rowSums(sign(expected_loo) == sign(expected_r), na.rm = TRUE) / valid_loo
  sign_fraction[valid_loo == 0L] <- NA_real_
  same(edges$loo_valid, valid_loo, "Every edge valid omission count")
  same(edges$loo_sign_fraction, sign_fraction, "Every edge finite-omission sign fraction with explicit denominator")
  if (!is.null(a$edge_loo) && is.matrix(a$edge_loo)) {
    same(a$edge_loo, expected_loo, "Complete edge-by-omission coefficient bank")
  }
  if ("display" %in% names(edges)) {
    active <- paste(edges$gene, edges$sex, edges$diet, edges$zone, sep = "|") %in% key(ce)
    stopifnot(identical(as.logical(edges$active_candidate), active),
              identical(as.logical(edges$display), primary & active & is.finite(expected_r) & abs(expected_r) >= 0.65),
              identical(as.logical(edges$stable_display), edges$display & valid_loo == 7L & sign_fraction == 1))
    record("Exploratory display threshold retains both positive and negative associations")
  }

  # Class denominators include all measured quality-eligible CF453 species,
  # including species with no selected association. Each species appears once.
  membership <- function(quality = "keep_primary", threshold = 0.65) {
    ledger <- ll[ll[[quality]], , drop = FALSE]
    context_order <- c("female_PP", "female_PC", "male_PP", "male_PC")
    result <- list()
    for (diet in c("SD", "WD")) {
      bit <- matrix(FALSE, nrow(ledger), 4L,
                    dimnames = list(ledger$feature_id, context_order))
      select <- which(edges$diet == diet & edges$active_candidate &
                        is.finite(expected_r) & abs(expected_r) >= threshold)
      species <- match(edges$lipid[select], rownames(bit))
      context <- match(paste(edges$sex[select], edges$zone[select], sep = "_"), context_order)
      keep <- !is.na(species) & !is.na(context)
      if (any(keep)) bit[cbind(species[keep], context[keep])] <- TRUE
      n <- rowSums(bit)
      category <- rep("unassigned", nrow(bit))
      category[n > 1L] <- "shared"
      if (any(n == 1L)) category[n == 1L] <- context_order[max.col(bit[n == 1L, , drop = FALSE], ties.method = "first")]
      result[[diet]] <- data.frame(lipid = ledger$feature_id,
        lipid_class = ledger$lipid_class, diet = diet, category = category,
        membership_bits = apply(bit, 1L, function(v) paste(as.integer(v), collapse = "")),
        n_contexts = n, stringsAsFactors = FALSE)
    }
    do.call(rbind, result)
  }
  em <- membership()
  membership_key <- function(d) paste(d$lipid, d$diet, sep = "|")
  observed_membership <- a$membership[match(membership_key(em), membership_key(a$membership)), , drop = FALSE]
  stopifnot(nrow(em) == nrow(a$membership),
            identical(observed_membership$category, em$category),
            identical(observed_membership$membership_bits, em$membership_bits),
            identical(as.integer(observed_membership$n_contexts), as.integer(em$n_contexts)))
  class_key <- function(d) paste(d$diet, d$lipid_class, d$category, sep = "|")
  cs <- a$class_summary
  expected_n <- expected_total <- numeric(nrow(cs))
  for (j in seq_len(nrow(cs))) {
    in_class <- em$diet == cs$diet[j] & em$lipid_class == cs$lipid_class[j]
    expected_total[j] <- sum(in_class)
    expected_n[j] <- sum(in_class & em$category == cs$category[j])
  }
  stopifnot(!anyDuplicated(class_key(cs)), nrow(cs) == 12L * length(unique(em$lipid_class)))
  same(cs$n, expected_n, "All exclusive/shared/unassigned class membership counts")
  same(cs$total, expected_total, "All measured primary-quality class denominators")
  same(cs$proportion, expected_n / expected_total, "Every class membership proportion")
  for (quality in c("keep_primary", "keep_broad", "keep_fill")) {
    qm <- membership(quality)
    rows <- a$quality_sensitivity[a$quality_sensitivity$quality == quality, , drop = FALSE]
    for (j in seq_len(nrow(rows))) {
      k <- qm$diet == rows$diet[j] & qm$lipid_class == rows$lipid_class[j]
      stopifnot(rows$total[j] == sum(k), rows$assigned[j] == sum(k & qm$n_contexts > 0L))
      same(rows$proportion[j], sum(k & qm$n_contexts > 0L) / sum(k),
           paste("Lipid quality sensitivity", quality, rows$diet[j], rows$lipid_class[j]))
    }
  }
  for (threshold in unique(a$thresholds$threshold)) {
    tm <- membership(threshold = threshold)
    rows <- a$thresholds[a$thresholds$threshold == threshold, , drop = FALSE]
    for (j in seq_len(nrow(rows))) {
      k <- tm$diet == rows$diet[j] & tm$lipid_class == rows$lipid_class[j]
      stopifnot(rows$total[j] == sum(k), rows$assigned[j] == sum(k & tm$n_contexts > 0L))
      same(rows$proportion[j], sum(k & tm$n_contexts > 0L) / sum(k),
           paste("Display threshold sensitivity", threshold, rows$diet[j], rows$lipid_class[j]))
    }
  }
  record("Complete class assignment and quality/threshold sensitivities use measured CF453 denominators")

  # Reconstruct every ranking row before checking the displayed winners. This
  # prevents a few convenient examples from silently replacing the stated rule.
  pe <- edges[primary, , drop = FALSE]
  ranking <- unique(pe[, c("sex", "zone", "gene"), drop = FALSE])
  rank_key <- function(d) paste(d$sex, d$zone, d$gene, sep = "|")
  stable_pairs <- unique(pe[pe$stable_display, c("sex", "zone", "gene", "lipid"), drop = FALSE])
  degrees <- table(rank_key(stable_pairs))
  ranking$stable_degree <- as.integer(degrees[rank_key(ranking)])
  ranking$stable_degree[is.na(ranking$stable_degree)] <- 0L
  interaction_min <- tapply(ce$q_interaction, rank_key(ce), min)
  ranking$q_interaction <- as.numeric(interaction_min[rank_key(ranking)])
  ranking <- ranking[order(ranking$sex, ranking$zone, -ranking$stable_degree,
                           ranking$q_interaction, ranking$gene), , drop = FALSE]
  ranking$rank <- ave(ranking$stable_degree, paste(ranking$sex, ranking$zone), FUN = seq_along)
  stopifnot(identical(rank_key(a$selection_ledger), rank_key(ranking)))
  same(a$selection_ledger$stable_degree, ranking$stable_degree, "Complete representative stable-degree ranking")
  same(a$selection_ledger$q_interaction, ranking$q_interaction, "Complete representative RNA interaction tie-break")
  expected_reps <- ranking[ranking$rank == 1L & ranking$stable_degree > 0, , drop = FALSE]
  stopifnot(identical(rank_key(a$representatives), rank_key(expected_reps)))
  selected_edge_keys <- character()
  for (j in seq_len(nrow(expected_reps))) {
    rep <- expected_reps[j, ]
    ee <- pe[pe$sex == rep$sex & pe$zone == rep$zone & pe$gene == rep$gene, , drop = FALSE]
    stable <- unique(ee$lipid[ee$stable_display])
    strength <- vapply(stable, function(species) max(abs(ee$r[ee$lipid == species])), numeric(1L))
    partners <- head(stable[order(-strength, stable)], 6L)
    selected <- ee[ee$lipid %in% partners, , drop = FALSE]
    selected_edge_keys <- c(selected_edge_keys, edge_key(selected))
    actual <- a$representative_edges[rank_key(a$representative_edges) == rank_key(rep), , drop = FALSE]
    stopifnot(setequal(edge_key(actual), edge_key(selected)), nrow(actual) == 2L * length(partners))
    same(actual$partner_rank, match(actual$lipid, partners), paste("Representative partner rule", rank_key(rep)))
    example <- a$examples[rank_key(a$examples) == rank_key(rep), , drop = FALSE]
    stopifnot(nrow(example) == 14L, all(example$lipid == partners[1L]),
              all(table(example$diet) == 7L), !anyDuplicated(example$mouse_id))
    profiles <- pm$profile_id[match(paste(example$mouse_id, rep$zone), paste(pm$mouse_id, pm$zone))]
    same(example$rna, y[rep$gene, profiles], paste("Representative individual-mouse RNA values", rank_key(rep)))
    same(example$lipid_log2, lipid[partners[1L], example$mouse_id],
         paste("Representative individual-mouse lipid values", rank_key(rep)))
  }
  stopifnot(setequal(selected_edge_keys, edge_key(a$representative_edges)))
  record("Every representative and up-to-six partner set independently reproduces the CF453-only rule")

  if (!is.null(a$diet_edge_changes)) {
    de <- a$diet_edge_changes
    key_no_diet <- function(d) paste(d$sex, d$zone, d$gene, d$lipid, sep = "|")
    sd <- pe[pe$diet == "SD", , drop = FALSE]; wd <- pe[pe$diet == "WD", , drop = FALSE]
    sd <- sd[match(key_no_diet(de), key_no_diet(sd)), , drop = FALSE]
    wd <- wd[match(key_no_diet(de), key_no_diet(wd)), , drop = FALSE]
    same(de$r_SD, sd$r, "All matched-pair SD association values")
    same(de$r_WD, wd$r, "All matched-pair WD association values")
    same(de$delta_r, wd$r - sd$r, "Every descriptive WD-minus-SD association delta")
  }

  # The fixed-selection paired bank is shared across genes, lipids and zones.
  # Replay all draws for all representative edges plus deterministic extra
  # edges throughout each context; production streams its all-edge summaries.
  bank <- a$bootstrap_indices
  stopifnot(length(dim(bank)) == 3L, dim(bank)[1L] == 7L,
            dim(bank)[2L] == a$settings$bootstraps, dim(bank)[3L] == 4L,
            setequal(dimnames(bank)[[3L]], c("female_SD", "female_WD", "male_SD", "male_WD")),
            !anyNA(bank), all(bank == trunc(bank)), all(bank >= 1L & bank <= 7L))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(a$settings$seed)
  expected_bank <- array(NA_integer_, dim(bank), dimnames = dimnames(bank))
  for (group in dimnames(bank)[[3L]]) {
    expected_bank[, , group] <- replicate(a$settings$bootstraps, sample.int(7L, 7L, TRUE))
  }
  if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv)
  stopifnot(identical(bank, expected_bank))
  record("All four paired bootstrap index banks exactly reproduce the recorded seed")
  bootstrap_checked <- 0L
  for (context in unique(context_keys)) {
    rows <- which(context_keys == context & primary)
    if (!length(rows)) next
    chosen <- unique(c(rows[unique(round(seq(1L, length(rows), length.out = min(12L, length(rows)))))],
                       rows[edge_key(edges[rows, , drop = FALSE]) %in% selected_edge_keys]))
    ee <- edges[chosen, , drop = FALSE]
    mice <- m$mouse_id[m$sex == ee$sex[1L] & m$diet == ee$diet[1L]]
    profile <- pm$profile_id[match(paste(mice, ee$zone[1L]), paste(pm$mouse_id, pm$zone))]
    genes <- unique(ee$gene); lipids <- unique(ee$lipid)
    xx <- y[genes, profile, drop = FALSE]; yy <- lipid[lipids, mice, drop = FALSE]
    index <- cbind(match(ee$gene, genes), match(ee$lipid, lipids))
    draws <- bank[, , paste(ee$sex[1L], ee$diet[1L], sep = "_")]
    rr <- vapply(seq_len(ncol(draws)), function(b) {
      z <- suppressWarnings(stats::cor(t(xx[, draws[, b], drop = FALSE]),
                                       t(yy[, draws[, b], drop = FALSE])))
      z[index]
    }, numeric(nrow(ee)))
    if (is.null(dim(rr))) rr <- matrix(rr, nrow(ee))
    valid <- rowSums(is.finite(rr))
    interval <- t(apply(rr, 1L, stats::quantile, probs = c(0.025, 0.975),
                        names = FALSE, na.rm = TRUE))
    fraction <- rowMeans(sign(rr) == sign(ee$r), na.rm = TRUE)
    same(ee$bootstrap_valid, valid, paste(context, "all bootstrap validity counts for representative/independent sample"))
    same(ee$bootstrap_low, interval[, 1L], paste(context, "independent bootstrap lower percentiles"))
    same(ee$bootstrap_high, interval[, 2L], paste(context, "independent bootstrap upper percentiles"))
    same(ee$bootstrap_sign_fraction, fraction, paste(context, "independent bootstrap sign fractions"))
    if (!is.null(a$representative_bootstrap)) {
      br <- a$representative_bootstrap
      br <- br[br$sex == ee$sex[1L] & br$zone == ee$zone[1L] & br$diet == ee$diet[1L], , drop = FALSE]
      expected_pairs <- a$representative_edges
      expected_pairs <- expected_pairs[expected_pairs$sex == ee$sex[1L] &
        expected_pairs$zone == ee$zone[1L] & expected_pairs$diet == ee$diet[1L], , drop = FALSE]
      stopifnot(nrow(br) == nrow(expected_pairs) * a$settings$bootstraps,
        !anyDuplicated(paste(br$gene, br$lipid, br$draw)),
        all(br$draw %in% seq_len(a$settings$bootstraps)))
      br_index <- cbind(match(paste(br$gene, br$lipid), paste(ee$gene, ee$lipid)), br$draw)
      stopifnot(!anyNA(br_index))
      same(br$r, rr[br_index], paste(context, "complete exported representative bootstrap draw bank"))
    }
    bootstrap_checked <- bootstrap_checked + nrow(ee) * ncol(draws)
  }
  record("Paired bootstrap uses complete-mouse positions shared across both RNA zones and lipids",
         detail = paste(bootstrap_checked, "independently replayed edge/draw correlations"))

  if (!is.null(a$representative_bootstrap)) {
    br <- a$representative_bootstrap
    stopifnot(nrow(br) == nrow(a$representative_edges) * a$settings$bootstraps)
    sd <- br[br$diet == "SD", , drop = FALSE]; wd <- br[br$diet == "WD", , drop = FALSE]
    draw_key <- function(d) paste(d$sex, d$zone, d$gene, d$lipid, d$draw, sep = "|")
    wd <- wd[match(draw_key(sd), draw_key(wd)), , drop = FALSE]
    delta <- a$representative_delta_draws
    stopifnot(identical(draw_key(delta), draw_key(sd)))
    same(delta$delta_r, wd$r - sd$r, "Every representative independent-SD/WD paired bootstrap delta")
    ds <- a$representative_delta_summary
    edge_no_diet_key <- function(d) paste(d$sex, d$zone, d$gene, d$lipid, sep = "|")
    stopifnot(!anyDuplicated(edge_no_diet_key(ds)),
              setequal(edge_no_diet_key(ds), unique(edge_no_diet_key(delta))))
    for (j in seq_len(nrow(ds))) {
      v <- delta$delta_r[edge_no_diet_key(delta) == edge_no_diet_key(ds[j, , drop = FALSE])]
      observed <- a$diet_edge_changes$delta_r[
        edge_no_diet_key(a$diet_edge_changes) == edge_no_diet_key(ds[j, , drop = FALSE])]
      same(ds$delta_r[j], observed, paste("Representative observed delta", j))
      same(ds$low[j], stats::quantile(v, 0.025, na.rm = TRUE, names = FALSE),
           paste("Representative delta lower percentile", j))
      same(ds$high[j], stats::quantile(v, 0.975, na.rm = TRUE, names = FALSE),
           paste("Representative delta upper percentile", j))
      same(ds$positive_fraction[j], mean(v > 0, na.rm = TRUE),
           paste("Representative delta positive fraction", j))
      stopifnot(ds$valid[j] == sum(is.finite(v)), ds$requested[j] == length(v))
    }
    record("Complete representative delta intervals retain requested and undefined draw counts")
  }

  finish("complete standalone actual-data replay", expected_r, expected_loo)
}

qa_cf453_class_density_actual <- function(analysis_file,
    qa_dir = "tests/outputs/cf453_standalone_actual") {
  a <- readRDS(analysis_file)
  stopifnot(is.data.frame(a$class_pair_density), is.matrix(a$edge_loo),
            nrow(a$edge_loo) == nrow(a$edges), ncol(a$edge_loo) == 7L)
  d <- a$class_pair_density; e <- a$edges
  checks <- list()
  record <- function(label, observed, expected) {
    stopifnot(length(observed) == length(expected),
              identical(is.na(as.numeric(observed)), is.na(as.numeric(expected))))
    difference <- if (length(observed)) max(abs(observed - expected), na.rm = TRUE) else 0
    if (!is.finite(difference) || difference > 1e-8) {
      stop("Independent CF453 density QA failed: ", label, "; difference=", difference)
    }
    checks[[length(checks) + 1L]] <<- data.frame(check = label,
      max_abs_difference = difference, stringsAsFactors = FALSE)
  }
  fields <- c("diet", "sex", "zone", "lipid_class")
  key <- function(x) do.call(paste, c(x[, fields, drop = FALSE], sep = "|"))
  stopifnot(!anyDuplicated(key(d)))
  primary <- which(e$primary_lipid)
  primary_keys <- key(e[primary, , drop = FALSE])
  stopifnot(setequal(unique(primary_keys), key(d)))
  groups <- split(primary, primary_keys)
  valid <- rowSums(is.finite(a$edge_loo))
  signs <- rowSums(sign(a$edge_loo) == sign(e$r), na.rm = TRUE)
  stable <- valid == 7L & signs == 7L
  pair_key <- paste(e$sex, e$zone, e$gene, e$lipid, sep = "|")
  for (j in seq_len(nrow(d))) {
    rows <- groups[[key(d[j, , drop = FALSE])]]
    opposite <- d[j, , drop = FALSE]
    opposite$diet <- if (opposite$diet == "SD") "WD" else "SD"
    other <- groups[[key(opposite)]]
    stopifnot(length(rows) == length(other), setequal(pair_key[rows], pair_key[other]),
              length(rows) > 0L, !anyDuplicated(pair_key[rows]))
    ar <- abs(e$r[rows]); qualified <- ar >= 0.65
    strong_stable <- qualified & stable[rows]
    expected <- c(total_pairs = length(rows), display_pairs = sum(qualified),
      stable_pairs = sum(strong_stable), fdr05_pairs = sum(e$q[rows] <= .05, na.rm = TRUE),
      mean_abs_r = mean(ar), median_abs_r = stats::median(ar),
      display_fraction = sum(qualified) / length(rows),
      stable_fraction = sum(strong_stable) / length(rows),
      fdr05_fraction = sum(e$q[rows] <= .05, na.rm = TRUE) / length(rows))
    for (field in names(expected)) {
      record(paste(key(d[j, , drop = FALSE]), field), d[[field]][j], expected[[field]])
    }
    # The common species/gene Cartesian product is a measured/tested-pair
    # denominator, rather than a thresholded subset or an extra mouse count.
    record(paste(key(d[j, , drop = FALSE]), "complete common gene-by-species denominator"),
           length(rows), length(unique(e$gene[rows])) * length(unique(e$lipid[rows])))
  }
  full_qa <- file.path(qa_dir, "independent_checks.csv")
  if (!file.exists(full_qa)) {
    stop("Run complete actual-data QA before the focused density-table audit.")
  }
  prior <- utils::read.csv(full_qa, stringsAsFactors = FALSE)
  stopifnot(any(prior$check == "All eight-context observed Pearson coefficients"),
            any(prior$check == "Complete edge-by-omission coefficient bank"))
  summary <- do.call(rbind, checks)
  utils::write.csv(summary, file.path(qa_dir, "independent_class_pair_density.csv"), row.names = FALSE)
  saveRDS(list(checks = summary, source_analysis = normalizePath(analysis_file),
               groups_checked = nrow(d), all_previous_edge_checks = full_qa),
          file.path(qa_dir, "independent_class_pair_density.rds"))
  cat("Independent CF453 class-pair density QA passed:", nrow(summary),
      "checks across", nrow(d), "class/context/diet groups.\n")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || identical(args, "--help")) {
    cat("Usage: Rscript tests/qa_cf453_standalone_actual.R ANALYSIS_RDS [QA_DIR]\n")
    cat("       Rscript tests/qa_cf453_standalone_actual.R --class-density-only ANALYSIS_RDS [QA_DIR]\n")
  } else if (args[1L] == "--class-density-only") {
    stopifnot(length(args) >= 2L)
    qa_cf453_class_density_actual(args[2L],
      if (length(args) >= 3L) args[3L] else "tests/outputs/cf453_standalone_actual")
  } else {
    qa_cf453_standalone_actual(args[1L],
      if (length(args) >= 2L) args[2L] else "tests/outputs/cf453_standalone_actual")
  }
}
