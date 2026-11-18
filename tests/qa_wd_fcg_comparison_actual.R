source("scripts/config/analysis_release.R")
# Independent actual-data replay for the WP3 original-male-cohort extension.
# Run: Rscript tests/qa_wd_fcg_comparison_actual.R COMPARISON_RDS CONTEXTS_RDS
#      ORIGINAL_ANALYSIS_RDS [QA_DIR]
# Source loading defines the verifier only. It does not change an analysis.

qa_wd_fcg_comparison_actual <- function(comparison_file, contexts_file,
                                       original_analysis_file,
                                       qa_dir = "tests/outputs/wp3_fcg_actual_qa",
                                       root = ".") {
  stopifnot(file.exists(comparison_file), file.exists(contexts_file),
            file.exists(original_analysis_file))
  fit <- readRDS(comparison_file)
  contexts <- readRDS(contexts_file)
  old <- readRDS(original_analysis_file)
  if (!is.null(contexts$contexts)) contexts <- contexts$contexts
  expected <- c("male_periportal", "male_pericentral")
  stopifnot(identical(names(contexts), expected),
            all(c("edges", "inference", "FCG_paired_index_matrix",
                  "original_WP3_index_matrices", "exact_CF459_indices",
                  "omission_ledger", "paired_draws") %in% names(fit)))
  checks <- list()
  record <- function(check, difference = 0, detail = "") {
    if (!is.finite(difference) || difference > 1e-9) {
      stop("Independent WP3 original-cohort QA failed: ", check,
           "; difference=", difference)
    }
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, max_abs_difference = difference, detail = detail,
      stringsAsFactors = FALSE)
  }
  compare <- function(a, b, check) {
    a <- as.numeric(a); b <- as.numeric(b)
    stopifnot(length(a) == length(b))
    if (!identical(is.na(a), is.na(b))) {
      stop("Independent QA missingness mismatch: ", check, "; positions ",
           paste(head(which(is.na(a) != is.na(b)), 12L), collapse = ","))
    }
    valid <- is.finite(a) & is.finite(b)
    stopifnot(all(valid | (is.na(a) & is.na(b)) |
                    (is.infinite(a) & a == b)))
    record(check, if (any(valid)) max(abs(a[valid] - b[valid])) else 0)
  }
  read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                           stringsAsFactors = FALSE)
  # Reconstruct FCG from original CSVs, without the analysis/reference helpers.
  rt <- read_csv(analysis_path("historical_fcg_zonation_csv",root))
  meta <- read_csv(analysis_path("historical_fcg_zonation_metadata",root))
  gt <- read_csv(file.path(root, "results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv"))
  lt <- read_csv(file.path(root, "data/metabolomics/Lipid_data_file_JPC_cleaned.csv"))
  stopifnot(names(rt)[1L] == "gene", nrow(meta) == 86L,
            nrow(gt) == 43L, ncol(lt) == 1255L,
            !anyDuplicated(rt$gene), !anyDuplicated(meta$sample))
  raw_lipid <- t(as.matrix(lt[, -1L, drop = FALSE]))
  lipid_mouse <- gsub("-", "_", gsub("-0", "-", lt[[1L]], fixed = TRUE))
  dimnames(raw_lipid) <- list(names(lt)[-1L], lipid_mouse)
  stopifnot(all(is.finite(raw_lipid)), all(raw_lipid > 0),
            !anyDuplicated(lipid_mouse))
  gt$mouse_id <- gsub("-", "_", gt$sequencing_name)
  fcg <- list(); first_ids <- first_groups <- NULL
  for (context in expected) {
    zone <- if (context == "male_periportal") "Periportal" else "Pericentral"
    columns <- names(rt)[grepl(paste0("_", zone, "_Testes$"), names(rt))]
    m <- meta[match(columns, meta$sample), , drop = FALSE]
    ids <- m$sample_name
    g <- gt[match(ids, gt$mouse_id), , drop = FALSE]
    stopifnot(length(columns) == 22L, !anyNA(m$sample),
              !anyDuplicated(ids), identical(m$gonadal_type, rep("Testes", 22L)),
              setequal(g$group, c("XXT", "XYT")),
              all(table(g$group) == 11L), !anyNA(g$mouse_id),
              identical(ids, sub(paste0("_", zone, "_Testes$"), "", columns)),
              all(ids %in% lipid_mouse))
    if (is.null(first_ids)) {
      first_ids <- ids; first_groups <- g$group
    } else stopifnot(identical(ids, first_ids), identical(g$group, first_groups))
    e <- contexts[[context]]$edges
    genes <- unique(e$gene); lipids <- unique(e$reference_lipid)
    stopifnot(length(lipids) == 610L, all(genes %in% rt$gene),
              all(lipids %in% rownames(raw_lipid)))
    rna <- as.matrix(rt[match(genes, rt$gene), columns, drop = FALSE])
    dimnames(rna) <- list(genes, ids)
    lipid <- log2(raw_lipid[lipids, ids, drop = FALSE])
    fcg[[context]] <- list(rna = rna, lipid = lipid, mouse_id = ids,
                           group = g$group)
    record(paste(context, "independent original FCG metadata and pairing"),
           detail = "22 paired testes-bearing mice; 11 XXT plus 11 XYT; no WT-only substitution")
  }
  rm(rt, lt, raw_lipid)
  # This replay deliberately uses base cor and lm, without WP2/WP3 helpers.
  edge_cor <- function(g, l, e, gi = NULL, li = gi) {
    if (is.null(gi)) gi <- seq_len(ncol(g))
    if (is.null(li)) li <- seq_len(ncol(l))
    cm <- suppressWarnings(stats::cor(t(g[, gi, drop = FALSE]),
                                      t(l[, li, drop = FALSE])))
    as.numeric(cm[cbind(match(e$gene, rownames(cm)),
                         match(e$reference_lipid, colnames(cm)))])
  }
  wd_cor <- function(x, cohort, e, gi = NULL, li = gi) {
    edge_cor(x[[paste0("rna_", cohort)]], x[[paste0("lipid_", cohort)]],
             e, gi, li)
  }
  vec_cor <- function(a, b) {
    if (length(a) < 3L || any(!is.finite(a)) || any(!is.finite(b))) return(NA_real_)
    # Repeated-mouse draws can make every edge exactly +1 or -1. The base
    # coefficient matrix can differ by a few ulps even then. Treat relative
    # roundoff around a constant vector as zero variance; this does not impose
    # an absolute amplitude threshold on weak, genuinely varying coefficients.
    constant_roundoff <- function(v) {
      diff(range(v)) <= 32 * .Machine$double.eps * max(abs(v))
    }
    if (constant_roundoff(a) || constant_roundoff(b)) return(NA_real_)
    suppressWarnings(stats::cor(a, b))
  }
  triplet <- function(a, b, c) {
    r1 <- vec_cor(a, c); r2 <- vec_cor(b, c)
    c(r_FCG_CF459 = r1, r_CF453_CF459 = r2,
      r_FCG_CF453 = vec_cor(a, b), delta_r = r1 - r2)
  }
  scope_indices <- function(e, rec) {
    keep <- e$eligible
    if (rec$universe == "displayed") keep <- keep & e$primary_edge
    if (rec$scope == "gene") keep <- keep & e$gene == rec$gene
    which(keep)
  }
  ix <- fit$original_WP3_index_matrices
  stopifnot(identical(ix, old$resampling_index_matrices))
  record("Original WP3 mouse-index banks remain exactly unchanged")
  # An exactly factorized empirical source-ID contingency proves covariance
  # zero for every possible pair of source features. This certifies undefined
  # coefficient-vector correlations, rather than correlating floating noise.
  independent_margin <- vapply(seq_len(nrow(ix$CF459_null_rna)), function(draw) {
    tab <- table(factor(ix$CF459_null_rna[draw, ], levels = 1:4),
                 factor(ix$CF459_null_lipid[draw, ], levels = 1:4))
    all(ncol(ix$CF459_null_rna) * tab == outer(rowSums(tab), colSums(tab)))
  }, logical(1L))
  record("Independent integer contingency identifies exactly zero-covariance null draws",
         detail = paste(sum(independent_margin), "of", length(independent_margin),
                        "empirical source-ID contingency tables factorize exactly"))
  stopifnot(all(c("historical_null_correction_summary", "historical_null_corrected_draws") %in% names(fit)))
  for (k in seq_len(nrow(old$inference))) {
    rec <- old$inference[k, ]; key <- rec$specification
    original <- old$bootstrap_null_draws[[key]]
    corrected <- original; corrected[independent_margin] <- NA_real_
    row <- fit$historical_null_correction_summary[
      fit$historical_null_correction_summary$specification == key, , drop = FALSE]
    stopifnot(nrow(row) == 1L)
    compare(corrected, fit$historical_null_corrected_draws[[key]],
            paste(key, "complete corrected historical CF453-anchor bank"))
    compare(row$observed_CF453_CF459_r, rec$coefficient_vector_r,
            paste(key, "historical observed coefficient remains unchanged"))
    compare(row$original_null_bootstrap_p, rec$null_bootstrap_p,
            paste(key, "historical raw p remains recorded"))
    valid <- is.finite(corrected); n <- length(corrected)
    exceed <- sum(abs(corrected[valid]) >= abs(rec$coefficient_vector_r) - 1e-12)
    p <- if (any(valid)) (1 + exceed) / (1 + sum(valid)) else NA_real_
    compare(row$corrected_null_bootstrap_p, p,
            paste(key, "independent corrected historical bootstrap tail"))
    compare(row$corrected_completion_lower, (1 + exceed) / (1 + n),
            paste(key, "corrected historical lower completion"))
    compare(row$corrected_completion_upper, (1 + exceed + sum(!valid)) / (1 + n),
            paste(key, "corrected historical upper completion"))
    stopifnot(row$original_valid_draws == sum(is.finite(original)),
              row$original_invalid_draws == sum(!is.finite(original)),
              row$corrected_valid_draws == sum(valid),
              row$corrected_invalid_draws == sum(!valid),
              row$newly_invalid_draws == sum(is.finite(original) & !valid))
  }
  pi <- fit$exact_CF459_indices
  stopifnot(nrow(pi) == 24L * 4L)
  permutations <- matrix(pi$source_position, ncol = 4L, byrow = TRUE)
  stopifnot(identical(as.integer(permutations[1L, ]), 1:4),
            !anyDuplicated(apply(permutations, 1L, paste, collapse = ",")),
            all(apply(permutations, 1L, function(z) identical(sort(as.integer(z)), 1:4))))
  record("CF459 exact bank contains all 24 unique pairings including identity")
  fi <- fit$FCG_paired_index_matrix
  stopifnot(is.matrix(fi), ncol(fi) == 22L,
            nrow(fi) == nrow(ix$CF459_paired),
            !anyNA(fi), all(fi >= 1L & fi <= 22L), all(fi == trunc(fi)))
  for (group in c("XXT", "XYT")) {
    columns <- which(first_groups == group)
    stopifnot(all(first_groups[fi[, columns]] == group))
  }
  ledger <- fit$FCG_paired_indices
  stopifnot(identical(as.integer(ledger$source_position), as.integer(as.vector(t(fi)))),
            identical(ledger$source_mouse, first_ids[ledger$source_position]),
            identical(ledger$source_group, first_groups[ledger$source_position]),
            all(ledger$source_group == ledger$position_group))
  record("FCG paired bank preserves both 11-mouse chromosome strata and shared zone order")
  if (!is.null(fit$metadata$FCG_paired_seed)) {
    stopifnot(fit$metadata$FCG_paired_seed != old$metadata$seed)
    # Independent seeded regeneration, with caller RNG state restored.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(fit$metadata$FCG_paired_seed)
    replica <- matrix(NA_integer_, nrow(fi), 22L)
    for (group in c("XXT", "XYT")) {
      positions <- which(first_groups == group)
      replica[, positions] <- matrix(positions[sample.int(length(positions),
         length(positions) * nrow(fi), replace = TRUE)], nrow(fi), length(positions))
    }
    if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv)
    stopifnot(identical(fi, replica))
    record("FCG paired bank independently replays its distinct recorded seed")
  }
  omission <- fit$omission_ledger
  stopifnot(nrow(omission) == 61L,
            sum(omission$method == "joint_omission") == 28L,
            sum(omission$method == "CF453_only") == 7L,
            sum(omission$method == "CF459_only") == 4L,
            sum(omission$method == "FCG_only") == 22L,
            !anyDuplicated(omission$draw))
  record("Omission bank includes original 39 configurations plus every FCG-only omission")
  for (context in expected) {
    x <- contexts[[context]]; f <- fcg[[context]]
    e <- fit$edges[fit$edges$context == context, , drop = FALSE]
    stopifnot(identical(e$edge_id, x$edges$edge_id), all(e$eligible),
              identical(e$primary_edge, x$edges$primary_edge),
              identical(e$reference_lipid, x$edges$reference_lipid),
              all(e$n_FCG == 22L))
    a <- edge_cor(f$rna, f$lipid, e)
    b <- wd_cor(x, "CF453", e); c <- wd_cor(x, "CF459", e)
    compare(a, e$r_FCG_SD, paste(context, "all original FCG observed edge coefficients"))
    compare(a, e$r_reference, paste(context, "complete archived FCG replay"))
    compare(b, e$r_CF453_WD, paste(context, "all unchanged CF453 edge coefficients"))
    compare(c, e$r_CF459_WD, paste(context, "all unchanged CF459 edge coefficients"))
    old_edges <- old$edges[old$edges$context == context, , drop = FALSE]
    stopifnot(identical(e$edge_id, old_edges$edge_id))
    compare(b, old_edges$r_CF453_WD, paste(context, "old WP3 CF453 preservation"))
    compare(c, old_edges$r_CF459_WD, paste(context, "old WP3 CF459 preservation"))
    rows <- fit$inference[fit$inference$context == context, , drop = FALSE]
    specs <- lapply(seq_len(nrow(rows)), function(k) scope_indices(e, rows[k, , drop = FALSE]))
    names(specs) <- rows$specification
    for (k in seq_len(nrow(rows))) {
      rec <- rows[k, ]; j <- specs[[k]]; key <- rec$specification
      z <- triplet(a[j], b[j], c[j])
      for (field in names(z)) compare(rec[[field]], z[field], paste(key, "observed", field))
      model <- summary(stats::lm(c[j] ~ a[j]))$coefficients[2L, ]
      compare(rec$ULM_FCG_CF459, model["t value"], paste(key, "independent OLS ULM"))
      compare(rec$ULM_p_FCG_CF459, model["Pr(>|t|)"], paste(key, "independent OLS nominal p"))
      stopifnot(rec$n_edges == length(j), rec$n_FCG == 22L,
                rec$n_CF453 == 7L, rec$n_CF459 == 4L)
    }
    # Replay every exact and null draw on the same frozen edges.
    for (method in c("permutation_CF459", "bootstrap_null")) {
      n <- nrow(fit[[paste0(method, "_draws")]])
      direct <- matrix(NA_real_, n, length(specs), dimnames = list(NULL, names(specs)))
      for (draw in seq_len(n)) {
        v <- if (method == "permutation_CF459") wd_cor(x, "CF459", e, 1:4, permutations[draw, ]) else wd_cor(x, "CF459", e, ix$CF459_null_rna[draw, ], ix$CF459_null_lipid[draw, ])
        if (method == "bootstrap_null" && independent_margin[draw]) {
          v[is.finite(v)] <- 0
        }
        for (k in seq_along(specs)) direct[draw, k] <- vec_cor(a[specs[[k]]], v[specs[[k]]])
      }
      for (k in seq_len(nrow(rows))) {
        rec <- rows[k, ]; key <- rec$specification; v <- direct[, k]
        compare(v, fit[[paste0(method, "_draws")]][[key]], paste(key, method, "every direct Pearson bank draw"))
        ulm <- v * sqrt((length(specs[[k]]) - 2) / (1 - v * v))
        compare(ulm, fit[[paste0(method, "_ULM_draws")]][[key]],
                paste(key, method, "complete-bank ULM identity"))
        valid <- is.finite(v); exceed <- sum(abs(v[valid]) >= abs(rec$r_FCG_CF459) - 1e-12)
        plus <- if (method == "permutation_CF459") 0 else 1
        p <- if (any(valid)) (plus + exceed) / (plus + sum(valid)) else NA_real_
        pfield <- if (method == "permutation_CF459") "pairing_CF459_given_FCG_p" else "null_bootstrap_p"
        compare(rec[[pfield]], p, paste(key, method, "independent complete-bank tail"))
        if (method == "bootstrap_null") {
          stopifnot(rec$bootstrap_valid_draws == sum(valid), rec$bootstrap_invalid_draws == sum(!valid))
          compare(rec$bootstrap_completion_lower, (1 + exceed) / (1 + n), paste(key, "undefined-mass lower completion"))
          compare(rec$bootstrap_completion_upper, (1 + exceed + sum(!valid)) / (1 + n), paste(key, "undefined-mass upper completion"))
        } else stopifnot(rec$permutation_CF459_valid_draws == 24L)
      }
      message("Independent FCG QA replayed ", context, " ", method, " complete bank.")
    }
    # One shared CF459 draw enters both anchor comparisons; reconstruct all draws.
    pb <- nrow(fi)
    direct_paired <- array(NA_real_, c(pb, length(specs), 4L),
      dimnames = list(NULL, names(specs), names(triplet(a, b, c))))
    for (draw in seq_len(pb)) {
      aa <- edge_cor(f$rna, f$lipid, e, fi[draw, ])
      bb <- wd_cor(x, "CF453", e, ix$CF453_paired[draw, ])
      cc <- wd_cor(x, "CF459", e, ix$CF459_paired[draw, ])
      for (k in seq_along(specs)) {
        j <- specs[[k]]; direct_paired[draw, k, ] <- triplet(aa[j], bb[j], cc[j])
      }
    }
    for (k in seq_len(nrow(rows))) {
      rec <- rows[k, ]; key <- rec$specification
      p <- fit$paired_draws[fit$paired_draws$specification == key, , drop = FALSE]
      p <- p[order(p$draw), , drop = FALSE]
      stopifnot(identical(as.integer(p$draw), seq_len(pb)), nrow(p) == pb)
      for (field in dimnames(direct_paired)[[3L]]) {
        compare(direct_paired[, k, field], p[[field]], paste(key, "every independent paired draw", field))
      }
      v <- direct_paired[, k, "delta_r"]; valid <- is.finite(v)
      stopifnot(rec$paired_delta_valid_draws == sum(valid), rec$paired_delta_invalid_draws == sum(!valid))
      compare(rec$paired_delta_positive_fraction, if (any(valid)) mean(v[valid] > 0) else NA_real_, paste(key, "delta positive fraction on defined draws"))
      compare(direct_paired[, k, "r_CF453_CF459"], old$paired_stability_draws[[key]], paste(key, "original WP3 paired-bank preservation"))
    }
    message("Independent FCG QA replayed ", context, " complete paired bank.")
    for (draw in seq_len(nrow(omission))) {
      o <- omission[draw, ]
      fa <- if (is.na(o$omit_FCG)) 1:22 else setdiff(1:22, o$omit_FCG)
      fb <- if (is.na(o$omit_CF453)) 1:7 else setdiff(1:7, o$omit_CF453)
      fc <- if (is.na(o$omit_CF459)) 1:4 else setdiff(1:4, o$omit_CF459)
      for (cohort in c("FCG", "CF453", "CF459")) {
        kept <- switch(cohort, FCG = fa, CF453 = fb, CF459 = fc)
        ids <- if (cohort == "FCG") f$mouse_id else colnames(x[[paste0("rna_", cohort)]])
        oi <- fit$omission_indices[fit$omission_indices$draw == o$draw &
                                   fit$omission_indices$cohort == cohort, , drop = FALSE]
        stopifnot(identical(as.integer(oi$source_position), as.integer(kept)),
                  identical(oi$source_mouse, ids[kept]))
      }
      aa <- edge_cor(f$rna, f$lipid, e, fa)
      bb <- wd_cor(x, "CF453", e, fb); cc <- wd_cor(x, "CF459", e, fc)
      for (k in seq_len(nrow(rows))) {
        rec <- rows[k, ]; key <- rec$specification; j <- specs[[k]]
        p <- fit$omission_draws[fit$omission_draws$specification == key & fit$omission_draws$draw == o$draw, , drop = FALSE]
        stopifnot(nrow(p) == 1L, p$n_FCG == length(fa), p$n_CF453 == length(fb), p$n_CF459 == length(fc))
        z <- triplet(aa[j], bb[j], cc[j])
        for (field in names(z)) compare(z[field], p[[field]], paste(key, "independent omission", o$draw, field))
      }
    }
    record(paste(context, "all 61 shared omission configurations replayed"))
  }
  result <- do.call(rbind, checks)
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(qa_dir, "independent_actual_replay_checks.csv"), row.names = FALSE)
  writeLines(c("Independent WP3 original-male-cohort actual-data QA passed.",
               paste("Checks:", nrow(result)),
               paste("Maximum absolute numerical difference:", format(max(result$max_abs_difference), digits = 10L)),
               "All exact, null, paired and 61 omission banks replayed using base cor/lm.",
               paste("Comparison:", normalizePath(comparison_file)),
               paste("Contexts:", normalizePath(contexts_file)),
               paste("Original analysis:", normalizePath(original_analysis_file))), file.path(qa_dir, "README.txt"))
  message("Independent WP3 original-cohort QA passed: ", nrow(result), " checks.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% c(3L, 4L)) {
    stop("Usage: Rscript tests/qa_wd_fcg_comparison_actual.R COMPARISON_RDS CONTEXTS_RDS ORIGINAL_ANALYSIS_RDS [QA_DIR]")
  }
  qa_wd_fcg_comparison_actual(args[1L], args[2L], args[3L],
      if (length(args) == 4L) args[4L] else "tests/outputs/wp3_fcg_actual_qa")
}
