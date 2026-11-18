# Independent replay of completed WP3 results against stored paired matrices.
# Run: Rscript tests/qa_wd_crosscohort_actual.R ANALYSIS_RDS CONTEXTS_RDS [QA_DIR]
# Sourcing defines the verifier only; it never starts or modifies an analysis.

qa_wd_crosscohort_actual <- function(analysis_file, contexts_file,
                                   qa_dir = "tests/outputs/wp3_actual_qa") {
  stopifnot(file.exists(analysis_file), file.exists(contexts_file))
  fit <- readRDS(analysis_file)
  contexts <- readRDS(contexts_file)
  if (!is.null(contexts$contexts)) contexts <- contexts$contexts
  stopifnot(all(c("male_periportal", "male_pericentral") %in% names(contexts)),
            all(c("edges", "inference", "pairing_permutation_indices",
                  "resampling_index_matrices", "loco_omission_ledger") %in% names(fit)))
  checks <- list()
  record <- function(check, difference = 0, detail = "") {
    if (!is.finite(difference) || difference > 1e-9) {
      stop("Independent WP3 QA failed: ", check, "; difference=", difference)
    }
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, max_abs_difference = difference, detail = detail,
      stringsAsFactors = FALSE)
  }
  compare <- function(a, b, check) {
    a <- as.numeric(a); b <- as.numeric(b)
    stopifnot(length(a) == length(b), identical(is.na(a), is.na(b)))
    valid <- is.finite(a) & is.finite(b)
    stopifnot(all(valid | (is.na(a) & is.na(b)) | (is.infinite(a) & a == b)))
    record(check, if (any(valid)) max(abs(a[valid] - b[valid])) else 0)
  }
  # This verifier deliberately uses base cor/lm rather than WP2/WP3 helpers.
  edges_direct <- function(x, cohort, e, rna_ix = NULL, lipid_ix = rna_ix) {
    g <- x[[paste0("rna_", cohort)]]
    l <- x[[paste0("lipid_", cohort)]]
    if (is.null(rna_ix)) rna_ix <- seq_len(ncol(g))
    if (is.null(lipid_ix)) lipid_ix <- seq_len(ncol(l))
    cmat <- suppressWarnings(stats::cor(t(g[, rna_ix, drop = FALSE]),
                                       t(l[, lipid_ix, drop = FALSE])))
    as.numeric(cmat[cbind(match(e$gene, rownames(cmat)),
                         match(e$feature_id, colnames(cmat)))])
  }
  vector_cor <- function(a, b) {
    if (any(!is.finite(a)) || any(!is.finite(b)) || length(a) < 3L) return(NA_real_)
    suppressWarnings(stats::cor(a, b))
  }
  scope_indices <- function(e, rec) {
    keep <- e$eligible
    if (rec$universe == "displayed") keep <- keep & e$primary_edge
    if (rec$scope == "gene") keep <- keep & e$gene == rec$gene
    which(keep)
  }
  exact_ledger <- fit$pairing_permutation_indices
  permutation_banks <- lapply(c("CF459", "CF453"), function(cohort) {
    n <- if (cohort == "CF459") 4L else 7L
    ledger <- exact_ledger[exact_ledger$cohort == cohort, , drop = FALSE]
    stopifnot(nrow(ledger) == factorial(n) * n)
    bank <- matrix(ledger$source_position, byrow = TRUE, ncol = n)
    stopifnot(identical(as.integer(bank[1L, ]), seq_len(n)),
              !anyDuplicated(apply(bank, 1L, paste, collapse = ",")),
              all(apply(bank, 1L, function(z) identical(sort(as.integer(z)), seq_len(n)))))
    ids <- colnames(contexts[[1L]][[paste0("rna_", cohort)]])
    stopifnot(identical(as.character(ledger$source_mouse), ids[ledger$source_position]))
    record(paste(cohort, "complete unique permutation ledger"), detail = paste(factorial(n), "pairings"))
    bank
  })
  names(permutation_banks) <- c("CF459", "CF453")

  for (context in names(contexts)) {
    x <- contexts[[context]]
    e <- fit$edges[fit$edges$context == context, , drop = FALSE]
    stopifnot(identical(e$edge_id, x$edges$edge_id), all(e$eligible))
    for (cohort in c("CF453", "CF459")) {
      ids <- colnames(x[[paste0("rna_", cohort)]])
      stopifnot(identical(ids, colnames(x[[paste0("lipid_", cohort)]])),
                identical(ids, colnames(contexts[[1L]][[paste0("rna_", cohort)]])))
      direct <- edges_direct(x, cohort, e)
      compare(direct, e[[paste0("r_", cohort, "_WD")]], paste(context, cohort, "all observed coefficients"))
      rev <- rev(seq_along(ids))
      compare(direct, edges_direct(x, cohort, e, rev, rev),
              paste(context, cohort, "joint mouse reorder invariance"))
    }
    a <- edges_direct(x, "CF453", e); b <- edges_direct(x, "CF459", e)
    rows <- fit$inference[fit$inference$context == context, , drop = FALSE]
    specs <- lapply(seq_len(nrow(rows)), function(k) scope_indices(e, rows[k, , drop = FALSE]))
    names(specs) <- rows$specification
    for (k in seq_len(nrow(rows))) {
      rec <- rows[k, , drop = FALSE]; ix <- specs[[k]]
      model <- summary(stats::lm(b[ix] ~ a[ix]))$coefficients[2L, ]
      compare(as.numeric(rec$coefficient_vector_r), vector_cor(a[ix], b[ix]), paste(rec$specification, "observed r"))
      compare(as.numeric(rec$ULM_score), as.numeric(model["t value"]), paste(rec$specification, "independent OLS t"))
      compare(as.numeric(rec$ULM_p_value), as.numeric(model["Pr(>|t|)"]), paste(rec$specification, "independent OLS nominal p"))
      stopifnot(rec$n_edges == length(ix), rec$n_CF453 == 7L, rec$n_CF459 == 4L)
    }
    for (cohort in c("CF459", "CF453")) {
      indices <- permutation_banks[[cohort]]
      direct <- matrix(NA_real_, nrow(indices), length(specs), dimnames = list(NULL, names(specs)))
      for (j in seq_len(nrow(indices))) {
        v <- edges_direct(x, cohort, e, seq_len(ncol(indices)), indices[j, ])
        for (k in seq_along(specs)) {
          ix <- specs[[k]]
          direct[j, k] <- if (cohort == "CF459") vector_cor(a[ix], v[ix]) else vector_cor(v[ix], b[ix])
        }
      }
      field <- paste0("permutation_", cohort, "_draws")
      pfield <- paste0("pairing_", cohort, "_given_", if (cohort == "CF459") "CF453" else "CF459", "_p")
      for (k in seq_len(nrow(rows))) {
        rec <- rows[k, , drop = FALSE]; key <- rec$specification
        compare(direct[, k], as.numeric(fit[[field]][[key]]), paste(key, cohort, "complete null-bank direct Pearson replay"))
        valid <- is.finite(direct[, k])
        p <- mean(abs(direct[valid, k]) >= abs(rec$coefficient_vector_r) - 1e-12)
        compare(as.numeric(rec[[pfield]]), p, paste(key, cohort, "independent exhaustive tail"))
        stopifnot(p >= 1 / nrow(indices))
      }
      message("Independent WP3 QA replayed ", context, " ", cohort, " complete exact bank.")
    }
    indices <- fit$resampling_index_matrices
    B <- nrow(indices$CF459_null_rna)
    paired_B <- nrow(indices$CF453_paired)
    for (k in seq_len(nrow(rows))) {
      rec <- rows[k, , drop = FALSE]; key <- rec$specification
      vals <- fit$bootstrap_null_draws[[key]]; valid <- is.finite(vals)
      exceed <- sum(abs(vals[valid]) >= abs(rec$coefficient_vector_r) - 1e-12)
      compare(as.numeric(rec$null_bootstrap_p), (1 + exceed) / (1 + sum(valid)), paste(key, "null bootstrap tail"))
      compare(as.numeric(rec$bootstrap_completion_lower), (1 + exceed) / (1 + B), paste(key, "undefined mass lower bound"))
      compare(as.numeric(rec$bootstrap_completion_upper), (1 + exceed + sum(!valid)) / (1 + B), paste(key, "undefined mass upper bound"))
      stopifnot(rec$bootstrap_valid_draws == sum(valid), rec$bootstrap_invalid_draws == sum(!valid))
    }
    chosen <- unique(c(1L, min(17L, B), B,
                       which(!is.finite(fit$bootstrap_null_draws[[names(specs)[1L]]]))[1L]))
    chosen <- chosen[!is.na(chosen)]
    for (j in chosen) {
      v <- edges_direct(x, "CF459", e, indices$CF459_null_rna[j, ], indices$CF459_null_lipid[j, ])
      for (k in seq_along(specs)) {
        ix <- specs[[k]]; key <- names(specs)[k]
        compare(vector_cor(a[ix], v[ix]), fit$bootstrap_null_draws[[key]][j], paste(key, "direct null draw", j))
      }
    }
    for (j in unique(c(1L, min(17L, paired_B), paired_B))) {
      aa <- edges_direct(x, "CF453", e, indices$CF453_paired[j, ])
      bb <- edges_direct(x, "CF459", e, indices$CF459_paired[j, ])
      for (k in seq_along(specs)) {
        ix <- specs[[k]]; key <- names(specs)[k]
        compare(vector_cor(aa[ix], bb[ix]), fit$paired_stability_draws[[key]][j], paste(key, "direct paired draw", j))
      }
    }
    lc <- fit$loco_omission_ledger
    stopifnot(nrow(lc) == 39L, sum(lc$method == "joint_omission") == 28L,
              sum(lc$method == "CF453_only") == 7L, sum(lc$method == "CF459_only") == 4L)
    for (j in seq_len(nrow(lc))) {
      ia <- if (is.na(lc$omit_CF453[j])) 1:7 else setdiff(1:7, lc$omit_CF453[j])
      ib <- if (is.na(lc$omit_CF459[j])) 1:4 else setdiff(1:4, lc$omit_CF459[j])
      aa <- edges_direct(x, "CF453", e, ia); bb <- edges_direct(x, "CF459", e, ib)
      compare(aa, fit$loco_matrices[[context]]$r_CF453_WD[, j], paste(context, "all-edge CF453 omission", j))
      compare(bb, fit$loco_matrices[[context]]$r_CF459_WD[, j], paste(context, "all-edge CF459 omission", j))
      for (k in seq_along(specs)) {
        ix <- specs[[k]]; key <- names(specs)[k]
        rr <- fit$loco_draws$coefficient_vector_r[fit$loco_draws$specification == key & fit$loco_draws$draw == j]
        compare(vector_cor(aa[ix], bb[ix]), rr, paste(key, "direct omission", j))
      }
    }
  }
  result <- do.call(rbind, checks)
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(qa_dir, "independent_actual_replay_checks.csv"), row.names = FALSE)
  writeLines(c("Independent WP3 actual-data QA passed.",
               paste("Checks:", nrow(result)),
               paste("Maximum absolute numerical difference:", format(max(result$max_abs_difference), digits = 10L)),
               paste("Analysis:", normalizePath(analysis_file)),
               paste("Contexts:", normalizePath(contexts_file))), file.path(qa_dir, "README.txt"))
  message("Independent WP3 actual-data QA passed: ", nrow(result), " checks.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% c(2L, 3L)) {
    stop("Usage: Rscript tests/qa_wd_crosscohort_actual.R ANALYSIS_RDS CONTEXTS_RDS [QA_DIR]")
  }
  qa_wd_crosscohort_actual(args[1L], args[2L],
                         if (length(args) == 3L) args[3L] else "tests/outputs/wp3_actual_qa")
}
