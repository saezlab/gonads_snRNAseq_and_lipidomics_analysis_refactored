source("scripts/config/analysis_release.R")
# Independent numerical audit of the completed pooled Figure 3 follow-up.
# This auditor uses scalar Pearson/OLS and saved draw indices, not the pooled
# statistics helper or its p-value functions. It never changes analysis files.
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/cyp2c38_reference.R")
source("scripts/core/wdsd_intake.R")
source("scripts/core/wdsd_pooled_reference.R")

qa_pooled_equal <- function(actual, expected, label, tolerance = 1e-9) {
  if (length(actual) != length(expected)) stop(label, ": length mismatch.")
  defined <- is.finite(actual) & is.finite(expected)
  same_missing <- is.na(actual) & is.na(expected)
  same_infinite <- is.infinite(actual) & is.infinite(expected) & sign(actual) == sign(expected)
  if (any(!(defined | same_missing | same_infinite))) stop(label, ": defined/undefined mismatch.")
  if (any(defined) && any(abs(actual[defined] - expected[defined]) >
                           tolerance * pmax(1, abs(expected[defined])))) {
    stop(label, ": numerical discrepancy exceeds tolerance.")
  }
  if (any(defined)) max(abs(actual[defined] - expected[defined])) else 0
}

qa_pooled_direct_edges <- function(rna, lipids, edges) {
  if (!identical(colnames(rna), colnames(lipids))) stop("Independent correlation matrices have different mouse order.")
  # stats::cor uses its own covariance calculation, independently of the
  # production helper's row-unit normalization and tcrossprod implementation.
  matrix_r <- suppressWarnings(stats::cor(t(rna), t(lipids)))
  index <- cbind(match(edges$gene, rownames(matrix_r)), match(edges$feature_id, colnames(matrix_r)))
  if (anyNA(index)) stop("Independent coefficient replay has an unmeasured frozen feature.")
  unname(matrix_r[index])
}

qa_pooled_center <- function(x, sex) {
  answer <- x
  # Scalar model residuals independently reproduce within-sex mean removal.
  for (i in seq_len(nrow(x))) answer[i, ] <- stats::residuals(stats::lm(x[i, ] ~ factor(sex)))
  answer
}

qa_pooled_actual <- function(out = file.path(analysis_path("historical_wf12"),"followup"), strict_only = FALSE) {
  required <- file.path(out, c("pooled_analysis.rds", "paired_network_inputs.rds", "audit.rds"))
  if (any(!file.exists(required))) stop("Actual outputs are not ready: ", paste(required[!file.exists(required)], collapse = ", "))
  qa_dir <- file.path(out, "QA")
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- wp2_load_inputs()
  original <- wp2_pooled_reference()
  audit <- readRDS(file.path(out, "audit.rds"))
  analyses <- list(strict = readRDS(file.path(out, "pooled_analysis.rds")))
  contexts <- list(strict = readRDS(file.path(out, "paired_network_inputs.rds")))
  paths <- c(strict = "pooled_analysis.rds")
  if (!isTRUE(strict_only)) {
    paths <- c(paths, conditional_compatibility = "conditional_compatibility_analysis.rds",
               QC_CV_le30 = "QC_CV_le30_analysis.rds", supplied_fill_ge50 = "supplied_fill_ge50_analysis.rds")
    if (any(!file.exists(file.path(out, c(unname(paths), "pooled_contexts.rds"))))) {
      stop("Sensitivity outputs are not ready; --strict-only can audit the main result immediately.")
    }
    analyses <- lapply(file.path(out, paths), readRDS)
    names(analyses) <- names(paths)
    contexts <- readRDS(file.path(out, "pooled_contexts.rds"))
    stopifnot(setequal(names(analyses), names(contexts)))
  }
  numerical <- null_tails <- reconstruction <- coverage <- list()
  evidence <- function(analysis, check, value, tolerance, notes = "") {
    data.frame(analysis = analysis, check = check, max_abs_error = value,
               tolerance = tolerance, passed = TRUE, notes = notes, stringsAsFactors = FALSE)
  }
  sd_metadata <- inputs$metadata[inputs$metadata$diet == "SD", , drop = FALSE]
  stopifnot(nrow(sd_metadata) == 14L, all(table(sd_metadata$sex) == 7L),
            nrow(original$metadata) == 43L,
            identical(as.integer(table(factor(original$metadata$sex, levels = c("female", "male")))), c(21L, 22L)),
            nrow(original$figure3_edges) == 192L, nrow(original$edges) == 8896L,
            !anyDuplicated(audit$safe_map$reference_lipid), !anyDuplicated(audit$safe_map$feature_id),
            !anyDuplicated(audit$conditional_compatibility_map$reference_lipid),
            !anyDuplicated(audit$conditional_compatibility_map$feature_id))
  original_key <- paste(original$edges$gene, original$edges$reference_lipid, sep = "\r")
  for (name in names(analyses)) {
    fit <- analyses[[name]]; context <- contexts[[name]]; e <- fit$edges
    stopifnot(identical(as.character(context$metadata_SD$mouse_id), colnames(context$rna_SD)),
              identical(colnames(context$rna_SD), colnames(context$lipid_SD)),
              identical(as.character(context$reference_metadata$mouse_id), colnames(context$reference_rna)),
              identical(colnames(context$reference_rna), colnames(context$reference_lipid)),
              identical(context$metadata_SD$mouse_id, sd_metadata$mouse_id),
              identical(context$metadata_SD$sex, sd_metadata$sex),
              identical(context$reference_metadata$mouse_id, original$metadata$mouse_id),
              identical(context$reference_metadata$group, original$metadata$group),
              !anyDuplicated(e$edge_id), !anyDuplicated(paste(e$gene, e$feature_id, sep = "\r")),
              all(e$eligible), nrow(context$edges) == nrow(e))
    mapping <- unique(e[, c("reference_lipid", "feature_id")])
    stopifnot(!anyDuplicated(mapping$reference_lipid), !anyDuplicated(mapping$feature_id))
    map_source <- if (name == "conditional_compatibility") audit$conditional_compatibility_map else audit$safe_map
    map_index <- match(mapping$reference_lipid, map_source$reference_lipid)
    stopifnot(!anyNA(map_index), identical(mapping$feature_id, map_source$feature_id[map_index]))
    if (name != "conditional_compatibility") {
      ledger_index <- match(mapping$feature_id, inputs$lipid_ledger$feature_id)
      stopifnot(identical(wp1_canonical_lipid(mapping$reference_lipid), inputs$lipid_ledger$canonical_lipid[ledger_index]))
    }
    if (name == "QC_CV_le30") stopifnot(all(e$feature_id %in% inputs$lipid_ledger$feature_id[inputs$lipid_ledger$keep_qc30]))
    if (name == "supplied_fill_ge50") {
      fi <- match(e$feature_id, inputs$lipid_ledger$feature_id)
      stopifnot(all(is.finite(inputs$lipid_ledger[["Fill %"]][fi])), all(inputs$lipid_ledger[["Fill %"]][fi] >= 0.5))
    }
    archived_index <- match(paste(e$gene, e$reference_lipid, sep = "\r"), original_key)
    stopifnot(!anyNA(archived_index), identical(e$primary_edge, original$edges$primary_edge[archived_index]))
    err <- qa_pooled_equal(e$r_reference, original$edges$r_reference[archived_index], paste(name, "archived reference"), 1e-12)
    numerical[[paste(name, "archive")]] <- evidence(name, "frozen archival identity/coefficient replay", err, 1e-12)
    # Reload cohort matrices independently; context inputs must equal the
    # supplied RNA layer and log2 lipid matrix in the recorded mouse order.
    expected_rna <- inputs$rna_hepatocytes[rownames(context$rna_SD), context$metadata_SD$mouse_id, drop = FALSE]
    expected_lipid <- inputs$lipid[rownames(context$lipid_SD), context$metadata_SD$mouse_id, drop = FALSE]
    err <- qa_pooled_equal(as.vector(context$rna_SD), as.vector(expected_rna), paste(name, "supplied RNA"), 1e-12)
    numerical[[paste(name, "rna")]] <- evidence(name, "independently reloaded WDSD RNA", err, 1e-12)
    err <- qa_pooled_equal(as.vector(context$lipid_SD), as.vector(expected_lipid), paste(name, "supplied lipid"), 1e-12)
    numerical[[paste(name, "lipid")]] <- evidence(name, "independently reloaded log2 WDSD lipids", err, 1e-12)
    expected_reference_rna <- original$rna[rownames(context$reference_rna), context$reference_metadata$mouse_id, drop = FALSE]
    reverse_map <- map_source$reference_lipid[match(rownames(context$reference_lipid), map_source$feature_id)]
    stopifnot(!anyNA(reverse_map))
    expected_reference_lipid <- original$lipids[reverse_map, context$reference_metadata$mouse_id, drop = FALSE]
    err <- qa_pooled_equal(as.vector(context$reference_rna), as.vector(expected_reference_rna), paste(name, "FCG RNA"), 1e-12)
    numerical[[paste(name, "reference_rna")]] <- evidence(name, "independently reloaded FCG RNA", err, 1e-12)
    err <- qa_pooled_equal(as.vector(context$reference_lipid), as.vector(expected_reference_lipid), paste(name, "FCG lipids"), 1e-12)
    numerical[[paste(name, "reference_lipid")]] <- evidence(name, "independently reloaded FCG log2 lipids", err, 1e-12)
    direct_SD <- qa_pooled_direct_edges(context$rna_SD, context$lipid_SD, e)
    direct_REF <- qa_pooled_direct_edges(context$reference_rna, context$reference_lipid, e)
    SD_centered <- qa_pooled_center(context$rna_SD, context$metadata_SD$sex)
    lipid_centered <- qa_pooled_center(context$lipid_SD, context$metadata_SD$sex)
    REF_centered <- qa_pooled_center(context$reference_rna, context$reference_metadata$sex)
    REF_lipid_centered <- qa_pooled_center(context$reference_lipid, context$reference_metadata$sex)
    direct_SD_centered <- qa_pooled_direct_edges(SD_centered, lipid_centered, e)
    direct_REF_centered <- qa_pooled_direct_edges(REF_centered, REF_lipid_centered, e)
    for (field in c("r_SD", "r_reference", "r_SD_withinsex", "r_reference_withinsex")) {
      expected <- switch(field, r_SD = direct_SD, r_reference = direct_REF,
                         r_SD_withinsex = direct_SD_centered, r_reference_withinsex = direct_REF_centered)
      err <- qa_pooled_equal(e[[field]], expected, paste(name, field), 1e-10)
      numerical[[paste(name, field)]] <- evidence(name, paste("independent Pearson", field), err, 1e-10)
    }
    coverage[[name]] <- data.frame(analysis = name, full_edges = nrow(e),
      primary_edges = sum(e$primary_edge), full_genes = length(unique(e$gene)),
      primary_genes = length(unique(e$gene[e$primary_edge])), matched_lipids = length(unique(e$feature_id)),
      n_SD = ncol(context$rna_SD), n_reference = ncol(context$reference_rna),
      mapping_one_to_one = TRUE, stringsAsFactors = FALSE)
    # Check the complete unthresholded sparse precursor, rather than accepting
    # a partial or post-selected complement. All 139 genes cross each matched
    # species; every matched displayed species retains its original three genes.
    stopifnot(length(unique(e$gene)) == 139L,
              nrow(e) == 139L * length(unique(e$feature_id)),
              sum(e$primary_edge) == 3L * length(unique(e$feature_id)),
              all(as.integer(table(e$feature_id[e$primary_edge])) == 3L))
    if (name == "strict") stopifnot(nrow(e) == 5004L, sum(e$primary_edge) == 108L,
                                     length(unique(e$feature_id)) == 36L, length(unique(e$gene[e$primary_edge])) == 15L)
    # Scalar lm and cor.test independently validate every primary/complement,
    # predefined gene and centered inference row.
    for (i in seq_len(nrow(fit$inference))) {
      row <- fit$inference[i, ]
      index <- which(e$eligible & (row$universe == "full" | e$primary_edge))
      if (row$scope == "gene") index <- index[e$gene[index] == row$gene]
      stopifnot(length(index) == row$n_edges, row$n_SD == 14L, row$n_reference == 43L)
      if (row$analysis == "pooled") {
        x <- direct_REF[index]; y <- direct_SD[index]
      } else {
        x <- direct_REF_centered[index]; y <- direct_SD_centered[index]
        stopifnot(all(is.na(unlist(row[grep("_(absolute|upper|lower|equal_tail)_p$", names(row))]))))
      }
      if (any(!is.finite(x)) || any(!is.finite(y)) || stats::sd(x) == 0 || stats::sd(y) == 0) {
        expected <- c(coefficient_vector_r = NA_real_, ULM_score = NA_real_, ULM_p_value = NA_real_, slope = NA_real_)
      } else {
        test <- stats::cor.test(x, y, method = "pearson")
        ols <- summary(stats::lm(y ~ x))$coefficients[2L, ]
        expected <- c(coefficient_vector_r = unname(test$estimate), ULM_score = unname(ols["t value"]),
                      ULM_p_value = unname(ols["Pr(>|t|)"]), slope = unname(ols["Estimate"]))
        qa_pooled_equal(expected["ULM_score"], unname(test$statistic), "OLS/cor.test equivalence", 1e-9)
        qa_pooled_equal(expected["ULM_p_value"], test$p.value, "OLS/cor.test nominal p equivalence", 1e-9)
      }
      errors <- vapply(names(expected), function(field) qa_pooled_equal(row[[field]], expected[field],
                        paste(name, row$analysis, row$specification, field), 1e-9), numeric(1L))
      numerical[[paste(name, row$analysis, row$specification)]] <- evidence(name,
        paste("scalar OLS/cor.test", row$analysis, row$specification), max(errors), 1e-9,
        if (is.finite(row$ULM_score) && row$ULM_p_value == 0) "Nominal p underflows machine arithmetic; interpret as a threshold, not exact zero." else "")
      if (row$analysis != "pooled") next
      for (method in c("withinsex_permutation", "unrestricted_permutation", "withinsex_bootstrap", "unrestricted_bootstrap")) {
        bank <- fit[[switch(method, withinsex_permutation = "permutation_withinsex_draws",
                             unrestricted_permutation = "permutation_unrestricted_draws",
                             withinsex_bootstrap = "bootstrap_withinsex_draws",
                             unrestricted_bootstrap = "bootstrap_unrestricted_draws")]]
        if (grepl("permutation", method)) {
          stopifnot(bank$draw[1L] == 0L, nrow(bank) == fit$metadata$random_permutations + 1L,
                    all(bank$draw[-1L] == seq_len(fit$metadata$random_permutations)))
          qa_pooled_equal(bank[[row$specification]][1L], row$coefficient_vector_r, paste(name, "identity"), 1e-12)
          null <- bank[[row$specification]][bank$draw > 0L]
        } else null <- bank[[row$specification]]
        valid <- is.finite(null); n_valid <- sum(valid); null <- null[valid]
        observed <- row$coefficient_vector_r
        absolute <- (1 + sum(abs(null) >= abs(observed) - 1e-12)) / (1 + n_valid)
        upper <- (1 + sum(null >= observed - 1e-12)) / (1 + n_valid)
        lower <- (1 + sum(null <= observed + 1e-12)) / (1 + n_valid)
        equal <- min(1, 2 * min(upper, lower))
        expected_tail <- c(absolute_p = absolute, upper_p = upper, lower_p = lower, equal_tail_p = equal,
                           valid_draws = n_valid, invalid_draws = sum(!valid), null_median = stats::median(null))
        error <- vapply(names(expected_tail), function(field) qa_pooled_equal(row[[paste(method, field, sep = "_")]],
                          expected_tail[field], paste(name, row$specification, method, field), 1e-12), numeric(1L))
        null_tails[[paste(name, row$specification, method)]] <- data.frame(
          analysis = name, specification = row$specification, method = method,
          n_valid_random_draws = n_valid, n_invalid_random_draws = sum(!valid),
          identity_audit_excluded = grepl("permutation", method), max_abs_error = max(error),
          absolute_p = absolute, upper_p = upper, lower_p = lower, equal_tail_p = equal,
          passed = TRUE, stringsAsFactors = FALSE)
      }
    }
    # Saved banks are shared whole-mouse positions. Their long-form ledger
    # must reproduce every position, mouse identity and sex label exactly.
    for (method in names(fit$index_banks)) {
      bank <- fit$index_banks[[method]]
      ledger <- fit$resampling_indices[fit$resampling_indices$method == method, , drop = FALSE]
      meta <- if (method == "paired_reference") context$reference_metadata else context$metadata_SD
      positions <- rep(seq_len(ncol(bank)), times = nrow(bank))
      source <- as.vector(t(bank))
      stopifnot(nrow(ledger) == length(source), identical(as.integer(ledger$position), positions),
                identical(as.integer(ledger$source_position), as.integer(source)),
                identical(as.character(ledger$position_mouse), as.character(meta$mouse_id[positions])),
                identical(as.character(ledger$source_mouse), as.character(meta$mouse_id[source])),
                identical(as.character(ledger$position_sex), as.character(meta$sex[positions])),
                identical(as.character(ledger$source_sex), as.character(meta$sex[source])))
      if (grepl("withinsex", method) || method %in% c("paired_SD", "paired_reference")) {
        stopifnot(all(meta$sex[positions] == meta$sex[source]))
      }
      if (method == "paired_reference") stopifnot(all(meta$group[positions] == meta$group[source]))
      if (grepl("permutation", method)) {
        stopifnot(all(apply(bank, 1L, function(z) identical(sort(as.integer(z)), seq_len(ncol(bank))))))
      }
    }
    stopifnot(!identical(fit$index_banks$bootstrap_withinsex_rna, fit$index_banks$bootstrap_withinsex_lipid),
              !identical(fit$index_banks$bootstrap_unrestricted_rna, fit$index_banks$bootstrap_unrestricted_lipid),
              fit$metadata$identity_audit_draw == 0L, !fit$metadata$identity_in_null_draw_count)
    # Reconstruct first, middle and last null/stability draws directly from the
    # saved mouse banks. Recompute all inference scopes from those same whole
    # vectors to detect feature-wise or edge-wise resampling mistakes.
    for (method in c("withinsex_permutation", "unrestricted_permutation", "withinsex_bootstrap", "unrestricted_bootstrap", "paired", "paired_reference")) {
      output_bank <- fit[[switch(method, withinsex_permutation = "permutation_withinsex_draws",
        unrestricted_permutation = "permutation_unrestricted_draws", withinsex_bootstrap = "bootstrap_withinsex_draws",
        unrestricted_bootstrap = "bootstrap_unrestricted_draws", paired = "paired_stability_draws",
        paired_reference = "paired_reference_stability_draws")]]
      count <- if (grepl("permutation", method)) fit$metadata$random_permutations else nrow(output_bank)
      for (draw in unique(c(1L, as.integer(ceiling(count / 2)), count))) {
        if (grepl("permutation", method)) {
          ri <- seq_len(14L)
          li <- fit$index_banks[[paste0("permutation_", sub("_permutation$", "", method))]][draw, ]
        } else if (grepl("bootstrap", method)) {
          kind <- sub("_bootstrap$", "", method)
          ri <- fit$index_banks[[paste0("bootstrap_", kind, "_rna")]][draw, ]
          li <- fit$index_banks[[paste0("bootstrap_", kind, "_lipid")]][draw, ]
        } else ri <- li <- fit$index_banks$paired_SD[draw, ]
        rna_draw <- context$rna_SD[, ri, drop = FALSE]
        lipid_draw <- context$lipid_SD[, li, drop = FALSE]
        # Null draws deliberately pair different source mice. Their resulting
        # columns denote common resampling slots, while source mouse identities
        # remain independently verified in the saved index/ID ledger above.
        colnames(rna_draw) <- colnames(lipid_draw) <- paste0("draw_slot_", seq_len(ncol(rna_draw)))
        values <- qa_pooled_direct_edges(rna_draw, lipid_draw, e)
        reference_values <- if (method == "paired_reference") {
          ix <- fit$index_banks$paired_reference[draw, ]
          qa_pooled_direct_edges(context$reference_rna[, ix, drop = FALSE], context$reference_lipid[, ix, drop = FALSE], e)
        } else e$r_reference
        for (i in which(fit$inference$analysis == "pooled")) {
          row <- fit$inference[i, ]
          index <- which(e$eligible & (row$universe == "full" | e$primary_edge))
          if (row$scope == "gene") index <- index[e$gene[index] == row$gene]
          expected <- if (any(!is.finite(reference_values[index])) || any(!is.finite(values[index]))) NA_real_ else {
            suppressWarnings(stats::cor(reference_values[index], values[index]))
          }
          output_row <- match(draw, output_bank$draw)
          error <- qa_pooled_equal(output_bank[[row$specification]][output_row], expected,
                                   paste(name, method, draw, row$specification), 1e-9)
          reconstruction[[paste(name, method, draw, row$specification)]] <- data.frame(
            analysis = name, method = method, draw = draw, specification = row$specification,
            expected_coefficient_vector_r = expected,
            saved_coefficient_vector_r = output_bank[[row$specification]][output_row],
            max_abs_error = error, passed = TRUE, stringsAsFactors = FALSE)
        }
      }
    }
    message("Independent actual-output numerical audit passed: ", name,
            " (", nrow(e), " full / ", sum(e$primary_edge), " displayed pairs).")
  }
  numerical <- do.call(rbind, numerical); null_tails <- do.call(rbind, null_tails)
  reconstruction <- do.call(rbind, reconstruction); coverage <- do.call(rbind, coverage)
  suffix <- if (strict_only) "_strict" else ""
  utils::write.csv(numerical, file.path(qa_dir, paste0("pooled_independent_numerical_audit", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(null_tails, file.path(qa_dir, paste0("pooled_independent_null_tail_audit", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(reconstruction, file.path(qa_dir, paste0("pooled_independent_draw_reconstruction", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(coverage, file.path(qa_dir, paste0("pooled_independent_context_coverage", suffix, ".csv")), row.names = FALSE)
  note <- c("Independent actual-output audit of pooled Figure 3 follow-up: PASS.",
    paste("Audited mappings/sensitivities:", paste(names(analyses), collapse = ", ")),
    paste("UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Scalar native ULM/Pearson checks:", nrow(numerical), "; maximum absolute numerical discrepancy:", format(max(numerical$max_abs_error), digits = 6L)),
    paste("Independently reconstructed complete null-tail summaries:", nrow(null_tails)),
    paste("Independent whole-mouse draw/scoped-statistic reconstructions:", nrow(reconstruction)),
    "Verified all saved whole-mouse index/ID/sex ledgers, within-sex and four-group strata, identity exclusion and +1 tails.",
    "Verified one-to-one strict/conditional recorded mapping, all unchanged archival primary memberships, and complete 139-gene sparse precursor complements.",
    "Native ULM p-values remain nominal under dependent edges; unrestricted and sex-conditional animal nulls answer distinct questions.",
    "Zero numeric nominal p means machine underflow, not a literal probability of zero.")
  writeLines(note, file.path(qa_dir, paste0("pooled_independent_numerical_audit", suffix, ".txt")))
  invisible(list(numerical = numerical, null_tails = null_tails, reconstruction = reconstruction, coverage = coverage))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--validate-only" %in% args) {
    cat("Independent pooled actual-output QA definitions loaded successfully.\n")
  } else {
    output <- args[startsWith(args, "--output=")]
    out <- if (length(output)) substring(output[1L], 10L) else file.path(analysis_path("historical_wf12"),"followup")
    qa_pooled_actual(out, strict_only = "--strict-only" %in% args)
  }
}
