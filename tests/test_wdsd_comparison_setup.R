# Run from repository root: Rscript tests/test_wdsd_comparison_setup.R
source("scripts/support_functions.R")
source("scripts/core/cyp2c38_reference.R")
source("scripts/core/wdsd_intake.R")
source("scripts/core/wdsd_reference.R")
source("scripts/core/wdsd_comparison_setup.R")
source("scripts/core/wdsd_network_statistics.R")

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected error matching: ", pattern)
  }
}
quiet_setup <- function(data, reference) {
  result <- NULL
  invisible(utils::capture.output(result <- wp2_prepare_contexts(data, reference)))
  result
}
data <- wp2_load_inputs()
reference <- wp2_load_reference()
prepared <- wp2_prepare_contexts(data, reference)
stopifnot(nrow(prepared$edge_eligibility_ledger) == 97812L,
          sum(prepared$edge_eligibility_ledger$eligible) == 58656L,
          identical(prepared$coverage$n_eligible_full_edges, c(18800L, 22560L, 9024L, 8272L)),
          identical(prepared$coverage$n_eligible_primary_edges, c(86L, 249L, 52L, 47L)),
          nrow(prepared$lipid_crosswalk) == 1254L,
          sum(prepared$lipid_crosswalk$match_status == "unique_exact") == 752L,
          sum(prepared$lipid_crosswalk$match_status == "missing_exact_molecular_annotation") == 502L,
          !anyDuplicated(prepared$edge_eligibility_ledger$edge_id))

# Independently join canonical identities and replay the original matrix values.
manual_matches <- match(wp1_canonical_lipid(prepared$lipid_crosswalk$reference_lipid),
                         data$lipid_ledger$canonical_lipid)
stopifnot(identical(prepared$lipid_crosswalk$feature_id, data$lipid_ledger$feature_id[manual_matches]))
for (name in names(prepared$contexts)) {
  x <- prepared$contexts[[name]]
  stopifnot(identical(colnames(x$rna_SD), x$metadata_SD$mouse_id),
            identical(colnames(x$lipid_SD), x$metadata_SD$mouse_id),
            identical(colnames(x$rna_WD), x$metadata_WD$mouse_id),
            identical(colnames(x$lipid_WD), x$metadata_WD$mouse_id),
            ncol(x$rna_SD) == 7L, ncol(x$rna_WD) == 7L,
            all(x$edges$eligible), all(x$metadata_SD$diet == "SD"), all(x$metadata_WD$diet == "WD"))
  direct <- stats::cor(t(x$reference_rna), t(x$reference_lipid))
  edge_lookup <- cbind(match(x$edges$gene, rownames(direct)), match(x$edges$feature_id, colnames(direct)))
  stopifnot(max(abs(direct[edge_lookup] - x$edges$r_reference)) < 1e-10)
  # Statistics must accept our exact field contract and frozen reference IDs.
  normalized <- wp2_normalize_context(x, name)
  stopifnot(all(normalized$edges$eligible),
            identical(normalized$edges$edge_id, x$edges$edge_id),
            all(is.finite(normalized$edges$r_SD)), all(is.finite(normalized$edges$r_WD)))
  selected <- unique(round(seq(1, nrow(x$edges), length.out = 13L)))
  for (i in selected) {
    gene <- x$edges$gene[i]; feature <- x$edges$feature_id[i]
    direct_SD <- stats::cor(data$rna[[x$zone]][gene, x$metadata_SD$mouse_id],
                            data$lipid[feature, x$metadata_SD$mouse_id])
    stopifnot(abs(direct_SD - normalized$edges$r_SD[i]) < 1e-10)
  }
}

# One-modality reorder and metadata relabeling must fail before computing edges.
bad <- data
bad$lipid <- bad$lipid[, rev(seq_len(ncol(bad$lipid))), drop = FALSE]
expect_error(quiet_setup(bad, reference), "ordered metadata")
bad <- data
bad$rna$Periportal <- bad$rna$Periportal[, rev(seq_len(28L)), drop = FALSE]
expect_error(quiet_setup(bad, reference), "same mice in identical order")
bad <- data
bad$metadata$mouse_id[1:2] <- rev(bad$metadata$mouse_id[1:2])
expect_error(quiet_setup(bad, reference), "ordered metadata")
bad_reference <- reference
bad_reference[[1L]]$lipids <- bad_reference[[1L]]$lipids[, rev(seq_len(ncol(bad_reference[[1L]]$lipids))), drop = FALSE]
expect_error(quiet_setup(data, bad_reference), "pairing is inconsistent")
bad <- data
bad$lipid_ledger$canonical_lipid[2L] <- bad$lipid_ledger$canonical_lipid[1L]
expect_error(quiet_setup(bad, reference), "unique feature/canonical")
bad_reference <- reference
bad_reference[[1L]]$edges$reference_lipid[2L] <- paste0(bad_reference[[1L]]$edges$reference_lipid[1L], " ")
expect_error(quiet_setup(data, bad_reference), "ambiguous after formatting")

# Conservation and diet-change stages must share a BOTH-diet eligible universe.
context <- "male_periportal"
gene <- prepared$contexts[[context]]$edges$gene[1L]
feature <- prepared$contexts[[context]]$edges$feature_id[1L]
SD <- which(data$metadata$sex == "male" & data$metadata$diet == "SD")
WD <- which(data$metadata$sex == "male" & data$metadata$diet == "WD")
bad <- data
bad$rna$Periportal[gene, SD] <- 3
out <- quiet_setup(bad, reference)
rows <- out$edge_eligibility_ledger$context == context & out$edge_eligibility_ledger$gene == gene
stopifnot(all(out$edge_eligibility_ledger$invariant_gene_SD[rows]),
          all(!out$edge_eligibility_ledger$invariant_gene_WD[rows]),
          all(!out$edge_eligibility_ledger$eligible[rows]),
          !gene %in% out$contexts[[context]]$edges$gene)
bad <- data
bad$lipid[feature, WD] <- 2
out <- quiet_setup(bad, reference)
rows <- out$edge_eligibility_ledger$feature_id == feature & out$edge_eligibility_ledger$sex == "male"
rows[is.na(rows)] <- FALSE
stopifnot(all(out$edge_eligibility_ledger$invariant_lipid_WD[rows]),
          all(!out$edge_eligibility_ledger$eligible[rows]))
bad <- data
bad$rna$Periportal[gene, WD[1L]] <- NA_real_
out <- quiet_setup(bad, reference)
rows <- out$edge_eligibility_ledger$context == context & out$edge_eligibility_ledger$gene == gene
stopifnot(all(out$edge_eligibility_ledger$nonfinite_gene_WD[rows]), all(!out$edge_eligibility_ledger$eligible[rows]))
bad <- data
bad$rna$Periportal <- bad$rna$Periportal[setdiff(rownames(bad$rna$Periportal), gene), , drop = FALSE]
out <- quiet_setup(bad, reference)
rows <- out$edge_eligibility_ledger$context == context & out$edge_eligibility_ledger$gene == gene
stopifnot(all(out$edge_eligibility_ledger$missing_gene[rows]), all(!out$edge_eligibility_ledger$eligible[rows]))
bad <- data
keep <- bad$lipid_ledger$feature_id == feature
bad$lipid_ledger$keep_primary[keep] <- FALSE
bad$lipid_ledger$keep_qc30[keep] <- FALSE
bad$lipid_ledger$exclusion_reason[keep] <- "test_analytical_exclusion"
bad$lipid <- bad$lipid[rownames(bad$lipid) != feature, , drop = FALSE]
out <- quiet_setup(bad, reference)
rows <- out$edge_eligibility_ledger$feature_id == feature
rows[is.na(rows)] <- FALSE
stopifnot(all(out$edge_eligibility_ledger$excluded_lipid[rows]),
          all(!out$edge_eligibility_ledger$missing_lipid[rows]), all(!out$edge_eligibility_ledger$eligible[rows]))

qc <- wp2_contexts_qc30(prepared$contexts, data$lipid_ledger)
stopifnot(identical(vapply(qc, function(x) sum(x$edges$primary_edge), integer(1L)),
                    c(male_periportal = 86L, male_pericentral = 245L,
                      female_periportal = 50L, female_pericentral = 45L)))
alternate <- wp2_contexts_data_layer(prepared$contexts, data)
for (name in names(alternate)) {
  original <- prepared$contexts[[name]]
  stopifnot(identical(alternate[[name]]$edges$edge_id, original$edges$edge_id),
            max(abs(wp2_cor_edges(alternate[[name]]$rna_SD, original$lipid_SD, original$edges) -
                    wp2_cor_edges(original$rna_SD, original$lipid_SD, original$edges))) < 1e-10)
}
cat("WDSD strict comparison setup tests passed.\n")
