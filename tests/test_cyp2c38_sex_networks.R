source("scripts/support_functions.R")
source("scripts/core/cyp2c38_reference.R")
source("scripts/core/cyp2c38_sex_networks.R")
fixed <- wp1_load_reference()
full <- wp1_full_reference()
stopifnot(nrow(full$lipids) == 1254L, ncol(full$lipids) == 43L,
          sum(full$metadata$gonadal_type == "Ovaries") == 21L,
          sum(full$metadata$gonadal_type == "Testes") == 22L,
          identical(names(full$gene), colnames(full$lipids)),
          identical(full$metadata$sample_id, colnames(full$lipids)),
          sum(full$fixed_reference) == 828L,
          identical(names(full$fixed_reference), rownames(full$lipids)))
idx <- match(fixed$coefficients$reference_lipid, full$coefficients$reference_lipid)
stopifnot(!anyNA(idx))
for (column in c("r_pooled", "r_female", "r_male")) {
  stopifnot(max(abs(full$coefficients[[column]][idx] - fixed$coefficients[[column]])) < 1e-12)
}
stopifnot(full$max_archive_difference < 1e-12)
for (group in c("Ovaries", "Testes")) {
  keep <- full$metadata$gonadal_type == group
  column <- if (group == "Ovaries") "r_female" else "r_male"
  for (j in c(1L, 252L, 1254L)) {
    direct <- stats::cor(full$gene[keep], full$lipids[j, keep])
    stopifnot(abs(full$coefficients[[column]][j] - direct) < 1e-12)
    order <- rev(which(keep))
    reordered <- stats::cor(full$gene[order], full$lipids[j, order])
    stopifnot(abs(direct - reordered) < 1e-12)
  }
}
out <- "/private/tmp/wp1_sex_network_audit"
manifest <- wp1_export_sex_networks(out, full)
stopifnot(nrow(manifest) == 6L,
          all(manifest$n_edges[manifest$universe == "full"] == 1254L),
          all(manifest$n_edges[manifest$universe == "fixed828"] == 828L))
for (path in manifest$file) {
  network <- utils::read.csv(path, check.names = FALSE)
  stopifnot(identical(names(network)[seq_len(3L)], c("node1", "node2", "mor")),
            !anyDuplicated(network$node2), all(network$node1 == "Cyp2c38"),
            all(is.finite(network$mor)))
}
malformed <- full
malformed$coefficients$fixed_reference[1L] <- NA
stopifnot(inherits(try(wp1_export_sex_networks(out, malformed), silent = TRUE), "try-error"))
ledger <- utils::read.csv("results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/feature_annotation_and_QC_ledger.csv",
                         check.names = FALSE)
endogenous <- ledger[ledger$keep_primary, , drop = FALSE]
stopifnot(!anyDuplicated(endogenous$canonical_lipid))
matched <- match(full$coefficients$reference_lipid, endogenous$canonical_lipid)
full_match <- full$coefficients
full_match$feature_id <- endogenous$feature_id[matched]
full_match$match_status <- ifelse(is.na(matched), "unmatched", "matched")
utils::write.csv(full_match, file.path(out, "full_network_strict_KO_matching.csv"), row.names = FALSE)
utils::write.csv(full$coefficients, file.path(out, "full_network_coefficients.csv"), row.names = FALSE)
utils::write.csv(manifest, file.path(out, "network_manifest.csv"), row.names = FALSE)
cat("PASS: full feature universe, mouse pairing, direct sex-specific Pearson replay,\n")
cat("archive and historical-subset replay, order invariance and six exports.\n")
cat("Full strict KO matches:", sum(!is.na(matched)), "out of", nrow(full_match), "\n")
cat("Fixed historical strict KO matches:", sum(!is.na(matched) & full$fixed_reference), "\n")
cat("New strict KO matches:", sum(!is.na(matched) & !full$fixed_reference), "\n")
print(full$fixed_subset_max_difference)
print(manifest)
