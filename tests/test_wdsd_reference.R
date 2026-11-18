# Run from the repository root: Rscript tests/test_wdsd_reference.R
source("scripts/support_functions.R")
source("scripts/core/wdsd_reference.R")

reference <- wp2_load_reference(".")
expected_contexts <- c("male_periportal", "male_pericentral",
                       "female_periportal", "female_pericentral")
stopifnot(identical(names(reference), expected_contexts))
summary <- attr(reference, "summary")
stopifnot(identical(summary$n_mice, c(22L, 22L, 21L, 21L)),
          identical(summary$n_genes, c(25L, 30L, 12L, 11L)),
          identical(summary$n_primary_edges, c(155L, 387L, 89L, 77L)),
          all(summary$n_lipids == 1254L),
          all(summary$log2_max_abs_archive_difference < 1e-10),
          all(summary$raw_max_abs_archive_difference > 0.1),
          all(summary$lipid_scale == "log2"))
for (context in names(reference)) {
  ref <- reference[[context]]
  stopifnot(identical(ref$metadata$mouse_id, colnames(ref$rna)),
            identical(colnames(ref$rna), colnames(ref$lipids)),
            !anyDuplicated(ref$metadata$mouse_id),
            all(ref$edges$gene %in% rownames(ref$rna)),
            all(ref$edges$reference_lipid %in% rownames(ref$lipids)),
            !anyDuplicated(paste(ref$edges$gene, ref$edges$reference_lipid)),
            nrow(ref$edges) == nrow(ref$rna) * nrow(ref$lipids))
  # Every archived coefficient is checked against the direct matrix Pearson
  # calculation, independent of the edge presentation order.
  direct <- stats::cor(t(ref$rna), t(ref$lipids), method = "pearson")
  lookup <- cbind(match(ref$edges$gene, rownames(direct)),
                  match(ref$edges$reference_lipid, colnames(direct)))
  stopifnot(max(abs(direct[lookup] - ref$edges$r_reference)) < 1e-10)
  selected <- unique(round(seq(1, nrow(ref$edges), length.out = 17L)))
  for (i in selected) {
    value <- stats::cor(ref$rna[ref$edges$gene[i], ],
                        ref$lipids[ref$edges$reference_lipid[i], ])
    stopifnot(abs(value - ref$edges$r_reference[i]) < 1e-10)
  }
  archived <- utils::read.csv(ref$source_paths[["displayed_SIF"]], check.names = FALSE)
  archived <- archived[archived$node2 %in% rownames(ref$lipids), , drop = FALSE]
  displayed <- ref$edges[ref$edges$primary_edge, , drop = FALSE]
  stopifnot(setequal(paste(archived$node1, archived$node2, sep = "\r"),
                     paste(displayed$gene, displayed$reference_lipid, sep = "\r")))
  # Jointly reordering mice retains coefficients; reordering only one modality
  # substantially changes them, verifying the pairing is consequential.
  permutation <- rev(seq_len(ncol(ref$rna)))
  jointly_reordered <- stats::cor(t(ref$rna[, permutation, drop = FALSE]),
                                  t(ref$lipids[, permutation, drop = FALSE]))
  unpaired <- stats::cor(t(ref$rna[, permutation, drop = FALSE]), t(ref$lipids))
  stopifnot(max(abs(jointly_reordered - direct)) < 1e-10,
            max(abs(unpaired - direct)) > 0.1)
}
# Pitx3 is an RNA gene with a TF annotation; its lipid branch is kept while
# the 559 male-pericentral TF-to-gene edges are independently excluded.
stopifnot(identical(reference$male_pericentral$RNA_genes_labeled_TF, "Pitx3"),
          sum(reference$male_pericentral$edges$gene == "Pitx3") == 1254L,
          reference$male_pericentral$replay_metrics$n_excluded_TF_gene_edges == 559L,
          summary$n_Mtmr11_full_edges[1L] == 1254L,
          summary$n_Mtmr11_primary_edges[1L] == 9L,
          all(summary$n_Mtmr11_full_edges[-1L] == 0L))
cat("WDSD original Figure 5 reference reconstruction tests passed.\n")
