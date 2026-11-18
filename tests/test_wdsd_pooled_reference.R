# Run from the repository root: Rscript tests/test_wdsd_pooled_reference.R
source("scripts/support_functions.R")
source("scripts/core/wdsd_pooled_reference.R")

reference <- wp2_pooled_reference(".")
stopifnot(identical(normalizePath(reference$provenance$figure_pdf, mustWork = TRUE),
                    normalizePath("docs/reference_figures/FCGL_Fig3.ai", mustWork = TRUE)),
          identical(reference$provenance$figure_pdf_md5,
                    unname(tools::md5sum("docs/reference_figures/FCGL_Fig3.ai"))))
stopifnot(identical(dim(reference$rna), c(139L, 43L)),
          identical(dim(reference$lipids), c(64L, 43L)),
          identical(colnames(reference$rna), colnames(reference$lipids)),
          identical(colnames(reference$rna), reference$metadata$mouse_id),
          sum(reference$metadata$sex == "female") == 21L,
          sum(reference$metadata$sex == "male") == 22L,
          nrow(reference$edges) == 8896L,
          sum(reference$edges$primary_edge) == 192L,
          sum(reference$edges$display_gene_edge) == 1600L,
          nrow(reference$figure3_edges) == 192L,
          length(unique(reference$figure3_edges$gene)) == 25L,
          all(as.integer(table(reference$figure3_edges$reference_lipid)) == 3L),
          min(abs(reference$figure3_edges$r_reference)) < 0.65,
          nrow(reference$compressed_gene_lipid_edges) == 75L,
          reference$provenance$n_excluded_TF_gene_edges == 75L,
          reference$permissive_inventory$n_edges == 636732L,
          !reference$permissive_inventory$included_in_conservation,
          all(reference$replay_metrics$log2_max_abs_archive_difference < 1e-10),
          reference$provenance$raw_max_abs_archive_difference > 0.1)

# Independent scalar calculations cover every displayed species-level pair,
# without relying on the vectorized correlation lookup in the loader.
for (i in seq_len(nrow(reference$figure3_edges))) {
  edge <- reference$figure3_edges[i, ]
  direct <- stats::cor(reference$rna[edge$gene, ],
                       reference$lipids[edge$reference_lipid, ])
  stopifnot(abs(direct - edge$r_reference) < 1e-10)
}
keys <- function(edges) paste(edges$gene, edges$reference_lipid, sep = "\r")
stopifnot(setequal(keys(reference$edges[reference$edges$primary_edge, ]),
                   keys(reference$figure3_edges)))
archived <- utils::read.csv(reference$source_paths[["figure3_species_SIF"]],
                           check.names = FALSE)
is_lipid <- archived$node1 %in% rownames(reference$lipids)
stopifnot(sum(is_lipid) == 192L,
          setequal(paste(archived$node2[is_lipid], archived$node1[is_lipid], sep = "\r"),
                   keys(reference$figure3_edges)),
          !any(archived$node1[!is_lipid] %in% rownames(reference$lipids)))

# Jointly reordering mice preserves every coefficient. A one-modality-only
# permutation with unchanged labels must fail numerical archival replay; a
# mismatched column roster must fail the explicit pairing guard.
reorder <- rev(seq_len(ncol(reference$rna)))
joint_delta <- wp2_pooled_replay(reference$rna[, reorder],
                                 reference$lipids[, reorder], reference$edges)
stopifnot(joint_delta < 1e-10)
expect_error <- function(expr, text) {
  result <- tryCatch({ force(expr); NULL }, error = function(e) e$message)
  stopifnot(is.character(result), grepl(text, result, fixed = TRUE))
}
wrong_pair <- reference$lipids[, reorder]
colnames(wrong_pair) <- colnames(reference$rna)
expect_error(wp2_pooled_replay(reference$rna, wrong_pair, reference$edges),
             "fail archival replay")
expect_error(wp2_pooled_replay(reference$rna, reference$lipids[, reorder], reference$edges),
             "identically ordered unique mice")
duplicate <- reference$edges[c(1L, 1L), ]
expect_error(wp2_pooled_replay(reference$rna, reference$lipids, duplicate),
             "unique and finite")
missing <- reference$figure3_edges[1L, ]
missing$gene <- "Definitely_not_a_measured_gene"
expect_error(wp2_pooled_replay(reference$rna, reference$lipids, missing),
             "requires a measured gene and lipid")

# Class compression is a mean of original species correlations for each frozen
# class/gene membership, not a replacement abundance correlation.
species <- reference$figure3_edges
species$lipid_class <- sub("[(].*$", "", species$reference_lipid)
for (i in seq_len(nrow(reference$compressed_gene_lipid_edges))) {
  edge <- reference$compressed_gene_lipid_edges[i, ]
  selected <- species$gene == edge$gene & species$lipid_class == edge$lipid_class
  stopifnot(sum(selected) >= 1L,
            abs(mean(species$r_reference[selected]) - edge$r_reference) < 1e-10)
}
cat("WDSD pooled Figure 3 frozen-reference tests passed.\n")
