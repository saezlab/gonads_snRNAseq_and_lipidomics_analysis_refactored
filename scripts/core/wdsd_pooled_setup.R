# Fixed whole-hepatocyte SD inputs; do not select edges on new correlations.
wp2_prepare_pooled_context <- function(data, reference, map, mapping = "strict") {
  stopifnot(all(c("reference_lipid", "feature_id", "mapping_rule") %in% names(map)),
            !anyDuplicated(map$reference_lipid), !anyDuplicated(map$feature_id))
  meta <- data$metadata[data$metadata$diet == "SD", , drop = FALSE]
  stopifnot(nrow(meta) == 14L, !anyDuplicated(meta$mouse_id),
            identical(as.integer(table(factor(meta$sex, c("female", "male")))), c(7L, 7L)))
  all <- reference$edges
  stopifnot(!anyDuplicated(paste(all$gene, all$reference_lipid, sep = "::")))
  all$edge_id <- paste(all$gene, all$reference_lipid, sep = "::")
  mi <- match(all$reference_lipid, map$reference_lipid)
  all$feature_id <- map$feature_id[mi]
  all$mapping_rule <- map$mapping_rule[mi]
  all$mapping <- mapping
  all$exclusion_reason <- ifelse(is.na(mi), "no_mapping_at_reported_resolution", "")
  all$exclusion_reason[!all$gene %in% rownames(data$rna_hepatocytes)] <- "gene_not_measured"
  all$exclusion_reason[!is.na(mi) & !all$feature_id %in% rownames(data$lipid)] <- "lipid_not_retained"
  mapped <- all[all$exclusion_reason == "", , drop = FALSE]
  lipmap <- unique(mapped[, c("reference_lipid", "feature_id")])
  stopifnot(!anyDuplicated(lipmap$reference_lipid), !anyDuplicated(lipmap$feature_id))
  genes <- unique(mapped$gene)
  rna <- data$rna_hepatocytes[genes, meta$mouse_id, drop = FALSE]
  lip <- data$lipid[lipmap$feature_id, meta$mouse_id, drop = FALSE]
  ref_rna <- reference$rna[genes, reference$metadata$mouse_id, drop = FALSE]
  ref_lip <- reference$lipids[lipmap$reference_lipid, reference$metadata$mouse_id, drop = FALSE]
  rownames(ref_lip) <- lipmap$feature_id
  context <- list(rna_SD = rna, lipid_SD = lip, metadata_SD = meta,
                  reference_rna = ref_rna, reference_lipid = ref_lip,
                  reference_metadata = reference$metadata, edges = mapped, mapping = mapping)
  checked <- wp2_pooled_context(context)
  jj <- match(all$edge_id, checked$edges$edge_id)
  all$r_SD <- checked$edges$r_SD[jj]
  all$eligible <- !is.na(jj) & checked$edges$eligible[jj]
  all$exclusion_reason[!is.na(jj)] <- checked$edges$exclusion_reason[jj[!is.na(jj)]]
  coverage <- do.call(rbind, lapply(c("primary", "full"), function(universe) {
    z <- all[universe == "full" | all$primary_edge, , drop = FALSE]
    keep <- z[z$eligible, , drop = FALSE]
    data.frame(mapping = mapping, universe = universe,
               original_edges = nrow(z), matched_edges = nrow(keep),
               original_lipids = length(unique(z$reference_lipid)),
               matched_lipids = length(unique(keep$reference_lipid)),
               original_genes = length(unique(z$gene)), matched_genes = length(unique(keep$gene)))
  }))
  list(context = context, coverage = coverage, ledger = all, lipid_map = lipmap)
}

wp2_pooled_annotate_edges <- function(edges, ledger) {
  ix <- match(edges$feature_id, ledger$feature_id)
  edges$canonical_lipid <- ledger$canonical_lipid[ix]
  edges$lipid_class <- ledger$lipid_class[ix]
  edges$qc_cv_percent <- ledger$qc_cv_percent[ix]
  edges$supplied_fill <- ledger[["Fill %"]][ix]
  edges
}
