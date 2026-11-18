# Fixed WP2 molecular joins and edge eligibility. Source support functions first.
# Source cyp2c38_reference.R for the shared formatting-only canonicalizer.

wp2_setup_named_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) || any(!nzchar(rownames(x))) ||
      any(!nzchar(colnames(x))) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(label, " must be a numeric feature-by-mouse matrix with unique, nonblank names.")
  }
  invisible(TRUE)
}

# Match the statistics helper's relative numerical constancy tolerance exactly.
# Report nonfinite measurements separately from finite constant measurements.
wp2_setup_row_flags <- function(x) {
  finite <- rowSums(!is.finite(x)) == 0L
  magnitude <- apply(abs(x), 1L, max)
  safe_magnitude <- magnitude
  safe_magnitude[!is.finite(safe_magnitude) | safe_magnitude == 0] <- 1
  scaled <- x / safe_magnitude
  centered <- scaled - rowMeans(scaled)
  ss <- rowSums(centered^2)
  invariant <- finite & (ss <= ncol(x) * .Machine$double.eps)
  data.frame(feature = rownames(x), nonfinite = !finite,
             invariant = invariant, stringsAsFactors = FALSE)
}

wp2_prepare_contexts <- function(data, reference) {
  if (!exists("wp1_canonical_lipid", mode = "function")) {
    stop("Source scripts/core/cyp2c38_reference.R before comparison setup.")
  }
  needed <- c("metadata", "lipid", "lipid_ledger", "rna")
  if (!all(needed %in% names(data))) stop("WDSD intake list lacks required named fields.")
  m <- data$metadata
  if (!all(c("mouse_id", "sex", "diet", "group") %in% names(m)) ||
      nrow(m) != 28L || anyNA(m[, c("mouse_id", "sex", "diet", "group")]) ||
      anyDuplicated(m$mouse_id) || !setequal(m$mouse_id, paste0("cf453-", seq_len(28L)))) {
    stop("WDSD metadata must contain 28 distinct verified CF453 mouse IDs and complete sex/diet labels.")
  }
  expected_groups <- rep(c("female-sd", "female-wd", "male-sd", "male-wd"), each = 7L)
  names(expected_groups) <- paste0("cf453-", seq_len(28L))
  if (!identical(as.character(m$group), unname(expected_groups[m$mouse_id])) ||
      !identical(as.character(m$group), paste(m$sex, tolower(m$diet), sep = "-")) ||
      !all(m$sex %in% c("female", "male")) || !all(m$diet %in% c("SD", "WD"))) {
    stop("WDSD metadata sex/diet labels disagree with the verified CF453 roster.")
  }
  wp2_setup_named_matrix(data$lipid, "WDSD lipid")
  if (!identical(colnames(data$lipid), as.character(m$mouse_id))) {
    stop("WDSD lipid columns must be identical to the ordered metadata mouse IDs.")
  }
  if (!all(c("Periportal", "Pericentral") %in% names(data$rna))) stop("Both WDSD RNA zones are required.")
  for (zone in c("Periportal", "Pericentral")) {
    wp2_setup_named_matrix(data$rna[[zone]], paste("WDSD RNA", zone))
    if (!identical(colnames(data$rna[[zone]]), colnames(data$lipid))) {
      stop("WDSD RNA-lipid columns must identify the same mice in identical order in each zone.")
    }
  }
  ledger <- data$lipid_ledger
  fields <- c("feature_id", "canonical_lipid", "keep_primary", "keep_qc30", "exclusion_reason")
  if (!all(fields %in% names(ledger)) || anyNA(ledger[, fields]) ||
      anyDuplicated(ledger$feature_id) || anyDuplicated(ledger$canonical_lipid) ||
      any(!nzchar(ledger$canonical_lipid)) || !is.logical(ledger$keep_primary) ||
      !is.logical(ledger$keep_qc30) || any(ledger$keep_qc30 & !ledger$keep_primary)) {
    stop("Lipid ledger must have unique feature/canonical identities and complete logical eligibility fields.")
  }
  if (!setequal(rownames(data$lipid), ledger$feature_id[ledger$keep_primary])) {
    stop("WDSD lipid matrix rows do not match the retained annotation ledger exactly.")
  }
  context_names <- c("male_periportal", "male_pericentral", "female_periportal", "female_pericentral")
  if (!is.list(reference) || !identical(names(reference), context_names)) {
    stop("Expected the four named frozen Figure 5 reference contexts in documented order.")
  }
  ref_labels <- unique(unlist(lapply(reference, function(x) as.character(x$edges$reference_lipid)), use.names = FALSE))
  ref_canonical <- wp1_canonical_lipid(ref_labels)
  if (anyNA(ref_canonical) || any(!nzchar(ref_canonical)) || anyDuplicated(ref_canonical)) {
    stop("Reference lipid annotations are ambiguous after formatting-only canonicalization.")
  }
  join <- match(ref_canonical, ledger$canonical_lipid)
  known <- !is.na(join)
  crosswalk <- data.frame(reference_lipid = ref_labels, canonical_lipid = ref_canonical,
                          feature_id = ledger$feature_id[join],
                          matched_annotation = if ("Metabolite name" %in% names(ledger)) ledger[["Metabolite name"]][join] else NA_character_,
                          matched_primary = ifelse(known, ledger$keep_primary[join], FALSE),
                          matched_qc30 = ifelse(known, ledger$keep_qc30[join], FALSE),
                          input_exclusion_reason = ledger$exclusion_reason[join],
                          stringsAsFactors = FALSE)
  crosswalk$match_status <- ifelse(!known, "missing_exact_molecular_annotation",
                                  ifelse(crosswalk$matched_primary, "unique_exact", "detected_but_excluded_by_intake"))
  crosswalk$matching_policy <- "last_pipe_annotation_and_class_formatting_only"
  crosswalk$input_exclusion_reason[!known] <- ""
  contexts <- stats::setNames(vector("list", length(reference)), names(reference))
  all_ledgers <- vector("list", length(reference))
  coverage <- vector("list", length(reference))
  for (i in seq_along(reference)) {
    name <- names(reference)[i]
    ref <- reference[[i]]
    sex <- sub("_.*$", "", name)
    zone <- if (endsWith(name, "periportal")) "Periportal" else "Pericentral"
    if (!all(c("edges", "rna", "lipids", "metadata") %in% names(ref))) stop("Incomplete reference context: ", name)
    e <- ref$edges
    edge_fields <- c("gene", "reference_lipid", "r_reference", "primary_edge")
    if (!all(edge_fields %in% names(e)) || anyNA(e[, c("gene", "reference_lipid", "primary_edge")]) ||
        !is.logical(e$primary_edge) || anyDuplicated(paste(e$gene, e$reference_lipid, sep = "\r"))) {
      stop("Reference edges require unique gene-lipid identities and complete displayed membership.")
    }
    wp2_setup_named_matrix(ref$rna, paste("Reference RNA", name))
    wp2_setup_named_matrix(ref$lipids, paste("Reference lipids", name))
    if (!all(c("mouse_id", "group") %in% names(ref$metadata)) || anyNA(ref$metadata$mouse_id) ||
        anyNA(ref$metadata$group) || anyDuplicated(ref$metadata$mouse_id) ||
        !identical(colnames(ref$rna), colnames(ref$lipids)) ||
        !identical(colnames(ref$rna), as.character(ref$metadata$mouse_id))) {
      stop("Original reference metadata/RNA/lipid pairing is inconsistent in ", name, ".")
    }
    if (!all(e$gene %in% rownames(ref$rna)) || !all(e$reference_lipid %in% rownames(ref$lipids))) {
      stop("Reference edges lack measured source features in ", name, ".")
    }
    if ("sex" %in% names(ref$metadata) && any(ref$metadata$sex != sex)) stop("Reference sex labels disagree with context.")
    if ("zone" %in% names(ref$metadata) && any(ref$metadata$zone != zone)) stop("Reference zone labels disagree with context.")
    reference_groups <- if (sex == "male") c("XXT", "XYT") else c("XXO", "XYO")
    reference_group_counts <- table(factor(ref$metadata$group, levels = reference_groups))
    if (any(!ref$metadata$group %in% reference_groups) ||
        !identical(as.integer(reference_group_counts), if (sex == "male") c(11L, 11L) else c(10L, 11L))) {
      stop("Reference gonadal/chromosome groups differ from the frozen original cohort in ", name, ".")
    }
    if (anyDuplicated(wp1_canonical_lipid(rownames(ref$lipids)))) stop("Ambiguous canonical original matrix lipids.")
    SD <- which(m$sex == sex & m$diet == "SD")
    WD <- which(m$sex == sex & m$diet == "WD")
    if (length(SD) != 7L || length(WD) != 7L) stop("Each WDSD sex/diet stratum must contain seven mice.")
    rna_SD <- data$rna[[zone]][, SD, drop = FALSE]
    rna_WD <- data$rna[[zone]][, WD, drop = FALSE]
    lipid_SD <- data$lipid[, SD, drop = FALSE]
    lipid_WD <- data$lipid[, WD, drop = FALSE]
    gSD <- wp2_setup_row_flags(rna_SD); gWD <- wp2_setup_row_flags(rna_WD)
    lSD <- wp2_setup_row_flags(lipid_SD); lWD <- wp2_setup_row_flags(lipid_WD)
    map <- crosswalk[match(e$reference_lipid, crosswalk$reference_lipid), , drop = FALSE]
    e$context <- name; e$sex <- sex; e$zone <- zone
    e$canonical_lipid <- map$canonical_lipid
    e$feature_id <- map$feature_id
    e$edge_id <- paste(name, e$gene, e$reference_lipid, sep = "::")
    e$missing_gene <- !e$gene %in% rownames(rna_SD)
    e$missing_lipid <- map$match_status == "missing_exact_molecular_annotation"
    e$excluded_lipid <- map$match_status == "detected_but_excluded_by_intake"
    e$input_lipid_exclusion_reason <- map$input_exclusion_reason
    e$nonfinite_reference <- !is.finite(e$r_reference)
    lookup_flag <- function(table, ids, key) {
      result <- table[[key]][match(ids, table$feature)]
      result[is.na(result)] <- FALSE
      result
    }
    for (diet in c("SD", "WD")) {
      gf <- if (diet == "SD") gSD else gWD
      lf <- if (diet == "SD") lSD else lWD
      e[[paste0("nonfinite_gene_", diet)]] <- lookup_flag(gf, e$gene, "nonfinite")
      e[[paste0("invariant_gene_", diet)]] <- lookup_flag(gf, e$gene, "invariant")
      e[[paste0("nonfinite_lipid_", diet)]] <- lookup_flag(lf, e$feature_id, "nonfinite")
      e[[paste0("invariant_lipid_", diet)]] <- lookup_flag(lf, e$feature_id, "invariant")
    }
    flags <- c("missing_gene", "missing_lipid", "excluded_lipid", "nonfinite_reference",
               "nonfinite_gene_SD", "invariant_gene_SD", "nonfinite_lipid_SD", "invariant_lipid_SD",
               "nonfinite_gene_WD", "invariant_gene_WD", "nonfinite_lipid_WD", "invariant_lipid_WD")
    e$eligible <- rowSums(as.matrix(e[, flags, drop = FALSE])) == 0L
    e$exclusion_reason <- ""
    for (flag in flags) {
      hit <- e[[flag]]
      e$exclusion_reason[hit] <- ifelse(nzchar(e$exclusion_reason[hit]),
                                      paste(e$exclusion_reason[hit], flag, sep = ";"), flag)
    }
    e$keep_qc30 <- map$matched_qc30
    e$n_SD <- length(SD); e$n_WD <- length(WD)
    all_ledgers[[i]] <- e
    use <- e[e$eligible, , drop = FALSE]
    genes <- unique(use$gene); features <- unique(use$feature_id)
    ref_lipids <- ref$lipids[use$reference_lipid[match(features, use$feature_id)], , drop = FALSE]
    rownames(ref_lipids) <- features
    contexts[[i]] <- list(
      sex = sex, zone = zone,
      rna_SD = rna_SD[genes, , drop = FALSE], rna_WD = rna_WD[genes, , drop = FALSE],
      lipid_SD = lipid_SD[features, , drop = FALSE], lipid_WD = lipid_WD[features, , drop = FALSE],
      metadata_SD = m[SD, , drop = FALSE], metadata_WD = m[WD, , drop = FALSE],
      edges = use,
      reference_rna = ref$rna[genes, , drop = FALSE], reference_lipid = ref_lipids,
      reference_metadata = ref$metadata,
      eligibility_policy = "exact identity; finite and varying features in both diets; fixed for stages A and B"
    )
    primary <- e$primary_edge
    eligible_primary <- e$eligible & primary
    count_features <- function(mask, field) length(unique(e[[field]][mask]))
    d <- data.frame(context = name, sex = sex, zone = zone,
                    n_reference_mice = ncol(ref$rna), n_SD = length(SD), n_WD = length(WD),
                    n_original_full_edges = nrow(e), n_original_primary_edges = sum(primary),
                    n_original_full_genes = count_features(rep(TRUE, nrow(e)), "gene"),
                    n_original_full_lipids = count_features(rep(TRUE, nrow(e)), "reference_lipid"),
                    n_original_primary_genes = count_features(primary, "gene"),
                    n_original_primary_lipids = count_features(primary, "reference_lipid"),
                    n_eligible_full_edges = sum(e$eligible), n_eligible_primary_edges = sum(eligible_primary),
                    n_eligible_full_genes = count_features(e$eligible, "gene"),
                    n_eligible_full_lipids = count_features(e$eligible, "reference_lipid"),
                    n_eligible_primary_genes = count_features(eligible_primary, "gene"),
                    n_eligible_primary_lipids = count_features(eligible_primary, "reference_lipid"),
                    n_eligible_QC30_full_edges = sum(e$eligible & e$keep_qc30),
                    n_eligible_QC30_primary_edges = sum(eligible_primary & e$keep_qc30),
                    stringsAsFactors = FALSE)
    for (flag in flags) d[[paste0("n_", flag, "_edges")]] <- sum(e[[flag]])
    coverage[[i]] <- d
    cat("WDSD ", name, ": ", nrow(use), "/", nrow(e), " full and ",
        sum(use$primary_edge), "/", sum(primary), " displayed eligible edges; ",
        length(genes), " genes, ", length(features), " lipids; 7 SD + 7 WD mice.\n", sep = "")
  }
  for (sex in c("male", "female")) {
    portal <- contexts[[paste0(sex, "_periportal")]]
    central <- contexts[[paste0(sex, "_pericentral")]]
    if (!identical(portal$reference_metadata$mouse_id, central$reference_metadata$mouse_id) ||
        !identical(portal$reference_metadata$group, central$reference_metadata$group)) {
      stop("Original reference zones must share the same ordered mouse roster and group strata.")
    }
  }
  list(contexts = contexts, lipid_crosswalk = crosswalk,
       edge_eligibility_ledger = do.call(rbind, all_ledgers), coverage = do.call(rbind, coverage),
       eligibility_policy = "Freeze exact measurable, finite, varying edges in BOTH diets before stages A/B; no WDSD correlation threshold")
}

wp2_contexts_qc30 <- function(contexts, lipid_ledger) {
  if (!all(c("feature_id", "keep_qc30") %in% names(lipid_ledger)) ||
      anyDuplicated(lipid_ledger$feature_id) || anyNA(lipid_ledger$keep_qc30)) {
    stop("QC30 sensitivity requires unique feature IDs and complete QC flags.")
  }
  retained <- lipid_ledger$feature_id[lipid_ledger$keep_qc30]
  lapply(contexts, function(x) {
    x$edges <- x$edges[x$edges$feature_id %in% retained, , drop = FALSE]
    x$eligibility_policy <- paste(x$eligibility_policy, "; prespecified QC-CV <=30% subset")
    x
  })
}

wp2_contexts_data_layer <- function(contexts, data) {
  if (is.null(data$rna_data)) stop("Stored data-layer sensitivity requires intake rna_data matrices.")
  lapply(contexts, function(x) {
    alternate <- data$rna_data[[x$zone]]
    wp2_setup_named_matrix(alternate, "Stored RNA data-layer sensitivity")
    if (!identical(colnames(alternate), data$metadata$mouse_id)) stop("Data-layer matrix is misaligned with metadata.")
    genes <- rownames(x$rna_SD)
    if (!all(genes %in% rownames(alternate))) stop("Data-layer sensitivity lacks frozen eligible genes.")
    x$rna_SD <- alternate[genes, x$metadata_SD$mouse_id, drop = FALSE]
    x$rna_WD <- alternate[genes, x$metadata_WD$mouse_id, drop = FALSE]
    x$RNA_sensitivity <- "stored_data_layer; frozen primary feature/edge universe"
    x
  })
}
