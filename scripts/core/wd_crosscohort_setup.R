source("scripts/config/analysis_release.R")
# WP3 strict molecular join and fixed male WD context setup.
# Source support functions and the WP2 reference/setup/statistics helpers first.
wp3_prepare_contexts <- function(cf453, cf459, reference, root = ".") {
  required <- c("metadata", "lipid", "lipid_ledger", "rna")
  if (!all(required %in% names(cf453)) || !all(required %in% names(cf459))) {
    stop("Both cohorts require metadata, lipid, lipid_ledger and zonal rna.")
  }
  for (cohort in c("CF453", "CF459")) {
    x <- if (cohort == "CF453") cf453 else cf459
    m <- x$metadata
    if (!all(c("mouse_id", "sex", "diet") %in% names(m)) ||
        anyNA(m[, c("mouse_id", "sex", "diet")]) || anyDuplicated(m$mouse_id) ||
        !all(m$sex %in% c("male", "female")) || !all(m$diet %in% c("WD", "SD"))) {
      stop("Invalid cohort mouse/sex/diet metadata: ", cohort)
    }
    wp2_setup_named_matrix(x$lipid, paste(cohort, "lipid"))
    if (!identical(colnames(x$lipid), m$mouse_id)) stop("Lipid mouse order disagrees with metadata: ", cohort)
    for (zone in c("Periportal", "Pericentral")) {
      if (!zone %in% names(x$rna)) stop("Missing RNA zone: ", zone)
      wp2_setup_named_matrix(x$rna[[zone]], paste(cohort, zone))
      if (!identical(colnames(x$rna[[zone]]), m$mouse_id)) stop("RNA-lipid mouse order disagrees: ", cohort)
    }
    l <- x$lipid_ledger
    if (!all(c("feature_id", "canonical_lipid", "keep_primary", "exclusion_reason") %in% names(l)) ||
        anyNA(l[, c("feature_id", "canonical_lipid", "keep_primary")]) ||
        anyDuplicated(l$feature_id) || !is.logical(l$keep_primary) ||
        !setequal(rownames(x$lipid), l$feature_id[l$keep_primary])) stop("Invalid feature ledger: ", cohort)
  }
  expected <- c("male_periportal", "male_pericentral")
  if (!all(expected %in% names(reference))) stop("Both frozen male Figure 5 references are required.")
  labels <- unique(unlist(lapply(reference[expected], function(x) x$edges$reference_lipid), use.names = FALSE))
  canonical <- wp1_canonical_lipid(labels)
  if (anyNA(canonical) || any(!nzchar(canonical)) || anyDuplicated(canonical)) stop("Ambiguous frozen reference molecular identities.")
  crosswalk <- data.frame(reference_lipid = labels, canonical_lipid = canonical, stringsAsFactors = FALSE)
  for (cohort in c("CF453", "CF459")) {
    l <- if (cohort == "CF453") cf453$lipid_ledger else cf459$lipid_ledger
    counts <- table(l$canonical_lipid)
    n <- as.integer(counts[canonical]); n[is.na(n)] <- 0L
    index <- match(canonical, l$canonical_lipid)
    unique <- n == 1L
    feature <- l$feature_id[index]; feature[!unique] <- NA_character_
    kept <- unique & !is.na(index) & l$keep_primary[index]
    status <- ifelse(n == 0L, "missing_exact_annotation", ifelse(n > 1L, "ambiguous_exact_annotation", ifelse(kept, "unique_exact_retained", "detected_but_excluded")))
    crosswalk[[paste0("feature_id_", cohort)]] <- feature
    crosswalk[[paste0("annotation_", cohort)]] <- as.character(l[["Metabolite name"]][index])
    crosswalk[[paste0("multiplicity_", cohort)]] <- n
    crosswalk[[paste0("status_", cohort)]] <- status
    crosswalk[[paste0("kept_", cohort)]] <- kept
    crosswalk[[paste0("exclusion_", cohort)]] <- ifelse(is.na(index), "", l$exclusion_reason[index])
  }
  crosswalk$shared_primary <- crosswalk$kept_CF453 & crosswalk$kept_CF459
  a <- match(crosswalk$feature_id_CF453, cf453$lipid_ledger$feature_id)
  crosswalk$CF453_QC30 <- !is.na(a) & cf453$lipid_ledger$keep_qc30[a]
  fill <- cf453$lipid_ledger[["Fill %"]][a]
  crosswalk$CF453_Fill50 <- is.finite(fill) & fill >= 0.5
  crosswalk$matching_policy <- "unique strict final-pipe/class formatting; preserve molecular resolution; no Alignment-ID join"
  indices453 <- which(cf453$metadata$sex == "male" & cf453$metadata$diet == "WD")
  indices459 <- which(cf459$metadata$sex == "male" & cf459$metadata$diet == "WD")
  if (length(indices453) != 7L || length(indices459) != 4L ||
      !setequal(cf453$metadata$mouse_id[indices453], paste0("cf453-", 22:28)) ||
      !setequal(cf459$metadata$mouse_id[indices459], paste0("cf459-", 13:16))) stop("Expected verified seven CF453/four CF459 male WD mice.")
  contexts <- stats::setNames(vector("list", 2L), expected)
  ledgers <- coverage <- vector("list", 2L)
  for (i in seq_along(expected)) {
    name <- expected[i]; ref <- reference[[name]]
    zone <- if (i == 1L) "Periportal" else "Pericentral"
    e <- ref$edges
    if (!all(c("gene", "reference_lipid", "r_reference", "primary_edge") %in% names(e)) ||
        anyNA(e[, c("gene", "reference_lipid", "primary_edge")]) || !is.logical(e$primary_edge) ||
        anyDuplicated(paste(e$gene, e$reference_lipid, sep = "\r"))) stop("Invalid frozen edge selection: ", name)
    map <- crosswalk[match(e$reference_lipid, crosswalk$reference_lipid), , drop = FALSE]
    e$context <- name; e$sex <- "male"; e$zone <- zone
    e$feature_id <- e$reference_lipid
    e$feature_id_CF453 <- map$feature_id_CF453; e$feature_id_CF459 <- map$feature_id_CF459
    e$CF453_QC30 <- map$CF453_QC30; e$CF453_Fill50 <- map$CF453_Fill50
    e$lipid_class <- sub("\\(.*$", "", map$canonical_lipid)
    e$edge_id <- paste(name, e$gene, e$reference_lipid, sep = "::")
    e$not_strictly_shared <- !map$shared_primary
    e$nonfinite_reference <- !is.finite(e$r_reference)
    cohort_rna <- cohort_lipid <- list()
    for (cohort in c("CF453", "CF459")) {
      x <- if (cohort == "CF453") cf453 else cf459
      ix <- if (cohort == "CF453") indices453 else indices459
      rn <- x$rna[[zone]][, ix, drop = FALSE]; lp <- x$lipid[, ix, drop = FALSE]
      gf <- wp2_setup_row_flags(rn); lf <- wp2_setup_row_flags(lp)
      e[[paste0("missing_gene_", cohort)]] <- !e$gene %in% rownames(rn)
      lookup <- function(flag_table, ids, field) { z <- flag_table[[field]][match(ids, flag_table$feature)]; z[is.na(z)] <- FALSE; z }
      for (field in c("nonfinite", "invariant")) {
        e[[paste0(field, "_gene_", cohort)]] <- lookup(gf, e$gene, field)
        e[[paste0(field, "_lipid_", cohort)]] <- lookup(lf, e[[paste0("feature_id_", cohort)]], field)
      }
      cohort_rna[[cohort]] <- rn; cohort_lipid[[cohort]] <- lp
    }
    flags <- c("not_strictly_shared", "nonfinite_reference", unlist(lapply(c("CF453", "CF459"), function(x) paste0(c("missing_gene_", "nonfinite_gene_", "invariant_gene_", "nonfinite_lipid_", "invariant_lipid_"), x))))
    e$eligible <- rowSums(as.matrix(e[, flags, drop = FALSE])) == 0L
    e$exclusion_reason <- ""
    for (flag in flags) e$exclusion_reason[e[[flag]]] <- ifelse(nzchar(e$exclusion_reason[e[[flag]]]), paste(e$exclusion_reason[e[[flag]]], flag, sep = ";"), flag)
    use <- e[e$eligible, , drop = FALSE]
    if (nrow(use) < 3L || sum(use$primary_edge) < 3L) stop("Insufficient fixed eligible full/displayed edges: ", name)
    genes <- unique(use$gene); features <- unique(use$feature_id)
    mat <- function(cohort) {
      fid <- use[[paste0("feature_id_", cohort)]][match(features, use$feature_id)]
      ans <- cohort_lipid[[cohort]][fid, , drop = FALSE]; rownames(ans) <- features; ans
    }
    contexts[[name]] <- list(sex = "male", zone = zone, edges = use,
      rna_CF453 = cohort_rna$CF453[genes, , drop = FALSE], lipid_CF453 = mat("CF453"),
      rna_CF459 = cohort_rna$CF459[genes, , drop = FALSE], lipid_CF459 = mat("CF459"),
      metadata_CF453 = cf453$metadata[indices453, , drop = FALSE],
      metadata_CF459 = cf459$metadata[indices459, , drop = FALSE],
      eligibility_policy = "fixed FCG full/displayed membership; strict shared finite/varying male WD features only; no SD eligibility requirement")
    ledgers[[i]] <- e
    count <- function(mask, column) length(unique(e[[column]][mask]))
    coverage[[i]] <- data.frame(context = name, zone = zone, n_CF453 = 7L, n_CF459 = 4L,
      original_full_edges = nrow(e), original_displayed_edges = sum(e$primary_edge),
      eligible_full_edges = sum(e$eligible), eligible_displayed_edges = sum(e$eligible & e$primary_edge),
      eligible_genes = count(e$eligible, "gene"), eligible_lipids = count(e$eligible, "reference_lipid"),
      eligible_displayed_genes = count(e$eligible & e$primary_edge, "gene"),
      eligible_displayed_lipids = count(e$eligible & e$primary_edge, "reference_lipid"),
      nonshared_edge_exclusions = sum(e$not_strictly_shared), other_edge_exclusions = sum(!e$eligible & !e$not_strictly_shared))
    cat("WP3 ", name, ": ", nrow(use), "/", nrow(e), " full; ", sum(use$primary_edge), "/", sum(e$primary_edge), " displayed pairs; ", length(features), " strict shared lipids.\n", sep = "")
  }
  for (cohort in c("CF453", "CF459")) if (!identical(colnames(contexts[[1]][[paste0("rna_", cohort)]]), colnames(contexts[[2]][[paste0("rna_", cohort)]]))) stop("Zones must share identical mouse order.")
  # Independently reconstruct CF453 WD coefficients and compare existing WP2
  # wherever that earlier SD+WD eligible universe overlaps the new WD-only one.
  prior_path <- file.path(analysis_path("historical_wf12",root),"network_edges.csv")
  replay <- vector("list", 2L)
  old <- if (file.exists(prior_path)) utils::read.csv(prior_path, check.names = FALSE) else NULL
  for (i in seq_along(contexts)) {
    x <- contexts[[i]]; e <- x$edges
    observed <- wp2_cor_edges(x$rna_CF453, x$lipid_CF453, e)
    if (!is.null(old)) {
      p <- old[old$context == names(contexts)[i], , drop = FALSE]
      ix <- match(paste(e$gene, e$reference_lipid, sep = "\r"), paste(p$gene, p$reference_lipid, sep = "\r"))
      ok <- !is.na(ix); diff <- if (any(ok)) max(abs(observed[ok] - p$r_WD[ix[ok]])) else NA_real_
      if (any(ok) && (!is.finite(diff) || diff > 1e-10)) stop("CF453 WD does not replay WP2: ", names(contexts)[i])
      replay[[i]] <- data.frame(context = names(contexts)[i], WP2_available = TRUE, overlapping_edges = sum(ok), new_WD_only_eligible_edges = sum(!ok), max_abs_difference = diff)
    } else replay[[i]] <- data.frame(context = names(contexts)[i], WP2_available = FALSE, overlapping_edges = 0L, new_WD_only_eligible_edges = nrow(e), max_abs_difference = NA_real_)
  }
  list(contexts = contexts, lipid_crosswalk = crosswalk,
       edge_eligibility_ledger = do.call(rbind, ledgers), coverage = do.call(rbind, coverage),
       WP2_replay = do.call(rbind, replay), WP2_replay_path = if (is.null(old)) character() else normalizePath(prior_path))
}

wp3_context_sensitivity <- function(contexts, cf453, cf459, type) {
  if (!type %in% c("CF453_QC30", "CF453_Fill50", "RNA_data")) stop("Unknown prespecified WP3 sensitivity.")
  lapply(contexts, function(x) {
    if (type == "RNA_data") {
      x$rna_CF453 <- cf453$rna_data[[x$zone]][rownames(x$rna_CF453), x$metadata_CF453$mouse_id, drop = FALSE]
      x$rna_CF459 <- cf459$rna_data[[x$zone]][rownames(x$rna_CF459), x$metadata_CF459$mouse_id, drop = FALSE]
    } else x$edges <- x$edges[x$edges[[type]], , drop = FALSE]
    x$sensitivity <- type; x
  })
}
