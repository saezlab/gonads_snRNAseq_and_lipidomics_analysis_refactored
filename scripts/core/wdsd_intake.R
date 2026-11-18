source("scripts/config/analysis_release.R",local=TRUE)
# WP2 WDSD intake and audit. Source scripts/support_functions.R before this file.
# Functions only read the supplied inputs; no normalization, filtering by outcome,
# RNA re-annotation or nucleus-count reconstruction is performed.

wp2_lipid_sample_metadata <- function(columns) {
  samples <- grep("^CF453-[0-9]+_[FM]-(SD|WD)$", columns, value = TRUE)
  if (length(samples) != 28L || anyDuplicated(samples)) {
    stop("Expected exactly 28 uniquely named CF453 biological lipid columns.")
  }
  metadata <- data.frame(
    mouse_id = tolower(sub("_.*$", "", samples)),
    sex = ifelse(grepl("_F-", samples, fixed = TRUE), "female", "male"),
    diet = sub("^.*-", "", samples), lipid_sample = samples,
    stringsAsFactors = FALSE
  )
  metadata$group <- paste(metadata$sex, tolower(metadata$diet), sep = "-")
  expected_ids <- paste0("cf453-", seq_len(28L))
  if (anyDuplicated(metadata$mouse_id) || !setequal(metadata$mouse_id, expected_ids)) {
    stop("Biological lipid columns must identify cf453-1 through cf453-28 exactly once.")
  }
  metadata <- metadata[match(expected_ids, metadata$mouse_id), , drop = FALSE]
  rownames(metadata) <- NULL
  expected_groups <- rep(c("female-sd", "female-wd", "male-sd", "male-wd"), each = 7L)
  if (!identical(metadata$group, expected_groups)) {
    stop("Lipid-column sex/diet labels disagree with the verified CF453 mouse roster.")
  }
  metadata
}

# Scale the exported data layer across its original complete profile roster.
# Constant genes stay zero, matching the supplied ScaleData output. This is an
# explicit sensitivity adapter, not a reconstruction of upstream normalization.
wp2_row_scale <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) < 2L || any(!is.finite(x))) {
    stop("Row scaling requires a finite numeric matrix with at least two profiles.")
  }
  centered <- x - rowMeans(x)
  row_sd <- sqrt(rowSums(centered^2) / (ncol(x) - 1L))
  result <- centered
  nonconstant <- row_sd > 0
  result[nonconstant, ] <- centered[nonconstant, , drop = FALSE] / row_sd[nonconstant]
  result[!nonconstant, ] <- 0
  dimnames(result) <- dimnames(x)
  result
}

wp2_align_rna_zone <- function(layer, profile_metadata, mouse_ids, zone) {
  required <- c("profile_id", "mouse_id", "zonation_type")
  if (!all(required %in% names(profile_metadata)) ||
      !zone %in% c("Periportal", "Pericentral")) {
    stop("Zone alignment requires named profile/mouse/zonation metadata.")
  }
  if (anyDuplicated(rownames(layer)) || anyDuplicated(colnames(layer)) ||
      !setequal(colnames(layer), profile_metadata$profile_id)) {
    stop("RNA layer profile identifiers are missing, duplicated or unmatched.")
  }
  rows <- profile_metadata[profile_metadata$zonation_type == zone, , drop = FALSE]
  if (anyDuplicated(rows$mouse_id) || !setequal(rows$mouse_id, mouse_ids)) {
    stop("Each mouse must contribute exactly one profile for each zone.")
  }
  ordered_profiles <- rows$profile_id[match(mouse_ids, rows$mouse_id)]
  matrix <- as.matrix(layer[, match(ordered_profiles, colnames(layer)), drop = FALSE])
  colnames(matrix) <- mouse_ids
  matrix
}

wp2_load_inputs <- function(root = ".", input_zonation = NULL, input_hepatocytes = NULL) {
  if (length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop("root must identify one existing analysis repository directory.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  source(file.path(root,"scripts/core/pseudobulk_csv.R"),local=TRUE)
  required_packages <- "Matrix"
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  }
  input_paths <- c(
    lipidomics = file.path(root, "data/metabolomics/wdsd/lipidomics_cf453.csv"),
    hepatocytes = analysis_path("wdsd_hepatocytes_transport",root),
    zonation = analysis_path("current_wdsd_zonation_transport",root)
  )
  if (!is.null(input_zonation)) {
    if (length(input_zonation) != 1L || is.na(input_zonation) || !nzchar(input_zonation))
      stop("input_zonation must identify one saved zonation pseudobulk file.")
    input_paths[["zonation"]] <- input_zonation
  }
  if(!is.null(input_hepatocytes))input_paths[["hepatocytes"]]<-input_hepatocytes
  absent <- input_paths[!file.exists(input_paths)]
  if (length(absent)) stop("Missing WDSD input(s): ", paste(absent, collapse = ", "))
  input_paths <- vapply(input_paths, normalizePath, character(1L), mustWork = TRUE)
  canonical_helpers <- new.env(parent = baseenv())
  sys.source(file.path(root, "scripts/core/cyp2c38_reference.R"), envir = canonical_helpers)

  raw <- utils::read.csv(input_paths[["lipidomics"]], check.names = FALSE,
                         stringsAsFactors = FALSE)
  if (anyDuplicated(names(raw)) || !all(c("Metabolite name", "Alignment ID",
                                         "Adduct type", "Ontology") %in% names(raw))) {
    stop("Lipidomics requires unique column names and the named annotation fields.")
  }
  metadata <- wp2_lipid_sample_metadata(names(raw))
  qc_columns <- grep("^QC[0-9]+$", names(raw), value = TRUE)
  if (length(qc_columns) != 7L) stop("Expected seven WDSD analytical QC injections.")
  lipid_samples <- metadata$lipid_sample
  x_all <- as.matrix(raw[, lipid_samples, drop = FALSE])
  qc_matrix <- as.matrix(raw[, qc_columns, drop = FALSE])
  if (!is.numeric(x_all) || !is.numeric(qc_matrix) || any(!is.finite(x_all)) ||
      any(!is.finite(qc_matrix)) || any(x_all <= 0) || any(qc_matrix <= 0)) {
    stop("Supplied biological and QC intensities must be numeric, finite and positive.")
  }
  names_raw <- raw[["Metabolite name"]]
  if (anyNA(names_raw) || any(!nzchar(trimws(names_raw)))) stop("Empty lipid annotation.")
  ledger <- raw[, setdiff(names(raw), c(lipid_samples, qc_columns)), drop = FALSE]
  ledger$feature_id <- sprintf("WDSD_%04d", seq_len(nrow(raw)))
  ledger$canonical_lipid <- canonical_helpers$wp1_canonical_lipid(names_raw)
  ledger$lipid_class <- sub("\\(.*$", "", ledger$canonical_lipid)
  ledger$internal_standard <- grepl("^IS:", names_raw) | raw$Ontology == "IS"
  ledger$isotope_labelled <- grepl("\\([dD][0-9]+\\)|\\b[0-9]+[Cc]\\b", names_raw)
  excluded_classes <- c("BA", "DGGA", "HBMP", "LPE-N", "NATau", "PEtOH", "PMeOH", "SL", "CAR", "SPB", "VAE")
  ledger$excluded_class <- ledger$lipid_class %in% excluded_classes
  ledger$exclusion_reason <- ifelse(ledger$internal_standard, "internal_standard",
                                   ifelse(ledger$isotope_labelled, "isotope_labelled_feature",
                                          ifelse(ledger$excluded_class, "repository_excluded_class", "")))
  ledger$keep_primary <- ledger$exclusion_reason == ""
  ledger$qc_cv_percent <- 100 * apply(qc_matrix, 1L, stats::sd) / rowMeans(qc_matrix)
  ledger$keep_qc30 <- ledger$keep_primary & ledger$qc_cv_percent <= 30
  multiplicity <- function(value) as.integer(table(value)[as.character(value)])
  ledger$alignment_id_multiplicity <- multiplicity(ledger[["Alignment ID"]])
  ledger$alignment_adduct_multiplicity <- multiplicity(paste(ledger[["Alignment ID"]], ledger[["Adduct type"]], sep = "__"))
  ledger$original_annotation_multiplicity <- multiplicity(names_raw)
  ledger$canonical_annotation_multiplicity <- multiplicity(ledger$canonical_lipid)
  ledger$pipe_delimited_annotation <- grepl("|", names_raw, fixed = TRUE)
  ledger$annotation_join_status <- ifelse(ledger$canonical_annotation_multiplicity == 1L,
                                           "unique_canonical_annotation", "ambiguous_canonical_annotation")
  ledger$identity_policy <- "last_pipe_annotation_and_class_formatting_only"
  lipid <- log2(x_all[ledger$keep_primary, , drop = FALSE])
  dimnames(lipid) <- list(ledger$feature_id[ledger$keep_primary], metadata$mouse_id)
  lipid_supplied <- x_all[ledger$keep_primary, , drop = FALSE]
  dimnames(lipid_supplied) <- dimnames(lipid)
  if (anyDuplicated(rownames(lipid)) || anyDuplicated(ledger$canonical_lipid[ledger$keep_primary])) {
    stop("Retained lipid canonical annotations must be unique; ambiguous rows cannot be forced into reference joins.")
  }
  metadata$total_supplied_lipid_intensity <- colSums(x_all)
  metadata$total_primary_lipid_intensity <- colSums(lipid_supplied)
  metadata$median_primary_log2_intensity <- apply(lipid, 2L, stats::median)
  analytical_qc <- data.frame(qc_injection = qc_columns,
                              total_supplied_intensity = colSums(qc_matrix),
                              total_primary_intensity = colSums(qc_matrix[ledger$keep_primary, , drop = FALSE]),
                              stringsAsFactors = FALSE)

  read_rna <- function(path, kind) {
    object <- pb_components(path)
    expected_layers <- c("counts", "data", "scale.data")
    if (!all(expected_layers %in% names(object$layers))) stop("Missing RNA layer.")
    m <- object$meta.data
    if (!all(c("orig.ident", "cf453_id", "group") %in% names(m)) ||
        anyNA(m$cf453_id) || anyDuplicated(rownames(m)) ||
        !identical(rownames(m), as.character(m$orig.ident)) ||
        !identical(rownames(m), object$cells)) stop("Inconsistent RNA profile metadata.")
    m$profile_id <- rownames(m)
    m$mouse_id <- as.character(m$cf453_id)
    ix <- match(m$mouse_id, metadata$mouse_id)
    if (anyNA(ix) || !identical(as.character(m$group), metadata$group[ix])) {
      stop("RNA and lipid mouse identities or sex/diet groups disagree.")
    }
    m$sex <- metadata$sex[ix]
    m$diet <- metadata$diet[ix]
    if (kind == "hepatocytes") {
      if (nrow(m) != 28L || anyDuplicated(m$mouse_id) ||
          !setequal(m$mouse_id, metadata$mouse_id)) stop("Expected one hepatocyte profile per mouse.")
      m$zonation_type <- "Hepatocytes"
    } else {
      if (!"zonation_type" %in% names(m) || nrow(m) != 56L || anyNA(m$zonation_type) ||
          !all(m$zonation_type %in% c("Periportal", "Pericentral")) ||
          anyDuplicated(paste(m$mouse_id, m$zonation_type))) stop("Expected two distinct zonation profiles per mouse.")
      counts_by_mouse <- table(factor(m$mouse_id, levels = metadata$mouse_id),
                               factor(m$zonation_type, levels = c("Periportal", "Pericentral")))
      if (any(counts_by_mouse != 1L)) stop("Missing or repeated mouse-zone profile.")
    }
    expected_profiles <- if (kind == "hepatocytes") paste(m$mouse_id, m$group, sep = "_") else
      paste(m$mouse_id, m$zonation_type, m$group, sep = "_")
    identities <- object$identities
    if (!identical(m$profile_id, expected_profiles) ||
        !identical(names(identities), m$profile_id) ||
        !identical(as.character(identities), m$profile_id)) {
      stop("RNA cf453_id, profile labels or active identity names disagree.")
    }
    layers <- lapply(expected_layers, function(layer) {
      value <- object$layers[[layer]]
      if (!is.numeric(value) && !inherits(value, "Matrix")) stop("Non-numeric RNA layer.")
      if (is.null(rownames(value)) || is.null(colnames(value)) || anyDuplicated(rownames(value)) ||
          anyDuplicated(colnames(value)) || anyNA(rownames(value)) ||
          any(!nzchar(rownames(value))) || any(!is.finite(value)) ||
          !setequal(colnames(value), m$profile_id)) stop("Invalid RNA gene/profile names or values.")
      value <- value[, match(m$profile_id, colnames(value)), drop = FALSE]
      if (layer != "scale.data" && any(value < 0)) stop("Negative values in RNA counts/data layer.")
      if (layer == "counts" && any(value != round(value))) stop("Non-integer RNA counts.")
      value
    })
    names(layers) <- expected_layers
    if (!identical(rownames(layers$counts), rownames(layers$data)) ||
        !identical(rownames(layers$data), rownames(layers$scale.data))) stop("RNA layer gene rosters differ.")
    qc <- data.frame(m[, c("profile_id", "mouse_id", "group", "sex", "diet", "zonation_type")],
                     object_kind = kind,
                     supplied_count_library_size = as.numeric(Matrix::colSums(layers$counts)),
                     supplied_count_detected_genes = as.numeric(Matrix::colSums(layers$counts > 0)),
                     supplied_data_sum = colSums(layers$data),
                     supplied_data_max = apply(layers$data, 2L, max),
                     stringsAsFactors = FALSE, row.names = NULL)
    inventory <- do.call(rbind, lapply(expected_layers, function(layer) {
      v <- layers[[layer]]
      data.frame(object_kind = kind, layer = layer,
                 matrix_class = paste(object$layer_classes[[layer]], collapse = "/"),
                 transport_matrix_class = paste(class(v), collapse = "/"), n_genes = nrow(v),
                 n_profiles = ncol(v), minimum = min(v), maximum = max(v),
                 nonfinite_values = sum(!is.finite(v)), stringsAsFactors = FALSE)
    }))
    reconstructed_scale <- wp2_row_scale(as.matrix(layers$data))
    scale_replay_max_difference <- max(abs(reconstructed_scale - layers$scale.data))
    list(metadata = m, layers = layers, profile_qc = qc, layer_inventory = inventory,
         data_scaled = reconstructed_scale, scale_replay_max_difference = scale_replay_max_difference,
         recorded_commands = object$recorded_commands, object_version = object$object_version,
         source_paths = object$source_paths, source_identity = object$source)
  }
  hepatocytes <- read_rna(input_paths[["hepatocytes"]], "hepatocytes")
  zonation <- read_rna(input_paths[["zonation"]], "zonation")
  payloads<-setdiff(unique(c(hepatocytes$source_paths,zonation$source_paths)),unname(input_paths))
  names(payloads)<-sprintf("RNA_payload_%03d",seq_along(payloads))
  input_paths<-c(input_paths,payloads)
  if (!identical(rownames(hepatocytes$layers$data), rownames(zonation$layers$data))) {
    stop("Hepatocyte and zonation objects do not have the same ordered gene roster.")
  }
  make_zones <- function(layer) {
    result <- lapply(c("Periportal", "Pericentral"), function(zone) {
      wp2_align_rna_zone(layer, zonation$metadata, metadata$mouse_id, zone)
    })
    names(result) <- c("Periportal", "Pericentral")
    result
  }
  rna <- make_zones(zonation$layers$scale.data)
  rna_data <- make_zones(zonation$layers$data)
  rna_data_scaled <- make_zones(zonation$data_scaled)
  rna_counts <- make_zones(zonation$layers$counts)
  ix_hep <- match(metadata$mouse_id, hepatocytes$metadata$mouse_id)
  for (kind in c("hepatocyte", "periportal", "pericentral")) {
    p <- if (kind == "hepatocyte") hepatocytes$metadata else
      zonation$metadata[zonation$metadata$zonation_type == if (kind == "periportal") "Periportal" else "Pericentral", , drop = FALSE]
    metadata[[paste0(kind, "_profile")]] <- p$profile_id[match(metadata$mouse_id, p$mouse_id)]
  }
  hepa <- function(layer) {
    result <- as.matrix(layer[, ix_hep, drop = FALSE])
    colnames(result) <- metadata$mouse_id
    result
  }
  marker_genes <- c("Cyp2e1", "Cyp2f2", "Glul", "Gls2", "Hal", "Sds", "Mtmr11", "Cyp2c38")
  absent_markers <- setdiff(marker_genes, rownames(rna$Periportal))
  if (length(absent_markers)) stop("Missing zonation/focus markers: ", paste(absent_markers, collapse = ", "))
  marker_by_profile <- do.call(rbind, lapply(c("Periportal", "Pericentral"), function(zone) {
    do.call(rbind, lapply(marker_genes, function(gene) {
      data.frame(mouse_id = metadata$mouse_id, sex = metadata$sex, diet = metadata$diet,
                 group = metadata$group, zonation_type = zone, gene = gene,
                 supplied_data_expression = unname(rna_data[[zone]][gene, ]),
                 stored_scaled_expression = unname(rna[[zone]][gene, ]),
                 supplied_counts = unname(rna_counts[[zone]][gene, ]), stringsAsFactors = FALSE)
    }))
  }))
  marker_summary <- stats::aggregate(
    cbind(supplied_data_expression, stored_scaled_expression, supplied_counts) ~
      gene + sex + diet + zonation_type, data = marker_by_profile, FUN = mean)
  # The conventional nucleus-level threshold is diagnostic only here: applying
  # it to aggregate profiles cannot reconstruct how nuclei were originally gated.
  zonation_marker_diagnostic <- do.call(rbind, lapply(c("Periportal", "Pericentral"), function(zone) {
    e1 <- rna_data[[zone]]["Cyp2e1", ]; f2 <- rna_data[[zone]]["Cyp2f2", ]
    expected <- if (zone == "Pericentral") e1 > f2 else f2 > e1
    threshold <- if (zone == "Pericentral") e1 > 1 & f2 < 1 else e1 < 1 & f2 > 1
    data.frame(mouse_id = metadata$mouse_id, sex = metadata$sex, diet = metadata$diet,
               zonation_type = zone, Cyp2e1_data = unname(e1), Cyp2f2_data = unname(f2),
               expected_marker_order = unname(expected), aggregate_data_threshold_rule = unname(threshold),
               stringsAsFactors = FALSE)
  }))
  profiles <- rbind(hepatocytes$profile_qc, zonation$profile_qc)
  inventory <- rbind(hepatocytes$layer_inventory, zonation$layer_inventory)
  metric <- function(key, value, unit = "count", note = "") {
    data.frame(metric = key, value = as.character(value), unit = unit, note = note,
               stringsAsFactors = FALSE)
  }
  qc_summary <- rbind(
    metric("biological_mice", nrow(metadata)),
    metric("mice_per_sex_diet_group", 7L),
    metric("analytical_QC_injections", length(qc_columns), note = "Not biological replicates"),
    metric("supplied_lipid_feature_rows", nrow(raw)),
    metric("internal_standard_exclusions", sum(ledger$exclusion_reason == "internal_standard")),
    metric("additional_isotope_label_exclusions", sum(ledger$exclusion_reason == "isotope_labelled_feature")),
    metric("repository_class_exclusions", sum(ledger$exclusion_reason == "repository_excluded_class")),
    metric("primary_lipid_features", nrow(lipid)),
    metric("QC_CV_le30_lipid_features", sum(ledger$keep_qc30), note = "Sensitivity set; primary retains all otherwise eligible features"),
    metric("primary_features_QC_CV_gt30", sum(ledger$keep_primary & !ledger$keep_qc30)),
    metric("primary_QC_CV_median", stats::median(ledger$qc_cv_percent[ledger$keep_primary]), "percent"),
    metric("primary_QC_CV_maximum", max(ledger$qc_cv_percent[ledger$keep_primary]), "percent"),
    metric("duplicate_alignment_ID_excess_rows", sum(duplicated(ledger[["Alignment ID"]]))),
    metric("duplicate_alignment_adduct_excess_rows", sum(duplicated(paste(ledger[["Alignment ID"]], ledger[["Adduct type"]], sep = "__")))),
    metric("duplicate_original_annotation_excess_rows", sum(duplicated(names_raw))),
    metric("duplicate_canonical_annotation_excess_rows", sum(duplicated(ledger$canonical_lipid))),
    metric("pipe_delimited_annotation_rows", sum(ledger$pipe_delimited_annotation), note = "Hierarchical annotations retained; final segment used for strict matching"),
    metric("hepatocyte_profiles", nrow(hepatocytes$metadata)),
    metric("zonation_profiles", nrow(zonation$metadata)),
    metric("RNA_genes", nrow(rna$Periportal)),
    metric("hepatocyte_count_library_minimum", min(hepatocytes$profile_qc$supplied_count_library_size), "supplied counts"),
    metric("zonation_count_library_minimum", min(zonation$profile_qc$supplied_count_library_size), "supplied counts"),
    metric("RNA_data_scale_replay_max_difference_hepatocyte", hepatocytes$scale_replay_max_difference, "absolute difference"),
    metric("RNA_data_scale_replay_max_difference_zonation", zonation$scale_replay_max_difference, "absolute difference"),
    metric("profiles_expected_marker_order", sum(zonation_marker_diagnostic$expected_marker_order)),
    metric("profiles_aggregate_marker_threshold", sum(zonation_marker_diagnostic$aggregate_data_threshold_rule), note = "Aggregate diagnostic cannot validate the original nucleus-level gate"),
    metric("mice_with_missing_nucleus_counts", nrow(metadata), note = "Nucleus contribution is unavailable in supplied metadata; count library size is not a nucleus count")
  )
  limitations <- c(
    "The supplied files are processed pseudobulks; per-mouse and per-zone nucleus counts, doublet/QC summaries and composition metadata are absent.",
    "Only ScaleData.RNA is recorded in each command history. The upstream normalization, missing-value treatment, imputation and nuclear zonation gate cannot be established from these objects.",
    "Row scaling the stored data layer reproduces scale.data numerically; it supplies a transparent scaling sensitivity but no independent normalization validation.",
    "Lipid intensity positivity and QC CV describe the supplied export, without establishing raw-data detection, batch correction or normalization provenance.",
    "No animal or lipid is removed on the basis of a network association, apparent biological outlier or favorable result. The QC-CV <=30% set is a sensitivity only.",
    "Strict canonical labels preserve chain order, separators, ether prefixes and oxygenation. Formatting agreement is not independent chemical identification; ambiguous annotations must not be forced into matches.",
    "Each mouse contributes two RNA zones but one whole-liver lipid profile. Twenty-eight mice remain the biological sample count."
  )
  cat("WDSD intake: ", nrow(metadata), " paired mice; ", nrow(lipid), " primary lipid features; ",
      nrow(rna$Periportal), " genes x 28 mice in each zone.\n", sep = "")
  cat("Supplied RNA metadata: ", nrow(hepatocytes$metadata), " hepatocyte and ", nrow(zonation$metadata),
      " zonation profiles; QC-CV <=30%: ", sum(ledger$keep_qc30), " lipid features.\n", sep = "")
  print(metadata[, c("mouse_id", "sex", "diet", "group")][seq_len(3L), , drop = FALSE])
  list(metadata = metadata, lipid = lipid, lipid_supplied = lipid_supplied,
       lipid_ledger = ledger, qc_summary = qc_summary, analytical_qc = analytical_qc,
       rna = rna, rna_data = rna_data, rna_data_scaled = rna_data_scaled, rna_counts = rna_counts,
       rna_hepatocytes = hepa(hepatocytes$layers$scale.data),
       rna_hepatocytes_data = hepa(hepatocytes$layers$data),
       rna_hepatocytes_data_scaled = hepa(hepatocytes$data_scaled),
       rna_profile_metadata = rbind(hepatocytes$metadata, zonation$metadata),
       rna_profile_qc = profiles, layer_inventory = inventory,
       marker_by_profile = marker_by_profile, marker_summary = marker_summary,
       zonation_marker_diagnostic = zonation_marker_diagnostic,
       recorded_commands = list(hepatocytes = hepatocytes$recorded_commands,
                                zonation = zonation$recorded_commands),
       object_versions = c(hepatocytes = hepatocytes$object_version, zonation = zonation$object_version),
       excluded_lipid_classes = excluded_classes, input_paths = input_paths,
       input_source_identity = list(hepatocytes=hepatocytes$source_identity,zonation=zonation$source_identity),
       limitations = limitations,
       preprocessing = list(lipid_primary = "log2_supplied_positive_intensities",
                            RNA_primary = "stored_scale.data_layer",
                            RNA_sensitivity = "row_scaled_stored_data_layer_full_original_profile_roster",
                            canonical_lipid = "wp1_canonical_lipid; final pipe segment; class formatting only"))
}
