# WP3 CF459 intake: processed WD-only cohort, not a within-cohort diet contrast.
# Source scripts/support_functions.R first. Existing layers are read unchanged.

wd_lipid_sample_metadata <- function(columns) {
  samples <- grep("^cf459-[0-9]+$", columns, value = TRUE)
  expected <- paste0("cf459-", c(1L:4L, 13L:16L))
  if (length(samples) != 8L || anyDuplicated(samples) ||
      !setequal(samples, expected)) {
    stop("Expected eight unique CF459 mouse IDs: 1-4 and 13-16.")
  }
  data.frame(mouse_id = expected, sex = rep(c("female", "male"), each = 4L),
             diet = "WD", group = rep(c("female-wd", "male-wd"), each = 4L),
             original_rna_group = rep(c("f-ctrl", "m-ctrl"), each = 4L),
             lipid_sample = samples[match(expected, samples)],
             stringsAsFactors = FALSE)
}

# Validate explicit mouse identities and the supplied control-label sex map.
# The WD assignment comes from the meeting, not from the ambiguous word ctrl.
wd_validate_rna_metadata <- function(m, object_columns, lipid_metadata, kind) {
  if (!kind %in% c("cell_types", "zonation")) stop("Unknown RNA object kind.")
  required <- c("orig.ident", "sample_name", "group",
                if (kind == "cell_types") "cell_type" else "zonation_type")
  if (!all(required %in% names(m)) || anyNA(m[, required, drop = FALSE]) ||
      anyDuplicated(rownames(m)) ||
      !identical(rownames(m), as.character(m$orig.ident)) ||
      !identical(rownames(m), object_columns)) {
    stop("Inconsistent CF459 RNA profile metadata or object columns.")
  }
  m$profile_id <- rownames(m)
  m$mouse_id <- as.character(m$sample_name)
  ix <- match(m$mouse_id, lipid_metadata$mouse_id)
  if (anyNA(ix) ||
      !identical(as.character(m$group), lipid_metadata$original_rna_group[ix])) {
    stop("CF459 RNA and lipid mouse identities or sex assignments disagree.")
  }
  m$original_group <- as.character(m$group)
  m$group <- lipid_metadata$group[ix]
  m$sex <- lipid_metadata$sex[ix]
  m$diet <- lipid_metadata$diet[ix]
  levels <- if (kind == "cell_types")
    c("Hepatocytes", "Endothelial cells", "Immune cells", "Stellate cells", "Cholangiocytes") else
      c("Periportal", "Pericentral")
  label <- as.character(m[[if (kind == "cell_types") "cell_type" else "zonation_type"]])
  if (nrow(m) != 8L * length(levels) || !all(label %in% levels) ||
      anyDuplicated(paste(m$mouse_id, label))) {
    stop("Missing or repeated CF459 mouse-cell-type/zone profile.")
  }
  roster <- table(factor(m$mouse_id, levels = lipid_metadata$mouse_id),
                  factor(label, levels = levels))
  if (any(roster != 1L)) stop("Each mouse must contribute exactly one profile per cell type/zone.")
  expected_profiles <- paste(m$mouse_id, label, m$original_group, sep = "_")
  if (!identical(m$profile_id, expected_profiles)) {
    stop("RNA sample_name, cell type/zone and profile labels disagree.")
  }
  if (kind == "cell_types") m$zonation_type <- as.character(m$cell_type) else m$cell_type <- "Hepatocytes"
  m
}

wd_load_inputs <- function(root = ".",input_cell_types=NULL,input_zonation=NULL) {
  if (length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop("root must identify one existing analysis repository directory.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  source(file.path(root,"scripts/config/analysis_release.R"),local=TRUE)
  source(file.path(root,"scripts/core/pseudobulk_csv.R"),local=TRUE)
  for (package in "Matrix") {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  }
  input_paths <- c(lipidomics = file.path(root, "data/metabolomics/wd/lipidomics_wd.csv"),
                   cell_types = analysis_path("wd_cell_types_transport",root),
                   zonation = analysis_path("wd_zonation_transport",root))
  if(!is.null(input_cell_types))input_paths[["cell_types"]]<-input_cell_types
  if(!is.null(input_zonation))input_paths[["zonation"]]<-input_zonation
  absent <- input_paths[!file.exists(input_paths)]
  if (length(absent)) stop("Missing CF459 input(s): ", paste(absent, collapse = ", "))
  input_paths <- vapply(input_paths, normalizePath, character(1L), mustWork = TRUE)
  helpers <- new.env(parent = baseenv())
  sys.source(file.path(root, "scripts/core/cyp2c38_reference.R"), envir = helpers)
  sys.source(file.path(root, "scripts/core/wdsd_intake.R"), envir = helpers)

  # read.csv2 is intentional: semicolons, decimal commas, scientific notation,
  # quoted annotation semicolons and a UTF-8 BOM all occur in the original file.
  raw <- utils::read.csv2(input_paths[["lipidomics"]], check.names = FALSE,
                          stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  if (anyDuplicated(names(raw)) ||
      !all(c("Alignment ID", "Metabolite name") %in% names(raw))) {
    stop("CF459 lipidomics requires unique annotation/sample column names.")
  }
  metadata <- wd_lipid_sample_metadata(names(raw))
  extra <- setdiff(names(raw), c("Alignment ID", "Metabolite name", metadata$lipid_sample))
  if (length(extra)) stop("Unexpected CF459 columns; review changed export: ", paste(extra, collapse = ", "))
  x_all <- as.matrix(raw[, metadata$lipid_sample, drop = FALSE])
  if (!is.numeric(x_all) || any(!is.finite(x_all)) || any(x_all <= 0)) {
    stop("Supplied CF459 intensities must parse as numeric, finite and positive.")
  }
  annotation <- raw[["Metabolite name"]]
  if (anyNA(annotation) || any(!nzchar(trimws(annotation))) || anyNA(raw[["Alignment ID"]])) {
    stop("Empty CF459 lipid annotation or alignment ID.")
  }
  ledger <- raw[, c("Alignment ID", "Metabolite name"), drop = FALSE]
  ledger$feature_id <- sprintf("WD_%04d", seq_len(nrow(raw)))
  ledger$canonical_lipid <- helpers$wp1_canonical_lipid(annotation)
  ledger$lipid_class <- sub("\\(.*$", "", ledger$canonical_lipid)
  ledger$internal_standard <- grepl("^IS:", annotation)
  ledger$isotope_labelled <- grepl("\\([dD][0-9]+\\)|\\b[0-9]+[Cc]\\b", annotation)
  excluded_classes <- c("BA", "DGGA", "HBMP", "LPE-N", "NATau", "PEtOH", "PMeOH", "SL", "CAR", "SPB", "VAE")
  ledger$excluded_class <- ledger$lipid_class %in% excluded_classes
  ledger$exclusion_reason <- ifelse(ledger$internal_standard, "internal_standard",
    ifelse(ledger$isotope_labelled, "isotope_labelled_feature",
      ifelse(ledger$excluded_class, "repository_excluded_class", "")))
  ledger$keep_primary <- ledger$exclusion_reason == ""
  multiplicity <- function(x) as.integer(table(x)[as.character(x)])
  ledger$alignment_id_multiplicity <- multiplicity(ledger[["Alignment ID"]])
  ledger$original_annotation_multiplicity <- multiplicity(annotation)
  ledger$canonical_annotation_multiplicity <- multiplicity(ledger$canonical_lipid)
  ledger$pipe_delimited_annotation <- grepl("|", annotation, fixed = TRUE)
  ledger$annotation_join_status <- ifelse(ledger$canonical_annotation_multiplicity == 1L,
    "unique_canonical_annotation", "ambiguous_canonical_annotation")
  ledger$identity_policy <- "last_pipe_annotation_and_class_formatting_only"
  # Missing analytical controls and annotation fields are explicitly unavailable.
  ledger$qc_cv_percent <- NA_real_
  ledger$keep_qc30 <- NA
  ledger$analytical_qc_available <- FALSE
  ledger$adduct_available <- FALSE
  ledger$fill_percentage_available <- FALSE
  if (anyDuplicated(ledger$feature_id) ||
      anyDuplicated(ledger$canonical_lipid[ledger$keep_primary])) {
    stop("Retained CF459 canonical annotations must be unique; ambiguous joins are prohibited.")
  }
  lipid_supplied <- x_all[ledger$keep_primary, , drop = FALSE]
  dimnames(lipid_supplied) <- list(ledger$feature_id[ledger$keep_primary], metadata$mouse_id)
  lipid <- log2(lipid_supplied)
  metadata$total_supplied_lipid_intensity <- colSums(x_all)
  metadata$total_primary_lipid_intensity <- colSums(lipid_supplied)
  metadata$median_primary_log2_intensity <- apply(lipid, 2L, stats::median)
  metadata$nuclei_contributing <- NA_integer_

  read_rna <- function(path, kind) {
    object <- pb_components(path)
    expected_layers <- c("counts", "data", "scale.data")
    if (!all(expected_layers %in% names(object$layers))) stop("Missing CF459 RNA layer.")
    m <- wd_validate_rna_metadata(object$meta.data, object$cells, metadata, kind)
    identities <- object$identities
    if (!identical(names(identities), m$profile_id) ||
        !identical(as.character(identities), m$profile_id)) stop("RNA active identities disagree with profile names.")
    layers <- lapply(expected_layers, function(layer) {
      v <- object$layers[[layer]]
      if ((!is.numeric(v) && !inherits(v, "Matrix")) ||
          is.null(rownames(v)) || is.null(colnames(v)) ||
          anyNA(rownames(v)) || any(!nzchar(rownames(v))) ||
          anyDuplicated(rownames(v)) || anyDuplicated(colnames(v)) ||
          any(!is.finite(v)) || !setequal(colnames(v), m$profile_id)) {
        stop("Invalid CF459 RNA layer names, types or values.")
      }
      v <- v[, match(m$profile_id, colnames(v)), drop = FALSE]
      if (layer != "scale.data" && any(v < 0)) stop("Negative RNA counts/data.")
      if (layer == "counts" && any(v != round(v))) stop("Non-integer RNA counts.")
      v
    })
    names(layers) <- expected_layers
    if (!identical(rownames(layers$counts), rownames(layers$data)) ||
        !identical(rownames(layers$data), rownames(layers$scale.data))) stop("RNA layer gene rosters differ.")
    qc <- data.frame(m[, c("profile_id", "mouse_id", "group", "original_group", "sex", "diet", "cell_type", "zonation_type")],
      object_kind = kind, supplied_count_library_size = as.numeric(Matrix::colSums(layers$counts)),
      supplied_count_detected_genes = as.numeric(Matrix::colSums(layers$counts > 0)),
      supplied_data_sum = colSums(layers$data), supplied_data_max = apply(layers$data, 2L, max),
      nuclei_contributing = NA_integer_, stringsAsFactors = FALSE, row.names = NULL)
    inventory <- do.call(rbind, lapply(expected_layers, function(layer) {
      v <- layers[[layer]]
      data.frame(object_kind = kind, layer = layer, matrix_class = paste(object$layer_classes[[layer]], collapse = "/"),
                 transport_matrix_class = paste(class(v), collapse = "/"),
        n_genes = nrow(v), n_profiles = ncol(v), minimum = min(v), maximum = max(v),
        nonfinite_values = sum(!is.finite(v)), stringsAsFactors = FALSE)
    }))
    reconstructed <- helpers$wp2_row_scale(as.matrix(layers$data))
    scale_replay_difference <- max(abs(reconstructed - layers$scale.data))
    list(metadata = m, layers = layers, profile_qc = qc, layer_inventory = inventory,
      data_scaled = reconstructed, scale_replay_max_difference = scale_replay_difference,
      recorded_commands = object$recorded_commands, object_version = object$object_version,
      source_paths = object$source_paths, source_identity = object$source)
  }
  cell_types <- read_rna(input_paths[["cell_types"]], "cell_types")
  zonation <- read_rna(input_paths[["zonation"]], "zonation")
  payloads<-setdiff(unique(c(cell_types$source_paths,zonation$source_paths)),unname(input_paths))
  names(payloads)<-sprintf("RNA_payload_%03d",seq_along(payloads))
  input_paths<-c(input_paths,payloads)
  if (!identical(rownames(cell_types$layers$data), rownames(zonation$layers$data))) stop("Cell-type and zonation gene rosters differ.")
  zones <- c("Periportal", "Pericentral")
  align_zones <- function(layer) {
    result <- lapply(zones, function(zone) helpers$wp2_align_rna_zone(layer, zonation$metadata, metadata$mouse_id, zone))
    names(result) <- zones
    result
  }
  rna <- align_zones(zonation$layers$scale.data)
  rna_data <- align_zones(zonation$layers$data)
  rna_data_scaled <- align_zones(zonation$data_scaled)
  rna_counts <- align_zones(zonation$layers$counts)
  hepa_metadata <- cell_types$metadata[cell_types$metadata$cell_type == "Hepatocytes", , drop = FALSE]
  hepa_profiles <- hepa_metadata$profile_id[match(metadata$mouse_id, hepa_metadata$mouse_id)]
  hepa <- function(layer) {
    value <- as.matrix(layer[, match(hepa_profiles, colnames(layer)), drop = FALSE])
    colnames(value) <- metadata$mouse_id
    value
  }
  metadata$hepatocyte_profile <- hepa_profiles
  for (zone in zones) {
    m <- zonation$metadata[zonation$metadata$zonation_type == zone, , drop = FALSE]
    metadata[[paste0(tolower(zone), "_profile")]] <- m$profile_id[match(metadata$mouse_id, m$mouse_id)]
  }
  marker_genes <- c("Cyp2e1", "Cyp2f2", "Glul", "Gls2", "Hal", "Sds", "Mtmr11", "Cyp2c38")
  if (length(setdiff(marker_genes, rownames(rna$Periportal)))) stop("Missing zonation/focus marker genes.")
  marker_by_profile <- do.call(rbind, lapply(zones, function(zone) {
    do.call(rbind, lapply(marker_genes, function(gene) {
      data.frame(mouse_id = metadata$mouse_id, sex = metadata$sex, diet = metadata$diet,
        group = metadata$group, zonation_type = zone, gene = gene,
        supplied_data_expression = unname(rna_data[[zone]][gene, ]),
        stored_scaled_expression = unname(rna[[zone]][gene, ]),
        supplied_counts = unname(rna_counts[[zone]][gene, ]), stringsAsFactors = FALSE)
    }))
  }))
  marker_summary <- stats::aggregate(cbind(supplied_data_expression, stored_scaled_expression, supplied_counts) ~
    gene + sex + diet + zonation_type, data = marker_by_profile, FUN = mean)
  zonation_marker_diagnostic <- do.call(rbind, lapply(zones, function(zone) {
    e1 <- rna_data[[zone]]["Cyp2e1", ]; f2 <- rna_data[[zone]]["Cyp2f2", ]
    expected <- if (zone == "Pericentral") e1 > f2 else f2 > e1
    threshold <- if (zone == "Pericentral") e1 > 1 & f2 < 1 else e1 < 1 & f2 > 1
    data.frame(mouse_id = metadata$mouse_id, sex = metadata$sex, diet = metadata$diet,
      zonation_type = zone, Cyp2e1_data = unname(e1), Cyp2f2_data = unname(f2),
      expected_marker_order = unname(expected), aggregate_data_threshold_rule = unname(threshold),
      stringsAsFactors = FALSE)
  }))
  metric <- function(key, value, unit = "count", note = "") {
    data.frame(metric = key, value = as.character(value), unit = unit, note = note, stringsAsFactors = FALSE)
  }
  qc_summary <- rbind(metric("biological_mice", nrow(metadata)), metric("male_WD_mice_primary", 4L),
    metric("female_WD_mice_intake_only", 4L), metric("analytical_QC_injections", 0L, note = "Unavailable in supplied export; not zero analytical variability"),
    metric("supplied_lipid_feature_rows", nrow(raw)), metric("primary_lipid_features", nrow(lipid)),
    metric("internal_standard_exclusions", sum(ledger$internal_standard)),
    metric("additional_isotope_label_exclusions", sum(ledger$exclusion_reason == "isotope_labelled_feature")),
    metric("repository_class_exclusions", sum(ledger$exclusion_reason == "repository_excluded_class")),
    metric("duplicate_alignment_ID_excess_rows", sum(duplicated(ledger[["Alignment ID"]]))),
    metric("duplicate_original_annotation_excess_rows", sum(duplicated(annotation))),
    metric("duplicate_canonical_annotation_excess_rows", sum(duplicated(ledger$canonical_lipid))),
    metric("pipe_delimited_annotation_rows", sum(ledger$pipe_delimited_annotation)),
    metric("cell_type_profiles", nrow(cell_types$metadata)), metric("hepatocyte_profiles", nrow(hepa_metadata)),
    metric("zonation_profiles", nrow(zonation$metadata)), metric("RNA_genes", nrow(rna$Periportal)),
    metric("RNA_data_scale_replay_max_difference_cell_types", cell_types$scale_replay_max_difference, "absolute difference", "Across all 40 original cell-type profiles, before hepatocyte subset"),
    metric("RNA_data_scale_replay_max_difference_zonation", zonation$scale_replay_max_difference, "absolute difference", "Across all 16 original zone profiles, before male subset"),
    metric("profiles_expected_marker_order", sum(zonation_marker_diagnostic$expected_marker_order)),
    metric("profiles_aggregate_marker_threshold", sum(zonation_marker_diagnostic$aggregate_data_threshold_rule), note = "Aggregate diagnostic cannot validate original nucleus-level gating"),
    metric("mice_with_missing_nucleus_counts", nrow(metadata), note = "Count-library size is not a nucleus count"))
  limitations <- c(
    "CF459 is WD-only according to the meeting; f-ctrl/m-ctrl metadata labels do not identify standard-diet controls. There is no within-CF459 diet contrast.",
    "Four CF459 male mice are the primary cross-cohort group; the four females are retained only for input/sex-map QC and original scaling-roster diagnostics.",
    "Exact ages, genetic background, diet formulation/exposure, litter and batch details remain unavailable; the meeting described qualitative cohort differences only.",
    "No analytical QC columns, adduct, retention time/mass-to-charge or Fill percentage are supplied for CF459; QC-CV and annotation-field sensitivities cannot be inferred.",
    "Only ScaleData.RNA is recorded. Upstream normalization, raw detection, missing-value treatment/imputation and nucleus-level QC/gating cannot be reconstructed.",
    "RNA primary uses existing scale.data unchanged. Row scaling stored data across each complete original object is diagnostic and does not independently validate normalization.",
    "Nucleus contributions, composition and doublet/QC summaries are unavailable. Supplied count totals and detected genes are aggregate profile diagnostics only.",
    "Strict canonical labels preserve chain order, separators, ether prefixes and oxygenation. Formatting agreement does not independently confirm molecular identity; alignment IDs are not cross-cohort identifiers.",
    "Each mouse contributes two RNA zones and one bulk liver lipid profile; zone profiles are not independent animals. No outcome-based animal/feature exclusions are applied.")
  cat("CF459 intake: ", nrow(metadata), " paired WD mice (4 male primary); ", nrow(lipid),
      " eligible lipid features; ", nrow(rna$Periportal), " genes x 8 mice per zone.\n", sep = "")
  print(metadata[, c("mouse_id", "sex", "diet", "original_rna_group")])
  list(metadata = metadata, lipid = lipid, lipid_supplied = lipid_supplied, lipid_ledger = ledger,
    qc_summary = qc_summary, analytical_qc = data.frame(qc_injection = character(), total_supplied_intensity = numeric(), total_primary_intensity = numeric()),
    rna = rna, rna_data = rna_data, rna_data_scaled = rna_data_scaled, rna_counts = rna_counts,
    rna_hepatocytes = hepa(cell_types$layers$scale.data), rna_hepatocytes_data = hepa(cell_types$layers$data),
    rna_hepatocytes_data_scaled = hepa(cell_types$data_scaled), rna_hepatocytes_counts = hepa(cell_types$layers$counts),
    rna_profile_metadata = rbind(cell_types$metadata, zonation$metadata[, names(cell_types$metadata), drop = FALSE]),
    rna_profile_qc = rbind(cell_types$profile_qc, zonation$profile_qc),
    layer_inventory = rbind(cell_types$layer_inventory, zonation$layer_inventory),
    marker_by_profile = marker_by_profile, marker_summary = marker_summary,
    zonation_marker_diagnostic = zonation_marker_diagnostic,
    recorded_commands = list(cell_types = cell_types$recorded_commands, zonation = zonation$recorded_commands),
    object_versions = c(cell_types = cell_types$object_version, zonation = zonation$object_version),
    excluded_lipid_classes = excluded_classes, input_paths = input_paths, limitations = limitations,
    input_source_identity=list(cell_types=cell_types$source_identity,zonation=zonation$source_identity),
    preprocessing = list(lipid_primary = "log2_supplied_positive_intensities", RNA_primary = "stored_scale.data_layer",
      RNA_sensitivity = "row_scaled_stored_data_layer_full_original_profile_roster",
      canonical_lipid = "wp1_canonical_lipid; final pipe segment; class formatting only",
      diet_assignment = "WD_only_from_2026_10_05_meeting; ctrl_is_not_SD"))
}
