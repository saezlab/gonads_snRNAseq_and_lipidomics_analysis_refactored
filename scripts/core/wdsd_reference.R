source("scripts/config/analysis_release.R")
# Frozen original Figure 5 FCG gene-lipid references for WDSD replication.
# Source scripts/support_functions.R before using this module. Sourcing this file
# defines helpers only and does not write or overwrite archived result files.

#' Replay four archived context-specific FCG gene-lipid networks.
#'
#' The returned list has four named contexts, each with ordered metadata, a
#' gene-by-mouse scaled RNA matrix, a lipid-by-the-same-mouse log2 abundance
#' matrix, and archived gene-lipid edges. primary_edge is exact membership in
#' the original standard SIF, never a newly applied cutoff in the new cohort.
#'
#' All lipid annotations retain their original structural resolution. Branch
#' identity is determined by an original measured lipid target and a measured
#' RNA source; an ATT label of TF alone does not exclude a gene-lipid edge when
#' that RNA gene is also present in the TF regulon. TF-to-gene edges are excluded.
#' Both raw and log2 lipid coefficients are audited, and only an exact log2
#' archive replay is accepted. No silent selection or transformation changes.
#'
#' @param root Existing analysis repository root.
#' @return Four context lists; attributes summary and provenance record counts,
#'   original selection details, input paths, and transformation replay checks.
wp2_load_reference <- function(root = ".") {
  if (length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop("root must identify one existing analysis repository directory.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  definition <- data.frame(
    context = c("male_periportal", "male_pericentral",
                "female_periportal", "female_pericentral"),
    archive_stem = c("testes_portal", "testes_central",
                     "Ovaries_portal", "ova_central"),
    gonadal_type = c("Testes", "Testes", "Ovaries", "Ovaries"),
    zone = c("Periportal", "Pericentral", "Periportal", "Pericentral"),
    expected_mice = c(22L, 22L, 21L, 21L),
    expected_genes = c(25L, 30L, 12L, 11L),
    expected_primary_edges = c(155L, 387L, 89L, 77L),
    stringsAsFactors = FALSE
  )
  rna_file <- analysis_path("historical_fcg_zonation_csv",root)
  metadata_file <- analysis_path("historical_fcg_zonation_metadata",root)
  lipid_file <- file.path(root, "data/metabolomics/Lipid_data_file_JPC_cleaned.csv")
  genotype_file <- file.path(root, "results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv")
  network_dir <- file.path(analysis_path("historical_fcg05",root),"networks")
  network_paths <- unlist(lapply(definition$archive_stem, function(stem) {
    file.path(network_dir, paste0("cormat_", stem,
                                 c("_full_SIF.csv", "_SIF.csv", "_full_ATT.csv")))
  }), use.names = FALSE)
  required_paths <- c(rna_file, metadata_file, lipid_file, genotype_file, network_paths)
  if (any(!file.exists(required_paths))) {
    stop("Missing Figure 5 reference input(s): ",
         paste(required_paths[!file.exists(required_paths)], collapse = ", "))
  }
  read_table <- function(path) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  require_columns <- function(x, columns, label) {
    if (!all(columns %in% names(x))) {
      stop(label, " requires columns: ", paste(columns, collapse = ", "))
    }
  }
  assert_unique <- function(x, label) {
    if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
      stop(label, " must be unique, nonempty and nonmissing.")
    }
  }
  metadata <- read_table(metadata_file)
  require_columns(metadata, c("sample", "sample_name", "zonation_type", "gonadal_type"),
                  "Zonation metadata")
  assert_unique(metadata$sample, "Zonation RNA sample names")
  if (nrow(metadata) != 86L ||
      !all(metadata$zonation_type %in% c("Periportal", "Pericentral")) ||
      !all(metadata$gonadal_type %in% c("Testes", "Ovaries"))) {
    stop("Expected 86 zonation pseudobulks from 43 FCG mice with two zones each.")
  }
  if (!identical(metadata$sample,
                 paste(metadata$sample_name, metadata$zonation_type,
                       metadata$gonadal_type, sep = "_"))) {
    stop("Zonation metadata sample suffixes disagree with mouse, zone or gonad.")
  }
  mouse_zone_table <- table(metadata$sample_name, metadata$zonation_type)
  if (nrow(mouse_zone_table) != 43L || any(mouse_zone_table != 1L)) {
    stop("Every original mouse must contribute exactly one pseudobulk per zone.")
  }
  rna_table <- read_table(rna_file)
  if (!identical(names(rna_table)[1L], "gene")) stop("Zonation RNA must have a first gene column.")
  assert_unique(rna_table$gene, "Zonation RNA gene symbols")
  assert_unique(names(rna_table)[-1L], "Zonation RNA matrix sample names")
  if (!setequal(names(rna_table)[-1L], metadata$sample)) {
    stop("Zonation RNA and metadata sample rosters differ.")
  }
  rna <- as.matrix(rna_table[, -1L, drop = FALSE])
  rownames(rna) <- rna_table$gene
  rm(rna_table)
  if (!is.numeric(rna) || any(!is.finite(rna))) {
    stop("Original scaled zonation RNA must be finite numeric measurements.")
  }
  lipid_table <- read_table(lipid_file)
  if (!identical(names(lipid_table)[1L], "Unnamed: 0") || ncol(lipid_table) != 1255L) {
    stop("Expected an original mouse-ID column and 1,254 measured lipidomics features.")
  }
  assert_unique(names(lipid_table)[-1L], "Original lipid annotations")
  lipid_mouse <- gsub("-", "_", gsub("-0", "-", lipid_table[[1L]], fixed = TRUE))
  assert_unique(lipid_mouse, "Normalized original lipid mouse IDs")
  lipid_raw <- t(as.matrix(lipid_table[, -1L, drop = FALSE]))
  dimnames(lipid_raw) <- list(names(lipid_table)[-1L], lipid_mouse)
  rm(lipid_table)
  if (!is.numeric(lipid_raw) || any(!is.finite(lipid_raw)) || any(lipid_raw <= 0)) {
    stop("Original lipid intensities must be finite positive numeric measurements.")
  }
  genotype <- read_table(genotype_file)
  require_columns(genotype, c("sequencing_name", "group", "gonadal_type", "allosome_type"),
                  "Original FCG genotype metadata")
  genotype$mouse_id <- gsub("-", "_", genotype$sequencing_name)
  assert_unique(genotype$mouse_id, "Original genotype mouse IDs")
  if (nrow(genotype) != 43L || !setequal(genotype$mouse_id, metadata$sample_name)) {
    stop("Original zonation and recorded FCG genotype mouse rosters differ.")
  }
  contexts <- stats::setNames(vector("list", nrow(definition)), definition$context)
  summaries <- vector("list", nrow(definition))
  for (i in seq_len(nrow(definition))) {
    def <- definition[i, , drop = FALSE]
    stem <- def$archive_stem
    full_path <- file.path(network_dir, paste0("cormat_", stem, "_full_SIF.csv"))
    primary_path <- file.path(network_dir, paste0("cormat_", stem, "_SIF.csv"))
    att_path <- file.path(network_dir, paste0("cormat_", stem, "_full_ATT.csv"))
    full <- read_table(full_path)
    primary <- read_table(primary_path)
    att <- read_table(att_path)
    require_columns(full, c("node1", "node2", "weight"), "Full archived SIF")
    require_columns(primary, c("node1", "node2", "weight"), "Primary archived SIF")
    require_columns(att, c("gene", "type"), "Archived ATT")
    assert_unique(att$gene, paste(def$context, "ATT node identities"))
    if (anyNA(att$type) || !all(att$type %in% c("gene", "lipid", "TF"))) {
      stop("Unknown archived ATT node type in ", def$context)
    }
    lipid_branch <- full$node2 %in% rownames(lipid_raw)
    edges <- full[lipid_branch, , drop = FALSE]
    genes <- unique(edges$node1)
    lipid_ids <- unique(edges$node2)
    if (anyNA(full$node1) || anyNA(full$node2) ||
        !is.numeric(full$weight) || any(!is.finite(full$weight)) ||
        length(genes) != def$expected_genes ||
        !all(genes %in% rownames(rna)) ||
        !setequal(lipid_ids, rownames(lipid_raw))) {
      stop("Archived gene-lipid universe or source RNA identities changed in ", def$context)
    }
    edge_key <- paste(edges$node1, edges$node2, sep = "\r")
    assert_unique(edge_key, paste(def$context, "archived gene-lipid edge identities"))
    if (nrow(edges) != length(genes) * length(lipid_ids)) {
      stop("Full archived gene-lipid network must retain every selected gene-by-lipid pair.")
    }
    # ATT TF precedence can label an RNA gene as a TF. The measured lipid target
    # identifies the gene-lipid branch independently of that display annotation.
    gene_types <- att$type[match(genes, att$gene)]
    lipid_types <- att$type[match(lipid_ids, att$gene)]
    excluded <- full[!lipid_branch, , drop = FALSE]
    excluded_source_type <- att$type[match(excluded$node1, att$gene)]
    if (anyNA(gene_types) || !all(gene_types %in% c("gene", "TF")) ||
        anyNA(lipid_types) || !all(lipid_types == "lipid") ||
        anyNA(excluded_source_type) || !all(excluded_source_type == "TF") ||
        !all(excluded$node2 %in% genes)) {
      stop("Archived ATT and measured-feature branch audit failed in ", def$context)
    }
    primary <- primary[primary$node2 %in% rownames(lipid_raw), , drop = FALSE]
    primary_key <- paste(primary$node1, primary$node2, sep = "\r")
    assert_unique(primary_key, paste(def$context, "displayed edge identities"))
    primary_index <- match(primary_key, edge_key)
    if (nrow(primary) != def$expected_primary_edges || anyNA(primary_index) ||
        max(abs(primary$weight - edges$weight[primary_index])) > 1e-10) {
      stop("Displayed archived edge membership or weights changed in ", def$context)
    }
    primary_mask <- edge_key %in% primary_key
    if (!identical(primary_mask, abs(edges$weight) >= 0.65)) {
      stop("Archived standard network does not reproduce its documented 0.65 cutoff.")
    }
    sample_columns <- colnames(rna)[grepl(paste0("_", def$zone, "_", def$gonadal_type, "$"),
                                            colnames(rna))]
    meta <- metadata[match(sample_columns, metadata$sample), , drop = FALSE]
    mouse <- meta$sample_name
    if (length(mouse) != def$expected_mice || anyNA(mouse) || anyDuplicated(mouse) ||
        !identical(mouse, sub(paste0("_", def$zone, "_", def$gonadal_type, "$"),
                              "", sample_columns))) {
      stop("RNA column, gonad, zone and mouse pairing failed in ", def$context)
    }
    g <- genotype[match(mouse, genotype$mouse_id), , drop = FALSE]
    if (anyNA(g$mouse_id) || !identical(g$gonadal_type, meta$gonadal_type) ||
        !all(g$group %in% c("XXO", "XYO", "XXT", "XYT")) ||
        !all(g$allosome_type %in% c("XX", "XY")) ||
        any(substr(g$group, 1L, 2L) != g$allosome_type) ||
        any(substr(g$group, 3L, 3L) != ifelse(meta$gonadal_type == "Testes", "T", "O"))) {
      stop("Recorded FCG gonad, group and sex-chromosome identities disagree.")
    }
    lipid_index <- match(mouse, colnames(lipid_raw))
    if (anyNA(lipid_index)) stop("Missing paired original lipid mouse in ", def$context)
    context_rna <- rna[genes, sample_columns, drop = FALSE]
    colnames(context_rna) <- mouse
    context_raw <- lipid_raw[lipid_ids, lipid_index, drop = FALSE]
    context_lipids <- log2(context_raw)
    if (!identical(colnames(context_rna), colnames(context_lipids))) {
      stop("Original RNA-lipid mouse order does not agree exactly.")
    }
    r_log2 <- stats::cor(t(context_rna), t(context_lipids), method = "pearson")
    r_raw <- stats::cor(t(context_rna), t(context_raw), method = "pearson")
    lookup <- cbind(match(edges$node1, rownames(r_log2)),
                    match(edges$node2, colnames(r_log2)))
    log_difference <- max(abs(r_log2[lookup] - edges$weight))
    raw_difference <- max(abs(r_raw[lookup] - edges$weight))
    if (any(!is.finite(r_log2)) || !is.finite(log_difference) || log_difference > 1e-10) {
      stop("Log2 lipid intensities fail exact archived Pearson replay in ", def$context)
    }
    context_metadata <- data.frame(
      mouse_id = mouse, sample_id = mouse, rna_sample = sample_columns,
      gonadal_type = meta$gonadal_type,
      sex = ifelse(meta$gonadal_type == "Testes", "male", "female"),
      zone = meta$zonation_type, group = g$group, allosome_type = g$allosome_type,
      stringsAsFactors = FALSE
    )
    context_edges <- data.frame(
      gene = edges$node1, reference_lipid = edges$node2,
      r_reference = edges$weight, primary_edge = primary_mask,
      stringsAsFactors = FALSE
    )
    replay <- data.frame(context = def$context, n_mice = length(mouse),
                         n_genes = length(genes), n_lipids = length(lipid_ids),
                         n_full_edges = nrow(context_edges),
                         n_primary_edges = sum(primary_mask),
                         n_excluded_TF_gene_edges = nrow(excluded),
                         n_RNA_genes_labeled_TF = sum(gene_types == "TF"),
                         n_Mtmr11_full_edges = sum(context_edges$gene == "Mtmr11"),
                         n_Mtmr11_primary_edges = sum(context_edges$gene == "Mtmr11" & primary_mask),
                         log2_max_abs_archive_difference = log_difference,
                         raw_max_abs_archive_difference = raw_difference,
                         lipid_scale = "log2", stringsAsFactors = FALSE)
    contexts[[def$context]] <- list(
      metadata = context_metadata, rna = context_rna, lipids = context_lipids,
      edges = context_edges, genes = genes, lipid_ids = lipid_ids,
      replay_metrics = replay, lipid_scale = "log2",
      RNA_genes_labeled_TF = genes[gene_types == "TF"],
      source_paths = c(full_SIF = full_path, displayed_SIF = primary_path, full_ATT = att_path)
    )
    summaries[[i]] <- replay
    cat("Figure 5 ", def$context, ": ", length(mouse), " mice, ",
        length(genes), " genes, ", nrow(context_edges), " full / ",
        sum(primary_mask), " displayed gene-lipid edges; log2 replay error ",
        format(log_difference, digits = 3L), "; raw error ",
        format(raw_difference, digits = 3L), "\n", sep = "")
  }
  for (sex in c("male", "female")) {
    portal <- contexts[[paste0(sex, "_periportal")]]$metadata
    central <- contexts[[paste0(sex, "_pericentral")]]$metadata
    if (!identical(portal$mouse_id, central$mouse_id) || !identical(portal$group, central$group)) {
      stop("Original zones must use the same ordered mice within each gonadal stratum.")
    }
  }
  attr(contexts, "summary") <- do.call(rbind, summaries)
  attr(contexts, "provenance") <- list(
    inputs = required_paths,
    lipid_scale = "log2",
    scale_precedent = paste("Workflow 05 section 8 calls load_lipidomic_data();",
                            "its default is log_transform = TRUE; section 13 uses Pearson cor().",
                            "All four full archives are independently replayed against both scales."),
    gene_selection = paste("Genes frozen by original full SIF membership, originally selected",
                           "by sex-specific zonation/interaction statistics in Workflow 05 sections 11-12.",
                           "No WDSD outcomes select genes or pairs."),
    primary_selection = paste("Exact original standard SIF membership (documented |r| >= 0.65);",
                              "no new-cohort correlation threshold is applied."),
    branch_rule = paste("Retain an original measured lipid target and original measured RNA source;",
                       "exclude archived TF-to-gene edges. A gene can also carry ATT type TF."),
    interpretation = paste("Bulk liver lipid profiles paired with zonal RNA; lipids are not",
                           "measured separately in the two zones. Original edge selection makes",
                           "displayed-edge conservation conditional on original strong associations.")
  )
  contexts
}
