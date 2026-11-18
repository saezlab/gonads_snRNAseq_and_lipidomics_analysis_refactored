# Original paired FCG Cyp2c38 reference for workflow 07b.
# Sourcing this file defines helpers without running an analysis. Shared support
# functions are the caller's responsibility; these helpers use base R and stats.

#' Convert annotation formatting without changing lipid structural identity.
#'
#' Select the final pipe-delimited annotation and change "CLASS structure" to
#' "CLASS(structure)". Chain order, slash/underscore, ether prefixes and oxygenation
#' remain intact. Exclude internal standards from their ORIGINAL annotation before
#' calling this function: a standard's final annotation may omit its IS prefix.
#' Formatting agreement does not independently establish chemical identity.
wp1_canonical_lipid <- function(raw_name) {
  label <- trimws(as.character(raw_name))
  label <- trimws(sub("^.*\\|", "", label))
  sub("^([^[:space:]()]+)[[:space:]]+(.+)$", "\\1(\\2)", label)
}

#' Load fixed Cyp2c38 predictions and their paired FCG source measurements.
#'
#' @param root Analysis repository root, independent of the working directory.
#' @return A list with coefficients, named gene expression, a log2 lipid-by-mouse
#'   matrix, ordered metadata and the maximum pooled archive discrepancy.
#'
#' RNA is already scaled. Original lipid abundances receive the log2 transformation
#' and mouse-name normalization used by load_lipidomic_data(). All 828 features are
#' fixed by the existing full MOFA export, independent of any KO response. Female
#' and male labels mean Ovaries and Testes in this FCG reference. Genotype metadata
#' comes from the standard WF02 hepatocyte PCA export, whose 43 group/allosome rows
#' were independently checked against Suppl_Table_1_OverviewData.xlsx on 2026-10-06.
#' A changed sample roster, feature universe, pairing or pooled prediction fails
#' explicitly rather than silently producing a different validation reference.
wp1_load_reference <- function(root = ".") {
  if (length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop("root must identify one existing analysis repository directory.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  reference_file <- file.path(root, "results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv")
  rna_files <- file.path(root, "data/RNA", paste0("scaled_RNA_", seq_len(3L), ".csv"))
  lipid_file <- file.path(root, "data/metabolomics/Lipid_data_file_JPC_cleaned.csv")
  metadata_file <- file.path(root, "data/RNA/meta_data_df.csv")
  genotype_file <- file.path(root, "results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv")
  required_files <- c(reference_file, rna_files, lipid_file, metadata_file, genotype_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop("Missing Cyp2c38 reference input(s): ", paste(missing_files, collapse = ", "))
  }

  archive <- utils::read.csv(reference_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("node1", "node2", "mor") %in% names(archive))) {
    stop("Archived reference requires named columns node1, node2 and mor.")
  }
  archive <- archive[!is.na(archive$node1) & archive$node1 == "Cyp2c38", , drop = FALSE]
  if (nrow(archive) != 828L || anyNA(archive$node2) ||
      any(!nzchar(archive$node2)) || anyDuplicated(archive$node2)) {
    stop("Expected 828 unique, nonempty archived Cyp2c38 lipid predictions.")
  }
  if (!is.numeric(archive$mor) || any(!is.finite(archive$mor))) {
    stop("Archived Cyp2c38 correlation coefficients must be finite numeric values.")
  }

  # Only keep Cyp2c38 while reading the large numerical CSVs in bounded chunks.
  # The documented source schema has one gene record per line, with gene first.
  read_gene_row <- function(path) {
    connection <- file(path, open = "rt")
    on.exit(close(connection))
    header <- readLines(connection, n = 1L, warn = FALSE)
    hit <- character(0L)
    repeat {
      chunk <- readLines(connection, n = 1024L, warn = FALSE)
      if (!length(chunk)) break
      hit <- c(hit, chunk[grepl('^("Cyp2c38"|Cyp2c38),', chunk)])
    }
    if (length(header) != 1L || length(hit) != 1L) {
      stop("Expected exactly one Cyp2c38 RNA row in ", path)
    }
    row <- utils::read.csv(text = paste(c(header, hit), collapse = "\n"),
                           check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(names(row)[1L], "gene") || row$gene != "Cyp2c38" ||
        anyDuplicated(names(row))) {
      stop("Unexpected gene/sample schema in ", path)
    }
    row
  }
  pieces <- lapply(rna_files, read_gene_row)
  sample_columns <- unlist(lapply(pieces, function(x) names(x)[-1L]), use.names = FALSE)
  if (anyDuplicated(sample_columns)) stop("RNA split files contain duplicate sample columns.")
  rna <- Reduce(function(x, y) merge(x, y, by = "gene", sort = FALSE), pieces)
  hepatocyte_columns <- names(rna)[grepl("_Hepatocytes_", names(rna), fixed = TRUE)]
  if (length(hepatocyte_columns) != 43L) {
    stop("Expected 43 hepatocyte RNA samples in the frozen FCG reference.")
  }
  gene <- unlist(rna[1L, hepatocyte_columns, drop = FALSE], use.names = FALSE)
  if (!is.numeric(gene) || any(!is.finite(gene))) {
    stop("Cyp2c38 RNA measurements must be finite numeric values.")
  }
  sample_ids <- gsub("-", "_", gsub("_Hepatocytes|_Testes|_Ovaries", "", hepatocyte_columns))
  if (anyDuplicated(sample_ids)) stop("RNA normalization creates duplicate mouse IDs.")
  names(gene) <- sample_ids

  # Require explicitly named metadata fields rather than positional columns.
  rna_metadata <- utils::read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
  metadata_columns <- c("sample", "sample_name", "cell_type", "gonadal_type")
  if (!all(metadata_columns %in% names(rna_metadata))) {
    stop("RNA metadata requires columns: ", paste(metadata_columns, collapse = ", "))
  }
  hepa_metadata <- rna_metadata[!is.na(rna_metadata$cell_type) &
                                  rna_metadata$cell_type == "Hepatocytes", , drop = FALSE]
  if (anyNA(hepa_metadata$sample) || anyDuplicated(hepa_metadata$sample)) {
    stop("Hepatocyte metadata contains missing or duplicated sample identifiers.")
  }
  metadata_index <- match(hepatocyte_columns, hepa_metadata$sample)
  if (anyNA(metadata_index)) stop("Some hepatocyte RNA columns have no matching metadata.")
  hepa_metadata <- hepa_metadata[metadata_index, , drop = FALSE]
  metadata <- data.frame(
    sample_id = gsub("-", "_", hepa_metadata$sample_name),
    gonadal_type = hepa_metadata$gonadal_type,
    stringsAsFactors = FALSE
  )
  if (!identical(metadata$sample_id, sample_ids)) {
    stop("RNA column-derived and metadata-derived mouse IDs do not agree.")
  }
  if (anyNA(metadata$gonadal_type) ||
      sum(metadata$gonadal_type == "Testes") != 22L ||
      sum(metadata$gonadal_type == "Ovaries") != 21L) {
    stop("Expected 22 Testes and 21 Ovaries mice in the frozen reference.")
  }

  # The PCA export preserves overview metadata. Validate its pairing and
  # categorical identities; do not infer sex chromosomes from gonadal type.
  genotype <- utils::read.csv(genotype_file, check.names = FALSE, stringsAsFactors = FALSE)
  genotype_columns <- c("sequencing_name", "group", "gonadal_type", "allosome_type")
  if (!all(genotype_columns %in% names(genotype))) {
    stop("Genotype metadata requires columns: ", paste(genotype_columns, collapse = ", "))
  }
  genotype_ids <- gsub("-", "_", genotype$sequencing_name)
  if (anyNA(genotype_ids) || anyDuplicated(genotype_ids)) {
    stop("Genotype metadata has missing or duplicated mouse IDs.")
  }
  genotype_index <- match(sample_ids, genotype_ids)
  if (anyNA(genotype_index)) stop("Some paired FCG mice have no recorded genotype metadata.")
  genotype <- genotype[genotype_index, , drop = FALSE]
  if (!identical(genotype$gonadal_type, metadata$gonadal_type) ||
      anyNA(genotype$group) || anyNA(genotype$allosome_type) ||
      !all(genotype$group %in% c("XXO", "XYO", "XXT", "XYT")) ||
      !all(genotype$allosome_type %in% c("XX", "XY")) ||
      any(substr(genotype$group, 1L, 2L) != genotype$allosome_type) ||
      any(substr(genotype$group, 3L, 3L) !=
            ifelse(metadata$gonadal_type == "Ovaries", "O", "T"))) {
    stop("Recorded FCG group, gonad and chromosome identities are inconsistent.")
  }
  group_counts <- table(factor(genotype$group, levels = c("XXO", "XYO", "XXT", "XYT")))
  if (!identical(as.integer(group_counts), c(10L, 11L, 11L, 11L))) {
    stop("Expected original FCG groups XXO=10, XYO=11, XXT=11, XYT=11.")
  }
  metadata$group <- genotype$group
  metadata$allosome_type <- genotype$allosome_type
  rownames(metadata) <- NULL

  lipid_table <- utils::read.csv(lipid_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Unnamed: 0" %in% names(lipid_table) || anyDuplicated(names(lipid_table))) {
    stop("Original lipid table requires unique feature names and Unnamed: 0 mouse IDs.")
  }
  lipid_sample_ids <- gsub("-", "_", gsub("-0", "-", lipid_table[["Unnamed: 0"]], fixed = TRUE))
  if (anyNA(lipid_sample_ids) || anyDuplicated(lipid_sample_ids)) {
    stop("Original lipid sample IDs are missing or duplicated after normalization.")
  }
  mouse_index <- match(sample_ids, lipid_sample_ids)
  feature_index <- match(archive$node2, names(lipid_table))
  if (anyNA(mouse_index)) stop("Some FCG hepatocyte mice have no paired lipid sample.")
  if (anyNA(feature_index)) stop("Some archived reference lipids are absent from original measurements.")
  paired_values <- as.matrix(lipid_table[mouse_index, feature_index, drop = FALSE])
  if (!is.numeric(paired_values) || any(!is.finite(paired_values)) || any(paired_values <= 0)) {
    stop("Paired original lipid values must be finite and positive before log2 transformation.")
  }
  lipids <- t(log2(paired_values))
  dimnames(lipids) <- list(archive$node2, sample_ids)

  strata <- list(pooled = rep(TRUE, length(gene)),
                 female = metadata$gonadal_type == "Ovaries",
                 male = metadata$gonadal_type == "Testes")
  coefficients <- data.frame(reference_lipid = archive$node2, stringsAsFactors = FALSE)
  for (stratum in names(strata)) {
    keep <- strata[[stratum]]
    correlations <- vapply(seq_len(nrow(lipids)), function(i) {
      stats::cor(gene[keep], lipids[i, keep], method = "pearson", use = "pairwise.complete.obs")
    }, numeric(1L))
    if (any(!is.finite(correlations))) {
      stop("Undefined Cyp2c38 correlation(s) in the ", stratum, " reference stratum.")
    }
    coefficients[[paste0("r_", stratum)]] <- correlations
  }
  for (stratum in names(strata)) {
    coefficients[[paste0("n_", stratum)]] <- sum(strata[[stratum]])
  }
  coefficients$r_archived <- archive$mor
  max_archive_difference <- max(abs(coefficients$r_pooled - coefficients$r_archived))
  if (!is.finite(max_archive_difference) || max_archive_difference > 1e-12) {
    stop("Reconstructed Cyp2c38 pooled coefficients fail the archived check: max difference = ",
         format(max_archive_difference, digits = 16L), ".")
  }
  cat("Cyp2c38 reference: ", nrow(lipids), " lipids x ", ncol(lipids),
      " paired mice (21 Ovaries, 22 Testes); max archive difference = ",
      format(max_archive_difference, digits = 4L), "\n", sep = "")
  list(coefficients = coefficients, gene = gene, lipids = lipids,
       metadata = metadata, max_archive_difference = max_archive_difference)
}
