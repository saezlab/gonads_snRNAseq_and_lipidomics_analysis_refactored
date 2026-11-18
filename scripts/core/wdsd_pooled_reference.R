# Frozen, non-zonated Figure 3 FCG gene-lipid references.
# Source scripts/support_functions.R before calling this helper. Sourcing this
# file defines functions only; it never writes or replaces archived outputs.

#' Check archived pooled coefficients against paired source measurements.
#' @param rna Gene-by-mouse numeric matrix.
#' @param lipids Log2 lipid-by-mouse numeric matrix, in the same mouse order.
#' @param edges Data frame with gene, reference_lipid and r_reference.
#' @param tolerance Maximum accepted absolute replay discrepancy.
wp2_pooled_replay <- function(rna, lipids, edges, tolerance = 1e-10) {
  if (!is.matrix(rna) || !is.numeric(rna) || !is.matrix(lipids) ||
      !is.numeric(lipids) || any(!is.finite(rna)) ||
      any(!is.finite(lipids)) || is.null(colnames(rna)) ||
      !identical(colnames(rna), colnames(lipids)) ||
      anyDuplicated(colnames(rna)) || ncol(rna) < 3L) {
    stop("Pooled RNA and log2 lipids require finite, identically ordered unique mice.")
  }
  required <- c("gene", "reference_lipid", "r_reference")
  if (!is.data.frame(edges) || !all(required %in% names(edges)) ||
      anyNA(edges[, required]) || anyDuplicated(paste(edges$gene, edges$reference_lipid,
                                                     sep = "\r")) ||
      !is.numeric(edges$r_reference) || any(!is.finite(edges$r_reference))) {
    stop("Pooled edge identities and archived coefficients must be unique and finite.")
  }
  lookup <- cbind(match(edges$gene, rownames(rna)),
                  match(edges$reference_lipid, rownames(lipids)))
  if (anyNA(lookup)) stop("Every frozen pooled edge requires a measured gene and lipid.")
  observed <- stats::cor(t(rna), t(lipids), method = "pearson")
  delta <- max(abs(observed[lookup] - edges$r_reference))
  if (!is.finite(delta) || delta > tolerance) {
    stop("Pooled log2 reference coefficients fail archival replay; maximum error ",
         format(delta, digits = 16L), ".")
  }
  delta
}

#' Load the fixed Figure 3 species network and its paired 43-mouse FCG inputs.
#'
#' The current revised Figure 3 B labels the sparse species-level network as
#' "FULL SPECIES-LEVEL NETWORK" and Figure 3 A specifies the top three correlated
#' genes per lipid. This is the 192 gene-lipid branch of the integrated 267-edge
#' export, rather than the separate permissive 769-by-828, 33%-tail reference.
#' TF-to-gene and class-compressed edges do not enter species-level conservation.
#' The larger reference is returned separately with explicit provenance. Original
#' feature selection remains fixed and does not depend on WDSD observations.
#'
#' @return A list with rna (139 genes), lipids (64 Figure 3 species), ordered
#'   metadata, edges (8,896 sparse-precursor pairs with primary_edge marking the
#'   192 displayed pairs), figure3_edges, sparse_edges, compressed_gene_lipid_edges,
#'   permissive_inventory (the separate full33% network is inventoried, not part
#'   of this conservation comparison), replay_metrics, provenance and source_paths.
#'   All numerical lipid measurements are log2.
wp2_pooled_reference <- function(root = ".") {
  if (length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop("root must identify one existing analysis repository directory.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  paths <- c(
    full33_SIF = file.path(root, "results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv"),
    figure3_species_SIF = file.path(root, "results/mofa/networks/hepatocyte_TF_gene_lipid_network_full_SIF.csv"),
    figure3_species_ATT = file.path(root, "results/mofa/networks/hepatocyte_TF_gene_lipid_network_full_ATT.csv"),
    figure3_compressed_SIF = file.path(root, "results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv"),
    sparse_matrix = file.path(root, "results/mofa/networks/gene_lipid_correlations_factor1.csv"),
    RNA_1 = file.path(root, "data/RNA/scaled_RNA_1.csv"),
    RNA_2 = file.path(root, "data/RNA/scaled_RNA_2.csv"),
    RNA_3 = file.path(root, "data/RNA/scaled_RNA_3.csv"),
    original_lipids = file.path(root, "data/metabolomics/Lipid_data_file_JPC_cleaned.csv"),
    RNA_metadata = file.path(root, "data/RNA/meta_data_df.csv"),
    genotype_metadata = file.path(root, "results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv")
  )
  if (any(!file.exists(paths))) {
    stop("Missing frozen pooled source files: ", paste(paths[!file.exists(paths)], collapse = ", "))
  }
  read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                            stringsAsFactors = FALSE)
  assert_ids <- function(x, label) {
    if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
      stop(label, " must be unique, nonempty identifiers.")
    }
  }
  full <- read_csv(paths[["full33_SIF"]])
  integrated <- read_csv(paths[["figure3_species_SIF"]])
  att <- read_csv(paths[["figure3_species_ATT"]])
  compressed <- read_csv(paths[["figure3_compressed_SIF"]])
  if (!all(c("node1", "node2", "mor") %in% names(full)) ||
      !all(c("node1", "node2", "mor") %in% names(integrated)) ||
      !all(c("node", "type") %in% names(att))) {
    stop("Frozen pooled network archives require their recorded named columns.")
  }
  raw <- read_csv(paths[["original_lipids"]])
  if (!identical(names(raw)[1L], "Unnamed: 0") || ncol(raw) != 1255L) {
    stop("Expected the original mouse-by-1,254-species lipid table.")
  }
  original_lipid_ids <- names(raw)[-1L]
  assert_ids(original_lipid_ids, "Original lipid identities")
  # In this archive, sparse gene-lipid edges run lipid -> gene; TF edges run
  # TF -> gene. Identify the branch by measured species, independently of ATT.
  displayed <- integrated[integrated$node1 %in% original_lipid_ids, , drop = FALSE]
  excluded <- integrated[!integrated$node1 %in% original_lipid_ids, , drop = FALSE]
  if (nrow(integrated) != 267L || nrow(displayed) != 192L || nrow(excluded) != 75L ||
      length(unique(displayed$node2)) != 25L ||
      length(unique(displayed$node1)) != 64L ||
      any(as.integer(table(displayed$node1)) != 3L)) {
    stop("The fixed Figure 3 network must retain 192 species edges and 75 TF edges.")
  }
  assert_ids(att$node, "Integrated ATT nodes")
  if (!all(att$type[match(displayed$node1, att$node)] == "lipid") ||
      !all(att$type[match(excluded$node1, att$node)] == "TF") ||
      !all(excluded$node2 %in% displayed$node2)) {
    stop("Measured-lipid and TF-to-gene branch identities disagree with the archive.")
  }
  full_genes <- unique(full$node1)
  full_lipids <- unique(full$node2)
  if (nrow(full) != 636732L || length(full_genes) != 769L ||
      length(full_lipids) != 828L || !all(full_lipids %in% original_lipid_ids)) {
    stop("Expected the frozen 769-gene by 828-lipid permissive reference.")
  }
  edges <- data.frame(gene = full$node1, reference_lipid = full$node2,
                      r_reference = full$mor, stringsAsFactors = FALSE)
  assert_ids(paste(edges$gene, edges$reference_lipid, sep = "\r"), "Full33% edge identities")
  figure3_edges <- data.frame(gene = displayed$node2,
                              reference_lipid = displayed$node1,
                              r_reference = displayed$mor,
                              primary_edge = TRUE,
                              network_universe = "figure3_displayed_species",
                              stringsAsFactors = FALSE)
  assert_ids(paste(figure3_edges$gene, figure3_edges$reference_lipid, sep = "\r"),
             "Figure 3 displayed edge identities")
  primary_key <- paste(figure3_edges$gene, figure3_edges$reference_lipid, sep = "\r")
  full_key <- paste(edges$gene, edges$reference_lipid, sep = "\r")
  primary_index <- match(primary_key, full_key)
  if (anyNA(primary_index) ||
      max(abs(edges$r_reference[primary_index] - figure3_edges$r_reference)) > 1e-10) {
    stop("The sparse species branch must be an unchanged subset of the permissive archive.")
  }
  edges$primary_edge <- full_key %in% primary_key
  edges$network_universe <- "permissive_33percent_all_pairs"
  sparse_mat <- as.matrix(utils::read.csv(paths[["sparse_matrix"]], row.names = 1L,
                                          check.names = FALSE))
  if (!is.numeric(sparse_mat) || any(!is.finite(sparse_mat)) ||
      !identical(dim(sparse_mat), c(139L, 64L)) ||
      !all(rownames(sparse_mat) %in% full_genes) ||
      !setequal(colnames(sparse_mat), unique(displayed$node1))) {
    stop("Expected the frozen 139-gene by 64-species sparse precursor matrix.")
  }
  sparse_edges <- data.frame(gene = rep(rownames(sparse_mat), times = ncol(sparse_mat)),
                             reference_lipid = rep(colnames(sparse_mat), each = nrow(sparse_mat)),
                             r_reference = as.vector(sparse_mat),
                             stringsAsFactors = FALSE)
  sparse_edges$primary_edge <- paste(sparse_edges$gene, sparse_edges$reference_lipid,
                                     sep = "\r") %in% primary_key
  sparse_edges$display_gene_edge <- sparse_edges$gene %in% figure3_edges$gene
  sparse_edges$network_universe <- "sparse_selection_unthresholded_precursor"
  # The archived sparse network used top-three absolute correlations, allowing
  # the recorded values to establish membership rather than adding an r cutoff.
  expected_top <- unlist(lapply(seq_len(ncol(sparse_mat)), function(j) {
    cut <- sort(abs(sparse_mat[, j]), decreasing = TRUE)[3L]
    paste(rownames(sparse_mat)[abs(sparse_mat[, j]) >= cut],
          colnames(sparse_mat)[j], sep = "\r")
  }), use.names = FALSE)
  if (!setequal(expected_top, primary_key)) {
    stop("Figure 3 sparse species membership no longer replays the top-three rule.")
  }
  genes <- rownames(sparse_mat)
  # Read only requested gene rows from each of the original split RNA tables.
  read_rows <- function(path) {
    con <- file(path, "rt")
    on.exit(close(con))
    header <- readLines(con, n = 1L, warn = FALSE)
    selected <- character()
    repeat {
      chunk <- readLines(con, n = 1024L, warn = FALSE)
      if (!length(chunk)) break
      ids <- sub('^"?([^",]+)"?,.*$', "\\1", chunk)
      selected <- c(selected, chunk[ids %in% genes])
    }
    tab <- utils::read.csv(text = paste(c(header, selected), collapse = "\n"),
                           check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(names(tab)[1L], "gene") || anyDuplicated(names(tab)) ||
        nrow(tab) != length(genes) || !setequal(tab$gene, genes)) {
      stop("Selected pooled genes must be uniquely measured in every RNA split file.")
    }
    tab[match(genes, tab$gene), , drop = FALSE]
  }
  pieces <- lapply(paths[c("RNA_1", "RNA_2", "RNA_3")], read_rows)
  sample_columns <- unlist(lapply(pieces, function(p) names(p)[-1L]), use.names = FALSE)
  assert_ids(sample_columns, "Original RNA sample columns")
  rna <- do.call(cbind, lapply(pieces, function(p) as.matrix(p[, -1L, drop = FALSE])))
  rownames(rna) <- genes
  hepatocyte_columns <- colnames(rna)[grepl("_Hepatocytes_", colnames(rna), fixed = TRUE)]
  if (length(hepatocyte_columns) != 43L) stop("Expected 43 paired original hepatocyte mice.")
  rna <- rna[, hepatocyte_columns, drop = FALSE]
  rna_meta <- read_csv(paths[["RNA_metadata"]])
  required_meta <- c("sample", "sample_name", "cell_type", "gonadal_type")
  if (!all(required_meta %in% names(rna_meta))) stop("Missing original RNA metadata fields.")
  hepa <- rna_meta[rna_meta$cell_type == "Hepatocytes", , drop = FALSE]
  assert_ids(hepa$sample, "Hepatocyte metadata samples")
  index <- match(hepatocyte_columns, hepa$sample)
  if (anyNA(index)) stop("Every original RNA column requires hepatocyte metadata.")
  hepa <- hepa[index, , drop = FALSE]
  mouse_ids <- gsub("-", "_", hepa$sample_name)
  from_columns <- gsub("-", "_", gsub("_Hepatocytes|_Testes|_Ovaries", "", hepatocyte_columns))
  assert_ids(mouse_ids, "Paired original mouse IDs")
  if (!identical(mouse_ids, from_columns) ||
      sum(hepa$gonadal_type == "Ovaries") != 21L ||
      sum(hepa$gonadal_type == "Testes") != 22L) {
    stop("Original pooled hepatocyte sample and gonadal identities disagree.")
  }
  genotype <- read_csv(paths[["genotype_metadata"]])
  geno_fields <- c("sequencing_name", "group", "gonadal_type", "allosome_type")
  if (!all(geno_fields %in% names(genotype))) stop("Missing recorded FCG genotype fields.")
  genotype_ids <- gsub("-", "_", genotype$sequencing_name)
  assert_ids(genotype_ids, "Original genotype mouse IDs")
  geno_index <- match(mouse_ids, genotype_ids)
  if (anyNA(geno_index)) stop("Every original pooled mouse requires recorded FCG genotype.")
  genotype <- genotype[geno_index, , drop = FALSE]
  if (!identical(genotype$gonadal_type, hepa$gonadal_type) ||
      !all(genotype$group %in% c("XXO", "XYO", "XXT", "XYT")) ||
      !all(genotype$allosome_type %in% c("XX", "XY")) ||
      any(substr(genotype$group, 1L, 2L) != genotype$allosome_type) ||
      any(substr(genotype$group, 3L, 3L) != ifelse(hepa$gonadal_type == "Ovaries", "O", "T")) ||
      !identical(as.integer(table(factor(genotype$group,
                                          levels = c("XXO", "XYO", "XXT", "XYT")))),
                  c(10L, 11L, 11L, 11L))) {
    stop("Recorded original group, chromosome and gonadal identities are inconsistent.")
  }
  raw_ids <- gsub("-", "_", gsub("-0", "-", raw[["Unnamed: 0"]], fixed = TRUE))
  assert_ids(raw_ids, "Original lipid mouse IDs")
  lipid_index <- match(mouse_ids, raw_ids)
  if (anyNA(lipid_index)) stop("Every original RNA mouse requires its paired lipid sample.")
  paired_raw <- t(as.matrix(raw[lipid_index, -1L, drop = FALSE]))
  if (!is.numeric(paired_raw) || any(!is.finite(paired_raw)) || any(paired_raw <= 0)) {
    stop("Original lipid measurements must be finite positive values before log2 transformation.")
  }
  dimnames(paired_raw) <- list(original_lipid_ids, mouse_ids)
  lipids <- log2(paired_raw[colnames(sparse_mat), , drop = FALSE])
  colnames(rna) <- mouse_ids
  metadata <- data.frame(mouse_id = mouse_ids, sample_id = mouse_ids,
                         rna_sample = hepatocyte_columns,
                         gonadal_type = hepa$gonadal_type,
                         sex = ifelse(hepa$gonadal_type == "Ovaries", "female", "male"),
                         group = genotype$group,
                         allosome_type = genotype$allosome_type,
                         stringsAsFactors = FALSE)
  displayed_delta <- wp2_pooled_replay(rna, lipids, figure3_edges)
  sparse_delta <- wp2_pooled_replay(rna, lipids, sparse_edges)
  raw_cor <- stats::cor(t(rna), t(paired_raw[colnames(sparse_mat), , drop = FALSE]),
                        method = "pearson")
  sparse_lookup <- cbind(match(sparse_edges$gene, rownames(raw_cor)),
                         match(sparse_edges$reference_lipid, colnames(raw_cor)))
  raw_delta <- max(abs(raw_cor[sparse_lookup] - sparse_edges$r_reference))
  displayed_classes <- figure3_edges
  displayed_classes$lipid_class <- sub("[(].*$", "", displayed_classes$reference_lipid)
  class_edges <- stats::aggregate(r_reference ~ gene + lipid_class,
                                  data = displayed_classes, FUN = mean)
  classes <- unique(displayed_classes$lipid_class)
  compressed_gene_lipid <- compressed[compressed$node1 %in% classes, , drop = FALSE]
  ck <- paste(class_edges$gene, class_edges$lipid_class, sep = "\r")
  lookup_class <- match(paste(compressed_gene_lipid$node2,
                              compressed_gene_lipid$node1, sep = "\r"), ck)
  if (nrow(class_edges) != 75L || nrow(compressed_gene_lipid) != 75L ||
      anyNA(lookup_class) ||
      max(abs(compressed_gene_lipid$mor - class_edges$r_reference[lookup_class])) > 1e-10) {
    stop("Compressed display does not reproduce means of frozen Figure 3 species coefficients.")
  }
  class_edges$network_universe <- "display_class_means_of_sparse_species"
  replay_metrics <- data.frame(
    network_universe = c("figure3_displayed_species", "sparse_selection_unthresholded_precursor"),
    n_mice = 43L, n_genes = c(25L, 139L), n_lipids = c(64L, 64L),
    n_edges = c(192L, 8896L),
    log2_max_abs_archive_difference = c(displayed_delta, sparse_delta),
    lipid_scale = "log2", stringsAsFactors = FALSE)
  provenance <- list(
    primary_network = paste("Current Figure 3 A states top 3 correlated genes per lipid;",
                            "Figure 3 B labels the sparse species-level display FULL SPECIES-LEVEL NETWORK.",
                            "Integrated export contains 192 species pairs and 75 TF-gene pairs."),
    primary_selection = paste("MOFA Factor1 top 2.5% in each tail, lipid-metabolism gene filter,",
                              "then top 3 absolute gene correlations per lipid; no additional 0.65 cutoff."),
    full_selection = paste("Separate permissive archive: Factor1 top 33% in each tail,",
                           "lipid-metabolism gene filter, all 769-by-828 pairs without correlation filtering."),
    compressed_selection = paste("75 class-gene coefficients are arithmetic means of the original sparse",
                                 "species coefficients, not correlations with class-summed lipid abundance."),
    pairing = "43 original FCG mice: 21 ovaries and 22 testes; chromosomes and gonads retained in metadata.",
    sources = paths,
    raw_max_abs_archive_difference = raw_delta,
    n_excluded_TF_gene_edges = 75L,
    n_class_gene_edges = 75L,
    figure_pdf = file.path(root, "docs/reference_figures/FCGL_Fig3.ai")
  )
  if (file.exists(provenance$figure_pdf)) {
    provenance$figure_pdf_md5 <- unname(tools::md5sum(provenance$figure_pdf))
  }
  permissive_inventory <- data.frame(
    network_universe = "permissive_33percent_all_pairs",
    n_mice = 43L, n_genes = 769L, n_lipids = 828L, n_edges = 636732L,
    included_in_conservation = FALSE,
    source = paths[["full33_SIF"]], stringsAsFactors = FALSE)
  cat("Figure 3 pooled reference: 43 paired FCG mice, 192 displayed species pairs; ",
      "139 x 64 precursor and 769 x 828 permissive universes; maximum log2 replay error ",
      format(max(replay_metrics$log2_max_abs_archive_difference), digits = 4L), "\n", sep = "")
  list(rna = rna, lipids = lipids, metadata = metadata,
       edges = sparse_edges, figure3_edges = figure3_edges, sparse_edges = sparse_edges,
       compressed_gene_lipid_edges = class_edges,
       permissive_inventory = permissive_inventory,
       replay_metrics = replay_metrics, provenance = provenance,
       source_paths = paths)
}
