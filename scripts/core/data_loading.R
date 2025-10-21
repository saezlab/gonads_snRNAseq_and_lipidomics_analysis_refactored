#' Core Data Loading Functions
#'
#' This module provides unified data loading functions for RNA-seq and lipidomic data.
#' It consolidates the data loading logic from multiple scripts (seurat_to_df*.R,
#' pseudobulk_DA*.R, lipidomic_DA*.R) into reusable, parameterized functions.

library(readr)
library(readxl)
library(dplyr)

#' Load RNA Expression Data
#'
#' Loads scaled RNA expression data, handling split files if necessary.
#' The RNA data may be split across multiple CSV files due to GitHub size constraints.
#'
#' @param data_dir Character. Path to directory containing RNA data files.
#'                 Default: "data/RNA"
#' @param split_files Logical. Whether data is split across multiple files.
#'                    Default: TRUE
#' @param file_prefix Character. Prefix for split files. Default: "scaled_RNA"
#' @param n_splits Integer. Number of split files. Default: 3
#' @param single_file Character. Name of single file if split_files=FALSE.
#'                    Default: "scaled_RNA.csv"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Matrix with genes as rows and samples as columns. Row names are gene symbols.
#'
#' @details
#' The function handles two scenarios:
#' - Split files: Reads multiple CSVs and merges by "gene" column
#' - Single file: Reads one CSV file
#'
#' Split files are expected to have:
#' - File 1: gene column + samples 1 to n
#' - File 2: gene column + samples n+1 to m
#' - File 3: gene column + samples m+1 to end
#'
#' @examples
#' rna_data <- load_rna_data()  # Default: load from data/RNA, 3 split files
#' rna_data <- load_rna_data(split_files = FALSE, single_file = "my_rna.csv")
load_rna_data <- function(data_dir = "data/RNA",
                          split_files = TRUE,
                          file_prefix = "scaled_RNA",
                          n_splits = 3,
                          single_file = "scaled_RNA.csv",
                          verbose = TRUE) {

  if (split_files) {
    if (verbose) cat("Loading split RNA data files...\n")

    # Load first file
    file_path_1 <- file.path(data_dir, paste0(file_prefix, "_1.csv"))
    if (verbose) cat("  Loading:", file_path_1, "\n")
    rna_data <- as.data.frame(read_csv(file_path_1, show_col_types = FALSE))

    # Load and merge remaining files
    for (i in 2:n_splits) {
      file_path <- file.path(data_dir, paste0(file_prefix, "_", i, ".csv"))
      if (verbose) cat("  Loading:", file_path, "\n")
      rna_part <- as.data.frame(read_csv(file_path, show_col_types = FALSE))
      rna_data <- merge(rna_data, rna_part, by = "gene")
      rm(rna_part)  # Free memory
    }

  } else {
    # Load single file
    file_path <- file.path(data_dir, single_file)
    if (verbose) cat("Loading RNA data from:", file_path, "\n")
    rna_data <- as.data.frame(read_csv(file_path, show_col_types = FALSE))
  }

  # Set gene names as row names and remove gene column
  row.names(rna_data) <- rna_data$gene
  rna_data <- rna_data[, -1]

  if (verbose) {
    cat("RNA data loaded successfully:\n")
    cat("  Genes:", nrow(rna_data), "\n")
    cat("  Samples:", ncol(rna_data), "\n")
  }

  return(rna_data)
}


#' Load Lipidomic Data
#'
#' Loads and preprocesses lipidomic abundance data.
#'
#' @param file_path Character. Path to lipidomic data file.
#'                  Default: "data/metabolomics/Lipid_data_file_JPC_cleaned.csv"
#' @param log_transform Logical. Apply log2 transformation. Default: TRUE
#' @param sample_subset Character vector. Sample names to keep. If NULL, keep all.
#'                      Default: NULL
#' @param clean_names Logical. Clean sample names (replace "-0" with "-", "-" with "_").
#'                    Default: TRUE
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Matrix with lipids as rows and samples as columns.
#'
#' @details
#' Processing steps:
#' 1. Load CSV file (lipids in rows, samples in columns after transpose)
#' 2. Transpose to lipids as rows, samples as columns
#' 3. Optional: Log2 transform
#' 4. Optional: Clean sample names (remove "-0", replace "-" with "_")
#' 5. Optional: Subset to specified samples
#'
#' @examples
#' lipid_data <- load_lipidomic_data()
#' lipid_data <- load_lipidomic_data(log_transform = FALSE)
load_lipidomic_data <- function(file_path = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
                                 log_transform = TRUE,
                                 sample_subset = NULL,
                                 clean_names = TRUE,
                                 verbose = TRUE) {

  if (verbose) cat("Loading lipidomic data from:", file_path, "\n")

  # Load and transpose
  lipidomic <- as.data.frame(read_csv(file_path, show_col_types = FALSE))
  row.names(lipidomic) <- lipidomic$`Unnamed: 0`
  lipidomic <- as.data.frame(t(lipidomic[, -1]))

  # Log2 transform if requested
  if (log_transform) {
    if (verbose) cat("  Applying log2 transformation...\n")
    lipidomic <- log2(lipidomic)
  }

  # Clean sample names if requested
  if (clean_names) {
    if (verbose) cat("  Cleaning sample names...\n")
    names(lipidomic) <- gsub("-0", "-", names(lipidomic))
    names(lipidomic) <- gsub("-", "_", names(lipidomic))
  }

  # Subset samples if requested
  if (!is.null(sample_subset)) {
    if (verbose) cat("  Subsetting to", length(sample_subset), "samples...\n")
    lipidomic <- lipidomic[, sample_subset]
  }

  if (verbose) {
    cat("Lipidomic data loaded successfully:\n")
    cat("  Lipids:", nrow(lipidomic), "\n")
    cat("  Samples:", ncol(lipidomic), "\n")
  }

  return(lipidomic)
}


#' Load Sample Metadata
#'
#' Loads sample metadata and optionally merges with RNA metadata.
#'
#' @param overview_file Character. Path to overview data Excel file.
#'                      Default: "data/Suppl_Table_1_OverviewData.xlsx"
#' @param rna_metadata_file Character. Path to RNA metadata CSV file. If NULL,
#'                          only loads overview data. Default: NULL
#' @param exclude_wt Logical. Exclude wild-type samples (group contains "-WT").
#'                   Default: FALSE
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame with sample metadata.
#'
#' @details
#' When rna_metadata_file is provided:
#' - Loads RNA metadata from CSV
#' - Merges with overview data by sample_name (RNA) and sequencing_name (overview)
#' - Removes 4th column from overview data (column without header in Excel)
#'
#' When rna_metadata_file is NULL:
#' - Only loads overview data
#' - Renames first column to "sample"
#'
#' @examples
#' # For pseudobulk RNA analysis (with RNA metadata)
#' metadata <- load_metadata(rna_metadata_file = "data/RNA/meta_data_df.csv")
#'
#' # For lipidomic analysis (overview only)
#' metadata <- load_metadata()
load_metadata <- function(overview_file = "data/Suppl_Table_1_OverviewData.xlsx",
                          rna_metadata_file = NULL,
                          exclude_wt = FALSE,
                          verbose = TRUE) {

  if (verbose) cat("Loading sample metadata...\n")

  # Load overview data
  if (verbose) cat("  Loading overview data from:", overview_file, "\n")
  overview_data <- as.data.frame(read_excel(overview_file))

  # If RNA metadata provided, merge
  if (!is.null(rna_metadata_file)) {
    # For RNA analysis: remove 4th column before merging
    overview_data <- overview_data[, -4]  # Remove unnamed 4th column

    if (verbose) cat("  Loading RNA metadata from:", rna_metadata_file, "\n")
    rna_metadata <- as.data.frame(read_csv(rna_metadata_file, show_col_types = FALSE))

    # Merge by sample names
    metadata <- merge(rna_metadata, overview_data,
                     by.x = "sample_name", by.y = "sequencing_name")
  } else {
    # For lipidomic analysis: select columns directly (don't remove column 4 first)
    metadata <- overview_data[, c(2, 3, 4, 5)]  # Select relevant columns
    names(metadata)[1] <- "sample"
  }

  # Exclude wild-type samples if requested
  if (exclude_wt) {
    if (verbose) cat("  Excluding wild-type samples...\n")
    metadata <- metadata[!grepl("-WT", metadata$group), ]
  }

  if (verbose) {
    cat("Metadata loaded successfully:\n")
    cat("  Samples:", nrow(metadata), "\n")
    if ("cell_type" %in% names(metadata)) {
      cat("  Cell types:", paste(unique(metadata$cell_type), collapse = ", "), "\n")
    }
  }

  return(metadata)
}


#' Prepare Pseudobulk Data for a Specific Cell Type
#'
#' Subsets RNA expression data to samples of a specific cell type.
#'
#' @param rna_data Matrix. RNA expression data (genes x samples).
#' @param metadata Data frame. Sample metadata with columns "sample" and "cell_type".
#' @param cell_type Character. Cell type to extract (e.g., "Hepatocytes").
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Matrix with genes as rows and samples of the specified cell type as columns.
#'
#' @examples
#' hepatocyte_data <- prepare_pseudobulk(rna_data, metadata, "Hepatocytes")
prepare_pseudobulk <- function(rna_data, metadata, cell_type, verbose = TRUE) {

  if (verbose) cat("Preparing pseudobulk data for:", cell_type, "\n")

  # Get samples for this cell type
  samples_in_celltype <- metadata %>%
    filter(cell_type == !!cell_type) %>%
    pull(sample)

  # Subset RNA data
  pseudobulk_data <- rna_data[, colnames(rna_data) %in% samples_in_celltype]

  if (verbose) {
    cat("  Samples in", cell_type, ":", ncol(pseudobulk_data), "\n")
  }

  return(pseudobulk_data)
}


#' Filter Samples by Criteria
#'
#' Filters sample data based on metadata criteria.
#'
#' @param data Matrix or data frame. Expression/abundance data with samples as columns.
#' @param metadata Data frame. Sample metadata.
#' @param criteria Named list. Filtering criteria (e.g., list(cell_type = "Hepatocytes")).
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Filtered data matrix.
#'
#' @examples
#' # Filter to hepatocytes from testes
#' filtered <- filter_samples(rna_data, metadata,
#'                            list(cell_type = "Hepatocytes", gonadal_type = "Testes"))
filter_samples <- function(data, metadata, criteria, verbose = TRUE) {

  if (verbose) cat("Filtering samples by criteria...\n")

  # Apply each criterion
  filtered_metadata <- metadata
  for (column in names(criteria)) {
    value <- criteria[[column]]
    if (verbose) cat("  ", column, "=", value, "\n")
    filtered_metadata <- filtered_metadata[filtered_metadata[[column]] == value, ]
  }

  # Subset data
  filtered_data <- data[, colnames(data) %in% filtered_metadata$sample]

  if (verbose) {
    cat("Samples after filtering:", ncol(filtered_data), "\n")
  }

  return(filtered_data)
}
