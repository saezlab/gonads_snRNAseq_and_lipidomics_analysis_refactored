################################################################################
# Workflow 1: Data Preparation and Quality Control
################################################################################
#
# This workflow prepares the raw data for downstream analysis and performs
# comprehensive quality control checks.
#
# Author: Refactored codebase
# Date: 2025
#
# Replaces legacy scripts:
#   - scripts/seurat_to_df.R
#   - Initial QC steps scattered across multiple scripts
#
# Key steps:
#   1. Load Seurat object with snRNA-seq data
#   2. Extract and validate metadata
#   3. Perform quality control (cells per sample, genes detected, etc.)
#   4. Extract scaled expression data
#   5. Create data files for downstream workflows
#   6. Generate initial exploratory visualizations
#   7. Validate data integrity
#
# Data generated:
#   - data/RNA/scaled_RNA_1.csv (first third of genes)
#   - data/RNA/scaled_RNA_2.csv (second third of genes)
#   - data/RNA/scaled_RNA_3.csv (third third of genes)
#   - data/RNA/meta_data_df.csv (cell-level metadata)
#   - results/qc/data_quality_report.txt
#
# Dataset scale:
#   - 466,986 high-quality nuclei
#   - 48 samples (43 FCG + 5 WT males)
#   - 5 major liver cell types
#   - 20,123 genes detected
#
################################################################################

# Clear workspace
rm(list = ls())

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(readxl)

################################################################################
# 1. CONFIGURATION
################################################################################

cat(strrep("=", 80), "\n")
cat("WORKFLOW 1: DATA PREPARATION AND QUALITY CONTROL\n")
cat(strrep("=", 80), "\n\n")

# Load configuration
source("scripts/config/analysis_config.R")
source("scripts/config/datasets.R")

# Load support functions
source("scripts/support_functions.R")

# Load compatibility functions for Seurat v4/v5
source("scripts/seurat_compatibility.R")

# Load visualization functions
source("scripts/visualization/dimred_plots.R")

# Define output directories
qc_dir <- "results/qc"
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

data_rna_dir <- "data/RNA"
if (!dir.exists(data_rna_dir)) {
  dir.create(data_rna_dir, recursive = TRUE)
}

################################################################################
# 2. LOAD SEURAT OBJECT
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("2. LOADING SEURAT OBJECT\n")
cat(strrep("-", 80), "\n\n")

seurat_path <- "data/RNA/fcg_pseudobulks.rds"

if (!file.exists(seurat_path)) {
  cat("ERROR: Seurat object not found at:", seurat_path, "\n")
  cat("Please ensure the Seurat object is in the correct location.\n")
  stop("Seurat object not found!")
}

cat("Loading Seurat object from:", seurat_path, "\n")
cat("This may take several minutes for large datasets...\n\n")

seurat_obj <- readRDS(seurat_path)

cat("Seurat object loaded successfully!\n\n")

# Display basic information
cat("Seurat object information:\n")
cat("  NOTE: This is a PSEUDOBULK object (not single-cell data)\n")
cat("  Pseudobulk samples (sample × cell type):", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")
cat("  Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")

# Check available data slots (compatible with Seurat v4 and v5)
cat("\n  Available data slots:\n")
if ("RNA" %in% names(seurat_obj@assays)) {
  # Check if Seurat v5 (Assay5 class) or v4 (Assay class)
  assay_class <- class(seurat_obj@assays$RNA)[1]
  cat("    Assay class:", assay_class, "\n")

  if (assay_class == "Assay5") {
    # Seurat v5 uses layers
    tryCatch({
      layers <- Layers(seurat_obj, assay = "RNA")
      cat("    Available layers:", paste(layers, collapse = ", "), "\n")
      cat("    - counts:", ifelse("counts" %in% layers, "✓", "✗"), "\n")
      cat("    - data:", ifelse("data" %in% layers, "✓", "✗"), "\n")
      cat("    - scale.data:", ifelse("scale.data" %in% layers, "✓", "✗"), "\n")
    }, error = function(e) {
      cat("    (Unable to check layers)\n")
    })
  } else {
    # Seurat v4 uses direct slots
    tryCatch({
      cat("    - counts:", ifelse(!is.null(seurat_obj@assays$RNA@counts) &&
                                  length(seurat_obj@assays$RNA@counts) > 0, "✓", "✗"), "\n")
      cat("    - data:", ifelse(!is.null(seurat_obj@assays$RNA@data) &&
                                length(seurat_obj@assays$RNA@data) > 0, "✓", "✗"), "\n")
      cat("    - scale.data:", ifelse(!is.null(seurat_obj@assays$RNA@scale.data) &&
                                      length(seurat_obj@assays$RNA@scale.data) > 0, "✓", "✗"), "\n")
    }, error = function(e) {
      cat("    (Unable to check slots)\n")
    })
  }
}

# Check metadata
cat("\n  Metadata columns:\n")
cat("   ", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")

################################################################################
# 3. EXTRACT AND VALIDATE METADATA
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("3. EXTRACTING METADATA\n")
cat(strrep("-", 80), "\n\n")

# Extract metadata
metadata <- seurat_obj@meta.data

cat("Metadata extracted:\n")
cat("  NOTE: Each row is a pseudobulk SAMPLE (not individual cells)\n")
cat("  Rows (pseudobulk samples):", nrow(metadata), "\n")
cat("  Columns:", ncol(metadata), "\n\n")

# Check for expected columns in raw Seurat metadata
# NOTE: The Seurat metadata only contains basic columns at this stage.
# Additional columns (group, allosome_type) are added later by merging with overview data.
expected_cols <- c("cell_type", "gonadal_type")
available_cols <- colnames(metadata)

cat("Checking for expected metadata columns in raw Seurat object:\n")
for (col in expected_cols) {
  if (col %in% available_cols) {
    cat("  ✓", col, "found\n")
  } else {
    cat("  ✗", col, "NOT FOUND\n")

    # Try to infer missing columns
    if (col == "cell_type" && "celltype" %in% available_cols) {
      cat("    → Using 'celltype' instead\n")
      metadata$cell_type <- metadata$celltype
    } else if (col == "gonadal_type" && "gonad_type" %in% available_cols) {
      cat("    → Using 'gonad_type' instead\n")
      metadata$gonadal_type <- metadata$gonad_type
    }
  }
}

# Check for sample identifier column (usually the rownames or first column)
if (!("sample" %in% available_cols)) {
  cat("  ✗ sample column NOT FOUND\n")
  cat("    → Creating 'sample' column from rownames\n")
  metadata$sample <- rownames(metadata)
}

# Note about additional columns
cat("\nNOTE: Additional columns (group, allosome_type, etc.) will be added\n")
cat("      by merging with overview data in downstream workflows.\n")

# Display metadata summary
cat("\nPseudobulk sample summary:\n")
cat("  Total pseudobulk samples:", nrow(metadata), "\n\n")

if ("cell_type" %in% colnames(metadata)) {
  cat("  Cell types represented:\n")
  cell_type_counts <- table(metadata$cell_type)
  print(cell_type_counts)
  cat("\n")
}

if ("sample_name" %in% colnames(metadata)) {
  # sample_name contains the biological sample ID (mouse ID)
  # Each mouse has multiple pseudobulk samples (one per cell type)
  n_biological_samples <- length(unique(metadata$sample_name))
  cat("  Unique biological samples (mice):", n_biological_samples, "\n")
  if ("cell_type" %in% colnames(metadata)) {
    n_cell_types <- length(unique(metadata$cell_type))
    cat("  Expected pseudobulk samples:", n_biological_samples, "×", n_cell_types,
        "=", n_biological_samples * n_cell_types, "\n")
    cat("  Actual pseudobulk samples:", nrow(metadata), "\n\n")
  }
}

if ("gonadal_type" %in% colnames(metadata)) {
  cat("  Gonadal types:\n")
  print(table(metadata$gonadal_type))
  cat("\n")
}

# NOTE: The raw Seurat metadata does NOT contain "group" or "allosome_type".
# These are added by merging with overview data in downstream workflows.
# At this stage, we only have: cell_type, gonadal_type, sample_name, and sample.

# Save metadata
cat("\nSaving metadata to:", file.path(data_rna_dir, "meta_data_df.csv"), "\n")
write.csv(metadata, file.path(data_rna_dir, "meta_data_df.csv"), row.names = TRUE)

################################################################################
# 4. QUALITY CONTROL METRICS
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("4. QUALITY CONTROL METRICS\n")
cat(strrep("-", 80), "\n\n")

# Calculate QC metrics for PSEUDOBULK data
cat("Computing quality control metrics for pseudobulk data...\n\n")

# 1. Pseudobulk samples per biological sample (mouse)
if ("sample_name" %in% colnames(metadata)) {
  samples_per_bio <- metadata %>%
    group_by(sample_name) %>%
    summarize(n_pseudobulk = n(), .groups = "drop") %>%
    arrange(desc(n_pseudobulk))

  cat("Pseudobulk samples per biological sample:\n")
  cat("  Total biological samples:", nrow(samples_per_bio), "\n")
  cat("  Mean pseudobulk samples per bio sample:", round(mean(samples_per_bio$n_pseudobulk), 1), "\n")
  cat("  Median:", median(samples_per_bio$n_pseudobulk), "\n")
  cat("  Range:", min(samples_per_bio$n_pseudobulk), "-", max(samples_per_bio$n_pseudobulk), "\n\n")

  # Plot distribution
  p_samples <- ggplot(samples_per_bio, aes(x = reorder(sample_name, n_pseudobulk), y = n_pseudobulk)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Biological Sample (Mouse ID)", y = "Number of Pseudobulk Samples",
         title = "Pseudobulk Samples per Biological Sample\n(usually 1 per cell type = 5 total)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6))

  ggsave(file.path(qc_dir, "pseudobulk_per_biosample.pdf"), p_samples, width = 8, height = 10)
  cat("  Saved: results/qc/pseudobulk_per_biosample.pdf\n")
}

# 2. Pseudobulk samples per cell type
if ("cell_type" %in% colnames(metadata)) {
  samples_per_type <- metadata %>%
    group_by(cell_type) %>%
    summarize(n_pseudobulk = n(), .groups = "drop") %>%
    arrange(desc(n_pseudobulk))

  cat("\nPseudobulk samples per cell type:\n")
  print(samples_per_type)

  # Plot
  p_types <- ggplot(samples_per_type, aes(x = reorder(cell_type, n_pseudobulk), y = n_pseudobulk)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    coord_flip() +
    labs(x = "Cell Type", y = "Number of Pseudobulk Samples",
         title = "Pseudobulk Samples per Cell Type\n(usually = number of biological samples)") +
    theme_bw()

  ggsave(file.path(qc_dir, "pseudobulk_per_celltype.pdf"), p_types, width = 8, height = 6)
  cat("\n  Saved: results/qc/pseudobulk_per_celltype.pdf\n")
}

# 3. Sample × cell type design matrix
if ("sample_name" %in% colnames(metadata) && "cell_type" %in% colnames(metadata)) {
  sample_celltype_design <- metadata %>%
    group_by(sample_name, cell_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = cell_type, values_from = n, values_fill = 0)

  cat("\nBiological sample (mouse) × cell type design matrix (first 10 rows):\n")
  print(head(sample_celltype_design, 10))

  # Save full table
  write.csv(sample_celltype_design,
            file.path(qc_dir, "sample_celltype_design.csv"),
            row.names = FALSE)

  cat("\n  Saved: results/qc/sample_celltype_design.csv\n")
}

# 4. Experimental design balance
if (all(c("group", "cell_type") %in% colnames(metadata))) {
  samples_per_group_type <- metadata %>%
    group_by(group, cell_type) %>%
    summarize(n_pseudobulk = n(), .groups = "drop") %>%
    pivot_wider(names_from = cell_type, values_from = n_pseudobulk, values_fill = 0)

  cat("\nExperimental design (pseudobulk samples per group × cell type):\n")
  print(samples_per_group_type)

  write.csv(samples_per_group_type,
            file.path(qc_dir, "experimental_design.csv"),
            row.names = FALSE)
}

################################################################################
# 5. EXTRACT SCALED EXPRESSION DATA
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("5. EXTRACTING SCALED EXPRESSION DATA\n")
cat(strrep("-", 80), "\n\n")

# Check if scaled data exists (compatible with Seurat v4 and v5)
if (!has_layer_data(seurat_obj, "RNA", "scale.data")) {
  cat("No scaled data found. Scaling data now...\n")
  cat("This may take several minutes...\n\n")

  # Scale all genes
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  cat("Scaling complete.\n\n")
}

# Extract scaled data (compatible with Seurat v4 and v5)
cat("Extracting scaled expression matrix...\n")
scaled_data <- get_layer_data(seurat_obj, "RNA", "scale.data")

cat("  Genes:", nrow(scaled_data), "\n")
cat("  Pseudobulk samples:", ncol(scaled_data), "\n")
cat("  Matrix size:", format(object.size(scaled_data), units = "MB"), "\n\n")

# Convert to data frame
cat("Converting to data frame...\n")
pseudobulk_df <- as.data.frame(scaled_data)

# Clean up column names (sample identifiers)
# Ensure column names match metadata rownames
if (!all(colnames(pseudobulk_df) %in% rownames(metadata))) {
  cat("  WARNING: Some sample names don't match metadata. Attempting to fix...\n")

  # For now, assume they match by position if lengths are equal
  if (ncol(pseudobulk_df) == nrow(metadata)) {
    cat("  Setting column names to match metadata rownames\n")
    colnames(pseudobulk_df) <- rownames(metadata)
  }
}

# Verify match
matching_samples <- sum(colnames(pseudobulk_df) %in% rownames(metadata))
cat("  Samples matching metadata:", matching_samples, "/", ncol(pseudobulk_df), "\n")

################################################################################
# 6. VERIFY PSEUDOBULK DATA STRUCTURE
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("6. VERIFYING PSEUDOBULK DATA STRUCTURE\n")
cat(strrep("-", 80), "\n\n")

cat("NOTE: Data is already in pseudobulk format (no aggregation needed).\n")
cat("Each column represents a biological sample × cell type combination.\n\n")

# The "sample" column already contains the full pseudobulk sample ID
# Format: [MouseID]_[CellType]_[GonadType]
# We can use it directly as the pseudobulk sample identifier
if ("sample" %in% colnames(metadata)) {
  metadata$sample_celltype <- metadata$sample
  cat("Pseudobulk sample naming format: [MouseID]_[CellType]_[GonadType]\n")
  cat("Example:", head(metadata$sample_celltype, 3), "\n\n")
}

cat("Pseudobulk data summary:\n")
cat("  Genes:", nrow(pseudobulk_df), "\n")
cat("  Pseudobulk samples:", ncol(pseudobulk_df), "\n")

if ("cell_type" %in% colnames(metadata)) {
  n_cell_types <- length(unique(metadata$cell_type))
  cat("  Cell types:", n_cell_types, "\n")
  if ("sample_name" %in% colnames(metadata)) {
    n_bio_samples <- length(unique(metadata$sample_name))
    cat("  Biological samples (mice):", n_bio_samples, "\n")
    cat("  Expected total:", n_bio_samples, "×", n_cell_types, "=", n_bio_samples * n_cell_types, "\n")
  }
}

################################################################################
# 7. SPLIT AND SAVE DATA FILES
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("7. SAVING DATA FILES\n")
cat(strrep("-", 80), "\n\n")

cat("Splitting data into 3 files (GitHub size constraint)...\n\n")

# IMPORTANT: Split by SAMPLES (columns), NOT by genes (rows)
# Each file should have ALL genes but different samples
# This allows proper merging in load_rna_data() function

n_samples <- ncol(pseudobulk_df)

# Add gene column before splitting
pseudobulk_df$gene <- rownames(pseudobulk_df)

# Reorder to put gene column first
pseudobulk_df <- pseudobulk_df[, c("gene", setdiff(colnames(pseudobulk_df), "gene"))]

# Split samples into 3 parts (keeping gene column in each)
# File 1: gene + samples 1-79
# File 2: gene + samples 80-159
# File 3: gene + samples 160-215

sample_cols_1 <- 2:80  # Columns 2-80 (after gene column)
sample_cols_2 <- 81:160  # Columns 81-160
sample_cols_3 <- 161:ncol(pseudobulk_df)  # Remaining columns

scaled_RNA_1 <- pseudobulk_df[, c(1, sample_cols_1), drop = FALSE]  # gene + samples 1-79
scaled_RNA_2 <- pseudobulk_df[, c(1, sample_cols_2), drop = FALSE]  # gene + samples 80-159
scaled_RNA_3 <- pseudobulk_df[, c(1, sample_cols_3), drop = FALSE]  # gene + remaining samples

cat("File 1:", nrow(scaled_RNA_1), "genes x", ncol(scaled_RNA_1) - 1, "samples\n")
cat("File 2:", nrow(scaled_RNA_2), "genes x", ncol(scaled_RNA_2) - 1, "samples\n")
cat("File 3:", nrow(scaled_RNA_3), "genes x", ncol(scaled_RNA_3) - 1, "samples\n\n")

# Save files
cat("Saving data files...\n")

write.csv(scaled_RNA_1,
          file.path(data_rna_dir, "scaled_RNA_1.csv"),
          row.names = FALSE)
cat("  Saved: data/RNA/scaled_RNA_1.csv (",
    format(file.size(file.path(data_rna_dir, "scaled_RNA_1.csv")) / 1024^2, digits = 2),
    "MB)\n")

write.csv(scaled_RNA_2,
          file.path(data_rna_dir, "scaled_RNA_2.csv"),
          row.names = FALSE)
cat("  Saved: data/RNA/scaled_RNA_2.csv (",
    format(file.size(file.path(data_rna_dir, "scaled_RNA_2.csv")) / 1024^2, digits = 2),
    "MB)\n")

write.csv(scaled_RNA_3,
          file.path(data_rna_dir, "scaled_RNA_3.csv"),
          row.names = FALSE)
cat("  Saved: data/RNA/scaled_RNA_3.csv (",
    format(file.size(file.path(data_rna_dir, "scaled_RNA_3.csv")) / 1024^2, digits = 2),
    "MB)\n")

################################################################################
# 8. INITIAL EXPLORATORY ANALYSIS
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("8. EXPLORATORY ANALYSIS\n")
cat(strrep("-", 80), "\n\n")

# Create overall PCA of ALL samples colored by cell type
# This shows that cell type drives clustering more than biological sample of origin
cat("Creating overall PCA of all samples colored by cell type...\n\n")

# Load data EXACTLY as in pseudobulk_DA.R to ensure same gene set and format
cat("  Loading scaled RNA data (matching pseudobulk_DA.R approach)...\n")

# Load the three parts and merge (matching pseudobulk_DA.R lines 24-32)
scaled_RNA_1_pca <- as.data.frame(read_csv("data/RNA/scaled_RNA_1.csv"))
scaled_RNA_2_pca <- as.data.frame(read_csv("data/RNA/scaled_RNA_2.csv"))
scaled_RNA_3_pca <- as.data.frame(read_csv("data/RNA/scaled_RNA_3.csv"))

scaled_RNA_pca <- merge(scaled_RNA_1_pca, scaled_RNA_2_pca, by = "gene")
scaled_RNA_pca <- merge(scaled_RNA_pca, scaled_RNA_3_pca, by = "gene")

row.names(scaled_RNA_pca) <- scaled_RNA_pca$gene
scaled_RNA_pca <- scaled_RNA_pca[, -1]

rm(scaled_RNA_1_pca)
rm(scaled_RNA_2_pca)
rm(scaled_RNA_3_pca)

cat("  Loaded", nrow(scaled_RNA_pca), "genes x", ncol(scaled_RNA_pca), "samples\n")

# Load metadata and merge with overview data (matching pseudobulk_DA.R lines 38-40)
meta_data_df_pca <- as.data.frame(read_csv("data/RNA/meta_data_df.csv"))
Suppl_Table_1 <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
Suppl_Table_1 <- Suppl_Table_1[, -4]
meta_data_df_pca$sample_name <- gsub("-","_",meta_data_df_pca$sample_name)
meta_data_df_pca <- merge(meta_data_df_pca, Suppl_Table_1, by.x = "sample_name", by.y = "sequencing_name")

cat("  Loaded metadata for", nrow(meta_data_df_pca), "samples\n\n")

# Define cell type color mapping (matching pseudobulk_DA.R lines 42-46)
color_mapping <- c("Hepatocytes" = "#ca3659",
                   "Endothelial cells" = "#8f7b2f",
                   "Immune cells" = "#0d934d",
                   "Stellate cells" = "#1892b1",
                   "Cholangiocytes" = "#834f9f")

# Create PCA plot EXACTLY as in pseudobulk_DA.R (lines 48-50)
cat("  Generating PCA plot...\n")

pdf(file.path(qc_dir, "PCA_all_samples_by_celltype.pdf"),
    width = 10, height = 10)
grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
grid::grid.draw(nicePCA(scaled_RNA_pca, meta_data_df_pca[, c(2, 4)],
                        no_label = TRUE, color_mapping = color_mapping))
dev.off()

cat("  Saved: results/qc/PCA_all_samples_by_celltype.pdf\n")
cat("  This plot shows that cell type is the primary driver of variation,\n")
cat("  with samples clustering by cell type rather than by biological sample.\n\n")

# Clean up temporary objects
rm(scaled_RNA_pca)
rm(meta_data_df_pca)
rm(Suppl_Table_1)

# Create overall PCA for each cell type
cat("Creating individual PCA plots for each cell type...\n\n")

cell_types <- unique(metadata$cell_type)

for (cell_type in cell_types) {
  cat("  Processing", cell_type, "...\n")

  # Get samples for this cell type
  cell_type_samples <- metadata %>%
    filter(cell_type == !!cell_type) %>%
    pull(sample_celltype) %>%
    unique()

  # Filter to samples present in pseudobulk data
  cell_type_samples <- intersect(cell_type_samples, colnames(pseudobulk_df))

  if (length(cell_type_samples) < 5) {
    cat("    Skipping (too few samples)\n")
    next
  }

  # Get data for this cell type
  cell_type_data <- pseudobulk_df[, cell_type_samples, drop = FALSE]

  # Remove genes with zero variance
  gene_vars <- apply(cell_type_data, 1, var, na.rm = TRUE)
  cell_type_data <- cell_type_data[gene_vars > 0, , drop = FALSE]

  # Create targets dataframe for nicePCA
  # nicePCA expects: sample column matching data colnames + condition column(s)
  cell_type_meta <- metadata %>%
    filter(sample_celltype %in% cell_type_samples)

  # Build targets with correct structure for nicePCA
  cell_type_targets <- data.frame(
    sample = cell_type_meta$sample_celltype,
    stringsAsFactors = FALSE
  )

  # Add condition column (gonadal_type if available)
  if ("gonadal_type" %in% colnames(metadata)) {
    cell_type_targets$condition = cell_type_meta$gonadal_type
  }

  # Set rownames to match sample IDs
  rownames(cell_type_targets) <- cell_type_targets$sample

  # Align data and targets
  common_samples <- intersect(colnames(cell_type_data), rownames(cell_type_targets))
  cell_type_data <- cell_type_data[, common_samples, drop = FALSE]
  cell_type_targets <- cell_type_targets[common_samples, , drop = FALSE]

  # Create color mapping for gonadal type
  color_mapping <- setNames(c("#e54141", "#4173d8"),
                            c("Ovaries", "Testes"))

  # Create PCA plot
  if (ncol(cell_type_data) >= 3 && nrow(cell_type_data) >= 3) {
    pca_plot <- nicePCA(
      df = cell_type_data,
      targets = cell_type_targets,
      components = c(1, 2, 3),
      centering = TRUE,
      scaling = FALSE,
      pointSize = 3,
      color_mapping = color_mapping
    )

    # Save
    pdf(file.path(qc_dir, paste0("PCA_", gsub(" ", "_", cell_type), ".pdf")),
        width = 10, height = 10)
    grid::grid.draw(pca_plot)
    dev.off()

    cat("    Saved: results/qc/PCA_", gsub(" ", "_", cell_type), ".pdf\n", sep = "")
  }
}

cat("\nExploratory analysis complete.\n")

################################################################################
# 9. DATA VALIDATION
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("9. DATA VALIDATION\n")
cat(strrep("-", 80), "\n\n")

cat("Performing data validation checks...\n\n")

# NOTE: pseudobulk_df now has a 'gene' column, so we need to select only numeric columns
# for validation operations that require numeric data
numeric_cols <- sapply(pseudobulk_df, is.numeric)
pseudobulk_numeric <- pseudobulk_df[, numeric_cols, drop = FALSE]

validation_results <- list()

# 1. Check for missing values
pct_missing <- sum(is.na(pseudobulk_numeric)) / (nrow(pseudobulk_numeric) * ncol(pseudobulk_numeric)) * 100
validation_results$missing_values <- pct_missing
cat("  Missing values:", round(pct_missing, 3), "%",
    ifelse(pct_missing < 1, "✓", "⚠"), "\n")

# 2. Check for infinite values
n_infinite <- sum(!is.finite(as.matrix(pseudobulk_numeric)))
validation_results$infinite_values <- n_infinite
cat("  Infinite values:", n_infinite, ifelse(n_infinite == 0, "✓", "⚠"), "\n")

# 3. Check gene variance
gene_variances <- apply(pseudobulk_numeric, 1, var, na.rm = TRUE)
n_zero_var <- sum(gene_variances == 0, na.rm = TRUE)
validation_results$zero_variance_genes <- n_zero_var
cat("  Zero-variance genes:", n_zero_var, "/", nrow(pseudobulk_numeric),
    paste0("(", round(100 * n_zero_var / nrow(pseudobulk_numeric), 2), "%)"), "\n")

# 4. Check sample completeness
samples_complete <- colSums(!is.na(pseudobulk_numeric)) == nrow(pseudobulk_numeric)
n_complete_samples <- sum(samples_complete)
validation_results$complete_samples <- n_complete_samples
cat("  Complete samples (no NAs):", n_complete_samples, "/", ncol(pseudobulk_numeric), "\n")

# 5. Check data range
data_range <- range(pseudobulk_numeric, na.rm = TRUE)
validation_results$data_range <- data_range
cat("  Data range:", round(data_range[1], 3), "to", round(data_range[2], 3), "\n")

# 6. Check for duplicated genes
n_dup_genes <- sum(duplicated(rownames(pseudobulk_df)))
validation_results$duplicated_genes <- n_dup_genes
cat("  Duplicated gene names:", n_dup_genes, ifelse(n_dup_genes == 0, "✓", "⚠"), "\n")

# 7. Check for duplicated samples
n_dup_samples <- sum(duplicated(colnames(pseudobulk_df)))
validation_results$duplicated_samples <- n_dup_samples
cat("  Duplicated sample names:", n_dup_samples, ifelse(n_dup_samples == 0, "✓", "⚠"), "\n")

cat("\n")

# Overall validation status
all_checks_pass <- (pct_missing < 5) && (n_infinite == 0) &&
                   (n_dup_genes == 0) && (n_dup_samples == 0)

if (all_checks_pass) {
  cat("  ✓✓✓ ALL VALIDATION CHECKS PASSED ✓✓✓\n")
} else {
  cat("  ⚠ SOME VALIDATION CHECKS FAILED - REVIEW RESULTS ⚠\n")
}

################################################################################
# 10. GENERATE QC REPORT
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("10. GENERATING QC REPORT\n")
cat(strrep("-", 80), "\n\n")

# Create comprehensive report
report_file <- file.path(qc_dir, "data_quality_report.txt")

sink(report_file)

cat(strrep("=", 80), "\n")
cat("DATA QUALITY CONTROL REPORT\n")
cat(strrep("=", 80), "\n\n")

cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Workflow: 01_data_preparation_qc.R\n\n")

cat("DATASET OVERVIEW\n")
cat(strrep("-", 80), "\n")
cat("NOTE: This is PSEUDOBULK data (sample × cell type combinations)\n")
cat("Total pseudobulk samples:", nrow(metadata), "\n")
cat("Total genes:", nrow(pseudobulk_df), "\n")
if ("sample_name" %in% colnames(metadata)) {
  cat("Biological samples (mice):", length(unique(metadata$sample_name)), "\n")
}
cat("Cell types:", length(unique(metadata$cell_type)), "\n")
if ("group" %in% colnames(metadata)) {
  cat("Experimental groups:", length(unique(metadata$group)), "\n")
} else {
  cat("Experimental groups: (not available in raw metadata)\n")
}
cat("\n")

cat("PSEUDOBULK SAMPLES BY CELL TYPE\n")
cat(strrep("-", 80), "\n")
if ("cell_type" %in% colnames(metadata)) {
  ct_summary <- metadata %>%
    group_by(cell_type) %>%
    summarize(n_pseudobulk = n(), pct = 100 * n() / nrow(metadata), .groups = "drop") %>%
    arrange(desc(n_pseudobulk))
  print(ct_summary)
}

cat("\n")
cat("PSEUDOBULK SAMPLES BY EXPERIMENTAL GROUP\n")
cat(strrep("-", 80), "\n")
if ("group" %in% colnames(metadata)) {
  group_summary <- metadata %>%
    group_by(group) %>%
    summarize(n_pseudobulk = n(), pct = 100 * n() / nrow(metadata), .groups = "drop") %>%
    arrange(desc(n_pseudobulk))
  print(group_summary)
}

cat("\n")
cat("VALIDATION RESULTS\n")
cat(strrep("-", 80), "\n")
cat("Missing values:        ", round(validation_results$missing_values, 3), "%\n")
cat("Infinite values:       ", validation_results$infinite_values, "\n")
cat("Zero-variance genes:   ", validation_results$zero_variance_genes, "\n")
cat("Complete samples:      ", validation_results$complete_samples, "/", ncol(pseudobulk_df), "\n")
cat("Data range:            ", round(validation_results$data_range[1], 3), "to",
    round(validation_results$data_range[2], 3), "\n")
cat("Duplicated genes:      ", validation_results$duplicated_genes, "\n")
cat("Duplicated samples:    ", validation_results$duplicated_samples, "\n")

cat("\n")
cat("FILES GENERATED\n")
cat(strrep("-", 80), "\n")
cat("data/RNA/scaled_RNA_1.csv\n")
cat("data/RNA/scaled_RNA_2.csv\n")
cat("data/RNA/scaled_RNA_3.csv\n")
cat("data/RNA/meta_data_df.csv\n")
cat("results/qc/pseudobulk_per_biosample.pdf\n")
cat("results/qc/pseudobulk_per_celltype.pdf\n")
cat("results/qc/sample_celltype_design.csv\n")
cat("results/qc/experimental_design.csv\n")
cat("results/qc/PCA_all_samples_by_celltype.pdf\n")
cat("results/qc/PCA_[celltype].pdf (for each cell type)\n")

cat("\n")
cat(strrep("=", 80), "\n")
cat("END OF REPORT\n")
cat(strrep("=", 80), "\n")

sink()

cat("QC report saved to:", report_file, "\n")

################################################################################
# WORKFLOW COMPLETE
################################################################################

################################################################################
# 11. ZONATION PSEUDOBULK: SEURAT TO CSV CONVERSION
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("11. ZONATION PSEUDOBULK: SEURAT TO CSV CONVERSION\n")
cat(strrep("-", 80), "\n\n")

cat("Converting zonation-annotated hepatocyte pseudobulk from Seurat to CSV.\n")
cat("These CSVs are consumed by workflow 05 (zonation analysis).\n\n")

zonation_rds_path <- "data/RNA/fcg_zonation_pseudobulks.rds"

if (!file.exists(zonation_rds_path)) {
  cat("WARNING: Zonation Seurat object not found at:", zonation_rds_path, "\n")
  cat("  Skipping zonation CSV generation.\n")
  cat("  Workflow 05 will use existing CSVs if available.\n\n")
} else {
  cat("Loading zonation Seurat object:", zonation_rds_path, "\n")
  zonation_obj <- readRDS(zonation_rds_path)

  cat("  Pseudobulk samples:", ncol(zonation_obj), "\n")
  cat("  Genes:", nrow(zonation_obj), "\n\n")

  # Extract scaled data (compatible with Seurat v4 and v5)
  cat("Extracting scaled expression data...\n")
  zonation_scaled <- as.data.frame(
    get_layer_data(zonation_obj, assay = "RNA", layer = "scale.data")
  )
  rownames(zonation_scaled) <- Features(zonation_obj[["RNA"]])

  # Extract metadata
  zonation_meta <- as.data.frame(zonation_obj@meta.data)
  names(zonation_meta)[1] <- "sample"

  # Normalize sample name separators (hyphens -> underscores)
  zonation_meta$sample <- gsub("-", "_", zonation_meta$sample)
  if ("sample_name" %in% colnames(zonation_meta)) {
    zonation_meta$sample_name <- gsub("-", "_", zonation_meta$sample_name)
  }

  # Set column names to match metadata sample names
  colnames(zonation_scaled) <- zonation_meta$sample

  cat("  Scaled data:", nrow(zonation_scaled), "genes x", ncol(zonation_scaled), "samples\n")
  cat("  Metadata:", nrow(zonation_meta), "rows x", ncol(zonation_meta), "columns\n")
  cat("  Metadata columns:", paste(colnames(zonation_meta), collapse = ", "), "\n\n")

  # Add gene column and reorder (gene first, then samples)
  zonation_scaled$gene <- rownames(zonation_scaled)
  zonation_scaled <- zonation_scaled[, c("gene", setdiff(colnames(zonation_scaled), "gene"))]

  # Save CSVs
  cat("Saving zonation CSVs...\n")

  readr::write_csv(zonation_scaled, file.path(data_rna_dir, "scaled_RNA_zonation.csv"))
  cat("  Saved:", file.path(data_rna_dir, "scaled_RNA_zonation.csv"),
      "(", format(file.size(file.path(data_rna_dir, "scaled_RNA_zonation.csv")) / 1024^2, digits = 2), "MB)\n")

  rownames(zonation_meta) <- NULL
  readr::write_csv(zonation_meta, file.path(data_rna_dir, "meta_data_df_zonation.csv"))
  cat("  Saved:", file.path(data_rna_dir, "meta_data_df_zonation.csv"), "\n\n")

  # Verify output
  cat("Verification:\n")
  cat("  Sample names (first 4):", paste(head(zonation_meta$sample, 4), collapse = ", "), "\n")
  if ("zonation_type" %in% colnames(zonation_meta)) {
    cat("  Zonation types:", paste(unique(zonation_meta$zonation_type), collapse = ", "), "\n")
  }
  if ("gonadal_type" %in% colnames(zonation_meta)) {
    cat("  Gonadal types:", paste(unique(zonation_meta$gonadal_type), collapse = ", "), "\n")
  }
  cat("\n")

  # Clean up
  rm(zonation_obj, zonation_scaled, zonation_meta)
  gc(verbose = FALSE)
}

################################################################################
# 12. GDX PSEUDOBULK: SEURAT TO CSV CONVERSION
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("12. GDX PSEUDOBULK: SEURAT TO CSV CONVERSION\n")
cat(strrep("-", 80), "\n\n")

cat("Converting GDX (gonadectomy) pseudobulk from Seurat to split CSVs.\n")
cat("These CSVs are consumed by workflow 02 (step 3B) and workflow 03.\n\n")

gdx_rds_path <- "data/RNA/gdx_pseudobulks.rds"

if (!file.exists(gdx_rds_path)) {
  cat("WARNING: GDX Seurat object not found at:", gdx_rds_path, "\n")
  cat("  Skipping GDX CSV generation.\n")
  cat("  Workflows 02/03 will use existing CSVs if available.\n\n")
} else {
  cat("Loading GDX Seurat object:", gdx_rds_path, "\n")
  gdx_obj <- readRDS(gdx_rds_path)

  cat("  Pseudobulk samples:", ncol(gdx_obj), "\n")
  cat("  Genes:", nrow(gdx_obj), "\n\n")

  # Extract scaled data (compatible with Seurat v4 and v5)
  cat("Extracting scaled expression data...\n")
  gdx_scaled <- as.data.frame(
    get_layer_data(gdx_obj, assay = "RNA", layer = "scale.data")
  )
  rownames(gdx_scaled) <- Features(gdx_obj[["RNA"]])

  # Extract metadata
  gdx_meta <- as.data.frame(gdx_obj@meta.data)
  names(gdx_meta)[1] <- "sample"

  # Normalize sample name separators (hyphens -> underscores)
  gdx_meta$sample <- gsub("-", "_", gdx_meta$sample)
  if ("sample_name" %in% colnames(gdx_meta)) {
    gdx_meta$sample_name <- gsub("-", "_", gdx_meta$sample_name)
  }

  # Set column names to match metadata sample names
  colnames(gdx_scaled) <- gdx_meta$sample

  cat("  Scaled data:", nrow(gdx_scaled), "genes x", ncol(gdx_scaled), "samples\n")
  cat("  Metadata:", nrow(gdx_meta), "rows x", ncol(gdx_meta), "columns\n")
  cat("  Metadata columns:", paste(colnames(gdx_meta), collapse = ", "), "\n\n")

  # Add gene column and reorder
  gdx_scaled$gene <- rownames(gdx_scaled)
  gdx_scaled <- gdx_scaled[, c("gene", setdiff(colnames(gdx_scaled), "gene"))]

  # Split into 3 files (same pattern as FCG, for GitHub size constraint)
  n_gdx_samples <- ncol(gdx_scaled) - 1  # exclude gene column
  split_size <- ceiling(n_gdx_samples / 3)
  sample_cols_1 <- 2:(1 + split_size)
  sample_cols_2 <- (2 + split_size):(1 + 2 * split_size)
  sample_cols_3 <- (2 + 2 * split_size):ncol(gdx_scaled)

  gdx_rna_1 <- gdx_scaled[, c(1, sample_cols_1), drop = FALSE]
  gdx_rna_2 <- gdx_scaled[, c(1, sample_cols_2), drop = FALSE]
  gdx_rna_3 <- gdx_scaled[, c(1, sample_cols_3), drop = FALSE]

  cat("Splitting GDX data into 3 files...\n")
  cat("  File 1:", nrow(gdx_rna_1), "genes x", ncol(gdx_rna_1) - 1, "samples\n")
  cat("  File 2:", nrow(gdx_rna_2), "genes x", ncol(gdx_rna_2) - 1, "samples\n")
  cat("  File 3:", nrow(gdx_rna_3), "genes x", ncol(gdx_rna_3) - 1, "samples\n\n")

  # Save split CSVs
  cat("Saving GDX CSVs...\n")

  readr::write_csv(gdx_rna_1, file.path(data_rna_dir, "scaled_RNA_gdxhr_1.csv"))
  cat("  Saved:", file.path(data_rna_dir, "scaled_RNA_gdxhr_1.csv"),
      "(", format(file.size(file.path(data_rna_dir, "scaled_RNA_gdxhr_1.csv")) / 1024^2, digits = 2), "MB)\n")

  readr::write_csv(gdx_rna_2, file.path(data_rna_dir, "scaled_RNA_gdxhr_2.csv"))
  cat("  Saved:", file.path(data_rna_dir, "scaled_RNA_gdxhr_2.csv"),
      "(", format(file.size(file.path(data_rna_dir, "scaled_RNA_gdxhr_2.csv")) / 1024^2, digits = 2), "MB)\n")

  readr::write_csv(gdx_rna_3, file.path(data_rna_dir, "scaled_RNA_gdxhr_3.csv"))
  cat("  Saved:", file.path(data_rna_dir, "scaled_RNA_gdxhr_3.csv"),
      "(", format(file.size(file.path(data_rna_dir, "scaled_RNA_gdxhr_3.csv")) / 1024^2, digits = 2), "MB)\n")

  # Save metadata
  rownames(gdx_meta) <- NULL
  readr::write_csv(gdx_meta, file.path(data_rna_dir, "meta_data_df_gdxhr.csv"))
  cat("  Saved:", file.path(data_rna_dir, "meta_data_df_gdxhr.csv"), "\n\n")

  # Clean up split objects
  rm(gdx_rna_1, gdx_rna_2, gdx_rna_3)

  # Verify output
  cat("Verification:\n")
  cat("  Sample names (first 4):", paste(head(gdx_meta$sample, 4), collapse = ", "), "\n")
  if ("cell_type" %in% colnames(gdx_meta)) {
    cat("  Cell types:", paste(sort(unique(gdx_meta$cell_type)), collapse = ", "), "\n")
  }
  if ("group" %in% colnames(gdx_meta)) {
    cat("  Groups:", paste(sort(unique(gdx_meta$group)), collapse = ", "), "\n")
  }
  cat("\n")

  # Clean up
  rm(gdx_obj, gdx_scaled, gdx_meta)
  gc(verbose = FALSE)
}

################################################################################
# WORKFLOW COMPLETE
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("DATA PREPARATION AND QC WORKFLOW COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("Summary:\n")
cat("  ✓ Pseudobulk Seurat object loaded and processed\n")
cat("  ✓ Metadata extracted and saved\n")
cat("  ✓ Pseudobulk data verified (already aggregated)\n")
cat("  ✓ Data split into 3 files for GitHub\n")
cat("  ✓ Quality control plots generated\n")
cat("  ✓ Data validation completed\n")
cat("  ✓ QC report generated\n")
cat("  ✓ Zonation pseudobulk CSVs generated\n")
cat("  ✓ GDX pseudobulk CSVs generated\n\n")

cat("Output files:\n")
cat("  - data/RNA/scaled_RNA_1.csv\n")
cat("  - data/RNA/scaled_RNA_2.csv\n")
cat("  - data/RNA/scaled_RNA_3.csv\n")
cat("  - data/RNA/meta_data_df.csv\n")
cat("  - data/RNA/scaled_RNA_zonation.csv\n")
cat("  - data/RNA/meta_data_df_zonation.csv\n")
cat("  - data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv\n")
cat("  - data/RNA/meta_data_df_gdxhr.csv\n")
cat("  - results/qc/ (QC plots and reports)\n\n")

cat("Next steps:\n")
cat("  Run workflow 02: Pseudobulk differential analysis\n")
cat("  source('workflows/02_fcg_pseudobulk_analysis.R')\n\n")

cat("This workflow provides the foundation data for all subsequent analyses.\n")

cat("\n")
cat(strrep("=", 80), "\n\n")
