#' Analysis Configuration
#'
#' This file defines general analysis settings used across all workflows.
#' It provides default parameters for statistical tests, thresholds, and
#' general analysis options.
#'
#' @description
#' Centralized configuration for:
#' - Statistical significance thresholds
#' - Default analysis parameters
#' - Output directory structure
#' - General settings used across multiple analyses

################################################################################
# STATISTICAL THRESHOLDS
################################################################################

# Differential expression/abundance significance thresholds
SIGNIFICANCE_THRESHOLDS <- list(
  # P-value and FDR thresholds
  p_value = 0.05,
  adj_p_value = 0.05,  # FDR threshold

  # Log fold change thresholds
  logFC = 1,  # Default for volcano plots

  # Volcano plot asymptotic curve parameters
  volcano_hAss = 0.05,   # Horizontal asymptote (p-value)
  volcano_vAss = 1       # Vertical asymptote (logFC)
)

################################################################################
# ENRICHMENT ANALYSIS SETTINGS
################################################################################

ENRICHMENT_SETTINGS <- list(
  # Minimum gene set size for pathway enrichment
  min_geneset_size = 10,
  max_geneset_size = 500,

  # decoupleR settings
  method = "ulm",  # Univariate linear model
  minsize = 15,    # Minimum overlap for enrichment

  # Enrichment score thresholds
  score_threshold = 11  # For heatmap filtering
)

################################################################################
# PCA SETTINGS
################################################################################

PCA_SETTINGS <- list(
  # Default components to display
  components = c(1, 2, 3),

  # Scaling and centering
  center = TRUE,
  scale = FALSE,  # Don't scale for already-scaled data

  # Visualization
  point_size = 4,

  # Variance explained threshold for reporting
  var_explained_threshold = 0.05  # Report PCs explaining >5% variance
)

################################################################################
# LIMMA DIFFERENTIAL ANALYSIS SETTINGS
################################################################################

LIMMA_SETTINGS <- list(
  # Trend for eBayes
  trend = FALSE,

  # Robust empirical Bayes
  robust = FALSE,

  # Top table settings
  number = Inf,  # Return all genes
  adjust.method = "BH",  # Benjamini-Hochberg FDR

  # Filtering
  min_samples = 3  # Minimum samples per group
)

################################################################################
# MOFA SETTINGS
################################################################################

MOFA_SETTINGS <- list(
  # Factor range to test
  min_factors = 5,
  max_factors = 15,
  optimal_factors = 12,  # Based on model selection

  # Training parameters
  max_iter = 1000,
  convergence_mode = "fast",

  # Feature selection
  top_percent = 0.5,  # Top 50% of features by |weight|
  top_percent_network = 0.025,  # Top 2.5% for network analysis

  # Correlation thresholds
  cor_threshold = 0.65  # For gene-lipid networks
)

################################################################################
# STRUCTURAL ANALYSIS SETTINGS
################################################################################

STRUCTURAL_SETTINGS <- list(
  # Regression parameters
  min_lipids_per_class = 10,  # Minimum lipids for class-specific regression

  # Significance threshold
  p_threshold = 0.05,

  # Model formula components
  predictors = c("lipid_class", "logfc", "chain_length", "db_position"),
  interactions = c("logfc:chain_length", "logfc:db_position")
)

################################################################################
# ZONATION ANALYSIS SETTINGS
################################################################################

ZONATION_SETTINGS <- list(
  # Zonation markers (Halpern et al. 2017)
  pericentral_markers = c("Cyp2e1", "Glul"),
  periportal_markers = c("Cyp2f2", "Gls2", "Hal", "Sds"),

  # Zone assignment threshold
  score_threshold = 0.5,  # Difference in zonation scores for assignment

  # Minimum cells per zone per sample
  min_cells_per_zone = 10
)

################################################################################
# VISUALIZATION SETTINGS
################################################################################

VISUALIZATION_SETTINGS <- list(
  # Default colors
  colors = list(
    gonadal = c(Ovaries = "red", Testes = "royalblue"),
    allosome = c(XX = "red", XY = "royalblue"),
    cell_types = c(
      Hepatocytes = "#ca3659",
      `Endothelial cells` = "#8f7b2f",
      `Immune cells` = "#0d934d",
      `Stellate cells` = "#1892b1",
      Cholangiocytes = "#834f9f"
    ),
    heatmap_diverging = c("red", "white", "royalblue")
  ),

  # Figure dimensions (inches)
  figure_width = 8,
  figure_height = 7,

  # Font sizes
  title_size = 14,
  axis_title_size = 12,
  axis_text_size = 10,
  legend_text_size = 10,

  # Volcano plot settings
  volcano_nlabels = 20,  # Number of genes to label
  volcano_label_size = 3,

  # Heatmap settings
  heatmap_cluster_rows = TRUE,
  heatmap_cluster_cols = TRUE,
  heatmap_show_rownames = FALSE
)

################################################################################
# OUTPUT DIRECTORY STRUCTURE
################################################################################

OUTPUT_DIRS <- list(
  results = "results",
  qc = "results/qc",
  pseudobulk = "results/pseudobulk",
  lipidomic = "results/lipidomic",
  mofa = "results/mofa",
  zonation = "results/zonation",
  structural = "results/multi-omic/structural_model",
  multiomics = "results/multi-omic",
  figures = "results/figures"
)

################################################################################
# FILE NAME PATTERNS
################################################################################

FILE_PATTERNS <- list(
  # Results files
  ttop_results = "ttop_full_list.RData",
  tval_matrix = "tval_full_df.RData",

  # Plots
  volcano_prefix = "volcano_",
  pca_prefix = "PCA_",
  heatmap_prefix = "heatmap_",

  # Data files
  rna_prefix = "scaled_RNA_",
  metadata_file = "meta_data_df.csv"
)

################################################################################
# SESSION SETTINGS
################################################################################

SESSION_SETTINGS <- list(
  # Random seed for reproducibility
  seed = 42,

  # Parallel processing
  use_parallel = FALSE,  # Enable for large datasets
  n_cores = 4,

  # Memory management
  gc_frequency = "high",  # Run garbage collection frequently

  # Progress reporting
  verbose = TRUE
)

################################################################################
# DATA QUALITY THRESHOLDS
################################################################################

QC_THRESHOLDS <- list(
  # Missing data
  max_missing_percent = 5,  # Maximum % missing values allowed

  # Gene filtering
  min_detection_rate = 0.1,  # Minimum % samples with detected expression

  # Sample filtering
  min_genes_detected = 1000,  # Minimum genes detected per sample

  # Variance filtering
  min_variance = 0  # Minimum variance for inclusion
)

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Get output directory path
#'
#' @param dir_name Name of directory from OUTPUT_DIRS
#' @return Full path to directory
get_output_dir <- function(dir_name) {
  if (dir_name %in% names(OUTPUT_DIRS)) {
    return(OUTPUT_DIRS[[dir_name]])
  } else {
    stop("Unknown output directory: ", dir_name)
  }
}

#' Get default colors for a variable
#'
#' @param var_name Variable name (e.g., "gonadal", "cell_types")
#' @return Named vector of colors
get_default_colors <- function(var_name) {
  if (var_name %in% names(VISUALIZATION_SETTINGS$colors)) {
    return(VISUALIZATION_SETTINGS$colors[[var_name]])
  } else {
    warning("No default colors defined for: ", var_name)
    return(NULL)
  }
}

#' Set random seed for reproducibility
set_analysis_seed <- function() {
  set.seed(SESSION_SETTINGS$seed)
  cat("Random seed set to:", SESSION_SETTINGS$seed, "\n")
}

################################################################################
# INITIALIZATION
################################################################################

# Print configuration loaded message
if (SESSION_SETTINGS$verbose) {
  cat("Analysis configuration loaded successfully.\n")
  cat("  - Statistical thresholds: adj.P.Val <", SIGNIFICANCE_THRESHOLDS$adj_p_value, "\n")
  cat("  - MOFA optimal factors:", MOFA_SETTINGS$optimal_factors, "\n")
  cat("  - Random seed:", SESSION_SETTINGS$seed, "\n")
}

################################################################################
# END OF CONFIGURATION
################################################################################
