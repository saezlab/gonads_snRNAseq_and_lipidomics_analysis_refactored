#' Core Enrichment Analysis Functions
#'
#' This module provides functions for pathway and transcription factor enrichment analysis.
#' It wraps decoupleR functions with a cleaner interface and consolidates patterns from:
#' - pseudobulk_pathways.R
#' - pseudobulk_TFs.R
#' - lipidomic_DA.R (lipid class enrichment)
#'
#' All enrichment uses decoupleR::run_ulm() (Univariate Linear Model method).

library(decoupleR)
library(GSEABase)
library(reshape2)
library(dplyr)

# Source support functions for import_gmt
source("scripts/support_functions.R")


#' Run Pathway Enrichment Analysis
#'
#' Performs pathway enrichment on differential analysis results using decoupleR.
#'
#' @param stat_matrix Matrix. Statistics matrix (genes x comparisons). Typically
#'                    t-statistics from limma. Row names should be gene symbols.
#' @param pathway_file Character. Path to GMT file with pathway gene sets.
#'                     Default: "support/m5.all.v2023.2.Mm.symbols.gmt" (MSigDB GO BP)
#' @param minsize Integer. Minimum gene set size. Default: 10
#' @param method Character. decoupleR method to use. Default: "ulm"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame with pathway enrichment scores (long format) or
#'         matrix (wide format) if return_matrix = TRUE.
#'
#' @details
#' The function:
#' 1. Loads pathway gene sets from GMT file
#' 2. Formats as decoupleR network (target, source, mor)
#' 3. Runs enrichment with run_ulm()
#' 4. Returns scores in long or wide format
#'
#' The score represents the pathway activity (positive = pathway up-regulated,
#' negative = pathway down-regulated).
#'
#' @examples
#' # Run on t-statistics matrix
#' pathway_scores <- run_pathway_enrichment(t_matrix)
#'
#' # Use custom GMT file
#' pathway_scores <- run_pathway_enrichment(t_matrix, pathway_file = "my_pathways.gmt")
run_pathway_enrichment <- function(stat_matrix,
                                   pathway_file = "support/m5.all.v2023.2.Mm.symbols.gmt",
                                   minsize = 10,
                                   method = "ulm",
                                   verbose = TRUE) {

  if (verbose) {
    cat("Running pathway enrichment analysis...\n")
    cat("  Pathway file:", pathway_file, "\n")
    cat("  Statistics matrix:", nrow(stat_matrix), "genes x",
        ncol(stat_matrix), "comparisons\n")
  }

  # Load pathway gene sets
  if (verbose) cat("  Loading pathways from GMT file...\n")
  pathways <- import_gmt(pathway_file)
  names(pathways) <- c("target", "source")
  pathways$mor <- 1

  n_pathways <- length(unique(pathways$source))
  if (verbose) cat("  Loaded", n_pathways, "pathways\n")

  # Run enrichment
  if (verbose) cat("  Running decoupleR::run_ulm()...\n")
  pathway_scores <- decoupleR::run_ulm(
    mat = stat_matrix,
    network = pathways,
    minsize = minsize
  )

  if (verbose) {
    n_tested <- length(unique(pathway_scores$source))
    cat("  ✓ Enrichment complete:", n_tested, "pathways tested\n")
  }

  return(pathway_scores)
}


#' Run Transcription Factor Activity Inference
#'
#' Infers transcription factor activity from differential analysis results.
#'
#' @param stat_matrix Matrix. Statistics matrix (genes x comparisons).
#' @param regulon_file Character. Path to RData file with TF regulons.
#'                     Default: "support/TF_regulons.RData" (CollecTRI regulons)
#' @param minsize Integer. Minimum regulon size. Default: 20
#' @param method Character. decoupleR method to use. Default: "ulm"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame with TF activity scores (long format).
#'
#' @details
#' The function:
#' 1. Loads TF regulons (source, target, mor format)
#' 2. Filters to complete cases
#' 3. Runs activity inference with run_ulm()
#' 4. Returns TF activity scores
#'
#' The score represents TF activity (positive = TF more active,
#' negative = TF less active).
#'
#' Regulons can be from:
#' - CollecTRI (default): decoupleR::get_collectri(organism = "mouse")
#' - Custom regulons: Must have columns "source" (TF), "target" (gene), "mor" (mode)
#'
#' @examples
#' # Run TF activity inference
#' tf_scores <- run_tf_activity(t_matrix)
#'
#' # Use custom regulons
#' tf_scores <- run_tf_activity(t_matrix, regulon_file = "my_regulons.RData")
run_tf_activity <- function(stat_matrix,
                             regulon_file = "support/TF_regulons.RData",
                             minsize = 20,
                             method = "ulm",
                             verbose = TRUE) {

  if (verbose) {
    cat("Running TF activity inference...\n")
    cat("  Regulon file:", regulon_file, "\n")
    cat("  Statistics matrix:", nrow(stat_matrix), "genes x",
        ncol(stat_matrix), "comparisons\n")
  }

  # Load TF regulons
  if (verbose) cat("  Loading TF regulons...\n")
  load(regulon_file)  # Should load object named "TF_regulons"
  TF_regulons <- TF_regulons[complete.cases(TF_regulons), ]

  n_tfs <- length(unique(TF_regulons$source))
  if (verbose) cat("  Loaded", n_tfs, "transcription factors\n")

  # Run activity inference
  if (verbose) cat("  Running decoupleR::run_ulm()...\n")
  tf_scores <- decoupleR::run_ulm(
    mat = stat_matrix,
    network = TF_regulons,
    minsize = minsize
  )

  if (verbose) {
    n_tested <- length(unique(tf_scores$source))
    cat("  ✓ Inference complete:", n_tested, "TFs tested\n")
  }

  return(tf_scores)
}


#' Run Lipid Class Enrichment
#'
#' Tests for enrichment of lipid classes in differential lipidomic results.
#'
#' @param stat_matrix Matrix. Statistics matrix (lipids x comparisons).
#' @param extract_class_function Function. Function to extract lipid class from
#'                               lipid names. Default: extracts text before first "(".
#' @param minsize Integer. Minimum class size. Default: 5
#' @param method Character. decoupleR method to use. Default: "ulm"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame with lipid class enrichment scores (long format).
#'
#' @details
#' The function:
#' 1. Extracts lipid classes from lipid names (e.g., "PC(16:0/18:1)" → "PC")
#' 2. Creates a network mapping lipids to their classes
#' 3. Runs enrichment with run_ulm()
#' 4. Returns class-level scores
#'
#' Positive scores indicate lipid classes that are generally more abundant,
#' negative scores indicate classes that are less abundant.
#'
#' @examples
#' # Run lipid class enrichment
#' class_scores <- run_lipid_class_enrichment(lipid_t_matrix)
#'
#' # Custom class extraction (for different nomenclature)
#' class_scores <- run_lipid_class_enrichment(
#'   lipid_t_matrix,
#'   extract_class_function = function(x) gsub("_.*", "", x)
#' )
run_lipid_class_enrichment <- function(stat_matrix,
                                       extract_class_function = NULL,
                                       minsize = 5,
                                       method = "ulm",
                                       verbose = TRUE) {

  if (verbose) {
    cat("Running lipid class enrichment...\n")
    cat("  Statistics matrix:", nrow(stat_matrix), "lipids x",
        ncol(stat_matrix), "comparisons\n")
  }

  # Default class extraction: text before "("
  if (is.null(extract_class_function)) {
    extract_class_function <- function(lipid_name) {
      gsub("\\(.*", "", lipid_name)
    }
  }

  # Extract lipid classes
  lipid_names <- rownames(stat_matrix)
  lipid_classes <- sapply(lipid_names, extract_class_function, USE.NAMES = FALSE)

  # Create network (lipid → class)
  lipid_network <- data.frame(
    target = lipid_names,
    source = lipid_classes,
    mor = 1,
    stringsAsFactors = FALSE
  )

  n_classes <- length(unique(lipid_classes))
  if (verbose) cat("  Detected", n_classes, "lipid classes\n")

  # Run enrichment
  if (verbose) cat("  Running decoupleR::run_ulm()...\n")
  class_scores <- decoupleR::run_ulm(
    mat = stat_matrix,
    network = lipid_network,
    minsize = minsize
  )

  if (verbose) {
    n_tested <- length(unique(class_scores$source))
    cat("  ✓ Enrichment complete:", n_tested, "classes tested\n")
  }

  return(class_scores)
}


#' Convert Enrichment Results to Wide Format
#'
#' Converts long-format decoupleR results to wide matrix format for heatmaps.
#'
#' @param enrichment_scores Data frame. Output from run_pathway_enrichment(),
#'                          run_tf_activity(), or run_lipid_class_enrichment().
#' @param value_column Character. Column to use for values. Default: "score"
#'
#' @return Data frame in wide format (sources x conditions).
#'
#' @examples
#' pathway_scores_wide <- enrichment_to_wide(pathway_scores)
#' tf_scores_wide <- enrichment_to_wide(tf_scores)
enrichment_to_wide <- function(enrichment_scores, value_column = "score") {

  scores_wide <- dcast(
    enrichment_scores,
    source ~ condition,
    value.var = value_column
  )

  return(scores_wide)
}


#' Filter Top Enriched Terms
#'
#' Filters enrichment results to top terms by score threshold.
#'
#' @param enrichment_scores Data frame. Long-format enrichment results.
#' @param score_threshold Numeric. Minimum absolute score to keep.
#'                        Default: 10 for pathways, 4 for TFs, 3 for lipid classes.
#' @param return_wide Logical. Return in wide format? Default: FALSE
#'
#' @return Data frame with top enriched terms (long or wide format).
#'
#' @details
#' A term is kept if its maximum absolute score across all conditions exceeds
#' the threshold. This ensures we get terms that show strong enrichment in at
#' least one comparison.
#'
#' @examples
#' top_pathways <- filter_top_enriched(pathway_scores, score_threshold = 12)
#' top_tfs <- filter_top_enriched(tf_scores, score_threshold = 4)
filter_top_enriched <- function(enrichment_scores,
                                score_threshold = 10,
                                return_wide = FALSE) {

  # Convert to wide format for filtering
  scores_wide <- enrichment_to_wide(enrichment_scores)

  # Calculate max absolute score per source
  score_cols <- setdiff(names(scores_wide), "source")
  max_abs_scores <- apply(
    scores_wide[, score_cols],
    1,
    function(x) max(abs(x), na.rm = TRUE)
  )

  # Filter
  top_indices <- which(max_abs_scores > score_threshold)
  top_scores <- scores_wide[top_indices, ]

  if (return_wide) {
    return(top_scores)
  } else {
    # Convert back to long format
    top_scores_long <- melt(
      top_scores,
      id.vars = "source",
      variable.name = "condition",
      value.name = "score"
    )
    return(top_scores_long)
  }
}


#' Subset Enrichment Scores by Comparison Pattern
#'
#' Subsets enrichment scores to specific comparisons (e.g., only gonadal comparisons).
#'
#' @param scores_wide Data frame. Wide-format enrichment scores.
#' @param pattern Character. Pattern to match in column names (regex).
#' @param keep_source Logical. Keep the "source" column. Default: TRUE
#'
#' @return Data frame with filtered columns.
#'
#' @examples
#' # Get only gonadal comparisons
#' gonadal_scores <- subset_by_pattern(pathway_scores_wide, pattern = "testes")
#'
#' # Get only allosome comparisons
#' allosome_scores <- subset_by_pattern(pathway_scores_wide, pattern = "XY_vs")
subset_by_pattern <- function(scores_wide, pattern, keep_source = TRUE) {

  # Find matching columns
  matching_cols <- grep(pattern, names(scores_wide), value = TRUE)

  if (length(matching_cols) == 0) {
    warning("No columns matched pattern: ", pattern)
    return(NULL)
  }

  # Subset
  if (keep_source && "source" %in% names(scores_wide)) {
    subset_scores <- scores_wide[, c("source", matching_cols)]
  } else {
    subset_scores <- scores_wide[, matching_cols]
  }

  return(subset_scores)
}


#' Build T-statistic Matrix from Differential Analysis Results
#'
#' Combines toptables from multiple cell types and comparisons into a single
#' t-statistic matrix suitable for enrichment analysis.
#'
#' @param results_list List. Results from run_differential_analysis_multi_celltype().
#'                     Format: list(CellType1 = list(Comparison1 = ttop, ...))
#' @param stat_column Character. Column to extract. Default: "t"
#'
#' @return Matrix with genes as rows and "comparison_celltype" as columns.
#'
#' @details
#' This function is the modular equivalent of the code in pseudobulk_pathways.R
#' that creates tval_full_df.
#'
#' Column names follow pattern: "comparison_stat_celltype"
#' (e.g., "gonadal_type_tval_Hepatocytes")
#'
#' @examples
#' # Build from multi-cell-type results
#' t_matrix <- build_stat_matrix_from_results(results_all)
build_stat_matrix_from_results <- function(results_list, stat_column = "t") {

  # Storage
  all_stat_dfs <- list()

  # Loop through grouping types (if results has that structure)
  # For now, assume results_list is cell types directly
  for (cell_type in names(results_list)) {
    comparisons <- results_list[[cell_type]]

    # Extract stat columns from each comparison
    stat_list <- lapply(names(comparisons), function(comp_name) {
      ttop <- comparisons[[comp_name]]
      stat_df <- ttop[, c("ID", stat_column)]
      names(stat_df)[2] <- paste(comp_name, "tval", cell_type, sep = "_")
      return(stat_df)
    })

    # Merge all comparisons for this cell type
    stat_df_merged <- Reduce(function(x, y) merge(x, y, by = "ID"), stat_list)
    all_stat_dfs[[cell_type]] <- stat_df_merged
  }

  # Merge all cell types
  stat_matrix_df <- Reduce(function(x, y) merge(x, y, by = "ID"), all_stat_dfs)

  # Convert to matrix
  rownames(stat_matrix_df) <- stat_matrix_df$ID
  stat_matrix <- as.matrix(stat_matrix_df[, -1])

  return(stat_matrix)
}


#' Save Enrichment Results
#'
#' Saves enrichment scores to RData and optionally CSV.
#'
#' @param enrichment_scores Data frame. Enrichment results (long or wide format).
#' @param output_file Character. Path for output file (without extension).
#' @param save_csv Logical. Also save as CSV. Default: TRUE
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @examples
#' save_enrichment_results(pathway_scores, "results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/pathways/pathway_scores_all")
save_enrichment_results <- function(enrichment_scores,
                                   output_file,
                                   save_csv = TRUE,
                                   verbose = TRUE) {

  # Save RData
  rdata_file <- paste0(output_file, ".RData")
  save(enrichment_scores, file = rdata_file)
  if (verbose) cat("Saved RData:", rdata_file, "\n")

  # Save CSV
  if (save_csv) {
    csv_file <- paste0(output_file, ".csv")
    write.csv(enrichment_scores, file = csv_file, row.names = FALSE)
    if (verbose) cat("Saved CSV:", csv_file, "\n")
  }
}


#' Summarize Enrichment Results
#'
#' Prints a summary of enrichment results.
#'
#' @param enrichment_scores Data frame. Enrichment results.
#' @param score_thresholds Numeric vector. Thresholds to report. Default: c(5, 10, 15)
#' @param enrichment_type Character. Type of enrichment for labeling. Default: "Terms"
#'
#' @examples
#' summarize_enrichment(pathway_scores, enrichment_type = "Pathways")
#' summarize_enrichment(tf_scores, score_thresholds = c(2, 4, 6), enrichment_type = "TFs")
summarize_enrichment <- function(enrichment_scores,
                                 score_thresholds = c(5, 10, 15),
                                 enrichment_type = "Terms") {

  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Enrichment Analysis Summary:", enrichment_type, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  n_sources <- length(unique(enrichment_scores$source))
  n_conditions <- length(unique(enrichment_scores$condition))

  cat("Total", tolower(enrichment_type), "tested:", n_sources, "\n")
  cat("Comparisons:", n_conditions, "\n\n")

  cat("Number of", tolower(enrichment_type), "above score thresholds:\n")
  for (threshold in score_thresholds) {
    scores_wide <- enrichment_to_wide(enrichment_scores)
    score_cols <- setdiff(names(scores_wide), "source")
    max_abs_scores <- apply(
      scores_wide[, score_cols],
      1,
      function(x) max(abs(x), na.rm = TRUE)
    )
    n_above <- sum(max_abs_scores > threshold, na.rm = TRUE)
    cat("  |score| >", threshold, ":", n_above, "\n")
  }

  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}
