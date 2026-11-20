#' Core Differential Analysis Functions
#'
#' This module provides high-level functions for differential analysis (RNA and lipid).
#' It wraps the existing limma functions from support_functions.R with a cleaner,
#' more user-friendly interface.
#'
#' The module builds on top of existing functions:
#' - runLimma() - core limma wrapper
#' - limma_res_to_ttop_list() - extract toptables
#' - makeContrastsAlt() - create contrast matrices
#' - checkInputs() - validate inputs
#'
#' These are sourced from support_functions.R and not duplicated.

library(limma)
library(dplyr)

# Source existing support functions
source("scripts/support_functions.R")


#' Run Differential Analysis for Multiple Cell Types
#'
#' High-level function to run differential analysis across multiple cell types.
#' This is the main function you should use in workflows.
#'
#' @param data Matrix or data frame. Expression/abundance data (features x samples).
#' @param metadata Data frame. Sample metadata with columns "sample", "cell_type",
#'                 and grouping variable.
#' @param group_variable Character. Column name in metadata to use for grouping
#'                       (e.g., "gonadal_type", "allosome_type", "group").
#' @param comparisons Named list. Contrasts to test (e.g., list("Testes_vs_Ovaries" = c(-1, 2))).
#' @param cell_types Character vector. Cell types to analyze. If NULL, analyzes all
#'                   cell types in metadata. Default: NULL
#' @param adjust_method Character. Method for p-value adjustment. Default: "fdr"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Named list (by cell type) of named lists (by comparison) of toptables.
#'         Each toptable is a data frame with columns: ID, logFC, AveExpr, t, P.Value, adj.P.Val, B
#'
#' @details
#' This function:
#' 1. Loops through each cell type
#' 2. Subsets data to samples from that cell type
#' 3. Runs limma differential analysis for each comparison
#' 4. Extracts toptables with specified number of features
#' 5. Returns results organized by cell type and comparison
#'
#' The contrast vectors follow the pattern:
#' - c(-1, 2): Condition 2 - Condition 1
#' - c(1, -4): Condition 1 - Condition 4
#' - c(2, 2, -1, -3): (Cond 2 + Cond 2) - (Cond 1 + Cond 3) / 2
#'
#' Conditions are ordered by unique(targets$condition) after subsetting.
#'
#' @examples
#' # Analyze gonadal type effect across all cell types
#' results <- run_differential_analysis_multi_celltype(
#'   data = rna_data,
#'   metadata = metadata,
#'   group_variable = "gonadal_type",
#'   comparisons = list("Testes_vs_Ovaries" = c(-1, 2))
#' )
#'
#' # Analyze specific cell types only
#' results <- run_differential_analysis_multi_celltype(
#'   data = rna_data,
#'   metadata = metadata,
#'   group_variable = "allosome_type",
#'   comparisons = list("XY_vs_XX" = c(-1, 2)),
#'   cell_types = c("Hepatocytes", "Endothelial cells")
#' )
run_differential_analysis_multi_celltype <- function(data,
                                                      metadata,
                                                      group_variable,
                                                      comparisons,
                                                      cell_types = NULL,
                                                      adjust_method = "fdr",
                                                      verbose = TRUE) {

  # Determine which cell types to analyze
  if (is.null(cell_types)) {
    cell_types <- unique(metadata$cell_type)
  }

  if (verbose) {
    cat("Running differential analysis for", length(cell_types), "cell type(s):\n")
    cat(" ", paste(cell_types, collapse = ", "), "\n")
    cat("Group variable:", group_variable, "\n")
    cat("Comparisons:", paste(names(comparisons), collapse = ", "), "\n\n")
  }

  # Store results for all cell types
  results_all <- list()

  # Loop through each cell type
  for (cell_type in cell_types) {
    if (verbose) cat("Analyzing", cell_type, "...\n")

    # Subset metadata to this cell type
    targets <- metadata[metadata$cell_type == cell_type, ]

    # Subset data to samples from this cell type
    data_subset <- data[, targets$sample]

    # Prepare targets for limma (sample + condition columns)
    targets_limma <- targets[, c("sample", group_variable)]
    names(targets_limma)[2] <- "condition"

    if (verbose) {
      cat("  Conditions:", paste(unique(targets_limma$condition), collapse = ", "), "\n")
      cat("  Samples:", ncol(data_subset), "\n")
    }

    # Run differential analysis
    ttop_list <- run_differential_analysis(
      data = data_subset,
      targets = targets_limma,
      comparisons = comparisons,
      adjust_method = adjust_method,
      verbose = FALSE
    )

    # Rename t column to include cell type
    for (i in seq_along(ttop_list)) {
      names(ttop_list[[i]])[4] <- paste("t", cell_type, sep = "_")
    }

    # Store results
    results_all[[cell_type]] <- ttop_list

    if (verbose) cat("  ✓ Complete\n\n")
  }

  if (verbose) cat("Differential analysis complete!\n")

  return(results_all)
}


#' Run Differential Analysis (Single Dataset)
#'
#' Runs differential analysis for a single dataset (e.g., one cell type).
#' Lower-level function - most users should use run_differential_analysis_multi_celltype().
#'
#' @param data Matrix or data frame. Expression/abundance data (features x samples).
#' @param targets Data frame. Must have columns "sample" and "condition".
#' @param comparisons Named list. Contrasts to test (e.g., list("Test_vs_Control" = c(-1, 2))).
#' @param regress_out Character vector. Variables to regress out (batch effects). Default: NULL
#' @param adjust_method Character. Method for p-value adjustment. Default: "fdr"
#' @param number Integer. Number of top features to return in toptables. Default: all features
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Named list of toptables (one per comparison).
#'
#' @examples
#' # Simple usage
#' results <- run_differential_analysis(
#'   data = hepatocyte_data,
#'   targets = targets,  # must have "sample" and "condition" columns
#'   comparisons = list("Treatment_vs_Control" = c(-1, 2))
#' )
#'
#' # With batch effect regression
#' results <- run_differential_analysis(
#'   data = hepatocyte_data,
#'   targets = targets,
#'   comparisons = list("Treatment_vs_Control" = c(-1, 2)),
#'   regress_out = c("batch")
#' )
run_differential_analysis <- function(data,
                                      targets,
                                      comparisons,
                                      regress_out = NULL,
                                      adjust_method = "fdr",
                                      number = NULL,
                                      verbose = TRUE) {

  # Convert to data.frame if matrix
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # Default number: all features
  if (is.null(number)) {
    number <- nrow(data)
  }

  if (verbose) {
    cat("Running limma differential analysis...\n")
    cat("  Features:", nrow(data), "\n")
    cat("  Samples:", ncol(data), "\n")
    cat("  Conditions:", paste(unique(targets$condition), collapse = ", "), "\n")
  }

  # Run limma using existing support function
  limma_res <- runLimma(
    measurements = data,
    targets = targets,
    comparisons = comparisons,
    regress_out = regress_out
  )

  # Extract toptables using existing support function
  ttop_list <- limma_res_to_ttop_list(
    limma_res = limma_res,
    comp_names = names(comparisons),
    number = number,
    adjust.method = adjust_method
  )

  if (verbose) {
    cat("  ✓ Analysis complete\n")
    for (comp_name in names(ttop_list)) {
      n_sig <- sum(ttop_list[[comp_name]]$adj.P.Val < 0.05, na.rm = TRUE)
      cat("  ", comp_name, ":", n_sig, "features at FDR < 0.05\n")
    }
  }

  return(ttop_list)
}


#' Filter Significant Features from Toptable
#'
#' Extracts features meeting specified significance criteria.
#'
#' @param ttop Data frame. Toptable from limma (output of topTable()).
#' @param fdr_cutoff Numeric. FDR cutoff. Default: 0.05
#' @param fc_cutoff Numeric. Absolute log fold-change cutoff. Default: 0 (no filtering)
#' @param direction Character. Filter by direction: "up", "down", or "both". Default: "both"
#'
#' @return Filtered toptable data frame.
#'
#' @examples
#' sig_features <- filter_significant_features(ttop, fdr_cutoff = 0.05, fc_cutoff = 1)
#' up_features <- filter_significant_features(ttop, fdr_cutoff = 0.01, direction = "up")
filter_significant_features <- function(ttop,
                                        fdr_cutoff = 0.05,
                                        fc_cutoff = 0,
                                        direction = "both") {

  # Filter by FDR
  sig_features <- ttop[ttop$adj.P.Val < fdr_cutoff, ]

  # Filter by fold change
  if (fc_cutoff > 0) {
    sig_features <- sig_features[abs(sig_features$logFC) >= fc_cutoff, ]
  }

  # Filter by direction
  if (direction == "up") {
    sig_features <- sig_features[sig_features$logFC > 0, ]
  } else if (direction == "down") {
    sig_features <- sig_features[sig_features$logFC < 0, ]
  }

  return(sig_features)
}


#' Save Differential Analysis Results
#'
#' Saves toptables to RData and optionally to CSV files.
#'
#' @param results List. Results from run_differential_analysis() or
#'                run_differential_analysis_multi_celltype().
#' @param output_dir Character. Directory to save results. Default: "results/"
#' @param prefix Character. Prefix for output files. Default: ""
#' @param save_csv Logical. Also save as CSV files. Default: TRUE
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @examples
#' save_da_results(results, output_dir = "results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/", prefix = "fcg_main")
save_da_results <- function(results,
                            output_dir = "results/",
                            prefix = "",
                            save_csv = TRUE,
                            verbose = TRUE) {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created directory:", output_dir, "\n")
  }

  # Determine if results are multi-cell-type or single
  is_multi_celltype <- all(sapply(results, is.list))

  # Save as RData
  rdata_file <- file.path(output_dir, paste0(prefix, "ttop_full_list.RData"))
  save(results, file = rdata_file)
  if (verbose) cat("Saved RData:", rdata_file, "\n")

  # Optionally save as CSV
  if (save_csv) {
    if (is_multi_celltype) {
      # Multi-cell-type results: save each cell type x comparison
      for (cell_type in names(results)) {
        for (comparison in names(results[[cell_type]])) {
          csv_file <- file.path(
            output_dir,
            paste0(prefix, cell_type, "_", comparison, ".csv")
          )
          write.csv(results[[cell_type]][[comparison]], file = csv_file, row.names = FALSE)
          if (verbose) cat("  Saved CSV:", csv_file, "\n")
        }
      }
    } else {
      # Single dataset results: save each comparison
      for (comparison in names(results)) {
        csv_file <- file.path(output_dir, paste0(prefix, comparison, ".csv"))
        write.csv(results[[comparison]], file = csv_file, row.names = FALSE)
        if (verbose) cat("  Saved CSV:", csv_file, "\n")
      }
    }
  }

  if (verbose) cat("✓ Results saved successfully\n")
}


#' Convert Toptable List to T-statistic Matrix
#'
#' Extracts t-statistics from multiple toptables and combines into a matrix.
#' Useful for pathway/TF enrichment analysis with decoupleR.
#'
#' @param ttop_list Named list of toptables (from run_differential_analysis()).
#' @param stat_column Character. Column name containing statistics. Default: "t"
#'                    (can also use "logFC", "B", etc.)
#'
#' @return Matrix with features as rows and comparisons as columns.
#'
#' @examples
#' t_matrix <- ttop_list_to_stat_matrix(results$Hepatocytes)
#' fc_matrix <- ttop_list_to_stat_matrix(results$Hepatocytes, stat_column = "logFC")
ttop_list_to_stat_matrix <- function(ttop_list, stat_column = "t") {

  # Get all unique feature IDs across all comparisons
  all_features <- unique(unlist(lapply(ttop_list, function(x) x$ID)))

  # Create empty matrix
  stat_matrix <- matrix(NA, nrow = length(all_features), ncol = length(ttop_list))
  rownames(stat_matrix) <- all_features
  colnames(stat_matrix) <- names(ttop_list)

  # Fill in statistics
  for (comparison in names(ttop_list)) {
    ttop <- ttop_list[[comparison]]
    stat_matrix[ttop$ID, comparison] <- ttop[[stat_column]]
  }

  return(stat_matrix)
}


#' Summarize Differential Analysis Results
#'
#' Prints a summary of significant features across comparisons and cell types.
#'
#' @param results List. Results from differential analysis.
#' @param fdr_cutoff Numeric. FDR cutoff for significance. Default: 0.05
#'
#' @examples
#' summarize_da_results(results)
#' summarize_da_results(results, fdr_cutoff = 0.01)
summarize_da_results <- function(results, fdr_cutoff = 0.05) {

  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Differential Analysis Results Summary\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("FDR cutoff:", fdr_cutoff, "\n\n")

  # Determine if multi-cell-type
  is_multi_celltype <- all(sapply(results, is.list))

  if (is_multi_celltype) {
    # Multi-cell-type results
    for (cell_type in names(results)) {
      cat("Cell Type:", cell_type, "\n")
      cat(paste(rep("-", 80), collapse = ""), "\n")

      for (comparison in names(results[[cell_type]])) {
        ttop <- results[[cell_type]][[comparison]]
        n_sig <- sum(ttop$adj.P.Val < fdr_cutoff, na.rm = TRUE)
        n_up <- sum(ttop$adj.P.Val < fdr_cutoff & ttop$logFC > 0, na.rm = TRUE)
        n_down <- sum(ttop$adj.P.Val < fdr_cutoff & ttop$logFC < 0, na.rm = TRUE)

        cat("  ", comparison, ":\n")
        cat("    Total significant:", n_sig, "\n")
        cat("    Up-regulated:     ", n_up, "\n")
        cat("    Down-regulated:   ", n_down, "\n")
      }
      cat("\n")
    }
  } else {
    # Single dataset results
    for (comparison in names(results)) {
      ttop <- results[[comparison]]
      n_sig <- sum(ttop$adj.P.Val < fdr_cutoff, na.rm = TRUE)
      n_up <- sum(ttop$adj.P.Val < fdr_cutoff & ttop$logFC > 0, na.rm = TRUE)
      n_down <- sum(ttop$adj.P.Val < fdr_cutoff & ttop$logFC < 0, na.rm = TRUE)

      cat(comparison, ":\n")
      cat("  Total significant:", n_sig, "\n")
      cat("  Up-regulated:     ", n_up, "\n")
      cat("  Down-regulated:   ", n_down, "\n")
    }
  }

  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}
