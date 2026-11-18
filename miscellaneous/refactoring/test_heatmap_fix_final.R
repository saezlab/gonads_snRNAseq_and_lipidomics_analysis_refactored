# Test script to verify heatmap fix works
# This should create a valid, non-corrupted PDF
#
# Run this in R or RStudio after sourcing the fixed heatmaps.R

library(dplyr)
library(tidyr)
library(tibble)

cat("\n========================================\n")
cat("TESTING HEATMAP FIX\n")
cat("========================================\n\n")

# Source the fixed visualization module
cat("Loading visualization module...\n")
source("scripts/visualization/heatmaps.R")

# Load pathway results
cat("Loading pathway results...\n")
load("results/pseudobulk/pathways/pathway_scores_df_all.RData")

cat("Pathway results loaded:\n")
cat("  Total rows:", nrow(pathway_results), "\n")
cat("  Significant (|score| > 11):", sum(abs(pathway_results$score) > 11, na.rm = TRUE), "\n\n")

# Filter and prepare matrix
cat("Preparing pathway matrix...\n")
pathway_results_sig <- pathway_results %>%
  filter(abs(score) > 11)

if (nrow(pathway_results_sig) > 0) {
  # Convert to wide format
  pathway_matrix <- pathway_results_sig %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score, values_fill = 0) %>%
    tibble::column_to_rownames("source")

  # Select top 50
  if (nrow(pathway_matrix) > 50) {
    mean_scores <- rowMeans(abs(pathway_matrix), na.rm = TRUE)
    top_pathways <- names(sort(mean_scores, decreasing = TRUE)[1:50])
    pathway_matrix <- pathway_matrix[top_pathways, , drop = FALSE]
  }

  cat("Matrix dimensions:", nrow(pathway_matrix), "x", ncol(pathway_matrix), "\n\n")

  # Test heatmap creation
  cat("Creating heatmap with fixed function...\n")
  test_filename <- "results/pseudobulk/pathways/TEST_pathway_enrichment_heatmap.pdf"

  # Remove old test file if exists
  if (file.exists(test_filename)) {
    unlink(test_filename)
    cat("Removed old test file\n")
  }

  tryCatch({
    # Create heatmap using the fixed create_enrichment_heatmap function
    heatmap_result <- create_enrichment_heatmap(
      enrichment_matrix = as.matrix(pathway_matrix),
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = "Top Enriched Pathways (GO Biological Process) - TEST",
      angle_col = 315,
      filename = test_filename,
      width = 11,
      height = 7.5
    )

    cat("\n✓ Heatmap function completed without error\n")

    # Check if file was created
    if (file.exists(test_filename)) {
      file_info <- file.info(test_filename)
      file_size <- file_info$size

      cat("\n========================================\n")
      cat("SUCCESS!\n")
      cat("========================================\n")
      cat("File created:", test_filename, "\n")
      cat("File size:", file_size, "bytes\n")

      if (file_size > 5000) {
        cat("\n✓ File size looks good (> 5KB)\n")
        cat("✓ PDF is likely NOT corrupted\n")
        cat("\nPlease open the PDF file to verify it displays correctly.\n")
        cat("If it opens and shows the heatmap, the fix is successful!\n")
      } else if (file_size > 1000) {
        cat("\n⚠ File size is small but might be OK\n")
        cat("Please try opening the PDF to check if it's valid\n")
      } else {
        cat("\n✗ File size is very small - may be corrupted\n")
        cat("Size:", file_size, "bytes\n")
      }

      cat("\nTo use this as the actual output, run:\n")
      cat("  file.rename('", test_filename, "',\n")
      cat("              'results/pseudobulk/pathways/pathway_enrichment_heatmap.pdf')\n")

    } else {
      cat("\n========================================\n")
      cat("FAILURE\n")
      cat("========================================\n")
      cat("✗ File was not created\n")
      cat("This indicates a problem with file creation\n")
    }

  }, error = function(e) {
    cat("\n========================================\n")
    cat("ERROR\n")
    cat("========================================\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("\nStack trace:\n")
    print(e)
  })

} else {
  cat("ERROR: No significant pathways found!\n")
  cat("Cannot create heatmap with no data.\n")
}

cat("\n========================================\n")
cat("TEST COMPLETE\n")
cat("========================================\n\n")
