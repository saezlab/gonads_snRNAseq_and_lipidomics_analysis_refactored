# Test script to debug heatmap generation issue
# This will run pathway analysis and test heatmap creation

library(dplyr)
library(tidyr)
library(tibble)

# Source required modules
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/differential_analysis.R")
source("scripts/core/enrichment_analysis.R")
source("scripts/config/contrasts.R")
source("scripts/visualization/heatmaps.R")

cat("Step 1: Loading or creating t-statistic matrix...\n")

# Check if we already have the tval matrix
if (file.exists("results/pseudobulk/tval_full_df.RData")) {
  cat("Loading existing tval_full_df...\n")
  load("results/pseudobulk/tval_full_df.RData")
  cat("Matrix dimensions:", nrow(tval_full_df), "x", ncol(tval_full_df), "\n")
} else {
  cat("ERROR: Please run workflow 02 first to generate tval_full_df.RData\n")
  stop("Missing tval_full_df.RData")
}

cat("\nStep 2: Running pathway enrichment...\n")

# Run pathway enrichment
pathway_results <- run_pathway_enrichment(
  stat_matrix = tval_full_df,
  pathway_file = "support/m5.all.v2023.2.Mm.symbols.gmt",
  minsize = 15,
  method = "ulm",
  verbose = TRUE
)

cat("\nPathway enrichment complete!\n")
cat("Total results:", nrow(pathway_results), "\n")
cat("Significant (|score| > 11):", sum(abs(pathway_results$score) > 11, na.rm = TRUE), "\n")

# Save for future debugging
save(pathway_results, file = "test_pathway_results.RData")
cat("Saved pathway_results to test_pathway_results.RData\n")

cat("\nStep 3: Filtering and preparing matrix...\n")

# Filter to significant pathways
pathway_results_sig <- pathway_results %>%
  filter(abs(score) > 11)

cat("Significant pathways:", nrow(pathway_results_sig), "\n")

if (nrow(pathway_results_sig) == 0) {
  cat("WARNING: No significant pathways found!\n")
  cat("Trying lower threshold (|score| > 5)...\n")
  pathway_results_sig <- pathway_results %>%
    filter(abs(score) > 5)
  cat("Pathways with |score| > 5:", nrow(pathway_results_sig), "\n")
}

if (nrow(pathway_results_sig) > 0) {
  cat("\nStep 4: Converting to wide format...\n")

  # Check structure
  cat("Columns in pathway_results_sig:", paste(colnames(pathway_results_sig), collapse = ", "), "\n")
  cat("First few rows:\n")
  print(head(pathway_results_sig))

  # Convert to wide format
  pathway_matrix <- pathway_results_sig %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score, values_fill = 0) %>%
    tibble::column_to_rownames("source")

  cat("\nMatrix dimensions after pivot_wider:", nrow(pathway_matrix), "x", ncol(pathway_matrix), "\n")
  cat("Column names:", paste(colnames(pathway_matrix), collapse = ", "), "\n")
  cat("First few row names:", paste(head(rownames(pathway_matrix)), collapse = ", "), "\n")

  # Select top 50 pathways
  if (nrow(pathway_matrix) > 50) {
    mean_scores <- rowMeans(abs(pathway_matrix), na.rm = TRUE)
    top_pathways <- names(sort(mean_scores, decreasing = TRUE)[1:50])
    pathway_matrix <- pathway_matrix[top_pathways, , drop = FALSE]
    cat("Selected top 50 pathways\n")
  }

  cat("Final matrix dimensions:", nrow(pathway_matrix), "x", ncol(pathway_matrix), "\n")

  # Save matrix for debugging
  save(pathway_matrix, file = "test_pathway_matrix.RData")
  cat("Saved pathway_matrix to test_pathway_matrix.RData\n")

  cat("\nStep 5: Testing heatmap generation...\n")

  # Test 1: Try without file output first
  cat("\nTest 1: Creating heatmap without file output...\n")
  tryCatch({
    p1 <- create_enrichment_heatmap(
      enrichment_matrix = as.matrix(pathway_matrix),
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = "Top Enriched Pathways (GO Biological Process)",
      angle_col = 315
    )
    cat("SUCCESS: Heatmap created without file output\n")
  }, error = function(e) {
    cat("ERROR creating heatmap without file:\n")
    cat("  ", conditionMessage(e), "\n")
  })

  # Test 2: Try with file output
  cat("\nTest 2: Creating heatmap with file output...\n")
  test_file <- "test_heatmap_output.pdf"
  if (file.exists(test_file)) unlink(test_file)

  tryCatch({
    p2 <- create_enrichment_heatmap(
      enrichment_matrix = as.matrix(pathway_matrix),
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = "Top Enriched Pathways (GO Biological Process)",
      angle_col = 315,
      filename = test_file,
      width = 11,
      height = 7.5
    )

    if (file.exists(test_file)) {
      file_size <- file.info(test_file)$size
      cat("SUCCESS: PDF created!\n")
      cat("  File:", test_file, "\n")
      cat("  Size:", file_size, "bytes\n")

      # Check if file is corrupted by trying to read it
      if (file_size < 1000) {
        cat("  WARNING: File size is suspiciously small - may be corrupted\n")
      } else {
        cat("  File appears valid (size > 1000 bytes)\n")
      }
    } else {
      cat("ERROR: File was not created\n")
    }
  }, error = function(e) {
    cat("ERROR creating heatmap with file:\n")
    cat("  ", conditionMessage(e), "\n")

    if (file.exists(test_file)) {
      cat("  Corrupted file was created - deleting\n")
      unlink(test_file)
    }
  })

  # Test 3: Try using base pheatmap directly
  cat("\nTest 3: Using base pheatmap directly...\n")
  test_file2 <- "test_heatmap_base.pdf"
  if (file.exists(test_file2)) unlink(test_file2)

  tryCatch({
    library(pheatmap)

    # Create color palette
    color_palette <- colorRampPalette(c("red", "grey90", "royalblue3"))(100)
    max_val <- max(abs(pathway_matrix), na.rm = TRUE)
    breaks <- seq(-max_val, max_val, length.out = 101)

    pheatmap(
      pathway_matrix,
      color = color_palette,
      breaks = breaks,
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize = 10,
      border_color = NA,
      main = "Top Enriched Pathways (Test)",
      filename = test_file2,
      width = 11,
      height = 7.5
    )

    if (file.exists(test_file2)) {
      file_size <- file.info(test_file2)$size
      cat("SUCCESS: Base pheatmap PDF created!\n")
      cat("  File:", test_file2, "\n")
      cat("  Size:", file_size, "bytes\n")
    }
  }, error = function(e) {
    cat("ERROR with base pheatmap:\n")
    cat("  ", conditionMessage(e), "\n")
  })

} else {
  cat("ERROR: No pathways to plot even with lower threshold\n")
}

cat("\n=== DEBUG COMPLETE ===\n")
