# Simple debug script - load existing pathway results and create heatmap
# Run this in RStudio or R console

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

cat("Loading pathway results...\n")
load("results/pseudobulk/pathways/pathway_scores_df_all.RData")

cat("Pathway results loaded:\n")
cat("  Rows:", nrow(pathway_results), "\n")
cat("  Columns:", paste(colnames(pathway_results), collapse = ", "), "\n")
cat("  Score range:", min(pathway_results$score, na.rm = TRUE), "to",
    max(pathway_results$score, na.rm = TRUE), "\n")
cat("  Significant (|score| > 11):", sum(abs(pathway_results$score) > 11, na.rm = TRUE), "\n\n")

# Filter to significant
pathway_results_sig <- pathway_results %>%
  filter(abs(score) > 11)

cat("Significant pathways:", nrow(pathway_results_sig), "\n")

if (nrow(pathway_results_sig) > 0) {
  # Show sample of data
  cat("\nFirst few rows:\n")
  print(head(pathway_results_sig))

  # Convert to wide format
  cat("\nConverting to wide format...\n")
  pathway_matrix <- pathway_results_sig %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score, values_fill = 0) %>%
    tibble::column_to_rownames("source")

  cat("Matrix dimensions:", nrow(pathway_matrix), "x", ncol(pathway_matrix), "\n")
  cat("Matrix class:", class(pathway_matrix), "\n")

  # Select top 50
  if (nrow(pathway_matrix) > 50) {
    mean_scores <- rowMeans(abs(pathway_matrix), na.rm = TRUE)
    top_pathways <- names(sort(mean_scores, decreasing = TRUE)[1:50])
    pathway_matrix <- pathway_matrix[top_pathways, , drop = FALSE]
    cat("Selected top 50 pathways\n")
  }

  cat("Final matrix dimensions:", nrow(pathway_matrix), "x", ncol(pathway_matrix), "\n")

  # Check for any issues
  cat("\nChecking matrix:\n")
  cat("  Any NA values:", any(is.na(pathway_matrix)), "\n")
  cat("  Any Inf values:", any(is.infinite(as.matrix(pathway_matrix))), "\n")
  cat("  Value range:", min(pathway_matrix, na.rm = TRUE), "to",
      max(pathway_matrix, na.rm = TRUE), "\n")

  # Create color palette and breaks
  cat("\nCreating color palette...\n")
  color_palette <- colorRampPalette(c("red", "grey90", "royalblue3"))(100)
  max_val <- max(abs(pathway_matrix), na.rm = TRUE)
  breaks <- seq(-max_val, max_val, length.out = 101)
  cat("  Color breaks range:", min(breaks), "to", max(breaks), "\n")

  # Try creating heatmap
  cat("\nAttempt 1: Basic pheatmap to screen...\n")
  tryCatch({
    p <- pheatmap(
      as.matrix(pathway_matrix),
      color = color_palette,
      breaks = breaks,
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize = 10,
      border_color = NA,
      main = "Test Heatmap"
    )
    cat("SUCCESS: Heatmap displayed\n")
  }, error = function(e) {
    cat("ERROR:\n")
    print(e)
  })

  # Try saving to file
  cat("\nAttempt 2: Save directly with pheatmap...\n")
  test_file <- "results/pseudobulk/pathways/test_direct_save.pdf"
  if (file.exists(test_file)) unlink(test_file)

  tryCatch({
    pheatmap(
      as.matrix(pathway_matrix),
      color = color_palette,
      breaks = breaks,
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize = 10,
      border_color = NA,
      main = "Top Enriched Pathways (GO Biological Process)",
      filename = test_file,
      width = 11,
      height = 7.5
    )

    if (file.exists(test_file)) {
      size <- file.info(test_file)$size
      cat("SUCCESS: PDF created\n")
      cat("  File:", test_file, "\n")
      cat("  Size:", size, "bytes\n")

      if (size < 1000) {
        cat("  WARNING: File is very small, may be corrupted\n")
      }
    } else {
      cat("ERROR: File was not created\n")
    }
  }, error = function(e) {
    cat("ERROR:\n")
    print(e)
  })

  # Try manual PDF approach
  cat("\nAttempt 3: Manual PDF device...\n")
  test_file2 <- "results/pseudobulk/pathways/test_manual_pdf.pdf"
  if (file.exists(test_file2)) unlink(test_file2)

  tryCatch({
    pdf(test_file2, width = 11, height = 7.5)

    pheatmap(
      as.matrix(pathway_matrix),
      color = color_palette,
      breaks = breaks,
      display_numbers = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize = 10,
      border_color = NA,
      main = "Top Enriched Pathways (Manual PDF)"
    )

    dev.off()

    if (file.exists(test_file2)) {
      size <- file.info(test_file2)$size
      cat("SUCCESS: Manual PDF created\n")
      cat("  File:", test_file2, "\n")
      cat("  Size:", size, "bytes\n")
    }
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("ERROR:\n")
    print(e)
    if (file.exists(test_file2)) unlink(test_file2)
  })

} else {
  cat("No significant pathways to plot\n")
}

cat("\n=== DEBUG COMPLETE ===\n")
cat("Please run this script in R/RStudio and report results\n")
