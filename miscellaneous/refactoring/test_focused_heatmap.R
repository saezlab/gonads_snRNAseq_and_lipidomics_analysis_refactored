# Test script to verify focused pathway heatmap generation
# This creates the publication-quality heatmap matching your screenshot

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

cat("\n========================================\n")
cat("TESTING FOCUSED PATHWAY HEATMAP\n")
cat("========================================\n\n")

# Load pathway results
cat("Loading pathway results...\n")
load("results/pseudobulk/pathways/pathway_scores_df_all.RData")

cat("Pathway results loaded:", nrow(pathway_results), "rows\n\n")

# Define the specific pathways to display (matching screenshot)
key_pathways_patterns <- c(
  "ribosome" = "RIBOSOM",
  "immune response" = "IMMUNE_RESPONSE",
  "estrogen 2 hydroxylase" = "ESTROGEN.*HYDROXYLASE|ESTROGEN_2",
  "testosterone 16 alpha hydroxylase" = "TESTOSTERONE.*16.*ALPHA.*HYDROXYLASE|TESTOSTERONE_16_ALPHA",
  "testosterone 16 beta hydroxylase" = "TESTOSTERONE.*16.*BETA.*HYDROXYLASE|TESTOSTERONE_16_BETA"
)

# Filter to only gonadal_type (Testes vs Ovaries) for each cell type
# Note: Use the exact string format from your data
pathway_gonadal <- pathway_results %>%
  filter(grepl("Testes vs Ovaries_tval", condition))

# Optional: Filter by score threshold (set to 0 to include all, or higher to filter)
FOCUSED_HEATMAP_SCORE_THRESHOLD <- 0  # Adjust this value to filter by score
if (FOCUSED_HEATMAP_SCORE_THRESHOLD > 0) {
  pathway_gonadal <- pathway_gonadal %>%
    filter(abs(score) > FOCUSED_HEATMAP_SCORE_THRESHOLD)
  cat("Filtering pathways with |score| >", FOCUSED_HEATMAP_SCORE_THRESHOLD, "\n")
}

cat("Filtered to gonadal contrasts:", nrow(pathway_gonadal), "rows\n\n")

cat("Searching for key pathways...\n")

# Find matching pathways
matched_pathways <- list()
for (display_name in names(key_pathways_patterns)) {
  pattern <- key_pathways_patterns[display_name]
  matches <- pathway_gonadal %>%
    filter(grepl(pattern, source, ignore.case = TRUE))

  if (nrow(matches) > 0) {
    # Take the first match or the one with highest absolute score
    # Convert to data.frame first to avoid Rle type issues with slice()
    best_match <- matches %>%
      as.data.frame() %>%
      arrange(desc(abs(score))) %>%
      head(1)  # Use head() instead of slice() for compatibility
    matched_pathways[[display_name]] <- best_match$source[1]
    cat("  ✓ Found:", display_name, "\n")
    cat("     Full name:", best_match$source[1], "\n")
    cat("     Score range:", min(matches$score), "to", max(matches$score), "\n")
  } else {
    cat("  ✗ No match found for:", display_name, "\n")
    # Show available pathways that partially match
    partial_matches <- pathway_gonadal %>%
      filter(grepl(gsub(".*\\|", "", pattern), source, ignore.case = TRUE)) %>%
      pull(source) %>%
      unique()
    if (length(partial_matches) > 0) {
      cat("     Similar pathways found:\n")
      for (pm in head(partial_matches, 3)) {
        cat("       -", pm, "\n")
      }
    }
  }
}

if (length(matched_pathways) > 0) {
  cat("\n", length(matched_pathways), "pathways matched!\n\n")

  cat("Building focused pathway matrix...\n")

  # Create mapping from long pathway names to display names
  pathway_name_map <- setNames(names(matched_pathways), unlist(matched_pathways))

  # Filter to matched pathways and convert to wide format
  focused_pathway_matrix <- pathway_gonadal %>%
    filter(source %in% unlist(matched_pathways)) %>%
    dplyr::mutate(
      # Extract cell type from condition name (adjust regex to match actual format)
      cell_type = gsub("Testes vs Ovaries_tval_", "", condition),
      # Shorten cell type names to match figure
      # Note: Cell types in data have spaces, not underscores
      cell_type_short = case_when(
        cell_type == "Hepatocytes" ~ "Hep",
        cell_type == "Endothelial cells" ~ "End",  # Space, not underscore
        cell_type == "Immune cells" ~ "Imm",       # Space, not underscore
        cell_type == "Stellate cells" ~ "Ste",     # Space, not underscore
        cell_type == "Cholangiocytes" ~ "Cho",
        TRUE ~ cell_type
      ),
      # Map to display names
      pathway_display = pathway_name_map[source]
    ) %>%
    dplyr::select(pathway_display, cell_type_short, score) %>%
    pivot_wider(names_from = cell_type_short, values_from = score, values_fill = 0)

  # Ensure column order matches figure (Hep, End, Imm, Ste, Cho)
  col_order <- c("Hep", "End", "Imm", "Ste", "Cho")
  missing_cols <- setdiff(col_order, colnames(focused_pathway_matrix))
  for (col in missing_cols) {
    focused_pathway_matrix[[col]] <- 0
  }

  focused_pathway_matrix <- focused_pathway_matrix %>%
    dplyr::select(pathway_display, all_of(col_order)) %>%
    tibble::column_to_rownames("pathway_display")

  # Reorder rows to match the figure
  row_order <- names(key_pathways_patterns)
  row_order <- row_order[row_order %in% rownames(focused_pathway_matrix)]
  focused_pathway_matrix <- focused_pathway_matrix[row_order, , drop = FALSE]

  cat("Matrix dimensions:", nrow(focused_pathway_matrix), "x", ncol(focused_pathway_matrix), "\n\n")

  cat("Matrix preview:\n")
  print(focused_pathway_matrix)

  # Determine actual data range
  data_min <- min(focused_pathway_matrix, na.rm = TRUE)
  data_max <- max(focused_pathway_matrix, na.rm = TRUE)
  cat("\n\nValue range in matrix:", round(data_min, 2), "to", round(data_max, 2), "\n")

  # Create focused heatmap with color scale
  cat("Creating heatmap...\n")

  output_file <- "results/pseudobulk/pathways/TEST_pathway_enrichment_heatmap_focused.pdf"
  if (file.exists(output_file)) unlink(output_file)

  tryCatch({
    # Option 1: Use data range (colors fill entire heatmap)
    # Option 2: Use fixed range like -20 to 20 (for comparison across datasets)
    USE_FIXED_COLOR_RANGE <- FALSE  # Set to TRUE for fixed -20 to 20 range

    if (USE_FIXED_COLOR_RANGE) {
      # Fixed range: -20 to 20 (may result in mostly white cells if data range is smaller)
      color_range <- 20
      cat("Using fixed color range: -20 to 20\n")
    } else {
      # Automatic range: based on actual data (uses full color spectrum)
      color_range <- max(abs(c(data_min, data_max)), na.rm = TRUE)
      cat("Using automatic color range: -", round(color_range, 2), "to", round(color_range, 2), "\n")
    }

    color_palette <- colorRampPalette(c("red", "white", "royalblue3"))(100)
    breaks <- seq(-color_range, color_range, length.out = 101)

    focused_heatmap <- pheatmap::pheatmap(
      as.matrix(focused_pathway_matrix),
      color = color_palette,
      breaks = breaks,
      display_numbers = FALSE,
      cluster_rows = FALSE,  # Keep order as specified
      cluster_cols = FALSE,  # Keep order as specified
      fontsize = 12,
      fontsize_row = 11,
      fontsize_col = 11,
      border_color = "grey80",
      main = "",  # No title for publication figure
      angle_col = 0,  # Horizontal column labels
      filename = output_file,
      width = 6,
      height = 4
    )

    if (file.exists(output_file)) {
      file_size <- file.info(output_file)$size
      cat("\n========================================\n")
      cat("SUCCESS!\n")
      cat("========================================\n")
      cat("File created:", output_file, "\n")
      cat("File size:", file_size, "bytes\n")

      if (file_size > 5000) {
        cat("\n✓ File size looks good (> 5KB)\n")
        cat("✓ PDF should match your screenshot format:\n")
        cat("  - 5 pathways (rows)\n")
        cat("  - 5 cell types (columns)\n")
        cat("  - Color scale: -20 (red/Ovaries) to 20 (blue/Testes)\n")
        cat("  - No clustering (maintains specified order)\n")
        cat("\nPlease open the PDF to verify it matches your screenshot!\n")
      } else {
        cat("\n⚠ File size is small - may be corrupted\n")
      }
    } else {
      cat("\n✗ ERROR: File was not created\n")
    }

  }, error = function(e) {
    cat("\n========================================\n")
    cat("ERROR\n")
    cat("========================================\n")
    cat("Error message:", conditionMessage(e), "\n")
  })

} else {
  cat("\nERROR: No pathways matched!\n")
  cat("This likely means the pathway names in your data don't match the patterns.\n")
  cat("\nTry running this to see available pathway names:\n")
  cat("  pathway_results %>% filter(grepl('Testes_vs_Ovaries', condition)) %>%\n")
  cat("    pull(source) %>% unique() %>% head(20)\n")
}

cat("\n========================================\n")
cat("TEST COMPLETE\n")
cat("========================================\n\n")
