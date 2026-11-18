## Quick script to regenerate volcano plots from saved ttop_full_list.RData
## Uses the same logic as workflow 02, step 4

library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")
source("scripts/visualization/volcano_plots.R")

# Load saved DE results
load("results/pseudobulk/ttop_full_list.RData")

dir.create("results/pseudobulk/volcanos", recursive = TRUE, showWarnings = FALSE)

# Define manual labels for specific contrast + cell type combinations
manual_label_map <- list(
  "XY vs XX__Hepatocytes" = c("Xist", "Tsix", "Kdm6a", "Ddx3y", "Uty", "Kdm5d"),
  "Testes vs Ovaries__Hepatocytes" = c("Cyp2b9", "Cyp2b13", "Stambpl1", "Cyp2d9", "Mup20", "Slco1a1")
)

# Create volcano plots for each contrast and cell type
for (contrast_name in names(ttop_full_list)) {
  cat("Processing:", contrast_name, "\n")

  for (cell_type in names(ttop_full_list[[contrast_name]])) {
    cat("  ", cell_type, "...\n", sep = "")

    tryCatch({
      ttop <- ttop_full_list[[contrast_name]][[cell_type]]

      ttop$P.Value <- ifelse(-log10(ttop$P.Value) > 25, 10^-25, ttop$P.Value)

      # Check if this combination has manual labels
      map_key <- paste0(contrast_name, "__", cell_type)
      has_labels <- map_key %in% names(manual_label_map)

      # Create volcano plot
      p <- volcano_nice(
        df = ttop,
        FCIndex = which(names(ttop) == "logFC"),
        pValIndex = which(names(ttop) == "P.Value"),
        IDIndex = which(names(ttop) == "ID"),
        label = has_labels,
        nlabels = 0,
        manual_labels = if (has_labels) manual_label_map[[map_key]] else NA
      )

      # Add labels
      p <- customize_volcano(
        p,
        title = paste(cell_type, "-", contrast_name),
        xlab = "log2 Fold Change",
        ylab = "-log10(p-value)"
      )

      # Save plot - square aspect ratio
      filename <- paste0("volcano_", gsub(" ", "_", contrast_name), "_",
                        gsub(" ", "_", cell_type), ".pdf")
      filepath <- file.path("results/pseudobulk/volcanos", filename)
      ggsave(filepath, plot = p, width = 5, height = 5)
    }, error = function(e) {
      cat("    WARNING: Error generating volcano plot for", contrast_name, "-", cell_type, "\n")
      cat("    Error message:", conditionMessage(e), "\n")
    })
  }
  cat("\n")
}

cat("Volcano plots complete!\n")
