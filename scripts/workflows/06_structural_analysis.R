################################################################################
# Workflow 6: Lipid Structural Analysis
################################################################################
#
# This workflow demonstrates the structural analysis pipeline that identifies
# genes whose association with lipids depends on lipid structural features
# (chain length, double bond position/saturation).
#
# Author: Refactored codebase
# Date: 2025
#
# Replaces legacy scripts:
#   - scripts/lipid_rna_correaltion_for_structural_analysis.R (note: typo in filename)
#   - scripts/lipid_structure_analysis.R
#   - scripts/lipid_structure_visualisation.R
#
# Key biological insights:
#   - Identifies substrate specificity: which genes/enzymes preferentially
#     process lipids with specific chain lengths or saturation levels
#   - Example: Cyp2c38 shows structural preferences in lipid associations
#   - Reveals coordination between enzyme expression and lipid composition
#
# Structural features analyzed:
#   - Chain length (C): Total number of carbons in fatty acid chains
#   - Double bond position (DB): Number of double bonds (saturation level)
#   - Interactions: How abundance × structure affects gene-lipid correlations
#
# Regression model:
#   gene_correlation ~ lipid_class + lipid_abundance + chain_length + db_position +
#                      lipid_abundance:chain_length + lipid_abundance:db_position
#
################################################################################

# Clear workspace
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(readr)
  library(patchwork)
})

################################################################################
# 1. CONFIGURATION
################################################################################

cat(strrep("=", 80), "\n")
cat("WORKFLOW 6: LIPID STRUCTURAL ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

# Load configuration
source("scripts/config/analysis_config.R")
source("scripts/config/datasets.R")

# Load core functions
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/structural_analysis.R")

# Define output directory
structural_dir <- "results/lipidomic/structural_model"
if (!dir.exists(structural_dir)) {
  dir.create(structural_dir, recursive = TRUE)
}

################################################################################
# 2. LOAD LIPID STRUCTURAL FEATURES (FROM LIMMA RESULTS)
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("2. LOADING LIPID STRUCTURAL FEATURES\n")
cat(strrep("-", 80), "\n\n")

cat("Generating structural data from limma results (FULL distribution, not filtered).\n")
cat("This ensures we have ALL lipids, not just those with |logFC| > 1.\n\n")

# Load limma results for lipidomics
cat("Loading limma results...\n")
load("results/lipidomic/ttop_full_list.RData")

# Extract gonadal type comparison (Testes vs Ovaries)
ttop_gonadal <- ttop_full_list[["Testes vs Ovaries"]][["Testes vs Ovaries"]]

cat("  Total lipids in limma results:", nrow(ttop_gonadal), "\n")
cat("  LogFC range:", round(min(ttop_gonadal$logFC), 3), "to",
    round(max(ttop_gonadal$logFC), 3), "\n")
cat("  Lipids with |logFC| >= 1:", sum(abs(ttop_gonadal$logFC) >= 1), "\n")
cat("  Lipids with |logFC| < 1:", sum(abs(ttop_gonadal$logFC) < 1), "\n\n")

# Generate structural data by parsing lipid names
cat("Parsing lipid names to extract chain length (C) and double bonds (DB)...\n")
lipid_structure <- generate_lipid_structural_data(ttop_gonadal, verbose = TRUE)

cat("Lipid structural features generated:\n")
cat("  Total lipids with valid structure:", nrow(lipid_structure), "\n")
cat("  Lipid classes:", length(unique(lipid_structure$lipid_class)), "\n")

# Display structure
cat("\nStructural features available:\n")
cat("  - lipid_name: Full lipid identifier\n")
cat("  - lipid_class: Category (PC, PE, TG, etc.)\n")
cat("  - chain_length: Total carbons in fatty acid chains\n")
cat("  - db_position: Number of double bonds\n")
cat("  - logfc: Log fold change (gonadal type comparison)\n")

# Verify we have the full distribution
cat("\nVerification - LogFC distribution covers full range:\n")
cat("  Min |logFC|:", round(min(abs(lipid_structure$logfc)), 4), "\n")
cat("  Max |logFC|:", round(max(abs(lipid_structure$logfc)), 4), "\n")

# Save structural features
write.csv(
  lipid_structure,
  file = file.path(structural_dir, "lipids_with_structure.csv"),
  row.names = FALSE
)

cat("\nSaved: results/lipidomic/structural_model/lipids_with_structure.csv\n")

# Clean up
rm(ttop_full_list, ttop_gonadal)
gc()

################################################################################
# 5. LIPID-ONLY STRUCTURAL ANALYSIS
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("3. LIPID-ONLY STRUCTURAL ANALYSIS\n")
cat(strrep("-", 80), "\n\n")

cat("Analyzing how lipid structural features predict differential abundance.\n")
cat("This simplified analysis focuses on lipids only, before gene-lipid correlations.\n\n")
cat("Model for each lipid class:\n")
cat("  logfc ~ chain_length + db_position + chain_length:db_position\n\n")

# Run class-specific lipid-only regression
lipid_only_results <- run_lipid_only_class_regression(
  lipid_structure = lipid_structure,
  min_lipids_per_class = 10,
  verbose = TRUE
)

# Extract significant associations
lipid_only_sig <- extract_lipid_structural_associations(
  class_results = lipid_only_results,
  p_cutoff = 0.05
)

# Save lipid-only results
write.csv(
  lipid_only_results$coefficients,
  file = file.path(structural_dir, "lipid_only_class_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  lipid_only_results$class_summaries,
  file = file.path(structural_dir, "lipid_only_class_summaries.csv"),
  row.names = FALSE
)

if (nrow(lipid_only_sig$classes_with_chain_effect) > 0) {
  write.csv(
    lipid_only_sig$classes_with_chain_effect,
    file = file.path(structural_dir, "lipid_only_chain_length_effects.csv"),
    row.names = FALSE
  )
}

if (nrow(lipid_only_sig$classes_with_db_effect) > 0) {
  write.csv(
    lipid_only_sig$classes_with_db_effect,
    file = file.path(structural_dir, "lipid_only_db_position_effects.csv"),
    row.names = FALSE
  )
}

if (nrow(lipid_only_sig$classes_with_interaction) > 0) {
  write.csv(
    lipid_only_sig$classes_with_interaction,
    file = file.path(structural_dir, "lipid_only_interaction_effects.csv"),
    row.names = FALSE
  )
}

cat("\nSaved lipid-only structural analysis results:\n")
cat("  - lipid_only_class_coefficients.csv\n")
cat("  - lipid_only_class_summaries.csv\n")
cat("  - lipid_only_chain_length_effects.csv (if any)\n")
cat("  - lipid_only_db_position_effects.csv (if any)\n")
cat("  - lipid_only_interaction_effects.csv (if any)\n")

#-------------------------------------------------------------------------------
# 5b. VISUALIZE LIPID-ONLY STRUCTURAL EFFECTS
#-------------------------------------------------------------------------------

cat("\n")
cat("Creating visualizations for lipid-only structural analysis...\n\n")

lipid_only_plots_dir <- file.path(structural_dir, "lipid_only_plots")
if (!dir.exists(lipid_only_plots_dir)) {

  dir.create(lipid_only_plots_dir, recursive = TRUE)
}

#-------------------------------------------------------------------------------
# Define FOCUS LIPID CLASSES for visualization
# Only these classes will be plotted, and axis ranges will be calculated from them
#-------------------------------------------------------------------------------
focus_lipid_classes <- c("SM", "BMP", "PG", "TG")

cat("Focus lipid classes for visualization:", paste(focus_lipid_classes, collapse = ", "), "\n\n")

# Get classes with sufficient data, filtered to focus classes
class_counts <- lipid_structure %>%
  dplyr::group_by(lipid_class) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::filter(n >= 10) %>%
  dplyr::filter(lipid_class %in% focus_lipid_classes)

cat("Classes with sufficient data (n >= 10):\n")
for (i in seq_len(nrow(class_counts))) {
  cat("  ", class_counts$lipid_class[i], ": n =", class_counts$n[i], "\n")
}
cat("\n")

#-------------------------------------------------------------------------------
# Calculate GLOBAL axis ranges using ONLY focus lipid classes
#-------------------------------------------------------------------------------
cat("Calculating global axis ranges from focus lipid classes only...\n")

# Filter to valid data across focus classes only
all_valid_data <- lipid_structure %>%
  dplyr::filter(lipid_class %in% focus_lipid_classes) %>%
  dplyr::filter(!is.na(chain_length) & !is.na(db_position) & !is.na(logfc))

# Calculate global ranges with 5% padding
global_chain_range <- range(all_valid_data$chain_length, na.rm = TRUE)
global_chain_padding <- diff(global_chain_range) * 0.05
global_chain_limits <- c(global_chain_range[1] - global_chain_padding,
                          global_chain_range[2] + global_chain_padding)

global_db_range <- range(all_valid_data$db_position, na.rm = TRUE)
global_db_padding <- diff(global_db_range) * 0.05
global_db_limits <- c(global_db_range[1] - global_db_padding,
                       global_db_range[2] + global_db_padding)

global_logfc_range <- range(all_valid_data$logfc, na.rm = TRUE)
global_logfc_padding <- diff(global_logfc_range) * 0.05
global_logfc_limits <- c(global_logfc_range[1] - global_logfc_padding,
                          global_logfc_range[2] + global_logfc_padding)

cat("  Chain Length (C) range:", round(global_chain_limits[1], 1), "to", round(global_chain_limits[2], 1), "\n")
cat("  Double Bonds (DB) range:", round(global_db_limits[1], 1), "to", round(global_db_limits[2], 1), "\n")
cat("  LogFC axis range:", round(global_logfc_limits[1], 2), "to", round(global_logfc_limits[2], 2), "\n")

#-------------------------------------------------------------------------------
# Fixed COLOR SCALE limits for the C × DB interaction plot
# This ensures colors are vivid and comparable across lipid classes
# Values beyond these limits will be capped at the extreme colors
#-------------------------------------------------------------------------------
logfc_color_limits <- c(-1, 1)

cat("  LogFC COLOR scale (capped):", logfc_color_limits[1], "to", logfc_color_limits[2], "\n\n")

#-------------------------------------------------------------------------------
# Helper function to create structural effects plot
#-------------------------------------------------------------------------------
create_structural_plot <- function(class_data, cls, lipid_only_results,
                                    chain_limits, db_limits, logfc_limits,
                                    logfc_color_limits,
                                    filter_label = "") {

  # Get significance info for this class from multivariate regression
  class_sig <- lipid_only_results$coefficients %>%
    dplyr::filter(lipid_class == cls)

  # Extract p-values from multivariate model
  chain_pval <- class_sig %>%
    dplyr::filter(term == "chain_length") %>%
    dplyr::pull(p_value)
  chain_pval <- if(length(chain_pval) > 0) chain_pval else NA

  db_pval <- class_sig %>%
    dplyr::filter(term == "db_position") %>%
    dplyr::pull(p_value)
  db_pval <- if(length(db_pval) > 0) db_pval else NA

  interaction_pval <- class_sig %>%
    dplyr::filter(term == "chain_length:db_position") %>%
    dplyr::pull(p_value)
  interaction_pval <- if(length(interaction_pval) > 0) interaction_pval else NA

  chain_sig <- !is.na(chain_pval) && chain_pval < 0.05
  db_sig <- !is.na(db_pval) && db_pval < 0.05
  interaction_sig <- !is.na(interaction_pval) && interaction_pval < 0.05

  # Format p-values for display
  format_pval <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return(sprintf("%.2e", p))
    if (p < 0.01) return(sprintf("%.3f", p))
    return(sprintf("%.2f", p))
  }

  chain_pval_text <- format_pval(chain_pval)
  db_pval_text <- format_pval(db_pval)
  interaction_pval_text <- format_pval(interaction_pval)

  # Color coding based on significance
  chain_color <- if(chain_sig) "#00BA38" else "grey50"
  db_color <- if(db_sig) "#00BA38" else "grey50"
  interaction_color <- if(interaction_sig) "#00BA38" else "grey50"

  # Get model R-squared
  class_rsq <- lipid_only_results$class_summaries %>%
    dplyr::filter(lipid_class == cls) %>%
    dplyr::pull(r_squared)
  class_rsq <- if(length(class_rsq) > 0) round(class_rsq, 3) else NA

  # Clean class name for display
  cls_clean <- gsub("_with_name", "", cls)

  # Plot 1: LogFC vs Chain Length with consistent axes
  p1 <- ggplot(class_data, aes(x = chain_length, y = logfc)) +
    geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB", shape = 15) +
    geom_smooth(method = "lm", se = TRUE, color = chain_color, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(xlim = chain_limits, ylim = logfc_limits) +
    labs(
      title = paste0(cls_clean, " (n=", nrow(class_data), ", R²=", class_rsq, ")", filter_label),
      subtitle = paste0("Chain Length: p = ", chain_pval_text, " *"),
      x = "Chain Length (C)",
      y = "LogFC (Testes vs Ovaries)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = chain_color, face = "bold", size = 11),
      axis.title.x = element_text(color = chain_color, face = if(chain_sig) "bold" else "plain")
    )

  # Plot 2: LogFC vs DB Position with consistent axes
  p2 <- ggplot(class_data, aes(x = db_position, y = logfc)) +
    geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB", shape = 15) +
    geom_smooth(method = "lm", se = TRUE, color = db_color, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(xlim = db_limits, ylim = logfc_limits) +
    labs(
      title = "",
      subtitle = paste0("Double Bonds: p = ", db_pval_text, " *"),
      x = "Double Bonds (DB)",
      y = "LogFC (Testes vs Ovaries)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.subtitle = element_text(color = db_color, face = "bold", size = 11),
      axis.title.x = element_text(color = db_color, face = if(db_sig) "bold" else "plain")
    )

  # Plot 3: 2D plot showing C × DB interaction (colored by logFC) with consistent axes
  # Uses CAPPED color limits for vivid colors; values beyond limits are squished to extremes
  p3 <- ggplot(class_data, aes(x = chain_length, y = db_position, color = logfc)) +
    geom_point(size = 3, alpha = 1, shape = 15) +
    scale_color_gradient2(low = "#e54141", mid = "white", high = "#4173d8",
                          midpoint = 0, name = "LogFC",
                          limits = logfc_color_limits,
                          oob = scales::squish) +
    coord_cartesian(xlim = chain_limits, ylim = db_limits) +
    labs(
      title = "",
      subtitle = paste0("C × DB Interaction: p = ", interaction_pval_text, " *"),
      x = "Chain Length (C)",
      y = "Double Bonds (DB)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.subtitle = element_text(color = interaction_color, face = "bold", size = 11),
      legend.position = "right"
    )

  # Combine plots
  combined <- p1 + p2 + p3 +
    plot_layout(ncol = 3) +
    plot_annotation(
      title = paste0("Lipid Structural Effects: ", cls_clean, filter_label),
      subtitle = "* p-values from multivariate regression: logfc ~ chain_length + db_position + chain_length:db_position\n  (Green = significant at p < 0.05; regression lines show bivariate trends for visualization)",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30")
      )
    )

  return(combined)
}

#-------------------------------------------------------------------------------
# Create scatter plots for each lipid class (two versions: all data and |LFC| > 1)
#-------------------------------------------------------------------------------
cat("Generating plots for each lipid class (two versions each)...\n\n")

for (cls in class_counts$lipid_class) {

  # Get all data for this class
  class_data_all <- lipid_structure %>%
    dplyr::filter(lipid_class == cls) %>%
    dplyr::filter(!is.na(chain_length) & !is.na(db_position) & !is.na(logfc))

  if (nrow(class_data_all) < 10) next

  cls_clean <- gsub("_with_name", "", cls)

  # VERSION 1: All data (original plot)
  combined_all <- create_structural_plot(
    class_data = class_data_all,
    cls = cls,
    lipid_only_results = lipid_only_results,
    chain_limits = global_chain_limits,
    db_limits = global_db_limits,
    logfc_limits = global_logfc_limits,
    logfc_color_limits = logfc_color_limits,
    filter_label = ""
  )

  ggsave(
    filename = file.path(lipid_only_plots_dir, paste0(cls_clean, "_structural_effects.pdf")),
    plot = combined_all,
    width = 14,
    height = 5
  )

  cat("  ", cls_clean, ": n=", nrow(class_data_all), "\n", sep = "")
}

cat("\n  Saved individual class plots to: results/lipidomic/structural_model/lipid_only_plots/\n")

#-------------------------------------------------------------------------------
# 5c. SUMMARY FIGURE FOR LIPID-ONLY ANALYSIS
#-------------------------------------------------------------------------------

cat("\nCreating summary heatmap of structural effects across lipid classes...\n")

# Create summary heatmap
if (nrow(lipid_only_results$coefficients) > 0) {

  # Prepare data for heatmap
  heatmap_data <- lipid_only_results$coefficients %>%
    dplyr::mutate(
      lipid_class = gsub("_with_name", "", lipid_class),
      # Use signed -log10(p) for visualization
      signed_log_p = -log10(p_value) * sign(estimate)
    ) %>%
    dplyr::select(lipid_class, term, signed_log_p) %>%
    tidyr::pivot_wider(names_from = term, values_from = signed_log_p) %>%
    tibble::column_to_rownames("lipid_class")

  # Reorder columns
  col_order <- c("chain_length", "db_position", "chain_length:db_position")
  col_order <- col_order[col_order %in% colnames(heatmap_data)]
  heatmap_data <- heatmap_data[, col_order, drop = FALSE]

  # Order rows by chain_length effect
  if ("chain_length" %in% colnames(heatmap_data)) {
    row_order <- order(heatmap_data$chain_length, na.last = TRUE, decreasing = TRUE)
    heatmap_data <- heatmap_data[row_order, , drop = FALSE]
  }

  # Create heatmap
  suppressPackageStartupMessages(library(pheatmap))

  max_val <- max(abs(heatmap_data), na.rm = TRUE)
  if (is.finite(max_val) && max_val > 0) {
    breaks <- seq(-max_val, max_val, length.out = 101)
  } else {
    breaks <- seq(-3, 3, length.out = 101)
  }

  pdf(file.path(structural_dir, "lipid_only_summary_heatmap.pdf"), width = 8, height = 10)
  pheatmap(
    as.matrix(heatmap_data),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = FALSE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = breaks,
    na_col = "grey90",
    main = "Lipid Structural Effects on Differential Abundance\n(signed -log10 p-value)",
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 11,
    angle_col = "45"
  )
  dev.off()

  cat("  Saved: lipid_only_summary_heatmap.pdf\n")
}

cat("\nLipid-only structural analysis complete!\n")
cat("This provides a baseline understanding of how structural features relate to\n")
cat("lipid abundance changes before examining gene-lipid correlations.\n")

################################################################################
# WORKFLOW COMPLETE
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("LIPID STRUCTURAL ANALYSIS WORKFLOW COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("Output files:\n")
cat("  - results/lipidomic/structural_model/lipids_with_structure.csv\n")
cat("  - results/lipidomic/structural_model/lipid_only_class_coefficients.csv\n")
cat("  - results/lipidomic/structural_model/lipid_only_class_summaries.csv\n")
cat("  - results/lipidomic/structural_model/lipid_only_chain_length_effects.csv\n")
cat("  - results/lipidomic/structural_model/lipid_only_db_position_effects.csv\n")
cat("  - results/lipidomic/structural_model/lipid_only_interaction_effects.csv\n")
cat("  - results/lipidomic/structural_model/lipid_only_summary_heatmap.pdf\n")
cat("  - results/lipidomic/structural_model/lipid_only_plots/[CLASS]_structural_effects.pdf\n")

cat("\n")
cat(strrep("=", 80), "\n\n")

# NOTE: Gene-lipid correlation structural analysis (sections 6-14) is
# commented out for now. It can be re-enabled when needed by uncommenting
# the code below. That analysis requires loading MOFA weights and RNA data
# (sections 2-3 in the original workflow) which are also commented out above.
