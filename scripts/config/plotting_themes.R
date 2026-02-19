#' Plotting Themes and Color Schemes
#'
#' This file centralizes all color schemes, plot dimensions, and themes
#' used across the FCG liver analysis project.
#'
#' Consolidates plotting configuration from:
#' - pseudobulk_DA.R
#' - lipidomic_DA.R
#' - Various figure scripts

library(ggplot2)

# Cell Type Colors
# ================
# Standard colors used for all cell type visualizations

CELLTYPE_COLORS <- c(
  "Hepatocytes" = "#ca3659",
  "Endothelial cells" = "#8f7b2f",
  "Immune cells" = "#0d934d",
  "Stellate cells" = "#1892b1",
  "Cholangiocytes" = "#834f9f"
)


# Condition Colors
# ================
# Colors for main experimental conditions

# Gonadal Type Colors
GONADAL_COLORS <- c(
  "Ovaries" = "#e54141",
  "Testes" = "#4173d8"
)

# Allosome Type Colors
ALLOSOME_COLORS <- c(
  "XX" = "#e54141",
  "XY" = "#4173d8"
)

# Alternative allosome colors (for specific comparisons)
ALLOSOME_COLORS_ALT <- c(
  "XX" = "brown",
  "XY" = "cyan"
)

# FCG Group Colors
GROUP_COLORS <- c(
  "XXO" = "#E74C3C",    # Red
  "XXT" = "#3498DB",    # Blue
  "XYO" = "#F39C12",    # Orange
  "XYT" = "#9B59B6"     # Purple
)

# Zonation Colors (simple)
ZONATION_COLORS <- c(
  "Periportal" = "#8B4513",     # Brown
  "Pericentral" = "#4169E1"     # Royal blue
)

# Zonation x Gonad Colors (harmonized scheme)
# Uses darker shades for pericentral, lighter shades for periportal
# Red tones for Ovaries, Blue tones for Testes
ZONATION_GONAD_COLORS <- c(
  # Testes (blue tones)
  "TPC" = "#2b64d6",            # Testes Pericentral - Dark blue
  "TPP" = "#5f8dff",            # Testes Periportal - Light blue
  # Ovaries (red tones)
  "OPC" = "#d32f2f",            # Ovaries Pericentral - Dark red
  "OPP" = "#ff5555"             # Ovaries Periportal - Light red
)

# Full names version
ZONATION_GONAD_COLORS_FULL <- c(
  "Testes_Pericentral" = "#2b64d6",
  "Testes_Periportal" = "#5f8dff",
  "Ovaries_Pericentral" = "#d32f2f",
  "Ovaries_Periportal" = "#ff5555"
)

# Lipid-Zonation Association Colors (by gonad x zonation)
# For plots showing lipid associations with zonation-biased genes
# Combines gonad (red=Ovaries, blue=Testes) with zonation (dark=pericentral, light=periportal)
LIPID_ZONATION_ASSOC_COLORS <- c(
  # Testes (blue tones)
  "TPC" = "#2b64d6",            # Testes Pericentral - Dark blue
  "TPP" = "#5f8dff",            # Testes Periportal - Light blue
  # Ovaries (red tones)
  "OPC" = "#d32f2f",            # Ovaries Pericentral - Dark red
  "OPP" = "#ff5555",            # Ovaries Periportal - Light red
  # Mixed/None
  "both" = "#808080",           # Grey for both pericentral and periportal
  "no_association" = "#FFFFFF"  # White for no association
)

# Border colors for lipid-zonation associations (darker versions)
LIPID_ZONATION_ASSOC_BORDER_COLORS <- c(
  "TPC" = "#1c4da8",            # Darker blue
  "TPP" = "#3d6fd6",            # Medium blue
  "OPC" = "#a62424",            # Darker red
  "OPP" = "#d63939",            # Medium red
  "both" = "#505050",           # Darker grey
  "no_association" = "#cccccc"  # Light grey border for white fill
)


# Diverging Color Palettes
# =========================
# For heatmaps with negative and positive values

# Red-White-Blue palette (standard)
DIVERGING_PALETTE_RWB <- c("#e54141", "white", "#4173d8")

# Red-Grey-Blue palette (alternative)
DIVERGING_PALETTE_RGB <- c("#e54141", "grey90", "#4173d8")


# Plot Dimensions
# ===============
# Standard dimensions for different plot types

PLOT_DIMS <- list(
  volcano = list(width = 9, height = 9),
  heatmap_small = list(width = 8, height = 6),
  heatmap_medium = list(width = 10, height = 8),
  heatmap_large = list(width = 11, height = 7.5),
  pca = list(width = 7, height = 6),
  umap = list(width = 7, height = 6),
  bubble = list(width = 10, height = 8),
  pathway_heatmap = list(width = 11, height = 7.5),
  tf_heatmap = list(width = 11, height = 7.5)
)


# ggplot2 Themes
# ==============

#' Standard theme for publication-quality plots
#'
#' @param base_size Numeric. Base font size. Default: 12
#' @param base_family Character. Font family. Default: ""
#'
#' @return ggplot2 theme object
theme_publication <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Grid
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),

      # Axis
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 2, face = "bold"),

      # Legend
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1, face = "bold"),

      # Plot title
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 2, hjust = 0.5),

      # Strip (for facets)
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = base_size, face = "bold")
    )
}


#' Minimal theme for clean plots
#'
#' @param base_size Numeric. Base font size. Default: 12
#'
#' @return ggplot2 theme object
theme_clean <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey95", size = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.3),
      legend.position = "right"
    )
}


# Helper Functions
# ================

#' Get color palette for a specific condition type
#'
#' @param condition_type Character. Type of condition: "gonadal", "allosome",
#'                       "group", "celltype", "zonation", "zonation_gonad",
#'                       "lipid_zonation_assoc"
#' @param alternative Logical. Use alternative color scheme? Default: FALSE
#'
#' @return Named vector of colors
#'
#' @examples
#' get_color_palette("gonadal")
#' get_color_palette("allosome", alternative = TRUE)
#' get_color_palette("zonation_gonad")  # TPC, TPP, OPC, OPP
#' get_color_palette("lipid_zonation_assoc")  # pericentral, periportal, both, no_association
get_color_palette <- function(condition_type, alternative = FALSE) {

  color_palette <- switch(
    condition_type,
    "gonadal" = GONADAL_COLORS,
    "allosome" = if (alternative) ALLOSOME_COLORS_ALT else ALLOSOME_COLORS,
    "group" = GROUP_COLORS,
    "celltype" = CELLTYPE_COLORS,
    "zonation" = ZONATION_COLORS,
    "zonation_gonad" = if (alternative) ZONATION_GONAD_COLORS_FULL else ZONATION_GONAD_COLORS,
    "lipid_zonation_assoc" = LIPID_ZONATION_ASSOC_COLORS,
    "lipid_zonation_assoc_border" = LIPID_ZONATION_ASSOC_BORDER_COLORS,
    stop("Unknown condition type: ", condition_type)
  )

  return(color_palette)
}


#' Get plot dimensions for a specific plot type
#'
#' @param plot_type Character. Type of plot: "volcano", "heatmap_small",
#'                  "pca", "umap", etc.
#'
#' @return List with width and height
#'
#' @examples
#' dims <- get_plot_dims("volcano")
#' ggsave("plot.pdf", width = dims$width, height = dims$height)
get_plot_dims <- function(plot_type) {

  if (!plot_type %in% names(PLOT_DIMS)) {
    warning("Unknown plot type: ", plot_type, ". Using default 8x6")
    return(list(width = 8, height = 6))
  }

  return(PLOT_DIMS[[plot_type]])
}


#' Create diverging color palette for heatmaps
#'
#' Wrapper around createLinearColors from support_functions.R with
#' standardized palette.
#'
#' @param values Numeric vector. Values to create palette for
#' @param palette Character vector. Colors to use. Default: DIVERGING_PALETTE_RWB
#' @param withZero Logical. Include zero? Default: FALSE
#' @param maximum Numeric. Maximum value for scaling. Default: auto-calculated
#'
#' @return Character vector of colors
#'
#' @examples
#' palette <- create_heatmap_palette(my_matrix)
create_heatmap_palette <- function(values,
                                   palette = DIVERGING_PALETTE_RWB,
                                   withZero = FALSE,
                                   maximum = NULL) {

  # Source support functions if not already loaded
  if (!exists("createLinearColors")) {
    source("scripts/support_functions.R")
  }

  # Calculate maximum if not provided
  if (is.null(maximum)) {
    maximum <- max(abs(values), na.rm = TRUE) * 10
  }

  # Separate negative and positive values
  neg_values <- values[values < 0]
  pos_values <- values[values > 0]

  # Create palettes
  if (length(neg_values) > 0) {
    palette1 <- createLinearColors(
      neg_values,
      withZero = withZero,
      maximum = abs(min(values, na.rm = TRUE)) * 10,
      my_colors = palette[c(1, 2)]
    )
  } else {
    palette1 <- NULL
  }

  if (length(pos_values) > 0) {
    palette2 <- createLinearColors(
      pos_values,
      withZero = withZero,
      maximum = abs(max(values, na.rm = TRUE)) * 10,
      my_colors = palette[c(2, 3)]
    )
  } else {
    palette2 <- NULL
  }

  # Combine
  full_palette <- c(palette1, palette2)

  return(full_palette)
}


#' Apply standard ggplot2 theme
#'
#' Convenience function to quickly apply publication theme and colors.
#'
#' @param plot ggplot2 object
#' @param condition_type Character. Type of condition for colors.
#' @param add_theme Logical. Add theme_publication()? Default: TRUE
#'
#' @return Modified ggplot2 object
#'
#' @examples
#' p <- ggplot(data, aes(x = x, y = y, color = gonadal_type)) + geom_point()
#' p <- apply_standard_theme(p, "gonadal")
apply_standard_theme <- function(plot, condition_type = NULL, add_theme = TRUE) {

  if (add_theme) {
    plot <- plot + theme_publication()
  }

  if (!is.null(condition_type)) {
    colors <- get_color_palette(condition_type)
    plot <- plot + scale_color_manual(values = colors) +
                  scale_fill_manual(values = colors)
  }

  return(plot)
}


#' Create color mapping for cell types
#'
#' Returns a named vector of colors for specified cell types.
#' If cell types not in standard list, generates colors.
#'
#' @param cell_types Character vector. Cell type names.
#'
#' @return Named vector of colors
#'
#' @examples
#' colors <- create_celltype_colors(c("Hepatocytes", "Endothelial cells"))
create_celltype_colors <- function(cell_types) {

  # Use standard colors for known cell types
  colors <- CELLTYPE_COLORS[cell_types]

  # Generate colors for unknown cell types
  unknown <- cell_types[!cell_types %in% names(CELLTYPE_COLORS)]
  if (length(unknown) > 0) {
    library(RColorBrewer)
    new_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unknown))
    names(new_colors) <- unknown
    colors <- c(colors[!is.na(colors)], new_colors)
  }

  return(colors)
}


#' Clean pathway names for plotting
#'
#' Removes prefixes and cleans up pathway names for better display.
#'
#' @param pathway_names Character vector. Pathway names.
#'
#' @return Character vector of cleaned names
#'
#' @examples
#' clean_names <- clean_pathway_names(c("GOBP_LIPID_METABOLIC_PROCESS"))
clean_pathway_names <- function(pathway_names) {

  cleaned <- pathway_names

  # Remove GO prefixes
  cleaned <- gsub("GO[CMB][CFP]_", "", cleaned)

  # Replace underscores with spaces
  cleaned <- gsub("_", " ", cleaned)

  # Convert to lowercase
  cleaned <- tolower(cleaned)

  # Capitalize first letter
  cleaned <- paste0(toupper(substring(cleaned, 1, 1)), substring(cleaned, 2))

  return(cleaned)
}


# Export color schemes for easy access
# =====================================

#' List all available color schemes
#'
#' Prints all defined color schemes for reference.
list_color_schemes <- function() {
  cat("Available Color Schemes:\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  cat("Cell Type Colors:\n")
  print(CELLTYPE_COLORS)
  cat("\n")

  cat("Gonadal Type Colors:\n")
  print(GONADAL_COLORS)
  cat("\n")

  cat("Allosome Type Colors:\n")
  print(ALLOSOME_COLORS)
  cat("\n")

  cat("FCG Group Colors:\n")
  print(GROUP_COLORS)
  cat("\n")

  cat("Zonation Colors (simple):\n")
  print(ZONATION_COLORS)
  cat("\n")

  cat("Zonation x Gonad Colors (harmonized):\n")
  cat("  Dark shades = Pericentral, Light shades = Periportal\n")
  cat("  Blue tones = Testes, Red tones = Ovaries\n")
  print(ZONATION_GONAD_COLORS)
  cat("\n")

  cat("Lipid-Zonation Association Colors:\n")
  print(LIPID_ZONATION_ASSOC_COLORS)
  cat("\n")

  cat(paste(rep("=", 80), collapse = ""), "\n")
}
