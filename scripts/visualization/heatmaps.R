#' Heatmap Visualization Module
#'
#' This module provides functions for creating publication-quality heatmaps
#' including:
#' - Differential expression heatmaps
#' - Pathway/TF enrichment heatmaps
#' - Correlation heatmaps
#' - Custom color palettes (diverging, sequential)
#'
#' These functions reproduce the exact plotting behavior from support_functions.R
#' and are used to generate manuscript figures.

library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)


#' Create Linear Diverging Color Palette
#'
#' Creates a diverging color palette for heatmaps based on data values.
#' Automatically handles positive-only, negative-only, or mixed values.
#'
#' This function replicates createLinearColors() from support_functions.R.
#'
#' @param numbers Numeric vector. Values to create palette for.
#' @param withZero Logical. Include zero in palette? (default: TRUE)
#' @param maximum Numeric. Maximum palette resolution (default: 100).
#' @param my_colors Character vector of 3 colors: c(low, mid, high).
#'   Default: c("red", "white", "royalblue3")
#'
#' @return Character vector of colors
#'
#' @examples
#' # For t-statistics
#' palette <- createLinearColors(t_values)
#'
#' # Custom colors
#' palette <- createLinearColors(logFC, my_colors = c("blue", "white", "red"))
#'
#' @export
createLinearColors <- function(numbers, withZero = TRUE, maximum = 100,
                               my_colors = c("red", "white", "royalblue3")) {

  if (min(numbers, na.rm = TRUE) > 0) {
    # All positive values
    if (withZero) {
      numbers <- c(0, numbers)
    }
    myPalette <- colorRampPalette(my_colors[c(2, 3)])
    myColors <- myPalette(maximum)
  } else {
    if (max(numbers, na.rm = TRUE) < 0) {
      # All negative values
      if (withZero) {
        numbers <- c(0, numbers)
      }
      myPalette <- colorRampPalette(my_colors[c(1, 2)])
      myColors <- myPalette(maximum)
    } else {
      # Mixed positive and negative values
      myPalette_pos <- colorRampPalette(c("white", "royalblue3"))
      myPalette_neg <- colorRampPalette(c("red", "white"))

      npos <- length(numbers[numbers >= 0]) + 1
      nneg <- length(numbers[numbers <= 0]) + 1

      myColors_pos <- myPalette_pos(npos)
      myColors_neg <- myPalette_neg(nneg)

      myColors <- c(myColors_neg[-(nneg)], myColors_pos[-1])
    }
  }

  return(myColors)
}


#' Create Expression Heatmap
#'
#' Creates a heatmap of gene expression or other measurements with clustering.
#'
#' @param data Matrix or data frame. Features x Samples.
#' @param annotation_col Data frame. Column annotations (samples).
#' @param annotation_row Data frame. Row annotations (features).
#' @param annotation_colors List. Custom colors for annotations.
#' @param scale Character. Scale rows ("row"), columns ("column"), or "none" (default: "row").
#' @param cluster_rows Logical. Cluster rows? (default: TRUE)
#' @param cluster_cols Logical. Cluster columns? (default: TRUE)
#' @param show_rownames Logical. Show row names? (default: FALSE)
#' @param show_colnames Logical. Show column names? (default: TRUE)
#' @param color_palette Character vector. Colors for heatmap.
#'   If NULL, uses red-white-blue diverging palette.
#' @param breaks Numeric vector. Color breaks. If NULL, auto-calculated.
#' @param main Character. Plot title.
#' @param fontsize Numeric. Base font size (default: 10).
#' @param filename Character. If provided, saves plot to file.
#' @param width Numeric. Plot width in inches (default: 8).
#' @param height Numeric. Plot height in inches (default: 10).
#' @param ... Additional arguments passed to pheatmap().
#'
#' @return pheatmap object
#'
#' @examples
#' # Basic expression heatmap
#' create_expression_heatmap(scaled_RNA[top_genes, ])
#'
#' # With annotations
#' annot <- data.frame(
#'   gonadal_type = metadata$gonadal_type,
#'   row.names = metadata$sample
#' )
#' colors <- list(gonadal_type = c("Ovaries" = "red", "Testes" = "blue"))
#' create_expression_heatmap(data, annotation_col = annot,
#'                          annotation_colors = colors)
#'
#' @export
create_expression_heatmap <- function(data, annotation_col = NULL,
                                      annotation_row = NULL,
                                      annotation_colors = NULL,
                                      scale = "row",
                                      cluster_rows = TRUE,
                                      cluster_cols = TRUE,
                                      show_rownames = FALSE,
                                      show_colnames = TRUE,
                                      color_palette = NULL,
                                      breaks = NULL,
                                      main = "",
                                      fontsize = 10,
                                      filename = NULL,
                                      width = 8,
                                      height = 10,
                                      ...) {

  # Create color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("red", "white", "royalblue3"))(100)
  }

  # Create breaks if not provided
  if (is.null(breaks)) {
    if (scale == "row") {
      max_val <- max(abs(data), na.rm = TRUE)
      breaks <- seq(-max_val, max_val, length.out = length(color_palette) + 1)
    } else {
      breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE),
                    length.out = length(color_palette) + 1)
    }
  }

  # Create heatmap using explicit namespace to avoid ComplexHeatmap masking
  p <- pheatmap::pheatmap(
    data,
    color = color_palette,
    breaks = breaks,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    main = main,
    fontsize = fontsize,
    border_color = NA,
    filename = filename,
    width = width,
    height = height,
    ...
  )

  return(p)
}


#' Create Enrichment Heatmap
#'
#' Creates a heatmap of pathway/TF enrichment scores with display_numbers option.
#'
#' @param enrichment_matrix Matrix or data frame. Pathways/TFs x Contrasts.
#' @param display_numbers Logical or matrix. Display values in cells? (default: FALSE)
#' @param number_format Character. Format string for numbers (default: "%.1f").
#' @param number_color Character. Color for numbers (default: "black").
#' @param cluster_rows Logical. Cluster rows? (default: TRUE)
#' @param cluster_cols Logical. Cluster columns? (default: FALSE)
#' @param color_palette Character vector. Colors. If NULL, red-grey-blue.
#' @param breaks Numeric vector. Color breaks.
#' @param angle_col Numeric. Angle for column labels (default: 45).
#' @param fontsize Numeric. Base font size (default: 10).
#' @param fontsize_number Numeric. Number font size (default: 8).
#' @param main Character. Plot title.
#' @param filename Character. Output file path.
#' @param width Numeric. Width in inches (default: 11).
#' @param height Numeric. Height in inches (default: 7.5).
#' @param ... Additional arguments for pheatmap().
#'
#' @return pheatmap object
#'
#' @examples
#' # Pathway enrichment heatmap
#' create_enrichment_heatmap(pathway_scores, display_numbers = TRUE,
#'                          main = "Pathway Enrichment")
#'
#' # TF activity heatmap
#' create_enrichment_heatmap(tf_scores, display_numbers = TRUE,
#'                          number_format = "%.2f")
#'
#' @export
create_enrichment_heatmap <- function(enrichment_matrix, display_numbers = FALSE,
                                      number_format = "%.1f",
                                      number_color = "black",
                                      cluster_rows = TRUE,
                                      cluster_cols = FALSE,
                                      color_palette = NULL,
                                      breaks = NULL,
                                      angle_col = 45,
                                      fontsize = 10,
                                      fontsize_number = 8,
                                      main = "",
                                      filename = NULL,
                                      width = 11,
                                      height = 7.5,
                                      ...) {

  # Create color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("red", "grey90", "royalblue3"))(100)
  }

  # Create breaks if not provided
  if (is.null(breaks)) {
    max_val <- max(abs(enrichment_matrix), na.rm = TRUE)
    breaks <- seq(-max_val, max_val, length.out = length(color_palette) + 1)
  }

  # Check if using ComplexHeatmap's pheatmap (it doesn't support filename directly)
  # ComplexHeatmap's pheatmap also has stricter requirements for angle_col
  using_complexheatmap <- "ComplexHeatmap" %in% loadedNamespaces()

  # Validate input matrix
  if (nrow(enrichment_matrix) == 0 || ncol(enrichment_matrix) == 0) {
    warning("Cannot create heatmap: input matrix is empty (",
            nrow(enrichment_matrix), " rows x ",
            ncol(enrichment_matrix), " cols)")
    return(NULL)
  }

  # IMPORTANT: Always use pheatmap:: to avoid ComplexHeatmap masking the function
  # ComplexHeatmap::pheatmap behaves differently and doesn't save files properly

  # Create directory if needed
  if (!is.null(filename)) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  }

  # Use original pheatmap with explicit namespace
  # This works correctly regardless of what packages are loaded
  tryCatch({
    p <- pheatmap::pheatmap(
      enrichment_matrix,
      color = color_palette,
      breaks = breaks,
      display_numbers = display_numbers,
      number_format = number_format,
      number_color = number_color,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      angle_col = angle_col,
      fontsize = fontsize,
      fontsize_number = fontsize_number,
      border_color = NA,
      main = main,
      filename = filename,
      width = width,
      height = height,
      ...
    )

    # Verify file was created successfully if filename provided
    if (!is.null(filename)) {
      if (file.exists(filename)) {
        file_size <- file.info(filename)$size
        if (file_size < 1000) {
          warning("Heatmap file created but is suspiciously small (", file_size, " bytes)")
        }
      } else {
        warning("Heatmap file was not created: ", filename)
      }
    }

    return(p)
  }, error = function(e) {
    if (!is.null(filename) && file.exists(filename)) {
      unlink(filename)  # Remove corrupted file
    }
    stop("Error creating enrichment heatmap: ", conditionMessage(e))
  })
}


#' Create Correlation Heatmap
#'
#' Creates a correlation heatmap with optional hierarchical clustering.
#'
#' @param data Matrix or data frame. Features x Samples or Samples x Features.
#' @param method Character. Correlation method: "pearson", "spearman", "kendall" (default: "pearson").
#' @param cluster_rows Logical. Cluster rows? (default: TRUE)
#' @param cluster_cols Logical. Cluster columns? (default: TRUE)
#' @param show_rownames Logical. Show row names? (default: TRUE)
#' @param show_colnames Logical. Show column names? (default: TRUE)
#' @param color_palette Character vector. Colors. If NULL, uses blue-white-red.
#' @param main Character. Plot title.
#' @param fontsize Numeric. Base font size (default: 10).
#' @param filename Character. Output file path.
#' @param width Numeric. Width in inches (default: 8).
#' @param height Numeric. Height in inches (default: 8).
#' @param ... Additional arguments for pheatmap().
#'
#' @return pheatmap object
#'
#' @examples
#' # Sample correlation heatmap
#' create_correlation_heatmap(scaled_RNA, main = "Sample Correlations")
#'
#' # Spearman correlation
#' create_correlation_heatmap(lipidomic, method = "spearman")
#'
#' @export
create_correlation_heatmap <- function(data, method = "pearson",
                                       cluster_rows = TRUE,
                                       cluster_cols = TRUE,
                                       show_rownames = TRUE,
                                       show_colnames = TRUE,
                                       color_palette = NULL,
                                       main = "Correlation",
                                       fontsize = 10,
                                       filename = NULL,
                                       width = 8,
                                       height = 8,
                                       ...) {

  # Compute correlation matrix
  cor_matrix <- cor(data, method = method, use = "pairwise.complete.obs")

  # Create color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("royalblue3", "white", "red"))(100)
  }

  # Create breaks
  breaks <- seq(-1, 1, length.out = length(color_palette) + 1)

  # Create heatmap using explicit namespace to avoid ComplexHeatmap masking
  p <- pheatmap::pheatmap(
    cor_matrix,
    color = color_palette,
    breaks = breaks,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    display_numbers = FALSE,
    border_color = NA,
    main = main,
    fontsize = fontsize,
    filename = filename,
    width = width,
    height = height,
    ...
  )

  return(p)
}


#' Create ComplexHeatmap with Advanced Features
#'
#' Creates a ComplexHeatmap object with multiple annotations and legends.
#'
#' @param data Matrix. Features x Samples.
#' @param top_annotation HeatmapAnnotation object. Column annotations.
#' @param left_annotation HeatmapAnnotation object. Row annotations.
#' @param col_fun colorRamp2 function. Color mapping function.
#' @param cluster_rows Logical or hclust object. Cluster rows? (default: TRUE)
#' @param cluster_columns Logical or hclust object. Cluster columns? (default: TRUE)
#' @param show_row_names Logical. Show row names? (default: FALSE)
#' @param show_column_names Logical. Show column names? (default: TRUE)
#' @param row_names_gp gpar object. Row name graphics parameters.
#' @param column_names_gp gpar object. Column name graphics parameters.
#' @param heatmap_legend_param List. Legend parameters.
#' @param column_title Character. Column title.
#' @param row_title Character. Row title.
#' @param ... Additional arguments for Heatmap().
#'
#' @return Heatmap object
#'
#' @examples
#' # Create column annotation
#' col_annot <- HeatmapAnnotation(
#'   gonadal_type = metadata$gonadal_type,
#'   col = list(gonadal_type = c("Ovaries" = "red", "Testes" = "blue"))
#' )
#'
#' # Create heatmap
#' ht <- create_complex_heatmap(
#'   data = scaled_RNA[top_genes, ],
#'   top_annotation = col_annot,
#'   column_title = "Gonad-Dependent Genes"
#' )
#' draw(ht)
#'
#' @export
create_complex_heatmap <- function(data, top_annotation = NULL,
                                   left_annotation = NULL,
                                   col_fun = NULL,
                                   cluster_rows = TRUE,
                                   cluster_columns = TRUE,
                                   show_row_names = FALSE,
                                   show_column_names = TRUE,
                                   row_names_gp = gpar(fontsize = 8),
                                   column_names_gp = gpar(fontsize = 8),
                                   heatmap_legend_param = list(),
                                   column_title = "",
                                   row_title = "",
                                   ...) {

  # Create color function if not provided
  if (is.null(col_fun)) {
    max_val <- max(abs(data), na.rm = TRUE)
    col_fun <- colorRamp2(
      c(-max_val, 0, max_val),
      c("red", "white", "royalblue3")
    )
  }

  # Create heatmap
  ht <- Heatmap(
    data,
    col = col_fun,
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    heatmap_legend_param = heatmap_legend_param,
    column_title = column_title,
    row_title = row_title,
    border = TRUE,
    ...
  )

  return(ht)
}


#' Create Bubble Heatmap
#'
#' Creates a bubble heatmap where bubble size represents significance
#' and color represents effect size (commonly used for enrichment results).
#'
#' @param data Data frame with columns: term, contrast, score, p_value.
#' @param size_column Character. Column name for bubble size (default: "p_value").
#' @param color_column Character. Column name for bubble color (default: "score").
#' @param max_pvalue Numeric. Maximum p-value to display (default: 0.05).
#' @param size_range Numeric vector. c(min, max) for bubble size (default: c(2, 10)).
#' @param color_palette Character vector. Colors for gradient.
#' @param term_column Character. Column name for terms (default: "term").
#' @param contrast_column Character. Column name for contrasts (default: "contrast").
#' @param title Character. Plot title.
#'
#' @return ggplot object
#'
#' @examples
#' # Pathway enrichment bubble plot
#' bubble_plot <- create_bubble_heatmap(
#'   enrichment_results,
#'   max_pvalue = 0.01,
#'   title = "Top Enriched Pathways"
#' )
#'
#' @export
create_bubble_heatmap <- function(data, size_column = "p_value",
                                  color_column = "score",
                                  max_pvalue = 0.05,
                                  size_range = c(2, 10),
                                  color_palette = NULL,
                                  term_column = "term",
                                  contrast_column = "contrast",
                                  title = "") {

  # Filter by p-value
  data <- data[data[[size_column]] <= max_pvalue, ]

  # Transform p-value for size (larger bubbles = more significant)
  data$neg_log_p <- -log10(data[[size_column]])

  # Create color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- c("red", "white", "royalblue3")
  }

  # Create plot
  p <- ggplot(data, aes_string(
    x = contrast_column,
    y = term_column,
    size = "neg_log_p",
    color = color_column
  )) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(
      name = "-log10(p-value)",
      range = size_range
    ) +
    scale_color_gradient2(
      low = color_palette[1],
      mid = color_palette[2],
      high = color_palette[3],
      midpoint = 0,
      name = "Score"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = title)

  return(p)
}


#' Save Heatmap to File
#'
#' Saves a pheatmap or ComplexHeatmap object to a file.
#'
#' @param heatmap pheatmap or Heatmap object.
#' @param filename Character. Output file path.
#' @param width Numeric. Width in inches.
#' @param height Numeric. Height in inches.
#' @param device Character. Device: "pdf", "png", "svg" (default: "pdf").
#'
#' @return Invisibly returns filename
#'
#' @examples
#' hm <- create_expression_heatmap(data)
#' save_heatmap(hm, "results/heatmap.pdf", width = 8, height = 10)
#'
#' @export
save_heatmap <- function(heatmap, filename, width, height, device = "pdf") {

  # Determine heatmap type
  is_complex <- inherits(heatmap, "Heatmap")

  if (is_complex) {
    # ComplexHeatmap
    if (device == "pdf") {
      pdf(filename, width = width, height = height)
    } else if (device == "png") {
      png(filename, width = width * 100, height = height * 100, res = 100)
    } else if (device == "svg") {
      svg(filename, width = width, height = height)
    }

    draw(heatmap)
    dev.off()
  } else {
    # pheatmap - already saves to file if filename provided in pheatmap()
    # This function is just for consistency
    cat("Note: For pheatmap objects, use filename parameter in pheatmap() directly\n")
  }

  cat("Heatmap saved to:", filename, "\n")

  invisible(filename)
}


#' Prepare Data for Heatmap
#'
#' Prepares data for heatmap visualization by filtering, scaling, and selecting top features.
#'
#' @param data Matrix or data frame. Features x Samples.
#' @param top_n Integer. Select top N features by variance/MAD.
#' @param selection_method Character. "variance", "mad", or "custom" (default: "variance").
#' @param custom_features Character vector. Custom feature selection (for selection_method = "custom").
#' @param scale_data Character. Scale "row", "column", or "none" (default: "row").
#' @param center Logical. Center data? (default: TRUE)
#'
#' @return Matrix ready for heatmap
#'
#' @examples
#' # Select top 100 most variable genes
#' heatmap_data <- prepare_data_for_heatmap(scaled_RNA, top_n = 100)
#'
#' # Custom gene selection
#' heatmap_data <- prepare_data_for_heatmap(
#'   scaled_RNA,
#'   selection_method = "custom",
#'   custom_features = c("Cyp2b9", "Cyp2c38", "Esr1")
#' )
#'
#' @export
prepare_data_for_heatmap <- function(data, top_n = 50,
                                     selection_method = "variance",
                                     custom_features = NULL,
                                     scale_data = "row",
                                     center = TRUE) {

  # Remove incomplete cases
  data <- data[complete.cases(data), ]

  # Feature selection
  if (selection_method == "custom") {
    if (is.null(custom_features)) {
      stop("custom_features must be provided when selection_method = 'custom'")
    }
    selected_features <- custom_features[custom_features %in% rownames(data)]
    data_selected <- data[selected_features, , drop = FALSE]

  } else if (selection_method == "variance") {
    # Select by variance
    variances <- apply(data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE)[1:min(top_n, length(variances))])
    data_selected <- data[top_features, , drop = FALSE]

  } else if (selection_method == "mad") {
    # Select by median absolute deviation
    mads <- apply(data, 1, mad, na.rm = TRUE)
    top_features <- names(sort(mads, decreasing = TRUE)[1:min(top_n, length(mads))])
    data_selected <- data[top_features, , drop = FALSE]

  } else {
    stop("selection_method must be 'variance', 'mad', or 'custom'")
  }

  # Scaling
  if (scale_data == "row") {
    data_selected <- t(scale(t(data_selected), center = center, scale = TRUE))
  } else if (scale_data == "column") {
    data_selected <- scale(data_selected, center = center, scale = TRUE)
  }

  return(data_selected)
}
