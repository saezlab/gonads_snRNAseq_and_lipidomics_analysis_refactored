#' Dimensionality Reduction Plots Module
#'
#' This module provides functions for PCA, UMAP, and other dimensionality reduction
#' visualizations.
#'
#' The nicePCA() function creates a comprehensive 3x3 grid showing:
#' - Component density distributions
#' - Pairwise component scatter plots
#' - Scree plot (explained variance)
#' - Legend
#'
#' These functions reproduce the exact plotting behavior from support_functions.R
#' and are used to generate manuscript figures.

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(cowplot)


#' Create Comprehensive PCA Visualization Grid
#'
#' Creates a 3x3 grid of PCA plots including:
#' - Density distributions for each component
#' - Pairwise scatter plots (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3)
#' - Scree plot showing explained variance
#' - Legend
#'
#' This function replicates the exact behavior of nicePCA() from support_functions.R.
#'
#' @param df Data frame or matrix. Features x Samples (rows x columns).
#' @param targets Data frame with 2 columns: sample, condition (or sample, color).
#'   If "color" column present, uses it for coloring; otherwise uses "condition".
#' @param components Integer vector of length 3. Components to plot (default: c(1, 2, 3)).
#' @param centering Logical. Center data before PCA? (default: TRUE)
#' @param scaling Logical. Scale data before PCA? (default: FALSE)
#' @param pointSize Numeric. Point size for scatter plots (default: 4).
#' @param no_label Logical. Omit sample labels? (default: FALSE)
#' @param color_mapping Named vector. Custom color mapping for conditions.
#'
#' @return arrangeGrob object (3x3 grid) that can be displayed with grid.draw() or saved.
#'
#' @examples
#' # Basic PCA plot
#' pca_grid <- nicePCA(scaled_RNA, targets)
#' grid.draw(pca_grid)
#'
#' # Without labels (cleaner for many samples)
#' pca_grid <- nicePCA(scaled_RNA, targets, no_label = TRUE)
#' grid.draw(pca_grid)
#'
#' # With custom colors
#' colors <- c("Ovaries" = "red", "Testes" = "royalblue")
#' pca_grid <- nicePCA(lipidomic, targets, color_mapping = colors)
#' grid.draw(pca_grid)
#'
#' # Specific components
#' pca_grid <- nicePCA(data, targets, components = c(1, 3, 4))
#' grid.draw(pca_grid)
#'
#' @export
nicePCA <- function(df, targets, components = c(1, 2, 3), centering = TRUE,
                    scaling = FALSE, pointSize = 4, no_label = FALSE,
                    color_mapping = NULL) {

  # Match and align df and targets
  if (!is.null(targets)) {
    df_and_targets <- make_df_and_targets_great_again(df, targets)
    df <- df_and_targets[[1]]
    targets <- df_and_targets[[2]]
  }

  # Perform PCA
  sample <- names(df)
  pca <- prcomp(t(df[complete.cases(df), ]), center = centering, scale. = scaling)
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)

  # Component labels with explained variance
  xCompLab <- paste("PC", components[1], " (", sep = "")
  yCompLab <- paste("PC", components[2], " (", sep = "")
  zCompLab <- paste("PC", components[3], " (", sep = "")

  # Prepare data for plotting
  data.to.plot <- as.data.frame(cbind(
    pca$x[, components[1]],
    pca$x[, components[2]],
    pca$x[, components[3]]
  ))

  # Use "color" column if present, otherwise "condition"
  if ("color" %in% names(targets)) {
    data.to.plot$condition <- targets$color
  } else {
    data.to.plot$condition <- targets$condition
  }
  data.to.plot$sample <- targets$sample
  colnames(data.to.plot) <- c("pc1", "pc2", "pc3", "condition", "sample")

  # Print for debugging
  print(head(data.to.plot))

  # Create density plots for each component
  # ========================================

  comp1 <- ggplot(data.to.plot, aes(x = pc1, color = condition)) +
    geom_density() +
    theme_classic() +
    xlab(paste(xCompLab, round(100 * explained[components[1]], digits = 2), '%)', sep = '')) +
    theme(legend.position = "none") +
    scale_x_continuous(
      position = "top",
      limits = c(-max(abs(pca$x[, components[1]])), max(abs(pca$x[, components[1]])))
    ) +
    scale_y_continuous(position = "left")

  comp2 <- ggplot(data.to.plot, aes(x = pc2, color = condition)) +
    geom_density() +
    theme_classic() +
    xlab(paste(yCompLab, round(100 * explained[components[2]], digits = 2), '%)', sep = '')) +
    theme(legend.position = "none") +
    scale_x_continuous(
      position = "top",
      limits = c(-max(abs(pca$x[, components[2]])), max(abs(pca$x[, components[2]])))
    ) +
    scale_y_continuous(position = "right")

  comp3 <- ggplot(data.to.plot, aes(x = pc3, color = condition)) +
    geom_density() +
    theme_classic() +
    xlab(paste(zCompLab, round(100 * explained[components[3]], digits = 2), '%)', sep = '')) +
    theme(legend.position = "none") +
    scale_x_continuous(
      position = "top",
      limits = c(-max(abs(pca$x[, components[3]])), max(abs(pca$x[, components[3]])))
    ) +
    scale_y_continuous(position = "right")

  # Extract legend
  # ==============

  just_for_legend <- ggplot(data.to.plot, aes(x = pc1, y = pc2, color = condition)) +
    geom_point(size = pointSize, alpha = 0.3) +
    geom_text(aes(label = sample), nudge_y = 1.0, size = 6) +
    scale_alpha_discrete(range = c(0.3, 1.0)) +
    theme_classic() +
    xlab(paste(xCompLab, round(100 * explained[components[1]], digits = 2), '%)', sep = '')) +
    ylab(paste(yCompLab, round(100 * explained[components[2]], digits = 2), '%)', sep = '')) +
    xlim(c(-max(abs(pca$x[, components[1]])), max(abs(pca$x[, components[1]])))) +
    ylim(c(-max(abs(pca$x[, components[2]])), max(abs(pca$x[, components[2]]))))

  if (!is.null(color_mapping)) {
    just_for_legend <- just_for_legend + scale_color_manual(values = color_mapping)
  }

  legend <- get_legend(just_for_legend)
  rm(just_for_legend)

  # Create scatter plots
  # ====================

  if (no_label) {
    # Without labels (cleaner for many samples)
    comp1vcomp2 <- ggplot(data.to.plot, aes(x = pc1, y = pc2, color = condition)) +
      geom_point(size = pointSize, alpha = 0.6) +
      scale_alpha_discrete(range = c(0.3, 1.0)) +
      theme_classic() +
      xlab("") +
      ylab(paste(yCompLab, round(100 * explained[components[2]], digits = 2), '%)', sep = '')) +
      ylim(c(-max(abs(pca$x[, components[3]])), max(abs(pca$x[, components[3]])))) +
      # scale_x_continuous(limits = c(-30, 30), position = "bottom") +
      # scale_y_continuous(limits = c(-30, 30), position = "left") +
      theme(legend.position = "none")

    comp1vcomp3 <- ggplot(data.to.plot, aes(x = pc1, y = pc3, color = condition)) +
      geom_point(size = pointSize, alpha = 0.3) +
      scale_alpha_discrete(range = c(0.3, 1.0)) +
      theme_classic() +
      xlab("") +
      ylab(paste(zCompLab, round(100 * explained[components[3]], digits = 2), '%)', sep = '')) +
      ylim(c(-max(abs(pca$x[, components[3]])), max(abs(pca$x[, components[3]])))) +
      theme(legend.position = "none") +
      scale_x_continuous(
        position = "top",
        limits = c(-max(abs(pca$x[, components[1]])), max(abs(pca$x[, components[1]])))
      )

    comp2vcomp3 <- ggplot(data.to.plot, aes(x = pc2, y = pc3, color = condition)) +
      geom_point(size = pointSize, alpha = 0.3) +
      scale_alpha_discrete(range = c(0.3, 1.0)) +
      theme_classic() +
      xlab("") +
      ylab("") +
      scale_y_continuous(
        position = "right",
        limits = c(-max(abs(pca$x[, components[3]])), max(abs(pca$x[, components[3]])))
      ) +
      scale_x_continuous(
        position = "top",
        limits = c(-max(abs(pca$x[, components[2]])), max(abs(pca$x[, components[2]])))
      ) +
      theme(legend.position = "none")
  } else {
    # With labels
    comp1vcomp2 <- ggplot(data.to.plot, aes(x = pc1, y = pc2, color = condition)) +
      geom_point(size = pointSize, alpha = 0.3) +
      geom_text_repel(aes(label = sample), size = 2) +
      scale_alpha_discrete(range = c(0.3, 1.0)) +
      theme_classic() +
      xlab("") +
      ylab(paste(yCompLab, round(100 * explained[components[2]], digits = 2), '%)', sep = '')) +
      ylim(c(-max(abs(pca$x[, components[2]])), max(abs(pca$x[, components[2]])))) +
      theme(legend.position = "none") +
      scale_x_continuous(
        position = "top",
        limits = c(-max(abs(pca$x[, components[1]])), max(abs(pca$x[, components[1]])))
      )

    comp1vcomp3 <- ggplot(data.to.plot, aes(x = pc1, y = pc3, color = condition)) +
      geom_point(size = pointSize, alpha = 0.3) +
      geom_text_repel(aes(label = sample), size = 2) +
      scale_alpha_discrete(range = c(0.3, 1.0)) +
      theme_classic() +
      xlab("") +
      ylab(paste(zCompLab, round(100 * explained[components[3]], digits = 2), '%)', sep = '')) +
      ylim(c(-max(abs(pca$x[, components[3]])), max(abs(pca$x[, components[3]])))) +
      theme(legend.position = "none") +
      scale_x_continuous(
        position = "top",
        limits = c(-max(abs(pca$x[, components[1]])), max(abs(pca$x[, components[1]])))
      )

    comp2vcomp3 <- ggplot(data.to.plot, aes(x = pc2, y = pc3, color = condition)) +
      geom_point(size = pointSize, alpha = 0.3) +
      geom_text_repel(aes(label = sample), size = 2) +
      scale_alpha_discrete(range = c(0.3, 1.0)) +
      theme_classic() +
      xlab("") +
      ylab("") +
      scale_y_continuous(
        position = "right",
        limits = c(-max(abs(pca$x[, components[3]])), max(abs(pca$x[, components[3]])))
      ) +
      scale_x_continuous(
        position = "top",
        limits = c(-max(abs(pca$x[, components[2]])), max(abs(pca$x[, components[2]])))
      ) +
      theme(legend.position = "none")
  }

  # Create scree plot (boulder plot)
  # =================================

  explained_df <- as.data.frame(explained)
  explained_df$Components <- as.numeric(rownames(explained_df))
  names(explained_df)[1] <- "% of explained variance"

  boulder <- ggplot(explained_df, aes(x = Components, y = `% of explained variance`, fill = "red")) +
    geom_col() +
    theme_classic() +
    scale_y_continuous(position = "right") +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(position = "top", breaks = seq_along(explained_df[, 1]))

  # Apply color mapping if provided
  # ================================

  if (!is.null(color_mapping)) {
    comp1 <- comp1 + scale_color_manual(values = color_mapping)
    comp2 <- comp2 + scale_color_manual(values = color_mapping)
    comp3 <- comp3 + scale_color_manual(values = color_mapping)
    comp1vcomp2 <- comp1vcomp2 + scale_color_manual(values = color_mapping)
    comp1vcomp3 <- comp1vcomp3 + scale_color_manual(values = color_mapping)
    comp2vcomp3 <- comp2vcomp3 + scale_color_manual(values = color_mapping)
  }

  # Arrange in 3x3 grid
  # ===================

  return(arrangeGrob(
    comp1, legend, boulder,
    comp1vcomp2, comp2, legend,
    comp1vcomp3, comp2vcomp3, comp3,
    ncol = 3
  ))
}


#' Match and Align Data and Targets
#'
#' Ensures df and targets have matching samples in the same order.
#' This is a helper function used by nicePCA().
#'
#' @param df Data frame. Features x Samples.
#' @param targets Data frame with columns: sample, condition.
#'
#' @return List with aligned df and targets.
#'
#' @keywords internal
make_df_and_targets_great_again <- function(df, targets)
{
  names(targets)[c(1,2)] <- c("sample","condition")
  bad_samples <- c(NA)
  i <- 1
  for (sample in names(df))
  {
    if (!(sample %in% targets[,1]))
    {
      bad_samples <- c(bad_samples, sample)
    }
    else
    {
      if (i == 1)
      {
        clean_targets <- as.data.frame(matrix(NA, 1, length(targets[1,])))
        names(clean_targets) <- names(targets)
        clean_targets[1,] <- targets[targets$sample == sample,]
        i <- i+1
      }
      else
      {
        clean_targets <- as.data.frame(rbind(clean_targets,targets[targets$sample == sample,]))
        i <- i+1
      }
    }
  }
  targets <- clean_targets
  if (length(bad_samples) > 1)
  {
    bad_samples <- bad_samples[-1]
    print(paste("These samples were not found : ", bad_samples, sep = ""))
    print(bad_samples)
    df <- df[,!(names(df) %in% bad_samples)]
  }
  return(list(df,targets))
}


#' Create Simple PCA Plot (Single Panel)
#'
#' Creates a simple 2D PCA plot (PC1 vs PC2) with custom styling.
#'
#' @param df Data frame or matrix. Features x Samples.
#' @param targets Data frame with columns: sample, condition.
#' @param components Integer vector of length 2. Components to plot (default: c(1, 2)).
#' @param centering Logical. Center data? (default: TRUE)
#' @param scaling Logical. Scale data? (default: FALSE)
#' @param pointSize Numeric. Point size (default: 3).
#' @param add_labels Logical. Add sample labels? (default: FALSE)
#' @param color_mapping Named vector. Custom colors.
#' @param title Character. Plot title.
#'
#' @return ggplot object
#'
#' @examples
#' # Simple PCA plot
#' p <- simple_pca(scaled_RNA, targets, title = "Hepatocytes")
#'
#' # With custom colors
#' colors <- c("Ovaries" = "red", "Testes" = "royalblue")
#' p <- simple_pca(data, targets, color_mapping = colors)
#'
#' @export
simple_pca <- function(df, targets, components = c(1, 2), centering = TRUE,
                       scaling = FALSE, pointSize = 3, add_labels = FALSE,
                       color_mapping = NULL, title = NULL) {

  # Match and align
  if (!is.null(targets)) {
    df_and_targets <- make_df_and_targets_great_again(df, targets)
    df <- df_and_targets[[1]]
    targets <- df_and_targets[[2]]
  }

  # Perform PCA
  pca <- prcomp(t(df[complete.cases(df), ]), center = centering, scale. = scaling)
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)

  # Prepare data
  data.to.plot <- as.data.frame(cbind(
    pca$x[, components[1]],
    pca$x[, components[2]]
  ))

  if ("color" %in% names(targets)) {
    data.to.plot$condition <- targets$color
  } else {
    data.to.plot$condition <- targets$condition
  }
  data.to.plot$sample <- targets$sample
  colnames(data.to.plot) <- c("pc1", "pc2", "condition", "sample")

  # Create plot
  p <- ggplot(data.to.plot, aes(x = pc1, y = pc2, color = condition)) +
    geom_point(size = pointSize, alpha = 0.6) +
    theme_classic() +
    xlab(paste0("PC", components[1], " (",
                round(100 * explained[components[1]], 1), "%)")) +
    ylab(paste0("PC", components[2], " (",
                round(100 * explained[components[2]], 1), "%)")) +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  if (add_labels) {
    p <- p + geom_text_repel(aes(label = sample), size = 3)
  }

  if (!is.null(color_mapping)) {
    p <- p + scale_color_manual(values = color_mapping)
  }

  return(p)
}


#' Save PCA Grid to File
#'
#' Saves a nicePCA grid to a file.
#'
#' @param pca_grid arrangeGrob object from nicePCA().
#' @param filename Character. Output filename.
#' @param width Numeric. Width in inches (default: 12).
#' @param height Numeric. Height in inches (default: 12).
#'
#' @return Invisibly returns filename
#'
#' @examples
#' pca_grid <- nicePCA(scaled_RNA, targets)
#' save_pca_grid(pca_grid, "results/pca_hepatocytes.pdf")
#'
#' @export
save_pca_grid <- function(pca_grid, filename, width = 12, height = 12) {
  ggsave(
    filename = filename,
    plot = pca_grid,
    width = width,
    height = height
  )

  cat("PCA grid saved to:", filename, "\n")

  invisible(filename)
}


#' Extract PCA Loadings
#'
#' Extracts feature loadings from PCA for a specific component.
#'
#' @param df Data frame or matrix. Features x Samples.
#' @param component Integer. Component number (default: 1).
#' @param centering Logical. Center data? (default: TRUE)
#' @param scaling Logical. Scale data? (default: FALSE)
#' @param top_n Integer. Return only top N features by absolute loading?
#'
#' @return Data frame with feature names and loadings, sorted by absolute loading.
#'
#' @examples
#' # Get top 50 features driving PC1
#' loadings <- extract_pca_loadings(scaled_RNA, component = 1, top_n = 50)
#'
#' @export
extract_pca_loadings <- function(df, component = 1, centering = TRUE,
                                  scaling = FALSE, top_n = NULL) {

  # Perform PCA
  pca <- prcomp(t(df[complete.cases(df), ]), center = centering, scale. = scaling)

  # Extract loadings (rotation matrix)
  loadings <- pca$rotation[, component]

  # Create data frame
  loadings_df <- data.frame(
    feature = names(loadings),
    loading = loadings,
    abs_loading = abs(loadings),
    stringsAsFactors = FALSE
  )

  # Sort by absolute loading
  loadings_df <- loadings_df[order(loadings_df$abs_loading, decreasing = TRUE), ]

  # Filter to top N if requested
  if (!is.null(top_n)) {
    loadings_df <- head(loadings_df, top_n)
  }

  rownames(loadings_df) <- NULL

  return(loadings_df)
}


#' Calculate Variance Explained by Factors
#'
#' Tests which experimental factors (e.g., gonadal type, sex chromosome)
#' explain variance in PCA components.
#'
#' @param df Data frame. Features x Samples.
#' @param metadata Data frame with sample metadata including factors to test.
#' @param factors Character vector. Column names in metadata to test.
#' @param n_components Integer. Number of PCs to test (default: 10).
#' @param centering Logical. Center data? (default: TRUE)
#' @param scaling Logical. Scale data? (default: FALSE)
#'
#' @return Data frame with R² values for each factor x component combination.
#'
#' @examples
#' # Test which factors explain PC variance
#' var_explained <- calculate_factor_variance(
#'   scaled_RNA, metadata,
#'   factors = c("gonadal_type", "allosome_type", "group")
#' )
#'
#' @export
calculate_factor_variance <- function(df, metadata, factors, n_components = 10,
                                      centering = TRUE, scaling = FALSE) {

  # Perform PCA
  pca <- prcomp(t(df[complete.cases(df), ]), center = centering, scale. = scaling)

  # Match samples
  common_samples <- intersect(colnames(df), metadata$sample)
  metadata <- metadata[metadata$sample %in% common_samples, ]
  metadata <- metadata[match(common_samples, metadata$sample), ]

  # Extract PC scores
  pc_scores <- pca$x[, 1:n_components]

  # Test each factor
  results <- list()

  for (factor_name in factors) {
    if (!factor_name %in% names(metadata)) {
      warning("Factor '", factor_name, "' not found in metadata. Skipping.")
      next
    }

    factor_values <- metadata[[factor_name]]

    # Calculate R² for each component
    r_squared <- sapply(1:n_components, function(i) {
      model <- lm(pc_scores[, i] ~ factor_values)
      summary(model)$r.squared
    })

    results[[factor_name]] <- r_squared
  }

  # Combine into data frame
  results_df <- as.data.frame(do.call(cbind, results))
  results_df$Component <- paste0("PC", 1:n_components)
  results_df <- results_df[, c("Component", factors)]

  return(results_df)
}
