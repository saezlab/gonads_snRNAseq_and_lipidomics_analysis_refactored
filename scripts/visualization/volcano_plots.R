#' Volcano Plot Visualization Module
#'
#' This module provides functions for creating publication-quality volcano plots
#' with asymptotic significance curves.
#'
#' These functions reproduce the exact plotting behavior from support_functions.R
#' and are used to generate manuscript figures.
#'
#' Key features:
#' - Asymptotic curves for significance thresholds (non-straight)
#' - Automatic labeling of top hits
#' - Manual label override capability
#' - Customizable colors for up/down/neutral
#' - Exact reproduction of legacy plotting behavior

library(ggplot2)
library(ggrepel)


#' Create Volcano Plot with Asymptotic Curves
#'
#' Creates a publication-quality volcano plot with asymptotic significance curves.
#' This function replicates the exact behavior of volcano_nice() from support_functions.R.
#'
#' The asymptotic curves use the formula:
#'   For negative logFC: |hAss - 1 + x/(x + vAss)|
#'   For positive logFC: |hAss - 1 + x/(x - vAss)|
#'
#' This creates curves that approach the significance threshold asymptotically,
#' providing more stringent significance criteria for larger fold changes.
#'
#' @param df Data frame with differential analysis results. Should contain:
#'   - logFC (or equivalent fold change column)
#'   - adj.P.Val or P.Value (adjusted p-values)
#'   - ID column (gene/feature names)
#' @param hAss Numeric. Horizontal FDR threshold (default: 0.05).
#' @param FCIndex Integer. Column index for log fold change.
#' @param pValIndex Integer. Column index for adjusted p-value.
#' @param IDIndex Integer. Column index for feature IDs/names.
#' @param vAss Numeric. Vertical fold change threshold. If NULL, set to xlimAbs/10.
#' @param label Logical. Add labels to significant points? (default: FALSE)
#' @param straight Logical. Use straight threshold lines instead of curves? (default: FALSE)
#' @param nlabels Integer. Number of top hits to label (default: 20).
#' @param manual_labels Character vector. Specific labels to include regardless of ranking.
#' @param col_high_neut_low Character vector of 3 colors for c(high, neutral, low).
#'   Default: c("#e54141", "grey30", "#4173d8")
#'
#' @return ggplot object
#'
#' @examples
#' # Basic volcano plot
#' volcano_nice(ttop, FCIndex = 2, pValIndex = 5, IDIndex = 1)
#'
#' # With labels
#' volcano_nice(ttop, FCIndex = 2, pValIndex = 5, IDIndex = 1,
#'              label = TRUE, nlabels = 40)
#'
#' # With manual labels
#' volcano_nice(ttop, FCIndex = 2, pValIndex = 5, IDIndex = 1,
#'              label = TRUE, nlabels = 20,
#'              manual_labels = c("Cyp2b9", "Cyp2c38"))
#'
#' # Straight thresholds (classical volcano plot)
#' volcano_nice(ttop, FCIndex = 2, pValIndex = 5, IDIndex = 1,
#'              straight = TRUE, vAss = 0.5, hAss = 0.05)
#'
#' @export
volcano_nice <- function(df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                         label = FALSE, straight = FALSE, nlabels = 20,
                         manual_labels = NA,
                         col_high_neut_low = c("#e54141", "grey30", "#4173d8")) {

  # Remove incomplete cases
  df <- df[complete.cases(df), ]

  # Ensure first column is named consistently
  names(df)[1] <- "X1"

  # Store original hAss for comparison
  hAssOri <- hAss
  hAss <- -log10(hAss)

  # Standardize column names
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"

  # Calculate plot limits
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log10(df[, pValIndex]))))
  } else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log10(df[, pValIndex])))
  }

  # Set vertical threshold if not provided
  if (is.null(vAss)) {
    vAss <- xlimAbs / 10
  }

  # Define asymptotic curve functions
  xneg <- function(x) abs(hAss - 1 + x / (x + vAss))
  xpos <- function(x) abs(hAss - 1 + x / (x - vAss))

  # Test function to determine significance
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log10(y)) {
        return("1")
      } else {
        return("0")
      }
    } else {
      if (x > vAss) {
        if (xpos(x) < -log10(y)) {
          return("1")
        } else {
          return("0")
        }
      } else {
        return("0")
      }
    }
  }

  # Assign colors based on significance
  if (straight) {
    # Straight threshold: simple FC and p-value cutoffs
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <= hAssOri,
                         "1", "0")
  } else {
    # Asymptotic curves: more stringent
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) {
      test(as.numeric(x[FCIndex]), as.numeric(x[pValIndex]), vAss)
    })
  }

  # Sort by p-value for label selection
  df <- df[order(df$adj.P.Val, decreasing = FALSE), ]

  # Prepare labels
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA

  # Keep only top n labels + manual labels
  labels_to_keep <- c(df[1:nlabels, "condLabel"], manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA

  # Differentiate up-regulated (high) from down-regulated (low)
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2", df$couleur)

  # Print color distribution for diagnostics
  print(table(df$couleur))

  # Create plot
  if (label) {
    # With labels
    a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) +
      geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) +
      ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) +
      stat_function(fun = xpos, xlim = c(vAss, xlimAbs),
                    color = "black", alpha = 0.7) +
      scale_colour_manual(values = c(col_high_neut_low[2],
                                     col_high_neut_low[3],
                                     col_high_neut_low[1])) +
      theme_minimal() +
      theme(legend.position = "none")
  } else {
    # Without labels
    if (straight) {
      # Straight thresholds
      a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
        geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = col_high_neut_low[3]) +
        geom_vline(xintercept = vAss, color = col_high_neut_low[3]) +
        ylim(c(0, ylimAbs)) +
        xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = col_high_neut_low[1]) +
        scale_colour_manual(values = c(col_high_neut_low[2],
                                       col_high_neut_low[1],
                                       col_high_neut_low[3])) +
        theme_minimal() +
        theme(legend.position = "none")
    } else {
      # Asymptotic curves
      a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
        geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                      color = "black", alpha = 0.7) +
        ylim(c(0, ylimAbs)) +
        xlim(c(-xlimAbs, xlimAbs)) +
        stat_function(fun = xpos, xlim = c(vAss, xlimAbs),
                      color = "black", alpha = 0.7) +
        scale_colour_manual(values = c(col_high_neut_low[2],
                                       col_high_neut_low[1],
                                       col_high_neut_low[3])) +
        theme_minimal() +
        theme(legend.position = "none")
    }
  }

  return(a)
}


#' Create Multiple Volcano Plots
#'
#' Creates volcano plots for multiple comparisons/contrasts.
#'
#' @param ttop_list Named list of toptables (output from limma).
#' @param FC_column Character or integer. Column name or index for log fold change (default: "logFC").
#' @param pval_column Character or integer. Column name or index for adjusted p-value (default: "adj.P.Val").
#' @param id_column Character or integer. Column name or index for feature IDs (default: "ID").
#' @param label Logical. Add labels? (default: FALSE)
#' @param nlabels Integer. Number of top hits to label (default: 20).
#' @param save_plots Logical. Save plots to files? (default: TRUE)
#' @param output_dir Character. Output directory (default: "results/volcanos/").
#' @param width Numeric. Plot width in inches (default: 9).
#' @param height Numeric. Plot height in inches (default: 9).
#' @param verbose Logical. Print progress messages?
#'
#' @return List of ggplot objects (named by comparison)
#'
#' @examples
#' # Create volcano plots for all contrasts
#' plots <- create_volcano_plots(ttop_list, label = TRUE, nlabels = 40)
#'
#' # Save plots
#' create_volcano_plots(ttop_list, save_plots = TRUE,
#'                      output_dir = "results/pseudobulk/volcanos/")
#'
#' @export
create_volcano_plots <- function(ttop_list,
                                  FC_column = "logFC",
                                  pval_column = "adj.P.Val",
                                  id_column = "ID",
                                  label = FALSE,
                                  nlabels = 20,
                                  save_plots = TRUE,
                                  output_dir = "results/volcanos/",
                                  width = 9,
                                  height = 9,
                                  verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Creating Volcano Plots\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Create output directory
  if (save_plots) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  plot_list <- list()

  for (comparison_name in names(ttop_list)) {
    if (verbose) {
      cat("  Processing:", comparison_name, "\n")
    }

    ttop <- ttop_list[[comparison_name]]

    # Determine column indices
    if (is.character(FC_column)) {
      FCIndex <- which(names(ttop) == FC_column)
    } else {
      FCIndex <- FC_column
    }

    if (is.character(pval_column)) {
      pValIndex <- which(names(ttop) == pval_column)
    } else {
      pValIndex <- pval_column
    }

    if (is.character(id_column)) {
      IDIndex <- which(names(ttop) == id_column)
    } else {
      IDIndex <- id_column
    }

    # Create plot
    p <- volcano_nice(
      df = ttop,
      FCIndex = FCIndex,
      pValIndex = pValIndex,
      IDIndex = IDIndex,
      label = label,
      nlabels = nlabels
    )

    plot_list[[comparison_name]] <- p

    # Save plot
    if (save_plots) {
      # Clean comparison name for filename
      filename <- gsub(" ", "_", comparison_name)
      filename <- gsub("[^A-Za-z0-9_-]", "", filename)
      filepath <- file.path(output_dir, paste0("volcano_", filename, ".pdf"))

      ggsave(filepath, plot = p, width = width, height = height)

      if (verbose) {
        cat("    Saved:", filepath, "\n")
      }
    }
  }

  if (verbose) {
    cat("\n")
    cat("Volcano plots complete!\n")
    cat("  Total plots created:", length(plot_list), "\n")
    if (save_plots) {
      cat("  Saved to:", output_dir, "\n")
    }
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(plot_list)
}


#' Extract Significant Features from Volcano Plot Criteria
#'
#' Extracts features that meet the asymptotic curve significance criteria.
#'
#' @param df Data frame with differential analysis results.
#' @param hAss Numeric. FDR threshold (default: 0.05).
#' @param vAss Numeric. Fold change threshold. If NULL, calculated automatically.
#' @param straight Logical. Use straight thresholds? (default: FALSE)
#'
#' @return Data frame with significant features
#'
#' @examples
#' sig_features <- extract_volcano_significant(ttop, hAss = 0.05, vAss = 0.5)
#'
#' @export
extract_volcano_significant <- function(df, hAss = 0.05, vAss = NULL, straight = FALSE) {

  df <- df[complete.cases(df), ]

  # Ensure logFC and adj.P.Val columns exist
  if (!"logFC" %in% names(df)) {
    stop("Column 'logFC' not found in data frame")
  }
  if (!"adj.P.Val" %in% names(df)) {
    stop("Column 'adj.P.Val' not found in data frame")
  }

  hAssOri <- hAss
  hAss <- -log10(hAss)

  # Calculate vAss if not provided
  if (is.null(vAss)) {
    if (max(abs(df$logFC)) >= 1) {
      xlimAbs <- ceiling(max(abs(df$logFC)))
    } else {
      xlimAbs <- max(abs(df$logFC))
    }
    vAss <- xlimAbs / 10
  }

  if (straight) {
    # Straight threshold
    sig_df <- df[abs(df$logFC) >= vAss & df$adj.P.Val <= hAssOri, ]
  } else {
    # Asymptotic curves
    xneg <- function(x) abs(hAss - 1 + x / (x + vAss))
    xpos <- function(x) abs(hAss - 1 + x / (x - vAss))

    test <- function(x, y, vAss) {
      if (x < -vAss) {
        if (xneg(x) < -log10(y)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      } else {
        if (x > vAss) {
          if (xpos(x) < -log10(y)) {
            return(TRUE)
          } else {
            return(FALSE)
          }
        } else {
          return(FALSE)
        }
      }
    }

    df$is_sig <- apply(df, 1, FUN = function(x) {
      test(as.numeric(x["logFC"]), as.numeric(x["adj.P.Val"]), vAss)
    })

    sig_df <- df[df$is_sig, ]
    sig_df$is_sig <- NULL
  }

  return(sig_df)
}


#' Customize Volcano Plot Appearance
#'
#' Applies custom theme and labels to an existing volcano plot.
#'
#' @param plot ggplot object from volcano_nice().
#' @param title Character. Plot title.
#' @param xlab Character. X-axis label (default: "log2 Fold Change").
#' @param ylab Character. Y-axis label (default: "-log10(FDR)").
#' @param base_size Numeric. Base font size (default: 12).
#'
#' @return Modified ggplot object
#'
#' @examples
#' p <- volcano_nice(ttop, FCIndex = 2, pValIndex = 5, IDIndex = 1)
#' p <- customize_volcano(p, title = "Testes vs Ovaries",
#'                        xlab = "log2 FC", ylab = "-log10(adj P)")
#'
#' @export
customize_volcano <- function(plot, title = NULL, xlab = "log2 Fold Change",
                               ylab = "-log10(FDR)", base_size = 12) {

  plot <- plot +
    labs(x = xlab, y = ylab, title = title) +
    theme(
      text = element_text(size = base_size),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5)
    )

  return(plot)
}


#' Create Xenium Validation Volcano Plot
#'
#' Creates a volcano plot specifically for Xenium spatial transcriptomics validation.
#' This variant of volcano_nice() uses manual label selection only, making it ideal
#' for highlighting specific validated genes from the Xenium 100-gene panel.
#'
#' Key differences from volcano_nice():
#' - No automatic top N label selection - only manual_labels are shown
#' - Designed for pre-filtered gene sets (e.g., Xenium panel genes)
#' - Optimized for rotated plots (use with coord_flip())
#'
#' @param df Data frame with differential analysis results. Should contain:
#'   - logFC (or equivalent fold change column)
#'   - adj.P.Val or P.Value (adjusted p-values)
#'   - ID column (gene/feature names)
#' @param hAss Numeric. Horizontal FDR threshold (default: 0.05).
#' @param FCIndex Integer. Column index for log fold change.
#' @param pValIndex Integer. Column index for adjusted p-value.
#' @param IDIndex Integer. Column index for feature IDs/names.
#' @param vAss Numeric. Vertical fold change threshold. If NULL, set to xlimAbs/10.
#' @param label Logical. Add labels to points? (default: FALSE)
#' @param straight Logical. Use straight threshold lines instead of curves? (default: FALSE)
#' @param nlabels Integer. Unused in this function (kept for compatibility).
#' @param manual_labels Character vector. Specific labels to include. Only these will be shown.
#' @param col_high_neut_low Character vector of 3 colors for c(high, neutral, low).
#'   Default: c("#e54141", "grey30", "#4173d8")
#'
#' @return ggplot object
#'
#' @examples
#' # Load Xenium panel and filter to validation genes
#' xenium_panel <- read.csv("support/FCGL_Xenium_Panel.csv")
#' xenium_genes <- xenium_panel$Gene[xenium_panel$Reason != "cell type marker"]
#'
#' # Filter DE results to Xenium genes
#' ttop_xenium <- ttop[ttop$ID %in% xenium_genes, ]
#'
#' # Define key genes to label
#' key_genes <- c("Cyp2b9", "C8a", "Slco1a1", "Nipal1")
#'
#' # Create rotated volcano plot
#' p <- volcano_nice_xenium(ttop_xenium,
#'                          FCIndex = 2, pValIndex = 5, IDIndex = 1,
#'                          label = TRUE, manual_labels = key_genes) +
#'      coord_flip()
#'
#' @export
volcano_nice_xenium <- function(df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                                label = FALSE, straight = FALSE, nlabels = 100,
                                manual_labels = NA,
                                col_high_neut_low = c("#e54141", "grey30", "#4173d8")) {

  # Remove incomplete cases
  df <- df[complete.cases(df), ]

  # Ensure first column is named consistently
  names(df)[1] <- "X1"

  # Store original hAss for comparison
  hAssOri <- hAss
  hAss <- -log10(hAss)

  # Standardize column names
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"

  # Calculate plot limits
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log10(df[, pValIndex]))))
  } else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log10(df[, pValIndex])))
  }

  # Set vertical threshold if not provided
  if (is.null(vAss)) {
    vAss <- xlimAbs / 10
  }

  # Define asymptotic curve functions
  xneg <- function(x) abs(hAss - 1 + x / (x + vAss))
  xpos <- function(x) abs(hAss - 1 + x / (x - vAss))

  # Test function to determine significance
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log10(y)) {
        return("1")
      } else {
        return("0")
      }
    } else {
      if (x > vAss) {
        if (xpos(x) < -log10(y)) {
          return("1")
        } else {
          return("0")
        }
      } else {
        return("0")
      }
    }
  }

  # Assign colors based on significance
  if (straight) {
    # Straight threshold: simple FC and p-value cutoffs
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <= hAssOri,
                         "1", "0")
  } else {
    # Asymptotic curves: more stringent
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) {
      test(as.numeric(x[FCIndex]), as.numeric(x[pValIndex]), vAss)
    })
  }

  # Sort by p-value
  df <- df[order(df$adj.P.Val, decreasing = FALSE), ]

  # Prepare labels - ONLY use manual labels (no automatic selection)
  df$condLabel <- df[, IDIndex]
  labels_to_keep <- c(manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA

  # Differentiate up-regulated (high) from down-regulated (low)
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2", df$couleur)

  # Create plot
  if (label) {
    # With labels
    a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) +
      geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) +
      ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) +
      stat_function(fun = xpos, xlim = c(vAss, xlimAbs),
                    color = "black", alpha = 0.7) +
      scale_colour_manual(values = c(col_high_neut_low[2],
                                     col_high_neut_low[3],
                                     col_high_neut_low[1])) +
      theme_minimal() +
      theme(legend.position = "none")
  } else {
    # Without labels
    if (straight) {
      # Straight thresholds
      a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
        geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = col_high_neut_low[3]) +
        geom_vline(xintercept = vAss, color = col_high_neut_low[3]) +
        ylim(c(0, ylimAbs)) +
        xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = col_high_neut_low[1]) +
        scale_colour_manual(values = c(col_high_neut_low[2],
                                       col_high_neut_low[1],
                                       col_high_neut_low[3])) +
        theme_minimal() +
        theme(legend.position = "none")
    } else {
      # Asymptotic curves
      a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
        geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                      color = "black", alpha = 0.7) +
        ylim(c(0, ylimAbs)) +
        xlim(c(-xlimAbs, xlimAbs)) +
        stat_function(fun = xpos, xlim = c(vAss, xlimAbs),
                      color = "black", alpha = 0.7) +
        scale_colour_manual(values = c(col_high_neut_low[2],
                                       col_high_neut_low[1],
                                       col_high_neut_low[3])) +
        theme_minimal() +
        theme(legend.position = "none")
    }
  }

  return(a)
}
