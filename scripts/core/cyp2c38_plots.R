# Scientific QC and result plots for Workflow 07b.
# Base graphics keep these figures independent of additional plotting packages.
# Concordance p-values below are supplied mouse-level permutation results;
# no independent-lipid regression p-value is calculated or displayed.

wp1_plot_pdf <- function(path, width, height, draw) {
  grDevices::pdf(path, width = width, height = height, onefile = TRUE,
                 family = "Helvetica", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
  invisible(path)
}

wp1_plot_range <- function(x, symmetric = FALSE, extra = 0.07) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("A plot axis contains no finite values.")
  if (symmetric) {
    radius <- max(abs(x))
    if (radius == 0) radius <- 1
    return(c(-radius, radius) * (1 + extra))
  }
  limits <- range(x)
  span <- diff(limits)
  if (span == 0) span <- max(abs(limits), 1)
  limits + c(-1, 1) * span * extra
}

wp1_plot_number <- function(x, digits = 3L) {
  if (length(x) != 1L || !is.finite(x)) return("not available")
  formatC(x, digits = digits, format = "f")
}

wp1_plot_probability <- function(x) {
  if (length(x) != 1L || !is.finite(x)) return("not available")
  format.pval(x, digits = 3L, eps = 0.0001)
}

wp1_make_plots <- function(out, y, meta, feature_ledger, fit_primary,
                           differential, consistency, concordance) {
  required_meta <- c("sample_id", "sex", "genotype", "group",
                     "precision_weight")
  required_differential <- c("feature_id", "contrast", "logFC", "t",
                             "adj.P.Val", "canonical_lipid", "lipid_class")
  required_consistency <- c("feature_id", "contrast", "reference_set",
                            "reference_r", "t", "logFC", "consistency",
                            "lipid_class", "canonical_lipid")
  required_concordance <- c("contrast", "reference_set", "pearson_r",
                            "concordance_score", "permutation_p_two_sided",
                            "permutation_fdr", "bootstrap_ci_low",
                            "bootstrap_ci_high")
  if (!is.matrix(y) || !is.numeric(y) || !all(is.finite(y)) ||
      ncol(y) != nrow(meta) || nrow(y) < 2L || ncol(y) < 4L) {
    stop("Plots require finite feature-by-mouse y and aligned metadata.")
  }
  if (!all(required_meta %in% names(meta)) ||
      !identical(colnames(y), as.character(meta$sample_id)) ||
      !all(required_differential %in% names(differential)) ||
      !all(required_consistency %in% names(consistency)) ||
      !all(required_concordance %in% names(concordance)) ||
      !all(c("qc_cv_percent", "keep_primary") %in% names(feature_ledger))) {
    stop("Plot input schemas or sample order are inconsistent.")
  }
  if (!inherits(fit_primary, "MArrayLM")) {
    stop("fit_primary must be the moderated limma MArrayLM fit.")
  }
  group_levels <- c("F_WT", "F_KO", "M_WT", "M_KO")
  groups <- factor(meta$group, levels = group_levels)
  sex <- toupper(as.character(meta$sex))
  if (anyNA(groups) || !all(sex %in% c("F", "M")) ||
      any(!is.finite(meta$precision_weight)) ||
      any(meta$precision_weight <= 0)) {
    stop("Plot metadata contains invalid groups, sex or precision weights.")
  }
  contrast_names <- c("pooled_KO_vs_WT", "female_KO_vs_WT",
                      "male_KO_vs_WT", "interaction_M_minus_F")
  contrast_titles <- c("Pooled KO - WT (equal sexes)", "Female KO - WT",
                       "Male KO - WT", "Interaction: male - female effect")
  if (!all(contrast_names %in% differential$contrast)) {
    stop("Differential plots require all four prespecified contrasts.")
  }
  plot_dir <- file.path(out, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  palette <- c(F_WT = "#e58c8c", F_KO = "#c73333",
               M_WT = "#88a5df", M_KO = "#315da8")
  point_colors <- unname(palette[as.character(groups)])
  point_symbols <- ifelse(toupper(meta$genotype) == "KO", 17L, 16L)
  animal_labels <- sub("_.*$", "", as.character(meta$sample_id))
  group_labels <- c("Female WT", "Female KO", "Male WT", "Male KO")
  n_group <- table(groups)
  group_legend <- paste0(group_labels, " (n = ", as.integer(n_group), ")")
  point_alpha <- grDevices::adjustcolor("#44596b", alpha.f = 0.42)

  qc_file <- file.path(plot_dir, "QC_PCA.pdf")
  wp1_plot_pdf(qc_file, 11, 7.5, function() {
    # PCA uses the supplied log2 features, centered and not scaled.
    pca <- stats::prcomp(t(y), center = TRUE, scale. = FALSE)
    explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    graphics::par(mar = c(5, 5, 4, 2), las = 1)
    graphics::plot(pca$x[, 1], pca$x[, 2], col = point_colors,
                   pch = point_symbols, cex = 1.25,
                   xlim = wp1_plot_range(pca$x[, 1], extra = 0.18),
                   ylim = wp1_plot_range(pca$x[, 2], extra = 0.2),
                   xlab = sprintf("PC1 (%.1f%% variance)", explained[1]),
                   ylab = sprintf("PC2 (%.1f%% variance)", explained[2]),
                   main = "Cyp2c38 lipidomics: biological sample PCA")
    graphics::abline(h = 0, v = 0, col = "#e1e6ea", lty = 3)
    graphics::points(pca$x[, 1], pca$x[, 2], col = point_colors,
                     pch = point_symbols, cex = 1.25)
    # Identify at most two distant animals per group for readable labels.
    # All animals remain plotted and are retained in the sample manifest.
    labelled <- unlist(lapply(group_levels, function(group_name) {
      selected <- which(groups == group_name)
      coordinates <- pca$x[selected, 1:2, drop = FALSE]
      distance <- rowSums(sweep(coordinates, 2L, colMeans(coordinates), "-")^2)
      selected[order(distance, decreasing = TRUE)[seq_len(min(2L, length(selected)))]]
    }), use.names = FALSE)
    graphics::text(pca$x[labelled, 1], pca$x[labelled, 2],
                   labels = animal_labels[labelled],
                   pos = ifelse(point_symbols[labelled] == 17L, 1L, 3L),
                   cex = 0.72, col = "#263d50", offset = 0.6)
    graphics::legend("topright", legend = group_legend,
                     col = palette, pch = c(16, 17, 16, 17),
                     bg = "white", box.col = "#d9e0e6", cex = 0.85)
    graphics::mtext(sprintf("%d primary features; log2, centered without scaling; two distant animals labelled per group",
                            nrow(y)), side = 3, line = 0.3, cex = 0.82)

    graphics::par(mfrow = c(1, 2), mar = c(5, 5, 4, 1.5), las = 1)
    qc_values <- feature_ledger$qc_cv_percent[
      feature_ledger$keep_primary %in% TRUE &
        is.finite(feature_ledger$qc_cv_percent)]
    if (length(qc_values)) {
      graphics::hist(qc_values, breaks = "FD", col = "#88a5b9",
                     border = "white", xlab = "CV across nine QC injections (%)",
                     ylab = "Primary features", main = "Analytical precision")
      graphics::abline(v = 30, col = "#a04f21", lwd = 2, lty = 2)
      graphics::legend("topright", legend = c(
        sprintf("Median CV = %.2f%%", stats::median(qc_values)),
        sprintf("%d exceed 30%%", sum(qc_values > 30)),
        "30% cutoff: sensitivity only"), bty = "n", cex = 0.8)
    } else {
      graphics::plot.new()
      graphics::title("Analytical precision")
      graphics::text(0.5, 0.5, "No finite QC CV values supplied.")
    }
    residuals <- y
    for (group_name in group_levels) {
      selected <- groups == group_name
      residuals[, selected] <- sweep(y[, selected, drop = FALSE], 1L,
                                     rowMeans(y[, selected, drop = FALSE]), "-")
    }
    residual_variances <- lapply(c("F", "M"), function(sex_name) {
      selected <- sex == sex_name
      rowSums(residuals[, selected, drop = FALSE]^2) / (sum(selected) - 2L)
    })
    names(residual_variances) <- c("Female", "Male")
    log_variances <- lapply(residual_variances, function(x) log2(x[x > 0]))
    graphics::boxplot(log_variances, outline = FALSE,
                      col = c("#e58c8c", "#88a5df"),
                      ylab = "log2 residual variance", main = "Within-group residual scale")
    graphics::mtext("Residuals of log2 intensity after four group means", side = 3,
                    line = 0.3, cex = 0.75)

    graphics::par(mfrow = c(1, 2), mar = c(8, 5, 4, 1.5), las = 1)
    graphics::barplot(meta$precision_weight, names.arg = animal_labels,
                      col = point_colors, border = NA, las = 2, cex.names = 0.62,
                      ylab = "Sample precision weight",
                      main = "Prespecified sex-group precision weights")
    graphics::abline(h = 1, col = "#6c7a86", lty = 2)
    graphics::mtext("These weights affect precision, not the equal-sex pooled contrast",
                    side = 3, line = 0.2, cex = 0.72)
    graphics::par(mar = c(5, 5, 4, 1.5))
    limma::plotSA(fit_primary, main = "Moderated mean-variance diagnostic")

    graphics::par(mfrow = c(1, 1), mar = c(8, 5, 4, 2), las = 1)
    graphics::boxplot(lapply(seq_len(ncol(residuals)), function(i) residuals[, i]),
                      names = animal_labels, outline = FALSE, col = point_colors,
                      border = "#526a7e", las = 2, cex.axis = 0.7,
                      ylab = "Residual log2 intensity",
                      main = "Biological sample residual distributions")
    graphics::abline(h = 0, col = "#526a7e", lty = 2)
    graphics::mtext("Per-feature residuals after sex/genotype group means; outlier points omitted for readability",
                    side = 3, line = 0.3, cex = 0.78)
  })

  differential_file <- file.path(plot_dir, "differential.pdf")
  wp1_plot_pdf(differential_file, 11, 8.5, function() {
    selected_tables <- lapply(contrast_names, function(name) {
      differential[differential$contrast == name, , drop = FALSE]
    })
    x_limits <- wp1_plot_range(unlist(lapply(selected_tables, function(x) x$logFC)),
                               symmetric = TRUE)
    y_values <- lapply(selected_tables, function(x) {
      -log10(pmax(x$adj.P.Val, .Machine$double.xmin))
    })
    y_limits <- c(0, max(c(unlist(y_values), -log10(0.05)), na.rm = TRUE) * 1.08 + 0.2)
    graphics::par(mfrow = c(2, 2), mar = c(4.6, 4.7, 3.5, 1), las = 1)
    for (i in seq_along(selected_tables)) {
      tab <- selected_tables[[i]]
      significant <- is.finite(tab$adj.P.Val) & tab$adj.P.Val <= 0.05
      colors <- rep(grDevices::adjustcolor("#6c7a86", alpha.f = 0.4), nrow(tab))
      colors[significant & tab$logFC > 0] <- "#593d8a"
      colors[significant & tab$logFC < 0] <- "#087f82"
      graphics::plot(tab$logFC, y_values[[i]], pch = 16, cex = 0.65,
                     col = colors, xlim = x_limits, ylim = y_limits,
                     xlab = if (i == 4L) "Male - female difference in log2 effects" else "KO - WT log2 fold change",
                     ylab = "-log10 BH-adjusted P value",
                     main = contrast_titles[i])
      graphics::abline(h = -log10(0.05), col = "#a04f21", lty = 2)
      graphics::abline(v = 0, col = "#d9e0e6", lty = 3)
      graphics::mtext(sprintf("%d / %d features at FDR <= 0.05",
                              sum(significant), nrow(tab)), side = 3,
                      line = 0.25, cex = 0.78)
      if (i == 1L) {
        graphics::legend("topright", legend = c("Increase: FDR <= 0.05",
                         "Decrease: FDR <= 0.05", "FDR > 0.05"),
                         col = c("#593d8a", "#087f82", "#6c7a86"),
                         pch = 16, bty = "n", cex = 0.7)
      }
    }
    female <- selected_tables[[2]][, c("feature_id", "logFC")]
    male <- selected_tables[[3]][, c("feature_id", "logFC")]
    effects <- merge(female, male, by = "feature_id", suffixes = c("_F", "_M"))
    limits <- wp1_plot_range(c(effects$logFC_F, effects$logFC_M), symmetric = TRUE)
    graphics::par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), las = 1)
    graphics::plot(effects$logFC_F, effects$logFC_M, pch = 16, cex = 0.75,
                   col = point_alpha, xlim = limits, ylim = limits, asp = 1,
                   xlab = "Female KO - WT log2 fold change",
                   ylab = "Male KO - WT log2 fold change",
                   main = "Cyp2c38 perturbation effects across sexes")
    graphics::abline(h = 0, v = 0, col = "#d9e0e6", lty = 3)
    graphics::abline(a = 0, b = 1, col = "#526a7e", lty = 2)
    graphics::mtext(sprintf("%d features; dashed line denotes equal effects", nrow(effects)),
                    side = 3, line = 0.3, cex = 0.82)
  })

  concordance_file <- file.path(plot_dir, "concordance.pdf")
  main_contrasts <- contrast_names[1:3]
  main_references <- c("pooled", "female", "male")
  main_titles <- c("Pooled prediction and pooled KO effect",
                   "Female prediction and female KO effect",
                   "Male prediction and male KO effect")
  main_edges <- lapply(seq_along(main_contrasts), function(i) {
    consistency[consistency$contrast == main_contrasts[i] &
                  consistency$reference_set == main_references[i], , drop = FALSE]
  })
  if (any(vapply(main_edges, nrow, integer(1)) < 3L)) {
    stop("Each main concordance panel needs at least three matched lipid features.")
  }
  t_limits <- wp1_plot_range(unlist(lapply(main_edges, function(x) x$t)))
  draw_concordance_panel <- function(edges, contrast, reference, title, calibrated) {
    finite <- is.finite(edges$reference_r) & is.finite(edges$t)
    edges <- edges[finite, , drop = FALSE]
    stats_row <- concordance[concordance$contrast == contrast &
                               concordance$reference_set == reference, , drop = FALSE]
    if (nrow(stats_row) != 1L || nrow(edges) < 3L) {
      stop("Concordance panel requires one summary row and finite matched edges.")
    }
    graphics::plot(edges$reference_r, edges$t, pch = 16, cex = 0.65,
                   col = point_alpha, xlim = c(-1, 1), ylim = t_limits,
                   xlab = "Reference Cyp2c38-lipid correlation (r)",
                   ylab = "KO - WT lipid effect (moderated t)", main = title,
                   cex.main = 0.94, cex.lab = 0.92)
    graphics::abline(h = 0, v = 0, col = "#d9e0e6", lty = 3)
    if (stats::sd(edges$reference_r) > 0) {
      guide <- stats::lm(t ~ reference_r, data = edges)
      graphics::abline(guide, col = "#a04f21", lwd = 1.8)
    }
    annotation <- c(
      sprintf("Matched lipids: %d", nrow(edges)),
      paste0("Pearson r: ", wp1_plot_number(stats_row$pearson_r)),
      paste0("Signed concordance (-r): ",
             wp1_plot_number(stats_row$concordance_score)),
      paste0("95% combined-bootstrap score CI: [",
             wp1_plot_number(stats_row$bootstrap_ci_low), ", ",
             wp1_plot_number(stats_row$bootstrap_ci_high), "]")
    )
    if (calibrated) {
      annotation <- c(annotation,
        paste0("Genotype-permutation P (two-sided): ",
               wp1_plot_probability(stats_row$permutation_p_two_sided)),
        paste0("BH FDR across three main tests: ",
               wp1_plot_probability(stats_row$permutation_fdr)))
    }
    graphics::mtext(annotation, side = 1,
                    line = seq(4.1, by = 1.05, length.out = length(annotation)),
                    adj = 0, cex = 0.76)
  }
  wp1_plot_pdf(concordance_file, 15, 7.5, function() {
    graphics::par(mfrow = c(1, 3), cex = 1, mar = c(11.2, 4.6, 4.2, 1.1),
                  oma = c(0, 0, 2.7, 0), las = 1)
    for (i in seq_along(main_contrasts)) {
      draw_concordance_panel(main_edges[[i]], main_contrasts[i],
                             main_references[i], main_titles[i], TRUE)
    }
    graphics::mtext("Cyp2c38 in vivo validation: prespecified main comparisons",
                    side = 3, outer = TRUE, line = 1.3, font = 2, cex = 1.1)
    graphics::mtext("Negative r and positive signed concordance support the predicted perturbation direction; the fitted line is descriptive.",
                    side = 3, outer = TRUE, line = 0, cex = 0.83)
    secondary <- lapply(2:3, function(i) {
      consistency[consistency$contrast == main_contrasts[i] &
                    consistency$reference_set == "pooled", , drop = FALSE]
    })
    if (all(vapply(secondary, nrow, integer(1)) >= 3L)) {
      t_limits <- wp1_plot_range(unlist(lapply(secondary, function(x) x$t)))
      graphics::par(mfrow = c(1, 2), mar = c(9, 5, 4, 2),
                    oma = c(0, 0, 2.7, 0), las = 1)
      for (i in seq_along(secondary)) {
        draw_concordance_panel(secondary[[i]], main_contrasts[i + 1L], "pooled",
                               c("Pooled prediction and female KO effect",
                                 "Pooled prediction and male KO effect")[i], FALSE)
      }
      graphics::mtext("Secondary descriptions: sex-specific effects against the pooled reference",
                      side = 3, outer = TRUE, line = 1.3, font = 2, cex = 1.1)
      graphics::mtext("These panels are outside the three-test main permutation family.",
                      side = 3, outer = TRUE, line = 0, cex = 0.83)
    }
  })
  cat("  Plot files:", qc_file, differential_file, concordance_file, sep = "\n")
  invisible(c(QC_PCA = qc_file, differential = differential_file,
              concordance = concordance_file))
}
