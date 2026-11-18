# Four independent panels for the direct CF453 diet-interaction analysis.
# Sourcing defines functions only. Models, tests and example selection are inputs.

cfdi_theme <- function(size = 11) {
  ggplot2::theme_classic(base_size = size, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "#4C4C4C", size = size - 1,
                                           margin = ggplot2::margin(b = 9)),
      plot.caption = ggplot2::element_text(color = "#555555", size = size - 2,
                                          hjust = 0, lineheight = 1.2),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.text = ggplot2::element_text(color = "#303030"),
      axis.line = ggplot2::element_line(color = "#333333", linewidth = .35),
      axis.ticks = ggplot2::element_line(linewidth = .3),
      strip.background = ggplot2::element_rect(fill = "#F3F3F3", color = NA),
      strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(7, 5, 7, 5)),
      panel.spacing = grid::unit(1.2, "lines"), legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = size - 1),
      legend.text = ggplot2::element_text(size = size - 1),
      plot.margin = ggplot2::margin(10, 13, 10, 10)
    )
}

cfdi_q <- function(x) {
  vapply(x, function(z) if (!is.finite(z)) "NA" else if (z < .001)
    formatC(z, digits = 1, format = "e") else formatC(z, digits = 2, format = "g"), character(1L))
}

cfdi_bind <- function(...) {
  xs <- list(...)
  columns <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, xs)
}

cfdi_empty <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = .5, y = .5, label = message, size = 4,
                      color = "#555555", family = "Helvetica", lineheight = 1.3) +
    ggplot2::labs(title = title) + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    cfdi_theme() + ggplot2::theme(axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank())
}

cfdi_plot_overview <- function(data) {
  labels <- c(diet_sex = "Diet x sex\n(mean over zones)",
              diet_zone = "Diet x zone\n(mean over sexes)",
              diet_sex_zone = "Diet x sex x zone")
  if (any(!data$contrast %in% names(labels))) stop("Unknown primary interaction in overview.")
  d <- data
  d$contrast <- factor(d$contrast, levels = names(labels))
  d$neglog10q <- -log10(pmax(d$q, .Machine$double.xmin))
  d$significant <- d$q <= .05
  d$label <- FALSE
  for (contrast in names(labels)) {
    selected <- which(d$contrast == contrast & d$significant & is.finite(d$q))
    selected <- selected[order(d$q[selected], d$p[selected], d$gene[selected])]
    d$label[head(selected, 2)] <- TRUE
  }
  counts <- stats::aggregate(significant ~ contrast, data = d, FUN = sum)
  counts$label <- paste0(format(counts$significant, big.mark = ",", trim = TRUE), " genes at FDR <= 0.05")
  p <- ggplot2::ggplot(d, ggplot2::aes(effect, neglog10q)) +
    ggplot2::geom_hline(yintercept = -log10(.05), color = "#AAAAAA", linetype = "22", linewidth = .35) +
    ggplot2::geom_vline(xintercept = 0, color = "#BBBBBB", linewidth = .35) +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = .65, alpha = .65) +
    ggplot2::geom_text(data = counts, ggplot2::aes(x = -Inf, y = Inf, label = label),
                       inherit.aes = FALSE, hjust = -.04, vjust = 1.3,
                       size = 3, color = "#444444", family = "Helvetica") +
    ggplot2::facet_wrap(~contrast, nrow = 1, scales = "free_x", labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#C8C8C8", `TRUE` = "#76517B"),
      labels = c(`FALSE` = "FDR > 0.05", `TRUE` = "FDR <= 0.05"), name = "Direct interaction") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.03, .18))) +
    ggplot2::labs(title = "Direct diet effects on sex and zonation differences",
      subtitle = "WD-minus-SD change in each expression contrast",
      x = "Interaction effect (log-normalized RNA)", y = "-log10 BH-adjusted P value") +
    cfdi_theme(10)
  if (any(d$label)) p <- p + ggrepel::geom_label_repel(
    data = d[d$label, , drop = FALSE], ggplot2::aes(label = gene),
    seed = 1506, size = 3, color = "#53365A", fill = "white", label.size = .1,
    label.padding = grid::unit(.1, "lines"), point.padding = .4, box.padding = .4,
    min.segment.length = 0, segment.size = .25, max.overlaps = Inf, family = "Helvetica")
  list(plot = p, source = d)
}

cfdi_plot_heatmap <- function(data) {
  title <- "How diet changes sex and zonation effects"
  if (!nrow(data)) return(list(plot = cfdi_empty(title, "No genes pass the joint diet-interaction FDR threshold"), source = data))
  columns <- c(diet_sex_PC = "Sex gap\nPC", diet_sex_PP = "Sex gap\nPP",
               diet_zone_female = "Zone gap\nFemale", diet_zone_male = "Zone gap\nMale")
  if (any(!data$contrast %in% names(columns))) stop("Unknown conditional interaction in heatmap.")
  d <- data
  gene_order <- unique(d$gene[order(d$rank, d$gene)])
  if (length(gene_order) > 20L) stop("Heatmap is restricted to at most twenty input-selected genes.")
  d$gene <- factor(d$gene, levels = rev(gene_order))
  d$x <- match(d$contrast, names(columns))
  max_effect <- max(abs(d$effect), na.rm = TRUE)
  max_effect <- max(max_effect, .01)
  d$marker <- ifelse(d$q <= .05, "*", "")
  d$text_color <- ifelse(abs(d$effect) >= .55 * max_effect, "white", "#303030")
  qtab <- unique(d[, c("gene", "joint_q", "threeway_q")])
  if (anyDuplicated(qtab$gene)) stop("Gene-level heatmap evidence is inconsistent across columns.")
  qtab <- rbind(data.frame(gene = qtab$gene, x = 5.35, text = cfdi_q(qtab$joint_q)),
                data.frame(gene = qtab$gene, x = 6.5, text = cfdi_q(qtab$threeway_q)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, gene)) +
    ggplot2::geom_tile(ggplot2::aes(fill = effect), width = .95, height = .95, color = "white", linewidth = .3) +
    ggplot2::geom_text(ggplot2::aes(label = marker, color = text_color), size = 4.3, vjust = .7, family = "Helvetica") +
    ggplot2::geom_vline(xintercept = 4.65, color = "#DDDDDD", linewidth = .4) +
    ggplot2::geom_text(data = qtab, ggplot2::aes(x = x, y = gene, label = text),
                       inherit.aes = FALSE, size = 2.9, color = "#333333", family = "Helvetica") +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(low = "#C64E54", mid = "#FAF9F7", high = "#467CBC", midpoint = 0,
      limits = c(-max_effect, max_effect), name = "Change in effect\n(log-normalized RNA)",
      guide = ggplot2::guide_colourbar(display = "rectangles", barheight = grid::unit(32, "mm"))) +
    ggplot2::scale_x_continuous(breaks = c(seq_len(4), 5.35, 6.5),
      labels = c(unname(columns), "Joint\nq", "Three-way\nq"), limits = c(.45, 7.05), expand = c(0, 0)) +
    ggplot2::labs(title = title,
      subtitle = "Sex gap: male - female; zone gap: PP - PC. Values show WD-minus-SD change.",
      x = NULL, y = NULL,
      caption = "* Conditional-interaction FDR <= 0.05. Rows are ranked by the joint 3-df diet-interaction test.") +
    cfdi_theme(10) + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "italic"), axis.text.x = ggplot2::element_text(size = 9),
      legend.position = "right", panel.grid = ggplot2::element_blank())
  list(plot = p, source = d)
}

cfdi_plot_examples <- function(points, means) {
  title <- "Diet-dependent expression patterns in individual mice"
  if (!nrow(points)) return(list(plot = cfdi_empty(title, "No supported examples at the primary FDR threshold"), source = points))
  p <- points
  m <- means
  contrast_labels <- c(diet_sex = "Diet x sex", diet_zone = "Diet x zone", diet_sex_zone = "Three-way")
  definitions <- unique(p[, c("gene", "metric", "example_rank", "selected_contrast", "q_primary")])
  definitions <- definitions[order(definitions$example_rank, definitions$gene), , drop = FALSE]
  if (nrow(definitions) > 4L || anyDuplicated(definitions$gene)) stop("Examples require at most four distinct input-selected genes.")
  definitions$facet <- paste0(definitions$gene, " | ", definitions$metric, "\n",
    unname(contrast_labels[definitions$selected_contrast]), " q = ", cfdi_q(definitions$q_primary))
  p$facet <- factor(definitions$facet[match(p$gene, definitions$gene)], levels = definitions$facet)
  m$facet <- factor(definitions$facet[match(m$gene, definitions$gene)], levels = definitions$facet)
  p$sex <- factor(tolower(p$sex), levels = c("female", "male"))
  m$sex <- factor(tolower(m$sex), levels = c("female", "male"))
  p$diet <- factor(toupper(p$diet), levels = c("SD", "WD"))
  m$diet <- factor(toupper(m$diet), levels = c("SD", "WD"))
  dodge <- ggplot2::position_dodge(width = .4)
  plot <- ggplot2::ggplot(p, ggplot2::aes(diet, value, color = sex, group = sex)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .08, dodge.width = .4, seed = 1506),
                        size = 1.6, alpha = .65) +
    ggplot2::geom_line(data = m, ggplot2::aes(y = estimate), position = dodge, linewidth = .65) +
    ggplot2::geom_linerange(data = m, ggplot2::aes(y = estimate, ymin = low, ymax = high),
                            position = dodge, linewidth = .6) +
    ggplot2::geom_point(data = m, ggplot2::aes(y = estimate), position = dodge,
                        shape = 21, fill = "white", size = 2.7, stroke = .9) +
    ggplot2::facet_wrap(~facet, ncol = 2, scales = "free_y", drop = FALSE) +
    ggplot2::scale_color_manual(values = c(female = "#CC4549", male = "#3C74BC"),
      labels = c(female = "Female", male = "Male"), name = NULL, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.10, .12))) +
    ggplot2::labs(title = title,
      subtitle = "Points: mice; lines and bars: fitted group means with nominal 95% CIs",
      x = NULL, y = "Within-mouse expression summary (log-normalized RNA)",
      caption = "Different mice received SD and WD; lines join group means, not paired animals.") +
    cfdi_theme(10) + ggplot2::theme(strip.text = ggplot2::element_text(size = 9.2))
  p$row_type <- "individual_mouse"
  m$row_type <- "fitted_group_mean"
  list(plot = plot, source = cfdi_bind(p, m))
}

cfdi_plot_lipids <- function(data) {
  d <- data
  d$significant <- d$q_interaction <= .05
  d$label <- FALSE
  selected <- which(d$significant & is.finite(d$q_interaction))
  selected <- selected[order(d$q_interaction[selected], d$p_interaction[selected], d$lipid[selected])]
  d$label[head(selected, 6)] <- TRUE
  extent <- max(abs(c(d$effect_female, d$effect_male)), na.rm = TRUE)
  limits <- c(-1, 1) * max(.5, extent * 1.14)
  count <- paste0(sum(d$significant), " / ", format(nrow(d), big.mark = ",", trim = TRUE),
    " lipids: diet x sex FDR <= 0.05")
  plot <- ggplot2::ggplot(d, ggplot2::aes(effect_female, effect_male)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = .5, color = "#AFAFAF") +
    ggplot2::geom_hline(yintercept = 0, linewidth = .3, color = "#D0D0D0") +
    ggplot2::geom_vline(xintercept = 0, linewidth = .3, color = "#D0D0D0") +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 1.4, alpha = .7) +
    ggplot2::annotate("text", x = limits[1], y = limits[2], label = count,
      hjust = -.02, vjust = 1.25, size = 3.1, color = "#444444", family = "Helvetica") +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#BFBFBF", `TRUE` = "#76517B"),
      labels = c(`FALSE` = "FDR > 0.05", `TRUE` = "FDR <= 0.05"), name = "Diet x sex interaction") +
    ggplot2::coord_equal(xlim = limits, ylim = limits, clip = "off") +
    ggplot2::labs(title = "Diet effects on whole-liver lipids",
      subtitle = "Deviation from the identity line indicates a sex difference in diet response",
      x = "Female WD-minus-SD lipid effect (log2)", y = "Male WD-minus-SD lipid effect (log2)") +
    cfdi_theme(10) + ggplot2::theme(axis.title.x = ggplot2::element_text(color = "#A43F44"),
                                   axis.title.y = ggplot2::element_text(color = "#305F9C"))
  if (any(d$label)) plot <- plot + ggrepel::geom_label_repel(
    data = d[d$label, , drop = FALSE], ggplot2::aes(label = lipid_label),
    seed = 1506, size = 2.8, color = "#53365A", fill = "white", label.size = .1,
    label.padding = grid::unit(.1, "lines"), point.padding = .35, box.padding = .4,
    min.segment.length = 0, segment.size = .25, max.overlaps = Inf, family = "Helvetica")
  list(plot = plot, source = d)
}

cfdi_save_panel <- function(result, output, stem, width, height, family) {
  panel_dir <- file.path(output, "panels")
  source_dir <- file.path(output, "source_data")
  if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
  if (!dir.exists(source_dir)) dir.create(source_dir, recursive = TRUE)
  source <- file.path(source_dir, paste0(stem, ".csv"))
  utils::write.csv(result$source, source, row.names = FALSE, na = "")
  files <- file.path(panel_dir, paste0(stem, c(".pdf", ".svg", ".png")))
  ggplot2::ggsave(files[1], plot = result$plot, width = width, height = height, units = "in",
    device = grDevices::pdf, family = "Helvetica", bg = "white", useDingbats = FALSE)
  ggplot2::ggsave(files[2], plot = result$plot, width = width, height = height, units = "in",
    device = svglite::svglite, bg = "white")
  ggplot2::ggsave(files[3], plot = result$plot, width = width, height = height, units = "in",
    device = ragg::agg_png, dpi = 300, bg = "white")
  data.frame(panel = stem, family = family, file = normalizePath(files, mustWork = TRUE),
    format = c("PDF", "SVG", "PNG"), width_in = width, height_in = height,
    source_table = normalizePath(source, mustWork = TRUE), stringsAsFactors = FALSE)
}

cfdi_make_panels <- function(a, output) {
  absent <- setdiff(c("overview", "heatmap", "example_points", "example_means", "lipid_scatter"), names(a))
  if (length(absent)) stop("Missing compact panel inputs: ", paste(absent, collapse = ", "))
  a_plot <- cfdi_plot_overview(a$overview)
  b_plot <- cfdi_plot_heatmap(a$heatmap)
  c_plot <- cfdi_plot_examples(a$example_points, a$example_means)
  d_plot <- cfdi_plot_lipids(a$lipid_scatter)
  n_genes <- if (nrow(a$heatmap)) length(unique(a$heatmap$gene)) else 1L
  manifest <- rbind(
    cfdi_save_panel(a_plot, output, "Fig6A_direct_RNA_interactions", 10.2, 4.7, "A"),
    cfdi_save_panel(b_plot, output, "Fig6B_conditional_diet_effects", 9.3, max(4.7, 2.2 + .23 * n_genes), "B"),
    cfdi_save_panel(c_plot, output, "Fig6C_mouse_examples", 8.4, 6.8, "C"),
    cfdi_save_panel(d_plot, output, "Fig6D_lipid_diet_effects", 6.6, 6.4, "D")
  )
  rownames(manifest) <- NULL
  if (nrow(manifest) != 12L || length(unique(manifest$panel)) != 4L) stop("Compact output must have exactly four panel groups.")
  utils::write.csv(manifest, file.path(output, "panel_manifest.csv"), row.names = FALSE)
  manifest
}
