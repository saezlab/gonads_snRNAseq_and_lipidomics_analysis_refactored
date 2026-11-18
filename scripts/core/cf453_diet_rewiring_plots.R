# Three independent Figure 5-style panels from standalone CF453 analysis tables.
# No model fitting, gene selection or network selection is performed here.

cfdr_colors <- c(female_PC = "#B83338", female_PP = "#ED8585",
                 male_PC = "#285CA8", male_PP = "#80A8E0")
cfdr_context_labels <- c(female_PC = "Female pericentral", female_PP = "Female periportal",
                         male_PC = "Male pericentral", male_PP = "Male periportal")

cfdr_context <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[ :.-]+", "_", x)
  x <- sub("_(periportal|portal|pp)$", "_PP", x)
  sub("_(pericentral|central|pc)$", "_PC", x)
}

cfdr_theme <- function(size = 11) {
  ggplot2::theme_classic(base_size = size, base_family = "Helvetica") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = size - 1, color = "#4B4B4B", margin = ggplot2::margin(b = 9)),
      plot.caption = ggplot2::element_text(size = size - 2, color = "#555555", hjust = 0, lineheight = 1.25),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.text = ggplot2::element_text(color = "#303030"),
      axis.line = ggplot2::element_line(linewidth = .35, color = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = .3),
      strip.background = ggplot2::element_rect(fill = "#F3F3F3", color = NA),
      strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(7, 5, 7, 5)),
      legend.title = ggplot2::element_text(face = "bold", size = size - 1),
      legend.text = ggplot2::element_text(size = size - 1),
      legend.position = "bottom", panel.spacing = grid::unit(1.1, "lines"),
      plot.margin = ggplot2::margin(10, 14, 10, 10))
}

cfdr_bind <- function(...) {
  xs <- list(...)
  fields <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) {
    for (field in setdiff(fields, names(x))) x[[field]] <- rep(NA, nrow(x))
    x[, fields, drop = FALSE]
  })
  do.call(rbind, xs)
}

cfdr_plot_anchors <- function(data) {
  d <- data
  d$context <- cfdr_context(d$context)
  if (any(!d$context %in% names(cfdr_colors))) stop("Unrecognised anchor context.")
  d$anchor_diet <- toupper(d$anchor_diet)
  if (any(!d$anchor_diet %in% c("SD", "WD"))) stop("Anchor diet must be SD or WD.")
  d$x_start <- ifelse(d$anchor_diet == "SD", d$x_SD, d$x_WD)
  d$y_start <- ifelse(d$anchor_diet == "SD", d$y_SD, d$y_WD)
  d$x_end <- ifelse(d$anchor_diet == "SD", d$x_WD, d$x_SD)
  d$y_end <- ifelse(d$anchor_diet == "SD", d$y_WD, d$y_SD)
  d$direct_test <- factor(ifelse(d$threeway_q <= .05, "supported", "other"), levels = c("supported", "other"))
  d$anchor_diet <- factor(d$anchor_diet, levels = c("SD", "WD"))
  counts <- table(d$anchor_diet)
  facet_labels <- c(SD = paste0("SD-selected (n = ", counts[["SD"]], "): SD to WD"),
                     WD = paste0("WD-selected (n = ", counts[["WD"]], "): WD to SD"))
  sd <- d; sd$x <- sd$x_SD; sd$y <- sd$y_SD; sd$endpoint <- "SD"
  wd <- d; wd$x <- wd$x_WD; wd$y <- wd$y_WD; wd$endpoint <- "WD"
  endpoints <- rbind(sd, wd)
  extent <- max(abs(c(d$x_SD, d$y_SD, d$x_WD, d$y_WD)), na.rm = TRUE)
  limits <- c(-1, 1) * max(1, extent * 1.20)
  fills <- vapply(cfdr_colors, grDevices::adjustcolor, character(1), alpha.f = .12)
  plot <- ggplot2::ggplot(d) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#E0E0E0", linewidth = .35) +
    ggplot2::geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = .35) +
    ggplot2::geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = .35) +
    ggplot2::geom_segment(ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                                      color = context, linetype = direct_test),
                           linewidth = .55, alpha = .8,
                           arrow = grid::arrow(length = grid::unit(2.4, "mm"), type = "closed")) +
    ggplot2::geom_point(data = endpoints, ggplot2::aes(x, y, color = context, shape = endpoint),
                         size = 2.5, stroke = .7) +
    ggrepel::geom_label_repel(ggplot2::aes(x_start, y_start, label = gene, fill = context),
      color = "#303030", size = 2.8, seed = 1606, box.padding = .35, point.padding = .45,
      min.segment.length = 0, segment.size = .25, segment.color = "#888888", label.size = .1,
      label.padding = grid::unit(.08, "lines"), label.r = grid::unit(.06, "lines"),
      force = 2, max.overlaps = Inf, max.time = 3, max.iter = 10000,
      show.legend = FALSE, family = "Helvetica") +
    ggplot2::scale_color_manual(values = cfdr_colors, labels = cfdr_context_labels,
                                breaks = names(cfdr_colors), name = "Anchor programme", drop = FALSE) +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_shape_manual(values = c(SD = 1, WD = 16), name = "Endpoints") +
    ggplot2::scale_linetype_manual(values = c(supported = "solid", other = "22"),
      labels = c(supported = "q <= 0.05", other = "q > 0.05"), name = "Direct three-way test") +
    ggplot2::facet_wrap(~anchor_diet, nrow = 1, labeller = ggplot2::as_labeller(facet_labels)) +
    ggplot2::coord_equal(xlim = limits, ylim = limits, clip = "off") +
    ggplot2::labs(title = "Diet shifts in selected sex-specific zonation genes",
      subtitle = "Arrows start in the selection diet; positive t = periportal, negative t = pericentral",
      x = "Female zonation bias (t-statistic)", y = "Male zonation bias (t-statistic)",
      caption = "SD, standard diet; WD, Western diet. Arrow style reports the direct three-way model test.") +
    cfdr_theme(10) + ggplot2::theme(legend.box = "vertical") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, order = 1,
                      override.aes = list(linetype = "solid", shape = 16, alpha = 1)),
      shape = ggplot2::guide_legend(nrow = 1, order = 2, override.aes = list(color = "#333333")),
      linetype = ggplot2::guide_legend(nrow = 1, order = 3, override.aes = list(color = "#333333")))
  list(plot = plot, source = d)
}

cfdr_plot_class_shift <- function(data) {
  d <- data
  contexts <- c("female_PC", "female_PP", "male_PC", "male_PP")
  d$context <- cfdr_context(d$context)
  if (any(!d$context %in% contexts)) stop("Unrecognised class-summary context.")
  if (anyDuplicated(paste(d$context, d$lipid_class))) stop("Class-summary cells must be unique.")
  classes <- rev(sort(unique(as.character(d$lipid_class))))
  d$x <- match(d$context, contexts)
  d$y <- match(d$lipid_class, classes)
  d$delta_percentage_points <- 100 * d$delta_E
  d$stability_marker <- is.finite(d$bootstrap_direction) & d$bootstrap_direction >= .90 & d$bootstrap_valid > 0
  d$label <- ifelse(is.finite(d$delta_percentage_points), sprintf("%+.1f", d$delta_percentage_points), "NA")
  extent <- max(abs(d$delta_percentage_points), na.rm = TRUE)
  extent <- max(extent, .1)
  d$text_color <- ifelse(abs(d$delta_percentage_points) >= .6 * extent, "white", "#303030")
  denominators <- unique(d[, c("lipid_class", "n_lipids")])
  if (anyDuplicated(denominators$lipid_class)) stop("Lipid-class species counts differ across fixed contexts.")
  row_labels <- paste0(classes, "  (n = ", denominators$n_lipids[match(classes, denominators$lipid_class)], ")")
  plot <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = delta_percentage_points), width = .96, height = .96,
                        color = "white", linewidth = .3) +
    ggplot2::geom_text(ggplot2::aes(label = label, color = text_color), size = 2.8, family = "Helvetica") +
    ggplot2::geom_point(data = d[d$stability_marker, , drop = FALSE],
      ggplot2::aes(x = x + .33, y = y + .29), shape = 21, size = 1.6, stroke = .45,
      fill = "white", color = "#303030") +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(low = "#C64E54", mid = "#FAF9F7", high = "#467CBC", midpoint = 0,
      limits = c(-extent, extent), name = "WD - SD change in\nclass excess (pp)", na.value = "#DDDDDD",
      guide = ggplot2::guide_colourbar(display = "rectangles", barheight = grid::unit(35, "mm"))) +
    ggplot2::scale_x_continuous(breaks = seq_along(contexts),
      labels = c("Female PC", "Female PP", "Male PC", "Male PP"), expand = c(.015, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(classes), labels = row_labels, expand = c(.015, 0)) +
    ggplot2::labs(title = "Diet change in lipid-class associations",
      subtitle = "Class-minus-other-class pair-hit rates; fixed anchor-selected gene sets",
      x = NULL, y = NULL,
      caption = "Circle: >=90% bootstrap directional stability among defined draws. Values are percentage points.") +
    cfdr_theme(10) + ggplot2::theme(axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), legend.position = "right")
  list(plot = plot, source = d)
}

cfdr_lipid_label <- function(x) {
  x <- as.character(x)
  long <- nchar(x) > 20
  x[long] <- gsub("([_/])", "\\1\n", x[long])
  x
}

cfdr_plot_networks <- function(edges, genes) {
  order <- c("female_PP", "male_PP", "female_PC", "male_PC")
  facets <- data.frame(context = factor(order, levels = order), stringsAsFactors = FALSE)
  e <- edges
  e$context <- cfdr_context(e$context)
  g <- genes
  g$context <- cfdr_context(g$context)
  g$context <- factor(g$context, levels = order)
  e$context <- factor(e$context, levels = order)
  e$plot_direction <- tolower(as.character(e$direction))
  e$plot_direction[e$plot_direction %in% c("gain", "gained")] <- "gained"
  e$plot_direction[e$plot_direction %in% c("loss", "lost")] <- "lost"
  if (any(!e$plot_direction %in% c("gained", "lost"))) stop("Network examples require gained or lost edges.")
  e$slot <- rep(NA_integer_, nrow(e))
  for (context in order) {
    indices <- which(e$context == context)
    if (length(indices) > 6L) stop("At most six selected partners may be drawn per context.")
    if (length(unique(e$gene[indices])) > 1L) stop("Each network context must have one selected gene.")
    for (direction in c("gained", "lost")) {
      ii <- indices[e$plot_direction[indices] == direction]
      ii <- ii[order(-abs(e$delta_abs[ii]), e$lipid[ii])]
      ii <- head(ii, 3L)
      e$slot[ii] <- if (direction == "gained") seq_along(ii) else 3L + seq_along(ii)
    }
    remaining <- indices[is.na(e$slot[indices])]
    remaining <- remaining[order(-abs(e$delta_abs[remaining]), e$lipid[remaining])]
    e$slot[remaining] <- head(setdiff(seq_len(6), e$slot[indices]), length(remaining))
  }
  e$node_x <- c(.17, .5, .83, .17, .5, .83)[e$slot]
  e$node_y <- ifelse(e$slot <= 3, .5, .12)
  sd_label <- formatC(e$r_SD, digits = 2, format = "f")
  wd_label <- formatC(e$r_WD, digits = 2, format = "f")
  sd_label[sd_label == "-0.00"] <- "0.00"
  wd_label[wd_label == "-0.00"] <- "0.00"
  e$edge_label <- if (nrow(e)) paste0(sd_label, " -> ", wd_label) else character(0)
  e$node_label <- cfdr_lipid_label(e$lipid_label)
  e$node_fill <- vapply(as.character(e$context), function(z) grDevices::adjustcolor(cfdr_colors[[z]], alpha.f = .14), character(1))
  e$edge_width <- abs(e$delta_abs)
  path <- lapply(seq_len(nrow(e)), function(i) {
    k <- e$slot[i]
    coordinates <- if (k <= 3L) {
      data.frame(path_x = c(c(.46, .5, .54)[k], e$node_x[i]), path_y = c(.84, .61))
    } else if (k == 4L) {
      data.frame(path_x = c(.43, -.015, -.015, .17), path_y = c(.84, .72, .29, .23))
    } else if (k == 5L) {
      data.frame(path_x = c(.52, .665, .665, .5), path_y = c(.84, .72, .29, .23))
    } else {
      data.frame(path_x = c(.57, 1.015, 1.015, .83), path_y = c(.84, .72, .29, .23))
    }
    cbind(e[rep(i, nrow(coordinates)), , drop = FALSE], coordinates, edge_id = i)
  })
  paths <- if (length(path)) do.call(rbind, path) else data.frame()
  has_gene <- !is.na(g$gene) & nzchar(g$gene) & g$context %in% e$context
  selected_genes <- g[has_gene, , drop = FALSE]
  selected_genes$node_fill <- unname(cfdr_colors[as.character(selected_genes$context)])
  empty <- facets[!facets$context %in% e$context, , drop = FALSE]
  empty$message <- rep("No supported gain/loss example", nrow(empty))
  max_width <- if (nrow(e)) max(e$edge_width, na.rm = TRUE) else 1
  plot <- ggplot2::ggplot(facets) +
    ggplot2::geom_blank(ggplot2::aes(x = .5, y = .5))
  if (nrow(e)) {
    plot <- plot +
      ggplot2::geom_path(data = paths, ggplot2::aes(path_x, path_y, group = edge_id,
        color = plot_direction, linewidth = edge_width), lineend = "round", linejoin = "round") +
      ggplot2::geom_rect(data = e, ggplot2::aes(xmin = node_x - .145, xmax = node_x + .145,
        ymin = node_y - .09, ymax = node_y + .09, fill = node_fill, color = context), linewidth = .35) +
      ggplot2::geom_text(data = e, ggplot2::aes(node_x, node_y, label = node_label),
        size = 2.7, color = "#252525", lineheight = 1.05, family = "Helvetica") +
      ggplot2::geom_label(data = e, ggplot2::aes(node_x, node_y + .16, label = edge_label),
        fill = "white", color = "#333333", linewidth = 0, size = 2.6,
        label.padding = grid::unit(.04, "lines"), family = "Helvetica") +
      ggplot2::geom_label(data = selected_genes, ggplot2::aes(x = .5, y = .91, label = gene, fill = node_fill),
        color = "white", fontface = "italic", linewidth = 0, size = 3.8,
        label.r = grid::unit(5, "mm"), label.padding = grid::unit(3.2, "mm"), family = "Helvetica")
  }
  if (nrow(empty)) plot <- plot + ggplot2::geom_text(data = empty,
    ggplot2::aes(x = .5, y = .5, label = message), color = "#666666", size = 3.5, family = "Helvetica")
  plot <- plot +
    ggplot2::scale_color_manual(values = c(cfdr_colors, gained = "#238D7B", lost = "#C7793B"),
      breaks = c("gained", "lost"), labels = c("Gained", "Lost"), name = "Association") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_linewidth_continuous(range = c(.1, 1.6), limits = c(0, max(max_width, .1)),
      name = "Change in |r|", breaks = pretty(c(0, max(max_width, .1)), n = 2)) +
    ggplot2::facet_wrap(~context, ncol = 2, drop = FALSE,
      labeller = ggplot2::as_labeller(cfdr_context_labels)) +
    ggplot2::coord_cartesian(xlim = c(-.07, 1.07), ylim = c(-.015, 1.02), clip = "off") +
    ggplot2::labs(title = "Diet-associated gains and losses in gene-lipid links",
      subtitle = "Selected examples within the same sex and hepatocyte zone",
      caption = "Exploratory examples with bootstrap and leave-one-mouse-out support.\nEdge labels: r(SD) -> r(WD); line width: absolute change in |r|.") +
    cfdr_theme(10) + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      legend.key.width = grid::unit(9, "mm"), panel.spacing = grid::unit(1.3, "lines"))
  e$row_type <- rep("selected_edge", nrow(e))
  g$row_type <- rep("context_representative", nrow(g))
  list(plot = plot, source = cfdr_bind(e, g))
}

cfdr_save_panel <- function(result, output, stem, width, height, family) {
  panel_dir <- file.path(output, "panels"); source_dir <- file.path(output, "source_data")
  if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
  if (!dir.exists(source_dir)) dir.create(source_dir, recursive = TRUE)
  source <- file.path(source_dir, paste0(stem, ".csv"))
  utils::write.csv(result$source, source, row.names = FALSE, na = "")
  files <- file.path(panel_dir, paste0(stem, c(".pdf", ".svg", ".png")))
  ggplot2::ggsave(files[1], plot = result$plot, width = width, height = height, units = "in",
    device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE, bg = "white")
  ggplot2::ggsave(files[2], plot = result$plot, width = width, height = height, units = "in",
    device = svglite::svglite, bg = "white")
  ggplot2::ggsave(files[3], plot = result$plot, width = width, height = height, units = "in",
    device = ragg::agg_png, dpi = 300, bg = "white")
  data.frame(panel = stem, family = family, file = normalizePath(files, mustWork = TRUE),
    format = c("PDF", "SVG", "PNG"), width_in = width, height_in = height,
    source_table = normalizePath(source, mustWork = TRUE), stringsAsFactors = FALSE)
}

cfdr_make_panels <- function(a, output) {
  required <- c("anchors", "class_summary", "network_genes", "network_edges")
  if (!all(required %in% names(a))) stop("Missing rewiring panel inputs: ", paste(setdiff(required, names(a)), collapse = ", "))
  result_a <- cfdr_plot_anchors(a$anchors)
  result_b <- cfdr_plot_class_shift(a$class_summary)
  result_c <- cfdr_plot_networks(a$network_edges, a$network_genes)
  n_classes <- length(unique(a$class_summary$lipid_class))
  manifest <- rbind(
    cfdr_save_panel(result_a, output, "Fig6A_anchor_gene_shifts", 11.4, 7.4, "A"),
    cfdr_save_panel(result_b, output, "Fig6B_lipid_class_shift", 7.8, max(5.5, 2.0 + .26 * n_classes), "B"),
    cfdr_save_panel(result_c, output, "Fig6C_network_gain_loss", 10.0, 8.7, "C"))
  rownames(manifest) <- NULL
  if (nrow(manifest) != 9L || length(unique(manifest$panel)) != 3L) stop("Exactly three panel groups are required.")
  utils::write.csv(manifest, file.path(output, "panel_manifest.csv"), row.names = FALSE)
  manifest
}
