# Figure 6A paired-dot alternative. Source the common layout helpers first.
# This module creates plot objects only; it does not fit models or export files.
cf6a_dots <- function(d) {
  required <- c("gene", "context", "sex", "zone", "selection", "display_rank",
                "effect_female_SD", "effect_female_WD", "effect_male_SD", "effect_male_WD",
                "threeway_q", "direct_supported")
  if (length(setdiff(required, names(d)))) stop("Missing Figure 6A dot-layout fields.")
  if (!exists("cf6a_long", mode = "function") || !exists("cf6a_theme", mode = "function"))
    stop("Source scripts/core/cf453_figure6a_layouts.R before calling cf6a_dots().")
  contexts <- c("female_PP", "male_PP", "female_PC", "male_PC")
  headings <- c(female_PP = "Selected female periportal genes",
                male_PP = "Selected male periportal genes",
                female_PC = "Selected female pericentral genes",
                male_PC = "Selected male pericentral genes")
  stopifnot(nrow(d) > 0L, !anyDuplicated(d$gene),
            all(d$context %in% contexts), all(d$selection %in% c("SD", "WD", "Both")),
            all(d$direct_supported == (d$threeway_q <= .05)))
  effect_columns <- c("effect_female_SD", "effect_female_WD", "effect_male_SD", "effect_male_WD")
  effect_values <- unlist(d[effect_columns], use.names = FALSE)
  if (any(!is.finite(effect_values))) stop("Dot layout requires four finite effects per gene.")
  effect_range <- range(effect_values)
  span <- diff(effect_range)
  if (span <= 0) span <- max(abs(effect_range), 1)
  point_limits <- effect_range + c(-1, 1) * .045 * span
  # A shared label gutter keeps gene names and selection badges outside the
  # numeric point field, while every block retains the same effect axis.
  gene_x <- effect_range[1L] - .62 * span
  badge_x <- effect_range[1L] - .155 * span
  plot_limits <- c(gene_x - .025 * span, point_limits[2L])
  d$programme <- factor(d$context, levels = contexts)
  d$programme_rows <- as.integer(ave(d$display_rank, d$context, FUN = max))
  for (context in unique(d$context)) {
    ranks <- d$display_rank[d$context == context]
    if (!identical(sort(as.integer(ranks)), seq_along(ranks)))
      stop("Display ranks must be consecutive and unique within each programme.")
  }
  d$y_gene <- d$programme_rows + 1 - d$display_rank
  d$gene_label <- ifelse(d$direct_supported, paste0(d$gene, " *"), d$gene)
  long <- cf6a_long(d)
  long$sex <- factor(long$sex, levels = c("female", "male"))
  long$diet <- factor(long$diet, levels = c("SD", "WD"))
  long$track_offset <- ifelse(long$sex == "female", .145, -.145)
  long$y <- long$y_gene + long$track_offset
  long$q <- long$threeway_q
  long$common_effect_min <- effect_range[1L]
  long$common_effect_max <- effect_range[2L]
  rownames(long) <- NULL
  segments <- long[long$diet == "SD", , drop = FALSE]
  segments$effect_end <- ifelse(segments$sex == "female", segments$effect_female_WD,
                               segments$effect_male_WD)
  direction <- d[!duplicated(d$context), c("context", "programme"), drop = FALSE]
  direction$y <- .16
  left <- direction; left$x <- point_limits[1L]; left$label <- "PC-biased"; left$hjust <- 0
  right <- direction; right$x <- point_limits[2L]; right$label <- "PP-biased"; right$hjust <- 1
  direction <- rbind(left, right)
  ticks <- pretty(effect_range, n = 5L)
  ticks <- ticks[ticks >= point_limits[1L] & ticks <= point_limits[2L]]
  if (!any(ticks == 0) && point_limits[1L] < 0 && point_limits[2L] > 0) ticks <- sort(c(ticks, 0))
  colors <- c(female = "#C9464C", male = "#346DB3")
  plot <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linewidth = .35, color = "#C5C5C5") +
    ggplot2::geom_segment(data = segments,
      ggplot2::aes(x = effect, xend = effect_end, y = y, yend = y, color = sex),
      linewidth = .55, alpha = .8, show.legend = FALSE) +
    ggplot2::geom_point(data = long,
      ggplot2::aes(x = effect, y = y, color = sex, shape = diet), size = 2.05, stroke = .65) +
    ggplot2::geom_text(data = d,
      ggplot2::aes(x = gene_x, y = y_gene, label = gene_label), hjust = 0,
      fontface = "italic", size = 2.85, color = "#252525", family = "Helvetica") +
    ggplot2::geom_label(data = d,
      ggplot2::aes(x = badge_x, y = y_gene, label = selection),
      size = 2.15, label.padding = grid::unit(.55, "mm"),
      label.r = grid::unit(.45, "mm"), linewidth = .2,
      fill = "#F5F5F5", color = "#555555", family = "Helvetica") +
    ggplot2::geom_text(data = direction,
      ggplot2::aes(x = x, y = y, label = label, hjust = hjust),
      size = 2.35, color = "#6B6B6B", family = "Helvetica") +
    ggplot2::scale_color_manual(values = colors, breaks = names(colors),
      labels = c(female = "Female", male = "Male"), name = "RNA effect") +
    ggplot2::scale_shape_manual(values = c(SD = 1, WD = 16), name = "Diet") +
    ggplot2::scale_x_continuous(limits = plot_limits, breaks = ticks,
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(add = c(.35, .55))) +
    ggplot2::facet_wrap(~programme, ncol = 2L, scales = "free_y", drop = FALSE,
      labeller = ggplot2::as_labeller(headings), axes = "all_x") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "Diet shifts in selected zonation programmes",
      subtitle = "Both sexes shown for every gene; open SD to filled WD on a common effect scale",
      x = "Periportal - pericentral effect (log-normalized expression)", y = NULL,
      caption = "Selection badge: SD, WD or Both. * Direct sex-by-diet-by-zone q <= 0.05.\nLines join fitted SD and WD effects; diet groups contain different mice.") +
    cf6a_theme(10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8.5),
      axis.title.x = ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 9)),
      axis.line.x = ggplot2::element_line(linewidth = .3, color = "#9C9C9C"),
      strip.background = ggplot2::element_rect(fill = "#F1F2F3", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 10, hjust = 0,
        margin = ggplot2::margin(6, 5, 6, 5)),
      panel.spacing.x = grid::unit(9, "mm"), panel.spacing.y = grid::unit(7, "mm"),
      legend.box = "horizontal", legend.spacing.x = grid::unit(7, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      plot.caption = ggplot2::element_text(size = 8.2, hjust = 0, lineheight = 1.2)) +
    ggplot2::guides(color = ggplot2::guide_legend(order = 1, override.aes = list(shape = 16)),
      shape = ggplot2::guide_legend(order = 2, override.aes = list(color = "#444444")))
  list(plot = plot, source = long)
}
