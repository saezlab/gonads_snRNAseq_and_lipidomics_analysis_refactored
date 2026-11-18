# Read-only layout alternative for the saved Workflow 16 Figure 6A effects.
# The shared cf453_figure6a_layouts.R module supplies ordering and source data.
cf6a_heatmap <- function(d) {
  required <- c("gene", "context", "selection", "display_rank", "global_rank",
    "effect_female_SD", "effect_female_WD", "effect_male_SD", "effect_male_WD",
    "threeway_q", "direct_supported")
  if (!all(required %in% names(d))) stop("Missing Figure 6A heatmap fields.")
  if (anyDuplicated(d$gene) || any(!d$context %in% cf6a_contexts) ||
      any(!d$selection %in% c("SD", "WD", "Both"))) {
    stop("The heatmap requires one row per gene and valid programme/selection labels.")
  }
  d <- d[order(match(d$context, cf6a_contexts), d$display_rank, d$gene), , drop = FALSE]
  if (!identical(as.integer(d$global_rank), seq_len(nrow(d)))) {
    stop("Heatmap input order differs from the shared gene ledger.")
  }
  # Extra space between programmes leaves an aligned heading above each block.
  d$display_y <- -d$global_rank - 1.8 * (match(d$context, cf6a_contexts) - 1L)
  source <- cf6a_long(d)
  state_positions <- c(female_SD = 0, female_WD = 1, male_SD = 2.35, male_WD = 3.35)
  source$display_x <- unname(state_positions[source$state])
  if (nrow(source) != 4L * nrow(d) || any(!is.finite(source$effect))) {
    stop("Expected four finite fitted zonation effects per gene.")
  }
  colour_limit <- max(abs(source$effect))
  if (!is.finite(colour_limit) || colour_limit <= 0) stop("Invalid heatmap colour range.")
  source$colour_limit <- colour_limit
  source$effect_units <- "Difference in log-normalized RNA expression (PP - PC)"

  programme <- do.call(rbind, lapply(cf6a_contexts, function(context) {
    z <- d[d$context == context, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    data.frame(context = context, heading_y = max(z$display_y) + 0.94,
      top = max(z$display_y) + 0.43, bottom = min(z$display_y) - 0.43,
      label = unname(cf6a_context_labels[context]), stringsAsFactors = FALSE)
  }))
  programme_colours <- c(female_PP = "#ED8585", female_PC = "#B83338",
    male_PP = "#80A8E0", male_PC = "#285CA8")
  supported <- d[!is.na(d$direct_supported) & d$direct_supported, , drop = FALSE]
  p <- ggplot2::ggplot(source, ggplot2::aes(x = display_x, y = display_y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = effect), width = 0.94, height = 0.87,
      colour = "white", linewidth = 0.25) +
    ggplot2::geom_segment(data = programme,
      ggplot2::aes(x = -2.40, xend = -2.40, y = bottom, yend = top, colour = context),
      linewidth = 2.0, inherit.aes = FALSE, lineend = "butt") +
    ggplot2::geom_text(data = programme,
      ggplot2::aes(x = -2.22, y = heading_y, label = label),
      size = 3.0, fontface = "bold", hjust = 0, colour = "#303030",
      inherit.aes = FALSE) +
    ggplot2::geom_text(data = d,
      ggplot2::aes(x = -2.22, y = display_y, label = gene),
      size = 3.0, fontface = "italic", hjust = 0, colour = "#252525", inherit.aes = FALSE) +
    ggplot2::geom_text(data = d,
      ggplot2::aes(x = 4.55, y = display_y, label = selection),
      size = 2.8, colour = "#505050", inherit.aes = FALSE) +
    ggplot2::geom_point(data = supported,
      ggplot2::aes(x = 5.45, y = display_y), shape = 16, size = 1.85,
      colour = "#252525", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = c(0.5, 2.85), y = 2.10,
      label = c("Female", "Male"), fontface = "bold", size = 3.3, colour = "#303030") +
    ggplot2::annotate("text", x = unname(state_positions), y = 1.12,
      label = c("SD", "WD", "SD", "WD"), size = 3.0, colour = "#303030") +
    ggplot2::annotate("text", x = c(4.55, 5.45), y = 2.10,
      label = c("Selected", "Direct q"), fontface = "bold", size = 2.8, colour = "#303030") +
    ggplot2::annotate("text", x = c(4.55, 5.45), y = 1.12,
      label = c("in diet", "<= 0.05"), size = 2.8, colour = "#505050") +
    ggplot2::scale_fill_gradient2(low = "#70508C", mid = "#FAF9F5", high = "#C99A2E",
      midpoint = 0, limits = c(-colour_limit, colour_limit),
      breaks = scales::breaks_pretty(n = 5),
      name = "PP - PC expression difference (log-normalized units)",
      guide = ggplot2::guide_colourbar(title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(86, "mm"), barheight = grid::unit(3.6, "mm"),
        display = "rectangles", ticks = TRUE)) +
    ggplot2::scale_colour_manual(values = programme_colours, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-2.55, 5.94), expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(min(d$display_y) - 0.85, 2.85),
      expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(title = "Sex-specific zonation across diets",
      subtitle = paste(nrow(d), "genes selected for sex-specific zonation in either diet"),
      x = NULL, y = NULL,
      caption = paste(
        "SD, standard diet; WD, Western diet. PP, periportal; PC, pericentral.",
        "Positive effects favour PP; negative effects favour PC. Common colour scale; no row scaling.",
        "Dot: direct diet x sex x zonation interaction q <= 0.05.",
        "Selected: anchor-diet membership.",
        sep = "\n")) +
    cf6a_theme(11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#555555",
        margin = ggplot2::margin(b = 9)),
      plot.caption = ggplot2::element_text(size = 8.1, colour = "#555555", hjust = 0,
        lineheight = 1.15, margin = ggplot2::margin(t = 8)),
      legend.title = ggplot2::element_text(size = 9, face = "plain"),
      legend.text = ggplot2::element_text(size = 8.5),
      legend.margin = ggplot2::margin(t = 5, b = 0),
      legend.box.spacing = grid::unit(0, "mm"),
      plot.margin = ggplot2::margin(10, 12, 10, 12))
  list(plot = p, source = source)
}
