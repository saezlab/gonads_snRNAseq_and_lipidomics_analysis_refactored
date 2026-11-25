# Manuscript Figure 6B and 6E from saved CF453 results.
# Figure numbers here refer to the manuscript, and upstream CSV names retain Workflow 16 numbering.

cf453_manuscript_lipid_contrasts <- function(tests) {
  required <- c("lipid", "lipid_label", "lipid_class", "contrast", "family", "effect", "t", "p", "q")
  stopifnot(all(required %in% names(tests)), !anyDuplicated(tests[c("lipid", "contrast")]))
  interaction <- subset(tests, contrast == "sex_diet")
  sex_contrasts <- lapply(c("diet_female", "diet_male"), function(contrast_name) {
    selected <- tests[tests$contrast == contrast_name, ]
    stopifnot(setequal(selected$lipid, interaction$lipid))
    selected[match(interaction$lipid, selected$lipid), ]
  })
  female <- sex_contrasts[[1L]]
  male <- sex_contrasts[[2L]]
  stopifnot(
    nrow(interaction) > 0L, all(interaction$family == "lipid_interaction"),
    identical(female$lipid_label, interaction$lipid_label),
    identical(male$lipid_label, interaction$lipid_label),
    identical(female$lipid_class, interaction$lipid_class),
    identical(male$lipid_class, interaction$lipid_class),
    all(is.finite(c(female$t, male$t, interaction$t))),
    all(is.finite(c(female$effect, male$effect, interaction$effect))),
    max(abs(interaction$effect - (male$effect - female$effect))) < 1e-12
  )
  data.frame(
    lipid = interaction$lipid, lipid_label = interaction$lipid_label,
    lipid_class = interaction$lipid_class,
    female_log2_WD_SD = female$effect, female_t = female$t, female_BH_q = female$q,
    male_log2_WD_SD = male$effect, male_t = male$t, male_BH_q = male$q,
    diet_sex_male_minus_female = interaction$effect, interaction_t = interaction$t,
    interaction_p = interaction$p, interaction_BH_q = interaction$q
  )
}

cf453_manuscript_figure6b <- function(lipids) {
  # Fixed display thresholds transferred from the expanded-marker Figure 5A rule.
  interaction_cutoff <- 4.4
  other_cutoff <- 1
  interaction_gate <- abs(lipids$interaction_t) > interaction_cutoff
  female_highlight <- interaction_gate & abs(lipids$male_t) < other_cutoff & lipids$female_t != 0
  male_highlight <- interaction_gate & abs(lipids$female_t) < other_cutoff & lipids$male_t != 0
  stopifnot(!any(female_highlight & male_highlight))

  categories <- c("female_lower", "female_higher", "male_lower", "male_higher")
  category <- rep("other", nrow(lipids))
  category[female_highlight] <- ifelse(lipids$female_t[female_highlight] < 0, "female_lower", "female_higher")
  category[male_highlight] <- ifelse(lipids$male_t[male_highlight] < 0, "male_lower", "male_higher")
  lipids$interaction_abs_t_cutoff <- interaction_cutoff
  lipids$other_sex_abs_t_cutoff <- other_cutoff
  lipids$category <- factor(category, levels = c(categories, "other"))

  counts <- table(lipids$category)
  legend_labels <- sprintf("%s (%d)", c("Female lower", "Female higher",
    "Male lower", "Male higher"), counts[categories])
  category_colors <- stats::setNames(unname(ZONATION_GONAD_COLORS[c("OPC", "OPP", "TPC", "TPP")]), categories)

  # Approved label placement at 64 mm. Use a representative per populated group on other inputs.
  label_positions <- data.frame(
    lipid = c("WDSD_0371", "WDSD_0300", "WDSD_0964"),
    text_x = c(-5.7, 3.8, -2.65), text_y = c(5.4, -7.7, 26)
  )
  label_indices <- match(label_positions$lipid, lipids$lipid)
  manual_labels <- identical(as.character(lipids$category[label_indices]),
      c("female_lower", "female_higher", "male_higher")) && counts[["male_lower"]] == 0L
  labels <- cbind(
    label_positions[rep(manual_labels, nrow(label_positions)), ],
    lipids[label_indices[rep(manual_labels, length(label_indices))], setdiff(names(lipids), "lipid")]
  )
  automatic_labels <- lipids[FALSE, ]
  if (!manual_labels) {
    highlighted <- subset(lipids, category != "other")
    strongest_first <- order(-pmax(abs(highlighted$female_t), abs(highlighted$male_t)), highlighted$lipid)
    highlighted <- highlighted[strongest_first, ]
    automatic_labels <- highlighted[!duplicated(highlighted$category), ]
  }
  lipids$example_labelled <- lipids$lipid %in% c(labels$lipid, automatic_labels$lipid)

  # Asinh(t/2) retains the full range while making the near-zero bands visible at print width.
  axis_transform <- scales::transform_pseudo_log(sigma = 1, base = exp(1))
  labels$lead_y <- axis_transform$inverse(
    axis_transform$transform(labels$text_y) - .27 * sign(labels$text_y - labels$male_t)
  )
  limit <- ceiling(max(abs(c(lipids$female_t, lipids$male_t))) / 10) * 10
  limit <- max(30, limit)

  plot <- ggplot2::ggplot(lipids, ggplot2::aes(female_t, male_t)) +
    ggplot2::annotate("rect", xmin = -limit, xmax = limit, ymin = -other_cutoff,
                      ymax = other_cutoff, fill = "#FCE8E8") +
    ggplot2::annotate("rect", xmin = -other_cutoff, xmax = other_cutoff, ymin = -limit,
                      ymax = limit, fill = "#DFEAFD", alpha = .78) +
    ggplot2::geom_hline(yintercept = 0, color = "#AAA1A4", linewidth = .23) +
    ggplot2::geom_vline(xintercept = 0, color = "#A0A8B2", linewidth = .23) +
    ggplot2::geom_point(data = subset(lipids, category == "other"),
                        color = "#9DA3AA", alpha = .48, size = .81) +
    ggplot2::geom_point(
      data = subset(lipids, category != "other"), ggplot2::aes(fill = category),
      shape = 21, color = "white", stroke = .24, size = 1.58, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = category_colors, limits = categories, labels = legend_labels, name = NULL,
      guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
                                    override.aes = list(size = 1.64, stroke = .15))
    ) +
    ggplot2::geom_segment(
      data = labels, ggplot2::aes(xend = text_x, yend = lead_y),
      color = "#57575A", linewidth = .29
    ) +
    ggplot2::geom_text(
      data = labels, ggplot2::aes(x = text_x, y = text_y, label = lipid_label),
      color = "#26272A", family = "Helvetica", size = 8.3, size.unit = "pt"
    ) +
    ggplot2::scale_x_continuous(transform = axis_transform, breaks = c(-30, -5, 0, 5, 30)) +
    ggplot2::scale_y_continuous(transform = axis_transform, breaks = c(-30, -5, 0, 5, 30)) +
    ggplot2::coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit), expand = FALSE) +
    ggplot2::labs(
      title = "B   Sex-polarised lipid response",
      x = expression("Female WD-SD (" * italic(t) * ")"),
      y = expression("Male WD-SD (" * italic(t) * ")"),
      caption = "Shading marks |t| < 1. Asinh axes."
    ) +
    cfdi_theme(9) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(3, 2, 2, 3),
      plot.title = ggplot2::element_text(size = 10.1, margin = ggplot2::margin(b = 4)),
      axis.text = ggplot2::element_text(size = 8.2),
      axis.title.x = ggplot2::element_text(size = 9, color = "#A43F44", margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(size = 9, color = "#305F9C", margin = ggplot2::margin(r = 3)),
      axis.ticks.length = grid::unit(2, "pt"),
      legend.position = "bottom", legend.location = "plot", legend.justification = "center",
      legend.text = ggplot2::element_text(size = 8.0),
      legend.key.width = grid::unit(9, "pt"), legend.key.height = grid::unit(8.5, "pt"),
      legend.spacing.x = grid::unit(2, "pt"), legend.spacing.y = grid::unit(1, "pt"),
      legend.margin = ggplot2::margin(t = 2, b = 1),
      plot.caption = ggplot2::element_text(size = 7.7, margin = ggplot2::margin(t = 2))
    )
  if (nrow(automatic_labels) > 0L) {
    plot <- plot + ggrepel::geom_text_repel(
      data = automatic_labels, ggplot2::aes(label = lipid_label),
      size = 2.9, family = "Helvetica", seed = 1602, max.overlaps = Inf
    )
  }
  list(plot = plot, source = lipids, stem = "Fig6B_diet_specificity",
       width_mm = 64, height_mm = 69, dpi = 600)
}

cf453_manuscript_figure6e <- function(lipids, network) {
  required <- c("row_type", "context", "sex", "zone", "gene", "lipid", "lipid_label",
    "lipid_class", "r_SD", "r_WD", "display_rank")
  stopifnot(all(required %in% names(network)))
  contexts <- c("female_PP", "female_PC", "male_PP", "male_PC")
  edges <- subset(network, row_type == "selected_edge")
  stopifnot(
    all(edges$context %in% contexts),
    all(edges$context == paste(edges$sex, edges$zone, sep = "_")),
    !anyDuplicated(edges[c("context", "lipid")]),
    !anyDuplicated(unique(edges[c("context", "gene")])$context),
    all(table(edges$context) <= 6L)
  )
  edges$context <- factor(edges$context, levels = contexts)
  edges <- edges[order(edges$context, edges$display_rank), ]
  effects <- lipids[match(edges$lipid, lipids$lipid), ]
  stopifnot(
    identical(edges$lipid, effects$lipid),
    identical(edges$lipid_label, effects$lipid_label),
    identical(edges$lipid_class, effects$lipid_class),
    all(is.finite(c(edges$r_SD, edges$r_WD))),
    all(abs(c(edges$r_SD, edges$r_WD)) <= 1)
  )
  log2_WD_over_SD <- ifelse(edges$sex == "female", effects$female_log2_WD_SD, effects$male_log2_WD_SD)
  plot_data <- data.frame(
    context = edges$context, gene = edges$gene, sex = edges$sex, zone = edges$zone,
    lipid = edges$lipid, lipid_label = edges$lipid_label, lipid_class = edges$lipid_class,
    log2_WD_over_SD = log2_WD_over_SD, r_SD = edges$r_SD, r_WD = edges$r_WD
  )
  positions <- plot_data[c("context", "log2_WD_over_SD")]
  points <- rbind(
    transform(positions, r = plot_data$r_SD, diet = rep("SD", nrow(plot_data)),
              point_fill = rep("white", nrow(plot_data))),
    transform(positions, r = plot_data$r_WD, diet = rep("WD", nrow(plot_data)),
              point_fill = unname(cfdr_colors[as.character(plot_data$context)]))
  )

  # Approved leaders move only the text. The SD and WD points retain their exact coordinates.
  label_placement <- utils::read.table(header = TRUE, text = "
  context    lipid     endpoint text_x text_y hjust
  female_PP  WDSD_1023 WD        -2.06   1.110  0.5
  female_PP  WDSD_0389 WD        -2.51  -0.270  0.5
  female_PP  WDSD_0802 SD        -1.40  -0.590  0.5
  female_PP  WDSD_0805 SD        -1.01  -0.180  0
  female_PP  WDSD_1080 WD         1.29   1.110  0.5
  female_PP  WDSD_0322 WD         0.67   0.750  0
  female_PC  WDSD_0404 WD        -1.02   1.075  1
  female_PC  WDSD_0909 SD        -2.01  -0.962  0.5
  female_PC  WDSD_0302 WD         0.02  -1.110  0.5
  female_PC  WDSD_0599 SD         0.64   0.540  0.5
  female_PC  WDSD_0688 SD         2.42  -1.110  0.5
  female_PC  WDSD_1114 SD         2.11   1.075  1
  male_PP    WDSD_0822 WD        -2.04   0.850  0
  male_PP    WDSD_0598 WD        -0.34  -0.335  1
  male_PP    WDSD_0151 WD         0.98   0.895  0
  male_PP    WDSD_0591 SD         1.12  -1.080  1
  male_PP    WDSD_0062 WD         2.41   1.090  1
  male_PP    WDSD_0463 SD         2.44  -1.087  0.5
  male_PC    WDSD_0179 SD        -1.90   0.420  0.5
  male_PC    WDSD_0175 WD        -2.20  -1.104  0.5
  male_PC    WDSD_0852 WD         0.10   0.510  0
  male_PC    WDSD_0382 WD         0.62  -1.104  0.5
  male_PC    WDSD_0997 SD         2.93   1.045  1
  male_PC    WDSD_0590 SD         2.95  -0.690  1
  ")

  label_key <- paste(label_placement$context, label_placement$lipid)
  plot_key <- paste(plot_data$context, plot_data$lipid)
  label_indices <- match(plot_key, label_key)
  max_abundance <- max(abs(plot_data$log2_WD_over_SD), 0)
  manual_labels <- setequal(plot_key, label_key) && max_abundance <= 3.4
  labels <- label_placement[label_indices[rep(manual_labels, length(label_indices))], ]
  labelled_data <- plot_data[rep(manual_labels, nrow(plot_data)), ]
  labels$context <- labelled_data$context
  labels$lipid_label <- labelled_data$lipid_label
  labels$lead_x <- labelled_data$log2_WD_over_SD
  labels$lead_y <- ifelse(labels$endpoint == "SD", labelled_data$r_SD, labelled_data$r_WD)
  labels$lead_end_x <- labels$text_x + 0.08 * sign(labels$hjust - 0.5)
  labels$lead_end_y <- labels$text_y
  centered <- labels$hjust == 0.5
  labels$lead_end_y[centered] <- labels$text_y[centered] +
    0.063 * sign(labels$lead_y[centered] - labels$text_y[centered])

  representatives <- unique(plot_data[c("context", "gene")])
  genes <- representatives$gene[match(contexts, representatives$context)]
  strip_titles <- stats::setNames(ifelse(
    is.na(genes), encodeString(cfdr_context_labels[contexts], quote = '"'),
    sprintf("bolditalic(%s)~~%s", encodeString(genes, quote = '"'),
            encodeString(cfdr_context_labels[contexts], quote = '"'))
  ), contexts)
  outline_colors <- cfdr_colors[contexts]
  outline_colors[c("female_PP", "male_PP")] <- c("#B86269", "#5480BD")
  x_limit <- if (max_abundance <= 3.4) 3.65 else ceiling(max_abundance * 1.10)
  x_breaks <- if (x_limit == 3.65) -3:3 else pretty(c(-x_limit, x_limit), n = 7)
  caption <- if (all(table(plot_data$context) == 6L)) {
    "Six selected pairs per gene. The SD and WD groups each contain seven mice per sex."
  } else {
    "Selected pairs shown. The SD and WD groups each contain seven mice per sex."
  }

  plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = "#AEAEAE", linetype = "dashed", linewidth = 0.23) +
    ggplot2::annotate(
      "segment", x = 0, xend = 0, y = -1, yend = 1,
      color = "#AEAEAE", linetype = "dashed", linewidth = 0.23
    ) +
    ggplot2::geom_segment(
      data = labels, ggplot2::aes(x = lead_x, y = lead_y, xend = lead_end_x, yend = lead_end_y),
      color = "#777777", linewidth = 0.25, lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(x = log2_WD_over_SD, xend = log2_WD_over_SD, y = r_SD,
                   yend = r_WD - 0.031 * sign(r_WD - r_SD), color = context),
      linewidth = 0.48, lineend = "round",
      arrow = grid::arrow(angle = 24, length = grid::unit(1.42, "mm"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = points, ggplot2::aes(x = log2_WD_over_SD, y = r, shape = diet,
                                  color = context, fill = point_fill),
      size = 1.62, stroke = 0.45
    ) +
    ggplot2::geom_text(
      data = labels, ggplot2::aes(x = text_x, y = text_y, label = lipid_label, hjust = hjust),
      family = "Helvetica", color = "#242424", size = 8.8, size.unit = "pt"
    ) +
    ggplot2::scale_color_manual(values = outline_colors, guide = "none") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_manual(
      values = c(SD = 21, WD = 21), limits = c("SD", "WD"), breaks = c("SD", "WD"),
      labels = c("Standard diet (SD)", "Western diet (WD)"), name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(
        color = "#3D454B", fill = c("white", "#3D454B"), size = 1.85
      ))
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
    ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(-1.19, 1.19), expand = FALSE) +
    ggplot2::facet_wrap(
      ggplot2::vars(context), ncol = 2, axes = "all", axis.labels = "all", drop = FALSE,
      labeller = ggplot2::as_labeller(strip_titles, default = ggplot2::label_parsed)
    ) +
    ggplot2::labs(
      title = expression(bold("E") ~~ bold("Gene-lipid correlations by diet")),
      x = expression("Whole-liver lipid abundance, " * log[2] * "(WD/SD)"),
      y = expression("Pearson correlation (" * italic(r) * ")"),
      caption = caption
    ) +
    cfdr_theme(10) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(4, 5, 4, 4),
      panel.spacing.x = grid::unit(9, "pt"), panel.spacing.y = grid::unit(5, "pt"),
      plot.title = ggplot2::element_text(size = 11.6, margin = ggplot2::margin(b = 2)),
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10.5, hjust = 0, margin = ggplot2::margin(b = 3)),
      axis.text = ggplot2::element_text(size = 8.6),
      axis.title.x = ggplot2::element_text(size = 10.0, margin = ggplot2::margin(t = 5)),
      axis.title.y = ggplot2::element_text(size = 10.0, margin = ggplot2::margin(r = 4)),
      legend.position = "top", legend.justification = "right",
      legend.text = ggplot2::element_text(size = 9.0),
      legend.key.width = grid::unit(3.9, "mm"), legend.key.height = grid::unit(3.5, "mm"),
      legend.spacing.x = grid::unit(1.4, "mm"), legend.margin = ggplot2::margin(b = 3),
      plot.caption = ggplot2::element_text(size = 8.4, margin = ggplot2::margin(t = 5))
    )
  if (!manual_labels && nrow(plot_data) > 0L) {
    plot <- plot + ggrepel::geom_text_repel(
      data = plot_data, ggplot2::aes(x = log2_WD_over_SD, y = r_WD, label = lipid_label),
      size = 2.9, family = "Helvetica", seed = 1605, max.overlaps = Inf
    )
  }
  empty_contexts <- data.frame(context = factor(contexts[is.na(genes)], levels = contexts))
  if (nrow(empty_contexts) > 0L) {
    plot <- plot + ggplot2::geom_text(
      data = empty_contexts, x = 0, y = 0.35, label = "No qualifying pairs",
      size = 9, size.unit = "pt", family = "Helvetica", color = "#454545"
    )
  }
  list(plot = plot, source = plot_data, stem = "Fig6E_four_examples",
       width_mm = 180, height_mm = 190, dpi = 450)
}
