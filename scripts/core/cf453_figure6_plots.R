# Publication panels for the CF453 standalone analysis.
# Sourcing this module defines functions only; analysis is performed elsewhere.

cf6_plot_colors <- c(
  female_pericentral = "#B83338", female_periportal = "#ED8585",
  male_pericentral = "#285CA8", male_periportal = "#80A8E0",
  shared = "#929292", no_association = "#F3F3F3", other = "#D5D5D5"
)

cf6_plot_category_labels <- c(
  female_pericentral = "Female PC", female_periportal = "Female PP",
  male_pericentral = "Male PC", male_periportal = "Male PP",
  shared = "Shared", no_association = "No association", other = "Other genes"
)

cf6_plot_context <- function(sex, zone) {
  zone <- tolower(as.character(zone))
  zone[zone %in% c("pp", "portal")] <- "periportal"
  zone[zone %in% c("pc", "central")] <- "pericentral"
  paste(tolower(as.character(sex)), zone, sep = "_")
}

cf6_plot_category <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("_pp$", "_periportal", x)
  x <- gsub("_pc$", "_pericentral", x)
  x[x %in% c("unselected", "none")] <- "other"
  x[x %in% c("unassigned", "no_network_data")] <- "no_association"
  x
}

cf6_plot_write_source <- function(x, output, stem) {
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  utils::write.csv(x, file.path(output, paste0(stem, ".csv")), row.names = FALSE, na = "")
}

cf6_plot_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, color = "#4A4A4A", margin = ggplot2::margin(b = 9)),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.title = ggplot2::element_text(color = "#252525"),
      axis.text = ggplot2::element_text(color = "#303030"),
      axis.line = ggplot2::element_line(linewidth = .35, color = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = .3),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.height = grid::unit(4.2, "mm"),
      strip.background = ggplot2::element_rect(fill = "#F4F4F4", color = NA),
      strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(6, 5, 6, 5)),
      panel.spacing = grid::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(10, 13, 10, 10)
    )
}

cf6_plot_stat <- function(x, digits = 2) {
  out <- rep("NA", length(x))
  finite <- is.finite(x)
  out[finite] <- formatC(x[finite], digits = digits, format = "f")
  out
}

cf6_plot_save <- function(plot, output, stem, width, height, family,
                          diet = NA_character_, context = NA_character_, source_table = NA_character_) {
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  stem <- gsub("[^A-Za-z0-9_.-]", "_", stem)
  files <- file.path(output, paste0(stem, c(".pdf", ".svg", ".png")))
  ggplot2::ggsave(files[1], plot = plot, width = width, height = height, units = "in",
                 device = grDevices::pdf, family = "Helvetica", bg = "white", useDingbats = FALSE)
  ggplot2::ggsave(files[2], plot = plot, width = width, height = height, units = "in",
                 device = svglite::svglite, bg = "white")
  ggplot2::ggsave(files[3], plot = plot, width = width, height = height, units = "in",
                 device = ragg::agg_png, dpi = 300, bg = "white")
  data.frame(panel = stem, family = family, diet = diet, context = context,
             file = normalizePath(files, mustWork = TRUE), format = c("PDF", "SVG", "PNG"),
             width_in = width, height_in = height, source_table = source_table,
             stringsAsFactors = FALSE)
}

cf6_plot_schematic <- function(metadata) {
  if (!all(c("sex", "diet") %in% names(metadata))) stop("Schematic metadata requires sex and diet.")
  n_sex <- function(diet, sex) sum(tolower(metadata$sex) == sex & toupper(metadata$diet) == diet)
  diet_boxes <- data.frame(
    x = .17, y = c(.70, .28), diet = c("SD", "WD"),
    label = c("Standard diet", "Western diet")
  )
  sex_boxes <- expand.grid(sex = c("female", "male"), diet = c("SD", "WD"), stringsAsFactors = FALSE)
  sex_boxes$x <- rep(c(.35, .53), 2)
  sex_boxes$y <- rep(c(.70, .28), each = 2)
  sex_boxes$label <- paste0(tools::toTitleCase(sex_boxes$sex), "\nn = ",
                            mapply(n_sex, sex_boxes$diet, sex_boxes$sex))
  sex_boxes$color <- ifelse(sex_boxes$sex == "female", "#CC4549", "#3C74BC")
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = .03, y = .96, label = "CF453", hjust = 0, fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = .03, y = .89,
                      label = paste0(nrow(metadata), " mice with paired transcriptomics and lipidomics"),
                      hjust = 0, color = "#4A4A4A", size = 3.8) +
    ggplot2::geom_label(data = diet_boxes, ggplot2::aes(x, y, label = label),
                        linewidth = .35, label.padding = grid::unit(4, "mm"),
                        fill = "#F5F5F5", size = 3.6, lineheight = 1.1) +
    ggplot2::geom_label(data = sex_boxes, ggplot2::aes(x, y, label = label, color = color),
                        linewidth = .45, label.padding = grid::unit(4, "mm"),
                        fill = "white", size = 3.7, lineheight = 1.25) +
    ggplot2::scale_color_identity() +
    ggplot2::annotate("segment", x = .64, xend = .72, y = .70, yend = .70,
                      linewidth = .45, arrow = grid::arrow(length = grid::unit(2, "mm"))) +
    ggplot2::annotate("segment", x = .64, xend = .72, y = .28, yend = .28,
                      linewidth = .45, arrow = grid::arrow(length = grid::unit(2, "mm"))) +
    ggplot2::annotate("label", x = .855, y = .70, label = "Periportal RNA\nPericentral RNA\nWhole-liver lipids",
                      fill = "white", color = "#303030", linewidth = .35,
                      label.padding = grid::unit(4, "mm"), size = 3.7, lineheight = 1.3) +
    ggplot2::annotate("label", x = .855, y = .28, label = "Periportal RNA\nPericentral RNA\nWhole-liver lipids",
                      fill = "white", color = "#303030", linewidth = .35,
                      label.padding = grid::unit(4, "mm"), size = 3.7, lineheight = 1.3) +
    ggplot2::annotate("text", x = .03, y = .025,
                      label = "Gene-lipid associations are estimated within each sex and diet.",
                      hjust = 0, size = 3.7, color = "#404040") +
    ggplot2::coord_cartesian(xlim = c(0, 1.04), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 15, 12, 10), legend.position = "none")
}

cf6_plot_zonation <- function(data, diet, representatives = NULL, limits = NULL) {
  d <- data[toupper(data$diet) == diet, , drop = FALSE]
  d$plot_category <- cf6_plot_category(d$category)
  d$plot_category[!d$plot_category %in% names(cf6_plot_colors)] <- "other"
  d$plot_category <- factor(d$plot_category, levels = c("other", "shared", "female_pericentral", "female_periportal", "male_pericentral", "male_periportal"))
  d <- d[order(d$plot_category), , drop = FALSE]
  d$label <- FALSE
  selected <- which(d$plot_category != "other" & is.finite(d$q_interaction))
  if (length(selected)) {
    selected <- selected[order(d$q_interaction[selected], d$gene[selected])]
    d$label[head(selected, 8)] <- TRUE
  }
  if (!is.null(representatives) && nrow(representatives)) d$label[d$gene %in% representatives$gene] <- TRUE
  if (is.null(limits)) {
    extent <- max(abs(c(data$t_female, data$t_male)), na.rm = TRUE)
    limits <- c(-1, 1) * max(1, extent * 1.08)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(t_female, t_male)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#C5C5C5", linewidth = .4) +
    ggplot2::geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = .3) +
    ggplot2::geom_point(ggplot2::aes(color = plot_category), alpha = .65, size = .85) +
    ggplot2::scale_color_manual(values = cf6_plot_colors,
                               labels = cf6_plot_category_labels[c("female_pericentral", "female_periportal", "male_pericentral", "male_periportal", "shared")],
                               name = "Zonation program", drop = FALSE,
                               breaks = c("female_pericentral", "female_periportal", "male_pericentral", "male_periportal", "shared")) +
    ggplot2::coord_equal(xlim = limits, ylim = limits, clip = "off") +
    ggplot2::labs(title = paste0("Hepatocyte zonation | ", diet),
                  subtitle = "Positive = periportal; negative = pericentral",
                  x = "Female zonation effect (t-statistic)", y = "Male zonation effect (t-statistic)") +
    cf6_plot_theme() + ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.3, alpha = 1)))
  if (any(d$label)) p <- p + ggrepel::geom_label_repel(
    data = d[d$label, , drop = FALSE], ggplot2::aes(label = gene, color = plot_category),
    seed = 1406, size = 3.1, max.overlaps = Inf, min.segment.length = 0,
    segment.size = .25, segment.color = "#777777", box.padding = .6,
    point.padding = .6, force = 2, fill = "white", label.size = .1,
    label.padding = grid::unit(.10, "lines"), label.r = grid::unit(.08, "lines"),
    show.legend = FALSE, family = "Helvetica")
  list(plot = p, source = d)
}

cf6_plot_class_membership <- function(data, diet, classes = NULL) {
  d <- data[toupper(data$diet) == diet, , drop = FALSE]
  d$plot_category <- cf6_plot_category(d$category)
  lev <- c("female_pericentral", "female_periportal", "male_pericentral", "male_periportal", "shared", "no_association")
  if (any(!d$plot_category %in% lev)) stop("Unknown lipid membership category.")
  if (is.null(classes)) classes <- sort(unique(as.character(data$lipid_class)))
  d$lipid_class <- factor(d$lipid_class, levels = rev(classes))
  d$plot_category <- factor(d$plot_category, levels = lev)
  sums <- stats::aggregate(proportion ~ lipid_class, data = d, FUN = sum)
  if (any(abs(sums$proportion - 1) > 1e-8)) stop("Within-class membership proportions must sum to one.")
  denominators <- unique(d[, c("lipid_class", "total")])
  class_labels <- stats::setNames(paste0(denominators$lipid_class, "  (n = ", denominators$total, ")"), denominators$lipid_class)
  p <- ggplot2::ggplot(d, ggplot2::aes(lipid_class, proportion, fill = plot_category)) +
    ggplot2::geom_col(width = .78, position = ggplot2::position_stack(reverse = TRUE), color = "white", linewidth = .2) +
    ggplot2::scale_fill_manual(values = cf6_plot_colors[lev], labels = cf6_plot_category_labels[lev],
                              name = "Associated gene program", drop = FALSE) +
    ggplot2::scale_x_discrete(labels = class_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, .01))) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste0("Lipid associations by class | ", diet),
                  subtitle = "Shared = associations with more than one sex/zone context",
                  x = NULL, y = "Proportion of eligible lipid species") +
    cf6_plot_theme() + ggplot2::theme(legend.position = "bottom", axis.ticks.y = ggplot2::element_blank(),
                                    axis.line.y = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3, byrow = TRUE))
  list(plot = p, source = d)
}

cf6_plot_pair_density <- function(data, diet, limits = NULL) {
  d <- data[toupper(data$diet) == diet, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$context <- factor(cf6_plot_context(d$sex, d$zone),
                       levels = c("female_pericentral", "female_periportal", "male_pericentral", "male_periportal"))
  d$lipid_class <- factor(d$lipid_class, levels = rev(sort(unique(as.character(data$lipid_class)))))
  if (is.null(limits)) limits <- c(0, max(data$display_fraction, na.rm = TRUE))
  limits[2] <- max(limits[2], .01)
  d$label <- scales::label_percent(accuracy = .1)(d$display_fraction)
  d$label_color <- ifelse(d$display_fraction > .60 * limits[2], "white", "#303030")
  p <- ggplot2::ggplot(d, ggplot2::aes(context, lipid_class, fill = display_fraction)) +
    ggplot2::geom_tile(color = "white", linewidth = .45, width = .96, height = .96) +
    ggplot2::geom_text(ggplot2::aes(label = label, color = label_color), size = 2.7, family = "Helvetica") +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(low = "#F5F1F6", high = "#664071", limits = limits,
                                 labels = scales::label_percent(accuracy = 1),
                                 name = "Fraction of pairs\n|r| >= 0.65",
                                 guide = ggplot2::guide_colourbar(display = "rectangles")) +
    ggplot2::scale_x_discrete(labels = cf6_plot_category_labels, drop = FALSE) +
    ggplot2::labs(title = paste0("Gene-lipid association density | ", diet),
                  subtitle = "Exploratory associations in a fixed candidate universe",
                  x = NULL, y = NULL) +
    cf6_plot_theme(10) + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                      legend.position = "right", panel.grid = ggplot2::element_blank(),
                                      axis.text.x = ggplot2::element_text(size = 8.7))
  list(plot = p, source = d)
}

cf6_plot_empty <- function(title, message) {
  ggplot2::ggplot() + ggplot2::annotate("text", x = .5, y = .5, label = message,
                                       size = 4, lineheight = 1.3, color = "#555555", family = "Helvetica") +
    ggplot2::labs(title = title) + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    cf6_plot_theme() + ggplot2::theme(axis.line = ggplot2::element_blank(),
                                     axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                     axis.title = ggplot2::element_blank())
}

cf6_plot_network <- function(edges, sex, zone, diet, representative = NULL) {
  context <- cf6_plot_context(sex, zone)
  z_label <- ifelse(grepl("periportal", context), "periportal", "pericentral")
  title <- paste0(tools::toTitleCase(sex), " ", z_label, " | ", diet)
  keep <- cf6_plot_context(edges$sex, edges$zone) == context & toupper(edges$diet) == diet
  d <- edges[keep, , drop = FALSE]
  if (!is.null(representative) && length(representative)) d <- d[d$gene == representative[1], , drop = FALSE]
  if (nrow(d) == 0L) return(list(plot = cf6_plot_empty(title, "No eligible representative\ngene-lipid associations"), source = d))
  d <- d[order(d$partner_rank, d$lipid), , drop = FALSE]
  d <- d[!duplicated(d$lipid), , drop = FALSE]
  if (nrow(d) > 6L) stop("Representative network has more than six distinct partners.")
  d$y <- seq(.84, .16, length.out = nrow(d))
  if (nrow(d) == 1L) d$y <- .5
  d$line_type <- ifelse(d$r >= 0, "Positive", "Negative")
  d$plot_lipid <- if ("lipid_label" %in% names(d)) as.character(d$lipid_label) else as.character(d$lipid)
  d$edge_label <- paste0("r = ", cf6_plot_stat(d$r))
  color <- unname(cf6_plot_colors[context])
  if (is.na(color)) color <- "#555555"
  gene <- unique(d$gene)
  if (length(gene) != 1L) stop("A representative network must contain exactly one gene.")
  nomination <- if ("active_candidate" %in% names(d) && !all(d$active_candidate))
    "RNA nomination in the other diet" else "RNA-selected in this diet"
  d$edge_opacity <- ifelse(abs(d$r) >= .65, .85, .3)
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_segment(ggplot2::aes(x = .16, xend = .68, y = .5, yend = y,
                                     linewidth = abs(r), linetype = line_type, alpha = edge_opacity),
                          color = color) +
    ggplot2::geom_label(ggplot2::aes(x = .69, y = y, label = plot_lipid), hjust = 0,
                        fill = "white", color = "#202020", linewidth = .3,
                        label.padding = grid::unit(2.4, "mm"), label.r = grid::unit(.7, "mm"),
                        family = "Helvetica", size = 3.35) +
    ggplot2::geom_text(ggplot2::aes(x = .52, y = .5 + (y - .5) * .73 + .035, label = edge_label),
                       size = 3.1, color = "#303030", family = "Helvetica") +
    ggplot2::annotate("label", x = .075, y = .5, label = gene,
                      fill = color, color = "white", fontface = "italic", size = 4.5,
                      linewidth = 0, label.padding = grid::unit(4, "mm"),
                      label.r = grid::unit(5, "mm"), family = "Helvetica") +
    ggplot2::scale_linewidth_continuous(range = c(.65, 2.2), limits = c(0, 1), guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_linetype_manual(values = c(Positive = "solid", Negative = "22"), name = "Association") +
    ggplot2::coord_cartesian(xlim = c(-.08, 1.45), ylim = c(.02, .98), clip = "off") +
    ggplot2::labs(title = title, subtitle = paste0(gene, " | n = ", unique(d$n)[1], " mice | ", nomination),
                  caption = "Exploratory Pearson associations. Same partners shown in both diets; faint edges: |r| < 0.65.") +
    cf6_plot_theme() + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                                     axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                                     legend.position = "bottom", legend.key.width = grid::unit(10, "mm"),
                                     plot.caption = ggplot2::element_text(size = 8.4, color = "#555555", hjust = 0))
  list(plot = p, source = d)
}

cf6_plot_diet_zonation <- function(data, representatives = NULL, primary_interactions = NULL) {
  d <- data
  if (nrow(d) == 0L) return(list(plot = cf6_plot_empty("Diet modification of zonation", "No estimable zonation contrasts"), source = d))
  d$sex <- tolower(d$sex)
  genes <- if (!is.null(representatives) && nrow(representatives)) unique(as.character(representatives$gene)) else character()
  primary_q <- if (!is.null(primary_interactions)) primary_interactions[primary_interactions$contrast == "sex_diet_zone", , drop = FALSE] else NULL
  rank <- if (!is.null(primary_q) && nrow(primary_q)) primary_q[, c("gene", "q")] else stats::aggregate(q ~ gene, data = d, FUN = min)
  rank <- rank[order(rank$q, rank$gene), , drop = FALSE]
  genes <- unique(c(genes, head(rank$gene, 16)))
  genes <- genes[genes %in% d$gene]
  d <- d[d$gene %in% genes, , drop = FALSE]
  d$gene <- factor(d$gene, levels = rev(genes))
  d$sex <- factor(d$sex, levels = c("female", "male"))
  pd <- ggplot2::position_dodge(width = .55)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = effect, y = gene, color = sex)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = .4, color = "#AAAAAA") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = low, xmax = high), position = pd, linewidth = .6) +
    ggplot2::geom_point(ggplot2::aes(shape = q <= .05), position = pd, size = 2.5, stroke = .8) +
    ggplot2::scale_color_manual(values = c(female = "#CC4549", male = "#3C74BC"),
                               labels = c(female = "Female", male = "Male"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                               labels = c(`TRUE` = "q <= 0.05", `FALSE` = "q > 0.05"), name = "Within-sex diet effect") +
    ggplot2::labs(title = "Diet modification of hepatocyte zonation",
                  subtitle = "WD minus SD change in the periportal-minus-pericentral effect",
                  x = "Change in zonation effect (log-normalized expression)", y = NULL) +
    cf6_plot_theme() + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(face = "italic"))
  if (!is.null(primary_interactions) && nrow(primary_interactions)) {
    q <- primary_interactions[primary_interactions$contrast == "sex_diet_zone", , drop = FALSE]
    if (nrow(q) && !anyDuplicated(q$gene)) {
      q <- q[match(levels(d$gene), q$gene), , drop = FALSE]
      q$gene <- factor(q$gene, levels = levels(d$gene))
      q <- q[!is.na(q$gene), , drop = FALSE]
      q$label <- paste0("3-way q = ", formatC(q$q, digits = 2, format = "g"))
      effect_range <- range(c(d$low, d$high), finite = TRUE)
      span <- max(diff(effect_range), .1)
      q$x <- max(effect_range) + .15 * span
      p <- p + ggplot2::geom_text(data = q, ggplot2::aes(x = x, y = gene, label = label),
                                  inherit.aes = FALSE, hjust = 0, size = 2.7, color = "#4A4A4A", family = "Helvetica") +
        ggplot2::expand_limits(x = max(effect_range) + .56 * span)
      d$q_three_way <- primary_q$q[match(d$gene, primary_q$gene)]
    }
  }
  list(plot = p, source = d)
}

cf6_plot_examples <- function(data, sex, zone) {
  context <- cf6_plot_context(sex, zone)
  d <- data[cf6_plot_context(data$sex, data$zone) == context, , drop = FALSE]
  title <- paste0(tools::toTitleCase(sex), " ", ifelse(grepl("periportal", context), "periportal", "pericentral"), " associations")
  if (!nrow(d)) return(list(plot = cf6_plot_empty(title, "No eligible representative\ngene-lipid examples"), source = d))
  d$diet <- factor(toupper(d$diet), levels = c("SD", "WD"))
  d$example <- paste0(d$gene, " : ", d$lipid_label)
  summaries <- lapply(split(d, interaction(d$example, d$diet, drop = TRUE)), function(z) {
    z <- z[is.finite(z$rna) & is.finite(z$lipid_log2), , drop = FALSE]
    omitted <- vapply(seq_len(nrow(z)), function(j) stats::cor(z$rna[-j], z$lipid_log2[-j]), numeric(1L))
    data.frame(example = z$example[1], diet = as.character(z$diet[1]), n_plot = nrow(z),
               r_plot = stats::cor(z$rna, z$lipid_log2),
               loo_min_plot = min(omitted, na.rm = TRUE), loo_max_plot = max(omitted, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  summary <- do.call(rbind, summaries)
  summary$label <- paste0(summary$diet, ": r = ", cf6_plot_stat(summary$r_plot), " (n = ", summary$n_plot,
                           "); leave-one-mouse-out r: ", cf6_plot_stat(summary$loo_min_plot), " to ", cf6_plot_stat(summary$loo_max_plot))
  key <- paste(d$example, d$diet)
  si <- match(key, paste(summary$example, summary$diet))
  for (nm in c("n_plot", "r_plot", "loo_min_plot", "loo_max_plot")) d[[nm]] <- summary[[nm]][si]
  caption <- if (length(unique(d$example)) == 1L) paste(summary$label, collapse = "\n") else
    "Lines are within-diet linear fits; correlation and omission ranges are supplied with source data."
  color <- unname(cf6_plot_colors[context])
  p <- ggplot2::ggplot(d, ggplot2::aes(rna, lipid_log2, shape = diet, linetype = diet)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .6, color = color) +
    ggplot2::geom_point(color = color, size = 2.2, alpha = .9) +
    ggplot2::facet_wrap(~example, ncol = 2, scales = "free") +
    ggplot2::scale_shape_manual(values = c(SD = 16, WD = 17), name = "Diet") +
    ggplot2::scale_linetype_manual(values = c(SD = "solid", WD = "22"), name = "Diet") +
    ggplot2::labs(title = title, subtitle = "Exploratory selected examples | each point represents one mouse",
                  x = "Log-normalized RNA expression", y = "Whole-liver lipid abundance (log2)",
                  caption = caption) +
    cf6_plot_theme(10) + ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(size = 8.5),
                                      plot.caption = ggplot2::element_text(size = 8.5, color = "#555555", hjust = 0, lineheight = 1.25))
  list(plot = p, source = d)
}

cf6_plot_pca <- function(data, kind) {
  d <- data[tolower(data$kind) == tolower(kind), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$sex <- factor(tolower(d$sex), levels = c("female", "male"))
  d$diet <- factor(toupper(d$diet), levels = c("SD", "WD"))
  d$zone <- as.character(d$zone)
  is_rna <- tolower(kind) == "rna"
  if (is_rna) {
    d$context <- cf6_plot_context(d$sex, d$zone)
    p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, color = context, shape = diet)) +
      ggplot2::geom_point(size = 2.6, alpha = .85) +
      ggplot2::scale_color_manual(values = cf6_plot_colors, labels = cf6_plot_category_labels, name = "Sex and zone")
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, color = sex, shape = diet)) +
      ggplot2::geom_point(size = 2.6, alpha = .85) +
      ggplot2::scale_color_manual(values = c(female = "#CC4549", male = "#3C74BC"),
                                 labels = c(female = "Female", male = "Male"), name = "Sex")
  }
  v1 <- unique(d$variance_PC1[is.finite(d$variance_PC1)])[1]
  v2 <- unique(d$variance_PC2[is.finite(d$variance_PC2)])[1]
  if (v1 <= 1) v1 <- 100 * v1
  if (v2 <= 1) v2 <- 100 * v2
  p <- p + ggplot2::scale_shape_manual(values = c(SD = 16, WD = 17), name = "Diet") +
    ggplot2::labs(title = paste0(ifelse(is_rna, "Zonated hepatocyte RNA", "Whole-liver lipids"), " | PCA"),
                  x = paste0("PC1 (", cf6_plot_stat(v1, 1), "%)"),
                  y = paste0("PC2 (", cf6_plot_stat(v2, 1), "%)")) +
    cf6_plot_theme() + ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = ifelse(is_rna, 2, 1), byrow = TRUE),
                    shape = ggplot2::guide_legend(nrow = 1))
  list(plot = p, source = d)
}

cf6_plot_lipid_effect <- function(data, contrast) {
  d <- data[data$contrast == contrast, , drop = FALSE]
  d$significance <- ifelse(d$q <= .05, "q <= 0.05", "q > 0.05")
  d$neglog10q <- -log10(pmax(d$q, .Machine$double.xmin))
  d$plot_lipid <- if ("lipid_label" %in% names(d)) d$lipid_label else d$lipid
  titles <- c(sex_SD = "Lipid sex difference | SD", sex_WD = "Lipid sex difference | WD",
               diet_female = "Lipid diet effect | Female", diet_male = "Lipid diet effect | Male",
               sex_diet = "Sex difference in lipid diet response")
  axes <- c(sex_SD = "Male-minus-female lipid effect (log2)", sex_WD = "Male-minus-female lipid effect (log2)",
             diet_female = "WD-minus-SD lipid effect (log2)", diet_male = "WD-minus-SD lipid effect (log2)",
             sex_diet = "Male (WD-SD) minus female (WD-SD) lipid effect (log2)")
  title <- if (contrast %in% names(titles)) titles[[contrast]] else paste0("Lipid abundance | ", gsub("_", " ", contrast))
  axis <- if (contrast %in% names(axes)) axes[[contrast]] else "Lipid abundance effect (log2)"
  labels <- head(order(d$q, -abs(d$effect), d$lipid), 8)
  p <- ggplot2::ggplot(d, ggplot2::aes(effect, neglog10q)) +
    ggplot2::geom_hline(yintercept = -log10(.05), color = "#AAAAAA", linetype = "22", linewidth = .4) +
    ggplot2::geom_vline(xintercept = 0, color = "#AAAAAA", linewidth = .3) +
    ggplot2::geom_point(ggplot2::aes(color = significance), alpha = .75, size = 1.4) +
    ggplot2::scale_color_manual(values = c(`q <= 0.05` = "#75607E", `q > 0.05` = "#BEBEBE"), name = NULL) +
    ggplot2::labs(title = title, x = axis, y = "-log10 adjusted P value") +
    cf6_plot_theme() + ggplot2::theme(legend.position = "bottom")
  if (length(labels)) p <- p + ggrepel::geom_text_repel(
    data = d[labels, , drop = FALSE], ggplot2::aes(label = plot_lipid),
    seed = 1406, size = 2.8, max.overlaps = Inf, min.segment.length = 0,
    segment.size = .25, segment.color = "#888888", family = "Helvetica")
  list(plot = p, source = d)
}

cf6_plot_thresholds <- function(data) {
  d <- stats::aggregate(cbind(assigned, total) ~ diet + threshold, data = data, FUN = sum)
  d$proportion <- d$assigned / d$total
  d$diet <- factor(toupper(d$diet), levels = c("SD", "WD"))
  p <- ggplot2::ggplot(d, ggplot2::aes(threshold, proportion, linetype = diet)) +
    ggplot2::geom_line(color = "#5C5369", linewidth = .85) +
    ggplot2::geom_point(color = "#5C5369", size = 1.8) +
    ggplot2::geom_vline(xintercept = .65, color = "#AAAAAA", linetype = "22", linewidth = .45) +
    ggplot2::scale_linetype_manual(values = c(SD = "solid", WD = "22"), name = "Diet") +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(title = "Association display threshold sensitivity",
                  x = "Minimum absolute Pearson correlation", y = "Eligible lipids with at least one association") +
    cf6_plot_theme() + ggplot2::theme(legend.position = "bottom")
  list(plot = p, source = d)
}

cf6_plot_quality <- function(data) {
  d <- stats::aggregate(cbind(assigned, total) ~ diet + quality, data = data, FUN = sum)
  d$proportion <- d$assigned / d$total
  d$quality <- factor(d$quality, levels = c("keep_broad", "keep_primary", "keep_fill"),
                       labels = c("All eligible", "QC CV <= 30%", "QC + Fill >= 50%"))
  d$diet <- factor(toupper(d$diet), levels = c("SD", "WD"))
  d$label <- paste0(d$assigned, " / ", d$total)
  p <- ggplot2::ggplot(d, ggplot2::aes(quality, proportion)) +
    ggplot2::geom_col(fill = "#88768F", width = .64) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -.55, size = 3.2, family = "Helvetica") +
    ggplot2::facet_wrap(~diet, nrow = 1) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::labs(title = "Lipid quality sensitivity", subtitle = "Labels show associated / eligible lipid species",
                  x = NULL, y = "Lipids with at least one exploratory association") +
    cf6_plot_theme(10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  list(plot = p, source = d)
}

cf6_plot_markers <- function(data) {
  markers <- c("Cyp2e1", "Cyp2f2", "Glul", "Gls2", "Hal", "Sds")
  d <- data[data$gene %in% markers, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$gene <- factor(d$gene, levels = markers)
  d$diet <- factor(toupper(d$diet), levels = c("SD", "WD"))
  d$sex <- factor(tolower(d$sex), levels = c("female", "male"))
  d$zone <- factor(ifelse(tolower(d$zonation_type) == "periportal", "PP", "PC"), levels = c("PC", "PP"))
  p <- ggplot2::ggplot(d, ggplot2::aes(zone, supplied_data_expression, color = sex, fill = sex)) +
    ggplot2::geom_boxplot(width = .65, alpha = .1, outlier.shape = NA, linewidth = .35,
                          position = ggplot2::position_dodge(width = .7)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .08, dodge.width = .7, seed = 1406),
                        size = 1.4, alpha = .8) +
    ggplot2::facet_grid(gene ~ diet, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(female = "#CC4549", male = "#3C74BC"),
                               labels = c(female = "Female", male = "Male"), name = NULL) +
    ggplot2::scale_fill_manual(values = c(female = "#CC4549", male = "#3C74BC"), guide = "none") +
    ggplot2::labs(title = "Hepatocyte zonation markers", x = "Hepatocyte zone",
                  y = "Log-normalized RNA expression") +
    cf6_plot_theme(10) + ggplot2::theme(legend.position = "bottom", strip.text.y = ggplot2::element_text(face = "italic"))
  list(plot = p, source = d)
}

cf6_make_panels <- function(analysis, output) {
  required <- c("metadata", "zonation", "class_summary", "representatives", "representative_edges",
                "examples", "zonation_diet_effects", "lipid_results", "qc_pca", "thresholds")
  absent <- setdiff(required, names(analysis))
  if (length(absent)) stop("Analysis object lacks plotting inputs: ", paste(absent, collapse = ", "))
  panel_dir <- file.path(output, "panels")
  source_dir <- file.path(output, "source_data")
  if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
  if (!dir.exists(source_dir)) dir.create(source_dir, recursive = TRUE)
  manifest <- list()
  add <- function(result, stem, width, height, family, diet = NA_character_, context = NA_character_) {
    if (is.null(result)) return(invisible(NULL))
    source_file <- file.path(source_dir, paste0(stem, ".csv"))
    cf6_plot_write_source(result$source, source_dir, stem)
    manifest[[length(manifest) + 1L]] <<- cf6_plot_save(
      result$plot, panel_dir, stem, width, height, family, diet, context,
      normalizePath(source_file, mustWork = TRUE))
    invisible(NULL)
  }
  for (diet in c("SD", "WD")) {
    add(cf6_plot_zonation(analysis$zonation, diet, analysis$representatives),
        paste0("Fig6A_", diet), 6.3, 6.5, "A", diet)
  }
  add(list(plot = cf6_plot_schematic(analysis$metadata), source = analysis$metadata),
      "Fig6B_design", 8.4, 4.4, "B")
  classes <- sort(unique(as.character(analysis$class_summary$lipid_class)))
  for (diet in c("SD", "WD")) {
    add(cf6_plot_class_membership(analysis$class_summary, diet, classes),
        paste0("Fig6C_", diet), 7.0, max(5.5, 2.2 + .26 * length(classes)), "C", diet)
    if (!is.null(analysis$class_pair_density)) add(cf6_plot_pair_density(analysis$class_pair_density, diet),
        paste0("Fig6C_pair_density_", diet), 7.2, max(5.5, 1.5 + .26 * length(classes)), "C_pair_density", diet)
  }
  for (sex in c("female", "male")) for (zone in c("pericentral", "periportal")) {
    context <- cf6_plot_context(sex, zone)
    reps <- analysis$representatives[
      cf6_plot_context(analysis$representatives$sex, analysis$representatives$zone) == context, , drop = FALSE]
    gene <- if (nrow(reps)) reps$gene[1] else NULL
    for (diet in c("SD", "WD")) {
      add(cf6_plot_network(analysis$representative_edges, sex, zone, diet, gene),
          paste0("Fig6D_", context, "_", diet), 7.4, 5.4, "D", diet, context)
    }
  }
  e <- cf6_plot_diet_zonation(analysis$zonation_diet_effects, analysis$representatives,
                              analysis$primary_interactions)
  n_genes <- if (nrow(e$source)) length(unique(e$source$gene)) else 1L
  add(e, "Fig6E_diet_zonation", 7.6, max(4.8, 2.3 + .28 * n_genes), "E")
  for (sex in c("female", "male")) for (zone in c("pericentral", "periportal")) {
    context <- cf6_plot_context(sex, zone)
    f <- cf6_plot_examples(analysis$examples, sex, zone)
    n_pairs <- if (nrow(f$source)) length(unique(f$source$example)) else 1L
    add(f, paste0("Fig6F_", context), ifelse(n_pairs == 1L, 6.3, 8.5),
        ifelse(n_pairs == 1L, 4.8, max(4.4, 1.8 + 2.5 * ceiling(n_pairs / 2))), "F", context = context)
  }
  for (kind in unique(analysis$qc_pca$kind)) add(cf6_plot_pca(analysis$qc_pca, kind),
    paste0("FigS6_PCA_", kind), 6.8, 5.4, "supplementary_QC")
  for (contrast in unique(analysis$lipid_results$contrast)) add(cf6_plot_lipid_effect(analysis$lipid_results, contrast),
    paste0("FigS6_lipid_", contrast), 6.8, 5.4, "supplementary_lipid")
  if (nrow(analysis$thresholds)) add(cf6_plot_thresholds(analysis$thresholds),
    "FigS6_display_threshold", 6.8, 4.9, "supplementary_threshold")
  if (!is.null(analysis$quality_sensitivity) && nrow(analysis$quality_sensitivity))
    add(cf6_plot_quality(analysis$quality_sensitivity), "FigS6_lipid_quality", 7.8, 5.0, "supplementary_quality")
  marker_file <- file.path(output, "QA", "provenance", "marker_expression.csv")
  markers <- if (!is.null(analysis$marker_expression)) analysis$marker_expression else
    if (file.exists(marker_file)) utils::read.csv(marker_file, stringsAsFactors = FALSE) else NULL
  if (is.null(markers) && !is.null(analysis$input_qc$marker_by_profile))
    markers <- analysis$input_qc$marker_by_profile
  if (!is.null(markers)) add(cf6_plot_markers(markers), "FigS6_zonation_markers", 7.2, 9.2, "supplementary_QC")
  out <- do.call(rbind, manifest)
  rownames(out) <- NULL
  utils::write.csv(out, file.path(output, "panel_manifest.csv"), row.names = FALSE, na = "")
  out
}
