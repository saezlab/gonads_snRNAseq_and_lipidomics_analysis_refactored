# Figure 3F–G uses matched lipids and the original per-lipid limma statistics.
fig3fg_genes <- c(F = "Cyp2d9", G = "Cyp2c38")

fig3fg_prepare <- function(t_table, network) {
  if (!all(c("ID", tolower(fig3fg_genes)) %in% names(t_table)) ||
      anyNA(t_table$ID) || any(!nzchar(t_table$ID)) || anyDuplicated(t_table$ID)) {
    stop("Figure 3F–G requires unique lipid IDs and both knockdown t-statistic columns.")
  }
  if (!all(c("source", "target", "mor") %in% names(network))) {
    stop("Figure 3F–G requires network source, target and mor columns.")
  }

  panels <- lapply(names(fig3fg_genes), function(panel) {
    gene <- fig3fg_genes[[panel]]
    condition <- tolower(gene)
    edges <- network[which(tolower(network$source) == condition &
                           grepl("[(]", network$target)), c("source", "target", "mor")]
    if (anyDuplicated(edges$target)) stop("Multiple reference weights per lipid for ", gene)
    edge_index <- match(t_table$ID, edges$target)
    matched <- !is.na(edge_index)
    lipid <- t_table$ID[matched]
    r <- edges$mor[edge_index[matched]]
    knockdown_t <- t_table[[condition]][matched]
    if (length(lipid) < 3L || !is.numeric(r) || !is.numeric(knockdown_t) ||
        any(!is.finite(r)) || any(abs(r) > 1) || any(!is.finite(knockdown_t)) ||
        length(unique(r)) < 2L || length(unique(knockdown_t)) < 2L) {
      stop("Figure 3", panel, " requires at least three matched lipids and finite, varying t and r values.")
    }

    mat <- matrix(knockdown_t, ncol = 1L, dimnames = list(lipid, condition))
    ulm <- decoupleR::run_ulm(
      mat, edges[edge_index[matched], , drop = FALSE],
      .source = "source", .target = "target", .mor = "mor", minsize = 2L
    )
    stopifnot(nrow(ulm) == 1L, is.finite(ulm$score), is.finite(ulm$p_value))
    display_t <- pmax(-10, pmin(10, knockdown_t))
    list(
      points = data.frame(panel, gene, lipid, pearson_r = r, knockdown_t, display_t),
      summary = data.frame(panel, gene, n_matched_lipids = length(lipid),
                           n_capped_for_display = sum(display_t != knockdown_t),
                           ULM_score = ulm$score, nominal_ULM_p = ulm$p_value)
    )
  })
  stats::setNames(panels, names(fig3fg_genes))
}

fig3fg_write <- function(panels, output) {
  source("scripts/config/plotting_themes.R", local = TRUE)
  gonadal_colors <- get_color_palette("gonadal")
  colors <- c(F = gonadal_colors[["Testes"]], G = gonadal_colors[["Ovaries"]])
  labelled_lipids <- c("PS(16:0_20:3)", "PI(36:3)", "PI(38:3)",
                       "SM(18:1;O2/16:0)", "HexCer(18:0;O2/24:1)", "SM(18:0;O2/16:0)")
  note <- "Fit and ULM use the original t-statistics. Points beyond ±10 are displayed at ±10."

  plots <- lapply(names(panels), function(panel) {
    points <- panels[[panel]]$points
    scores <- panels[[panel]]$summary
    points$label <- ""
    if (panel == "G") {
      label <- points$lipid %in% labelled_lipids & points$pearson_r * points$knockdown_t < 0
      points$label[label] <- points$lipid[label]
    }
    color <- colors[[panel]]
    gene_symbol <- as.name(fig3fg_genes[[panel]])
    ggplot2::ggplot(points, ggplot2::aes(x = pearson_r, y = knockdown_t)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey84", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "grey84", linewidth = 0.35) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, level = 0.95,
                           color = "grey35", fill = "grey65", alpha = 0.38, linewidth = 0.65) +
      ggplot2::geom_point(ggplot2::aes(y = display_t), color = color, alpha = 0.76, size = 1.8) +
      ggrepel::geom_text_repel(
        ggplot2::aes(y = display_t, label = label), color = color, seed = 1124,
        size = 3.1, box.padding = 0.35, point.padding = 0.20,
        max.overlaps = Inf, max.iter = 12000L, max.time = 2,
        segment.color = color, segment.size = 0.28, min.segment.length = 0
      ) +
      ggplot2::scale_x_continuous(breaks = seq(-1, 1, by = 0.5)) +
      ggplot2::scale_y_continuous(breaks = seq(-10, 10, by = 5)) +
      ggplot2::coord_cartesian(xlim = c(-1.04, 1.04), ylim = c(-10.7, 10.7), expand = FALSE) +
      ggplot2::labs(
        tag = panel, title = bquote(italic(.(gene_symbol))~"knockdown"),
        subtitle = sprintf("%d matched lipids  |  ULM %.2f", scores$n_matched_lipids, scores$ULM_score),
        x = "In vivo gene-lipid association (Pearson r)",
        y = "Knockdown - control (limma t)"
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, color = "#202735", margin = ggplot2::margin(b = 4)),
        plot.subtitle = ggplot2::element_text(size = 10.5, color = "#49566a", margin = ggplot2::margin(b = 10)),
        plot.tag = ggplot2::element_text(size = 17, face = "bold"),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 7)),
        plot.margin = ggplot2::margin(10, 10, 8, 8)
      )
  })
  names(plots) <- names(panels)
  for (panel in names(plots)) {
    filename <- paste0(tolower(fig3fg_genes[[panel]]), "_validation_scatter.pdf")
    ggplot2::ggsave(
      file.path(output, filename),
      plots[[panel]] + ggplot2::labs(caption = note) +
        ggplot2::theme(plot.caption = ggplot2::element_text(size = 8.5, hjust = 0)),
      width = 5.65, height = 5.35, device = "pdf", useDingbats = FALSE, bg = "white"
    )
  }
  combined <- patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(
      caption = note,
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 8.5, hjust = 0))
    )
  ggplot2::ggsave(file.path(output, "Fig3F-G_matched.pdf"), combined,
                 width = 11.3, height = 5.35, device = "pdf", useDingbats = FALSE, bg = "white")
  utils::write.csv(do.call(rbind, lapply(panels, "[[", "points")),
                   file.path(output, "Fig3FG_matched_lipids.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, lapply(panels, "[[", "summary")),
                   file.path(output, "Fig3FG_matched_ULM.csv"), row.names = FALSE)
  invisible(output)
}
