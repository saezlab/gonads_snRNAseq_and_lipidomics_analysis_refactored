#!/usr/bin/env Rscript
# Base-R publication/export plots. Exact IDs remain in the joined data tables.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript scripts/moon/plot_factor1_network.R RESULTS EVIDENCE_STEM PATH_STEM")
}
output <- args[[1L]]
evidence_stem <- args[[2L]]
path_stem <- args[[3L]]
extraction_label <- if (grepl("solution_catalyst_bounded6", path_stem, fixed = TRUE)) "Exploratory: catalyst routes, <=6 edges" else if (grepl("relaxed", path_stem, fixed = TRUE)) "relaxed exploratory extraction" else if (grepl("intermediate", path_stem, fixed = TRUE)) "intermediate exploratory extraction" else if (grepl("strict", path_stem, fixed = TRUE)) "strict extraction" else if (grepl("solution_p", path_stem, fixed = TRUE)) "parameter-search extraction" else "primary extraction"
if (grepl("uncompressed_sensitivity", output, fixed = TRUE)) extraction_label <- paste0("uncompressed sensitivity, ", if (grepl("strict", path_stem, fixed = TRUE)) "strict thresholds" else if (grepl("relaxed", path_stem, fixed = TRUE)) "relaxed thresholds" else "primary thresholds")
audit_path <- file.path(output, "extraction_audit.json")
if (file.exists(audit_path)) {
  audit <- jsonlite::read_json(audit_path, simplifyVector = TRUE)
  scoring_depth <- sub(".*depth([0-9]+).*", "\\1", audit$run_id)
  topology_radius <- sub("radius([0-9]+).*", "\\1", audit$run_id)
  extraction_label <- paste0("compressed, depth ", scoring_depth, ", radius ", topology_radius,
                            "; cutoffs ", audit$primary_thresh, " / ", audit$secondary_thresh)
}
catalyst_only <- grepl("catalyst", path_stem, fixed = TRUE)
nodes <- utils::read.csv(file.path(output, paste0(evidence_stem, "_nodes.csv")), check.names = FALSE, stringsAsFactors = FALSE)
edges <- utils::read.csv(file.path(output, paste0(evidence_stem, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
paths <- utils::read.csv(file.path(output, paste0(path_stem, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
illustrated <- utils::read.csv(file.path(output, paste0(path_stem, "_illustrated.csv")), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(nodes$node_id))
figures <- file.path(output, "figures")
dir.create(figures, showWarnings = FALSE)

node_color <- function(score) ifelse(score > 0, "#B5403E", ifelse(score < 0, "#3563A7", "#737373"))
short_label <- function(id) {
  row <- nodes[match(id, nodes$node_id), , drop = FALSE]
  if (row$node_type == "enzyme_reaction_instance") {
    catalyst <- gsub("Gene__", "", strsplit(row$catalysts, "|", fixed = TRUE)[[1L]], fixed = TRUE)
    catalyst_label <- paste(utils::head(catalyst, 3L), collapse = "/")
    if (length(catalyst) > 3L) catalyst_label <- paste0(catalyst_label, "/+", length(catalyst) - 3L)
    reaction_label <- if (grepl("LIPID_RULE:", row$reaction_id, fixed = TRUE)) "species-rule reaction" else row$reaction_id
    return(paste(catalyst_label, reaction_label, sep = "\n"))
  }
  paste(strwrap(row$label, width = 23L), collapse = "\n")
}
text_utf8_safe <- function(x) gsub("\u2192", "->", x, fixed = TRUE)

draw_paths <- function() {
  graphics::par(mar = c(1, 1, 1, 1), family = "sans")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  graphics::text(0.015, 0.977, paste0("Factor 1 MOON solution: ", extraction_label), adj = c(0, 1), cex = 1.40, font = 2)
  graphics::text(0.015, 0.934, "Arrows follow PKN direction. Colors show score sign; dashed diamonds mark measured-species reaction hypotheses.", adj = c(0, 1), cex = 0.85, col = "#424242")
  if (catalyst_only) graphics::text(0.015, 0.906, "Route audit: protein regulatory intermediates, then one terminal catalyst/reaction/product step; at most six edges.", adj = c(0, 1), cex = 0.85, col = "#424242")
  if (!nrow(illustrated)) {
    graphics::text(0.5, 0.5, "No upstream TF-to-measured-lipid paths pass this extraction.")
    return(invisible(NULL))
  }
  y_positions <- seq(0.79, 0.26, length.out = max(3L, nrow(illustrated)))
  for (i in seq_len(nrow(illustrated))) {
    path <- illustrated[i, , drop = FALSE]
    ids <- strsplit(path$path_node_ids, "|", fixed = TRUE)[[1L]]
    x <- seq(0.06, 0.94, length.out = length(ids))
    y <- y_positions[[i]]
    graphics::text(0.025, y + 0.09,
                   paste0(i, ". ", path$upstream_label, " to ", path$downstream_label,
                          "   |   ", path$edge_count, " edges; minimum |inferred score| = ", formatC(path$minimum_absolute_inferred_score, digits = 2L, format = "f")),
                   adj = c(0, 0.5), cex = 0.9, font = 2)
    for (j in seq_len(length(ids) - 1L)) {
      edge <- edges[edges$source == ids[[j]] & edges$target == ids[[j + 1L]], , drop = FALSE]
      stopifnot(nrow(edge) == 1L)
      graphics::arrows(x[[j]] + 0.026, y, x[[j + 1L]] - 0.026, y, length = 0.065, lwd = 1.3, col = "#555555")
      graphics::text(mean(x[c(j, j + 1L)]), y + 0.014, ifelse(edge$interaction > 0, "+", "-"), cex = 1, col = "#555555")
    }
    for (j in seq_along(ids)) {
      row <- nodes[match(ids[[j]], nodes$node_id), , drop = FALSE]
      color <- node_color(row$score)
      if (row$node_type == "enzyme_reaction_instance") {
        graphics::polygon(x[[j]] + c(-0.024, 0, 0.024, 0), y + c(0, 0.02, 0, -0.02),
                          col = "white", border = color, lwd = 2,
                          lty = ifelse(grepl("LIPID_RULE:", row$reaction_id, fixed = TRUE), 2L, 1L))
      } else if (row$level == 0) {
        graphics::polygon(x[[j]] + c(-0.022, -0.011, 0.011, 0.022, 0.011, -0.011),
                          y + c(0, 0.018, 0.018, 0, -0.018, -0.018), col = color, border = color)
      } else if (row$type == "upstream_input") {
        graphics::points(x[[j]], y, pch = 21L, bg = color, col = color, cex = 2.7)
      } else {
        graphics::rect(x[[j]] - 0.019, y - 0.014, x[[j]] + 0.019, y + 0.014, col = color, border = color)
      }
      graphics::text(x[[j]], y - 0.031, short_label(ids[[j]]), adj = c(0.5, 1), cex = max(0.53, 0.85 - 0.035 * (length(ids) - 4L)))
      value_label <- if (row$level == 0) paste0("W = ", formatC(row$score, digits = 3L, format = "g")) else paste0("MOON = ", formatC(row$score, digits = 2L, format = "f"))
      graphics::text(x[[j]], y + 0.032, value_label, cex = 0.68, col = color)
    }
  }
  graphics::text(0.015, 0.096, "Selection: up to three distinct observed lipids with largest |Factor 1 weight|; strongest bottleneck path per endpoint, then shortest path.", adj = c(0, 1), cex = 0.77)
  graphics::text(0.015, 0.065, "Each TF/lipid pair uses a shortest route; ties maximize the minimum inferred score. All pairs and exact IDs are exported separately.", adj = c(0, 1), cex = 0.77)
  graphics::text(0.015, 0.034, "Hypotheses only: class-to-species substrate specificity is assumed. Positive reaction traversal does not establish flux or metabolite production.", adj = c(0, 1), cex = 0.77, col = "#444444")
}

draw_reachability <- function() {
  if (!nrow(paths)) {
    graphics::plot.new(); graphics::text(0.5, 0.5, "No TF-to-observation paths.")
    return(invisible(NULL))
  }
  tf_counts <- sort(table(paths$upstream_tf), decreasing = TRUE)
  tfs <- names(utils::head(tf_counts, 30L))
  targets <- unique(paths$downstream_node[order(-abs(paths$downstream_weight), paths$downstream_node)])
  targets <- utils::head(targets, 25L)
  values <- matrix(NA_real_, nrow = length(targets), ncol = length(tfs), dimnames = list(targets, tfs))
  for (i in seq_len(nrow(paths))) {
    row <- paths[i, , drop = FALSE]
    if (row$upstream_tf %in% tfs && row$downstream_node %in% targets) values[row$downstream_node, row$upstream_tf] <- row$minimum_absolute_inferred_score
  }
  graphics::par(mar = c(11, 15, 5, 3), family = "sans")
  limit <- max(values, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) limit <- 1
  graphics::image(seq_along(tfs), seq_along(targets), t(values), xlab = "", ylab = "", axes = FALSE,
                  col = grDevices::colorRampPalette(c("#E2EDF4", "#185A7D"))(64L), zlim = c(0, limit),
                  main = paste0("TF-to-lipid reachability: ", extraction_label))
  graphics::axis(1, seq_along(tfs), gsub("Gene__", "", tfs, fixed = TRUE), las = 2, cex.axis = 0.75, tick = FALSE)
  graphics::axis(2, seq_along(targets), vapply(targets, short_label, character(1)), las = 2, cex.axis = 0.72, tick = FALSE)
  for (x in seq_along(tfs)) for (y in seq_along(targets)) {
    if (is.finite(values[y, x])) graphics::text(x, y, formatC(values[y, x], format = "f", digits = 1L), cex = 0.57, col = ifelse(values[y, x] > 0.55 * limit, "white", "#222222"))
  }
  graphics::mtext(paste0("Cell: minimum |inferred MOON score| along chosen shortest path. Blank: no retained path.\nShowing ", length(tfs), "/", length(tf_counts), " TFs (most reachable endpoints) and ", length(targets), "/", length(unique(paths$downstream_node)), " lipids (largest |W|)."), side = 3, line = 0.5, cex = 0.75)
}

pdf_file <- file.path(figures, paste0(path_stem, "_network.pdf"))
grDevices::pdf(pdf_file, width = 14, height = 10, useDingbats = FALSE)
draw_paths()
draw_reachability()
grDevices::dev.off()
for (kind in c("illustrative", "reachability")) {
  grDevices::png(file.path(figures, paste0(path_stem, "_", kind, ".png")), width = 2100, height = 1500, res = 150)
  if (kind == "illustrative") draw_paths() else draw_reachability()
  grDevices::dev.off()
}
cat("Wrote figures:", pdf_file, "\n")
