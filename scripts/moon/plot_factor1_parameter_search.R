#!/usr/bin/env Rscript
# Comparison figures only; no network scores or grid artifacts are modified.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/moon/plot_factor1_parameter_search.R SUMMARY_DIRECTORY")
output <- args[[1L]]
grid <- utils::read.csv(file.path(output, "grid_assessment.csv"), check.names = FALSE, stringsAsFactors = FALSE)
shortlist <- utils::read.csv(file.path(output, "configuration_shortlist.csv"), check.names = FALSE, stringsAsFactors = FALSE)
stability <- utils::read.csv(file.path(output, "eligibility_stability_counts.csv"), check.names = FALSE, stringsAsFactors = FALSE)
truth <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")
grid$eligible <- truth(grid$eligible)
grid$interpretable_eligible <- truth(grid$interpretable_eligible)
if (nrow(shortlist)) shortlist$interpretable_eligible <- truth(shortlist$interpretable_eligible)
figures <- file.path(output, "figures")
dir.create(figures, showWarnings = FALSE)

draw_shortlist <- function() {
  graphics::par(mar = c(7, 13, 4, 2), family = "sans")
  if (!nrow(shortlist)) {
    graphics::plot.new(); graphics::text(0.5, 0.5, "No networks meet the declared eligibility criteria.")
    return(invisible(NULL))
  }
  config <- paste(shortlist$compression, shortlist$n_layers, shortlist$topology_radius, sep = " / ")
  # Prefer the explicitly interpretable shortlist within each configuration for
  # this panel; retain configurations without an interpretable option in gray.
  pick <- unlist(lapply(split(seq_len(nrow(shortlist)), config), function(ids) {
    readable <- ids[shortlist$interpretable_eligible[ids]]
    if (length(readable)) readable[[1L]] else ids[[1L]]
  }), use.names = FALSE)
  selected <- shortlist[pick, , drop = FALSE]
  selected <- selected[order(selected$n_layers, selected$compression, selected$topology_radius), , drop = FALSE]
  if (nrow(selected) > 24L) selected <- selected[seq_len(24L), , drop = FALSE]
  labels <- paste0(selected$compression, " | depth ", selected$n_layers, " | radius ", selected$topology_radius)
  heights <- selected$secondary_threshold
  span <- max(heights, 1, na.rm = TRUE)
  y <- graphics::barplot(heights, horiz = TRUE, names.arg = labels, las = 1,
                         col = ifelse(selected$interpretable_eligible, "#286991", "#A8A8A8"), border = NA,
                         xlim = c(0, span * 1.90), xlab = "Secondary cutoff", cex.names = min(0.9, 16 / nrow(selected)),
                         main = "Configuration-specific threshold-priority candidates")
  values <- paste0("primary ", selected$primary_threshold, " | ", selected$upstream_tfs, " TF / ", selected$downstream_nodes, " lipid | ", selected$edges, " edges")
  graphics::text(heights + span * 0.025, y, values, adj = c(0, 0.5), cex = min(0.75, 15 / nrow(selected)))
  graphics::mtext("Blue: >=2 TFs and >=2 lipids. Gray: only the base criteria pass. Depth/compression configurations are separate scoring problems.", side = 1, line = 5.2, cex = 0.75)
  if (length(pick) > 24L) graphics::mtext("First 24 configurations in depth/compression/radius order; full table is configuration_shortlist.csv.", side = 3, line = 0.4, cex = 0.75)
}

draw_stability <- function() {
  graphics::par(mar = c(8, 5, 4, 3), family = "sans")
  if (!nrow(stability)) {
    graphics::plot.new(); graphics::text(0.5, 0.5, "No finite cutoff pairs are present.")
    return(invisible(NULL))
  }
  primary <- sort(unique(stability$primary_threshold))
  secondary <- sort(unique(stability$secondary_threshold))
  values <- matrix(NA_real_, length(primary), length(secondary))
  numerator <- denominator <- values
  for (i in seq_len(nrow(stability))) {
    p <- match(stability$primary_threshold[[i]], primary)
    s <- match(stability$secondary_threshold[[i]], secondary)
    values[p, s] <- stability$interpretable_fraction_of_tested[[i]]
    numerator[p, s] <- stability$interpretable_configurations[[i]]
    denominator[p, s] <- stability$tested_configurations[[i]]
  }
  # Integer coordinates keep irregular coarse/refined cutoff spacing explicit
  # in labels, without implying untested cells exist between sampled cutoffs.
  graphics::image(seq_along(primary), seq_along(secondary), values,
                  col = grDevices::colorRampPalette(c("#F0F0F0", "#23668F"))(100), zlim = c(0, 1),
                  axes = FALSE, xlab = "Primary cutoff (sampled values)", ylab = "Secondary cutoff (sampled values)",
                  main = "Eligibility across tested configurations at matched cutoffs")
  px <- unique(round(seq(1, length(primary), length.out = min(12, length(primary)))))
  sy <- unique(round(seq(1, length(secondary), length.out = min(12, length(secondary)))))
  graphics::axis(1, at = px, labels = primary[px], las = 2)
  graphics::axis(2, at = sy, labels = secondary[sy], las = 1)
  if (length(values) <= 120L) {
    for (p in seq_along(primary)) for (s in seq_along(secondary)) {
      if (is.finite(values[p, s])) graphics::text(p, s, paste0(numerator[p, s], "/", denominator[p, s]), cex = 0.7, col = ifelse(values[p, s] >= 0.55, "white", "#333333"))
    }
  }
  graphics::mtext("Color: fraction of tested configurations retaining >=2 TFs and >=2 lipids (light=0, dark=1). White cells were not tested.", side = 1, line = 5.0, cex = 0.75)
  graphics::mtext("Counts describe eligibility, not agreement in edges, endpoints or mechanisms. Coarse/refined coverage changes denominators.", side = 1, line = 6.2, cex = 0.75)
}

draw_coverage <- function() {
  graphics::par(mar = c(8, 5, 4, 2), family = "sans")
  eligible <- grid[grid$eligible, , drop = FALSE]
  if (!nrow(eligible)) {
    graphics::plot.new(); graphics::text(0.5, 0.5, "No eligible network sizes to compare.")
    return(invisible(NULL))
  }
  compressions <- sort(unique(eligible$compression))
  palette <- grDevices::hcl.colors(length(compressions), palette = "Dark 3")
  colors <- palette[match(eligible$compression, compressions)]
  size <- 0.6 + 1.3 * eligible$secondary_threshold / max(c(eligible$secondary_threshold, 0.01))
  graphics::plot(eligible$edges, eligible$downstream_nodes, pch = 21L,
                  bg = grDevices::adjustcolor(colors, alpha.f = 0.5), col = colors, cex = size,
                  xlab = "Edges in the extracted network", ylab = "Measured lipid endpoints",
                  main = "Eligible network size and endpoint coverage")
  graphics::legend("topright", legend = compressions, pt.bg = palette, col = palette, pch = 21L, bty = "n", title = "Compression")
  graphics::mtext("Larger circles indicate higher secondary cutoffs. Each point is one configuration/cutoff row; overplotting is possible.", side = 1, line = 5.0, cex = 0.75)
  graphics::mtext("Within-configuration Pareto tables retain cutoff stringency, TF coverage, lipid coverage and edge count as separate objectives.", side = 1, line = 6.2, cex = 0.75)
}

drawers <- list(shortlist = draw_shortlist, eligibility = draw_stability, coverage = draw_coverage)
grDevices::pdf(file.path(figures, "parameter_search_comparison.pdf"), width = 14, height = 10, useDingbats = FALSE)
for (draw in drawers) draw()
grDevices::dev.off()
for (name in names(drawers)) {
  grDevices::png(file.path(figures, paste0("parameter_search_", name, ".png")), width = 2100, height = 1500, res = 150)
  drawers[[name]]()
  grDevices::dev.off()
}
cat("Wrote parameter-search comparison figures in", figures, "\n")
