#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript plot_factor1_candidate_audit.R AUDIT_DIRECTORY")
output <- args[[1L]]
summary <- utils::read.csv(file.path(output, "candidate_summary.csv"), stringsAsFactors = FALSE)
scores <- utils::read.csv(file.path(output, "candidate_node_scores.csv"), stringsAsFactors = FALSE)
endpoints <- utils::read.csv(file.path(output, "endpoint_identities.csv"), stringsAsFactors = FALSE)
overlap <- utils::read.csv(file.path(output, "pairwise_overlap.csv"), stringsAsFactors = FALSE)
candidates <- summary$candidate
labels <- paste0(candidates, ": ", summary$compression, ", depth ", summary$n_layers)
figures <- file.path(output, "figures")
dir.create(figures, showWarnings = FALSE)

draw_scores <- function() {
  graphics::par(mar = c(7, 5, 5, 2), family = "sans")
  groups <- lapply(candidates, function(candidate) abs(scores$score[scores$candidate == candidate & scores$score_type == "inferred_MOON"]))
  graphics::boxplot(groups, names = labels, col = c("#93B8CD", "#93B8CD", "#D9B67A"),
                    ylab = "Absolute inferred MOON score", main = "Refined candidates: inferred scores and extraction cutoffs", ylim = c(0, max(unlist(groups)) * 1.05))
  for (i in seq_along(candidates)) {
    graphics::segments(i - 0.3, summary$secondary_threshold[[i]], i + 0.3, summary$secondary_threshold[[i]], col = "#A33333", lwd = 2, lty = 2)
    graphics::segments(i - 0.3, summary$primary_threshold[[i]], i + 0.3, summary$primary_threshold[[i]], col = "#333333", lwd = 2, lty = 3)
  }
  graphics::legend("topright", legend = c("Secondary cutoff", "Primary cutoff"), col = c("#A33333", "#333333"), lty = c(2, 3), lwd = 2, bty = "n")
  graphics::mtext("Level-0 lipid weights are excluded. These are distributions within selected networks, not comparable effect-size or significance distributions.", side = 1, line = 4.3, cex = 0.75)
  graphics::mtext("U4/U5 share identical scores on common inferred nodes; U5 adds nodes and also uses a different primary threshold.", side = 1, line = 5.5, cex = 0.75)
}

draw_overlap <- function() {
  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 4, 2), oma = c(5.2, 0, 3, 0), family = "sans")
  for (metric in c("nodes", "edges", "upstream")) {
    values <- diag(1, length(candidates))
    for (i in seq_len(nrow(overlap))) {
      left <- match(overlap$left[[i]], candidates); right <- match(overlap$right[[i]], candidates)
      values[left, right] <- values[right, left] <- overlap[[paste0(metric, "_jaccard")]][[i]]
    }
    graphics::image(seq_along(candidates), seq_along(candidates), values, zlim = c(0, 1),
                    col = grDevices::colorRampPalette(c("#F0F0F0", "#285E7E"))(100), axes = FALSE, xlab = "", ylab = "",
                    main = switch(metric, nodes = "Node overlap", edges = "Signed edge overlap", upstream = "TF endpoint overlap"))
    graphics::axis(1, seq_along(candidates), candidates)
    graphics::axis(2, seq_along(candidates), candidates, las = 1)
    for (x in seq_along(candidates)) for (y in seq_along(candidates)) graphics::text(x, y, formatC(values[x, y], digits = 2, format = "f"), col = ifelse(values[x, y] > 0.5, "white", "#333333"))
  }
  graphics::mtext("Refined candidate overlap: exact identities and signs", outer = TRUE, side = 3, line = 1, cex = 1.2, font = 2)
  graphics::mtext("Jaccard = intersection / union. Both lipid endpoints are shared by all three candidates; no upstream TF endpoint is shared by all three.", outer = TRUE, side = 1, line = 1.5, cex = 0.8)
  graphics::mtext("U4: uncompressed depth 4, primary/secondary 1.90/1.60. U5: uncompressed depth 5, 2.05/1.60. C4: compressed depth 4, 2.00/0.85.", outer = TRUE, side = 1, line = 3.1, cex = 0.75)
  graphics::par(mfrow = c(1, 1))
}

draw_endpoints <- function() {
  graphics::par(mar = c(7, 8, 5, 2), family = "sans")
  tf <- endpoints[endpoints$role == "upstream_tf", , drop = FALSE]
  names <- sort(unique(tf$node_id))
  values <- matrix(0, length(candidates), length(names))
  labels_score <- matrix("", length(candidates), length(names))
  for (i in seq_len(nrow(tf))) {
    x <- match(tf$candidate[[i]], candidates); y <- match(tf$node_id[[i]], names)
    values[x, y] <- sign(tf$inferred_moon_score[[i]])
    labels_score[x, y] <- formatC(tf$inferred_moon_score[[i]], digits = 2, format = "f")
  }
  graphics::image(seq_along(candidates), seq_along(names), values, zlim = c(-1, 1), col = c("#4477AA", "#F0F0F0", "#AD4848"),
                  axes = FALSE, xlab = "", ylab = "", main = "Upstream TF endpoints retained by each candidate")
  graphics::axis(1, seq_along(candidates), labels)
  graphics::axis(2, seq_along(names), sub("Gene__", "", names, fixed = TRUE), las = 1)
  for (x in seq_along(candidates)) for (y in seq_along(names)) graphics::text(x, y, ifelse(values[x, y] == 0, "not retained", labels_score[x, y]), col = ifelse(values[x, y] == 0, "#666666", "white"))
  graphics::mtext("Numbers are inferred MOON scores. Gray means absent as an upstream endpoint in the selected network, not necessarily unscored.", side = 1, line = 4.2, cex = 0.8)
  graphics::mtext("Shared downstream observations: DG(24:0_18:2) and DG(24:1_18:2). Their species-specific Smpd1 reaction rules remain hypotheses.", side = 1, line = 5.4, cex = 0.8)
}

drawers <- list(score_distributions = draw_scores, overlap = draw_overlap, tf_endpoints = draw_endpoints)
grDevices::pdf(file.path(figures, "refined_candidate_comparison.pdf"), width = 14, height = 9, useDingbats = FALSE)
for (draw in drawers) draw()
grDevices::dev.off()
for (name in names(drawers)) {
  grDevices::png(file.path(figures, paste0("refined_candidate_", name, ".png")), width = 2100, height = 1350, res = 150)
  drawers[[name]]()
  grDevices::dev.off()
}
cat("Wrote refined-candidate comparison figures in", figures, "\n")
